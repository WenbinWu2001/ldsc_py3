"""Workflow-layer LD-score orchestration and normalized result objects.

This module is the public boundary for LD-score path handling and result
normalization. Callers may provide exact paths, glob patterns, or explicit
chromosome-suite tokens using ``@``. The workflow resolves those tokens into
concrete per-chromosome files and asks ``ref_panel.prepare_chromosome`` for
aligned reference data, annotations, window bounds, and an owned reader. The
numerical kernel consumes that prepared state without reopening the panel.

For ordinary unpartitioned LD-score runs, callers may omit both baseline and
query annotations. In that case the workflow constructs an all-ones baseline
annotation named exactly ``base`` over the retained reference-panel metadata and
continues through the same calculator and output writer used by partitioned
runs. Query annotations remain partitioned-LDSC inputs: they are accepted only
when explicit baseline annotations are supplied.

Parsed workflow entry points write ``diagnostics/ldscore.log`` under
``output_dir`` after preflighting the complete LD-score artifact family. Direct
``LDScoreCalculator.run(...)`` calls remain data-oriented and do not create log
files.
"""

from __future__ import annotations

import argparse
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import dataclass, field, replace as dataclass_replace
import logging
import multiprocessing as mp
import os
from pathlib import Path
import sys
from typing import Any, Sequence
import warnings

import numpy as np
import pandas as pd

from ._chr_sampler import sample_frame_from_chr_pattern
from ._kernel.snp_identity import empty_identity_drop_frame, identity_mode_family, is_allele_aware_mode
from .column_inference import normalize_genome_build, normalize_snp_identifier_mode
from .config import (
    GlobalConfig,
    AnnotationBuildConfig,
    LDScoreConfig,
    RefPanelConfig,
    get_global_config,
    print_global_config_banner,
    validate_config_compatibility,
)
from .hm3 import packaged_hm3_curated_map_path
from .genome_build_inference import resolve_genome_build
from .outputs import (
    LDScoreDirectoryWriter,
    LDScoreOutputConfig,
    REGRESSION_LD_SCORE_COLUMN,
)
from .path_resolution import (
    ensure_output_directory,
    normalize_optional_path_token,
    normalize_path_token,
    resolve_scalar_path,
    split_cli_path_tokens,
)
from ._logging import log_inputs, log_outputs, materializing_overwrite_guard, workflow_logging
from ._kernel import ldscore as kernel_ldscore
from ._kernel import regions as kernel_regions
from ._kernel.overlap import OverlapContribution, sum_overlap_contributions
from ._kernel.regions import REGR_SNPS_EXCLUDE_REGIONS_CHOICES, regr_snps_exclude_regions_choice_to_presets
from ._kernel.identifiers import build_snp_id_series, read_snp_restriction_keys
from ._kernel.snp_identity import RestrictionIdentityKeys
from ._row_alignment import assert_same_snp_rows
from .errors import LDSCConfigError, LDSCInputError, LDSCInternalError, LDSCUsageError
from .query_annotations import (
    QueryAnnotationStatus, finalize_query_statuses, gene_query_statuses, gene_viability_errors,
    _log_gene_list_rejections, _log_gene_list_snp_support,
    _log_query_annotation_statuses, _all_query_annotations_skipped_message,
)
from .annotation_semantics import (
    require_unique_annotation_names,
)


LOGGER = logging.getLogger("LDSC.ldscore_calculator")
MAX_CONSOLE_GENE_ISSUES = 10


_LDSCORE_SUFFIX_COLUMNS = ("CHR", "SNP", "POS", "BP", "CM", "MAF")
_QUERY_REQUIRES_BASELINE_MESSAGE = (
    "ldscore cannot run query annotations without baseline annotations. Most likely "
    "a prebuilt, BED, or gene-list query source was supplied without `--baseline-annot-sources`. "
    "Pass matching baseline annotations, or create an explicit all-ones `base` "
    "baseline annotation over the query annotation universe before running "
    "partitioned LDSC."
)
_LDSCORE_INTERSECTION_DOC = (
    "docs/troubleshooting.md#ldscore-no-annotation-snps-remain-after-reference-panel-intersection"
)
_LDSCORE_PARQUET_DOC = "docs/troubleshooting.md#ldscore-parquet-r2-input-is-incompatible"
_LDSCORE_BUILD_DOC = "docs/troubleshooting.md#ldscore-genome-build-could-not-be-resolved-consistently"




@dataclass(frozen=True)
class ChromLDScoreResult:
    """Normalized per-chromosome LD-score output in split baseline/query form.

    Parameters
    ----------
    chrom : str
        Chromosome label.
    baseline_table : pandas.DataFrame
        Table with ``CHR``, ``SNP``, ``POS``, ``regression_ld_scores``, and baseline
        LD-score columns.
    query_table : pandas.DataFrame or None
        Optional table with ``CHR``, ``SNP``, ``POS``, and query LD-score
        columns.
    count_records : list of dict
        Manifest-ready count records keyed by annotation column.
    count_config : dict, optional
        Manifest-ready metadata describing how optional common-SNP counts were
        computed, including the inclusive common-MAF threshold.
    baseline_columns, query_columns : list of str
        Ordered annotation LD-score columns in the baseline and query tables.
    ld_reference_snps, ld_regression_snps : frozenset of str
        SNP universes used for count records and regression rows.
    annotation_types : dict of str to str, optional
        Per-column ``binary`` or ``quantitative`` classification used only for
        interpretation messages and persisted provenance.
    """
    chrom: str
    baseline_table: pd.DataFrame
    query_table: pd.DataFrame | None
    count_records: list[dict[str, Any]]
    baseline_columns: list[str]
    query_columns: list[str]
    ld_reference_snps: frozenset[str]
    ld_regression_snps: frozenset[str]
    snp_count_totals: dict[str, np.ndarray] = field(default_factory=dict, repr=False)
    count_config: dict[str, Any] = field(default_factory=dict)
    output_paths: dict[str, str] = field(default_factory=dict)
    config_snapshot: GlobalConfig | None = None
    overlap: OverlapContribution | None = field(default=None, repr=False)
    reference_snp_count: int = 0
    regression_selected_snp_count: int = 0
    regression_region_removed_snp_count: int = 0
    identity_drops: pd.DataFrame = field(default_factory=empty_identity_drop_frame, repr=False)
    annotation_types: dict[str, str] = field(default_factory=dict)

    def validate(self) -> None:
        """Check the normalized public contract for chromosome-level results."""
        required = {"CHR", "SNP", "POS", REGRESSION_LD_SCORE_COLUMN, *self.baseline_columns}
        missing = required - set(self.baseline_table.columns)
        if missing:
            raise LDSCInternalError(
                "LD-score result validation failed in ChromLDScoreResult.validate(): "
                f"baseline_table is missing required columns {sorted(missing)}. "
                "Most likely an internal result assembly step dropped metadata or LD-score columns. "
                "Re-run with DEBUG logging and report the traceback."
            )
        if self.query_columns and self.query_table is None:
            raise LDSCInternalError(
                "LD-score result validation failed in ChromLDScoreResult.validate(): "
                "query columns are present but query_table is missing. Most likely an internal "
                "split-table assembly step lost the query table. Re-run with DEBUG logging and "
                "report the traceback."
            )
        if self.query_table is not None:
            missing_query = {"CHR", "SNP", "POS", *self.query_columns} - set(self.query_table.columns)
            if missing_query:
                raise LDSCInternalError(
                    "LD-score result validation failed in ChromLDScoreResult.validate(): "
                    f"query_table is missing required columns {sorted(missing_query)}. "
                    "Most likely an internal query-table assembly step dropped metadata or "
                    "annotation columns. Re-run with DEBUG logging and report the traceback."
                )
            assert_same_snp_rows(
                self.baseline_table,
                self.query_table,
                context="query rows must match baseline rows on CHR/SNP/POS",
                snp_identifier=getattr(self.config_snapshot, "snp_identifier", "chr_pos_allele_aware"),
            )

    def summary(self) -> dict[str, Any]:
        """Return a compact summary of chromosome-level retained rows."""
        return {
            "chrom": self.chrom,
            "n_rows": len(self.baseline_table),
            "count_columns": [record["column"] for record in self.count_records],
        }


@dataclass(frozen=True)
class LDScoreResult:
    """Aggregated cross-chromosome LD-score result in split persisted form.

    Parameters
    ----------
    baseline_table : pandas.DataFrame
        Cross-chromosome table persisted as ``ldscore.baseline.parquet`` when
        the result is written. Required columns are ``CHR``, ``SNP``, ``POS``,
        ``regression_ld_scores``, and every entry in ``baseline_columns``.
    query_table : pandas.DataFrame or None
        Optional cross-chromosome table persisted as ``ldscore.query.parquet``.
        Required columns are ``CHR``, ``SNP``, ``POS``, and every entry in
        ``query_columns``.
    count_records : list of dict
        Manifest-ready count records. Each record names an annotation column and
        its all-SNP and optional common-SNP counts.
    count_config : dict, optional
        Manifest-ready metadata describing how optional common-SNP counts were
        computed, including the inclusive common-MAF threshold.
    baseline_columns, query_columns : list of str
        Ordered LD-score columns available to regression workflows.
    ld_reference_snps, ld_regression_snps : frozenset of str
        SNP universes used for count records and persisted regression rows.
    chromosome_results : list of ChromLDScoreResult
        Per-chromosome components used to build the aggregate result.
    output_paths : dict, optional
        Paths written by ``LDScoreDirectoryWriter`` when output was requested.
    config_snapshot : GlobalConfig or None, optional
        Shared configuration active when the result was computed.
    query_statuses : tuple of QueryAnnotationStatus, optional
        Ordered per-source outcomes for BED and gene-list query annotations.
    gene_list_batch : GeneListBatchResolution or None, optional
        Complete gene-list audit, source summary, policy, and selected genes.
    chromosome_scope : dict, optional
        Validated input chromosome coverage and effective analysis scope.
        Empty for workflows that do not perform query coverage preflight.
    legacy_ldsc2_import : dict or None, optional
        Explicit converter provenance for an LDSC2 LD-score suite. ``None``
        for natively computed LDSC3 results.
    annotation_types : dict of str to str, optional
        Per-column semantic classification used for interpretation only.
    """
    baseline_table: pd.DataFrame
    query_table: pd.DataFrame | None
    count_records: list[dict[str, Any]]
    baseline_columns: list[str]
    query_columns: list[str]
    ld_reference_snps: frozenset[str]
    ld_regression_snps: frozenset[str]
    chromosome_results: list[ChromLDScoreResult]
    output_paths: dict[str, str] = field(default_factory=dict)
    count_config: dict[str, Any] = field(default_factory=dict)
    config_snapshot: GlobalConfig | None = None
    overlap: "LDScoreOverlap | None" = field(default=None, repr=False)
    query_statuses: tuple[QueryAnnotationStatus, ...] = ()
    gene_list_batch: Any | None = field(default=None, repr=False)
    snp_universe_policy: dict[str, Any] | None = None
    index_provenance: dict[str, str] | None = None
    legacy_ldsc2_import: dict[str, Any] | None = None
    identity_drops_by_chrom: dict[str, pd.DataFrame] = field(default_factory=dict, repr=False)
    annotation_types: dict[str, str] = field(default_factory=dict)
    chromosome_scope: dict[str, Any] = field(default_factory=dict)

    def validate(self, *, require_query_alignment: bool = True) -> None:
        """Check the normalized public contract for aggregated results."""
        require_unique_annotation_names(self.baseline_columns, self.query_columns)
        required = {"CHR", "SNP", "POS", REGRESSION_LD_SCORE_COLUMN, *self.baseline_columns}
        missing = required - set(self.baseline_table.columns)
        if missing:
            raise LDSCInternalError(
                "LD-score result validation failed in LDScoreResult.validate(): "
                f"baseline_table is missing required columns {sorted(missing)}. "
                "Most likely chromosome aggregation dropped metadata or LD-score columns. "
                "Re-run with DEBUG logging and report the traceback."
            )
        if self.query_columns and self.query_table is None:
            raise LDSCInternalError(
                "LD-score result validation failed in LDScoreResult.validate(): "
                "query columns are present but query_table is missing. Most likely chromosome "
                "aggregation lost the query table. Re-run with DEBUG logging and report the traceback."
            )
        if not self.query_columns and self.query_table is not None:
            raise LDSCInternalError(
                "LD-score result validation failed in LDScoreResult.validate(): "
                "query_table was provided but query_columns is empty. Most likely an internal "
                "result assembly step preserved an unexpected query table. Re-run with DEBUG "
                "logging and report the traceback."
            )
        if self.query_table is not None:
            missing_query = {"CHR", "SNP", "POS", *self.query_columns} - set(self.query_table.columns)
            if missing_query:
                raise LDSCInternalError(
                    "LD-score result validation failed in LDScoreResult.validate(): "
                    f"query_table is missing required columns {sorted(missing_query)}. "
                    "Most likely chromosome aggregation dropped metadata or query annotation "
                    "columns. Re-run with DEBUG logging and report the traceback."
                )
            if require_query_alignment:
                assert_same_snp_rows(
                    self.baseline_table,
                    self.query_table,
                    context="query rows must match baseline rows on CHR/SNP/POS",
                    snp_identifier=getattr(self.config_snapshot, "snp_identifier", "chr_pos_allele_aware"),
                )

    def summary(self) -> dict[str, Any]:
        """Return a compact cross-chromosome summary."""
        return {
            "n_rows": len(self.baseline_table),
            "chromosomes": [result.chrom for result in self.chromosome_results],
            "count_columns": [record["column"] for record in self.count_records],
        }


class LDScoreCalculator:
    """Orchestrate chromosome-wise LD-score calculation.

    This service assembles annotation and reference-panel inputs, delegates the
    heavy computation to the internal LD-score kernel, aggregates chromosome
    outputs, and optionally hands the result to the output layer. The calculator
    delegates reference filtering, annotation alignment, LD windows, and reader
    policy to ``ref_panel.prepare_chromosome``. Reference-panel CM and MAF are
    authoritative. The prepared reader is closed at the chromosome boundary,
    including when numerical computation fails.
    """

    def __init__(self, output_writer: LDScoreDirectoryWriter | None = None) -> None:
        """Initialize the calculator with the directory writer used for LD-score outputs."""
        self.output_writer = output_writer or LDScoreDirectoryWriter()

    def run(
        self,
        annotation_bundle,
        ref_panel,
        ldscore_config: LDScoreConfig,
        global_config: GlobalConfig,
        output_config: LDScoreOutputConfig | None = None,
        regression_snps: set[str] | RestrictionIdentityKeys | None = None,
        regression_regions: kernel_regions.RegionIntervals | None = None,
        config_snapshot: dict[str, Any] | None = None,
    ) -> LDScoreResult:
        """Compute and aggregate LD scores across all chromosomes.

        Chromosomes are computed independently and then aggregated in input
        order. ``ldscore_config.threads`` controls cross-chromosome parallelism
        (joblib ``n_jobs`` convention): ``1`` (default) runs sequentially
        in-process, while ``-1`` (all cores), ``-2`` (all but one), or any
        positive ``N`` fan out over a spawn ``ProcessPoolExecutor``. The
        aggregated output is identical regardless of this setting.

        Parameters
        ----------
        annotation_bundle : AnnotationBundle
            Aligned SNP-level baseline and query annotations.
        ref_panel : RefPanel
            Reference-panel adapter that prepares the aligned SNP universe,
            annotation matrix, LD windows, and reader. Preparation applies the
            configured SNP/sample/MAF restrictions and supplies reference CM
            and MAF independently of annotation-file metadata.
        ldscore_config : LDScoreConfig
            LD-window and retained-SNP settings.
        global_config : GlobalConfig
            Shared SNP identifier, genome-build, and logging settings.
        output_config : LDScoreOutputConfig or None, optional
            If provided, write the canonical LD-score result directory after
            the aggregate result is built. Existing canonical files are refused
            unless ``output_config.overwrite`` is true.
            Default is ``None``, which keeps the result in memory only.
            ``ldscore_config.export_ref_metadata=True`` requires an output
            configuration because the exported metadata is a file artifact.
        regression_snps : set of str, RestrictionIdentityKeys, or None, optional
            Optional regression SNP universe used to define persisted rows and
            regression-weight contributions. The CLI always provides its
            bundled HM3 default; ``None`` is retained for direct API callers.
        config_snapshot : dict or None, optional
            Optional run metadata forwarded to the output layer. Default is
            ``None``.

        Returns
        -------
        LDScoreResult
            Aggregated cross-chromosome result with aligned metadata and output
            paths if writing was requested.
        """
        if ldscore_config.export_ref_metadata and output_config is None:
            raise LDSCUsageError(
                "LDScoreCalculator.run() requires output_config when export_ref_metadata=True. "
                "Reference metadata is a filesystem artifact, so pass LDScoreOutputConfig(output_dir=...)."
            )
        print_global_config_banner(type(self).__name__, global_config)
        if annotation_bundle.config_snapshot is not None:
            validate_config_compatibility(
                annotation_bundle.config_snapshot,
                global_config,
                context="AnnotationBundle and LDScoreCalculator runtime config",
            )
        chromosomes = _chromosomes_from_bundle(annotation_bundle)
        LOGGER.info(_format_ldscore_start_message(annotation_bundle, len(chromosomes)))
        worker_count = _resolve_worker_count(ldscore_config.threads, len(chromosomes))
        cm_source, maf_source = _reference_metadata_sources(ref_panel)
        LOGGER.info(
            f"Reference-panel CM source: {cm_source}; MAF source: {maf_source} "
            "(annotation CM/MAF are ignored)."
        )
        export_dir = None
        if output_config is not None and ldscore_config.export_ref_metadata and cm_source != "parquet_sidecar":
            export_dir = str(output_config.output_dir)
        outcomes = self._run_chromosomes(
            chromosomes=chromosomes,
            annotation_bundle=annotation_bundle,
            ref_panel=ref_panel,
            ldscore_config=ldscore_config,
            global_config=global_config,
            regression_snps=regression_snps,
            regression_regions=regression_regions,
            worker_count=worker_count,
            export_dir=export_dir,
        )
        chromosome_results: list[ChromLDScoreResult] = []
        for chrom in chromosomes:
            outcome = outcomes[chrom]
            if outcome.skipped:
                warnings.warn(f"Skipping chromosome {chrom}: {outcome.skip_message}", UserWarning, stacklevel=2)
                continue
            chromosome_results.append(outcome.result)
        if not chromosome_results:
            raise LDSCInputError(
                "ldscore could not compute any chromosome results after intersecting annotations "
                "with the reference panel. Most likely the annotation SNP IDs, genome build, "
                "or allele-aware identifier mode do not match the reference panel. Use matching "
                "annotation and reference-panel artifacts, or rerun with the correct "
                "`--snp-identifier` and `--genome-build`. "
                f"Other causes & fixes: {_LDSCORE_INTERSECTION_DOC}"
            )
        result = self._aggregate_chromosome_results(
            chromosome_results,
            global_config=global_config,
            count_config=_count_config_from_ldscore_config(ldscore_config),
        )
        binary = [name for name, kind in result.annotation_types.items() if kind == "binary"]
        quantitative = [name for name, kind in result.annotation_types.items() if kind == "quantitative"]
        LOGGER.info(
            "Annotation classification (advisory only): binary=%s; quantitative=%s. "
            "Classification does not change LD scores or downstream fitted coefficients.",
            ", ".join(binary) or "none",
            ", ".join(quantitative) or "none",
        )
        result = dataclass_replace(
            result,
            gene_list_batch=getattr(annotation_bundle, "gene_list_batch", None),
            chromosome_scope=dict(annotation_bundle.source_summary.get("chromosome_scope", {})),
            snp_universe_policy=_snp_universe_policy(
                ref_panel=ref_panel,
                ldscore_config=ldscore_config,
                regression_regions=regression_regions,
                chromosome_results=chromosome_results,
            ),
        )
        result = finalize_query_statuses(
            result,
            tuple(getattr(annotation_bundle, "query_statuses", ())),
        )
        if output_config is not None:
            errors = gene_viability_errors(result.gene_list_batch, result.query_statuses, result)
            if errors:
                self.output_writer.write_query_diagnostics(result, output_config)
                raise LDSCInputError("; ".join(errors))
            output_paths = self.output_writer.write(result, output_config)
            result = _replace_result_output_paths(result, output_paths)
            LOGGER.info(f"Wrote LD-score result directory to '{output_config.output_dir}'.")
        LOGGER.info(
            f"Computed LD scores for {len(chromosome_results)} chromosomes "
            f"and {len(result.baseline_table)} retained SNP rows."
        )
        return result

    def _run_chromosomes(
        self,
        chromosomes,
        annotation_bundle,
        ref_panel,
        ldscore_config: LDScoreConfig,
        global_config: GlobalConfig,
        regression_snps,
        worker_count: int,
        regression_regions: kernel_regions.RegionIntervals | None = None,
        export_dir: str | None = None,
    ) -> dict[str, _ChromOutcome]:
        """Compute every chromosome, sequentially or via a spawn process pool.

        Returns outcomes keyed by chromosome. ``worker_count == 1`` runs inline
        on this calculator and the passed ``ref_panel`` (the exact pre-parallel
        path); larger counts fan out over a spawn ``ProcessPoolExecutor`` whose
        workers rebuild the panel from ``ref_panel.spec`` because a live reader
        cannot cross the process boundary. Callers re-order by the original
        chromosome list for deterministic aggregation, so completion order does
        not matter.
        """
        if worker_count == 1:
            outcomes: dict[str, _ChromOutcome] = {}
            for chrom in chromosomes:
                chrom_bundle = _slice_annotation_bundle(annotation_bundle, chrom)
                try:
                    compute_kwargs = dict(
                        chrom=chrom,
                        annotation_bundle=chrom_bundle,
                        ref_panel=ref_panel,
                        ldscore_config=ldscore_config,
                        global_config=global_config,
                        regression_snps=regression_snps,
                        export_dir=export_dir,
                    )
                    if regression_regions is not None:
                        compute_kwargs["regression_regions"] = regression_regions
                    result = self.compute_chromosome(**compute_kwargs)
                except (ValueError, LDSCInputError) as exc:
                    if _is_empty_intersection(exc, chrom):
                        outcomes[chrom] = _ChromOutcome(chrom=chrom, result=None, skipped=True, skip_message=str(exc))
                        continue
                    raise
                outcomes[chrom] = _ChromOutcome(chrom=chrom, result=result, skipped=False)
            return outcomes

        ref_panel_spec = ref_panel.spec
        ctx = mp.get_context("spawn")
        outcomes = {}
        with ProcessPoolExecutor(
            max_workers=worker_count,
            mp_context=ctx,
            initializer=_init_worker,
            initargs=(regression_snps, regression_regions, global_config.log_level),
        ) as pool:
            remaining = iter(chromosomes)
            futures = set()
            def submit_next():
                chrom = next(remaining,None)
                if chrom is None:
                    return
                futures.add(pool.submit(_compute_one_chromosome,chrom,
                    _slice_annotation_bundle(annotation_bundle,chrom),ref_panel_spec,
                    ldscore_config,global_config,export_dir))
            for _ in range(worker_count):
                submit_next()
            try:
                while futures:
                    completed, _ = wait(futures,return_when=FIRST_COMPLETED)
                    for future in completed:
                        futures.remove(future)
                        outcome = future.result()
                        outcomes[outcome.chrom] = outcome
                        submit_next()
                    completed.clear()
            except mp.ProcessError as exc:  # BrokenProcessPool subclasses this
                raise LDSCInternalError(
                    "ldscore parallel chromosome computation aborted: a worker process "
                    "crashed (most likely out-of-memory or a native library fault). "
                    "Re-run with a lower `--threads`, or `--threads 1` to isolate "
                    "the failing chromosome."
                ) from exc
        return outcomes

    def compute_chromosome(
        self,
        chrom: str,
        annotation_bundle,
        ref_panel,
        ldscore_config: LDScoreConfig,
        global_config: GlobalConfig,
        regression_snps: set[str] | RestrictionIdentityKeys | None = None,
        regression_regions: kernel_regions.RegionIntervals | None = None,
        export_dir: str | None = None,
    ) -> ChromLDScoreResult:
        """Compute normalized LD-score outputs for one chromosome.

        ``ref_panel.prepare_chromosome`` aligns annotations to the filtered
        reference universe ``B_chrom ∩ A'_chrom`` and owns the reader until the
        kernel returns or raises. The same aligned annotation values feed
        numerical projection, counts, overlap, and result provenance.
        ``regression_snps`` defines weight contributions and later restricts
        the normalized output rows.
        """
        backend = getattr(getattr(ref_panel, "spec", None), "backend", None)
        legacy_bundle = _kernel_annotation_bundle(annotation_bundle,chrom)
        LOGGER.info(f"Computing chromosome {chrom} LD scores with backend '{backend or 'unknown'}' and {len(legacy_bundle.metadata)} annotation rows.")
        with ref_panel.prepare_chromosome(chrom, legacy_bundle, ldscore_config) as prepared:
            legacy_result = kernel_ldscore.compute_chromosome(
                chrom, prepared, snp_identifier=global_config.snp_identifier,
                snp_batch_size=ldscore_config.snp_batch_size,
                common_maf_min=ldscore_config.common_maf_min,
                regression_keys=regression_snps, regression_regions=regression_regions,
                query_batch_size=ldscore_config.query_batch_size,
            )
            if export_dir is not None:
                _write_one_ref_metadata_sidecar(prepared.metadata, chrom, export_dir)
        result = self._wrap_legacy_chrom_result(
            legacy_result,
            global_config=global_config,
            regression_snps=regression_snps,
            regression_regions=regression_regions,
            annotation_types=legacy_result.annotation_types,
        )
        LOGGER.info(f"Finished chromosome {chrom} with {len(result.baseline_table)} retained SNP rows.")
        return result

    def _wrap_legacy_chrom_result(
        self,
        legacy_result: kernel_ldscore.ChromComputationResult,
        global_config: GlobalConfig,
        regression_snps: set[str] | RestrictionIdentityKeys | None = None,
        regression_regions: kernel_regions.RegionIntervals | None = None,
        annotation_types: dict[str, str] | None = None,
    ) -> ChromLDScoreResult:
        """Convert one kernel chromosome result into the typed public result."""
        reference_metadata = legacy_result.metadata.reset_index(drop=True).copy()
        ld_scores = pd.DataFrame(legacy_result.ld_scores, columns=list(legacy_result.ldscore_columns))
        regression_keep = pd.Series(True, index=reference_metadata.index)
        ld_regression_snps = frozenset(
            build_snp_id_series(
                reference_metadata.loc[regression_keep],
                global_config.snp_identifier,
            )
        )
        pos_column = "POS" if "POS" in reference_metadata.columns else "BP"
        metadata_columns = ["CHR", "SNP", pos_column, *[column for column in ("A1", "A2") if column in reference_metadata.columns]]
        regression_weights = np.asarray(legacy_result.w_ld, dtype=np.float32).reshape(-1)
        ldscore_table = pd.concat(
            [
                reference_metadata.loc[regression_keep, metadata_columns].rename(columns={pos_column: "POS"}).reset_index(drop=True),
                ld_scores.loc[regression_keep].reset_index(drop=True),
                pd.DataFrame({REGRESSION_LD_SCORE_COLUMN: regression_weights[regression_keep.to_numpy()]}).reset_index(drop=True),
            ],
            axis=1,
        )
        ldscore_table = kernel_ldscore.sort_frame_by_genomic_position(ldscore_table)
        count_map = {"all_reference_snp_counts": np.asarray(legacy_result.M, dtype=np.float64)}
        if legacy_result.M_5_50 is not None:
            count_map["common_reference_snp_counts"] = np.asarray(legacy_result.M_5_50, dtype=np.float64)
        baseline_table, query_table = _split_ldscore_table(
            ldscore_table,
            baseline_columns=list(legacy_result.baseline_columns),
            query_columns=list(legacy_result.query_columns),
            snp_identifier=global_config.snp_identifier,
        )
        result = ChromLDScoreResult(
            chrom=str(legacy_result.chrom),
            baseline_table=baseline_table,
            query_table=query_table,
            count_records=_count_records_from_totals(
                baseline_columns=list(legacy_result.baseline_columns),
                query_columns=list(legacy_result.query_columns),
                count_totals=count_map,
            ),
            baseline_columns=list(legacy_result.baseline_columns),
            query_columns=list(legacy_result.query_columns),
            ld_reference_snps=frozenset(),
            ld_regression_snps=ld_regression_snps,
            snp_count_totals=count_map,
            count_config={},
            config_snapshot=global_config,
            overlap=getattr(legacy_result, "overlap", None),
            reference_snp_count=legacy_result.reference_snp_count,
            regression_selected_snp_count=legacy_result.regression_selected_snp_count,
            regression_region_removed_snp_count=legacy_result.regression_region_removed_snp_count,
            identity_drops=getattr(legacy_result, "identity_drops", empty_identity_drop_frame()),
            annotation_types=dict(annotation_types or {}),
        )
        result.validate()
        return result

    def _aggregate_chromosome_results(
        self,
        chromosome_results: Sequence[ChromLDScoreResult],
        global_config: GlobalConfig,
        count_config: dict[str, Any] | None = None,
    ) -> LDScoreResult:
        """Concatenate and sum per-chromosome results into one aggregate object.

        The annotation overlap matrix is aggregated only for partitioned runs (two
        or more LD-score columns). Single-annotation runs (e.g. the synthetic
        unpartitioned ``base``) leave ``LDScoreResult.overlap`` as ``None`` so the
        downstream writer emits no redundant overlap artifact.
        """
        if not chromosome_results:
            raise LDSCInternalError(
                "LD-score aggregation failed in LDScoreCalculator._aggregate_chromosome_results(): "
                "no chromosome results were supplied. Most likely all chromosomes were skipped "
                "before aggregation. Re-run with DEBUG logging and report the traceback."
            )
        snapshots = [result.config_snapshot for result in chromosome_results if result.config_snapshot is not None]
        if snapshots:
            for snapshot in snapshots[1:]:
                validate_config_compatibility(
                    snapshots[0],
                    snapshot,
                    context="ChromLDScoreResult aggregation",
                )

        count_keys = sorted(set.intersection(*(set(result.snp_count_totals) for result in chromosome_results)))
        count_totals = {
            key: np.sum(np.vstack([result.snp_count_totals[key] for result in chromosome_results]), axis=0)
            for key in count_keys
        }
        from .overlap_matrix import LDScoreOverlap

        # The overlap (Gram) matrix is only consumed by partitioned-h2 and the
        # design-matrix collinearity check, both of which need >=2 LD-score
        # columns. A single-column run (e.g. the synthetic unpartitioned `base`)
        # yields a degenerate 1x1 overlap equal to a SNP count already in
        # metadata, so the artifact is suppressed rather than written redundantly.
        n_ld_columns = len(chromosome_results[0].baseline_columns) + len(chromosome_results[0].query_columns)
        overlaps = [result.overlap for result in chromosome_results if result.overlap is not None]
        aggregated_overlap = None
        if n_ld_columns >= 2 and len(overlaps) == len(chromosome_results):
            aggregated_overlap = LDScoreOverlap.from_contribution(
                sum_overlap_contributions(overlaps),
                list(chromosome_results[0].baseline_columns),
                list(chromosome_results[0].query_columns),
            )
        merged_table = pd.concat(
            [
                _join_split_tables(
                    result.baseline_table,
                    result.query_table,
                    result.query_columns,
                    snp_identifier=getattr(result.config_snapshot, "snp_identifier", global_config.snp_identifier),
                )
                for result in chromosome_results
            ],
            axis=0,
            ignore_index=True,
        )
        merged_table = kernel_ldscore.sort_frame_by_genomic_position(merged_table)
        _assert_canonical_maf(merged_table)
        baseline_table, query_table = _split_ldscore_table(
            merged_table,
            baseline_columns=list(chromosome_results[0].baseline_columns),
            query_columns=list(chromosome_results[0].query_columns),
            snp_identifier=global_config.snp_identifier,
        )
        annotation_names = [
            *chromosome_results[0].baseline_columns,
            *chromosome_results[0].query_columns,
        ]
        annotation_types = {
            name: (
                "quantitative"
                if any(chrom_result.annotation_types.get(name) == "quantitative" for chrom_result in chromosome_results)
                else "binary"
            )
            for name in annotation_names
        }
        result = LDScoreResult(
            baseline_table=baseline_table,
            query_table=query_table,
            count_records=_count_records_from_totals(
                baseline_columns=list(chromosome_results[0].baseline_columns),
                query_columns=list(chromosome_results[0].query_columns),
                count_totals=count_totals,
            ),
            baseline_columns=list(chromosome_results[0].baseline_columns),
            query_columns=list(chromosome_results[0].query_columns),
            ld_reference_snps=frozenset(),
            ld_regression_snps=frozenset().union(*(result.ld_regression_snps for result in chromosome_results)),
            chromosome_results=list(chromosome_results),
            count_config=dict(count_config or {}),
            config_snapshot=snapshots[0] if snapshots else None,
            overlap=aggregated_overlap,
            identity_drops_by_chrom={
                result.chrom: result.identity_drops.copy()
                for result in chromosome_results
            },
            annotation_types=annotation_types,
        )
        result.validate()
        return result

    def write_outputs(
        self,
        result: LDScoreResult,
        output_config: LDScoreOutputConfig,
        config_snapshot: dict[str, Any] | None = None,
    ):
        """Write a previously computed result as a canonical LD-score directory.

        Parameters
        ----------
        result : LDScoreResult
            Aggregate LD-score result to serialize.
        output_config : LDScoreOutputConfig
            Directory path and overwrite/compression controls. Existing
            canonical files are refused unless ``output_config.overwrite`` is
            true.
        config_snapshot : dict or None, optional
            Optional metadata recorded in the emitted run summary. Default is
            ``None``.

        Returns
        -------
        dict
            Resolved output paths keyed by artifact name.
        """
        del config_snapshot
        return self.output_writer.write(result, output_config)


def _assert_canonical_maf(metadata: pd.DataFrame) -> None:
    """Verify the canonical invariant that MAF = freq(A1) is the minor allele.

    LD-score artifacts inherit A1/A2/MAF from the reference panel. A panel built
    before allele-orientation canonicalization can carry MAF > 0.5 (folded MAF not
    tied to A1); fail loudly rather than silently mis-tie MAF to A1.
    """
    if "MAF" not in metadata.columns:
        return
    maf = pd.to_numeric(metadata["MAF"], errors="coerce")
    over = maf[maf > 0.5 + 1e-9]
    if len(over):
        raise LDSCInternalError(
            "LD-score metadata carries MAF > 0.5, violating the canonical "
            f"A1=minor invariant ({len(over)} rows; first={float(over.iloc[0]):.4f}). "
            "Most likely the reference panel was built before allele-orientation "
            "canonicalization. Rebuild the reference panel with the current package."
        )


def _split_ldscore_table(
    ldscore_table: pd.DataFrame,
    *,
    baseline_columns: list[str],
    query_columns: list[str],
    snp_identifier: str = "chr_pos_allele_aware",
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Split a merged LD-score table into baseline and optional query tables."""
    if is_allele_aware_mode(snp_identifier) and not {"A1", "A2"}.issubset(ldscore_table.columns):
        raise LDSCInternalError(
            "LD-score split-table assembly failed in _split_ldscore_table(): allele-aware "
            f"SNP identity mode '{snp_identifier}' requires A1/A2 columns, but the LD-score "
            "table does not contain both. Most likely an upstream reference-panel or annotation "
            "assembly step dropped allele columns. Re-run with DEBUG logging and report the traceback."
        )
    metadata_columns = ["CHR", "SNP", "POS", *[column for column in ("A1", "A2") if column in ldscore_table.columns]]
    baseline_order = [*metadata_columns, REGRESSION_LD_SCORE_COLUMN, *baseline_columns]
    query_order = [*metadata_columns, *query_columns]
    missing_baseline = [column for column in baseline_order if column not in ldscore_table.columns]
    if missing_baseline:
        raise LDSCInternalError(
            "LD-score split-table assembly failed in _split_ldscore_table(): "
            f"the merged table is missing baseline columns {missing_baseline}. Most likely "
            "a kernel result omitted expected baseline annotation scores. Re-run with DEBUG "
            "logging and report the traceback."
        )
    baseline_table = ldscore_table.loc[:, baseline_order].reset_index(drop=True).copy()
    query_table = None
    if query_columns:
        missing_query = [column for column in query_order if column not in ldscore_table.columns]
        if missing_query:
            raise LDSCInternalError(
                "LD-score split-table assembly failed in _split_ldscore_table(): "
                f"the merged table is missing query columns {missing_query}. Most likely "
                "a kernel result omitted expected query annotation scores. Re-run with DEBUG "
                "logging and report the traceback."
            )
        query_table = ldscore_table.loc[:, query_order].reset_index(drop=True).copy()
        assert_same_snp_rows(
            baseline_table,
            query_table,
            context="query rows must match baseline rows on CHR/SNP/POS",
            snp_identifier=snp_identifier,
        )
    return baseline_table, query_table


def _join_split_tables(
    baseline_table: pd.DataFrame,
    query_table: pd.DataFrame | None,
    query_columns: Sequence[str],
    *,
    snp_identifier: str = "chr_pos_allele_aware",
) -> pd.DataFrame:
    """Join split LD-score tables for sorting or regression assembly."""
    if query_table is None:
        return baseline_table.copy()
    assert_same_snp_rows(
        baseline_table,
        query_table,
        context="query rows must match baseline rows on CHR/SNP/POS",
        snp_identifier=snp_identifier,
    )
    query_values = query_table.loc[:, list(query_columns)].reset_index(drop=True)
    return pd.concat([baseline_table.reset_index(drop=True), query_values], axis=1)


def _count_records_from_totals(
    *,
    baseline_columns: list[str],
    query_columns: list[str],
    count_totals: dict[str, np.ndarray],
) -> list[dict[str, Any]]:
    """Convert positional count vectors into metadata-friendly column records."""
    columns = [*baseline_columns, *query_columns]
    groups = ["baseline"] * len(baseline_columns) + ["query"] * len(query_columns)
    all_counts = np.asarray(count_totals.get("all_reference_snp_counts"), dtype=np.float64)
    if all_counts.size != len(columns):
        raise LDSCInternalError(
            "LD-score count assembly failed in _count_records_from_totals(): "
            f"all_reference_snp_counts has length {all_counts.size}, but there are "
            f"{len(columns)} annotation columns. Most likely the kernel returned count "
            "vectors for a different annotation matrix. Re-run with DEBUG logging and "
            "report the traceback."
        )
    common_raw = count_totals.get("common_reference_snp_counts")
    common_counts = None if common_raw is None else np.asarray(common_raw, dtype=np.float64)
    if common_counts is not None and common_counts.size != len(columns):
        raise LDSCInternalError(
            "LD-score count assembly failed in _count_records_from_totals(): "
            f"common_reference_snp_counts has length {common_counts.size}, but there are "
            f"{len(columns)} annotation columns. Most likely the kernel returned common-SNP "
            "count vectors for a different annotation matrix. Re-run with DEBUG logging and "
            "report the traceback."
        )
    records: list[dict[str, Any]] = []
    for idx, (group, column) in enumerate(zip(groups, columns)):
        record: dict[str, Any] = {
            "group": group,
            "column": column,
            "all_reference_snp_count": float(all_counts[idx]),
        }
        if common_counts is not None:
            record["common_reference_snp_count"] = float(common_counts[idx])
        records.append(record)
    return records


def _count_config_from_ldscore_config(ldscore_config: LDScoreConfig) -> dict[str, Any]:
    """Return count metadata for common-SNP count vectors."""
    return {
        "common_reference_snp_maf_min": float(ldscore_config.common_maf_min),
        "common_reference_snp_maf_operator": ">=",
    }


def build_parser() -> argparse.ArgumentParser:
    """Build the feature parser for LD-score calculation."""
    parser = argparse.ArgumentParser(
        description="Estimate LDSC-compatible LD scores from SNP-level annotation files using PLINK or sorted parquet R2 input.",
        allow_abbrev=False,
    )
    parser.add_argument("--output-dir", required=True, help="Output directory for the canonical LD-score result.")
    parser.add_argument(
        "--gene-ldscore-index-dir",
        default=None,
        help="Explicit complete gene LD-score index directory. Requires gene-list queries and forbids live reference/baseline inputs.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Replace LD-score output artifacts and remove stale owned siblings.",
    )
    query_group = parser.add_mutually_exclusive_group()
    query_group.add_argument(
        "--query-annot-sources",
        default=None,
        help="Comma-separated query annotation path tokens: exact paths, globs, or explicit @ suite tokens. Requires --baseline-annot-sources.",
    )
    query_group.add_argument(
        "--query-annot-bed-sources",
        default=None,
        help="Comma-separated BED file path tokens projected in memory as query annotations. Requires --baseline-annot-sources.",
    )
    query_group.add_argument(
        "--query-annot-gene-list-sources",
        default=None,
        help=(
            "Comma-separated one-column gene-list exact paths or glob patterns projected in memory "
            "as query annotations. Direct mode requires --baseline-annot-sources, "
            "--gene-coordinate-file, and an explicit --padding-bp."
        ),
    )
    parser.add_argument(
        "--padding-bp",
        type=int,
        default=None,
        help=(
            "Base pairs to add to both sides of each live BED- or gene-list-derived interval "
            "before in-memory projection. Omission means 0 for BED mode; live gene-list mode "
            "requires an explicit value, including --padding-bp 0 for gene bodies. Forbidden "
            "with prebuilt annotations, no query, or --gene-ldscore-index-dir."
        ),
    )
    parser.add_argument(
        "--gene-coordinate-file",
        default=None,
        help=(
            "Required headered TSV/TSV.GZ coordinate catalog for live gene-list mode. "
            "Coordinates are one-based inclusive and the catalog is the sole gene-resolution authority."
        ),
    )
    parser.add_argument(
        "--gene-list-resolution-policy",
        choices=("strict", "resolved-only"),
        default="strict",
        help=(
            "Gene identifier policy. 'strict' stops on any rejected identifier; 'resolved-only' "
            "explicitly continues with the audited usable subset. Default: strict."
        ),
    )
    parser.add_argument(
        "--gene-exclude-regions",
        choices=("none", "mhc"),
        default="none",
        help="Gene regions excluded before padding in gene-list workflows. Default: none.",
    )
    parser.add_argument(
        "--control-gene-list-file",
        default=None,
        help="Optional fixed-control gene-list file. Omit to add no gene control.",
    )
    parser.add_argument(
        "--baseline-annot-sources",
        default=None,
        help="Comma-separated baseline annotation path tokens. If omitted with no query inputs, an all-ones `base` annotation is synthesized.",
    )
    parser.add_argument(
        "--plink-prefix",
        default=None,
        help=(
            "PLINK prefix shared by a complete .bed/.bim/.fam trio, or a plain stem that discovers "
            "chromosome-coded complete trios. Globs and @ chromosome patterns are also supported."
        ),
    )
    parser.add_argument(
        "--r2-dir",
        dest="r2_dir",
        default=None,
        help="Build-specific R2 directory containing chr*_r2.parquet and optional chr*_meta.tsv.gz sidecars.",
    )
    parser.add_argument(
        "--snp-identifier",
        default="chr_pos_allele_aware",
        choices=("rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"),
        help="Identifier mode used to match annotations to the reference panel.",
    )
    parser.add_argument(
        "--genome-build",
        choices=("auto", "hg19", "hg37", "GRCh37", "hg38", "GRCh38"),
        default=None,
        help=(
            "Genome build for chr_pos-family inputs and gene-list interval projection. "
            "Required when --snp-identifier is a chr_pos mode; gene-list runs default "
            "to 'auto'. Use 'auto' to infer hg19/hg38 from baseline/reference-panel "
            "evidence. In rsid-family gene-list runs, the resolved build selects both "
            "gene projection intervals and named regression-region presets."
        ),
    )
    parser.add_argument(
        "--ref-panel-snps-file",
        default=None,
        help=(
            "Optional identity-only SNP list defining the retained reference-panel universe A'. "
            "Duplicate restriction keys collapse to one retained key; non-identity columns such as CM or MAF are ignored. "
            "The workflow intersects each chromosome annotation bundle with this prepared panel before LD computation."
        ),
    )
    parser.add_argument(
        "--regr-snps-exclude-regions",
        choices=REGR_SNPS_EXCLUDE_REGIONS_CHOICES,
        default="mhc-and-centromeres",
        help="Curated region presets subtracted from regression/output SNPs after selecting bundled HM3 or --regr-snps-file. Baseline/query LD-score contributors and M/overlap counts remain unchanged. Defaults to mhc-and-centromeres; use 'none' to keep all regions.",
    )
    parser.add_argument(
        "--exclude-regions",
        dest="regr_snps_exclude_regions",
        choices=REGR_SNPS_EXCLUDE_REGIONS_CHOICES,
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--regr-snps-file",
        default=None,
        help=(
            "Optional identity-only SNP list defining the regression SNP set and the written LD-score row set. "
            "Duplicate restriction keys collapse to one retained key; non-identity columns such as CM or MAF are ignored."
        ),
    )
    parser.add_argument(
        "--keep-indivs-file",
        default=None,
        help="File with individuals to include in LD Score estimation. The file should contain one IID per row.",
    )
    parser.add_argument("--ld-wind-snps", default=None, type=int, help="LD window size in SNPs.")
    parser.add_argument("--ld-wind-kb", default=None, type=float, help="LD window size in kilobases.")
    parser.add_argument("--ld-wind-cm", default=None, type=float, help="LD window size in centiMorgans.")
    parser.add_argument("--maf-min", default=None, type=float, help="Optional MAF filter for retained reference-panel SNPs when MAF is available.")
    parser.add_argument("--common-maf-min", default=0.05, type=float, help="MAF threshold used only for common-SNP annotation count vectors.")
    parser.add_argument("--genetic-map-hg19-sources", default=None, help="Genetic map (hg19) used to derive CM for cM windows when a PLINK .bim CM column is uninformative.")
    parser.add_argument("--genetic-map-hg38-sources", default=None, help="Genetic map (hg38) used to derive CM for cM windows when a PLINK .bim CM column is uninformative.")
    parser.add_argument("--export-ref-metadata", default=False, action="store_true", help="Write a chrN_meta.tsv.gz reference-metadata sidecar next to the LD-score output (PLINK backend only).")
    parser.add_argument("--snp-batch-size", default=128, type=int, help="Genotype batch size for the PLINK reference-panel backend; ignored by the parquet-R2 backend, which streams stored pairs. Defaults to 128.")
    parser.add_argument("--query-batch-size", default=1000, type=int, help="Maximum active focal query columns per chromosome projection. Default: 1000; independent of chromosome workers.")
    parser.add_argument("--threads", default=1, type=int, help="Worker processes for cross-chromosome parallelism (joblib n_jobs convention): 1=sequential (default), N=N workers, -1=all cores, -2=all but one. Respects CPU affinity; capped at the chromosome count.")
    parser.add_argument("--yes-really", default=False, action="store_true", help="Allow whole-chromosome LD windows.")
    parser.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"), help="Logging verbosity.")
    return parser


@materializing_overwrite_guard(
    lambda args: (
        (getattr(args, "output_dir"), getattr(args, "overwrite", False), "RUN_FAILED.txt")
        if getattr(args, "output_dir", None)
        else None
    ),
    command="run_ldscore_from_args(...)",
)
def run_ldscore_from_args(args: argparse.Namespace) -> LDScoreResult:
    """Run LD-score calculation from a parsed CLI namespace.

    The workflow resolves unified path tokens for baseline annotations, optional
    prebuilt, BED, or gene-list queries, and the reference panel. If no baseline
    or query annotations are supplied, it synthesizes an all-ones ``base``
    annotation over the retained reference-panel metadata. Before calculation
    it preflights ``metadata.json``, ``ldscore.baseline.parquet``, optional
    ``ldscore.query.parquet``, conditional ``ldscore.overlap.parquet``, query
    diagnostics when applicable, and ``diagnostics/ldscore.log`` under
    ``output_dir``. With
    overwrite enabled, successful baseline-only runs remove stale query parquet
    siblings. Each chromosome uses ``ref_panel.prepare_chromosome`` to bind
    aligned annotations to one configured reader before numerical projection.
    The workflow returns the normalized public ``LDScoreResult`` with split
    baseline/query tables.
    The result ``output_paths`` mapping contains data artifacts only.
    """
    _validate_padding_usage(args)
    _validate_gene_list_mode_args(args)
    if getattr(args, "gene_ldscore_index_dir", None) is not None:
        return _run_explicit_indexed_ldscore(args)

    from .gene_list_resolver import GeneCatalog
    from ._gene_query_storage import resolve_gene_lists_staged, persistent_gene_diagnostics
    from ._annotation_storage import AnnotationWorkspace
    from ._direct_annotation import prepare_direct_annotations, prepare_synthetic_base

    has_queries = any(_has_cli_tokens(getattr(args, name, None)) for name in (
        "query_annot_gene_list_sources", "query_annot_bed_sources", "query_annot_sources",
    ))
    if has_queries and not _has_cli_tokens(getattr(args, "baseline_annot_sources", None)):
        raise LDSCUsageError(_QUERY_REQUIRES_BASELINE_MESSAGE)
    normalized_args = global_config = None
    if not _has_cli_tokens(getattr(args,"query_annot_gene_list_sources",None)):
        normalized_args, global_config = _normalize_run_args(args)
        _validate_run_args(normalized_args)
    output_config = _output_config_from_args(args)
    output_dir = ensure_output_directory(output_config.output_dir, label="LD-score output directory")
    log_path = output_dir / "diagnostics" / "ldscore.log"
    LDScoreDirectoryWriter.artifact_family(output_dir).preflight(
        overwrite=output_config.overwrite, additional_paths=[log_path])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with AnnotationWorkspace(output_dir) as workspace:
        with workflow_logging("ldscore", log_path, log_level=getattr(args,"log_level","INFO")):
            batch = catalog = None
            if _has_cli_tokens(getattr(args,"query_annot_gene_list_sources",None)):
                catalog = GeneCatalog.load(args.gene_coordinate_file)
                _log_live_gene_catalog_issues(catalog)
                batch = resolve_gene_lists_staged(split_cli_path_tokens(args.query_annot_gene_list_sources),catalog,workspace,
                    control_path=getattr(args,"control_gene_list_file",None),
                    resolution_policy=getattr(args,"gene_list_resolution_policy","strict"),
                    gene_exclude_regions=getattr(args,"gene_exclude_regions","none"))
                _log_gene_list_rejections(batch)
            if normalized_args is None:
                normalized_args, global_config = _normalize_run_args(args,gene_catalog=catalog)
            print_global_config_banner("run_ldscore_from_args",global_config)
            _validate_run_args(normalized_args)
            ldscore_config = _ldscore_config_from_args(normalized_args)
            regression_snps = _load_regression_snps(_regr_snps_file_from_config(ldscore_config),global_config,
                label="packaged HM3 regression SNP map" if ldscore_config.regr_snps_file is None else "regression SNP list")
            regression_regions = _regression_region_intervals(normalized_args,global_config)
            scope = {}
            if _has_cli_tokens(normalized_args.baseline_annot_sources):
                annotation_bundle, scope = prepare_direct_annotations(normalized_args,global_config,
                    _annotation_build_config_from_args(normalized_args),workspace,output_config,batch=batch)
                ref_panel = _ref_panel_from_args(normalized_args,global_config,
                    chromosome_prefixes=scope.get("reference_prefixes_by_chrom"))
            else:
                ref_panel = _ref_panel_from_args(normalized_args,global_config)
                annotation_bundle = prepare_synthetic_base(ref_panel,global_config,workspace)
            if scope:
                LOGGER.info("Chromosomes resolved and entering the analysis: %s.", ", ".join(scope['chromosomes']))
            log_inputs(output_dir=str(output_dir),reference_mode=ref_panel.spec.backend,
                snp_identifier=global_config.snp_identifier,genome_build=global_config.genome_build,
                ref_panel_snps_file=normalized_args.ref_panel_snps_file or "none",
                regr_snps_file=ldscore_config.regr_snps_file or "packaged_hm3_default",
                regression_region_presets=", ".join(regression_regions.source_labels) if regression_regions else "none")
            if getattr(normalized_args,"gene_catalog_build",None) is not None:
                LOGGER.info("Gene-list catalog projection build: %s (selected from --genome-build and baseline/reference-panel evidence).",normalized_args.gene_catalog_build)
            calculator = LDScoreCalculator()
            if annotation_bundle.gene_list_batch is not None:
                annotation_bundle, gate_b_error = _apply_direct_gene_gate_b(annotation_bundle,ref_panel,ldscore_config,global_config)
                _log_gene_list_snp_support(annotation_bundle.gene_list_batch)
                if gate_b_error:
                    log_outputs(**calculator.output_writer.write_query_diagnostics(annotation_bundle,output_config))
                    raise LDSCInputError(gate_b_error)
            statuses = annotation_bundle.query_statuses
            if statuses and not annotation_bundle.query_columns:
                _log_query_annotation_statuses(statuses)
                log_outputs(**calculator.output_writer.write_query_diagnostics(annotation_bundle,output_config))
                raise LDSCInputError(_all_query_annotations_skipped_message(statuses))
            result = calculator.run(annotation_bundle,ref_panel,ldscore_config,global_config,
                output_config=output_config,regression_snps=regression_snps,regression_regions=regression_regions)
            _log_query_annotation_statuses(result.query_statuses)
            log_outputs(**result.output_paths)
            if bool(getattr(args,"_emit_gene_console_notices",False)):
                _emit_gene_gate_b_notice(result.query_statuses,result.gene_list_batch)
                _emit_resolved_only_notice(result.gene_list_batch)
            return dataclass_replace(result,gene_list_batch=persistent_gene_diagnostics(result.gene_list_batch,result.output_paths))


def _run_explicit_indexed_ldscore(args: argparse.Namespace) -> LDScoreResult:
    """Validate and dispatch the closed indexed gene-list mode."""
    explicit_options = set(getattr(args, "_explicit_cli_options", ()))
    live_options = {
        "--baseline-annot-sources",
        "--query-annot-sources",
        "--query-annot-bed-sources",
        "--padding-bp",
        "--gene-coordinate-file",
        "--gene-exclude-regions",
        "--plink-prefix",
        "--r2-dir",
        "--snp-identifier",
        "--genome-build",
        "--ref-panel-snps-file",
        "--regr-snps-file",
        "--regr-snps-exclude-regions",
        "--exclude-regions",
        "--keep-indivs-file",
        "--ld-wind-snps",
        "--ld-wind-kb",
        "--ld-wind-cm",
        "--maf-min",
        "--common-maf-min",
        "--genetic-map-hg19-sources",
        "--genetic-map-hg38-sources",
        "--export-ref-metadata",
        "--snp-batch-size",
        "--threads",
        "--yes-really",
    }
    forbidden = {
        "--baseline-annot-sources": getattr(args, "baseline_annot_sources", None),
        "--query-annot-sources": getattr(args, "query_annot_sources", None),
        "--query-annot-bed-sources": getattr(args, "query_annot_bed_sources", None),
        "--gene-coordinate-file": getattr(args, "gene_coordinate_file", None),
        "--plink-prefix": getattr(args, "plink_prefix", None),
        "--r2-dir": getattr(args, "r2_dir", None),
        "--ref-panel-snps-file": getattr(args, "ref_panel_snps_file", None),
        "--regr-snps-file": getattr(args, "regr_snps_file", None),
        "--keep-indivs-file": getattr(args, "keep_indivs_file", None),
        "--genetic-map-hg19-sources": getattr(args, "genetic_map_hg19_sources", None),
        "--genetic-map-hg38-sources": getattr(args, "genetic_map_hg38_sources", None),
        "--ld-wind-snps": getattr(args, "ld_wind_snps", None),
        "--ld-wind-kb": getattr(args, "ld_wind_kb", None),
        "--ld-wind-cm": getattr(args, "ld_wind_cm", None),
        "--maf-min": getattr(args, "maf_min", None),
    }
    if getattr(args, "gene_exclude_regions", "none") != "none":
        forbidden["--gene-exclude-regions"] = getattr(args, "gene_exclude_regions")
    if "--genome-build" in explicit_options or getattr(args, "genome_build", None) is not None:
        forbidden["--genome-build"] = getattr(args, "genome_build")
    if (
        "--snp-identifier" in explicit_options
        or getattr(args, "snp_identifier", "chr_pos_allele_aware") != "chr_pos_allele_aware"
    ):
        forbidden["--snp-identifier"] = getattr(args, "snp_identifier")
    if getattr(args, "regr_snps_exclude_regions", "mhc-and-centromeres") != "mhc-and-centromeres":
        forbidden["--regr-snps-exclude-regions"] = getattr(args, "regr_snps_exclude_regions")
    if getattr(args, "common_maf_min", 0.05) != 0.05:
        forbidden["--common-maf-min"] = getattr(args, "common_maf_min")
    if getattr(args, "snp_batch_size", 128) != 128:
        forbidden["--snp-batch-size"] = getattr(args, "snp_batch_size")
    if getattr(args, "threads", 1) != 1:
        forbidden["--threads"] = getattr(args, "threads")
    if bool(getattr(args, "export_ref_metadata", False)):
        forbidden["--export-ref-metadata"] = True
    if bool(getattr(args, "yes_really", False)):
        forbidden["--yes-really"] = True
    supplied = sorted(
        (explicit_options & live_options)
        | {
            option
            for option, value in forbidden.items()
            if value not in {None, ""}
        }
    )
    if supplied:
        raise LDSCInputError(
            "ldscore indexed mode inherits immutable scientific inputs, SNP identity, and genome build "
            "from the index and cannot accept live overrides: "
            + ", ".join(supplied)
            + ". Remove these options, or remove --gene-ldscore-index-dir to run direct mode."
        )
    gene_lists = split_cli_path_tokens(getattr(args, "query_annot_gene_list_sources", None))
    if not gene_lists:
        raise LDSCInputError(
            "ldscore indexed mode requires --query-annot-gene-list-sources."
        )
    from .gene_ldscore_index import run_indexed_ldscore
    output_dir = Path(args.output_dir)
    log_path = output_dir / "diagnostics" / "ldscore.log"
    LDScoreDirectoryWriter.artifact_family(output_dir).preflight(
        overwrite=bool(getattr(args, "overwrite", False)), additional_paths=[log_path],
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with workflow_logging(
        "ldscore",
        log_path,
        log_level=getattr(args, "log_level", "INFO"),
    ):
        log_inputs(
            output_dir=str(output_dir),
            reference_mode="gene_ldscore_index",
            gene_ldscore_index=str(args.gene_ldscore_index_dir),
            query_gene_list_count=len(gene_lists),
            control_gene_list_file=getattr(args, "control_gene_list_file", None),
        )
        result = run_indexed_ldscore(
            args.gene_ldscore_index_dir,
            query_gene_list_sources=tuple(gene_lists),
            control_gene_list_file=getattr(args, "control_gene_list_file", None),
            gene_list_resolution_policy=getattr(args, "gene_list_resolution_policy", "strict"),
            output_dir=args.output_dir,
            overwrite=bool(getattr(args, "overwrite", False)),
        )
        if bool(getattr(args, "_emit_gene_console_notices", False)):
            _emit_gene_gate_b_notice(
                getattr(result, "query_statuses", ()),
                getattr(result, "gene_list_batch", None),
            )
            _emit_resolved_only_notice(getattr(result, "gene_list_batch", None))
        output_paths = getattr(result, "output_paths", None)
        if output_paths:
            log_outputs(**output_paths)
        return result


def _validate_padding_usage(args: argparse.Namespace) -> None:
    """Reject explicit padding outside live BED or gene-list projection.

    The parser's ``None`` sentinel distinguishes omission from an explicit
    zero; normalization converts an allowed omission to the effective default
    only after this mode check.
    """
    has_gene_lists = _has_cli_tokens(getattr(args, "query_annot_gene_list_sources", None))
    uses_index = getattr(args, "gene_ldscore_index_dir", None) is not None
    if getattr(args, "padding_bp", None) is None:
        if has_gene_lists and not uses_index:
            raise LDSCUsageError(
                "Live gene-list mode requires an explicit --padding-bp value; it must be chosen deliberately. "
                "Use --padding-bp 0 for gene bodies or a positive value for padded intervals."
            )
        return
    uses_live_intervals = (
        _has_cli_tokens(getattr(args, "query_annot_bed_sources", None))
        or _has_cli_tokens(getattr(args, "query_annot_gene_list_sources", None))
    )
    if uses_live_intervals and not uses_index:
        return
    raise LDSCUsageError(
        "ldscore accepts `--padding-bp` only with live `--query-annot-bed-sources` "
        "or `--query-annot-gene-list-sources`. It cannot be used with prebuilt "
        "annotation queries, without a query, or in indexed mode. Remove "
        "`--padding-bp` from the command (or remove the Python `padding_bp` keyword)."
    )


def _validate_gene_list_mode_args(args: argparse.Namespace) -> None:
    """Reject gene-only options outside their direct or indexed input mode."""
    has_gene_lists = _has_cli_tokens(getattr(args, "query_annot_gene_list_sources", None))
    uses_index = getattr(args, "gene_ldscore_index_dir", None) is not None
    explicit_options = set(getattr(args, "_explicit_cli_options", ()))
    explicitly_supplied_gene_options = explicit_options & {
        "--gene-coordinate-file",
        "--gene-exclude-regions",
        "--control-gene-list-file",
        "--gene-list-resolution-policy",
    }
    if not has_gene_lists:
        if (
            explicitly_supplied_gene_options
            or getattr(args, "gene_exclude_regions", "none") != "none"
            or getattr(args, "control_gene_list_file", None) is not None
            or getattr(args, "gene_coordinate_file", None) is not None
            or getattr(args, "gene_list_resolution_policy", "strict") != "strict"
        ):
            raise LDSCUsageError(
                "Gene coordinate, exclusion, control-list, and resolution policy options are valid only with "
                "gene-list query annotations. Remove the gene-specific option or use gene-list mode."
            )
        return
    if not uses_index and getattr(args, "gene_coordinate_file", None) is None:
        raise LDSCUsageError(
            "Live gene-list mode requires --gene-coordinate-file. Supply a one-build, one-based "
            "gene-coordinate TSV/TSV.GZ; there is no packaged catalog fallback."
        )


def _annotation_build_config_from_args(args: argparse.Namespace) -> AnnotationBuildConfig:
    """Build the annotation source config after the public mode matrix passes."""
    return AnnotationBuildConfig(
        baseline_annot_sources=tuple(split_cli_path_tokens(args.baseline_annot_sources)),
        query_annot_sources=tuple(split_cli_path_tokens(args.query_annot_sources)),
        query_annot_bed_sources=tuple(split_cli_path_tokens(getattr(args, "query_annot_bed_sources", None))),
        query_annot_gene_list_sources=tuple(
            split_cli_path_tokens(getattr(args, "query_annot_gene_list_sources", None))
        ),
        gene_coordinate_file=getattr(args, "gene_coordinate_file", None),
        control_gene_list_file=getattr(args, "control_gene_list_file", None),
        gene_list_resolution_policy=getattr(args, "gene_list_resolution_policy", "strict"),
        gene_exclude_regions=getattr(args, "gene_exclude_regions", "none"),
        padding_bp=getattr(args, "padding_bp", 0),
    )


def _log_live_gene_catalog_issues(catalog: Any) -> None:
    """Log every nonstructural catalog defect discovered in live mode."""
    for row in catalog.issues.itertuples(index=False):
        LOGGER.warning(
            "Gene-coordinate catalog defect: source=%s line=%s field=%s reason=%s observed=%r",
            row.source,
            row.catalog_line,
            row.field,
            row.reason,
            row.observed_value,
        )


def _apply_direct_gene_gate_b(annotation_bundle, ref_panel, ldscore_config, global_config):
    """Measure interval support on the prepared computational SNP universe."""
    from .errors import EmptyReferenceSNPs, LDSCUserError

    batch = annotation_bundle.gene_list_batch
    padding_bp = int(annotation_bundle.source_summary.get("padding_bp", 0))
    intervals_by_index = {}
    for selection in batch.selections:
        intervals_by_index.update(zip(selection.catalog_indices,selection.intervals))
    interval_frame = pd.DataFrame({'catalog_index':list(intervals_by_index), 'interval':list(intervals_by_index.values())})
    interval_coordinates = pd.DataFrame(
        interval_frame.pop("interval").tolist(),
        columns=["chrom", "start0", "end"],
        index=interval_frame.index,
    )
    interval_frame = pd.concat([interval_frame, interval_coordinates], axis=1)
    interval_frame["catalog_index"] = interval_frame["catalog_index"].astype(int)
    interval_frame["start0"] = (
        interval_frame["start0"].astype(np.int64) - padding_bp
    ).clip(lower=0)
    interval_frame["end"] = interval_frame["end"].astype(np.int64) + padding_bp
    support = pd.Series(pd.NA, index=interval_frame["catalog_index"], dtype="Int64")
    errors = []
    input_issues = []
    for chrom in _chromosomes_from_bundle(annotation_bundle):
        genes = interval_frame[interval_frame["chrom"].astype(str).eq(str(chrom))]
        chrom_bundle = _slice_annotation_bundle(annotation_bundle, chrom)
        legacy_bundle = _kernel_annotation_bundle(chrom_bundle,chrom)
        try:
            with ref_panel.prepare_chromosome(chrom, legacy_bundle, ldscore_config) as prepared:
                positions0 = np.sort(prepared.metadata["POS"].to_numpy(dtype=np.int64) - 1)
        except EmptyReferenceSNPs:
            positions0 = np.array([], dtype=np.int64)
        except (OSError, EOFError, ValueError, LDSCUserError) as exc:
            errors.append(f"Chromosome {chrom} reference preparation failed: {exc}")
            input_issues.append({"input_role": "reference", "source": f"chromosome {chrom}", "chrom": str(chrom), "reason": "invalid_required_input",
                                 "details": str(exc), "repair": "Repair reference inputs and runtime filters before retrying."})
            continue
        counts = np.searchsorted(positions0, genes["end"].to_numpy(dtype=np.int64), side="left") - np.searchsorted(
            positions0, genes["start0"].to_numpy(dtype=np.int64), side="left"
        )
        support.loc[genes["catalog_index"].astype(int).to_numpy()] = counts
    updated_batch = batch.with_snp_support(support)
    statuses = gene_query_statuses(updated_batch)
    retained_queries = [status.query for status in statuses if status.status in {"ok", "warning"}]
    errors.extend(gene_viability_errors(updated_batch, statuses))
    changes = dict(query_columns=retained_queries, query_statuses=statuses, gene_list_batch=updated_batch,
                   input_issues=pd.DataFrame(input_issues) if input_issues else None)
    if not hasattr(annotation_bundle,'shards'):
        changes['query_annotations'] = annotation_bundle.query_annotations.loc[:,retained_queries].copy()
    return dataclass_replace(annotation_bundle,**changes), "; ".join(errors) if errors else None


def _emit_resolved_only_notice(batch: Any | None) -> None:
    """Emit the one bounded successful-run warning that must reach the console."""
    if batch is None or batch.resolution_policy != "resolved-only":
        return
    affected = batch.summary[batch.summary["rejected_rows"].fillna(0).gt(0)]
    if affected.empty:
        return
    total = int(batch.summary["nonblank_input_rows"].fillna(0).sum())
    rejected = int(batch.summary["rejected_rows"].fillna(0).sum())
    sources = ", ".join(
        f"{row.source} ({int(row.rejected_rows)}/{int(row.nonblank_input_rows)} omitted)"
        for row in affected.head(MAX_CONSOLE_GENE_ISSUES).itertuples(index=False)
    )
    extra = max(0, len(affected) - MAX_CONSOLE_GENE_ISSUES)
    message = (
        f"WARNING: --gene-list-resolution-policy resolved-only completed with {rejected} of {total} "
        f"submitted row(s) omitted: {sources}"
        + (f", plus {extra} more affected source(s)" if extra else "")
        + ". See diagnostics/gene_list_resolution_summary.tsv and diagnostics/gene_list_audit.tsv.gz."
    )
    LOGGER.warning(message.removeprefix("WARNING: "))
    print(message, file=sys.stderr)


def _emit_gene_gate_b_notice(
    statuses: Sequence[QueryAnnotationStatus],
    batch: Any | None,
) -> None:
    """Emit one bounded successful-run notice for gene support/query viability."""
    if batch is None:
        return
    unsupported_rows = sum(int(frame.disposition.eq("unsupported").sum()) for frame in batch.audit_frames())
    gate_b_reasons = {
        "partial_snp_support",
        "zero_annotation_snps",
        "zero_variance_ld_scores",
    }
    affected = [
        status
        for status in statuses
        if status.reason in gate_b_reasons
        or "zero reference-SNP support" in (status.details or "")
    ]
    if unsupported_rows == 0 and not affected:
        return
    displayed = affected[:MAX_CONSOLE_GENE_ISSUES]
    outcomes = ", ".join(f"{status.query}={status.reason}" for status in displayed)
    extra = max(0, len(affected) - MAX_CONSOLE_GENE_ISSUES)
    message = (
        f"WARNING: Gene-list SNP-universe preflight found zero retained-SNP support for "
        f"{unsupported_rows} gene-list row(s)"
        + (f"; affected query outcomes: {outcomes}" if outcomes else "")
        + (f", plus {extra} more affected query(s)" if extra else "")
        + ". See diagnostics/query_annotation_status.tsv, diagnostics/gene_list_resolution_summary.tsv, "
        "and diagnostics/gene_list_audit.tsv.gz."
    )
    LOGGER.warning(message.removeprefix("WARNING: "))
    print(message, file=sys.stderr)


def _validate_run_args(args: argparse.Namespace) -> None:
    """Validate public LD-score workflow arguments before loading inputs.

    This validator intentionally lives in the workflow layer rather than in the
    kernel because optional baseline synthesis is a public orchestration rule:
    the numerical kernels still receive an explicit annotation bundle.
    """
    has_gene_lists = _has_cli_tokens(getattr(args, "query_annot_gene_list_sources", None))
    if not has_gene_lists and (
        getattr(args, "gene_exclude_regions", "none") != "none"
        or getattr(args, "control_gene_list_file", None) is not None
        or getattr(args, "gene_coordinate_file", None) is not None
        or getattr(args, "gene_list_resolution_policy", "strict") != "strict"
    ):
        raise LDSCUsageError(
            "Gene coordinate, exclusion, control-list, and resolution policy options are valid only with "
            "gene-list query annotations. Remove the gene-specific option or use gene-list mode."
        )
    if has_gene_lists and getattr(args, "gene_coordinate_file", None) is None:
        raise LDSCUsageError(
            "Live gene-list mode requires --gene-coordinate-file. Supply a one-build, one-based "
            "gene-coordinate TSV/TSV.GZ; there is no packaged catalog fallback."
        )
    if not _has_cli_tokens(args.baseline_annot_sources) and (
        _has_cli_tokens(args.query_annot_sources)
        or _has_cli_tokens(getattr(args, "query_annot_bed_sources", None))
        or _has_cli_tokens(getattr(args, "query_annot_gene_list_sources", None))
    ):
        raise LDSCUsageError(_QUERY_REQUIRES_BASELINE_MESSAGE)
    keep = getattr(args, "keep", None)
    if _uses_parquet_reference(args) == bool(args.bfile):
        raise LDSCUsageError(
            "ldscore could not choose a reference-panel backend. Most likely both parquet "
            "R2 input and PLINK input were supplied, or neither was supplied. Pass exactly "
            "one of `--r2-dir` for parquet mode or `--plink-prefix` for PLINK mode."
        )
    if _uses_parquet_reference(args):
        if keep:
            raise LDSCUsageError(
                "ldscore cannot apply `--keep-indivs-file` in parquet R2 mode. Most likely "
                "`--keep-indivs-file` was combined with `--r2-dir`, but individual-level "
                "filtering only exists before PLINK genotype LD calculation. Remove "
                "`--keep-indivs-file`, or rerun in PLINK mode with `--plink-prefix`."
            )
        if identity_mode_family(args.snp_identifier) == "chr_pos" and args.genome_build is None:
            raise LDSCUsageError(
                "ldscore cannot run parquet R2 mode with chr_pos-family SNP identifiers "
                "without a genome build. Most likely `--snp-identifier chr_pos` or an "
                "allele-aware chr_pos mode was supplied without `--genome-build`. Pass "
                "`--genome-build hg19`, `--genome-build hg38`, or `--genome-build auto`."
            )
    if args.ld_wind_cm is not None and args.ld_wind_cm <= 0:
        raise LDSCConfigError(
            f"ldscore received invalid `--ld-wind-cm={args.ld_wind_cm}`. Most likely "
            "the LD window was set to zero or a negative value. Pass a positive "
            "centimorgan window, for example `--ld-wind-cm 1.0`."
        )
    if args.ld_wind_kb is not None and args.ld_wind_kb <= 0:
        raise LDSCConfigError(
            f"ldscore received invalid `--ld-wind-kb={args.ld_wind_kb}`. Most likely "
            "the LD window was set to zero or a negative value. Pass a positive "
            "kilobase window."
        )
    if args.ld_wind_snps is not None and args.ld_wind_snps <= 0:
        raise LDSCConfigError(
            f"ldscore received invalid `--ld-wind-snps={args.ld_wind_snps}`. Most likely "
            "the LD window was set to zero or a negative SNP count. Pass a positive "
            "SNP-count window."
        )
    if getattr(args, "maf_min", None) is not None and not 0 <= args.maf_min <= 0.5:
        raise LDSCConfigError(
            f"ldscore received invalid `--maf-min={args.maf_min}`. Most likely the "
            "minor-allele frequency threshold was entered outside the valid [0, 0.5] "
            "range. Pass a value between 0 and 0.5, or omit the option."
        )
    if not 0 <= getattr(args, "common_maf_min", 0.05) <= 0.5:
        raise LDSCConfigError(
            f"ldscore received invalid `--common-maf-min={args.common_maf_min}`. Most "
            "likely the common-SNP MAF threshold was entered outside the valid [0, 0.5] "
            "range. Pass a value between 0 and 0.5."
        )
    if args.snp_batch_size <= 0:
        raise LDSCConfigError(
            f"ldscore received invalid `--snp-batch-size={args.snp_batch_size}`. Most "
            "likely the parquet query batch size was set to zero or a negative value. "
            "Pass a positive integer batch size."
        )


def _has_cli_tokens(value: str | Sequence[str] | None) -> bool:
    """Return whether a CLI path field contains at least one non-empty token."""
    return bool(split_cli_path_tokens(value))


def _format_ldscore_start_message(annotation_bundle, n_chromosomes: int) -> str:
    """Return the workflow-level LD-score start message for a bundle."""
    source_summary = getattr(annotation_bundle, "source_summary", {}) or {}
    baseline_source = str(source_summary.get("baseline", ""))
    if (
        list(getattr(annotation_bundle, "baseline_columns", [])) == ["base"]
        and not list(getattr(annotation_bundle, "query_columns", []))
        and baseline_source.startswith("synthetic all-ones base annotation")
    ):
        return (
            f"Starting LD-score calculation for {n_chromosomes} chromosomes "
            "with synthetic base annotation and no query annotations."
        )
    return (
        f"Starting LD-score calculation for {n_chromosomes} chromosomes "
        f"with {len(annotation_bundle.baseline_columns)} baseline columns "
        f"and {len(annotation_bundle.query_columns)} query columns."
    )


def _uses_parquet_reference(args: argparse.Namespace) -> bool:
    """Return whether normalized args select parquet reference-panel mode."""
    return bool(_r2_dir_from_args(args))


def _r2_dir_from_args(args: argparse.Namespace) -> str | None:
    """Return the R2 directory supplied to the LD-score workflow."""
    return normalize_optional_path_token(getattr(args, "r2_dir", None))


def _pseudo_base_annotation_bundle_from_ref_panel(ref_panel, global_config: GlobalConfig):
    """Build an all-ones ``base`` bundle from retained reference-panel metadata.

    The reference-panel adapter has already applied retained-panel SNP filters
    when ``load_metadata(chrom)`` returns. Later runtime regression-SNP
    restriction remains in the normal LD-score compute path.
    """
    from .annotation_builder import AnnotationBundle

    metadata_frames = []
    chromosomes = [str(chrom) for chrom in ref_panel.available_chromosomes()]
    for chrom in chromosomes:
        metadata = ref_panel.load_metadata(chrom).copy()
        if len(metadata) == 0:
            continue
        if "POS" not in metadata.columns and "BP" in metadata.columns:
            metadata = metadata.rename(columns={"BP": "POS"})
        metadata_columns = ["CHR", "SNP", "CM", "POS", *[column for column in ("A1", "A2") if column in metadata.columns]]
        metadata_frames.append(metadata.loc[:, metadata_columns].reset_index(drop=True))
    if not metadata_frames:
        raise LDSCInputError(
            "ldscore could not build the synthetic `base` annotation from the reference panel: "
            "no retained reference-panel SNP metadata rows were available. Most likely the "
            "reference panel was filtered to zero SNPs by `--ref-panel-snps-file` or the "
            "selected chromosomes are absent. Use a reference panel and SNP restriction file "
            "with overlapping SNPs, or pass explicit baseline annotations."
        )
    metadata = pd.concat(metadata_frames, axis=0, ignore_index=True)
    baseline = pd.DataFrame({"base": np.ones(len(metadata), dtype=np.float32)})
    query = pd.DataFrame(index=metadata.index)
    bundle = AnnotationBundle(
        metadata=metadata,
        baseline_annotations=baseline,
        query_annotations=query,
        baseline_columns=["base"],
        query_columns=[],
        chromosomes=chromosomes,
        source_summary={"baseline": "synthetic all-ones base annotation from retained reference-panel metadata"},
        config_snapshot=global_config,
    )
    bundle.validate()
    return bundle


def run_ldscore(**kwargs) -> LDScoreResult:
    """Run LD-score calculation from Python using public CLI-style names.

    Keyword arguments are interpreted as CLI-equivalent option names without
    leading ``--``; for example ``baseline_annot_sources``, ``query_annot_sources``,
    ``query_annot_bed_sources``, ``query_annot_gene_list_sources``, ``plink_prefix``, ``r2_dir``,
    ``keep_indivs_file``, ``snp_batch_size``, ``common_maf_min``, and
    ``output_dir``. Shared runtime assumptions such as ``snp_identifier`` and
    ``genome_build`` must be supplied through ``set_global_config(...)`` first,
    while per-run controls such as ``ref_panel_snps_file`` and
    ``regr_snps_file`` remain ordinary keyword arguments here. ``output_dir``
    is required because this public workflow always writes canonical artifacts.

    When ``baseline_annot_sources`` and query inputs are omitted, the workflow
    builds a synthetic all-ones baseline column named ``base`` from retained
    reference-panel metadata. ``query_annot_sources`` and
    ``query_annot_bed_sources`` and ``query_annot_gene_list_sources`` require
    explicit baseline annotations because query columns are interpreted
    relative to that baseline SNP universe. Gene-list inputs are projected
    from the required one-based coordinate catalog for the resolved build.
    ``padding_bp`` defaults to zero for live BED input but must be chosen
    explicitly for live gene lists; passing it with a prebuilt annotation
    query, no query, or an exact index raises ``LDSCUsageError`` even when its
    value is zero.

    Returns
    -------
    LDScoreResult
        Aggregated result with ``baseline_table``, optional ``query_table``,
        count records, and canonical output paths.
    """
    forbidden = sorted({"snp_identifier", "genome_build", "log_level"} & set(kwargs))
    if forbidden:
        joined = ", ".join(forbidden)
        raise LDSCUsageError(
            f"Python run_ldscore() cannot accept shared runtime option(s): {joined}. "
            "Most likely this call still passes pre-restructure keyword arguments for "
            "SNP identity, genome build, or logging. Call set_global_config(...) first, "
            "then pass only LD-score run-specific options to run_ldscore()."
        )
    removed = sorted(
        {
            "out",
            "baseline_annot",
            "query_annot",
            "query_annot_bed",
            "bfile",
            "r2_table",
            "frqfile",
            "keep",
            "baseline_annot_paths",
            "query_annot_paths",
            "query_annot_bed_paths",
            "plink_path",
            "r2_paths",
            "metadata_paths",
            "r2_ref_panel_dir",
            "ref_panel_dir",
            "r2_sources",
            "metadata_sources",
            "ref_panel_snps_path",
            "regression_snps_path",
            "keep_indivs_path",
            "chunk_size",
            "maf",
            "control_gene_list_source",
        }
        & set(kwargs)
    )
    if removed:
        joined = ", ".join(removed)
        raise LDSCUsageError(
            f"Python run_ldscore() cannot accept removed IO argument(s): {joined}. "
            "Most likely this call still uses legacy keyword names from the old LDSC API. "
            "Use CLI-style names such as `baseline_annot_sources`, `query_annot_sources`, "
            "`plink_prefix`, `r2_dir`, and `output_dir`."
        )
    if not kwargs.get("output_dir"):
        raise LDSCUsageError(
            "run_ldscore() requires output_dir. Pass an explicit directory for the canonical LD-score artifacts."
        )
    parser = build_parser()
    defaults = {
        action.dest: action.default
        for action in parser._actions
        if action.dest != "help" and action.default is not argparse.SUPPRESS
    }
    global_config = get_global_config()
    defaults["snp_identifier"] = global_config.snp_identifier
    defaults["genome_build"] = global_config.genome_build
    defaults["log_level"] = global_config.log_level
    defaults.update(kwargs)
    args = argparse.Namespace(**defaults)
    args._explicit_cli_options = frozenset(
        "--" + key.replace("_", "-") for key in kwargs
    )
    args._emit_gene_console_notices = False
    return run_ldscore_from_args(args)


def main(argv: Sequence[str] | None = None) -> LDScoreResult:
    """Command-line entry point for the LD-score workflow."""
    parser = build_parser()
    argv_list = list(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(argv_list)
    args._explicit_cli_options = frozenset(
        token.split("=", 1)[0] for token in argv_list if token.startswith("--")
    )
    args._emit_gene_console_notices = True
    return run_ldscore_from_args(args)


def _normalize_run_args(
    args: argparse.Namespace,
    *,
    gene_catalog: Any | None = None,
) -> tuple[argparse.Namespace, GlobalConfig]:
    """Normalize CLI-style args and derive the shared ``GlobalConfig`` object."""
    normalized_mode = normalize_snp_identifier_mode(args.snp_identifier)
    normalized_args = argparse.Namespace(**vars(args))
    for attr in ("query_annot_chr", "baseline_annot_chr", "bfile_chr", "r2_table_chr", "frqfile_chr"):
        if not hasattr(normalized_args, attr):
            setattr(normalized_args, attr, None)
    for attr in (
        "query_annot_sources",
        "baseline_annot_sources",
        "plink_prefix",
        "r2_dir",
        "query_annot_bed_sources",
        "query_annot_gene_list_sources",
        "gene_coordinate_file",
        "keep_indivs_file",
    ):
        if not hasattr(normalized_args, attr):
            setattr(normalized_args, attr, None)
    for attr in ("ref_panel_snps_file", "regr_snps_file"):
        if not hasattr(normalized_args, attr):
            setattr(normalized_args, attr, None)
    for attr in ("regr_snps_exclude_regions",):
        if not hasattr(normalized_args, attr):
            setattr(normalized_args, attr, None)
    if not hasattr(normalized_args, "control_gene_list_file"):
        normalized_args.control_gene_list_file = None
    if not hasattr(normalized_args, "gene_list_resolution_policy"):
        normalized_args.gene_list_resolution_policy = "strict"
    if not hasattr(normalized_args, "gene_exclude_regions"):
        normalized_args.gene_exclude_regions = "none"
    if not hasattr(normalized_args, "maf_min"):
        normalized_args.maf_min = None
    if not hasattr(normalized_args, "common_maf_min"):
        normalized_args.common_maf_min = 0.05
    if not hasattr(normalized_args, "padding_bp") or normalized_args.padding_bp is None:
        normalized_args.padding_bp = 0
    if not hasattr(normalized_args, "snp_batch_size"):
        normalized_args.snp_batch_size = 128
    normalized_args.snp_identifier = normalized_mode
    normalized_args.output_dir = normalize_path_token(args.output_dir)
    normalized_args.r2_dir = _r2_dir_from_args(normalized_args)
    normalized_args.keep_indivs_file = normalize_optional_path_token(getattr(args, "keep_indivs_file", None))
    normalized_args.gene_coordinate_file = normalize_optional_path_token(getattr(args, "gene_coordinate_file", None))
    normalized_args.ref_panel_snps_file = normalize_optional_path_token(getattr(args, "ref_panel_snps_file", None))
    normalized_args.regr_snps_file = normalize_optional_path_token(getattr(args, "regr_snps_file", None))
    normalized_args.regr_snps_exclude_regions = getattr(args, "regr_snps_exclude_regions", None)
    # The numerical kernel still consumes the historical namespace shape.
    normalized_args.query_annot = normalized_args.query_annot_sources
    normalized_args.baseline_annot = normalized_args.baseline_annot_sources
    normalized_args.bfile = normalized_args.plink_prefix
    normalized_args.r2_table = None
    normalized_args.frqfile = None
    normalized_args.keep = normalized_args.keep_indivs_file
    has_gene_lists = _has_cli_tokens(normalized_args.query_annot_gene_list_sources)
    normalized_args.gene_catalog_build = None
    if identity_mode_family(normalized_mode) == "rsid":
        requested_build = normalize_genome_build(getattr(args, "genome_build", None))
        if has_gene_lists:
            if requested_build in {None, "auto"}:
                normalized_args.gene_catalog_build = _resolve_ldscore_chr_pos_genome_build(
                    normalized_args,
                    "auto",
                    gene_catalog_build=getattr(gene_catalog, "genome_build", None),
                )
            else:
                normalized_args.gene_catalog_build = requested_build
            global_config = GlobalConfig(
                snp_identifier=normalized_mode,
                genome_build=None,
                log_level=getattr(args, "log_level", "INFO"),
            )
            normalized_args.genome_build = None
        else:
            global_config = GlobalConfig(
                snp_identifier=normalized_mode,
                genome_build=requested_build,
                log_level=getattr(args, "log_level", "INFO"),
            )
            normalized_args.genome_build = global_config.genome_build
    else:
        requested_build = getattr(args, "genome_build", None)
        if has_gene_lists and requested_build is None:
            requested_build = "auto"
        resolved_genome_build = _resolve_ldscore_chr_pos_genome_build(
            normalized_args,
            requested_build,
            gene_catalog_build=getattr(gene_catalog, "genome_build", None),
        )
        normalized_args.genome_build = resolved_genome_build
        if has_gene_lists:
            normalized_args.gene_catalog_build = resolved_genome_build
        global_config = GlobalConfig(
            snp_identifier=normalized_mode,
            genome_build=resolved_genome_build,
            log_level=getattr(args, "log_level", "INFO"),
        )
    return normalized_args, global_config


def _resolve_ldscore_chr_pos_genome_build(
    args: argparse.Namespace,
    genome_build: str | None,
    *,
    gene_catalog_build: str | None = None,
) -> str:
    normalized = normalize_genome_build(genome_build)
    if normalized is None:
        raise LDSCUsageError(
            "ldscore cannot resolve chr_pos-family SNP identifiers without a genome build. "
            "Most likely `--snp-identifier` was set to a chr_pos mode but `--genome-build` "
            "was omitted. Pass `--genome-build auto`, `--genome-build hg19`, or "
            "`--genome-build hg38`."
        )
    if normalized != "auto":
        return normalized

    resolved: list[tuple[str, str]] = []
    normalized_catalog_build = normalize_genome_build(gene_catalog_build)
    if normalized_catalog_build in {"hg19", "hg38"}:
        resolved.append(("gene-coordinate catalog", normalized_catalog_build))
    annotation_tokens = split_cli_path_tokens(getattr(args, "baseline_annot_sources", None))
    if annotation_tokens:
        frame, sampled_path = sample_frame_from_chr_pattern(
            annotation_tokens,
            context="LD-score annotation inputs",
        )
        resolved.append(
            (
                "annotation",
                resolve_genome_build(
                    "auto",
                    "chr_pos",
                    frame,
                    context="LD-score annotation inputs",
                    logger=LOGGER,
                ),
            )
        )
        LOGGER.info(f"Resolved LD-score annotation genome build from '{sampled_path}'.")
    r2_dir = _r2_dir_from_args(args)
    if r2_dir is not None:
        ref_panel_build = _infer_r2_dir_genome_build(r2_dir)
        if ref_panel_build is not None:
            resolved.append(("reference panel", ref_panel_build))
            LOGGER.info(f"Resolved LD-score reference-panel genome build from parquet schema metadata in '{r2_dir}'.")
    if not resolved:
        raise LDSCInputError(
            "ldscore could not infer the genome build for chr_pos inputs. Most likely "
            "no chromosome-suite annotation sample or R2 parquet build metadata was "
            "available to inspect. Pass `--genome-build hg19` or `--genome-build hg38`, "
            "or provide LDSC-generated parquet R2 files with build metadata. "
            f"Other causes & fixes: {_LDSCORE_BUILD_DOC}"
        )
    builds = {build for _label, build in resolved}
    if len(builds) != 1:
        details = ", ".join(f"{label}={build}" for label, build in resolved)
        raise LDSCInputError(
            f"ldscore found conflicting genome-build evidence: {details}. Most likely "
            "the annotation files and parquet reference panel were generated on different "
            "builds. Regenerate one input on the same build, or rerun with a matching "
            "reference panel and `--genome-build`. "
            f"Other causes & fixes: {_LDSCORE_BUILD_DOC}"
        )
    return resolved[0][1]


def _load_regression_snps(
    path: str | None,
    global_config: GlobalConfig,
    *,
    label: str = "regression SNP list",
) -> RestrictionIdentityKeys | None:
    """Load ``LDScoreConfig.regr_snps_file`` using the active identifier mode."""
    if not path:
        return None
    return read_snp_restriction_keys(
        resolve_scalar_path(path, label=label),
        global_config.snp_identifier,
        genome_build=global_config.genome_build,
    )


def _regr_snps_file_from_config(config: LDScoreConfig) -> str:
    """Return the explicit regression list or the canonical bundled HM3 default."""
    return config.regr_snps_file or packaged_hm3_curated_map_path()


def _infer_r2_dir_genome_build(r2_dir: str) -> str | None:
    """Infer hg19/hg38 from R2 parquet schema metadata, never path names."""
    builds_by_path: dict[str, str] = {}
    for path in _candidate_r2_schema_paths(r2_dir):
        build = _read_r2_sorted_by_build(path)
        if build is not None:
            builds_by_path[str(path)] = build

    builds = set(builds_by_path.values())
    if len(builds) > 1:
        details = ", ".join(f"{path}={build}" for path, build in sorted(builds_by_path.items()))
        raise LDSCInputError(
            f"ldscore found conflicting R2 parquet genome-build metadata in '{r2_dir}': "
            f"{details}. Most likely the R2 directory mixes files from different reference-panel "
            "builds. Keep only one build in the R2 directory or regenerate the panel as a single "
            "consistent artifact. "
            f"Other causes & fixes: {_LDSCORE_BUILD_DOC}"
        )
    return next(iter(builds), None)


def _candidate_r2_schema_paths(r2_dir: str) -> list[Path]:
    """Return candidate R2 parquet files whose schema metadata can identify a panel build."""
    root = Path(r2_dir)
    if not root.is_dir():
        return []

    direct = sorted(root.glob("chr*_r2.parquet"))
    if direct:
        return direct

    paths: list[Path] = []
    for build in ("hg19", "hg38"):
        child = root / build
        if child.is_dir():
            paths.extend(sorted(child.glob("chr*_r2.parquet")))
    return paths


def _read_r2_sorted_by_build(path: Path) -> str | None:
    """Read ``ldsc:sorted_by_build`` from one R2 parquet schema."""
    import pyarrow.parquet as pq

    raw_meta = pq.read_schema(str(path)).metadata or {}
    raw_build = raw_meta.get(b"ldsc:sorted_by_build")
    if raw_build is None:
        return None
    build = normalize_genome_build(raw_build.decode("utf-8"))
    return build if build in {"hg19", "hg38"} else None


def _resolve_regression_region_build(
    args: argparse.Namespace, global_config: GlobalConfig, presets: tuple[str, ...]
) -> str | None:
    """Resolve the genome build used to select region-exclusion preset BEDs.

    Coordinate-family runs reuse the resolved panel build. Gene-list runs may
    reuse their concrete inferred projection build while keeping rsID identity
    metadata build-independent. Other rsID-family runs need an explicit build
    because rsIDs themselves carry no coordinate-build information.
    """
    if not presets:
        return None
    if global_config.genome_build in {"hg19", "hg38"}:
        return global_config.genome_build
    requested = normalize_genome_build(getattr(args, "genome_build", None))
    if requested in {"hg19", "hg38"}:
        return requested
    gene_catalog_build = normalize_genome_build(getattr(args, "gene_catalog_build", None))
    if gene_catalog_build in {"hg19", "hg38"}:
        return gene_catalog_build
    raise LDSCUsageError(
        "ldscore cannot select named regression exclusion regions without a concrete genome build. "
        "Pass `--genome-build hg19` or `--genome-build hg38`, or use `--regr-snps-exclude-regions none`."
    )


def _regression_region_intervals(
    args: argparse.Namespace, global_config: GlobalConfig
) -> kernel_regions.RegionIntervals | None:
    """Load named intervals applied only after selecting regression SNPs."""
    presets = regr_snps_exclude_regions_choice_to_presets(getattr(args, "regr_snps_exclude_regions", None) or "none")
    build = _resolve_regression_region_build(args, global_config, presets)
    if not presets:
        return None
    assert build is not None
    return kernel_regions.load_preset_intervals(presets, build)


def _ref_panel_from_args(args: argparse.Namespace, global_config: GlobalConfig, *, chromosome_prefixes=None):
    """Build the reference-panel adapter that owns the ``A -> A'`` restriction."""
    from ._kernel.ref_panel import RefPanelLoader

    ref_panel_snps_file = normalize_optional_path_token(getattr(args, "ref_panel_snps_file", None))
    r2_dir = _r2_dir_from_args(args)
    if r2_dir is not None:
        spec = RefPanelConfig(
            backend="parquet_r2",
            r2_dir=r2_dir,
            ref_panel_snps_file=ref_panel_snps_file,
            maf_min=getattr(args, "maf_min", None),
            keep_indivs_file=getattr(args, "keep_indivs_file", None),
            genetic_map_hg19_sources=getattr(args, "genetic_map_hg19_sources", None),
            genetic_map_hg38_sources=getattr(args, "genetic_map_hg38_sources", None),
        )
    else:
        spec = RefPanelConfig(
            backend="plink",
            plink_prefix=getattr(args, "plink_prefix", None),
            ref_panel_snps_file=ref_panel_snps_file,
            maf_min=getattr(args, "maf_min", None),
            keep_indivs_file=getattr(args, "keep_indivs_file", None),
            genetic_map_hg19_sources=getattr(args, "genetic_map_hg19_sources", None),
            genetic_map_hg38_sources=getattr(args, "genetic_map_hg38_sources", None),
        )
    if spec.backend == "plink" and chromosome_prefixes is not None:
        from ._kernel.ref_panel import PlinkRefPanel

        return PlinkRefPanel(global_config, spec, chromosome_prefixes=chromosome_prefixes)
    return RefPanelLoader(global_config).load(spec)


def _ldscore_config_from_args(args: argparse.Namespace) -> LDScoreConfig:
    """Build an ``LDScoreConfig`` from normalized LD-score arguments."""
    return LDScoreConfig(
        ld_wind_snps=getattr(args, "ld_wind_snps", None),
        ld_wind_kb=getattr(args, "ld_wind_kb", None),
        ld_wind_cm=getattr(args, "ld_wind_cm", None),
        regr_snps_file=getattr(args, "regr_snps_file", None),
        snp_batch_size=getattr(args, "snp_batch_size", 128),
        query_batch_size=getattr(args, "query_batch_size", 1000),
        common_maf_min=getattr(args, "common_maf_min", 0.05),
        whole_chromosome_ok=getattr(args, "yes_really", False),
        export_ref_metadata=getattr(args, "export_ref_metadata", False),
        threads=getattr(args, "threads", 1),
    )


def _reference_metadata_sources(ref_panel) -> tuple[str, str]:
    """Return ``(cm_source, maf_source)`` provenance labels for the reference panel."""
    spec = getattr(ref_panel, "spec", None)
    backend = getattr(spec, "backend", None)
    if backend == "parquet_r2":
        return "parquet_sidecar", "parquet_sidecar"
    has_map = bool(getattr(spec, "genetic_map_hg19_sources", None) or getattr(spec, "genetic_map_hg38_sources", None))
    return ("genetic_map" if has_map else "plink_bim"), "plink_genotypes"


def _snp_universe_policy(
    *,
    ref_panel,
    ldscore_config: LDScoreConfig,
    regression_regions: kernel_regions.RegionIntervals | None,
    chromosome_results: Sequence[ChromLDScoreResult],
) -> dict[str, Any]:
    """Describe the independently selected LD-reference and regression universes.

    This is persisted in ``metadata.json`` so an LD-score artifact cannot make
    a regression-row exclusion look like a change to its LD-reference
    estimand. Counts are after retained-reference/annotation alignment, before
    and after the regression-region subtraction respectively.
    """
    ref_panel_snps_file = getattr(getattr(ref_panel, "spec", None), "ref_panel_snps_file", None)
    custom_regression_file = ldscore_config.regr_snps_file
    return {
        "ld_reference_universe": {
            "selection": "full_retained_reference_panel"
            if ref_panel_snps_file is None
            else "explicit_ref_panel_snps_file",
            "ref_panel_snps_file": None if ref_panel_snps_file is None else str(ref_panel_snps_file),
            "retained_snp_count": int(sum(result.reference_snp_count for result in chromosome_results)),
        },
        "regression_rows_and_weights": {
            "selection": "bundled_hm3_default" if custom_regression_file is None else "explicit_regr_snps_file",
            "regr_snps_file": None if custom_regression_file is None else str(custom_regression_file),
            "region_exclusion_sources": [] if regression_regions is None else list(regression_regions.source_labels),
            "selected_snp_count_before_region_exclusion": int(
                sum(result.regression_selected_snp_count for result in chromosome_results)
            ),
            "region_excluded_snp_count": int(
                sum(result.regression_region_removed_snp_count for result in chromosome_results)
            ),
            "written_snp_count": int(sum(len(result.baseline_table) for result in chromosome_results)),
            "weight_contributors": "same_filtered_regression_set",
        },
    }


def _write_one_ref_metadata_sidecar(metadata: pd.DataFrame, chrom: str, output_dir: str) -> None:
    """Write one opt-in reference-metadata sidecar (PLINK backend).

    The format matches the parquet panel sidecar: a gzip TSV with
    ``CHR POS SNP A1 A2 CM MAF`` columns, so an exported sidecar can seed a future
    ``build-ref-panel`` or a QC diff. Written from inside the per-chromosome
    worker so the lean cross-process result does not carry per-SNP metadata.
    """
    columns = [c for c in ("CHR", "POS", "SNP", "A1", "A2", "CM", "MAF") if c in metadata.columns]
    export_dir = Path(output_dir) / "ref_metadata"
    export_dir.mkdir(parents=True, exist_ok=True)
    path = export_dir / f"chr{chrom}_meta.tsv.gz"
    metadata[columns].to_csv(path, sep="\t", index=False, na_rep="NA", float_format="%.6g", compression="gzip")
    LOGGER.info(f"Wrote reference-metadata sidecar '{path}'.")


def _output_config_from_args(args: argparse.Namespace) -> LDScoreOutputConfig:
    """Translate LD-score CLI arguments into the canonical directory output config."""
    return LDScoreOutputConfig(
        output_dir=normalize_path_token(args.output_dir),
        overwrite=getattr(args, "overwrite", False),
    )


def _replace_result_output_paths(result: LDScoreResult, output_paths: dict[str, str]) -> LDScoreResult:
    """Return ``result`` with updated artifact-path metadata after writing outputs.

    Uses ``dataclasses.replace`` so every other field (including ``overlap``) is
    carried over; a manual reconstruction silently drops fields added later.
    """
    return dataclass_replace(result, output_paths=dict(output_paths))


def _available_cpu_count() -> int:
    """Return the number of CPUs available to this process.

    Prefers ``os.sched_getaffinity`` (Linux) so SLURM/cgroup/cpuset/Docker CPU
    allocations are respected rather than the whole machine's core count; falls
    back to ``os.cpu_count()`` on platforms without affinity support.
    """
    getaffinity = getattr(os, "sched_getaffinity", None)
    if getaffinity is not None:
        try:
            return len(getaffinity(0)) or 1
        except OSError:
            pass
    return os.cpu_count() or 1


def _resolve_worker_count(threads: int, n_chromosomes: int) -> int:
    """Resolve the effective worker-process count from the ``threads`` setting.

    Uses the joblib ``n_jobs`` convention: ``1`` is sequential, a positive ``N``
    requests ``N`` workers, ``-1`` requests all available cores, ``-2`` all but
    one, and any negative ``-k`` requests ``n_cpus + 1 - k``. Core counts respect
    CPU affinity (see :func:`_available_cpu_count`). The result is capped at
    ``n_chromosomes`` and floored at ``1`` so a single chromosome never spawns a
    pool.
    """
    if n_chromosomes <= 0:
        return 1
    if threads < 0:
        resolved = _available_cpu_count() + 1 + threads
    else:
        resolved = threads
    return max(1, min(resolved, n_chromosomes))


def _chromosomes_from_bundle(annotation_bundle) -> list[str]:
    """Infer the chromosome processing order from an annotation bundle."""
    chromosomes = getattr(annotation_bundle, "chromosomes", None)
    if chromosomes:
        return list(chromosomes)
    return sorted(annotation_bundle.metadata["CHR"].astype(str).unique().tolist())


def _slice_annotation_bundle(annotation_bundle, chrom: str):
    """Return the per-chromosome view of an annotation bundle."""
    if hasattr(annotation_bundle,'shards'):
        return dataclass_replace(annotation_bundle,shards={str(chrom):annotation_bundle.shard(chrom)},gene_list_batch=None)
    keep = annotation_bundle.metadata["CHR"].astype(str) == str(chrom)
    return type(annotation_bundle)(
        metadata=annotation_bundle.metadata.loc[keep].reset_index(drop=True),
        baseline_annotations=annotation_bundle.baseline_annotations.loc[keep].reset_index(drop=True),
        query_annotations=annotation_bundle.query_annotations.loc[keep].reset_index(drop=True),
        baseline_columns=list(annotation_bundle.baseline_columns),
        query_columns=list(annotation_bundle.query_columns),
        chromosomes=[str(chrom)],
        source_summary=dict(getattr(annotation_bundle, "source_summary", {})),
        config_snapshot=getattr(annotation_bundle, "config_snapshot", None),
        query_statuses=tuple(getattr(annotation_bundle, "query_statuses", ())),
        gene_list_batch=getattr(annotation_bundle, "gene_list_batch", None),
    )


def _kernel_annotation_bundle(bundle, chrom):
    """Borrow one chromosome's metadata and selected annotation reader."""
    if hasattr(bundle,'shards'):
        from ._kernel.ldscore_projection import MappedAnnotations
        shard = bundle.shard(chrom)
        metadata = shard.metadata()
        values = MappedAnnotations(shard,np.arange(shard.n_rows),tuple(bundle.baseline_columns+bundle.query_columns))
    else:
        metadata,values = bundle.metadata.copy(),_float32_annotation_frame(bundle)
    return kernel_ldscore.AnnotationBundle(metadata,values,list(bundle.baseline_columns),list(bundle.query_columns))


def _float32_annotation_frame(annotation_bundle) -> pd.DataFrame:
    """Convert one chromosome's annotation blocks to one float32 numerical frame.

    The destination is allocated once. Contiguous same-dtype column runs are
    copied with vectorized NumPy casting, avoiding a second full mixed-type
    staging matrix.
    """
    frames = (annotation_bundle.baseline_annotations, annotation_bundle.query_annotations)
    columns = [*annotation_bundle.baseline_columns, *annotation_bundle.query_columns]
    values = np.empty((len(annotation_bundle.metadata), len(columns)), dtype=np.float32)
    destination_start = 0
    for frame in frames:
        source_start = 0
        dtypes = frame.dtypes.tolist()
        while source_start < len(dtypes):
            source_end = source_start + 1
            while source_end < len(dtypes) and dtypes[source_end] == dtypes[source_start]:
                source_end += 1
            values[
                :,
                destination_start + source_start : destination_start + source_end,
            ] = frame.iloc[:, source_start:source_end].to_numpy(copy=False)
            source_start = source_end
        destination_start += frame.shape[1]
    return pd.DataFrame(values, columns=columns, copy=False)


def _is_empty_intersection(error: Exception, chrom: str) -> bool:
    """Return whether ``error`` is the recoverable empty-intersection case."""
    from .errors import EmptyReferenceSNPs

    if isinstance(error, EmptyReferenceSNPs):
        return True
    message = str(error)
    old_shape = message.startswith("No retained annotation SNPs remain on chromosome ")
    new_shape = message.startswith("ldscore retained no annotation SNPs on chromosome ")
    if not (old_shape or new_shape):
        return False
    if old_shape and not any(suffix in message for suffix in (" after parquet intersection.", " after PLINK intersection.")):
        return False
    if new_shape and "after intersecting with the " not in message:
        return False
    return True


def _warn_and_skip_empty_intersection(error: Exception, chrom: str) -> bool:
    """Warn and signal skip when a chromosome loses all SNPs after reference intersection."""
    if not _is_empty_intersection(error, chrom):
        return False
    warnings.warn(f"Skipping chromosome {chrom}: {error}", UserWarning, stacklevel=3)
    return True


# Worker state shared within a pool worker process. ``regression_snps`` is set
# once by the pool initializer so a large SNP set is transferred per worker, not
# per task; the inline (single-worker) caller passes it explicitly instead.
_WORKER_UNSET = "__use_worker_global__"
_WORKER_STATE: dict[str, Any] = {}


@dataclass(frozen=True)
class _ChromOutcome:
    """Tagged result of one chromosome's LD-score computation in a worker."""

    chrom: str
    result: "ChromLDScoreResult | None"
    skipped: bool
    skip_message: str | None = None


def _compute_one_chromosome(
    chrom: str,
    chrom_bundle,
    ref_panel_spec,
    ldscore_config: LDScoreConfig,
    global_config: GlobalConfig,
    export_dir: str | None = None,
    regression_snps=_WORKER_UNSET,
    regression_regions=_WORKER_UNSET,
) -> _ChromOutcome:
    """Compute one chromosome end-to-end and return a tagged outcome.

    Rebuilds the reference-panel adapter locally from ``ref_panel_spec`` so no
    live reader/file handle crosses the process boundary. ``regression_snps``
    defaults to the pool-initializer-provided global; the inline single-worker
    caller passes it explicitly. Recoverable empty-intersection errors return a
    skipped outcome rather than raising across the pool boundary.
    """
    from ._kernel.ref_panel import RefPanelLoader

    if regression_snps == _WORKER_UNSET:
        regression_snps = _WORKER_STATE.get("regression_snps")
    if regression_regions == _WORKER_UNSET:
        regression_regions = _WORKER_STATE.get("regression_regions")
    prefixes = chrom_bundle.source_summary.get('chromosome_scope',{}).get('reference_prefixes_by_chrom')
    if ref_panel_spec.backend == 'plink' and prefixes:
        from ._kernel.ref_panel import PlinkRefPanel
        ref_panel = PlinkRefPanel(global_config,ref_panel_spec,chromosome_prefixes=prefixes)
    else:
        ref_panel = RefPanelLoader(global_config).load(ref_panel_spec)
    calculator = LDScoreCalculator()
    try:
        result = calculator.compute_chromosome(
            chrom=chrom,
            annotation_bundle=chrom_bundle,
            ref_panel=ref_panel,
            ldscore_config=ldscore_config,
            global_config=global_config,
            regression_snps=regression_snps,
            regression_regions=regression_regions,
            export_dir=export_dir,
        )
    except (ValueError, LDSCInputError) as exc:
        if _is_empty_intersection(exc, chrom):
            return _ChromOutcome(chrom=chrom, result=None, skipped=True, skip_message=str(exc))
        raise
    return _ChromOutcome(chrom=chrom, result=result, skipped=False)


def _init_worker(regression_snps, regression_regions=None, log_level: str = "INFO") -> None:
    """Initialize a pool worker: shared regression keys, logging, BLAS threads.

    Pins BLAS thread counts to 1 unless the user already set them, so ``W``
    worker processes do not oversubscribe ``W x BLAS_threads`` cores.
    """
    from ._logging import configure_package_logging

    _WORKER_STATE["regression_snps"] = regression_snps
    _WORKER_STATE["regression_regions"] = regression_regions
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(var, "1")
    configure_package_logging(log_level)
