"""Workflow-layer LD-score orchestration and normalized result objects.

This module is the public boundary for LD-score path handling and result
normalization. Callers may provide exact paths, glob patterns, or explicit
chromosome-suite tokens using ``@``. The workflow resolves those tokens into
concrete per-chromosome files, aligns each chromosome-local annotation bundle
to the prepared reference-panel metadata returned by ``ref_panel.load_metadata``,
and only then dispatches to the primitive-only kernel.

For ordinary unpartitioned LD-score runs, callers may omit both baseline and
query annotations. In that case the workflow constructs an all-ones baseline
annotation named exactly ``base`` over the retained reference-panel metadata and
continues through the same calculator and output writer used by partitioned
runs. Query annotations remain partitioned-LDSC inputs: they are accepted only
when explicit baseline annotations are supplied.

Parsed workflow entry points write ``ldscore.log`` under ``output_dir`` after
preflighting all deterministic outputs. Direct ``LDScoreCalculator.run(...)``
calls remain data-oriented and do not create log files.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import logging
from pathlib import Path
from typing import Any, Sequence
import warnings

import numpy as np
import pandas as pd

from ._chr_sampler import sample_frame_from_chr_pattern
from .column_inference import normalize_genome_build, normalize_snp_identifier_mode
from .config import (
    ConfigMismatchError,
    GlobalConfig,
    AnnotationBuildConfig,
    LDScoreConfig,
    RefPanelConfig,
    get_global_config,
    print_global_config_banner,
    validate_config_compatibility,
)
from .genome_build_inference import resolve_genome_build
from .outputs import LDScoreDirectoryWriter, LDScoreOutputConfig
from .path_resolution import (
    ensure_output_directory,
    ensure_output_paths_available,
    normalize_optional_path_token,
    normalize_path_token,
    resolve_plink_prefix,
    resolve_scalar_path,
    split_cli_path_tokens,
)
from ._logging import log_inputs, log_outputs, workflow_logging
from ._kernel import ldscore as kernel_ldscore
from ._kernel.identifiers import build_snp_id_series, read_global_snp_restriction
from ._row_alignment import assert_same_snp_rows


LOGGER = logging.getLogger("LDSC.ldscore_calculator")
_LDSCORE_SUFFIX_COLUMNS = ("CHR", "SNP", "POS", "BP", "CM", "MAF")
_QUERY_REQUIRES_BASELINE_MESSAGE = (
    "Query annotations require baseline annotations. If you intentionally want to test query annotations "
    "against an all-ones baseline, create an explicit all-ones `base` baseline annotation over the query "
    "annotation universe, then run the partitioned LDSC workflow with both baseline and query annotations."
)


@dataclass(frozen=True)
class _LegacyChromResult:
    """Compatibility record for chromosome outputs produced by the kernel."""
    chrom: str
    metadata: pd.DataFrame
    ld_scores: np.ndarray
    w_ld: np.ndarray
    M: np.ndarray
    M_5_50: np.ndarray | None
    ldscore_columns: list[str]
    baseline_columns: list[str]
    query_columns: list[str]


@dataclass(frozen=True)
class ChromLDScoreResult:
    """Normalized per-chromosome LD-score output in split baseline/query form.

    Parameters
    ----------
    chrom : str
        Chromosome label.
    baseline_table : pandas.DataFrame
        Table with ``CHR``, ``SNP``, ``POS``, ``regr_weight``, and baseline
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

    def validate(self) -> None:
        """Check the normalized public contract for chromosome-level results."""
        required = {"CHR", "SNP", "POS", "regr_weight", *self.baseline_columns}
        missing = required - set(self.baseline_table.columns)
        if missing:
            raise ValueError(f"baseline_table is missing required columns: {sorted(missing)}")
        if self.query_columns and self.query_table is None:
            raise ValueError("query_table is required when query_columns are present.")
        if self.query_table is not None:
            missing_query = {"CHR", "SNP", "POS", *self.query_columns} - set(self.query_table.columns)
            if missing_query:
                raise ValueError(f"query_table is missing required columns: {sorted(missing_query)}")
            assert_same_snp_rows(
                self.baseline_table,
                self.query_table,
                context="query rows must match baseline rows on CHR/SNP/POS",
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
        ``regr_weight``, and every entry in ``baseline_columns``.
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

    def validate(self, *, require_query_alignment: bool = True) -> None:
        """Check the normalized public contract for aggregated results."""
        required = {"CHR", "SNP", "POS", "regr_weight", *self.baseline_columns}
        missing = required - set(self.baseline_table.columns)
        if missing:
            raise ValueError(f"baseline_table is missing required columns: {sorted(missing)}")
        if self.query_columns and self.query_table is None:
            raise ValueError("query_table is required when query_columns are present.")
        if not self.query_columns and self.query_table is not None:
            raise ValueError("query_table was provided but query_columns is empty.")
        if self.query_table is not None:
            missing_query = {"CHR", "SNP", "POS", *self.query_columns} - set(self.query_table.columns)
            if missing_query:
                raise ValueError(f"query_table is missing required columns: {sorted(missing_query)}")
            if require_query_alignment:
                assert_same_snp_rows(
                    self.baseline_table,
                    self.query_table,
                    context="query rows must match baseline rows on CHR/SNP/POS",
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
    treats ``ref_panel`` as the owner of the reference-panel SNP universe:
    chromosome bundles are intersected with ``ref_panel.load_metadata(chrom)``
    before the kernel runs so ``RefPanelConfig.ref_panel_snps_file`` is honored
    without leaking that setting into the calculator interface. For LD-score
    calculation, the annotation file's ``CM`` is the first source. The sidecar
    metadata only fills missing ``CM`` values.
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
        regression_snps: set[str] | None = None,
        config_snapshot: dict[str, Any] | None = None,
    ) -> LDScoreResult:
        """Compute and aggregate LD scores across all chromosomes.

        Parameters
        ----------
        annotation_bundle : AnnotationBundle
            Aligned SNP-level baseline and query annotations.
        ref_panel : RefPanel
            Reference-panel adapter that supplies chromosome readers and
            metadata. Any ``RefPanelConfig.ref_panel_snps_file`` restriction is
            already applied when the workflow calls ``ref_panel.load_metadata``.
            Reference-panel sidecar ``CM`` values are used only to fill missing
            annotation ``CM`` values during LD-score calculation.
        ldscore_config : LDScoreConfig
            LD-window and retained-SNP settings.
        global_config : GlobalConfig
            Shared SNP identifier, genome-build, and logging settings.
        output_config : LDScoreOutputConfig or None, optional
            If provided, write the canonical LD-score result directory after
            the aggregate result is built. Existing canonical files are refused
            unless ``output_config.overwrite`` is true.
            Default is ``None``, which keeps the result in memory only.
        regression_snps : set of str or None, optional
            Optional regression SNP universe used to define the weight table.
            Default is ``None``, which uses the retained reference SNP universe.
        config_snapshot : dict or None, optional
            Optional run metadata forwarded to the output layer. Default is
            ``None``.

        Returns
        -------
        LDScoreResult
            Aggregated cross-chromosome result with aligned metadata and output
            paths if writing was requested.
        """
        print_global_config_banner(type(self).__name__, global_config)
        if annotation_bundle.config_snapshot is not None:
            validate_config_compatibility(
                annotation_bundle.config_snapshot,
                global_config,
                context="AnnotationBundle and LDScoreCalculator runtime config",
            )
        chromosomes = _chromosomes_from_bundle(annotation_bundle)
        LOGGER.info(_format_ldscore_start_message(annotation_bundle, len(chromosomes)))
        chromosome_results: list[ChromLDScoreResult] = []
        for chrom in chromosomes:
            chrom_bundle = _slice_annotation_bundle(annotation_bundle, chrom)
            try:
                chromosome_result = self.compute_chromosome(
                    chrom=chrom,
                    annotation_bundle=chrom_bundle,
                    ref_panel=ref_panel,
                    ldscore_config=ldscore_config,
                    global_config=global_config,
                    regression_snps=regression_snps,
                )
            except ValueError as exc:
                if _warn_and_skip_empty_intersection(exc, chrom):
                    continue
                raise
            chromosome_results.append(chromosome_result)
        if not chromosome_results:
            raise ValueError("No chromosome results were produced after intersecting annotations with the reference panel.")
        result = self._aggregate_chromosome_results(
            chromosome_results,
            global_config=global_config,
            count_config=_count_config_from_ldscore_config(ldscore_config),
        )
        if output_config is not None:
            output_paths = self.output_writer.write(result, output_config)
            result = _replace_result_output_paths(result, output_paths)
            LOGGER.info(f"Wrote LD-score result directory to '{output_config.output_dir}'.")
        LOGGER.info(
            f"Computed LD scores for {len(chromosome_results)} chromosomes "
            f"and {len(result.baseline_table)} retained SNP rows."
        )
        return result

    def compute_chromosome(
        self,
        chrom: str,
        annotation_bundle,
        ref_panel,
        ldscore_config: LDScoreConfig,
        global_config: GlobalConfig,
        regression_snps: set[str] | None = None,
    ) -> ChromLDScoreResult:
        """Compute normalized LD-score outputs for one chromosome.

        The workflow first asks ``ref_panel`` for the prepared chromosome
        metadata, which already reflects any ``ref_panel_snps_file`` restriction.
        It then restricts the chromosome-local ``AnnotationBundle`` to the same
        identifier set so the legacy kernel sees ``B_chrom ∩ A'_chrom`` rather
        than the raw annotation universe ``B_chrom``. ``regression_snps`` is
        applied later when the normalized row table is formed.
        """
        backend = getattr(getattr(ref_panel, "spec", None), "backend", None)
        LOGGER.info(
            f"Computing chromosome {chrom} LD scores with backend '{backend or 'unknown'}' "
            f"and {len(annotation_bundle.metadata)} annotation rows."
        )
        annotation_bundle = _align_annotation_bundle_to_ref_panel(
            annotation_bundle=annotation_bundle,
            ref_panel=ref_panel,
            chrom=chrom,
            global_config=global_config,
        )
        args = _namespace_from_configs(
            chrom=chrom,
            ref_panel=ref_panel,
            ldscore_config=ldscore_config,
            global_config=global_config,
        )
        legacy_bundle = kernel_ldscore.AnnotationBundle(
            metadata=annotation_bundle.metadata.copy(),
            annotations=annotation_bundle.annotation_matrix(include_query=True).copy(),
            baseline_columns=list(annotation_bundle.baseline_columns),
            query_columns=list(annotation_bundle.query_columns),
        )
        if getattr(ref_panel.spec, "backend", None) == "parquet_r2":
            legacy_result = kernel_ldscore.compute_chrom_from_parquet(chrom, legacy_bundle, args, regression_snps)
        else:
            legacy_result = kernel_ldscore.compute_chrom_from_plink(chrom, legacy_bundle, args, regression_snps)
        result = self._wrap_legacy_chrom_result(legacy_result, global_config=global_config, regression_snps=regression_snps)
        LOGGER.info(f"Finished chromosome {chrom} with {len(result.baseline_table)} retained SNP rows.")
        return result

    def _wrap_legacy_chrom_result(
        self,
        legacy_result: _LegacyChromResult | Any,
        global_config: GlobalConfig,
        regression_snps: set[str] | None = None,
    ) -> ChromLDScoreResult:
        """Convert one kernel chromosome result into the typed public result."""
        reference_metadata = legacy_result.metadata.reset_index(drop=True).copy()
        ld_scores = pd.DataFrame(legacy_result.ld_scores, columns=list(legacy_result.ldscore_columns))
        reference_ids = frozenset(build_snp_id_series(reference_metadata, global_config.snp_identifier))
        retained_regression_snps = (
            reference_ids if regression_snps is None else frozenset(reference_ids.intersection(regression_snps))
        )
        regression_keep = build_snp_id_series(reference_metadata, global_config.snp_identifier).isin(retained_regression_snps)
        pos_column = "POS" if "POS" in reference_metadata.columns else "BP"
        regression_weights = np.asarray(legacy_result.w_ld, dtype=np.float32).reshape(-1)
        ldscore_table = pd.concat(
            [
                reference_metadata.loc[regression_keep, ["CHR", "SNP", pos_column]].rename(columns={pos_column: "POS"}).reset_index(drop=True),
                ld_scores.loc[regression_keep].reset_index(drop=True),
                pd.DataFrame({"regr_weight": regression_weights[regression_keep.to_numpy()]}).reset_index(drop=True),
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
            ld_regression_snps=frozenset(build_snp_id_series(baseline_table, global_config.snp_identifier)),
            snp_count_totals=count_map,
            count_config={},
            config_snapshot=global_config,
        )
        result.validate()
        return result

    def _aggregate_chromosome_results(
        self,
        chromosome_results: Sequence[ChromLDScoreResult],
        global_config: GlobalConfig,
        count_config: dict[str, Any] | None = None,
    ) -> LDScoreResult:
        """Concatenate and sum per-chromosome results into one aggregate object."""
        if not chromosome_results:
            raise ValueError("At least one chromosome result is required.")
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
        merged_table = pd.concat(
            [_join_split_tables(result.baseline_table, result.query_table, result.query_columns) for result in chromosome_results],
            axis=0,
            ignore_index=True,
        )
        merged_table = kernel_ldscore.sort_frame_by_genomic_position(merged_table)
        baseline_table, query_table = _split_ldscore_table(
            merged_table,
            baseline_columns=list(chromosome_results[0].baseline_columns),
            query_columns=list(chromosome_results[0].query_columns),
        )
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


def _split_ldscore_table(
    ldscore_table: pd.DataFrame,
    *,
    baseline_columns: list[str],
    query_columns: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Split a merged LD-score table into baseline and optional query tables."""
    metadata_columns = ["CHR", "SNP", "POS"]
    baseline_order = [*metadata_columns, "regr_weight", *baseline_columns]
    query_order = [*metadata_columns, *query_columns]
    missing_baseline = [column for column in baseline_order if column not in ldscore_table.columns]
    if missing_baseline:
        raise ValueError(f"LD-score table is missing baseline columns: {missing_baseline}")
    baseline_table = ldscore_table.loc[:, baseline_order].reset_index(drop=True).copy()
    query_table = None
    if query_columns:
        missing_query = [column for column in query_order if column not in ldscore_table.columns]
        if missing_query:
            raise ValueError(f"LD-score table is missing query columns: {missing_query}")
        query_table = ldscore_table.loc[:, query_order].reset_index(drop=True).copy()
        assert_same_snp_rows(
            baseline_table,
            query_table,
            context="query rows must match baseline rows on CHR/SNP/POS",
        )
    return baseline_table, query_table


def _join_split_tables(
    baseline_table: pd.DataFrame,
    query_table: pd.DataFrame | None,
    query_columns: Sequence[str],
) -> pd.DataFrame:
    """Join split LD-score tables for sorting or regression assembly."""
    if query_table is None:
        return baseline_table.copy()
    assert_same_snp_rows(
        baseline_table,
        query_table,
        context="query rows must match baseline rows on CHR/SNP/POS",
    )
    query_values = query_table.loc[:, list(query_columns)].reset_index(drop=True)
    return pd.concat([baseline_table.reset_index(drop=True), query_values], axis=1)


def _count_records_from_totals(
    *,
    baseline_columns: list[str],
    query_columns: list[str],
    count_totals: dict[str, np.ndarray],
) -> list[dict[str, Any]]:
    """Convert positional count vectors into manifest-friendly column records."""
    columns = [*baseline_columns, *query_columns]
    groups = ["baseline"] * len(baseline_columns) + ["query"] * len(query_columns)
    all_counts = np.asarray(count_totals.get("all_reference_snp_counts"), dtype=np.float64)
    if all_counts.size != len(columns):
        raise ValueError("all_reference_snp_counts length does not match annotation columns.")
    common_raw = count_totals.get("common_reference_snp_counts")
    common_counts = None if common_raw is None else np.asarray(common_raw, dtype=np.float64)
    if common_counts is not None and common_counts.size != len(columns):
        raise ValueError("common_reference_snp_counts length does not match annotation columns.")
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
    """Return manifest count metadata for common-SNP count vectors."""
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
    parser.add_argument("--overwrite", action="store_true", default=False, help="Replace existing fixed output files.")
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
    parser.add_argument(
        "--baseline-annot-sources",
        default=None,
        help="Comma-separated baseline annotation path tokens. If omitted with no query inputs, an all-ones `base` annotation is synthesized.",
    )
    parser.add_argument("--plink-prefix", default=None, help="PLINK prefix token for the reference panel.")
    parser.add_argument(
        "--r2-dir",
        dest="r2_dir",
        default=None,
        help="Build-specific R2 directory containing chr*_r2.parquet and optional chr*_meta.tsv.gz sidecars.",
    )
    parser.add_argument("--snp-identifier", default="chr_pos", help="Identifier mode used to match annotations to the reference panel.")
    parser.add_argument(
        "--genome-build",
        choices=("auto", "hg19", "hg37", "GRCh37", "hg38", "GRCh38"),
        default=None,
        help=(
            "Genome build for chr_pos inputs. Required when --snp-identifier chr_pos "
            "(the default). Use 'auto' to infer hg19/hg38 and 0-based/1-based coordinates "
            "from data. Not used when --snp-identifier rsid."
        ),
    )
    parser.add_argument("--r2-bias-mode", choices=("raw", "unbiased"), default="unbiased", help="Whether parquet R2 values are raw or already unbiased.")
    parser.add_argument("--r2-sample-size", default=None, type=float, help="LD reference sample size used to correct raw parquet R2 values.")
    parser.add_argument(
        "--ref-panel-snps-file",
        default=None,
        help="Optional SNP list defining the retained reference-panel universe A'; the workflow intersects each chromosome annotation bundle with this prepared panel before LD computation.",
    )
    parser.add_argument("--regression-snps-file", default=None, help="Optional SNP list defining the regression SNP set and the written LD-score row set.")
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
    parser.add_argument("--snp-batch-size", default=128, type=int, help="Number of SNPs processed per LD-score sliding batch.")
    parser.add_argument("--yes-really", default=False, action="store_true", help="Allow whole-chromosome LD windows.")
    parser.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"), help="Logging verbosity.")
    return parser


def run_ldscore_from_args(args: argparse.Namespace) -> LDScoreResult:
    """Run LD-score calculation from a parsed CLI namespace.

    The workflow resolves unified path tokens for baseline annotations, optional
    query annotations or BED files, and the reference panel. If no baseline or
    query annotations are supplied, it synthesizes an all-ones ``base``
    annotation over the retained reference-panel metadata. Before calculation
    it preflights ``manifest.json``, ``ldscore.baseline.parquet``, optional
    ``ldscore.query.parquet``, and ``ldscore.log`` under ``output_dir``. For
    each chromosome it intersects annotation rows with
    ``ref_panel.load_metadata(chrom)`` before calling the kernel, then returns
    the normalized public ``LDScoreResult`` with split baseline/query tables.
    The result ``output_paths`` mapping contains data artifacts only.
    """
    from .annotation_builder import AnnotationBuilder

    normalized_args, global_config = _normalize_run_args(args)
    print_global_config_banner("run_ldscore_from_args", global_config)
    _validate_run_args(normalized_args)
    ldscore_config = _ldscore_config_from_args(normalized_args)
    regression_snps = _load_regression_snps(ldscore_config.regression_snps_file, global_config)
    ref_mode = "parquet" if _uses_parquet_reference(normalized_args) else "plink"
    LOGGER.info(
        f"Starting LD-score workflow with reference mode '{ref_mode}', "
        f"output directory '{normalized_args.output_dir}', "
        f"snp_identifier='{global_config.snp_identifier}', genome_build='{global_config.genome_build}'."
    )
    if _has_cli_tokens(normalized_args.baseline_annot_sources):
        source_spec = AnnotationBuildConfig(
            baseline_annot_sources=tuple(split_cli_path_tokens(normalized_args.baseline_annot_sources)),
            query_annot_sources=tuple(split_cli_path_tokens(normalized_args.query_annot_sources)),
            query_annot_bed_sources=tuple(split_cli_path_tokens(getattr(normalized_args, "query_annot_bed_sources", None))),
        )
        annotation_bundle = AnnotationBuilder(global_config).run(source_spec)
        ref_panel = _ref_panel_from_args(normalized_args, global_config)
    else:
        ref_panel = _ref_panel_from_args(normalized_args, global_config)
        annotation_bundle = _pseudo_base_annotation_bundle_from_ref_panel(ref_panel, global_config)
    output_config = _output_config_from_args(normalized_args)
    output_dir = ensure_output_directory(output_config.output_dir, label="LD-score output directory")
    log_path = output_dir / "ldscore.log"
    ensure_output_paths_available(
        [*_expected_ldscore_output_paths(output_dir, bool(annotation_bundle.query_columns)), log_path],
        overwrite=output_config.overwrite,
        label="LD-score output artifact",
    )
    calculator = LDScoreCalculator()
    with workflow_logging("ldscore", log_path, log_level=global_config.log_level):
        log_inputs(
            output_dir=str(output_dir),
            reference_mode=ref_mode,
            snp_identifier=global_config.snp_identifier,
            genome_build=global_config.genome_build,
        )
        result = calculator.run(
            annotation_bundle=annotation_bundle,
            ref_panel=ref_panel,
            ldscore_config=ldscore_config,
            global_config=global_config,
            output_config=output_config,
            regression_snps=regression_snps,
        )
        log_outputs(**result.output_paths)
    return result


def _validate_run_args(args: argparse.Namespace) -> None:
    """Validate public LD-score workflow arguments before loading inputs.

    This validator intentionally lives in the workflow layer rather than in the
    kernel because optional baseline synthesis is a public orchestration rule:
    the numerical kernels still receive an explicit annotation bundle.
    """
    if not _has_cli_tokens(args.baseline_annot_sources) and (
        _has_cli_tokens(args.query_annot_sources) or _has_cli_tokens(getattr(args, "query_annot_bed_sources", None))
    ):
        raise ValueError(_QUERY_REQUIRES_BASELINE_MESSAGE)
    keep = getattr(args, "keep", None)
    if _uses_parquet_reference(args) == bool(args.bfile):
        raise ValueError("Specify exactly one reference-panel mode: parquet or PLINK.")
    if _uses_parquet_reference(args):
        if keep:
            raise ValueError("--keep-indivs-file is only supported in PLINK mode.")
        if args.r2_bias_mode is None:
            args.r2_bias_mode = "unbiased"
        if args.r2_bias_mode == "raw" and args.r2_sample_size is None:
            raise ValueError("--r2-sample-size is required when --r2-bias-mode raw.")
        if args.snp_identifier == "chr_pos" and args.genome_build is None:
            raise ValueError("--genome-build is required in parquet mode when --snp-identifier chr_pos.")
    if args.ld_wind_cm is not None and args.ld_wind_cm <= 0:
        raise ValueError("--ld-wind-cm must be positive.")
    if args.ld_wind_kb is not None and args.ld_wind_kb <= 0:
        raise ValueError("--ld-wind-kb must be positive.")
    if args.ld_wind_snps is not None and args.ld_wind_snps <= 0:
        raise ValueError("--ld-wind-snps must be positive.")
    if getattr(args, "maf_min", None) is not None and not 0 <= args.maf_min <= 0.5:
        raise ValueError("--maf-min must lie in [0, 0.5].")
    if not 0 <= getattr(args, "common_maf_min", 0.05) <= 0.5:
        raise ValueError("--common-maf-min must lie in [0, 0.5].")
    if args.snp_batch_size <= 0:
        raise ValueError("--snp-batch-size must be positive.")


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
        metadata_frames.append(metadata.loc[:, ["CHR", "SNP", "CM", "POS"]].reset_index(drop=True))
    if not metadata_frames:
        raise ValueError("No reference-panel SNP metadata rows are available for pseudo `base` annotation.")
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
    bundle.validate(global_config.snp_identifier)
    return bundle


def run_ldscore(**kwargs) -> LDScoreResult:
    """Run LD-score calculation from Python using public CLI-style names.

    Keyword arguments are interpreted as CLI-equivalent option names without
    leading ``--``; for example ``baseline_annot_sources``, ``query_annot_sources``,
    ``query_annot_bed_sources``, ``plink_prefix``, ``r2_dir``,
    ``keep_indivs_file``, ``snp_batch_size``, ``common_maf_min``, and
    ``output_dir``. Shared runtime assumptions such as ``snp_identifier`` and
    ``genome_build`` must be supplied through ``set_global_config(...)`` first,
    while per-run controls such as ``ref_panel_snps_file`` and
    ``regression_snps_file`` remain ordinary keyword arguments here.

    When ``baseline_annot_sources`` and query inputs are omitted, the workflow
    builds a synthetic all-ones baseline column named ``base`` from retained
    reference-panel metadata. ``query_annot_sources`` and
    ``query_annot_bed_sources`` require explicit baseline annotations because
    query columns are interpreted relative to that baseline SNP universe.

    Returns
    -------
    LDScoreResult
        Aggregated result with ``baseline_table``, optional ``query_table``,
        manifest count records, and canonical output paths.
    """
    forbidden = sorted({"snp_identifier", "genome_build", "log_level"} & set(kwargs))
    if forbidden:
        joined = ", ".join(forbidden)
        raise ValueError(f"Python run_ldscore() no longer accepts {joined}; call set_global_config(...) first.")
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
        }
        & set(kwargs)
    )
    if removed:
        joined = ", ".join(removed)
        raise ValueError(f"Python run_ldscore() no longer accepts removed IO argument(s): {joined}.")
    parser = build_parser()
    defaults = vars(parser.parse_args(["--output-dir", "placeholder"]))
    global_config = get_global_config()
    defaults["snp_identifier"] = global_config.snp_identifier
    defaults["genome_build"] = global_config.genome_build
    defaults["log_level"] = global_config.log_level
    defaults.update(kwargs)
    args = argparse.Namespace(**defaults)
    return run_ldscore_from_args(args)


def main(argv: Sequence[str] | None = None) -> LDScoreResult:
    """Command-line entry point for the LD-score workflow."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return run_ldscore_from_args(args)


def _normalize_run_args(args: argparse.Namespace) -> tuple[argparse.Namespace, GlobalConfig]:
    """Normalize CLI-style args and derive the shared ``GlobalConfig`` object."""
    normalized_mode = normalize_snp_identifier_mode(args.snp_identifier)
    normalized_args = argparse.Namespace(**vars(args))
    for attr in ("query_annot_chr", "baseline_annot_chr", "bfile_chr", "r2_table_chr", "frqfile_chr"):
        if not hasattr(normalized_args, attr):
            setattr(normalized_args, attr, None)
    for attr in ("query_annot_sources", "baseline_annot_sources", "plink_prefix", "r2_dir", "query_annot_bed_sources", "keep_indivs_file"):
        if not hasattr(normalized_args, attr):
            setattr(normalized_args, attr, None)
    for attr in ("ref_panel_snps_file", "regression_snps_file"):
        if not hasattr(normalized_args, attr):
            setattr(normalized_args, attr, None)
    if not hasattr(normalized_args, "maf_min"):
        normalized_args.maf_min = None
    if not hasattr(normalized_args, "common_maf_min"):
        normalized_args.common_maf_min = 0.05
    if not hasattr(normalized_args, "snp_batch_size"):
        normalized_args.snp_batch_size = 128
    normalized_args.snp_identifier = normalized_mode
    normalized_args.output_dir = normalize_path_token(args.output_dir)
    normalized_args.r2_dir = _r2_dir_from_args(normalized_args)
    normalized_args.keep_indivs_file = normalize_optional_path_token(getattr(args, "keep_indivs_file", None))
    normalized_args.ref_panel_snps_file = normalize_optional_path_token(getattr(args, "ref_panel_snps_file", None))
    normalized_args.regression_snps_file = normalize_optional_path_token(getattr(args, "regression_snps_file", None))
    # The numerical kernel still consumes the historical namespace shape.
    normalized_args.query_annot = normalized_args.query_annot_sources
    normalized_args.baseline_annot = normalized_args.baseline_annot_sources
    normalized_args.bfile = normalized_args.plink_prefix
    normalized_args.r2_table = None
    normalized_args.frqfile = None
    normalized_args.keep = normalized_args.keep_indivs_file
    if normalized_mode == "rsid":
        global_config = GlobalConfig(
            snp_identifier=normalized_mode,
            genome_build=normalize_genome_build(getattr(args, "genome_build", None)),
            log_level=getattr(args, "log_level", "INFO"),
        )
        normalized_args.genome_build = global_config.genome_build
    else:
        resolved_genome_build = _resolve_ldscore_chr_pos_genome_build(
            normalized_args,
            getattr(args, "genome_build", None),
        )
        normalized_args.genome_build = resolved_genome_build
        global_config = GlobalConfig(
            snp_identifier=normalized_mode,
            genome_build=resolved_genome_build,
            log_level=getattr(args, "log_level", "INFO"),
        )
    return normalized_args, global_config


def _resolve_ldscore_chr_pos_genome_build(args: argparse.Namespace, genome_build: str | None) -> str:
    normalized = normalize_genome_build(genome_build)
    if normalized is None:
        raise ValueError(
            "genome_build is required when snp_identifier='chr_pos'. "
            "Pass --genome-build auto, --genome-build hg19, or --genome-build hg38."
        )
    if normalized != "auto":
        return normalized

    resolved: list[tuple[str, str]] = []
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
        raise ValueError(
            "Cannot infer --genome-build for LD-score chr_pos inputs because no chromosome-suite "
            "annotation or R2 parquet build metadata was available."
        )
    builds = {build for _label, build in resolved}
    if len(builds) != 1:
        details = ", ".join(f"{label}={build}" for label, build in resolved)
        raise ValueError(f"LD-score input genome build sources disagree: {details}.")
    return resolved[0][1]


def _load_regression_snps(path: str | None, global_config: GlobalConfig) -> set[str] | None:
    """Load ``LDScoreConfig.regression_snps_file`` using the active identifier mode."""
    if not path:
        return None
    return read_global_snp_restriction(
        resolve_scalar_path(path, label="regression SNP list"),
        global_config.snp_identifier,
        genome_build=global_config.genome_build,
    )


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
        raise ValueError(f"Conflicting R2 parquet genome-build metadata in '{r2_dir}': {details}.")
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


def _ref_panel_from_args(args: argparse.Namespace, global_config: GlobalConfig):
    """Build the reference-panel adapter that owns the ``A -> A'`` restriction."""
    from ._kernel.ref_panel import RefPanelLoader

    ref_panel_snps_file = normalize_optional_path_token(getattr(args, "ref_panel_snps_file", None))
    r2_dir = _r2_dir_from_args(args)
    if r2_dir is not None:
        spec = RefPanelConfig(
            backend="parquet_r2",
            r2_dir=r2_dir,
            r2_bias_mode=getattr(args, "r2_bias_mode", None),
            sample_size=getattr(args, "r2_sample_size", None),
            ref_panel_snps_file=ref_panel_snps_file,
            maf_min=getattr(args, "maf_min", None),
            keep_indivs_file=getattr(args, "keep_indivs_file", None),
        )
    else:
        spec = RefPanelConfig(
            backend="plink",
            plink_prefix=getattr(args, "plink_prefix", None),
            ref_panel_snps_file=ref_panel_snps_file,
            maf_min=getattr(args, "maf_min", None),
            keep_indivs_file=getattr(args, "keep_indivs_file", None),
        )
    return RefPanelLoader(global_config).load(spec)


def _ldscore_config_from_args(args: argparse.Namespace) -> LDScoreConfig:
    """Build an ``LDScoreConfig`` from normalized LD-score arguments."""
    return LDScoreConfig(
        ld_wind_snps=getattr(args, "ld_wind_snps", None),
        ld_wind_kb=getattr(args, "ld_wind_kb", None),
        ld_wind_cm=getattr(args, "ld_wind_cm", None),
        regression_snps_file=getattr(args, "regression_snps_file", None),
        snp_batch_size=getattr(args, "snp_batch_size", 128),
        common_maf_min=getattr(args, "common_maf_min", 0.05),
        whole_chromosome_ok=getattr(args, "yes_really", False),
    )


def _output_config_from_args(args: argparse.Namespace) -> LDScoreOutputConfig:
    """Translate LD-score CLI arguments into the canonical directory output config."""
    return LDScoreOutputConfig(
        output_dir=normalize_path_token(args.output_dir),
        overwrite=getattr(args, "overwrite", False),
    )


def _expected_ldscore_output_paths(output_dir: Path, has_query: bool) -> list[Path]:
    """Return canonical LD-score output paths written by the directory writer."""
    paths = [output_dir / "manifest.json", output_dir / "ldscore.baseline.parquet"]
    if has_query:
        paths.append(output_dir / "ldscore.query.parquet")
    return paths


def _replace_result_output_paths(result: LDScoreResult, output_paths: dict[str, str]) -> LDScoreResult:
    """Return ``result`` with updated artifact-path metadata after writing outputs."""
    return LDScoreResult(
        baseline_table=result.baseline_table,
        query_table=result.query_table,
        count_records=result.count_records,
        baseline_columns=result.baseline_columns,
        query_columns=result.query_columns,
        ld_reference_snps=result.ld_reference_snps,
        ld_regression_snps=result.ld_regression_snps,
        chromosome_results=result.chromosome_results,
        output_paths=dict(output_paths),
        count_config=dict(result.count_config),
        config_snapshot=result.config_snapshot,
    )


def _chromosomes_from_bundle(annotation_bundle) -> list[str]:
    """Infer the chromosome processing order from an annotation bundle."""
    chromosomes = getattr(annotation_bundle, "chromosomes", None)
    if chromosomes:
        return list(chromosomes)
    return sorted(annotation_bundle.metadata["CHR"].astype(str).unique().tolist())


def _slice_annotation_bundle(annotation_bundle, chrom: str):
    """Return the per-chromosome view of an annotation bundle."""
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
    )


def _align_annotation_bundle_to_ref_panel(annotation_bundle, ref_panel, chrom: str, global_config: GlobalConfig):
    """Restrict one chromosome bundle to the prepared reference-panel universe.

    ``AnnotationBuilder`` defines the annotation universe ``B``. The reference
    panel owns the optional ``A -> A'`` restriction through
    ``RefPanelConfig.ref_panel_snps_file``. This helper materializes the intended
    compute-time universe ``B_chrom ∩ A'_chrom`` immediately before the legacy
    kernel call. It restricts annotation rows but does not replace annotation
    metadata; for LD-score calculation, the annotation file's ``CM`` is the
    first source, and sidecar metadata only fills missing ``CM`` values later in
    the kernel.
    """
    reference_metadata = ref_panel.load_metadata(chrom)
    reference_ids = set(build_snp_id_series(reference_metadata, global_config.snp_identifier))
    keep = build_snp_id_series(annotation_bundle.metadata, global_config.snp_identifier).isin(reference_ids)
    if not bool(keep.any()):
        backend = getattr(getattr(ref_panel, "spec", None), "backend", None)
        intersection = "parquet" if backend == "parquet_r2" else "PLINK"
        raise ValueError(f"No retained annotation SNPs remain on chromosome {chrom} after {intersection} intersection.")
    if bool(keep.all()):
        return annotation_bundle
    return type(annotation_bundle)(
        metadata=annotation_bundle.metadata.loc[keep].reset_index(drop=True),
        baseline_annotations=annotation_bundle.baseline_annotations.loc[keep].reset_index(drop=True),
        query_annotations=annotation_bundle.query_annotations.loc[keep].reset_index(drop=True),
        baseline_columns=list(annotation_bundle.baseline_columns),
        query_columns=list(annotation_bundle.query_columns),
        chromosomes=list(getattr(annotation_bundle, "chromosomes", [str(chrom)])),
        source_summary=dict(getattr(annotation_bundle, "source_summary", {})),
        config_snapshot=getattr(annotation_bundle, "config_snapshot", None),
    )


def _warn_and_skip_empty_intersection(error: ValueError, chrom: str) -> bool:
    """Warn and signal skip when a chromosome loses all SNPs after reference intersection."""
    message = str(error)
    if not message.startswith("No retained annotation SNPs remain on chromosome "):
        return False
    if not any(suffix in message for suffix in (" after parquet intersection.", " after PLINK intersection.")):
        return False
    warnings.warn(f"Skipping chromosome {chrom}: {message}", UserWarning, stacklevel=3)
    return True


def _namespace_from_configs(chrom: str, ref_panel, ldscore_config: LDScoreConfig, global_config: GlobalConfig) -> argparse.Namespace:
    """Build the legacy-kernel namespace expected by the chromosome backends."""
    spec = getattr(ref_panel, "spec", None)
    backend = getattr(spec, "backend", None)
    bfile = None
    r2_table = None
    frqfile = None
    if backend == "plink" and getattr(spec, "plink_prefix", None) is not None:
        bfile = resolve_plink_prefix(spec.plink_prefix, chrom=chrom)
    if backend == "parquet_r2":
        if hasattr(ref_panel, "resolve_r2_paths"):
            r2_table = ",".join(ref_panel.resolve_r2_paths(chrom, required=False)) or None
    if hasattr(ref_panel, "resolve_metadata_paths"):
        frqfile = ",".join(ref_panel.resolve_metadata_paths(chrom)) or None
    return argparse.Namespace(
        out=None,
        query_annot=None,
        query_annot_chr=None,
        baseline_annot=None,
        baseline_annot_chr=None,
        bfile=bfile,
        bfile_chr=None,
        r2_table=r2_table,
        r2_table_chr=None,
        snp_identifier=global_config.snp_identifier,
        genome_build=global_config.genome_build,
        r2_bias_mode=getattr(spec, "r2_bias_mode", None),
        r2_sample_size=getattr(spec, "sample_size", None),
        frqfile=frqfile,
        frqfile_chr=None,
        keep=getattr(spec, "keep_indivs_file", None),
        ld_wind_snps=ldscore_config.ld_wind_snps,
        ld_wind_kb=ldscore_config.ld_wind_kb,
        ld_wind_cm=ldscore_config.ld_wind_cm,
        maf=None,
        maf_min=None,
        common_maf_min=ldscore_config.common_maf_min,
        snp_batch_size=ldscore_config.snp_batch_size,
        per_chr_output=False,
        yes_really=ldscore_config.whole_chromosome_ok,
        log_level=global_config.log_level,
    )
