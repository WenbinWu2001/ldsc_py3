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
runs. Query annotations remain partitioned-heritability inputs: they are accepted only
when explicit baseline annotations are supplied.

Writing workflows create ``diagnostics/ldscore.log`` and detailed
``diagnostics/chromosome-logs/chrN.log`` files after output preflight.
``LDScoreCalculator.run(..., output_config=None)`` creates no files.
Every chromosome executes in a worker with one BLAS thread, including
serial execution; workers are reused across query batches.
"""

from __future__ import annotations

import argparse
from contextlib import ExitStack, contextmanager
from copy import copy
from tempfile import TemporaryDirectory
from dataclasses import dataclass, field, replace as dataclass_replace
import logging
from pathlib import Path
import sys
from typing import Any, Sequence

import numpy as np
import pandas as pd

from .._runtime.serial_calculation import serial_calculation
from .diagnostics import (
    chromosome_log,
    emit_reference_snp_funnel_notice,
    fail_gene_list_resolution,
    fail_unusable_queries,
    ldscore_logs,
    log_query_outcomes,
    log_reference_snp_funnel_warning,
    reference_snp_funnel_path,
)
from .batch_execution import finalize_batch, run_chromosome_batch, run_query_batches
from .direct_batches import DirectChromosomeFragment, stage_direct_chromosome, stage_query_table
from .fragments import StagedLDScoreBatch
from .._core.chromosomes import chrom_sort_key
from .._runtime.cli_help import CLIHelpFormatter, CHROMOSOME_PATH_HELP, SCALAR_PATH_HELP, output_directory_argument
from .._runtime.workflow_logging import args_output_scope, LOG_LEVEL_HELP
from .._runtime.parallelism import (
    WORKERS_HELP, chromosome_pool, parse_workers, resolve_worker_count, validate_query_batch_size, validate_workers,
)
from .._runtime.progress import report_phase, PhaseProgress
from .._annotation.saved_results import AnnotationResult
from .._annotation.storage import TsvDiagnostics, AnnotationWorkspace
from .source import LDScoreSource, validate_ldscore_schemas
from .._core.snp_identity import coerce_identity_drop_frame, IDENTITY_DROP_COLUMNS, identity_mode_family, is_allele_aware_mode
from .._core.columns import normalize_genome_build
from .._core.snp_identity import emit_dropped_snps_notice, normalize_snp_identifier_mode
from ..config import (
    DEFAULT_SNP_BATCH_SIZE,
    GlobalConfig,
    AnnotationBuildConfig,
    LDScoreConfig,
    RefPanelConfig,
    get_global_config,
    print_global_config_banner,
    validate_config_compatibility,
)
from .._genome.hm3 import packaged_hm3_curated_map_path
from .outputs import (
    LDSCORE_LOG,
    LDScoreDirectoryWriter,
    LDScoreOutputConfig,
    QueryDiagnostics,
    REGRESSION_LD_SCORE_COLUMN,
    common_count_config,
    ref_metadata_file,
    snp_universe_policy_record,
)
from .._runtime.paths import (
    ensure_output_directory,
    normalize_optional_path_token,
    normalize_path_token,
    output_directory_path,
    resolve_scalar_path,
    split_cli_path_tokens,
)
from .._runtime.workflow_logging import (
    log_inputs, log_outputs, materializing_overwrite_guard, overwrite_failure_marker, workflow_logging,
)
from . import compute as kernel_ldscore
from .._refpanel.funnel import (
    REFERENCE_PANEL,
    annotation_snp_grid_emptied,
    concat_funnels,
    funnel_summary,
    read_reference_snp_funnel,
)
from .._refpanel.prepared import ChromosomeAnnotations, MappedAnnotations
from .._core import regions as kernel_regions
from .._annotation.overlap import OverlapContribution, sum_overlap_contributions
from .._annotation.semantics import combine_annotation_types
from .._core.regions import REGR_SNPS_EXCLUDE_REGIONS_CHOICES, regr_snps_exclude_regions_choice_to_presets
from .._refpanel.r2_format import SIDECAR_FLOAT_FORMAT
from .._core.identifiers import read_snp_restriction_keys
from .._core.snp_identity import RestrictionIdentityKeys
from .._core.row_alignment import assert_same_snp_rows, concat_chromosome_tables, sort_frame_by_genomic_position
from ..errors import (
    EmptyReferenceSNPs,
    LDSCConfigError,
    LDSCInputError,
    LDSCInternalError,
    LDSCUsageError,
    LDSCUserError,
)
from .._annotation.query_outcomes import (
    QueryAnnotationStatus, gene_query_statuses, gene_viability_errors, all_query_annotations_skipped_message,
    emit_bed_query_notice, emit_gene_gate_b_notice, emit_resolved_only_notice,
)
from .._genome.build_policy import GenomeBuildCheck
from .._refpanel.adapters import PlinkRefPanel, R2MatrixRefPanel
from .._annotation.interval_queries import execution_query_bundle
from .._annotation.gene_lists import GeneCatalog
from .._annotation.gene_query_storage import resolve_gene_lists_staged, persistent_gene_diagnostics
from .._annotation.preflight import (
    annotation_columns,
    require_usable_bed_query_names,
    validate_gene_query_names,
    inspect_annotation_inputs,
)
from .._runtime.preflight import InputGate
from .gene_index.scoring import run_indexed_ldscore
from .._runtime.worker_runtime import initialize_worker


LOGGER = logging.getLogger("LDSC.ldscore_calculator")


_QUERY_REQUIRES_BASELINE_MESSAGE = (
    "ldscore cannot run query annotations without baseline annotations. Most likely "
    "a prebuilt, BED, or gene-list query source was supplied without `--baseline-annot-sources`. "
    "Pass matching baseline annotations, or create an explicit all-ones `base` "
    "baseline annotation over the same SNP rows as the query annotations before estimating "
    "partitioned heritability."
)
_LDSCORE_INTERSECTION_DOC = (
    "docs/user-guide/troubleshooting.md#ldscore-no-reference-snps-remain-on-the-annotation-snp-grid"
)


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
    baseline_columns, query_columns : list of str
        Ordered annotation LD-score columns in the baseline and query tables.
    snp_count_totals : dict of str to numpy.ndarray, optional
        All-SNP (``all_reference_snp_counts``) and common-SNP
        (``common_reference_snp_counts``) count vectors over the reference
        SNPs, one entry per baseline column followed by each query column.
    annotation_types : dict of str to str, optional
        Per-column ``binary`` or ``quantitative`` classification used only for
        interpretation messages and persisted provenance.
    reference_snp_funnel : pandas.DataFrame or None, optional
        SNPs removed by each reference-preparation step on this chromosome;
        see ``ldsc._refpanel.funnel``.
    """
    chrom: str
    baseline_table: pd.DataFrame
    query_table: pd.DataFrame | None
    baseline_columns: list[str]
    query_columns: list[str]
    snp_count_totals: dict[str, np.ndarray] = field(default_factory=dict, repr=False)
    config_snapshot: GlobalConfig | None = None
    overlap: OverlapContribution | None = field(default=None, repr=False)
    reference_snp_count: int = 0
    regression_selected_snp_count: int = 0
    regression_region_removed_snp_count: int = 0
    identity_drops: TsvDiagnostics | None = field(default=None, repr=False)
    annotation_types: dict[str, str] = field(default_factory=dict)
    reference_snp_funnel: pd.DataFrame | None = field(default=None, repr=False)


@dataclass(frozen=True)
class LDScoreResult:
    """Materialized LD scores for one cross-chromosome calculation batch.

    ``LDScoreCalculator.run`` returns this object for a single prepared-input
    batch without an output configuration. Writing workflows return
    ``LDScoreSource`` so completed query tables can be released.

    Parameters
    ----------
    baseline_table : pandas.DataFrame
        Cross-chromosome table persisted as ``ldscore.baseline.parquet`` when
        the result is written. Required columns are ``CHR``, ``SNP``, ``POS``,
        ``regression_ld_scores``, and every entry in ``baseline_columns``.
    query_table : pandas.DataFrame or None
        Optional materialized cross-chromosome query table. A standalone write
        saves it as ``ldscore.query.parquet``; sequential output uses the
        applicable numbered batch filename. Required columns are ``CHR``,
        ``SNP``, ``POS``, and every entry in ``query_columns``.
    count_records : list of dict
        Manifest-ready count records. Each record names an annotation column and
        its all-SNP and optional common-SNP counts.
    count_config : dict, optional
        Manifest-ready metadata describing how optional common-SNP counts were
        computed, including the inclusive common-MAF threshold.
    baseline_columns, query_columns : list of str
        Ordered LD-score columns available to regression workflows.
    output_paths : dict, optional
        Paths written by ``LDScoreDirectoryWriter`` when output was requested.
    config_snapshot : GlobalConfig or None, optional
        Shared configuration active when the result was computed.
    query_statuses : tuple of QueryAnnotationStatus, optional
        Ordered per-source outcomes for BED and gene-list query annotations.
    gene_list_batch : object or None, optional
        Compact gene-list summaries and a streamed audit. Completed command
        results reference persistent diagnostics.
    identity_drops_by_chrom : dict, optional
        Chromosome diagnostic artifacts with explicit ``path`` and ``frames()``
        access. Published results reference persistent files. A calculator run
        without output configuration owns complete in-memory diagnostics,
        which remain available after the input bundle is closed.
    chromosome_scope : dict, optional
        Validated input chromosome coverage and effective analysis scope.
        Empty for workflows that do not perform query coverage preflight.
    legacy_ldscore_conversion : dict or None, optional
        Explicit converter provenance for a legacy LD-score suite. ``None``
        for natively computed LDSC3 results.
    ld_window : dict or None, optional
        The run's LD window and whole-chromosome override, recorded in
        ``metadata.json``: ``unit``, ``value`` and ``whole_chromosome_ok``.
        Indexed runs take the index's window. ``None`` for converted legacy
        LD scores, whose window is unknown.
    annotation_types : dict of str to str, optional
        Per-column semantic classification used for interpretation only.
    ld_window_snp_counts : pandas.DataFrame or None, optional
        For an R2 matrix with ``--ld-wind-snps``, one row per chromosome
        with ``CHR``, ``panel_snps`` (sidecar rows), ``kept_snps`` (SNPs this
        run kept), ``ld_wind_snps``, and ``counting_basis="panel_rows"``.
        ``None`` for kb/cM windows and PLINK reference panels.
    reference_snp_funnel : pandas.DataFrame or None, optional
        Reference-SNP funnel of every chromosome, including skipped ones: the
        SNPs each reference-preparation step removed, in the backend's order.
        Columns are ``funnel`` (``reference_panel`` or ``annotation_snp_grid``),
        ``chrom``, ``backend``, ``step_index``, ``step``, ``removed`` and
        ``remaining``; ``removed`` is missing for a step that did not run.
        Writing runs save it as ``diagnostics/reference_snp_funnel.tsv``.
        ``None`` when the reference panel records no funnel.
    """
    baseline_table: pd.DataFrame
    query_table: pd.DataFrame | None
    count_records: list[dict[str, Any]]
    baseline_columns: list[str]
    query_columns: list[str]
    output_paths: dict[str, str] = field(default_factory=dict)
    count_config: dict[str, Any] = field(default_factory=dict)
    config_snapshot: GlobalConfig | None = None
    overlap: "LDScoreOverlap | None" = field(default=None, repr=False)
    query_statuses: tuple[QueryAnnotationStatus, ...] = ()
    gene_list_batch: Any | None = field(default=None, repr=False)
    snp_universe_policy: dict[str, Any] | None = None
    index_provenance: dict[str, str] | None = None
    legacy_ldscore_conversion: dict[str, Any] | None = None
    ld_window: dict[str, Any] | None = None
    identity_drops_by_chrom: dict[str, TsvDiagnostics] = field(default_factory=dict, repr=False)
    annotation_types: dict[str, str] = field(default_factory=dict)
    chromosome_scope: dict[str, Any] = field(default_factory=dict)
    ld_window_snp_counts: pd.DataFrame | None = field(default=None, repr=False)
    reference_snp_funnel: pd.DataFrame | None = field(default=None, repr=False)

    def validate(self) -> None:
        """Check the normalized public contract for aggregated results."""
        validate_ldscore_schemas(
            self.baseline_columns, self.query_columns, self.baseline_table.columns,
            None if self.query_table is None else self.query_table.columns,
        )
        if self.query_table is not None:
            assert_same_snp_rows(
                self.baseline_table,
                self.query_table,
                context="query rows must match baseline rows on CHR/SNP/POS",
                snp_identifier=getattr(self.config_snapshot, "snp_identifier", "chr_pos_allele_aware"),
            )

    def read_queries(self, columns):
        """Select named queries from this materialized result without file I/O.

        Parameters
        ----------
        columns : sequence of str
            Query names in the desired output column order.

        Returns
        -------
        pandas.DataFrame or None
            SNP metadata followed by selected query values, in existing row
            order with a reset index. An empty selection returns metadata
            only, or ``None`` when this result has no query table.

        Raises
        ------
        LDSCInputError
            Nonempty query selection is requested from a baseline-only result.
        KeyError
            A selected name is absent from the materialized table.
        """
        if self.query_table is None:
            if columns:
                raise LDSCInputError("LD-score result has no query table.")
            return None
        metadata = [name for name in ("CHR", "SNP", "POS", "A1", "A2") if name in self.query_table]
        selected = self.query_table.loc[:, metadata + list(columns)]
        if not isinstance(selected.index, pd.RangeIndex) or selected.index.start != 0:
            selected = selected.reset_index(drop=True)
        return selected


def _log_annotation_classification(result):
    """Report whole-run type counts; full fitted column identities stay in metadata."""
    counts = pd.Series(list(result.annotation_types.values()), dtype="string").value_counts()
    LOGGER.info("Annotation classification: binary=%d; quantitative=%d. Classification does not change fitted coefficients; see annotation_types in LD-score metadata.",
                counts.get("binary", 0), counts.get("quantitative", 0))
    if LOGGER.isEnabledFor(logging.DEBUG):
        LOGGER.debug("Annotation types: %s", result.annotation_types)


@dataclass(frozen=True)
class _RunContext:
    """What a direct LD-score run settled before ``LDScoreCalculator`` computes it.

    The command-line workflow settles everything first and says so;
    ``LDScoreCalculator.run`` settles nothing, so the calculation does each
    remaining step once and repeats none.

    Attributes
    ----------
    global_config : GlobalConfig
        The run's SNP identity, genome build and log level.
    ref_panel : RefPanel
        The reference-panel adapter, built with ``global_config``.
    log_path : Path or None
        The open workflow log ``diagnostics/ldscore.log``, whose output
        preflight already ran; ``None`` when the calculation preflights its
        outputs and opens the log itself.
    build_resolved : bool
        The genome build in ``global_config`` and ``ref_panel`` was resolved
        against every input.
    reference_checked : bool
        Every chromosome's R2 matrix footer and sidecar were checked.
    """

    global_config: GlobalConfig
    ref_panel: Any
    log_path: Path | None = None
    build_resolved: bool = False
    reference_checked: bool = False


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

    def run(
        self, annotation_bundle, ref_panel, ldscore_config: LDScoreConfig,
        global_config: GlobalConfig, output_config: LDScoreOutputConfig | None = None,
        regression_snps=None, regression_regions=None,
    ):
        """Calculate sequential query batches over one fixed set of reference SNPs.

        Parameters
        ----------
        annotation_bundle : AnnotationBundle or AnnotationResult
            Prepared baseline and query inputs. ``AnnotationBundle.from_frames``
            supplies a small in-memory dataset without filesystem staging. An
            annotation result (``AnnotationResult``) from ``run_annotate``
            needs ``output_config``; it is prepared inside the output
            directory and the preparation is removed when the run ends.
        ref_panel : RefPanel
            Resolved PLINK or R2 matrix adapter; input files are read without mutation.
        ldscore_config : LDScoreConfig
            Scientific settings, query execution width, and chromosome worker
            count. Chromosome workers run within each active query batch.
        global_config : GlobalConfig
            Shared identifier, genome-build, and logging configuration.
        output_config : LDScoreOutputConfig or None, optional
            Required for multiple query batches. If supplied, all batch files
            stay private until success. If omitted, one materialized batch and
            complete diagnostics are returned without filesystem writes.
        regression_snps : set of str or RestrictionIdentityKeys or None, optional
            Already-resolved regression SNP candidates. ``None`` selects all
            retained reference rows before region exclusions. This low-level
            method does not load ``ldscore_config.regr_snps_file`` or apply the
            public workflow's packaged HapMap3 default.
        regression_regions : RegionIntervals or None, optional
            Already-resolved region exclusions for the regression SNPs.
            ``None`` applies no region exclusions. Neither this nor
            ``regression_snps`` changes the reference SNPs used for
            annotation counts and LD scores.

        Returns
        -------
        LDScoreSource or LDScoreResult
            Saved-artifact handle when writing; otherwise an in-memory result.
            Explicit ``read_queries(names)`` on a saved handle has no width cap.

        Raises
        ------
        LDSCUsageError
            Output is omitted for multiple batches or an annotation result
            (``AnnotationResult``), annotations still need disk preparation, or
            reference-metadata export is requested without an output directory.
        LDSCInputError
            An R2 matrix has an invalid chromosome footer or sidecar.
            Every chromosome is checked before any is computed, and one error
            lists all issues. Also raised, naming every such chromosome, when
            no reference-panel SNP left on a chromosome has an annotation row.

        Notes
        -----
        The annotation rows define the candidate SNPs: reference-panel SNPs
        without an annotation row are excluded, counted in the result's
        ``reference_snp_funnel``, and reported once as a WARNING in the log.
        A chromosome emptied by an explicit restriction or by genotype QC is
        skipped with a WARNING.
        Each direct batch repeats reference traversal. Writing workers stage
        their post-kernel float32 query scores in private chromosome fragments;
        the parent combines scientific summaries and streams those fragments
        without constructing a genome-wide query matrix. Query workspace scales
        with active batch width and chromosome workers. Counts and LD
        accumulation use float64; canonical output retains its usual precision.
        The calculator borrows inputs and never closes the caller's bundle.
        Each worker uses 1 BLAS thread, including
        serial execution. Writing runs retain chromosome logs under
        ``diagnostics/chromosome-logs``; no-output calls create no files.
        """
        return self._run(annotation_bundle, ldscore_config, output_config, regression_snps, regression_regions,
                         context=_RunContext(global_config, ref_panel))

    @report_phase(LOGGER, "computation", "LD-score calculation")
    def _run(self, annotation_bundle, ldscore_config, output_config=None, regression_snps=None,
             regression_regions=None, *, context, skipped_funnels=()):
        """Implement :meth:`run`; ``context`` (``_RunContext``) says what the caller already settled.

        An annotation result (``AnnotationResult``) is prepared inside this
        run's output directory, so the annotate output directory, an input
        here, is never written; the preparation is removed when the run ends.
        ``skipped_funnels`` are the reference-SNP funnels of chromosomes the
        synthetic ``base`` preparation already skipped.
        """
        arguments = (ldscore_config, output_config, regression_snps, regression_regions)
        options = dict(context=context, skipped_funnels=skipped_funnels)
        if not isinstance(annotation_bundle, AnnotationResult):
            return self._run_prepared(annotation_bundle, *arguments, **options)
        if output_config is None:
            raise LDSCUsageError(_annotation_result_message("LDScoreCalculator.run() without output_config") +
                                 " Alternatively pass output_config, and run() prepares it inside the output directory.")
        with annotation_bundle.prepared_in(output_config.output_dir) as bundle:
            return self._run_prepared(bundle, *arguments, **options)

    def _run_prepared(self, annotation_bundle, ldscore_config, output_config, regression_snps, regression_regions,
                      *, context, skipped_funnels):
        """Run the LD-score calculation for an annotation bundle that needs no further preparation.

        Resolves the genome build, checks the R2 matrix and opens the workflow
        log only where ``context`` says the caller has not.
        """
        ref_panel = context.ref_panel
        if ldscore_config.export_ref_metadata and output_config is None:
            raise LDSCUsageError(
                "LDScoreCalculator.run() requires output_config when export_ref_metadata=True. "
                "Reference metadata is a filesystem artifact, so pass LDScoreOutputConfig(output_dir=...)."
            )
        if ldscore_config.export_ref_metadata and _reference_metadata_sources(ref_panel)[0] == "parquet_sidecar":
            raise LDSCUsageError(
                "ldscore cannot export reference metadata from an R2 matrix: --export-ref-metadata "
                "(LDScoreConfig.export_ref_metadata) applies to a PLINK reference panel (--plink-prefix). "
                "The R2 matrix sidecars chr<N>_meta.tsv.gz already hold the same CHR, POS, SNP, A1, A2, CM "
                "and MAF metadata; use them directly and drop the option."
            )
        global_config = context.global_config
        if not context.build_resolved:
            global_config, ref_panel = _resolve_calculator_build(annotation_bundle, ref_panel, global_config)
        chromosomes = _chromosomes_from_bundle(annotation_bundle)
        if not context.reference_checked:
            _check_parquet_reference(ref_panel, chromosomes, ldscore_config, global_config)
        with ExitStack() as stack:
            output_dir = None if output_config is None else output_config.output_dir
            log_paths = stack.enter_context(ldscore_logs(
                output_dir, chromosomes, overwrite=bool(output_config and output_config.overwrite),
                log_level=global_config.log_level, main_log=context.log_path))
            worker_count = resolve_worker_count(ldscore_config.workers, len(chromosomes))
            export_root = (stack.enter_context(TemporaryDirectory(prefix=".ldsc-ref-", dir=output_dir))
                           if ldscore_config.export_ref_metadata else None)
            pool = stack.enter_context(chromosome_pool(
                worker_count, operation="ldscore chromosome computation", initializer=_init_worker,
                initargs=(regression_snps, regression_regions, global_config.log_level)))
            LOGGER.info(_format_ldscore_start_message(annotation_bundle, len(annotation_bundle.chromosomes)))
            cm_source, maf_source = _reference_metadata_sources(ref_panel)
            LOGGER.info("Reference-panel CM source: %s; MAF source: %s (annotation CM/MAF are ignored).",
                        cm_source, maf_source)
            producer = _DirectProducer(self, annotation_bundle, ref_panel, ldscore_config, global_config,
                                       regression_regions, log_paths, output_dir=output_dir, export_dir=export_root,
                                       skipped_funnels=skipped_funnels)
            if output_config is None:
                if len(annotation_bundle.query_columns) > ldscore_config.query_batch_size:
                    raise LDSCUsageError("Multiple query execution batches require output_config with an output directory. Use one batch for a small in-memory calculation.")
                if annotation_bundle.query_preparation is not None:
                    raise LDSCUsageError("In-memory LD-score calculation requires already-prepared annotations. Use AnnotationBundle.from_frames() to prepare without disk staging.")
                batch = producer.prepare(dataclass_replace(annotation_bundle, diagnostics_in_memory=True))
                outcomes = run_chromosome_batch(producer, pool, worker_count, batch, "batch 1/1")
                result, _ = producer.reduce(outcomes, batch)
                result = finalize_batch(result, tuple(batch.bundle.query_statuses), ())
                _warn_ld_window_snp_subset(result.ld_window_snp_counts, "LDScoreResult.ld_window_snp_counts")
                log_reference_snp_funnel_warning(result.reference_snp_funnel, "LDScoreResult.reference_snp_funnel")
                _log_annotation_classification(result)
                return result

            results = run_query_batches(
                producer, annotation_bundle.query_columns, tuple(annotation_bundle.query_statuses),
                width=ldscore_config.query_batch_size, pool=pool, workers=worker_count,
                scratch_root=output_config.output_dir)
            result = LDScoreDirectoryWriter().write_batches(results, output_config,
                query_status_order=[status.query for status in annotation_bundle.query_statuses],
                ref_metadata_dir=export_root)
            _log_annotation_classification(result)
            return result

    @serial_calculation
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

        ``ref_panel.prepare_chromosome`` aligns annotations to the retained
        reference SNPs ``B_chrom & A'_chrom`` and owns the reader until the
        kernel returns or raises. The same aligned annotation values feed
        numerical projection, counts, overlap, and result provenance.
        ``regression_snps`` defines weight contributions and later restricts
        the normalized output rows.
        Standalone calls use one worker with 1 BLAS thread; calls inside an
        existing calculation worker reuse that process.
        """
        if isinstance(annotation_bundle, AnnotationResult):
            raise LDSCUsageError(_annotation_result_message("compute_chromosome()"))
        backend = getattr(getattr(ref_panel, "spec", None), "backend", None)
        chrom_annotations = _chromosome_annotations(annotation_bundle,chrom)
        LOGGER.info(f"Computing chromosome {chrom} LD scores with backend '{backend or 'unknown'}' and {len(chrom_annotations.metadata)} annotation rows.")
        with ref_panel.prepare_chromosome(chrom, chrom_annotations, ldscore_config) as prepared:
            funnel = prepared.reference_snp_funnel
            if funnel is not None:
                LOGGER.info("Chromosome %s reference-SNP funnel: %s.", chrom, funnel_summary(funnel))
            chrom_result = kernel_ldscore.compute_chromosome(
                chrom, prepared, snp_identifier=global_config.snp_identifier,
                snp_batch_size=ldscore_config.snp_batch_size,
                common_maf_min=ldscore_config.common_maf_min,
                regression_keys=regression_snps, regression_regions=regression_regions,
                query_batch_size=ldscore_config.query_batch_size,
            )
            if export_dir is not None:
                _write_one_ref_metadata_sidecar(prepared.metadata, chrom, export_dir)
        drops = _stage_chromosome_drops(annotation_bundle, chrom, chrom_result.identity_drops)
        result = self._wrap_chrom_result(
            chrom_result,
            global_config=global_config,
            annotation_types=chrom_result.annotation_types,
            identity_drops=drops,
        )
        LOGGER.info(f"Finished chromosome {chrom} with {len(result.baseline_table)} retained SNP rows.")
        return dataclass_replace(result, reference_snp_funnel=funnel)

    def _wrap_chrom_result(
        self,
        chrom_result: kernel_ldscore.ChromComputationResult,
        global_config: GlobalConfig,
        annotation_types: dict[str, str] | None = None,
        identity_drops: TsvDiagnostics | None = None,
    ) -> ChromLDScoreResult:
        """Convert one kernel chromosome result into the typed public result.

        The kernel's metadata is already a private frame with a RangeIndex. The
        baseline and query tables are column subsets of one sorted table, so
        their rows match by construction and need no row comparison.
        """
        reference_metadata = chrom_result.metadata
        ld_scores = pd.DataFrame(chrom_result.ld_scores, columns=list(chrom_result.ldscore_columns))
        metadata_columns = ["CHR", "SNP", "POS", *[column for column in ("A1", "A2") if column in reference_metadata.columns]]
        regression_weights = np.asarray(chrom_result.w_ld, dtype=np.float32).reshape(-1)
        ldscore_table = pd.concat(
            [
                reference_metadata.loc[:, metadata_columns],
                ld_scores,
                pd.DataFrame({REGRESSION_LD_SCORE_COLUMN: regression_weights}, copy=False),
            ],
            axis=1,
            copy=False,
        )
        ldscore_table = sort_frame_by_genomic_position(ldscore_table)
        count_map = {
            "all_reference_snp_counts": np.asarray(chrom_result.M, dtype=np.float64),
            "common_reference_snp_counts": np.asarray(chrom_result.M_5_50, dtype=np.float64),
        }
        baseline_table, query_table = _split_ldscore_table(
            ldscore_table,
            baseline_columns=list(chrom_result.baseline_columns),
            query_columns=list(chrom_result.query_columns),
            snp_identifier=global_config.snp_identifier,
        )
        result = ChromLDScoreResult(
            chrom=str(chrom_result.chrom),
            baseline_table=baseline_table,
            query_table=query_table,
            baseline_columns=list(chrom_result.baseline_columns),
            query_columns=list(chrom_result.query_columns),
            snp_count_totals=count_map,
            config_snapshot=global_config,
            overlap=getattr(chrom_result, "overlap", None),
            reference_snp_count=chrom_result.reference_snp_count,
            regression_selected_snp_count=chrom_result.regression_selected_snp_count,
            regression_region_removed_snp_count=chrom_result.regression_region_removed_snp_count,
            identity_drops=identity_drops,
            annotation_types=dict(annotation_types or {}),
        )
        return result

    def _aggregate_chromosome_results(
        self,
        chromosome_results: Sequence[ChromLDScoreResult | DirectChromosomeFragment],
        count_config: dict[str, Any] | None = None,
    ) -> LDScoreResult | StagedLDScoreBatch:
        """Combine chromosome statistics, retaining query paths for writing runs.

        The annotation overlap matrix is aggregated only for partitioned runs (two
        or more LD-score columns). Single-annotation runs (e.g. the synthetic
        unpartitioned ``base``) leave ``LDScoreResult.overlap`` as ``None`` so the
        downstream writer emits no redundant overlap artifact.
        """
        if not chromosome_results:
            raise LDSCInternalError(
                "LD-score calculation could not assemble its chromosome results because none were returned. "
                "This is an unexpected internal state; its cause has not been established. "
                "Report the command, package version, and saved traceback. "
                "--log-level DEBUG adds developer context if the run needs to be reproduced."
            )
        snapshots = [result.config_snapshot for result in chromosome_results if result.config_snapshot is not None]
        if snapshots:
            for snapshot in snapshots[1:]:
                validate_config_compatibility(
                    snapshots[0],
                    snapshot,
                    context="ChromLDScoreResult aggregation",
                )

        # Later direct batches carry query statistics only; batch 1 supplied the baseline.
        include_baseline = chromosome_results[0].baseline_table is not None
        baseline_columns = list(chromosome_results[0].baseline_columns)
        query_columns = list(chromosome_results[0].query_columns)
        counted_baseline = baseline_columns if include_baseline else []
        count_keys = sorted(set.intersection(*(set(result.snp_count_totals) for result in chromosome_results)))
        count_totals = {
            key: np.sum(np.vstack([result.snp_count_totals[key] for result in chromosome_results]), axis=0)
            for key in count_keys
        }
        from .overlap_table import LDScoreOverlap

        # The overlap (Gram) matrix is only consumed by partitioned-h2 and the
        # design-matrix collinearity check, both of which need >=2 LD-score
        # columns. A single-column run (e.g. the synthetic unpartitioned `base`)
        # yields a degenerate 1x1 overlap equal to a SNP count already in
        # metadata, so the artifact is suppressed rather than written redundantly.
        n_ld_columns = len(baseline_columns) + len(query_columns)
        overlaps = [result.overlap for result in chromosome_results if result.overlap is not None]
        aggregated_overlap = None
        if n_ld_columns >= 2 and len(overlaps) == len(chromosome_results):
            aggregated_overlap = LDScoreOverlap.from_contribution(
                sum_overlap_contributions(overlaps), baseline_columns, query_columns,
                include_baseline=include_baseline,
            )
        annotation_types = combine_annotation_types(
            [*counted_baseline, *query_columns], [result.annotation_types for result in chromosome_results])
        shared = dict(
            baseline_columns=baseline_columns,
            query_columns=query_columns,
            count_records=_count_records_from_totals(
                baseline_columns=counted_baseline,
                query_columns=query_columns,
                count_totals=count_totals,
            ),
            count_config=dict(count_config or {}),
            config_snapshot=snapshots[0] if snapshots else None,
            overlap=aggregated_overlap,
            identity_drops_by_chrom={
                result.chrom: result.identity_drops
                for result in chromosome_results if result.identity_drops is not None
            },
            annotation_types=annotation_types,
        )
        if isinstance(chromosome_results[0], DirectChromosomeFragment):
            baseline_table = None
            if include_baseline:
                # Each worker sorted its chromosome; placing chromosomes in order
                # gives the genomic order that the query fragments below follow.
                baseline_table = concat_chromosome_tables(
                    [(result.chrom, result.baseline_table) for result in chromosome_results])
            ordered = sorted(chromosome_results, key=lambda result: chrom_sort_key(result.chrom))
            return StagedLDScoreBatch(
                baseline_table=baseline_table,
                query_fragments=tuple(result.query_fragment for result in ordered if result.query_fragment is not None),
                query_statuses=(), gene_list_batch=None, chromosome_scope={}, index_provenance={}, **shared,
            )
        baseline_table = concat_chromosome_tables([(result.chrom, result.baseline_table) for result in chromosome_results])
        query_table = None
        if query_columns:
            query_table = concat_chromosome_tables([(result.chrom, result.query_table) for result in chromosome_results])
        result = LDScoreResult(
            baseline_table=baseline_table,
            query_table=query_table,
            **shared,
        )
        result.validate()
        return result


@dataclass(frozen=True)
class _DirectBatch:
    """One direct query batch: its prepared annotations and private fragment directory."""

    bundle: Any
    source_drops: dict
    ordinal: int = 0
    fragment_dir: str | None = None


class _DirectProducer:
    """Direct LD-score producer: each task computes one chromosome from the reference panel.

    A batch projects its active query annotations, workers compute and stage
    each chromosome's LD scores, and ``reduce`` combines them with the run's
    direct provenance. The no-output path reuses ``prepare``, the tasks, and
    ``reduce`` for its single materialized batch. ``skipped_funnels`` belong to
    chromosomes skipped before the bundle was built, which have no task.
    """

    batch_name, phase, task = "batch", "LD-score calculation", "LD-score computation"

    def __init__(self, calculator, annotation_bundle, ref_panel, ldscore_config, global_config,
                 regression_regions, log_paths, *, output_dir=None, export_dir=None, skipped_funnels=()):
        self.calculator = calculator
        self.annotation_bundle = annotation_bundle
        self.ref_panel = ref_panel
        self.ldscore_config = ldscore_config
        self.global_config = global_config
        self.regression_regions = regression_regions
        self.log_paths = log_paths or {}
        self.logger = LOGGER
        self.output_dir = output_dir
        self.export_dir = export_dir
        self.skipped_funnels = list(skipped_funnels)
        self.chromosomes = _chromosomes_from_bundle(annotation_bundle)
        self.prefixes = ref_panel._resolved_prefixes() if isinstance(ref_panel, PlinkRefPanel) else None

    def prepare(self, bundle, *, ordinal=0, fragment_dir=None) -> _DirectBatch:
        """Check the bundle against the run's identity and partition its source identity drops."""
        global_config = self.global_config
        if bundle.config_snapshot is not None:
            validate_config_compatibility(
                (dataclass_replace(bundle.config_snapshot, genome_build=global_config.genome_build)
                 if global_config.genome_build in {"hg19", "hg38"} else bundle.config_snapshot),
                global_config,
                context="AnnotationBundle and LDScoreCalculator runtime config",
            )
        # Identity drops are published from the first batch only.
        source_drops = _partition_annotation_drops(bundle) if ordinal == 0 else {}
        return _DirectBatch(dataclass_replace(bundle, chromosome_identity_drops=source_drops), source_drops,
                            ordinal, fragment_dir)

    @contextmanager
    def open_batch(self, ordinal, query_columns, statuses, scratch):
        """Project the batch's BED or gene-list queries and release them after consumption."""
        with execution_query_bundle(self.annotation_bundle, query_columns, statuses, self.output_dir) as active:
            yield self.prepare(active, ordinal=ordinal, fragment_dir=str(scratch))

    def submit(self, pool, chrom, batch, label):
        """Submit one chromosome's computation with its slice of the batch's annotations.

        Only the first batch exports reference metadata, and only the first
        batch of a staged (writing) run computes the baseline; later staged
        batches compute their query scores and statistics alone.
        """
        first = batch.ordinal == 0
        return pool.submit(_compute_one_chromosome, chrom,
            _slice_annotation_bundle(batch.bundle, chrom), self.ref_panel.spec,
            self.ldscore_config, self.global_config, self.export_dir if first else None,
            chromosome_prefixes=self.prefixes, log_path=self.log_paths.get(chrom), batch_label=label,
            fragment_dir=batch.fragment_dir, include_baseline=first or batch.fragment_dir is None)

    def reduce(self, outcomes, batch):
        """Combine chromosome outcomes; skipped chromosomes contribute only their funnels.

        A chromosome emptied by the annotation SNP grid fails the run; one
        emptied by an explicit restriction or genotype QC is skipped with a
        WARNING. Returns the staged or materialized batch and, for staged
        batches, each computed chromosome's query variation.
        """
        ref_panel, ldscore_config = self.ref_panel, self.ldscore_config
        chromosome_results, funnels, grid_emptied = [], list(self.skipped_funnels), []
        for chrom in self.chromosomes:
            outcome = outcomes[chrom]
            funnel = outcome.reference_snp_funnel if outcome.skipped else outcome.result.reference_snp_funnel
            funnels.append(funnel)
            if outcome.skipped:
                if funnel is not None and annotation_snp_grid_emptied(funnel):
                    grid_emptied.append(funnel)
                else:
                    LOGGER.warning("Skipping chromosome %s: %s", chrom, outcome.skip_message)
                continue
            chromosome_results.append(outcome.result)
        if grid_emptied:
            raise LDSCInputError(_grid_emptied_message(grid_emptied))
        if not chromosome_results:
            raise LDSCInputError(
                "ldscore could not compute any chromosome results after intersecting annotations "
                "with the reference panel. Most likely the annotation SNP identities, genome build, "
                "or allele-aware identifier mode do not match the reference panel. Use matching "
                "annotation and reference-panel artifacts, or rerun with the correct "
                "`--snp-identifier` and `--genome-build`. "
                f"Other causes & fixes: {_LDSCORE_INTERSECTION_DOC}"
            )
        result = self.calculator._aggregate_chromosome_results(
            chromosome_results,
            count_config=common_count_config(ldscore_config.common_maf_min),
        )
        result = dataclass_replace(
            result,
            identity_drops_by_chrom={**batch.source_drops, **result.identity_drops_by_chrom},
            gene_list_batch=getattr(batch.bundle, "gene_list_batch", None),
            chromosome_scope=dict(batch.bundle.source_summary.get("chromosome_scope", {})),
            ld_window={**ldscore_config.ld_window.record(), "whole_chromosome_ok": ldscore_config.whole_chromosome_ok},
        )
        if batch.ordinal == 0:  # every batch shares the reference and regression SNPs; batch 1 records them
            result = dataclass_replace(
                result,
                snp_universe_policy=_snp_universe_policy(
                    ref_panel=ref_panel,
                    ldscore_config=ldscore_config,
                    regression_regions=self.regression_regions,
                    chromosome_results=chromosome_results,
                ),
                ld_window_snp_counts=_ld_window_snp_counts(ref_panel, ldscore_config, chromosome_results),
                reference_snp_funnel=concat_funnels(funnels),
            )
        if batch.ordinal == 0 and self.output_dir is not None:  # every batch shares chromosomes and counts
            _warn_ld_window_snp_subset(
                result.ld_window_snp_counts,
                str(Path(self.output_dir) / "diagnostics" / "ld_window_snp_counts.tsv"))
            log_reference_snp_funnel_warning(result.reference_snp_funnel, reference_snp_funnel_path(self.output_dir))
        LOGGER.debug(f"Computed LD scores for {len(chromosome_results)} chromosomes.")
        staged = isinstance(result, StagedLDScoreBatch)
        return result, [item.variation for item in chromosome_results] if staged else []


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
            f"SNP identifier mode '{snp_identifier}' requires A1/A2 columns, but the LD-score "
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
    if not isinstance(ldscore_table.index, pd.RangeIndex) or ldscore_table.index.start != 0:
        ldscore_table = ldscore_table.reset_index(drop=True)
    # Both tables are column selections (copies) of one frame, so their rows match.
    baseline_table = ldscore_table.loc[:, baseline_order]
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
        query_table = ldscore_table.loc[:, query_order]
    return baseline_table, query_table


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


def build_parser() -> argparse.ArgumentParser:
    """Build the feature parser for LD-score calculation."""
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.prog = 'ldsc ldscore'
    parser.formatter_class = CLIHelpFormatter
    parser.description = (
        'Compute LD scores from a reference panel by following Steps 1-6, or assemble gene-list LD scores from a '
        "prebuilt gene LD-score index (see 'Alternative: gene LD-score index')."
    )
    parser.usage = (
        '%(prog)s [--snp-identifier MODE] --genome-build BUILD\n'
        '    [--baseline-annot-sources SOURCES\n'
        '      [--query-annot-bed-sources SOURCES | --query-annot-gene-list-sources SOURCES]]\n'
        '    (--r2-matrix-dir DIR | --plink-prefix PREFIX)\n'
        '    (--ld-wind-cm CM | --ld-wind-kb KB | --ld-wind-snps N)\n'
        '    --output-dir DIR [options]\n'
        '   or: %(prog)s --gene-ldscore-index-dir DIR --query-annot-gene-list-sources SOURCES\n'
        '    --output-dir DIR [options]'
    )
    identity =parser.add_argument_group('Step 1: SNP identity and genome build')
    baseline = parser.add_argument_group('Step 2: Baseline annotations', description='Required when Step 3 adds query annotations. Omit Steps 2 and 3 to compute unpartitioned LD scores.')
    queries = parser.add_argument_group('Step 3: Query annotations (optional; choose at most one source type)')
    panel = parser.add_argument_group('Step 4: Reference panel (choose one: R2 matrix or PLINK)')
    plink = parser.add_argument_group('Step 4 continued: PLINK reference panel only (optional)')
    window = parser.add_argument_group('Step 5: LD window (choose one)', description='Specify exactly one window; it must fit within the window of any R2 matrix.')
    output = parser.add_argument_group('Step 6: Output')
    snps = parser.add_argument_group('Optional: regression SNPs and reference SNPs', description='By default, the regression SNPs are the HM3 map SNPs outside MHC and centromeres.')
    indexed = parser.add_argument_group('Alternative: gene LD-score index', description='Use instead of Steps 1, 2, 4, and 5 and the optional SNP settings, which the index stores. Supply --query-annot-gene-list-sources (Step 3) and --output-dir (Step 6); omit --padding-bp, --gene-coordinate-file, and --gene-exclude-regions.')
    computing = parser.add_argument_group('Computing resources')
    logs = parser.add_argument_group('Logging')

    identity.add_argument(
        '--snp-identifier', default='chr_pos_allele_aware', choices=('rsid', 'rsid_allele_aware', 'chr_pos', 'chr_pos_allele_aware'),
        help=(
            'Match SNPs by rsID (rsid) or chromosome and position (chr_pos); allele-aware variants also use '
            'A1/A2. Default: chr_pos_allele_aware; usable alleles are required. With '
            '--gene-ldscore-index-dir, inherit the stored identity and omit this flag.'
        ),
    )
    identity.add_argument(
        '--genome-build', choices=('auto', 'hg19', 'hg37', 'GRCh37', 'hg38', 'GRCh38'), default=None,
        help=(
            'Genome build for coordinate matching and gene intervals. Required for coordinate-family modes; direct '
            'gene-list runs default to auto. Auto checks all coordinate inputs and shards; concrete builds '
            'must agree with input declarations. Other rsid runs need '
            'no build. Gene-list runs also use this build for --regr-snps-exclude-regions. With '
            '--gene-ldscore-index-dir, inherit its build and omit this flag.'
        ),
    )

    baseline.add_argument(
        '--baseline-annot-sources', default=None, metavar='SOURCES',
        help=(
            "Baseline annotation tables: comma-separated paths, '*' or '@' patterns. Required with direct query "
            'inputs. If omitted without queries, use a single all-SNP baseline. Cannot be combined with '
            '--gene-ldscore-index-dir. '
            + CHROMOSOME_PATH_HELP
        ),
    )

    query_group = queries.add_mutually_exclusive_group()
    query_group.add_argument(
        '--query-annot-bed-sources', default=None, metavar='SOURCES',
        help=(
            "BED files whose intervals define query SNPs; accepts comma-separated paths or '*' patterns. Requires "
            '--baseline-annot-sources; cannot be combined with --query-annot-sources, '
            "--query-annot-gene-list-sources, or --gene-ldscore-index-dir. If omitted, add no BED queries. '@' is "
            'not expanded; see --baseline-annot-sources for wildcard and quoting rules.'
        ),
    )
    query_group.add_argument(
        '--query-annot-gene-list-sources', default=None, metavar='SOURCES',
        help=(
            'One-column gene-list files, separated by commas. Direct mode requires --baseline-annot-sources, '
            '--gene-coordinate-file, and --padding-bp; indexed mode requires --gene-ldscore-index-dir. Cannot be '
            'combined with --query-annot-sources or --query-annot-bed-sources. If omitted, add no gene-list '
            "queries. Accepts '*' patterns (see --baseline-annot-sources); '@' is not expanded."
        ),
    )
    query_group.add_argument(
        '--query-annot-sources', default=None, metavar='SOURCES',
        help=(
            "Query annotation tables: comma-separated paths, '*' or '@' patterns (see --baseline-annot-sources for "
            'pattern rules). Requires --baseline-annot-sources; cannot be combined with --query-annot-bed-sources, '
            '--query-annot-gene-list-sources, or --gene-ldscore-index-dir. If omitted, add no table-based '
            'queries.'
        ),
    )
    queries.add_argument(
        '--padding-bp', type=int, default=None, metavar='BP',
        help=(
            'Extend BED intervals or gene boundaries by this many base pairs at both ends, then find covered '
            'SNPs. Omission means 0 for BED input; direct gene-list input requires an explicit value, including '
            '0. Cannot be used with prebuilt annotations, no query, or --gene-ldscore-index-dir. Set to 0 if your '
            'BED intervals are already padded to avoid double padding.'
        ),
    )
    queries.add_argument(
        '--gene-coordinate-file', default=None, metavar='FILE',
        help=(
            'Headered TSV/TSV.GZ of gene coordinates, using one-based inclusive intervals. Required with direct '
            '--query-annot-gene-list-sources; indexed runs use their stored catalog and reject this flag. If '
            'omitted outside gene-list mode, no catalog is used. Supply an exact file path; patterns are not '
            'expanded.'
        ),
    )
    queries.add_argument(
        '--control-gene-list-file', default=None, metavar='FILE',
        help=(
            'One-column gene list added as a fixed control alongside each gene-list query. Requires '
            '--query-annot-gene-list-sources. If omitted, add no gene control. Supply an exact file path; '
            'patterns are not expanded.'
        ),
    )
    queries.add_argument(
        '--gene-exclude-regions', choices=('none', 'mhc'), default='none',
        help=(
            'Exclude genes overlapping MHC before extending their intervals; none excludes no genes. '
            'Default: none for direct gene-list input. With --gene-ldscore-index-dir, use the stored policy '
            'and omit this flag.'
        ),
    )
    queries.add_argument(
        '--allow-unresolved-genes', action='store_const', const='resolved-only',
        dest='gene_list_resolution_policy', default='strict',
        help=(
            'Omit rejected gene identifiers and continue with the usable subset, recording omissions '
            'in diagnostics. Default: off (strict resolution; rejected identifiers stop the run). '
            'Requires --query-annot-gene-list-sources. Other validation errors remain fatal.'
        ),
    )

    panel.add_argument(
        '--r2-matrix-dir', dest='r2_matrix_dir', default=None, metavar='DIR',
        help=(
            'R2 matrix directory written by build-r2-matrix, containing chrN_r2.parquet and required chrN_meta.tsv.gz '
            'files. The LD window must use the unit the R2 matrix was built with and be no wider. Direct mode '
            'requires exactly one of --r2-matrix-dir or --plink-prefix. Cannot be combined with '
            '--gene-ldscore-index-dir; no default.'
        ),
    )
    panel.add_argument(
        '--plink-prefix', default=None, metavar='PREFIX',
        help=(
            'PLINK .bed/.bim/.fam filename prefix or chromosome suite stem. Direct mode requires exactly one of '
            '--plink-prefix or --r2-matrix-dir. Cannot be combined with --gene-ldscore-index-dir; no default. '
            "Direct query runs require every declared chromosome member. Accepts '*' and '@' patterns; see "
            '--baseline-annot-sources for pattern rules.'
        ),
    )

    plink.add_argument(
        '--keep-indivs-file', default=None, metavar='FILE',
        help=(
            'File listing one individual ID (IID) per row to retain from PLINK genotypes. Requires '
            '--plink-prefix; cannot be used with --r2-matrix-dir or --gene-ldscore-index-dir. If omitted, retain all '
            'individuals. Supply an exact file path; patterns are not expanded.'
        ),
    )
    plink.add_argument(
        '--genetic-map-hg19-sources', default=None, metavar='SOURCES',
        help=(
            'hg19 genetic-map files used to derive centiMorgan coordinates for PLINK SNPs. Accepts '
            "comma-separated exact file paths; '*' and '@' are not expanded. If omitted, use informative PLINK CM "
            'values; cM windows require one of these sources. Direct PLINK mode only.'
        ),
    )
    plink.add_argument(
        '--genetic-map-hg38-sources', default=None, metavar='SOURCES',
        help=(
            'hg38 genetic-map files used to derive centiMorgan coordinates for PLINK SNPs. Accepts '
            "comma-separated exact file paths; '*' and '@' are not expanded. If omitted, use informative PLINK CM "
            'values; cM windows require one of these sources. Direct PLINK mode only.'
        ),
    )
    plink.add_argument(
        '--export-ref-metadata', default=False, action='store_true',
        help=(
            'Write reference-panel SNP metadata, with the minor allele as A1, to ref_metadata/chrN_meta.tsv.gz '
            'under --output-dir. Applies to --plink-prefix only; rejected with --r2-matrix-dir, whose '
            'chrN_meta.tsv.gz sidecars hold the same metadata, and with --gene-ldscore-index-dir. Default: off; '
            'these additional files are not written.'
        ),
    )

    # Not required: indexed mode takes its window from the index.
    window_choice = window.add_mutually_exclusive_group()
    window_choice.add_argument(
        '--ld-wind-cm', default=None, type=float, metavar='CM',
        help=(
            'Maximum LD-window distance on either side of each SNP, in centiMorgans; requires usable '
            'genetic-map coordinates. No default; specify exactly one of --ld-wind-snps, --ld-wind-kb, or '
            '--ld-wind-cm. Direct mode only.'
        ),
    )
    window_choice.add_argument(
        '--ld-wind-kb', default=None, type=float, metavar='KB',
        help=(
            'Maximum LD-window distance on either side of each SNP, in kilobases. No default; specify '
            'exactly one of --ld-wind-snps, --ld-wind-kb, or --ld-wind-cm. Direct mode only.'
        ),
    )
    window_choice.add_argument(
        '--ld-wind-snps', default=None, type=int, metavar='N',
        help=(
            'Maximum LD-window distance on either side of each SNP, counted in retained SNPs with '
            '--plink-prefix and in R2 matrix rows with --r2-matrix-dir. No default; specify exactly one of '
            '--ld-wind-snps, --ld-wind-kb, or --ld-wind-cm. Direct mode only.'
        ),
    )
    window.add_argument(
        '--yes-really', default=False, action='store_true',
        help=(
            'Allow an LD window covering an entire chromosome, bypassing the usual guard. Default: off; '
            'reject whole-chromosome windows. Direct mode only.'
        ),
    )

    output.add_argument(
        '--output-dir', required=True, metavar='DIR', type=output_directory_argument,
        help=(
            'Required destination for LD-score tables, SNP counts, and diagnostics.'
        ),
    )
    output.add_argument(
        '--overwrite', action='store_true', default=False,
        help=(
            "Replace this command's existing output files and remove obsolete outputs from an earlier run. "
            'Default: off; stop if output files already exist.'
        ),
    )

    snps.add_argument(
        '--regr-snps-file', default=None, metavar='FILE',
        help=(
            'Headered SNP identity table selecting the regression SNPs: stored output rows and regression weights. If '
            'omitted, use bundled HapMap3 SNPs. --regr-snps-exclude-regions then applies. Direct mode only; '
            'rejected with --gene-ldscore-index-dir. '
            + SCALAR_PATH_HELP
        ),
    )
    snps.add_argument(
        '--regr-snps-exclude-regions', choices=REGR_SNPS_EXCLUDE_REGIONS_CHOICES, default='mhc-and-centromeres',
        help=(
            'Exclude MHC, centromeres, both, or neither from regression SNPs after --regr-snps-file '
            'or HapMap3 selection. Default: mhc-and-centromeres. Reference SNPs and annotation '
            'counts are unchanged. Direct mode only; indexed runs inherit this setting.'
        ),
    )
    snps.add_argument(
        '--ref-panel-snps-file', default=None, metavar='FILE',
        help=(
            'Headered SNP identity table restricting the reference SNPs contributing to LD scores. If omitted, '
            'use all otherwise eligible reference SNPs. Extra columns are ignored and duplicate keys collapse. '
            'Direct mode only; rejected with --gene-ldscore-index-dir. Supply an exact file path; patterns are '
            'not expanded.'
        ),
    )
    snps.add_argument(
        '--maf-min', default=None, type=float, metavar='VALUE',
        help=(
            'Keep reference SNPs with MAF at or above this value, in [0, 0.5]. If omitted, apply no '
            'additional MAF threshold. Direct mode only; indexed runs inherit their reference selection.'
        ),
    )
    snps.add_argument(
        '--common-maf-min', default=0.05, type=float, metavar='VALUE',
        help=(
            'Minimum MAF for common-SNP annotation counts, including the boundary. Default: 0.05; does not '
            'filter reference SNPs. Direct mode only; indexed runs inherit their count threshold.'
        ),
    )

    indexed.add_argument(
        '--gene-ldscore-index-dir', default=None, metavar='DIR',
        help=(
            'Use a complete prebuilt gene LD-score index instead of a live reference panel. Requires '
            '--query-annot-gene-list-sources; cannot be combined with direct-input or direct-computation '
            'flags. If omitted, compute LD scores directly.'
        ),
    )

    computing.add_argument(
        '--workers', default=1, type=parse_workers, metavar='N',
        help='Maximum simultaneous chromosomes in direct or indexed mode. ' + WORKERS_HELP,
    )
    computing.add_argument(
        '--query-batch-size', default=100, type=int, metavar='N',
        help=(
            'Maximum number of query annotations processed at once. Default: 100; use a smaller positive '
            'value to reduce memory. Each completed batch has its own query file; direct mode repeats reference traversal.'
        ),
    )
    computing.add_argument(
        '--snp-batch-size', default=DEFAULT_SNP_BATCH_SIZE, type=int, metavar='N',
        help=(
            f'Number of reference SNPs processed per PLINK genotype batch. Default: {DEFAULT_SNP_BATCH_SIZE}; '
            'larger batches use more memory. LD scores do not depend on it. Ignored with --r2-matrix-dir; '
            'rejected with --gene-ldscore-index-dir.'
        ),
    )

    logs.add_argument(
        '--log-level', default='INFO', choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'),
        help=LOG_LEVEL_HELP,
    )
    return parser


@dataclass(frozen=True)
class _DirectConfigs:
    """The configuration of one direct ``ldscore`` run, built before anything is written.

    ``global_config`` holds the requested genome build, which the run resolves
    against its inputs (``_resolve_build``): ``"auto"`` or a concrete build in
    coordinate-family modes. ``genome_build`` is ``--genome-build`` as given,
    normalized, or ``None``. rsID-family identity records no build, but
    gene-list projection still needs one.
    """

    global_config: GlobalConfig
    genome_build: str | None
    annotation_config: AnnotationBuildConfig
    ref_panel_config: RefPanelConfig
    ldscore_config: LDScoreConfig
    output_config: LDScoreOutputConfig
    regr_snps_exclude_regions: str
    emit_console_notices: bool

    @property
    def has_queries(self) -> bool:
        """Return whether the run adds prebuilt, BED, or gene-list query annotations."""
        config = self.annotation_config
        return bool(config.query_annot_sources or config.query_annot_bed_sources
                    or config.query_annot_gene_list_sources)

    @property
    def has_gene_lists(self) -> bool:
        """Return whether the run projects gene-list query annotations."""
        return bool(self.annotation_config.query_annot_gene_list_sources)


@dataclass(frozen=True)
class _IndexedConfigs:
    """One indexed ``ldscore`` run, checked before anything is written.

    The gene LD-score index supplies the SNP identity, genome build, reference
    and regression SNPs, and LD window; the run adds its gene lists and
    execution settings.
    """

    index_dir: str
    query_gene_list_sources: tuple[str, ...]
    control_gene_list_file: str | None
    gene_list_resolution_policy: str
    query_batch_size: int
    workers: int
    output_dir: str
    overwrite: bool
    log_level: str
    emit_console_notices: bool


@materializing_overwrite_guard(
    args_output_scope,
    command="run_ldscore_from_args(...)",
)
def run_ldscore_from_args(args: argparse.Namespace) -> "LDScoreSource":
    """Run LD-score calculation from a parsed CLI namespace.

    The workflow resolves unified path tokens for baseline annotations, optional
    prebuilt, BED, or gene-list queries, and the reference panel. If no baseline
    or query annotations are supplied, it synthesizes an all-ones ``base``
    annotation over the retained reference-panel metadata. Before calculation
    it preflights the complete baseline/query/overlap artifact family, root
    metadata, applicable diagnostics, and ``diagnostics/ldscore.log`` under
    ``output_dir``. Every run writes ``diagnostics/reference_snp_funnel.tsv``;
    reference-panel SNPs without an annotation row are excluded and reported
    as one log WARNING, which a successful command-line run also prints to
    stderr. Each execution batch is privately written and released;
    successful publication records ordered files and chromosome row groups in
    ``metadata.json.query_batches``. A single batch uses
    ``ldscore.query.parquet``; multiple batches use numbered query files.
    Successful authorized overwrites remove stale owned batch files.

    Direct batches repeat reference traversal. Indexed runs reuse one
    chromosome operator per worker through its query batches. Both modes cap
    chromosome workers at the chromosome count. Scientific restrictions and
    count/overlap definitions are independent of these runtime settings.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed LD-score command arguments, including an output directory and
        the record of explicitly supplied options that ``parse_command_line``
        and ``namespace_from_kwargs`` keep (``_explicit_cli_options``).

    Returns
    -------
    LDScoreSource
        Validated saved-artifact handle with shared baseline values and
        explicit ``read_queries(names)`` access. No query LD-score tables
        remain resident in the returned object. ``output_paths`` contains
        scientific artifacts only.
    """
    configs = _configs_from_args(args)
    if isinstance(configs, _IndexedConfigs):
        return _run_explicit_indexed_ldscore(configs)
    return _run_direct_ldscore(configs)


# Range errors that a configuration words by field name, worded with the option a command-line user typed.
_OPTION_RANGE_ERRORS = {
    "maf_min": (
        "ldscore received invalid `--maf-min={value}`. Most likely the minor-allele frequency threshold was "
        "entered outside the valid [0, 0.5] range. Pass a value between 0 and 0.5, or omit the option."
    ),
    "common_maf_min": (
        "ldscore received invalid `--common-maf-min={value}`. Most likely the common-SNP MAF threshold was "
        "entered outside the valid [0, 0.5] range. Pass a value between 0 and 0.5."
    ),
    "snp_batch_size": (
        "ldscore received invalid `--snp-batch-size={value}`. Most likely the parquet query batch size was set "
        "to zero or a negative value. Pass a positive integer batch size."
    ),
}


def _config_from_options(config_type, **fields):
    """Build ``config_type(**fields)``; a range error on a field in ``_OPTION_RANGE_ERRORS`` names its option."""
    try:
        return config_type(**fields)
    except LDSCConfigError as exc:
        for name, message in _OPTION_RANGE_ERRORS.items():
            if name in fields and str(exc).startswith(f"Could not construct {config_type.__name__}: {name}="):
                raise LDSCConfigError(message.format(value=fields[name])) from None
        raise


def _configs_from_args(args: argparse.Namespace) -> "_DirectConfigs | _IndexedConfigs":
    """Check the input-mode rules and build every configuration of one ``ldscore`` run.

    Runs before anything is written, for the command line and for
    ``run_ldscore()`` alike. An option counts as supplied when the caller
    named it (``_explicit_cli_options``), whatever its value.
    """
    explicit_options = args._explicit_cli_options
    _validate_padding_usage(args)
    _validate_gene_list_mode_args(args, explicit_options)
    emit_console_notices = bool(getattr(args, "_emit_query_console_notices", False))
    if args.gene_ldscore_index_dir is not None:
        return _indexed_configs(args, explicit_options, emit_console_notices)
    gene_lists = tuple(split_cli_path_tokens(args.query_annot_gene_list_sources))
    if not _has_cli_tokens(args.baseline_annot_sources) and (
        gene_lists or _has_cli_tokens(args.query_annot_sources) or _has_cli_tokens(args.query_annot_bed_sources)
    ):
        raise LDSCUsageError(_QUERY_REQUIRES_BASELINE_MESSAGE)
    r2_matrix_dir = normalize_optional_path_token(args.r2_matrix_dir)
    if bool(r2_matrix_dir) == bool(args.plink_prefix):
        raise LDSCUsageError(
            "ldscore could not choose a reference-panel backend. Most likely both R2 matrix "
            "input and PLINK input were supplied, or neither was supplied. Pass exactly "
            "one of `--r2-matrix-dir` for R2 matrix mode or `--plink-prefix` for PLINK mode."
        )
    snp_identifier = normalize_snp_identifier_mode(args.snp_identifier)
    genome_build = normalize_genome_build(args.genome_build)
    if r2_matrix_dir:
        if args.keep_indivs_file:
            raise LDSCUsageError(
                "ldscore cannot apply `--keep-indivs-file` in R2 matrix mode. Most likely "
                "`--keep-indivs-file` was combined with `--r2-matrix-dir`, but individual-level "
                "filtering only exists before PLINK genotype LD calculation. Remove "
                "`--keep-indivs-file`, or rerun in PLINK mode with `--plink-prefix`."
            )
        if identity_mode_family(snp_identifier) == "chr_pos" and genome_build is None:
            raise LDSCUsageError(
                "ldscore cannot run R2 matrix mode with a coordinate-family SNP identifier mode "
                "without a genome build. Most likely `--snp-identifier chr_pos` or "
                "`--snp-identifier chr_pos_allele_aware` was supplied without `--genome-build`. Pass "
                "`--genome-build hg19`, `--genome-build hg38`, or `--genome-build auto`."
            )
    ref_panel_config = _config_from_options(
        RefPanelConfig,
        backend="r2_matrix" if r2_matrix_dir else "plink",
        plink_prefix=None if r2_matrix_dir else args.plink_prefix,
        r2_matrix_dir=r2_matrix_dir or None,
        ref_panel_snps_file=args.ref_panel_snps_file,
        maf_min=args.maf_min,
        keep_indivs_file=args.keep_indivs_file,
        genetic_map_hg19_sources=args.genetic_map_hg19_sources,
        genetic_map_hg38_sources=args.genetic_map_hg38_sources,
    )
    ldscore_config = _config_from_options(
        LDScoreConfig,
        ld_wind_snps=args.ld_wind_snps,
        ld_wind_kb=args.ld_wind_kb,
        ld_wind_cm=args.ld_wind_cm,
        regr_snps_file=args.regr_snps_file,
        snp_batch_size=args.snp_batch_size,
        query_batch_size=args.query_batch_size,
        common_maf_min=args.common_maf_min,
        whole_chromosome_ok=args.yes_really,
        export_ref_metadata=args.export_ref_metadata,
        workers=args.workers,
    )
    global_config = _requested_global_config(
        snp_identifier, genome_build, has_gene_lists=bool(gene_lists), log_level=args.log_level)
    annotation_config = AnnotationBuildConfig(
        baseline_annot_sources=tuple(split_cli_path_tokens(args.baseline_annot_sources)),
        query_annot_sources=tuple(split_cli_path_tokens(args.query_annot_sources)),
        query_annot_bed_sources=tuple(split_cli_path_tokens(args.query_annot_bed_sources)),
        query_annot_gene_list_sources=gene_lists,
        gene_coordinate_file=args.gene_coordinate_file,
        control_gene_list_file=args.control_gene_list_file,
        gene_list_resolution_policy=args.gene_list_resolution_policy,
        gene_exclude_regions=args.gene_exclude_regions,
        padding_bp=0 if args.padding_bp is None else args.padding_bp,
    )
    return _DirectConfigs(
        global_config=global_config,
        genome_build=genome_build,
        annotation_config=annotation_config,
        ref_panel_config=ref_panel_config,
        ldscore_config=ldscore_config,
        output_config=LDScoreOutputConfig(output_dir=normalize_path_token(args.output_dir), overwrite=args.overwrite),
        regr_snps_exclude_regions=args.regr_snps_exclude_regions,
        emit_console_notices=emit_console_notices,
    )


def _requested_global_config(snp_identifier: str, genome_build: str | None, *, has_gene_lists: bool,
                             log_level: str) -> GlobalConfig:
    """Return the run's ``GlobalConfig`` with the requested genome build, before it is resolved.

    Coordinate-family modes need a build; a gene-list run without one infers
    it (``auto``). rsID-family identity records no build in gene-list runs,
    where ``--genome-build`` selects only the gene projection build.
    """
    if identity_mode_family(snp_identifier) == "chr_pos":
        if genome_build is None and not has_gene_lists:
            raise LDSCUsageError(
                "ldscore cannot resolve a coordinate-family SNP identifier mode without a genome build. "
                "Most likely `--snp-identifier` was set to a coordinate-family mode but `--genome-build` "
                "was omitted. Pass `--genome-build auto`, `--genome-build hg19`, or "
                "`--genome-build hg38`."
            )
        return GlobalConfig(snp_identifier=snp_identifier, genome_build=genome_build or "auto", log_level=log_level)
    return GlobalConfig(snp_identifier=snp_identifier, genome_build=None if has_gene_lists else genome_build,
                        log_level=log_level)


def _run_direct_ldscore(configs: _DirectConfigs) -> "LDScoreSource":
    """Run one direct ``ldscore`` request: output preflight, the input gates, then the calculation."""
    from .direct_annotation import prepare_direct_annotations, prepare_synthetic_base

    annotation_config, panel_config = configs.annotation_config, configs.ref_panel_config
    ldscore_config, output_config = configs.ldscore_config, configs.output_config
    output_dir = ensure_output_directory(output_config.output_dir, label="LD-score output directory")
    log_path = output_dir / LDSCORE_LOG
    LDScoreDirectoryWriter.artifact_family(output_dir, include_logs=True).preflight(
        overwrite=output_config.overwrite, additional_paths=[log_path])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    issues_path = output_dir / "diagnostics/input_issues.tsv"
    with AnnotationWorkspace(output_dir) as workspace:
        with workflow_logging("ldscore", log_path, log_level=configs.global_config.log_level):
            from .._refpanel.preflight import inspect_r2_headers
            from .preflight import inspect_reference_inputs
            r2_matrix_dir = panel_config.r2_matrix_dir
            declared_inputs = inspect_annotation_inputs(
                baseline=annotation_config.baseline_annot_sources, query=annotation_config.query_annot_sources,
                files=[("query_annot_bed_sources", annotation_config.query_annot_bed_sources),
                       ("query_annot_gene_list_sources", annotation_config.query_annot_gene_list_sources)],
                scalar_files=[("control_gene_list_file", annotation_config.control_gene_list_file),
                              ("gene_coordinate_file", annotation_config.gene_coordinate_file),
                              ("regr_snps_file", ldscore_config.regr_snps_file),
                              ("ref_panel_snps_file", panel_config.ref_panel_snps_file),
                              ("keep_indivs_file", panel_config.keep_indivs_file)],
                plink=panel_config.plink_prefix, require_complete_suite=configs.has_queries,
                mode=configs.global_config.snp_identifier, issues_path=issues_path,
                checks=[("reference", r2_matrix_dir, lambda: inspect_r2_headers(
                    r2_matrix_dir, configs.genome_build, command="ldscore", ld_window=ldscore_config.ld_window))]
                       if r2_matrix_dir else ())
            require_usable_bed_query_names(declared_inputs.files["query_annot_bed_sources"],
                                           annotation_columns(declared_inputs), issues_path)
            batch = catalog = None
            if configs.has_gene_lists:
                catalog = GeneCatalog.load(annotation_config.gene_coordinate_file)
                batch = resolve_gene_lists_staged(list(annotation_config.query_annot_gene_list_sources), catalog,
                    workspace, control_path=annotation_config.control_gene_list_file,
                    resolution_policy=annotation_config.gene_list_resolution_policy,
                    gene_exclude_regions=annotation_config.gene_exclude_regions)
                batch = validate_gene_query_names(batch, annotation_columns(declared_inputs))
                if batch.has_fatal_gate_a_issues:
                    fail_gene_list_resolution(batch, output_config)
            global_config, gene_catalog_build = _resolve_build(configs, gene_catalog=catalog)
            print_global_config_banner("run_ldscore_from_args", global_config)
            with PhaseProgress(LOGGER, "validation", "reference content and compatibility"):
                reference_inputs = inspect_reference_inputs(
                    panel_config, global_config, require_autosomes=configs.has_queries,
                    plink_members=declared_inputs.files.get("plink_members"))
                gate = InputGate("Reference/control content validation", issues_path,
                                 reference_inputs.issues.to_dict("records"))
                regression_path = _regr_snps_file_from_config(ldscore_config)
                regression_snps = gate.check("regression SNPs", regression_path,
                    lambda: _load_regression_snps(regression_path, global_config,
                        label="packaged HM3 regression SNP map" if ldscore_config.regr_snps_file is None
                        else "regression SNP list"))
                regression_regions = gate.check("regression regions", configs.regr_snps_exclude_regions,
                    lambda: _regression_region_intervals(
                        configs.regr_snps_exclude_regions, _regression_region_build(global_config, gene_catalog_build)))
                gate.finish()
            chromosome_prefixes = reference_inputs.scope.get("reference_prefixes_by_chrom")
            scope, skipped_funnels = {}, []
            if annotation_config.baseline_annot_sources:
                annotation_bundle, scope = prepare_direct_annotations(
                    annotation_config, global_config, workspace, output_config, batch=batch,
                    gene_catalog_build=gene_catalog_build, declared_inputs=declared_inputs,
                    reference_inputs=reference_inputs)
                ref_panel = _ref_panel_from_config(panel_config, global_config, chromosome_prefixes=chromosome_prefixes)
            else:
                ref_panel = _ref_panel_from_config(panel_config, global_config, chromosome_prefixes=chromosome_prefixes)
                annotation_bundle, skipped_funnels = prepare_synthetic_base(ref_panel, global_config, workspace)
            if scope.get("chromosomes"):
                LOGGER.info("Chromosomes resolved and entering the analysis: %s.", ", ".join(scope['chromosomes']))
            log_inputs(output_dir=str(output_dir), reference_mode=ref_panel.spec.backend,
                reference_panel_source=ref_panel.spec.plink_prefix or ref_panel.spec.r2_matrix_dir,
                snp_identifier=global_config.snp_identifier, genome_build=global_config.genome_build,
                ref_panel_snps_file=panel_config.ref_panel_snps_file or "none",
                regr_snps_file=ldscore_config.regr_snps_file or "packaged_hm3_default",
                regression_region_presets=", ".join(regression_regions.source_labels) if regression_regions else "none")
            if gene_catalog_build is not None:
                LOGGER.info("Gene-list catalog projection build: %s (selected from --genome-build and baseline/reference-panel evidence).", gene_catalog_build)
            if annotation_bundle.gene_list_batch is not None:
                annotation_bundle, gate_b_error, gate_b_issues = _apply_direct_gene_gate_b(
                    annotation_bundle, ref_panel, ldscore_config)
                if gate_b_error:
                    fail_unusable_queries(_bundle_diagnostics(annotation_bundle, gate_b_issues), gate_b_error, output_config)
            statuses = annotation_bundle.query_statuses
            if statuses and not annotation_bundle.query_columns:
                fail_unusable_queries(_bundle_diagnostics(annotation_bundle),
                                      all_query_annotations_skipped_message(statuses), output_config)
            context = _RunContext(global_config, ref_panel, log_path=log_path, build_resolved=True,
                                  reference_checked=True)
            result = LDScoreCalculator()._run(annotation_bundle, ldscore_config, output_config=output_config,
                regression_snps=regression_snps, regression_regions=regression_regions, context=context,
                skipped_funnels=skipped_funnels)
            log_query_outcomes(result.gene_list_batch, result.query_statuses)
            log_outputs(**result.output_paths)
            if configs.emit_console_notices:
                emit_gene_gate_b_notice(result.query_statuses, result.gene_list_batch, command="ldscore")
                emit_resolved_only_notice(result.gene_list_batch)
                emit_bed_query_notice(result.query_statuses, command="ldscore")
                _emit_reference_snp_notices(result.output_paths)
            return dataclass_replace(result, gene_list_batch=persistent_gene_diagnostics(result.gene_list_batch,
                                                                                         result.output_paths))


def _resolve_build(configs: _DirectConfigs, *, gene_catalog=None) -> tuple[GlobalConfig, str | None]:
    """Resolve the run's genome build once against every coordinate input.

    Returns the run's ``GlobalConfig`` and the gene-list projection build
    (``None`` without gene lists). A coordinate-family run carries the
    resolved build. rsID-family identity records none: a gene-list run
    resolves its projection build (``auto`` when ``--genome-build`` is
    omitted), and another run checks a concrete ``--genome-build`` against
    the input declarations.
    """
    global_config = configs.global_config
    catalog_build = getattr(gene_catalog, "genome_build", None)
    if identity_mode_family(global_config.snp_identifier) == "chr_pos":
        build = _resolve_ldscore_chr_pos_genome_build(configs, global_config.genome_build,
                                                      gene_catalog_build=catalog_build)
        return dataclass_replace(global_config, genome_build=build), build if configs.has_gene_lists else None
    if configs.has_gene_lists:
        return global_config, _resolve_ldscore_chr_pos_genome_build(
            configs, configs.genome_build or "auto", gene_catalog_build=catalog_build)
    if configs.genome_build in {"hg19", "hg38"}:
        _resolve_ldscore_chr_pos_genome_build(configs, configs.genome_build)
    return global_config, None


# Options of a direct run; an indexed run takes their values from its gene LD-score index.
_INDEX_OWNED_OPTIONS = frozenset({
    "--baseline-annot-sources",
    "--query-annot-sources",
    "--query-annot-bed-sources",
    "--padding-bp",
    "--gene-coordinate-file",
    "--gene-exclude-regions",
    "--plink-prefix",
    "--r2-matrix-dir",
    "--snp-identifier",
    "--genome-build",
    "--ref-panel-snps-file",
    "--regr-snps-file",
    "--regr-snps-exclude-regions",
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
    "--yes-really",
})


def _indexed_configs(args: argparse.Namespace, explicit_options, emit_console_notices: bool) -> _IndexedConfigs:
    """Check the closed indexed gene-list mode: no direct-run option, gene lists, valid execution settings."""
    supplied = sorted(explicit_options & _INDEX_OWNED_OPTIONS)
    if supplied:
        raise LDSCInputError(
            "ldscore indexed mode inherits immutable scientific inputs, SNP identity, and genome build "
            "from the index and cannot accept live overrides: "
            + ", ".join(supplied)
            + ". Remove these options, or remove --gene-ldscore-index-dir to run direct mode."
        )
    gene_lists = split_cli_path_tokens(args.query_annot_gene_list_sources)
    if not gene_lists:
        raise LDSCInputError(
            "ldscore indexed mode requires --query-annot-gene-list-sources."
        )
    return _IndexedConfigs(
        index_dir=args.gene_ldscore_index_dir,
        query_gene_list_sources=tuple(gene_lists),
        control_gene_list_file=args.control_gene_list_file,
        gene_list_resolution_policy=args.gene_list_resolution_policy,
        query_batch_size=validate_query_batch_size(args.query_batch_size),
        workers=validate_workers(args.workers, error_type=LDSCConfigError),
        output_dir=args.output_dir,
        overwrite=args.overwrite,
        log_level=args.log_level,
        emit_console_notices=emit_console_notices,
    )


def _run_explicit_indexed_ldscore(configs: _IndexedConfigs) -> "LDScoreSource":
    """Run one indexed ``ldscore`` request under its own workflow log."""
    output_dir = output_directory_path(configs.output_dir)
    log_path = output_dir / LDSCORE_LOG
    LDScoreDirectoryWriter.artifact_family(output_dir, include_logs=True).preflight(
        overwrite=configs.overwrite, additional_paths=[log_path],
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with workflow_logging("ldscore", log_path, log_level=configs.log_level):
        log_inputs(
            output_dir=str(output_dir),
            reference_mode="gene_ldscore_index",
            gene_ldscore_index=str(configs.index_dir),
            query_gene_list_count=len(configs.query_gene_list_sources),
            control_gene_list_file=configs.control_gene_list_file,
        )
        result = run_indexed_ldscore(
            configs.index_dir,
            query_gene_list_sources=configs.query_gene_list_sources,
            control_gene_list_file=configs.control_gene_list_file,
            gene_list_resolution_policy=configs.gene_list_resolution_policy,
            query_batch_size=configs.query_batch_size,
            workers=configs.workers,
            output_dir=str(output_dir),
            overwrite=configs.overwrite,
            main_log=log_path,
        )
        if configs.emit_console_notices:
            emit_gene_gate_b_notice(
                getattr(result, "query_statuses", ()),
                getattr(result, "gene_list_batch", None),
                command="ldscore",
            )
            emit_resolved_only_notice(getattr(result, "gene_list_batch", None))
            _emit_reference_snp_notices(getattr(result, "output_paths", None) or {})
        output_paths = getattr(result, "output_paths", None)
        if output_paths:
            log_outputs(**output_paths)
        return result


def _validate_padding_usage(args: argparse.Namespace) -> None:
    """Reject explicit padding outside live BED or gene-list projection.

    The parser's ``None`` sentinel distinguishes omission from an explicit
    zero; the direct configuration turns an allowed omission into 0.
    """
    has_gene_lists = _has_cli_tokens(args.query_annot_gene_list_sources)
    uses_index = args.gene_ldscore_index_dir is not None
    if args.padding_bp is None:
        if has_gene_lists and not uses_index:
            raise LDSCUsageError(
                "Live gene-list mode requires an explicit --padding-bp value; it must be chosen deliberately. "
                "Use --padding-bp 0 for gene bodies or a positive value for padded intervals."
            )
        return
    if (has_gene_lists or _has_cli_tokens(args.query_annot_bed_sources)) and not uses_index:
        return
    raise LDSCUsageError(
        "ldscore accepts `--padding-bp` only with live `--query-annot-bed-sources` "
        "or `--query-annot-gene-list-sources`. It cannot be used with prebuilt "
        "annotation queries, without a query, or in indexed mode. Remove "
        "`--padding-bp` from the command (or remove the Python `padding_bp` keyword)."
    )


# Options of gene-list mode; --gene-list-resolution-policy is the keyword route's --allow-unresolved-genes.
_GENE_LIST_OPTIONS = frozenset({
    "--gene-coordinate-file",
    "--gene-exclude-regions",
    "--control-gene-list-file",
    "--allow-unresolved-genes",
    "--gene-list-resolution-policy",
})


def _validate_gene_list_mode_args(args: argparse.Namespace, explicit_options) -> None:
    """Reject gene-only options outside their direct or indexed input mode."""
    if not _has_cli_tokens(args.query_annot_gene_list_sources):
        if explicit_options & _GENE_LIST_OPTIONS:
            raise LDSCUsageError(
                "Gene coordinate, exclusion, control-list, and resolution policy options are valid only with "
                "gene-list query annotations. Remove the gene-specific option or use gene-list mode."
            )
        return
    if args.gene_ldscore_index_dir is None and args.gene_coordinate_file is None:
        raise LDSCUsageError(
            "Live gene-list mode requires --gene-coordinate-file. Supply a one-build, one-based "
            "gene-coordinate TSV/TSV.GZ; there is no packaged catalog fallback."
        )


def _apply_direct_gene_gate_b(annotation_bundle, ref_panel, ldscore_config):
    """Measure interval support on the prepared reference SNPs.

    Return the updated bundle, the Gate B error message or ``None``, and the
    chromosomes whose reference preparation failed as an input-issue table or
    ``None``.
    """
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
    grid_emptied = []
    for chrom in _chromosomes_from_bundle(annotation_bundle):
        genes = interval_frame[interval_frame["chrom"].astype(str).eq(str(chrom))]
        chrom_bundle = _slice_annotation_bundle(annotation_bundle, chrom)
        chrom_annotations = _chromosome_annotations(chrom_bundle,chrom)
        try:
            with ref_panel.prepare_chromosome(chrom, chrom_annotations, ldscore_config) as prepared:
                positions0 = np.sort(prepared.metadata["POS"].to_numpy(dtype=np.int64) - 1)
        except EmptyReferenceSNPs as exc:
            funnel = getattr(exc, "reference_snp_funnel", None)
            if funnel is not None and annotation_snp_grid_emptied(funnel):
                grid_emptied.append(funnel)
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
    if grid_emptied:
        errors.append(_grid_emptied_message(grid_emptied))
    updated_batch = batch.with_snp_support(support)
    statuses = gene_query_statuses(updated_batch)
    retained_queries = [status.query for status in statuses if status.status in {"ok", "warning"}]
    errors.extend(gene_viability_errors(updated_batch, statuses))
    changes = dict(query_columns=retained_queries, query_statuses=statuses, gene_list_batch=updated_batch)
    return (dataclass_replace(annotation_bundle, **changes), "; ".join(errors) if errors else None,
            pd.DataFrame(input_issues) if input_issues else None)


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


def run_ldscore(**kwargs) -> "LDScoreSource":
    """Run LD-score calculation from Python using public CLI-style names.

    Keyword arguments are interpreted as CLI-equivalent option names without
    leading ``--``; for example ``baseline_annot_sources``, ``query_annot_sources``,
    ``query_annot_bed_sources``, ``query_annot_gene_list_sources``, ``plink_prefix``, ``r2_matrix_dir``,
    ``keep_indivs_file``, ``snp_batch_size``, ``common_maf_min``, ``workers``, and
    ``output_dir``. Shared runtime assumptions such as ``snp_identifier`` and
    ``genome_build`` must be supplied through ``set_global_config(...)`` first,
    while per-run controls such as ``ref_panel_snps_file`` and
    ``regr_snps_file`` remain ordinary keyword arguments here. ``output_dir``
    is required because this public workflow always writes canonical artifacts.
    ``workers`` defaults to one worker; each worker uses 1 BLAS thread
    in serial and parallel execution. Python scripts must launch workflow
    calls under ``if __name__ == "__main__":``.

    When ``baseline_annot_sources`` and query inputs are omitted, the workflow
    builds a synthetic all-ones baseline column named ``base`` from retained
    reference-panel metadata. ``query_annot_sources`` and
    ``query_annot_bed_sources`` and ``query_annot_gene_list_sources`` require
    explicit baseline annotations because query columns are interpreted
    relative to the baseline annotation rows. Gene-list inputs are projected
    from the required one-based gene coordinate catalog for the resolved build.
    ``padding_bp`` defaults to zero for live BED input but must be chosen
    explicitly for live gene lists; passing it with a prebuilt annotation
    query, no query, or a gene LD-score index raises ``LDSCUsageError`` even when its
    value is zero.

    With ``gene_ldscore_index_dir``, the gene LD-score index supplies the SNP
    identifier mode and genome build, as on the command line; the registered
    ``GlobalConfig.snp_identifier`` and ``genome_build`` are ignored, and the
    returned ``config_snapshot`` records the index values.

    Gene-list resolution defaults to ``gene_list_resolution_policy="strict"``.
    For deliberate subset analysis, pass
    ``gene_list_resolution_policy="resolved-only"`` to omit allowlisted rejected
    rows with diagnostics. Structural input, chromosome-coverage, and requested
    control validation still apply. This Python keyword retains its policy
    strings; the CLI uses the value-free ``--allow-unresolved-genes`` flag to
    select ``"resolved-only"`` and rejects ``--gene-list-resolution-policy``.

    Returns
    -------
    LDScoreSource
        Saved artifacts with shared baseline, counts, overlap, diagnostics, and
        canonical paths. Use ``read_queries(names)`` to load selected queries
        across batch files. No query LD tables remain resident after return.
        For a small calculation without writes, use prepared inputs with
        ``LDScoreCalculator.run(..., output_config=None)`` instead.

    Raises
    ------
    LDSCUsageError
        Before anything is written, when a keyword is not an ``ldscore``
        option name, ``output_dir`` is missing, a value is outside the
        command's choices, a flag is not a bool, or ``snp_identifier``,
        ``genome_build`` or ``log_level`` is passed. Also before anything is
        written, when the keywords combine input modes the command rejects;
        a passed keyword counts as supplied even at its default value, so
        ``gene_coordinate_file``, ``gene_exclude_regions``,
        ``control_gene_list_file`` or ``gene_list_resolution_policy``
        without gene lists is rejected.
    LDSCConfigError
        Before anything is written, when a value is outside its valid range.
    """
    from .._runtime.workflow_kwargs import GLOBAL_CONFIG_RUNTIME_OPTIONS, namespace_from_kwargs

    # Indexed mode ignores the registered identity: its live-override check reads explicit options only.
    args = namespace_from_kwargs(
        "run_ldscore", build_parser(), kwargs, rejected=GLOBAL_CONFIG_RUNTIME_OPTIONS, registered=get_global_config(),
    )
    args._emit_query_console_notices = False
    # The outer marker scope makes a failed overwrite name this function, not run_ldscore_from_args.
    with overwrite_failure_marker(args.output_dir, overwrite=args.overwrite, command="Python API: run_ldscore(...)"):
        return run_ldscore_from_args(args)


def main(argv: Sequence[str] | None = None) -> LDScoreSource:
    """Command-line entry point for the LD-score workflow."""
    from .._runtime.workflow_kwargs import parse_command_line

    args = parse_command_line(build_parser(), list(sys.argv[1:] if argv is None else argv))
    args._emit_query_console_notices = True
    return run_ldscore_from_args(args)


def _resolve_ldscore_chr_pos_genome_build(
    configs: _DirectConfigs,
    genome_build: str,
    *,
    gene_catalog_build: str | None = None,
) -> str:
    """Collect every coordinate input, then resolve or report all build evidence.

    A concrete build checks declarations without fresh inference. Auto trusts declarations, otherwise
    samples each file/chromosome; generic coordinate restrictions must resolve
    independently and BED inherits the result.
    """
    annotations, panel = configs.annotation_config, configs.ref_panel_config
    check = GenomeBuildCheck(genome_build, context="ldscore", logger=LOGGER)
    if gene_catalog_build is not None:
        check.declared("--gene-coordinate-file", annotations.gene_coordinate_file or "gene coordinate catalog",
                       gene_catalog_build)
    check.paths("--baseline-annot-sources", annotations.baseline_annot_sources)
    check.paths("--query-annot-sources", annotations.query_annot_sources)
    check.plink("--plink-prefix", panel.plink_prefix)
    if panel.r2_matrix_dir is not None:
        check.r2("--r2-matrix-dir", panel.r2_matrix_dir)
    mode = configs.global_config.snp_identifier
    check.paths("--ref-panel-snps-file", panel.ref_panel_snps_file, restriction_mode=mode)
    check.paths("--regr-snps-file", configs.ldscore_config.regr_snps_file, restriction_mode=mode)
    check.unchecked("--query-annot-bed-sources", annotations.query_annot_bed_sources,
                    "ordinary BED; inherits workflow build")
    resolved = check.resolve()
    if panel.plink_prefix:
        map_sources = getattr(panel, f"genetic_map_{resolved}_sources")
        if map_sources:
            maps = GenomeBuildCheck(resolved, context="ldscore genetic map", logger=LOGGER)
            maps.paths(f"--genetic-map-{resolved}-sources", map_sources)
            maps.resolve()
    return resolved


@report_phase(LOGGER, 'validation', 'regression SNP restriction')
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


def _regression_region_build(global_config: GlobalConfig, gene_catalog_build: str | None) -> str | None:
    """Return the genome build that selects region-exclusion preset BEDs, or ``None`` without one.

    Coordinate-family runs use their resolved build. rsID-family identity
    carries no build, so gene-list runs use their projection build.
    """
    for build in (global_config.genome_build, gene_catalog_build):
        if build in {"hg19", "hg38"}:
            return build
    return None


def _regression_region_intervals(
    choice: str, build: str | None
) -> kernel_regions.RegionIntervals | None:
    """Load the named intervals of ``choice`` in ``build``, applied only after selecting regression SNPs."""
    presets = regr_snps_exclude_regions_choice_to_presets(choice)
    if not presets:
        return None
    if build is None:
        raise LDSCUsageError(
            "ldscore cannot select named regression exclusion regions without a concrete genome build. "
            "Pass `--genome-build hg19` or `--genome-build hg38`, or use `--regr-snps-exclude-regions none`."
        )
    return kernel_regions.load_preset_intervals(presets, build)


def _ref_panel_from_config(config: RefPanelConfig, global_config: GlobalConfig, *, chromosome_prefixes=None):
    """Build the reference-panel adapter that owns the ``A -> A'`` restriction.

    ``global_config`` carries the build the run already resolved against the
    panel, so unlike ``RefPanelLoader.load`` this checks no build again.
    """
    if config.backend == "plink":
        return PlinkRefPanel(global_config, config, chromosome_prefixes=chromosome_prefixes)
    return R2MatrixRefPanel(global_config, config)


def _reference_metadata_sources(ref_panel) -> tuple[str, str]:
    """Return ``(cm_source, maf_source)`` provenance labels for the reference panel."""
    spec = getattr(ref_panel, "spec", None)
    backend = getattr(spec, "backend", None)
    if backend == "r2_matrix":
        return "parquet_sidecar", "parquet_sidecar"
    has_map = bool(getattr(spec, "genetic_map_hg19_sources", None) or getattr(spec, "genetic_map_hg38_sources", None))
    return ("genetic_map" if has_map else "plink_bim"), "plink_genotypes"


def _snp_universe_policy(
    *,
    ref_panel,
    ldscore_config: LDScoreConfig,
    regression_regions: kernel_regions.RegionIntervals | None,
    chromosome_results: Sequence[ChromLDScoreResult | DirectChromosomeFragment],
) -> dict[str, Any]:
    """Describe the computed chromosomes' reference SNPs and regression SNPs.

    See ``snp_universe_policy_record`` for the persisted fields.
    """
    return snp_universe_policy_record(
        ref_panel_snps_file=getattr(getattr(ref_panel, "spec", None), "ref_panel_snps_file", None),
        regr_snps_file=ldscore_config.regr_snps_file,
        regression_regions=regression_regions,
        retained_snp_count=sum(result.reference_snp_count for result in chromosome_results),
        selected_snp_count_before_region_exclusion=sum(
            result.regression_selected_snp_count for result in chromosome_results),
        region_excluded_snp_count=sum(result.regression_region_removed_snp_count for result in chromosome_results),
        written_snp_count=sum(len(result.baseline_table) for result in chromosome_results),
    )


def _ld_window_snp_counts(ref_panel, ldscore_config: LDScoreConfig, chromosome_results) -> pd.DataFrame | None:
    """Per-chromosome R2 matrix and kept SNP counts behind an R2 matrix ``--ld-wind-snps`` window.

    On an R2 matrix the window counts R2 matrix (sidecar) rows, so a chromosome
    keeping fewer SNPs than the R2 matrix gives its SNPs fewer than K reference SNPs
    per side. ``None`` for kb/cM windows and PLINK reference panels.
    """
    if ldscore_config.ld_wind_snps is None or getattr(getattr(ref_panel, "spec", None), "backend", None) != "r2_matrix":
        return None
    rows = [
        {
            "CHR": str(result.chrom),
            "panel_snps": ref_panel.open_chromosome(result.chrom).contract.n_snps,
            "kept_snps": int(result.reference_snp_count),
            "ld_wind_snps": int(ldscore_config.ld_wind_snps),
            "counting_basis": "panel_rows",
        }
        for result in chromosome_results
    ]
    return pd.DataFrame(rows, columns=["CHR", "panel_snps", "kept_snps", "ld_wind_snps", "counting_basis"])


def _resolve_calculator_build(annotation_bundle, ref_panel, global_config: GlobalConfig):
    """Resolve the build of a prepared-input calculation against its annotations and reference panel.

    ``LDScoreCalculator.run`` callers supply inputs the workflow has not
    checked. Returns the ``GlobalConfig`` and the reference panel, both with
    a concrete build in coordinate-family modes.
    """
    if identity_mode_family(global_config.snp_identifier) == "chr_pos":
        check = GenomeBuildCheck(global_config.genome_build, context='LDScoreCalculator', logger=LOGGER)
        for role, obj in (("annotations", annotation_bundle), ("reference panel", ref_panel)):
            snapshot = getattr(obj, 'config_snapshot', None) or getattr(obj, 'global_config', None)
            build = getattr(snapshot, 'genome_build', None)
            if build in {'hg19', 'hg38'}:
                check.declared(role, role, build)
            if role == 'annotations':
                for chrom in obj.chromosomes:
                    check.frames(role, f'annotations [chr{chrom}]', [obj.metadata_for_chromosome(chrom)],
                                 infer_unlabeled=build not in {'hg19', 'hg38'})
            if role == 'reference panel':
                if obj.spec.backend == 'plink':
                    check.plink('--plink-prefix', obj.spec.plink_prefix)
                else:
                    check.r2('--r2-matrix-dir', obj.spec.r2_matrix_dir)
                check.paths('--ref-panel-snps-file', obj.spec.ref_panel_snps_file, restriction_mode=global_config.snp_identifier)
        global_config = dataclass_replace(global_config, genome_build=check.resolve())
    if hasattr(ref_panel, 'global_config') and global_config.genome_build in {'hg19', 'hg38'}:
        ref_panel = copy(ref_panel)
        ref_panel.global_config = dataclass_replace(ref_panel.global_config, genome_build=global_config.genome_build)
    return global_config, ref_panel


def _check_parquet_reference(ref_panel, chromosomes, ldscore_config: LDScoreConfig, global_config: GlobalConfig) -> None:
    """Check every chromosome's R2 matrix footer and sidecar before any is computed.

    The CLI runs the same checks in its input preflight. Python API callers would
    otherwise learn of a defect only from the first chromosome worker to reach it.
    """
    if getattr(getattr(ref_panel, "spec", None), "backend", None) != "r2_matrix":
        return
    from .._refpanel.preflight import inspect_r2_headers

    def check_sidecar(chrom):
        chromosome = ref_panel.open_chromosome(chrom)
        chromosome.pair_reader(chromosome.sidecar(), global_config.snp_identifier, global_config.genome_build).close()

    r2_matrix_dir = ref_panel.spec.r2_matrix_dir
    gate = InputGate("R2 matrix validation")
    gate.check("reference", r2_matrix_dir, lambda: inspect_r2_headers(
        r2_matrix_dir, global_config.genome_build, command="ldscore", ld_window=ldscore_config.ld_window))
    for chrom in chromosomes:
        gate.check("reference", r2_matrix_dir, lambda chrom=chrom: check_sidecar(chrom), chrom=chrom)
    gate.finish()


def _warn_ld_window_snp_subset(counts: pd.DataFrame | None, where: str) -> None:
    """Warn once per run when a parquet ``--ld-wind-snps`` window spans dropped R2 matrix SNPs."""
    if counts is None:
        return
    subset = counts.loc[counts["kept_snps"] < counts["panel_snps"]]
    if subset.empty:
        return
    first = subset.iloc[0]
    window = int(first["ld_wind_snps"])
    LOGGER.warning(
        f"ldscore `--ld-wind-snps {window}` counts positions in the R2 matrix; this run kept fewer "
        f"SNPs than the R2 matrix on {len(subset)} chromosome(s) (e.g. chr{first['CHR']}: {first['kept_snps']} of "
        f"{first['panel_snps']}), so those SNPs have fewer than {window} reference SNPs per side. "
        f"Per-chromosome counts: {where}."
    )


def _emit_reference_snp_notices(output_paths: dict) -> None:
    """Print the funnel and dropped-SNP WARNINGs to stderr after a successful command-line run."""
    path = output_paths.get("reference_snp_funnel")
    if path is not None:
        emit_reference_snp_funnel_notice(read_reference_snp_funnel(path), path)
    emit_dropped_snps_notice([path for key, path in output_paths.items() if key.startswith("dropped_snps")],
                             command="ldscore")


def _grid_emptied_message(funnels: Sequence[pd.DataFrame]) -> str:
    """Describe every chromosome whose annotation SNP grid matched no remaining panel SNP."""
    details = []
    for funnel in funnels:
        grid = funnel.loc[funnel["step"].eq("annotation_snp_grid")].set_index("funnel")
        chrom = funnel["chrom"].iloc[0]
        panel_snps = int(grid.loc[REFERENCE_PANEL, "removed"])
        grid_rows = int(grid.loc["annotation_snp_grid", "removed"])
        details.append(f"chromosome {chrom} ({grid_rows} annotation rows, {panel_snps} reference-panel SNPs)")
    return (
        "ldscore found no reference SNPs on the annotation SNP grid for " + "; ".join(details) + ": "
        "after the panel's filters, none of the panel SNPs on these chromosomes has an annotation row. "
        "Most likely the annotation SNP identifiers, genome build, or allele-aware identifier mode do not "
        "match the reference panel. Use annotation and reference-panel artifacts built with the same SNP "
        "identifier mode and genome build, or leave a chromosome out of the annotation files to exclude it on purpose. "
        f"Other causes & fixes: {_LDSCORE_INTERSECTION_DOC}"
    )


def _annotation_result_message(caller: str) -> str:
    """Explain how to read an annotation result where a bundle is required."""
    return (f"{caller} cannot read an annotation result (AnnotationResult) from run_annotate directly: "
            "it reads no files until you prepare it. Call result.prepared_in(directory), which prepares the "
            "exported annotation files in a private folder of a directory you choose, and pass the returned "
            "AnnotationBundle; "
            "close that bundle when you are done.")


def _bundle_diagnostics(bundle, input_issues=None) -> QueryDiagnostics:
    """Return the query diagnostics an annotation bundle holds, plus any Gate B input issues, for the writer."""
    return QueryDiagnostics(dict(bundle.source_summary.get("chromosome_scope", {})), tuple(bundle.query_statuses),
                            bundle.gene_list_batch, input_issues)


def _write_one_ref_metadata_sidecar(metadata: pd.DataFrame, chrom: str, output_dir: str) -> None:
    """Write one opt-in reference-metadata file (PLINK backend).

    The file is a gzip TSV with the ``CHR POS SNP A1 A2 CM MAF`` columns,
    the minor allele in the retained reference individuals as ``A1``, and
    the float format of an R2 matrix sidecar. It is not an R2 matrix
    sidecar: it has no ``# ldsc:`` metadata lines, and ``CM`` is the value
    this run used, interpolated from a supplied genetic map for a cM window
    and otherwise the ``.bim`` value (often 0), where ``build-r2-matrix``
    writes ``NA`` without a genetic map. Written from
    inside the per-chromosome worker so the lean cross-process result does
    not carry per-SNP metadata.
    """
    columns = [c for c in ("CHR", "POS", "SNP", "A1", "A2", "CM", "MAF") if c in metadata.columns]
    path = Path(output_dir) / ref_metadata_file(chrom)
    path.parent.mkdir(parents=True, exist_ok=True)
    metadata[columns].to_csv(path, sep="\t", index=False, na_rep="NA", float_format=SIDECAR_FLOAT_FORMAT,
                             compression="gzip")
    LOGGER.info(f"Wrote reference-metadata sidecar '{ref_metadata_file(chrom)}'.")


def _chromosomes_from_bundle(annotation_bundle) -> list[str]:
    """Return declared chromosome order without reading annotation metadata."""
    return list(annotation_bundle.chromosomes)


def _slice_annotation_bundle(annotation_bundle, chrom: str):
    """Borrow one chromosome descriptor without copying annotation payloads."""
    drops = (annotation_bundle.identity_drops if annotation_bundle.chromosome_identity_drops is None
             else annotation_bundle.chromosome_identity_drops.get(str(chrom)))
    return dataclass_replace(annotation_bundle, shards={str(chrom): annotation_bundle.shard(chrom)},
                             gene_list_batch=None, identity_drops=drops, chromosome_identity_drops=None)


def _partition_annotation_drops(bundle):
    """Partition source diagnostics once, including fully removed chromosomes."""
    from uuid import uuid4

    artifacts = {}
    if bundle.identity_drops is None:
        return artifacts
    if bundle.diagnostics_in_memory:
        from .._annotation.in_memory import MemoryDiagnostics
        frames = list(bundle.identity_drops.frames())
        if frames:
            for chrom, rows in pd.concat(frames, ignore_index=True).groupby("CHR", sort=False):
                artifacts[str(chrom)] = MemoryDiagnostics(rows.reset_index(drop=True))
        return artifacts
    prefix = uuid4().hex
    for frame in bundle.identity_drops.frames():
        for chrom, rows in frame.groupby("CHR", sort=False):
            chrom = str(chrom)
            first = chrom not in artifacts
            if first:
                artifacts[chrom] = TsvDiagnostics(bundle.workspace.path / f"source-drops-{prefix}-{chrom}.tsv.gz", bundle.workspace)
            rows.to_csv(artifacts[chrom].path, sep="\t", index=False, header=first, mode="a", compression="gzip", na_rep="")
    return artifacts


def _stage_chromosome_drops(bundle, chrom, reference_drops):
    """Stream source cleanup followed by reference-stage drops into one artifact."""
    import gzip
    from uuid import uuid4

    if bundle.diagnostics_in_memory:
        from .._annotation.in_memory import MemoryDiagnostics
        frames = [frame.loc[frame.CHR.astype(str).eq(str(chrom))] for frame in bundle.identity_drops.frames()] if bundle.identity_drops is not None else []
        frames.append(coerce_identity_drop_frame(reference_drops))
        return MemoryDiagnostics(coerce_identity_drop_frame(pd.concat(frames, ignore_index=True)))

    path = bundle.workspace.path / f"ldscore-drops-{chrom}-{uuid4().hex}.tsv.gz"
    with gzip.open(path, "wt") as stream:
        pd.DataFrame(columns=IDENTITY_DROP_COLUMNS).to_csv(stream, sep="\t", index=False)
        if bundle.identity_drops is not None:
            for frame in bundle.identity_drops.frames():
                frame = frame.loc[frame.CHR.astype(str).eq(str(chrom))]
                if not frame.empty:
                    frame.to_csv(stream, sep="\t", index=False, header=False, na_rep="", chunksize=65536)
        coerce_identity_drop_frame(reference_drops).to_csv(stream, sep="\t", index=False, header=False,
                                                         na_rep="", chunksize=65536)
    return TsvDiagnostics(path, bundle.workspace)


def _chromosome_annotations(bundle, chrom):
    """Borrow one chromosome's metadata and selected annotation reader."""
    shard = bundle.shard(chrom)
    metadata = shard.metadata()
    values = MappedAnnotations(shard, np.arange(shard.n_rows), tuple(bundle.baseline_columns + bundle.query_columns))
    return ChromosomeAnnotations(metadata, values, list(bundle.baseline_columns), list(bundle.query_columns))


# Worker state shared within a pool worker process. ``regression_snps`` is set
# once by the pool initializer so a large SNP set is transferred per worker, not
# per task; direct task probes can pass it explicitly instead.
_WORKER_UNSET = "__use_worker_global__"
_WORKER_STATE: dict[str, Any] = {}


@dataclass(frozen=True)
class _ChromOutcome:
    """Tagged result of one chromosome's LD-score computation in a worker."""

    chrom: str
    result: "ChromLDScoreResult | DirectChromosomeFragment | None"
    skipped: bool
    skip_message: str | None = None
    reference_snp_funnel: pd.DataFrame | None = None


def _compute_one_chromosome(
    chrom: str,
    chrom_bundle,
    ref_panel_spec,
    ldscore_config: LDScoreConfig,
    global_config: GlobalConfig,
    export_dir: str | None = None,
    regression_snps=_WORKER_UNSET,
    regression_regions=_WORKER_UNSET,
    chromosome_prefixes=None,
    log_path=None, batch_label="batch 1/1", fragment_dir=None, include_baseline=True,
) -> _ChromOutcome:
    """Compute one chromosome end-to-end and return a tagged outcome.

    Rebuilds the reference-panel adapter locally from ``ref_panel_spec`` so no
    live reader/file handle crosses the process boundary; ``global_config``
    carries the build the parent resolved, which is not checked again.
    ``regression_snps``
    defaults to the pool-initializer-provided global; direct task callers may
    pass it explicitly. Only ``EmptyReferenceSNPs`` returns a skipped outcome,
    carrying the error's reference-SNP funnel so the parent can tell a
    chromosome emptied by the annotation SNP grid from one emptied by an
    explicit restriction; other failures propagate across the pool boundary.
    With ``fragment_dir``, the outcome carries a private query path and compact
    summaries; otherwise it retains the materialized chromosome result.
    ``include_baseline=False`` (later batches of a writing run) computes only
    the query scores and statistics and needs ``fragment_dir``.
    """
    with chromosome_log(log_path, chrom, batch_label):
        if regression_snps == _WORKER_UNSET:
            regression_snps = _WORKER_STATE.get("regression_snps")
        if regression_regions == _WORKER_UNSET:
            regression_regions = _WORKER_STATE.get("regression_regions")
        prefixes = chromosome_prefixes
        if prefixes is None:
            prefixes = chrom_bundle.source_summary.get('chromosome_scope',{}).get('reference_prefixes_by_chrom')
        ref_panel = _ref_panel_from_config(ref_panel_spec, global_config, chromosome_prefixes=prefixes or None)
        try:
            if not include_baseline:
                result = _compute_query_chromosome(chrom, chrom_bundle, ref_panel, ldscore_config, global_config,
                                                   regression_snps, regression_regions, fragment_dir)
                return _ChromOutcome(chrom=chrom, result=result, skipped=False)
            result = LDScoreCalculator().compute_chromosome(
                chrom=chrom,
                annotation_bundle=chrom_bundle,
                ref_panel=ref_panel,
                ldscore_config=ldscore_config,
                global_config=global_config,
                regression_snps=regression_snps,
                regression_regions=regression_regions,
                export_dir=export_dir,
            )
        except EmptyReferenceSNPs as exc:
            return _ChromOutcome(chrom=chrom, result=None, skipped=True, skip_message=str(exc),
                                 reference_snp_funnel=getattr(exc, "reference_snp_funnel", None))
        if fragment_dir is not None:
            result = stage_direct_chromosome(result, fragment_dir)
        return _ChromOutcome(chrom=chrom, result=result, skipped=False)


def _compute_query_chromosome(chrom, chrom_bundle, ref_panel, ldscore_config, global_config,
                              regression_snps, regression_regions, fragment_dir) -> DirectChromosomeFragment:
    """Compute one chromosome of a later query batch and stage its query scores.

    The first batch already returned this chromosome's baseline table, counts,
    overlap and identity drops. This batch projects only its queries; the rows
    are sorted exactly as the first batch sorted its combined table, and the
    parent checks them against the shared first-batch baseline.
    """
    chrom_annotations = _chromosome_annotations(chrom_bundle, chrom)
    LOGGER.info(f"Computing chromosome {chrom} query LD scores with backend '{ref_panel.spec.backend}' "
                f"and {len(chrom_annotations.metadata)} annotation rows.")
    with ref_panel.prepare_chromosome(chrom, chrom_annotations, ldscore_config) as prepared:
        funnel = prepared.reference_snp_funnel
        chrom_result = kernel_ldscore.compute_chromosome(
            chrom, prepared, snp_identifier=global_config.snp_identifier,
            snp_batch_size=ldscore_config.snp_batch_size,
            common_maf_min=ldscore_config.common_maf_min,
            regression_keys=regression_snps, regression_regions=regression_regions,
            query_batch_size=ldscore_config.query_batch_size, include_baseline=False,
        )
    metadata = chrom_result.metadata
    metadata_columns = ["CHR", "SNP", "POS", *[column for column in ("A1", "A2") if column in metadata.columns]]
    query_table = sort_frame_by_genomic_position(pd.concat(
        [metadata.loc[:, metadata_columns],
         pd.DataFrame(chrom_result.ld_scores, columns=list(chrom_result.query_columns))], axis=1))
    path, variation = stage_query_table(chrom_result.chrom, query_table, chrom_result.query_columns, fragment_dir)
    LOGGER.info(f"Finished chromosome {chrom} with {len(query_table)} retained SNP rows.")
    return DirectChromosomeFragment(
        chrom=str(chrom_result.chrom), baseline_table=None, query_fragment=path,
        baseline_columns=list(chrom_result.baseline_columns), query_columns=list(chrom_result.query_columns),
        snp_count_totals={
            "all_reference_snp_counts": np.asarray(chrom_result.M, dtype=np.float64),
            "common_reference_snp_counts": np.asarray(chrom_result.M_5_50, dtype=np.float64),
        },
        config_snapshot=global_config, overlap=chrom_result.overlap,
        reference_snp_count=chrom_result.reference_snp_count,
        regression_selected_snp_count=chrom_result.regression_selected_snp_count,
        regression_region_removed_snp_count=chrom_result.regression_region_removed_snp_count,
        identity_drops=None, annotation_types=dict(chrom_result.annotation_types),
        variation=variation, reference_snp_funnel=funnel,
    )


def _init_worker(regression_snps, regression_regions=None, log_level: str = "INFO") -> None:
    """Initialize shared keys and fixed, worker-owned native thread limits."""
    _WORKER_STATE["regression_snps"] = regression_snps
    _WORKER_STATE["regression_regions"] = regression_regions
    initialize_worker(log_level)
