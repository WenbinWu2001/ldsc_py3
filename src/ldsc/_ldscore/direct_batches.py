"""The direct LD-score producer: one chromosome of the reference panel per task.

``_ldscore.batch_execution`` runs one query batch at a time through the run's
persistent pool. Every task, in the first and in later query batches, runs
``compute_direct_chromosome``: it prepares the chromosome through
``ref_panel.prepare_chromosome``, projects the batch's annotations with the
LD-score kernel once, and returns one record built directly from the kernel's
result. A writing run's task writes its post-kernel float32 query scores to a
private Parquet fragment and returns a ``DirectChromosomeFragment`` with the
fragment path and compact summaries; the no-output path returns a
``ChromLDScoreResult``. Only the first batch computes the baseline table and
the identity drops; later batches compute their queries alone, and the writer
checks their rows against the shared first-batch baseline. ``_DirectProducer``
combines the chromosome records into one staged or materialized batch. The
caller owns the batch's fragment directory through consumption or pool
shutdown.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, replace
import gzip
import logging
from pathlib import Path
from typing import Any, Sequence
from uuid import uuid4

import numpy as np
import pandas as pd

from .._annotation.interval_queries import execution_query_bundle
from .._annotation.overlap import sum_overlap_contributions
from .._annotation.semantics import combine_annotation_types
from .._annotation.storage import TsvDiagnostics
from .._core.chromosomes import chrom_sort_key
from .._core.row_alignment import concat_chromosome_tables, sort_frame_by_genomic_position
from .._core.snp_identity import coerce_identity_drop_frame, IDENTITY_DROP_COLUMNS, is_allele_aware_mode
from .._refpanel.adapters import PlinkRefPanel, resolved_ref_panel
from .._refpanel.funnel import REFERENCE_PANEL, annotation_snp_grid_emptied, concat_funnels, funnel_summary
from .._refpanel.prepared import ChromosomeAnnotations, MappedAnnotations
from .._refpanel.r2_format import SIDECAR_FLOAT_FORMAT
from .._runtime.worker_runtime import initialize_worker
from ..config import GlobalConfig, LDScoreConfig, validate_config_compatibility
from ..errors import EmptyReferenceSNPs, LDSCInputError, LDSCInternalError
from . import compute as kernel_ldscore
from .diagnostics import chromosome_log, log_reference_snp_funnel_warning, reference_snp_funnel_path
from .outputs import LD_WINDOW_SNP_COUNTS_FILE, common_count_config, ref_metadata_file, snp_universe_policy_record
from .source import REGRESSION_LD_SCORE_COLUMN, count_record
from .query_variation import summarize_variation
from .results import ChromLDScoreResult, DirectChromosomeFragment, LDScoreResult, StagedLDScoreBatch


LOGGER = logging.getLogger("LDSC.ldscore_calculator")

_LDSCORE_INTERSECTION_DOC = (
    "docs/user-guide/troubleshooting.md#ldscore-no-reference-snps-remain-on-the-annotation-snp-grid"
)


def compute_direct_chromosome(
    chrom: str,
    annotation_bundle,
    ref_panel,
    ldscore_config: LDScoreConfig,
    global_config: GlobalConfig,
    regression_snps=None,
    regression_regions=None,
    *,
    include_baseline: bool = True,
    fragment_dir: str | None = None,
    export_dir: str | None = None,
) -> ChromLDScoreResult | DirectChromosomeFragment:
    """Compute one chromosome's LD scores from the reference panel.

    ``ref_panel.prepare_chromosome`` aligns annotations to the retained
    reference SNPs ``B_chrom & A'_chrom`` and owns the reader until the kernel
    returns or raises; the kernel projects the annotations once. Returns a
    ``ChromLDScoreResult``, or with ``fragment_dir`` a
    ``DirectChromosomeFragment`` whose query scores are staged in
    ``fragment_dir``. ``include_baseline=False`` (later query batches of a
    writing run) computes only the query scores and statistics, stages no
    identity drops, and needs ``fragment_dir``. With ``export_dir`` the
    chromosome's reference metadata is written there. Raises
    ``EmptyReferenceSNPs`` when preparation leaves no reference SNP, and
    ``LDSCInternalError`` for ``include_baseline=False`` without
    ``fragment_dir``.
    """
    if not include_baseline and fragment_dir is None:
        raise LDSCInternalError(
            "LD-score chromosome task failed in compute_direct_chromosome(): a query-only task "
            "(include_baseline=False) returns no baseline table and needs fragment_dir to stage its "
            "query scores. Most likely a caller other than the direct LD-score producer's later query "
            "batches requested it. Re-run with DEBUG logging and report the traceback."
        )
    backend =getattr(getattr(ref_panel, "spec", None), "backend", None)
    chrom_annotations = _chromosome_annotations(annotation_bundle, chrom)
    LOGGER.info(f"Computing chromosome {chrom} {'LD scores' if include_baseline else 'query LD scores'} "
                f"with backend '{backend or 'unknown'}' and {len(chrom_annotations.metadata)} annotation rows.")
    with ref_panel.prepare_chromosome(chrom, chrom_annotations, ldscore_config) as prepared:
        funnel = prepared.reference_snp_funnel
        # Later batches prepare the chromosome exactly as the first did, which logged its funnel.
        if include_baseline and funnel is not None:
            LOGGER.info("Chromosome %s reference-SNP funnel: %s.", chrom, funnel_summary(funnel))
        computed = kernel_ldscore.compute_chromosome(
            chrom, prepared, snp_identifier=global_config.snp_identifier,
            snp_batch_size=ldscore_config.snp_batch_size,
            common_maf_min=ldscore_config.common_maf_min,
            regression_keys=regression_snps, regression_regions=regression_regions,
            query_batch_size=ldscore_config.query_batch_size, include_baseline=include_baseline,
        )
        if export_dir is not None:
            _write_one_ref_metadata_sidecar(prepared.metadata, chrom, export_dir)
    drops = (_stage_chromosome_drops(annotation_bundle, chrom, computed.identity_drops)
             if include_baseline else None)
    result = _chromosome_record(computed, global_config, include_baseline=include_baseline, identity_drops=drops,
                                reference_snp_funnel=funnel, fragment_dir=fragment_dir)
    LOGGER.info(f"Finished chromosome {chrom} with {len(computed.metadata)} retained SNP rows.")
    return result


def _chromosome_record(computed: kernel_ldscore.ChromComputationResult, global_config: GlobalConfig, *,
                       include_baseline: bool = True, identity_drops: TsvDiagnostics | None = None,
                       reference_snp_funnel: pd.DataFrame | None = None, fragment_dir: str | None = None):
    """Turn one kernel chromosome result into the chromosome's record.

    Returns a ``ChromLDScoreResult`` holding the query table, or with
    ``fragment_dir`` a ``DirectChromosomeFragment`` holding the path of the
    query table staged there and its variation summary. Both records take
    every other field from the same values.
    """
    baseline_table, query_table = _ldscore_tables(computed, global_config.snp_identifier, include_baseline)
    fields = dict(
        chrom=str(computed.chrom),
        baseline_table=baseline_table,
        baseline_columns=list(computed.baseline_columns),
        query_columns=list(computed.query_columns),
        snp_count_totals={
            "all_reference_snp_counts": np.asarray(computed.M, dtype=np.float64),
            "common_reference_snp_counts": np.asarray(computed.M_5_50, dtype=np.float64),
        },
        config_snapshot=global_config,
        overlap=computed.overlap,
        reference_snp_count=computed.reference_snp_count,
        regression_selected_snp_count=computed.regression_selected_snp_count,
        regression_region_removed_snp_count=computed.regression_region_removed_snp_count,
        identity_drops=identity_drops,
        annotation_types=dict(computed.annotation_types),
        reference_snp_funnel=reference_snp_funnel,
    )
    if fragment_dir is None:
        return ChromLDScoreResult(query_table=query_table, **fields)
    # Persist the existing float32 scores without recomputation or further casting.
    variation = summarize_variation(query_table[name].to_numpy(copy=False) for name in computed.query_columns)
    path = None
    if query_table is not None:
        path = Path(fragment_dir) / f"chr{computed.chrom}.parquet"
        query_table.to_parquet(path, index=False, compression="zstd", use_dictionary=False)
    return DirectChromosomeFragment(query_fragment=path, variation=variation, **fields)


def _ldscore_tables(computed, snp_identifier: str, include_baseline: bool = True):
    """Sort one kernel chromosome result genomically and split it into baseline and query tables.

    The kernel's metadata is already a private frame with a RangeIndex. The
    baseline and query tables are column subsets of one sorted table, so their
    rows match by construction and need no row comparison. Without the
    baseline, the rows sort exactly as the first batch's combined table did.
    """
    metadata = computed.metadata
    metadata_columns = ["CHR", "SNP", "POS", *[column for column in ("A1", "A2") if column in metadata.columns]]
    pieces = [metadata.loc[:, metadata_columns], pd.DataFrame(computed.ld_scores, columns=list(computed.ldscore_columns))]
    if include_baseline:
        regression_weights = np.asarray(computed.w_ld, dtype=np.float32).reshape(-1)
        pieces.append(pd.DataFrame({REGRESSION_LD_SCORE_COLUMN: regression_weights}, copy=False))
    ldscore_table = sort_frame_by_genomic_position(pd.concat(pieces, axis=1, copy=False))
    return _split_ldscore_table(
        ldscore_table,
        baseline_columns=list(computed.baseline_columns) if include_baseline else None,
        query_columns=list(computed.query_columns),
        snp_identifier=snp_identifier,
    )


def _split_ldscore_table(
    ldscore_table: pd.DataFrame,
    *,
    baseline_columns: list[str] | None,
    query_columns: list[str],
    snp_identifier: str = "chr_pos_allele_aware",
) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    """Split a merged LD-score table into optional baseline and query tables.

    ``baseline_columns=None`` (later query batches) returns no baseline table.
    """
    if is_allele_aware_mode(snp_identifier) and not {"A1", "A2"}.issubset(ldscore_table.columns):
        raise LDSCInternalError(
            "LD-score split-table assembly failed in _split_ldscore_table(): allele-aware "
            f"SNP identifier mode '{snp_identifier}' requires A1/A2 columns, but the LD-score "
            "table does not contain both. Most likely an upstream reference-panel or annotation "
            "assembly step dropped allele columns. Re-run with DEBUG logging and report the traceback."
        )
    metadata_columns = ["CHR", "SNP", "POS", *[column for column in ("A1", "A2") if column in ldscore_table.columns]]
    query_order = [*metadata_columns, *query_columns]
    if not isinstance(ldscore_table.index, pd.RangeIndex) or ldscore_table.index.start != 0:
        ldscore_table = ldscore_table.reset_index(drop=True)
    # Both tables are column selections (copies) of one frame, so their rows match.
    baseline_table = None
    if baseline_columns is not None:
        baseline_order = [*metadata_columns, REGRESSION_LD_SCORE_COLUMN, *baseline_columns]
        missing_baseline = [column for column in baseline_order if column not in ldscore_table.columns]
        if missing_baseline:
            raise LDSCInternalError(
                "LD-score split-table assembly failed in _split_ldscore_table(): "
                f"the merged table is missing baseline columns {missing_baseline}. Most likely "
                "a kernel result omitted expected baseline annotation scores. Re-run with DEBUG "
                "logging and report the traceback."
            )
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


def _aggregate_chromosome_results(
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
    return [count_record(group, column, all_counts[idx], None if common_counts is None else common_counts[idx])
            for idx, (group, column) in enumerate(zip(groups, columns))]


@dataclass(frozen=True)
class _DirectBatch:
    """One direct query batch: its prepared annotations and private fragment directory."""

    bundle: Any
    source_drops: dict
    ordinal: int = 0
    fragment_dir: str | None = None


@dataclass(frozen=True)
class _SkippedChromosome:
    """A chromosome task whose reference preparation left no SNP (``EmptyReferenceSNPs``).

    ``reference_snp_funnel`` tells a chromosome emptied by the annotation SNP
    grid from one emptied by an explicit restriction or genotype QC.
    """

    message: str
    reference_snp_funnel: pd.DataFrame | None = None


class _DirectProducer:
    """Direct LD-score producer: each task computes one chromosome from the reference panel.

    A batch projects its active query annotations, workers compute and stage
    each chromosome's LD scores, and ``reduce`` combines them with the run's
    direct provenance. The no-output path reuses ``prepare``, the tasks, and
    ``reduce`` for its single materialized batch. ``skipped_funnels`` belong to
    chromosomes skipped before the bundle was built, which have no task.
    """

    batch_name, phase, task = "batch", "LD-score calculation", "LD-score computation"

    def __init__(self, annotation_bundle, ref_panel, ldscore_config, global_config, regression_regions, log_paths,
                 *, output_dir=None, export_dir=None, skipped_funnels=()):
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
        self.chromosomes = list(annotation_bundle.chromosomes)
        # Tasks use the panel rebuilt once with the run's configuration, which
        # carries the resolved build; its PLINK chromosome prefixes are
        # resolved here, and its SNP restriction and keep list at the first
        # submission, rather than in every task.
        self.task_panel = resolved_ref_panel(
            global_config, ref_panel.spec,
            chromosome_prefixes=ref_panel._resolved_prefixes() if isinstance(ref_panel, PlinkRefPanel) else None)

    def prepare(self, bundle, *, ordinal=0, fragment_dir=None) -> _DirectBatch:
        """Check the bundle against the run's identity and partition its source identity drops."""
        global_config = self.global_config
        if bundle.config_snapshot is not None:
            validate_config_compatibility(
                (replace(bundle.config_snapshot, genome_build=global_config.genome_build)
                 if global_config.genome_build in {"hg19", "hg38"} else bundle.config_snapshot),
                global_config,
                context="AnnotationBundle and LDScoreCalculator runtime config",
            )
        # Identity drops are published from the first batch only.
        source_drops = _partition_annotation_drops(bundle) if ordinal == 0 else {}
        return _DirectBatch(replace(bundle, chromosome_identity_drops=source_drops), source_drops,
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
        self.task_panel._resolve_run_inputs()  # read once in this process; each task's copy carries them
        return pool.submit(_compute_one_chromosome, chrom, _slice_annotation_bundle(batch.bundle, chrom),
                           self.task_panel, self.ldscore_config, self.global_config,
                           include_baseline=first or batch.fragment_dir is None, fragment_dir=batch.fragment_dir,
                           export_dir=self.export_dir if first else None,
                           log_path=self.log_paths.get(chrom), batch_label=label)

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
            funnels.append(outcome.reference_snp_funnel)
            if isinstance(outcome, _SkippedChromosome):
                funnel = outcome.reference_snp_funnel
                if funnel is not None and annotation_snp_grid_emptied(funnel):
                    grid_emptied.append(funnel)
                else:
                    LOGGER.warning("Skipping chromosome %s: %s", chrom, outcome.message)
                continue
            chromosome_results.append(outcome)
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
        result = _aggregate_chromosome_results(
            chromosome_results,
            count_config=common_count_config(ldscore_config.common_maf_min),
        )
        result = replace(
            result,
            identity_drops_by_chrom={**batch.source_drops, **result.identity_drops_by_chrom},
            gene_list_batch=getattr(batch.bundle, "gene_list_batch", None),
            chromosome_scope=dict(batch.bundle.source_summary.get("chromosome_scope", {})),
            ld_window={**ldscore_config.ld_window.record(), "whole_chromosome_ok": ldscore_config.whole_chromosome_ok},
        )
        if batch.ordinal == 0:  # every batch shares the reference and regression SNPs; batch 1 records them
            result = replace(
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
                str(Path(self.output_dir) / LD_WINDOW_SNP_COUNTS_FILE))
            log_reference_snp_funnel_warning(result.reference_snp_funnel, reference_snp_funnel_path(self.output_dir))
        LOGGER.debug(f"Computed LD scores for {len(chromosome_results)} chromosomes.")
        staged = isinstance(result, StagedLDScoreBatch)
        return result, [item.variation for item in chromosome_results] if staged else []


def _snp_universe_policy(
    *,
    ref_panel,
    ldscore_config: LDScoreConfig,
    regression_regions,
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


def _slice_annotation_bundle(annotation_bundle, chrom: str):
    """Borrow one chromosome descriptor without copying annotation payloads."""
    drops = (annotation_bundle.identity_drops if annotation_bundle.chromosome_identity_drops is None
             else annotation_bundle.chromosome_identity_drops.get(str(chrom)))
    return replace(annotation_bundle, shards={str(chrom): annotation_bundle.shard(chrom)},
                   gene_list_batch=None, identity_drops=drops, chromosome_identity_drops=None)


def _partition_annotation_drops(bundle):
    """Partition source diagnostics once, including fully removed chromosomes."""
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


def _compute_one_chromosome(
    chrom: str,
    chrom_bundle,
    ref_panel,
    ldscore_config: LDScoreConfig,
    global_config: GlobalConfig,
    *,
    include_baseline: bool = True,
    fragment_dir: str | None = None,
    export_dir: str | None = None,
    regression_snps=_WORKER_UNSET,
    regression_regions=_WORKER_UNSET,
    log_path=None,
    batch_label: str = "batch 1/1",
) -> ChromLDScoreResult | DirectChromosomeFragment | _SkippedChromosome:
    """Run one chromosome task of any query batch in a pool worker.

    ``ref_panel`` holds only configuration and resolved PLINK prefixes, so no
    live reader or file handle crosses the process boundary; ``global_config``
    carries the build the parent resolved, which is not checked again.
    ``regression_snps`` and ``regression_regions`` default to the
    pool-initializer-provided globals; direct task callers may pass them
    explicitly. Only ``EmptyReferenceSNPs`` returns a ``_SkippedChromosome``,
    carrying the error's reference-SNP funnel so the parent can tell a
    chromosome emptied by the annotation SNP grid from one emptied by an
    explicit restriction; other failures propagate across the pool boundary.
    """
    with chromosome_log(log_path, chrom, batch_label):
        if regression_snps == _WORKER_UNSET:
            regression_snps = _WORKER_STATE.get("regression_snps")
        if regression_regions == _WORKER_UNSET:
            regression_regions = _WORKER_STATE.get("regression_regions")
        try:
            return compute_direct_chromosome(
                chrom, chrom_bundle, ref_panel, ldscore_config, global_config, regression_snps, regression_regions,
                include_baseline=include_baseline, fragment_dir=fragment_dir, export_dir=export_dir,
            )
        except EmptyReferenceSNPs as exc:
            return _SkippedChromosome(str(exc), getattr(exc, "reference_snp_funnel", None))


def _init_worker(regression_snps, regression_regions=None, log_level: str = "INFO") -> None:
    """Initialize shared keys and fixed, worker-owned native thread limits."""
    _WORKER_STATE["regression_snps"] = regression_snps
    _WORKER_STATE["regression_regions"] = regression_regions
    initialize_worker(log_level)
