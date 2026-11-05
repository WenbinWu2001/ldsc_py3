"""Run a partitioned-h2 batch: one trait, shared preparation, one model per query.

``estimate_partitioned_h2_batch`` prepares the trait and the baseline LD
scores once, then fits the complete baseline model for a baseline-only
LD-score directory, or one baseline-plus-query model per query annotation,
loading query columns in bounded batches and fitting each model in a query
worker (``partitioned_parallel``) through ``fitting.fit_partitioned_query``.
Query statuses follow the requested query order; the final summary sort
fixes the published row and folder order. The batch publishes its result
directory through ``PartitionedH2DirectoryWriter``. This module also logs
the analysis regime of a run and sorts the summary.
"""

from __future__ import annotations

from dataclasses import replace
import logging
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from .._annotation.storage import AnnotationWorkspace
from .._ldscore.results import LDScoreResult
from .._ldscore.source import LDScoreSource
from .._runtime.parallelism import resolve_worker_count, validate_query_batch_size, validate_workers
from .._runtime.paths import output_directory_path
from .._runtime.progress import PhaseProgress
from .._runtime.serial_calculation import serial_calculation
from ..errors import LDSCInputError, LDSCInternalError, LDSCUsageError
from . import partitioned_parallel
from .outputs import PartitionedH2DirectoryWriter, PartitionedH2OutputConfig, published_query_artifacts
from .results import PartitionedH2BatchResult
from .alignment import (
    _PartitionedModelInputs,
    _PreparedLDScores,
    _count_totals_for_columns,
    _dataset_from_prepared,
    _log_regression_identity,
    _prepare_h2_inputs,
    _prepare_trait_sumstats,
    _read_query_batch,
    _select_count_key,
    _trait_label,
    _validate_partitioned_query_columns,
)
from .fitting import (
    PARTITIONED_H2_REQUIRES_OVERLAP_MESSAGE,
    _absolute_ldscore_dir,
    _coefficient_delete_frame,
    _config_prevalence,
    _fit_h2_dataset,
    _fit_metadata,
    _genome_build_of,
    _intercept_policies,
    _ldscore_dir_of,
    _log_phase_timing,
    fit_partitioned_query,
    summarize_partitioned_h2,
)


LOGGER = logging.getLogger("LDSC.regression_runner")


PARTITIONED_H2_SUMMARY_SORT_COLUMNS = {
    "category": "category",
    "prop-snps": "prop_snps",
    "prop-h2": "prop_h2",
    "enrichment": "enrichment",
    "enrichment-p": "enrichment_p",
    "coefficient": "coefficient",
    "coefficient-p": "coefficient_p",
}
PARTITIONED_H2_SUMMARY_ASCENDING_SORTS = {
    "enrichment-p",
    "coefficient-p",
}


@serial_calculation
def _fit_functional_model(dataset, config):
    """Fit and summarize the complete baseline-only model in one child."""
    outcome = _fit_h2_dataset(dataset, config)
    samp_prev, pop_prev = _config_prevalence(config)
    summary = summarize_partitioned_h2(
        outcome.estimator, outcome.dataset, outcome.dataset.retained_ld_columns,
        samp_prev=samp_prev, pop_prev=pop_prev,
    )
    return outcome, summary


def estimate_partitioned_h2_batch(
    sumstats_table, ldscore_result, config, global_config, *, output_dir, query_columns=None,
    query_batch_size=100, overwrite=False, summary_sort_by="auto", metadata=None,
    continue_on_query_error=False, workers=1,
) -> PartitionedH2BatchResult:
    """Fit and write separate pathway models; see ``RegressionRunner.estimate_partitioned_h2_batch``.

    ``config`` holds the statistical settings and ``global_config`` is the
    SNP identifier fallback for an input without recorded provenance. The
    run still publishes its result directory, as the CLI and the facade
    method do.
    """
    validate_workers(workers)
    validate_query_batch_size(query_batch_size, error_type=LDSCUsageError)
    if not isinstance(continue_on_query_error, bool):
        raise ValueError("continue_on_query_error must be a boolean.")
    if ldscore_result.overlap is None:
        raise LDSCInputError(PARTITIONED_H2_REQUIRES_OVERLAP_MESSAGE)
    queries = list(ldscore_result.query_columns if query_columns is None else query_columns)
    queries = _validate_partitioned_query_columns(ldscore_result, queries)
    # Check requested counts before fitting: baseline counts are shared by
    # every model; query counts too unless continuation was requested, which
    # records a query lacking its counts as a failed query instead.
    model_columns = [*ldscore_result.baseline_columns, *([] if continue_on_query_error else queries)]
    count_key = _select_count_key(_count_totals_for_columns(ldscore_result.count_records, model_columns),
                                  config.use_common_counts, model_columns)
    worker_count = resolve_worker_count(workers, max(1, len(queries)),
                                         task_unit="query-plus-baseline model", batch_size=query_batch_size)
    sort_by = _resolve_summary_sort(summary_sort_by, has_queries=bool(queries))
    output_dir = str(output_directory_path(output_dir))
    writer = PartitionedH2DirectoryWriter()
    writer.artifact_family(output_dir, write_per_query_results=bool(queries),
                           coefficient_delete_values=not queries, query_status=bool(queries)).preflight(overwrite=overwrite)
    metadata = dict(metadata or {})
    defaults = {
        "trait_name": _trait_label(sumstats_table),
        "count_kind": "common" if config.use_common_counts else "all",
        "analysis_type": "cell_type_specific" if queries else "functional_category",
        "headline_metric": "coefficient" if queries else "enrichment",
        "enrichment_p_test": "two_sided_t", "coefficient_p_test": "one_sided_greater",
        "annotation_types": dict(ldscore_result.annotation_types),
    }
    ldscore_dir = _ldscore_dir_of(ldscore_result)
    if ldscore_dir:
        defaults["ldscore_dir"] = ldscore_dir
    for name, value in defaults.items():
        metadata.setdefault(name, value)
    if metadata.get("ldscore_dir"):
        metadata["ldscore_dir"] = _absolute_ldscore_dir(metadata["ldscore_dir"])
    with AnnotationWorkspace(output_dir) as workspace:
        with _log_phase_timing("shared regression preparation"):
            trait = _prepare_trait_sumstats(sumstats_table, ldscore_result)
            prepared = _prepare_h2_inputs(trait, _PreparedLDScores(ldscore_result), config, global_config)
            _log_regression_identity([sumstats_table], ldscore_result, global_config,
                                     [(prepared.dataset.effective_snp_identifier,
                                       prepared.dataset.identity_downgrade_applied)])
            # The trait's provenance, which h2 metadata records as well.
            metadata.update(
                sumstats_file=sumstats_table.source_path, genome_build=_genome_build_of(prepared.dataset),
                effective_snp_identifier=prepared.dataset.effective_snp_identifier,
                identity_downgrade_applied=prepared.dataset.identity_downgrade_applied,
                count_key_used_for_regression=count_key,
                dropped_nonfinite_statistic_rows=trait.dropped_nonfinite_statistic_rows,
                intercept_h2_policy=_intercept_policies(config)["intercept_h2_policy"],
            )
            # Fits need only the aligned rows; release the trait's keyed copy.
            legacy_drops, trait = trait.legacy_drops, None
        if not queries:
            dataset = _dataset_from_prepared(prepared, config)
            if len(dataset.retained_ld_columns) == 1:
                LOGGER.warning(
                    "partitioned-h2 functional regime retained a single annotation (%s); the fit is "
                    "degenerate and collapses to single-annotation h2 behavior (two-step estimator, no "
                    "default chi-square cap). Most likely the LD-score directory has only one baseline "
                    "annotation. Provide multiple baseline annotations for a meaningful partitioned analysis.",
                    dataset.retained_ld_columns[0],
                )
            outcome, summary = _fit_functional_model(dataset, config)
            summary = _sort_partitioned_h2_summary(summary, sort_by)
            with _log_phase_timing("output writing"):
                paths = writer.write(
                    summary, PartitionedH2OutputConfig(output_dir=output_dir, overwrite=overwrite),
                    metadata={**metadata, **_fit_metadata(outcome)},
                    coefficient_delete_values=_coefficient_delete_frame(outcome.estimator, outcome.dataset.retained_ld_columns),
                    legacy_sumstats_drops=legacy_drops,
                )
            return PartitionedH2BatchResult(summary, paths, legacy_sumstats_drops=legacy_drops)

        rows, private_artifacts, statuses = [], {}, []
        first_error = None
        metadata.update(query_workers_requested=workers, query_workers_effective=worker_count,
                        query_worker_native_threads=1)
        # Workers read the model's numeric columns, SNP and the panel coordinates only.
        fit_columns = {*prepared.dataset.ref_ld_columns, prepared.dataset.weight_column, "Z", "N", "SNP",
                       "_ldsc_panel_CHR", "_ldsc_panel_POS"}
        shared = prepared.dataset.merged
        models = _PartitionedModelInputs(
            replace(prepared.dataset, merged=pd.DataFrame(
                {name: shared[name] for name in shared.columns if name in fit_columns}, copy=False)),
            prepared.counts_by_column, prepared.ld_variances)
        del shared
        with _log_phase_timing("query worker preparation"):
            worker_pool = partitioned_parallel._QueryWorkers(
                workspace.path, models, fit_partitioned_query, config, worker_count)
        total_batches = (len(queries) + query_batch_size - 1) // query_batch_size
        del models
        # Loading a batch needs only the LD-score source and row map, not the aligned frame.
        prepared = replace(prepared, dataset=replace(prepared.dataset, merged=pd.DataFrame()))
        with worker_pool as pool, PhaseProgress(LOGGER, "computation", "query fitting", len(queries)) as progress:
            for start in range(0, len(queries), query_batch_size):
                ordinal = start // query_batch_size + 1
                with PhaseProgress(LOGGER, "batch", f"query fitting {ordinal}/{total_batches}",
                                   min(query_batch_size, len(queries) - start)) as batch_progress:
                    batch_columns = queries[start:start + query_batch_size]
                    with _log_phase_timing("query batch loading", level=logging.DEBUG):
                        batch = _read_query_batch(prepared, batch_columns)
                    with TemporaryDirectory(prefix="query-batch-", dir=workspace.path) as batch_root:
                        with _log_phase_timing("query batch mapping", level=logging.DEBUG):
                            mapped = partitioned_parallel._stage_query_batch(batch, Path(batch_root) / "values")
                        batch = None
                        tasks = [(query, mapped[query], workspace.path / f"fit-{start + index}")
                                 for index, query in enumerate(batch_columns)]
                        with _log_phase_timing("query batch execution", level=logging.DEBUG):
                            outcomes = pool.map(tasks, progress=progress)
                    for focal, artifacts, status, error in outcomes:
                        statuses.append(status)
                        if error is not None:
                            first_error = first_error or error
                        else:
                            rows.append(focal)
                            private_artifacts[status["query_annotation"]] = artifacts
                    del outcomes, batch
                    batch_progress.completed = len(batch_columns)
        del prepared
        query_status = pd.DataFrame(statuses)
        n_failed = len(queries) - len(rows)
        metadata.update(query_error_policy="continue" if continue_on_query_error else "raise",
                        n_queries_requested=len(queries), n_queries_successful=len(rows), n_queries_failed=n_failed)
        if first_error is not None:
            LOGGER.warning("Partitioned-h2 query scan: %d successful, %d failed. See diagnostics/query_status.tsv "
                           "and diagnostics/partitioned-h2.log.", len(rows), n_failed)
            if not continue_on_query_error or not rows:
                writer.write_query_status(output_dir, query_status)
                first_error.add_note(
                    f"{n_failed}/{len(queries)} queries failed; no results published. "
                    + ("No query fit succeeded." if not rows else
                       "Use --continue-on-query-error to publish successful query fits.")
                )
                raise first_error
        summary = _sort_partitioned_h2_summary(pd.concat(rows, ignore_index=True), sort_by)
        with _log_phase_timing("output writing"):
            paths = writer.write(summary, PartitionedH2OutputConfig(output_dir=output_dir, overwrite=overwrite,
                                                                   write_per_query_results=True),
                                 per_query_artifacts=private_artifacts, metadata=metadata, query_status=query_status,
                                 legacy_sumstats_drops=legacy_drops)
        return PartitionedH2BatchResult(summary, paths, published_query_artifacts(output_dir, summary), query_status,
                                        legacy_drops)


def _resolve_summary_sort(sort_by: str, *, has_queries: bool) -> str:
    """Resolve the ``auto`` summary sort to a regime-appropriate column key.

    Cell-type runs (query annotations present) surface the most significant query
    first by ``coefficient-p``; functional runs preserve baseline order with
    ``category``. Any explicit choice is returned unchanged.
    """
    if sort_by != "auto":
        return sort_by
    return "coefficient-p" if has_queries else "category"


def _log_partitioned_h2_regime(ldscore_result: LDScoreSource | LDScoreResult, has_queries: bool) -> None:
    """Log a one-block banner naming the regime and the column to focus on."""
    if has_queries:
        LOGGER.info(
            "Cell-type-specific regime: baseline + one query per model (%d queries). "
            "Focus on `coefficient` together with the one-sided `coefficient_p` "
            "(a one-sided test of whether `coefficient` > 0): a positive `coefficient` with a small "
            "`coefficient_p` means the query annotation contributes heritability beyond the baseline "
            "annotations. The `enrichment` column is confounded by the query annotation's overlap with "
            "the baseline annotations, so use `coefficient_p` to judge whether the additional "
            "contribution is significant.",
            len(ldscore_result.query_columns),
        )
    else:
        LOGGER.info(
            "Functional-category regime: joint fit of %d baseline annotations. "
            "Focus on `enrichment` and `enrichment_p`: `enrichment > 1` means the annotation's SNPs "
            "explain a larger share of heritability than their share of SNPs (`< 1` means a smaller "
            "share); a small `enrichment_p` indicates the enrichment is significantly different from 1, "
            "i.e. significantly larger or smaller.",
            len(ldscore_result.baseline_columns),
        )


def _log_quantitative_annotation_interpretation(ldscore_result: LDScoreResult) -> None:
    """Warn once when a fitted model contains continuous annotations."""
    quantitative = [
        name for name, annotation_type in ldscore_result.annotation_types.items()
        if annotation_type == "quantitative"
    ]
    if not quantitative:
        return
    LOGGER.warning(
        "Continuous annotation(s) detected: count=%d; names are recorded in LD-score annotation_types metadata. Numerical prop_snps, prop_h2, and enrichment summaries "
        "remain visible for compatibility, but they are weighted summaries and do not have the ordinary "
        "binary-annotation interpretation. Interpret the coefficient (`coefficient`), its SE, and its coefficient zero-test; "
        "use `ldsc quantile-h2` for post-fit quantile enrichment. To estimate distinct conditional coefficients "
        "for bins, construct binary quantile annotations and refit partitioned-h2.",
        len(quantitative),
    )


def _sort_partitioned_h2_summary(summary: pd.DataFrame, sort_by: str = "category") -> pd.DataFrame:
    """Return ``summary`` ordered by a public partitioned-h2 summary column."""
    if sort_by == "category":
        return summary.reset_index(drop=True)
    column = PARTITIONED_H2_SUMMARY_SORT_COLUMNS.get(sort_by)
    if column is None:
        valid = ", ".join(PARTITIONED_H2_SUMMARY_SORT_COLUMNS)
        raise LDSCUsageError(
            f"partitioned-h2 summary sort key {sort_by!r} is not supported. "
            f"Choose one of: {valid}."
        )
    if column not in summary.columns:
        raise LDSCInternalError(
            f"partitioned-h2 summary cannot be sorted by {sort_by!r} because required column {column!r} "
            "is missing. Most likely the regression workflow returned a result table with the wrong schema. "
            "Re-run with DEBUG logging and report the traceback."
        )
    return summary.sort_values(
        by=column,
        ascending=sort_by in PARTITIONED_H2_SUMMARY_ASCENDING_SORTS,
        na_position="last",
        kind="mergesort",
    ).reset_index(drop=True)
