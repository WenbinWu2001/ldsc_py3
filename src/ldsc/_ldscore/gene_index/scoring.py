"""Scoring gene lists with a published gene LD-score index.

``run_indexed_ldscore`` resolves query and control gene lists against the
index's catalog and writes their exact LD scores without re-reading the
reference panel. Loading and verifying the index is
``ldsc._ldscore.gene_index.artifact.load_gene_ldscore_index``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Sequence

import pandas as pd

from ..._annotation.gene_lists import GeneCatalog
from ..._runtime.scratch import ScratchWorkspace
from ..outputs import INPUT_ISSUES_FILE, LDScoreDirectoryWriter, QueryDiagnostics, LDScoreOutputConfig
from .artifact import inspect_index_paths, load_gene_ldscore_index
from ..._runtime.parallelism import validate_query_batch_size, validate_workers
from ..._runtime.paths import output_directory_path
from ...errors import LDSCInputError
from ..._runtime.preflight import inspect_declared_inputs
from ..._annotation.gene_query_storage import resolve_gene_lists_staged
from ..._annotation.query_outcomes import gene_query_statuses, gene_control_errors, gene_viability_errors
from ..diagnostics import (
    fail_gene_list_resolution,
    fail_unusable_queries,
    log_query_outcomes,
    ldscore_logs,
    log_reference_snp_funnel_warning,
)

LOGGER = logging.getLogger("LDSC.gene_ldscore_index")

# Chromosome scope recorded when the index or its inputs fail validation.
_FAILED_INDEX_SCOPE = {"selection": "validated_immutable_index", "chromosomes": None,
                       "analysis_chromosomes": [], "validation_status": "failed"}


def run_indexed_ldscore(
    index_dir: str | Path,
    *,
    query_gene_list_sources: Sequence[str | Path],
    control_gene_list_file: str | Path | None = None,
    gene_list_resolution_policy: str = "strict",
    query_batch_size: int = 100,
    workers: int = 1,
    output_dir: str | Path,
    overwrite: bool = False,
    main_log_open: bool = False,
):
    """Assemble exact gene-list LD scores from one explicit index directory.

    Parameters
    ----------
    index_dir : path-like
        Exact validated single-index directory.
    query_gene_list_sources : sequence of path-like
        One-column gene-list files. Exact authoritative gene IDs and
        case-sensitive gene names are resolved against the embedded catalog in
        source order.
    control_gene_list_file : path-like or None, optional
        Optional fixed-control gene-list file. When omitted, no
        ``gene_control`` annotation is added.
    gene_list_resolution_policy : {"strict", "resolved-only"}, optional
        Default ``"strict"`` stops on rejected identifiers after auditing all
        query gene lists and the optional control. ``"resolved-only"`` omits
        allowlisted rejected rows with diagnostics; invalid indexes, incomplete
        chromosome coverage, and an unusable requested control remain fatal.
        This Python parameter retains its string values. In ``ldsc ldscore``,
        omit ``--allow-unresolved-genes`` for strict resolution or supply that
        value-free flag for ``"resolved-only"``.
    query_batch_size : int, optional
        Positive maximum query columns per operator multiplication, default 100.
        Bounds outstanding chromosome fragments and saved query-file width.
        One batch is consumed before the next is computed; the parent writes
        one chromosome at a time without a genome-wide query array.
        Chromosome operators are reloaded for each batch in a persistent pool.
        Each pathway remains a separate annotation and regression model.
    workers : int, optional
        Maximum concurrent workers, default 1, including on HPC.
        Use joblib.cpu_count() as the detected available CPU budget.
        Positive requests above that budget fail before capping by task
        count; negative -k uses max(1, budget + 1 - k). Zero, booleans,
        and non-integers are invalid. Even one worker runs in its own process.
        Every worker uses 1 BLAS thread in serial and parallel execution.
        Numerical OpenMP and Arrow CPU/I/O limits are also one; unsupported
        numerical limits are logged as unverified and computation continues.
        This does not detect idle CPUs, guarantee scheduler allocation
        compliance, or ensure sufficient memory.
    output_dir : path-like
        Destination for the canonical self-contained LD-score directory.
        User and environment tokens are expanded for both scratch and results.
    overwrite : bool, optional
        Replace workflow-owned output artifacts. Default is ``False``.
    main_log_open : bool, optional
        ``True`` when the calling command already preflighted the outputs and
        opened ``diagnostics/ldscore.log``. Default ``False``: this function
        opens that log itself while chromosomes are computed.

    Returns
    -------
    LDScoreSource
        Shared baseline/control values, counts, overlap, statuses, provenance,
        and saved paths. ``read_queries(names)`` reads explicit selections
        across the ordered query batch files; no query values are cached.

    Raises
    ------
    LDSCInputError
        If the index is invalid, a control is unusable, or every requested
        query is skipped. The all-skipped case writes diagnostics but no
        scientific LD-score tables.
    ValueError
        Query batch width is not a positive integer or the worker request is
        not a nonzero integer.
    FileExistsError
        An owned output artifact exists and overwrite is disabled.

    Notes
    -----
    Overlapping, nested, duplicated, and alias-selected genes are combined by
    Boolean union. A persistent worker pool processes one active query batch;
    each chromosome task reloads its operator and releases it after projection.
    Private fragments preserve float64 through global query-status checks.
    Assembly reads one chromosome at a time in index order, narrows to the
    public float32 format, and removes consumed fragments before the next batch.
    Completed batches remain private until the complete run succeeds.
    One query file is named ``ldscore.query.parquet``; multiple files use
    numbered names in the required ``query_batches`` manifest. This API writes
    outputs even for a single batch; batch resumption is not supported.
    """
    validate_query_batch_size(query_batch_size)
    validate_workers(workers)
    output_dir = str(output_directory_path(output_dir))
    LDScoreDirectoryWriter.artifact_family(output_dir).preflight(overwrite=overwrite)
    try:
        inspect_declared_inputs(files=[('query genes', query_gene_list_sources), ('control genes', control_gene_list_file)],
            checks=[('index', index_dir, lambda: inspect_index_paths(index_dir))],
            issues_path=Path(output_dir)/INPUT_ISSUES_FILE)
    except LDSCInputError as exc:
        LDScoreDirectoryWriter.write_input_preflight_diagnostics(
            QueryDiagnostics(_FAILED_INDEX_SCOPE, input_issues=exc.input_issues), output_dir)
        raise
    with ScratchWorkspace(output_dir) as workspace:
        return _run_indexed_ldscore_in_workspace(
            index_dir, query_gene_list_sources=query_gene_list_sources,
            control_gene_list_file=control_gene_list_file,
            gene_list_resolution_policy=gene_list_resolution_policy,
            query_batch_size=query_batch_size, workers=workers, output_dir=output_dir,
            overwrite=overwrite, workspace=workspace, main_log_open=main_log_open,
        )


def _run_indexed_ldscore_in_workspace(
    index_dir, *, query_gene_list_sources, control_gene_list_file,
    gene_list_resolution_policy, query_batch_size, workers, output_dir, overwrite,
    workspace, main_log_open=False,
):
    """Own staged gene selections until the canonical result has been written."""
    try:
        index = load_gene_ldscore_index(index_dir)
    except (OSError, ValueError, LDSCInputError) as exc:
        from ..preflight import ISSUE_COLUMNS

        issues = getattr(exc, "input_issues", pd.DataFrame([{
            "input_role": "index", "source": str(index_dir), "chrom": "", "reason": "invalid_required_input",
            "details": str(exc), "repair": "Restore or rebuild the complete immutable autosomes 1-22 index.",
        }], columns=ISSUE_COLUMNS))
        LDScoreDirectoryWriter().write_query_diagnostics(
            QueryDiagnostics(_FAILED_INDEX_SCOPE, input_issues=issues),
            LDScoreOutputConfig(output_dir=output_dir, overwrite=overwrite))
        raise
    projection_build = str(index.index_identity.get("projection_build"))
    gene_policy = str(index.index_identity.get("gene_exclude_regions", "none"))
    catalog = GeneCatalog.from_embedded_frame(index.gene_catalog)
    if catalog.genome_build != projection_build:
        raise LDSCInputError("Embedded gene coordinate catalog build disagrees with the index projection build.")
    batch = resolve_gene_lists_staged(
        query_gene_list_sources,
        catalog,
        workspace,
        existing_columns=index.baseline_columns,
        control_path=control_gene_list_file,
        resolution_policy=gene_list_resolution_policy,
        gene_exclude_regions=gene_policy,
    )
    output_config = LDScoreOutputConfig(output_dir=output_dir, overwrite=overwrite)
    if batch.has_fatal_gate_a_issues:
        fail_gene_list_resolution(batch, output_config)
    batch, coverage_errors = batch.with_coverage(index.chromosomes)
    scope = {"chromosomes": list(index.chromosomes), "analysis_chromosomes": list(index.chromosomes),
             "selection": "validated_immutable_index", "index_id": index.index_id,
             "validation_status": "passed"}
    LOGGER.info("Chromosomes resolved and entering the analysis: %s.", ", ".join(index.chromosomes))
    coverage_errors.extend(gene_control_errors(batch))
    if coverage_errors:
        LDScoreDirectoryWriter().write_query_diagnostics(QueryDiagnostics(
            {**scope, "validation_status": "failed", "analysis_chromosomes": []},
            query_statuses=gene_query_statuses(batch), gene_list_batch=batch), output_config)
        raise LDSCInputError("; ".join(coverage_errors))
    batch = batch.with_snp_support(index.gene_support)
    statuses = gene_query_statuses(batch)
    usable = [
        declaration
        for declaration, status in zip(
            (item for item in batch.declarations if item['input_role'] == "focal"),
            statuses,
            strict=True,
        )
        if status.status in {"ok", "warning"}
    ]
    control_resolution = batch.selection("control", 0) if control_gene_list_file is not None else None
    control_errors = gene_viability_errors(batch, statuses)
    if control_errors:
        fail_unusable_queries(QueryDiagnostics(scope, query_statuses=statuses, gene_list_batch=batch),
                              "; ".join(control_errors), output_config)
    if "gene_control" in index.baseline_columns:
        raise LDSCInputError(
            "Index supplied baseline columns contain reserved name 'gene_control'. "
            "Rebuild the index with a distinct baseline annotation name; gene_control is reserved for the requested gene control."
        )
    from .batches import indexed_results

    with ldscore_logs(output_dir, index.chromosomes, overwrite=overwrite,
                      log_level=logging.getLevelName(LOGGER.getEffectiveLevel()),
                      main_log_open=main_log_open) as log_paths:
        results = indexed_results(index, batch, usable, control_resolution, query_batch_size,
                                  workspace.path, workers, statuses, scope, log_paths)
        result = LDScoreDirectoryWriter().write_batches(
            results, output_config, query_status_order=[status.query for status in statuses])
        log_reference_snp_funnel_warning(index.reference_snp_funnel, result.output_paths["reference_snp_funnel"])
        log_query_outcomes(result.gene_list_batch, result.query_statuses)
        return result
