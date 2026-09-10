"""Shared gene-query coverage, viability, and aligned result finalization.

Resolution selects canonical intervals; backend adapters measure SNP support.
This module applies the same focal/control policy to both direct and indexed
measurements, and prunes scientific columns with their counts and provenance.
It does not load reference artifacts or implement numerical LD-score kernels.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, replace
from typing import Any, Sequence

import numpy as np

import pandas as pd

from ._kernel.overlap import OverlapContribution
from .errors import LDSCInternalError

LOGGER = logging.getLogger("LDSC.query_annotations")
MAX_CONSOLE_GENE_ISSUES = 10


@dataclass(frozen=True)
class QueryAnnotationStatus:
    """Status of one requested BED or gene-list query annotation."""

    query: str
    source: str
    input_type: str
    status: str
    reason: str = ""
    n_annotation_snps: float | None = None
    details: str | None = None

    def updated(self, **changes: Any) -> "QueryAnnotationStatus":
        """Return a copy with final workflow status fields replaced."""
        return replace(self, **changes)

    def as_dict(self) -> dict[str, Any]:
        """Return the fixed query-status manifest row schema."""
        return {
            "query": self.query,
            "source": self.source,
            "input_type": self.input_type,
            "status": self.status,
            "reason": self.reason,
            "n_annotation_snps": self.n_annotation_snps,
            "details": self.details,
        }


def assess_gene_coverage(batch, chromosomes):
    """Assess unique post-exclusion selections without altering resolver truth.

    ``chromosomes`` must come from validated inputs. Support is deliberately
    left unknown here, including for genes outside the evaluated scope.
    Return the enriched batch and every focal/control coverage failure.
    """
    scope = set(map(str, chromosomes))
    summary, audit = batch.summary.copy(), batch.audit.copy()
    for column in ("coverage_status", "missing_chromosomes", "uncovered_gene_ids"):
        summary[column] = summary[column].astype(object)
    audit["coverage_status"] = audit["coverage_status"].astype(object)
    errors = []
    for selection in batch.selections:
        mask = summary.input_role.eq(selection.input_role) & summary.source_ordinal.eq(selection.source_ordinal)
        missing = [(gene, interval[0]) for gene, interval in zip(selection.canonical_gene_ids, selection.intervals, strict=True) if interval[0] not in scope]
        selected = len(selection.canonical_gene_ids)
        covered = selected - len(missing)
        status = "empty" if not selected else "full" if not missing else "none" if not covered else "partial"
        missing_chroms = ",".join(sorted({chrom for _, chrom in missing}, key=int))
        summary.loc[mask, "coverage_status"] = status
        summary.loc[mask, "selected_genes"] = selected
        summary.loc[mask, "covered_genes"] = covered
        summary.loc[mask, "uncovered_genes"] = len(missing)
        summary.loc[mask, "missing_chromosomes"] = missing_chroms
        summary.loc[mask, "uncovered_gene_ids"] = ",".join(gene for gene, _ in missing)
        rows = audit.input_role.eq(selection.input_role) & audit.source_ordinal.eq(selection.source_ordinal) & audit.canonical_gene_id.isin(selection.canonical_gene_ids)
        audit.loc[rows, "coverage_status"] = audit.loc[rows, "chrom"].map(lambda chrom: "covered" if str(chrom) in scope else "uncovered")
        if missing:
            errors.append(f"{selection.input_role} {selection.query!r}: incomplete chromosome coverage ({covered}/{selected} selected genes covered); missing chromosomes {missing_chroms}; affected genes {', '.join(gene for gene, _ in missing)}. Supply matching baseline/reference inputs covering these genes or explicitly revise the submitted list.")
    return replace(batch, summary=summary, audit=audit), errors


def gene_query_statuses(batch):
    """Apply shared resolution, coverage, and measured-support precedence."""
    statuses = []
    for selection in batch.selections:
        if selection.input_role != "focal":
            continue
        row = batch.summary.loc[batch.summary.input_role.eq("focal") & batch.summary.source_ordinal.eq(selection.source_ordinal)].iloc[0]
        rejected = int(row.rejected_rows)
        unsupported, supported = row.zero_support_genes, row.genes_with_snp_support
        coverage = row.get("coverage_status")
        details = []
        if rejected:
            details.append(f"{rejected} submitted row(s) were rejected during catalog resolution.")
        if pd.notna(unsupported) and unsupported:
            universe = "annotation-grid" if getattr(batch, "support_kind", "reference") == "annotation" else "retained reference-SNP"
            details.append(f"{int(unsupported)} resolved gene(s) have zero {universe} support.")
        if not selection.canonical_gene_ids:
            status, reason = "skipped", "empty_gene_list" if not row.nonblank_input_rows else "zero_resolved_genes"
        elif pd.notna(coverage) and coverage in {"partial", "none"}:
            status, reason = "error", "incomplete_chromosome_coverage"
            details.append(f"{int(row.covered_genes)}/{int(row.selected_genes)} selected genes covered; missing chromosomes {row.missing_chromosomes}; affected genes {row.uncovered_gene_ids}.")
        elif pd.notna(supported) and supported == 0:
            status, reason = "skipped", "zero_annotation_snps"
        elif rejected:
            status, reason = "warning", "partial_gene_resolution"
        elif pd.notna(unsupported) and unsupported:
            status, reason = "warning", "partial_snp_support"
        else:
            status, reason = "ok", ""
        if status != "ok":
            details.append("See diagnostics/gene_list_resolution_summary.tsv and diagnostics/gene_list_audit.tsv.gz.")
        statuses.append(QueryAnnotationStatus(selection.query, selection.source, "gene_list", status, reason,
                                             0.0 if reason == "zero_annotation_snps" else None,
                                             " ".join(details) or None))
    return tuple(statuses)


def gene_control_errors(batch, result=None):
    """Return required-control viability failures without dropping the control."""
    if batch is None:
        return []
    errors = []
    for selection in batch.selections:
        if selection.input_role != "control":
            continue
        row = batch.summary.loc[batch.summary.input_role.eq("control")].iloc[0]
        if not selection.canonical_gene_ids:
            errors.append("The requested control gene list has zero selected genes.")
        elif pd.notna(row.genes_with_snp_support) and row.genes_with_snp_support == 0:
            errors.append("The requested control gene list has zero annotation SNPs (zero retained-SNP support).")
        elif result is not None:
            if "gene_control" not in result.baseline_table:
                raise LDSCInternalError("The requested gene control is missing from the computed baseline LD-score table.")
            if pd.to_numeric(result.baseline_table["gene_control"], errors="coerce").nunique(dropna=False) <= 1:
                errors.append("The requested control gene list produced zero-variance LD scores on regression SNP rows.")
    if errors:
        errors.append("The conditioning model cannot be fitted. Inspect the gene-list diagnostics and supply a usable control.")
    return errors


def gene_viability_errors(batch, statuses, result=None):
    """Collect control and all-focal-skipped failures at one viability gate."""
    errors = gene_control_errors(batch, result)
    if (statuses or batch is not None) and all(status.status == "skipped" for status in statuses):
        errors.append(_all_query_annotations_skipped_message(statuses))
    return errors


def finalize_query_statuses(
    result: Any,
    statuses: Sequence[QueryAnnotationStatus],
) -> Any:
    """Attach final counts and prune zero-hit or constant query LD scores."""
    if not statuses:
        return result
    count_by_column = {
        str(record["column"]): float(record["all_reference_snp_count"])
        for record in result.count_records
        if record.get("group") == "query"
    }
    finalized: list[QueryAnnotationStatus] = []
    retained_queries: list[str] = []
    for status in statuses:
        if status.status not in {"ok", "warning"}:
            finalized.append(status)
            continue
        if status.query not in count_by_column or result.query_table is None:
            raise LDSCInternalError(
                "LD-score query-status finalization could not find the scientific column for "
                f"usable query {status.query!r}. Re-run with DEBUG logging and report the traceback."
            )
        count = count_by_column[status.query]
        if count <= 0:
            finalized.append(
                status.updated(
                    status="skipped",
                    reason="zero_annotation_snps",
                    n_annotation_snps=count,
                    details="query overlaps no retained LD-reference SNPs",
                )
            )
            continue
        values = pd.to_numeric(result.query_table[status.query], errors="coerce")
        if values.nunique(dropna=False) <= 1:
            finalized.append(
                status.updated(
                    status="skipped",
                    reason="zero_variance_ld_scores",
                    n_annotation_snps=count,
                    details="query LD scores have zero variance on regression SNP rows",
                )
            )
            continue
        finalized.append(status.updated(n_annotation_snps=count))
        retained_queries.append(status.query)

    original_queries = list(result.query_columns)
    query_table = None
    if retained_queries:
        metadata_columns = [column for column in result.query_table.columns if column not in original_queries]
        query_table = result.query_table.loc[:, [*metadata_columns, *retained_queries]].copy()
    count_records = [
        record
        for record in result.count_records
        if record.get("group") != "query" or record.get("column") in retained_queries
    ]
    overlap = result.overlap
    if overlap is not None:
        overlap = overlap.select_queries(retained_queries)
    chromosome_results = [
        _select_chromosome_result_queries(chromosome_result, retained_queries)
        for chromosome_result in result.chromosome_results
    ]
    finalized_result = replace(
        result,
        query_table=query_table,
        query_columns=retained_queries,
        count_records=count_records,
        overlap=overlap,
        chromosome_results=chromosome_results,
        query_statuses=tuple(finalized),
        annotation_types={name: kind for name, kind in result.annotation_types.items()
                          if name in [*result.baseline_columns, *retained_queries]},
    )
    finalized_result.validate()
    return finalized_result



def _select_chromosome_result_queries(
    result: Any,
    retained_queries: Sequence[str],
) -> Any:
    """Return one chromosome result restricted to final usable query columns."""
    retained = list(retained_queries)
    original = list(result.query_columns)
    query_positions = [original.index(column) for column in retained]
    query_table = None
    if retained and result.query_table is not None:
        metadata_columns = [column for column in result.query_table.columns if column not in original]
        query_table = result.query_table.loc[:, [*metadata_columns, *retained]].copy()
    keep_positions = [*range(len(result.baseline_columns)), *[len(result.baseline_columns) + pos for pos in query_positions]]
    snp_count_totals = {
        key: np.asarray(values)[keep_positions]
        for key, values in result.snp_count_totals.items()
    }
    overlap = result.overlap
    if overlap is not None:
        overlap = OverlapContribution(
            baseline_block_all=np.asarray(overlap.baseline_block_all)[:, keep_positions],
            baseline_block_common=(
                None
                if overlap.baseline_block_common is None
                else np.asarray(overlap.baseline_block_common)[:, keep_positions]
            ),
            query_diagonal_all=np.asarray(overlap.query_diagonal_all)[query_positions],
            query_diagonal_common=(
                None
                if overlap.query_diagonal_common is None
                else np.asarray(overlap.query_diagonal_common)[query_positions]
            ),
            n_all=overlap.n_all,
            n_common=overlap.n_common,
        )
    selected = replace(
        result,
        query_table=query_table,
        query_columns=retained,
        count_records=[
            record
            for record in result.count_records
            if record.get("group") != "query" or record.get("column") in retained
        ],
        snp_count_totals=snp_count_totals,
        overlap=overlap,
        annotation_types={name: kind for name, kind in result.annotation_types.items()
                          if name in [*result.baseline_columns, *retained]},
    )
    selected.validate()
    return selected


def _log_gene_list_rejections(batch: Any) -> None:
    """Log every rejected row and intentional gene-region exclusion."""
    source_errors = batch.summary[batch.summary["source_status"].eq("error")]
    for row in source_errors.itertuples(index=False):
        LOGGER.warning(
            "Gene-list source rejected: role=%s source=%s reason=%s",
            row.input_role,
            row.source,
            row.source_reasons,
        )
    rejected = batch.audit[batch.audit["disposition"].eq("rejected")]
    for row in rejected.itertuples(index=False):
        LOGGER.warning(
            "Gene-list row rejected: role=%s source=%s line=%s input_gene=%r reason=%s",
            row.input_role,
            row.source,
            row.line,
            row.input_gene,
            row.reason,
        )
    excluded = batch.audit[batch.audit["disposition"].eq("excluded")]
    for row in excluded.itertuples(index=False):
        LOGGER.info(
            "Gene intentionally excluded by region policy: role=%s source=%s line=%s "
            "input_gene=%r canonical_gene_id=%s reason=%s",
            row.input_role,
            row.source,
            row.line,
            row.input_gene,
            row.canonical_gene_id,
            row.reason,
        )


def _log_gene_list_snp_support(batch: Any) -> None:
    """Write every zero-support gene to the complete file log."""
    unsupported = batch.audit[batch.audit["disposition"].eq("unsupported")]
    for row in unsupported.itertuples(index=False):
        LOGGER.warning(
            "Gene has zero retained reference-SNP support: role=%s source=%s line=%s input_gene=%r canonical_gene_id=%s",
            row.input_role,
            row.source,
            row.line,
            row.input_gene,
            row.canonical_gene_id,
        )


def _gene_list_gate_a_message(batch: Any) -> str:
    """Return one bounded, actionable catalog-preflight failure message."""
    rejected = batch.audit[batch.audit["disposition"].eq("rejected")]
    submitted = int(batch.summary["nonblank_input_rows"].fillna(0).sum())
    breakdown = "; ".join(
        f"{row.source}: {int(row.rejected_rows or 0)}/{int(row.nonblank_input_rows or 0)} rejected"
        for row in batch.summary.itertuples(index=False)
        if pd.notna(row.nonblank_input_rows)
    )
    sample = "; ".join(
        f"line {row.line}, {row.input_gene!r}, {row.reason}"
        for row in rejected.head(MAX_CONSOLE_GENE_ISSUES).itertuples(index=False)
    )
    omitted = max(0, len(rejected) - MAX_CONSOLE_GENE_ISSUES)
    omitted_text = f"; {omitted} more rejected row(s) are in the audit" if omitted else ""
    source_errors = batch.summary[batch.summary["source_status"].eq("error")]
    def source_reason_text(reason: str, source: str, role: str) -> str:
        if any(token in source for token in ("*", "?", "[", "]")):
            return (
                "glob patterns are not allowed"
                if role == "control"
                else "glob pattern matched no files"
            )
        return reason

    source_error_text = "; ".join(
        (
            f"control gene-list file {row.source!r}: {source_reason_text(row.source_reasons, row.source, row.input_role)}"
            if row.input_role == "control"
            else f"query gene-list source {row.source!r}: {source_reason_text(row.source_reasons, row.source, row.input_role)}"
        )
        for row in source_errors.itertuples(index=False)
    )
    issue_text = "; ".join(value for value in (source_error_text, sample) if value)
    zero_usable_sources = batch.summary[
        batch.summary["source_status"].eq("ok")
        & batch.summary["unique_resolved_genes"].fillna(0).eq(0)
    ]
    zero_usable_text = "; ".join(
        f"{row.input_role} source {row.source!r} has zero usable genes"
        for row in zero_usable_sources.itertuples(index=False)
    )
    issue_text = "; ".join(value for value in (issue_text, zero_usable_text) if value)
    if not rejected.empty:
        first_rejected = rejected.iloc[0]
        source_filter = (
            f"source == {first_rejected['source']!r} & disposition == 'rejected'"
        )
    elif not source_errors.empty:
        first_source_error = source_errors.iloc[0]
        source_filter = (
            f"source == {first_source_error['source']!r} & source_status == 'error' "
            "(in gene_list_resolution_summary.tsv)"
        )
    else:
        source_filter = "disposition == 'rejected'"
    return (
        f"Gene-list catalog preflight rejected {len(rejected)} of {submitted} submitted row(s) "
        f"({breakdown}). First issues: {issue_text}{omitted_text}. "
        "Complete diagnostics: diagnostics/gene_list_audit.tsv.gz and "
        "diagnostics/gene_list_resolution_summary.tsv. "
        f"Suggested audit filter: {source_filter}. Repair the listed gene-list rows first, then rerun "
        "with --overwrite into this owned output directory."
    )


def _log_query_annotation_statuses(statuses: Sequence[QueryAnnotationStatus]) -> None:
    """Log every non-ok query and one compact batch-status summary."""
    if not statuses:
        return
    counts: dict[str, int] = {"ok": 0, "warning": 0, "skipped": 0}
    for status in statuses:
        counts[status.status] = counts.get(status.status, 0) + 1
        if status.status != "ok":
            diagnostic = (
                "diagnostics/gene_list_resolution_summary.tsv and diagnostics/gene_list_audit.tsv.gz"
                if status.input_type == "gene_list"
                else "diagnostics/query_annotation_status.tsv"
            )
            LOGGER.warning(
                f"Query annotation '{status.query}' status={status.status}, reason={status.reason}; "
                f"see {diagnostic}."
            )
    LOGGER.info(
        "Query annotation status summary: "
        f"ok={counts.get('ok', 0)}, warning={counts.get('warning', 0)}, skipped={counts.get('skipped', 0)}."
    )


def _all_query_annotations_skipped_message(statuses: Sequence[QueryAnnotationStatus]) -> str:
    """Return the consolidated user-facing error for an unusable batch."""
    reasons = ", ".join(f"{status.query}={status.reason}" for status in statuses)
    return (
        f"ldscore could not continue because all {len(statuses)} requested query annotations were skipped "
        f"({reasons}). Check diagnostics/query_annotation_status.tsv for details, correct the inputs, "
        "and rerun with --overwrite when reusing this output directory."
    )
