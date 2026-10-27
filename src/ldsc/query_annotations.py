"""Shared gene-query coverage, viability, and aligned result finalization.

Resolution selects canonical intervals; backend adapters measure SNP support.
This module applies the same focal/control policy to both direct and indexed
measurements, and prunes scientific columns with their counts and provenance.
It does not load reference artifacts or implement numerical LD-score kernels.
"""

from __future__ import annotations

import logging
import sys
from dataclasses import dataclass, replace
from collections import Counter
from typing import Any, Sequence

import numpy as np

import pandas as pd

from ._kernel.overlap import OverlapContribution
from .errors import LDSCInternalError

LOGGER = logging.getLogger("LDSC.query_annotations")
MAX_CONSOLE_GENE_ISSUES = 10
_QUERY_REASON_EXPLANATIONS = {
    "empty_gene_list": "the supplied gene list is empty",
    "zero_resolved_genes": "no gene identifiers resolved to usable catalog entries",
    "unresolved_genes": "some gene identifiers could not be resolved",
    "partial_snp_support": "some selected genes have no retained reference-SNP support",
    "zero_annotation_snps": "the query covers no SNP rows in the evaluated annotation/reference set",
    "zero_variance_ld_scores": "the LD scores are constant on the SNP rows used for regression",
    "empty_input": "the BED file contains no intervals",
}
# Successful-run notice wording: checks performed and where skipped queries are omitted.
_NOTICE_WORDING = {
    "ldscore": ("Gene-list SNP-support and LD-score checks", "LD-score outputs"),
    "annotate": ("Gene-list SNP-support checks", "annotation outputs"),
}


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
    return batch.with_coverage(chromosomes)


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
    if batch is None or not batch.summary.input_role.eq("control").any():
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


def finalize_query_measurements(statuses, count_by_column, varying):
    """Finalize statuses from counts and exact pre-publication variation.

    ``varying[name]`` is true iff the complete regression-row column
    contains more than one distinct value, treating NaNs as one value just as
    ``Series.nunique(dropna=False)`` does. Writing workflows combine chromosome
    summaries in direct post-kernel float32 or indexed float64 precision.
    Materialized direct callers inspect the same float32 score columns.
    """
    finalized: list[QueryAnnotationStatus] = []
    retained_queries: list[str] = []
    for status in statuses:
        if status.status not in {"ok", "warning"}:
            finalized.append(status)
            continue
        if status.query not in count_by_column or status.query not in varying:
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
        if not varying[status.query]:
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
    return tuple(finalized), retained_queries


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
    varying = {} if result.query_table is None else {
        status.query: (count_by_column[status.query] > 0 and
                       pd.to_numeric(result.query_table[status.query], errors="coerce").nunique(dropna=False) > 1)
        for status in statuses if status.status in {"ok", "warning"} and status.query in count_by_column
    }
    finalized, retained_queries = finalize_query_measurements(statuses, count_by_column, varying)

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


def _log_gene_catalog_issues(issues: pd.DataFrame) -> None:
    """Summarize catalog defects after their complete diagnostic table is saved."""
    if not issues.empty:
        LOGGER.warning("Gene-coordinate catalog defects: records=%d. See diagnostics/gene_catalog_issues.tsv.", len(issues))


def _log_gene_list_rejections(batch: Any) -> None:
    """Summarize source failures, rejected rows, and intentional exclusions by role.

    Call after publishing diagnostics; physical rows and gene identities remain
    in the complete gene-list audit rather than being replayed into INFO logs.
    """
    for role, rows in batch.summary.groupby("input_role", sort=False):
        source_errors = int(rows.source_status.eq("error").sum())
        if source_errors:
            LOGGER.warning("Gene-list sources rejected: role=%s; sources=%d. See diagnostics/gene_list_resolution_summary.tsv.",
                           role, source_errors)
        rejected = rows.rejected_rows.fillna(0)
        if rejected.sum():
            LOGGER.warning("Gene-list rows rejected: role=%s; sources=%d; input rows=%d. See diagnostics/gene_list_audit.tsv.gz.",
                           role, int(rejected.gt(0).sum()), int(rejected.sum()))
        excluded = rows.excluded_genes.fillna(0)
        if excluded.sum():
            LOGGER.info("Genes intentionally excluded by region policy: role=%s; sources=%d; input rows=%d. See diagnostics/gene_list_audit.tsv.gz.",
                        role, int(excluded.gt(0).sum()), int(excluded.sum()))


def _log_gene_list_snp_support(batch: Any) -> None:
    """Summarize unsupported gene selections after publishing their row audit."""
    for role, rows in batch.summary.groupby("input_role", sort=False):
        unsupported = pd.to_numeric(rows.zero_support_genes).fillna(0)
        if unsupported.sum():
            LOGGER.warning("Genes with zero retained SNP support: role=%s; sources=%d; gene selections=%d. See diagnostics/gene_list_audit.tsv.gz.",
                           role, int(unsupported.gt(0).sum()), int(unsupported.sum()))


def _gene_list_gate_a_message(batch: Any) -> str:
    """Return one bounded, actionable catalog-preflight failure message."""
    examples, n_rejected = [], 0
    for frame in batch.audit_frames():
        rejected_rows = frame.loc[frame.disposition.eq('rejected')]
        n_rejected += len(rejected_rows)
        examples.extend(rejected_rows.head(max(0,MAX_CONSOLE_GENE_ISSUES-len(examples))).to_dict('records'))
    rejected = pd.DataFrame(examples)
    submitted = int(batch.summary["nonblank_input_rows"].fillna(0).sum())
    breakdown = "; ".join(
        f"{row.source}: {int(row.rejected_rows or 0)}/{int(row.nonblank_input_rows or 0)} rejected"
        for row in batch.summary.loc[batch.summary.rejected_rows.fillna(0).gt(0)].head(MAX_CONSOLE_GENE_ISSUES).itertuples(index=False)
        if pd.notna(row.nonblank_input_rows)
    )
    sample = "; ".join(
        f"line {row.line}, {row.input_gene!r}, {row.reason}"
        for row in rejected.head(MAX_CONSOLE_GENE_ISSUES).itertuples(index=False)
    )
    omitted = max(0, n_rejected - MAX_CONSOLE_GENE_ISSUES)
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
        for row in source_errors.head(MAX_CONSOLE_GENE_ISSUES).itertuples(index=False)
    )
    issue_text = "; ".join(value for value in (source_error_text, sample) if value)
    zero_usable_sources = batch.summary[
        batch.summary["source_status"].eq("ok")
        & batch.summary["unique_resolved_genes"].fillna(0).eq(0)
    ]
    zero_usable_text = "; ".join(
        f"{row.input_role} source {row.source!r} has zero usable genes"
        for row in zero_usable_sources.head(MAX_CONSOLE_GENE_ISSUES).itertuples(index=False)
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
        f"Gene-list catalog preflight rejected {n_rejected} of {submitted} submitted row(s) "
        f"({breakdown}). First issues: {issue_text}{omitted_text}. "
        "Complete diagnostics: diagnostics/gene_list_audit.tsv.gz and "
        "diagnostics/gene_list_resolution_summary.tsv. "
        f"Suggested audit filter: {source_filter}. Repair the listed gene-list rows first, then rerun "
        "with --overwrite into this owned output directory."
    )


def _log_query_annotation_statuses(statuses: Sequence[QueryAnnotationStatus]) -> None:
    """Log whole-run counts and reason summaries after diagnostics are written."""
    if not statuses:
        return
    counts = Counter(status.status for status in statuses)
    for (status, reason), count in Counter((item.status, item.reason) for item in statuses if item.status != "ok").items():
        consequence = (" These queries are omitted from outputs." if status == "skipped" else
                       " These queries remain in the outputs; review the recorded limitation." if status == "warning" else
                       " These queries cannot be used; correct the affected inputs.")
        LOGGER.warning(
            "Query annotation outcomes: status=%s; reason=%s; queries=%d; %s.%s "
            "See diagnostics/query_annotation_status.tsv for affected queries and reasons.",
            status, reason, count, _QUERY_REASON_EXPLANATIONS.get(reason, "the recorded query check did not pass"), consequence,
        )
    LOGGER.info("Query annotation status summary: total=%d, ok=%d, warning=%d, skipped=%d, error=%d.",
                len(statuses), counts["ok"], counts["warning"], counts["skipped"], counts["error"])


def _all_query_annotations_skipped_message(statuses: Sequence[QueryAnnotationStatus]) -> str:
    """Return the consolidated user-facing error for an unusable batch."""
    reasons = "; ".join(
        f"{status.query}: {_QUERY_REASON_EXPLANATIONS.get(status.reason, 'the recorded query check did not pass')} "
        f"({status.reason})"
        for status in statuses[:MAX_CONSOLE_GENE_ISSUES]
    )
    if len(statuses) > MAX_CONSOLE_GENE_ISSUES:
        reasons += f"; {len(statuses) - MAX_CONSOLE_GENE_ISSUES} more in the status table"
    return (
        f"The run could not continue because all {len(statuses)} requested query annotations were skipped "
        f"({reasons}). Check diagnostics/query_annotation_status.tsv for each query's skip reason. "
        "Correct unintended input problems before rerunning. "
        "Use a fresh output directory, or --overwrite if replacing this run's owned outputs is intended."
    )


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
        f"WARNING: --allow-unresolved-genes completed with {rejected} of {total} "
        f"submitted row(s) omitted: {sources}"
        + (f", plus {extra} more affected source(s)" if extra else "")
        + ". See diagnostics/gene_list_resolution_summary.tsv and diagnostics/gene_list_audit.tsv.gz."
    )
    LOGGER.warning(message.removeprefix("WARNING: "))
    print(message, file=sys.stderr)


def _emit_gene_gate_b_notice(
    statuses: Sequence[QueryAnnotationStatus],
    batch: Any | None,
    *,
    command: str,
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
    checks, outputs = _NOTICE_WORDING[command]
    displayed = affected[:MAX_CONSOLE_GENE_ISSUES]
    outcomes = ", ".join(f"{status.query}={status.reason}" for status in displayed)
    extra = max(0, len(affected) - MAX_CONSOLE_GENE_ISSUES)
    summaries = []
    if unsupported_rows:
        summaries.append(f"{unsupported_rows} gene-list row(s) have no retained reference-SNP support")
    for reason, explanation in (
        ("partial_snp_support", "retained despite some genes having no retained reference-SNP support"),
        ("zero_annotation_snps", "skipped because they overlap no retained reference SNPs"),
        ("zero_variance_ld_scores", "skipped because LD scores are constant on regression SNP rows"),
    ):
        count = sum(status.reason == reason for status in affected)
        if count:
            summaries.append(f"{count} query annotation(s) {explanation} ({reason})")
    message = (
        f"WARNING: {checks}"
        + (": " + "; ".join(summaries) if summaries else "")
        + "."
        + (f" Skipped queries are omitted from {outputs}." if any(status.status == "skipped" for status in affected) else "")
        + (f" Affected queries: {outcomes}" if outcomes else "")
        + (f", plus {extra} more affected query(s)" if extra else "")
        + ("." if outcomes else "")
        + " See diagnostics/query_annotation_status.tsv, diagnostics/gene_list_resolution_summary.tsv, "
        "and diagnostics/gene_list_audit.tsv.gz."
    )
    LOGGER.warning(message.removeprefix("WARNING: "))
    print(message, file=sys.stderr)


def _emit_bed_query_notice(statuses: Sequence[QueryAnnotationStatus], *, command: str) -> None:
    """Emit one bounded successful-run notice for skipped BED queries."""
    bed = [status for status in statuses if status.input_type == "bed"]
    skipped = [status for status in bed if status.status == "skipped"]
    if not skipped:
        return
    outcomes = ", ".join(f"{status.query}={status.reason}" for status in skipped[:MAX_CONSOLE_GENE_ISSUES])
    extra = max(0, len(skipped) - MAX_CONSOLE_GENE_ISSUES)
    reasons = "; ".join(
        f"{reason}: {_QUERY_REASON_EXPLANATIONS.get(reason, 'the recorded query check did not pass')}"
        for reason in dict.fromkeys(status.reason for status in skipped)
    )
    message = (
        f"WARNING: {len(skipped)} of {len(bed)} BED query annotation(s) were skipped and are omitted from "
        f"{_NOTICE_WORDING[command][1]}: {outcomes}"
        + (f", plus {extra} more" if extra else "")
        + f". Reasons: {reasons}. See diagnostics/query_annotation_status.tsv."
    )
    LOGGER.warning(message.removeprefix("WARNING: "))
    print(message, file=sys.stderr)
