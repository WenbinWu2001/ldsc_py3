"""Public-layer overlap-matrix container, (de)serialization, and assembly.

The labeled in-memory form keeps the baseline-rows block as a DataFrame indexed
by baseline annotation name with all annotation columns, plus a query
self-overlap Series. Persistence is a long/tidy parquet so the artifact scales
linearly with the query-annotation count and aligns by name.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from ._kernel import regression as reg
from ._kernel.overlap import OverlapContribution


@dataclass(frozen=True)
class LDScoreOverlap:
    """Labeled overlap blocks for both SNP universes.

    ``baseline_block_*`` is indexed by baseline annotation name with columns over
    all annotations (``baseline + query``). ``query_diagonal_*`` is indexed by
    query annotation name. ``total_*_reference_snps`` are the ``M_tot`` scalars.
    """
    baseline_block_all: pd.DataFrame
    baseline_block_common: pd.DataFrame | None
    query_diagonal_all: pd.Series
    query_diagonal_common: pd.Series | None
    total_all_reference_snps: float
    total_common_reference_snps: float | None

    @classmethod
    def from_contribution(
        cls,
        contribution: OverlapContribution,
        baseline_columns: list[str],
        query_columns: list[str],
    ) -> "LDScoreOverlap":
        """Label a raw :class:`OverlapContribution` with annotation names."""
        all_columns = [*baseline_columns, *query_columns]
        block_all = pd.DataFrame(
            np.asarray(contribution.baseline_block_all, dtype=np.float64),
            index=baseline_columns,
            columns=all_columns,
        )
        qdiag_all = pd.Series(np.asarray(contribution.query_diagonal_all, dtype=np.float64), index=query_columns)
        block_common = qdiag_common = total_common = None
        if contribution.baseline_block_common is not None:
            block_common = pd.DataFrame(
                np.asarray(contribution.baseline_block_common, dtype=np.float64),
                index=baseline_columns,
                columns=all_columns,
            )
            qdiag_common = pd.Series(
                np.asarray(contribution.query_diagonal_common, dtype=np.float64), index=query_columns
            )
            total_common = float(contribution.n_common)
        return cls(block_all, block_common, qdiag_all, qdiag_common, float(contribution.n_all), total_common)


def overlap_to_long_frame(overlap: LDScoreOverlap) -> pd.DataFrame:
    """Flatten the baseline-rows block plus query diagonal into long/tidy rows."""
    rows = []
    block_all = overlap.baseline_block_all
    block_common = overlap.baseline_block_common
    for row_name in block_all.index:
        for col_name in block_all.columns:
            rows.append({
                "row_annotation": row_name,
                "col_annotation": col_name,
                "overlap_all_snps": float(block_all.at[row_name, col_name]),
                "overlap_common_snps": (
                    None if block_common is None else float(block_common.at[row_name, col_name])
                ),
            })
    for q in overlap.query_diagonal_all.index:
        rows.append({
            "row_annotation": q,
            "col_annotation": q,
            "overlap_all_snps": float(overlap.query_diagonal_all.at[q]),
            "overlap_common_snps": (
                None if overlap.query_diagonal_common is None else float(overlap.query_diagonal_common.at[q])
            ),
        })
    return pd.DataFrame(rows, columns=["row_annotation", "col_annotation", "overlap_all_snps", "overlap_common_snps"])


def _clean_block(pivoted: pd.DataFrame, baseline_columns: list[str], all_columns: list[str]) -> pd.DataFrame:
    """Reindex a pivoted block to canonical order and drop axis names."""
    return pivoted.reindex(index=baseline_columns, columns=all_columns).rename_axis(index=None, columns=None)


def _clean_diagonal(series: pd.Series, query_columns: list[str]) -> pd.Series:
    """Reindex a diagonal series to canonical order and drop the index/name."""
    return series.reindex(query_columns).rename_axis(None).rename(None)


def overlap_from_long_frame(
    frame: pd.DataFrame,
    baseline_columns: list[str],
    query_columns: list[str],
    total_all_reference_snps: float,
    total_common_reference_snps: float | None,
) -> LDScoreOverlap:
    """Rebuild labeled blocks from the long/tidy parquet frame."""
    all_columns = [*baseline_columns, *query_columns]
    baseline_set = set(baseline_columns)
    block_rows = frame[frame["row_annotation"].isin(baseline_set)]
    diag_rows = frame[~frame["row_annotation"].isin(baseline_set)].set_index("row_annotation")
    block_all = _clean_block(
        block_rows.pivot(index="row_annotation", columns="col_annotation", values="overlap_all_snps"),
        baseline_columns, all_columns,
    )
    qdiag_all = _clean_diagonal(diag_rows["overlap_all_snps"], query_columns)
    has_common = total_common_reference_snps is not None and frame["overlap_common_snps"].notna().any()
    block_common = qdiag_common = None
    if has_common:
        block_common = _clean_block(
            block_rows.pivot(index="row_annotation", columns="col_annotation", values="overlap_common_snps"),
            baseline_columns, all_columns,
        )
        qdiag_common = _clean_diagonal(diag_rows["overlap_common_snps"], query_columns)
    return LDScoreOverlap(
        block_all, block_common, qdiag_all, qdiag_common,
        float(total_all_reference_snps),
        None if total_common_reference_snps is None else float(total_common_reference_snps),
    )


def overlap_aware_category_table(hsq, overlap_matrix, m_annot, m_tot, category_names) -> pd.DataFrame:
    """Legacy overlap-aware category table plus LDSC3 augmentation columns.

    The overlap-aware columns are produced verbatim by the ported
    ``Hsq._overlap_output``. This function then renames the z-score column,
    appends the one-sided ``Coefficient_p`` (``H1: coefficient > 0``), the
    conditional ``Category_h2`` / ``Category_h2_std_error`` (from ``hsq.cat`` /
    ``hsq.cat_se``), and a per-model ``overlap_aware`` flag.
    """
    overlap_matrix = np.asarray(overlap_matrix, dtype=np.float64)
    table = reg.Hsq._overlap_output(
        hsq,
        list(category_names),
        overlap_matrix,
        np.asarray(m_annot, dtype=np.float64).reshape(1, -1),
        float(m_tot),
        True,
    ).reset_index(drop=True)
    table = table.rename(columns={"Coefficient_z-score": "Coefficient_z"})
    # Legacy emits the string 'NA' for a category that contains every SNP; use a numeric NaN.
    table["Enrichment_p"] = pd.to_numeric(table["Enrichment_p"], errors="coerce")
    coef = np.ravel(hsq.coef).astype(np.float64)
    coef_se = np.ravel(hsq.coef_se).astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(coef_se > 0, coef / coef_se, np.nan)
    table["Coefficient_p"] = stats.norm.sf(z)
    table["Category_h2"] = np.ravel(hsq.cat).astype(np.float64)
    table["Category_h2_std_error"] = np.ravel(hsq.cat_se).astype(np.float64)
    off_diagonal = overlap_matrix.copy()
    np.fill_diagonal(off_diagonal, 0.0)
    table["overlap_aware"] = bool(np.any(off_diagonal > 0))
    return table


def assemble_model_overlap(
    overlap: LDScoreOverlap, retained_columns: list[str], use_common: bool
) -> np.ndarray:
    """Build the ``K×K`` model overlap matrix for ``retained_columns`` by name.

    Baseline-baseline and baseline-query entries come from the baseline-rows
    block (using symmetry for query-baseline); the query-query diagonal comes
    from the stored self-overlap. ``use_common`` selects the common-SNP universe.
    """
    block = overlap.baseline_block_common if use_common else overlap.baseline_block_all
    qdiag = overlap.query_diagonal_common if use_common else overlap.query_diagonal_all
    if block is None:
        raise ValueError("common-universe overlap requested but not stored")
    baseline_index = set(block.index)
    k = len(retained_columns)
    out = np.zeros((k, k), dtype=np.float64)
    for a, ca in enumerate(retained_columns):
        for b, cb in enumerate(retained_columns):
            if ca in baseline_index:
                out[a, b] = float(block.at[ca, cb])
            elif cb in baseline_index:
                out[a, b] = float(block.at[cb, ca])
            else:  # both query -> only the diagonal ca == cb is reachable per model
                out[a, b] = float(qdiag.at[ca])
    return out
