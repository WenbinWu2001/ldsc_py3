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

from .._annotation.overlap import OverlapContribution


@dataclass(frozen=True)
class LDScoreOverlap:
    """Labeled overlap blocks for all and common reference SNPs.

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

    def select_queries(self, query_columns: list[str]) -> "LDScoreOverlap":
        """Return the same labeled overlap data restricted to selected queries.

        Baseline columns are kept as present: all of them for a complete block,
        none for a later batch's baseline-by-query block.
        """
        queries = set(self.query_diagonal_all.index)
        all_columns = [*[name for name in self.baseline_block_all.columns if name not in queries], *query_columns]
        block_all = self.baseline_block_all.loc[:, all_columns].copy()
        block_common = None
        if self.baseline_block_common is not None:
            block_common = self.baseline_block_common.loc[:, all_columns].copy()
        diagonal_all = self.query_diagonal_all.reindex(query_columns).copy()
        diagonal_common = None
        if self.query_diagonal_common is not None:
            diagonal_common = self.query_diagonal_common.reindex(query_columns).copy()
        return LDScoreOverlap(
            block_all,
            block_common,
            diagonal_all,
            diagonal_common,
            self.total_all_reference_snps,
            self.total_common_reference_snps,
        )

    @classmethod
    def from_contribution(
        cls,
        contribution: OverlapContribution,
        baseline_columns: list[str],
        query_columns: list[str],
        *,
        include_baseline: bool = True,
    ) -> "LDScoreOverlap":
        """Label a raw :class:`OverlapContribution` with annotation names.

        ``include_baseline=False`` labels a baseline-by-query block, whose
        columns are the queries only.
        """
        all_columns = [*(baseline_columns if include_baseline else []), *query_columns]
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
    """Flatten the baseline-rows block plus query diagonal into long/tidy rows.

    Rows run over the block row by row, then over the query diagonal. A missing
    common-SNP part is ``None`` (an all-``None`` object column when no common
    blocks exist, ``NaN`` beside common values otherwise).
    """
    columns = ["row_annotation", "col_annotation", "overlap_all_snps", "overlap_common_snps"]
    block_all, block_common = overlap.baseline_block_all, overlap.baseline_block_common
    diagonal_all, diagonal_common = overlap.query_diagonal_all, overlap.query_diagonal_common
    n_rows, n_cols = block_all.shape
    if n_rows * n_cols + len(diagonal_all) == 0:
        return pd.DataFrame([], columns=columns)
    queries = diagonal_all.index.to_numpy(dtype=object)

    def common_part(values, size):
        return np.full(size, None, dtype=object) if values is None else np.asarray(values, dtype=np.float64).ravel()

    common = [common_part(None if block_common is None else block_common.to_numpy(dtype=np.float64), n_rows * n_cols),
              common_part(None if diagonal_common is None else diagonal_common.to_numpy(dtype=np.float64), len(queries))]
    common_values = (np.concatenate(common) if block_common is None and diagonal_common is None
                     else np.concatenate([np.asarray(part, dtype=np.float64) if part.dtype != object
                                          else np.full(len(part), np.nan) for part in common]))
    return pd.DataFrame({
        "row_annotation": np.concatenate([np.repeat(block_all.index.to_numpy(dtype=object), n_cols), queries]),
        "col_annotation": np.concatenate([np.tile(block_all.columns.to_numpy(dtype=object), n_rows), queries]),
        "overlap_all_snps": np.concatenate([block_all.to_numpy(dtype=np.float64).ravel(),
                                            diagonal_all.to_numpy(dtype=np.float64)]),
        "overlap_common_snps": common_values,
    }, columns=columns)


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


