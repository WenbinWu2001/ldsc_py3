"""Annotation overlap-matrix blocks computed during LD-score estimation.

The overlap matrix ``O = A.T @ A`` (``A`` the SNP-by-annotation matrix) drives
legacy overlap-aware partitioned-heritability summaries. Because partitioned heritability
fits ``baseline + one query`` per model, only the baseline-rows block
``A_B.T @ A`` (shape ``B x (B+Q)``) and each query's self-overlap ``sum_s A[s,q]**2``
are ever needed; cross-query overlaps are not. This module computes those blocks
for all reference SNPs and for the common reference SNPs, with no I/O.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class OverlapContribution:
    """Per-chromosome (or aggregated) overlap blocks for all and common reference SNPs.

    Attributes
    ----------
    baseline_block_all, baseline_block_common : numpy.ndarray or None
        ``A_B.T @ A`` over all or common reference SNPs, shape ``(B, B+Q)`` with columns
        ordered ``baseline + query``. The common block is ``None`` when MAF is
        unavailable.
    query_diagonal_all, query_diagonal_common : numpy.ndarray or None
        Query self-overlaps ``sum_s A[s,q]**2``, shape ``(Q,)``.
    n_all, n_common : int or None
        Reference-SNP counts (``M_tot``); ``n_common`` is ``None`` without MAF.
    """
    baseline_block_all: np.ndarray
    baseline_block_common: np.ndarray | None
    query_diagonal_all: np.ndarray
    query_diagonal_common: np.ndarray | None
    n_all: int
    n_common: int | None


def compute_overlap(
    metadata: pd.DataFrame,
    annotations: pd.DataFrame,
    n_baseline: int,
    common_maf_min: float = 0.05,
) -> OverlapContribution:
    """Compute baseline-rows overlap blocks and query self-overlaps.

    ``annotations`` columns are ordered ``baseline + query`` with ``n_baseline``
    leading baseline columns. The common reference SNPs have ``MAF >= common_maf_min``
    (inclusive, matching the common-count mask); their blocks are omitted when MAF metadata
    is absent or all-missing.
    """
    A = annotations.to_numpy(dtype=np.float64, copy=False)
    A_baseline = A[:, :n_baseline]
    block_all = A_baseline.T @ A
    query_diag_all = np.einsum("sq,sq->q", A[:, n_baseline:], A[:, n_baseline:])
    n_all = int(A.shape[0])
    if "MAF" not in metadata.columns or metadata["MAF"].isna().all():
        return OverlapContribution(block_all, None, query_diag_all, None, n_all, None)
    common = (metadata["MAF"] >= common_maf_min).to_numpy()
    A_common = A[common]
    block_common = A_common[:, :n_baseline].T @ A_common
    query_diag_common = np.einsum(
        "sq,sq->q", A_common[:, n_baseline:], A_common[:, n_baseline:]
    )
    return OverlapContribution(
        block_all, block_common, query_diag_all, query_diag_common, n_all, int(common.sum())
    )


def sum_overlap_contributions(contributions: list[OverlapContribution]) -> OverlapContribution:
    """Sum per-chromosome overlap contributions into one aggregate.

    The common-SNP blocks survive only if every contribution carries them; if any
    chromosome lacks MAF, the aggregate common blocks are ``None`` (mirroring how
    aggregated common-SNP counts are dropped when unavailable).
    """
    block_all = np.sum([c.baseline_block_all for c in contributions], axis=0)
    query_diag_all = np.sum([c.query_diagonal_all for c in contributions], axis=0)
    n_all = int(sum(c.n_all for c in contributions))
    has_common = all(c.baseline_block_common is not None for c in contributions)
    if not has_common:
        return OverlapContribution(block_all, None, query_diag_all, None, n_all, None)
    block_common = np.sum([c.baseline_block_common for c in contributions], axis=0)
    query_diag_common = np.sum([c.query_diagonal_common for c in contributions], axis=0)
    n_common = int(sum(c.n_common for c in contributions))
    return OverlapContribution(
        block_all, block_common, query_diag_all, query_diag_common, n_all, n_common
    )


def annotation_statistics(metadata, annotations, n_baseline, *, query_batch_size=100,
                          common_maf_min=0.05, read_budget_bytes=16*1024*1024, include_baseline=True):
    """Accumulate counts, overlap blocks, and classifications in bounded tiles.

    Counts and products accumulate in float64 so column width and row tiling
    do not amplify float32 cancellation. Only baseline-by-query overlaps and query diagonals are retained;
    different query annotations are never crossed. All rows are reference SNPs,
    independently of which SNPs receive output LD scores.

    ``include_baseline=False`` skips the baseline group: counts and
    classifications cover only the queries and the overlap block is
    baseline-by-query, ``(B, Q)``. Row tiles and query arithmetic are
    unchanged, so query statistics match a full call exactly.
    """
    from .semantics import classify_annotation_values

    n_rows, n_columns = annotations.shape
    n_query = n_columns-n_baseline
    common = None if 'MAF' not in metadata or metadata.MAF.isna().all() else (metadata.MAF >= common_maf_min).to_numpy()
    counts = np.zeros(n_columns, dtype=np.float64)
    counts_common = None if common is None else np.zeros(n_columns, dtype=np.float64)
    block = np.zeros((n_baseline,n_columns), dtype=np.float64)
    block_common = None if common is None else np.zeros_like(block)
    diagonal = np.zeros(n_query, dtype=np.float64)
    diagonal_common = None if common is None else np.zeros_like(diagonal)
    kinds = dict.fromkeys(annotations.columns, 'binary')
    batches = [(0,n_baseline)] if n_baseline and include_baseline else []
    batches.extend((start,min(start+query_batch_size,n_columns)) for start in range(n_baseline,n_columns,query_batch_size))
    tile = max(1, read_budget_bytes // (8*max(1,2*n_baseline+min(n_query,query_batch_size))))
    for first in range(0,n_rows,tile):
        last = min(first+tile,n_rows)
        rows = slice(first,last)
        baseline_read = annotations.read(rows=rows, columns=annotations.columns[:n_baseline])
        baseline = baseline_read.astype(np.float64)
        selected = None if common is None else common[rows]
        baseline_common = None if selected is None else baseline[selected]
        for start,stop in batches:
            names = annotations.columns[start:stop]
            # The baseline group reuses the tile's read: its float32 values
            # equal the float64 copy's narrowing exactly.
            values = baseline_read.astype(np.float32, copy=False) if start == 0 and stop == n_baseline else annotations.read(rows=rows,columns=names)
            counts[start:stop] += values.sum(axis=0,dtype=np.float64)
            numeric = values.astype(np.float64)
            block[:,start:stop] += baseline.T @ numeric
            if start >= n_baseline:
                diagonal[start-n_baseline:stop-n_baseline] += np.einsum('sq,sq->q',numeric,numeric)
            if selected is not None:
                counts_common[start:stop] += values[selected].sum(axis=0,dtype=np.float64)
                common_values = numeric[selected]
                block_common[:,start:stop] += baseline_common.T @ common_values
                if start >= n_baseline:
                    diagonal_common[start-n_baseline:stop-n_baseline] += np.einsum('sq,sq->q',common_values,common_values)
            if values.dtype == bool:
                continue  # Boolean values are binary by definition.
            for name,kind in classify_annotation_values(pd.DataFrame(values,columns=names)).items():
                if kind == 'quantitative':
                    kinds[name] = kind
    if not include_baseline:
        counts = counts[n_baseline:]
        counts_common = None if counts_common is None else counts_common[n_baseline:]
        block = block[:, n_baseline:]
        block_common = None if block_common is None else block_common[:, n_baseline:]
        kinds = {name: kind for name, kind in kinds.items() if name in annotations.columns[n_baseline:]}
    overlap = OverlapContribution(block,block_common,diagonal,diagonal_common,n_rows,
                                  None if common is None else int(common.sum()))
    return counts.astype(np.float64), None if counts_common is None else counts_common.astype(np.float64), overlap, kinds
