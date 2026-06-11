"""Annotation overlap-matrix blocks computed during LD-score estimation.

The overlap matrix ``O = AᵀA`` (``A`` the SNP-by-annotation matrix) drives
legacy overlap-aware partitioned-heritability summaries. Because partitioned-h2
fits ``baseline + one query`` per model, only the baseline-rows block
``A_Bᵀ·A`` (shape ``B × (B+Q)``) and each query's self-overlap ``Σ_s A[s,q]²``
are ever needed; cross-query overlaps are not. This module computes those blocks
for the all-SNP and common-SNP universes, with no I/O.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class OverlapContribution:
    """Per-chromosome (or aggregated) overlap blocks for both SNP universes.

    Attributes
    ----------
    baseline_block_all, baseline_block_common : numpy.ndarray or None
        ``A_Bᵀ·A`` over the all/common universe, shape ``(B, B+Q)`` with columns
        ordered ``baseline + query``. The common block is ``None`` when MAF is
        unavailable.
    query_diagonal_all, query_diagonal_common : numpy.ndarray or None
        Query self-overlaps ``Σ_s A[s,q]²``, shape ``(Q,)``.
    n_all, n_common : int or None
        SNP-universe sizes (``M_tot``); ``n_common`` is ``None`` without MAF.
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
    leading baseline columns. The common universe is ``MAF >= common_maf_min``
    (inclusive, matching the common-count mask); it is omitted when MAF metadata
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

    The common universe survives only if every contribution carries it; if any
    chromosome lacks MAF, the aggregate common blocks are ``None`` (mirroring how
    aggregated common counts are dropped when unavailable).
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
