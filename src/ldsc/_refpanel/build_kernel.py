"""Pairwise R2 emission for ``build-r2-matrix``.

``yield_pairwise_r2_rows`` streams standardized genotypes through the exact LD
window and emits the unbiased R2 and correlation sign of every SNP pair inside
it as ``PairColumns`` chunks; ``make_oriented_snp_getter`` re-orients each
genotype column so ``A1`` is the minor allele. Genetic maps are in
``genetic_map.py``, and the PLINK selection and metadata shared with
``ldscore`` are in ``plink_source.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator

import numpy as np

from ..errors import LDSCConfigError, LDSCInternalError
from .r2_format import PairColumns


def make_oriented_snp_getter(base_getter, flip_sign: np.ndarray):
    """Wrap a sequential standardized-SNP getter to apply ``flip_sign`` per column.

    The reference-panel windowed builder reads each retained SNP column exactly
    once, in increasing order, via successive ``base_getter(width)`` calls. This
    wrapper multiplies each returned batch by the matching slice of ``flip_sign``
    (``-1`` negates a standardized column, i.e. recodes it to the other allele),
    tracking the sequential position. Negation flips the correlation sign for
    pairs that include exactly one flipped SNP and leaves R2 (``sign**2``)
    unchanged.
    """
    flip_sign = np.asarray(flip_sign, dtype=np.float32)
    state = {"pos": 0}

    def getter(width):
        columns = base_getter(width)
        ncols = columns.shape[1]
        signs = flip_sign[state["pos"]:state["pos"] + ncols]
        state["pos"] += ncols
        if signs.size:
            columns = columns * signs  # (n, ncols) * (ncols,) broadcast over rows
        return columns

    return getter


def _unbiased_r2_array(correlation: np.ndarray, n_samples: int) -> np.ndarray:
    """Convert correlation coefficients to the unbiased :math:`R^2` estimate.

    Upper-clipped at ``1.0``: float roundoff at perfect LD can push raw
    :math:`r^2` a hair above 1, and the int16 quantization endpoint must map a
    true maximum of exactly ``1.0`` to ``32767``. Negative unbiased values are
    kept (not floored) so the bias correction stays visible downstream.
    """
    sq = correlation * correlation
    denom = n_samples - 2 if n_samples > 2 else n_samples
    return np.minimum(sq - (1.0 - sq) / denom, 1.0)


@dataclass
class _PendingPairs:
    """Per-left-index pair columns held until the left index is final."""

    j: list[np.ndarray]
    r2: list[np.ndarray]
    sign: list[np.ndarray]


def _emit_cross_block_pairs(
    correlation_matrix: np.ndarray,
    a_indices: np.ndarray,
    b_indices: np.ndarray,
    block_left: np.ndarray,
    n_samples: int,
    r2_min: float = 0.0,
) -> PairColumns:
    """Vectorize pair extraction between the carry-over block and the batch."""
    left = block_left[b_indices]
    valid = (a_indices[:, None] >= left[None, :]) & (a_indices[:, None] < b_indices[None, :])
    local_i, local_j = np.nonzero(valid)
    corr = correlation_matrix[local_i, local_j]
    i = a_indices[local_i].astype(np.int64, copy=False)
    j = b_indices[local_j].astype(np.int64, copy=False)
    return _finalize_pair_columns(i, j, corr, n_samples, r2_min)


def _emit_within_block_pairs(
    correlation_matrix: np.ndarray,
    b_indices: np.ndarray,
    block_left: np.ndarray,
    n_samples: int,
    r2_min: float = 0.0,
) -> PairColumns:
    """Vectorize pair extraction within the current SNP batch under block_left."""
    local_i, local_j = np.triu_indices(len(b_indices), k=1)
    global_i = b_indices[local_i]
    global_j = b_indices[local_j]
    valid = global_i >= block_left[global_j]
    local_i, local_j = local_i[valid], local_j[valid]
    corr = correlation_matrix[local_i, local_j]
    i = global_i[valid].astype(np.int64, copy=False)
    j = global_j[valid].astype(np.int64, copy=False)
    return _finalize_pair_columns(i, j, corr, n_samples, r2_min)


def _finalize_pair_columns(
    i: np.ndarray, j: np.ndarray, corr: np.ndarray, n_samples: int, r2_min: float = 0.0
) -> PairColumns:
    """Compute unbiased R2 (float32) and correlation sign (int8) for a batch.

    When ``r2_min > 0`` the batch is filtered to pairs whose unbiased R2 is at
    least ``r2_min``; ``r2_min <= 0`` disables filtering so negative unbiased R2
    pairs are retained (exact backward-compatible completeness).
    """
    r2 = _unbiased_r2_array(corr, n_samples).astype(np.float32, copy=False)
    sign = np.where(corr >= 0, np.int8(1), np.int8(-1))
    if r2_min > 0.0:
        keep = r2 >= r2_min
        i, j, r2, sign = i[keep], j[keep], r2[keep], sign[keep]
    return i, j, r2, sign


def _stash_pair_rows(pending: dict[int, _PendingPairs], columns: PairColumns) -> None:
    """Group an emitted pair batch by left SNP index until that index is final."""
    i, j, r2, sign = columns
    if i.size == 0:
        return
    # Emission yields left indices in order already; sorting is then the identity.
    order = np.arange(i.size) if np.all(i[1:] >= i[:-1]) else np.argsort(i, kind="stable")
    boundaries = np.flatnonzero(np.diff(i[order])) + 1
    for segment in np.split(order, boundaries):
        bucket = pending.setdefault(int(i[segment[0]]), _PendingPairs([], [], []))
        # Right indices are SNP rows of one chromosome, so int32 holds them while pending.
        bucket.j.append(j[segment].astype(np.int32))
        bucket.r2.append(r2[segment])
        bucket.sign.append(sign[segment])


def _pop_pair_rows_before(
    pending: dict[int, _PendingPairs],
    min_future_i: int,
) -> Iterator[PairColumns]:
    """Yield one column chunk per left index whose pairs are final.

    Each chunk is ``(i, j, r2, sign)`` for a single left index with ``j`` sorted
    ascending. Chunks are emitted in increasing left-index order, so concatenated
    ``i`` stays non-decreasing downstream.
    """
    flushable = sorted(i for i in pending if i < int(min_future_i))
    for i in flushable:
        bucket = pending.pop(i)
        j = np.concatenate(bucket.j).astype(np.int64)
        r2 = np.concatenate(bucket.r2)
        sign = np.concatenate(bucket.sign)
        if not np.all(j[1:] >= j[:-1]):
            order = np.argsort(j, kind="stable")
            j, r2, sign = j[order], r2[order], sign[order]
        i_arr = np.full(j.shape, i, dtype=np.int64)
        yield (i_arr, j, r2, sign)


def yield_pairwise_r2_rows(
    *,
    block_left: np.ndarray,
    snp_batch_size: int,
    standardized_snp_getter,
    m: int,
    n: int,
    r2_min: float = 0.0,
) -> Iterator[PairColumns]:
    """
    Yield retained SNP-pair column chunks inside the LD window.

    Each yielded item is one ``PairColumns`` chunk (left index, right index,
    unbiased R2 float32, correlation sign int8) for a single finalized left SNP.

    ``snp_batch_size`` controls the number of SNP genotype columns requested
    from ``standardized_snp_getter`` per computation batch. Window-spanning
    carry-over columns may be retained in memory in addition to the current
    batch.

    ``r2_min`` is an opt-in unbiased-R2 floor: ``r2_min <= 0`` (the default)
    emits every retained pair, while ``r2_min > 0`` drops pairs whose unbiased
    R2 is below the threshold. Filtering trades completeness for memory and
    output size, so callers that need exact completeness must keep the default.
    """

    block_left = np.asarray(block_left, dtype=int)
    if len(block_left) != m:
        raise LDSCInternalError(
            "build-r2-matrix pairwise R2 emission failed: block_left length does not "
            f"match SNP count ({len(block_left)} vs {m}). Most likely LD-window setup "
            "desynchronized from retained SNP metadata. Re-run with DEBUG logging and report the traceback."
        )
    if snp_batch_size <= 0:
        raise LDSCConfigError(
            f"build-r2-matrix received invalid snp_batch_size={snp_batch_size}. Most "
            "likely the batch size was set to zero or a negative value. Pass a positive "
            "integer batch size."
        )

    block_sizes = np.array(np.arange(m) - block_left)
    block_sizes = np.ceil(block_sizes / snp_batch_size) * snp_batch_size

    first_nonzero = np.nonzero(block_left > 0)[0]
    b = int(first_nonzero[0]) if len(first_nonzero) > 0 else m
    b = int(np.ceil(b / snp_batch_size) * snp_batch_size)
    if b > m:
        snp_batch_size = 1
        b = m

    # Flush in increasing-i order; non-decreasing IDX_1 holds only because
    # keep_snps was sorted by source position before reaching here.
    pending_rows: dict[int, _PendingPairs] = {}
    l_A = 0
    A = standardized_snp_getter(b)
    for l_B in range(0, b, snp_batch_size):
        width = min(snp_batch_size, b - l_B)
        yield from _pop_pair_rows_before(pending_rows, int(block_left[l_B]))
        B = A[:, l_B:l_B + width]
        correlation_matrix = np.dot(A.T, B / n)
        _stash_pair_rows(
            pending_rows,
            _emit_cross_block_pairs(
                correlation_matrix=correlation_matrix,
                a_indices=np.arange(l_A, l_A + b),
                b_indices=np.arange(l_B, l_B + width),
                block_left=block_left,
                n_samples=n,
                r2_min=r2_min,
            ),
        )

    b0 = b
    md = int(snp_batch_size * np.floor(m / snp_batch_size))
    end = md + 1 if md != m else md
    previous_chunk_width = snp_batch_size
    for l_B in range(b0, end, snp_batch_size):
        yield from _pop_pair_rows_before(pending_rows, int(block_left[l_B]))
        old_b = b
        b = int(block_sizes[l_B])
        if l_B > b0 and b > 0:
            A = np.hstack((A[:, old_b - b + previous_chunk_width:old_b], B))
            l_A += old_b - b + previous_chunk_width
        elif l_B == b0 and b > 0:
            A = A[:, b0 - b:b0]
            l_A = b0 - b
        elif b == 0:
            A = np.array(()).reshape((n, 0))
            l_A = l_B

        current_chunk_width = snp_batch_size
        if l_B == md:
            current_chunk_width = m - md
        B = standardized_snp_getter(current_chunk_width)
        b_indices = np.arange(l_B, l_B + current_chunk_width)

        if b > 0:
            correlation_matrix = np.dot(A.T, B / n)
            _stash_pair_rows(
                pending_rows,
                _emit_cross_block_pairs(
                    correlation_matrix=correlation_matrix,
                    a_indices=np.arange(l_A, l_A + b),
                    b_indices=b_indices,
                    block_left=block_left,
                    n_samples=n,
                    r2_min=r2_min,
                ),
            )

        within_block = np.dot(B.T, B / n)
        _stash_pair_rows(
            pending_rows,
            _emit_within_block_pairs(
                correlation_matrix=within_block,
                b_indices=b_indices,
                block_left=block_left,
                n_samples=n,
                r2_min=r2_min,
            ),
        )
        previous_chunk_width = current_chunk_width
    yield from _pop_pair_rows_before(pending_rows, m + 1)


