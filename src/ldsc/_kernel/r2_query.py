"""Point-lookup compute for index-format reference-panel R² parquet files.

Given a per-chromosome ``pyarrow.parquet.ParquetFile`` and canonical ``(i, j)``
sidecar-row index arrays (``i < j``), return the stored adjusted R² and the
panel-orientation sign for each pair. Pairs are matched by an int64 key
``i * n_snps + j`` against decoded row groups; absent pairs yield ``NaN`` R² and
sign ``0``. Two strategies share the match primitive: random-access (row-group
pruning on ``IDX_1`` footer statistics) and streaming (scan every row group).

This module is private kernel code. The user-facing surface is ``ldsc.r2_query``.
"""
from __future__ import annotations

import numpy as np

# Auto-strategy crossover: at or below this many queried pairs, prune+random-access
# row groups; above it, scan every row group sequentially (a large query touches
# most groups anyway, so the linear scan beats per-group seek + pruning overhead).
_RANDOM_STRATEGY_MAX_PAIRS = 50_000


def _arrow_to_numpy(column):
    """Convert an Arrow column to NumPy across PyArrow versions."""
    try:
        return column.to_numpy(zero_copy_only=False)
    except TypeError:
        return column.to_numpy()


def _idx1_column_index(parquet_file) -> int:
    """Return the parquet column index of ``IDX_1`` (for footer statistics)."""
    return parquet_file.schema_arrow.names.index("IDX_1")


def _prune_row_groups(parquet_file, left_sorted: np.ndarray, idx1_col: int):
    """Return row groups whose ``IDX_1`` [min,max] covers a needed left index.

    ``left_sorted`` is the sorted unique set of canonical left indices queried.
    A row group qualifies when at least one queried left index falls within its
    ``IDX_1`` footer range. Groups without usable statistics are kept (safe).
    """
    md = parquet_file.metadata
    keep: list[int] = []
    for rg in range(md.num_row_groups):
        stats = md.row_group(rg).column(idx1_col).statistics
        if stats is None or not stats.has_min_max:
            keep.append(rg)
            continue
        lo, hi = int(stats.min), int(stats.max)
        pos = int(np.searchsorted(left_sorted, lo))
        if pos < left_sorted.shape[0] and int(left_sorted[pos]) <= hi:
            keep.append(rg)
    return keep


def lookup_pairs_in_parquet(
    parquet_file,
    i: np.ndarray,
    j: np.ndarray,
    *,
    n_snps: int,
    r2_scale: float | None,
):
    """Return ``(r2, sign)`` for canonical pairs ``(i, j)`` from one parquet.

    ``r2`` is float32 with ``NaN`` where the pair is not stored; ``sign`` is int8
    with ``+1``/``-1`` for stored pairs (panel orientation) and ``0`` where the
    pair is not stored. When ``r2_scale`` is set the int16 R² column is divided
    by it (dequantization); a float R² column passes through unscaled.
    """
    i = np.asarray(i, dtype=np.int64)
    j = np.asarray(j, dtype=np.int64)
    n = i.shape[0]
    r2_out = np.full(n, np.nan, dtype=np.float32)
    sign_out = np.zeros(n, dtype=np.int8)
    if n == 0:
        return r2_out, sign_out

    n_snps = np.int64(n_snps)
    qkey = i * n_snps + j
    uniq_key = np.unique(qkey)
    uniq_r2 = np.full(uniq_key.shape[0], np.nan, dtype=np.float32)
    uniq_sign = np.zeros(uniq_key.shape[0], dtype=np.int8)

    # Auto strategy: prune+random-access row groups for small queries, full
    # sequential scan once the query is large enough to touch most groups anyway.
    chosen = "random" if n <= _RANDOM_STRATEGY_MAX_PAIRS else "stream"

    idx1_col = _idx1_column_index(parquet_file)
    if chosen == "random":
        row_groups = _prune_row_groups(parquet_file, np.unique(i), idx1_col)
    else:
        row_groups = range(parquet_file.metadata.num_row_groups)

    for rg in row_groups:
        table = parquet_file.read_row_group(int(rg), columns=["IDX_1", "IDX_2", "R2", "SIGN"])
        gi = _arrow_to_numpy(table.column("IDX_1")).astype(np.int64, copy=False)
        gj = _arrow_to_numpy(table.column("IDX_2")).astype(np.int64, copy=False)
        gkey = gi * n_snps + gj
        pos = np.searchsorted(uniq_key, gkey)
        in_range = pos < uniq_key.shape[0]
        matched = np.zeros(gkey.shape[0], dtype=bool)
        matched[in_range] = uniq_key[pos[in_range]] == gkey[in_range]
        if not matched.any():
            continue
        mpos = pos[matched]
        r2_vals = _arrow_to_numpy(table.column("R2")).astype(np.float32, copy=False)[matched]
        if r2_scale is not None:
            r2_vals = r2_vals / np.float32(r2_scale)
        sign_vals = _arrow_to_numpy(table.column("SIGN")).astype(bool, copy=False)[matched]
        uniq_r2[mpos] = r2_vals
        uniq_sign[mpos] = np.where(sign_vals, np.int8(1), np.int8(-1))

    qpos = np.searchsorted(uniq_key, qkey)
    r2_out[:] = uniq_r2[qpos]
    sign_out[:] = uniq_sign[qpos]
    return r2_out, sign_out
