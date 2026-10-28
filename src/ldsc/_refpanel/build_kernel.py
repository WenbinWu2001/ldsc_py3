"""Internal helpers for building R2 matrices."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from os import PathLike
from typing import Iterator, Sequence

import numpy as np
import pandas as pd

from .._core.coordinates import CHR_POS_KEY_COLUMN, build_chr_pos_key_frame
from .._core.chromosomes import chrom_sort_key, normalize_chromosome
from .._core.identifiers import open_text
from .._core.columns import (
    CHR_COLUMN_SPEC,
    POS_COLUMN_SPEC,
    ColumnSpec,
    resolve_required_column,
    resolve_restriction_chr_pos_columns,
    resolve_restriction_rsid_column,
)
from ..errors import LDSCConfigError, LDSCInputError, LDSCInternalError
from .._core.snp_identity import (
    RestrictionIdentityKeys,
    allele_set_series,
    base_key_series,
    identity_mode_family,
    restriction_membership_mask,
)

LOGGER = logging.getLogger("LDSC.r2_matrix_builder.kernel")
_BUILD_REF_PANEL_LIFTOVER_DOC = (
    "docs/user-guide/troubleshooting.md#build-r2-matrix-liftover-or-genetic-map-configuration-is-incomplete"
)


GENETIC_MAP_CM_SPEC = ColumnSpec(
    "CM",
    ("CM", "GENETIC_MAP_CM", "GENETICMAPCM", "GENETIC_MAP(CM)", "GENETICMAP(CM)"),
    "genetic map centiMorgan",
)

def _read_text_table(path: str | PathLike[str]) -> pd.DataFrame:
    """Read a whitespace-delimited plain-text or gzip-compressed table."""
    text_path = str(path)
    compression = "gzip" if text_path.endswith(".gz") else None
    return pd.read_csv(text_path, sep=r"\s+", compression=compression)


def _read_non_comment_lines(path: str | PathLike[str], limit: int = 5) -> list[str]:
    """Return up to ``limit`` non-empty, non-comment lines from ``path``."""
    lines: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            lines.append(stripped)
            if len(lines) >= limit:
                break
    return lines


def load_genetic_map(path: str | PathLike[str]) -> pd.DataFrame:
    """Load one genetic map table into canonical `CHR/POS/CM` columns."""

    df = _read_text_table(path)
    context = str(path)
    chr_col = resolve_required_column(df.columns, CHR_COLUMN_SPEC, context=context)
    pos_col = resolve_required_column(df.columns, POS_COLUMN_SPEC, context=context)
    cm_col = resolve_required_column(df.columns, GENETIC_MAP_CM_SPEC, context=context)

    out = pd.DataFrame(
        {
            "CHR": df[chr_col].map(lambda value: normalize_chromosome(value, context=context)),
            "POS": pd.to_numeric(df[pos_col], errors="raise").astype(np.int64),
            "CM": pd.to_numeric(df[cm_col], errors="raise").astype(float),
        }
    )
    duplicated = out.duplicated(subset=["CHR", "POS"], keep=False)
    if duplicated.any():
        raise LDSCInputError(
            f"build-r2-matrix could not load genetic map '{path}': duplicate CHR/POS "
            "rows were found. Most likely the map file contains repeated positions. "
            "Deduplicate the map or provide a clean genetic map. "
            f"Other causes & fixes: {_BUILD_REF_PANEL_LIFTOVER_DOC}"
        )

    chrom_order = out["CHR"].map(chrom_sort_key)
    if not chrom_order.equals(chrom_order.sort_values(ignore_index=True)):
        raise LDSCInputError(
            f"build-r2-matrix could not load genetic map '{path}': rows are not sorted "
            "by chromosome. Most likely map shards were concatenated out of order. Sort "
            "the file by chromosome and position, or pass correctly ordered map shards."
        )
    for chrom, chrom_frame in out.groupby("CHR", sort=False):
        positions = chrom_frame["POS"].reset_index(drop=True)
        if not positions.equals(positions.sort_values(ignore_index=True)):
            raise LDSCInputError(
                f"build-r2-matrix could not load genetic map '{path}': positions are not "
                f"sorted within chromosome {chrom}. Most likely this chromosome's map rows "
                "were shuffled. Sort the map by chromosome and position before rerunning."
            )

    return out.reset_index(drop=True)


def load_genetic_map_group(paths: Sequence[str | PathLike[str]]) -> pd.DataFrame:
    """Load and combine one or more genetic map files into one canonical table."""

    if not paths:
        raise LDSCInputError(
            "build-r2-matrix could not load a genetic map group because no map paths "
            "were supplied. Most likely an empty path token reached the kernel. Pass at "
            "least one genetic map file or omit cM-window output."
        )
    frames = [load_genetic_map(path) for path in paths]
    combined = pd.concat(frames, axis=0, ignore_index=True)
    combined = combined.sort_values(by=["CHR", "POS"], key=lambda col: col.map(chrom_sort_key) if col.name == "CHR" else col, kind="mergesort")
    duplicated = combined.duplicated(subset=["CHR", "POS"], keep=False)
    if duplicated.any():
        raise LDSCInputError(
            "build-r2-matrix could not combine genetic map inputs: duplicate CHR/POS "
            "rows were found across map files. Most likely overlapping map shards were "
            "supplied. Use non-overlapping shards or a single deduplicated genetic map. "
            f"Other causes & fixes: {_BUILD_REF_PANEL_LIFTOVER_DOC}"
        )
    return combined.reset_index(drop=True)


def interpolate_genetic_map_cm(
    chrom: str,
    positions: np.ndarray,
    genetic_map: pd.DataFrame,
) -> np.ndarray:
    """Interpolate cumulative cM values for `positions` on one chromosome."""

    chrom = normalize_chromosome(chrom)
    chrom_map = genetic_map.loc[genetic_map["CHR"] == chrom, ["POS", "CM"]].reset_index(drop=True)
    if len(chrom_map) == 0:
        raise LDSCInputError(
            f"build-r2-matrix could not interpolate cM positions for chromosome {chrom}: "
            "the genetic map contains no rows for that chromosome. Most likely the wrong "
            "map build or incomplete chromosome shard was supplied. Pass a map that covers "
            f"chromosome {chrom}. Other causes & fixes: {_BUILD_REF_PANEL_LIFTOVER_DOC}"
        )
    return np.interp(
        np.asarray(positions, dtype=np.int64),
        chrom_map["POS"].to_numpy(dtype=np.int64),
        chrom_map["CM"].to_numpy(dtype=float),
    )


def detect_restriction_identifier_mode(path: str | PathLike[str]) -> str:
    """Infer whether a SNP restriction file is keyed by rsID or coordinates."""

    lines = _read_non_comment_lines(path, limit=2)
    if not lines:
        return "rsid"

    header_fields = re.split(r"\t|,|\s+", lines[0])
    try:
        resolve_restriction_chr_pos_columns(header_fields, context=str(path))
        return "chr_pos"
    except (ValueError, LDSCInputError):
        # Not a chr_pos restriction file; try rsid before raising a combined
        # user-facing header error.
        pass
    try:
        resolve_restriction_rsid_column(header_fields, context=str(path))
        return "rsid"
    except (ValueError, LDSCInputError):
        # Neither supported restriction schema matched.
        pass

    raise LDSCInputError(
        f"Cannot detect SNP restriction mode for {path}: no recognizable SNP or CHR/POS header was found. "
        "Most likely the file has no header row or uses unsupported column names. "
        "Add a header with SNP or CHR/POS columns."
    )


def orientation_flip_sign(freq_values: Sequence[float]) -> np.ndarray:
    """Per-SNP sign that re-orients genotypes so ``A1`` becomes the minor allele.

    ``freq_values`` is the frequency of the ``.bim`` second allele (the package's
    ``A2``). A SNP needs flipping when ``A2`` is the minor allele (``freq < 0.5``);
    at ``freq == 0.5`` PLINK order is kept (no flip), so the enforced invariant is
    ``freq(A1) <= 0.5``. Returns ``+1.0`` (keep) / ``-1.0`` (flip) as ``float32``
    so it broadcasts against ``float32`` standardized genotype columns without
    upcasting.
    """
    freq = np.asarray(freq_values, dtype=np.float64)
    return np.where(freq < 0.5, np.float32(-1.0), np.float32(1.0)).astype(np.float32)


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


def build_plink_metadata_frame(
    *,
    bim,
    kept_snps: Sequence[int],
    maf_values: Sequence[float],
    freq_values: Sequence[float],
) -> pd.DataFrame:
    """Build the retained SNP metadata table after PLINK filtering.

    ``A1``/``A2`` are emitted in canonical orientation: ``A1`` is the minor allele
    (``freq(A1) <= 0.5``). ``freq_values`` is the ``.bim`` ``A2`` allele frequency;
    rows with ``freq < 0.5`` have ``A1``/``A2`` swapped. ``maf_values`` is the
    folded ``min(freq, 1-freq)`` and, after the swap, equals ``freq(A1)``.
    """

    kept_snps = np.asarray(kept_snps, dtype=int)
    maf_values = np.asarray(maf_values, dtype=float)
    freq_values = np.asarray(freq_values, dtype=float)
    if not (len(kept_snps) == len(maf_values) == len(freq_values)):
        raise LDSCInternalError(
            "build-r2-matrix metadata assembly failed: kept_snps, maf_values, and "
            f"freq_values have different lengths ({len(kept_snps)} vs "
            f"{len(maf_values)} vs {len(freq_values)}). Most likely PLINK filtering "
            "bookkeeping desynchronized. Re-run with DEBUG logging and report the traceback."
        )
    kept = bim.df.iloc[kept_snps].reset_index(drop=True).copy()
    a1 = kept["A1"].astype(str).to_numpy()
    a2 = kept["A2"].astype(str).to_numpy()
    flip = orientation_flip_sign(freq_values) < 0
    canonical_a1 = np.where(flip, a2, a1)
    canonical_a2 = np.where(flip, a1, a2)
    out = pd.DataFrame(
        {
            "CHR": kept["CHR"].map(normalize_chromosome),
            "SNP": kept["SNP"].astype(str),
            "CM": pd.to_numeric(kept["CM"], errors="coerce"),
            "POS": pd.to_numeric(kept["BP"], errors="raise").astype(np.int64),
            "A1": canonical_a1,
            "A2": canonical_a2,
            "MAF": maf_values.astype(float),
        }
    )
    return out.reset_index(drop=True)


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


# One emitted batch as parallel column arrays: (left index, right index,
# unbiased R2 as float32, sign of the correlation as int8 +1/-1).
PairColumns = tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]


def _emit_cross_block_pairs(
    correlation_matrix: np.ndarray,
    a_indices: np.ndarray,
    b_indices: np.ndarray,
    block_left: np.ndarray,
    n_samples: int,
    min_r2: float = 0.0,
) -> PairColumns:
    """Vectorize pair extraction between the carry-over block and the batch."""
    left = block_left[b_indices]
    valid = (a_indices[:, None] >= left[None, :]) & (a_indices[:, None] < b_indices[None, :])
    local_i, local_j = np.nonzero(valid)
    corr = correlation_matrix[local_i, local_j]
    i = a_indices[local_i].astype(np.int64)
    j = b_indices[local_j].astype(np.int64)
    return _finalize_pair_columns(i, j, corr, n_samples, min_r2)


def _emit_within_block_pairs(
    correlation_matrix: np.ndarray,
    b_indices: np.ndarray,
    block_left: np.ndarray,
    n_samples: int,
    min_r2: float = 0.0,
) -> PairColumns:
    """Vectorize pair extraction within the current SNP batch under block_left."""
    local_i, local_j = np.triu_indices(len(b_indices), k=1)
    global_i = b_indices[local_i]
    global_j = b_indices[local_j]
    valid = global_i >= block_left[global_j]
    local_i, local_j = local_i[valid], local_j[valid]
    corr = correlation_matrix[local_i, local_j]
    i = global_i[valid].astype(np.int64)
    j = global_j[valid].astype(np.int64)
    return _finalize_pair_columns(i, j, corr, n_samples, min_r2)


def _finalize_pair_columns(
    i: np.ndarray, j: np.ndarray, corr: np.ndarray, n_samples: int, min_r2: float = 0.0
) -> PairColumns:
    """Compute unbiased R2 (float32) and correlation sign (int8) for a batch.

    When ``min_r2 > 0`` the batch is filtered to pairs whose unbiased R2 is at
    least ``min_r2``; ``min_r2 <= 0`` disables filtering so negative unbiased R2
    pairs are retained (exact backward-compatible completeness).
    """
    r2 = _unbiased_r2_array(corr, n_samples).astype(np.float32)
    sign = np.where(corr >= 0, np.int8(1), np.int8(-1)).astype(np.int8)
    if min_r2 > 0.0:
        keep = r2 >= min_r2
        i, j, r2, sign = i[keep], j[keep], r2[keep], sign[keep]
    return i, j, r2, sign


def _stash_pair_rows(pending: dict[int, _PendingPairs], columns: PairColumns) -> None:
    """Group an emitted pair batch by left SNP index until that index is final."""
    i, j, r2, sign = columns
    if i.size == 0:
        return
    order = np.argsort(i, kind="stable")
    boundaries = np.flatnonzero(np.diff(i[order])) + 1
    for segment in np.split(order, boundaries):
        bucket = pending.setdefault(int(i[segment[0]]), _PendingPairs([], [], []))
        bucket.j.append(j[segment])
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
        j = np.concatenate(bucket.j)
        r2 = np.concatenate(bucket.r2)
        sign = np.concatenate(bucket.sign)
        order = np.argsort(j, kind="stable")
        j = j[order]
        i_arr = np.full(j.shape, i, dtype=np.int64)
        yield (i_arr, j, r2[order], sign[order])


def yield_pairwise_r2_rows(
    *,
    block_left: np.ndarray,
    snp_batch_size: int,
    standardized_snp_getter,
    m: int,
    n: int,
    min_r2: float = 0.0,
) -> Iterator[PairColumns]:
    """
    Yield retained SNP-pair column chunks inside the LD window.

    Each yielded item is one ``PairColumns`` chunk (left index, right index,
    unbiased R2 float32, correlation sign int8) for a single finalized left SNP.

    ``snp_batch_size`` controls the number of SNP genotype columns requested
    from ``standardized_snp_getter`` per computation batch. Window-spanning
    carry-over columns may be retained in memory in addition to the current
    batch.

    ``min_r2`` is an opt-in unbiased-R2 floor: ``min_r2 <= 0`` (the default)
    emits every retained pair, while ``min_r2 > 0`` drops pairs whose unbiased
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
                min_r2=min_r2,
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
                    min_r2=min_r2,
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
                min_r2=min_r2,
            ),
        )
        previous_chunk_width = current_chunk_width
    yield from _pop_pair_rows_before(pending_rows, m + 1)


def iter_pairwise_r2_rows(
    *,
    block_left: np.ndarray,
    snp_batch_size: int,
    standardized_snp_getter,
    m: int,
    n: int,
    min_r2: float = 0.0,
) -> list[PairColumns]:
    """Materialize :func:`yield_pairwise_r2_rows` as an in-memory list of chunks."""

    return list(
        yield_pairwise_r2_rows(
            block_left=block_left,
            snp_batch_size=snp_batch_size,
            standardized_snp_getter=standardized_snp_getter,
            m=m,
            n=n,
            min_r2=min_r2,
        )
    )


def build_restriction_mask(
    metadata: pd.DataFrame,
    restriction_values: set[str] | RestrictionIdentityKeys,
    mode: str,
) -> np.ndarray:
    """Build the retained-SNP boolean mask for one SNP restriction set."""

    if isinstance(restriction_values, RestrictionIdentityKeys):
        if restriction_values.match_kind == "identity":
            return _build_allele_aware_restriction_mask(
                metadata,
                restriction_values,
                mode,
                context="reference-panel SNP restriction matching",
            ).to_numpy(dtype=bool)
        return restriction_membership_mask(
            metadata,
            restriction_values,
            mode,
            context="reference-panel SNP restriction matching",
        ).to_numpy(dtype=bool)

    if identity_mode_family(mode) == "rsid":
        return metadata["SNP"].astype(str).isin(restriction_values).to_numpy(dtype=bool)
    keyed, _report = build_chr_pos_key_frame(
        metadata.loc[:, ["CHR", "POS"]].copy(),
        context="reference-panel SNP restriction matching",
        drop_missing=True,
        logger=LOGGER,
    )
    keep = pd.Series(False, index=metadata.index)
    keep.loc[keyed.index] = keyed[CHR_POS_KEY_COLUMN].isin(restriction_values).to_numpy(dtype=bool)
    return keep.to_numpy(dtype=bool)


def _build_allele_aware_restriction_mask(
    metadata: pd.DataFrame,
    restriction: RestrictionIdentityKeys,
    mode: str,
    *,
    context: str,
) -> pd.Series:
    """Match allele-aware restrictions without validating unrelated alleles."""
    base = base_key_series(metadata, mode, context=context)
    restricted_bases = _restriction_base_keys(restriction.keys)
    base_match = base.isin(restricted_bases)
    keep = pd.Series(False, index=metadata.index, dtype=bool)
    if not bool(base_match.any()):
        return keep

    allele_set, reasons = allele_set_series(metadata.loc[base_match], context=context)
    candidate_keys = pd.Series(pd.NA, index=allele_set.index, dtype=object)
    valid = base.loc[base_match].notna() & allele_set.notna()
    candidate_keys.loc[valid] = base.loc[base_match].loc[valid].astype(str) + ":" + allele_set.loc[valid].astype(str)
    keep.loc[base_match] = (candidate_keys.isin(restriction.keys) | reasons.notna()).to_numpy(dtype=bool)
    return keep


def _restriction_base_keys(keys: set[str]) -> set[str]:
    """Extract base keys from allele-aware identity keys."""
    base_keys: set[str] = set()
    for key in keys:
        parts = str(key).rsplit(":", 2)
        if len(parts) == 3:
            base_keys.add(parts[0])
    return base_keys
