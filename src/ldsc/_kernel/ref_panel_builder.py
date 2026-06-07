"""Internal helpers for building standard parquet reference panels."""

from __future__ import annotations

import gzip
import logging
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from os import PathLike
from typing import Iterable, Iterator, Sequence

import numpy as np
import pandas as pd

from .._coordinates import CHR_POS_KEY_COLUMN, build_chr_pos_key_frame
from ..chromosome_inference import chrom_sort_key, normalize_chromosome
from ..column_inference import (
    CHR_COLUMN_SPEC,
    POS_COLUMN_SPEC,
    ColumnSpec,
    resolve_required_column,
    resolve_restriction_chr_pos_columns,
    resolve_restriction_rsid_column,
)
from ..errors import LDSCConfigError, LDSCDependencyError, LDSCInputError, LDSCInternalError
from .liftover import LiftOverMappingResult, LiftOverTranslator
from .snp_identity import (
    RestrictionIdentityKeys,
    allele_set_series,
    base_key_series,
    identity_artifact_metadata,
    identity_mode_family,
    restriction_membership_mask,
    sidecar_identity_sha256,
)

LOGGER = logging.getLogger("LDSC.ref_panel_builder.kernel")
_BUILD_REF_PANEL_LIFTOVER_DOC = (
    "docs/troubleshooting.md#build-ref-panel-liftover-or-genetic-map-configuration-is-incomplete"
)


GENETIC_MAP_CM_SPEC = ColumnSpec(
    "CM",
    ("CM", "GENETIC_MAP_CM", "GENETICMAPCM", "GENETIC_MAP(CM)", "GENETICMAP(CM)"),
    "genetic map centiMorgan",
)
_INDEX_R2_COLUMNS = ["IDX_1", "IDX_2", "R2", "SIGN"]

def _open_text(path: str | PathLike[str]):
    """Open a plain-text or gzip-compressed text file for reading."""
    text_path = str(path)
    if text_path.endswith(".gz"):
        return gzip.open(text_path, "rt", encoding="utf-8")
    return open(text_path, "rt", encoding="utf-8")


def _read_text_table(path: str | PathLike[str]) -> pd.DataFrame:
    """Read a whitespace-delimited plain-text or gzip-compressed table."""
    text_path = str(path)
    compression = "gzip" if text_path.endswith(".gz") else None
    return pd.read_csv(text_path, sep=r"\s+", compression=compression)


def _read_non_comment_lines(path: str | PathLike[str], limit: int = 5) -> list[str]:
    """Return up to ``limit`` non-empty, non-comment lines from ``path``."""
    lines: list[str] = []
    with _open_text(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            lines.append(stripped)
            if len(lines) >= limit:
                break
    return lines


def _normalize_map_chromosome(value: object) -> str:
    """Normalize one genetic-map chromosome label to the package canon."""
    return normalize_chromosome(value)


def _chrom_sort_key(chrom: object) -> tuple[int, object]:
    """Return the stable chromosome sort key used for map-like inputs."""
    return chrom_sort_key(chrom)


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
            f"build-ref-panel could not load genetic map '{path}': duplicate CHR/POS "
            "rows were found. Most likely the map file contains repeated positions. "
            "Deduplicate the map or provide a clean genetic map. "
            f"Other causes & fixes: {_BUILD_REF_PANEL_LIFTOVER_DOC}"
        )

    chrom_order = out["CHR"].map(_chrom_sort_key)
    if not chrom_order.equals(chrom_order.sort_values(ignore_index=True)):
        raise LDSCInputError(
            f"build-ref-panel could not load genetic map '{path}': rows are not sorted "
            "by chromosome. Most likely map shards were concatenated out of order. Sort "
            "the file by chromosome and position, or pass correctly ordered map shards."
        )
    for chrom, chrom_frame in out.groupby("CHR", sort=False):
        positions = chrom_frame["POS"].reset_index(drop=True)
        if not positions.equals(positions.sort_values(ignore_index=True)):
            raise LDSCInputError(
                f"build-ref-panel could not load genetic map '{path}': positions are not "
                f"sorted within chromosome {chrom}. Most likely this chromosome's map rows "
                "were shuffled. Sort the map by chromosome and position before rerunning."
            )

    return out.reset_index(drop=True)


def load_genetic_map_group(paths: Sequence[str | PathLike[str]]) -> pd.DataFrame:
    """Load and combine one or more genetic map files into one canonical table."""

    if not paths:
        raise LDSCInputError(
            "build-ref-panel could not load a genetic map group because no map paths "
            "were supplied. Most likely an empty path token reached the kernel. Pass at "
            "least one genetic map file or omit cM-window output."
        )
    frames = [load_genetic_map(path) for path in paths]
    combined = pd.concat(frames, axis=0, ignore_index=True)
    combined = combined.sort_values(by=["CHR", "POS"], key=lambda col: col.map(_chrom_sort_key) if col.name == "CHR" else col, kind="mergesort")
    duplicated = combined.duplicated(subset=["CHR", "POS"], keep=False)
    if duplicated.any():
        raise LDSCInputError(
            "build-ref-panel could not combine genetic map inputs: duplicate CHR/POS "
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

    chrom = _normalize_map_chromosome(chrom)
    chrom_map = genetic_map.loc[genetic_map["CHR"] == chrom, ["POS", "CM"]].reset_index(drop=True)
    if len(chrom_map) == 0:
        raise LDSCInputError(
            f"build-ref-panel could not interpolate cM positions for chromosome {chrom}: "
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


def build_plink_metadata_frame(
    *,
    bim,
    kept_snps: Sequence[int],
    maf_values: Sequence[float],
) -> pd.DataFrame:
    """Build the retained SNP metadata table after PLINK filtering."""

    kept_snps = np.asarray(kept_snps, dtype=int)
    maf_values = np.asarray(maf_values, dtype=float)
    if len(kept_snps) != len(maf_values):
        raise LDSCInternalError(
            "build-ref-panel metadata assembly failed: kept_snps and maf_values have "
            f"different lengths ({len(kept_snps)} vs {len(maf_values)}). Most likely "
            "PLINK filtering bookkeeping desynchronized. Re-run with DEBUG logging and report the traceback."
        )
    kept = bim.df.iloc[kept_snps].reset_index(drop=True).copy()
    out = pd.DataFrame(
        {
            "CHR": kept["CHR"].map(_normalize_map_chromosome),
            "SNP": kept["SNP"].astype(str),
            "CM": pd.to_numeric(kept["CM"], errors="coerce"),
            "POS": pd.to_numeric(kept["BP"], errors="raise").astype(np.int64),
            "A1": kept["A1"].astype(str),
            "A2": kept["A2"].astype(str),
            "MAF": maf_values.astype(float),
        }
    )
    return out.reset_index(drop=True)


def build_window_coordinates(
    *,
    metadata: pd.DataFrame,
    cm_values: np.ndarray,
    ld_wind_snps: int | None,
    ld_wind_kb: float | None,
    ld_wind_cm: float | None,
) -> tuple[np.ndarray, float]:
    """Build the coordinate vector used to define block-left windows."""

    if ld_wind_snps is not None:
        return np.arange(len(metadata), dtype=float), float(ld_wind_snps)
    if ld_wind_kb is not None:
        return metadata["POS"].to_numpy(dtype=float), float(ld_wind_kb) * 1000.0
    cm_values = np.asarray(cm_values, dtype=float)
    if np.isnan(cm_values).any():
        raise LDSCInputError(
            "build-ref-panel cannot use cM LD windows because interpolated CM values "
            "are missing for retained SNPs. Most likely the genetic map does not cover "
            "all retained positions. Provide a complete genetic map or use an SNP/kb "
            f"LD window. Other causes & fixes: {_BUILD_REF_PANEL_LIFTOVER_DOC}"
        )
    return cm_values, float(ld_wind_cm)


def compute_block_left(coords: np.ndarray, max_dist: float) -> np.ndarray:
    """Compute the standard LDSC block-left window array."""

    coords = np.asarray(coords, dtype=float)
    m = len(coords)
    block_left = np.zeros(m, dtype=int)
    j = 0
    for i in range(m):
        while j < m and abs(coords[j] - coords[i]) > max_dist:
            j += 1
        block_left[i] = j
    return block_left


# int16 symmetric R2 quantization (refactor #11). Scale 32767 = int16 max, so
# the endpoint is exact: round(1.0 * 32767) = 32767, decoded 32767/32767 = 1.0.
R2_QUANT_SCALE = 32767
R2_ENCODING = "int16_symmetric"


def _quantize_r2(r2: np.ndarray) -> np.ndarray:
    """Symmetric int16 quantization of unbiased R2 (scale ``R2_QUANT_SCALE``).

    ``stored = clip(round(R2 * 32767), -32767, 32767)``. The clip is defensive:
    with the upstream ``min(R2, 1.0)`` and ``|negatives| ~ 1/N`` it never fires
    on real data. float64 multiply makes rounding deterministic at half-steps.
    """
    scaled = np.rint(r2.astype(np.float64) * R2_QUANT_SCALE)
    return np.clip(scaled, -R2_QUANT_SCALE, R2_QUANT_SCALE).astype(np.int16)


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
            "build-ref-panel pairwise R2 emission failed: block_left length does not "
            f"match SNP count ({len(block_left)} vs {m}). Most likely LD-window setup "
            "desynchronized from retained SNP metadata. Re-run with DEBUG logging and report the traceback."
        )
    if snp_batch_size <= 0:
        raise LDSCConfigError(
            f"build-ref-panel received invalid snp_batch_size={snp_batch_size}. Most "
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


def build_runtime_metadata_table(
    *,
    metadata: pd.DataFrame,
    positions: np.ndarray,
    cm_values: np.ndarray | None,
) -> pd.DataFrame:
    """Build the LDSC runtime metadata sidecar for one emitted build.

    When ``cm_values`` is ``None``, the sidecar keeps the SNP rows and writes a
    nullable ``CM`` column. This is used by SNP- and kb-window builds that omit
    genetic maps.
    """
    cm_column = (
        pd.Series(pd.array([pd.NA] * len(metadata), dtype="Float64"))
        if cm_values is None
        else pd.Series(np.asarray(cm_values, dtype=float))
    )

    values = {
        "CHR": metadata["CHR"].map(_normalize_map_chromosome).astype(str),
        "POS": np.asarray(positions, dtype=np.int64),
        "SNP": metadata["SNP"].astype(str),
        "A1": metadata["A1"].astype(str),
        "A2": metadata["A2"].astype(str),
        "CM": cm_column,
        "MAF": pd.to_numeric(metadata["MAF"], errors="coerce").astype(float),
    }
    return pd.DataFrame(values).reset_index(drop=True)


def _ensure_parent_dir(path: str | PathLike[str]) -> None:
    """Create the parent directory for ``path`` if it does not already exist."""
    parent = Path(path).parent
    parent.mkdir(parents=True, exist_ok=True)


def write_dataframe_to_parquet(df: pd.DataFrame, path: str | PathLike[str]) -> str:
    """Write one generic DataFrame to parquet without LDSC format guarantees."""

    _ensure_parent_dir(path)
    if "R2" in df.columns:
        df = df.copy()
        df["R2"] = df["R2"].astype(np.float32)
    try:
        df.to_parquet(path, index=False)
    except ImportError as exc:
        raise LDSCDependencyError(
            "build-ref-panel could not write reference-panel parquet artifacts because no parquet engine is installed. "
            "Most likely the environment is missing pyarrow or fastparquet. Install pyarrow, or install fastparquet "
            "for this legacy parquet writer path."
        ) from exc
    return str(path)


def _standard_r2_index_table(pa, schema, *, i, j, r2, sign):
    """Build one 4-column index ``pyarrow.Table`` from concatenated pair columns.

    ``i``/``j`` are sidecar-row indices (stored int32), ``r2`` is unbiased float32
    stored as symmetric int16 (scale ``R2_QUANT_SCALE``), and ``sign`` is int8
    ``+1``/``-1`` stored as the ``SIGN`` bool (``True`` when the correlation r >= 0).
    The indices are stored as-is; no reference-SNP join or identifier expansion.
    """
    if i.size == 0:
        return schema.empty_table()
    return pa.Table.from_arrays(
        [
            pa.array(i.astype(np.int32), type=pa.int32()),
            pa.array(j.astype(np.int32), type=pa.int32()),
            pa.array(_quantize_r2(r2), type=pa.int16()),
            pa.array(sign == 1, type=pa.bool_()),
        ],
        schema=schema,
    )


def write_r2_parquet(
    *,
    pair_chunks: Iterable[PairColumns],
    path: str | PathLike[str],
    genome_build: str,
    n_samples: int,
    snp_identifier: str,
    n_snps: int,
    sidecar_identity_sha256: str,
    batch_size: int = 100_000,
    row_group_size: int = 50_000,
    min_r2: float = 0.0,
    ld_window_mode: str | None = None,
    ld_window_value: float | None = None,
) -> str:
    """
    Write one canonical 4-column index R2 parquet with LDSC schema metadata.

    Columns are ``IDX_1``/``IDX_2`` (int32 sidecar-row indices), ``R2``
    (int16 symmetric quantization, scale ``R2_QUANT_SCALE``=32767, with
    ``BYTE_STREAM_SPLIT`` encoding; the endpoint ``1.0`` maps to exactly 32767),
    and ``SIGN`` (bit-packed bool, ``True`` when the Pearson correlation r >= 0).
    Pairs must arrive in non-decreasing ``IDX_1`` order because row-group pruning
    depends on monotonic footer statistics.

    The parquet is bound to its per-SNP sidecar by ``ldsc:n_snps`` (the index
    space size) and ``ldsc:sidecar_identity_sha256`` (a hash of the sidecar's
    ``CHR:POS:A1:A2`` identity). The stored indices are meaningless without the
    matching sidecar, so the reader hard-fails on a binding mismatch.

    ``min_r2`` records the unbiased-R2 floor applied upstream during pair
    emission; it is provenance only and does not itself filter rows here. The
    table is pairwise-complete only when ``ldsc:min_r2`` is ``"0.0"`` (or any
    non-positive value), because the read path treats absent pairs as R2=0.

    ``ld_window_mode`` and ``ld_window_value`` record the build-time LD window
    that bounded emitted pairs. LD-score calculation uses these fields to reject
    wider requested windows that would require pairs the parquet does not store.

    Column chunks are compressed with zstd level 9 when the codec is available
    (falling back to snappy otherwise). Compression is recorded per column-chunk
    in the footer and auto-detected on read.
    """

    _ensure_parent_dir(path)
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise LDSCDependencyError(
            "build-ref-panel could not write canonical R2 parquet artifacts because pyarrow is not installed. "
            "Most likely parquet reference-panel generation was requested in an environment missing pyarrow. "
            "Install pyarrow and rerun build-ref-panel."
        ) from exc

    if int(n_snps) >= 2**31:
        raise LDSCInputError(
            f"build-ref-panel cannot write R2 parquet with n_snps={n_snps}: the index "
            "space exceeds the int32 IDX column range. Most likely the retained reference "
            "panel is too large for the current parquet format. Split the panel or reduce "
            "the retained SNP universe."
        )
    if (ld_window_mode is None) != (ld_window_value is None):
        raise LDSCInternalError(
            "build-ref-panel could not write R2 parquet LD-window metadata because "
            "ld_window_mode and ld_window_value were not provided together. Most likely "
            "the reference-panel builder passed an incomplete build-window record. "
            "Re-run with DEBUG logging and report the traceback."
        )

    pa_meta = {
        **{
            f"ldsc:{key}".encode("utf-8"): str(value).encode("utf-8")
            for key, value in identity_artifact_metadata(
                artifact_type="ref_panel_r2",
                snp_identifier=snp_identifier,
                genome_build=genome_build,
            ).items()
        },
        b"ldsc:sorted_by_build": genome_build.encode("utf-8"),
        b"ldsc:row_group_size": str(row_group_size).encode("utf-8"),
        b"ldsc:n_samples": str(n_samples).encode("utf-8"),
        b"ldsc:r2_bias": b"unbiased",
        b"ldsc:min_r2": str(min_r2).encode("utf-8"),
        b"ldsc:n_snps": str(int(n_snps)).encode("utf-8"),
        b"ldsc:sidecar_identity_sha256": str(sidecar_identity_sha256).encode("utf-8"),
        b"ldsc:r2_encoding": R2_ENCODING.encode("utf-8"),
        b"ldsc:r2_scale": str(R2_QUANT_SCALE).encode("utf-8"),
    }
    if ld_window_mode is not None:
        pa_meta[b"ldsc:ld_window_mode"] = str(ld_window_mode).encode("utf-8")
        pa_meta[b"ldsc:ld_window_value"] = str(ld_window_value).encode("utf-8")
    schema = pa.schema(
        [
            ("IDX_1", pa.int32()),
            ("IDX_2", pa.int32()),
            ("R2", pa.int16()),
            ("SIGN", pa.bool_()),
        ]
    ).with_metadata(pa_meta)
    writer = None
    buf_i: list[np.ndarray] = []
    buf_j: list[np.ndarray] = []
    buf_r2: list[np.ndarray] = []
    buf_sign: list[np.ndarray] = []
    buffered = 0
    prev_idx1: int | None = None

    def _flush() -> None:
        nonlocal writer, prev_idx1, buf_i, buf_j, buf_r2, buf_sign, buffered
        i = np.concatenate(buf_i) if buf_i else np.empty(0, dtype=np.int64)
        j = np.concatenate(buf_j) if buf_j else np.empty(0, dtype=np.int64)
        r2 = np.concatenate(buf_r2) if buf_r2 else np.empty(0, dtype=np.float32)
        sign = np.concatenate(buf_sign) if buf_sign else np.empty(0, dtype=np.int8)
        table = _standard_r2_index_table(pa, schema, i=i, j=j, r2=r2, sign=sign)
        if i.size:
            sequence = i if prev_idx1 is None else np.concatenate(([prev_idx1], i))
            bad = np.flatnonzero(np.diff(sequence) < 0)
            if bad.size:
                k = int(bad[0])
                raise LDSCInternalError(
                    "build-ref-panel R2 parquet writing failed: pair rows arrived out of "
                    "non-decreasing IDX_1 order. "
                    f"Received IDX_1={int(sequence[k + 1])} after IDX_1={int(sequence[k])}. "
                    "Most likely retained SNP sorting or pair emission is inconsistent. "
                    "Re-run with DEBUG logging and report the traceback."
                )
            prev_idx1 = int(i[-1])
        if writer is None:
            # IDX_2: DELTA_BINARY_PACKED exploits sorted right-neighbors (median
            # forward-gap = 1), collapsing ~650 MB to ~6 MB vs the default
            # dictionary path. R2: BYTE_STREAM_SPLIT splits the int16 byte planes
            # so zstd compresses the low-entropy high byte (R2 clusters near 0)
            # apart from the noisy low byte. use_dictionary restricts auto-
            # dictionary to IDX_1 only so PyArrow cannot override either encoding.
            enc_kwargs = dict(
                column_encoding={"IDX_2": "DELTA_BINARY_PACKED", "R2": "BYTE_STREAM_SPLIT"},
                use_dictionary=["IDX_1"],
            )
            if pa.Codec.is_available("zstd"):
                writer = pq.ParquetWriter(
                    str(path), schema, compression="zstd", compression_level=9, **enc_kwargs
                )
            else:
                warnings.warn(
                    "zstd codec is not available in this pyarrow build; "
                    "falling back to snappy compression. "
                    "Install pyarrow from conda-forge or PyPI to enable zstd.",
                    UserWarning,
                    stacklevel=2,
                )
                writer = pq.ParquetWriter(str(path), schema, compression="snappy", **enc_kwargs)
        writer.write_table(table, row_group_size=row_group_size)
        buf_i, buf_j, buf_r2, buf_sign = [], [], [], []
        buffered = 0

    try:
        for chunk in pair_chunks:
            ci, cj, cr2, csign = chunk
            if ci.size == 0:
                continue
            buf_i.append(ci)
            buf_j.append(cj)
            buf_r2.append(cr2)
            buf_sign.append(csign)
            buffered += int(ci.size)
            if buffered >= batch_size:
                _flush()
        if buffered or writer is None:
            _flush()
    finally:
        if writer is not None:
            writer.close()
    return str(path)


def write_runtime_metadata_sidecar(
    df: pd.DataFrame,
    path: str | PathLike[str],
    *,
    genome_build: str,
    snp_identifier: str,
) -> str:
    """Write one gzip-compressed LDSC runtime metadata sidecar."""

    _ensure_parent_dir(path)
    identity_metadata = identity_artifact_metadata(
        artifact_type="ref_panel_metadata",
        snp_identifier=snp_identifier,
        genome_build=genome_build,
    )
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for key, value in identity_metadata.items():
            handle.write(f"# ldsc:{key}={value}\n")
        df.to_csv(handle, sep="\t", index=False, na_rep="NA", float_format="%.6g")
    return str(path)

def build_restriction_mask(
    metadata: pd.DataFrame,
    restriction_values: set[str] | RestrictionIdentityKeys,
    mode: str,
) -> np.ndarray:
    """Build the retained-SNP boolean mask for one restriction universe."""

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
    """Extract allele-blind base keys from allele-aware identity keys."""
    base_keys: set[str] = set()
    for key in keys:
        parts = str(key).rsplit(":", 2)
        if len(parts) == 3:
            base_keys.add(parts[0])
    return base_keys
