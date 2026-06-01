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
from ..errors import LDSCDependencyError
from .liftover import LiftOverMappingResult, LiftOverTranslator
from .snp_identity import (
    RestrictionIdentityKeys,
    allele_set_series,
    base_key_series,
    identity_artifact_metadata,
    identity_mode_family,
    restriction_membership_mask,
)

LOGGER = logging.getLogger("LDSC.ref_panel_builder.kernel")


GENETIC_MAP_CM_SPEC = ColumnSpec(
    "CM",
    ("CM", "GENETIC_MAP_CM", "GENETICMAPCM", "GENETIC_MAP(CM)", "GENETICMAP(CM)"),
    "genetic map centiMorgan",
)
_STANDARD_R2_COLUMNS = [
    "CHR",
    "POS_1",
    "POS_2",
    "SNP_1",
    "SNP_2",
    "A1_1",
    "A2_1",
    "A1_2",
    "A2_2",
    "R2",
]


def _empty_standard_r2_table() -> pd.DataFrame:
    """Return an empty canonical R2 table with stable physical dtypes."""
    return pd.DataFrame(
        {
            "CHR": pd.Series(dtype="string"),
            "POS_1": pd.Series(dtype=np.int64),
            "POS_2": pd.Series(dtype=np.int64),
            "SNP_1": pd.Series(dtype="string"),
            "SNP_2": pd.Series(dtype="string"),
            "A1_1": pd.Series(dtype="string"),
            "A2_1": pd.Series(dtype="string"),
            "A1_2": pd.Series(dtype="string"),
            "A2_2": pd.Series(dtype="string"),
            "R2": pd.Series(dtype=np.float32),
        },
        columns=_STANDARD_R2_COLUMNS,
    )

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
        raise ValueError(f"Genetic map {path} contains duplicate chromosome/position rows.")

    chrom_order = out["CHR"].map(_chrom_sort_key)
    if not chrom_order.equals(chrom_order.sort_values(ignore_index=True)):
        raise ValueError(f"Genetic map {path} must be sorted by chromosome.")
    for chrom, chrom_frame in out.groupby("CHR", sort=False):
        positions = chrom_frame["POS"].reset_index(drop=True)
        if not positions.equals(positions.sort_values(ignore_index=True)):
            raise ValueError(f"Genetic map {path} must be sorted by position within chromosome {chrom}.")

    return out.reset_index(drop=True)


def load_genetic_map_group(paths: Sequence[str | PathLike[str]]) -> pd.DataFrame:
    """Load and combine one or more genetic map files into one canonical table."""

    if not paths:
        raise ValueError("At least one genetic map path is required.")
    frames = [load_genetic_map(path) for path in paths]
    combined = pd.concat(frames, axis=0, ignore_index=True)
    combined = combined.sort_values(by=["CHR", "POS"], key=lambda col: col.map(_chrom_sort_key) if col.name == "CHR" else col, kind="mergesort")
    duplicated = combined.duplicated(subset=["CHR", "POS"], keep=False)
    if duplicated.any():
        raise ValueError("Combined genetic map inputs contain duplicate chromosome/position rows.")
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
        raise ValueError(f"Genetic map does not contain chromosome {chrom}.")
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
    except ValueError:
        # Not a chr_pos restriction file; try rsid before raising a combined
        # user-facing header error.
        pass
    try:
        resolve_restriction_rsid_column(header_fields, context=str(path))
        return "rsid"
    except ValueError:
        # Neither supported restriction schema matched.
        pass

    raise ValueError(
        f"Restriction file {path} must contain a header row with recognizable SNP or CHR/POS columns."
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
        raise ValueError("kept_snps and maf_values must have the same length.")
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
        raise ValueError("--ld-wind-cm requires non-missing interpolated CM values.")
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


def _unbiased_r2_array(correlation: np.ndarray, n_samples: int) -> np.ndarray:
    """Convert correlation coefficients to the unbiased :math:`R^2` estimate."""
    sq = correlation * correlation
    denom = n_samples - 2 if n_samples > 2 else n_samples
    return sq - (1.0 - sq) / denom


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
) -> Iterator[dict[str, float | int | str]]:
    """Yield pending rows whose left index cannot appear in future batches."""
    flushable = sorted(i for i in pending if i < int(min_future_i))
    for i in flushable:
        bucket = pending.pop(i)
        j = np.concatenate(bucket.j)
        r2 = np.concatenate(bucket.r2)
        sign = np.concatenate(bucket.sign)
        order = np.argsort(j, kind="stable")
        for k in order:
            yield {
                "i": i,
                "j": int(j[k]),
                "R2": float(r2[k]),
                "sign": "+" if sign[k] == 1 else "-",
            }


def yield_pairwise_r2_rows(
    *,
    block_left: np.ndarray,
    snp_batch_size: int,
    standardized_snp_getter,
    m: int,
    n: int,
    min_r2: float = 0.0,
) -> Iterator[dict[str, float | int | str]]:
    """
    Yield one unordered `R2` row per retained SNP pair inside the LD window.

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
        raise ValueError("block_left length must match the SNP count.")
    if snp_batch_size <= 0:
        raise ValueError("snp_batch_size must be positive.")

    block_sizes = np.array(np.arange(m) - block_left)
    block_sizes = np.ceil(block_sizes / snp_batch_size) * snp_batch_size

    first_nonzero = np.nonzero(block_left > 0)[0]
    b = int(first_nonzero[0]) if len(first_nonzero) > 0 else m
    b = int(np.ceil(b / snp_batch_size) * snp_batch_size)
    if b > m:
        snp_batch_size = 1
        b = m

    # Flush in increasing-i order; non-decreasing POS_1 holds only because
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
) -> list[dict[str, float | int | str]]:
    """Materialize :func:`yield_pairwise_r2_rows` as an in-memory list."""

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


def _build_unique_ids(chromosomes: pd.Series, positions: np.ndarray, a1: pd.Series, a2: pd.Series) -> pd.Series:
    """Build ``CHR:POS:A1:A2`` identifiers for one allele-orientation table."""
    return (
        chromosomes.astype(str)
        + ":"
        + pd.Series(np.asarray(positions, dtype=np.int64), index=chromosomes.index).astype(str)
        + ":"
        + a1.astype(str)
        + ":"
        + a2.astype(str)
    )


def _optional_position_series(positions: np.ndarray | None, index: pd.Index) -> pd.Series:
    """Return concrete positions or nullable missing values for one build."""
    if positions is None:
        return pd.Series(pd.array([pd.NA] * len(index), dtype="Int64"), index=index)
    return pd.Series(np.asarray(positions, dtype=np.int64), index=index)


def _optional_unique_ids(
    chromosomes: pd.Series,
    positions: np.ndarray | None,
    a1: pd.Series,
    a2: pd.Series,
) -> pd.Series:
    """Build unique IDs when positions exist, otherwise return missing strings."""
    if positions is None:
        return pd.Series(pd.array([pd.NA] * len(chromosomes), dtype="string"), index=chromosomes.index)
    return _build_unique_ids(chromosomes, positions, a1, a2)


def build_reference_snp_table(
    *,
    metadata: pd.DataFrame,
    hg19_positions: np.ndarray | None,
    hg38_positions: np.ndarray | None,
) -> pd.DataFrame:
    """Build the in-memory reference SNP table for one chromosome.

    ``hg19_positions`` or ``hg38_positions`` may be ``None`` in source-only
    builds. The returned table keeps both build-specific coordinate slots and
    fills unavailable coordinates and unique-ID fields with missing values.
    """

    chromosomes = metadata["CHR"].map(_normalize_map_chromosome)
    a1 = metadata["A1"].astype(str)
    a2 = metadata["A2"].astype(str)
    table = pd.DataFrame(
        {
            "chr": chromosomes.astype(str),
            "hg19_pos": _optional_position_series(hg19_positions, chromosomes.index),
            "hg38_pos": _optional_position_series(hg38_positions, chromosomes.index),
            "hg19_Uniq_ID": _optional_unique_ids(chromosomes, hg19_positions, a1, a2),
            "hg38_Uniq_ID": _optional_unique_ids(chromosomes, hg38_positions, a1, a2),
            "rsID": metadata["SNP"].astype(str),
            "MAF": pd.to_numeric(metadata["MAF"], errors="coerce").astype(float),
            "A1": a1,
            "A2": a2,
        }
    )
    return table.reset_index(drop=True)


def build_standard_r2_table(
    *,
    pair_rows: list[dict[str, float | int | str]],
    reference_snp_table: pd.DataFrame,
    genome_build: str,
) -> pd.DataFrame:
    """
    Build one canonical R2 table batch for a chromosome.

    The returned frame always uses the package-written parquet schema:
    string-valued ``CHR``/``SNP_1``/``SNP_2``/endpoint allele columns,
    ``int64`` positions, and ``float32`` R2 values. Empty batches keep the same
    dtypes so chromosomes with no emitted pairs still serialize as canonical R2
    parquet files.
    """

    if not pair_rows:
        return _empty_standard_r2_table()

    i = np.asarray([int(row["i"]) for row in pair_rows], dtype=np.int64)
    j = np.asarray([int(row["j"]) for row in pair_rows], dtype=np.int64)
    r2 = np.asarray([float(row["R2"]) for row in pair_rows], dtype=np.float32)
    left = reference_snp_table.iloc[i].reset_index(drop=True)
    right = reference_snp_table.iloc[j].reset_index(drop=True)
    pos_col = f"{genome_build}_pos"
    chr_col = resolve_required_column(left.columns, CHR_COLUMN_SPEC, context="standard R2 reference SNP table")
    return pd.DataFrame(
        {
            "CHR": left[chr_col].astype(str),
            "POS_1": left[pos_col].to_numpy(dtype=np.int64),
            "POS_2": right[pos_col].to_numpy(dtype=np.int64),
            "SNP_1": left["rsID"].astype(str),
            "SNP_2": right["rsID"].astype(str),
            "A1_1": left["A1"].astype(str),
            "A2_1": left["A2"].astype(str),
            "A1_2": right["A1"].astype(str),
            "A2_2": right["A2"].astype(str),
            "R2": r2,
        },
        columns=_STANDARD_R2_COLUMNS,
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
            "Writing reference-panel parquet artifacts requires pyarrow or fastparquet."
        ) from exc
    return str(path)


def _standard_r2_arrow_table(
    pa,
    schema,
    *,
    pair_rows: list[dict[str, float | int | str]],
    reference_snp_table: pd.DataFrame,
    genome_build: str,
):
    """Build one canonical R2 ``pyarrow.Table`` batch directly from numpy.

    Mirrors :func:`build_standard_r2_table` column-for-column but skips the
    pandas intermediate, so large builder batches avoid an extra full-width
    DataFrame copy on the way to the parquet writer. ``schema`` fixes the
    canonical Arrow column types (string/int64/float32) and carries the
    ``ldsc:*`` metadata.
    """
    if not pair_rows:
        return schema.empty_table()
    count = len(pair_rows)
    i = np.fromiter((int(row["i"]) for row in pair_rows), dtype=np.int64, count=count)
    j = np.fromiter((int(row["j"]) for row in pair_rows), dtype=np.int64, count=count)
    r2 = np.fromiter((float(row["R2"]) for row in pair_rows), dtype=np.float32, count=count)
    left = reference_snp_table.iloc[i]
    right = reference_snp_table.iloc[j]
    pos_col = f"{genome_build}_pos"
    chr_col = resolve_required_column(left.columns, CHR_COLUMN_SPEC, context="standard R2 reference SNP table")

    def _strings(series: pd.Series):
        return pa.array(series.astype(str).to_numpy(dtype=object), type=pa.string())

    columns = [
        _strings(left[chr_col]),
        pa.array(left[pos_col].to_numpy(dtype=np.int64), type=pa.int64()),
        pa.array(right[pos_col].to_numpy(dtype=np.int64), type=pa.int64()),
        _strings(left["rsID"]),
        _strings(right["rsID"]),
        _strings(left["A1"]),
        _strings(left["A2"]),
        _strings(right["A1"]),
        _strings(right["A2"]),
        pa.array(r2, type=pa.float32()),
    ]
    return pa.Table.from_arrays(columns, schema=schema)


def write_r2_parquet(
    *,
    pair_rows: Iterable[dict[str, float | int | str]],
    reference_snp_table: pd.DataFrame,
    path: str | PathLike[str],
    genome_build: str,
    n_samples: int,
    snp_identifier: str,
    batch_size: int = 100_000,
    row_group_size: int = 50_000,
    min_r2: float = 0.0,
) -> str:
    """
    Write one canonical R2 parquet table with LDSC schema metadata.

    The writer requires ``pyarrow`` because the canonical format depends on
    Arrow schema metadata and explicit row-group sizing. It writes exactly the
    canonical R2 columns and records ``ldsc:sorted_by_build``,
    ``ldsc:row_group_size``, ``ldsc:n_samples``, ``ldsc:r2_bias``, and
    ``ldsc:min_r2`` in the Arrow schema. Current package-built panels always
    store unbiased R2 values, so ``ldsc:r2_bias`` is written as ``"unbiased"``
    and ``n_samples`` captures the PLINK reader's ``geno.n`` for downstream
    provenance and future raw-R2 compatibility.

    ``min_r2`` records the unbiased-R2 floor applied upstream during pair
    emission; it is provenance only and does not itself filter rows here. The
    table is pairwise-complete only when ``ldsc:min_r2`` is ``"0.0"`` (or any
    non-positive value), because the read path treats absent pairs as R2=0.

    Each batch is assembled as a :class:`pyarrow.Table` directly from numpy
    column arrays via :func:`_standard_r2_arrow_table`, bypassing the pandas
    intermediate to keep peak memory low on deep panels. The on-disk schema is
    identical to the pandas path (string/int64/float32 columns).

    Column chunks are compressed with zstd level 9 when the codec is available
    in the local ``pyarrow`` build (falling back to snappy otherwise, which has
    no level knob). Compression is recorded per column-chunk in the parquet
    footer and auto-detected on read, so it is purely a storage concern and
    never affects decoded values.

    Incoming pair rows must already be sorted by non-decreasing ``POS_1``. The
    writer validates that invariant per batch (vectorized over ``POS_1`` and
    across the batch boundary) because row-group pruning depends on monotonic
    footer statistics.
    """

    _ensure_parent_dir(path)
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise LDSCDependencyError(
            "Writing canonical reference-panel R2 parquet artifacts requires pyarrow."
        ) from exc

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
    }
    schema = pa.schema(
        [
            ("CHR", pa.string()),
            ("POS_1", pa.int64()),
            ("POS_2", pa.int64()),
            ("SNP_1", pa.string()),
            ("SNP_2", pa.string()),
            ("A1_1", pa.string()),
            ("A2_1", pa.string()),
            ("A1_2", pa.string()),
            ("A2_2", pa.string()),
            ("R2", pa.float32()),
        ]
    ).with_metadata(pa_meta)
    writer = None
    batch: list[dict[str, float | int | str]] = []
    prev_pos1: int | None = None

    def _flush(rows: list[dict[str, float | int | str]]) -> None:
        nonlocal writer, prev_pos1
        table = _standard_r2_arrow_table(
            pa,
            schema,
            pair_rows=rows,
            reference_snp_table=reference_snp_table,
            genome_build=genome_build,
        )
        pos1 = table.column("POS_1").to_numpy()
        if pos1.size:
            sequence = pos1 if prev_pos1 is None else np.concatenate(([prev_pos1], pos1))
            bad = np.flatnonzero(np.diff(sequence) < 0)
            if bad.size:
                k = int(bad[0])
                raise ValueError(
                    "Pairs must arrive in non-decreasing POS_1 order. "
                    f"Received POS_1={int(sequence[k + 1])} after POS_1={int(sequence[k])}. "
                    "Verify that the reference panel builder traverses SNPs in ascending positional order."
                )
            prev_pos1 = int(pos1[-1])
        if writer is None:
            if pa.Codec.is_available("zstd"):
                # zstd level 9: ~3.5% smaller than the default level 1 with
                # identical read speed and memory (snappy has no level knob).
                writer = pq.ParquetWriter(str(path), schema, compression="zstd", compression_level=9)
            else:
                warnings.warn(
                    "zstd codec is not available in this pyarrow build; "
                    "falling back to snappy compression. "
                    "Install pyarrow from conda-forge or PyPI to enable zstd.",
                    UserWarning,
                    stacklevel=2,
                )
                writer = pq.ParquetWriter(str(path), schema, compression="snappy")
        writer.write_table(table, row_group_size=row_group_size)

    try:
        for row in pair_rows:
            batch.append(row)
            if len(batch) >= batch_size:
                _flush(batch)
                batch = []
        if batch or writer is None:
            _flush(batch)
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
