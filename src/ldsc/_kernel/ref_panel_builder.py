"""Internal helpers for building standard parquet reference panels."""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from os import PathLike
from typing import Iterable, Iterator, Sequence

import numpy as np
import pandas as pd

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
from .identifiers import build_snp_id_series


GENETIC_MAP_CM_SPEC = ColumnSpec(
    "CM",
    ("CM", "GENETIC_MAP_CM", "GENETICMAPCM", "GENETIC_MAP(CM)", "GENETICMAP(CM)"),
    "genetic map centiMorgan",
)
_STANDARD_R2_COLUMNS = [
    "CHR",
    "POS_1",
    "POS_2",
    "R2",
    "SNP_1",
    "SNP_2",
]


def _empty_standard_r2_table() -> pd.DataFrame:
    """Return an empty canonical R2 table with stable physical dtypes."""
    return pd.DataFrame(
        {
            "CHR": pd.Series(dtype="string"),
            "POS_1": pd.Series(dtype=np.int64),
            "POS_2": pd.Series(dtype=np.int64),
            "R2": pd.Series(dtype=np.float32),
            "SNP_1": pd.Series(dtype="string"),
            "SNP_2": pd.Series(dtype="string"),
        },
        columns=_STANDARD_R2_COLUMNS,
    )


@dataclass(frozen=True)
class LiftOverMappingResult:
    """Partial liftover result for one chromosome."""

    translated_positions: np.ndarray
    keep_mask: np.ndarray
    unmapped_count: int
    cross_chrom_count: int


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


def _try_normalize_map_chromosome(value: object) -> str | None:
    """Normalize a liftover hit chromosome, ignoring unsupported auxiliary contigs."""
    try:
        return _normalize_map_chromosome(value)
    except ValueError:
        return None


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
    """Infer whether a SNP restriction file is keyed by `rsid` or `chr_pos`."""

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


def _unbiased_r2_from_correlation(correlation: float, n_samples: int) -> float:
    """Convert a correlation coefficient to the unbiased :math:`R^2` estimate."""
    sq = correlation * correlation
    denom = n_samples - 2 if n_samples > 2 else n_samples
    return sq - (1.0 - sq) / denom


def _emit_cross_block_pairs(
    correlation_matrix: np.ndarray,
    a_indices: np.ndarray,
    b_indices: np.ndarray,
    block_left: np.ndarray,
    n_samples: int,
) -> list[dict[str, float | int | str]]:
    """Build pair rows between the carry-over block and the current chunk."""
    rows: list[dict[str, float | int | str]] = []
    for local_j, global_j in enumerate(b_indices):
        valid = (a_indices >= int(block_left[global_j])) & (a_indices < global_j)
        for local_i in np.flatnonzero(valid):
            corr = float(correlation_matrix[local_i, local_j])
            rows.append(
                {
                    "i": int(a_indices[local_i]),
                    "j": int(global_j),
                    "R2": _unbiased_r2_from_correlation(corr, n_samples),
                    "sign": "+" if corr >= 0 else "-",
                }
            )
    return rows


def _emit_within_block_pairs(
    correlation_matrix: np.ndarray,
    b_indices: np.ndarray,
    block_left: np.ndarray,
    n_samples: int,
) -> list[dict[str, float | int | str]]:
    """Build pair rows within the current chunk subject to ``block_left``."""
    rows: list[dict[str, float | int | str]] = []
    for local_j, global_j in enumerate(b_indices):
        for local_i in range(local_j):
            global_i = int(b_indices[local_i])
            if global_i < int(block_left[global_j]):
                continue
            corr = float(correlation_matrix[local_i, local_j])
            rows.append(
                {
                    "i": global_i,
                    "j": int(global_j),
                    "R2": _unbiased_r2_from_correlation(corr, n_samples),
                    "sign": "+" if corr >= 0 else "-",
                }
            )
    return rows


def _stash_pair_rows(
    pending: dict[int, list[dict[str, float | int | str]]],
    rows: Iterable[dict[str, float | int | str]],
) -> None:
    """Group emitted pair rows by left SNP index until that index is final."""
    for row in rows:
        pending.setdefault(int(row["i"]), []).append(row)


def _pop_pair_rows_before(
    pending: dict[int, list[dict[str, float | int | str]]],
    min_future_i: int,
) -> Iterator[dict[str, float | int | str]]:
    """Yield pending rows whose left index cannot appear in future chunks."""
    flushable = sorted(i for i in pending if i < int(min_future_i))
    for i in flushable:
        rows = pending.pop(i)
        rows.sort(key=lambda row: int(row["j"]))
        yield from rows


def yield_pairwise_r2_rows(
    *,
    block_left: np.ndarray,
    chunk_size: int,
    standardized_snp_getter,
    m: int,
    n: int,
) -> Iterator[dict[str, float | int | str]]:
    """Yield one unordered `R2` row per retained SNP pair inside the LD window."""

    block_left = np.asarray(block_left, dtype=int)
    if len(block_left) != m:
        raise ValueError("block_left length must match the SNP count.")
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive.")

    block_sizes = np.array(np.arange(m) - block_left)
    block_sizes = np.ceil(block_sizes / chunk_size) * chunk_size

    first_nonzero = np.nonzero(block_left > 0)[0]
    b = int(first_nonzero[0]) if len(first_nonzero) > 0 else m
    b = int(np.ceil(b / chunk_size) * chunk_size)
    if b > m:
        chunk_size = 1
        b = m

    # Flush in increasing-i order; non-decreasing POS_1 holds only because
    # keep_snps was sorted by source position before reaching here.
    pending_rows: dict[int, list[dict[str, float | int | str]]] = {}
    l_A = 0
    A = standardized_snp_getter(b)
    for l_B in range(0, b, chunk_size):
        width = min(chunk_size, b - l_B)
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
            ),
        )

    b0 = b
    md = int(chunk_size * np.floor(m / chunk_size))
    end = md + 1 if md != m else md
    previous_chunk_width = chunk_size
    for l_B in range(b0, end, chunk_size):
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

        current_chunk_width = chunk_size
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
            ),
        )
        previous_chunk_width = current_chunk_width
    yield from _pop_pair_rows_before(pending_rows, m + 1)


def iter_pairwise_r2_rows(
    *,
    block_left: np.ndarray,
    chunk_size: int,
    standardized_snp_getter,
    m: int,
    n: int,
) -> list[dict[str, float | int | str]]:
    """Materialize :func:`yield_pairwise_r2_rows` as an in-memory list."""

    return list(
        yield_pairwise_r2_rows(
            block_left=block_left,
            chunk_size=chunk_size,
            standardized_snp_getter=standardized_snp_getter,
            m=m,
            n=n,
        )
    )


def _build_unique_ids(chromosomes: pd.Series, positions: np.ndarray, ref: pd.Series, alt: pd.Series) -> pd.Series:
    """Build ``CHR:POS:REF:ALT`` identifiers for one allele-orientation table."""
    return (
        chromosomes.astype(str)
        + ":"
        + pd.Series(np.asarray(positions, dtype=np.int64), index=chromosomes.index).astype(str)
        + ":"
        + ref.astype(str)
        + ":"
        + alt.astype(str)
    )


def _optional_position_series(positions: np.ndarray | None, index: pd.Index) -> pd.Series:
    """Return concrete positions or nullable missing values for one build."""
    if positions is None:
        return pd.Series(pd.array([pd.NA] * len(index), dtype="Int64"), index=index)
    return pd.Series(np.asarray(positions, dtype=np.int64), index=index)


def _optional_unique_ids(
    chromosomes: pd.Series,
    positions: np.ndarray | None,
    ref: pd.Series,
    alt: pd.Series,
) -> pd.Series:
    """Build unique IDs when positions exist, otherwise return missing strings."""
    if positions is None:
        return pd.Series(pd.array([pd.NA] * len(chromosomes), dtype="string"), index=chromosomes.index)
    return _build_unique_ids(chromosomes, positions, ref, alt)


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
    ref = metadata["A1"].astype(str)
    alt = metadata["A2"].astype(str)
    table = pd.DataFrame(
        {
            "chr": chromosomes.astype(str),
            "hg19_pos": _optional_position_series(hg19_positions, chromosomes.index),
            "hg38_pos": _optional_position_series(hg38_positions, chromosomes.index),
            "hg19_Uniq_ID": _optional_unique_ids(chromosomes, hg19_positions, ref, alt),
            "hg38_Uniq_ID": _optional_unique_ids(chromosomes, hg38_positions, ref, alt),
            "rsID": metadata["SNP"].astype(str),
            "MAF": pd.to_numeric(metadata["MAF"], errors="coerce").astype(float),
            "REF": ref,
            "ALT": alt,
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
    Build one canonical six-column R2 table chunk for a chromosome.

    The returned frame always uses the package-written parquet schema:
    string-valued ``CHR``/``SNP_1``/``SNP_2``, ``int64`` positions, and
    ``float32`` R2 values. Empty chunks keep the same dtypes so chromosomes with
    no emitted pairs still serialize as canonical R2 parquet files.
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
            "R2": r2,
            "SNP_1": left["rsID"].astype(str),
            "SNP_2": right["rsID"].astype(str),
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

    return pd.DataFrame(
        {
            "CHR": metadata["CHR"].map(_normalize_map_chromosome).astype(str),
            "POS": np.asarray(positions, dtype=np.int64),
            "SNP": metadata["SNP"].astype(str),
            "CM": cm_column,
            "MAF": pd.to_numeric(metadata["MAF"], errors="coerce").astype(float),
        }
    ).reset_index(drop=True)


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


def write_r2_parquet(
    *,
    pair_rows: Iterable[dict[str, float | int | str]],
    reference_snp_table: pd.DataFrame,
    path: str | PathLike[str],
    genome_build: str,
    n_samples: int,
    batch_size: int = 100_000,
    row_group_size: int = 50_000,
) -> str:
    """
    Write a canonical R2 parquet table with row-group metadata.

    The writer requires ``pyarrow`` because the canonical format depends on
    Arrow schema metadata and explicit row-group sizing. It writes exactly the
    six canonical R2 columns, stores ``ldsc:sorted_by_build`` and
    ``ldsc:row_group_size`` in schema metadata, and validates that incoming
    pair rows are sorted by non-decreasing ``POS_1``.
    """

    _ensure_parent_dir(path)
    pos_col = f"{genome_build}_pos"
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise LDSCDependencyError(
            "Writing canonical reference-panel R2 parquet artifacts requires pyarrow."
        ) from exc

    pa_meta = {
        b"ldsc:sorted_by_build": genome_build.encode("utf-8"),
        b"ldsc:row_group_size": str(row_group_size).encode("utf-8"),
        b"ldsc:n_samples": str(n_samples).encode("utf-8"),
        b"ldsc:r2_bias": b"unbiased",
    }
    writer = None
    batch: list[dict[str, float | int | str]] = []
    prev_pos1: int | None = None
    try:
        for row in pair_rows:
            current_pos1 = int(reference_snp_table.iloc[int(row["i"])][pos_col])
            if prev_pos1 is not None and current_pos1 < prev_pos1:
                raise ValueError(
                    "Pairs must arrive in non-decreasing POS_1 order. "
                    f"Received POS_1={current_pos1} after POS_1={prev_pos1}. "
                    "Verify that the reference panel builder traverses SNPs in ascending positional order."
                )
            prev_pos1 = current_pos1
            batch.append(row)
            if len(batch) < batch_size:
                continue
            frame = build_standard_r2_table(
                pair_rows=batch,
                reference_snp_table=reference_snp_table,
                genome_build=genome_build,
            )
            table = pa.Table.from_pandas(frame, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(
                    str(path),
                    table.schema.with_metadata({**(table.schema.metadata or {}), **pa_meta}),
                )
            writer.write_table(table, row_group_size=row_group_size)
            batch = []

        if batch or writer is None:
            frame = build_standard_r2_table(
                pair_rows=batch,
                reference_snp_table=reference_snp_table,
                genome_build=genome_build,
            )
            table = pa.Table.from_pandas(frame, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(
                    str(path),
                    table.schema.with_metadata({**(table.schema.metadata or {}), **pa_meta}),
                )
            writer.write_table(table, row_group_size=row_group_size)
    finally:
        if writer is not None:
            writer.close()
    return str(path)


def write_runtime_metadata_sidecar(df: pd.DataFrame, path: str | PathLike[str]) -> str:
    """Write one gzip-compressed LDSC runtime metadata sidecar."""

    _ensure_parent_dir(path)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        df.to_csv(handle, sep="\t", index=False, na_rep="NA", float_format="%.6g")
    return str(path)


class LiftOverTranslator:
    """Translate one chromosome's 1-based positions between explicit genome-build chain files."""

    def __init__(
        self,
        *,
        source_build: str,
        target_build: str,
        chain_path: str | PathLike[str] | None = None,
    ) -> None:
        """Initialize a build-to-build position translator for one run."""
        self.source_build = source_build
        self.target_build = target_build
        self.chain_path = None if chain_path is None else Path(chain_path)
        if source_build == target_build:
            self._identity = True
            self._liftover = None
            return
        self._identity = False
        if self.chain_path is None:
            flag = (
                "--liftover-chain-hg19-to-hg38"
                if (source_build, target_build) == ("hg19", "hg38")
                else "--liftover-chain-hg38-to-hg19"
            )
            raise ValueError(
                f"An explicit liftover chain path is required for {source_build} -> {target_build}. "
                f"Provide it via {flag}."
            )
        try:
            from pyliftover import LiftOver
        except ImportError as exc:
            raise LDSCDependencyError(
                "Building reference panels across hg19/hg38 requires the optional dependency 'pyliftover'."
            ) from exc
        self.chain_path = Path(self.chain_path)
        self._liftover = LiftOver(str(self.chain_path))

    def map_positions(self, chrom: str, positions: Sequence[int] | np.ndarray) -> LiftOverMappingResult:
        """Translate 1-based positions and retain only same-chromosome mappings."""

        array = np.asarray(positions, dtype=np.int64)
        if self._identity:
            return LiftOverMappingResult(
                translated_positions=array.copy(),
                keep_mask=np.ones(len(array), dtype=bool),
                unmapped_count=0,
                cross_chrom_count=0,
            )
        chrom = _normalize_map_chromosome(chrom)
        query_chrom = f"chr{chrom}"
        translated = np.zeros(len(array), dtype=np.int64)
        keep_mask = np.zeros(len(array), dtype=bool)
        unmapped_count = 0
        cross_chrom_count = 0
        for idx, pos in enumerate(array):
            hits = self._liftover.convert_coordinate(query_chrom, int(pos) - 1)
            if not hits:
                unmapped_count += 1
                continue
            hit = next(
                (candidate for candidate in hits if _try_normalize_map_chromosome(candidate[0]) == chrom),
                None,
            )
            if hit is None:
                cross_chrom_count += 1
                continue
            translated[idx] = int(hit[1]) + 1
            keep_mask[idx] = True
        return LiftOverMappingResult(
            translated_positions=translated[keep_mask],
            keep_mask=keep_mask,
            unmapped_count=unmapped_count,
            cross_chrom_count=cross_chrom_count,
        )

    def translate_positions(self, chrom: str, positions: Sequence[int] | np.ndarray) -> np.ndarray:
        """Translate 1-based positions and require complete same-chromosome mapping."""

        result = self.map_positions(chrom, positions)
        if not result.keep_mask.all():
            failed = np.asarray(positions, dtype=np.int64)[~result.keep_mask]
            preview = ", ".join(str(int(value)) for value in failed[:5])
            raise ValueError(
                f"Failed to liftover {len(failed)} positions on chromosome {chrom} "
                f"from {self.source_build} to {self.target_build}: {preview}"
            )
        return result.translated_positions


def build_restriction_mask(metadata: pd.DataFrame, restriction_values: set[str], mode: str) -> np.ndarray:
    """Build the retained-SNP boolean mask for one restriction universe."""

    if mode == "rsid":
        return metadata["SNP"].astype(str).isin(restriction_values).to_numpy(dtype=bool)
    keys = build_snp_id_series(metadata.loc[:, ["CHR", "POS"]].copy(), "chr_pos")
    return keys.isin(restriction_values).to_numpy(dtype=bool)
