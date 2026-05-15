"""
identifiers.py

Core functionality:
    Parse user-supplied SNP restriction lists and delegate active identity
    semantics to ``snp_identity``.

Overview
--------
This module owns low-level restriction-file parsing, delimiter handling, and
header resolution. Exact SNP identity policy lives in
``ldsc._kernel.snp_identity``: four public modes, base/effective keys,
allele-aware cleanup, and restriction key collapse semantics. Callers pass
tables or file paths in, and receive normalized identifier series or
restriction sets back.

Key Functions
-------------
normalize_snp_identifier_mode :
    Validate one exact public SNP identifier mode.
build_snp_id_series :
    Construct the canonical SNP identifier series for a metadata table.
read_global_snp_restriction :
    Read a restriction list using the active identifier mode.
read_global_chr_pos_restriction_key_set :
    Stream a coordinate restriction file into packed uint64 CHR/POS keys for
    memory-sensitive chunk filtering.

Design Notes
------------
- Header inference is intentionally permissive to preserve legacy behavior for
  user-provided files with inconsistent capitalization or separators.
- Public ``snp_identifier`` mode values are exact; header aliases do not create
  mode aliases.
- Public and metadata-facing ``chr_pos`` identifiers remain normalized
  ``CHR:POS`` strings. Memory-sensitive kernels can instead use packed uint64
  CHR/POS keys that carry the same normalized coordinate identity.
"""

from __future__ import annotations

import csv
import gzip
import logging
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .._coordinates import (
    CHR_POS_KEY_COLUMN,
    build_chr_pos_key_frame,
    coordinate_missing_mask,
    normalize_chr_pos_frame,
    positive_int_position_series,
)
from ..chromosome_inference import normalize_chromosome, normalize_chromosome_series
from ..column_inference import (
    A1_COLUMN_SPEC,
    A2_COLUMN_SPEC,
    clean_header,
    infer_chr_bp_columns,
    infer_chr_pos_columns,
    infer_snp_column,
    normalize_genome_build,
    normalize_snp_identifier_mode,
    resolve_optional_column,
    resolve_restriction_chr_pos_columns,
    resolve_restriction_rsid_column,
)
from .snp_identity import (
    RestrictionIdentityKeys,
    collapse_restriction_identity_keys,
    effective_merge_key_series,
    empty_identity_drop_frame,
    identity_mode_family,
    is_allele_aware_mode,
)

LOGGER = logging.getLogger("LDSC.identifiers")
_CHROMOSOME_PACK_CODES = {
    **{str(chrom): chrom for chrom in range(1, 23)},
    "X": 23,
    "Y": 24,
    "M": 25,
    "MT": 26,
}
_MAX_PACKED_POS = (1 << 32) - 1


def build_chr_pos_snp_id(chrom: object, pos: object, *, context: str | None = None) -> str:
    """Build the canonical ``CHR:POS`` identifier used in ``chr_pos`` mode."""
    context = "CHR/POS identifier" if context is None else context
    chrom_norm = normalize_chromosome(chrom, context=context)
    pos_int = int(positive_int_position_series(pd.Series([pos]), context=context).iloc[0])
    return f"{chrom_norm}:{pos_int}"


def build_snp_id_series(df: pd.DataFrame, mode: str) -> pd.Series:
    """
    Build the canonical SNP identifier series for a metadata-like table.

    Parameters
    ----------
    df : pandas.DataFrame
        Table containing either an inferred SNP column or inferred chromosome
        and base-pair columns.
    mode : {"rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"}
        Exact SNP identifier mode. Base modes use only the base key; allele-aware
        modes require usable ``A1/A2`` and build effective keys with the
        normalized allele set.

    Returns
    -------
    pandas.Series
        One identifier per input row, aligned to ``df.index``.
    """
    mode = normalize_snp_identifier_mode(mode)
    try:
        return effective_merge_key_series(df, mode, context="metadata")
    except ValueError:
        if is_allele_aware_mode(mode):
            raise

    if mode == "rsid":
        snp_col = infer_snp_column(df.columns)
        return df[snp_col].astype(str)
    chr_col, pos_col = infer_chr_pos_columns(df.columns)
    return pd.Series(
        [build_chr_pos_snp_id(chrom, pos) for chrom, pos in zip(df[chr_col], df[pos_col])],
        index=df.index,
        name="snp_id",
    )


def build_packed_chr_pos_series(
    chrom: pd.Series,
    pos: pd.Series,
    *,
    context: str,
) -> pd.Series:
    """
    Return compact uint64 CHR/POS keys for memory-sensitive matching.

    Chromosomes are normalized to small integer codes and positions occupy the
    low 32 bits. The output is aligned to the input series index and is intended
    for internal set membership tests where constructing ``CHR:POS`` strings
    would add avoidable memory pressure.
    """
    normalized_chrom = normalize_chromosome_series(chrom, context=context)
    chrom_codes = normalized_chrom.map(_CHROMOSOME_PACK_CODES)
    if bool(chrom_codes.isna().any()):
        bad = normalized_chrom.loc[chrom_codes.isna()].iloc[0]
        raise ValueError(f"Unsupported chromosome label {bad!r} in {context}.")
    pos_values = positive_int_position_series(pos, context=context)
    too_large = pos_values > _MAX_PACKED_POS
    if bool(too_large.any()):
        bad = int(pos_values.loc[too_large].iloc[0])
        raise ValueError(f"POS values in {context} must be <= {_MAX_PACKED_POS}; got {bad}.")

    chrom_array = chrom_codes.to_numpy(dtype=np.uint64)
    pos_array = pos_values.to_numpy(dtype=np.uint64)
    packed = (chrom_array << np.uint64(32)) | pos_array
    return pd.Series(packed, index=chrom.index, dtype="uint64")


def validate_unique_snp_ids(df: pd.DataFrame, mode: str, context: str = "table") -> None:
    """Raise ``ValueError`` if a table does not have unique canonical SNP IDs."""
    mode = normalize_snp_identifier_mode(mode)
    snp_ids = build_snp_id_series(df, mode)
    duplicated = snp_ids[snp_ids.duplicated()].unique().tolist()
    if duplicated:
        preview = ", ".join(map(str, duplicated[:5]))
        details = _duplicate_identifier_details(df, snp_ids, duplicated[:5])
        message = f"{context} contains non-unique SNP identifiers in {mode} mode: {preview}"
        if details:
            message += f". Colliding rows: {details}"
        if mode == "chr_pos":
            message += (
                ". chr_pos mode requires one variant per chromosome/position; "
                "use rsid mode or prune duplicate-coordinate variants before building or running the reference panel."
            )
        raise ValueError(message)


def _duplicate_identifier_details(df: pd.DataFrame, snp_ids: pd.Series, duplicated: list[object]) -> str:
    """Return a compact description of duplicate keys when SNP labels are available."""
    try:
        snp_col = infer_snp_column(df.columns)
    except ValueError:
        return ""
    parts: list[str] = []
    for identifier in duplicated:
        rows = df.loc[snp_ids == identifier, snp_col].astype(str).dropna().unique().tolist()
        if not rows:
            continue
        parts.append(f"{identifier} has SNPs {', '.join(rows[:5])}")
    return "; ".join(parts)


def read_global_snp_restriction(
    path: str | Path,
    snp_identifier: str,
    genome_build: str | None = None,
    logger=None,
) -> set[str]:
    """
    Read a SNP-set definition file into the active canonical identifier set.

    The returned values are ready for direct set intersection with annotation,
    reference-panel, or regression SNP universes, depending on which config
    field supplied the path.
    """
    mode = normalize_snp_identifier_mode(snp_identifier)
    genome_build = normalize_genome_build(genome_build)
    assert genome_build in {"hg19", "hg38", None}, (
        f"genome_build reached kernel as {genome_build!r}; "
        "should have been resolved at workflow entry."
    )
    path = Path(path)
    if identity_mode_family(mode) == "rsid":
        return _read_rsid_restriction(path)
    return _read_chr_pos_restriction(path, genome_build=genome_build, logger=LOGGER if logger is None else logger)


def read_snp_restriction_keys(
    path: str | Path,
    snp_identifier: str,
    genome_build: str | None = None,
    logger=None,
) -> RestrictionIdentityKeys:
    """Read a restriction file into base or allele-aware identity keys."""
    mode = normalize_snp_identifier_mode(snp_identifier)
    genome_build = normalize_genome_build(genome_build)
    assert genome_build in {"hg19", "hg38", None}, (
        f"genome_build reached kernel as {genome_build!r}; "
        "should have been resolved at workflow entry."
    )
    path = Path(path)
    raw_frame = _read_snp_restriction_table(path)
    if len(raw_frame.columns) == 0:
        return RestrictionIdentityKeys(
            keys=set(),
            match_kind="base",
            dropped=empty_identity_drop_frame(),
            n_input_rows=0,
            n_retained_keys=0,
        )
    canonical_frame, has_allele_columns = _canonicalize_snp_restriction_table(
        raw_frame,
        snp_identifier=mode,
        genome_build=genome_build,
        context=str(path),
    )
    return collapse_restriction_identity_keys(
        canonical_frame,
        mode,
        context=f"SNP restriction file {path}",
        has_allele_columns=has_allele_columns,
        logger=LOGGER if logger is None else logger,
    )


def restriction_file_has_allele_columns(path: str | Path, *, context: str | None = None) -> bool:
    """Return whether a restriction file header has a complete A1/A2 pair."""
    path = Path(path)
    header = _read_restriction_header(path)
    if not header:
        return False
    _a1_col, _a2_col, has_allele_columns = _resolve_restriction_allele_columns(
        header,
        context=context or str(path),
    )
    return has_allele_columns


def read_global_chr_pos_restriction_key_set(
    path: str | Path,
    *,
    genome_build: str | None = None,
    logger=None,
) -> set[int]:
    """
    Stream a CHR/POS restriction file as compact packed coordinate keys.

    The reader resolves build-specific position columns in the same way as
    ``read_global_snp_restriction`` but emits uint64 keys instead of public
    ``CHR:POS`` strings. Missing or invalid coordinate rows are dropped with the
    shared coordinate-policy warning and are not included in the returned set.
    """
    genome_build = normalize_genome_build(genome_build)
    assert genome_build in {"hg19", "hg38", None}, (
        f"genome_build reached kernel as {genome_build!r}; "
        "should have been resolved at workflow entry."
    )
    path = Path(path)
    header, row_iter, _delimiter = _iter_restriction_rows(path)
    if not header:
        return set()
    chr_col, pos_col = resolve_restriction_chr_pos_columns(header, genome_build=genome_build, context=str(path))
    chr_idx = header.index(chr_col)
    pos_idx = header.index(pos_col)
    return _read_packed_chr_pos_restriction_rows(
        path=path,
        rows=row_iter,
        chr_idx=chr_idx,
        pos_idx=pos_idx,
        logger=LOGGER if logger is None else logger,
    )


def _open_text(path: Path):
    """Open a plain-text or gzipped restriction file for text iteration."""
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def _detect_delimiter(sample_line: str) -> str | None:
    """Infer a simple delimiter from one non-comment input line."""
    if "\t" in sample_line:
        return "\t"
    if "," in sample_line:
        return ","
    return None


def _non_comment_lines(path: Path) -> list[str]:
    """Read non-empty, non-comment lines from a restriction file."""
    with _open_text(path) as handle:
        lines: list[str] = []
        for line in handle:
            stripped = line.lstrip()
            if not stripped.strip():
                continue
            if stripped.startswith("#") and not stripped.upper().startswith("#CHROM"):
                continue
            lines.append(line.rstrip("\n"))
        return lines


def _iter_non_comment_lines(path: Path):
    """Yield non-empty, non-comment lines from a restriction file."""
    with _open_text(path) as handle:
        for line in handle:
            stripped = line.lstrip()
            if not stripped.strip():
                continue
            if stripped.startswith("#") and not stripped.upper().startswith("#CHROM"):
                continue
            yield line.rstrip("\n")


def _split_restriction_row(line: str, delimiter: str | None) -> list[str]:
    if delimiter is None:
        return re.split(r"\s+", line.strip())
    return next(csv.reader([line], delimiter=delimiter))


def _iter_restriction_rows(path: Path):
    """Return a header plus a streaming row iterator for one restriction file."""
    lines = _iter_non_comment_lines(path)
    try:
        first = next(lines)
    except StopIteration:
        return [], iter(()), None
    delimiter = _detect_delimiter(first)
    header = _split_restriction_row(first, delimiter)
    return header, (_split_restriction_row(line, delimiter) for line in lines), delimiter


def _read_restriction_header(path: Path) -> list[str]:
    """Read only the first non-comment restriction-file row as a header."""
    with _open_text(path) as handle:
        for line in handle:
            stripped = line.lstrip()
            if not stripped.strip():
                continue
            if stripped.startswith("#") and not stripped.upper().startswith("#CHROM"):
                continue
            return _split_restriction_row(line.rstrip("\n"), _detect_delimiter(line))
    return []


def _parse_restriction_rows(path: Path) -> tuple[list[str], list[list[str]], str | None]:
    """Parse one restriction file into a header row plus raw data rows."""
    lines = _non_comment_lines(path)
    if not lines:
        return [], [], None
    delimiter = _detect_delimiter(lines[0])
    if delimiter is None:
        header = re.split(r"\s+", lines[0].strip())
        rows = [re.split(r"\s+", line.strip()) for line in lines[1:]]
        return header, rows, None
    reader = csv.reader(lines, delimiter=delimiter)
    parsed = list(reader)
    return parsed[0], parsed[1:], delimiter


def _read_snp_restriction_table(path: Path) -> pd.DataFrame:
    """Read one restriction file into a raw string table."""
    header, rows, _delimiter = _parse_restriction_rows(path)
    if not header:
        return pd.DataFrame()
    return pd.DataFrame(rows, columns=header, dtype=object)


def _resolve_restriction_allele_columns(columns: Iterable[str], *, context: str | None = None) -> tuple[str | None, str | None, bool]:
    """Resolve optional restriction allele columns, requiring all-or-none."""
    columns = list(columns)
    a1_col = resolve_optional_column(columns, A1_COLUMN_SPEC, context=context)
    a2_col = resolve_optional_column(columns, A2_COLUMN_SPEC, context=context)
    if (a1_col is None) != (a2_col is None):
        present = "A1" if a1_col is not None else "A2"
        missing = "A2" if a1_col is not None else "A1"
        raise ValueError(
            f"Restriction file {context or ''} provides only one allele column ({present}); "
            f"both A1 and A2 allele columns are required for allele-aware restriction matching. "
            f"Missing {missing}."
        )
    return a1_col, a2_col, a1_col is not None and a2_col is not None


def _canonicalize_snp_restriction_table(
    raw_frame: pd.DataFrame,
    *,
    snp_identifier: str,
    genome_build: str | None = None,
    context: str | None = None,
) -> tuple[pd.DataFrame, bool]:
    """Canonicalize a raw restriction table to SNP or CHR/POS plus optional A1/A2."""
    mode = normalize_snp_identifier_mode(snp_identifier)
    columns = list(raw_frame.columns)
    data: dict[str, pd.Series] = {}
    if identity_mode_family(mode) == "rsid":
        snp_col = resolve_restriction_rsid_column(columns, context=context)
        data["SNP"] = raw_frame[snp_col]
    else:
        chr_col, pos_col = resolve_restriction_chr_pos_columns(columns, genome_build=genome_build, context=context)
        data["CHR"] = raw_frame[chr_col]
        data["POS"] = raw_frame[pos_col]

    has_allele_columns = False
    if is_allele_aware_mode(mode):
        a1_col, a2_col, has_allele_columns = _resolve_restriction_allele_columns(columns, context=context)
        if has_allele_columns:
            assert a1_col is not None and a2_col is not None
            data["A1"] = raw_frame[a1_col]
            data["A2"] = raw_frame[a2_col]
    return pd.DataFrame(data), has_allele_columns


def _read_rsid_restriction(path: Path) -> set[str]:
    """Read an rsid-style restriction file into a set of SNP IDs."""
    header, rows, _delimiter = _parse_restriction_rows(path)
    if not header:
        return set()
    snp_idx = header.index(resolve_restriction_rsid_column(header, context=str(path)))
    values: set[str] = set()
    for row in rows:
        if not row:
            continue
        if len(row) <= snp_idx:
            raise ValueError(f"Restriction row in {path} is missing the SNP column.")
        value = row[snp_idx].strip()
        if value:
            values.add(value)
    return values


def _finalize_chr_pos_restriction_frame(
    frame: pd.DataFrame,
    *,
    path: Path,
    logger,
) -> set[str]:
    """Normalize one restriction frame, dropping rows with missing or invalid CHR/POS."""
    frame = frame.loc[:, ["CHR", "POS"]].copy()
    keyed, _report = build_chr_pos_key_frame(
        frame,
        context=f"restriction file '{path}'",
        drop_missing=True,
        logger=logger,
        log_level=logging.WARNING,
    )
    return set(keyed[CHR_POS_KEY_COLUMN].astype(str))


def _read_chr_pos_restriction(path: Path, genome_build: str | None = None, logger=None) -> set[str]:
    """Read a ``CHR``/``POS`` restriction file into canonical ``CHR:POS`` IDs."""
    header, rows, _delimiter = _parse_restriction_rows(path)
    if not header:
        return set()
    chr_col, pos_col = resolve_restriction_chr_pos_columns(header, genome_build=genome_build, context=str(path))
    chr_idx = header.index(chr_col)
    pos_idx = header.index(pos_col)
    values: list[tuple[object, object]] = []
    for row in rows:
        if not row:
            continue
        if len(row) <= max(chr_idx, pos_idx):
            raise ValueError(f"Restriction row in {path} is missing the CHR or POS column.")
        values.append((row[chr_idx], row[pos_idx]))
    frame = pd.DataFrame(values, columns=["CHR", "POS"])
    return _finalize_chr_pos_restriction_frame(frame, path=path, logger=logger)


def _read_packed_chr_pos_restriction_rows(
    *,
    path: Path,
    rows: Iterable[list[str]],
    chr_idx: int,
    pos_idx: int,
    logger,
    chunk_size: int = 100_000,
) -> set[int]:
    """Read restriction rows in bounded chunks and pack retained CHR/POS pairs."""
    values: set[int] = set()
    chrom_chunk: list[object] = []
    pos_chunk: list[object] = []
    missing_dropped = 0
    retained = 0

    def flush() -> None:
        nonlocal chrom_chunk, pos_chunk, missing_dropped, retained
        if not chrom_chunk:
            return
        frame = pd.DataFrame({"CHR": chrom_chunk, "POS": pos_chunk})
        normalized, report = normalize_chr_pos_frame(
            frame,
            context=f"restriction file '{path}'",
            coordinate_policy="drop",
            logger=logger,
            log_level=logging.WARNING,
        )
        missing_dropped += int(report.n_dropped)
        retained += int(report.n_retained)
        if normalized.empty:
            chrom_chunk = []
            pos_chunk = []
            return
        keys = build_packed_chr_pos_series(
            normalized["CHR"].reset_index(drop=True),
            normalized["POS"].reset_index(drop=True),
            context=f"restriction file '{path}'",
        )
        values.update(int(value) for value in keys.to_numpy(dtype=np.uint64))
        chrom_chunk = []
        pos_chunk = []

    for row in rows:
        if not row:
            continue
        if len(row) <= max(chr_idx, pos_idx):
            raise ValueError(f"Restriction row in {path} is missing the CHR or POS column.")
        chrom_chunk.append(row[chr_idx])
        pos_chunk.append(row[pos_idx])
        if len(chrom_chunk) >= chunk_size:
            flush()
    flush()
    return values
