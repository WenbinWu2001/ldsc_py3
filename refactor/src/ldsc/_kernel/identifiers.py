"""
identifiers.py

Core functionality:
    Normalize SNP identifier conventions and load user-supplied restriction
    lists in either ``rsid`` or ``chr_pos`` mode.

Overview
--------
This module keeps all identifier-mode normalization in one place so the
annotation, LD-score, and regression layers can share the same SNP matching
rules. It is intentionally small and stateless: callers pass tables or file
paths in, and receive normalized identifier series or restriction sets back.

Key Functions
-------------
normalize_snp_identifier_mode :
    Collapse flexible user input such as ``rsID`` or ``snp_id`` into one
    internal identifier mode.
build_snp_id_series :
    Construct the canonical SNP identifier series for a metadata table.
read_global_snp_restriction :
    Read a restriction list using the active identifier mode.

Design Notes
------------
- Header inference is intentionally permissive to preserve legacy behavior for
  user-provided files with inconsistent capitalization or separators.
- ``chr_pos`` mode always uses normalized chromosome labels and 1-based base
  pair coordinates encoded as ``CHR:BP`` strings.
"""

from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path
from typing import Iterable

import pandas as pd


SNP_COLUMN_ALIASES = ("SNP", "SNPID", "SNP_ID", "RSID", "RS_ID", "RS", "MARKERNAME", "MARKER")
CHR_COLUMN_ALIASES = ("CHR", "CHROM", "CHROMOSOME")
BP_COLUMN_ALIASES = ("BP", "POS", "POSITION", "BASE_PAIR", "BASEPAIR")


def clean_header(header: str) -> str:
    """Normalize a raw column header into the alias-matching form."""
    return header.upper().replace("-", "_").replace(".", "_").replace("\n", "").strip()


def normalize_snp_identifier_mode(value: str) -> str:
    """
    Normalize flexible user-facing SNP identifier labels into one internal mode.

    Accepted spellings include variants such as ``rsid``, ``rsID``, ``snp``,
    ``snp_id``, and ``chr_pos``. The return value is always one of the internal
    modes used throughout the refactored codebase.
    """
    normalized = re.sub(r"[^a-z0-9]+", "", value.lower())
    if normalized in {"rsid", "rs", "snp", "snpid"}:
        return "rsid"
    if normalized in {"chrpos", "position", "chrompos", "chromosomeposition"}:
        return "chr_pos"
    raise ValueError(f"Unsupported snp_identifier mode: {value!r}")


def normalize_chromosome(value: object) -> str:
    """Return a canonical chromosome label without a leading ``chr`` prefix."""
    text = str(value).strip()
    if not text:
        raise ValueError("Encountered an empty chromosome label.")
    text = re.sub(r"^chr", "", text, flags=re.IGNORECASE)
    if text.isdigit():
        return str(int(text))
    return text.upper()


def build_chr_pos_snp_id(chrom: object, bp: object) -> str:
    """Build the canonical ``CHR:BP`` identifier used in ``chr_pos`` mode."""
    chrom_norm = normalize_chromosome(chrom)
    bp_int = int(bp)
    if bp_int <= 0:
        raise ValueError(f"Base-pair position must be positive; got {bp!r}.")
    return f"{chrom_norm}:{bp_int}"


def build_snp_id_series(df: pd.DataFrame, mode: str) -> pd.Series:
    """
    Build the canonical SNP identifier series for a metadata-like table.

    Parameters
    ----------
    df : pandas.DataFrame
        Table containing either an inferred SNP column or inferred chromosome
        and base-pair columns.
    mode : str
        SNP identifier mode. Flexible user spellings are accepted and
        normalized internally.

    Returns
    -------
    pandas.Series
        One identifier per input row, aligned to ``df.index``.
    """
    mode = normalize_snp_identifier_mode(mode)
    if mode == "rsid":
        snp_col = infer_snp_column(df.columns)
        return df[snp_col].astype(str)
    chr_col, bp_col = infer_chr_bp_columns(df.columns)
    return pd.Series(
        [build_chr_pos_snp_id(chrom, bp) for chrom, bp in zip(df[chr_col], df[bp_col])],
        index=df.index,
        name="snp_id",
    )


def infer_column(header: Iterable[str], aliases: Iterable[str], label: str) -> str:
    """Infer the first header entry matching the supplied alias family."""
    header = list(header)
    alias_set = {clean_header(alias) for alias in aliases}
    for column in header:
        normalized = clean_header(column)
        if normalized in alias_set or any(normalized.endswith(alias) for alias in alias_set):
            return column
    raise ValueError(f"Could not infer a {label} column from: {', '.join(header)}")


def infer_snp_column(header: Iterable[str]) -> str:
    """Infer the SNP identifier column from a table header."""
    return infer_column(header, SNP_COLUMN_ALIASES, "SNP identifier")


def infer_chr_bp_columns(header: Iterable[str]) -> tuple[str, str]:
    """Infer chromosome and base-pair columns from a table header."""
    return (
        infer_column(header, CHR_COLUMN_ALIASES, "chromosome"),
        infer_column(header, BP_COLUMN_ALIASES, "base-pair position"),
    )


def validate_unique_snp_ids(df: pd.DataFrame, mode: str, context: str = "table") -> None:
    """Raise ``ValueError`` if a table does not have unique canonical SNP IDs."""
    snp_ids = build_snp_id_series(df, mode)
    duplicated = snp_ids[snp_ids.duplicated()].unique().tolist()
    if duplicated:
        preview = ", ".join(map(str, duplicated[:5]))
        raise ValueError(f"{context} contains non-unique SNP identifiers: {preview}")


def read_global_snp_restriction(path: str | Path, snp_identifier: str) -> set[str]:
    """
    Read a global SNP restriction file into the active canonical identifier set.

    The returned values are ready for direct set intersection with annotation,
    reference-panel, or regression SNP universes.
    """
    mode = normalize_snp_identifier_mode(snp_identifier)
    path = Path(path)
    if mode == "rsid":
        return _read_rsid_restriction(path)
    return _read_chr_pos_restriction(path)


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
        return [line.rstrip("\n") for line in handle if line.strip() and not line.lstrip().startswith("#")]


def _read_rsid_restriction(path: Path) -> set[str]:
    """Read an rsid-style restriction file into a set of SNP IDs."""
    lines = _non_comment_lines(path)
    if not lines:
        return set()
    delimiter = _detect_delimiter(lines[0])
    if delimiter is None:
        return {line.split()[0] for line in lines}

    reader = csv.reader(lines, delimiter=delimiter)
    rows = list(reader)
    header = rows[0]
    try:
        snp_idx = header.index(infer_snp_column(header))
        data_rows = rows[1:]
    except ValueError:
        snp_idx = 0
        data_rows = rows
    return {row[snp_idx].strip() for row in data_rows if row and row[snp_idx].strip()}


def _read_chr_pos_restriction(path: Path) -> set[str]:
    """Read a ``CHR``/``BP`` restriction file into canonical ``CHR:BP`` IDs."""
    lines = _non_comment_lines(path)
    if not lines:
        return set()
    delimiter = _detect_delimiter(lines[0])
    if delimiter is None:
        values = set()
        for line in lines:
            fields = re.split(r"\s+", line.strip())
            if len(fields) < 2:
                raise ValueError(f"chr_pos restriction rows must contain CHR and BP: {line!r}")
            values.add(build_chr_pos_snp_id(fields[0], fields[1]))
        return values

    reader = csv.reader(lines, delimiter=delimiter)
    rows = list(reader)
    header = rows[0]
    try:
        chr_col, bp_col = infer_chr_bp_columns(header)
        chr_idx = header.index(chr_col)
        bp_idx = header.index(bp_col)
        data_rows = rows[1:]
    except ValueError:
        chr_idx, bp_idx = 0, 1
        data_rows = rows
    values = set()
    for row in data_rows:
        if not row:
            continue
        values.add(build_chr_pos_snp_id(row[chr_idx], row[bp_idx]))
    return values
