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
    Collapse flexible user input such as ``rsid`` or ``snp_id`` into one
    internal identifier mode.
build_snp_id_series :
    Construct the canonical SNP identifier series for a metadata table.
read_global_snp_restriction :
    Read a restriction list using the active identifier mode.

Design Notes
------------
- Header inference is intentionally permissive to preserve legacy behavior for
  user-provided files with inconsistent capitalization or separators.
- ``chr_pos`` mode always uses normalized chromosome labels and 1-based
  position coordinates encoded as ``CHR:POS`` strings.
"""

from __future__ import annotations

import csv
import gzip
import logging
import re
from pathlib import Path
from typing import Iterable

import pandas as pd

from ..chromosome_inference import normalize_chromosome
from ..column_inference import (
    clean_header,
    infer_chr_bp_columns,
    infer_chr_pos_columns,
    infer_snp_column,
    normalize_genome_build,
    normalize_snp_identifier_mode,
    resolve_restriction_chr_pos_columns,
    resolve_restriction_rsid_column,
)

LOGGER = logging.getLogger("LDSC.identifiers")


def build_chr_pos_snp_id(chrom: object, pos: object, *, context: str | None = None) -> str:
    """Build the canonical ``CHR:POS`` identifier used in ``chr_pos`` mode."""
    chrom_norm = normalize_chromosome(chrom, context=context)
    pos_int = int(pos)
    if pos_int <= 0:
        raise ValueError(f"Position must be positive; got {pos!r}.")
    return f"{chrom_norm}:{pos_int}"


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
    chr_col, pos_col = infer_chr_pos_columns(df.columns)
    return pd.Series(
        [build_chr_pos_snp_id(chrom, pos) for chrom, pos in zip(df[chr_col], df[pos_col])],
        index=df.index,
        name="snp_id",
    )


def validate_unique_snp_ids(df: pd.DataFrame, mode: str, context: str = "table") -> None:
    """Raise ``ValueError`` if a table does not have unique canonical SNP IDs."""
    snp_ids = build_snp_id_series(df, mode)
    duplicated = snp_ids[snp_ids.duplicated()].unique().tolist()
    if duplicated:
        preview = ", ".join(map(str, duplicated[:5]))
        raise ValueError(f"{context} contains non-unique SNP identifiers: {preview}")


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
    if mode == "rsid":
        return _read_rsid_restriction(path)
    return _read_chr_pos_restriction(path, genome_build=genome_build, logger=LOGGER if logger is None else logger)


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
    """Normalize one restriction frame, dropping rows missing CHR or POS."""
    frame = frame.loc[:, ["CHR", "POS"]].copy()
    chr_tokens = frame["CHR"].astype("string")
    chr_missing = chr_tokens.isna() | chr_tokens.str.strip().str.lower().isin({"", "na", "nan", "none"})

    pos_tokens = frame["POS"].astype("string")
    pos_missing = pos_tokens.isna() | pos_tokens.str.strip().str.lower().isin({"", "na", "nan", "none"})
    pos_numeric = pd.to_numeric(frame["POS"], errors="coerce")
    invalid_pos = (~pos_missing) & pos_numeric.isna()
    if invalid_pos.any():
        bad_value = pos_tokens.loc[invalid_pos].iloc[0]
        raise ValueError(f"Restriction POS values in {path} must be numeric; got {bad_value!r}.")

    keep = (~chr_missing) & (~pos_missing)
    dropped = int((~keep).sum())
    if dropped and logger is not None:
        logger.warning(f"Dropped {dropped} restriction rows with missing CHR or POS from '{path}'.")

    normalized = pd.DataFrame(
        {
            "CHR": chr_tokens.loc[keep].astype(str),
            "POS": pos_numeric.loc[keep].astype(int),
        }
    )
    return {
        build_chr_pos_snp_id(chrom, pos, context=str(path))
        for chrom, pos in normalized.itertuples(index=False, name=None)
    }


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
