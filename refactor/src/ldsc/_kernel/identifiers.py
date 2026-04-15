"""SNP identifier normalization and restriction-list helpers."""

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
    return header.upper().replace("-", "_").replace(".", "_").replace("\n", "").strip()


def normalize_snp_identifier_mode(value: str) -> str:
    normalized = re.sub(r"[^a-z0-9]+", "", value.lower())
    if normalized in {"rsid", "rs", "snp", "snpid"}:
        return "rsid"
    if normalized in {"chrpos", "position", "chrompos", "chromosomeposition"}:
        return "chr_pos"
    raise ValueError(f"Unsupported snp_identifier mode: {value!r}")


def normalize_chromosome(value: object) -> str:
    text = str(value).strip()
    if not text:
        raise ValueError("Encountered an empty chromosome label.")
    text = re.sub(r"^chr", "", text, flags=re.IGNORECASE)
    if text.isdigit():
        return str(int(text))
    return text.upper()


def build_chr_pos_snp_id(chrom: object, bp: object) -> str:
    chrom_norm = normalize_chromosome(chrom)
    bp_int = int(bp)
    if bp_int <= 0:
        raise ValueError(f"Base-pair position must be positive; got {bp!r}.")
    return f"{chrom_norm}:{bp_int}"


def build_snp_id_series(df: pd.DataFrame, mode: str) -> pd.Series:
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
    header = list(header)
    alias_set = {clean_header(alias) for alias in aliases}
    for column in header:
        normalized = clean_header(column)
        if normalized in alias_set or any(normalized.endswith(alias) for alias in alias_set):
            return column
    raise ValueError(f"Could not infer a {label} column from: {', '.join(header)}")


def infer_snp_column(header: Iterable[str]) -> str:
    return infer_column(header, SNP_COLUMN_ALIASES, "SNP identifier")


def infer_chr_bp_columns(header: Iterable[str]) -> tuple[str, str]:
    return (
        infer_column(header, CHR_COLUMN_ALIASES, "chromosome"),
        infer_column(header, BP_COLUMN_ALIASES, "base-pair position"),
    )


def validate_unique_snp_ids(df: pd.DataFrame, mode: str, context: str = "table") -> None:
    snp_ids = build_snp_id_series(df, mode)
    duplicated = snp_ids[snp_ids.duplicated()].unique().tolist()
    if duplicated:
        preview = ", ".join(map(str, duplicated[:5]))
        raise ValueError(f"{context} contains non-unique SNP identifiers: {preview}")


def read_global_snp_restriction(path: str | Path, snp_identifier: str) -> set[str]:
    mode = normalize_snp_identifier_mode(snp_identifier)
    path = Path(path)
    if mode == "rsid":
        return _read_rsid_restriction(path)
    return _read_chr_pos_restriction(path)


def _open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def _detect_delimiter(sample_line: str) -> str | None:
    if "\t" in sample_line:
        return "\t"
    if "," in sample_line:
        return ","
    return None


def _non_comment_lines(path: Path) -> list[str]:
    with _open_text(path) as handle:
        return [line.rstrip("\n") for line in handle if line.strip() and not line.lstrip().startswith("#")]


def _read_rsid_restriction(path: Path) -> set[str]:
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
