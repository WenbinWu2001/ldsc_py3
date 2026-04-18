"""Shared column-alias inference for workflow and kernel table parsing.

This module centralizes canonical field names, accepted aliases, and the
matching logic used to resolve user-provided headers and parquet schema names.
Callers should infer columns here and then rename to canonical names locally.
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
import re
import warnings
from typing import Iterable, Sequence


LOGGER = logging.getLogger("LDSC.columns")
_WARNED_INFERENCES: set[tuple[str, str, str]] = set()


def normalize_column_token(value: str) -> str:
    """Normalize a column-like token for alias matching."""
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


@dataclass(frozen=True)
class ColumnSpec:
    """Canonical field description used for alias-based inference."""

    canonical: str
    aliases: tuple[str, ...]
    label: str
    allow_suffix_match: bool = True


SNP_COLUMN_ALIASES = ("SNP", "SNPID", "SNP_ID", "RSID", "RS_ID", "RS", "MARKERNAME", "MARKER")
CHR_COLUMN_ALIASES = ("CHR", "CHROM", "CHROMOSOME")
POS_COLUMN_ALIASES = ("POS", "BP", "POSITION", "BASE_PAIR", "BASEPAIR")
CM_COLUMN_ALIASES = ("CM", "CMBP", "CENTIMORGAN")
MAF_COLUMN_ALIASES = ("MAF", "FRQ", "FREQ", "FREQUENCY")

CHR_COLUMN_SPEC = ColumnSpec("CHR", CHR_COLUMN_ALIASES, "chromosome")
POS_COLUMN_SPEC = ColumnSpec("POS", POS_COLUMN_ALIASES, "position")
SNP_COLUMN_SPEC = ColumnSpec("SNP", SNP_COLUMN_ALIASES, "SNP identifier")
CM_COLUMN_SPEC = ColumnSpec("CM", CM_COLUMN_ALIASES, "centiMorgan")
MAF_COLUMN_SPEC = ColumnSpec("MAF", MAF_COLUMN_ALIASES, "minor-allele frequency")

_HG19_BUILD_ALIASES = ("hg19", "hg37", "grch37")
_HG38_BUILD_ALIASES = ("hg38", "grch38")
_POS_TOKEN_ALIASES = ("pos", "bp", "position", "base_pair", "basepair")
_UNIQ_ID_TOKEN_ALIASES = ("uniq_id", "uniqid", "unique_id", "uniqueid")


def _indexed_aliases(base_aliases: Sequence[str], index: int) -> tuple[str, ...]:
    values: list[str] = []
    for alias in base_aliases:
        values.extend((f"{alias}_{index}", f"{alias}{index}"))
    return tuple(values)


def _build_position_aliases(build_aliases: Sequence[str], index: int) -> tuple[str, ...]:
    values: list[str] = []
    for build in build_aliases:
        for pos in _POS_TOKEN_ALIASES:
            values.extend((f"{build}_{pos}_{index}", f"{build}_{pos}{index}"))
    return tuple(values)


def _build_uniq_id_aliases(build_aliases: Sequence[str], index: int) -> tuple[str, ...]:
    values: list[str] = []
    for build in build_aliases:
        for token in _UNIQ_ID_TOKEN_ALIASES:
            values.extend((f"{build}_{token}_{index}", f"{build}_{token}{index}"))
    return tuple(values)


R2_SOURCE_COLUMN_SPECS = (
    ColumnSpec("chr", CHR_COLUMN_ALIASES, "R2 chromosome"),
    ColumnSpec("rsID_1", _indexed_aliases(SNP_COLUMN_ALIASES, 1), "R2 rsID_1"),
    ColumnSpec("rsID_2", _indexed_aliases(SNP_COLUMN_ALIASES, 2), "R2 rsID_2"),
    ColumnSpec("hg38_pos_1", _build_position_aliases(_HG38_BUILD_ALIASES, 1), "R2 hg38 position 1"),
    ColumnSpec("hg38_pos_2", _build_position_aliases(_HG38_BUILD_ALIASES, 2), "R2 hg38 position 2"),
    ColumnSpec("hg19_pos_1", _build_position_aliases(_HG19_BUILD_ALIASES, 1), "R2 hg19 position 1"),
    ColumnSpec("hg19_pos_2", _build_position_aliases(_HG19_BUILD_ALIASES, 2), "R2 hg19 position 2"),
    ColumnSpec("hg38_Uniq_ID_1", _build_uniq_id_aliases(_HG38_BUILD_ALIASES, 1), "R2 hg38 unique ID 1"),
    ColumnSpec("hg38_Uniq_ID_2", _build_uniq_id_aliases(_HG38_BUILD_ALIASES, 2), "R2 hg38 unique ID 2"),
    ColumnSpec("hg19_Uniq_ID_1", _build_uniq_id_aliases(_HG19_BUILD_ALIASES, 1), "R2 hg19 unique ID 1"),
    ColumnSpec("hg19_Uniq_ID_2", _build_uniq_id_aliases(_HG19_BUILD_ALIASES, 2), "R2 hg19 unique ID 2"),
    ColumnSpec("R2", ("R2",), "R2 value"),
    ColumnSpec("Dprime", ("DPRIME", "D_PRIME"), "D-prime"),
    ColumnSpec("+/-corr", ("+/-CORR", "CORR", "SIGNCORR", "SIGN_CORR"), "signed correlation"),
)
R2_SOURCE_COLUMN_SPEC_MAP = {spec.canonical: spec for spec in R2_SOURCE_COLUMN_SPECS}

R2_HELPER_COLUMN_SPECS = (
    ColumnSpec(
        "pos_1",
        ("POS_1", "POS1", "BP_1", "BP1", "POSITION_1", "POSITION1"),
        "normalized R2 position 1",
        allow_suffix_match=False,
    ),
    ColumnSpec(
        "pos_2",
        ("POS_2", "POS2", "BP_2", "BP2", "POSITION_2", "POSITION2"),
        "normalized R2 position 2",
        allow_suffix_match=False,
    ),
)
R2_HELPER_COLUMN_SPEC_MAP = {spec.canonical: spec for spec in R2_HELPER_COLUMN_SPECS}


def _best_matches(columns: Iterable[str], spec: ColumnSpec) -> list[str]:
    columns = list(columns)
    canonical_norm = normalize_column_token(spec.canonical)
    alias_norms = {normalize_column_token(alias) for alias in spec.aliases}
    ranked: dict[int, list[str]] = {}
    for column in columns:
        normalized = normalize_column_token(column)
        rank = None
        if normalized == canonical_norm:
            rank = 0
        elif normalized in alias_norms:
            rank = 1
        elif spec.allow_suffix_match and normalized.endswith(canonical_norm):
            rank = 2
        elif spec.allow_suffix_match and any(normalized.endswith(alias) for alias in alias_norms):
            rank = 3
        if rank is not None:
            ranked.setdefault(rank, []).append(column)
    if not ranked:
        return []
    return ranked[min(ranked)]


def _warn_if_inferred(canonical: str, actual: str, context: str | None) -> None:
    if actual == canonical:
        return
    key = (context or "", canonical, actual)
    if key in _WARNED_INFERENCES:
        return
    _WARNED_INFERENCES.add(key)
    suffix = "" if not context else f" in {context}"
    message = f"Inferred canonical field '{canonical}' from input column '{actual}'{suffix}."
    LOGGER.info(message)
    warnings.warn(message, UserWarning, stacklevel=3)


def resolve_required_column(
    columns: Iterable[str],
    spec: ColumnSpec,
    *,
    context: str | None = None,
) -> str:
    """Resolve one required column from ``columns`` using ``spec``."""
    matches = _best_matches(columns, spec)
    if not matches:
        header = ", ".join(str(column) for column in columns)
        raise ValueError(f"Could not infer a {spec.label} column for canonical field {spec.canonical!r} from: {header}")
    if len(matches) > 1:
        raise ValueError(
            f"Ambiguous {spec.label} columns for canonical field {spec.canonical!r}: "
            + ", ".join(matches)
        )
    actual = matches[0]
    _warn_if_inferred(spec.canonical, actual, context)
    return actual


def resolve_optional_column(
    columns: Iterable[str],
    spec: ColumnSpec,
    *,
    context: str | None = None,
) -> str | None:
    """Resolve one optional column from ``columns`` using ``spec``."""
    matches = _best_matches(columns, spec)
    if not matches:
        return None
    if len(matches) > 1:
        raise ValueError(
            f"Ambiguous {spec.label} columns for canonical field {spec.canonical!r}: "
            + ", ".join(matches)
        )
    actual = matches[0]
    _warn_if_inferred(spec.canonical, actual, context)
    return actual


def resolve_required_columns(
    columns: Iterable[str],
    specs: Sequence[ColumnSpec],
    *,
    context: str | None = None,
) -> dict[str, str]:
    """Resolve many required canonical fields from one header/schema."""
    return {spec.canonical: resolve_required_column(columns, spec, context=context) for spec in specs}
