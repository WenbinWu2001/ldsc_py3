# The clean_header helper and munging aliases derive from LDSC munge_sumstats.py.
# (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
# GPLv3 (LICENSE). Python 3 adaptation by Anthony Abrantes, 2024. Modified 2026:
# centralized aliases and added typed identifier and column normalization.
# See NOTICE for attribution and dated modifications.

"""Central naming and alias inference registry for LDSC3.

This module is the single source of truth for:

- column-header normalization
- column-alias families scoped by input context
- SNP identifier modes
- genome build aliases

Callers should infer columns here and then rename to canonical field names
locally. Internal package-written artifacts use canonical-only families so they
fail fast when their schema drifts.
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
import re
from typing import Iterable, Sequence

from ..errors import LDSCConfigError, LDSCInputError, LDSCInternalError

LOGGER = logging.getLogger("LDSC.columns")
_LOGGED_INFERENCES: set[tuple[str, str, str]] = set()


def clean_header(header: str) -> str:
    """Normalize a raw header using the legacy ``munge_sumstats`` rules."""
    return str(header).upper().replace("-", "_").replace(".", "_").replace("\n", "").strip()


def normalize_column_token(value: str) -> str:
    """Normalize a column-like token for permissive alias matching."""
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


_HG19_BUILD_ALIASES = ("hg19", "hg37", "grch37")
_HG38_BUILD_ALIASES = ("hg38", "grch38")
_AUTO_BUILD_ALIASES = ("auto",)


def normalize_genome_build(genome_build: str | None) -> str | None:
    """Normalize a flexible genome-build label to ``hg19``, ``hg38``, or ``auto``."""
    if genome_build is None:
        return None
    normalized = normalize_column_token(genome_build)
    if normalized in {normalize_column_token(alias) for alias in _AUTO_BUILD_ALIASES}:
        return "auto"
    if normalized in {normalize_column_token(alias) for alias in _HG19_BUILD_ALIASES}:
        return "hg19"
    if normalized in {normalize_column_token(alias) for alias in _HG38_BUILD_ALIASES}:
        return "hg38"
    raise LDSCConfigError(
        f"Could not normalize genome_build={genome_build!r}. "
        "Most likely the genome build label is misspelled or unsupported. "
        "Use None, 'auto', 'hg19', 'hg37', 'GRCh37', 'hg38', or 'GRCh38'."
    )


@dataclass(frozen=True)
class ColumnSpec:
    """Canonical field description used for alias-based inference."""

    canonical: str
    aliases: tuple[str, ...]
    label: str
    allow_suffix_match: bool = True


def _canonical_only_spec(canonical: str, label: str) -> ColumnSpec:
    """Return a strict spec that accepts only the canonical field name."""
    return ColumnSpec(canonical, (canonical,), label, allow_suffix_match=False)


_POS_TOKEN_ALIASES = ("pos", "bp", "position", "base_pair", "basepair")


# Shared external column families.
SNP_COLUMN_ALIASES = ("SNP", "SNPID", "SNP_ID", "RSID", "RS_ID", "RS", "ID", "MARKERNAME", "MARKER", "MARKERID")
CHR_COLUMN_ALIASES = ("CHR", "#CHROM", "CHROM", "CHROMOSOME")
POS_COLUMN_ALIASES = ("POS", "BP", "POSITION", "BASE_PAIR", "BASEPAIR")
CM_COLUMN_ALIASES = ("CM", "CMBP", "CENTIMORGAN")
MAF_COLUMN_ALIASES = ("MAF",)

CHR_COLUMN_SPEC = ColumnSpec("CHR", CHR_COLUMN_ALIASES, "chromosome")
POS_COLUMN_SPEC = ColumnSpec("POS", POS_COLUMN_ALIASES, "position")
SNP_COLUMN_SPEC = ColumnSpec("SNP", SNP_COLUMN_ALIASES, "SNP identifier")
CM_COLUMN_SPEC = ColumnSpec("CM", CM_COLUMN_ALIASES, "centiMorgan")
MAF_COLUMN_SPEC = ColumnSpec("MAF", MAF_COLUMN_ALIASES, "folded minor-allele frequency")
A1_COLUMN_SPEC = ColumnSpec(
    "A1",
    ("A1", "ALLELE1", "ALLELE_1", "EFFECT_ALLELE", "REFERENCE_ALLELE", "REF", "EA"),
    "allele 1",
)
A2_COLUMN_SPEC = ColumnSpec(
    "A2",
    ("A2", "ALLELE2", "ALLELE_2", "OTHER_ALLELE", "NON_EFFECT_ALLELE", "ALT", "NEA"),
    "allele 2",
)
RESTRICTION_ALLELE_SPECS = (A1_COLUMN_SPEC, A2_COLUMN_SPEC)


# Raw summary-stat alias families preserve legacy behavior and add current
# non-conflicting aliases where semantics match.
RAW_SUMSTATS_REQUIRED_OR_OPTIONAL_SPECS = (
    ColumnSpec(
        "SNP",
        ("SNP", "MARKERNAME", "MARKERID", "SNPID", "SNP_ID", "RS", "RSID", "RS_ID", "ID", "RS_NUMBER", "RS_NUMBERS", "MARKER"),
        "summary-stat SNP identifier",
    ),
    CHR_COLUMN_SPEC,
    POS_COLUMN_SPEC,
    ColumnSpec("NSTUDY", ("NSTUDY", "N_STUDY", "NSTUDIES", "N_STUDIES"), "summary-stat number of studies"),
    ColumnSpec("P", ("P", "PVALUE", "P_VALUE", "PVAL", "P_VAL", "GC_PVALUE"), "summary-stat p-value"),
    ColumnSpec(
        "A1",
        ("A1", "ALLELE1", "ALLELE_1", "EFFECT_ALLELE", "REFERENCE_ALLELE", "INC_ALLELE", "EA"),
        "signed-statistic allele",
    ),
    ColumnSpec(
        "A2",
        ("A2", "ALLELE2", "ALLELE_2", "OTHER_ALLELE", "NON_EFFECT_ALLELE", "DEC_ALLELE", "NEA"),
        "counterpart allele",
    ),
    ColumnSpec("N", ("N", "WEIGHT"), "summary-stat sample size"),
    ColumnSpec("N_CAS", ("NCASE", "CASES_N", "N_CASE", "N_CASES", "N_CAS", "NCAS", "Nca"), "summary-stat case count"),
    ColumnSpec("N_CON", ("N_CONTROLS", "N_CON", "NCONTROL", "CONTROLS_N", "N_CONTROL", "N_CTRL", "NCON", "Nco"), "summary-stat control count"),
    ColumnSpec("INFO", ("INFO", "IMPINFO"), "summary-stat INFO"),
    ColumnSpec("FRQ", ("EAF", "FRQ", "MAF", "FRQ_U", "F_U"), "summary-stat allele frequency"),
)

RAW_SUMSTATS_SIGNED_STAT_SPECS = (
    ColumnSpec("Z", ("ZSCORE", "Z-SCORE", "GC_ZSCORE", "Z"), "summary-stat Z-score"),
    ColumnSpec("OR", ("OR",), "summary-stat odds ratio"),
    ColumnSpec("BETA", ("B", "BETA", "EFFECTS", "EFFECT"), "summary-stat beta"),
    ColumnSpec("LOG_ODDS", ("LOG_ODDS",), "summary-stat log odds"),
    ColumnSpec("SIGNED_SUMSTAT", ("SIGNED_SUMSTAT",), "explicit signed summary statistic"),
)


# Context families for permissive external inputs.
ANNOTATION_METADATA_SPECS = (CHR_COLUMN_SPEC, POS_COLUMN_SPEC, SNP_COLUMN_SPEC, CM_COLUMN_SPEC, MAF_COLUMN_SPEC)
ANNOTATION_METADATA_SPEC_MAP = {spec.canonical: spec for spec in ANNOTATION_METADATA_SPECS}

REFERENCE_METADATA_SPECS = (CHR_COLUMN_SPEC, POS_COLUMN_SPEC, SNP_COLUMN_SPEC, CM_COLUMN_SPEC, MAF_COLUMN_SPEC)
REFERENCE_METADATA_SPEC_MAP = {spec.canonical: spec for spec in REFERENCE_METADATA_SPECS}

RESTRICTION_RSID_SPECS = (SNP_COLUMN_SPEC,)
RESTRICTION_RSID_SPEC_MAP = {spec.canonical: spec for spec in RESTRICTION_RSID_SPECS}

RESTRICTION_CHRPOS_SPECS = (CHR_COLUMN_SPEC, POS_COLUMN_SPEC)
RESTRICTION_CHRPOS_SPEC_MAP = {spec.canonical: spec for spec in RESTRICTION_CHRPOS_SPECS}
RESTRICTION_HG19_POS_SPEC = ColumnSpec("hg19_POS", tuple(f"{build}_{pos}" for build in _HG19_BUILD_ALIASES for pos in _POS_TOKEN_ALIASES), "hg19 position")
RESTRICTION_HG38_POS_SPEC = ColumnSpec("hg38_POS", tuple(f"{build}_{pos}" for build in _HG38_BUILD_ALIASES for pos in _POS_TOKEN_ALIASES), "hg38 position")


# Strict internal artifact families.
INTERNAL_SUMSTATS_ARTIFACT_SPECS = (
    _canonical_only_spec("CHR", "munged summary-statistics chromosome"),
    _canonical_only_spec("POS", "munged summary-statistics position"),
    _canonical_only_spec("SNP", "munged summary-statistics SNP identifier"),
    _canonical_only_spec("N", "munged summary-statistics sample size"),
    _canonical_only_spec("Z", "munged summary-statistics Z-score"),
    _canonical_only_spec("A1", "munged summary-statistics allele 1"),
    _canonical_only_spec("A2", "munged summary-statistics allele 2"),
    _canonical_only_spec("FRQ", "munged summary-statistics allele frequency"),
)
INTERNAL_SUMSTATS_ARTIFACT_SPEC_MAP = {
    spec.canonical: spec for spec in INTERNAL_SUMSTATS_ARTIFACT_SPECS
}


INTERNAL_LDSCORE_ARTIFACT_SPECS = (
    _canonical_only_spec("CHR", "internal LD-score chromosome"),
    _canonical_only_spec("POS", "internal LD-score position"),
    _canonical_only_spec("SNP", "internal LD-score SNP identifier"),
    _canonical_only_spec("CM", "internal LD-score centiMorgan"),
    _canonical_only_spec("MAF", "internal LD-score minor-allele frequency"),
)


def build_cleaned_alias_lookup(specs: Sequence[ColumnSpec]) -> dict[str, str]:
    """Build a legacy-style cleaned-header-to-canonical alias map."""
    alias_map: dict[str, str] = {}
    for spec in specs:
        values = (spec.canonical,) + tuple(spec.aliases)
        for value in values:
            key = clean_header(value)
            existing = alias_map.get(key)
            if existing is not None and existing != spec.canonical:
                raise LDSCInternalError(
                    f"Cannot build column-alias registry: alias {value!r} maps to both {existing!r} "
                    f"and {spec.canonical!r}. Most likely a new column alias collides with an existing field. "
                    "Remove the duplicate alias or split the conflicting spec."
                )
            alias_map[key] = spec.canonical
    return alias_map


def _best_matches(columns: Iterable[str], spec: ColumnSpec) -> list[str]:
    """Return the best-ranked header matches for one canonical field spec."""
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


def _log_if_inferred(canonical: str, actual: str, context: str | None) -> None:
    """Emit a one-time DEBUG record when alias inference renames a source column."""
    if actual == canonical:
        return
    key = (context or "", canonical, actual)
    if key in _LOGGED_INFERENCES:
        return
    _LOGGED_INFERENCES.add(key)
    suffix = "" if not context else f" in {context}"
    message = f"Inferred canonical field '{canonical}' from input column '{actual}'{suffix}."
    LOGGER.debug(message)


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
        context_text = f" in {context}" if context else ""
        raise LDSCInputError(
            f"Could not map the required {spec.label} column{context_text}: expected canonical field {spec.canonical!r}; "
            f"available columns: {header}. Most likely the header uses an unrecognized alias or the file delimiter is wrong. "
            "Pass the matching column explicitly when the workflow supports it, or rename the input header."
        )
    if len(matches) > 1:
        context_text = f" in {context}" if context else ""
        raise LDSCInputError(
            f"Could not map the {spec.label} column{context_text}: multiple columns match canonical field {spec.canonical!r}: "
            f"{', '.join(matches)}. Most likely the file contains duplicate aliases for the same concept. "
            "Rename or remove the extra matching columns so exactly one remains."
        )
    actual = matches[0]
    _log_if_inferred(spec.canonical, actual, context)
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
        context_text = f" in {context}" if context else ""
        raise LDSCInputError(
            f"Could not map the optional {spec.label} column{context_text}: multiple columns match canonical field {spec.canonical!r}: "
            f"{', '.join(matches)}. Most likely the file contains duplicate aliases for the same optional field. "
            "Rename or remove the extra matching columns so at most one remains."
        )
    actual = matches[0]
    _log_if_inferred(spec.canonical, actual, context)
    return actual


def resolve_required_columns(
    columns: Iterable[str],
    specs: Sequence[ColumnSpec],
    *,
    context: str | None = None,
) -> dict[str, str]:
    """Resolve many required canonical fields from one header/schema."""
    return {spec.canonical: resolve_required_column(columns, spec, context=context) for spec in specs}


def infer_snp_column(header: Iterable[str], *, context: str | None = None) -> str:
    """Infer the SNP identifier column from a restriction-like header."""
    return resolve_required_column(header, RESTRICTION_RSID_SPEC_MAP["SNP"], context=context)


def infer_chr_pos_columns(header: Iterable[str], *, context: str | None = None) -> tuple[str, str]:
    """Infer chromosome and position columns from a generic CHR/POS header."""
    return (
        resolve_required_column(header, RESTRICTION_CHRPOS_SPEC_MAP["CHR"], context=context),
        resolve_required_column(header, RESTRICTION_CHRPOS_SPEC_MAP["POS"], context=context),
    )


def _restriction_header_error(columns: Iterable[str], missing: Sequence[str]) -> LDSCInputError:
    """Build the shared error when a restriction header lacks required aliases."""
    header = ", ".join(str(column) for column in columns)
    joined = ", ".join(missing)
    return LDSCInputError(
        f"Could not read SNP restriction file header: missing recognizable {joined} column(s). "
        f"Available columns: {header}. Most likely the restriction file uses unsupported column names, "
        "has no header row, or was parsed with the wrong delimiter. Rename the header to recognizable "
        "SNP or CHR/POS column names."
    )


def resolve_restriction_rsid_column(
    columns: Iterable[str],
    *,
    context: str | None = None,
) -> str:
    """Resolve the single relevant column for an ``rsid`` restriction file."""
    columns = list(columns)
    try:
        return infer_snp_column(columns, context=context)
    except (ValueError, LDSCInputError) as exc:
        raise _restriction_header_error(columns, ("SNP",)) from exc


def resolve_restriction_chr_pos_columns(
    columns: Iterable[str],
    *,
    genome_build: str | None = None,
    context: str | None = None,
) -> tuple[str, str]:
    """Prefer the matching build-specific position; otherwise use generic POS.

    Build validation belongs to the workflow; this reader selects columns
    using its resolved build. Unused coordinate alternatives are ignored.
    """
    columns = list(columns)
    try:
        normalized_build = normalize_genome_build(genome_build)
    except LDSCConfigError as exc:
        raise LDSCConfigError(
            f"Could not resolve SNP restriction columns: genome_build={genome_build!r} is invalid. "
            "Most likely a misspelled build was passed while selecting build-specific position columns. "
            "Use genome_build='hg19', genome_build='hg38', genome_build='auto', or None."
        ) from exc

    try:
        chr_col = resolve_required_column(columns, RESTRICTION_CHRPOS_SPEC_MAP["CHR"], context=context)
    except (ValueError, LDSCInputError) as exc:
        raise _restriction_header_error(columns, ("CHR", "POS")) from exc

    # Select a matching build-specific representation before generic POS.
    specific = {
        build: resolve_optional_column(columns, spec, context=context)
        for build, spec in (("hg19", RESTRICTION_HG19_POS_SPEC),
                            ("hg38", RESTRICTION_HG38_POS_SPEC))
    }
    if normalized_build in specific and specific[normalized_build] is not None:
        return chr_col, specific[normalized_build]
    generic_columns = [column for column in columns if column not in specific.values()]
    pos_col = resolve_optional_column(generic_columns, RESTRICTION_CHRPOS_SPEC_MAP["POS"], context=context)
    if pos_col is not None:
        return chr_col, pos_col
    if normalized_build == "auto":
        matches = [value for value in specific.values() if value is not None]
        if len(matches) == 1:
            return chr_col, matches[0]
        if len(matches) > 1:
            raise LDSCInputError(
                f"Could not resolve SNP restriction position column in {context}: "
                "multiple build-specific position columns matched. Resolve the workflow "
                "genome build before selecting hg19_POS or hg38_POS."
            )
    raise _restriction_header_error(columns, ("POS",))
