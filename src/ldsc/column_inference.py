"""Central naming and alias inference registry for the refactored package.

This module is the single source of truth for:

- column-header normalization
- column-alias families scoped by input context
- SNP identifier mode aliases
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

LOGGER = logging.getLogger("LDSC.columns")
_LOGGED_INFERENCES: set[tuple[str, str, str]] = set()


def clean_header(header: str) -> str:
    """Normalize a raw header using the legacy ``munge_sumstats`` rules."""
    return str(header).upper().replace("-", "_").replace(".", "_").replace("\n", "").strip()


def normalize_column_token(value: str) -> str:
    """Normalize a column-like token for permissive alias matching."""
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


_RSID_IDENTIFIER_ALIASES = ("rsid", "rs_id", "rs", "snp", "snpid", "snp_id")
_CHR_POS_IDENTIFIER_ALIASES = ("chr_pos", "chrpos", "chrom_pos", "chromosome_position", "position")
_HG19_BUILD_ALIASES = ("hg19", "hg37", "grch37")
_HG38_BUILD_ALIASES = ("hg38", "grch38")
_AUTO_BUILD_ALIASES = ("auto",)


def normalize_snp_identifier_mode(value: str) -> str:
    """Normalize a flexible SNP identifier-mode label to ``rsid`` or ``chr_pos``."""
    normalized = normalize_column_token(value)
    if normalized in {normalize_column_token(alias) for alias in _RSID_IDENTIFIER_ALIASES}:
        return "rsid"
    if normalized in {normalize_column_token(alias) for alias in _CHR_POS_IDENTIFIER_ALIASES}:
        return "chr_pos"
    raise ValueError(f"Unsupported snp_identifier mode: {value!r}")


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
    raise ValueError("genome_build must be None, 'auto', 'hg19', 'hg37', 'GRCh37', 'hg38', or 'GRCh38'.")


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


def _index_spec(base_aliases: Sequence[str], index: int, canonical: str, label: str) -> ColumnSpec:
    """Build a numbered column spec such as ``rsID_1`` or ``rsID_2``."""
    values: list[str] = []
    for alias in base_aliases:
        values.extend((f"{alias}_{index}", f"{alias}{index}"))
    return ColumnSpec(canonical, tuple(values), label)


_POS_TOKEN_ALIASES = ("pos", "bp", "position", "base_pair", "basepair")
_UNIQ_ID_TOKEN_ALIASES = ("uniq_id", "uniqid", "unique_id", "uniqueid")


def _build_position_aliases(build_aliases: Sequence[str], index: int) -> tuple[str, ...]:
    """Build numbered position aliases for one genome-build family."""
    values: list[str] = []
    for build in build_aliases:
        for pos in _POS_TOKEN_ALIASES:
            values.extend((f"{build}_{pos}_{index}", f"{build}_{pos}{index}"))
    return tuple(values)


def _build_uniq_id_aliases(build_aliases: Sequence[str], index: int) -> tuple[str, ...]:
    """Build numbered unique-ID aliases for one genome-build family."""
    values: list[str] = []
    for build in build_aliases:
        for token in _UNIQ_ID_TOKEN_ALIASES:
            values.extend((f"{build}_{token}_{index}", f"{build}_{token}{index}"))
    return tuple(values)


# Shared external column families.
SNP_COLUMN_ALIASES = ("SNP", "SNPID", "SNP_ID", "RSID", "RS_ID", "RS", "ID", "MARKERNAME", "MARKER")
CHR_COLUMN_ALIASES = ("CHR", "#CHROM", "CHROM", "CHROMOSOME")
POS_COLUMN_ALIASES = ("POS", "BP", "POSITION", "BASE_PAIR", "BASEPAIR")
CM_COLUMN_ALIASES = ("CM", "CMBP", "CENTIMORGAN")
MAF_COLUMN_ALIASES = ("MAF", "FRQ", "FREQ", "FREQUENCY")

CHR_COLUMN_SPEC = ColumnSpec("CHR", CHR_COLUMN_ALIASES, "chromosome")
POS_COLUMN_SPEC = ColumnSpec("POS", POS_COLUMN_ALIASES, "position")
SNP_COLUMN_SPEC = ColumnSpec("SNP", SNP_COLUMN_ALIASES, "SNP identifier")
CM_COLUMN_SPEC = ColumnSpec("CM", CM_COLUMN_ALIASES, "centiMorgan")
MAF_COLUMN_SPEC = ColumnSpec("MAF", MAF_COLUMN_ALIASES, "minor-allele frequency")


# Raw summary-stat alias families preserve legacy behavior and add current
# non-conflicting aliases where semantics match.
RAW_SUMSTATS_REQUIRED_OR_OPTIONAL_SPECS = (
    ColumnSpec(
        "SNP",
        ("SNP", "MARKERNAME", "SNPID", "SNP_ID", "RS", "RSID", "RS_ID", "ID", "RS_NUMBER", "RS_NUMBERS", "MARKER"),
        "summary-stat SNP identifier",
    ),
    CHR_COLUMN_SPEC,
    POS_COLUMN_SPEC,
    ColumnSpec("NSTUDY", ("NSTUDY", "N_STUDY", "NSTUDIES", "N_STUDIES"), "summary-stat number of studies"),
    ColumnSpec("P", ("P", "PVALUE", "P_VALUE", "PVAL", "P_VAL", "GC_PVALUE"), "summary-stat p-value"),
    ColumnSpec("A1", ("A1", "ALLELE1", "ALLELE_1", "EFFECT_ALLELE", "REFERENCE_ALLELE", "INC_ALLELE", "EA"), "effect allele"),
    ColumnSpec("A2", ("A2", "ALLELE2", "ALLELE_2", "OTHER_ALLELE", "NON_EFFECT_ALLELE", "DEC_ALLELE", "NEA"), "non-effect allele"),
    ColumnSpec("N", ("N", "WEIGHT"), "summary-stat sample size"),
    ColumnSpec("N_CAS", ("NCASE", "CASES_N", "N_CASE", "N_CASES", "N_CAS", "NCAS"), "summary-stat case count"),
    ColumnSpec("N_CON", ("N_CONTROLS", "N_CON", "NCONTROL", "CONTROLS_N", "N_CONTROL", "NCON"), "summary-stat control count"),
    ColumnSpec("INFO", ("INFO",), "summary-stat INFO"),
    ColumnSpec("FRQ", ("EAF", "FRQ", "MAF", "FRQ_U", "F_U"), "summary-stat allele frequency"),
)
RAW_SUMSTATS_REQUIRED_OR_OPTIONAL_SPEC_MAP = {
    spec.canonical: spec for spec in RAW_SUMSTATS_REQUIRED_OR_OPTIONAL_SPECS
}

RAW_SUMSTATS_SIGNED_STAT_SPECS = (
    ColumnSpec("Z", ("ZSCORE", "Z-SCORE", "GC_ZSCORE", "Z"), "summary-stat Z-score"),
    ColumnSpec("OR", ("OR",), "summary-stat odds ratio"),
    ColumnSpec("BETA", ("B", "BETA", "EFFECTS", "EFFECT"), "summary-stat beta"),
    ColumnSpec("LOG_ODDS", ("LOG_ODDS",), "summary-stat log odds"),
    ColumnSpec("SIGNED_SUMSTAT", ("SIGNED_SUMSTAT",), "explicit signed summary statistic"),
)
RAW_SUMSTATS_SIGNED_STAT_SPEC_MAP = {
    spec.canonical: spec for spec in RAW_SUMSTATS_SIGNED_STAT_SPECS
}


# Context families for permissive external inputs.
ANNOTATION_METADATA_SPECS = (CHR_COLUMN_SPEC, POS_COLUMN_SPEC, SNP_COLUMN_SPEC, CM_COLUMN_SPEC, MAF_COLUMN_SPEC)
ANNOTATION_METADATA_SPEC_MAP = {spec.canonical: spec for spec in ANNOTATION_METADATA_SPECS}

REFERENCE_METADATA_SPECS = (CHR_COLUMN_SPEC, POS_COLUMN_SPEC, SNP_COLUMN_SPEC, CM_COLUMN_SPEC, MAF_COLUMN_SPEC)
REFERENCE_METADATA_SPEC_MAP = {spec.canonical: spec for spec in REFERENCE_METADATA_SPECS}

RESTRICTION_RSID_SPECS = (SNP_COLUMN_SPEC,)
RESTRICTION_RSID_SPEC_MAP = {spec.canonical: spec for spec in RESTRICTION_RSID_SPECS}

RESTRICTION_CHRPOS_SPECS = (CHR_COLUMN_SPEC, POS_COLUMN_SPEC)
RESTRICTION_CHRPOS_SPEC_MAP = {spec.canonical: spec for spec in RESTRICTION_CHRPOS_SPECS}
RESTRICTION_HG19_POS_SPEC = ColumnSpec("POS", tuple(f"{build}_{pos}" for build in _HG19_BUILD_ALIASES for pos in _POS_TOKEN_ALIASES), "hg19 position")
RESTRICTION_HG38_POS_SPEC = ColumnSpec("POS", tuple(f"{build}_{pos}" for build in _HG38_BUILD_ALIASES for pos in _POS_TOKEN_ALIASES), "hg38 position")


# Strict internal artifact families.
INTERNAL_SUMSTATS_ARTIFACT_SPECS = (
    _canonical_only_spec("SNP", "munged sumstats SNP identifier"),
    _canonical_only_spec("CHR", "munged sumstats chromosome"),
    _canonical_only_spec("POS", "munged sumstats position"),
    _canonical_only_spec("N", "munged sumstats sample size"),
    _canonical_only_spec("Z", "munged sumstats Z-score"),
    _canonical_only_spec("A1", "munged sumstats allele 1"),
    _canonical_only_spec("A2", "munged sumstats allele 2"),
    _canonical_only_spec("FRQ", "munged sumstats allele frequency"),
)
INTERNAL_SUMSTATS_ARTIFACT_SPEC_MAP = {
    spec.canonical: spec for spec in INTERNAL_SUMSTATS_ARTIFACT_SPECS
}

INTERNAL_ANNOT_ARTIFACT_SPECS = (
    _canonical_only_spec("CHR", "internal annotation chromosome"),
    _canonical_only_spec("POS", "internal annotation position"),
    _canonical_only_spec("SNP", "internal annotation SNP identifier"),
    _canonical_only_spec("CM", "internal annotation centiMorgan"),
    _canonical_only_spec("MAF", "internal annotation minor-allele frequency"),
)
INTERNAL_ANNOT_ARTIFACT_SPEC_MAP = {
    spec.canonical: spec for spec in INTERNAL_ANNOT_ARTIFACT_SPECS
}

INTERNAL_LDSCORE_ARTIFACT_SPECS = (
    _canonical_only_spec("CHR", "internal LD-score chromosome"),
    _canonical_only_spec("SNP", "internal LD-score SNP identifier"),
    _canonical_only_spec("POS", "internal LD-score position"),
    _canonical_only_spec("CM", "internal LD-score centiMorgan"),
    _canonical_only_spec("MAF", "internal LD-score minor-allele frequency"),
)
INTERNAL_LDSCORE_ARTIFACT_SPEC_MAP = {
    spec.canonical: spec for spec in INTERNAL_LDSCORE_ARTIFACT_SPECS
}


# Canonical parquet R2 schema inference.
PARQUET_R2_CANONICAL_SPECS = (
    ColumnSpec("CHR", CHR_COLUMN_ALIASES, "canonical parquet R2 chromosome", allow_suffix_match=False),
    ColumnSpec(
        "POS_1",
        ("POS_1", "POS1", "BP_1", "BP1", "POSITION_1", "POSITION1"),
        "canonical parquet R2 left position",
        allow_suffix_match=False,
    ),
    ColumnSpec(
        "POS_2",
        ("POS_2", "POS2", "BP_2", "BP2", "POSITION_2", "POSITION2"),
        "canonical parquet R2 right position",
        allow_suffix_match=False,
    ),
    _index_spec(SNP_COLUMN_ALIASES, 1, "SNP_1", "canonical parquet R2 left SNP"),
    _index_spec(SNP_COLUMN_ALIASES, 2, "SNP_2", "canonical parquet R2 right SNP"),
    ColumnSpec("R2", ("R2",), "canonical parquet R2 value", allow_suffix_match=False),
)
PARQUET_R2_CANONICAL_SPEC_MAP = {spec.canonical: spec for spec in PARQUET_R2_CANONICAL_SPECS}


# Legacy/raw parquet R2 schema inference.
R2_SOURCE_COLUMN_SPECS = (
    ColumnSpec("chr", CHR_COLUMN_ALIASES, "R2 chromosome"),
    _index_spec(SNP_COLUMN_ALIASES, 1, "rsID_1", "R2 rsID_1"),
    _index_spec(SNP_COLUMN_ALIASES, 2, "rsID_2", "R2 rsID_2"),
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


def build_cleaned_alias_lookup(specs: Sequence[ColumnSpec]) -> dict[str, str]:
    """Build a legacy-style cleaned-header-to-canonical alias map."""
    alias_map: dict[str, str] = {}
    for spec in specs:
        values = (spec.canonical,) + tuple(spec.aliases)
        for value in values:
            key = clean_header(value)
            existing = alias_map.get(key)
            if existing is not None and existing != spec.canonical:
                raise ValueError(
                    f"Alias {value!r} maps to both {existing!r} and {spec.canonical!r} in the central registry."
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
    """Emit a one-time info log when alias inference renames a source column."""
    if actual == canonical:
        return
    key = (context or "", canonical, actual)
    if key in _LOGGED_INFERENCES:
        return
    _LOGGED_INFERENCES.add(key)
    suffix = "" if not context else f" in {context}"
    message = f"Inferred canonical field '{canonical}' from input column '{actual}'{suffix}."
    LOGGER.info(message)


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
        raise ValueError(
            f"Ambiguous {spec.label} columns for canonical field {spec.canonical!r}: "
            + ", ".join(matches)
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


def infer_chr_bp_columns(header: Iterable[str], *, context: str | None = None) -> tuple[str, str]:
    """Backward-compatible alias for :func:`infer_chr_pos_columns`."""
    return infer_chr_pos_columns(header, context=context)


def _restriction_header_error(columns: Iterable[str], missing: Sequence[str]) -> ValueError:
    """Build the shared error when a restriction header lacks required aliases."""
    header = ", ".join(str(column) for column in columns)
    joined = ", ".join(missing)
    return ValueError(
        "Restriction files must contain a header row with recognizable "
        f"{joined} column(s). Available columns: {header}"
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
    except ValueError as exc:
        raise _restriction_header_error(columns, ("SNP",)) from exc


def resolve_restriction_chr_pos_columns(
    columns: Iterable[str],
    *,
    genome_build: str | None = None,
    context: str | None = None,
) -> tuple[str, str]:
    """Resolve the relevant columns for a ``chr_pos`` restriction file."""
    columns = list(columns)
    normalized_build = normalize_genome_build(genome_build)

    try:
        chr_col = resolve_required_column(columns, RESTRICTION_CHRPOS_SPEC_MAP["CHR"], context=context)
    except ValueError as exc:
        raise _restriction_header_error(columns, ("CHR", "POS")) from exc

    try:
        pos_col = resolve_required_column(columns, RESTRICTION_CHRPOS_SPEC_MAP["POS"], context=context)
        return chr_col, pos_col
    except ValueError as generic_error:
        # Fall through to build-specific position aliases before reporting the
        # generic POS failure to the caller.
        pass

    if normalized_build == "hg19":
        try:
            return chr_col, resolve_required_column(columns, RESTRICTION_HG19_POS_SPEC, context=context)
        except ValueError as exc:
            raise _restriction_header_error(columns, ("POS",)) from exc
    if normalized_build == "hg38":
        try:
            return chr_col, resolve_required_column(columns, RESTRICTION_HG38_POS_SPEC, context=context)
        except ValueError as exc:
            raise _restriction_header_error(columns, ("POS",)) from exc

    if normalized_build == "auto":
        matches: list[str] = []
        for spec in (RESTRICTION_HG19_POS_SPEC, RESTRICTION_HG38_POS_SPEC):
            match = resolve_optional_column(columns, spec, context=context)
            if match is not None and match not in matches:
                matches.append(match)
        if len(matches) == 1:
            return chr_col, matches[0]
        if len(matches) > 1:
            raise ValueError(
                "Restriction file exposes multiple build-specific position columns; "
                "set genome_build to 'hg19' or 'hg38', or rename the desired column to POS/BP."
            ) from generic_error

    raise _restriction_header_error(columns, ("POS",)) from generic_error
