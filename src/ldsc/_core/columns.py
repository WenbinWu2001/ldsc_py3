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


def _canonical_only_spec(canonical: str, label: str) -> ColumnSpec:
    """Return a spec that accepts only the canonical field name."""
    return ColumnSpec(canonical, (canonical,), label)


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


# Raw summary-stat alias families preserve legacy behavior and add current
# non-conflicting aliases where semantics match. Allele aliases are only names
# that say which allele the signed statistic refers to; see ORIENTATION_AMBIGUOUS_ALLELE_COLUMNS
# in _sumstats/munge_input.py for allele names that are never mapped.
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
        ("A1", "EFFECT_ALLELE", "INC_ALLELE", "EA"),
        "effect allele (A1), the allele the signed statistic refers to",
    ),
    ColumnSpec(
        "A2",
        ("A2", "OTHER_ALLELE", "NON_EFFECT_ALLELE", "DEC_ALLELE", "NEA"),
        "other allele (A2)",
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


# Context families for external inputs.
# Every family matches a header only by its canonical name or a listed alias.
# In annotation files every column that is not metadata is an annotation, so a
# looser rule would silently move an annotation such as lowMAF, eQTL_SNP or
# Enhancers (which ends in the alias RS) out of the model.
ANNOTATION_METADATA_SPECS = (
    CHR_COLUMN_SPEC,
    POS_COLUMN_SPEC,
    SNP_COLUMN_SPEC,
    CM_COLUMN_SPEC,
    MAF_COLUMN_SPEC,
    A1_COLUMN_SPEC,
    A2_COLUMN_SPEC,
)
ANNOTATION_METADATA_SPEC_MAP = {spec.canonical: spec for spec in ANNOTATION_METADATA_SPECS}

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
    _canonical_only_spec("A1", "munged summary-statistics A1 (effect allele)"),
    _canonical_only_spec("A2", "munged summary-statistics A2 (other allele)"),
    _canonical_only_spec("FRQ", "munged summary-statistics allele frequency"),
)
INTERNAL_SUMSTATS_ARTIFACT_SPEC_MAP = {
    spec.canonical: spec for spec in INTERNAL_SUMSTATS_ARTIFACT_SPECS
}


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
    """Return the headers named ``spec.canonical``, or failing that, the headers named by an alias.

    Headers are compared after :func:`normalize_column_token`, so case and
    punctuation never matter, but a header that only ends in a name (``EUR_MAF``)
    matches nothing.
    """
    canonical_norm = normalize_column_token(spec.canonical)
    alias_norms = {normalize_column_token(alias) for alias in spec.aliases}
    columns = list(columns)
    canonical = [column for column in columns if normalize_column_token(column) == canonical_norm]
    if canonical:
        return canonical
    return [column for column in columns if normalize_column_token(column) in alias_norms]


_MISSING_COLUMN_TROUBLESHOOTING = "docs/user-guide/troubleshooting.md#common-a-required-column-is-not-recognized"


def lookalike_hint(columns: Iterable[str], spec: ColumnSpec) -> str:
    """Return a sentence naming headers that end in a name of ``spec`` after a separator, or ``""``.

    The hint only helps a user repair a missing required column; such headers
    are never matched. A name that ends inside one word (``lowMAF``, ``Thyroid``)
    is not counted.
    """
    names = {normalize_column_token(name) for name in (spec.canonical, *spec.aliases)}
    found = []
    for column in columns:
        tokens = [token for token in re.split(r"[^0-9A-Za-z]+", str(column)) if token]
        if normalize_column_token(column) in names:
            continue
        if any("".join(tokens[start:]).lower() in names for start in range(1, len(tokens))):
            found.append(str(column))
    if not found:
        return ""
    listed = ", ".join(found)
    noun = "a header that ends" if len(found) == 1 else "headers that end"
    return (
        f" The header has {listed}, {noun} in a name for {spec.canonical} without being one; "
        f"if it holds the {spec.label}, rename it to {spec.canonical}."
    )


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
    columns = list(columns)
    matches = _best_matches(columns, spec)
    if not matches:
        header = ", ".join(str(column) for column in columns)
        context_text = f" in {context}" if context else ""
        raise LDSCInputError(
            f"Could not map the required {spec.label} column{context_text}: expected canonical field {spec.canonical!r}; "
            f"available columns: {header}. Most likely the header uses an unrecognized alias or the file delimiter is wrong. "
            "Pass the matching column explicitly when the workflow supports it, or rename the input header."
            + lookalike_hint(columns, spec)
            + f" Other causes & fixes: {_MISSING_COLUMN_TROUBLESHOOTING}"
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
    """Build the shared error when a restriction header lacks required aliases.

    A build-specific position column is never offered as the generic ``POS``:
    in a run on another build, renaming it would bring in the wrong coordinates.
    """
    columns = list(columns)
    header = ", ".join(str(column) for column in columns)
    joined = ", ".join(missing)
    build_specific = {column for spec in (RESTRICTION_HG19_POS_SPEC, RESTRICTION_HG38_POS_SPEC)
                      for column in _best_matches(columns, spec)}
    generic = [column for column in columns if column not in build_specific]
    specs = {"SNP": RESTRICTION_RSID_SPEC_MAP["SNP"], **RESTRICTION_CHRPOS_SPEC_MAP}
    hints = "".join(lookalike_hint(generic, specs[name]) for name in missing)
    return LDSCInputError(
        f"Could not read SNP restriction file header: missing recognizable {joined} column(s). "
        f"Available columns: {header}. Most likely the restriction file uses unsupported column names, "
        "has no header row, or was parsed with the wrong delimiter. Rename the header to recognizable "
        "SNP or CHR/POS column names." + hints
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
