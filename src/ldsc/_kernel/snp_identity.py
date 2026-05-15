"""Shared SNP identity policy for artifact and restriction matching."""

from __future__ import annotations

from dataclasses import dataclass
import logging
from typing import Literal

import pandas as pd

from .._coordinates import CHR_POS_KEY_COLUMN, build_chr_pos_key_frame

SNPIdentifierMode = Literal["rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"]
SNP_IDENTIFIER_MODES = ("rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware")
RSID_FAMILY_MODES = ("rsid", "rsid_allele_aware")
CHR_POS_FAMILY_MODES = ("chr_pos", "chr_pos_allele_aware")
ALLELE_AWARE_MODES = ("rsid_allele_aware", "chr_pos_allele_aware")
BASE_MODES = ("rsid", "chr_pos")
SCHEMA_VERSION = 1
ARTIFACT_TYPES = frozenset({"sumstats", "ref_panel_r2", "ref_panel_metadata", "ldscore"})
REGENERATE_ARTIFACT_MESSAGE = (
    "This artifact was not written with the current LDSC schema/provenance contract. "
    "Regenerate it with the current LDSC package."
)

IDENTITY_DROP_COLUMNS = [
    "CHR",
    "SNP",
    "source_pos",
    "target_pos",
    "reason",
    "base_key",
    "identity_key",
    "allele_set",
    "stage",
]
IDENTITY_DROP_REASONS = (
    "missing_allele",
    "invalid_allele",
    "strand_ambiguous_allele",
    "multi_allelic_base_key",
    "duplicate_identity",
)
COMPLEMENT = str.maketrans({"A": "T", "T": "A", "C": "G", "G": "C"})
_VALID_ALLELES = frozenset({"A", "C", "G", "T"})
_INTERNAL_COLUMNS = ["_ldsc_base_key", "_ldsc_allele_set", "_ldsc_identity_key"]


@dataclass(frozen=True)
class IdentityCleanupResult:
    """Artifact table after SNP identity cleanup plus rows removed by policy."""

    cleaned: pd.DataFrame
    dropped: pd.DataFrame


@dataclass(frozen=True)
class RestrictionIdentityKeys:
    """Collapsed restriction identity keys and rows rejected while building them."""

    keys: set[str]
    match_kind: Literal["base", "identity"]
    dropped: pd.DataFrame
    n_input_rows: int
    n_retained_keys: int


def normalize_snp_identifier_mode(value: str) -> str:
    """Validate one exact public SNP identifier mode."""
    if value in SNP_IDENTIFIER_MODES:
        return value
    allowed = ", ".join(repr(mode) for mode in SNP_IDENTIFIER_MODES)
    raise ValueError(f"Unsupported snp_identifier mode: {value!r}. Expected one of {allowed}.")


def identity_artifact_metadata(
    *,
    artifact_type: str,
    snp_identifier: str,
    genome_build: str | None,
) -> dict[str, object]:
    """Return the minimal identity metadata persisted with reloadable artifacts."""
    if artifact_type not in ARTIFACT_TYPES:
        raise ValueError(f"Unknown LDSC artifact_type {artifact_type!r}.")
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_type": artifact_type,
        "snp_identifier": normalize_snp_identifier_mode(snp_identifier),
        "genome_build": genome_build,
    }


def validate_identity_artifact_metadata(metadata: dict[str, object], *, expected_artifact_type: str) -> str:
    """Validate minimal identity metadata and return the normalized SNP identifier mode."""
    if metadata.get("schema_version") != SCHEMA_VERSION or metadata.get("artifact_type") != expected_artifact_type:
        raise ValueError(REGENERATE_ARTIFACT_MESSAGE)
    return normalize_snp_identifier_mode(str(metadata.get("snp_identifier")))


def identity_mode_family(mode: str) -> Literal["rsid", "chr_pos"]:
    """Return the base family for an exact SNP identifier mode."""
    mode = normalize_snp_identifier_mode(mode)
    if mode in RSID_FAMILY_MODES:
        return "rsid"
    return "chr_pos"


def identity_base_mode(mode: str) -> Literal["rsid", "chr_pos"]:
    """Return the allele-blind base mode for ``mode``."""
    return identity_mode_family(mode)


def is_allele_aware_mode(mode: str) -> bool:
    """Return whether ``mode`` includes normalized A1/A2 identity."""
    return normalize_snp_identifier_mode(mode) in ALLELE_AWARE_MODES


def _clean_allele(value: object) -> str | None:
    """Return an uppercase allele token, or ``None`` for missing/empty values."""
    if pd.isna(value):
        return None
    token = str(value).strip().upper()
    if not token:
        return None
    return token


def _allele_failure_reason(a1: object, a2: object) -> str | None:
    """Classify an allele pair failure using public identity drop reasons."""
    left = _clean_allele(a1)
    right = _clean_allele(a2)
    if left is None or right is None:
        return "missing_allele"
    if left not in _VALID_ALLELES or right not in _VALID_ALLELES:
        return "invalid_allele"
    if left == right:
        return "invalid_allele"
    if frozenset((left, right)) in {frozenset(("A", "T")), frozenset(("C", "G"))}:
        return "strand_ambiguous_allele"
    return None


def normalize_allele_set(a1: object, a2: object) -> object:
    """
    Normalize an unordered, strand-aware allele pair.

    Invalid, missing, identical, or strand-ambiguous pairs raise ``ValueError``
    with the public drop reason as the message.
    Valid pairs return the lexicographically smallest token among the observed
    unordered pair and its unordered complement pair.
    """
    reason = _allele_failure_reason(a1, a2)
    if reason is not None:
        raise ValueError(reason)
    left = _clean_allele(a1)
    right = _clean_allele(a2)
    assert left is not None and right is not None
    observed = ":".join(sorted((left, right)))
    complement = ":".join(sorted((left.translate(COMPLEMENT), right.translate(COMPLEMENT))))
    return min(observed, complement)


def _require_columns(frame: pd.DataFrame, columns: tuple[str, ...], *, context: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        joined = ", ".join(missing)
        raise ValueError(f"{context} requires canonical column(s): {joined}.")


def base_key_series(frame: pd.DataFrame, mode: str, *, context: str) -> pd.Series:
    """Build allele-blind base identity keys for an exact SNP identifier mode."""
    family = identity_mode_family(mode)
    if family == "rsid":
        _require_columns(frame, ("SNP",), context=context)
        return frame["SNP"].astype("string").astype(object)

    _require_columns(frame, ("CHR", "POS"), context=context)
    result = pd.Series(pd.NA, index=frame.index, dtype=object)
    keyed, _report = build_chr_pos_key_frame(
        frame.loc[:, ["CHR", "POS"]].copy(),
        context=context,
        drop_missing=True,
    )
    if not keyed.empty:
        result.loc[keyed.index] = keyed[CHR_POS_KEY_COLUMN].astype("string").astype(object)
    return result


def allele_set_series(frame: pd.DataFrame, *, context: str) -> tuple[pd.Series, pd.Series]:
    """Return normalized allele sets plus per-row failure reasons."""
    _require_columns(frame, ("A1", "A2"), context=context)
    values: list[object] = []
    reasons: list[object] = []
    for a1, a2 in zip(frame["A1"], frame["A2"]):
        reason = _allele_failure_reason(a1, a2)
        if reason is None:
            values.append(normalize_allele_set(a1, a2))
            reasons.append(pd.NA)
        else:
            values.append(pd.NA)
            reasons.append(reason)
    return (
        pd.Series(values, index=frame.index, dtype=object),
        pd.Series(reasons, index=frame.index, dtype=object),
    )


def effective_merge_key_series(frame: pd.DataFrame, mode: str, *, context: str = "table") -> pd.Series:
    """Build the active effective SNP merge key for a canonical table."""
    mode = normalize_snp_identifier_mode(mode)
    base = base_key_series(frame, mode, context=context)
    if not is_allele_aware_mode(mode):
        return base.rename("snp_id")

    allele_set, reasons = allele_set_series(frame, context=context)
    invalid = reasons.notna()
    if bool(invalid.any()):
        reason = str(reasons.loc[invalid].iloc[0])
        raise ValueError(f"{context} contains invalid allele identity rows: {reason}.")
    keys = pd.Series(pd.NA, index=frame.index, dtype=object)
    valid = base.notna() & allele_set.notna()
    keys.loc[valid] = base.loc[valid].astype(str) + ":" + allele_set.loc[valid].astype(str)
    return keys.rename("snp_id")


def empty_identity_drop_frame() -> pd.DataFrame:
    """Return an empty identity drop-report frame with the public schema."""
    return pd.DataFrame(columns=IDENTITY_DROP_COLUMNS)


def coerce_identity_drop_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Coerce ``frame`` to the identity drop-report schema and column order."""
    if frame is None or frame.empty:
        return empty_identity_drop_frame()
    coerced = frame.copy()
    for column in IDENTITY_DROP_COLUMNS:
        if column not in coerced.columns:
            coerced[column] = pd.NA
    coerced = coerced.loc[:, IDENTITY_DROP_COLUMNS].reset_index(drop=True)
    for column in ("source_pos", "target_pos"):
        original = coerced[column]
        numeric = pd.to_numeric(original, errors="coerce")
        if bool(original.isna().equals(numeric.isna())):
            coerced[column] = numeric.astype("Int64")
    return coerced


def _identity_drop_rows(frame: pd.DataFrame, *, reason: str, stage: str) -> pd.DataFrame:
    """Build drop-report rows for one policy reason."""
    if reason not in IDENTITY_DROP_REASONS:
        raise ValueError(f"Unsupported identity drop reason: {reason!r}.")
    if frame.empty:
        return empty_identity_drop_frame()

    def column_or_na(column: str) -> pd.Series:
        if column in frame.columns:
            return frame[column].reset_index(drop=True)
        return pd.Series([pd.NA] * len(frame), dtype=object)

    rows = pd.DataFrame(
        {
            "CHR": column_or_na("CHR"),
            "SNP": column_or_na("SNP"),
            "source_pos": column_or_na("POS") if "POS" in frame.columns else pd.Series(frame.index.to_list(), dtype=object),
            "target_pos": column_or_na("target_pos"),
            "reason": reason,
            "base_key": column_or_na("_ldsc_base_key"),
            "identity_key": column_or_na("_ldsc_identity_key"),
            "allele_set": column_or_na("_ldsc_allele_set"),
            "stage": stage,
        }
    )
    return coerce_identity_drop_frame(rows)


def clean_identity_artifact_table(
    frame: pd.DataFrame,
    mode: str,
    *,
    context: str,
    stage: str,
    logger: logging.Logger | None = None,
) -> IdentityCleanupResult:
    """Drop rows that cannot provide a unique effective SNP identity."""
    mode = normalize_snp_identifier_mode(mode)
    work = frame.copy()
    work["_ldsc_base_key"] = base_key_series(work, mode, context=context)
    drop_frames: list[pd.DataFrame] = []

    if is_allele_aware_mode(mode):
        allele_set, reasons = allele_set_series(work, context=context)
        work["_ldsc_allele_set"] = allele_set
        invalid_mask = reasons.notna()
        for reason in IDENTITY_DROP_REASONS[:3]:
            reason_mask = reasons == reason
            if bool(reason_mask.any()):
                drop_frames.append(_identity_drop_rows(work.loc[reason_mask], reason=reason, stage=stage))
        work = work.loc[~invalid_mask].copy()

        multi_mask = _multi_allelic_base_key_mask(work)
        if bool(multi_mask.any()):
            drop_frames.append(_identity_drop_rows(work.loc[multi_mask], reason="multi_allelic_base_key", stage=stage))
            work = work.loc[~multi_mask].copy()

        work["_ldsc_identity_key"] = _combine_base_and_allele_keys(work["_ldsc_base_key"], work["_ldsc_allele_set"])
    else:
        work["_ldsc_allele_set"] = pd.NA
        work["_ldsc_identity_key"] = work["_ldsc_base_key"]

    duplicate_mask = _duplicate_identity_mask(work["_ldsc_identity_key"])
    if bool(duplicate_mask.any()):
        drop_frames.append(_identity_drop_rows(work.loc[duplicate_mask], reason="duplicate_identity", stage=stage))
        work = work.loc[~duplicate_mask].copy()

    cleaned = work.drop(columns=[column for column in _INTERNAL_COLUMNS if column in work.columns]).reset_index(drop=True)
    dropped = _combine_drop_frames(drop_frames)
    _log_identity_drops(dropped, context=context, logger=logger)
    return IdentityCleanupResult(cleaned=cleaned, dropped=dropped)


def collapse_restriction_identity_keys(
    frame: pd.DataFrame,
    mode: str,
    *,
    context: str,
    has_allele_columns: bool | None = None,
    logger: logging.Logger | None = None,
) -> RestrictionIdentityKeys:
    """Build set-like restriction keys under the active identity mode."""
    mode = normalize_snp_identifier_mode(mode)
    n_input = int(len(frame))
    has_alleles = ("A1" in frame.columns and "A2" in frame.columns) if has_allele_columns is None else has_allele_columns
    base = base_key_series(frame, mode, context=context)

    if not is_allele_aware_mode(mode) or not has_alleles:
        keys = {str(value) for value in base.dropna().tolist()}
        return RestrictionIdentityKeys(
            keys=keys,
            match_kind="base",
            dropped=empty_identity_drop_frame(),
            n_input_rows=n_input,
            n_retained_keys=len(keys),
        )

    work = frame.copy()
    work["_ldsc_base_key"] = base
    allele_set, reasons = allele_set_series(work, context=context)
    work["_ldsc_allele_set"] = allele_set
    drop_frames: list[pd.DataFrame] = []

    invalid_mask = reasons.notna()
    for reason in IDENTITY_DROP_REASONS[:3]:
        reason_mask = reasons == reason
        if bool(reason_mask.any()):
            drop_frames.append(_identity_drop_rows(work.loc[reason_mask], reason=reason, stage="restriction"))
    work = work.loc[~invalid_mask].copy()

    multi_mask = _multi_allelic_base_key_mask(work)
    if bool(multi_mask.any()):
        drop_frames.append(_identity_drop_rows(work.loc[multi_mask], reason="multi_allelic_base_key", stage="restriction"))
        work = work.loc[~multi_mask].copy()

    work["_ldsc_identity_key"] = _combine_base_and_allele_keys(work["_ldsc_base_key"], work["_ldsc_allele_set"])
    keys = {str(value) for value in work["_ldsc_identity_key"].dropna().tolist()}
    dropped = _combine_drop_frames(drop_frames)
    _log_identity_drops(dropped, context=context, logger=logger)
    return RestrictionIdentityKeys(
        keys=keys,
        match_kind="identity",
        dropped=dropped,
        n_input_rows=n_input,
        n_retained_keys=len(keys),
    )


def _combine_base_and_allele_keys(base: pd.Series, allele_set: pd.Series) -> pd.Series:
    keys = pd.Series(pd.NA, index=base.index, dtype=object)
    valid = base.notna() & allele_set.notna()
    keys.loc[valid] = base.loc[valid].astype(str) + ":" + allele_set.loc[valid].astype(str)
    return keys


def _multi_allelic_base_key_mask(frame: pd.DataFrame) -> pd.Series:
    if frame.empty:
        return pd.Series(False, index=frame.index, dtype=bool)
    valid = frame["_ldsc_base_key"].notna() & frame["_ldsc_allele_set"].notna()
    counts = frame.loc[valid].groupby("_ldsc_base_key")["_ldsc_allele_set"].nunique()
    multi_keys = counts[counts > 1].index
    return frame["_ldsc_base_key"].isin(multi_keys)


def _duplicate_identity_mask(identity_key: pd.Series) -> pd.Series:
    valid = identity_key.notna()
    return valid & identity_key.duplicated(keep=False)


def _combine_drop_frames(drop_frames: list[pd.DataFrame]) -> pd.DataFrame:
    if not drop_frames:
        return empty_identity_drop_frame()
    return coerce_identity_drop_frame(pd.concat(drop_frames, ignore_index=True))


def _log_identity_drops(dropped: pd.DataFrame, *, context: str, logger: logging.Logger | None) -> None:
    if logger is None or dropped.empty:
        return
    counts = dropped["reason"].value_counts().to_dict()
    logger.warning("Dropped %d SNP identity rows in %s: %s.", len(dropped), context, counts)
