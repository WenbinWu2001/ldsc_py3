"""Advisory annotation classification and semantic fingerprints.

The helpers in this module describe annotation values without changing them.
They are shared by LD-score construction, regression diagnostics, and the
post-fit quantile workflow so all three boundaries use the same definitions.
"""

from __future__ import annotations

import hashlib
from collections.abc import Sequence

import numpy as np
import pandas as pd

from .errors import LDSCInputError


FINGERPRINT_ALGORITHM = "sha256"
FINGERPRINT_CANONICALIZATION = "ldsc_common_annotation_v1"


def require_unique_annotation_names(
    baseline_names: Sequence[str], query_names: Sequence[str] = ()
) -> None:
    """Require annotation names to be globally unique across fitted groups.

    Parameters
    ----------
    baseline_names, query_names : sequence of str
        Ordered names associated with the baseline and query matrices.

    Raises
    ------
    LDSCInputError
        If any name occurs more than once within or across the two groups.
    """
    names = [str(name) for name in (*baseline_names, *query_names)]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise LDSCInputError(
            "Annotation names must be globally unique within an LD-score artifact. "
            f"Duplicate annotation name(s): {', '.join(duplicates)}. Rename the conflicting "
            "baseline or query columns before computing LD scores."
        )


def classify_annotation_values(values: pd.DataFrame) -> dict[str, str]:
    """Classify finite annotation columns as binary or quantitative.

    A column is binary only when every retained value is exactly zero or one.
    The classification is descriptive: the input values are not modified and
    the result must not influence LD-score or regression calculations.

    Parameters
    ----------
    values : pandas.DataFrame
        Numeric SNP-by-annotation values. Missing and nonfinite values are not
        accepted.

    Returns
    -------
    dict of str to {"binary", "quantitative"}
        Classification in input column order.

    Raises
    ------
    LDSCInputError
        If a column contains a nonnumeric, missing, or infinite value.
    """
    classifications: dict[str, str] = {}
    for column in values.columns:
        numeric = pd.to_numeric(values[column], errors="coerce").to_numpy(dtype=np.float64)
        if not np.isfinite(numeric).all():
            raise LDSCInputError(
                f"Annotation '{column}' contains missing, nonnumeric, or infinite values. "
                "Fitted annotations must be finite numeric values for every retained SNP."
            )
        classifications[str(column)] = "binary" if np.isin(numeric, (0.0, 1.0)).all() else "quantitative"
    return classifications


def common_universe_fingerprint(effective_snp_ids: Sequence[object] | pd.Series) -> str:
    """Hash sorted effective identities in the common reference-SNP universe."""
    identities = _normalized_identities(effective_snp_ids)
    payload = "effective_snp_id\n" + "".join(f"{identity}\n" for identity in identities)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def annotation_value_fingerprint(
    effective_snp_ids: Sequence[object] | pd.Series,
    values: Sequence[object] | pd.Series,
) -> str:
    """Hash one annotation over the common reference-SNP universe.

    Values are normalized to float32, paired with effective identities, sorted
    by identity, formatted with 17 significant digits, and serialized as a
    headered UTF-8 TSV with ``\n`` line endings.
    """
    identities = pd.Series(effective_snp_ids, dtype="string").reset_index(drop=True)
    numeric = pd.to_numeric(pd.Series(values).reset_index(drop=True), errors="coerce").to_numpy(dtype=np.float64)
    if len(identities) != len(numeric):
        raise LDSCInputError("Cannot fingerprint annotation values because identity and value lengths differ.")
    if identities.isna().any() or not np.isfinite(numeric).all():
        raise LDSCInputError("Cannot fingerprint annotation values with missing identities or nonfinite values.")
    frame = pd.DataFrame(
        {
            "effective_snp_id": identities.astype(str),
            "annotation_value": numeric.astype(np.float32),
        }
    ).sort_values("effective_snp_id", kind="stable")
    if frame["effective_snp_id"].duplicated().any():
        raise LDSCInputError("Cannot fingerprint annotation values with duplicate effective SNP identities.")
    lines = ["effective_snp_id\tannotation_value\n"]
    lines.extend(
        f"{row.effective_snp_id}\t{format(float(row.annotation_value), '.17g')}\n"
        for row in frame.itertuples(index=False)
    )
    return hashlib.sha256("".join(lines).encode("utf-8")).hexdigest()


def build_annotation_fingerprint_metadata(
    effective_snp_ids: Sequence[object] | pd.Series,
    annotation_values: pd.DataFrame,
) -> dict[str, object]:
    """Build versioned common reference-SNP universe and annotation fingerprints."""
    return {
        "algorithm": FINGERPRINT_ALGORITHM,
        "canonicalization": FINGERPRINT_CANONICALIZATION,
        "common_reference_snp_universe": common_universe_fingerprint(effective_snp_ids),
        "annotation_values": {
            str(column): annotation_value_fingerprint(effective_snp_ids, annotation_values[column])
            for column in annotation_values.columns
        },
    }


def _normalized_identities(effective_snp_ids: Sequence[object] | pd.Series) -> list[str]:
    """Return unique nonmissing identities in deterministic sorted order."""
    identities = pd.Series(effective_snp_ids, dtype="string")
    if identities.isna().any():
        raise LDSCInputError("Cannot fingerprint a common reference-SNP universe with missing identities.")
    normalized = sorted(identities.astype(str).tolist())
    if len(set(normalized)) != len(normalized):
        raise LDSCInputError("Cannot fingerprint a common reference-SNP universe with duplicate identities.")
    return normalized


__all__ = [
    "FINGERPRINT_ALGORITHM",
    "FINGERPRINT_CANONICALIZATION",
    "annotation_value_fingerprint",
    "build_annotation_fingerprint_metadata",
    "classify_annotation_values",
    "common_universe_fingerprint",
    "require_unique_annotation_names",
]
