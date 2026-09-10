"""Annotation-name validation and advisory value classification.

The helpers in this module describe annotation values without changing them.
LD-score construction and regression diagnostics share these definitions
without changing the numerical calculations.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from .errors import LDSCInputError


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


__all__ = [
    "classify_annotation_values",
    "require_unique_annotation_names",
]
