"""Concrete annotation table parsing and numeric-value validation.

These helpers are shared by annotation preparation and LD-score readers.
Interval parsing lives in ``_core.regions`` and live BED/gene projection in
``_annotation_projection``. Workflow orchestration and output policy remain
with their public workflow and writer owners.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from ..errors import LDSCInputError

_ANNOTATION_VALUES_DOC = "docs/user-guide/troubleshooting.md#ldscore-annotation-values-are-malformed"


def _read_text_table(path: str) -> pd.DataFrame:
    """Read a whitespace-delimited annotation table with optional gzip compression."""
    compression = "gzip" if str(path).endswith(".gz") else None
    try:
        return pd.read_csv(path, sep=r"\s+", compression=compression)
    except (pd.errors.ParserError, UnicodeDecodeError, ValueError) as exc:
        raise LDSCInputError(_annotation_parse_error_message(path, details=str(exc))) from exc


def _validate_annotation_values(
    frame: pd.DataFrame,
    annotation_columns: Sequence[str],
    *,
    path: str | Path,
    dtype=np.float32,
) -> pd.DataFrame:
    """Validate numeric, non-missing annotations and apply the requested dtype.

    ``dtype=None`` preserves parsed precision so exact binary classification
    can precede the ordinary float32 normalization.
    """
    values = frame.loc[:, annotation_columns]
    numeric = values.apply(pd.to_numeric, errors="coerce")
    missing = values.isna()
    non_numeric = numeric.isna() & ~missing
    if missing.any().any() or non_numeric.any().any():
        problems: list[str] = []
        for column in annotation_columns:
            missing_rows = np.flatnonzero(missing[column].to_numpy()) + 1
            non_numeric_rows = np.flatnonzero(non_numeric[column].to_numpy()) + 1
            if len(missing_rows):
                problems.append(
                    f"'{column}': {_counted(len(missing_rows), 'missing value')} at data row(s) "
                    f"{_summarize_row_numbers(missing_rows)}"
                )
            if len(non_numeric_rows):
                problems.append(
                    f"'{column}': {_counted(len(non_numeric_rows), 'non-numeric value')} at data row(s) "
                    f"{_summarize_row_numbers(non_numeric_rows)}"
                )
        raise LDSCInputError(_annotation_parse_error_message(path, details="; ".join(problems)))
    return numeric if dtype is None else numeric.astype(dtype)


def _summarize_row_numbers(rows: np.ndarray, *, limit: int = 5) -> str:
    """Format a bounded set of one-based data-row numbers for an error."""
    shown = ", ".join(str(int(row)) for row in rows[:limit])
    return shown if len(rows) <= limit else f"{shown}, ..."


def _counted(count: int, noun: str) -> str:
    """Return a count with a minimally pluralized noun."""
    return f"{count} {noun}{'' if count == 1 else 's'}"


def _annotation_parse_error_message(path: str | Path, *, details: str) -> str:
    """Build the shared actionable annotation-table format error."""
    return (
        f"Could not parse annotation file '{path}': annotation value columns must be numeric "
        f"and non-missing. Affected values: {details}. Most likely a row has fewer fields "
        "than the header, delimiters are inconsistent, an empty placeholder was collapsed "
        "by whitespace parsing, a row/file is truncated, or an annotation value is not "
        "numeric. Ensure each row has one field per header column, use the package's `NA` "
        "token only for missing metadata values, and provide numeric, non-missing annotation values. "
        f"Other causes & fixes: {_ANNOTATION_VALUES_DOC}"
    )
