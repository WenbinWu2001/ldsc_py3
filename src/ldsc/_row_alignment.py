"""Shared SNP row-alignment checks for normalized LDSC tables."""

from __future__ import annotations

import numpy as np
import pandas as pd


SNP_ROW_COLUMNS = ("CHR", "SNP", "POS")
NUMERIC_METADATA_COLUMNS = ("CM", "MAF")
FLOAT_METADATA_TOLERANCE = float(np.finfo(np.float16).eps)


def assert_same_snp_rows(left: pd.DataFrame, right: pd.DataFrame, *, context: str) -> None:
    """Raise if two normalized tables do not share identical SNP rows."""
    _require_columns(left, side="left", context=context)
    _require_columns(right, side="right", context=context)
    if len(left) != len(right):
        raise ValueError(f"{context}: row count mismatch ({len(left)} != {len(right)}).")

    left_keys = _row_key_frame(left, side="left", context=context)
    right_keys = _row_key_frame(right, side="right", context=context)
    mismatched = left_keys.ne(right_keys).any(axis=1)
    if mismatched.any():
        row = int(np.flatnonzero(mismatched.to_numpy())[0])
        columns = left_keys.columns[left_keys.loc[row].ne(right_keys.loc[row])].tolist()
        raise ValueError(
            f"{context}: SNP row mismatch at row {row} in {columns}: "
            f"left={left_keys.loc[row].to_dict()}, right={right_keys.loc[row].to_dict()}."
        )

    for column in NUMERIC_METADATA_COLUMNS:
        if column in left.columns and column in right.columns:
            _assert_close_numeric_metadata(left, right, column=column, context=context)


def _require_columns(frame: pd.DataFrame, *, side: str, context: str) -> None:
    missing = [column for column in SNP_ROW_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(f"{context}: {side} table is missing required SNP row columns: {missing}")


def _row_key_frame(frame: pd.DataFrame, *, side: str, context: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "CHR": frame["CHR"].astype(str).reset_index(drop=True),
            "SNP": frame["SNP"].astype(str).reset_index(drop=True),
            "POS": _integer_pos_series(frame["POS"], side=side, context=context),
        }
    )


def _integer_pos_series(values: pd.Series, *, side: str, context: str) -> pd.Series:
    numeric = pd.to_numeric(values, errors="raise").reset_index(drop=True)
    missing = numeric.isna()
    if missing.any():
        row = int(np.flatnonzero(missing.to_numpy())[0])
        raise ValueError(f"{context}: {side} POS is missing at row {row}.")
    non_integral = (numeric % 1) != 0
    if non_integral.any():
        row = int(np.flatnonzero(non_integral.to_numpy())[0])
        raise ValueError(f"{context}: {side} POS must be integer-valued; got {numeric.iloc[row]!r} at row {row}.")
    return numeric.astype(np.int64)


def _numeric_array(frame: pd.DataFrame, column: str) -> np.ndarray:
    values = pd.to_numeric(frame[column], errors="raise").reset_index(drop=True)
    return values.to_numpy(dtype=np.float64, na_value=np.nan)


def _assert_close_numeric_metadata(
    left: pd.DataFrame,
    right: pd.DataFrame,
    *,
    column: str,
    context: str,
) -> None:
    left_values = _numeric_array(left, column)
    right_values = _numeric_array(right, column)
    close = np.isclose(
        left_values,
        right_values,
        rtol=FLOAT_METADATA_TOLERANCE,
        atol=FLOAT_METADATA_TOLERANCE,
        equal_nan=True,
    )
    if not bool(np.all(close)):
        row = int(np.flatnonzero(~close)[0])
        raise ValueError(
            f"{context}: numeric metadata column {column!r} differs at row {row}: "
            f"left={left_values[row]!r}, right={right_values[row]!r}."
        )
