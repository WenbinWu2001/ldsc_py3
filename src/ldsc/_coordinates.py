"""Shared CHR/POS normalization and missing-coordinate policy helpers."""

from __future__ import annotations

from dataclasses import dataclass
import logging
from typing import Sequence

import pandas as pd

from .chromosome_inference import normalize_chromosome_series

CHR_POS_KEY_COLUMN = "_ldsc_chr_pos_key"
NA_COORDINATE_TOKENS = frozenset({"", ".", "na", "nan", "none", "null"})


@dataclass(frozen=True)
class CoordinateDropReport:
    """Summary of rows dropped because CHR or POS was missing."""

    context: str
    n_input: int
    n_dropped: int
    n_retained: int
    n_missing_chr: int
    n_missing_pos: int
    examples: list[dict[str, object]]


def coordinate_missing_mask(values: pd.Series) -> pd.Series:
    """Return rows whose coordinate token is missing or package NA-like."""
    tokens = values.astype("string")
    return tokens.isna() | tokens.str.strip().str.lower().isin(NA_COORDINATE_TOKENS)


def complete_coordinate_mask(frame: pd.DataFrame, *, chr_col: str = "CHR", pos_col: str = "POS") -> pd.Series:
    """Return rows with both CHR and POS present under the shared NA-token policy."""
    return ~(coordinate_missing_mask(frame[chr_col]) | coordinate_missing_mask(frame[pos_col]))


def positive_int_position_series(values: pd.Series, *, context: str, label: str = "POS") -> pd.Series:
    """Convert non-missing POS values to positive int64 positions or raise."""
    tokens = values.astype("string")
    missing = coordinate_missing_mask(values)
    numeric = pd.to_numeric(values, errors="coerce")

    invalid = (~missing) & numeric.isna()
    if bool(invalid.any()):
        bad_value = tokens.loc[invalid].iloc[0]
        raise ValueError(f"{label} values in {context} must be numeric; got {bad_value!r}.")

    if bool(missing.any()):
        bad_row = int(missing.to_numpy().nonzero()[0][0])
        raise ValueError(f"{label} values in {context} are missing at row {bad_row}.")

    non_integral = (numeric % 1) != 0
    if bool(non_integral.any()):
        bad_value = values.loc[non_integral].iloc[0]
        raise ValueError(f"{label} values in {context} must be integer base-pair positions; got {bad_value!r}.")

    non_positive = numeric <= 0
    if bool(non_positive.any()):
        bad_value = values.loc[non_positive].iloc[0]
        raise ValueError(f"{label} values in {context} must be positive integer base-pair positions; got {bad_value!r}.")

    return numeric.astype("int64")


def normalize_chr_pos_frame(
    frame: pd.DataFrame,
    *,
    context: str,
    chr_col: str = "CHR",
    pos_col: str = "POS",
    drop_missing: bool = False,
    logger: logging.Logger | None = None,
    log_level: int = logging.INFO,
    example_columns: Sequence[str] | None = None,
) -> tuple[pd.DataFrame, CoordinateDropReport]:
    """
    Normalize CHR/POS columns and optionally drop rows missing either coordinate.

    Missing coordinates are blank or NA-like tokens and may be dropped at
    match/map boundaries. Invalid non-missing coordinates remain hard errors.
    """
    work = frame.copy()
    missing = ~complete_coordinate_mask(work, chr_col=chr_col, pos_col=pos_col)
    report = _drop_report(
        work,
        missing,
        context=context,
        chr_col=chr_col,
        pos_col=pos_col,
        example_columns=example_columns,
    )
    if report.n_dropped:
        if not drop_missing:
            first = report.examples[0] if report.examples else {}
            raise ValueError(f"{context} requires complete CHR/POS coordinates; first missing row: {first}.")
        _log_drop_report(report, logger=logger, level=log_level)
        work = work.loc[~missing].copy()

    if work.empty:
        return work, report

    work[chr_col] = normalize_chromosome_series(work[chr_col], context=context).astype(object)
    work[pos_col] = positive_int_position_series(work[pos_col], context=context, label=pos_col)
    return work, report


def build_chr_pos_key_frame(
    frame: pd.DataFrame,
    *,
    context: str,
    chr_col: str = "CHR",
    pos_col: str = "POS",
    key_col: str = CHR_POS_KEY_COLUMN,
    drop_missing: bool = False,
    logger: logging.Logger | None = None,
    log_level: int = logging.INFO,
    example_columns: Sequence[str] | None = None,
) -> tuple[pd.DataFrame, CoordinateDropReport]:
    """Return a normalized copy with a canonical private ``CHR:POS`` key."""
    keyed, report = normalize_chr_pos_frame(
        frame,
        context=context,
        chr_col=chr_col,
        pos_col=pos_col,
        drop_missing=drop_missing,
        logger=logger,
        log_level=log_level,
        example_columns=example_columns,
    )
    if keyed.empty:
        keyed[key_col] = pd.Series(dtype="string")
        return keyed, report
    keyed[key_col] = keyed[chr_col].astype(str) + ":" + keyed[pos_col].astype("int64").astype(str)
    return keyed, report


def _drop_report(
    frame: pd.DataFrame,
    missing: pd.Series,
    *,
    context: str,
    chr_col: str,
    pos_col: str,
    example_columns: Sequence[str] | None,
) -> CoordinateDropReport:
    missing_chr = coordinate_missing_mask(frame[chr_col])
    missing_pos = coordinate_missing_mask(frame[pos_col])
    examples = _example_rows(frame, missing, chr_col=chr_col, pos_col=pos_col, example_columns=example_columns)
    dropped = int(missing.sum())
    return CoordinateDropReport(
        context=context,
        n_input=int(len(frame)),
        n_dropped=dropped,
        n_retained=int(len(frame) - dropped),
        n_missing_chr=int(missing_chr.sum()),
        n_missing_pos=int(missing_pos.sum()),
        examples=examples,
    )


def _example_rows(
    frame: pd.DataFrame,
    missing: pd.Series,
    *,
    chr_col: str,
    pos_col: str,
    example_columns: Sequence[str] | None,
) -> list[dict[str, object]]:
    if not bool(missing.any()):
        return []
    if example_columns is None:
        columns = [column for column in ("SNP", chr_col, pos_col) if column in frame.columns]
    else:
        columns = [column for column in example_columns if column in frame.columns]
    if not columns:
        columns = [chr_col, pos_col]
    return frame.loc[missing, columns].head(5).to_dict(orient="records")


def _log_drop_report(report: CoordinateDropReport, *, logger: logging.Logger | None, level: int) -> None:
    if report.n_dropped == 0 or logger is None:
        return
    logger.log(
        level,
        "Dropped %d SNPs with missing CHR/POS in %s; %d rows remain.",
        report.n_dropped,
        report.context,
        report.n_retained,
    )
    if report.examples:
        logger.debug("Example rows dropped for missing CHR/POS in %s: %s", report.context, report.examples)
