"""Shared CHR/POS normalization and missing-coordinate policy helpers."""

from __future__ import annotations

from dataclasses import dataclass
import logging
from typing import Literal, Sequence

import pandas as pd

from .chromosome_inference import normalize_chromosome, normalize_chromosome_series

CHR_POS_KEY_COLUMN = "_ldsc_chr_pos_key"
NA_COORDINATE_TOKENS = frozenset({"", ".", "na", "nan", "none", "null"})
CoordinatePolicy = Literal["drop", "raise"]
DEFAULT_COORDINATE_POLICY: CoordinatePolicy = "drop"
STRICT_COORDINATE_POLICY: CoordinatePolicy = "raise"


@dataclass(frozen=True)
class CoordinateDropReport:
    """Summary of rows dropped because CHR or POS was invalid or missing."""

    context: str
    n_input: int
    n_dropped: int
    n_retained: int
    n_missing_chr: int
    n_missing_pos: int
    n_invalid_chr: int
    n_invalid_pos: int
    examples: list[dict[str, object]]

    @property
    def n_invalid_or_missing(self) -> int:
        """Return the total number of rows removed by the coordinate policy."""
        return self.n_dropped


def coordinate_missing_mask(values: pd.Series) -> pd.Series:
    """Return rows whose coordinate token is missing or package NA-like."""
    tokens = values.astype("string")
    return tokens.isna() | tokens.str.strip().str.lower().isin(NA_COORDINATE_TOKENS)


def complete_coordinate_mask(frame: pd.DataFrame, *, chr_col: str = "CHR", pos_col: str = "POS") -> pd.Series:
    """Return rows with both CHR and POS present under the shared NA-token policy."""
    return ~(coordinate_missing_mask(frame[chr_col]) | coordinate_missing_mask(frame[pos_col]))


def positive_int_position_series(values: pd.Series, *, context: str, label: str = "POS") -> pd.Series:
    """Convert non-missing POS values to positive int64 positions or raise."""
    return _integer_position_series(values, context=context, label=label, min_position=1)


def _integer_position_series(
    values: pd.Series,
    *,
    context: str,
    label: str,
    min_position: int,
) -> pd.Series:
    """Convert non-missing POS values to int64 positions with a lower bound."""
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

    below_minimum = numeric < min_position
    if bool(below_minimum.any()):
        bad_value = values.loc[below_minimum].iloc[0]
        requirement = "positive" if min_position == 1 else f">= {min_position}"
        raise ValueError(f"{label} values in {context} must be {requirement} integer base-pair positions; got {bad_value!r}.")

    return numeric.astype("int64")


def normalize_chr_pos_frame(
    frame: pd.DataFrame,
    *,
    context: str,
    chr_col: str = "CHR",
    pos_col: str = "POS",
    coordinate_policy: CoordinatePolicy | str | None = None,
    drop_missing: bool | None = None,
    logger: logging.Logger | None = None,
    log_level: int = logging.WARNING,
    example_columns: Sequence[str] | None = None,
    min_position: int = 1,
) -> tuple[pd.DataFrame, CoordinateDropReport]:
    """
    Normalize CHR/POS columns under an explicit invalid/missing row policy.

    The package-wide default drops invalid or missing coordinate rows with a
    warning report. Pass ``coordinate_policy="raise"`` where a strict file
    contract requires complete, valid coordinates.
    """
    work = frame.copy()
    policy = _resolve_coordinate_policy(coordinate_policy=coordinate_policy, drop_missing=drop_missing)
    masks = _coordinate_quality_masks(work, chr_col=chr_col, pos_col=pos_col, context=context, min_position=min_position)
    droppable = masks["droppable"]
    report = _drop_report(
        work,
        droppable,
        context=context,
        chr_col=chr_col,
        pos_col=pos_col,
        missing_chr=masks["missing_chr"],
        missing_pos=masks["missing_pos"],
        invalid_chr=masks["invalid_chr"],
        invalid_pos=masks["invalid_pos"],
        example_columns=example_columns,
    )
    if report.n_dropped:
        if policy == STRICT_COORDINATE_POLICY:
            _raise_coordinate_policy_error(work, masks, context=context, chr_col=chr_col, pos_col=pos_col, min_position=min_position)
        _log_drop_report(report, logger=logger, level=log_level)
        work = work.loc[~droppable].copy()

    if work.empty:
        return work, report

    work[chr_col] = normalize_chromosome_series(work[chr_col], context=context).astype(object)
    work[pos_col] = _integer_position_series(work[pos_col], context=context, label=pos_col, min_position=min_position)
    return work, report


def build_chr_pos_key_frame(
    frame: pd.DataFrame,
    *,
    context: str,
    chr_col: str = "CHR",
    pos_col: str = "POS",
    key_col: str = CHR_POS_KEY_COLUMN,
    coordinate_policy: CoordinatePolicy | str | None = None,
    drop_missing: bool | None = None,
    logger: logging.Logger | None = None,
    log_level: int = logging.WARNING,
    example_columns: Sequence[str] | None = None,
    min_position: int = 1,
) -> tuple[pd.DataFrame, CoordinateDropReport]:
    """Return a normalized copy with a canonical private ``CHR:POS`` key."""
    keyed, report = normalize_chr_pos_frame(
        frame,
        context=context,
        chr_col=chr_col,
        pos_col=pos_col,
        coordinate_policy=coordinate_policy,
        drop_missing=drop_missing,
        logger=logger,
        log_level=log_level,
        example_columns=example_columns,
        min_position=min_position,
    )
    if keyed.empty:
        keyed[key_col] = pd.Series(dtype="string")
        return keyed, report
    keyed[key_col] = keyed[chr_col].astype(str) + ":" + keyed[pos_col].astype("int64").astype(str)
    return keyed, report


def _drop_report(
    frame: pd.DataFrame,
    droppable: pd.Series,
    *,
    context: str,
    chr_col: str,
    pos_col: str,
    missing_chr: pd.Series,
    missing_pos: pd.Series,
    invalid_chr: pd.Series,
    invalid_pos: pd.Series,
    example_columns: Sequence[str] | None,
) -> CoordinateDropReport:
    examples = _example_rows(frame, droppable, chr_col=chr_col, pos_col=pos_col, example_columns=example_columns)
    dropped = int(droppable.sum())
    return CoordinateDropReport(
        context=context,
        n_input=int(len(frame)),
        n_dropped=dropped,
        n_retained=int(len(frame) - dropped),
        n_missing_chr=int(missing_chr.sum()),
        n_missing_pos=int(missing_pos.sum()),
        n_invalid_chr=int(invalid_chr.sum()),
        n_invalid_pos=int(invalid_pos.sum()),
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
    message = (
        f"Dropped {report.n_dropped} SNPs with invalid or missing CHR/POS in {report.context}; "
        f"{report.n_retained} rows remain "
        f"(missing CHR={report.n_missing_chr}, missing POS={report.n_missing_pos}, "
        f"invalid CHR={report.n_invalid_chr}, invalid POS={report.n_invalid_pos})."
    )
    logger.log(level, message)
    if report.examples:
        logger.debug(f"Example rows dropped for invalid or missing CHR/POS in {report.context}: {report.examples}")


def _resolve_coordinate_policy(
    *,
    coordinate_policy: CoordinatePolicy | str | None,
    drop_missing: bool | None,
) -> CoordinatePolicy:
    if coordinate_policy is None:
        if drop_missing is None:
            return DEFAULT_COORDINATE_POLICY
        return "drop" if drop_missing else "raise"
    if coordinate_policy not in {"drop", "raise"}:
        raise ValueError("coordinate_policy must be 'drop' or 'raise'.")
    return coordinate_policy  # type: ignore[return-value]


def _coordinate_quality_masks(
    frame: pd.DataFrame,
    *,
    chr_col: str,
    pos_col: str,
    context: str,
    min_position: int,
) -> dict[str, pd.Series]:
    missing_chr = coordinate_missing_mask(frame[chr_col])
    missing_pos = coordinate_missing_mask(frame[pos_col])
    invalid_chr = _invalid_chromosome_mask(frame[chr_col], missing=missing_chr, context=context)
    invalid_pos = _invalid_position_mask(frame[pos_col], missing=missing_pos, min_position=min_position)
    droppable = missing_chr | missing_pos | invalid_chr | invalid_pos
    return {
        "missing_chr": missing_chr,
        "missing_pos": missing_pos,
        "invalid_chr": invalid_chr,
        "invalid_pos": invalid_pos,
        "droppable": droppable,
    }


def _invalid_chromosome_mask(values: pd.Series, *, missing: pd.Series, context: str) -> pd.Series:
    invalid_tokens: set[object] = set()
    for value in pd.unique(values.loc[~missing]):
        try:
            normalize_chromosome(value, context=context)
        except ValueError:
            invalid_tokens.add(value)
    return (~missing) & values.isin(invalid_tokens)


def _invalid_position_mask(values: pd.Series, *, missing: pd.Series, min_position: int) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    non_numeric = (~missing) & numeric.isna()
    numeric_present = (~missing) & numeric.notna()
    non_integral = numeric_present & ((numeric % 1) != 0)
    below_minimum = numeric_present & (numeric < min_position)
    return non_numeric | non_integral | below_minimum


def _raise_coordinate_policy_error(
    frame: pd.DataFrame,
    masks: dict[str, pd.Series],
    *,
    context: str,
    chr_col: str,
    pos_col: str,
    min_position: int,
) -> None:
    if bool((masks["missing_chr"] | masks["missing_pos"]).any()):
        first = _example_rows(
            frame,
            masks["missing_chr"] | masks["missing_pos"],
            chr_col=chr_col,
            pos_col=pos_col,
            example_columns=None,
        )[0]
        raise ValueError(f"{context} requires valid CHR/POS coordinates; first missing row: {first}.")
    if bool(masks["invalid_pos"].any()):
        _integer_position_series(frame[pos_col], context=context, label=pos_col, min_position=min_position)
    if bool(masks["invalid_chr"].any()):
        normalize_chromosome_series(frame[chr_col], context=context)
    first = _example_rows(frame, masks["droppable"], chr_col=chr_col, pos_col=pos_col, example_columns=None)[0]
    raise ValueError(f"{context} requires valid CHR/POS coordinates; first invalid row: {first}.")
