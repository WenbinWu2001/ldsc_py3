"""Shared SNP row-alignment checks for normalized LDSC3 tables."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .coordinates import positive_int_position_series
from .chromosomes import chrom_sort_key
from .snp_identity import effective_merge_key_series, identity_mode_family, is_allele_aware_mode
from .columns import POS_COLUMN_SPEC, resolve_optional_column
from .snp_identity import normalize_snp_identifier_mode
from ..errors import LDSCInputError, LDSCInternalError


NUMERIC_METADATA_COLUMNS = ("CM", "MAF")
FLOAT_METADATA_TOLERANCE = float(np.finfo(np.float16).eps)


def assert_same_snp_rows(
    left: pd.DataFrame,
    right: pd.DataFrame,
    *,
    context: str,
    snp_identifier: str = "chr_pos_allele_aware",
) -> None:
    """Raise if two normalized tables do not share identical SNP rows."""
    mode = normalize_snp_identifier_mode(snp_identifier)
    _require_columns(left, side="left", context=context, snp_identifier=mode)
    _require_columns(right, side="right", context=context, snp_identifier=mode)
    if len(left) != len(right):
        raise LDSCInputError(
            f"{context}: row count mismatch ({len(left)} != {len(right)}). "
            "Most likely paired LDSC3 artifacts were generated from different SNP sets. "
            "Regenerate the paired artifacts together with the same SNP identifier mode."
        )

    extra_columns = _extra_row_metadata_columns(left, right, snp_identifier=mode)
    left_keys = _row_key_frame(left, side="left", context=context, snp_identifier=mode, extra_columns=extra_columns)
    right_keys = _row_key_frame(right, side="right", context=context, snp_identifier=mode, extra_columns=extra_columns)
    mismatched = left_keys.ne(right_keys).any(axis=1)
    if mismatched.any():
        row = int(np.flatnonzero(mismatched.to_numpy())[0])
        columns = left_keys.columns[left_keys.loc[row].ne(right_keys.loc[row])].tolist()
        raise LDSCInputError(
            f"{context}: SNP row mismatch at row {row} in {columns}: "
            f"left={left_keys.loc[row].to_dict()}, right={right_keys.loc[row].to_dict()}. "
            "Most likely paired LDSC3 artifacts were sorted or filtered differently. "
            "Regenerate the paired artifacts together from the same SNPs."
        )

    for column in NUMERIC_METADATA_COLUMNS:
        if column in left.columns and column in right.columns:
            _assert_close_numeric_metadata(left, right, column=column, context=context)


def _require_columns(frame: pd.DataFrame, *, side: str, context: str, snp_identifier: str) -> None:
    required = ["SNP"] if identity_mode_family(snp_identifier) == "rsid" else ["CHR", "POS"]
    if is_allele_aware_mode(snp_identifier):
        required.extend(["A1", "A2"])
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise LDSCInputError(
            f"{context}: {side} table is missing required SNP row columns: {missing}. "
            "Most likely the artifact was written by an older LDSC3 version or with a different "
            "SNP identifier mode. Regenerate the artifact with the current workflow."
        )


def _extra_row_metadata_columns(left: pd.DataFrame, right: pd.DataFrame, *, snp_identifier: str) -> list[str]:
    if identity_mode_family(snp_identifier) != "rsid":
        return []
    return [column for column in ("CHR", "POS") if column in left.columns and column in right.columns]


def _row_key_frame(
    frame: pd.DataFrame,
    *,
    side: str,
    context: str,
    snp_identifier: str,
    extra_columns: list[str],
) -> pd.DataFrame:
    positional = isinstance(frame.index, pd.RangeIndex) and frame.index.start == 0 and frame.index.step == 1
    key = effective_merge_key_series(
        frame if positional else frame.reset_index(drop=True),
        snp_identifier,
        context=f"{context}: {side}",
    ).reset_index(drop=True)
    keys = {"identity_key": key}
    if "CHR" in extra_columns:
        keys["CHR"] = frame["CHR"].astype(str).reset_index(drop=True)
    if "POS" in extra_columns:
        keys["POS"] = positive_int_position_series(
            frame["POS"].reset_index(drop=True),
            context=f"{context}: {side}",
            label="POS",
        )
    return pd.DataFrame(keys)


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
        raise LDSCInputError(
            f"{context}: numeric metadata column {column!r} differs at row {row}: "
            f"left={left_values[row]!r}, right={right_values[row]!r}. "
            "Most likely paired LDSC3 artifacts have inconsistent metadata sidecars. "
            "Regenerate the paired artifacts together."
        )


def sort_frame_by_genomic_position(df: pd.DataFrame) -> pd.DataFrame:
    """Sort by chromosome/position and use a present SNP label as a stable tie-breaker.

    Only the sort keys are sorted; the rows are then taken once by the
    resulting order. Already-ordered input is returned as a shallow copy with a
    fresh RangeIndex, so callers must not modify its values in place.
    """
    pos_col = resolve_optional_column(df.columns, POS_COLUMN_SPEC)
    if pos_col is None:
        raise KeyError("No POS-like column available for genomic sorting.")
    labels = df["CHR"].to_numpy()
    keys = {"_chrom_key": pd.Series(labels).map({label: chrom_sort_key(label) for label in pd.unique(labels)}),
            pos_col: df[pos_col].to_numpy()}
    if "SNP" in df.columns:
        keys["SNP"] = df["SNP"].to_numpy()
    order = pd.DataFrame(keys).sort_values(by=list(keys), kind="mergesort").index.to_numpy()
    if np.array_equal(order, np.arange(len(order))):
        ordered = df.copy(deep=False)
    else:
        ordered = df.take(order)
    ordered.index = pd.RangeIndex(len(ordered))
    return ordered


def first_identity_mismatch(left: pd.DataFrame, right: pd.DataFrame, columns) -> tuple[int, str] | None:
    """Return the first ``(row, column)`` where two tables' identity columns differ, else ``None``.

    Values are compared exactly and in row order (missing equals missing); a
    row-count difference is reported at the end of the shorter table. This is
    stricter than :func:`assert_same_snp_rows`, which compares identifier-mode keys.
    """
    for column in columns:
        expected, actual = left[column].to_numpy(), right[column].to_numpy()
        if len(expected) != len(actual):
            return min(len(expected), len(actual)), column
        same = (expected == actual) | (pd.isna(expected) & pd.isna(actual))
        if not same.all():
            return int(np.flatnonzero(~same)[0]), column
    return None


def concat_chromosome_tables(pieces) -> pd.DataFrame:
    """Concatenate ``(chrom, table)`` pieces, each already in genomic order, in chromosome order.

    The result equals ``sort_frame_by_genomic_position`` of the concatenation:
    a stable sort cannot move rows that are already ordered. Instead of sorting,
    each piece is checked to hold only its own chromosome, in nondecreasing
    position with a nondecreasing SNP label among equal positions; a piece that
    breaks this raises ``LDSCInternalError`` before any row is published.
    """
    ordered = sorted(pieces, key=lambda piece: chrom_sort_key(piece[0]))
    for chrom, table in ordered:
        if table.empty:
            continue
        if {chrom_sort_key(label) for label in table["CHR"].unique()} != {chrom_sort_key(chrom)}:
            raise LDSCInternalError(
                f"Chromosome {chrom} LD-score rows carry other chromosome labels; the per-chromosome "
                "tables cannot be concatenated in genomic order. Report the command and traceback."
            )
        pos_col = resolve_optional_column(table.columns, POS_COLUMN_SPEC)
        step = np.diff(table[pos_col].to_numpy())
        unordered = step < 0
        ties = np.flatnonzero(step == 0)
        if len(ties) and "SNP" in table.columns:
            snps = table["SNP"].to_numpy()
            unordered[ties] = snps[ties + 1] < snps[ties]
        if unordered.any():
            row = int(np.flatnonzero(unordered)[0]) + 1
            raise LDSCInternalError(
                f"Chromosome {chrom} LD-score rows are not in genomic order at row {row}; the "
                "per-chromosome tables cannot be concatenated. Report the command and traceback."
            )
    return pd.concat([table for _, table in ordered], ignore_index=True)
