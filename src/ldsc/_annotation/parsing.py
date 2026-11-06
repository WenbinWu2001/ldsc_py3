"""Shared annotation chunk normalization and numeric-value validation.

Annotation preparation, ``ldscore``, ``quantile-h2`` and ``AnnotationBuilder``
share these helpers; they read no sources and retain no matrices.
"""

import logging
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from .._core.snp_identity import identity_mode_family
from .._core.chromosomes import normalize_chromosome
from .._core.columns import ANNOTATION_METADATA_SPEC_MAP, resolve_optional_column, resolve_required_column
from ..errors import LDSCInputError

LOGGER = logging.getLogger("LDSC.annotation")


def normalize_annotation_chunk(df, path, snp_identifier, *, log_ignored_metadata=True, value_dtype=np.float32):
    """Normalize aligned metadata and numeric values, ordinarily to float32.

    Metadata columns (``CHR``, ``POS``, ``SNP``, ``CM``, ``MAF``, ``A1``,
    ``A2``) are recognized by canonical name or listed alias only, as for every
    input, so a header that merely ends in one of them (``lowMAF``,
    ``EUR_MAF``, ``eQTL_SNP``) stays an annotation value column.

    Source preparation enables ``log_ignored_metadata`` only on the first
    classification chunk, so the DEBUG CM/MAF notice appears once per source per
    preparation. The second, numeric-only pass does not call this normalizer.
    ``value_dtype=None`` preserves parsed precision for exact binary detection.
    """
    context = str(path)
    chr_col = resolve_required_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["CHR"], context=context)
    pos_col = resolve_required_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["POS"], context=context)
    # SNP is the identity key only in rsID-family modes; coordinate-family modes
    # key on CHR/POS, so SNP is optional there.
    if identity_mode_family(snp_identifier) == "rsid":
        snp_col = resolve_required_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["SNP"], context=context)
    else:
        snp_col = resolve_optional_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["SNP"], context=context)
    # CM/MAF are population-specific; the reference panel is authoritative. CM is
    # optional here and kept only as a placeholder for the legacy .annot layout
    # and the BED-query projection; ldscore ignores annotation CM/MAF.
    cm_col = resolve_optional_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["CM"], context=context)
    a1_col = resolve_optional_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["A1"], context=context)
    a2_col = resolve_optional_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["A2"], context=context)
    if (a1_col is None) ^ (a2_col is None):
        raise LDSCInputError(
            f"Could not read annotation file '{path}': exactly one allele column was found. "
            "Most likely the file has A1 without A2 or A2 without A1. Provide both allele columns or neither."
        )
    columns = {"CHR": df[chr_col], "POS": df[pos_col]}
    if snp_col is not None:
        columns["SNP"] = df[snp_col]
    metadata = pd.DataFrame(columns)
    metadata["CHR"] = metadata["CHR"].map(lambda value: normalize_chromosome(value, context=context))
    metadata["POS"] = pd.to_numeric(metadata["POS"], errors="raise").astype(np.int64)
    if "SNP" in metadata.columns:
        metadata["SNP"] = metadata["SNP"].astype(str)
    # CM is a population-agnostic placeholder (always NaN): ldscore sources CM from
    # the reference panel, and the legacy .annot layout only needs the column present.
    # Any CM value in the input is intentionally discarded.
    metadata["CM"] = np.nan
    if a1_col is not None and a2_col is not None:
        metadata["A1"] = df[a1_col]
        metadata["A2"] = df[a2_col]
    # MAF is population-specific and never carried into annotation metadata: resolve
    # the column only to keep it out of the annotation value columns.
    maf_col = resolve_optional_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["MAF"], context=context)
    if log_ignored_metadata and (cm_col is not None or maf_col is not None):
        LOGGER.debug(
            f"Annotation file '{path}' contains CM/MAF columns; these are ignored "
            "(the reference panel is authoritative for CM and MAF)."
        )
    if len(metadata) == 0:
        return metadata, pd.DataFrame(index=metadata.index)

    metadata_source_columns = {chr_col, pos_col, snp_col, cm_col, maf_col, a1_col, a2_col}
    annotation_columns = [column for column in df.columns if column not in metadata_source_columns]
    if not annotation_columns:
        raise LDSCInputError(
            f"Could not read annotation file '{path}': no annotation value columns remain after metadata columns. "
            "Most likely the file contains only CHR/POS/SNP/CM/allele metadata. Add at least one numeric annotation column."
        )
    annotations = validate_annotation_values(
        df,
        annotation_columns,
        path=path,
        dtype=value_dtype,
    ).reset_index(drop=True)
    metadata = metadata.reset_index(drop=True)
    return metadata, annotations


_ANNOTATION_VALUES_DOC = "docs/user-guide/troubleshooting.md#ldscore-annotation-values-are-malformed"


def validate_annotation_values(
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
            for mask, noun in ((missing, 'missing value'), (non_numeric, 'non-numeric value')):
                rows = np.flatnonzero(mask[column].to_numpy()) + 1
                if len(rows):
                    problems.append(
                        f"'{column}': {_counted(len(rows), noun)} at data row(s) {_summarize_row_numbers(rows)}"
                    )
        raise LDSCInputError(annotation_parse_error_message(path, details="; ".join(problems)))
    return numeric if dtype is None else numeric.astype(dtype)


def _summarize_row_numbers(rows: np.ndarray, *, limit: int = 5) -> str:
    """Format a bounded set of one-based data-row numbers for an error."""
    shown = ", ".join(str(int(row)) for row in rows[:limit])
    return shown if len(rows) <= limit else f"{shown}, ..."


def _counted(count: int, noun: str) -> str:
    """Return a count with a minimally pluralized noun."""
    return f"{count} {noun}{'' if count == 1 else 's'}"


def annotation_parse_error_message(path: str | Path, *, details: str) -> str:
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
