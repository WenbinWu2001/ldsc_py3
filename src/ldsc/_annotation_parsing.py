"""Shared annotation chunk normalization; no source reads or retained matrices."""

import logging
import numpy as np
import pandas as pd

from ._kernel import annotation as kernel_annotation
from ._kernel.snp_identity import identity_mode_family
from .chromosome_inference import normalize_chromosome
from .column_inference import (
    A1_COLUMN_SPEC, A2_COLUMN_SPEC, ANNOTATION_METADATA_SPEC_MAP,
    ColumnSpec, CHR_COLUMN_SPEC, CM_COLUMN_SPEC, POS_COLUMN_SPEC,
    SNP_COLUMN_SPEC, resolve_optional_column, resolve_required_column,
)
from .errors import LDSCInputError

LOGGER = logging.getLogger("LDSC.annotation")
_ANNOTATION_A1_COLUMN_SPEC = ColumnSpec(A1_COLUMN_SPEC.canonical, A1_COLUMN_SPEC.aliases, A1_COLUMN_SPEC.label, allow_suffix_match=False)
_ANNOTATION_A2_COLUMN_SPEC = ColumnSpec(A2_COLUMN_SPEC.canonical, A2_COLUMN_SPEC.aliases, A2_COLUMN_SPEC.label, allow_suffix_match=False)


def normalize_annotation_chunk(df, path, snp_identifier, chrom=None, *, log_ignored_metadata=True):
    """Normalize a bounded frame to aligned metadata and float32 values.

    Chunk readers enable ``log_ignored_metadata`` only on their first chunk,
    so the CM/MAF notice appears once per file read without process-wide state.
    """
    context = str(path)
    chr_col = resolve_required_column(df.columns, CHR_COLUMN_SPEC, context=context)
    pos_col = resolve_required_column(df.columns, POS_COLUMN_SPEC, context=context)
    # SNP is the identity key only in rsID-family modes; coordinate-family modes
    # key on CHR/POS, so SNP is optional there.
    if identity_mode_family(snp_identifier) == "rsid":
        snp_col = resolve_required_column(df.columns, SNP_COLUMN_SPEC, context=context)
    else:
        snp_col = resolve_optional_column(df.columns, SNP_COLUMN_SPEC, context=context)
    # CM/MAF are population-specific; the reference panel is authoritative. CM is
    # optional here and kept only as a placeholder for the legacy .annot layout
    # and the BED-query projection; ldscore ignores annotation CM/MAF.
    cm_col = resolve_optional_column(df.columns, CM_COLUMN_SPEC, context=context)
    a1_col = resolve_optional_column(df.columns, _ANNOTATION_A1_COLUMN_SPEC, context=context)
    a2_col = resolve_optional_column(df.columns, _ANNOTATION_A2_COLUMN_SPEC, context=context)
    if (a1_col is None) ^ (a2_col is None):
        raise LDSCInputError(
            f"annotate could not read annotation file '{path}': exactly one allele column was found. "
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
        LOGGER.info(
            f"Annotation file '{path}' contains CM/MAF columns; these are ignored "
            "(the reference panel is authoritative for CM and MAF)."
        )
    if chrom is not None:
        keep = metadata["CHR"] == normalize_chromosome(chrom, context=context)
        metadata = metadata.loc[keep].reset_index(drop=True)
        df = df.loc[keep].reset_index(drop=True)

    if len(metadata) == 0:
        return metadata, pd.DataFrame(index=metadata.index)

    metadata_source_columns = {chr_col, pos_col, snp_col, cm_col, maf_col, a1_col, a2_col}
    annotation_columns = [column for column in df.columns if column not in metadata_source_columns]
    if not annotation_columns:
        raise LDSCInputError(
            f"annotate could not read annotation file '{path}': no annotation value columns remain after metadata columns. "
            "Most likely the file contains only CHR/POS/SNP/CM/allele metadata. Add at least one numeric annotation column."
        )
    annotations = kernel_annotation._validate_annotation_values(
        df,
        annotation_columns,
        path=path,
    ).reset_index(drop=True)
    metadata = metadata.reset_index(drop=True)
    return metadata, annotations
