"""Packaged HapMap3 curated SNP map helpers."""

from __future__ import annotations

from functools import lru_cache
from importlib import resources
from os import PathLike

import pandas as pd

from .chromosome_inference import normalize_chromosome
from .column_inference import CHR_COLUMN_SPEC, SNP_COLUMN_SPEC, ColumnSpec, resolve_required_column
from ._coordinates import positive_int_position_series


HM3_HG19_POS_SPEC = ColumnSpec(
    "hg19_POS",
    (
        "hg19_POS",
        "hg19_BP",
        "hg19_POSITION",
        "hg37_POS",
        "hg37_BP",
        "GRCh37_POS",
        "GRCh37_BP",
        "GRCh37_POSITION",
    ),
    "hg19 position",
    allow_suffix_match=False,
)
HM3_HG38_POS_SPEC = ColumnSpec(
    "hg38_POS",
    (
        "hg38_POS",
        "hg38_BP",
        "hg38_POSITION",
        "GRCh38_POS",
        "GRCh38_BP",
        "GRCh38_POSITION",
    ),
    "hg38 position",
    allow_suffix_match=False,
)


def packaged_hm3_curated_map_path() -> str:
    """Return the installed path for the packaged curated dual-build HM3 map."""
    return str(resources.files("ldsc").joinpath("data", "hm3_curated_map.tsv.gz"))


def load_hm3_curated_map() -> pd.DataFrame:
    """Load the packaged curated dual-build HM3 map.

    Returns all packaged columns while canonicalizing ``CHR``, ``hg19_POS``,
    ``hg38_POS``, and ``SNP`` for reliable filtering and coordinate lookup.
    """
    return _load_hm3_curated_map_from_path(packaged_hm3_curated_map_path(), preserve_extra_columns=True)


def _load_hm3_curated_map_from_path(
    path: str | PathLike[str],
    *,
    preserve_extra_columns: bool,
) -> pd.DataFrame:
    """Load one HM3-compatible map from ``path`` with canonical core columns."""
    frame = _load_hm3_curated_map_cached(str(path), preserve_extra_columns)
    return frame.copy()


@lru_cache(maxsize=8)
def _load_hm3_curated_map_cached(path: str, preserve_extra_columns: bool) -> pd.DataFrame:
    """Cached implementation for HM3 map loading."""
    raw = pd.read_csv(path, sep=None, engine="python")
    context = str(path)
    chr_col = resolve_required_column(raw.columns, CHR_COLUMN_SPEC, context=context)
    hg19_col = resolve_required_column(raw.columns, HM3_HG19_POS_SPEC, context=context)
    hg38_col = resolve_required_column(raw.columns, HM3_HG38_POS_SPEC, context=context)
    snp_col = resolve_required_column(raw.columns, SNP_COLUMN_SPEC, context=context)
    core = pd.DataFrame(
        {
            "CHR": pd.Series(
                [normalize_chromosome(value, context=context) for value in raw[chr_col]],
                dtype="string",
            ),
            "hg19_POS": _positive_int_position(raw[hg19_col], label="hg19_POS", context=context),
            "hg38_POS": _positive_int_position(raw[hg38_col], label="hg38_POS", context=context),
            "SNP": raw[snp_col].astype("string"),
        }
    )
    _reject_duplicate_hm3_coordinates(core, build="hg19")
    _reject_duplicate_hm3_coordinates(core, build="hg38")
    if not preserve_extra_columns:
        return core

    out = raw.copy()
    out = out.drop(columns=[col for col in (chr_col, hg19_col, hg38_col, snp_col) if col in out.columns])
    return pd.concat([core, out.reset_index(drop=True)], axis=1)


def _positive_int_position(values, *, label: str, context: str) -> pd.Series:
    """Return positive integer positions with a clear HM3-specific label."""
    return positive_int_position_series(pd.Series(values), context=f"{context} {label}").astype("int64")


def _reject_duplicate_hm3_coordinates(frame: pd.DataFrame, *, build: str) -> None:
    """Reject HM3 maps that are not one-to-one by build-specific coordinates."""
    pos_col = f"{build}_POS"
    duplicate_mask = frame.duplicated(subset=["CHR", pos_col], keep=False)
    if not duplicate_mask.any():
        return
    duplicate_keys = (
        frame.loc[duplicate_mask, ["CHR", pos_col]]
        .drop_duplicates()
        .head(5)
        .apply(lambda row: f"{row['CHR']}:{int(row[pos_col])}", axis=1)
        .tolist()
    )
    raise ValueError(f"HM3 curated map contains duplicate {build} coordinates: {', '.join(duplicate_keys)}")


__all__ = ["load_hm3_curated_map"]
