"""Genetic maps: the header check, loading, and centiMorgan interpolation.

``ldscore`` (``PlinkRefPanel``), ``build-gene-ldscore-index`` and
``build-r2-matrix`` read genetic maps through these functions. Each command
keeps its own policy for which map files it accepts and when it needs them,
and passes its name as ``command``, so an error names the running command and
links its troubleshooting section (``GENETIC_MAP_TROUBLESHOOTING``). The
input gate checks map headers with ``GENETIC_MAP_HEADER_CHECKS``.
"""

from __future__ import annotations

from os import PathLike
from typing import Sequence

import numpy as np
import pandas as pd

from .._core.chromosomes import chrom_sort_key, normalize_chromosome
from .._core.columns import CHR_COLUMN_SPEC, POS_COLUMN_ALIASES, ColumnSpec, resolve_required_column
from ..errors import LDSCInputError

# Troubleshooting section for genetic-map errors, by the command that reads the map.
GENETIC_MAP_TROUBLESHOOTING = {
    "build-r2-matrix":
        "https://wenbinwu2001.github.io/ldsc_py3/user-guide/troubleshooting.html#build-r2-matrix-liftover-or-genetic-map-configuration-is-incomplete",
    "ldscore": "https://wenbinwu2001.github.io/ldsc_py3/user-guide/troubleshooting.html#ldscore-a-genetic-map-cannot-be-used",
    "build-gene-ldscore-index":
        "https://wenbinwu2001.github.io/ldsc_py3/user-guide/troubleshooting.html#build-gene-ldscore-index-a-genetic-map-cannot-be-used",
}


# Genetic maps also accept the HapMap header (Chromosome Position(bp) Rate(cM/Mb) Map(cM)).
GENETIC_MAP_POS_SPEC = ColumnSpec("POS", (*POS_COLUMN_ALIASES, "POSITION(BP)"), "genetic map position")
GENETIC_MAP_CM_SPEC = ColumnSpec(
    "CM",
    ("CM", "GENETIC_MAP_CM", "GENETICMAPCM", "GENETIC_MAP(CM)", "GENETICMAP(CM)", "MAP(CM)"),
    "genetic map centiMorgan",
)


def check_genetic_map_header(path):
    """Require chromosome, position and centimorgan columns in a genetic-map header."""
    columns = pd.read_csv(path, sep=r'\s+', nrows=0).columns
    for spec in (CHR_COLUMN_SPEC, GENETIC_MAP_POS_SPEC, GENETIC_MAP_CM_SPEC):
        resolve_required_column(columns, spec, context=str(path))


GENETIC_MAP_HEADER_CHECKS = {
    'genetic_map_hg19_sources': check_genetic_map_header,
    'genetic_map_hg38_sources': check_genetic_map_header,
}


def _read_text_table(path: str | PathLike[str]) -> pd.DataFrame:
    """Read a whitespace-delimited plain-text or gzip-compressed table."""
    text_path = str(path)
    compression = "gzip" if text_path.endswith(".gz") else None
    return pd.read_csv(text_path, sep=r"\s+", compression=compression)


def load_genetic_map(path: str | PathLike[str], *, command: str) -> pd.DataFrame:
    """Load one genetic map table into canonical `CHR/POS/CM` columns; errors name ``command``."""

    df = _read_text_table(path)
    context = str(path)
    chr_col = resolve_required_column(df.columns, CHR_COLUMN_SPEC, context=context)
    pos_col = resolve_required_column(df.columns, GENETIC_MAP_POS_SPEC, context=context)
    cm_col = resolve_required_column(df.columns, GENETIC_MAP_CM_SPEC, context=context)

    out = pd.DataFrame(
        {
            # One normalized label per distinct value, shared by its rows.
            "CHR": df[chr_col].map({value: normalize_chromosome(value, context=context)
                                    for value in pd.unique(df[chr_col])}),
            "POS": pd.to_numeric(df[pos_col], errors="raise").astype(np.int64),
            "CM": pd.to_numeric(df[cm_col], errors="raise").astype(float),
        }
    )
    duplicated = out.duplicated(subset=["CHR", "POS"], keep=False)
    if duplicated.any():
        raise LDSCInputError(
            f"{command} could not load genetic map '{path}': duplicate CHR/POS "
            "rows were found. Most likely the map file contains repeated positions. "
            "Deduplicate the map or provide a clean genetic map. "
            f"Other causes & fixes: {GENETIC_MAP_TROUBLESHOOTING[command]}"
        )

    chrom_order = out["CHR"].map(chrom_sort_key)
    if not chrom_order.equals(chrom_order.sort_values(ignore_index=True)):
        raise LDSCInputError(
            f"{command} could not load genetic map '{path}': rows are not sorted "
            "by chromosome. Most likely map shards were concatenated out of order. Sort "
            "the file by chromosome and position, or pass correctly ordered map shards."
        )
    for chrom, chrom_frame in out.groupby("CHR", sort=False):
        positions = chrom_frame["POS"].reset_index(drop=True)
        if not positions.equals(positions.sort_values(ignore_index=True)):
            raise LDSCInputError(
                f"{command} could not load genetic map '{path}': positions are not "
                f"sorted within chromosome {chrom}. Most likely this chromosome's map rows "
                "were shuffled. Sort the map by chromosome and position before rerunning."
            )

    return out.reset_index(drop=True)


def load_genetic_map_group(paths: Sequence[str | PathLike[str]], *, command: str) -> pd.DataFrame:
    """Load and combine one or more genetic map files into one canonical table; errors name ``command``."""

    if not paths:
        raise LDSCInputError(
            f"{command} could not load a genetic map group because no map paths "
            "were supplied. Most likely an empty path token reached the kernel. Pass at "
            "least one genetic map file or omit cM-window output."
        )
    frames = [load_genetic_map(path, command=command) for path in paths]
    combined = pd.concat(frames, axis=0, ignore_index=True)
    combined = combined.sort_values(by=["CHR", "POS"], key=lambda col: col.map(chrom_sort_key) if col.name == "CHR" else col, kind="mergesort")
    duplicated = combined.duplicated(subset=["CHR", "POS"], keep=False)
    if duplicated.any():
        raise LDSCInputError(
            f"{command} could not combine genetic map inputs: duplicate CHR/POS "
            "rows were found across map files. Most likely overlapping map shards were "
            "supplied. Use non-overlapping shards or a single deduplicated genetic map. "
            f"Other causes & fixes: {GENETIC_MAP_TROUBLESHOOTING[command]}"
        )
    return combined.reset_index(drop=True)


def interpolate_genetic_map_cm(
    chrom: str,
    positions: np.ndarray,
    genetic_map: pd.DataFrame,
    *,
    command: str,
) -> np.ndarray:
    """Interpolate cumulative cM values for `positions` on one chromosome; errors name ``command``."""

    chrom = normalize_chromosome(chrom)
    chrom_map = genetic_map.loc[genetic_map["CHR"] == chrom, ["POS", "CM"]].reset_index(drop=True)
    if len(chrom_map) == 0:
        raise LDSCInputError(
            f"{command} could not interpolate cM positions for chromosome {chrom}: "
            "the genetic map contains no rows for that chromosome. Most likely the wrong "
            "map build or incomplete chromosome shard was supplied. Pass a map that covers "
            f"chromosome {chrom}. Other causes & fixes: {GENETIC_MAP_TROUBLESHOOTING[command]}"
        )
    return np.interp(
        np.asarray(positions, dtype=np.int64),
        chrom_map["POS"].to_numpy(dtype=np.int64),
        chrom_map["CM"].to_numpy(dtype=float),
    )
