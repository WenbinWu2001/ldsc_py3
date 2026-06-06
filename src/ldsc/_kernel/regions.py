"""Region-exclusion interval loading and SNP keep-masking (private kernel).

Intervals use the standard BED convention: 0-based, half-open ``[start, end)``.
SNP positions in panel metadata are 1-based (PLINK convention). A 1-based
position ``p`` is excluded iff some interval contains it, i.e.
``start < p <= end``. See ``region_exclusion_keep_mask``.
"""
from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from os import PathLike
from typing import Sequence

import numpy as np
import pandas as pd

from ..chromosome_inference import normalize_chromosome
from ..errors import LDSCConfigError, LDSCInputError, LDSCInternalError, LDSCUsageError

REGION_PRESETS: frozenset[str] = frozenset({"mhc", "centromeres"})
_PRESET_BUILDS: frozenset[str] = frozenset({"hg19", "hg38"})


@dataclass(frozen=True)
class RegionIntervals:
    """Per-chromosome exclusion intervals in 0-based half-open BED coordinates.

    Parameters
    ----------
    intervals : dict of str to numpy.ndarray
        Maps normalized chromosome token to a sorted, disjoint ``(k, 2)`` int64
        array of ``[start, end)`` rows. An empty mapping excludes nothing.
    source_labels : tuple of str
        Human-readable provenance for logging, e.g. ``("preset:mhc[hg19]",)``.
    """

    intervals: dict[str, np.ndarray]
    source_labels: tuple[str, ...]


def _coalesce(rows: list[tuple[int, int]]) -> np.ndarray:
    """Sort and union possibly-overlapping intervals into disjoint sorted rows."""
    if not rows:
        return np.empty((0, 2), dtype=np.int64)
    merged: list[list[int]] = []
    for start, end in sorted(rows):
        # `<=` (not `<`) so adjacent half-open intervals like [10,20) and [20,30)
        # coalesce into one contiguous block under the `start < p <= end` rule.
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return np.asarray(merged, dtype=np.int64).reshape(-1, 2)


def _intervals_from_rows(rows: Sequence[tuple[str, int, int]]) -> dict[str, np.ndarray]:
    """Group ``(chrom, start, end)`` rows by normalized chromosome and coalesce."""
    grouped: dict[str, list[tuple[int, int]]] = {}
    for chrom, start, end in rows:
        grouped.setdefault(normalize_chromosome(chrom), []).append((int(start), int(end)))
    return {chrom: _coalesce(pairs) for chrom, pairs in grouped.items()}


def _parse_bed_text(text: str, label: str) -> list[tuple[str, int, int]]:
    """Parse standard BED text into ``(chrom, start, end)`` rows."""
    rows: list[tuple[str, int, int]] = []
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith(("#", "track", "browser")):
            continue
        fields = line.split("\t") if "\t" in line else line.split()
        if len(fields) < 3:
            raise LDSCInputError(
                f"Region BED parsing failed in '{label}' line {lineno}: expected at least 3 "
                f"columns (chrom start end), found {len(fields)}. Most likely the file is not "
                "tab/space-delimited BED. Provide a standard BED file."
            )
        chrom, start_token, end_token = fields[0], fields[1], fields[2]
        try:
            start, end = int(start_token), int(end_token)
        except ValueError as exc:
            raise LDSCInputError(
                f"Region BED parsing failed in '{label}' line {lineno}: start/end "
                f"('{start_token}', '{end_token}') are not integers. Provide 0-based "
                "half-open integer BED coordinates."
            ) from exc
        if start >= end:
            raise LDSCInputError(
                f"Region BED parsing failed in '{label}' line {lineno}: start {start} is not "
                f"less than end {end}. BED intervals are 0-based half-open with start < end."
            )
        rows.append((chrom, start, end))
    return rows


def load_preset_intervals(names: Sequence[str], build: str) -> RegionIntervals:
    """Read packaged preset BEDs for ``names`` in genome ``build``.

    Parameters
    ----------
    names : sequence of str
        Preset names; each must be in :data:`REGION_PRESETS`.
    build : str
        Genome build, either ``"hg19"`` or ``"hg38"``.

    Returns
    -------
    RegionIntervals
        Union of the requested presets' intervals.
    """
    if build not in _PRESET_BUILDS:
        raise LDSCConfigError(
            f"Region preset loading received unsupported genome build '{build}'. Most likely "
            "an invalid --exclude-regions-build was passed. Use 'hg19' or 'hg38'."
        )
    unknown = sorted(set(names) - REGION_PRESETS)
    if unknown:
        raise LDSCUsageError(
            f"Region preset loading received unknown region preset(s): {unknown}. Most likely "
            f"a name was misspelled. Valid presets are {sorted(REGION_PRESETS)}."
        )
    rows: list[tuple[str, int, int]] = []
    labels: list[str] = []
    for name in names:
        path = resources.files("ldsc").joinpath("data", "regions", f"{name}.{build}.bed")
        try:
            text = path.read_text(encoding="utf-8")
        except (FileNotFoundError, OSError) as exc:
            raise LDSCInternalError(
                f"Region preset loading could not read packaged BED for '{name}.{build}'. "
                "Most likely the package data was not installed. Reinstall the package or "
                "report the traceback."
            ) from exc
        rows.extend(_parse_bed_text(text, label=f"{name}.{build}.bed"))
        labels.append(f"preset:{name}[{build}]")
    return RegionIntervals(intervals=_intervals_from_rows(rows), source_labels=tuple(labels))


def load_bed_intervals(paths: Sequence[str | PathLike[str]]) -> RegionIntervals:
    """Read one or more user BED files, applied as-is on panel CHR/POS.

    Parameters
    ----------
    paths : sequence of str or os.PathLike[str]
        Standard BED files (0-based half-open). Chromosome tokens are
        normalized (``chr6`` -> ``6``); overlapping intervals are coalesced.

    Returns
    -------
    RegionIntervals
        Union of all intervals across the supplied files, labeled ``bed:<path>``.
    """
    rows: list[tuple[str, int, int]] = []
    labels: list[str] = []
    for path in paths:
        path_str = str(path)
        try:
            with open(path_str, encoding="utf-8") as handle:
                text = handle.read()
        except (FileNotFoundError, OSError) as exc:
            raise LDSCInputError(
                f"Region BED loading could not read user BED file '{path_str}'. Most likely "
                "the path is wrong or unreadable. Provide an existing standard BED file."
            ) from exc
        rows.extend(_parse_bed_text(text, label=path_str))
        labels.append(f"bed:{path_str}")
    return RegionIntervals(intervals=_intervals_from_rows(rows), source_labels=tuple(labels))


def merge_intervals(*groups: RegionIntervals) -> RegionIntervals:
    """Union several RegionIntervals into one (coalescing overlaps per chromosome)."""
    rows: list[tuple[str, int, int]] = []
    labels: list[str] = []
    for group in groups:
        for chrom, bounds in group.intervals.items():
            for start, end in bounds.tolist():
                rows.append((chrom, int(start), int(end)))
        labels.extend(group.source_labels)
    return RegionIntervals(intervals=_intervals_from_rows(rows), source_labels=tuple(labels))


def region_exclusion_keep_mask(
    metadata: pd.DataFrame,
    intervals: RegionIntervals,
    *,
    chr_col: str = "CHR",
    pos_col: str = "POS",
) -> np.ndarray:
    """Return a boolean KEEP mask over ``metadata`` rows.

    A 1-based SNP position ``p`` is excluded iff it falls within a BED interval:
    ``start < p <= end``, equivalently the 0-based position ``p - 1`` satisfies
    ``start <= p - 1 < end``. Intervals are assumed sorted and disjoint per
    chromosome (guaranteed by the loaders). Rows whose chromosome has no
    intervals are kept. Returns all-True when ``intervals`` is empty.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Must contain ``chr_col`` and ``pos_col``. Positions are 1-based.
    intervals : RegionIntervals
        Exclusion intervals to apply.
    chr_col, pos_col : str, optional
        Column names for chromosome and 1-based position. Defaults ``"CHR"`` /
        ``"POS"``.

    Returns
    -------
    numpy.ndarray
        Boolean array, ``True`` to keep the row.
    """
    n = len(metadata)
    keep = np.ones(n, dtype=bool)
    if not intervals.intervals or n == 0:
        return keep
    chrom_tokens = metadata[chr_col].astype(str).map(normalize_chromosome).to_numpy()
    pos = pd.to_numeric(metadata[pos_col], errors="raise").to_numpy()
    for chrom, bounds in intervals.intervals.items():
        sel = chrom_tokens == chrom
        if not sel.any():
            continue
        p = pos[sel]
        starts = bounds[:, 0]
        ends = bounds[:, 1]
        # Candidate interval for each p: the last one with start < p.
        idx = np.searchsorted(starts, p, side="left") - 1
        excluded = np.zeros(p.shape, dtype=bool)
        valid = idx >= 0
        excluded[valid] = ends[idx[valid]] >= p[valid]
        local = keep[sel]
        local[excluded] = False
        keep[sel] = local
    return keep
