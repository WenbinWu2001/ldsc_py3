"""Interval-union projection onto a supplied chromosome SNP grid."""

from typing import Sequence
import numpy as np
import pandas as pd
from .chromosome_inference import normalize_chromosome


class ChromosomeProjector:
    """Reuse a chromosome's sorted SNP positions across all interval queries."""

    def __init__(self, metadata):
        self.positions = pd.to_numeric(metadata['POS'], errors='raise').to_numpy(dtype=np.int64) - 1
        self.order = np.argsort(self.positions, kind='stable')
        self.sorted_positions = self.positions[self.order]

    def project(self, intervals, *, padding_bp=0):
        """Return binary union membership for 0-based half-open intervals."""
        difference = np.zeros(len(self.positions) + 1, dtype=np.int64)
        for start, end in intervals:
            left = np.searchsorted(self.sorted_positions, max(0, int(start) - padding_bp), side='left')
            right = np.searchsorted(self.sorted_positions, int(end) + padding_bp, side='left')
            difference[left] += 1
            difference[right] -= 1
        result = np.empty(len(self.positions), dtype=bool)
        result[self.order] = np.cumsum(difference[:-1]) > 0
        return result

    def counts(self, starts, ends, *, padding_bp=0):
        """Count grid rows per individual gene interval, including shared SNPs."""
        return np.searchsorted(self.sorted_positions, np.asarray(ends) + padding_bp, side='left') - np.searchsorted(self.sorted_positions, np.maximum(0, np.asarray(starts) - padding_bp), side='left')


def _project_intervals_to_metadata(
    metadata: pd.DataFrame,
    intervals: Sequence[tuple[str, int, int]],
    *,
    padding_bp: int,
) -> np.ndarray:
    """Project a union of 0-based half-open intervals onto 1-based SNP positions.

    Intervals are padded, sorted, and merged per chromosome. SNP positions are
    sorted only within each chromosome, then interval membership is assigned by
    binary-search bounds. Runtime therefore scales with SNP and interval counts
    rather than their product.
    """
    chrom = metadata["CHR"].map(normalize_chromosome).astype(str).to_numpy()
    pos0 = pd.to_numeric(metadata["POS"], errors="raise").astype(np.int64).to_numpy() - 1
    values = np.zeros(len(metadata), dtype=bool)
    intervals_by_chrom: dict[str, list[tuple[int, int]]] = {}
    for interval_chrom, start0, end in intervals:
        padded_start = max(0, int(start0) - padding_bp)
        padded_end = int(end) + padding_bp
        normalized_chrom = normalize_chromosome(interval_chrom)
        intervals_by_chrom.setdefault(normalized_chrom, []).append((padded_start, padded_end))
    for interval_chrom, chrom_intervals in intervals_by_chrom.items():
        row_indices = np.flatnonzero(chrom == interval_chrom)
        if not len(row_indices):
            continue
        order = np.argsort(pos0[row_indices], kind="stable")
        sorted_indices = row_indices[order]
        sorted_pos = pos0[sorted_indices]
        difference = np.zeros(len(sorted_pos) + 1, dtype=np.int32)
        for start, end in _merge_intervals(chrom_intervals):
            left = int(np.searchsorted(sorted_pos, start, side="left"))
            right = int(np.searchsorted(sorted_pos, end, side="left"))
            difference[left] += 1
            difference[right] -= 1
        values[sorted_indices[np.cumsum(difference[:-1]) > 0]] = True
    return values


def _merge_intervals(intervals: Sequence[tuple[int, int]]) -> list[tuple[int, int]]:
    """Return a sorted union of half-open intervals, joining touching bounds."""
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged
