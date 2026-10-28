"""Interval-union projection onto one chromosome of a supplied annotation SNP grid."""

import numpy as np
import pandas as pd


class ChromosomeProjector:
    """Reuse a chromosome's sorted SNP positions across all interval queries."""

    def __init__(self, metadata):
        self.positions = pd.to_numeric(metadata['POS'], errors='raise').to_numpy(dtype=np.int64) - 1
        self.order = np.argsort(self.positions, kind='stable')
        self.sorted_positions = self.positions[self.order]

    def project(self, intervals, *, padding_bp=0):
        """Return binary union membership for zero-based half-open intervals."""
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
