"""Interval-union projection onto one chromosome of a supplied annotation SNP grid."""

import numpy as np
import pandas as pd


class ChromosomeProjector:
    """Reuse a chromosome's sorted SNP positions across all interval queries."""

    def __init__(self, metadata):
        self.positions = pd.to_numeric(metadata['POS'], errors='raise').to_numpy(dtype=np.int64) - 1
        # Rows already in position order (the usual grid) need no permutation.
        if np.all(self.positions[1:] >= self.positions[:-1]):
            self.order, self.sorted_positions = None, self.positions
        else:
            self.order = np.argsort(self.positions, kind='stable')
            self.sorted_positions = self.positions[self.order]

    def project(self, intervals, *, padding_bp=0):
        """Return binary union membership for zero-based half-open intervals."""
        # Overlap depths never exceed one call's interval count, so int32
        # holds them; the running sum is taken in place.
        depth = np.zeros(len(self.positions) + 1, dtype=np.int32)
        for start, end in intervals:
            left = np.searchsorted(self.sorted_positions, max(0, int(start) - padding_bp), side='left')
            right = np.searchsorted(self.sorted_positions, int(end) + padding_bp, side='left')
            depth[left] += 1
            depth[right] -= 1
        np.cumsum(depth, out=depth)
        covered = depth[:-1] > 0
        if self.order is None:
            return covered
        result = np.empty(len(self.positions), dtype=bool)
        result[self.order] = covered
        return result

    def counts(self, starts, ends, *, padding_bp=0):
        """Count grid rows per individual gene interval, including shared SNPs."""
        return np.searchsorted(self.sorted_positions, np.asarray(ends) + padding_bp, side='left') - np.searchsorted(self.sorted_positions, np.maximum(0, np.asarray(starts) - padding_bp), side='left')
