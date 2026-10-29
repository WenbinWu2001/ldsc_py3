"""Exact query-column variation across independently computed chromosomes."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ScoreVariation:
    """First values and within-chromosome flags in the LD-score producer's precision."""

    first: np.ndarray | None
    varying: np.ndarray


def summarize_variation(columns):
    """Summarize 1-D score columns without copying their complete matrix.

    NaNs count as one distinct value, matching ``Series.nunique(dropna=False)``.
    Empty chromosomes have no first values. Input dtype is preserved so direct
    float32 and indexed float64 scores retain their own filtering boundaries.
    """
    first, varying = [], []
    for column in columns:
        if not len(column):
            varying.append(False)
            continue
        value = column[0]
        first.append(value)
        varying.append(not np.all((column == value) | (np.isnan(column) & np.isnan(value))))
    return ScoreVariation(np.asarray(first) if first else None, np.asarray(varying, dtype=bool))


def combine_variation(summaries, width):
    """Detect global variation, including constants that differ by chromosome."""
    first = None
    varying = np.zeros(width, dtype=bool)
    for summary in summaries:
        if summary.first is None:
            continue
        if first is None:
            first = summary.first
        same = (summary.first == first) | (np.isnan(summary.first) & np.isnan(first))
        varying |= summary.varying | ~same
    return varying
