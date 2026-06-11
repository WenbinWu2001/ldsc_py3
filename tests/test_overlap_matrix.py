"""Tests for the annotation overlap matrix: computation, storage, and assembly."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest


def test_compute_overlap_blocks_binary_and_continuous():
    from ldsc._kernel.overlap import compute_overlap
    # 4 SNPs, baseline = [base(all ones), catA], query = [q, cont]
    annotations = pd.DataFrame({
        "base": [1.0, 1.0, 1.0, 1.0],
        "catA": [1.0, 1.0, 0.0, 0.0],
        "q":    [1.0, 0.0, 1.0, 0.0],
        "cont": [0.5, 2.0, 1.0, 0.0],
    })
    metadata = pd.DataFrame({"MAF": [0.05, 0.06, 0.2, 0.3]})  # row 0 (MAF==0.05) is NOT common
    A = annotations.to_numpy(dtype=np.float64)

    contribution = compute_overlap(metadata, annotations, n_baseline=2, common_maf_min=0.05)

    np.testing.assert_allclose(contribution.baseline_block_all, A[:, :2].T @ A)
    np.testing.assert_allclose(contribution.query_diagonal_all, (A[:, 2:] ** 2).sum(axis=0))
    assert contribution.n_all == 4

    common = np.array([False, True, True, True])
    Ac = A[common]
    np.testing.assert_allclose(contribution.baseline_block_common, Ac[:, :2].T @ Ac)
    np.testing.assert_allclose(contribution.query_diagonal_common, (Ac[:, 2:] ** 2).sum(axis=0))
    assert contribution.n_common == 3


def test_compute_overlap_without_maf_has_no_common_universe():
    from ldsc._kernel.overlap import compute_overlap
    annotations = pd.DataFrame({"base": [1.0, 1.0], "q": [1.0, 0.0]})
    metadata = pd.DataFrame({"CHR": ["1", "1"]})  # no MAF
    contribution = compute_overlap(metadata, annotations, n_baseline=1, common_maf_min=0.05)
    assert contribution.baseline_block_common is None
    assert contribution.query_diagonal_common is None
    assert contribution.n_common is None


def test_sum_overlap_contributions_adds_blocks_and_sizes():
    from ldsc._kernel.overlap import OverlapContribution, sum_overlap_contributions
    a = OverlapContribution(np.ones((1, 2)), np.ones((1, 2)), np.array([2.0]), np.array([1.0]), 5, 3)
    b = OverlapContribution(np.full((1, 2), 3.0), np.full((1, 2), 2.0), np.array([4.0]), np.array([2.0]), 7, 4)
    total = sum_overlap_contributions([a, b])
    np.testing.assert_allclose(total.baseline_block_all, [[4.0, 4.0]])
    np.testing.assert_allclose(total.baseline_block_common, [[3.0, 3.0]])
    np.testing.assert_allclose(total.query_diagonal_all, [6.0])
    np.testing.assert_allclose(total.query_diagonal_common, [3.0])
    assert total.n_all == 12 and total.n_common == 7


def test_sum_overlap_contributions_drops_common_if_any_chrom_lacks_it():
    from ldsc._kernel.overlap import OverlapContribution, sum_overlap_contributions
    a = OverlapContribution(np.ones((1, 1)), np.ones((1, 1)), np.array([1.0]), np.array([1.0]), 3, 2)
    b = OverlapContribution(np.ones((1, 1)), None, np.array([1.0]), None, 2, None)
    total = sum_overlap_contributions([a, b])
    assert total.baseline_block_common is None and total.n_common is None
