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
    metadata = pd.DataFrame({"MAF": [0.05, 0.06, 0.2, 0.3]})  # row 0 (MAF==0.05) IS common (inclusive >=)
    A = annotations.to_numpy(dtype=np.float64)

    contribution = compute_overlap(metadata, annotations, n_baseline=2, common_maf_min=0.05)

    np.testing.assert_allclose(contribution.baseline_block_all, A[:, :2].T @ A)
    np.testing.assert_allclose(contribution.query_diagonal_all, (A[:, 2:] ** 2).sum(axis=0))
    assert contribution.n_all == 4

    common = np.array([True, True, True, True])
    Ac = A[common]
    np.testing.assert_allclose(contribution.baseline_block_common, Ac[:, :2].T @ Ac)
    np.testing.assert_allclose(contribution.query_diagonal_common, (Ac[:, 2:] ** 2).sum(axis=0))
    assert contribution.n_common == 4


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


def test_overlap_long_frame_round_trip():
    from ldsc._kernel.overlap import OverlapContribution
    from ldsc.overlap_matrix import (
        LDScoreOverlap, overlap_to_long_frame, overlap_from_long_frame,
    )
    contribution = OverlapContribution(
        baseline_block_all=np.array([[10.0, 4.0, 3.0]]),       # base x [base, q, cont]
        baseline_block_common=np.array([[8.0, 3.0, 2.0]]),
        query_diagonal_all=np.array([6.0, 9.0]),                # q, cont
        query_diagonal_common=np.array([5.0, 7.0]),
        n_all=10, n_common=8,
    )
    overlap = LDScoreOverlap.from_contribution(
        contribution, baseline_columns=["base"], query_columns=["q", "cont"]
    )
    frame = overlap_to_long_frame(overlap)
    assert set(frame.columns) == {"row_annotation", "col_annotation", "overlap_all_snps", "overlap_common_snps"}
    restored = overlap_from_long_frame(
        frame, baseline_columns=["base"], query_columns=["q", "cont"],
        total_all_reference_snps=10.0, total_common_reference_snps=8.0,
    )
    pd.testing.assert_frame_equal(restored.baseline_block_all, overlap.baseline_block_all)
    pd.testing.assert_frame_equal(restored.baseline_block_common, overlap.baseline_block_common)
    pd.testing.assert_series_equal(restored.query_diagonal_all, overlap.query_diagonal_all)
    pd.testing.assert_series_equal(restored.query_diagonal_common, overlap.query_diagonal_common)
    assert restored.total_common_reference_snps == 8.0


def _legacy_overlap_reference(prop, prop_cov, coef, coef_cov, n_blocks, O, M, M_tot):
    """Transcription of legacy ``Hsq._overlap_output`` math (pure-NumPy oracle)."""
    from scipy.stats import t as tdist
    K = O.shape[0]
    P = np.zeros((K, K))
    for i in range(K):
        P[i, :] = O[i, :] / M
    prop_h2 = (P @ prop.reshape(-1)).reshape(K)
    prop_h2_se = np.sqrt(np.maximum(0, np.diag(P @ prop_cov @ P.T)))
    prop_snps = (M / M_tot).reshape(K)
    enrich = prop_h2 / prop_snps
    enrich_se = prop_h2_se / prop_snps
    D = np.zeros((K, K))
    for i in range(K):
        if M_tot != M[i]:
            D[i, :] = O[i, :] / M[i] - (M - O[i, :]) / (M_tot - M[i])
    diff_est = D @ coef
    diff_se = np.sqrt(np.diag(D @ coef_cov @ D.T))
    enr_p = [np.nan if diff_se[i] == 0 else 2 * tdist.sf(abs(diff_est[i] / diff_se[i]), n_blocks) for i in range(K)]
    return prop_snps, prop_h2, prop_h2_se, enrich, enrich_se, enr_p


def test_overlap_aware_category_table_matches_legacy_and_augments():
    from types import SimpleNamespace
    from scipy import stats
    from ldsc.overlap_matrix import overlap_aware_category_table

    names = ["base", "catA", "q"]
    M = np.array([10.0, 4.0, 6.0])
    M_tot = 10.0
    # base overlaps everything fully; catA & q share 2 SNPs.
    O = np.array([[10.0, 4.0, 6.0],
                  [4.0, 4.0, 2.0],
                  [6.0, 2.0, 6.0]])
    coef = np.array([1e-7, 3e-7, 5e-7])
    coef_cov = np.diag([1e-15, 2e-15, 3e-15])
    prop = np.array([[0.5, 0.2, 0.3]])
    prop_cov = np.diag([1e-3, 2e-3, 3e-3])
    cat = M * coef
    cat_se = np.array([1e-7, 1e-7, 1e-7])
    hsq = SimpleNamespace(
        n_annot=3, prop=prop, prop_cov=prop_cov, coef=coef, coef_cov=coef_cov,
        coef_se=np.sqrt(np.diag(coef_cov)), n_blocks=200, cat=cat, cat_se=cat_se,
    )
    table = overlap_aware_category_table(hsq, O, M, M_tot, names)

    ps, ph, phse, en, ense, enp = _legacy_overlap_reference(prop, prop_cov, coef, coef_cov, 200, O, M, M_tot)
    np.testing.assert_allclose(table["Prop._SNPs"], ps)
    np.testing.assert_allclose(table["Prop._h2"], ph)
    np.testing.assert_allclose(table["Prop._h2_std_error"], phse)
    np.testing.assert_allclose(table["Enrichment"], en)
    np.testing.assert_allclose(table["Enrichment_std_error"], ense)
    np.testing.assert_allclose(np.asarray(table["Enrichment_p"], float), np.asarray(enp, float))
    assert np.isnan(float(table.loc[0, "Enrichment_p"]))  # base contains all SNPs
    np.testing.assert_allclose(table["Coefficient_p"], stats.norm.sf(coef / np.sqrt(np.diag(coef_cov))))
    np.testing.assert_allclose(table["Category_h2"], cat)
    assert table["overlap_aware"].all()
    assert "Coefficient_z" in table.columns and "Coefficient_z-score" not in table.columns


def test_overlap_aware_hand_example_disjoint_plus_base():
    from types import SimpleNamespace
    from ldsc.overlap_matrix import overlap_aware_category_table
    # base = all 6 SNPs; A = first 3; Bcat = last 3 (disjoint A,B but both overlap base)
    names = ["base", "A", "Bcat"]
    M = np.array([6.0, 3.0, 3.0]); M_tot = 6.0
    O = np.array([[6.0, 3.0, 3.0],
                  [3.0, 3.0, 0.0],
                  [3.0, 0.0, 3.0]])
    coef = np.array([0.0, 1.0, 2.0]); cat = M * coef
    hsq = SimpleNamespace(
        n_annot=3, prop=np.array([[0.0, 1 / 3, 2 / 3]]), prop_cov=np.eye(3) * 1e-6,
        coef=coef, coef_cov=np.eye(3) * 1e-6, coef_se=np.full(3, 1e-3),
        n_blocks=200, cat=cat, cat_se=np.full(3, 1e-3),
    )
    table = overlap_aware_category_table(hsq, O, M, M_tot, names).set_index("Category")
    # h2_tot = sum M*coef = 3*1 + 3*2 = 9 ; Prop._h2[A] = (3*0 + 3*1 + 0*2)/9 = 3/9
    assert abs(float(table.loc["A", "Prop._h2"]) - (3 / 9)) < 1e-9
    assert abs(float(table.loc["Bcat", "Prop._h2"]) - (6 / 9)) < 1e-9
    assert abs(float(table.loc["base", "Prop._h2"]) - 1.0) < 1e-9
    assert abs(float(table.loc["A", "Prop._SNPs"]) - 0.5) < 1e-12
    assert abs(float(table.loc["A", "Enrichment"]) - (3 / 9) / 0.5) < 1e-9


def test_assemble_model_overlap_reconstructs_submatrix():
    from ldsc._kernel.overlap import OverlapContribution
    from ldsc.overlap_matrix import LDScoreOverlap, assemble_model_overlap
    # Full A: baseline=[base, catA], query=[q1, q2]
    A = np.array([
        [1, 1, 1, 0],
        [1, 1, 0, 1],
        [1, 0, 1, 0],
        [1, 0, 0, 1],
    ], dtype=float)
    O_full = A.T @ A
    contribution = OverlapContribution(
        baseline_block_all=A[:, :2].T @ A,
        baseline_block_common=None,
        query_diagonal_all=(A[:, 2:] ** 2).sum(axis=0),
        query_diagonal_common=None,
        n_all=4, n_common=None,
    )
    overlap = LDScoreOverlap.from_contribution(contribution, ["base", "catA"], ["q1", "q2"])
    # model for q1: retained = [base, catA, q1] -> rows/cols 0,1,2 of O_full
    O_R = assemble_model_overlap(overlap, ["base", "catA", "q1"], use_common=False)
    np.testing.assert_allclose(O_R, O_full[np.ix_([0, 1, 2], [0, 1, 2])])
    np.testing.assert_allclose(O_R[2, :], O_R[:, 2])  # assembled query row/col symmetric


def test_model_collinearity_error_names_worst_pair():
    from ldsc.overlap_matrix import model_collinearity_error
    # near-duplicate LD-score columns -> high condition number
    X = np.array([[1.0, 1.0001, 2.0], [2.0, 2.0001, 1.0], [3.0, 3.0001, 0.5], [4.0, 4.0002, 0.1]])
    columns = ["catA", "catA_dup", "other"]
    O = np.array([[10.0, 9.99, 1.0], [9.99, 10.0, 1.0], [1.0, 1.0, 10.0]])
    msg = model_collinearity_error(X, columns, O, threshold=1e5)
    assert msg is not None and "catA" in msg and "catA_dup" in msg


def test_model_collinearity_error_none_when_well_conditioned():
    from ldsc.overlap_matrix import model_collinearity_error
    X = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    assert model_collinearity_error(X, ["a", "b"], np.eye(2), threshold=1e5) is None
