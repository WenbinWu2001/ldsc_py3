from types import SimpleNamespace

import numpy as np

from ldsc.outputs import H2_REGRESSION_BIN_COLUMNS
from ldsc.regression_runner import summarize_ld_score_regression_bins


def test_summarize_ld_score_regression_bins_uses_final_fit_values():
    hsq = SimpleNamespace(
        coef=np.array([0.0002]),
        tot=np.array([0.2]),
        intercept=np.array([1.0]),
    )

    result = summarize_ld_score_regression_bins(
        hsq,
        ld_score=np.array([4.0, 1.0, 3.0, 2.0]),
        chi_square=np.array([2.1, 1.1, 1.6, 1.4]),
        sample_size=np.full(4, 100.0),
        regression_ld_score=np.ones(4),
        reference_snp_count=1000.0,
        max_bins=2,
    )

    assert result.columns.tolist() == H2_REGRESSION_BIN_COLUMNS
    np.testing.assert_array_equal(result["bin"], [1, 2])
    np.testing.assert_array_equal(result["n_snps"], [2, 2])
    np.testing.assert_allclose(result["mean_ld_score"], [1.5, 3.5])
    np.testing.assert_allclose(result["mean_chi_square"], [1.25, 1.85])
    np.testing.assert_allclose(result["sd_chi_square"], [0.2121320343559642, 0.3535533905932738])
    np.testing.assert_allclose(result["mean_fitted_chi_square"], [1.03, 1.07])
    np.testing.assert_allclose(result["mean_regression_weight"], [0.4714312485639342, 0.43683381507900576])


def test_summarize_bins_preserves_fitted_row_order_for_tied_ld_scores():
    hsq = SimpleNamespace(coef=np.array([0.0002]), tot=np.array([0.2]), intercept=np.array([1.0]))

    result = summarize_ld_score_regression_bins(
        hsq,
        ld_score=np.array([1.0, 1.0, 1.0, 2.0]),
        chi_square=np.array([1.0, 2.0, 100.0, 4.0]),
        sample_size=np.full(4, 100.0),
        regression_ld_score=np.ones(4),
        reference_snp_count=1000.0,
        max_bins=2,
    )

    np.testing.assert_allclose(result["mean_chi_square"], [1.5, 52.0])
    np.testing.assert_array_equal(result["n_snps"], [2, 2])


def test_summarize_bins_uses_one_bin_per_snp_below_cap_and_fifty_at_cap():
    hsq = SimpleNamespace(coef=np.array([0.0002]), tot=np.array([0.2]), intercept=np.array([1.0]))

    small = summarize_ld_score_regression_bins(
        hsq,
        ld_score=np.arange(1.0, 4.0),
        chi_square=np.arange(1.0, 4.0),
        sample_size=np.full(3, 100.0),
        regression_ld_score=np.ones(3),
        reference_snp_count=1000.0,
    )
    large = summarize_ld_score_regression_bins(
        hsq,
        ld_score=np.arange(1.0, 102.0),
        chi_square=np.arange(1.0, 102.0),
        sample_size=np.linspace(100.0, 200.0, 101),
        regression_ld_score=np.ones(101),
        reference_snp_count=1000.0,
    )

    assert len(small) == 3
    assert len(large) == 50
    assert large["n_snps"].sum() == 101
    assert large["n_snps"].max() - large["n_snps"].min() <= 1
