"""Tests for the reference-panel R² pair-query API and R²→r converter."""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ldsc.r2_query import unbiased_r2_to_pearson_r


class TestUnbiasedR2ToPearsonR:
    def test_inverts_unbiased_correction_then_roots(self):
        # Build adjusted R2 from a known raw R2, then recover signed r.
        n = 1000
        r2_raw = 0.49
        r2_adj = r2_raw - (1.0 - r2_raw) / (n - 2)
        r = unbiased_r2_to_pearson_r(r2_adj, n, sign=-1)
        assert r == pytest.approx(-0.7, abs=1e-3)

    def test_negative_adjusted_value_still_yields_real_r(self):
        # Near-zero LD: adjusted R2 can be slightly negative; r must be real (not NaN).
        r = unbiased_r2_to_pearson_r(-0.0005, 1000, sign=1)
        assert np.isfinite(r)
        assert r >= 0.0

    def test_sign_none_returns_magnitude(self):
        r = unbiased_r2_to_pearson_r(0.25, 1000, sign=None)
        assert r >= 0.0

    def test_vectorized_arrays(self):
        r2 = np.array([0.25, 0.04], dtype=np.float64)
        sign = np.array([1, -1], dtype=np.int8)
        out = unbiased_r2_to_pearson_r(r2, 1000, sign=sign)
        assert out.shape == (2,)
        assert out[0] > 0 and out[1] < 0
