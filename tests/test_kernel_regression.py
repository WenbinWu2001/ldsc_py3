import math
from types import SimpleNamespace
import unittest
from unittest import mock

import numpy as np

from ldsc._kernel import regression as reg


class RGKernelTest(unittest.TestCase):
    def test_rg_converts_one_by_one_arrays_to_python_scalars(self):
        hsq1 = SimpleNamespace(
            tot=np.array([0.25]),
            intercept=np.array([1.0]),
            tot_delete_values=np.array([[0.20], [0.30]]),
        )
        hsq2 = SimpleNamespace(
            tot=np.array([1.0]),
            intercept=np.array([1.0]),
            tot_delete_values=np.array([[0.80], [1.20]]),
        )
        gencov = SimpleNamespace(
            tot=np.array([0.125]),
            tot_delete_values=np.array([[0.10], [0.15]]),
            N1=None,
            N2=None,
        )
        ratio_jackknife = SimpleNamespace(
            jknife_est=np.array([[0.49]]),
            jknife_se=np.array([[0.07]]),
        )

        with mock.patch.object(reg, "Hsq", side_effect=[hsq1, hsq2]), mock.patch.object(
            reg,
            "Gencov",
            return_value=gencov,
        ), mock.patch.object(
            reg.jk,
            "RatioJackknife",
            return_value=ratio_jackknife,
        ):
            result = reg.RG(
                np.ones((2, 1)),
                np.ones((2, 1)),
                np.ones((2, 1)),
                np.ones((2, 1)),
                np.ones((2, 1)),
                np.ones((2, 1)),
                np.array([[10.0]]),
            )

        self.assertIsInstance(result.rg_jknife, float)
        self.assertIsInstance(result.rg_se, float)
        self.assertIsInstance(result.rg_ratio, float)
        self.assertEqual(result.rg_jknife, 0.49)
        self.assertEqual(result.rg_se, 0.07)
        self.assertEqual(result.rg_ratio, 0.25)


class LiabilityConversionTest(unittest.TestCase):
    def test_known_value(self):
        # c(P=0.5, K=0.01) = K^2(1-K)^2 / (P(1-P) phi(isf(K))^2), verified via scipy.
        c = reg.liability_conversion_factor(0.5, 0.01)
        self.assertAlmostEqual(c, 0.5519073, places=6)

    def test_both_nan_is_one(self):
        self.assertEqual(reg.liability_conversion_factor(float("nan"), float("nan")), 1.0)

    def test_vectorizes_over_K(self):
        K = np.array([0.02, 0.05, 0.10, 0.20])
        c = reg.liability_conversion_factor(0.5, K)
        self.assertEqual(c.shape, K.shape)
        # Liability h2 = c * h2_obs is monotone increasing in K over this range.
        self.assertTrue(np.all(np.diff(c) > 0))

    def test_rejects_out_of_range(self):
        with self.assertRaises(Exception):
            reg.liability_conversion_factor(0.5, 1.0)
        with self.assertRaises(Exception):
            reg.liability_conversion_factor(0.0, 0.01)

    def test_rejects_half_specified(self):
        # One of (P, K) NaN while the other is finite is underspecified.
        with self.assertRaises(Exception):
            reg.liability_conversion_factor(float("nan"), 0.01)
        with self.assertRaises(Exception):
            reg.liability_conversion_factor(0.5, float("nan"))

    def test_h2_obs_to_liab_matches_primitive(self):
        self.assertAlmostEqual(
            reg.h2_obs_to_liab(0.3, 0.5, 0.01),
            0.3 * reg.liability_conversion_factor(0.5, 0.01),
            places=12,
        )

    def test_gencov_obs_to_liab_uses_sqrt_product(self):
        g = reg.gencov_obs_to_liab(0.1, 0.5, 0.5, 0.01, 0.02)
        expected = (
            0.1
            * math.sqrt(reg.liability_conversion_factor(0.5, 0.01))
            * math.sqrt(reg.liability_conversion_factor(0.5, 0.02))
        )
        self.assertAlmostEqual(g, expected, places=12)


if __name__ == "__main__":
    unittest.main()
