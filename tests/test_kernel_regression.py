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


if __name__ == "__main__":
    unittest.main()
