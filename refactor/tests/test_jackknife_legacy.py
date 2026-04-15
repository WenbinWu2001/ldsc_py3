from pathlib import Path
import sys
import unittest

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc._kernel import _jackknife as jk


class TestJackknife(unittest.TestCase):
    def test_separators(self):
        n = 20
        x = np.arange(n)
        for i in range(2, int(np.floor(n / 2))):
            s = jk.Jackknife.get_separators(n, i)
            lengths = [len(x[s[j] : s[j + 1]]) for j in range(len(s) - 2)]
            self.assertLessEqual(max(lengths) - min(lengths), 1)

    def test_jknife_1d(self):
        pseudovalues = np.atleast_2d(np.arange(10)).T
        est, var, se, cov = jk.Jackknife.jknife(pseudovalues)
        assert_array_almost_equal(var, [[0.91666667]])
        assert_array_almost_equal(est, [[4.5]])
        assert_array_almost_equal(cov, [[0.91666667]])
        assert_array_almost_equal(se**2, var)
        self.assertFalse(np.any(np.isnan(cov)))
        assert_array_equal(cov.shape, (1, 1))
        assert_array_equal(var.shape, (1, 1))
        assert_array_equal(est.shape, (1, 1))
        assert_array_equal(se.shape, (1, 1))

    def test_jknife_2d(self):
        pseudovalues = np.vstack([np.arange(10), np.arange(10)]).T
        est, var, se, cov = jk.Jackknife.jknife(pseudovalues)
        assert_array_almost_equal(var, np.array([[0.91666667, 0.91666667]]))
        assert_array_almost_equal(est, np.array([[4.5, 4.5]]))
        assert_array_almost_equal(
            cov,
            np.matrix([[0.91666667, 0.91666667], [0.91666667, 0.91666667]]),
        )
        assert_array_almost_equal(se**2, var)
        assert_array_equal(cov.shape, (2, 2))
        assert_array_equal(var.shape, (1, 2))
        assert_array_equal(est.shape, (1, 2))
        assert_array_equal(se.shape, (1, 2))

    def test_delete_to_pseudo(self):
        for dim in [1, 2]:
            est = np.ones((1, dim))
            delete_values = np.ones((20, dim))
            x = jk.Jackknife.delete_values_to_pseudovalues(delete_values, est)
            assert_array_equal(x, np.ones_like(delete_values))

        with self.assertRaises(ValueError):
            jk.Jackknife.delete_values_to_pseudovalues(delete_values, est.T)


class TestLstsqJackknifeSlow(unittest.TestCase):
    def test_delete_values_1d(self):
        func = lambda x, y: np.atleast_2d(np.sum(x + y))
        s = [0, 5, 10]
        x = np.atleast_2d(np.arange(10)).T
        y = np.atleast_2d(np.arange(10)).T
        p = jk.LstsqJackknifeSlow.delete_values(x, y, func, s)
        assert_array_equal(p.shape, (2, 1))
        assert_array_almost_equal(p, [[70], [20]])

    def test_delete_values_2d_1(self):
        func = lambda x, y: np.atleast_2d(np.sum(x + y, axis=0))
        s = [0, 2, 4, 6, 8, 10]
        x = np.vstack([np.arange(10), 2 * np.arange(10)]).T
        y = np.atleast_2d(np.arange(10)).T
        p = jk.LstsqJackknifeSlow.delete_values(x, y, func, s)
        assert_array_equal(p.shape, (5, 2))
        assert_array_almost_equal(p, [[88, 132], [80, 120], [72, 108], [64, 96], [56, 84]])

    def test_delete_values_2d_2(self):
        func = lambda x, y: np.atleast_2d(np.sum(x + y, axis=0))
        s = [0, 5, 10]
        x = np.vstack([np.arange(10), 2 * np.arange(10), 3 * np.arange(10)]).T
        y = np.atleast_2d(np.arange(10)).T
        p = jk.LstsqJackknifeSlow.delete_values(x, y, func, s)
        assert_array_equal(p.shape, (2, 3))
        assert_array_almost_equal(p, [[70, 105, 140], [20, 30, 40]])

    def test_lstsqjackknifeslow(self):
        x = np.atleast_2d(np.arange(10)).T
        y = np.atleast_2d(2 * np.arange(10)).T
        reg = jk.LstsqJackknifeSlow(x, y, n_blocks=10)
        regnn = jk.LstsqJackknifeSlow(x, y, n_blocks=10, nn=True)
        assert_array_almost_equal(reg.est, [[2.0]])
        assert_array_almost_equal(regnn.est, [[2.0]])

    def test_bad_data(self):
        x = np.arange(10)
        y = 2 * np.arange(9)
        with self.assertRaises(ValueError):
            jk.LstsqJackknifeSlow(x, y, 10)

    def test_too_many_blocks(self):
        x = np.arange(10)
        y = 2 * np.arange(10)
        with self.assertRaises(ValueError):
            jk.LstsqJackknifeSlow(x, y, 11)


class TestLstsqJackknifeFast(unittest.TestCase):
    def test_block_values_1d(self):
        x = np.arange(6).reshape((6, 1))
        y = np.arange(6).reshape((6, 1))
        s = [0, 2, 4, 6]
        xty, xtx = jk.LstsqJackknifeFast.block_values(x, y, s)
        assert_array_equal(xty.shape, (3, 1))
        assert_array_equal(xtx.shape, (3, 1, 1))
        assert_array_almost_equal(xty, [[1], [13], [41]])
        assert_array_almost_equal(xtx, [[[1]], [[13]], [[41]]])

    def test_block_values_2d(self):
        x = np.vstack([np.arange(6), 2 * np.arange(6)]).T
        y = np.arange(6).reshape((6, 1))
        s = [0, 2, 4, 6]
        xty, xtx = jk.LstsqJackknifeFast.block_values(x, y, s)
        assert_array_equal(xty.shape, (3, 2))
        assert_array_equal(xtx.shape, (3, 2, 2))
        assert_array_almost_equal(xty, [[1, 2], [13, 26], [41, 82]])
        assert_array_almost_equal(
            xtx,
            [[[1, 2], [2, 4]], [[13, 26], [26, 52]], [[41, 82], [82, 164]]],
        )

    def test_block_to_est_1d(self):
        x = np.arange(6).reshape((6, 1))
        y = np.arange(6).reshape((6, 1))
        for s in [[0, 3, 6], [0, 2, 4, 6], [0, 1, 5, 6]]:
            xty, xtx = jk.LstsqJackknifeFast.block_values(x, y, s)
            est = jk.LstsqJackknifeFast.block_values_to_est(xty, xtx)
            assert_array_equal(est.shape, (1, 1))
            assert_array_almost_equal(est, [[1.0]])

    def test_block_to_est_2d(self):
        x = np.vstack([np.arange(6), [1, 7, 6, 5, 2, 10]]).T
        y = np.atleast_2d(np.sum(x, axis=1)).T
        for s in [[0, 3, 6], [0, 2, 4, 6], [0, 1, 5, 6]]:
            xty, xtx = jk.LstsqJackknifeFast.block_values(x, y, s)
            est = jk.LstsqJackknifeFast.block_values_to_est(xty, xtx)
            assert_array_equal(est.shape, (1, 2))
            assert_array_almost_equal(est, [[1, 1]])

        with self.assertRaises(ValueError):
            jk.LstsqJackknifeFast.block_values_to_est(xty[0:2], xtx)
        with self.assertRaises(ValueError):
            jk.LstsqJackknifeFast.block_values_to_est(xty, xtx[:, :, 0:1])
        with self.assertRaises(ValueError):
            jk.LstsqJackknifeFast.block_values_to_est(xty, xtx[:, :, 0])

    def test_block_to_delete_1d(self):
        x = np.arange(6).reshape((6, 1))
        y = np.arange(6).reshape((6, 1))
        for s in [[0, 3, 6], [0, 2, 4, 6], [0, 1, 5, 6]]:
            xty, xtx = jk.LstsqJackknifeFast.block_values(x, y, s)
            delete = jk.LstsqJackknifeFast.block_values_to_delete_values(xty, xtx)
            assert_array_equal(delete.shape, (len(s) - 1, 1))
            assert_array_almost_equal(delete, np.ones_like(delete))

    def test_block_to_delete_2d(self):
        x = np.vstack([np.arange(6), [1, 7, 6, 5, 2, 10]]).T
        y = np.atleast_2d(np.sum(x, axis=1)).T
        for s in [[0, 3, 6], [0, 2, 4, 6], [0, 1, 5, 6]]:
            xty, xtx = jk.LstsqJackknifeFast.block_values(x, y, s)
            delete = jk.LstsqJackknifeFast.block_values_to_delete_values(xty, xtx)
            assert_array_equal(delete.shape, (len(s) - 1, 2))
            assert_array_almost_equal(delete, np.ones_like(delete))

    def test_eq_slow(self):
        x = np.atleast_2d(np.random.normal(size=(100, 2)))
        y = np.atleast_2d(np.random.normal(size=(100, 1)))
        for n_blocks in range(2, 49):
            b1 = jk.LstsqJackknifeFast(x, y, n_blocks=n_blocks).est
            b2 = jk.LstsqJackknifeSlow(x, y, n_blocks=n_blocks).est
            assert_array_almost_equal(b1, b2)

    def test_bad_data(self):
        x = np.arange(6).reshape((1, 6))
        with self.assertRaises(ValueError):
            jk.LstsqJackknifeFast(x, x, n_blocks=3)
        with self.assertRaises(ValueError):
            jk.LstsqJackknifeFast(x.T, x.T, n_blocks=8)
        with self.assertRaises(ValueError):
            jk.LstsqJackknifeFast(x.T, x.T, separators=list(range(10)))


class TestRatioJackknife(unittest.TestCase):
    def test_1d(self):
        numer_delete_values = np.matrix(np.arange(1, 11)).T
        denom_delete_values = -np.matrix(np.arange(1, 11)).T
        denom_delete_values[9, 0] += 1
        est = np.matrix(-1)
        jknife = jk.RatioJackknife(est, numer_delete_values, denom_delete_values)
        assert_array_almost_equal(jknife.est, est)
        assert_array_almost_equal(jknife.pseudovalues[0:9, :], -1)
        assert_array_almost_equal(jknife.pseudovalues[9, :], 0)
        assert_array_almost_equal(jknife.jknife_est[0, 0], -0.9)
        assert_array_almost_equal(jknife.jknife_se, 0.1)
        assert_array_almost_equal(jknife.jknife_var, 0.01)
        assert_array_almost_equal(jknife.jknife_cov, 0.01)

    def test_divide_by_zero_1d(self):
        est = np.ones((1, 1))
        numer_delete_vals = np.ones((10, 1))
        denom_delete_vals = np.ones((10, 1))
        denom_delete_vals[9, 0] = 0
        with self.assertRaises(FloatingPointError):
            jk.RatioJackknife(est, numer_delete_vals, denom_delete_vals)

    def test_2d(self):
        numer_delete_values = np.matrix(np.vstack((np.arange(1, 11), 2 * np.arange(1, 11)))).T
        x = -np.arange(1, 11)
        x[9] += 1
        denom_delete_values = np.vstack((x, 4 * x)).T
        est = np.matrix((-1, -0.5))
        jknife = jk.RatioJackknife(est, numer_delete_values, denom_delete_values)
        assert_array_almost_equal(jknife.est, est)
        self.assertEqual(jknife.est.shape, (1, 2))
        assert_array_almost_equal(jknife.pseudovalues[0:9, 0], -1)
        assert_array_almost_equal(jknife.pseudovalues[9, 0], 0)
        assert_array_almost_equal(jknife.pseudovalues[0:9, 1], -0.5)
        assert_array_almost_equal(jknife.pseudovalues[9, 1], 0)
        assert_array_almost_equal(jknife.jknife_est, [[-0.9, -0.45]])
        assert_array_almost_equal(jknife.jknife_se, [[0.1, 0.05]])
        assert_array_almost_equal(jknife.jknife_var, [[0.01, 0.0025]])
        assert_array_almost_equal(jknife.jknife_cov, np.matrix(((0.01, 0.005), (0.005, 0.0025))))

    def test_divide_by_zero_2d(self):
        est = np.ones((1, 2))
        numer_delete_vals = np.ones((10, 2))
        denom_delete_vals = np.ones((10, 2))
        denom_delete_vals[9, 0] = 0
        with self.assertRaises(FloatingPointError):
            jk.RatioJackknife(est, numer_delete_vals, denom_delete_vals)
