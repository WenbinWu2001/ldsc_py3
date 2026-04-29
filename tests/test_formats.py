from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd
from numpy.testing import assert_array_equal
from pandas.testing import assert_frame_equal

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc._kernel import formats as ps


FIXTURES = Path(__file__).resolve().parent / "fixtures"
PARSE_FIXTURES = FIXTURES / "formats"
PLINK_FIXTURES = FIXTURES / "plink"


class ParseHelpersTest(unittest.TestCase):
    def test_series_eq(self):
        x = pd.Series([1, 2, 3])
        y = pd.Series([1, 2])
        z = pd.Series([1, 2, 4])
        self.assertTrue(ps.series_eq(x, x))
        self.assertFalse(ps.series_eq(x, y))
        self.assertFalse(ps.series_eq(x, z))

    def test_get_compression(self):
        self.assertEqual(ps.get_compression("gz"), "gzip")
        self.assertEqual(ps.get_compression("bz2"), "bz2")
        self.assertIsNone(ps.get_compression("asdf"))

    def test_read_cts(self):
        match_snps = pd.Series(["rs1", "rs2", "rs3"])
        path = PARSE_FIXTURES / "test.cts"
        assert_array_equal(ps.read_cts(str(path), match_snps), [1, 2, 3])
        with self.assertRaises(ValueError):
            ps.read_cts(str(path), match_snps.iloc[0:2])

    def test_read_sumstats(self):
        x = ps.sumstats(str(PARSE_FIXTURES / "test.sumstats"), dropna=True, alleles=True)
        self.assertEqual(len(x), 1)
        assert_array_equal(x.SNP, ["rs1"])
        with self.assertRaises(ValueError):
            ps.sumstats(str(PARSE_FIXTURES / "test.l2.ldscore.gz"))

    def test_frq_parser(self):
        x = ps.frq_parser(str(PARSE_FIXTURES / "test1.frq"), compression=None)
        assert_array_equal(x.columns, ["SNP", "FRQ"])
        assert_array_equal(x.SNP, [f"rs_{i}" for i in range(8)])
        assert_array_equal(x.FRQ, [0.01, 0.1, 0.7, 0.2, 0.2, 0.2, 0.99, 0.03])

        y = ps.frq_parser(str(PARSE_FIXTURES / "test2.frq.gz"), compression="gzip")
        assert_array_equal(y.columns, ["SNP", "FRQ"])
        assert_array_equal(y.SNP, [f"rs_{i}" for i in range(8)])
        assert_array_equal(y.FRQ, [0.01, 0.1, 0.3, 0.2, 0.2, 0.2, 0.01, 0.03])


class LDScoreParsingTest(unittest.TestCase):
    def test_ldscore(self):
        x = ps.ldscore(str(PARSE_FIXTURES / "test"))
        self.assertEqual(list(x["SNP"]), [f"rs{i}" for i in range(1, 23)])
        self.assertEqual(list(x["AL2"]), list(range(1, 23)))
        self.assertEqual(list(x["BL2"]), list(range(2, 46, 2)))

    def test_ldscore_loop(self):
        x = ps.ldscore(str(PARSE_FIXTURES / "test"), 2)
        self.assertEqual(list(x["SNP"]), ["rs1", "rs2"])
        self.assertEqual(list(x["AL2"]), [1, 2])
        self.assertEqual(list(x["BL2"]), [2, 4])

    def test_ldscore_fromlist(self):
        fh = str(PARSE_FIXTURES / "test")
        x = ps.ldscore_fromlist([fh, fh])
        self.assertEqual(x.shape, (22, 5))
        y = ps.ldscore(fh)
        assert_array_equal(x.iloc[:, 0:3].to_numpy(), y.to_numpy())
        assert_array_equal(x.iloc[:, [0, 3, 4]].to_numpy(), y.to_numpy())
        with self.assertRaises(ValueError):
            ps.ldscore_fromlist([fh, str(PARSE_FIXTURES / "test2")])


class MParsingTest(unittest.TestCase):
    def test_bad_M(self):
        with self.assertRaises(ValueError):
            ps.M(str(PARSE_FIXTURES / "test_bad"))

    def test_M(self):
        x = ps.M(str(PARSE_FIXTURES / "test"))
        self.assertEqual(x.shape, (1, 3))
        assert_array_equal(x, [[1000, 2000, 3000]])

    def test_M_loop(self):
        x = ps.M(str(PARSE_FIXTURES / "test"), 2)
        self.assertEqual(x.shape, (1, 2))
        assert_array_equal(x, [[3, 6]])

    def test_M_fromlist(self):
        fh = str(PARSE_FIXTURES / "test")
        x = ps.M_fromlist([fh, fh])
        self.assertEqual(x.shape, (1, 6))
        assert_array_equal(x, np.hstack((ps.M(fh), ps.M(fh))))


class PlinkFilesTest(unittest.TestCase):
    def test_fam(self):
        fam = ps.PlinkFAMFile(str(PLINK_FIXTURES / "plink.fam"))
        self.assertEqual(fam.n, 5)
        correct = np.array(["per0", "per1", "per2", "per3", "per4"])
        assert_array_equal(fam.IDList.values.reshape((5,)), correct)

    def test_fam_bad_filename(self):
        with self.assertRaises(ValueError):
            ps.PlinkFAMFile(str(PLINK_FIXTURES / "plink.bim"))

    def test_bim(self):
        bim = ps.PlinkBIMFile(str(PLINK_FIXTURES / "plink.bim"))
        self.assertEqual(bim.n, 8)
        correct = np.array([f"rs_{i}" for i in range(8)])
        assert_array_equal(bim.IDList.values.reshape(8), correct)

    def test_bim_bad_filename(self):
        with self.assertRaises(ValueError):
            ps.PlinkBIMFile(str(PLINK_FIXTURES / "plink.fam"))
