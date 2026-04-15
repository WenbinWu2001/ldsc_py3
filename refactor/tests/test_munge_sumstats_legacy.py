import unittest

import numpy as np
import pandas as pd
from numpy.testing import assert_allclose, assert_array_equal
from pandas.testing import assert_frame_equal, assert_series_equal

try:
    from ldsc._kernel import sumstats_munger as munge
except ModuleNotFoundError:
    munge = None


class Mock(object):
    def log(self, _message):
        return None


@unittest.skipIf(munge is None, "ld-score runtime dependencies are not installed")
class PToZTest(unittest.TestCase):
    def setUp(self):
        self.N = pd.Series([1, 2, 3])
        self.P = pd.Series([0.1, 0.1, 0.1])
        self.Z = pd.Series([1.644854, 1.644854, 1.644854])

    def test_p_to_z(self):
        assert_allclose(munge.p_to_z(self.P, self.N), self.Z, atol=1e-5)


@unittest.skipIf(munge is None, "ld-score runtime dependencies are not installed")
class CheckMedianTest(unittest.TestCase):
    def setUp(self):
        self.x = pd.Series([1, 2, 3])

    def test_good_median(self):
        message = munge.check_median(self.x, 2, 0, "TEST")
        self.assertEqual(message, "Median value of TEST was 2.0, which seems sensible.")

    def test_bad_median(self):
        with self.assertRaises(ValueError):
            munge.check_median(self.x, 0, 0.1, "TEST")


@unittest.skipIf(munge is None, "ld-score runtime dependencies are not installed")
class ProcessNTest(unittest.TestCase):
    def setUp(self):
        self.log = Mock()
        self.dat = pd.DataFrame(["rs1", "rs2", "rs3"], columns=["SNP"])
        self.dat_filtered = pd.DataFrame(["rs2", "rs3"], columns=["SNP"])
        self.dat_filtered["N"] = [1234, 1234.0]
        self.dat_filtered9999 = pd.DataFrame(["rs2", "rs3"], columns=["SNP"])
        self.dat_filtered9999["N"] = [9999, 9999.0]
        self.args = munge.parser.parse_args("")
        self.args.N = 9999.0
        self.args.N_cas = 9999.0
        self.args.N_con = 9999.0
        self.N_const = pd.Series([1234, 1234, 1234.0])
        self.N = pd.Series([1, 1234, 1234.0])

    def test_n_col(self):
        self.dat["N"] = self.N
        dat = munge.process_n(self.dat, self.args, self.log)
        assert_frame_equal(dat, self.dat_filtered)

    def test_nstudy(self):
        self.dat["NSTUDY"] = self.N
        dat = munge.process_n(self.dat, self.args, self.log)
        assert_frame_equal(dat, self.dat_filtered9999)

    def test_n_cas_con_col(self):
        self.dat["N_CAS"] = self.N
        self.dat["N_CON"] = [0.0, 0, 0]
        dat = munge.process_n(self.dat, self.args, self.log)
        assert_frame_equal(dat, self.dat_filtered)

    def test_n_flag(self):
        self.args.N = 1234.0
        self.args.N_cas = None
        self.args.N_con = None
        dat = munge.process_n(self.dat, self.args, self.log)
        assert_series_equal(dat.N, self.N_const, check_names=False)

    def test_n_cas_con_flag(self):
        self.args.N = None
        self.args.N_cas = 1000.0
        self.args.N_con = 234.0
        dat = munge.process_n(self.dat, self.args, self.log)
        assert_series_equal(dat.N, self.N_const, check_names=False)


@unittest.skipIf(munge is None, "ld-score runtime dependencies are not installed")
class MungeFiltersTest(unittest.TestCase):
    def setUp(self):
        self.log = Mock()
        self.args = munge.parser.parse_args("")

    def test_filter_pvals(self):
        pvals = pd.Series([0, 0.1, 1, 2])
        assert_series_equal(
            munge.filter_pvals(pvals, self.log, self.args),
            pd.Series([False, True, True, False]),
        )

    def test_single_info(self):
        dat = pd.Series([0.8, 1, 1])
        assert_series_equal(
            munge.filter_info(dat, self.log, self.args),
            pd.Series([False, True, True]),
        )

    def test_multiple_info(self):
        info1 = pd.Series([0.8, 1, 1])
        info2 = pd.Series([1.01, 0.5, 9])
        dat = pd.concat([info1, info2], axis=1).reset_index(drop=True)
        dat.columns = ["INFO", "INFO"]
        assert_series_equal(
            munge.filter_info(dat, self.log, self.args),
            pd.Series([True, False, True]),
        )

    def test_filter_frq(self):
        frq = pd.Series([-1, 0, 0.005, 0.4, 0.6, 0.999, 1, 2])
        assert_series_equal(
            munge.filter_frq(frq, self.log, self.args),
            pd.Series([False, False, False, True, True, False, False, False]),
        )

    def test_filter_alleles(self):
        alleles = pd.Series(
            ["AC", "AG", "CA", "CT", "GA", "GT", "TC", "TG", "DI", "AAT", "RA"]
        )
        expected = pd.Series([i < 8 for i in range(11)])
        assert_series_equal(munge.filter_alleles(alleles), expected)


@unittest.skipIf(munge is None, "ld-score runtime dependencies are not installed")
class AlleleMergeTest(unittest.TestCase):
    def setUp(self):
        self.log = Mock()
        self.dat = pd.DataFrame(
            np.transpose([["a", "b", "c"], ["A", "T", "C"], ["C", "G", "A"]]),
            columns=["SNP", "A1", "A2"],
        )
        self.alleles = pd.DataFrame(
            np.transpose([["a", "extra", "b", "c"], ["AG", "TC", "AC", "AC"]]),
            columns=["SNP", "MA"],
        )

    def test_merge(self):
        merged = munge.allele_merge(self.dat, self.alleles, self.log)
        expected = pd.DataFrame(
            np.transpose([["a", "extra", "b", "c"], ["a", "a", "T", "C"], ["a", "a", "G", "A"]]),
            columns=["SNP", "A1", "A2"],
        )
        expected.loc[[0, 1], ["A1", "A2"]] = float("nan")
        assert_frame_equal(merged, expected)


@unittest.skipIf(munge is None, "ld-score runtime dependencies are not installed")
class ParseDatTest(unittest.TestCase):
    def setUp(self):
        self.log = Mock()
        dat = pd.DataFrame()
        dat["SNP"] = [f"rs{i}" for i in range(10)]
        dat["A1"] = ["A" for _ in range(10)]
        dat["A2"] = ["G" for _ in range(10)]
        dat["INFO"] = np.ones(10)
        dat["FRQ"] = np.ones(10) / 2
        dat["P"] = np.ones(10)
        self.dat = dat
        self.dat_gen = [dat.loc[0:4, :], dat.loc[5:9, :].reset_index(drop=True)]
        self.convert_colname = {x: x for x in self.dat_gen[0].columns}
        self.args = munge.parser.parse_args("")

    def test_no_alleles(self):
        dat = self.dat.drop(["A1", "A2"], axis=1)
        dat_gen = [dat.loc[0:4, :], dat.loc[5:9, :].reset_index(drop=True)]
        self.args.no_alleles = True
        parsed = munge.parse_dat(dat_gen, self.convert_colname, None, self.log, self.args)
        assert_frame_equal(parsed, self.dat.drop(["INFO", "FRQ", "A1", "A2"], axis=1))

    def test_merge_alleles(self):
        self.args.merge_alleles = True
        merge_alleles = pd.DataFrame(
            {"SNP": [f"rs{i}" for i in range(3)], "MA": ["AG", "AG", "AG"]}
        )
        parsed = munge.parse_dat(
            self.dat_gen, self.convert_colname, merge_alleles, self.log, self.args
        )
        assert_frame_equal(parsed, self.dat.loc[0:2, ["SNP", "A1", "A2", "P"]])

    def test_standard(self):
        parsed = munge.parse_dat(self.dat_gen, self.convert_colname, None, self.log, self.args)
        assert_frame_equal(parsed, self.dat.drop(["INFO", "FRQ"], axis=1))

    def test_na(self):
        self.dat.loc[0, "SNP"] = float("nan")
        self.dat.loc[1, "A2"] = float("nan")
        self.dat_gen = [self.dat.loc[0:4, :], self.dat.loc[5:9, :].reset_index(drop=True)]
        parsed = munge.parse_dat(self.dat_gen, self.convert_colname, None, self.log, self.args)
        expected = self.dat.loc[2:, ["SNP", "A1", "A2", "P"]].reset_index(drop=True)
        assert_frame_equal(parsed, expected)


@unittest.skipIf(munge is None, "ld-score runtime dependencies are not installed")
class ParseFlagsAndHeadersTest(unittest.TestCase):
    def setUp(self):
        self.log = Mock()
        self.args = munge.parser.parse_args("")

    def test_clean_header(self):
        self.assertEqual(munge.clean_header("foo-bar.foo_BaR"), "FOO_BAR_FOO_BAR")

    def test_get_compression(self):
        _open_func, compression = munge.get_compression("foo.gz")
        self.assertEqual(compression, "gzip")
        _open_func, compression = munge.get_compression("foo.bz2")
        self.assertEqual(compression, "bz2")
        _open_func, compression = munge.get_compression("foo.bar")
        self.assertIsNone(compression)

    def test_parse_flag_cnames(self):
        self.args.nstudy = "nstudy1"
        self.args.snp = "snp1"
        self.args.N_col = "n.col1"
        self.args.N_cas_col = "n-cas.col1"
        self.args.N_con_col = "n-con.col1"
        self.args.a1 = "a11"
        self.args.a2 = "a21"
        self.args.p = "p1"
        self.args.frq = "frq1"
        self.args.info = "info1"
        self.args.info_list = "info111,info222"
        self.args.signed_sumstats = "beta1,0"
        cname_map, null_value = munge.parse_flag_cnames(self.log, self.args)
        self.assertEqual(null_value, 0)
        self.assertEqual(cname_map["NSTUDY1"], "NSTUDY")
        self.assertEqual(cname_map["SNP1"], "SNP")
        self.assertEqual(cname_map["N_COL1"], "N")
        self.assertEqual(cname_map["N_CAS_COL1"], "N_CAS")
        self.assertEqual(cname_map["N_CON_COL1"], "N_CON")
        self.assertEqual(cname_map["A11"], "A1")
        self.assertEqual(cname_map["A21"], "A2")
        self.assertEqual(cname_map["P1"], "P")
        self.assertEqual(cname_map["FRQ1"], "FRQ")
        self.assertEqual(cname_map["INFO1"], "INFO")
        self.assertEqual(cname_map["INFO111"], "INFO")
        self.assertEqual(cname_map["INFO222"], "INFO")
        self.assertEqual(cname_map["BETA1"], "SIGNED_SUMSTAT")

    def test_get_cname_map(self):
        flag = {"SNP": "SNP"}
        default = {"SNPID": "SNP", "PVAL": "P"}
        ignore = ["pval"]
        cname_map = munge.get_cname_map(flag, default, ignore)
        self.assertEqual(cname_map, {"SNP": "SNP"})

    def test_filter_signed_sumstats_null(self):
        series = pd.Series([0.0, 1.0, -1.0])
        retained = munge.filter_signed_sumstats(series, 0.0, self.log)
        assert_array_equal(retained.values, [False, True, True])
