from pathlib import Path
import sys
import unittest

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc._kernel import formats as ps

try:
    import bitarray as ba
    from ldsc._kernel import ldscore as ld
except ModuleNotFoundError:
    ba = None
    ld = None


FIXTURES = Path(__file__).resolve().parent / "fixtures" / "legacy"
PLINK_FIXTURES = FIXTURES / "plink_test"


@unittest.skipIf(ld is None or ba is None, "bitarray dependency is not installed")
class LDScoreHelpersTest(unittest.TestCase):
    def test_get_block_lefts(self):
        cases = [
            (np.arange(1, 6), 5, np.zeros(5)),
            (np.arange(1, 6), 0, np.arange(0, 5)),
            ((1, 4, 6, 7, 7, 8), 2, (0, 1, 1, 2, 2, 2)),
        ]
        for coords, max_dist, expected in cases:
            assert_array_equal(ld.getBlockLefts(coords, max_dist), expected)

    def test_block_left_to_right(self):
        cases = [
            ((0, 0, 0, 0, 0), (5, 5, 5, 5, 5)),
            ((0, 1, 2, 3, 4, 5), (1, 2, 3, 4, 5, 6)),
            ((0, 0, 2, 2), (2, 2, 4, 4)),
        ]
        for block_left, expected in cases:
            assert_array_equal(ld.block_left_to_right(block_left), expected)


@unittest.skipIf(ld is None or ba is None, "bitarray dependency is not installed")
class PlinkBedFileTest(unittest.TestCase):
    def setUp(self):
        self.snp_count = 8
        self.sample_count = 5
        self.bim = ps.PlinkBIMFile(str(PLINK_FIXTURES / "plink.bim"))

    def test_bed(self):
        bed = ld.PlinkBEDFile(str(PLINK_FIXTURES / "plink.bed"), self.sample_count, self.bim)
        self.assertEqual(bed.m, 4)
        self.assertEqual(bed.n, self.sample_count)
        self.assertEqual(len(bed.geno), 64)
        expected = np.array([0.6, 0.6, 0.625, 0.625])
        assert_array_almost_equal(bed.freq, expected)

    def test_filter_snps(self):
        bed = ld.PlinkBEDFile(
            str(PLINK_FIXTURES / "plink.bed"),
            self.sample_count,
            self.bim,
            keep_snps=[1, 4],
        )
        self.assertEqual(bed.m, 1)
        self.assertEqual(bed.n, 5)
        self.assertEqual(bed.geno[0:10], ba.bitarray("0001011111"))

    def test_filter_indivs(self):
        bed = ld.PlinkBEDFile(
            str(PLINK_FIXTURES / "plink.bed"),
            self.sample_count,
            self.bim,
            keep_indivs=[0, 1],
        )
        self.assertEqual(bed.m, 2)
        self.assertEqual(bed.n, 2)
        self.assertEqual(bed.geno[0:4], ba.bitarray("0001"))
        self.assertEqual(bed.geno[8:12], ba.bitarray("0001"))

    def test_filter_indivs_and_snps(self):
        bed = ld.PlinkBEDFile(
            str(PLINK_FIXTURES / "plink.bed"),
            self.sample_count,
            self.bim,
            keep_snps=[1, 5],
            keep_indivs=[0, 1],
        )
        self.assertEqual(bed.m, 1)
        self.assertEqual(bed.n, 2)
        self.assertEqual(bed.geno[0:4], ba.bitarray("0001"))

    def test_bad_filename(self):
        with self.assertRaises(ValueError):
            ld.PlinkBEDFile(str(PLINK_FIXTURES / "plink.bim"), 9, self.bim)

    def test_next_snps_errors(self):
        bed = ld.PlinkBEDFile(str(PLINK_FIXTURES / "plink.bed"), self.sample_count, self.bim)
        with self.assertRaises(ValueError):
            bed.nextSNPs(0)
        with self.assertRaises(ValueError):
            bed.nextSNPs(5)

    def test_next_snps(self):
        for width in [1, 2, 3]:
            bed = ld.PlinkBEDFile(
                str(PLINK_FIXTURES / "plink.bed"), self.sample_count, self.bim
            )
            x = bed.nextSNPs(width)
            self.assertEqual(x.shape, (5, width))
            self.assertTrue(np.all(np.abs(np.mean(x, axis=0)) < 0.01))
            self.assertTrue(np.all(np.abs(np.std(x, axis=0) - 1) < 0.01))

    def test_next_snps_maf_ref(self):
        width = 4
        bed = ld.PlinkBEDFile(str(PLINK_FIXTURES / "plink.bed"), self.sample_count, self.bim)
        x = bed.nextSNPs(width)
        bed._currentSNP -= width
        y = bed.nextSNPs(width, minorRef=True)
        assert_array_equal(x, -y)
