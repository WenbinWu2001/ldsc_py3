from pathlib import Path
import sys
import unittest
import warnings

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc.chromosome_inference import STANDARD_CHROMOSOMES, chrom_sort_key, normalize_chromosome


class ChromosomeInferenceTest(unittest.TestCase):
    def test_standard_chromosomes_define_package_wide_order(self):
        self.assertEqual(STANDARD_CHROMOSOMES[:4], ("1", "2", "3", "4"))
        self.assertEqual(STANDARD_CHROMOSOMES[-4:], ("X", "Y", "MT", "M"))

    def test_normalize_chromosome_warns_once_per_context(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            self.assertEqual(normalize_chromosome("chr1.0", context="shared-context"), "1")
            self.assertEqual(normalize_chromosome("chr1.0", context="shared-context"), "1")

        self.assertEqual(len(caught), 1)
        self.assertIn("canonical chromosome '1'", str(caught[0].message))
        self.assertIn("input chromosome 'chr1.0'", str(caught[0].message))

    def test_normalize_chromosome_warns_on_sex_chromosome_inference(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            self.assertEqual(normalize_chromosome("24", context="sex-context"), "Y")

        self.assertEqual(len(caught), 1)
        self.assertIn("sex chromosome", str(caught[0].message).lower())
        self.assertIn("'24'", str(caught[0].message))
        self.assertIn("'Y'", str(caught[0].message))

    def test_chrom_sort_key_uses_standard_package_order(self):
        values = ["Y", "2", "MT", "1", "M", "X"]
        self.assertEqual(sorted(values, key=chrom_sort_key), ["1", "2", "X", "Y", "MT", "M"])
