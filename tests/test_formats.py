from pathlib import Path
import sys
import unittest

import numpy as np
from numpy.testing import assert_array_equal

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc._kernel import formats as ps
from ldsc.errors import LDSCInputError


FIXTURES = Path(__file__).resolve().parent / "fixtures"
PLINK_FIXTURES = FIXTURES / "plink"


class PlinkFilesTest(unittest.TestCase):
    def test_fam(self):
        fam = ps.PlinkFAMFile(str(PLINK_FIXTURES / "plink.fam"))
        self.assertEqual(fam.n, 5)
        correct = np.array(["per0", "per1", "per2", "per3", "per4"])
        assert_array_equal(fam.IDList.values.reshape((5,)), correct)

    def test_fam_bad_filename(self):
        with self.assertRaises(LDSCInputError):
            ps.PlinkFAMFile(str(PLINK_FIXTURES / "plink.bim"))

    def test_bim(self):
        bim = ps.PlinkBIMFile(str(PLINK_FIXTURES / "plink.bim"))
        self.assertEqual(bim.n, 8)
        correct = np.array([f"rs_{i}" for i in range(8)])
        assert_array_equal(bim.IDList.values.reshape(8), correct)

    def test_bim_bad_filename(self):
        with self.assertRaises(LDSCInputError):
            ps.PlinkBIMFile(str(PLINK_FIXTURES / "plink.fam"))
