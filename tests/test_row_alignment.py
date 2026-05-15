from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc._row_alignment import assert_same_snp_rows


class RowAlignmentTest(unittest.TestCase):
    def make_rows(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "CHR": ["1", "1"],
                "SNP": ["rs1", "rs2"],
                "POS": [10, 20],
                "A1": ["A", "C"],
                "A2": ["C", "T"],
                "CM": [0.1, np.nan],
                "MAF": [0.2, np.nan],
            }
        )

    def make_rows_without_alleles(self) -> pd.DataFrame:
        return self.make_rows().drop(columns=["A1", "A2"])

    def test_exact_match_passes(self):
        assert_same_snp_rows(self.make_rows(), self.make_rows().copy(), context="test rows")

    def test_reordered_rows_fail(self):
        right = self.make_rows().iloc[[1, 0]].reset_index(drop=True)

        with self.assertRaisesRegex(ValueError, "test rows.*row 0"):
            assert_same_snp_rows(self.make_rows(), right, context="test rows")

    def test_default_chr_pos_identity_column_mismatch_fails(self):
        for column, value in (("CHR", "2"), ("POS", 99)):
            right = self.make_rows()
            right.loc[0, column] = value
            with self.subTest(column=column):
                with self.assertRaisesRegex(ValueError, "test rows.*identity_key"):
                    assert_same_snp_rows(self.make_rows(), right, context="test rows")

    def test_default_chr_pos_identity_ignores_snp_label_mismatch(self):
        right = self.make_rows()
        right.loc[0, "SNP"] = "rs9"

        assert_same_snp_rows(self.make_rows(), right, context="test rows")

    def test_explicit_rsid_identity_rejects_snp_label_mismatch(self):
        right = self.make_rows()
        right.loc[0, "SNP"] = "rs9"

        with self.assertRaisesRegex(ValueError, "test rows.*SNP"):
            assert_same_snp_rows(self.make_rows(), right, context="test rows", snp_identifier="rsid")

    def test_chr_pos_mode_ignores_snp_label_mismatch(self):
        right = self.make_rows()
        right["SNP"] = ["ld_rs1", "ld_rs2"]

        assert_same_snp_rows(self.make_rows(), right, context="test rows", snp_identifier="chr_pos")

    def test_chr_pos_mode_ignores_allele_mismatch(self):
        right = self.make_rows()
        right.loc[0, "A1"] = "A"
        right.loc[0, "A2"] = "G"

        assert_same_snp_rows(self.make_rows(), right, context="test rows", snp_identifier="chr_pos")

    def test_chr_pos_mode_still_rejects_coordinate_mismatch(self):
        right = self.make_rows()
        right["SNP"] = ["ld_rs1", "ld_rs2"]
        right.loc[0, "POS"] = 99

        with self.assertRaisesRegex(ValueError, "test rows.*identity_key.*1:10.*1:99"):
            assert_same_snp_rows(self.make_rows(), right, context="test rows", snp_identifier="chr_pos")

    def test_default_allele_aware_mode_rejects_missing_allele_columns(self):
        with self.assertRaisesRegex(ValueError, "test rows.*missing required.*A1.*A2"):
            assert_same_snp_rows(self.make_rows_without_alleles(), self.make_rows(), context="test rows")

    def test_default_allele_aware_mode_rejects_allele_mismatch(self):
        right = self.make_rows()
        right.loc[0, "A1"] = "A"
        right.loc[0, "A2"] = "G"

        with self.assertRaisesRegex(ValueError, "test rows.*identity_key.*1:10:A:C.*1:10:A:G"):
            assert_same_snp_rows(self.make_rows(), right, context="test rows")

    def test_float_metadata_within_float16_epsilon_passes(self):
        right = self.make_rows()
        delta = np.finfo(np.float16).eps / 4
        right.loc[0, "CM"] += delta
        right.loc[0, "MAF"] += delta

        assert_same_snp_rows(self.make_rows(), right, context="test rows")

    def test_float_metadata_larger_than_float16_epsilon_fails(self):
        right = self.make_rows()
        right.loc[0, "CM"] += np.finfo(np.float16).eps * 4

        with self.assertRaisesRegex(ValueError, "test rows.*CM"):
            assert_same_snp_rows(self.make_rows(), right, context="test rows")

    def test_nan_numeric_metadata_matches_nan(self):
        left = self.make_rows()
        right = self.make_rows()
        left.loc[1, "CM"] = np.nan
        right.loc[1, "CM"] = np.nan

        assert_same_snp_rows(left, right, context="test rows")

    def test_missing_numeric_metadata_on_one_side_fails(self):
        left = self.make_rows()
        right = self.make_rows()
        left.loc[0, "CM"] = np.nan

        with self.assertRaisesRegex(ValueError, "test rows.*CM"):
            assert_same_snp_rows(left, right, context="test rows")


if __name__ == "__main__":
    unittest.main()
