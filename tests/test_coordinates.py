from __future__ import annotations

import logging
from pathlib import Path
import sys
import unittest

import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc._coordinates import (
    build_chr_pos_key_frame,
    coordinate_missing_mask,
    normalize_chr_pos_frame,
)


class CoordinateHelperTest(unittest.TestCase):
    def test_coordinate_missing_mask_recognizes_shared_na_tokens(self):
        values = pd.Series([pd.NA, None, "", " ", "NA", "nan", "None", ".", "null", "1"])

        mask = coordinate_missing_mask(values)

        self.assertEqual(mask.tolist(), [True, True, True, True, True, True, True, True, True, False])

    def test_normalize_chr_pos_frame_drops_missing_and_reports_examples(self):
        frame = pd.DataFrame(
            {
                "SNP": ["rs1", "missing_chr", "missing_pos", "bad_chr", "bad_pos"],
                "CHR": ["chr1", ".", "2", "chrUn", "3"],
                "POS": ["100", "200", "NA", "300", "abc"],
            }
        )
        logger = logging.getLogger("tests.coordinates")

        with self.assertLogs(logger, level="WARNING") as caught:
            normalized, report = normalize_chr_pos_frame(
                frame,
                context="coordinate helper test",
                logger=logger,
            )

        self.assertEqual(normalized["SNP"].tolist(), ["rs1"])
        self.assertEqual(normalized["CHR"].tolist(), ["1"])
        self.assertEqual(normalized["POS"].tolist(), [100])
        self.assertEqual(report.n_dropped, 4)
        self.assertEqual(report.n_missing_chr, 1)
        self.assertEqual(report.n_missing_pos, 1)
        self.assertEqual(report.n_invalid_chr, 1)
        self.assertEqual(report.n_invalid_pos, 1)
        self.assertEqual(report.examples[0]["SNP"], "missing_chr")
        self.assertIn("Dropped 4 SNPs with invalid or missing CHR/POS", "\n".join(caught.output))

    def test_normalize_chr_pos_frame_strict_policy_errors_on_invalid_or_missing_coordinates(self):
        missing = pd.DataFrame({"CHR": ["1"], "POS": [pd.NA]})
        with self.assertRaisesRegex(ValueError, "requires valid CHR/POS coordinates"):
            normalize_chr_pos_frame(missing, context="missing POS", coordinate_policy="raise")

        bad_pos = pd.DataFrame({"CHR": ["1"], "POS": ["abc"]})
        with self.assertRaisesRegex(ValueError, "POS values in bad POS must be numeric"):
            normalize_chr_pos_frame(bad_pos, context="bad POS", coordinate_policy="raise")

        zero_pos = pd.DataFrame({"CHR": ["1"], "POS": [0]})
        with self.assertRaisesRegex(ValueError, "positive"):
            normalize_chr_pos_frame(zero_pos, context="zero POS", coordinate_policy="raise")

        bad_chr = pd.DataFrame({"CHR": ["foo"], "POS": [10]})
        with self.assertRaisesRegex(ValueError, "Unsupported chromosome label"):
            normalize_chr_pos_frame(bad_chr, context="bad CHR", coordinate_policy="raise")

    def test_build_chr_pos_key_frame_adds_canonical_key_after_drop(self):
        frame = pd.DataFrame({"CHR": ["chr1", pd.NA, "2"], "POS": [10, 20, 30], "value": [1, 2, 3]})

        keyed, report = build_chr_pos_key_frame(frame, context="key test", drop_missing=True)

        self.assertEqual(keyed["_ldsc_chr_pos_key"].tolist(), ["1:10", "2:30"])
        self.assertEqual(keyed["value"].tolist(), [1, 3])
        self.assertEqual(report.n_dropped, 1)


if __name__ == "__main__":
    unittest.main()
