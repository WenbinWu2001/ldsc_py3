from pathlib import Path
import importlib.util
import logging
import sys
import unittest

import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def _load_module(testcase: unittest.TestCase):
    spec = importlib.util.find_spec("ldsc.genome_build_inference")
    testcase.assertIsNotNone(spec, "ldsc.genome_build_inference module is missing")
    return __import__("ldsc.genome_build_inference", fromlist=["dummy"])


def _build_reference_table(n_rows: int = 250) -> pd.DataFrame:
    hg19_positions = [1000 + (idx * 10) for idx in range(n_rows)]
    hg38_positions = [5000 + (idx * 10) for idx in range(n_rows)]
    return pd.DataFrame(
        {
            "CHR": ["1"] * n_rows,
            "hg19_POS": hg19_positions,
            "hg38_POS": hg38_positions,
        }
    )


class GenomeBuildInferenceTest(unittest.TestCase):
    def test_resolve_chr_pos_table_infers_hg19_one_based(self):
        module = _load_module(self)
        reference = _build_reference_table()
        raw = pd.DataFrame({"CHR": ["1"] * len(reference), "POS": reference["hg19_POS"]})

        normalized, result = module.resolve_chr_pos_table(
            raw,
            context="unit-hg19-1based",
            reference_table=reference,
        )

        self.assertEqual(result.genome_build, "hg19")
        self.assertEqual(result.coordinate_basis, "1-based")
        self.assertEqual(result.inspected_snp_count, len(reference))
        self.assertEqual(normalized["POS"].tolist()[:3], [1000, 1010, 1020])

    def test_resolve_chr_pos_table_infers_hg38_zero_based_and_shifts_positions(self):
        module = _load_module(self)
        reference = _build_reference_table()
        raw = pd.DataFrame({"CHR": ["1"] * len(reference), "POS": (reference["hg38_POS"] - 1).tolist()})

        logger = logging.getLogger("LDSC.test.genome_build_inference")
        with self.assertLogs(logger, level="WARNING") as captured:
            normalized, result = module.resolve_chr_pos_table(
                raw,
                context="unit-hg38-0based",
                reference_table=reference,
                logger=logger,
            )

        self.assertEqual(result.genome_build, "hg38")
        self.assertEqual(result.coordinate_basis, "0-based")
        self.assertEqual(normalized["POS"].tolist()[:3], [5000, 5010, 5020])
        self.assertIn("hg38", result.summary_message)
        self.assertIn("0-based", result.summary_message)
        self.assertTrue(any("converted positions to canonical 1-based coordinates" in message for message in captured.output))

    def test_resolve_chr_pos_table_rejects_mixed_evidence(self):
        module = _load_module(self)
        reference = _build_reference_table()
        half = len(reference) // 2
        raw = pd.DataFrame(
            {
                "CHR": ["1"] * len(reference),
                "POS": reference["hg19_POS"].tolist()[:half] + reference["hg38_POS"].tolist()[half:],
            }
        )

        with self.assertRaisesRegex(ValueError, "could not be inferred confidently"):
            module.resolve_chr_pos_table(
                raw,
                context="unit-mixed",
                reference_table=reference,
            )

    def test_resolve_chr_pos_table_requires_minimum_evidence(self):
        module = _load_module(self)
        reference = _build_reference_table(n_rows=150)
        raw = pd.DataFrame({"CHR": ["1"] * len(reference), "POS": reference["hg19_POS"]})

        with self.assertRaisesRegex(ValueError, "Insufficient evidence"):
            module.resolve_chr_pos_table(
                raw,
                context="unit-insufficient",
                reference_table=reference,
            )
