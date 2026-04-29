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


class TestResolveGenomeBuild(unittest.TestCase):
    def _make_hg38_frame(self, n=500):
        """Return a synthetic CHR/POS frame in hg38 coordinates."""
        module = _load_module(self)
        ref = module.load_packaged_reference_table()
        sample = ref.head(n).copy()
        return pd.DataFrame({"CHR": sample["CHR"], "POS": sample["hg38_POS"]})

    def _make_hg19_frame(self, n=500):
        module = _load_module(self)
        ref = module.load_packaged_reference_table()
        sample = ref.head(n).copy()
        return pd.DataFrame({"CHR": sample["CHR"], "POS": sample["hg19_POS"]})

    def _make_tiny_frame(self, n=5):
        """Frame too small to reach MIN_INSPECTED_REFERENCE_SNPS."""
        return pd.DataFrame({"CHR": ["1"] * n, "POS": list(range(n))})

    def test_rsid_auto_hint_returns_none(self):
        module = _load_module(self)
        self.assertIsNone(module.resolve_genome_build("auto", "rsid", None, context="test"))

    def test_rsid_hg38_hint_returns_none(self):
        module = _load_module(self)
        self.assertIsNone(module.resolve_genome_build("hg38", "rsid", None, context="test"))

    def test_rsid_none_hint_returns_none(self):
        module = _load_module(self)
        self.assertIsNone(module.resolve_genome_build(None, "rsid", None, context="test"))

    def test_chr_pos_hg38_hint_returns_hg38(self):
        module = _load_module(self)
        self.assertEqual(module.resolve_genome_build("hg38", "chr_pos", None, context="test"), "hg38")

    def test_chr_pos_hg19_hint_returns_hg19(self):
        module = _load_module(self)
        self.assertEqual(module.resolve_genome_build("hg19", "chr_pos", None, context="test"), "hg19")

    def test_chr_pos_alias_hg37_returns_hg19(self):
        module = _load_module(self)
        self.assertEqual(module.resolve_genome_build("hg37", "chr_pos", None, context="test"), "hg19")

    def test_auto_hg38_sample_infers_hg38(self):
        module = _load_module(self)
        frame = self._make_hg38_frame()
        result = module.resolve_genome_build("auto", "chr_pos", frame, context="test")
        self.assertEqual(result, "hg38")

    def test_auto_hg19_sample_infers_hg19(self):
        module = _load_module(self)
        frame = self._make_hg19_frame()
        result = module.resolve_genome_build("auto", "chr_pos", frame, context="test")
        self.assertEqual(result, "hg19")

    def test_auto_none_sample_frame_raises(self):
        module = _load_module(self)
        with self.assertRaises(ValueError) as ctx:
            module.resolve_genome_build("auto", "chr_pos", None, context="test")
        self.assertIn("--genome-build", str(ctx.exception))

    def test_auto_insufficient_overlap_raises(self):
        module = _load_module(self)
        frame = self._make_tiny_frame()
        with self.assertRaises(ValueError) as ctx:
            module.resolve_genome_build("auto", "chr_pos", frame, context="test")
        self.assertIn("--genome-build", str(ctx.exception))

    def test_auto_logs_inference_summary(self):
        module = _load_module(self)
        frame = self._make_hg38_frame()
        with self.assertLogs(level="INFO") as log:
            module.resolve_genome_build(
                "auto",
                "chr_pos",
                frame,
                context="test",
                logger=logging.getLogger("test"),
            )
        self.assertTrue(any("Inferred genome build" in m for m in log.output))
