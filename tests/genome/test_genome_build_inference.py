import importlib.util
import logging
import unittest

import pandas as pd


from ldsc.errors import LDSCInputError, LDSCUsageError


def _load_module(testcase: unittest.TestCase):
    spec = importlib.util.find_spec("ldsc._genome.build_inference")
    testcase.assertIsNotNone(spec, "ldsc._genome.build_inference module is missing")
    return __import__("ldsc._genome.build_inference", fromlist=["dummy"])


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
        raw = pd.DataFrame(
            {
                "CHR": ["1"] * len(reference) + ["chrUn", "1"],
                "POS": reference["hg19_POS"].tolist() + [12345, "bad"],
                "SNP": [f"rs{i}" for i in range(len(reference))] + ["bad_chr", "bad_pos"],
            }
        )
        logger = logging.getLogger("LDSC.test.genome_build_inference.drop")

        with self.assertLogs(logger, level="WARNING") as caught:
            normalized, result = module.resolve_chr_pos_table(
                raw,
                context="unit-hg19-1based",
                reference_table=reference,
                logger=logger,
            )

        self.assertEqual(result.genome_build, "hg19")
        self.assertEqual(result.coordinate_basis, "1-based")
        self.assertEqual(result.inspected_snp_count, len(reference))
        self.assertEqual(normalized["POS"].tolist()[:3], [1000, 1010, 1020])
        self.assertEqual(result.coordinate_report.n_dropped, 2)
        self.assertEqual(normalized["SNP"].tolist()[-1], f"rs{len(reference) - 1}")
        self.assertIn("Dropped 2 SNPs with invalid or missing CHR/POS", "\n".join(caught.output))

    def test_resolve_chr_pos_table_strict_policy_rejects_invalid_coordinates(self):
        module = _load_module(self)
        reference = _build_reference_table()
        raw = pd.DataFrame({"CHR": ["1", "chrUn"], "POS": [1000, 1010]})

        with self.assertRaisesRegex(LDSCInputError, "Unsupported chromosome label"):
            module.resolve_chr_pos_table(
                raw,
                context="unit-strict",
                reference_table=reference,
                coordinate_policy="raise",
            )

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

    def test_below_threshold_error_reports_all_hypothesis_support(self):
        module = _load_module(self)
        reference = _build_reference_table()
        half = len(reference) // 2
        raw = pd.DataFrame(
            {
                "CHR": ["1"] * len(reference),
                "POS": reference["hg19_POS"].tolist()[:half] + reference["hg38_POS"].tolist()[half:],
            }
        )

        with self.assertRaises(LDSCInputError) as ctx:
            module.resolve_chr_pos_table(raw, context="unit-mixed", reference_table=reference)

        message = str(ctx.exception)
        self.assertIn(
            "Hypothesis support: hg19/1-based=125/250 (50.0%), hg19/0-based=0/250 (0.0%), "
            "hg38/1-based=125/250 (50.0%), hg38/0-based=0/250 (0.0%)",
            message,
        )
        self.assertIn("below the 99% confidence threshold", message)

    def test_resolve_chr_pos_table_reports_auto_inference_cause_fix_and_reference(self):
        module = _load_module(self)
        reference = _build_reference_table(n_rows=150)
        raw = pd.DataFrame({"CHR": ["1"] * len(reference), "POS": reference["hg19_POS"]})

        with self.assertRaises(LDSCInputError) as ctx:
            module.resolve_chr_pos_table(
                raw,
                context="unit-insufficient",
                reference_table=reference,
            )

        message = str(ctx.exception)
        self.assertIn("Could not infer genome build for unit-insufficient", message)
        self.assertIn("Most likely", message)
        self.assertIn("Pass --genome-build hg19 or --genome-build hg38", message)
        self.assertIn(
            "https://wenbinwu2001.github.io/ldsc_py3/user-guide/troubleshooting.html#common-genome-build-could-not-be-inferred",
            message,
        )


class TestResolveGenomeBuild(unittest.TestCase):
    def _frame(self, build, n=500):
        """Return a synthetic CHR/POS frame of packaged reference SNPs in ``build`` coordinates."""
        sample = _load_module(self).load_packaged_reference_table().head(n)
        return pd.DataFrame({"CHR": sample["CHR"], "POS": sample[f"{build}_POS"]})

    def _make_tiny_frame(self, n=5):
        """Frame too small to reach MIN_INSPECTED_REFERENCE_SNPS."""
        return pd.DataFrame({"CHR": ["1"] * n, "POS": list(range(n))})

    def test_hint_without_sample(self):
        # rsID modes carry no build; coordinate modes take the hint, with hg37 an alias of hg19.
        module = _load_module(self)
        for hint, mode, expected in [
            ("auto", "rsid", None), ("hg38", "rsid", None), (None, "rsid", None),
            ("hg38", "chr_pos", "hg38"), ("hg19", "chr_pos", "hg19"), ("hg37", "chr_pos", "hg19"),
        ]:
            with self.subTest(hint=hint, mode=mode):
                self.assertEqual(module.resolve_genome_build(hint, mode, None, context="test"), expected)

    def test_auto_sample_infers_its_build(self):
        module = _load_module(self)
        for build in ("hg38", "hg19"):
            with self.subTest(build=build):
                self.assertEqual(module.resolve_genome_build("auto", "chr_pos", self._frame(build), context="test"), build)

    def test_auto_none_sample_frame_raises(self):
        module = _load_module(self)
        with self.assertRaises(LDSCUsageError) as ctx:
            module.resolve_genome_build("auto", "chr_pos", None, context="test")
        self.assertIn("--genome-build", str(ctx.exception))

    def test_auto_insufficient_overlap_raises(self):
        module = _load_module(self)
        frame = self._make_tiny_frame()
        with self.assertRaises(LDSCInputError) as ctx:
            module.resolve_genome_build("auto", "chr_pos", frame, context="test")
        self.assertIn("--genome-build", str(ctx.exception))
        self.assertIn("Most likely", str(ctx.exception))
        self.assertIn("https://wenbinwu2001.github.io/ldsc_py3/user-guide/troubleshooting.html#common-genome-build-could-not-be-inferred", str(ctx.exception))

    def test_auto_logs_inference_summary(self):
        module = _load_module(self)
        frame = self._frame("hg38")
        with self.assertLogs(level="INFO") as log:
            module.resolve_genome_build(
                "auto",
                "chr_pos",
                frame,
                context="test",
                logger=logging.getLogger("test"),
            )
        self.assertTrue(any("Inferred genome build" in m for m in log.output))
