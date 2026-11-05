"""Genetic correlation over trait pairs: pair order, relabeling, anchors, and failed or unestimable pairs."""

from dataclasses import replace
import unittest
from unittest import mock

import numpy as np

from ldsc.config import GlobalConfig, RegressionConfig
from ldsc.errors import LDSCInputError
from ldsc._regression import fitting, rg
from ldsc._regression.results import RgResultFamily
from ldsc._regression.runner import RegressionRunner
from tests.regression_fixtures import RegressionFixtures, stub_h2_kernel


class RgTest(RegressionFixtures, unittest.TestCase):
    def test_estimate_rg_pairs_lists_every_trait_without_alleles_before_fitting(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        traits = [self.allele_free_table("A"), self.make_sumstats_table(), self.allele_free_table("C")]

        with mock.patch.object(runner, "build_dataset", side_effect=AssertionError("fitting started")):
            with self.assertRaisesRegex(LDSCInputError, r"A \('A.parquet'\), C \('C.parquet'\)(.|\n)*munge-sumstats"):
                runner.estimate_rg_pairs(traits, self.make_ldscore_result())

    def test_estimate_rg_pairs_returns_explicit_result_family_in_input_order(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        traits = [
            replace(self.make_sumstats_table(), trait_name="A", source_path="A.sumstats.gz"),
            replace(self.make_sumstats_table(), trait_name="B", source_path="B.sumstats.gz"),
            replace(self.make_sumstats_table(), trait_name="C", source_path="C.sumstats.gz"),
        ]

        with stub_h2_kernel(self.make_rg_kernel_result().hsq1,
        ), mock.patch.object(
            fitting.reg,
            "RG",
            side_effect=[
                self.make_rg_kernel_result(rg=0.30, p=0.03),
                self.make_rg_kernel_result(rg=0.10, p=0.01),
                self.make_rg_kernel_result(rg=0.20, p=0.02),
            ],
        ):
            result = runner.estimate_rg_pairs(traits, self.make_ldscore_result())

        self.assertIsInstance(result, RgResultFamily)
        self.assertEqual(result.rg[["trait_1", "trait_2"]].values.tolist(), [["A", "B"], ["A", "C"], ["B", "C"]])
        self.assertEqual(result.rg_full["pair_kind"].tolist(), ["all_pairs", "all_pairs", "all_pairs"])
        self.assertEqual(result.rg_full["status"].tolist(), ["ok", "ok", "ok"])
        self.assertEqual(result.rg["n_snps_used"].tolist(), [3, 3, 3])
        self.assertEqual(result.h2_per_trait["trait_name"].tolist(), ["A", "B", "C"])
        self.assertEqual(len(result.per_pair_metadata), 3)
        self.assertEqual(result.per_pair_metadata[0]["effective_snp_identifier"], "rsid")
        self.assertFalse(result.per_pair_metadata[0]["identity_downgrade_applied"])

    def test_estimate_rg_pairs_relabels_shared_trait_names_like_the_cli(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        traits = [
            replace(self.make_sumstats_table(), trait_name="SCZ", source_path="pgc/v3/SCZ.parquet"),
            replace(self.make_sumstats_table(), trait_name="SCZ", source_path="ukb/v3/SCZ.parquet"),
        ]

        with stub_h2_kernel(self.make_rg_kernel_result().hsq1), mock.patch.object(
            fitting.reg, "RG", return_value=self.make_rg_kernel_result(rg=0.5, p=0.01)
        ):
            result = runner.estimate_rg_pairs(traits, self.make_ldscore_result())

        # Documented example in docs/user-guide/workflows/genetic-correlation.md.
        self.assertEqual(result.h2_per_trait["trait_name"].tolist(), ["SCZ@pgc/v3", "SCZ@ukb/v3"])
        self.assertEqual(result.rg[["trait_1", "trait_2"]].values.tolist(), [["SCZ@pgc/v3", "SCZ@ukb/v3"]])

    def test_estimate_rg_pairs_rejects_shared_trait_names_without_source_paths(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        traits = [replace(self.make_sumstats_table(), trait_name="SCZ", source_path=None) for _ in range(2)]

        with self.assertRaisesRegex(LDSCInputError, r"'SCZ'(.|\n)*distinct trait_name"):
            runner.estimate_rg_pairs(traits, self.make_ldscore_result())

    def test_estimate_rg_pairs_anchor_mode_uses_anchor_against_rest_in_input_order(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        traits = [
            replace(self.make_sumstats_table(), trait_name="A", source_path="A.sumstats.gz"),
            replace(self.make_sumstats_table(), trait_name="B", source_path="B.sumstats.gz"),
            replace(self.make_sumstats_table(), trait_name="C", source_path="C.sumstats.gz"),
        ]

        with stub_h2_kernel(self.make_rg_kernel_result().hsq1), mock.patch.object(
            fitting.reg,
            "RG",
            return_value=self.make_rg_kernel_result(),
        ):
            result = runner.estimate_rg_pairs(traits, self.make_ldscore_result(), anchor_index=1)

        self.assertEqual(result.rg[["trait_1", "trait_2"]].values.tolist(), [["B", "A"], ["B", "C"]])
        self.assertEqual(result.rg_full["pair_kind"].tolist(), ["anchor", "anchor"])

    def test_estimate_rg_pairs_records_pair_exceptions_and_continues(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        traits = [
            replace(self.make_sumstats_table(), trait_name="A", source_path="A.sumstats.gz"),
            replace(self.make_sumstats_table(), trait_name="B", source_path="B.sumstats.gz"),
            replace(self.make_sumstats_table(), trait_name="C", source_path="C.sumstats.gz"),
        ]

        with stub_h2_kernel(self.make_rg_kernel_result().hsq1), mock.patch.object(
            fitting.reg,
            "RG",
            side_effect=[
                self.make_rg_kernel_result(rg=0.30, p=0.03),
                RuntimeError("pair exploded"),
                self.make_rg_kernel_result(rg=0.20, p=0.02),
            ],
        ), self.assertLogs("LDSC.regression_runner", level="WARNING") as logs:
            result = runner.estimate_rg_pairs(traits, self.make_ldscore_result())

        self.assertEqual(result.rg_full["status"].tolist(), ["ok", "failed", "ok"])
        self.assertEqual(result.rg_full.loc[1, "error"], "RuntimeError: pair exploded")
        self.assertTrue(np.isnan(result.rg.loc[1, "rg"]))
        self.assertIn("Failed", result.rg.loc[1, "note"])
        self.assertTrue(any("pair exploded" in record for record in logs.output))

    def test_estimate_rg_pairs_marks_out_of_bounds_h2_unestimable_but_preserves_numeric_out_of_range_rg(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        traits = [
            replace(self.make_sumstats_table(), trait_name="A", source_path="A.sumstats.gz"),
            replace(self.make_sumstats_table(), trait_name="B", source_path="B.sumstats.gz"),
            replace(self.make_sumstats_table(), trait_name="C", source_path="C.sumstats.gz"),
        ]

        with stub_h2_kernel(self.make_rg_kernel_result().hsq1), mock.patch.object(
            fitting.reg,
            "RG",
            side_effect=[
                self.make_rg_kernel_result(rg="NA", rg_se="NA", z="NA", p="NA", h2_2=-0.01),
                self.make_rg_kernel_result(rg=1.5, rg_se=0.1, z=15.0, p=1e-6),
                self.make_rg_kernel_result(rg=0.2, rg_se=0.1, z=2.0, p=0.04),
            ],
        ):
            result = runner.estimate_rg_pairs(traits, self.make_ldscore_result())

        self.assertEqual(result.rg_full.loc[0, "status"], "unestimable")
        self.assertTrue(np.isnan(result.rg.loc[0, "rg"]))
        self.assertIn("h2_2_obs (trait 'B') is -0.01", result.rg_full.loc[0, "error"])
        self.assertEqual(result.rg_full.loc[0, "h2_2_obs"], -0.01)
        self.assertEqual(result.rg.loc[0, "note"], rg.H2_OUT_OF_BOUNDS_RG_NOTE)
        self.assertEqual(result.rg_full.loc[1, "status"], "ok")
        self.assertEqual(result.rg.loc[1, "rg"], 1.5)
