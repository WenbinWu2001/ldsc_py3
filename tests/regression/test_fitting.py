"""Fit policy, h2 and rg fits, and the result rows and metadata built from fitted estimators."""

from dataclasses import replace
from types import SimpleNamespace
import json
import math
import unittest
from unittest import mock

import numpy as np
import pandas as pd
from scipy import stats

from ldsc.config import GlobalConfig, RegressionConfig
from ldsc.errors import LDSCInputError
from ldsc._annotation.overlap import OverlapContribution
from ldsc._ldscore.overlap_table import LDScoreOverlap
from ldsc._regression import fitting
from ldsc._regression.outputs import H2_REGRESSION_BIN_COLUMNS, PARTITIONED_H2_COLUMNS
from ldsc._regression.results import RegressionDataset
from ldsc._regression.runner import RegressionRunner
from tests.regression_fixtures import RegressionFixtures, fake_partitioned_hsq, fake_total_h2_hsq, stub_h2_kernel


class FittingTest(RegressionFixtures, unittest.TestCase):
    def test_single_annotation_h2_applies_no_default_chisq_cap(self):
        # Legacy estimate_h2 leaves single-annotation h2 uncapped (relies on the
        # two-step estimator). All three SNPs must reach the kernel even though
        # rs1 has chi^2=144, which the partitioned default cap would have dropped.
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        sumstats = self.make_sumstats_table_from_frame(
            pd.DataFrame(
                {
                    "SNP": ["rs1", "rs2", "rs3"],
                    "Z": [12.0, 10.0, 2.0],
                    "N": [100000.0, 100000.0, 100000.0],
                    "A1": ["A", "C", "G"],
                    "A2": ["G", "T", "A"],
                }
            )
        )
        dataset = runner.build_dataset(sumstats, self.make_ldscore_result())
        fitted_hsq = mock.Mock(
            coef=np.array([0.001]),
            tot=np.array([0.1]),
            intercept=np.array([1.0]),
        )
        expected_bins = pd.DataFrame(
            [{column: 0 for column in H2_REGRESSION_BIN_COLUMNS}]
        )
        with mock.patch.object(
            fitting, "summarize_ld_score_regression_bins", return_value=expected_bins
        ) as summarize_bins, mock.patch.object(
            fitting.reg, "Hsq", return_value=fitted_hsq
        ) as patched:
            result = runner.estimate_h2(dataset)
        chisq_arg = patched.call_args.args[0]
        self.assertEqual(chisq_arg.shape[0], 3)
        summarize_bins.assert_called_once()
        self.assertIs(result.diagnostic_bins, expected_bins)

    def test_rg_applies_default_two_step_when_single_annotation_and_unset(self):
        # Legacy estimate_rg (sumstats.py:400) sets two_step=30 for the rg fit
        # when there is a single annotation, no explicit cutoff, and a free h2
        # intercept. The restructured rg fit must pass the same default through.
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        sumstats_1 = self.make_sumstats_table()
        sumstats_2 = replace(self.make_sumstats_table(), trait_name="trait2")
        with mock.patch.object(fitting.reg, "RG", return_value=self.make_rg_kernel_result()) as patched:
            runner.estimate_rg(runner.build_rg_dataset(sumstats_1, sumstats_2, self.make_ldscore_result()))
        self.assertEqual(patched.call_args.kwargs["twostep"], 30)

    def test_summarize_total_h2_reports_post_filter_count(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        dataset = runner.build_dataset(self._high_chisq_sumstats(), self.make_ldscore_result())
        hsq = fake_total_h2_hsq()
        summary = fitting.summarize_total_h2(hsq, dataset, trait_name="trait", n_snps_used=2)
        self.assertEqual(int(summary.loc[0, "n_snps"]), 2)

    def _fake_total_h2_hsq(self, tot=0.2, tot_se=0.03):
        return fake_total_h2_hsq(tot=np.array([tot]), tot_se=np.array([tot_se]))

    def test_summarize_total_h2_observed_when_no_prevalence(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        dataset = runner.build_dataset(self._high_chisq_sumstats(), self.make_ldscore_result())
        summary = fitting.summarize_total_h2(
            self._fake_total_h2_hsq(), dataset, trait_name="t", n_snps_used=3
        )
        row = summary.iloc[0]
        self.assertEqual(row["total_h2_obs"], 0.2)
        self.assertEqual(row["total_h2_obs_se"], 0.03)
        self.assertTrue(math.isnan(row["total_h2_liab"]))
        self.assertTrue(math.isnan(row["total_h2_liab_se"]))
        self.assertTrue(math.isnan(row["samp_prev"]))
        self.assertTrue(math.isnan(row["pop_prev"]))

    def test_summarize_total_h2_liability_scaled(self):
        from ldsc._regression.estimators import liability_conversion_factor

        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        dataset = runner.build_dataset(self._high_chisq_sumstats(), self.make_ldscore_result())
        summary = fitting.summarize_total_h2(
            self._fake_total_h2_hsq(), dataset, trait_name="t", n_snps_used=3,
            samp_prev=0.5, pop_prev=0.01,
        )
        c = liability_conversion_factor(0.5, 0.01)
        row = summary.iloc[0]
        self.assertAlmostEqual(row["total_h2_liab"], 0.2 * c, places=12)
        self.assertAlmostEqual(row["total_h2_liab_se"], 0.03 * c, places=12)
        self.assertEqual(row["total_h2_obs"], 0.2)
        self.assertEqual(row["samp_prev"], 0.5)
        self.assertEqual(row["pop_prev"], 0.01)

    def test_partitioned_per_query_metadata_reports_post_filter_count_and_cap(self):
        # The real fit preparation selects the SNP population before the kernel stub.
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig(chisq_max=100.0))
        ldscore_result = self.make_ldscore_result()
        annotation_bundle = self.make_annotation_bundle()
        fake_hsq = fake_partitioned_hsq()
        # The toy LD columns are collinear; exercise selection independently.
        with mock.patch.object(fitting, "_raise_on_model_collinearity"), stub_h2_kernel(fake_hsq):
            result = runner.estimate_partitioned_h2_batch(self._high_chisq_sumstats(), ldscore_result, query_columns=annotation_bundle.query_columns, output_dir=self.temporary_output_dir(), summary_sort_by='category')
        meta = json.loads(result.per_query_artifacts["query1"].metadata.read_text())
        self.assertEqual(meta["n_snps"], 2)
        self.assertAlmostEqual(meta["effective_chisq_max"], 100.0)

    def test_rg_pair_count_and_metadata_reflect_product_filter(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig(chisq_max=1.0))
        t1 = self.make_sumstats_table()
        t2 = replace(self.make_sumstats_table(), trait_name="trait2")
        with stub_h2_kernel(self.make_rg_kernel_result().hsq1
        ), mock.patch.object(fitting.reg, "RG", return_value=self.make_rg_kernel_result()):
            result = runner.estimate_rg_pairs([t1, t2], self.make_ldscore_result())
        meta = result.per_pair_metadata[0]
        self.assertEqual(meta["n_snps_used"], 2)
        self.assertEqual(meta["n_blocks_used"], 2)
        self.assertAlmostEqual(meta["effective_chisq_max"], 1.0)
        self.assertEqual(int(result.rg_full.loc[0, "n_snps_used"]), 2)

    def test_estimate_rg_preserves_legacy_rg_array_call_contract(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        sumstats_1 = self.make_sumstats_table()
        sumstats_2 = replace(
            self.make_sumstats_table(),
            trait_name="trait2",
            data=self.make_sumstats_table().data.assign(Z=[1.0, 1.5, 0.25], N=[900.0, 900.0, 900.0]),
        )

        kernel_result = self.make_rg_kernel_result()
        with mock.patch.object(fitting.reg, "RG", return_value=kernel_result) as patched:
            result = runner.estimate_rg(runner.build_rg_dataset(sumstats_1, sumstats_2, self.make_ldscore_result()))

        self.assertIs(result.estimator, kernel_result)
        np.testing.assert_allclose(patched.call_args.args[0], [[2.0], [1.0], [0.5]])
        np.testing.assert_allclose(patched.call_args.args[1], [[1.0], [1.5], [0.25]])
        np.testing.assert_allclose(patched.call_args.args[2], [[1.0], [2.0], [3.0]])
        np.testing.assert_allclose(patched.call_args.args[3], [[2.0], [2.0], [2.0]])
        np.testing.assert_allclose(patched.call_args.args[4], [[1000.0], [1000.0], [1000.0]])
        np.testing.assert_allclose(patched.call_args.args[5], [[900.0], [900.0], [900.0]])
        np.testing.assert_allclose(patched.call_args.args[6], [[8.0]])

    def _base_only_partitioned_fixture(self):
        """Single-base-annotation dataset + fitted-Hsq mock for partitioned tests."""
        dataset = RegressionDataset(
            merged=pd.DataFrame(
                {
                    "SNP": ["rs1", "rs2"],
                    "Z": [1.0, 2.0],
                    "N": [1000.0, 1000.0],
                    "regression_ld_scores": [1.0, 1.0],
                }
            ),
            ref_ld_columns=["base"],
            weight_column="regression_ld_scores",
            reference_snp_count_totals={"common_reference_snp_counts": np.array([10.0])},
            count_key_used_for_regression="common_reference_snp_counts",
            retained_ld_columns=["base"],
            dropped_zero_variance_ld_columns=[],
            trait_names=["trait"],
            ldscore_overlap=LDScoreOverlap.from_contribution(
                OverlapContribution(
                    baseline_block_all=np.array([[10.0]]), baseline_block_common=np.array([[10.0]]),
                    query_diagonal_all=np.array([]), query_diagonal_common=np.array([]),
                    n_all=10, n_common=10,
                ),
                ["base"], [],
            ),
        )
        hsq = fake_partitioned_hsq(
            n_annot=1,
            intercept=np.array([1.01]),
            intercept_se=0.02,
            mean_chisq=np.array([1.2]),
            lambda_gc=np.array([1.1]),
            ratio=0.05,
            ratio_se=0.01,
        )
        return dataset, hsq

    def test_h2_and_partitioned_summaries_share_one_fitted_hsq_result(self):
        dataset, hsq = self._base_only_partitioned_fixture()

        total = fitting.summarize_total_h2(hsq, dataset, trait_name="trait")
        partitioned = fitting.summarize_partitioned_h2(hsq, dataset, ["base"])

        self.assertEqual(total.loc[0, "total_h2_obs"], partitioned.loc[0, "category_h2_obs"])
        self.assertEqual(total.loc[0, "total_h2_obs_se"], partitioned.loc[0, "category_h2_obs_se"])
        # The partitioned table also carries the whole-model total, equal to the h2 summary's.
        self.assertEqual(partitioned.loc[0, "total_h2_obs"], total.loc[0, "total_h2_obs"])
        self.assertEqual(partitioned.loc[0, "total_h2_obs_se"], total.loc[0, "total_h2_obs_se"])
        self.assertEqual(partitioned.loc[0, "category"], "base")
        self.assertEqual(partitioned.loc[0, "prop_h2"], 1.0)
        self.assertEqual(partitioned.loc[0, "enrichment"], 1.0)
        self.assertTrue(np.isnan(partitioned.loc[0, "enrichment_p"]))

    def test_summarize_partitioned_h2_both_scales(self):
        from ldsc._regression.estimators import liability_conversion_factor

        dataset, hsq = self._base_only_partitioned_fixture()
        obs = fitting.summarize_partitioned_h2(hsq, dataset, ["base"])
        liab = fitting.summarize_partitioned_h2(
            hsq, dataset, ["base"], samp_prev=0.5, pop_prev=0.01
        )
        c = liability_conversion_factor(0.5, 0.01)
        # Observed preserved; liability is c * observed; proportions are scale-invariant.
        self.assertEqual(liab.loc[0, "category_h2_obs"], obs.loc[0, "category_h2_obs"])
        self.assertAlmostEqual(
            liab.loc[0, "category_h2_liab"], obs.loc[0, "category_h2_obs"] * c, places=12
        )
        self.assertAlmostEqual(
            liab.loc[0, "category_h2_liab_se"],
            obs.loc[0, "category_h2_obs_se"] * c,
            places=12,
        )
        self.assertEqual(liab.loc[0, "prop_h2"], obs.loc[0, "prop_h2"])
        self.assertTrue(math.isnan(obs.loc[0, "category_h2_liab"]))
        self.assertEqual(liab.loc[0, "samp_prev"], 0.5)
        self.assertEqual(liab.loc[0, "pop_prev"], 0.01)

    def test_summarize_partitioned_h2_reports_total_h2(self):
        from ldsc._regression.estimators import liability_conversion_factor

        dataset, hsq = self._base_only_partitioned_fixture()  # tot=0.25, tot_se=0.03
        obs = fitting.summarize_partitioned_h2(hsq, dataset, ["base"])
        # Observed-scale total mirrors the fitted Hsq; liability is NaN without prevalences.
        self.assertEqual(obs.loc[0, "total_h2_obs"], 0.25)
        self.assertEqual(obs.loc[0, "total_h2_obs_se"], 0.03)
        self.assertTrue(math.isnan(obs.loc[0, "total_h2_liab"]))
        self.assertTrue(math.isnan(obs.loc[0, "total_h2_liab_se"]))
        # Liability total uses the same conversion factor as the category columns.
        liab = fitting.summarize_partitioned_h2(
            hsq, dataset, ["base"], samp_prev=0.5, pop_prev=0.01
        )
        c = liability_conversion_factor(0.5, 0.01)
        self.assertAlmostEqual(liab.loc[0, "total_h2_liab"], 0.25 * c, places=12)
        self.assertAlmostEqual(liab.loc[0, "total_h2_liab_se"], 0.03 * c, places=12)

    def test_summarize_rg_pair_liability_and_invariant_ratio(self):
        from ldsc._regression.estimators import liability_conversion_factor, gencov_obs_to_liab

        rg_result = self.make_rg_kernel_result(rg=0.4, rg_se=0.05)
        dataset = SimpleNamespace(merged=[0, 1, 2])
        obs_row = fitting._summarize_rg_pair(
            rg_result, dataset, trait_1="A", trait_2="B", pair_kind="all_pairs"
        )
        liab_row = fitting._summarize_rg_pair(
            rg_result, dataset, trait_1="A", trait_2="B", pair_kind="all_pairs",
            samp_prev_1=0.5, pop_prev_1=0.01, samp_prev_2=0.5, pop_prev_2=0.02,
        )
        c1 = liability_conversion_factor(0.5, 0.01)
        c2 = liability_conversion_factor(0.5, 0.02)
        # rg ratio (and its SE/z/p) is scale-invariant.
        self.assertEqual(liab_row["rg"], obs_row["rg"])
        self.assertEqual(liab_row["rg_se"], obs_row["rg_se"])
        # Per-trait h2 and the pair's gencov go to the liability scale.
        self.assertAlmostEqual(liab_row["h2_1_liab"], 0.20 * c1, places=12)
        self.assertAlmostEqual(liab_row["h2_1_liab_se"], 0.02 * c1, places=12)
        self.assertAlmostEqual(liab_row["h2_2_liab"], 0.30 * c2, places=12)
        self.assertAlmostEqual(
            liab_row["gencov_liab"], gencov_obs_to_liab(0.10, 0.5, 0.5, 0.01, 0.02), places=12
        )
        self.assertEqual(liab_row["h2_1_obs"], 0.20)
        self.assertEqual(liab_row["gencov_obs"], 0.10)
        # Observed-only row leaves liability columns NaN.
        self.assertTrue(math.isnan(obs_row["h2_1_liab"]))
        self.assertTrue(math.isnan(obs_row["gencov_liab"]))
        self.assertEqual(liab_row["samp_prev_1"], 0.5)
        self.assertEqual(liab_row["pop_prev_2"], 0.02)

    def test_collinear_partitioned_model_raises_input_error(self):
        # Two near-duplicate baseline annotations produce a singular design
        # matrix; the fit must abort (legacy hard stop), not warn-and-proceed.
        overlap = LDScoreOverlap.from_contribution(
            OverlapContribution(
                baseline_block_all=np.array([[10.0, 9.999], [9.999, 10.0]]),
                baseline_block_common=np.array([[10.0, 9.999], [9.999, 10.0]]),
                query_diagonal_all=np.array([]),
                query_diagonal_common=np.array([]),
                n_all=10,
                n_common=10,
            ),
            ["catA", "catA_dup"],
            [],
        )
        dataset = RegressionDataset(
            merged=pd.DataFrame(
                {"SNP": ["rs1"], "Z": [1.0], "N": [1000.0], "regression_ld_scores": [1.0]}
            ),
            ref_ld_columns=["catA", "catA_dup"],
            weight_column="regression_ld_scores",
            reference_snp_count_totals={"common_reference_snp_counts": np.array([10.0, 10.0])},
            count_key_used_for_regression="common_reference_snp_counts",
            retained_ld_columns=["catA", "catA_dup"],
            dropped_zero_variance_ld_columns=[],
            trait_names=["trait"],
            ldscore_overlap=overlap,
        )
        x = np.array([[1.0, 1.0001], [2.0, 2.0001], [3.0, 3.0001], [4.0, 4.0002]])
        with self.assertRaises(LDSCInputError):
            fitting._raise_on_model_collinearity(dataset, x)

    def test_partitioned_full_summary_uses_canonical_column_order_and_enrichment_p(self):
        dataset = RegressionDataset(
            merged=pd.DataFrame(
                {
                    "SNP": ["rs1", "rs2", "rs3"],
                    "Z": [1.0, 1.5, 2.0],
                    "N": [1000.0, 1000.0, 1000.0],
                    "regression_ld_scores": [1.0, 1.0, 1.0],
                }
            ),
            ref_ld_columns=["base", "query"],
            weight_column="regression_ld_scores",
            reference_snp_count_totals={"common_reference_snp_counts": np.array([10.0, 30.0])},
            count_key_used_for_regression="common_reference_snp_counts",
            retained_ld_columns=["base", "query"],
            dropped_zero_variance_ld_columns=[],
            trait_names=["trait"],
            ldscore_overlap=LDScoreOverlap.from_contribution(
                OverlapContribution(
                    baseline_block_all=np.array([[10.0, 8.0]]), baseline_block_common=np.array([[10.0, 8.0]]),
                    query_diagonal_all=np.array([30.0]), query_diagonal_common=np.array([30.0]),
                    n_all=40, n_common=40,
                ),
                ["base"], ["query"],
            ),
        )
        hsq = mock.Mock(
            tot=np.array([0.8]),
            tot_se=np.array([0.08]),
            coef=np.array([0.2, 0.5]),
            coef_cov=np.diag([0.01, 0.04]),
            coef_se=np.array([0.1, 0.2]),
            cat=np.array([0.2, 0.6]),
            cat_se=np.array([0.02, 0.06]),
            prop=np.array([0.25, 0.75]),
            prop_cov=np.diag([0.0009, 0.0081]),
            prop_se=np.array([0.03, 0.09]),
            enrichment=np.array([1.0, 1.0]),
            n_blocks=200,
            n_annot=2,
        )

        summary = fitting.summarize_partitioned_h2(hsq, dataset, ["query"])

        self.assertEqual(summary.columns.tolist(), PARTITIONED_H2_COLUMNS)
        self.assertEqual(summary.loc[0, "category"], "query")
        # Overlap-aware Prop._SNPs = M_query / M_tot on the common universe (30 / 40).
        self.assertAlmostEqual(summary.loc[0, "prop_snps"], 30.0 / 40.0)
        # One-sided coefficient_p = sf(coef/se) = sf(0.5 / 0.2).
        self.assertAlmostEqual(summary.loc[0, "coefficient_p"], stats.norm.sf(2.5))
        self.assertTrue(bool(summary.loc[0, "overlap_annot"]))
        self.assertTrue(np.isfinite(summary.loc[0, "enrichment_p"]))
