"""Partitioned-h2 batches and the single-query facade method: shared preparation, regimes, summary sorting."""

from dataclasses import replace
import json
from pathlib import Path
import unittest
from unittest import mock

import numpy as np
import pandas as pd
from scipy import stats

from ldsc.config import GlobalConfig, RegressionConfig
from ldsc.errors import LDSCInputError
from ldsc._regression import fitting, partitioned
from ldsc._regression.outputs import PARTITIONED_H2_COLUMNS
from ldsc._regression.results import PartitionedH2BatchResult
from ldsc._regression.runner import RegressionRunner
from tests.regression_fixtures import RegressionFixtures, stub_h2_kernel


class PartitionedTest(RegressionFixtures, unittest.TestCase):
    def test_partitioned_functional_regime_reports_aggregate_count_and_cap(self):
        # Functional regime applies the same chi-square filter as h2; record the
        # post-filter count and effective cap in the aggregate metadata.
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig(chisq_max=4.0))
        ldscore_result = self._baseline_only_ldscore_result()
        fake_hsq = mock.Mock(
            tot=np.array([0.25]), tot_se=np.array([0.03]),
            coef=np.array([0.025]), coef_cov=np.array([[9e-6]]), coef_se=np.array([0.003]),
            cat=np.array([0.25]), cat_se=np.array([0.03]),
            prop=np.array([1.0]), prop_cov=np.array([[0.0]]), prop_se=np.array([0.0]),
            enrichment=np.array([1.0]), n_blocks=200, n_annot=1,
        )
        with stub_h2_kernel(fake_hsq):
            result = runner.estimate_partitioned_h2_batch(self._high_chisq_sumstats(), ldscore_result, query_columns=replace(self.make_annotation_bundle(), query_columns=[]).query_columns, output_dir=self.temporary_output_dir(), summary_sort_by='category')
        self.assertEqual(json.loads(Path(result.output_paths["metadata"]).read_text())["n_snps"], 1)
        self.assertAlmostEqual(json.loads(Path(result.output_paths["metadata"]).read_text())["effective_chisq_max"], 4.0)

    def test_partitioned_functional_single_annotation_warns_degenerate(self):
        # A functional fit over a single baseline annotation is degenerate; warn.
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        ldscore_result = self._baseline_only_ldscore_result()
        fake_hsq = mock.Mock(
            tot=np.array([0.25]), tot_se=np.array([0.03]),
            coef=np.array([0.025]), coef_cov=np.array([[9e-6]]), coef_se=np.array([0.003]),
            cat=np.array([0.25]), cat_se=np.array([0.03]),
            prop=np.array([1.0]), prop_cov=np.array([[0.0]]), prop_se=np.array([0.0]),
            enrichment=np.array([1.0]), n_blocks=200, n_annot=1,
        )
        with stub_h2_kernel(fake_hsq), self.assertLogs(
            "LDSC.regression_runner", level="WARNING"
        ) as logs:
            runner.estimate_partitioned_h2_batch(self.make_sumstats_table(), ldscore_result, query_columns=replace(self.make_annotation_bundle(), query_columns=[]).query_columns, output_dir=self.temporary_output_dir(), summary_sort_by='category')
        self.assertTrue(any("degenerate" in m.lower() for m in logs.output))

    def test_partitioned_cell_type_regime_does_not_warn_degenerate(self):
        # Cell-type models always fit baseline + one query (n_annot >= 2); no warning.
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        fake_hsq = mock.Mock(
            tot=np.array([0.3]), tot_se=np.array([0.03]),
            coef=np.array([0.0, 1.0]), coef_cov=np.diag([0.01, 0.04]), coef_se=np.array([0.1, 0.2]),
            cat=np.array([0.0, 0.3]), cat_se=np.array([0.01, 0.03]),
            prop=np.array([0.0, 0.4]), prop_cov=np.diag([0.0001, 0.0016]), prop_se=np.array([0.01, 0.04]),
            enrichment=np.array([0.0, 2.0]), n_blocks=200, n_annot=2,
        )
        # The toy LD columns are collinear; exercise selection independently.
        with (
            mock.patch.object(fitting, "_raise_on_model_collinearity"),
            mock.patch.object(partitioned.LOGGER, "warning") as warn,
            stub_h2_kernel(fake_hsq),
        ):
            runner.estimate_partitioned_h2_batch(self.make_sumstats_table(), self.make_ldscore_result(), query_columns=self.make_annotation_bundle().query_columns, output_dir=self.temporary_output_dir(), summary_sort_by='category')
        self.assertFalse(any("degenerate" in str(call).lower() for call in warn.call_args_list))

    def _baseline_only_ldscore_result(self):
        from ldsc._annotation.overlap import OverlapContribution
        from ldsc._ldscore.overlap_table import LDScoreOverlap

        overlap = LDScoreOverlap.from_contribution(
            OverlapContribution(
                baseline_block_all=np.array([[10.0]]),
                baseline_block_common=np.array([[8.0]]),
                query_diagonal_all=np.array([]),
                query_diagonal_common=np.array([]),
                n_all=10, n_common=8,
            ),
            baseline_columns=["base"], query_columns=[],
        )
        return replace(
            self.make_ldscore_result(),
            query_table=None,
            query_columns=[],
            count_records=[
                {"group": "baseline", "column": "base", "all_reference_snp_count": 10.0, "common_reference_snp_count": 8.0}
            ],
            overlap=overlap,
        )

    def test_partitioned_h2_without_overlap_sidecar_fails_fast(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        ldscore_result = replace(self.make_ldscore_result(), overlap=None)
        with self.assertRaisesRegex(LDSCInputError, "overlap matrix"):
            runner.estimate_partitioned_h2_batch(self.make_sumstats_table(), ldscore_result, query_columns=self.make_annotation_bundle().query_columns, output_dir=self.temporary_output_dir(), summary_sort_by='category')

    def test_summary_sort_auto_is_regime_aware(self):
        from ldsc._regression.partitioned import _resolve_summary_sort

        self.assertEqual(_resolve_summary_sort("auto", has_queries=True), "coefficient-p")
        self.assertEqual(_resolve_summary_sort("auto", has_queries=False), "category")
        self.assertEqual(_resolve_summary_sort("enrichment", has_queries=True), "enrichment")

    def test_estimate_partitioned_h2_requires_explicit_query_column(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())

        with self.assertRaises(TypeError):
            runner.estimate_partitioned_h2(
                self.make_sumstats_table(),
                self.make_ldscore_result(),
                self.make_annotation_bundle(),
            )

    def test_estimate_partitioned_h2_uses_requested_query_column(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        fake_hsq = mock.Mock(
            tot=np.array([0.3]),
            tot_se=np.array([0.03]),
            coef=np.array([0.0, 1.0]),
            coef_cov=np.diag([0.01, 0.04]),
            coef_se=np.array([0.1, 0.2]),
            cat=np.array([0.0, 0.3]),
            cat_se=np.array([0.01, 0.03]),
            prop=np.array([0.0, 0.4]),
            prop_cov=np.diag([0.0001, 0.0016]),
            prop_se=np.array([0.01, 0.04]),
            enrichment=np.array([0.0, 2.0]),
            n_blocks=200,
            n_annot=2,
        )

        # The toy LD columns are collinear; exercise selection independently.
        with mock.patch.object(fitting, "_raise_on_model_collinearity"), mock.patch.object(
            runner,
            "build_dataset",
            wraps=runner.build_dataset,
        ) as build_dataset, stub_h2_kernel(fake_hsq):
            result = runner.estimate_partitioned_h2(
                self.make_sumstats_table(),
                self.make_ldscore_result(),
                query_column="query1",
            )

        self.assertEqual(build_dataset.call_args.kwargs["query_columns"], ["query1"])
        self.assertEqual(result["category"].tolist(), ["query1"])

    def test_estimate_partitioned_h2_rejects_unknown_query_column(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())

        with self.assertRaisesRegex(LDSCInputError, "requested query annotation columns.*missing"):
            runner.estimate_partitioned_h2(
                self.make_sumstats_table(),
                self.make_ldscore_result(),
                query_column="missing",
            )

    def test_estimate_partitioned_h2_batch_loops_over_queries(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        table = self.make_sumstats_table()
        ldscore_result = self.make_ldscore_result()
        annotation_bundle = self.make_annotation_bundle()
        fake_hsq = mock.Mock(
            tot=np.array([0.3]),
            tot_se=np.array([0.03]),
            coef=np.array([0.0, 1.0]),
            coef_cov=np.diag([0.01, 0.04]),
            coef_se=np.array([0.1, 0.2]),
            cat=np.array([0.0, 0.3]),
            cat_se=np.array([0.01, 0.03]),
            prop=np.array([0.0, 0.4]),
            prop_cov=np.diag([0.0001, 0.0016]),
            prop_se=np.array([0.01, 0.04]),
            enrichment=np.array([0.0, 2.0]),
            n_blocks=200,
            n_annot=2,
        )
        # The toy LD columns are collinear; exercise selection independently.
        with (
            mock.patch.object(fitting, "_raise_on_model_collinearity"),
            self.assertLogs("LDSC.regression_runner", level="INFO") as captured,
            stub_h2_kernel(fake_hsq) as patched,
        ):
            result = runner.estimate_partitioned_h2_batch(table, ldscore_result, query_columns=annotation_bundle.query_columns, output_dir=self.temporary_output_dir(), summary_sort_by='category')
        result = result.summary
        self.assertEqual(patched.call_count, 2)
        phase_lines = "\n".join(captured.output)
        self.assertEqual(phase_lines.count("Phase timing: shared regression preparation completed"), 1)
        self.assertNotIn("Phase timing: estimator execution completed", phase_lines)
        self.assertEqual(result.columns.tolist(), PARTITIONED_H2_COLUMNS)
        self.assertEqual(result["category"].tolist(), ["query1", "query2"])
        # Overlap-aware Prop._SNPs = M_query1 / M_tot on the common universe (18 / 45).
        self.assertAlmostEqual(result.loc[0, "prop_snps"], 18.0 / 45.0)
        # One-sided Coefficient_p = P(Z > coef/se) = sf(1.0 / 0.2).
        self.assertAlmostEqual(result.loc[0, "coefficient_p"], stats.norm.sf(5.0))

    def test_estimate_partitioned_h2_batch_can_return_full_partitioned_tables(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        table = self.make_sumstats_table()
        ldscore_result = self.make_ldscore_result()
        annotation_bundle = self.make_annotation_bundle()
        fake_hsq = mock.Mock(
            tot=np.array([0.3]),
            tot_se=np.array([0.03]),
            coef=np.array([0.0, 1.0]),
            coef_cov=np.diag([0.01, 0.04]),
            coef_se=np.array([0.1, 0.2]),
            cat=np.array([0.0, 0.3]),
            cat_se=np.array([0.01, 0.03]),
            prop=np.array([0.0, 0.4]),
            prop_cov=np.diag([0.0001, 0.0016]),
            prop_se=np.array([0.01, 0.04]),
            enrichment=np.array([0.0, 2.0]),
            n_blocks=200,
            n_annot=2,
        )
        # The toy LD columns are collinear; exercise selection independently.
        with mock.patch.object(fitting, "_raise_on_model_collinearity"), stub_h2_kernel(fake_hsq):
            result = runner.estimate_partitioned_h2_batch(table, ldscore_result, query_columns=annotation_bundle.query_columns, output_dir=self.temporary_output_dir(), summary_sort_by='category')

        self.assertIsInstance(result, PartitionedH2BatchResult)
        self.assertEqual(result.summary["category"].tolist(), ["query1", "query2"])
        self.assertIn("query1", result.per_query_artifacts)
        self.assertEqual(json.loads(result.per_query_artifacts["query1"].metadata.read_text())["effective_snp_identifier"], "rsid")
        self.assertFalse(json.loads(result.per_query_artifacts["query1"].metadata.read_text())["identity_downgrade_applied"])
        # Per-query full table holds all baseline+query model rows in the single schema.
        self.assertEqual(pd.read_csv(result.per_query_artifacts["query1"].full, sep="\t")["category"].tolist(), ["base", "query1"])
        self.assertEqual(
            pd.read_csv(result.per_query_artifacts["query1"].full, sep="\t").columns.tolist(),
            PARTITIONED_H2_COLUMNS,
        )
        # One fitted model -> one total h2, broadcast across that model's rows.
        self.assertEqual(pd.read_csv(result.per_query_artifacts["query1"].full, sep="\t")["total_h2_obs"].tolist(), [0.3, 0.3])

    def test_sort_partitioned_h2_summary_preserves_category_order(self):
        summary = pd.DataFrame(
            [
                {"category": "query_b", "enrichment": 1.0},
                {"category": "query_a", "enrichment": 5.0},
                {"category": "query_c", "enrichment": 2.0},
            ]
        )

        sorted_summary = partitioned._sort_partitioned_h2_summary(summary, "category")

        self.assertEqual(sorted_summary["category"].tolist(), ["query_b", "query_a", "query_c"])

    def test_sort_partitioned_h2_summary_sorts_p_values_ascending_with_nan_last(self):
        summary = pd.DataFrame(
            [
                {"category": "query_b", "enrichment_p": 0.20},
                {"category": "query_a", "enrichment_p": np.nan},
                {"category": "query_c", "enrichment_p": 0.01},
                {"category": "query_d", "enrichment_p": 0.20},
            ]
        )

        sorted_summary = partitioned._sort_partitioned_h2_summary(summary, "enrichment-p")

        self.assertEqual(sorted_summary["category"].tolist(), ["query_c", "query_b", "query_d", "query_a"])

    def test_sort_partitioned_h2_summary_sorts_effects_descending_with_nan_last(self):
        summary = pd.DataFrame(
            [
                {"category": "query_b", "enrichment": 1.0, "coefficient": 0.3},
                {"category": "query_a", "enrichment": np.nan, "coefficient": 0.8},
                {"category": "query_c", "enrichment": 5.0, "coefficient": np.nan},
                {"category": "query_d", "enrichment": 5.0, "coefficient": 0.8},
            ]
        )

        enrichment = partitioned._sort_partitioned_h2_summary(summary, "enrichment")
        coefficient = partitioned._sort_partitioned_h2_summary(summary, "coefficient")

        self.assertEqual(enrichment["category"].tolist(), ["query_c", "query_d", "query_b", "query_a"])
        self.assertEqual(coefficient["category"].tolist(), ["query_a", "query_d", "query_b", "query_c"])
