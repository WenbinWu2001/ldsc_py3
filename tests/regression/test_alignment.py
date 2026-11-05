"""Aligning summary statistics with LD scores into h2 and rg regression datasets."""

from dataclasses import replace
from types import SimpleNamespace
import unittest
from unittest import mock

import numpy as np
import pandas as pd

from ldsc.config import GlobalConfig, RegressionConfig
from ldsc.errors import LDSCInputError, LDSCUsageError
from ldsc._regression import alignment
from ldsc._regression.runner import RegressionRunner
from ldsc._sumstats.table import SumstatsTable
from tests.regression_fixtures import RegressionFixtures


class AlignmentTest(RegressionFixtures, unittest.TestCase):
    def test_build_dataset_uses_baseline_only_for_h2_style_runs(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        dataset = runner.build_dataset(self.make_sumstats_table(), self.make_ldscore_result())
        self.assertEqual(dataset.count_key_used_for_regression, "common_reference_snp_counts")
        self.assertEqual(dataset.weight_column, "regression_ld_scores")
        self.assertIn("regression_ld_scores", dataset.merged.columns)
        self.assertEqual(dataset.retained_ld_columns, ["base"])
        self.assertEqual(dataset.dropped_zero_variance_ld_columns, [])
        np.testing.assert_allclose(
            dataset.reference_snp_count_totals["common_reference_snp_counts"],
            [8.0],
        )
        self.assertEqual(dataset.config_snapshot, GlobalConfig(snp_identifier="rsid"))

    def test_build_dataset_baseline_only_ignores_unused_misaligned_query_table(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        result = self.make_ldscore_result()
        bad_query = result.query_table.iloc[[1, 0, 2]].reset_index(drop=True)
        result = replace(result, query_table=bad_query)

        dataset = runner.build_dataset(self.make_sumstats_table(), result)

        self.assertEqual(dataset.retained_ld_columns, ["base"])
        self.assertEqual(dataset.merged["SNP"].tolist(), ["rs1", "rs2", "rs3"])

    def test_build_dataset_orders_merged_rows_by_ld_score_genomic_order(self):
        # Gap C: jackknife block contiguity depends on the merged row order. The
        # LD-score artifact is genomically sorted, so merged rows must follow the
        # LD-score order even when the sumstats rows are shuffled (e.g. an unsorted
        # legacy .sumstats.gz). The order originates entirely from the LD-score side.
        ld_frame = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "SNP": ["rs_a", "rs_b", "rs_c"],
                "POS": [10, 20, 30],
                "regression_ld_scores": [2.0, 2.0, 2.0],
                "base": [1.0, 2.0, 3.0],
            }
        )
        sumstats = self.make_sumstats_table_from_frame(
            pd.DataFrame(
                {
                    "SNP": ["rs_c", "rs_a", "rs_b"],  # deliberately not genomic order
                    "Z": [0.5, 2.0, 1.0],
                    "N": [1000.0, 1000.0, 1000.0],
                }
            )
        )
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        dataset = runner.build_dataset(sumstats, self.make_ldscore_result_from_frame(ld_frame))
        self.assertEqual(dataset.merged["SNP"].tolist(), ["rs_a", "rs_b", "rs_c"])

    def test_build_dataset_can_include_one_query_annotation_for_partitioned_h2(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        dataset = runner.build_dataset(self.make_sumstats_table(), self.make_ldscore_result(), query_columns=["query2"])
        self.assertEqual(dataset.retained_ld_columns, ["base", "query2"])
        np.testing.assert_allclose(
            dataset.reference_snp_count_totals["common_reference_snp_counts"],
            [8.0, 28.0],
        )

    def test_build_dataset_allele_aware_mode_merges_by_effective_identity(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid_allele_aware"), RegressionConfig())
        sumstats = self.make_sumstats_table_from_frame(
            pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "SNP": ["rs1", "rs2"],
                    "POS": [10, 20],
                    "A1": ["A", "A"],
                    "A2": ["C", "C"],
                    "Z": [2.0, 3.0],
                    "N": [1000.0, 1000.0],
                }
            ),
            snp_identifier="rsid_allele_aware",
        )
        ldscore_result = self.make_ldscore_result_from_frame(
            pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "SNP": ["rs1", "rs2"],
                    "POS": [10, 20],
                    "A1": ["A", "A"],
                    "A2": ["C", "G"],
                    "regression_ld_scores": [2.0, 2.0],
                    "base": [1.0, 2.0],
                }
            ),
            snp_identifier="rsid_allele_aware",
        )

        dataset = runner.build_dataset(sumstats, ldscore_result)

        self.assertEqual(dataset.merged["SNP"].tolist(), ["rs1"])
        self.assertEqual(dataset.effective_snp_identifier, "rsid_allele_aware")
        self.assertFalse(dataset.identity_downgrade_applied)

    def test_build_dataset_cross_family_modes_reject_even_with_identity_downgrade(self):
        runner = RegressionRunner(
            GlobalConfig(snp_identifier="rsid_allele_aware", genome_build="hg38"),
            RegressionConfig(allow_identity_downgrade=True),
        )
        sumstats = self.make_sumstats_table_from_frame(
            pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "POS": [10],
                    "A1": ["A"],
                    "A2": ["C"],
                    "Z": [2.0],
                    "N": [1000.0],
                }
            ),
            snp_identifier="rsid_allele_aware",
        )
        ldscore_result = self.make_ldscore_result_from_frame(
            pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "POS": [10],
                    "regression_ld_scores": [2.0],
                    "base": [1.0],
                }
            ),
            snp_identifier="chr_pos",
        )

        with self.assertRaisesRegex(LDSCUsageError, "different families"):
            runner.build_dataset(sumstats, ldscore_result)

    def test_build_dataset_rejects_same_family_mode_mismatch_by_default(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="chr_pos_allele_aware"), RegressionConfig())
        sumstats = self.make_sumstats_table_from_frame(
            pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "POS": [10],
                    "A1": ["A"],
                    "A2": ["C"],
                    "Z": [2.0],
                    "N": [1000.0],
                }
            ),
            snp_identifier="chr_pos_allele_aware",
        )
        ldscore_result = self.make_ldscore_result_from_frame(
            pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["ld1"],
                    "POS": [10],
                    "regression_ld_scores": [2.0],
                    "base": [1.0],
                }
            ),
            snp_identifier="chr_pos",
        )

        with self.assertRaisesRegex(LDSCUsageError, "allow-identity-downgrade"):
            runner.build_dataset(sumstats, ldscore_result)

    def test_build_dataset_allows_chr_pos_identity_downgrade_with_effective_base_mode(self):
        runner = RegressionRunner(
            GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"),
            RegressionConfig(allow_identity_downgrade=True),
        )
        sumstats = self.make_sumstats_table_from_frame(
            pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "SNP": ["sum1", "sum2"],
                    "POS": [10, 20],
                    "A1": ["A", "A"],
                    "A2": ["C", "G"],
                    "Z": [2.0, 3.0],
                    "N": [1000.0, 1000.0],
                }
            ),
            snp_identifier="chr_pos_allele_aware",
        )
        ldscore_result = self.make_ldscore_result_from_frame(
            pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "SNP": ["ld1", "ld2"],
                    "POS": [10, 20],
                    "regression_ld_scores": [2.0, 2.0],
                    "base": [1.0, 2.0],
                }
            ),
            snp_identifier="chr_pos",
        )

        dataset = runner.build_dataset(sumstats, ldscore_result)

        self.assertEqual(dataset.merged["SNP"].tolist(), ["sum1", "sum2"])
        self.assertEqual(dataset.effective_snp_identifier, "chr_pos")
        self.assertTrue(dataset.identity_downgrade_applied)

    def test_build_dataset_allows_rsid_identity_downgrade_with_effective_base_mode(self):
        runner = RegressionRunner(
            GlobalConfig(snp_identifier="rsid_allele_aware", genome_build="hg38"),
            RegressionConfig(allow_identity_downgrade=True),
        )
        sumstats = self.make_sumstats_table_from_frame(
            pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "SNP": ["rs1", "rs2"],
                    "POS": [10, 20],
                    "A1": ["A", "A"],
                    "A2": ["C", "G"],
                    "Z": [2.0, 3.0],
                    "N": [1000.0, 1000.0],
                }
            ),
            snp_identifier="rsid_allele_aware",
        )
        ldscore_result = self.make_ldscore_result_from_frame(
            pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "SNP": ["rs1", "rs2"],
                    "POS": [100, 200],
                    "regression_ld_scores": [2.0, 2.0],
                    "base": [1.0, 2.0],
                }
            ),
            snp_identifier="rsid",
        )

        dataset = runner.build_dataset(sumstats, ldscore_result)

        self.assertEqual(dataset.merged["SNP"].tolist(), ["rs1", "rs2"])
        self.assertEqual(dataset.effective_snp_identifier, "rsid")
        self.assertTrue(dataset.identity_downgrade_applied)

    def test_identity_downgrade_drops_duplicate_effective_keys_before_merge_and_logs(self):
        runner = RegressionRunner(
            GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"),
            RegressionConfig(allow_identity_downgrade=True),
        )
        sumstats = self.make_sumstats_table_from_frame(
            pd.DataFrame(
                {
                    "CHR": ["1", "1", "1"],
                    "SNP": ["dup_a", "dup_b", "keep"],
                    "POS": [10, 10, 20],
                    "A1": ["A", "A", "A"],
                    "A2": ["C", "G", "C"],
                    "Z": [9.0, 8.0, 2.0],
                    "N": [1000.0, 1000.0, 1000.0],
                }
            ),
            snp_identifier="chr_pos_allele_aware",
        )
        ldscore_result = self.make_ldscore_result_from_frame(
            pd.DataFrame(
                {
                    "CHR": ["1", "1", "1"],
                    "SNP": ["keep_ld", "dup_ld_a", "dup_ld_b"],
                    "POS": [20, 30, 30],
                    "regression_ld_scores": [2.0, 2.0, 2.0],
                    "base": [1.0, 3.0, 4.0],
                }
            ),
            snp_identifier="chr_pos",
        )

        with self.assertLogs("LDSC.regression_runner", level="WARNING") as caught:
            dataset = runner.build_dataset(sumstats, ldscore_result)

        self.assertEqual(dataset.merged["SNP"].tolist(), ["keep"])
        self.assertEqual(dataset.merged["base"].tolist(), [1.0])
        log_text = "\n".join(caught.output)
        self.assertIn("Identity downgrade enabled", log_text)
        self.assertIn("LD-score mode chr_pos", log_text)
        self.assertIn("summary-statistics mode chr_pos_allele_aware", log_text)
        self.assertIn("effective snp_identifier='chr_pos'", log_text)
        self.assertIn("Dropped 4 duplicate effective-key rows", log_text)

    def _log_identity_of_prepared_h2(self, sumstats, ldscore_result, runner_config, config):
        """Prepare an h2 dataset and log its identity, as ``ldsc h2`` does."""
        dataset = alignment._prepare_h2_dataset(sumstats, ldscore_result, config, runner_config).dataset
        with self.assertLogs("LDSC.regression_runner", level="INFO") as caught:
            alignment._log_regression_identity(
                [sumstats], ldscore_result, runner_config,
                [(dataset.effective_snp_identifier, dataset.identity_downgrade_applied)],
            )
        return "\n".join(caught.output)

    def _identity_inputs(self, *, ldscore_mode, sumstats_mode):
        ldscore_result = self.make_ldscore_result_from_frame(
            pd.DataFrame({"CHR": ["1"], "SNP": ["rs1"], "POS": [10], "regression_ld_scores": [2.0], "base": [1.0]}),
            snp_identifier=ldscore_mode,
        )
        sumstats = self.make_sumstats_table_from_frame(
            pd.DataFrame({"CHR": ["1"], "SNP": ["rs1"], "POS": [10], "A1": ["A"], "A2": ["C"],
                          "Z": [2.0], "N": [1000.0]}),
            snp_identifier=sumstats_mode,
        )
        return ldscore_result, sumstats

    def test_log_regression_identity_reports_provenance_over_default(self):
        # GlobalConfig default differs from what the artifacts actually recorded;
        # the log must show the provenance-derived effective mode, not the default.
        runner_config = GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38")
        ldscore_result, sumstats = self._identity_inputs(ldscore_mode="chr_pos", sumstats_mode="chr_pos")

        log_text = self._log_identity_of_prepared_h2(sumstats, ldscore_result, runner_config, RegressionConfig())

        self.assertIn("SNP identifier used for this regression: chr_pos.", log_text)
        self.assertIn("'chr_pos_allele_aware' shown above was NOT used", log_text)
        self.assertIn("LD-score=chr_pos (from provenance)", log_text)
        self.assertNotIn("downgrade", log_text.lower())

    def test_log_regression_identity_notes_allele_downgrade(self):
        runner_config = GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38")
        ldscore_result, sumstats = self._identity_inputs(ldscore_mode="chr_pos", sumstats_mode="chr_pos_allele_aware")

        log_text = self._log_identity_of_prepared_h2(
            sumstats, ldscore_result, runner_config, RegressionConfig(allow_identity_downgrade=True))

        self.assertIn("SNP identifier used for this regression: chr_pos. An allele-awareness downgrade", log_text)

    def test_log_regression_identity_names_the_fallback_each_input_without_provenance_takes(self):
        # LD scores without provenance take the GlobalConfig default; summary
        # statistics without provenance take the LD-score directory's mode.
        runner_config = GlobalConfig(snp_identifier="rsid")
        ldscore_result, sumstats = self._identity_inputs(ldscore_mode="rsid", sumstats_mode="rsid")
        ldscore_result = replace(ldscore_result, config_snapshot=None)
        sumstats = replace(sumstats, config_snapshot=None)

        log_text = self._log_identity_of_prepared_h2(sumstats, ldscore_result, runner_config, RegressionConfig())

        self.assertIn("SNP identifier used for this regression: rsid.", log_text)
        self.assertIn("'rsid' shown above was used only as a fallback", log_text)
        self.assertIn("LD-score=rsid (no provenance recorded; using GlobalConfig default)", log_text)
        self.assertIn("=rsid (no provenance recorded; using the LD-score directory's mode)", log_text)

    def test_log_regression_identity_does_not_say_when_an_unrecorded_pair_failed(self):
        # A pair can fail before alignment or after its fit, and its metadata
        # then lacks the realized mode either way, so the line must not claim
        # which stage failed.
        runner_config = GlobalConfig(snp_identifier="chr_pos")
        ldscore_result, sumstats = self._identity_inputs(ldscore_mode="chr_pos", sumstats_mode="chr_pos")

        with self.assertLogs("LDSC.regression_runner", level="INFO") as caught:
            alignment._log_regression_identity(
                [sumstats, sumstats], ldscore_result, runner_config, [("chr_pos", False)], unresolved_pairs=2)
        log_text = "\n".join(caught.output)

        self.assertIn("chr_pos for 1 pair; not recorded for 2 pairs that failed.", log_text)
        self.assertNotIn("before alignment", log_text)

    def test_h2_keeps_sumstats_z_in_the_effect_allele_orientation(self):
        # h2 fits chi^2 = Z^2, so the sign of Z is never aligned to LD-score alleles.
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid_allele_aware"), RegressionConfig())
        sumstats = self.make_sumstats_table_from_frame(
            pd.DataFrame(
                {
                    "CHR": ["1", "1", "1", "1"],
                    "SNP": ["rs1", "rs2", "rs3", "rs4"],
                    "POS": [10, 20, 30, 40],
                    "A1": ["A", "C", "A", "C"],
                    "A2": ["C", "A", "C", "A"],
                    "Z": [2.0, 3.0, 4.0, 5.0],
                    "N": [1000.0, 1000.0, 1000.0, 1000.0],
                }
            ),
            snp_identifier="rsid_allele_aware",
        )
        ldscore_result = self.make_ldscore_result_from_frame(
            pd.DataFrame(
                {
                    "CHR": ["1", "1", "1", "1"],
                    "SNP": ["rs1", "rs2", "rs3", "rs4"],
                    "POS": [10, 20, 30, 40],
                    "A1": ["A", "A", "T", "T"],
                    "A2": ["C", "C", "G", "G"],
                    "regression_ld_scores": [2.0, 2.0, 2.0, 2.0],
                    "base": [1.0, 2.0, 3.0, 4.0],
                }
            ),
            snp_identifier="rsid_allele_aware",
        )

        dataset = runner.build_dataset(sumstats, ldscore_result)

        self.assertEqual(dataset.merged["SNP"].tolist(), ["rs1", "rs2", "rs3", "rs4"])
        np.testing.assert_allclose(dataset.merged["Z"], [2.0, 3.0, 4.0, 5.0])

    def test_build_dataset_empty_merge_error_includes_active_config(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        sumstats = replace(
            self.make_sumstats_table(),
            data=self.make_sumstats_table().data.assign(SNP=["missing1", "missing2", "missing3"]),
        )

        with self.assertRaisesRegex(LDSCInputError, "retained no overlapping.*Other causes"):
            runner.build_dataset(sumstats, self.make_ldscore_result())

    def test_build_dataset_with_query_rejects_misaligned_query_table(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        result = self.make_ldscore_result()
        bad_query = result.query_table.copy()
        bad_query.loc[1, "POS"] = 999
        result = replace(result, query_table=bad_query)

        with self.assertRaisesRegex(LDSCInputError, "query rows must match baseline rows"):
            runner.build_dataset(self.make_sumstats_table(), result, query_columns=["query2"])

    def test_missing_requested_common_counts_fail_without_all_count_substitution(self):
        old_result = replace(
            self.make_ldscore_result(),
            count_records=[
                {
                    "group": "baseline",
                    "column": "base",
                    "all_reference_snp_count": 10.0,
                    "common_reference_snp_count_maf_gt_0_05": 8.0,
                }
            ],
            query_table=None,
            query_columns=[],
        )
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())

        # The old metadata key is not a common-SNP count, and all-SNP counts never stand in for it.
        with self.assertRaisesRegex(LDSCInputError, r"common-SNP count(.|\n)*--count-kind all"):
            runner.build_dataset(self.make_sumstats_table(), old_result)

        all_runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig(use_common_counts=False))
        dataset = all_runner.build_dataset(self.make_sumstats_table(), old_result)
        self.assertEqual(dataset.count_key_used_for_regression, "all_reference_snp_counts")

    def test_build_dataset_chr_pos_mode_merges_on_coordinates_not_snp(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"), RegressionConfig())
        ldscore_result = replace(
            self.make_ldscore_result(),
            baseline_table=self.make_ldscore_result().baseline_table.assign(SNP=["ld1", "ld2", "ld3"]),
            query_table=self.make_ldscore_result().query_table.assign(SNP=["ld1", "ld2", "ld3"]),
            config_snapshot=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
        )
        sumstats = SumstatsTable(
            data=pd.DataFrame(
                {
                    "SNP": ["raw1", "raw2", "missing"],
                    "CHR": ["1", "1", pd.NA],
                    "POS": [10, 20, 30],
                    "Z": [2.0, 1.0, 0.5],
                    "N": [1000.0, 1000.0, 1000.0],
                    "A1": ["A", "C", "G"],
                    "A2": ["G", "T", "A"],
                }
            ),
            has_alleles=True,
            source_path="sumstats.gz",
            trait_name="trait",
            provenance={},
            config_snapshot=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
        )

        with self.assertLogs("LDSC.regression_runner", level="INFO") as caught:
            dataset = runner.build_dataset(sumstats, ldscore_result)

        self.assertEqual(dataset.merged["SNP"].tolist(), ["raw1", "raw2"])
        self.assertEqual(dataset.merged["base"].tolist(), [1.0, 2.0])
        self.assertIn("_ldsc_chr_pos_key", dataset.merged.columns)
        log_text = "\n".join(caught.output)
        self.assertIn("Dropped 1 SNPs with invalid or missing CHR/POS", log_text)
        self.assertIn("missing CHR=1", log_text)
        self.assertIn("sumstats.gz", log_text)

    def test_build_dataset_chr_pos_mode_ignores_query_snp_label_mismatch(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"), RegressionConfig())
        ldscore_result = replace(
            self.make_ldscore_result(),
            baseline_table=self.make_ldscore_result().baseline_table.assign(SNP=["base1", "base2", "base3"]),
            query_table=self.make_ldscore_result().query_table.assign(SNP=["query1", "query2", "query3"]),
            config_snapshot=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
        )
        sumstats = SumstatsTable(
            data=pd.DataFrame(
                {
                    "SNP": ["trait1", "trait2", "trait3"],
                    "CHR": ["1", "1", "1"],
                    "POS": [10, 20, 30],
                    "Z": [2.0, 1.0, 0.5],
                    "N": [1000.0, 1000.0, 1000.0],
                    "A1": ["A", "C", "G"],
                    "A2": ["G", "T", "A"],
                }
            ),
            has_alleles=True,
            source_path="sumstats.gz",
            trait_name="trait",
            provenance={},
            config_snapshot=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
        )

        dataset = runner.build_dataset(sumstats, ldscore_result, query_columns=["query1"])

        self.assertEqual(dataset.merged["SNP"].tolist(), ["trait1", "trait2", "trait3"])
        self.assertEqual(dataset.retained_ld_columns, ["base", "query1"])

    def test_build_dataset_unknown_ldscore_snapshot_defaults_query_alignment_to_chr_pos(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"), RegressionConfig())
        ldscore_result = replace(
            self.make_ldscore_result(),
            baseline_table=self.make_ldscore_result().baseline_table.assign(SNP=["base1", "base2", "base3"]),
            query_table=self.make_ldscore_result().query_table.assign(SNP=["query1", "query2", "query3"]),
            config_snapshot=None,
        )
        sumstats = SumstatsTable(
            data=pd.DataFrame(
                {
                    "SNP": ["trait1", "trait2", "trait3"],
                    "CHR": ["1", "1", "1"],
                    "POS": [10, 20, 30],
                    "Z": [2.0, 1.0, 0.5],
                    "N": [1000.0, 1000.0, 1000.0],
                    "A1": ["A", "C", "G"],
                    "A2": ["G", "T", "A"],
                }
            ),
            has_alleles=True,
            source_path="sumstats.gz",
            trait_name="trait",
            provenance={},
            config_snapshot=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
        )

        with mock.patch('ldsc._genome.build_policy.infer_chr_pos_build',
                        return_value=SimpleNamespace(genome_build='hg38', coordinate_basis='1-based', inspected_snp_count=500)):
            dataset = runner.build_dataset(sumstats, ldscore_result, query_columns=["query1"])

        self.assertEqual(dataset.merged["SNP"].tolist(), ["trait1", "trait2", "trait3"])
        self.assertEqual(dataset.retained_ld_columns, ["base", "query1"])

    def test_build_rg_dataset_uses_final_three_way_snp_intersection(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        sumstats_1 = self.make_sumstats_table()
        sumstats_2 = replace(
            self.make_sumstats_table(),
            trait_name="trait2",
            data=self.make_sumstats_table().data.assign(
                SNP=["rs2", "rs3", "missing"],
                Z=[1.5, 0.25, 9.0],
                A1=["C", "G", "A"],
                A2=["T", "A", "C"],
            ),
        )

        dataset = runner.build_rg_dataset(sumstats_1, sumstats_2, self.make_ldscore_result())

        self.assertEqual(dataset.merged["SNP"].tolist(), ["rs2", "rs3"])
        np.testing.assert_allclose(dataset.merged["Z1"], [1.0, 0.5])
        np.testing.assert_allclose(dataset.merged["Z2"], [1.5, 0.25])

    def test_build_rg_dataset_harmonizes_alleles_and_flips_second_trait_z(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        sumstats_1 = replace(
            self.make_sumstats_table(),
            data=self.make_sumstats_table().data.assign(
                A1=["A", "A", "A"],
                A2=["C", "C", "C"],
            ),
        )
        sumstats_2 = replace(
            self.make_sumstats_table(),
            trait_name="trait2",
            data=self.make_sumstats_table().data.assign(
                Z=[1.0, 2.0, 3.0],
                A1=["C", "A", "A"],
                A2=["A", "C", "T"],
            ),
        )

        dataset = runner.build_rg_dataset(sumstats_1, sumstats_2, self.make_ldscore_result())

        self.assertEqual(dataset.merged["SNP"].tolist(), ["rs1", "rs2"])
        np.testing.assert_allclose(dataset.merged["Z2"], [-1.0, 2.0])

    def test_build_rg_dataset_allele_aware_mode_harmonizes_after_identity_merge(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid_allele_aware"), RegressionConfig())
        ldscore_result = self.make_ldscore_result_from_frame(
            pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "SNP": ["rs1", "rs2"],
                    "POS": [10, 20],
                    "A1": ["A", "A"],
                    "A2": ["C", "G"],
                    "regression_ld_scores": [2.0, 2.0],
                    "base": [1.0, 2.0],
                }
            ),
            snp_identifier="rsid_allele_aware",
        )
        sumstats_1 = self.make_sumstats_table_from_frame(
            pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "SNP": ["rs1", "rs2"],
                    "POS": [10, 20],
                    "A1": ["A", "A"],
                    "A2": ["C", "G"],
                    "Z": [2.0, 1.0],
                    "N": [1000.0, 1000.0],
                }
            ),
            snp_identifier="rsid_allele_aware",
            trait_name="trait1",
        )
        sumstats_2 = self.make_sumstats_table_from_frame(
            pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "SNP": ["rs1", "rs2"],
                    "POS": [10, 20],
                    "A1": ["C", "A"],
                    "A2": ["A", "C"],
                    "Z": [5.0, 9.0],
                    "N": [900.0, 900.0],
                }
            ),
            snp_identifier="rsid_allele_aware",
            trait_name="trait2",
        )

        dataset = runner.build_rg_dataset(sumstats_1, sumstats_2, ldscore_result)

        self.assertEqual(dataset.merged["SNP"].tolist(), ["rs1"])
        np.testing.assert_allclose(dataset.merged["Z2"], [-5.0])
        self.assertEqual(dataset.effective_snp_identifier, "rsid_allele_aware")
        self.assertFalse(dataset.identity_downgrade_applied)

    def test_build_rg_dataset_rejects_trait_without_alleles(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())

        with self.assertRaisesRegex(LDSCInputError, r"A1/A2(.|\n)*trait2 \('trait2.parquet'\)"):
            runner.build_rg_dataset(self.make_sumstats_table(), self.allele_free_table("trait2"),
                                    self.make_ldscore_result())

    def test_build_dataset_warns_when_dropping_zero_variance_ld_columns(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        base_result = self.make_ldscore_result()
        n_rows = len(base_result.baseline_table)
        ldscore_result = replace(
            base_result,
            baseline_table=base_result.baseline_table.assign(flat=[1.0] * n_rows),
            baseline_columns=["base", "flat"],
            count_records=[
                {
                    "group": "baseline",
                    "column": "base",
                    "all_reference_snp_count": 10.0,
                    "common_reference_snp_count": 8.0,
                },
                {
                    "group": "baseline",
                    "column": "flat",
                    "all_reference_snp_count": 30.0,
                    "common_reference_snp_count": 28.0,
                },
            ],
        )
        with self.assertLogs("LDSC.regression_runner", level="WARNING") as captured:
            dataset = runner.build_dataset(self.make_sumstats_table(), ldscore_result)

        self.assertEqual(dataset.dropped_zero_variance_ld_columns, ["flat"])
        self.assertEqual(dataset.retained_ld_columns, ["base"])
        self.assertTrue(any("Dropped 1 zero-variance" in message and "fitted-model metadata" in message for message in captured.output))
        self.assertFalse(any("flat" in message for message in captured.output))

    def test_build_rg_dataset_drops_zero_variance_ld_columns_after_final_merge(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        ldscore_result = replace(
            self.make_ldscore_result(),
            baseline_table=self.make_ldscore_result().baseline_table.assign(late_zero=[0.0, 1.0, 1.0]),
            baseline_columns=["base", "late_zero"],
            count_records=[
                {
                    "group": "baseline",
                    "column": "base",
                    "all_reference_snp_count": 10.0,
                    "common_reference_snp_count": 8.0,
                },
                {
                    "group": "baseline",
                    "column": "late_zero",
                    "all_reference_snp_count": 30.0,
                    "common_reference_snp_count": 28.0,
                },
            ],
        )
        sumstats_2 = replace(
            self.make_sumstats_table(),
            trait_name="trait2",
            data=self.make_sumstats_table().data.iloc[[1, 2]].reset_index(drop=True).assign(Z=[1.5, 0.25]),
        )

        dataset = runner.build_rg_dataset(self.make_sumstats_table(), sumstats_2, ldscore_result)

        self.assertEqual(dataset.merged["SNP"].tolist(), ["rs2", "rs3"])
        self.assertEqual(dataset.retained_ld_columns, ["base"])
        self.assertEqual(dataset.dropped_zero_variance_ld_columns, ["late_zero"])
        np.testing.assert_allclose(
            dataset.reference_snp_count_totals["common_reference_snp_counts"],
            [8.0],
        )

    def test_build_rg_dataset_chr_pos_mode_merges_all_sources_on_coordinates(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"), RegressionConfig())
        ldscore_result = replace(
            self.make_ldscore_result(),
            baseline_table=self.make_ldscore_result().baseline_table.assign(SNP=["ld1", "ld2", "ld3"]),
            query_table=self.make_ldscore_result().query_table.assign(SNP=["ld1", "ld2", "ld3"]),
            config_snapshot=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
        )
        sumstats_1 = SumstatsTable(
            data=pd.DataFrame(
                {
                    "SNP": ["trait1_a", "trait1_b", "trait1_missing"],
                    "CHR": ["1", "1", "1"],
                    "POS": [10, 20, 999],
                    "Z": [2.0, 1.0, 9.0],
                    "N": [1000.0, 1000.0, 1000.0],
                    "A1": ["A", "A", "A"],
                    "A2": ["C", "C", "C"],
                }
            ),
            has_alleles=True,
            source_path="sumstats1.gz",
            trait_name="trait1",
            provenance={},
            config_snapshot=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
        )
        sumstats_2 = replace(
            sumstats_1,
            data=pd.DataFrame(
                {
                    "SNP": ["trait2_a", "trait2_c"],
                    "CHR": ["1", "1"],
                    "POS": [10, 30],
                    "Z": [1.0, 1.5],
                    "N": [900.0, 900.0],
                    "A1": ["A", "A"],
                    "A2": ["C", "C"],
                }
            ),
            trait_name="trait2",
            source_path="sumstats2.gz",
        )

        dataset = runner.build_rg_dataset(sumstats_1, sumstats_2, ldscore_result)

        self.assertEqual(dataset.merged["SNP"].tolist(), ["trait1_a"])
        self.assertEqual(dataset.merged["_ldsc_chr_pos_key"].tolist(), ["1:10"])
        np.testing.assert_allclose(dataset.merged["Z2"], [1.0])

    def test_build_dataset_raises_on_mismatched_sumstats_and_ldscore_snapshots(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"), RegressionConfig())
        sumstats = replace(
            self.make_sumstats_table(),
            config_snapshot=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
        )
        ldscore_result = replace(
            self.make_ldscore_result(),
            config_snapshot=GlobalConfig(genome_build="hg19", snp_identifier="chr_pos"),
        )

        with self.assertRaisesRegex(LDSCInputError, "Conflicting genome-build evidence"):
            runner.build_dataset(sumstats, ldscore_result)

    def _chr_pos_ldscore(self, genome_build="hg38"):
        return replace(
            self.make_ldscore_result(),
            config_snapshot=GlobalConfig(snp_identifier="chr_pos", genome_build=genome_build),
        )

    def test_regression_genome_build_skips_rsid_family(self):
        sumstats = replace(self.make_sumstats_table(), config_snapshot=None)
        ldscore = replace(self.make_ldscore_result(), config_snapshot=GlobalConfig(snp_identifier="rsid_allele_aware"))
        # rsID-family run: no error even though the sumstats carries no build metadata.
        alignment._regression_genome_build(sumstats, ldscore, "rsid_allele_aware", context="sumstats")

    def test_regression_genome_build_passes_on_metadata_match(self):
        sumstats = replace(
            self.make_sumstats_table(),
            config_snapshot=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
        )
        alignment._regression_genome_build(sumstats, self._chr_pos_ldscore("hg38"), "chr_pos", context="sumstats")

    def test_regression_genome_build_errors_on_metadata_mismatch(self):
        sumstats = replace(
            self.make_sumstats_table(),
            config_snapshot=GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"),
        )
        with self.assertRaisesRegex(LDSCInputError, "Conflicting genome-build evidence"):
            alignment._regression_genome_build(sumstats, self._chr_pos_ldscore("hg38"), "chr_pos", context="sumstats")

    def test_regression_genome_build_infers_when_metadata_absent(self):
        frame = pd.DataFrame({"SNP": ["rs1"], "CHR": ["1"], "POS": [10], "Z": [1.0], "N": [100.0]})
        sumstats = replace(self.make_sumstats_table(), data=frame, config_snapshot=None)
        with mock.patch(
            "ldsc._genome.build_policy.infer_chr_pos_build", return_value=SimpleNamespace(genome_build="hg38", coordinate_basis="1-based", inspected_snp_count=500)
        ):
            alignment._regression_genome_build(sumstats, self._chr_pos_ldscore("hg38"), "chr_pos", context="sumstats")
        with mock.patch(
            "ldsc._genome.build_policy.infer_chr_pos_build", return_value=SimpleNamespace(genome_build="hg19", coordinate_basis="1-based", inspected_snp_count=500)
        ):
            with self.assertRaisesRegex(LDSCInputError, "Conflicting genome-build evidence"):
                alignment._regression_genome_build(sumstats, self._chr_pos_ldscore("hg38"), "chr_pos", context="sumstats")

    def test_regression_genome_build_errors_when_inference_inconclusive(self):
        frame = pd.DataFrame({"SNP": ["rs1"], "CHR": ["1"], "POS": [10], "Z": [1.0], "N": [100.0]})
        sumstats = replace(self.make_sumstats_table(), data=frame, config_snapshot=None)
        with mock.patch(
            "ldsc._genome.build_policy.infer_chr_pos_build", side_effect=LDSCInputError("could not infer")
        ):
            with self.assertRaisesRegex(LDSCInputError, "unresolved:"):
                alignment._regression_genome_build(sumstats, self._chr_pos_ldscore("hg38"), "chr_pos", context="sumstats")

    def test_regression_genome_build_errors_when_chr_pos_columns_missing(self):
        # Metadata-absent chr_pos run on a sumstats lacking CHR/POS (e.g. a legacy rsID file).
        sumstats = replace(self.make_sumstats_table(), config_snapshot=None)
        with self.assertRaisesRegex(LDSCInputError, "unresolved:"):
            alignment._regression_genome_build(sumstats, self._chr_pos_ldscore("hg38"), "chr_pos", context="sumstats")

    def test_build_dataset_skips_compatibility_check_for_legacy_none_snapshots(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        sumstats = replace(self.make_sumstats_table(), config_snapshot=None)
        ldscore_result = replace(self.make_ldscore_result(), config_snapshot=None)

        dataset = runner.build_dataset(sumstats, ldscore_result)

        self.assertIsNone(dataset.config_snapshot)

    def test_build_dataset_allows_unknown_sumstats_with_known_ldscore_snapshot(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        sumstats = replace(self.make_sumstats_table(), config_snapshot=None)

        dataset = runner.build_dataset(sumstats, self.make_ldscore_result())

        self.assertEqual(dataset.config_snapshot, GlobalConfig(snp_identifier="rsid"))
