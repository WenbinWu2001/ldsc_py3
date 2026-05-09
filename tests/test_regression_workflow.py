import argparse
from dataclasses import replace
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
import warnings
from unittest import mock

import numpy as np
import pandas as pd
from scipy import stats

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc.config import ConfigMismatchError, GlobalConfig, RegressionConfig, reset_global_config, set_global_config

try:
    from ldsc.ldscore_calculator import LDScoreResult
    from ldsc import regression_runner
    from ldsc.regression_runner import RegressionRunner
    from ldsc.sumstats_munger import SumstatsTable
    from ldsc.annotation_builder import AnnotationBundle
except ImportError:
    LDScoreResult = None
    RegressionRunner = None
    SumstatsTable = None
    AnnotationBundle = None


@unittest.skipIf(RegressionRunner is None, "regression_runner module is not available")
class RegressionWorkflowTest(unittest.TestCase):
    def setUp(self):
        reset_global_config()

    def tearDown(self):
        reset_global_config()

    def test_load_ldscore_from_dir_is_public(self):
        from ldsc import load_ldscore_from_dir

        self.assertTrue(callable(load_ldscore_from_dir))

    def test_load_ldscore_from_dir_reads_manifest_and_parquet_files(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "POS": [10],
                    "regr_weight": [3.0],
                    "base": [1.0],
                }
            ).to_parquet(tmpdir / "baseline.parquet", index=False)
            pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "POS": [10],
                    "query": [2.0],
                }
            ).to_parquet(tmpdir / "query.parquet", index=False)
            (tmpdir / "manifest.json").write_text(
                json.dumps(
                    {
                        "format": "ldsc.ldscore_result.v1",
                        "files": {"baseline": "baseline.parquet", "query": "query.parquet"},
                        "snp_identifier": "rsid",
                        "genome_build": "hg38",
                        "chromosomes": ["1"],
                        "baseline_columns": ["base"],
                        "query_columns": ["query"],
                        "counts": [
                            {
                                "group": "baseline",
                                "column": "base",
                                "all_reference_snp_count": 5.0,
                                "common_reference_snp_count": 4.0,
                            },
                            {
                                "group": "query",
                                "column": "query",
                                "all_reference_snp_count": 6.0,
                                "common_reference_snp_count": 5.0,
                            },
                        ],
                        "config_snapshot": {"snp_identifier": "rsid", "genome_build": "hg38", "log_level": "INFO"},
                    }
                ),
                encoding="utf-8",
            )
            result = load_ldscore_from_dir(str(tmpdir))

        self.assertEqual(result.baseline_table.columns.tolist(), ["CHR", "SNP", "POS", "regr_weight", "base"])
        self.assertEqual(result.query_table.columns.tolist(), ["CHR", "SNP", "POS", "query"])
        self.assertEqual(result.count_records[0]["column"], "base")
        self.assertEqual(result.ld_regression_snps, frozenset({"rs1"}))
        self.assertEqual(result.config_snapshot, GlobalConfig(snp_identifier="rsid"))

    def test_load_ldscore_from_dir_rejects_legacy_bp_coordinate_schema(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "BP": [10],
                    "regr_weight": [3.0],
                    "base": [1.0],
                }
            ).to_parquet(tmpdir / "baseline.parquet", index=False)
            (tmpdir / "manifest.json").write_text(
                json.dumps(
                    {
                        "format": "ldsc.ldscore_result.v1",
                        "files": {"baseline": "baseline.parquet"},
                        "snp_identifier": "rsid",
                        "genome_build": "hg38",
                        "chromosomes": ["1"],
                        "baseline_columns": ["base"],
                        "query_columns": [],
                        "counts": [
                            {
                                "group": "baseline",
                                "column": "base",
                                "all_reference_snp_count": 5.0,
                            }
                        ],
                        "config_snapshot": {"snp_identifier": "rsid", "genome_build": "hg38", "log_level": "INFO"},
                    }
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "baseline_table is missing required columns.*POS"):
                load_ldscore_from_dir(str(tmpdir))

    def test_load_ldscore_from_dir_warns_when_config_snapshot_is_missing(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self.write_ldscore_dir(tmpdir, include_query=False)
            manifest_path = tmpdir / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            del manifest["config_snapshot"]
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                result = load_ldscore_from_dir(str(tmpdir))

        self.assertIsNone(result.config_snapshot)
        self.assertTrue(any("GlobalConfig provenance is missing" in str(item.message) for item in caught))

    def test_load_ldscore_from_dir_warns_when_config_snapshot_is_invalid(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self.write_ldscore_dir(tmpdir, include_query=False)
            manifest_path = tmpdir / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["config_snapshot"] = {"snp_identifier": "bad", "genome_build": "hg38", "log_level": "INFO"}
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                result = load_ldscore_from_dir(str(tmpdir))

        self.assertIsNone(result.config_snapshot)
        self.assertTrue(any("GlobalConfig provenance is invalid" in str(item.message) for item in caught))

    def make_ldscore_result(self):
        baseline_table = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "SNP": ["rs1", "rs2", "rs3"],
                "POS": [10, 20, 30],
                "regr_weight": [2.0, 2.0, 2.0],
                "base": [1.0, 2.0, 3.0],
            }
        )
        query_table = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "SNP": ["rs1", "rs2", "rs3"],
                "POS": [10, 20, 30],
                "query1": [1.0, 2.0, 3.0],
                "query2": [0.5, 1.5, 2.5],
            }
        )
        return LDScoreResult(
            baseline_table=baseline_table,
            query_table=query_table,
            count_records=[
                {
                    "group": "baseline",
                    "column": "base",
                    "all_reference_snp_count": 10.0,
                    "common_reference_snp_count": 8.0,
                },
                {
                    "group": "query",
                    "column": "query1",
                    "all_reference_snp_count": 20.0,
                    "common_reference_snp_count": 18.0,
                },
                {
                    "group": "query",
                    "column": "query2",
                    "all_reference_snp_count": 30.0,
                    "common_reference_snp_count": 28.0,
                },
            ],
            baseline_columns=["base"],
            query_columns=["query1", "query2"],
            ld_reference_snps=frozenset(),
            ld_regression_snps=frozenset({"rs1", "rs2", "rs3"}),
            chromosome_results=[],
            config_snapshot=GlobalConfig(snp_identifier="rsid"),
        )

    def make_sumstats_table(self):
        return SumstatsTable(
            data=pd.DataFrame(
                {
                    "SNP": ["rs1", "rs2", "rs3"],
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
            config_snapshot=GlobalConfig(snp_identifier="rsid"),
        )

    def make_annotation_bundle(self):
        metadata = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "POS": [10, 20, 30],
                "SNP": ["rs1", "rs2", "rs3"],
                "CM": [0.1, 0.2, 0.3],
            }
        )
        return AnnotationBundle(
            metadata=metadata,
            baseline_annotations=pd.DataFrame({"base": [1, 1, 1]}),
            query_annotations=pd.DataFrame({"query1": [1, 0, 0], "query2": [0, 1, 0]}),
            baseline_columns=["base"],
            query_columns=["query1", "query2"],
            chromosomes=["1"],
            source_summary={},
        )

    def write_ldscore_dir(self, root: Path, *, include_query: bool = True) -> Path:
        root.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            {
                "CHR": ["1"],
                "SNP": ["rs1"],
                "POS": [10],
                "regr_weight": [2.0],
                "base": [1.0],
            }
        ).to_parquet(root / "baseline.parquet", index=False)
        files = {"baseline": "baseline.parquet"}
        query_columns = []
        counts = [
            {
                "group": "baseline",
                "column": "base",
                "all_reference_snp_count": 5.0,
                "common_reference_snp_count": 5.0,
            }
        ]
        if include_query:
            pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "POS": [10],
                    "query": [2.0],
                }
            ).to_parquet(root / "query.parquet", index=False)
            files["query"] = "query.parquet"
            query_columns = ["query"]
            counts.append(
                {
                    "group": "query",
                    "column": "query",
                    "all_reference_snp_count": 6.0,
                    "common_reference_snp_count": 6.0,
                }
            )
        (root / "manifest.json").write_text(
            json.dumps(
                {
                    "format": "ldsc.ldscore_result.v1",
                    "files": files,
                    "snp_identifier": "rsid",
                    "genome_build": "hg38",
                    "chromosomes": ["1"],
                    "baseline_columns": ["base"],
                    "query_columns": query_columns,
                    "counts": counts,
                    "config_snapshot": {"snp_identifier": "rsid", "genome_build": "hg38", "log_level": "INFO"},
                }
            ),
            encoding="utf-8",
        )
        return root

    def test_build_dataset_uses_baseline_only_for_h2_style_runs(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        dataset = runner.build_dataset(self.make_sumstats_table(), self.make_ldscore_result())
        self.assertEqual(dataset.count_key_used_for_regression, "common_reference_snp_counts")
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

    def test_build_dataset_can_include_one_query_annotation_for_partitioned_h2(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        dataset = runner.build_dataset(self.make_sumstats_table(), self.make_ldscore_result(), query_columns=["query2"])
        self.assertEqual(dataset.retained_ld_columns, ["base", "query2"])
        np.testing.assert_allclose(
            dataset.reference_snp_count_totals["common_reference_snp_counts"],
            [8.0, 28.0],
        )

    def test_build_dataset_empty_merge_error_includes_active_config(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        sumstats = replace(
            self.make_sumstats_table(),
            data=self.make_sumstats_table().data.assign(SNP=["missing1", "missing2", "missing3"]),
        )

        with self.assertRaisesRegex(ValueError, "Active config:"):
            runner.build_dataset(sumstats, self.make_ldscore_result())

    def test_build_dataset_with_query_rejects_misaligned_query_table(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        result = self.make_ldscore_result()
        bad_query = result.query_table.copy()
        bad_query.loc[1, "POS"] = 999
        result = replace(result, query_table=bad_query)

        with self.assertRaisesRegex(ValueError, "query rows must match baseline rows"):
            runner.build_dataset(self.make_sumstats_table(), result, query_columns=["query2"])

    def test_old_common_count_manifest_key_is_not_recognized(self):
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

        dataset = runner.build_dataset(self.make_sumstats_table(), old_result)

        self.assertEqual(dataset.count_key_used_for_regression, "all_reference_snp_counts")
        self.assertNotIn("common_reference_snp_counts", dataset.reference_snp_count_totals)

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

        dataset = runner.build_dataset(sumstats, ldscore_result)

        self.assertEqual(dataset.merged["SNP"].tolist(), ["raw1", "raw2"])
        self.assertEqual(dataset.merged["base"].tolist(), [1.0, 2.0])
        self.assertIn("_ldsc_chr_pos_key", dataset.merged.columns)

    def test_estimate_rg_uses_baseline_only_when_query_exists(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        sumstats_1 = self.make_sumstats_table()
        sumstats_2 = replace(
            self.make_sumstats_table(),
            trait_name="trait2",
            data=self.make_sumstats_table().data.assign(Z=[1.0, 1.5, 0.25]),
        )
        with mock.patch.object(regression_runner.reg, "RG", return_value=mock.sentinel.rg_result) as patched:
            result = runner.estimate_rg(sumstats_1, sumstats_2, self.make_ldscore_result())

        self.assertIs(result, mock.sentinel.rg_result)
        self.assertEqual(patched.call_args.args[2].shape[1], 1)

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

    def test_build_rg_dataset_skips_allele_harmonization_when_alleles_are_absent(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        sumstats_1 = replace(
            self.make_sumstats_table(),
            has_alleles=False,
            data=self.make_sumstats_table().data.drop(columns=["A1", "A2"]),
        )
        sumstats_2 = replace(
            self.make_sumstats_table(),
            has_alleles=False,
            trait_name="trait2",
            data=self.make_sumstats_table().data.drop(columns=["A1", "A2"]).assign(Z=[1.0, 2.0, 3.0]),
        )

        dataset = runner.build_rg_dataset(sumstats_1, sumstats_2, self.make_ldscore_result())

        self.assertEqual(dataset.merged["SNP"].tolist(), ["rs1", "rs2", "rs3"])
        self.assertFalse({"A1", "A2", "A1x", "A2x"} & set(dataset.merged.columns))
        np.testing.assert_allclose(dataset.merged["Z2"], [1.0, 2.0, 3.0])

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

    def test_estimate_rg_chr_pos_mode_merges_traits_on_coordinates(self):
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
                    "SNP": ["trait1_a", "trait1_b"],
                    "CHR": ["1", "1"],
                    "POS": [10, 20],
                    "Z": [2.0, 1.0],
                    "N": [1000.0, 1000.0],
                    "A1": ["A", "C"],
                    "A2": ["G", "T"],
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
                    "SNP": ["trait2_a", "trait2_b"],
                    "CHR": ["1", "1"],
                    "POS": [10, 20],
                    "Z": [1.0, 1.5],
                    "N": [900.0, 900.0],
                    "A1": ["A", "C"],
                    "A2": ["G", "T"],
                }
            ),
            trait_name="trait2",
            source_path="sumstats2.gz",
        )
        with mock.patch.object(regression_runner.reg, "RG", return_value=mock.sentinel.rg_result) as patched:
            result = runner.estimate_rg(sumstats_1, sumstats_2, ldscore_result)

        self.assertIs(result, mock.sentinel.rg_result)
        self.assertEqual(patched.call_args.args[0].shape[0], 2)
        self.assertEqual(patched.call_args.args[1].shape[0], 2)

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

    def test_estimate_rg_preserves_legacy_rg_array_call_contract(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        sumstats_1 = self.make_sumstats_table()
        sumstats_2 = replace(
            self.make_sumstats_table(),
            trait_name="trait2",
            data=self.make_sumstats_table().data.assign(Z=[1.0, 1.5, 0.25], N=[900.0, 900.0, 900.0]),
        )

        with mock.patch.object(regression_runner.reg, "RG", return_value=mock.sentinel.rg_result) as patched:
            result = runner.estimate_rg(sumstats_1, sumstats_2, self.make_ldscore_result())

        self.assertIs(result, mock.sentinel.rg_result)
        np.testing.assert_allclose(patched.call_args.args[0], [[2.0], [1.0], [0.5]])
        np.testing.assert_allclose(patched.call_args.args[1], [[1.0], [1.5], [0.25]])
        np.testing.assert_allclose(patched.call_args.args[2], [[1.0], [2.0], [3.0]])
        np.testing.assert_allclose(patched.call_args.args[3], [[2.0], [2.0], [2.0]])
        np.testing.assert_allclose(patched.call_args.args[4], [[1000.0], [1000.0], [1000.0]])
        np.testing.assert_allclose(patched.call_args.args[5], [[900.0], [900.0], [900.0]])
        np.testing.assert_allclose(patched.call_args.args[6], [[8.0]])

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

        with self.assertRaisesRegex(ConfigMismatchError, "genome_build mismatch"):
            runner.build_dataset(sumstats, ldscore_result)

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

    def test_estimate_partitioned_h2_requires_query_annotations(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        ldscore_result = replace(
            self.make_ldscore_result(),
            query_table=None,
            query_columns=[],
            count_records=[
                {
                    "group": "baseline",
                    "column": "base",
                    "all_reference_snp_count": 10.0,
                    "common_reference_snp_count": 8.0,
                }
            ],
        )

        with self.assertRaisesRegex(ValueError, "partitioned-h2 requires query annotations"):
            runner.estimate_partitioned_h2(
                self.make_sumstats_table(),
                ldscore_result,
                query_column="base",
            )

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
            coef=np.array([0.0, 1.0]),
            coef_cov=np.diag([0.01, 0.04]),
            coef_se=np.array([0.1, 0.2]),
            cat=np.array([0.0, 0.3]),
            cat_se=np.array([0.01, 0.03]),
            prop=np.array([0.0, 0.4]),
            prop_se=np.array([0.01, 0.04]),
            enrichment=np.array([0.0, 2.0]),
            n_blocks=200,
        )

        with mock.patch.object(
            runner,
            "build_dataset",
            wraps=runner.build_dataset,
        ) as build_dataset, mock.patch.object(runner, "estimate_h2", return_value=fake_hsq):
            result = runner.estimate_partitioned_h2(
                self.make_sumstats_table(),
                self.make_ldscore_result(),
                query_column="query1",
            )

        self.assertEqual(build_dataset.call_args.kwargs["query_columns"], ["query1"])
        self.assertEqual(result["Category"].tolist(), ["query1"])

    def test_estimate_partitioned_h2_rejects_unknown_query_column(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())

        with self.assertRaisesRegex(ValueError, "Unknown query annotation requested.*missing"):
            runner.estimate_partitioned_h2(
                self.make_sumstats_table(),
                self.make_ldscore_result(),
                query_column="missing",
            )

    def test_estimate_partitioned_h2_batch_rejects_empty_query_columns(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        annotation_bundle = replace(
            self.make_annotation_bundle(),
            query_annotations=pd.DataFrame(index=pd.RangeIndex(3)),
            query_columns=[],
        )

        with self.assertRaisesRegex(ValueError, "partitioned-h2 requires query annotations"):
            runner.estimate_partitioned_h2_batch(
                self.make_sumstats_table(),
                self.make_ldscore_result(),
                annotation_bundle,
            )

    def test_estimate_partitioned_h2_batch_loops_over_queries(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        table = self.make_sumstats_table()
        ldscore_result = self.make_ldscore_result()
        annotation_bundle = self.make_annotation_bundle()
        fake_hsq = mock.Mock(
            coef=np.array([0.0, 1.0]),
            coef_cov=np.diag([0.01, 0.04]),
            coef_se=np.array([0.1, 0.2]),
            cat=np.array([0.0, 0.3]),
            cat_se=np.array([0.01, 0.03]),
            prop=np.array([0.0, 0.4]),
            prop_se=np.array([0.01, 0.04]),
            enrichment=np.array([0.0, 2.0]),
            n_blocks=200,
        )
        with mock.patch.object(runner, "estimate_h2", return_value=fake_hsq) as patched:
            result = runner.estimate_partitioned_h2_batch(table, ldscore_result, annotation_bundle)
        self.assertEqual(patched.call_count, 2)
        self.assertEqual(
            result.columns.tolist(),
            [
                "Category",
                "Prop._SNPs",
                "Prop._h2",
                "Enrichment",
                "Enrichment_p",
                "Coefficient",
                "Coefficient_p",
            ],
        )
        self.assertEqual(result["Category"].tolist(), ["query1", "query2"])
        self.assertAlmostEqual(result.loc[0, "Prop._SNPs"], 18.0 / 26.0)
        self.assertAlmostEqual(result.loc[0, "Coefficient_p"], 2 * stats.norm.sf(5.0))

    def test_estimate_partitioned_h2_batch_can_return_full_partitioned_tables(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        table = self.make_sumstats_table()
        ldscore_result = self.make_ldscore_result()
        annotation_bundle = self.make_annotation_bundle()
        fake_hsq = mock.Mock(
            coef=np.array([0.0, 1.0]),
            coef_cov=np.diag([0.01, 0.04]),
            coef_se=np.array([0.1, 0.2]),
            cat=np.array([0.0, 0.3]),
            cat_se=np.array([0.01, 0.03]),
            prop=np.array([0.0, 0.4]),
            prop_se=np.array([0.01, 0.04]),
            enrichment=np.array([0.0, 2.0]),
            n_blocks=200,
        )
        with mock.patch.object(runner, "estimate_h2", return_value=fake_hsq):
            result = runner.estimate_partitioned_h2_batch(
                table,
                ldscore_result,
                annotation_bundle,
                include_full_partitioned_h2=True,
            )

        self.assertIsInstance(result, regression_runner.PartitionedH2BatchResult)
        self.assertEqual(result.summary["Category"].tolist(), ["query1", "query2"])
        self.assertIn("query1", result.per_query_category_tables)
        self.assertEqual(
            result.per_query_category_tables["query1"].columns.tolist(),
            [
                "Category",
                "Prop._SNPs",
                "Category_h2",
                "Category_h2_std_error",
                "Prop._h2",
                "Prop._h2_std_error",
                "Enrichment",
                "Enrichment_std_error",
                "Enrichment_p",
                "Coefficient",
                "Coefficient_std_error",
                "Coefficient_p",
            ],
        )

    def test_h2_and_partitioned_summaries_share_one_fitted_hsq_result(self):
        dataset = regression_runner.RegressionDataset(
            merged=pd.DataFrame(
                {
                    "SNP": ["rs1", "rs2"],
                    "Z": [1.0, 2.0],
                    "N": [1000.0, 1000.0],
                    "regr_weight": [1.0, 1.0],
                }
            ),
            ref_ld_columns=["base"],
            weight_column="regr_weight",
            reference_snp_count_totals={"common_reference_snp_counts": np.array([10.0])},
            count_key_used_for_regression="common_reference_snp_counts",
            retained_ld_columns=["base"],
            dropped_zero_variance_ld_columns=[],
            trait_names=["trait"],
            chromosomes_aggregated=[],
        )
        hsq = mock.Mock(
            tot=np.array([0.25]),
            tot_se=np.array([0.03]),
            intercept=np.array([1.01]),
            intercept_se=0.02,
            mean_chisq=np.array([1.2]),
            lambda_gc=np.array([1.1]),
            ratio=0.05,
            ratio_se=0.01,
            coef=np.array([0.025]),
            coef_cov=np.array([[0.000009]]),
            coef_se=np.array([0.003]),
            cat=np.array([0.25]),
            cat_se=np.array([0.03]),
            prop=np.array([1.0]),
            prop_se=np.array([0.0]),
            enrichment=np.array([1.0]),
            n_blocks=200,
        )

        total = regression_runner.summarize_total_h2(hsq, dataset, trait_name="trait")
        partitioned = regression_runner.summarize_partitioned_h2(
            hsq,
            dataset,
            ["base"],
            include_full_columns=True,
        )

        self.assertEqual(total.loc[0, "total_h2"], partitioned.loc[0, "Category_h2"])
        self.assertEqual(total.loc[0, "total_h2_se"], partitioned.loc[0, "Category_h2_std_error"])
        self.assertEqual(partitioned.loc[0, "Category"], "base")
        self.assertEqual(partitioned.loc[0, "Prop._h2"], 1.0)
        self.assertEqual(partitioned.loc[0, "Enrichment"], 1.0)
        self.assertTrue(np.isnan(partitioned.loc[0, "Enrichment_p"]))

    def test_partitioned_full_summary_uses_legacy_column_order_and_enrichment_p(self):
        dataset = regression_runner.RegressionDataset(
            merged=pd.DataFrame(
                {
                    "SNP": ["rs1", "rs2", "rs3"],
                    "Z": [1.0, 1.5, 2.0],
                    "N": [1000.0, 1000.0, 1000.0],
                    "regr_weight": [1.0, 1.0, 1.0],
                }
            ),
            ref_ld_columns=["base", "query"],
            weight_column="regr_weight",
            reference_snp_count_totals={"common_reference_snp_counts": np.array([10.0, 30.0])},
            count_key_used_for_regression="common_reference_snp_counts",
            retained_ld_columns=["base", "query"],
            dropped_zero_variance_ld_columns=[],
            trait_names=["trait"],
            chromosomes_aggregated=[],
        )
        hsq = mock.Mock(
            coef=np.array([0.2, 0.5]),
            coef_cov=np.diag([0.01, 0.04]),
            coef_se=np.array([0.1, 0.2]),
            cat=np.array([0.2, 0.6]),
            cat_se=np.array([0.02, 0.06]),
            prop=np.array([0.25, 0.75]),
            prop_se=np.array([0.03, 0.09]),
            enrichment=np.array([1.0, 1.0]),
            n_blocks=200,
        )

        summary = regression_runner.summarize_partitioned_h2(
            hsq,
            dataset,
            ["query"],
            include_full_columns=True,
        )

        self.assertEqual(
            summary.columns.tolist(),
            [
                "Category",
                "Prop._SNPs",
                "Category_h2",
                "Category_h2_std_error",
                "Prop._h2",
                "Prop._h2_std_error",
                "Enrichment",
                "Enrichment_std_error",
                "Enrichment_p",
                "Coefficient",
                "Coefficient_std_error",
                "Coefficient_p",
            ],
        )
        self.assertEqual(summary.loc[0, "Category"], "query")
        self.assertAlmostEqual(summary.loc[0, "Prop._SNPs"], 0.75)
        self.assertAlmostEqual(summary.loc[0, "Enrichment_std_error"], 0.09 / 0.75)
        self.assertAlmostEqual(summary.loc[0, "Coefficient_p"], 2 * stats.norm.sf(2.5))
        expected_enrichment_p = 2 * stats.t.sf(abs(0.3 / np.sqrt(0.05)), 200)
        self.assertAlmostEqual(summary.loc[0, "Enrichment_p"], expected_enrichment_p)

    def test_run_h2_from_args_uses_ldscore_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\nrs1\t1.0\t1000\n")
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            args = type(
                "Args",
                (),
                {
                    "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                    "trait_name": "trait",
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": str(tmpdir / "h2_out"),
                    "overwrite": False,
                    "log_level": "INFO",
                    "n_blocks": 200,
                    "no_intercept": False,
                    "intercept_h2": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                    "log_level": "INFO",
                },
            )()

            with mock.patch.object(
                regression_runner.RegressionRunner,
                "estimate_h2",
                return_value=mock.Mock(
                    tot=np.array([0.1]),
                    tot_se=np.array([0.01]),
                    intercept=np.array([1.0]),
                    intercept_se=0.01,
                    mean_chisq=np.array([1.1]),
                    lambda_gc=np.array([1.0]),
                    ratio=0.0,
                    ratio_se=0.0,
                ),
            ) as patched:
                summary = regression_runner.run_h2_from_args(args)

            patched.assert_called_once()
            self.assertEqual(patched.call_args.args[0].retained_ld_columns, ["base"])
            self.assertEqual(summary.loc[0, "trait_name"], "trait")
            self.assertTrue((tmpdir / "h2_out" / "h2.tsv").exists())
            self.assertTrue((tmpdir / "h2_out" / "h2.log").exists())

    def test_run_h2_from_args_without_output_dir_creates_no_log_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\nrs1\t1.0\t1000\n")
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            args = type(
                "Args",
                (),
                {
                    "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                    "trait_name": "trait",
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": None,
                    "overwrite": False,
                    "log_level": "INFO",
                    "n_blocks": 200,
                    "no_intercept": False,
                    "intercept_h2": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                },
            )()

            with mock.patch.object(
                regression_runner.RegressionRunner,
                "estimate_h2",
                return_value=mock.Mock(
                    tot=np.array([0.1]),
                    tot_se=np.array([0.01]),
                    intercept=np.array([1.0]),
                    intercept_se=0.01,
                    mean_chisq=np.array([1.1]),
                    lambda_gc=np.array([1.0]),
                    ratio=0.0,
                    ratio_se=0.0,
                ),
            ):
                regression_runner.run_h2_from_args(args)

            self.assertFalse(list(tmpdir.glob("*.log")))

    def test_common_regression_arguments_expose_only_ldscore_dir(self):
        parser = argparse.ArgumentParser()
        regression_runner.add_h2_arguments(parser)

        args = parser.parse_args(["--ldscore-dir", "ldscores", "--sumstats-file", "trait.sumstats.gz"])

        self.assertEqual(args.ldscore_dir, "ldscores")
        self.assertEqual(args.sumstats_file, "trait.sumstats.gz")
        self.assertEqual(args.count_kind, "common")
        self.assertEqual(args.log_level, "INFO")
        parsed = parser.parse_args(
            [
                "--ldscore-dir",
                "ldscores",
                "--sumstats-file",
                "trait.sumstats.gz",
                "--count-kind",
                "all",
                "--log-level",
                "DEBUG",
            ]
        )
        self.assertEqual(parsed.count_kind, "all")
        self.assertEqual(parsed.log_level, "DEBUG")
        with self.assertRaises(SystemExit):
            parser.parse_args(["--ldscore", "x", "--counts", "m", "--sumstats-file", "trait.sumstats.gz"])
        with self.assertRaises(SystemExit):
            parser.parse_args(["--ldscore-dir", "ldscores", "--sumstats", "trait.sumstats.gz"])
        with self.assertRaises(SystemExit):
            parser.parse_args(["--ldscore-dir", "ldscores", "--sumstats-file", "trait.sumstats.gz", "--count-kind", "m_5_50"])

    def test_partitioned_h2_arguments_accept_per_query_output_flag(self):
        parser = argparse.ArgumentParser()
        regression_runner.add_partitioned_h2_arguments(parser)

        args = parser.parse_args(
            [
                "--ldscore-dir",
                "ldscores",
                "--sumstats-file",
                "trait.sumstats.gz",
                "--write-per-query-results",
            ]
        )

        self.assertTrue(args.write_per_query_results)

    def test_run_partitioned_h2_from_args_uses_query_columns_from_ldscore_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\nrs1\t1.0\t1000\n")
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            args = type(
                "Args",
                (),
                {
                    "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                    "trait_name": "trait",
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": None,
                    "n_blocks": 200,
                    "no_intercept": False,
                    "intercept_h2": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                    "write_per_query_results": False,
                },
            )()

            with mock.patch.object(
                regression_runner.RegressionRunner,
                "estimate_partitioned_h2_batch",
                return_value=pd.DataFrame([{"Category": "query", "Coefficient": 1.0}]),
            ) as patched:
                summary = regression_runner.run_partitioned_h2_from_args(args)

        patched.assert_called_once()
        self.assertEqual(patched.call_args.args[2].query_columns, ["query"])
        self.assertEqual(summary.loc[0, "Category"], "query")

    def test_run_partitioned_h2_from_args_writes_with_partitioned_writer(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\nrs1\t1.0\t1000\n")
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            output_dir = tmpdir / "out"
            args = type(
                "Args",
                (),
                {
                    "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                    "trait_name": "trait",
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": str(output_dir),
                    "overwrite": False,
                    "n_blocks": 200,
                    "no_intercept": False,
                    "intercept_h2": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                    "write_per_query_results": True,
                },
            )()

            with mock.patch.object(
                regression_runner.RegressionRunner,
                "estimate_partitioned_h2_batch",
                return_value=pd.DataFrame([{"Category": "query", "Coefficient": 1.0}]),
            ), mock.patch.object(
                regression_runner.PartitionedH2DirectoryWriter,
                "write",
            ) as writer:
                summary = regression_runner.run_partitioned_h2_from_args(args)

        writer.assert_called_once()
        output_config = writer.call_args.args[1]
        self.assertEqual(str(output_config.output_dir), str(output_dir))
        self.assertTrue(output_config.write_per_query_results)
        self.assertEqual(writer.call_args.kwargs["metadata"]["count_kind"], "common")
        self.assertEqual(writer.call_args.kwargs["metadata"]["trait_name"], "trait")
        self.assertEqual(summary.loc[0, "Category"], "query")

    def test_run_partitioned_h2_from_args_refuses_stale_per_query_tree_without_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "out"
            stale = output_dir / "query_annotations" / "old"
            stale.mkdir(parents=True)
            (stale / "metadata.json").write_text("{}\n", encoding="utf-8")
            args = type(
                "Args",
                (),
                {
                    "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                    "trait_name": "trait",
                    "ldscore_dir": str(tmpdir / "ldscores"),
                    "count_kind": "common",
                    "output_dir": str(output_dir),
                    "overwrite": False,
                    "n_blocks": 200,
                    "no_intercept": False,
                    "intercept_h2": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                    "write_per_query_results": False,
                },
            )()

            with mock.patch.object(
                regression_runner.RegressionRunner,
                "estimate_partitioned_h2_batch",
                side_effect=AssertionError("estimation should not run after preflight failure"),
            ):
                with self.assertRaisesRegex(FileExistsError, "overwrite"):
                    regression_runner.run_partitioned_h2_from_args(args)

            self.assertTrue(stale.exists())
            self.assertFalse((output_dir / "partitioned_h2.tsv").exists())
            self.assertFalse((output_dir / "partitioned-h2.log").exists())

    def test_run_partitioned_h2_from_args_overwrite_removes_stale_per_query_tree(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\nrs1\t1.0\t1000\n")
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            output_dir = tmpdir / "out"
            stale = output_dir / "query_annotations" / "old"
            stale.mkdir(parents=True)
            (stale / "metadata.json").write_text("{}\n", encoding="utf-8")
            args = type(
                "Args",
                (),
                {
                    "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                    "trait_name": "trait",
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": str(output_dir),
                    "overwrite": True,
                    "n_blocks": 200,
                    "no_intercept": False,
                    "intercept_h2": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                    "write_per_query_results": False,
                },
            )()
            partitioned_summary = pd.DataFrame(
                [
                    {
                        "Category": "query",
                        "Prop._SNPs": 1.0,
                        "Prop._h2": 1.0,
                        "Enrichment": 1.0,
                        "Enrichment_p": 0.5,
                        "Coefficient": 1.0,
                        "Coefficient_p": 0.5,
                    }
                ]
            )

            with mock.patch.object(
                regression_runner.RegressionRunner,
                "estimate_partitioned_h2_batch",
                return_value=partitioned_summary,
            ):
                summary = regression_runner.run_partitioned_h2_from_args(args)

            self.assertEqual(summary.loc[0, "Category"], "query")
            self.assertTrue((output_dir / "partitioned_h2.tsv").exists())
            self.assertTrue((output_dir / "partitioned-h2.log").exists())
            self.assertFalse((output_dir / "query_annotations").exists())

    def test_regression_cli_writes_fixed_result_filename_under_output_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\nrs1\t1.0\t1000\n")
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            output_dir = tmpdir / "out"
            args = type(
                "Args",
                (),
                {
                    "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                    "trait_name": "trait",
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": str(output_dir),
                    "n_blocks": 200,
                    "no_intercept": False,
                    "intercept_h2": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                },
            )()

            with mock.patch.object(
                regression_runner.RegressionRunner,
                "estimate_h2",
                return_value=mock.Mock(
                    tot=np.array([0.1]),
                    tot_se=np.array([0.01]),
                    intercept=np.array([1.0]),
                    intercept_se=0.01,
                    mean_chisq=np.array([1.1]),
                    lambda_gc=np.array([1.0]),
                    ratio=0.0,
                    ratio_se=0.0,
                ),
            ):
                regression_runner.run_h2_from_args(args)

            self.assertTrue((output_dir / "h2.tsv").exists())

    def test_regression_cli_refuses_existing_result_file_by_default(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\nrs1\t1.0\t1000\n")
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            output_dir = tmpdir / "out"
            output_dir.mkdir()
            existing = output_dir / "h2.tsv"
            existing.write_text("existing\n", encoding="utf-8")
            args = type(
                "Args",
                (),
                {
                    "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                    "trait_name": "trait",
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": str(output_dir),
                    "overwrite": False,
                    "n_blocks": 200,
                    "no_intercept": False,
                    "intercept_h2": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                },
            )()

            with mock.patch.object(
                regression_runner.RegressionRunner,
                "estimate_h2",
                return_value=mock.Mock(
                    tot=np.array([0.1]),
                    tot_se=np.array([0.01]),
                    intercept=np.array([1.0]),
                    intercept_se=0.01,
                    mean_chisq=np.array([1.1]),
                    lambda_gc=np.array([1.0]),
                    ratio=0.0,
                    ratio_se=0.0,
                ),
            ):
                with self.assertRaisesRegex(FileExistsError, "overwrite"):
                    regression_runner.run_h2_from_args(args)

            self.assertEqual(existing.read_text(encoding="utf-8"), "existing\n")
            self.assertFalse((output_dir / "h2.log").exists())

    def test_regression_cli_allows_existing_result_file_with_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\nrs1\t1.0\t1000\n")
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            output_dir = tmpdir / "out"
            output_dir.mkdir()
            existing = output_dir / "h2.tsv"
            existing.write_text("existing\n", encoding="utf-8")
            args = type(
                "Args",
                (),
                {
                    "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                    "trait_name": "trait",
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": str(output_dir),
                    "overwrite": True,
                    "n_blocks": 200,
                    "no_intercept": False,
                    "intercept_h2": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                },
            )()

            with mock.patch.object(
                regression_runner.RegressionRunner,
                "estimate_h2",
                return_value=mock.Mock(
                    tot=np.array([0.1]),
                    tot_se=np.array([0.01]),
                    intercept=np.array([1.0]),
                    intercept_se=0.01,
                    mean_chisq=np.array([1.1]),
                    lambda_gc=np.array([1.0]),
                    ratio=0.0,
                    ratio_se=0.0,
                ),
            ):
                regression_runner.run_h2_from_args(args)

            self.assertIn("total_h2", existing.read_text(encoding="utf-8"))

    def test_regression_writer_refuses_each_fixed_summary_filename(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            summary = pd.DataFrame({"value": [1]})
            for filename in ("h2.tsv", "partitioned_h2.tsv", "rg.tsv"):
                with self.subTest(filename=filename):
                    existing = output_dir / filename
                    existing.write_text("existing\n", encoding="utf-8")

                    with self.assertRaisesRegex(FileExistsError, "overwrite"):
                        regression_runner._maybe_write_dataframe(summary, str(output_dir), filename)

                    self.assertEqual(existing.read_text(encoding="utf-8"), "existing\n")
