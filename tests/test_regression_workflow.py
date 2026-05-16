import argparse
import contextlib
from dataclasses import replace
import gzip
import io
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
    from ldsc import cli
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

    def write_sumstats_sidecar(
        self,
        path: Path,
        *,
        trait_name: str | None = "trait",
        snp_identifier: str = "rsid",
        genome_build: str | None = None,
    ) -> None:
        path.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "artifact_type": "sumstats",
                    "trait_name": trait_name,
                    "snp_identifier": snp_identifier,
                    "genome_build": genome_build,
                }
            ),
            encoding="utf-8",
        )

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
                    "regression_ld_scores": [3.0],
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
                        "schema_version": 1,
                        "artifact_type": "ldscore",
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
                    }
                ),
                encoding="utf-8",
            )
            result = load_ldscore_from_dir(str(tmpdir))

        self.assertEqual(result.baseline_table.columns.tolist(), ["CHR", "SNP", "POS", "regression_ld_scores", "base"])
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
                    "regression_ld_scores": [3.0],
                    "base": [1.0],
                }
            ).to_parquet(tmpdir / "baseline.parquet", index=False)
            (tmpdir / "manifest.json").write_text(
                json.dumps(
                    {
                        "format": "ldsc.ldscore_result.v1",
                        "schema_version": 1,
                        "artifact_type": "ldscore",
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
                    }
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "baseline_table is missing required columns.*POS"):
                load_ldscore_from_dir(str(tmpdir))

    def test_load_ldscore_from_dir_reads_minimal_identity_manifest(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self.write_ldscore_dir(tmpdir, include_query=False)
            manifest_path = tmpdir / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                result = load_ldscore_from_dir(str(tmpdir))

        self.assertEqual(result.config_snapshot, GlobalConfig(snp_identifier="rsid"))
        self.assertFalse(any("GlobalConfig provenance is missing" in str(item.message) for item in caught))

    def test_load_ldscore_from_dir_rejects_missing_snp_identifier(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self.write_ldscore_dir(tmpdir, include_query=False)
            manifest_path = tmpdir / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            del manifest["snp_identifier"]
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "Unsupported snp_identifier mode"):
                load_ldscore_from_dir(str(tmpdir))

    def test_load_ldscore_from_dir_rejects_invalid_minimal_identity(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self.write_ldscore_dir(tmpdir, include_query=False)
            manifest_path = tmpdir / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["snp_identifier"] = "bad"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "Unsupported snp_identifier mode"):
                load_ldscore_from_dir(str(tmpdir))

    def test_load_ldscore_from_dir_rejects_allele_aware_baseline_without_alleles(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self.write_ldscore_dir(tmpdir, include_query=True)
            pd.read_parquet(tmpdir / "baseline.parquet").drop(columns=["A1", "A2"]).to_parquet(
                tmpdir / "baseline.parquet",
                index=False,
            )
            manifest_path = tmpdir / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["snp_identifier"] = "rsid_allele_aware"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "baseline_table.*A1/A2"):
                load_ldscore_from_dir(str(tmpdir))

    def test_load_ldscore_from_dir_rejects_allele_aware_query_without_alleles(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self.write_ldscore_dir(tmpdir, include_query=True)
            pd.read_parquet(tmpdir / "query.parquet").drop(columns=["A1", "A2"]).to_parquet(
                tmpdir / "query.parquet",
                index=False,
            )
            manifest_path = tmpdir / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["snp_identifier"] = "rsid_allele_aware"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "query_table.*A1/A2"):
                load_ldscore_from_dir(str(tmpdir))

    def test_load_ldscore_from_dir_rejects_snp_identifier_override_mismatch(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self.write_ldscore_dir(tmpdir, include_query=False)

            with self.assertRaisesRegex(ConfigMismatchError, "snp_identifier mismatch"):
                load_ldscore_from_dir(str(tmpdir), snp_identifier="chr_pos")

    def make_ldscore_result(self):
        baseline_table = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "SNP": ["rs1", "rs2", "rs3"],
                "POS": [10, 20, 30],
                "regression_ld_scores": [2.0, 2.0, 2.0],
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

    def make_sumstats_table_from_frame(self, frame, *, snp_identifier="rsid", trait_name="trait"):
        return SumstatsTable(
            data=frame,
            has_alleles={"A1", "A2"}.issubset(frame.columns),
            source_path=f"{trait_name}.sumstats.gz",
            trait_name=trait_name,
            provenance={},
            config_snapshot=GlobalConfig(snp_identifier=snp_identifier, genome_build="hg38"),
        )

    def make_ldscore_result_from_frame(self, frame, *, snp_identifier="rsid"):
        return LDScoreResult(
            baseline_table=frame,
            query_table=None,
            count_records=[
                {
                    "group": "baseline",
                    "column": "base",
                    "all_reference_snp_count": 10.0,
                    "common_reference_snp_count": 8.0,
                }
            ],
            baseline_columns=["base"],
            query_columns=[],
            ld_reference_snps=frozenset(),
            ld_regression_snps=frozenset(frame["SNP"].astype(str)),
            chromosome_results=[],
            config_snapshot=GlobalConfig(snp_identifier=snp_identifier, genome_build="hg38"),
        )

    def make_rg_kernel_result(self, *, rg=0.25, rg_se=0.05, z=5.0, p=0.01):
        hsq1 = mock.Mock(
            tot=np.array([0.20]),
            tot_se=np.array([0.02]),
            intercept=np.array([1.01]),
            intercept_se=0.01,
            mean_chisq=np.array([1.20]),
            lambda_gc=np.array([1.05]),
            ratio=0.05,
            ratio_se=0.01,
        )
        hsq2 = mock.Mock(
            tot=np.array([0.30]),
            tot_se=np.array([0.03]),
            intercept=np.array([1.02]),
            intercept_se=0.02,
            mean_chisq=np.array([1.30]),
            lambda_gc=np.array([1.06]),
            ratio=0.06,
            ratio_se=0.02,
        )
        gencov = mock.Mock(
            tot=np.array([0.10]),
            tot_se=np.array([0.01]),
            intercept=np.array([0.001]),
            intercept_se=0.0001,
        )
        return mock.Mock(
            rg_ratio=rg,
            rg_se=rg_se,
            z=z,
            p=p,
            hsq1=hsq1,
            hsq2=hsq2,
            gencov=gencov,
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
                "A1": ["A"],
                "A2": ["C"],
                "regression_ld_scores": [2.0],
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
                    "A1": ["A"],
                    "A2": ["C"],
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
                    "schema_version": 1,
                    "artifact_type": "ldscore",
                    "files": files,
                    "snp_identifier": "rsid",
                    "genome_build": "hg38",
                    "chromosomes": ["1"],
                    "baseline_columns": ["base"],
                    "query_columns": query_columns,
                    "counts": counts,
                }
            ),
            encoding="utf-8",
        )
        return root

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

        with self.assertRaisesRegex(ValueError, "Cannot mix SNP identifier families"):
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

        with self.assertRaisesRegex(ValueError, "allow-identity-downgrade"):
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
        self.assertIn("sumstats mode chr_pos_allele_aware", log_text)
        self.assertIn("effective snp_identifier='chr_pos'", log_text)
        self.assertIn("Dropped 4 duplicate effective-key rows", log_text)

    def test_h2_allele_aware_mode_orients_sumstats_z_to_ldscore_alleles(self):
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
        np.testing.assert_allclose(dataset.merged["Z"], [2.0, -3.0, 4.0, -5.0])

    def test_h2_allele_orientation_drops_only_truly_incompatible_pairs(self):
        merged = pd.DataFrame(
            {
                "SNP": ["same", "complement_same", "literal_reversed", "complement_reversed", "bad"],
                "A1": ["A", "A", "C", "C", "A"],
                "A2": ["C", "C", "A", "A", "C"],
                "A1_ld": ["A", "T", "A", "T", "A"],
                "A2_ld": ["C", "G", "C", "G", "G"],
                "Z": [1.0, 2.0, 3.0, 4.0, 5.0],
            }
        )

        with self.assertLogs("LDSC.regression_runner", level="WARNING") as caught:
            oriented = regression_runner._orient_sumstats_z_to_reference_alleles(merged)

        self.assertEqual(oriented["SNP"].tolist(), ["same", "complement_same", "literal_reversed", "complement_reversed"])
        np.testing.assert_allclose(oriented["Z"], [1.0, 2.0, -3.0, -4.0])
        self.assertIn("Dropping 1 SNPs with incompatible sumstats and LD-score allele order", "\n".join(caught.output))

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

        dataset = runner.build_dataset(sumstats, ldscore_result, query_columns=["query1"])

        self.assertEqual(dataset.merged["SNP"].tolist(), ["trait1", "trait2", "trait3"])
        self.assertEqual(dataset.retained_ld_columns, ["base", "query1"])

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

    def test_estimate_rg_pairs_returns_explicit_result_family_in_input_order(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        traits = [
            replace(self.make_sumstats_table(), trait_name="A", source_path="A.sumstats.gz"),
            replace(self.make_sumstats_table(), trait_name="B", source_path="B.sumstats.gz"),
            replace(self.make_sumstats_table(), trait_name="C", source_path="C.sumstats.gz"),
        ]

        with mock.patch.object(
            runner,
            "estimate_h2",
            return_value=self.make_rg_kernel_result().hsq1,
        ), mock.patch.object(
            regression_runner.reg,
            "RG",
            side_effect=[
                self.make_rg_kernel_result(rg=0.30, p=0.03),
                self.make_rg_kernel_result(rg=0.10, p=0.01),
                self.make_rg_kernel_result(rg=0.20, p=0.02),
            ],
        ):
            result = runner.estimate_rg_pairs(traits, self.make_ldscore_result())

        self.assertIsInstance(result, regression_runner.RgResultFamily)
        self.assertEqual(result.rg[["trait_1", "trait_2"]].values.tolist(), [["A", "B"], ["A", "C"], ["B", "C"]])
        self.assertEqual(result.rg_full["pair_kind"].tolist(), ["all_pairs", "all_pairs", "all_pairs"])
        self.assertEqual(result.rg_full["status"].tolist(), ["ok", "ok", "ok"])
        self.assertEqual(result.rg["n_snps_used"].tolist(), [3, 3, 3])
        self.assertEqual(result.h2_per_trait["trait_name"].tolist(), ["A", "B", "C"])
        self.assertEqual(len(result.per_pair_metadata), 3)
        self.assertEqual(result.per_pair_metadata[0]["effective_snp_identifier"], "rsid")
        self.assertFalse(result.per_pair_metadata[0]["identity_downgrade_applied"])

    def test_estimate_rg_pairs_anchor_mode_uses_anchor_against_rest_in_input_order(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        traits = [
            replace(self.make_sumstats_table(), trait_name="A", source_path="A.sumstats.gz"),
            replace(self.make_sumstats_table(), trait_name="B", source_path="B.sumstats.gz"),
            replace(self.make_sumstats_table(), trait_name="C", source_path="C.sumstats.gz"),
        ]

        with mock.patch.object(runner, "estimate_h2", return_value=self.make_rg_kernel_result().hsq1), mock.patch.object(
            regression_runner.reg,
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

        with mock.patch.object(runner, "estimate_h2", return_value=self.make_rg_kernel_result().hsq1), mock.patch.object(
            regression_runner.reg,
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

    def test_estimate_rg_pairs_converts_kernel_na_to_failed_row_but_preserves_numeric_out_of_range_rg(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        traits = [
            replace(self.make_sumstats_table(), trait_name="A", source_path="A.sumstats.gz"),
            replace(self.make_sumstats_table(), trait_name="B", source_path="B.sumstats.gz"),
            replace(self.make_sumstats_table(), trait_name="C", source_path="C.sumstats.gz"),
        ]

        with mock.patch.object(runner, "estimate_h2", return_value=self.make_rg_kernel_result().hsq1), mock.patch.object(
            regression_runner.reg,
            "RG",
            side_effect=[
                self.make_rg_kernel_result(rg="NA", rg_se="NA", z="NA", p="NA"),
                self.make_rg_kernel_result(rg=1.5, rg_se=0.1, z=15.0, p=1e-6),
                self.make_rg_kernel_result(rg=0.2, rg_se=0.1, z=2.0, p=0.04),
            ],
        ):
            result = runner.estimate_rg_pairs(traits, self.make_ldscore_result())

        self.assertEqual(result.rg_full.loc[0, "status"], "failed")
        self.assertTrue(np.isnan(result.rg.loc[0, "rg"]))
        self.assertIn("non-numeric", result.rg_full.loc[0, "error"])
        self.assertEqual(result.rg_full.loc[1, "status"], "ok")
        self.assertEqual(result.rg.loc[1, "rg"], 1.5)

    def test_rg_multiple_testing_ignores_nan_p_values(self):
        rg = pd.DataFrame({"p": [0.03, np.nan, 0.01, 0.02]})
        rg_full = pd.DataFrame({"p": [0.03, np.nan, 0.01, 0.02]})

        regression_runner._apply_rg_multiple_testing(rg, rg_full)

        np.testing.assert_allclose(rg["p_fdr_bh"].iloc[[0, 2, 3]], [0.03, 0.03, 0.03])
        self.assertTrue(np.isnan(rg.loc[1, "p_fdr_bh"]))
        np.testing.assert_allclose(rg_full["p_bonferroni"].iloc[[0, 2, 3]], [0.09, 0.03, 0.06])
        self.assertTrue(np.isnan(rg_full.loc[1, "p_bonferroni"]))

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
        self.assertEqual(result.per_query_metadata["query1"]["effective_snp_identifier"], "rsid")
        self.assertFalse(result.per_query_metadata["query1"]["identity_downgrade_applied"])
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
            self.write_sumstats_sidecar(tmpdir / "trait.metadata.json", trait_name="trait")
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
            metadata = json.loads((tmpdir / "h2_out" / "h2.metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(
                metadata,
                {
                    "schema_version": 1,
                    "artifact_type": "h2_result",
                    "trait_name": "trait",
                    "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                    "ldscore_dir": str(ldscore_dir),
                    "effective_snp_identifier": "rsid",
                    "genome_build": None,
                    "identity_downgrade_applied": False,
                    "count_key_used_for_regression": "common_reference_snp_counts",
                    "retained_ld_columns": ["base"],
                    "dropped_zero_variance_ld_columns": [],
                    "n_snps": 1,
                },
            )
            self.assertTrue((tmpdir / "h2_out" / "h2.log").exists())

    def test_run_h2_from_args_without_output_dir_creates_no_log_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\nrs1\t1.0\t1000\n")
            self.write_sumstats_sidecar(tmpdir / "trait.metadata.json", trait_name="trait")
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
            self.assertFalse((tmpdir / "h2.metadata.json").exists())

    def test_run_h2_from_args_uses_metadata_trait_name_when_cli_label_is_omitted(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\nrs1\t1.0\t1000\n")
            (tmpdir / "trait.metadata.json").write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "artifact_type": "sumstats",
                        "trait_name": "MDD",
                        "snp_identifier": "rsid",
                        "genome_build": None,
                    }
                ),
                encoding="utf-8",
            )
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            args = type(
                "Args",
                (),
                {
                    "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                    "trait_name": None,
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
                summary = regression_runner.run_h2_from_args(args)

        self.assertEqual(summary.loc[0, "trait_name"], "MDD")

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

    def test_rg_arguments_use_sumstats_sources_and_reject_old_pairwise_flags(self):
        parser = argparse.ArgumentParser()
        regression_runner.add_rg_arguments(parser)

        args = parser.parse_args(
            [
                "--ldscore-dir",
                "ldscores",
                "--sumstats-sources",
                "a.sumstats.gz",
                "b.sumstats.gz",
                "--intercept-h2",
                "1.02",
                "--intercept-gencov",
                "0.01",
            ]
        )

        self.assertEqual(args.sumstats_sources, ["a.sumstats.gz", "b.sumstats.gz"])
        self.assertEqual(args.intercept_h2, 1.02)
        self.assertEqual(args.intercept_gencov, 0.01)
        anchored = parser.parse_args(
            [
                "--ldscore-dir",
                "ldscores",
                "--sumstats-sources",
                "a.sumstats.gz",
                "b.sumstats.gz",
                "--anchor-trait",
                "MDD",
            ]
        )
        self.assertEqual(anchored.anchor_trait, "MDD")
        with self.assertRaises(SystemExit):
            parser.parse_args(
                [
                    "--ldscore-dir",
                    "ldscores",
                    "--sumstats-1-file",
                    "a.sumstats.gz",
                    "--sumstats-2-file",
                    "b.sumstats.gz",
                ]
            )
        with self.assertRaises(SystemExit):
            parser.parse_args(
                [
                    "--ldscore-dir",
                    "ldscores",
                    "--sumstats-sources",
                    "a.sumstats.gz",
                    "b.sumstats.gz",
                    "--anchor-trait-file",
                    "a.sumstats.gz",
                ]
            )

    def test_run_rg_from_args_rejects_per_pair_detail_without_output_dir_before_loading_inputs(self):
        args = type(
            "Args",
            (),
            {
                "sumstats_sources": ["a.sumstats.gz", "b.sumstats.gz"],
                "anchor_trait": None,
                "ldscore_dir": "ldscores",
                "count_kind": "common",
                "output_dir": None,
                "overwrite": False,
                "write_per_pair_detail": True,
                "n_blocks": 200,
                "no_intercept": False,
                "intercept_h2": None,
                "intercept_gencov": None,
                "two_step_cutoff": None,
                "chisq_max": None,
                "log_level": "INFO",
            },
        )()

        with mock.patch.object(
            regression_runner,
            "_load_sumstats_table",
            side_effect=AssertionError("sumstats should not load"),
        ):
            with self.assertRaisesRegex(ValueError, "--write-per-pair-detail requires --output-dir"):
                regression_runner.run_rg_from_args(args)

    def test_run_rg_from_args_returns_result_family_without_printing_when_output_dir_is_absent(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            for name in ("a", "b"):
                with gzip.open(tmpdir / f"{name}.sumstats.gz", "wt", encoding="utf-8") as handle:
                    handle.write("SNP\tZ\tN\nrs1\t1.0\t1000\n")
                self.write_sumstats_sidecar(tmpdir / f"{name}.metadata.json", trait_name=name)
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            expected = regression_runner.RgResultFamily(
                rg=pd.DataFrame([{"trait_1": "a", "trait_2": "b", "n_snps_used": 1, "rg": 0.1, "rg_se": 0.01, "p": 0.02, "p_fdr_bh": 0.02, "note": ""}]),
                rg_full=pd.DataFrame([{"trait_1": "a", "trait_2": "b", "n_snps_used": 1, "rg": 0.1, "rg_se": 0.01, "z": 10.0, "p": 0.02, "p_fdr_bh": 0.02, "p_bonferroni": 0.02, "status": "ok", "error": ""}]),
                h2_per_trait=pd.DataFrame([{"trait_name": "a"}, {"trait_name": "b"}]),
                per_pair_metadata=[],
            )
            args = type(
                "Args",
                (),
                {
                    "sumstats_sources": [str(tmpdir / "a.sumstats.gz"), str(tmpdir / "b.sumstats.gz")],
                    "anchor_trait": None,
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": None,
                    "overwrite": False,
                    "write_per_pair_detail": False,
                    "n_blocks": 200,
                    "no_intercept": False,
                    "intercept_h2": None,
                    "intercept_gencov": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                    "log_level": "INFO",
                },
            )()

            stdout = io.StringIO()
            with mock.patch.object(RegressionRunner, "estimate_rg_pairs", return_value=expected), contextlib.redirect_stdout(stdout):
                result = regression_runner.run_rg_from_args(args)

        self.assertIs(result, expected)
        self.assertEqual(stdout.getvalue(), "")

    def test_run_rg_from_args_recovers_metadata_labels_and_resolves_anchor_by_trait_name(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            for stem, trait_name in (("a", "MDD"), ("b", "SCZ")):
                with gzip.open(tmpdir / f"{stem}.sumstats.gz", "wt", encoding="utf-8") as handle:
                    handle.write("SNP\tZ\tN\nrs1\t1.0\t1000\n")
                (tmpdir / f"{stem}.metadata.json").write_text(
                    json.dumps(
                        {
                            "schema_version": 1,
                            "artifact_type": "sumstats",
                            "trait_name": trait_name,
                            "snp_identifier": "rsid",
                            "genome_build": None,
                        }
                    ),
                    encoding="utf-8",
                )
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            expected = regression_runner.RgResultFamily(
                rg=pd.DataFrame(),
                rg_full=pd.DataFrame(),
                h2_per_trait=pd.DataFrame(),
                per_pair_metadata=[],
            )
            args = type(
                "Args",
                (),
                {
                    "sumstats_sources": [str(tmpdir / "a.sumstats.gz"), str(tmpdir / "b.sumstats.gz")],
                    "anchor_trait": "SCZ",
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": None,
                    "overwrite": False,
                    "write_per_pair_detail": False,
                    "n_blocks": 200,
                    "no_intercept": False,
                    "intercept_h2": None,
                    "intercept_gencov": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                    "log_level": "INFO",
                },
            )()

            with mock.patch.object(RegressionRunner, "estimate_rg_pairs", return_value=expected) as patched:
                result = regression_runner.run_rg_from_args(args)

        self.assertIs(result, expected)
        tables = patched.call_args.args[0]
        self.assertEqual([table.trait_name for table in tables], ["MDD", "SCZ"])
        self.assertEqual(patched.call_args.kwargs["anchor_index"], 1)

    def test_run_rg_from_args_resolves_anchor_by_path_when_no_trait_name_matches(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            for name in ("a", "b"):
                with gzip.open(tmpdir / f"{name}.sumstats.gz", "wt", encoding="utf-8") as handle:
                    handle.write("SNP\tZ\tN\nrs1\t1.0\t1000\n")
                self.write_sumstats_sidecar(tmpdir / f"{name}.metadata.json", trait_name=name)
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            expected = regression_runner.RgResultFamily(
                rg=pd.DataFrame(),
                rg_full=pd.DataFrame(),
                h2_per_trait=pd.DataFrame(),
                per_pair_metadata=[],
            )
            args = type(
                "Args",
                (),
                {
                    "sumstats_sources": [str(tmpdir / "a.sumstats.gz"), str(tmpdir / "b.sumstats.gz")],
                    "anchor_trait": str(tmpdir / "b.sumstats.gz"),
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": None,
                    "overwrite": False,
                    "write_per_pair_detail": False,
                    "n_blocks": 200,
                    "no_intercept": False,
                    "intercept_h2": None,
                    "intercept_gencov": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                    "log_level": "INFO",
                },
            )()

            with mock.patch.object(RegressionRunner, "estimate_rg_pairs", return_value=expected) as patched:
                regression_runner.run_rg_from_args(args)

        self.assertEqual(patched.call_args.kwargs["anchor_index"], 1)

    def test_cli_prints_concise_rg_table_only_without_output_dir(self):
        expected = type(
            "Result",
            (),
            {
                "rg": pd.DataFrame(
                    [
                        {
                            "trait_1": "a",
                            "trait_2": "b",
                            "n_snps_used": 0,
                            "rg": np.nan,
                            "rg_se": np.nan,
                            "p": np.nan,
                            "p_fdr_bh": np.nan,
                            "note": "Failed",
                        }
                    ]
                )
            },
        )()

        stdout = io.StringIO()
        with mock.patch.object(regression_runner, "run_rg_from_args", return_value=expected), contextlib.redirect_stdout(stdout):
            result = cli.main(
                [
                    "rg",
                    "--ldscore-dir",
                    "ldscores",
                    "--sumstats-sources",
                    "a.sumstats.gz",
                    "b.sumstats.gz",
                ]
            )

        self.assertIs(result, expected)
        text = stdout.getvalue()
        self.assertIn("trait_1\ttrait_2\tn_snps_used\trg\trg_se\tp\tp_fdr_bh\tnote", text)
        self.assertIn("NaN", text)

    def test_cli_does_not_print_rg_table_when_output_dir_is_supplied(self):
        expected = type("Result", (), {"rg": pd.DataFrame([{"trait_1": "a"}])})()

        stdout = io.StringIO()
        with mock.patch.object(regression_runner, "run_rg_from_args", return_value=expected), contextlib.redirect_stdout(stdout):
            result = cli.main(
                [
                    "rg",
                    "--ldscore-dir",
                    "ldscores",
                    "--sumstats-sources",
                    "a.sumstats.gz",
                    "b.sumstats.gz",
                    "--output-dir",
                    "out",
                ]
            )

        self.assertIs(result, expected)
        self.assertEqual(stdout.getvalue(), "")

    def test_run_partitioned_h2_from_args_uses_query_columns_from_ldscore_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\nrs1\t1.0\t1000\n")
            self.write_sumstats_sidecar(tmpdir / "trait.metadata.json", trait_name="trait")
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
            self.write_sumstats_sidecar(tmpdir / "trait.metadata.json", trait_name="trait")
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

    def test_run_partitioned_h2_from_args_uses_metadata_trait_name_when_cli_label_is_omitted(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\nrs1\t1.0\t1000\n")
            (tmpdir / "trait.metadata.json").write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "artifact_type": "sumstats",
                        "trait_name": "MDD",
                        "snp_identifier": "rsid",
                        "genome_build": None,
                    }
                ),
                encoding="utf-8",
            )
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            output_dir = tmpdir / "out"
            args = type(
                "Args",
                (),
                {
                    "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                    "trait_name": None,
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
                regression_runner.run_partitioned_h2_from_args(args)

        self.assertEqual(writer.call_args.kwargs["metadata"]["trait_name"], "MDD")

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
            self.write_sumstats_sidecar(tmpdir / "trait.metadata.json", trait_name="trait")
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
            self.write_sumstats_sidecar(tmpdir / "trait.metadata.json", trait_name="trait")
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
            self.write_sumstats_sidecar(tmpdir / "trait.metadata.json", trait_name="trait")
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

    def test_regression_cli_refuses_existing_h2_metadata_by_default_before_estimation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "out"
            output_dir.mkdir()
            existing = output_dir / "h2.metadata.json"
            existing.write_text("{}\n", encoding="utf-8")
            args = type(
                "Args",
                (),
                {
                    "sumstats_file": str(tmpdir / "missing.sumstats.gz"),
                    "trait_name": "trait",
                    "ldscore_dir": str(tmpdir / "missing_ldscores"),
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
                side_effect=AssertionError("estimation should not run after preflight failure"),
            ):
                with self.assertRaisesRegex(FileExistsError, "overwrite"):
                    regression_runner.run_h2_from_args(args)

            self.assertEqual(existing.read_text(encoding="utf-8"), "{}\n")
            self.assertFalse((output_dir / "h2.tsv").exists())
            self.assertFalse((output_dir / "h2.log").exists())

    def test_regression_cli_allows_existing_result_file_with_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\nrs1\t1.0\t1000\n")
            self.write_sumstats_sidecar(tmpdir / "trait.metadata.json", trait_name="trait")
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            output_dir = tmpdir / "out"
            output_dir.mkdir()
            existing = output_dir / "h2.tsv"
            existing.write_text("existing\n", encoding="utf-8")
            existing_metadata = output_dir / "h2.metadata.json"
            existing_metadata.write_text('{"old": true}\n', encoding="utf-8")
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
            metadata = json.loads(existing_metadata.read_text(encoding="utf-8"))
            self.assertNotIn("old", metadata)
            self.assertEqual(metadata["artifact_type"], "h2_result")

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
