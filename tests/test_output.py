from __future__ import annotations

from dataclasses import replace as dataclass_replace
from pathlib import Path
import json
import sys
import tempfile
import unittest
import warnings

import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc.config import GlobalConfig
from ldsc.ldscore_calculator import LDScoreResult
from ldsc.outputs import LDScoreDirectoryWriter, LDScoreOutputConfig


def make_split_ldscore_result(query: bool = True) -> LDScoreResult:
    baseline_table = pd.DataFrame(
        {
            "CHR": ["1", "1"],
            "SNP": ["rs1", "rs2"],
            "BP": [10, 20],
            "regr_weight": [1.5, 2.5],
            "base": [1.0, 2.0],
        }
    )
    query_table = None
    query_columns: list[str] = []
    if query:
        query_table = pd.DataFrame(
            {
                "CHR": ["1", "1"],
                "SNP": ["rs1", "rs2"],
                "BP": [10, 20],
                "query": [0.5, 0.7],
            }
        )
        query_columns = ["query"]
    count_records = [
        {
            "group": "baseline",
            "column": "base",
            "all_reference_snp_count": 10.0,
            "common_reference_snp_count_maf_gt_0_05": 8.0,
        }
    ]
    if query:
        count_records.append(
            {
                "group": "query",
                "column": "query",
                "all_reference_snp_count": 20.0,
                "common_reference_snp_count_maf_gt_0_05": 18.0,
            }
        )
    return LDScoreResult(
        baseline_table=baseline_table,
        query_table=query_table,
        count_records=count_records,
        baseline_columns=["base"],
        query_columns=query_columns,
        ld_reference_snps=frozenset(),
        ld_regression_snps=frozenset({"rs1", "rs2"}),
        chromosome_results=[],
        config_snapshot=GlobalConfig(genome_build="hg38", snp_identifier="rsid"),
    )


class LDScoreDirectoryWriterTest(unittest.TestCase):
    def test_outputs_module_does_not_expose_prefix_based_artifact_pipeline(self):
        import ldsc.outputs as outputs

        removed_names = [
            "ArtifactConfig",
            "ArtifactOutputConfig",
            "ArtifactProducer",
            "OutputManager",
            "PostProcessor",
            "ResultFormatter",
            "ResultWriter",
        ]
        for name in removed_names:
            self.assertFalse(hasattr(outputs, name), f"{name} should not be exposed by ldsc.outputs")

    def test_writes_manifest_baseline_and_query_parquet(self):
        result = make_split_ldscore_result(query=True)
        writer = LDScoreDirectoryWriter()
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "ldscores"
            output_paths = writer.write(result, LDScoreOutputConfig(output_dir=output_dir))

            self.assertEqual(set(output_paths), {"manifest", "baseline", "query"})
            manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["format"], "ldsc.ldscore_result.v1")
            self.assertEqual(manifest["files"], {"baseline": "baseline.parquet", "query": "query.parquet"})
            self.assertEqual(manifest["baseline_columns"], ["base"])
            self.assertEqual(manifest["query_columns"], ["query"])
            self.assertEqual(manifest["counts"][1]["column"], "query")

            baseline = pd.read_parquet(output_dir / "baseline.parquet")
            query = pd.read_parquet(output_dir / "query.parquet")
            self.assertEqual(baseline.columns.tolist(), ["CHR", "SNP", "BP", "regr_weight", "base"])
            self.assertEqual(query.columns.tolist(), ["CHR", "SNP", "BP", "query"])

    def test_omits_query_parquet_for_baseline_only_result(self):
        result = make_split_ldscore_result(query=False)
        writer = LDScoreDirectoryWriter()
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "ldscores"
            output_paths = writer.write(result, LDScoreOutputConfig(output_dir=output_dir))

            self.assertEqual(set(output_paths), {"manifest", "baseline"})
            self.assertFalse((output_dir / "query.parquet").exists())
            manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["files"], {"baseline": "baseline.parquet"})
            self.assertEqual(manifest["query_columns"], [])

    def test_rejects_query_table_with_mismatched_row_keys(self):
        result = make_split_ldscore_result(query=True)
        bad_query = result.query_table.copy()
        bad_query.loc[1, "BP"] = 999
        bad_result = dataclass_replace(result, query_table=bad_query)

        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(ValueError, "query rows must match baseline rows"):
                LDScoreDirectoryWriter().write(bad_result, LDScoreOutputConfig(output_dir=tmpdir))

    def test_refuses_overwrite_before_writing_any_file(self):
        result = make_split_ldscore_result(query=True)
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            (output_dir / "baseline.parquet").write_text("existing", encoding="utf-8")

            with self.assertRaisesRegex(FileExistsError, "overwrite=True"):
                LDScoreDirectoryWriter().write(result, LDScoreOutputConfig(output_dir=output_dir))

            self.assertFalse((output_dir / "manifest.json").exists())
            self.assertEqual((output_dir / "baseline.parquet").read_text(encoding="utf-8"), "existing")

    def test_omits_missing_common_count_values_from_manifest_records(self):
        result = make_split_ldscore_result(query=False)
        result = dataclass_replace(
            result,
            count_records=[
                {
                    "group": "baseline",
                    "column": "base",
                    "all_reference_snp_count": 10.0,
                }
            ],
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "ldscores"
            LDScoreDirectoryWriter().write(result, LDScoreOutputConfig(output_dir=output_dir))

            manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))

        self.assertEqual(manifest["counts"][0]["all_reference_snp_count"], 10.0)
        self.assertNotIn("common_reference_snp_count_maf_gt_0_05", manifest["counts"][0])


class FixedOutputDirectoryTest(unittest.TestCase):
    def test_missing_munge_output_dir_warns_and_is_created(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "nested" / "trait"
            munger = __import__("ldsc.sumstats_munger", fromlist=["SumstatsMunger"]).SumstatsMunger()
            table = __import__("ldsc.sumstats_munger", fromlist=["SumstatsTable"]).SumstatsTable(
                data=pd.DataFrame({"SNP": ["rs1"], "N": [1000.0], "Z": [1.5]}),
                has_alleles=False,
                source_path="source.sumstats.gz",
                trait_name="trait",
            )

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                written = munger.write_output(table, str(output_dir))

            self.assertTrue(output_dir.exists())
            self.assertEqual(len(caught), 1)
            self.assertIn("output directory", str(caught[0].message).lower())
            self.assertIn("created", str(caught[0].message).lower())
            self.assertTrue(Path(written).exists())
