from __future__ import annotations

from dataclasses import dataclass, replace as dataclass_replace
from pathlib import Path
import json
import os
import sys
import tempfile
import unittest
import warnings

import numpy as np
import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc.config import GlobalConfig
from ldsc.ldscore_calculator import LDScoreResult
from ldsc.outputs import (
    Artifact,
    ArtifactConfig,
    ArtifactOutputConfig,
    ArtifactProducer,
    LDScoreDirectoryWriter,
    LDScoreOutputConfig,
    OutputManager,
    PostProcessor,
)


@dataclass
class FakeChromResult:
    chrom: str
    ldscore_table: pd.DataFrame
    snp_count_totals: dict[str, np.ndarray]
    baseline_columns: list[str]
    query_columns: list[str]
    ld_reference_snps: frozenset[str]
    ld_regression_snps: frozenset[str]


@dataclass
class FakeResult:
    ldscore_table: pd.DataFrame
    snp_count_totals: dict[str, np.ndarray]
    baseline_columns: list[str]
    query_columns: list[str]
    ld_reference_snps: frozenset[str]
    ld_regression_snps: frozenset[str]
    chromosome_results: list[FakeChromResult]
    config_snapshot: GlobalConfig | None = None


class ExtraProducer(ArtifactProducer):
    name = "extra_note"

    def supports(self, result):
        return True

    def build(self, result, run_summary, output_spec, artifact_config=None):
        return [Artifact(self.name, "extra.txt", "extra", "text")]


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


def _count_map() -> dict[str, np.ndarray]:
    return {
        "all_reference_snp_counts": np.array([10, 20]),
        "common_reference_snp_counts_maf_gt_0_05": np.array([8, 18]),
    }


def make_fake_result() -> FakeResult:
    ldscore_table = pd.DataFrame(
        {
            "CHR": ["1", "1"],
            "SNP": ["rs1", "rs2"],
            "BP": [10, 20],
            "base": [1.0, 2.0],
            "query": [0.5, 0.7],
            "regr_weight": [1.5, 2.5],
        }
    )
    chrom_result = FakeChromResult(
        chrom="1",
        ldscore_table=ldscore_table.copy(),
        snp_count_totals=_count_map(),
        baseline_columns=["base"],
        query_columns=["query"],
        ld_reference_snps=frozenset({"rs1", "rs2"}),
        ld_regression_snps=frozenset({"rs1", "rs2"}),
    )
    return FakeResult(
        ldscore_table=ldscore_table,
        snp_count_totals=_count_map(),
        baseline_columns=["base"],
        query_columns=["query"],
        ld_reference_snps=frozenset({"rs1", "rs2"}),
        ld_regression_snps=frozenset({"rs1", "rs2"}),
        chromosome_results=[chrom_result],
    )


def make_multi_chrom_result() -> FakeResult:
    chrom1_table = pd.DataFrame(
        {
            "CHR": ["1", "1"],
            "SNP": ["rs1", "rs2"],
            "BP": [10, 20],
            "base": [1.0, 2.0],
            "query": [0.5, 0.7],
            "regr_weight": [1.5, 2.5],
        }
    )
    chrom2_table = pd.DataFrame(
        {
            "CHR": ["2"],
            "SNP": ["rs3"],
            "BP": [30],
            "base": [3.0],
            "query": [0.9],
            "regr_weight": [3.5],
        }
    )
    chromosome_results = [
        FakeChromResult(
            chrom="1",
            ldscore_table=chrom1_table,
            snp_count_totals=_count_map(),
            baseline_columns=["base"],
            query_columns=["query"],
            ld_reference_snps=frozenset({"rs1", "rs2"}),
            ld_regression_snps=frozenset({"rs1", "rs2"}),
        ),
        FakeChromResult(
            chrom="2",
            ldscore_table=chrom2_table,
            snp_count_totals=_count_map(),
            baseline_columns=["base"],
            query_columns=["query"],
            ld_reference_snps=frozenset({"rs3"}),
            ld_regression_snps=frozenset({"rs3"}),
        ),
    ]
    return FakeResult(
        ldscore_table=pd.concat([chrom1_table, chrom2_table], axis=0, ignore_index=True),
        snp_count_totals=_count_map(),
        baseline_columns=["base"],
        query_columns=["query"],
        ld_reference_snps=frozenset({"rs1", "rs2", "rs3"}),
        ld_regression_snps=frozenset({"rs1", "rs2", "rs3"}),
        chromosome_results=chromosome_results,
    )


def make_single_table_result() -> FakeResult:
    ldscore_table = pd.DataFrame(
        {
            "CHR": ["1", "2"],
            "SNP": ["rs1", "rs2"],
            "BP": [10, 20],
            "base": [1.0, 2.0],
            "query": [0.5, 0.7],
            "regr_weight": [1.5, 2.5],
        }
    )
    return FakeResult(
        ldscore_table=ldscore_table,
        snp_count_totals=_count_map(),
        baseline_columns=["base"],
        query_columns=["query"],
        ld_reference_snps=frozenset({"rs1", "rs2"}),
        ld_regression_snps=frozenset({"rs1", "rs2"}),
        chromosome_results=[],
    )


class LDScoreDirectoryWriterTest(unittest.TestCase):
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


class OutputManagerTest(unittest.TestCase):
    def test_artifact_output_config_normalizes_pathlike_fields(self):
        spec = ArtifactOutputConfig(
            out_prefix=Path("results") / "example",
            output_dir=Path("artifacts"),
            log_path=Path("logs") / "run.log",
        )
        self.assertEqual(spec.out_prefix, "results/example")
        self.assertEqual(spec.output_dir, "artifacts")
        self.assertEqual(spec.log_path, "logs/run.log")

    def test_artifact_output_config_expands_env_vars_but_does_not_glob_resolve(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            os.environ["LDSC_OUTPUT_ROOT"] = tmpdir
            spec = ArtifactOutputConfig(
                out_prefix="$LDSC_OUTPUT_ROOT/results/*.out",
                output_dir="$LDSC_OUTPUT_ROOT/artifacts",
            )
            self.assertEqual(spec.out_prefix, str(Path(tmpdir) / "results" / "*.out"))
            self.assertEqual(spec.output_dir, str(Path(tmpdir) / "artifacts"))

    def test_default_outputs(self):
        result = make_fake_result()
        manager = OutputManager()
        with tempfile.TemporaryDirectory() as tmpdir:
            spec = ArtifactOutputConfig(out_prefix="example", output_dir=tmpdir)
            summary = manager.write_outputs(result, spec, config_snapshot={"stage": "test"})
            self.assertIn("ldscore.chrom_1", summary.output_paths)
            self.assertNotIn("w_ld", summary.output_paths)
            self.assertIn("counts.all_reference_snp_counts", summary.output_paths)
            self.assertIn("summary_tsv", summary.output_paths)
            self.assertIn("summary_json", summary.output_paths)
            self.assertIn("run_metadata", summary.output_paths)
            self.assertTrue(Path(summary.output_paths["ldscore.chrom_1"]).exists())

    def test_missing_output_dir_warns_and_is_created(self):
        result = make_fake_result()
        manager = OutputManager()
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "missing" / "artifacts"
            spec = ArtifactOutputConfig(out_prefix="example", output_dir=output_dir)

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                summary = manager.write_outputs(result, spec)

            self.assertTrue(output_dir.exists())
            self.assertEqual(len(caught), 1)
            self.assertIn("output directory", str(caught[0].message).lower())
            self.assertIn("created", str(caught[0].message).lower())
            self.assertTrue(Path(summary.output_paths["ldscore.chrom_1"]).exists())

    def test_existing_output_dir_does_not_warn(self):
        result = make_fake_result()
        manager = OutputManager()
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "artifacts"
            output_dir.mkdir(parents=True, exist_ok=True)
            spec = ArtifactOutputConfig(out_prefix="example", output_dir=output_dir)

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                manager.write_outputs(result, spec)

            self.assertEqual(caught, [])

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

    def test_explicit_enabled_artifacts(self):
        result = make_fake_result()
        manager = OutputManager()
        with tempfile.TemporaryDirectory() as tmpdir:
            spec = ArtifactOutputConfig(
                out_prefix="example",
                output_dir=tmpdir,
                enabled_artifacts=["summary_json"],
                write_ldscore=True,
                write_summary_tsv=True,
            )
            summary = manager.write_outputs(result, spec)
            self.assertEqual(set(summary.output_paths.keys()), {"summary_json"})

    def test_register_custom_producer(self):
        result = make_fake_result()
        manager = OutputManager()
        manager.register_producer(ExtraProducer())
        with tempfile.TemporaryDirectory() as tmpdir:
            spec = ArtifactOutputConfig(out_prefix="example", output_dir=tmpdir, enabled_artifacts=["extra_note"])
            summary = manager.write_outputs(result, spec, artifact_config=ArtifactConfig())
            self.assertIn("extra_note", summary.output_paths)
            self.assertEqual(Path(summary.output_paths["extra_note"]).read_text(encoding="utf-8"), "extra")

    def test_per_chrom_output(self):
        result = make_fake_result()
        manager = OutputManager()
        with tempfile.TemporaryDirectory() as tmpdir:
            spec = ArtifactOutputConfig(
                out_prefix="example",
                output_dir=tmpdir,
                enabled_artifacts=["ldscore"],
                write_per_chrom=True,
                aggregate_across_chromosomes=False,
                artifact_layout="by_chrom",
            )
            summary = manager.write_outputs(result, spec)
            self.assertIn("ldscore.chrom_1", summary.output_paths)
            self.assertTrue(Path(summary.output_paths["ldscore.chrom_1"]).exists())
            self.assertIn(f"{os.sep}chr1{os.sep}", summary.output_paths["ldscore.chrom_1"])

    def test_single_table_result_writes_aggregate_ldscore_artifact(self):
        result = make_single_table_result()
        manager = OutputManager()
        with tempfile.TemporaryDirectory() as tmpdir:
            spec = ArtifactOutputConfig(
                out_prefix="example",
                output_dir=tmpdir,
                enabled_artifacts=["ldscore"],
            )

            summary = manager.write_outputs(result, spec)

            self.assertEqual(set(summary.output_paths.keys()), {"ldscore"})
            ldscore_df = pd.read_csv(summary.output_paths["ldscore"], sep="\t")
            self.assertEqual(ldscore_df["SNP"].tolist(), ["rs1", "rs2"])

    def test_multi_chrom_result_writes_one_artifact_per_chromosome(self):
        result = make_multi_chrom_result()
        manager = OutputManager()
        with tempfile.TemporaryDirectory() as tmpdir:
            spec = ArtifactOutputConfig(
                out_prefix="example",
                output_dir=tmpdir,
                enabled_artifacts=["ldscore"],
            )

            summary = manager.write_outputs(result, spec)

            self.assertEqual(set(summary.output_paths.keys()), {"ldscore.chrom_1", "ldscore.chrom_2"})
            chrom2_df = pd.read_csv(summary.output_paths["ldscore.chrom_2"], sep="\t")
            self.assertEqual(chrom2_df["SNP"].tolist(), ["rs3"])

    def test_existing_artifact_raises_before_writing_any_outputs(self):
        result = make_multi_chrom_result()
        manager = OutputManager()
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            conflict_path = output_dir / "example.2.l2.ldscore.gz"
            conflict_path.write_text("existing", encoding="utf-8")
            spec = ArtifactOutputConfig(
                out_prefix="example",
                output_dir=output_dir,
                enabled_artifacts=["ldscore"],
            )

            with self.assertRaisesRegex(FileExistsError, "overwrite=True"):
                manager.write_outputs(result, spec)

            self.assertFalse((output_dir / "example.1.l2.ldscore.gz").exists())
            self.assertEqual(conflict_path.read_text(encoding="utf-8"), "existing")

    def test_unknown_artifact_raises(self):
        manager = OutputManager()
        spec = ArtifactOutputConfig(out_prefix="example", enabled_artifacts=["does_not_exist"])
        with self.assertRaises(ValueError):
            manager.resolve_enabled_artifacts(spec)


class PostProcessorTest(unittest.TestCase):
    def test_registers_producer_on_manager(self):
        post_processor = PostProcessor()
        post_processor.register_producer(ExtraProducer())
        self.assertIn("extra_note", post_processor.output_manager.available_artifacts())
