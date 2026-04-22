from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
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
from ldsc.outputs import Artifact, ArtifactConfig, ArtifactProducer, OutputManager, OutputSpec, PostProcessor


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


class OutputManagerTest(unittest.TestCase):
    def test_output_spec_normalizes_pathlike_fields(self):
        spec = OutputSpec(
            out_prefix=Path("results") / "example",
            output_dir=Path("artifacts"),
            log_path=Path("logs") / "run.log",
        )
        self.assertEqual(spec.out_prefix, "results/example")
        self.assertEqual(spec.output_dir, "artifacts")
        self.assertEqual(spec.log_path, "logs/run.log")

    def test_output_spec_expands_env_vars_but_does_not_glob_resolve(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            os.environ["LDSC_OUTPUT_ROOT"] = tmpdir
            spec = OutputSpec(out_prefix="$LDSC_OUTPUT_ROOT/results/*.out", output_dir="$LDSC_OUTPUT_ROOT/artifacts")
            self.assertEqual(spec.out_prefix, str(Path(tmpdir) / "results" / "*.out"))
            self.assertEqual(spec.output_dir, str(Path(tmpdir) / "artifacts"))

    def test_default_outputs(self):
        result = make_fake_result()
        manager = OutputManager()
        with tempfile.TemporaryDirectory() as tmpdir:
            spec = OutputSpec(out_prefix="example", output_dir=tmpdir)
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
            spec = OutputSpec(out_prefix="example", output_dir=output_dir)

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
            spec = OutputSpec(out_prefix="example", output_dir=output_dir)

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                manager.write_outputs(result, spec)

            self.assertEqual(caught, [])

    def test_missing_out_prefix_parent_warns_and_is_created(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_prefix = Path(tmpdir) / "nested" / "trait"
            munger = __import__("ldsc.sumstats_munger", fromlist=["SumstatsMunger"]).SumstatsMunger()
            table = __import__("ldsc.sumstats_munger", fromlist=["SumstatsTable"]).SumstatsTable(
                data=pd.DataFrame({"SNP": ["rs1"], "N": [1000.0], "Z": [1.5]}),
                has_alleles=False,
                source_path="source.sumstats.gz",
                trait_name="trait",
            )

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                written = munger.write_output(table, str(out_prefix))

            self.assertTrue((Path(tmpdir) / "nested").exists())
            self.assertEqual(len(caught), 1)
            self.assertIn("parent directory", str(caught[0].message).lower())
            self.assertIn("created", str(caught[0].message).lower())
            self.assertTrue(Path(written).exists())

    def test_explicit_enabled_artifacts(self):
        result = make_fake_result()
        manager = OutputManager()
        with tempfile.TemporaryDirectory() as tmpdir:
            spec = OutputSpec(
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
            spec = OutputSpec(out_prefix="example", output_dir=tmpdir, enabled_artifacts=["extra_note"])
            summary = manager.write_outputs(result, spec, artifact_config=ArtifactConfig())
            self.assertIn("extra_note", summary.output_paths)
            self.assertEqual(Path(summary.output_paths["extra_note"]).read_text(encoding="utf-8"), "extra")

    def test_per_chrom_output(self):
        result = make_fake_result()
        manager = OutputManager()
        with tempfile.TemporaryDirectory() as tmpdir:
            spec = OutputSpec(
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
            spec = OutputSpec(
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
            spec = OutputSpec(
                out_prefix="example",
                output_dir=tmpdir,
                enabled_artifacts=["ldscore"],
            )

            summary = manager.write_outputs(result, spec)

            self.assertEqual(set(summary.output_paths.keys()), {"ldscore.chrom_1", "ldscore.chrom_2"})
            chrom2_df = pd.read_csv(summary.output_paths["ldscore.chrom_2"], sep="\t")
            self.assertEqual(chrom2_df["SNP"].tolist(), ["rs3"])

    def test_unknown_artifact_raises(self):
        manager = OutputManager()
        spec = OutputSpec(out_prefix="example", enabled_artifacts=["does_not_exist"])
        with self.assertRaises(ValueError):
            manager.resolve_enabled_artifacts(spec)


class PostProcessorTest(unittest.TestCase):
    def test_registers_producer_on_manager(self):
        post_processor = PostProcessor()
        post_processor.register_producer(ExtraProducer())
        self.assertIn("extra_note", post_processor.output_manager.available_artifacts())
