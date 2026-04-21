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
    reference_metadata: pd.DataFrame
    ld_scores: pd.DataFrame
    regression_metadata: pd.DataFrame
    w_ld: pd.DataFrame
    regression_snps: set[str]


@dataclass
class FakeResult:
    reference_metadata: pd.DataFrame
    ld_scores: pd.DataFrame
    regression_metadata: pd.DataFrame
    w_ld: pd.DataFrame
    snp_count_totals: dict[str, np.ndarray]
    baseline_columns: list[str]
    query_columns: list[str]
    regression_snps: set[str]
    chromosome_results: list[FakeChromResult]
    config_snapshot: GlobalConfig | None = None


class ExtraProducer(ArtifactProducer):
    name = "extra_note"

    def supports(self, result):
        return True

    def build(self, result, run_summary, output_spec, artifact_config=None):
        note = "extra"
        return [Artifact(self.name, "extra.txt", note, "text")]


def make_fake_result() -> FakeResult:
    ref_meta = pd.DataFrame({"CHR": [1, 1], "BP": [10, 20], "SNP": ["rs1", "rs2"], "CM": [0.1, 0.2]})
    ld_scores = pd.DataFrame({"base": [1.0, 2.0], "query": [0.5, 0.7]})
    reg_meta = ref_meta.iloc[[0]].reset_index(drop=True)
    w_ld = pd.DataFrame({"w_base": [1.5]})
    chrom_result = FakeChromResult("1", ref_meta, ld_scores, reg_meta, w_ld, {"rs1"})
    return FakeResult(
        reference_metadata=ref_meta,
        ld_scores=ld_scores,
        regression_metadata=reg_meta,
        w_ld=w_ld,
        snp_count_totals={
            "all_reference_snp_counts": np.array([10, 20]),
            "common_reference_snp_counts_maf_gt_0_05": np.array([8, 18]),
        },
        baseline_columns=["base"],
        query_columns=["query"],
        regression_snps={"rs1"},
        chromosome_results=[chrom_result],
    )


def make_multi_chrom_result() -> FakeResult:
    ref_meta = pd.DataFrame(
        {
            "CHR": [1, 1, 2],
            "BP": [10, 20, 30],
            "SNP": ["rs1", "rs2", "rs3"],
            "CM": [0.1, 0.2, 0.3],
        }
    )
    ld_scores = pd.DataFrame({"base": [1.0, 2.0, 3.0], "query": [0.5, 0.7, 0.9]})
    reg_meta = ref_meta.iloc[[0, 2]].reset_index(drop=True)
    w_ld = pd.DataFrame({"w_base": [1.5, 3.5]})
    chrom1 = FakeChromResult(
        "1",
        ref_meta.iloc[:2].reset_index(drop=True),
        ld_scores.iloc[:2].reset_index(drop=True),
        reg_meta.iloc[:1].reset_index(drop=True),
        w_ld.iloc[:1].reset_index(drop=True),
        {"rs1"},
    )
    chrom2 = FakeChromResult(
        "2",
        ref_meta.iloc[2:].reset_index(drop=True),
        ld_scores.iloc[2:].reset_index(drop=True),
        reg_meta.iloc[1:].reset_index(drop=True),
        w_ld.iloc[1:].reset_index(drop=True),
        {"rs3"},
    )
    return FakeResult(
        reference_metadata=ref_meta,
        ld_scores=ld_scores,
        regression_metadata=reg_meta,
        w_ld=w_ld,
        snp_count_totals={
            "all_reference_snp_counts": np.array([10, 20]),
            "common_reference_snp_counts_maf_gt_0_05": np.array([8, 18]),
        },
        baseline_columns=["base"],
        query_columns=["query"],
        regression_snps={"rs1", "rs3"},
        chromosome_results=[chrom1, chrom2],
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
            self.assertIn("ldscore", summary.output_paths)
            self.assertIn("w_ld", summary.output_paths)
            self.assertIn("counts.all_reference_snp_counts", summary.output_paths)
            self.assertIn("summary_tsv", summary.output_paths)
            self.assertIn("summary_json", summary.output_paths)
            self.assertIn("run_metadata", summary.output_paths)
            self.assertTrue(Path(summary.output_paths["ldscore"]).exists())

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
            self.assertTrue(Path(summary.output_paths["ldscore"]).exists())

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

    def test_regression_snps_path_filters_written_ldscore_and_weight_artifacts(self):
        result = make_multi_chrom_result()
        result.config_snapshot = GlobalConfig(snp_identifier="rsid", regression_snps_path="filters/hm3.txt")
        manager = OutputManager()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            spec = OutputSpec(
                out_prefix="example",
                output_dir=tmpdir,
                enabled_artifacts=["ldscore", "w_ld", "counts"],
            )

            summary = manager.write_outputs(result, spec)

            ldscore_df = pd.read_csv(summary.output_paths["ldscore"], sep="\t")
            self.assertEqual(ldscore_df["SNP"].tolist(), ["rs1", "rs3"])
            weight_df = pd.read_csv(summary.output_paths["w_ld"], sep="\t")
            self.assertEqual(weight_df["SNP"].tolist(), ["rs1", "rs3"])
            self.assertEqual(
                Path(summary.output_paths["counts.all_reference_snp_counts"]).read_text(encoding="utf-8"),
                "10\t20",
            )
            self.assertEqual(
                Path(summary.output_paths["counts.common_reference_snp_counts_maf_gt_0_05"]).read_text(encoding="utf-8"),
                "8\t18",
            )
            self.assertEqual(result.reference_metadata["SNP"].tolist(), ["rs1", "rs2", "rs3"])

    def test_regression_snps_path_filters_per_chrom_outputs_consistently(self):
        result = make_multi_chrom_result()
        result.config_snapshot = GlobalConfig(snp_identifier="rsid", regression_snps_path="filters/hm3.txt")
        manager = OutputManager()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            spec = OutputSpec(
                out_prefix="example",
                output_dir=tmpdir,
                enabled_artifacts=["ldscore", "w_ld"],
                write_per_chrom=True,
                aggregate_across_chromosomes=False,
                artifact_layout="by_chrom",
            )

            summary = manager.write_outputs(result, spec)

            self.assertIn("ldscore.chrom_1", summary.output_paths)
            self.assertIn("ldscore.chrom_2", summary.output_paths)
            self.assertIn("w_ld.chrom_1", summary.output_paths)
            self.assertIn("w_ld.chrom_2", summary.output_paths)
            ldscore_df = pd.read_csv(summary.output_paths["ldscore.chrom_1"], sep="\t")
            self.assertEqual(ldscore_df["SNP"].tolist(), ["rs1"])
            weight_df = pd.read_csv(summary.output_paths["w_ld.chrom_2"], sep="\t")
            self.assertEqual(weight_df["SNP"].tolist(), ["rs3"])

    def test_regression_snps_path_raises_when_no_written_rows_remain(self):
        result = make_multi_chrom_result()
        result.config_snapshot = GlobalConfig(snp_identifier="rsid", regression_snps_path="filters/hm3.txt")
        result.regression_snps = {"rs999"}
        for chrom_result in result.chromosome_results:
            chrom_result.regression_snps = {"rs999"}
        manager = OutputManager()
        with tempfile.TemporaryDirectory() as tmpdir:
            spec = OutputSpec(
                out_prefix="example",
                output_dir=tmpdir,
                enabled_artifacts=["ldscore"],
            )

            with self.assertRaisesRegex(ValueError, "After filtering to regression SNPs, no SNPs remain."):
                manager.write_outputs(result, spec)

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
