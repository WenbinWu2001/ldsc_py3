from dataclasses import dataclass
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from ldscore.output import Artifact, ArtifactConfig, ArtifactProducer, OutputManager, OutputSpec, PostProcessor


@dataclass
class FakeChromResult:
    chrom: str
    reference_metadata: pd.DataFrame
    ld_scores: pd.DataFrame
    regression_metadata: pd.DataFrame
    w_ld: pd.DataFrame


@dataclass
class FakeResult:
    reference_metadata: pd.DataFrame
    ld_scores: pd.DataFrame
    regression_metadata: pd.DataFrame
    w_ld: pd.DataFrame
    snp_count_totals: dict[str, np.ndarray]
    baseline_columns: list[str]
    query_columns: list[str]
    chromosome_results: list[FakeChromResult]


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
    chrom_result = FakeChromResult("1", ref_meta, ld_scores, reg_meta, w_ld)
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
        chromosome_results=[chrom_result],
    )


class OutputManagerTest(unittest.TestCase):
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
