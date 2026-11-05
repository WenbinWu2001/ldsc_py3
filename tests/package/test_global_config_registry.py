import inspect
from pathlib import Path
import re
import tempfile
import unittest
from unittest import mock

import pandas as pd
import pytest


import ldsc
from ldsc import GlobalConfig, set_global_config, reset_global_config
from ldsc.config import ConfigMismatchError, validate_config_compatibility
from ldsc._ldscore import calculator as ldscore_calculator
from ldsc._ldscore.calculator import LDScoreResult
from ldsc._runtime.workflow_kwargs import parse_command_line
from ldsc._refpanel import builder as r2_matrix_builder
from ldsc._refpanel.builder import R2MatrixBuildConfig, R2MatrixBuilder
from ldsc._regression import quantile_h2, runner as regression_runner
from ldsc._regression.runner import RegressionRunner
from ldsc._sumstats.table import SumstatsTable
from tests.runtime.test_path_resolution import _write_plink_trio


def _make_build_config(tmpdir: Path) -> R2MatrixBuildConfig:
    return R2MatrixBuildConfig(
        plink_prefix=tmpdir / "panel.@",
        source_genome_build="hg19",
        genetic_map_hg19_sources=tmpdir / "hg19.map",
        genetic_map_hg38_sources=tmpdir / "hg38.map",
        liftover_chain_hg19_to_hg38_file=tmpdir / "hg19ToHg38.over.chain",
        output_dir=tmpdir / "out",
        ld_wind_kb=1.0,
    )


def _make_sumstats_table() -> SumstatsTable:
    return SumstatsTable(
        data=pd.DataFrame(
            {
                "SNP": ["rs1", "rs2"],
                "Z": [2.0, 1.0],
                "N": [1000.0, 1000.0],
                "A1": ["A", "C"],
                "A2": ["G", "T"],
            }
        ),
        has_alleles=True,
        source_path="sumstats.gz",
        trait_name="trait",
        provenance={},
    )


def _make_ldscore_result() -> LDScoreResult:
    baseline_table = pd.DataFrame(
        {
            "CHR": ["1", "1"],
            "SNP": ["rs1", "rs2"],
            "BP": [10, 20],
            "regression_ld_scores": [3.0, 4.0],
            "base": [1.0, 2.0],
        }
    )
    return LDScoreResult(
        baseline_table=baseline_table,
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
    )


class GlobalConfigRegistryTest(unittest.TestCase):
    def tearDown(self) -> None:
        reset_global_config()

    def test_registered_global_config_defaults_to_chr_pos_allele_aware_auto(self):
        config = ldsc.get_global_config()

        self.assertEqual(config.snp_identifier, "chr_pos_allele_aware")
        self.assertEqual(config.genome_build, "auto")

    def test_run_annotate_uses_registered_global_config_and_logs_once(self):
        import tempfile
        from pathlib import Path

        set_global_config(GlobalConfig(snp_identifier="rsid", log_level="DEBUG"))
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            baseline, bed = root / "base.annot", root / "query.bed"
            baseline.write_text("CHR\tPOS\tSNP\tbase\n1\t10\trs1\t1\n")
            bed.write_text("1\t9\t10\n")
            with mock.patch("builtins.print") as patched_print, self.assertLogs("LDSC.config", level="INFO") as caught:
                result = ldsc.run_annotate(query_annot_bed_sources=[bed], ref_snp_annot_sources=[baseline], output_dir=root / "out")
                self.assertEqual(result.config_snapshot.snp_identifier, "rsid")
                self.assertIsNone(result.config_snapshot.genome_build)
                self.assertEqual(result.config_snapshot.log_level, "DEBUG")
            patched_print.assert_not_called()
            self.assertEqual(len(caught.records), 1)
            self.assertIn("GlobalConfig", caught.records[0].getMessage())

    def test_reference_panel_builder_uses_registered_global_config_and_logs_once(self):
        from tests.worker_fixtures import calculation_body
        self.enterContext(calculation_body())
        set_global_config(GlobalConfig(snp_identifier="rsid", log_level="DEBUG"))

        with tempfile.TemporaryDirectory() as tmpdir:
            config = _make_build_config(Path(tmpdir))
            Path(config.genetic_map_hg19_sources).write_text("CHR POS CM\n1 100 0.1\n")
            Path(config.genetic_map_hg38_sources).write_text("CHR POS CM\n1 100 0.1\n")
            builder = R2MatrixBuilder()
            _write_plink_trio(Path(tmpdir) / "panel.1", ["1"])

            with mock.patch.object(
                builder, "_prepare_build_state", return_value=mock.sentinel.state
            ), mock.patch.object(
                builder, "_check_chromosome_selections", side_effect=lambda chrom_sources, *_args: chrom_sources
            ), mock.patch.object(
                builder,
                "_build_chromosome",
                return_value={"r2_hg19": "r2-19", "r2_hg38": "r2-38", "meta_hg19": "m19", "meta_hg38": "m38"},
            ), mock.patch("builtins.print") as patched_print:
                with self.assertLogs("LDSC.config", level="INFO") as caught:
                    result = builder.run(config)

        self.assertEqual(builder.global_config, ldsc.get_global_config())
        self.assertEqual(result.chromosomes, ["1"])
        patched_print.assert_not_called()
        self.assertEqual(len(caught.records), 1)
        self.assertEqual(caught.records[0].levelname, "INFO")
        self.assertIn("GlobalConfig", caught.records[0].getMessage())

    def test_regression_runner_uses_registered_global_config_and_logs_once(self):
        set_global_config(GlobalConfig(snp_identifier="rsid"))
        runner = RegressionRunner()

        with mock.patch("builtins.print") as patched_print:
            with self.assertLogs("LDSC.config", level="INFO") as caught:
                dataset = runner.build_dataset(_make_sumstats_table(), _make_ldscore_result())

        self.assertEqual(runner.global_config, ldsc.get_global_config())
        self.assertEqual(dataset.merged["SNP"].tolist(), ["rs1", "rs2"])
        patched_print.assert_not_called()
        self.assertEqual(len(caught.records), 1)
        self.assertEqual(caught.records[0].levelname, "INFO")
        self.assertIn("GlobalConfig", caught.records[0].getMessage())

    def test_cli_normalization_ignores_registered_python_global_config(self):
        set_global_config(GlobalConfig(snp_identifier="chr_pos", genome_build="hg38", log_level="WARNING"))

        args = parse_command_line(ldscore_calculator.build_parser(), [
            "--snp-identifier", "rsid", "--log-level", "DEBUG", "--output-dir", "results/example",
            "--plink-prefix", "panel", "--ld-wind-kb", "1",
        ])

        global_config = ldscore_calculator._configs_from_args(args).global_config

        self.assertEqual(global_config.snp_identifier, "rsid")
        self.assertIsNone(global_config.genome_build)
        self.assertEqual(global_config.log_level, "DEBUG")


_ABSENT = "<absent>"
_REGISTERED_FIELDS = ("snp_identifier", "genome_build", "log_level")


def _identity_and_log_level(registered):
    return {name: getattr(registered, name) for name in _REGISTERED_FIELDS}


def _r2_matrix_global_config(registered):
    # build-r2-matrix records no registered build: "auto" in coordinate-family modes, None otherwise.
    chr_pos = registered.snp_identifier.startswith("chr_pos")
    return {**_identity_and_log_level(registered), "genome_build": "auto" if chr_pos else None}


def _log_level_only(registered):
    # Regression and post-fit workflows read SNP identity and build from their saved inputs.
    return {"snp_identifier": _ABSENT, "genome_build": _ABSENT, "log_level": registered.log_level}


# Every keyword wrapper: the workflow body it calls, the position of the namespace or GlobalConfig that
# body receives, plausible keywords (the body is replaced, so no input is read), and the registered
# values the body must receive.
_KEYWORD_WRAPPERS = {
    "run_ldscore": (ldscore_calculator, "run_ldscore_from_args", 0, {
        "baseline_annot_sources": "baseline.annot.gz", "plink_prefix": "panel", "ld_wind_cm": 1.0,
    }, _identity_and_log_level),
    "run_build_r2_matrix": (r2_matrix_builder, "_build", 1, {
        "plink_prefix": "panel.@", "ld_wind_kb": 1.0, "ref_panel_snps_file": "hm3.tsv",
    }, _r2_matrix_global_config),
    "run_h2": (regression_runner, "run_h2_from_args", 0, {
        "sumstats_file": "trait.parquet", "ldscore_dir": "ldscores",
    }, _log_level_only),
    "run_partitioned_h2": (regression_runner, "run_partitioned_h2_from_args", 0, {
        "sumstats_file": "trait.parquet", "ldscore_dir": "ldscores",
    }, _log_level_only),
    "run_rg": (regression_runner, "run_rg_from_args", 0, {
        "sumstats_sources": ["a.parquet", "b.parquet"], "ldscore_dir": "ldscores",
    }, _log_level_only),
    "run_quantile_h2": (quantile_h2, "run_quantile_h2_from_args", 0, {
        "partitioned_h2_result_dir": "fit", "target_annot_sources": ["target.annot.gz"],
        "target_annot_column": "score", "baseline_annot_sources": ["baseline.annot.gz"],
        "ref_metadata_sources": ["meta.tsv.gz"],
    }, _log_level_only),
}


@pytest.mark.parametrize("registered", [
    GlobalConfig(snp_identifier="rsid_allele_aware", log_level="DEBUG"),
    GlobalConfig(snp_identifier="chr_pos", genome_build="hg38", log_level="WARNING"),
], ids=["rsid_allele_aware", "chr_pos-hg38"])
@pytest.mark.parametrize("wrapper", sorted(_KEYWORD_WRAPPERS))
def test_keyword_wrappers_hand_their_workflow_the_registered_global_config(wrapper, registered, tmp_path):
    module, body, position, keywords, expected = _KEYWORD_WRAPPERS[wrapper]
    received = []
    set_global_config(registered)
    try:
        with mock.patch.object(module, body, side_effect=lambda *args, **kwargs: received.append(args[position])):
            getattr(ldsc, wrapper)(**keywords, output_dir=str(tmp_path / "out"))
    finally:
        reset_global_config()

    assert {name: getattr(received[0], name, _ABSENT) for name in _REGISTERED_FIELDS} == expected(registered)


def test_the_registry_test_lists_every_keyword_wrapper():
    wrappers = {
        name for name in ldsc.__all__
        if inspect.isfunction(getattr(ldsc, name))
        and re.search(r"\b(?:namespace|run)_from_kwargs\(", inspect.getsource(getattr(ldsc, name)))
    }
    assert wrappers == set(_KEYWORD_WRAPPERS)


class ConfigCompatibilityTest(unittest.TestCase):
    def test_validate_config_compatibility_rejects_genome_build_mismatch(self):
        left = GlobalConfig(genome_build="hg19", snp_identifier="chr_pos")
        right = GlobalConfig(genome_build="hg38", snp_identifier="chr_pos")

        with self.assertRaisesRegex(ConfigMismatchError, "genome_build mismatch"):
            validate_config_compatibility(left, right, context="test")

    def test_validate_config_compatibility_rejects_snp_identifier_mismatch(self):
        with self.assertWarns(UserWarning):
            cfg = GlobalConfig(genome_build="hg38", snp_identifier="rsid")
        self.assertIsNone(cfg.genome_build)
        left = cfg
        right = GlobalConfig(genome_build="hg38", snp_identifier="chr_pos")

        with self.assertRaisesRegex(ConfigMismatchError, "snp_identifier mismatch"):
            validate_config_compatibility(left, right, context="test")

    def test_validate_config_compatibility_accepts_matching_critical_fields(self):
        shared = GlobalConfig(genome_build="hg38", snp_identifier="chr_pos")

        validate_config_compatibility(shared, shared, context="test")
