"""One LD-window rule shared by the three configurations that take a window."""

import numpy as np
import pytest

from ldsc._ldscore.calculator import build_parser, run_ldscore_from_args
from ldsc._refpanel.builder import build_parser as r2_build_parser
from ldsc.config import GeneLDScoreIndexBuildConfig, LDScoreConfig, LDWindow, R2MatrixBuildConfig
from ldsc.errors import LDSCConfigError

CHOICE = ("Cannot configure the LD window: set exactly one of --ld-wind-snps, --ld-wind-kb or --ld-wind-cm "
          "(ld_wind_snps, ld_wind_kb or ld_wind_cm in Python); ")
VALUE = "Use a finite number greater than 0; --ld-wind-snps takes a whole number of SNPs."

CONFIGS = {
    "ldscore": lambda **window: LDScoreConfig(**window),
    "r2_matrix": lambda **window: R2MatrixBuildConfig(plink_prefix="panel", output_dir="out", **window),
    "gene_index": lambda **window: GeneLDScoreIndexBuildConfig(
        baseline_annot_sources=("baseline",), plink_prefix="panel", output_dir="out", gene_coordinate_file="genes",
        genome_build="hg38", snp_identifier="rsid", **window),
}


@pytest.mark.parametrize("config", CONFIGS)
@pytest.mark.parametrize(("window", "found"), [
    ({}, "none was set."),
    ({"ld_wind_kb": 1.0, "ld_wind_snps": 5}, "--ld-wind-snps and --ld-wind-kb were both set."),
    ({"ld_wind_cm": 1.0, "ld_wind_kb": 1.0, "ld_wind_snps": 5},
     "--ld-wind-snps, --ld-wind-kb and --ld-wind-cm were all set."),
])
def test_every_config_requires_exactly_one_window(config, window, found):
    with pytest.raises(LDSCConfigError) as error:
        CONFIGS[config](**window)
    assert str(error.value) == CHOICE + found


@pytest.mark.parametrize("config", CONFIGS)
@pytest.mark.parametrize(("window", "shown"), [
    ({"ld_wind_snps": 0}, "--ld-wind-snps (ld_wind_snps in Python) is 0"),
    ({"ld_wind_kb": -1.0}, "--ld-wind-kb (ld_wind_kb in Python) is -1.0"),
    ({"ld_wind_kb": float("nan")}, "--ld-wind-kb (ld_wind_kb in Python) is nan"),
    ({"ld_wind_cm": float("inf")}, "--ld-wind-cm (ld_wind_cm in Python) is inf"),
    ({"ld_wind_snps": 2.5}, "--ld-wind-snps (ld_wind_snps in Python) is 2.5"),
    ({"ld_wind_snps": 5.0}, "--ld-wind-snps (ld_wind_snps in Python) is 5.0"),
    ({"ld_wind_snps": True}, "--ld-wind-snps (ld_wind_snps in Python) is True"),
    ({"ld_wind_kb": True}, "--ld-wind-kb (ld_wind_kb in Python) is True"),
    ({"ld_wind_kb": "1"}, "--ld-wind-kb (ld_wind_kb in Python) is '1'"),
])
def test_every_config_rejects_the_same_window_values(config, window, shown):
    with pytest.raises(LDSCConfigError) as error:
        CONFIGS[config](**window)
    assert str(error.value) == f"Cannot configure the LD window: {shown}. {VALUE}"


@pytest.mark.parametrize("config", CONFIGS)
@pytest.mark.parametrize(("window", "field", "expected"), [
    ({"ld_wind_snps": np.int64(5)}, "ld_wind_snps", LDWindow("snps", 5)),
    ({"ld_wind_kb": 100}, "ld_wind_kb", LDWindow("kb", 100.0)),
    ({"ld_wind_cm": np.float32(0.5)}, "ld_wind_cm", LDWindow("cm", 0.5)),
])
def test_window_values_are_stored_as_int_counts_and_float_distances(config, window, field, expected):
    built = CONFIGS[config](**window)
    assert built.ld_window == expected
    assert type(getattr(built, field)) is type(expected.value)
    assert type(built.ld_window.value) is type(expected.value)


def test_window_spells_its_option_and_record():
    window = LDScoreConfig(ld_wind_kb=100).ld_window
    assert window.flag == "--ld-wind-kb"
    assert str(window) == "--ld-wind-kb 100.0"
    assert str(LDWindow("snps", 50)) == "--ld-wind-snps 50"
    assert window.record() == {"unit": "kb", "value": 100.0}


def test_gene_index_batch_size_error_names_only_the_batch_sizes():
    with pytest.raises(LDSCConfigError) as error:
        CONFIGS["gene_index"](ld_wind_kb=1.0, snp_batch_size=0)
    assert "ld_wind" not in str(error.value)
    assert "snp_batch_size=0" in str(error.value)


@pytest.mark.parametrize(("command", "argv", "message"), [
    ("ldscore", ["--ld-wind-kb", "1", "--ld-wind-snps", "5"],
     "argument --ld-wind-snps: not allowed with argument --ld-wind-kb"),
    ("build-r2-matrix", ["--ld-wind-cm", "1", "--ld-wind-kb", "5"],
     "argument --ld-wind-kb: not allowed with argument --ld-wind-cm"),
    ("build-r2-matrix", [], "one of the arguments --ld-wind-cm --ld-wind-kb --ld-wind-snps is required"),
])
def test_commands_reject_window_choices_while_parsing(tmp_path, capsys, command, argv, message):
    parser = build_parser() if command == "ldscore" else r2_build_parser()
    with pytest.raises(SystemExit) as stop:
        parser.parse_args(["--plink-prefix", "panel", "--output-dir", str(tmp_path / "out"), *argv])
    assert stop.value.code == 2
    assert message in capsys.readouterr().err


def test_ldscore_command_without_a_window_fails_before_creating_the_output_directory(tmp_path):
    args = build_parser().parse_args(["--plink-prefix", str(tmp_path / "panel"), "--output-dir", str(tmp_path / "out")])
    with pytest.raises(LDSCConfigError, match="none was set"):
        run_ldscore_from_args(args)
    assert not (tmp_path / "out").exists()


def test_ldscore_command_rejects_a_bad_window_before_creating_the_output_directory(tmp_path):
    args = build_parser().parse_args(["--plink-prefix", str(tmp_path / "panel"), "--ld-wind-kb", "nan",
                                      "--output-dir", str(tmp_path / "out")])
    with pytest.raises(LDSCConfigError, match=r"--ld-wind-kb \(ld_wind_kb in Python\) is nan"):
        run_ldscore_from_args(args)
    assert not (tmp_path / "out").exists()
