"""Keyword workflows reject what the command line would reject, before writing anything."""

import argparse
import re
from types import SimpleNamespace

import pytest

import ldsc
from ldsc._runtime.workflow_kwargs import GLOBAL_CONFIG_RUNTIME_OPTIONS, namespace_from_kwargs, parse_command_line
from ldsc.errors import LDSCUsageError

# Plausible inputs; validation fails before any of them is read.
_INPUTS = {
    "run_ldscore": {"plink_prefix": "panel", "ld_wind_cm": 1.0},
    "run_build_r2_matrix": {"plink_prefix": "panel.@", "ld_wind_kb": 1.0},
}
_SHARED = "set_global_config(GlobalConfig("
# Each case would otherwise run with a default or a truthy string in place of the intended value.
_INVALID_KEYWORDS = [
    ("run_ldscore", {"ld_wind_cn": 1.0}, "unknown keyword argument(s) ld_wind_cn"),
    ("run_ldscore", {"bfile": "panel"}, "unknown keyword argument(s) bfile"),
    ("run_ldscore", {"regr_snps_exclude_regions": "MHC"},
     "regr_snps_exclude_regions='MHC' is not one of 'none', 'mhc', 'centromeres', 'mhc-and-centromeres'"),
    ("run_ldscore", {"yes_really": "no"}, "yes_really='no' must be False or True"),
    ("run_ldscore", {"overwrite": "False"}, "overwrite='False' must be False or True"),
    ("run_ldscore", {"gene_list_resolution_policy": "resolved_only"},
     "gene_list_resolution_policy='resolved_only' must be 'strict' or 'resolved-only'"),
    ("run_ldscore", {"snp_identifier": "rsid"}, _SHARED + "snp_identifier="),
    ("run_ldscore", {"genome_build": "hg19"}, _SHARED + "genome_build="),
    ("run_ldscore", {"log_level": "DEBUG"}, _SHARED + "log_level="),
    ("run_build_r2_matrix", {"maf": 0.01}, "unknown keyword argument(s) maf"),
    ("run_build_r2_matrix", {"source_genome_build": "HG19"}, "source_genome_build='HG19' is not one of"),
    ("run_build_r2_matrix", {"overwrite": "False"}, "overwrite='False' must be False or True"),
    ("run_build_r2_matrix", {"snp_identifier": "rsid"}, _SHARED + "snp_identifier="),
    ("run_build_r2_matrix", {"log_level": "DEBUG"}, _SHARED + "log_level="),
    ("run_build_r2_matrix", {"genome_build": "hg19"}, "pass source_genome_build"),
]


@pytest.mark.parametrize(("function", "extra", "message"), _INVALID_KEYWORDS)
def test_keyword_workflows_reject_invalid_keywords_before_writing(tmp_path, function, extra, message):
    output_dir = tmp_path / "out"
    with pytest.raises(LDSCUsageError, match=re.escape(f"{function}() cannot run: ") + ".*" + re.escape(message)):
        getattr(ldsc, function)(**_INPUTS[function], output_dir=str(output_dir), **extra)
    assert not output_dir.exists()


@pytest.mark.parametrize(("function", "required"), [
    ("run_ldscore", ["output_dir"]),
    ("run_build_r2_matrix", ["output_dir", "plink_prefix"]),
])
def test_keyword_workflows_list_every_missing_required_keyword(function, required):
    with pytest.raises(LDSCUsageError, match=re.escape("missing required keyword argument(s)")) as raised:
        getattr(ldsc, function)()
    listed = str(raised.value).split("missing required keyword argument(s) ")[1].split(";")[0].rstrip(".").split(", ")
    assert sorted(listed) == required


def _command_parser():
    """A command parser with identity, log-level and ordinary options."""
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument("--snp-identifier", default="chr_pos_allele_aware")
    parser.add_argument("--genome-build", default=None)
    parser.add_argument("--log-level", default="INFO")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--maf-min", type=float, default=None)
    return parser


# Any object with the GlobalConfig attributes; the helper reads it without importing ldsc.config.
_REGISTERED = SimpleNamespace(snp_identifier="rsid_allele_aware", genome_build="hg19", log_level="DEBUG")


def test_rejected_options_take_the_registered_value_instead_of_the_parser_default():
    args = namespace_from_kwargs(
        "run_x", _command_parser(), {"output_dir": "out"}, rejected=GLOBAL_CONFIG_RUNTIME_OPTIONS,
        registered=_REGISTERED,
    )
    assert (args.snp_identifier, args.genome_build, args.log_level) == ("rsid_allele_aware", "hg19", "DEBUG")
    assert (args.output_dir, args.maf_min) == ("out", None)


def test_a_rejected_option_without_a_registered_value_is_absent():
    # Reading it fails instead of returning the parser default.
    args = namespace_from_kwargs(
        "run_x", _command_parser(), {"output_dir": "out"},
        rejected={**GLOBAL_CONFIG_RUNTIME_OPTIONS, "maf_min": "the threshold comes from elsewhere"},
        registered=_REGISTERED,
    )
    assert not hasattr(args, "maf_min")
    assert args.snp_identifier == "rsid_allele_aware"


def test_keyword_and_command_line_routes_record_the_same_explicit_options():
    parser = _command_parser()
    from_kwargs = namespace_from_kwargs(
        "run_x", parser, {"output_dir": "out", "maf_min": 0.01}, rejected=GLOBAL_CONFIG_RUNTIME_OPTIONS,
        registered=_REGISTERED,
    )
    from_argv = parse_command_line(parser, ["--output-dir", "out", "--maf-min=0.01"])
    assert from_kwargs._explicit_cli_options == from_argv._explicit_cli_options == {"--output-dir", "--maf-min"}
    assert from_argv.maf_min == 0.01


# Every command with a keyword wrapper, as `ldsc` parses it; the workflow is replaced by a recorder.
_COMMAND_LINES = [
    ("ldsc._ldscore.calculator", "run_ldscore_from_args", ["ldscore", "--plink-prefix", "panel", "--output-dir", "out"]),
    ("ldsc._refpanel.builder", "run_build_r2_matrix_from_args",
     ["build-r2-matrix", "--plink-prefix", "panel", "--ld-wind-kb=1", "--output-dir", "out"]),
    ("ldsc._regression.runner", "run_h2_from_args",
     ["h2", "--sumstats-file", "trait.parquet", "--ldscore-dir", "ld", "--output-dir", "out"]),
    ("ldsc._regression.runner", "run_partitioned_h2_from_args",
     ["partitioned-h2", "--sumstats-file", "trait.parquet", "--ldscore-dir", "ld", "--output-dir=out"]),
    ("ldsc._regression.runner", "run_rg_from_args",
     ["rg", "--sumstats-sources", "a.parquet", "b.parquet", "--ldscore-dir", "ld", "--output-dir", "out"]),
    ("ldsc._regression.quantile_h2", "run_quantile_h2_from_args",
     ["quantile-h2", "--partitioned-h2-result-dir", "fit", "--target-annot-sources", "target.annot",
      "--target-annot-column", "target", "--baseline-annot-sources", "base.annot", "--ref-metadata-sources",
      "meta.tsv", "--output-dir", "out"]),
]


@pytest.mark.parametrize(("module", "run", "argv"), _COMMAND_LINES)
def test_command_lines_record_explicit_options_like_their_keyword_wrappers(monkeypatch, module, run, argv):
    from importlib import import_module

    from ldsc import cli

    parsed = []
    monkeypatch.setattr(import_module(module), run, lambda args: parsed.append(args) or SimpleNamespace(output_paths={}))
    cli.main(argv)
    assert parsed[0]._explicit_cli_options == {token.split("=", 1)[0] for token in argv if token.startswith("--")}
