"""User-facing help must explain every option and prioritize common inputs."""

import argparse

import pytest

from ldsc import cli


def test_old_build_ref_panel_command_is_rejected(capsys):
    with pytest.raises(SystemExit) as exc:
        cli.build_parser().parse_args(["build-ref-panel", "--help"])
    assert exc.value.code == 2
    assert "invalid choice" in capsys.readouterr().err


@pytest.fixture(scope="module")
def commands():
    parser = cli.build_parser()
    return next(action for action in parser._actions
                if isinstance(action, argparse._SubParsersAction)).choices


def test_every_visible_option_has_help_and_logging_levels(commands):
    for name, parser in commands.items():
        for action in parser._actions:
            if action.help == argparse.SUPPRESS:
                continue
            assert action.help, (name, action.option_strings)
        log_help = parser._option_string_actions["--log-level"].help
        for meaning in ("DEBUG", "INFO", "WARNING", "ERROR", "Default: INFO", "workflow log"):
            assert meaning in log_help, (name, meaning)


@pytest.mark.parametrize("command, primary", [
    ("munge-sumstats", "--raw-sumstats-file"),
    ("h2", "--sumstats-file"),
    ("partitioned-h2", "--sumstats-file"),
    ("rg", "--sumstats-sources"),
])
def test_essential_inputs_precede_advanced_controls(commands, command, primary):
    help_text = commands[command].format_help()
    assert "Inputs and output:" in help_text
    body = help_text.split("Inputs and output:", 1)[1]
    assert body.index(primary) < body.index("--log-level")
    if command == "munge-sumstats":
        assert body.index("--snp-identifier") < body.index("--chunksize")
        assert "Column overrides:" in body


def test_munge_companion_flags_and_window_alternatives_are_explained(commands):
    munging = commands["munge-sumstats"]._option_string_actions
    for flag, companion in [("--N-cas-col", "--N-con-col"), ("--N-con-col", "--N-cas-col")]:
        assert companion in munging[flag].help
        assert "--N-col" in munging[flag].help
    ldscore = commands["ldscore"]._option_string_actions
    for flag in ("--ld-wind-snps", "--ld-wind-kb", "--ld-wind-cm"):
        for alternative in ("--ld-wind-snps", "--ld-wind-kb", "--ld-wind-cm"):
            assert alternative in ldscore[flag].help
    assert "ref_metadata/chrN_meta.tsv.gz" in ldscore["--export-ref-metadata"].help
    assert "Default: off" in ldscore["--export-ref-metadata"].help


def test_wrapped_help_keeps_companion_flag_names_intact(commands):
    help_text = commands["ldscore"].format_help()
    r2_description = help_text.split("  --r2-dir DIR", 1)[1].split(
        "\n  --gene-ldscore-index-dir", 1
    )[0]
    assert "--gene-ldscore-index-dir" in r2_description


@pytest.mark.parametrize("command", ["h2", "partitioned-h2", "rg"])
def test_regression_rejects_removed_no_intercept(commands, command, capsys):
    inputs = (["--sumstats-sources", "one.parquet", "two.parquet"] if command == "rg"
              else ["--sumstats-file", "one.parquet"])
    with pytest.raises(SystemExit) as error:
        commands[command].parse_args([
            *inputs, "--ldscore-dir", "ldscores", "--output-dir", "results", "--no-intercept",
        ])
    assert error.value.code == 2
    assert "unrecognized arguments: --no-intercept" in capsys.readouterr().err
    assert "--no-intercept" not in commands[command].format_help()


def test_path_help_distinguishes_suites_globs_and_single_files(commands):
    for command, flag in (
        ("annotate", "--baseline-annot-sources"), ("ldscore", "--plink-prefix"),
        ("build-r2-panel", "--plink-prefix"), ("build-gene-ldscore-index", "--baseline-annot-sources"),
        ("quantile-h2", "--target-annot-sources"),
    ):
        help_text = commands[command]._option_string_actions[flag].help
        for meaning in ("'*'", "'@'", "1-22", "subset", "Quote"):
            assert meaning in help_text, (command, flag, meaning)
    for command in ("annotate", "ldscore", "quantile-h2"):
        for flag in ("--query-annot-bed-sources", "--query-annot-gene-list-sources"):
            help_text = commands[command]._option_string_actions[flag].help
            assert "'*'" in help_text
            assert "'@' is not expanded" in help_text
    for command, flag in (("h2", "--sumstats-file"), ("partitioned-h2", "--sumstats-file"),
                          ("munge-sumstats", "--raw-sumstats-file")):
        help_text = commands[command]._option_string_actions[flag].help
        assert "exactly one file" in help_text
        assert "'*'" in help_text


def test_ldscore_padding_help_warns_against_double_padding(commands):
    help_text = commands["ldscore"]._option_string_actions["--padding-bp"].help
    assert "Set to 0 if your BED intervals are already padded" in help_text
    assert "double padding" in help_text
