"""User-facing help must explain every option and list flags in guided step order."""

import argparse
import re

import pytest

from ldsc import cli


# Small commands keep plain section headings instead of numbered steps.
UNNUMBERED_COMMANDS = {"infer-genome-build", "convert-h2-scale", "plot"}

LDSCORE_SECTIONS = [
    ("Step 1: SNP identity and genome build", ["--snp-identifier", "--genome-build"]),
    ("Step 2: Baseline annotations", ["--baseline-annot-sources"]),
    ("Step 3: Query annotations (optional; choose at most one source type)", [
        "--query-annot-bed-sources", "--query-annot-gene-list-sources", "--query-annot-sources",
        "--padding-bp", "--gene-coordinate-file", "--control-gene-list-file",
        "--gene-exclude-regions", "--allow-unresolved-genes",
    ]),
    ("Step 4: Reference panel (choose one: R2 matrix or PLINK)", ["--r2-matrix-dir", "--plink-prefix"]),
    ("Step 4 continued: PLINK reference panel only (optional)", [
        "--keep-indivs-file", "--genetic-map-hg19-sources", "--genetic-map-hg38-sources",
        "--export-ref-metadata",
    ]),
    ("Step 5: LD window (choose one)", ["--ld-wind-cm", "--ld-wind-kb", "--ld-wind-snps", "--yes-really"]),
    ("Step 6: Output", ["--output-dir", "--overwrite"]),
    ("Optional: regression SNPs and reference SNPs", [
        "--regr-snps-file", "--regr-snps-exclude-regions", "--ref-panel-snps-file", "--maf-min",
        "--common-maf-min",
    ]),
    ("Alternative: gene LD-score index", ["--gene-ldscore-index-dir"]),
    ("Computing resources", ["--workers", "--query-batch-size", "--snp-batch-size"]),
    ("Logging", ["--log-level"]),
]

MUNGE_SECTIONS = [
    ("Step 1: SNP identity and genome build", [
        "--snp-identifier", "--source-genome-build", "--output-genome-build", "--liftover-chain-file",
    ]),
    ("Step 2: Raw summary statistics", ["--raw-sumstats-file", "--trait-name", "--input-format"]),
    ("Step 3: Sample size (only if not detected from the file)", [
        "--N-col", "--N-cas-col", "--N-con-col", "--N", "--N-cas", "--N-con",
    ]),
    ("Step 4: Output", ["--output-dir", "--overwrite", "--output-format", "--infer-only"]),
    ("Optional: SNP selection and quality filters", [
        "--sumstats-snps-file", "--no-snp-restriction", "--N-min", "--nstudy-min", "--info-min", "--maf-min",
    ]),
    ("Optional: column overrides", [
        "--snp", "--chr", "--pos", "--a1", "--a2", "--p", "--signed-sumstats", "--a1-inc", "--frq",
        "--info", "--info-list", "--nstudy", "--ignore",
    ]),
    ("Computing resources", ["--chunksize"]),
    ("Logging", ["--log-level"]),
]


def _help_sections(parser):
    """Return rendered help sections as (title, flags) pairs in display order."""
    lines = parser.format_help().splitlines()
    sections = []
    for line, following in zip(lines, [*lines[1:], ""]):
        # A heading starts at column 0, ends with a colon, and is followed by
        # indented section content; wrapped command-description text is not.
        if (line.endswith(":") and not line[:1].isspace() and not line.startswith("usage:")
                and following.startswith("  ")):
            sections.append((line[:-1], []))
        elif sections:
            sections[-1][1].append(line)
    result = []
    for title, content in sections:
        paragraphs = [text for text in "\n".join(content).split("\n\n") if text.strip()]
        # A group description is its own paragraph before the options, and a
        # wrapped description line may itself start with a flag name.
        options = paragraphs[1:] if len(paragraphs) > 1 else paragraphs
        flags = [match.group(1) for text in options for match in re.finditer(r"^  (-[-\w]+)", text, re.M)]
        if title not in ("options", "positional arguments"):
            result.append((title, flags))
    return result


@pytest.fixture(scope="module")
def commands():
    return {name: cli._command_parser(name) for name in cli._SUBCOMMAND_HELP}


@pytest.mark.parametrize("command, inputs, default", [
    ("annotate", ["--ref-snp-annot-sources", "baseline.annot.gz",
                  "--query-annot-gene-list-sources", "genes.txt"], None),
    ("ldscore", [], "strict"),
    ("quantile-h2", ["--partitioned-h2-result-dir", "results",
                     "--target-annot-sources", "target.annot.gz",
                     "--target-annot-column", "target",
                     "--baseline-annot-sources", "baseline.annot.gz",
                     "--ref-metadata-sources", "ref.tsv.gz"], "strict"),
])
def test_allow_unresolved_genes_is_an_opt_in_flag(commands, command, inputs, default, capsys):
    parser = commands[command]
    argv = [*inputs, "--output-dir", "out"]
    assert parser.parse_args(argv).gene_list_resolution_policy == default
    assert parser.parse_args([*argv, "--allow-unresolved-genes"]).gene_list_resolution_policy == "resolved-only"
    with pytest.raises(SystemExit) as error:
        parser.parse_args([*argv, "--allow-unresolved-genes", "resolved-only"])
    assert error.value.code == 2
    assert "unrecognized arguments" in capsys.readouterr().err
    assert "--allow-unresolved-genes" in parser.format_help()


def test_every_visible_option_has_help_and_logging_levels(commands):
    for name, parser in commands.items():
        for action in parser._actions:
            if action.help == argparse.SUPPRESS:
                continue
            assert action.help, (name, action.option_strings)
        if name == "infer-genome-build":
            assert set(parser._option_string_actions) == {"-h", "--help"}
            continue
        log_help = parser._option_string_actions["--log-level"].help
        for meaning in ("DEBUG", "INFO", "WARNING", "ERROR", "Default: INFO", "workflow log"):
            assert meaning in log_help, (name, meaning)


@pytest.mark.parametrize("command, expected", [
    ("ldscore", LDSCORE_SECTIONS),
    ("munge-sumstats", MUNGE_SECTIONS),
])
def test_help_lists_flags_in_guided_step_order(commands, command, expected):
    assert _help_sections(commands[command]) == expected


def test_help_puts_numbered_steps_first_and_logging_last(commands):
    for name, parser in commands.items():
        if name == "infer-genome-build":
            continue
        sections = _help_sections(parser)
        titles = [title for title, _ in sections]
        section_of = {flag: title for title, flags in sections for flag in flags}
        assert sections[-1] == ("Logging", ["--log-level"]), name
        if "Computing resources" in titles:
            assert titles[-2] == "Computing resources", name
        output_section = section_of.get("--output-dir", "Output")
        assert section_of["--overwrite"] == output_section, name
        if name in UNNUMBERED_COMMANDS:
            continue
        steps = [title for title in titles if title.startswith("Step ")]
        assert titles[:len(steps)] == steps, (name, titles)
        numbers = [int(re.match(r"Step (\d+)", title).group(1)) for title in steps]
        assert numbers == sorted(numbers), (name, steps)
        first_titles = [n for title, n in zip(steps, numbers) if " continued: " not in title]
        assert first_titles == list(range(1, len(first_titles) + 1)), (name, steps)
        for title in titles[len(steps):]:
            assert title.startswith(("Optional: ", "Alternative: ")) or title in (
                "Computing resources", "Logging",
            ), (name, title)
        for action in parser._actions:
            if action.required and action.option_strings:
                assert section_of[action.option_strings[0]].startswith("Step "), (name, action.option_strings)


@pytest.mark.parametrize("command", ["ldscore", "munge-sumstats", "annotate", "build-gene-ldscore-index"])
def test_usage_shows_the_typical_command(commands, command):
    parser = commands[command]
    usage = parser.format_usage()
    assert usage.startswith(f"usage: ldsc {command} ")
    assert usage.rstrip().endswith("[options]")
    named = set(re.findall(r"--[A-Za-z][\w-]*", usage))
    assert named <= set(parser._option_string_actions), named - set(parser._option_string_actions)
    required = {action.option_strings[0] for action in parser._actions
                if action.required and action.option_strings}
    assert required <= named, required - named
    if command == "ldscore":
        assert "or: ldsc ldscore --gene-ldscore-index-dir DIR --query-annot-gene-list-sources" in usage


def test_help_cross_references_point_to_earlier_flags(commands):
    for name, parser in commands.items():
        order = [flag for _, flags in _help_sections(parser) for flag in flags]
        for action in parser._actions:
            if not action.option_strings or not action.help:
                continue
            for target in re.findall(r"(?:see|rules as) (--[\w-]+)", action.help):
                assert order.index(target) < order.index(action.option_strings[0]), (
                    name, action.option_strings[0], target,
                )


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
    r2_description = help_text.split("  --r2-matrix-dir DIR", 1)[1].split("\n  --", 1)[0]
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
        ("annotate", "--ref-snp-annot-sources"), ("ldscore", "--baseline-annot-sources"),
        ("build-r2-matrix", "--plink-prefix"), ("build-gene-ldscore-index", "--baseline-annot-sources"),
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


def test_munge_restriction_options_are_mutually_exclusive(commands):
    with pytest.raises(SystemExit) as error:
        commands["munge-sumstats"].parse_args([
            "--raw-sumstats-file", "raw", "--output-dir", "out",
            "--sumstats-snps-file", "keep", "--no-snp-restriction",
        ])
    assert error.value.code == 2
