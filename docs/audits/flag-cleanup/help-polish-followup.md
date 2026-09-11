# CLI help cleanup implementation

Last updated on: 2026-09-11

Applied the [agreed CLI help guidelines](../../current/cli-help-guidelines.md) across all 13 public commands and 215 visible command-option entries. The [original audit](help-text-audit.md) remains the historical snapshot, including its 50 entries without descriptions.

The subsequent [intercept consolidation](intercept-consolidation.md) removed three occurrences of `--no-intercept`, leaving 212 visible command-option entries. The verification below records the help-polishing phase before that removal.

## Changes

- Every visible flag has a description covering its purpose, relevant input details, effective default or omission behavior, and dependencies. Mode-dependent defaults and ignored versus rejected options are described according to their workflows.
- Help groups prioritize essential inputs and output destinations, then analysis settings and specialized uses, with advanced controls, performance, and logging later. Munging separates sample sizes, filters, coordinate conversion, and column overrides; LD-score help distinguishes direct inputs/settings from indexed gene-list use.
- Required companions and mutually exclusive options are named explicitly, including case/control columns, prevalence pairs, query source alternatives, liftover choices, intercept controls, and LD-window choices.
- Reference-panel SNP counts, interval padding, frequency preservation, fallback sample sizes, and reference-metadata destinations are explained in user terms. Logging descriptions share the four level meanings and `INFO` default without implying that progress is printed to the terminal.
- Usage lines name the actual `ldsc COMMAND`. Argument placeholders identify files, directories, columns, and numerical units. `CLIHelpFormatter` wraps at spaces while preserving flag names and paths.
- Unified parser construction preserves argument groups, descriptions, and placeholders from workflow parsers. Existing parsing rules remain unchanged, including the separate parser construction paths' pre-existing mutual-exclusion handling.

Sources: workflow `build_parser` / `add_*_arguments` functions, [`cli._copy_actions`](../../../src/ldsc/cli.py), [`CLIHelpFormatter`](../../../src/ldsc/_cli_help.py), and [`LOG_LEVEL_HELP`](../../../src/ldsc/_logging.py). Contributor navigation is in [code structure](../../current/code-structure.md).

## Verification

Before edits, help tests failed on missing descriptions, ungrouped ordering, and incomplete dependency text. A separate rendered-help test reproduced hyphen splitting of a companion flag name before the formatter fix. All seven help tests now pass.

An independent before/after comparison covers 21 parser surfaces: 13 unified subparsers and eight workflow parser factories. Option strings, destinations, defaults, required status, argument counts, constants, choices, types, action classes, abbreviation policy, and existing mutual-exclusion groups match exactly. Presentation fields are deliberately excluded from that comparison. Snapshots and command output are retained locally under `/tmp/ldsc-help-polish/` for this run.

Final verification:

- Full test suite: 1,535 passed, one skipped, and 132 subtests passed in 95.20 seconds. The suite reported 189 warnings; the skip applies when the optional `pyarrow` dependency is installed.
- All 26 help invocations passed: both `ldsc COMMAND --help` and `python -m ldsc COMMAND --help` for every public command. The two entry points produced identical text, and all 215 visible flag entries had descriptions.
- All 21 parser contracts matched the saved baseline. Edited Python files compiled successfully, 70 documentation links and anchors passed validation, and `git diff --check` passed.

This cleanup changes help presentation only. Previously approved flag removals, retained spellings, numerical behavior, and file formats are unchanged. Unrelated edits in the shared working tree remain outside this task.
