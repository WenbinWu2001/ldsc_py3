"""Unified CLI for the refactored LDSC package.

Core functionality:
    Provide one command surface for annotation building, parquet
    reference-panel generation, LD-score calculation, summary-statistics
    munging, and regression workflows.

Overview
--------
The refactored package intentionally exposes a single command, ``ldsc``, with
subcommands grouped by user task rather than by historical script name. This
module owns only argument parsing and dispatch. Scientific work remains in the
public workflow modules, including the standalone parquet reference-panel
builder added during the restructuring pass.
"""

from __future__ import annotations

import argparse
import logging
import sys
from typing import Sequence

from . import annotation_builder, ldscore_calculator, r2_query, ref_panel_builder
from ._logging import (
    install_cli_console_handler,
    last_workflow_log_path,
    remove_cli_console_handler,
    reset_workflow_log_path,
)
from .errors import LDSCError, LDSCUserError

LOGGER = logging.getLogger("LDSC.cli")
_USER_ERROR_TYPES = (
    LDSCUserError,
    ValueError,
    FileNotFoundError,
    FileExistsError,
    NotADirectoryError,
    ImportError,
)

# Single source of truth for the subcommand list shown in `ldsc --help`. Used by
# both the full parser and the lightweight help parser so the two cannot drift.
_SUBCOMMAND_HELP = {
    "annotate": "Project BED files to SNP-level query annotations.",
    "ldscore": "Compute LD scores.",
    "build-ref-panel": "Build standard parquet reference panels.",
    "munge-sumstats": "Munge GWAS summary statistics.",
    "h2": "Estimate heritability from munged sumstats and LD scores.",
    "partitioned-h2": "Estimate partitioned heritability by looping over query annotations.",
    "rg": "Estimate genetic correlation.",
    "query-r2": "Query R2 for SNP pairs from a reference panel.",
}


class _NoAbbrevArgumentParser(argparse.ArgumentParser):
    """ArgumentParser variant that disables long-option abbreviation."""

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("allow_abbrev", False)
        super().__init__(*args, **kwargs)


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level ``ldsc`` command parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser whose subcommands mirror the supported public workflows.
    """
    parser = _NoAbbrevArgumentParser(prog="ldsc", description="LDSC command line interface.")
    subparsers = parser.add_subparsers(dest="command", required=True, parser_class=_NoAbbrevArgumentParser)

    annotate_parser = subparsers.add_parser("annotate", help=_SUBCOMMAND_HELP["annotate"])
    annotation_builder.add_annotate_arguments(annotate_parser)

    ldscore_parser = subparsers.add_parser("ldscore", help=_SUBCOMMAND_HELP["ldscore"])
    _copy_actions(ldscore_parser, ldscore_calculator.build_parser())

    ref_panel_parser = subparsers.add_parser("build-ref-panel", help=_SUBCOMMAND_HELP["build-ref-panel"])
    _copy_actions(ref_panel_parser, ref_panel_builder.build_parser())

    sumstats_munger = _load_sumstats_munger()
    munge_parser = subparsers.add_parser("munge-sumstats", help=_SUBCOMMAND_HELP["munge-sumstats"])
    _copy_actions(munge_parser, sumstats_munger.kernel_parser())

    regression_runner = _load_regression_runner()
    h2_parser = subparsers.add_parser("h2", help=_SUBCOMMAND_HELP["h2"])
    regression_runner.add_h2_arguments(h2_parser)

    partitioned_parser = subparsers.add_parser("partitioned-h2", help=_SUBCOMMAND_HELP["partitioned-h2"])
    regression_runner.add_partitioned_h2_arguments(partitioned_parser)

    rg_parser = subparsers.add_parser("rg", help=_SUBCOMMAND_HELP["rg"])
    regression_runner.add_rg_arguments(rg_parser)

    query_r2_parser = subparsers.add_parser("query-r2", help=_SUBCOMMAND_HELP["query-r2"])
    _copy_actions(query_r2_parser, r2_query.build_parser())
    return parser


def _build_top_level_parser() -> argparse.ArgumentParser:
    """Build a parser that only lists subcommands, for help and usage messages.

    The runtime help, no-argument, and unknown-command paths use this instead of
    :func:`build_parser` so that ``ldsc --help`` does not import the SciPy-heavy
    regression and munging workflows just to print the subcommand list.

    Returns
    -------
    argparse.ArgumentParser
        Parser exposing the subcommand names and one-line help only.
    """
    parser = _NoAbbrevArgumentParser(prog="ldsc", description="LDSC command line interface.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    for name, help_text in _SUBCOMMAND_HELP.items():
        subparsers.add_parser(name, help=help_text)
    return parser


def main(argv: Sequence[str] | None = None):
    """Run the unified command-line interface.

    Parameters
    ----------
    argv : sequence of str or None, optional
        Explicit argument vector. Default is ``None``, which reads arguments
        from ``sys.argv``.

    Returns
    -------
    object
        Workflow-specific result object returned by the dispatched subcommand.
    """
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv:
        command = argv[0]
        subargv = argv[1:]
        if command == "annotate":
            return annotation_builder.main(subargv)
        if command == "ldscore":
            return ldscore_calculator.main(subargv)
        if command == "build-ref-panel":
            return ref_panel_builder.main(subargv)
        if command == "query-r2":
            return r2_query.main(subargv)
        if command == "munge-sumstats":
            sumstats_munger = _load_sumstats_munger()
            return sumstats_munger.main(subargv)
        if command == "h2":
            regression_runner = _load_regression_runner()
            parser = _NoAbbrevArgumentParser(prog="ldsc h2", description="Estimate heritability from munged sumstats and LD scores.")
            regression_runner.add_h2_arguments(parser)
            parsed = parser.parse_args(subargv)
            result = regression_runner.run_h2_from_args(parsed)
            _print_table_stdout_if_needed(parsed, result)
            return result
        if command == "partitioned-h2":
            regression_runner = _load_regression_runner()
            parser = _NoAbbrevArgumentParser(
                prog="ldsc partitioned-h2",
                description="Estimate partitioned heritability by looping over query annotations.",
            )
            regression_runner.add_partitioned_h2_arguments(parser)
            parsed = parser.parse_args(subargv)
            result = regression_runner.run_partitioned_h2_from_args(parsed)
            _print_table_stdout_if_needed(parsed, result)
            return result
        if command == "rg":
            regression_runner = _load_regression_runner()
            parser = _NoAbbrevArgumentParser(prog="ldsc rg", description="Estimate genetic correlation.")
            regression_runner.add_rg_arguments(parser)
            parsed = parser.parse_args(subargv)
            result = regression_runner.run_rg_from_args(parsed)
            _print_rg_stdout_if_needed(parsed, result)
            return result
    # Only top-level help, no arguments, or an unknown command reach this point;
    # every valid subcommand is dispatched above without importing the heavy
    # workflow modules. The lightweight parser keeps this path fast: it always
    # exits here, printing help or an argparse error.
    _build_top_level_parser().parse_args(argv)


def run_cli(argv: Sequence[str] | None = None) -> int:
    """Run the CLI with a clean user-error boundary and traceback logging."""
    cli_mode = argv is None
    if cli_mode:
        reset_workflow_log_path()
        install_cli_console_handler()
    try:
        main(argv)
    except SystemExit as exc:
        if exc.code is None:
            return 0
        if isinstance(exc.code, int):
            return exc.code
        LOGGER.error(f"Error: {exc.code}")
        return 1
    except _USER_ERROR_TYPES as exc:
        LOGGER.error(f"Error: {exc}")
        return 1
    except LDSCError as exc:
        _log_internal_error(exc)
        return 2
    except Exception as exc:
        _log_internal_error(exc)
        return 2
    finally:
        if cli_mode:
            remove_cli_console_handler()
    return 0


def _log_internal_error(exc: BaseException) -> None:
    """Emit a concise internal-error line; full traceback is in the run log file."""
    log_path = last_workflow_log_path()
    pointer = f" See {log_path} for the full traceback." if log_path else ""
    LOGGER.error(f"Internal error while running ldsc: {exc}.{pointer}")


def _run_annotate(args: argparse.Namespace):
    """Dispatch the ``annotate`` subcommand."""
    return annotation_builder.run_annotate_from_args(args)


def _print_rg_stdout_if_needed(args: argparse.Namespace, result) -> None:
    """Print concise rg output for no-output-dir CLI runs."""
    if getattr(args, "output_dir", None):
        return
    rg = getattr(result, "rg", None)
    if rg is not None:
        _print_table_stdout_if_needed(args, rg)


def _print_table_stdout_if_needed(args: argparse.Namespace, table) -> None:
    """Print a compact TSV table for no-output-dir regression CLI runs."""
    if getattr(args, "output_dir", None):
        return
    if hasattr(table, "to_csv"):
        table.to_csv(sys.stdout, sep="\t", index=False, na_rep="NaN")


def _copy_actions(target: argparse.ArgumentParser, source: argparse.ArgumentParser) -> None:
    """Clone option actions from ``source`` onto ``target``.

    This keeps the unified CLI aligned with the feature parsers defined in the
    workflow modules without duplicating every flag definition in two places.
    """
    for action in source._actions:
        if action.dest == "help":
            continue
        option_strings = list(action.option_strings)
        kwargs = {
            "dest": action.dest,
            "default": action.default,
            "required": action.required,
            "help": action.help,
        }
        if getattr(action, "choices", None) is not None:
            kwargs["choices"] = action.choices
        if getattr(action, "type", None) is not None:
            kwargs["type"] = action.type
        if getattr(action, "nargs", None) is not None:
            kwargs["nargs"] = action.nargs
        if action.const is not None:
            kwargs["const"] = action.const
        if action.__class__.__name__ == "_StoreTrueAction":
            target.add_argument(*option_strings, action="store_true", default=action.default, help=action.help)
        elif action.__class__.__name__ == "_StoreFalseAction":
            target.add_argument(*option_strings, action="store_false", default=action.default, help=action.help)
        else:
            target.add_argument(*option_strings, **kwargs)


def _load_regression_runner():
    """Import the regression workflow lazily for SciPy-light entry paths."""
    from . import regression_runner

    return regression_runner


def _load_sumstats_munger():
    """Import the munging workflow lazily for SciPy-light entry paths."""
    from . import sumstats_munger

    return sumstats_munger
