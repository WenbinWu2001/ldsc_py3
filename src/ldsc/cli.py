"""Unified CLI for the refactored LDSC package.

Core functionality:
    Provide one command surface for annotation building, parquet
    reference-panel generation, LD-score calculation, summary-statistics
    munging, regression, result plotting, and liability-scale post-processing.

Overview
--------
The refactored package intentionally exposes a single command, ``ldsc``, with
subcommands grouped by user task rather than by historical script name. This
module owns only argument parsing and dispatch. Scientific work remains in the
public workflow modules; optional plotting dependencies remain behind lazy
workflow imports.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import re
import shlex
import sys
from typing import Sequence

from . import annotation_builder, gene_ldscore_index, ldscore_calculator, r2_query, ref_panel_builder
from ._logging import (
    install_cli_console_handler,
    last_workflow_log_path,
    remove_cli_console_handler,
    reset_workflow_log_path,
    overwrite_failure_marker,
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
    "build-gene-ldscore-index": "Build an exact disjoint-atom gene LD-score index.",
    "convert-ldsc2-ldscores": "Convert selected LDSC2 LD-score suites to LDSC3 format.",
    "munge-sumstats": "Munge GWAS summary statistics.",
    "h2": "Estimate heritability from munged sumstats and LD scores.",
    "partitioned-h2": "Estimate partitioned heritability by looping over query annotations.",
    "quantile-h2": "Project a fitted partitioned-LDSC model onto annotation quantiles.",
    "rg": "Estimate genetic correlation.",
    "query-r2": "Query R2 for SNP pairs from a reference panel.",
    "convert-h2-scale": "Convert a saved observed-scale h2 estimate to liability scale.",
    "plot": "Create the approved plot for a canonical LDSC result directory.",
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

    gene_index_parser = subparsers.add_parser(
        "build-gene-ldscore-index", help=_SUBCOMMAND_HELP["build-gene-ldscore-index"]
    )
    _copy_actions(gene_index_parser, gene_ldscore_index.build_parser())

    converter_parser = subparsers.add_parser(
        "convert-ldsc2-ldscores", help=_SUBCOMMAND_HELP["convert-ldsc2-ldscores"]
    )
    legacy_ldscore_converter = _load_legacy_ldscore_converter()
    _copy_actions(converter_parser, legacy_ldscore_converter.build_parser())

    sumstats_munger = _load_sumstats_munger()
    munge_parser = subparsers.add_parser("munge-sumstats", help=_SUBCOMMAND_HELP["munge-sumstats"])
    _copy_actions(munge_parser, sumstats_munger.kernel_parser())

    regression_runner = _load_regression_runner()
    h2_parser = subparsers.add_parser("h2", help=_SUBCOMMAND_HELP["h2"])
    regression_runner.add_h2_arguments(h2_parser)

    partitioned_parser = subparsers.add_parser("partitioned-h2", help=_SUBCOMMAND_HELP["partitioned-h2"])
    regression_runner.add_partitioned_h2_arguments(partitioned_parser)

    quantile_h2 = _load_quantile_h2()
    quantile_parser = subparsers.add_parser("quantile-h2", help=_SUBCOMMAND_HELP["quantile-h2"])
    quantile_h2.add_quantile_h2_arguments(quantile_parser)

    rg_parser = subparsers.add_parser("rg", help=_SUBCOMMAND_HELP["rg"])
    regression_runner.add_rg_arguments(rg_parser)

    query_r2_parser = subparsers.add_parser("query-r2", help=_SUBCOMMAND_HELP["query-r2"])
    _copy_actions(query_r2_parser, r2_query.build_parser())

    h2_scale = _load_h2_scale()
    conversion_parser = subparsers.add_parser(
        "convert-h2-scale", help=_SUBCOMMAND_HELP["convert-h2-scale"]
    )
    h2_scale.add_convert_h2_scale_arguments(conversion_parser)

    plotting = _load_plotting()
    plot_parser = subparsers.add_parser("plot", help=_SUBCOMMAND_HELP["plot"])
    plotting.add_plot_arguments(plot_parser)
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
        if command == "build-gene-ldscore-index":
            return gene_ldscore_index.main(subargv)
        if command == "convert-ldsc2-ldscores":
            legacy_ldscore_converter = _load_legacy_ldscore_converter()
            return legacy_ldscore_converter.main(subargv)
        if command == "query-r2":
            return r2_query.main(subargv)
        if command == "convert-h2-scale":
            h2_scale = _load_h2_scale()
            parser = _NoAbbrevArgumentParser(
                prog="ldsc convert-h2-scale",
                description="Convert a saved observed-scale h2 estimate to liability scale.",
            )
            h2_scale.add_convert_h2_scale_arguments(parser)
            return h2_scale.run_convert_h2_scale_from_args(parser.parse_args(subargv))
        if command == "plot":
            plotting = _load_plotting()
            parser = _NoAbbrevArgumentParser(
                prog="ldsc plot",
                description="Create the approved plot for a canonical LDSC result directory.",
            )
            plotting.add_plot_arguments(parser)
            return plotting.run_plot_from_args(parser.parse_args(subargv))
        if command == "munge-sumstats":
            sumstats_munger = _load_sumstats_munger()
            return sumstats_munger.main(subargv)
        if command == "h2":
            regression_runner = _load_regression_runner()
            parser = _NoAbbrevArgumentParser(prog="ldsc h2", description="Estimate heritability from munged sumstats and LD scores.")
            regression_runner.add_h2_arguments(parser)
            parsed = parser.parse_args(subargv)
            return regression_runner.run_h2_from_args(parsed)
        if command == "partitioned-h2":
            regression_runner = _load_regression_runner()
            parser = _NoAbbrevArgumentParser(
                prog="ldsc partitioned-h2",
                description="Estimate partitioned heritability by looping over query annotations.",
            )
            regression_runner.add_partitioned_h2_arguments(parser)
            parsed = parser.parse_args(subargv)
            return regression_runner.run_partitioned_h2_from_args(parsed)
        if command == "quantile-h2":
            quantile_h2 = _load_quantile_h2()
            parser = _NoAbbrevArgumentParser(
                prog="ldsc quantile-h2",
                description="Project a fitted partitioned-LDSC model onto continuous-annotation quantiles.",
            )
            quantile_h2.add_quantile_h2_arguments(parser)
            return quantile_h2.run_quantile_h2_from_args(parser.parse_args(subargv))
        if command == "rg":
            regression_runner = _load_regression_runner()
            parser = _NoAbbrevArgumentParser(prog="ldsc rg", description="Estimate genetic correlation.")
            regression_runner.add_rg_arguments(parser)
            parsed = parser.parse_args(subargv)
            return regression_runner.run_rg_from_args(parsed)
    # Only top-level help, no arguments, or an unknown command reach this point;
    # every valid subcommand is dispatched above without importing the heavy
    # workflow modules. The lightweight parser keeps this path fast: it always
    # exits here, printing help or an argparse error.
    _build_top_level_parser().parse_args(argv)


def run_cli(argv: Sequence[str] | None = None) -> int:
    """Run the CLI with a clean user-error boundary and traceback logging."""
    cli_mode = argv is None
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    if cli_mode:
        reset_workflow_log_path()
        install_cli_console_handler()
    try:
        marker_scope = _cli_failure_marker_scope(raw_argv)
        if marker_scope is None:
            main(argv)
        else:
            output_dir, marker_name = marker_scope
            with overwrite_failure_marker(
                output_dir,
                overwrite=True,
                command=shlex.join(["ldsc", *raw_argv]),
                marker_name=marker_name,
            ):
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


def _cli_failure_marker_scope(argv: Sequence[str]) -> tuple[Path, str] | None:
    """Return the authorized overwrite marker scope encoded by raw CLI args."""
    if not argv or "--overwrite" not in argv:
        return None
    command = argv[0]
    if command not in _SUBCOMMAND_HELP:
        return None
    if command == "munge-sumstats" and "--infer-only" in argv:
        return None
    if command == "plot":
        source = _raw_option_value(argv, "--result-dir")
        return (Path(source) / "plots", "RUN_FAILED.txt") if source is not None else None
    if command == "convert-h2-scale":
        source = _raw_option_value(argv, "--h2-result-dir")
        if source is None:
            return None
        return Path(source) / "postprocessing" / "liability-scale", "RUN_FAILED.txt"
    output_dir = _raw_option_value(argv, "--output-dir")
    if output_dir is None:
        return None
    if command == "build-ref-panel":
        prefix = _raw_option_value(argv, "--plink-prefix")
        chromosome = _chromosome_from_concrete_prefix(prefix)
        if chromosome is not None:
            return Path(output_dir), f"RUN_FAILED.chr{chromosome}.txt"
    return Path(output_dir), "RUN_FAILED.txt"


def _raw_option_value(argv: Sequence[str], option: str) -> str | None:
    """Extract one raw ``--option value`` or ``--option=value`` token."""
    for index, token in enumerate(argv):
        if token == option and index + 1 < len(argv):
            return argv[index + 1]
        prefix = option + "="
        if token.startswith(prefix):
            return token[len(prefix) :]
    return None


def _chromosome_from_concrete_prefix(prefix: str | None) -> str | None:
    """Infer an explicit chromosome token without resolving or reading inputs."""
    if not prefix or "@" in prefix:
        return None
    matches = re.findall(
        r"(?:^|[._-])(?:chr)?(1[0-9]|2[0-2]|[1-9]|X|Y)(?=$|[._-])",
        Path(prefix).name,
        flags=re.IGNORECASE,
    )
    normalized = {match.upper() for match in matches}
    return next(iter(normalized)) if len(normalized) == 1 else None


def _log_internal_error(exc: BaseException) -> None:
    """Emit a concise internal-error line; full traceback is in the run log file."""
    log_path = last_workflow_log_path()
    pointer = f" See {log_path} for the full traceback." if log_path else ""
    LOGGER.error(f"Internal error while running ldsc: {exc}.{pointer}")


def _run_annotate(args: argparse.Namespace):
    """Dispatch the ``annotate`` subcommand."""
    return annotation_builder.run_annotate_from_args(args)


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


def _load_quantile_h2():
    """Import the continuous-annotation post-fit workflow lazily."""
    from . import quantile_h2

    return quantile_h2


def _load_legacy_ldscore_converter():
    """Import the legacy converter lazily for lightweight top-level help."""
    from . import legacy_ldscore_converter

    return legacy_ldscore_converter


def _load_h2_scale():
    """Import post-fit h2 conversion lazily."""
    from . import h2_scale

    return h2_scale


def _load_plotting():
    """Import the optional plotting workflow lazily."""
    from . import plotting

    return plotting
