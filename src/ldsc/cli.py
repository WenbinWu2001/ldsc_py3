"""Unified CLI for the LDSC3 package.

Core functionality:
    Provide one command surface for annotation building, R2 matrix
    generation, LD-score calculation, summary-statistics
    munging, regression, result plotting, and liability-scale post-processing.

Overview
--------
LDSC3 intentionally exposes a single command, ``ldsc``, with
subcommands grouped by user task rather than by historical script name. This
module owns only argument parsing and dispatch. Scientific work remains in the
public workflow modules; the plotting runtime remains behind lazy workflow
imports.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import shlex
import sys
from typing import Sequence

from ._runtime.cli_help import CLIHelpFormatter
from ._annotation import builder as annotation_builder
from ._ldscore.gene_index import build as gene_ldscore_index
from ._ldscore import calculator as ldscore_calculator
from ._refpanel import query as r2_matrix_query, builder as r2_matrix_builder
from ._runtime.workflow_logging import (
    FAILURE_MARKER,
    install_cli_console_handler,
    last_workflow_log_path,
    overwrite_failure_marker,
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

# Single source of truth for the subcommand list shown in `ldsc --help`.
_SUBCOMMAND_HELP = {
    "infer-genome-build": "Inspect an artifact's genome build without writing files.",
    "annotate": "Find SNPs covered by BED intervals or gene boundaries.",
    "ldscore": "Compute LD scores.",
    "build-r2-matrix": "Build pairwise R2 matrices from PLINK genotypes.",
    "build-gene-ldscore-index": "Build a reusable gene LD-score index.",
    "convert-legacy-ldscores": "Convert selected legacy LD-score suites to LDSC3 format.",
    "munge-sumstats": "Munge GWAS summary statistics.",
    "h2": "Estimate heritability from munged summary statistics and LD scores.",
    "partitioned-h2": "Estimate partitioned heritability with independent query models.",
    "quantile-h2": "Summarize fitted heritability across annotation quantiles.",
    "rg": "Estimate genetic correlation.",
    "query-r2-matrix": "Look up R2 for SNP pairs in an R2 matrix.",
    "convert-h2-scale": "Convert a saved observed-scale SNP heritability estimate to liability scale.",
    "plot": "Plot a saved LDSC3 result.",
}


class _NoAbbrevArgumentParser(argparse.ArgumentParser):
    """ArgumentParser variant that disables long-option abbreviation."""

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("allow_abbrev", False)
        kwargs.setdefault("formatter_class", CLIHelpFormatter)
        super().__init__(*args, **kwargs)


def _command_parser(command: str) -> argparse.ArgumentParser:
    """Return the parser that ``ldsc <command>`` parses its arguments with.

    Workflows with their own ``build_parser`` supply it; the others get a
    no-abbreviation parser filled by their ``add_*_arguments`` function.
    Workflow modules are imported only for the requested command.
    """
    if command == "annotate":
        return annotation_builder.build_parser()
    if command == "ldscore":
        return ldscore_calculator.build_parser()
    if command == "build-r2-matrix":
        return r2_matrix_builder.build_parser()
    if command == "build-gene-ldscore-index":
        return gene_ldscore_index.build_parser()
    if command == "convert-legacy-ldscores":
        return _load_legacy_ldscore_converter().build_parser()
    if command == "query-r2-matrix":
        return r2_matrix_query.build_parser()
    if command == "munge-sumstats":
        return _load_sumstats_munger().build_parser()
    parser = _NoAbbrevArgumentParser(prog=f"ldsc {command}")
    if command == "infer-genome-build":
        _add_build_inspection_arguments(parser)
    elif command == "convert-h2-scale":
        _load_h2_scale().add_convert_h2_scale_arguments(parser)
    elif command == "plot":
        _load_plotting().add_plot_arguments(parser)
    elif command == "quantile-h2":
        _load_quantile_h2().add_quantile_h2_arguments(parser)
    else:
        regression_runner = _load_regression_runner()
        {"h2": regression_runner.add_h2_arguments,
         "partitioned-h2": regression_runner.add_partitioned_h2_arguments,
         "rg": regression_runner.add_rg_arguments}[command](parser)
    return parser


def _build_top_level_parser() -> argparse.ArgumentParser:
    """Build a parser that only lists subcommands, for help and usage messages.

    The runtime help, no-argument, and unknown-command paths use this so that
    ``ldsc --help`` does not import the SciPy-heavy regression and munging
    workflows just to print the subcommand list.

    Returns
    -------
    argparse.ArgumentParser
        Parser exposing the subcommand names and one-line help only.
    """
    parser = _NoAbbrevArgumentParser(prog="ldsc", description="LDSC3 command line interface.")
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
        if command == "infer-genome-build":
            from ._genome.build_inspection import inspect_genome_build

            report = inspect_genome_build(_command_parser(command).parse_args(subargv).path)
            print(report)
            return report
        if command == "annotate":
            return annotation_builder.main(subargv)
        if command == "ldscore":
            return ldscore_calculator.main(subargv)
        if command == "build-r2-matrix":
            return r2_matrix_builder.main(subargv)
        if command == "build-gene-ldscore-index":
            return gene_ldscore_index.main(subargv)
        if command == "convert-legacy-ldscores":
            return _load_legacy_ldscore_converter().main(subargv)
        if command == "query-r2-matrix":
            return r2_matrix_query.main(subargv)
        if command == "munge-sumstats":
            return _load_sumstats_munger().main(subargv)
        if command == "convert-h2-scale":
            return _load_h2_scale().run_convert_h2_scale_from_args(_command_parser(command).parse_args(subargv))
        if command == "plot":
            return _load_plotting().run_plot_from_args(_command_parser(command).parse_args(subargv))
        if command == "quantile-h2":
            return _load_quantile_h2().run_quantile_h2_from_args(_command_parser(command).parse_args(subargv))
        if command in {"h2", "partitioned-h2", "rg"}:
            regression_runner = _load_regression_runner()
            run = {"h2": regression_runner.run_h2_from_args,
                   "partitioned-h2": regression_runner.run_partitioned_h2_from_args,
                   "rg": regression_runner.run_rg_from_args}[command]
            return run(_command_parser(command).parse_args(subargv))
    # Only top-level help, no arguments, or an unknown command reach this point;
    # every valid subcommand is dispatched above without importing the heavy
    # workflow modules. The lightweight parser keeps this path fast: it always
    # exits here, printing help or an argparse error.
    _build_top_level_parser().parse_args(argv)


def _add_build_inspection_arguments(parser: argparse.ArgumentParser) -> None:
    parser.description = _SUBCOMMAND_HELP["infer-genome-build"]
    parser.add_argument(
        "path", metavar="PATH",
        help=("Required exact file or recognized LDSC3 artifact directory. Read build declarations "
              "or infer hg19/hg38 from SNP coordinates. No patterns or multiple paths. "
              "Print unknown when evidence is insufficient; unknown exits successfully. "
              "Creates no files or logs."),
    )


def run_cli(argv: Sequence[str] | None = None) -> int:
    """Run the CLI with a clean user-error boundary and traceback logging."""
    cli_mode = argv is None
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    if cli_mode:
        reset_workflow_log_path()
        install_cli_console_handler()
    result = None
    try:
        marker_scope = _cli_failure_marker_scope(raw_argv)
        if marker_scope is None:
            result = main(argv)
        else:
            output_dir, marker_name = marker_scope
            with overwrite_failure_marker(
                output_dir,
                overwrite=True,
                command=shlex.join(["ldsc", *raw_argv]),
                marker_name=marker_name,
            ):
                result = main(argv)
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
        _log_internal_error(exc, workflow_advice=raw_argv[:1] != ["infer-genome-build"])
        return 2
    except Exception as exc:
        _log_internal_error(exc, workflow_advice=raw_argv[:1] != ["infer-genome-build"])
        return 2
    finally:
        close = getattr(result, "close", None)
        if callable(close):
            close()
        if cli_mode:
            remove_cli_console_handler()
    return 0


def _cli_failure_marker_scope(argv: Sequence[str]) -> tuple[Path, str] | None:
    """Return the authorized overwrite marker scope encoded by raw CLI args."""
    if not argv or "--overwrite" not in argv:
        return None
    if "--help" in argv or "-h" in argv:
        return None
    command = argv[0]
    if command not in _SUBCOMMAND_HELP or command == "infer-genome-build":
        return None
    if command == "munge-sumstats" and "--infer-only" in argv:
        return None
    if command == "plot":
        from ._plotting.dispatch import default_plot_dir
        source = _raw_option_value(argv, "--result-dir")
        return (default_plot_dir(source), FAILURE_MARKER) if source is not None else None
    if command == "convert-h2-scale":
        from ._regression.h2_scale import default_liability_scale_dir
        source = _raw_option_value(argv, "--h2-result-dir")
        return (default_liability_scale_dir(source), FAILURE_MARKER) if source is not None else None
    output_dir = _raw_option_value(argv, "--output-dir")
    if output_dir is None:
        return None
    if command == "build-r2-matrix":
        return Path(output_dir), r2_matrix_builder.failure_marker_name(_raw_option_value(argv, "--plink-prefix"))
    return Path(output_dir), FAILURE_MARKER


def _raw_option_value(argv: Sequence[str], option: str) -> str | None:
    """Extract the last option value, matching argparse's store action."""
    value = None
    for index, token in enumerate(argv):
        if token == option and index + 1 < len(argv):
            value = argv[index + 1]
        prefix = option + "="
        if token.startswith(prefix):
            value = token[len(prefix) :]
    return value


def _log_internal_error(exc: BaseException, *, workflow_advice: bool = True) -> None:
    """Emit a concise internal-error line; full traceback is in the run log file."""
    log_path = last_workflow_log_path() if workflow_advice else None
    pointer = f" See {log_path} for the full traceback." if log_path else ""
    hint = " --log-level DEBUG adds developer context if the run needs to be reproduced." if workflow_advice else ""
    LOGGER.error(
        f"Internal error while running ldsc: {exc}.{pointer} "
        f"Report the command and package version with this error and any saved traceback.{hint}"
    )


def _load_regression_runner():
    """Import the regression workflow lazily for SciPy-light entry paths."""
    from ._regression import runner as regression_runner

    return regression_runner


def _load_sumstats_munger():
    """Import the munging workflow lazily for SciPy-light entry paths."""
    from ._sumstats import munger as sumstats_munger

    return sumstats_munger


def _load_quantile_h2():
    """Import the continuous-annotation post-fit workflow lazily."""
    from ._regression import quantile_h2

    return quantile_h2


def _load_legacy_ldscore_converter():
    """Import the legacy converter lazily for lightweight top-level help."""
    from ._ldscore import legacy_converter as legacy_ldscore_converter

    return legacy_ldscore_converter


def _load_h2_scale():
    """Import the post-fit heritability scale conversion lazily."""
    from ._regression import h2_scale

    return h2_scale


def _load_plotting():
    """Import the plotting workflow lazily."""
    from ._plotting import dispatch as plotting

    return plotting
