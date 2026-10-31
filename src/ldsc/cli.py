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
from dataclasses import dataclass
from importlib import import_module
import logging
from pathlib import Path
import shlex
import sys
from typing import Sequence

from ._runtime.cli_help import CLIHelpFormatter
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

@dataclass(frozen=True)
class _Command:
    """One ``ldsc`` subcommand.

    ``module`` is a module path relative to this package, imported only when
    the command's parser is built or the command runs, so ``ldsc --help`` does
    not import the SciPy-heavy workflows. A module that parses its own
    arguments provides ``build_parser()`` and ``main(argv)``; otherwise
    ``add_arguments`` fills a no-abbreviation parser and ``run`` takes the
    parsed namespace. ``read_only`` commands write nothing, so they stay
    outside the failure-marker lifecycle and their errors give no
    workflow-log advice.
    """

    help: str
    module: str
    add_arguments: str | None = None
    run: str | None = None
    read_only: bool = False


# The subcommands, in `ldsc --help` order.
_COMMANDS = {
    "infer-genome-build": _Command(
        "Inspect an artifact's genome build without writing files.",
        "._genome.build_inspection", "add_inspection_arguments", "run_inspection_from_args", read_only=True,
    ),
    "annotate": _Command("Find SNPs covered by BED intervals or gene boundaries.", "._annotation.builder"),
    "ldscore": _Command("Compute LD scores.", "._ldscore.calculator"),
    "build-r2-matrix": _Command("Build pairwise R2 matrices from PLINK genotypes.", "._refpanel.builder"),
    "build-gene-ldscore-index": _Command("Build a reusable gene LD-score index.", "._ldscore.gene_index.build"),
    "convert-legacy-ldscores": _Command(
        "Convert selected legacy LD-score suites to LDSC3 format.", "._ldscore.legacy_converter",
    ),
    "munge-sumstats": _Command("Munge GWAS summary statistics.", "._sumstats.munger"),
    "h2": _Command(
        "Estimate heritability from munged summary statistics and LD scores.",
        "._regression.runner", "add_h2_arguments", "run_h2_from_args",
    ),
    "partitioned-h2": _Command(
        "Estimate partitioned heritability with independent query models.",
        "._regression.runner", "add_partitioned_h2_arguments", "run_partitioned_h2_from_args",
    ),
    "quantile-h2": _Command(
        "Summarize fitted heritability across annotation quantiles.",
        "._regression.quantile_h2", "add_quantile_h2_arguments", "run_quantile_h2_from_args",
    ),
    "rg": _Command("Estimate genetic correlation.", "._regression.runner", "add_rg_arguments", "run_rg_from_args"),
    "query-r2-matrix": _Command("Look up R2 for SNP pairs in an R2 matrix.", "._refpanel.query"),
    "convert-h2-scale": _Command(
        "Convert a saved observed-scale SNP heritability estimate to liability scale.",
        "._regression.h2_scale", "add_convert_h2_scale_arguments", "run_convert_h2_scale_from_args",
    ),
    "plot": _Command("Plot a saved LDSC3 result.", "._plotting.dispatch", "add_plot_arguments", "run_plot_from_args"),
}


class _NoAbbrevArgumentParser(argparse.ArgumentParser):
    """ArgumentParser variant that disables long-option abbreviation."""

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("allow_abbrev", False)
        kwargs.setdefault("formatter_class", CLIHelpFormatter)
        super().__init__(*args, **kwargs)


def _command_parser(name: str) -> argparse.ArgumentParser:
    """Return the parser that ``ldsc <name>`` parses its arguments with.

    Only the command's own workflow module is imported.
    """
    command = _COMMANDS[name]
    module = import_module(command.module, __package__)
    if command.add_arguments is None:
        return module.build_parser()
    parser = _NoAbbrevArgumentParser(prog=f"ldsc {name}")
    getattr(module, command.add_arguments)(parser)
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
    for name, command in _COMMANDS.items():
        subparsers.add_parser(name, help=command.help)
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
    command = _COMMANDS.get(argv[0]) if argv else None
    if command is None:
        # Top-level help, no arguments, or an unknown command: the lightweight
        # parser prints help or an argparse error and exits without importing
        # any workflow module.
        _build_top_level_parser().parse_args(argv)
        return None
    module = import_module(command.module, __package__)
    if command.add_arguments is None:
        return module.main(argv[1:])
    return getattr(module, command.run)(_command_parser(argv[0]).parse_args(argv[1:]))


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
        _log_internal_error(exc, workflow_advice=not _read_only(raw_argv))
        return 2
    except Exception as exc:
        _log_internal_error(exc, workflow_advice=not _read_only(raw_argv))
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
    name = argv[0]
    if name not in _COMMANDS or _read_only(argv):
        return None
    if name == "munge-sumstats" and "--infer-only" in argv:
        return None
    if name == "plot":
        from ._plotting.dispatch import default_plot_dir
        source = _raw_option_value(argv, "--result-dir")
        return (default_plot_dir(source), FAILURE_MARKER) if source is not None else None
    if name == "convert-h2-scale":
        from ._regression.h2_scale import default_liability_scale_dir
        source = _raw_option_value(argv, "--h2-result-dir")
        return (default_liability_scale_dir(source), FAILURE_MARKER) if source is not None else None
    output_dir = _raw_option_value(argv, "--output-dir")
    if output_dir is None:
        return None
    if name == "build-r2-matrix":
        from ._refpanel.builder import failure_marker_name
        return Path(output_dir), failure_marker_name(_raw_option_value(argv, "--plink-prefix"))
    return Path(output_dir), FAILURE_MARKER


def _read_only(argv: Sequence[str]) -> bool:
    """Return whether ``argv`` names a subcommand that writes nothing."""
    command = _COMMANDS.get(argv[0]) if argv else None
    return command is not None and command.read_only


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
