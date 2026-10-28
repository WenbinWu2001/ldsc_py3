"""Keyword entry points for Python workflows that mirror one command.

Functions such as ``run_h2(**kwargs)`` and ``run_ldscore(**kwargs)`` take the
command's option names without the leading ``--`` and run the same
``*_from_args`` workflow as the CLI. A keyword call never passes through
argparse, so this module repeats the checks argparse makes on a command line --
known option names, required options, declared choices and on/off flags --
before any directory, file or log is created. A misspelled or invalid keyword
would otherwise fall back to a default and silently change the analysis.
"""

from __future__ import annotations

import argparse
from typing import Any, Callable, Mapping

from .workflow_logging import overwrite_failure_marker
from ..errors import LDSCUsageError

_LOG_LEVEL_REASON = "set the log level with set_global_config(GlobalConfig(log_level=...))"
# Regression and post-fit workflows read SNP identity and build from saved inputs.
SAVED_INPUT_RUNTIME_OPTIONS = {
    "log_level": _LOG_LEVEL_REASON,
    "snp_identifier": "the SNP identifier mode comes from the saved inputs",
    "genome_build": "the genome build comes from the saved inputs",
}
# Workflows that read identity and build from the registered GlobalConfig.
GLOBAL_CONFIG_RUNTIME_OPTIONS = {
    "log_level": _LOG_LEVEL_REASON,
    "snp_identifier": "set the SNP identifier mode with set_global_config(GlobalConfig(snp_identifier=...))",
    "genome_build": "set the genome build with set_global_config(GlobalConfig(genome_build=...))",
}


def command_parser(add_arguments: Callable[[argparse.ArgumentParser], None]) -> argparse.ArgumentParser:
    """Return a parser holding one command's options from its argument registrar."""
    parser = argparse.ArgumentParser(add_help=False)
    add_arguments(parser)
    return parser


def namespace_from_kwargs(
    function: str,
    parser: argparse.ArgumentParser,
    kwargs: dict[str, Any],
    *,
    rejected: Mapping[str, str],
    log_level: str,
) -> argparse.Namespace:
    """Return the command namespace for validated Python keyword arguments.

    Parameters
    ----------
    function : str
        Public function name used in error messages, for example ``"run_h2"``.
    parser : argparse.ArgumentParser
        The command's parser; its options define the valid keyword names.
    kwargs : dict
        Option names without the leading ``--``, with Python values.
    rejected : mapping of str to str
        Option names that are not keyword arguments of this function, each
        with the reason shown in the error, such as where the value comes from.
    log_level : str
        The log level of the registered ``GlobalConfig``, supplied by the
        calling workflow.

    Returns
    -------
    argparse.Namespace
        Parser defaults updated with ``kwargs``, with ``log_level`` set to the
        ``log_level`` argument.

    Raises
    ------
    LDSCUsageError
        Listing every problem found: rejected names, unknown names, missing
        required options, values outside the declared choices, and on/off
        options given a value other than ``True`` or ``False`` (or, for a
        value-free policy flag, its two policy strings).
    """
    actions = {
        action.dest: action for action in parser._actions
        if action.dest not in rejected and action.dest not in {"help", "log_level"}
    }
    problems = [f"{name} is not a keyword argument ({reason})" for name, reason in rejected.items() if name in kwargs]
    unknown = sorted(set(kwargs) - set(actions) - set(rejected))
    if unknown:
        problems.append(
            f"unknown keyword argument(s) {', '.join(unknown)}; valid names are {', '.join(sorted(actions))}"
        )
    missing = [dest for dest, action in actions.items() if action.required and not kwargs.get(dest)]
    if missing:
        problems.append(f"missing required keyword argument(s) {', '.join(missing)}")
    for name, value in kwargs.items():
        action = actions.get(name)
        if action is None or (value is None and action.default is None):
            continue
        if action.choices is not None and value not in action.choices:
            problems.append(f"{name}={value!r} is not one of {', '.join(map(repr, action.choices))}")
        elif action.nargs == 0 and value not in (action.default, action.const):
            problems.append(f"{name}={value!r} must be {action.default!r} or {action.const!r}")
    if problems:
        raise LDSCUsageError(f"{function}() cannot run: " + "; ".join(problems) + ".")
    defaults = {action.dest: action.default for action in parser._actions if action.dest != "help"}
    defaults["log_level"] = log_level
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


def run_from_kwargs(
    function: str,
    parser: argparse.ArgumentParser,
    run_from_args: Callable[[argparse.Namespace], Any],
    kwargs: dict[str, Any],
    *,
    log_level: str,
) -> Any:
    """Validate ``kwargs``, then run the command workflow under its failure marker.

    Used by the regression and post-fit workflows, whose SNP identity and build
    come from saved inputs. The marker names the Python function, so a failed
    authorized overwrite records ``run_h2(...)`` rather than the namespace
    entry point it calls.
    """
    args = namespace_from_kwargs(function, parser, kwargs, rejected=SAVED_INPUT_RUNTIME_OPTIONS, log_level=log_level)
    with overwrite_failure_marker(args.output_dir, overwrite=args.overwrite, command=f"Python API: {function}(...)"):
        return run_from_args(args)
