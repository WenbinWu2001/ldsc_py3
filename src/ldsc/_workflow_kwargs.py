"""Keyword entry points for Python workflows that mirror one command.

Functions such as ``run_h2(**kwargs)`` take the command's option names without
the leading ``--`` and run the same ``*_from_args`` workflow as the CLI. A
keyword call never passes through argparse, so this module repeats the checks
argparse makes on a command line -- known option names, required options,
declared choices and on/off flags -- before any directory, file or log is
created. A misspelled or invalid keyword would otherwise fall back to a
default and silently change the analysis.
"""

from __future__ import annotations

import argparse
from typing import Any, Callable

from ._logging import overwrite_failure_marker
from .config import get_global_config
from .errors import LDSCUsageError

_SHARED_RUNTIME_OPTIONS = {
    "log_level": "set the log level with set_global_config(GlobalConfig(log_level=...))",
    "snp_identifier": "the SNP identifier mode comes from the saved inputs",
    "genome_build": "the genome build comes from the saved inputs",
}


def namespace_from_kwargs(
    function: str,
    add_arguments: Callable[[argparse.ArgumentParser], None],
    kwargs: dict[str, Any],
) -> argparse.Namespace:
    """Return the command namespace for validated Python keyword arguments.

    Parameters
    ----------
    function : str
        Public function name used in error messages, for example ``"run_h2"``.
    add_arguments : callable
        The command's argument registrar, such as ``add_h2_arguments``.
    kwargs : dict
        Option names without the leading ``--``, with Python values.

    Returns
    -------
    argparse.Namespace
        Parser defaults updated with ``kwargs``; ``log_level`` comes from the
        registered ``GlobalConfig``.

    Raises
    ------
    LDSCUsageError
        Listing every problem found: shared runtime options passed as keywords,
        unknown names, missing required options, values outside the declared
        choices, and on/off options given a value other than ``True`` or
        ``False`` (or, for a value-free policy flag, its two policy strings).
    """
    parser = argparse.ArgumentParser(add_help=False)
    add_arguments(parser)
    actions = {action.dest: action for action in parser._actions if action.dest != "log_level"}
    problems = [
        f"{name} is not a keyword argument ({reason})"
        for name, reason in _SHARED_RUNTIME_OPTIONS.items()
        if name in kwargs
    ]
    unknown = sorted(set(kwargs) - set(actions) - set(_SHARED_RUNTIME_OPTIONS))
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
    defaults = {action.dest: action.default for action in parser._actions}
    defaults["log_level"] = get_global_config().log_level
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


def run_from_kwargs(
    function: str,
    add_arguments: Callable[[argparse.ArgumentParser], None],
    run_from_args: Callable[[argparse.Namespace], Any],
    kwargs: dict[str, Any],
) -> Any:
    """Validate ``kwargs``, then run the command workflow under its failure marker.

    The marker names the Python function, so a failed authorized overwrite
    records ``run_h2(...)`` rather than the namespace entry point it calls.
    """
    args = namespace_from_kwargs(function, add_arguments, kwargs)
    with overwrite_failure_marker(args.output_dir, overwrite=args.overwrite, command=f"Python API: {function}(...)"):
        return run_from_args(args)
