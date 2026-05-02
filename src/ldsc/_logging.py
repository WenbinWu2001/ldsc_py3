"""Shared per-run logging helpers for public LDSC workflows.

Overview
--------
This module centralizes the small amount of logging orchestration that belongs
at workflow boundaries: optional per-run file handlers, stable lifecycle audit
lines, and a shared LDSC package logger threshold. It intentionally avoids
decorating computational class methods so direct APIs stay focused on data
execution unless a public wrapper supplies a workflow log path.

Design Notes
------------
- File handlers attach to the ``LDSC`` ancestor logger so workflow and kernel
  child records are captured together.
- Scientific output paths and log paths should be preflighted by callers before
  entering :func:`workflow_logging`.
- Lifecycle audit lines are file-only; user-facing error reporting remains the
  responsibility of the CLI boundary.
"""

from __future__ import annotations

import logging
from pathlib import Path
import shlex
import sys
import threading
import time
from datetime import datetime
from os import PathLike
from types import TracebackType
from typing import Any


_LDSC_LOGGER_NAME = "LDSC"
_CURRENT: threading.local = threading.local()
_BORDER = "=" * 51
_FLOAT_FORMAT = ".4g"


def format_float(value: float) -> str:
    """Return the compact float representation used for log text only."""
    return format(value, _FLOAT_FORMAT)


def configure_package_logging(level: str = "INFO") -> None:
    """Set the ``LDSC`` logger threshold for console-style workflow logging."""
    logging.getLogger(_LDSC_LOGGER_NAME).setLevel(_level_number(level))


def workflow_logging(
    workflow_name: str,
    log_path: str | PathLike[str] | None,
    *,
    log_level: str = "INFO",
) -> "_WorkflowLoggingContext":
    """Create a workflow-scoped logging context.

    Parameters
    ----------
    workflow_name : str
        Human-readable workflow label written in lifecycle audit lines.
    log_path : str, path-like, or None
        Destination for the per-run log file. When ``None``, the context only
        applies the package logger threshold and writes no file.
    log_level : {"INFO", "DEBUG", "WARNING", "ERROR"}, optional
        Minimum level for module log records captured by the temporary file
        handler and emitted by existing console handlers. Default is ``"INFO"``.

    Returns
    -------
    context : _WorkflowLoggingContext
        Context manager that restores the previous ``LDSC`` logger level and
        removes its file handler on exit.

    Notes
    -----
    Lifecycle audit lines are written directly to the file, so ``Started``,
    ``Finished``/``Failed``, and elapsed-time lines are present even when
    ``log_level`` suppresses ordinary INFO records. Escaping exceptions are not
    logged here; callers and the CLI boundary decide how to present errors.
    """
    return _WorkflowLoggingContext(workflow_name, log_path, log_level=log_level)


def log_inputs(**items: Any) -> None:
    """Write an ``Inputs:`` audit section to the active workflow log, if any."""
    context = _active_context()
    if context is not None:
        context.log_inputs(**items)


def log_outputs(**items: Any) -> None:
    """Write an ``Outputs:`` audit section to the active workflow log, if any."""
    context = _active_context()
    if context is not None:
        context.log_outputs(**items)


def _active_context() -> "_WorkflowLoggingContext | None":
    return getattr(_CURRENT, "context", None)


def _level_number(level: str) -> int:
    normalized = str(level).upper()
    value = logging.getLevelName(normalized)
    if not isinstance(value, int):
        raise ValueError(f"log_level must be one of DEBUG, INFO, WARNING, ERROR; got {level!r}.")
    if normalized not in {"DEBUG", "INFO", "WARNING", "ERROR"}:
        raise ValueError(f"log_level must be one of DEBUG, INFO, WARNING, ERROR; got {level!r}.")
    return value


class _WorkflowLoggingContext:
    """Concrete context manager returned by :func:`workflow_logging`."""

    def __init__(
        self,
        workflow_name: str,
        log_path: str | PathLike[str] | None,
        *,
        log_level: str,
    ) -> None:
        self._workflow_name = workflow_name
        self._log_path = Path(log_path) if log_path is not None else None
        self._level = _level_number(log_level)
        self._handler: logging.FileHandler | None = None
        self._saved_level: int | None = None
        self._start_time: float | None = None
        self._nested = False

    def __enter__(self) -> "_WorkflowLoggingContext":
        if _active_context() is not None:
            self._nested = True
            return self

        logger = logging.getLogger(_LDSC_LOGGER_NAME)
        self._saved_level = logger.level
        logger.setLevel(self._level)
        self._start_time = time.monotonic()

        try:
            if self._log_path is not None:
                self._log_path.parent.mkdir(parents=True, exist_ok=True)
                self._handler = logging.FileHandler(self._log_path, mode="w", encoding="utf-8")
                self._handler.setLevel(self._level)
                self._handler.setFormatter(logging.Formatter("%(message)s"))
                logger.addHandler(self._handler)

            _CURRENT.context = self
            self._write_header()
        except Exception:
            if self._handler is not None:
                logger.removeHandler(self._handler)
                self._handler.close()
                self._handler = None
            logger.setLevel(self._saved_level)
            if _active_context() is self:
                _CURRENT.context = None
            raise
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> bool:
        if self._nested:
            return False

        status = "Failed" if exc_type is not None else "Finished"
        self._write_footer(status)

        logger = logging.getLogger(_LDSC_LOGGER_NAME)
        if self._handler is not None:
            logger.removeHandler(self._handler)
            self._handler.close()
        if self._saved_level is not None:
            logger.setLevel(self._saved_level)
        if _active_context() is self:
            _CURRENT.context = None
        return False

    def log_inputs(self, **items: Any) -> None:
        self._write_items("Inputs", items)

    def log_outputs(self, **items: Any) -> None:
        self._write_items("Outputs", items)

    def _write_header(self) -> None:
        start_dt = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self._write_line(_BORDER)
        self._write_line(f"LDSC {self._workflow_name} Started {start_dt}")
        self._write_line(_BORDER)
        self._write_line("")
        self._write_line("Invocation:")
        self._write_line(f"  {shlex.join(sys.argv)}")
        self._write_line("")

    def _write_footer(self, status: str) -> None:
        end_dt = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        elapsed = 0.0 if self._start_time is None else time.monotonic() - self._start_time
        minutes = int(elapsed // 60)
        seconds = int(round(elapsed - minutes * 60))
        if seconds == 60:
            minutes += 1
            seconds = 0
        self._write_line("")
        self._write_line(f"{status} {end_dt}")
        self._write_line(f"Elapsed: {minutes}:{seconds:02d}")
        self._write_line(_BORDER)

    def _write_items(self, title: str, items: dict[str, Any]) -> None:
        self._write_line(f"{title}:")
        for label, value in items.items():
            self._write_line(f"  {label:<18} {value}")
        self._write_line("")

    def _write_line(self, line: str) -> None:
        if self._handler is None:
            return
        self._handler.acquire()
        try:
            self._handler.stream.write(f"{line}\n")
            self._handler.flush()
        finally:
            self._handler.release()


__all__ = [
    "configure_package_logging",
    "format_float",
    "log_inputs",
    "log_outputs",
    "workflow_logging",
]
