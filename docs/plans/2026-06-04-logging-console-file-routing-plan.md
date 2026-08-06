# Logging Console/File Routing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the per-run `.log` file the authoritative sink for all workflow messages when an output directory is provided, keep the console (stderr) silent except for errors, fall back to console INFO for no-output quick runs, capture run-aborting tracebacks into the `.log` file, and never emit console output from the Python API.

**Architecture:** Console logging becomes a CLI-only concern owned by `run_cli`, which installs a single stderr `StreamHandler` on the `LDSC` logger at `ERROR` level (removed in a `finally`). `workflow_logging` keeps the `FileHandler` on `LDSC`, lowers the CLI console handler to `INFO` only when there is no log file (quick run), and writes the full traceback into the file on a failed exit before the handler is torn down. The root logger never receives a handler again — the three `basicConfig`/root-handler sources are removed — so `LDSC.propagate` can stay `True` and pytest `caplog` keeps working.

**Tech Stack:** Python stdlib `logging`, `traceback`; pytest + unittest.

---

## File Structure

- `src/ldsc/_logging.py` — gains CLI console-handler install/remove/accessor helpers, a thread-local "last workflow log path" for the error-boundary pointer, console-level toggling when `log_path is None`, and traceback-into-file capture on failed exit. Renames the `__exit__` traceback parameter to avoid shadowing the new `import traceback`.
- `src/ldsc/cli.py` — `run_cli` replaces `logging.basicConfig(...)` with the CLI console handler lifecycle and emits concise, traceback-free console messages (internal errors gain a logfile pointer).
- `src/ldsc/_kernel/ldscore.py` — `configure_logging` stops calling `basicConfig` (which installs a root handler) and only sets the `LDSC` logger threshold.
- `tests/test_logging_refactor.py` — two existing behaviors are intentionally reversed (traceback now in file; console boundary record now concise) plus new coverage for console routing.

---

## Task 1: Capture the run-aborting traceback into the `.log` file

**Files:**
- Modify: `src/ldsc/_logging.py` (imports ~26-35; `__exit__` 174-194; `_write_footer` 211-222)
- Test: `tests/test_logging_refactor.py`

- [ ] **Step 1: Update the failing-footer test to expect the traceback in the file**

Replace the existing `test_workflow_logging_failed_footer_does_not_record_exception_message` (lines 115-128) with:

```python
    def test_workflow_logging_failed_footer_records_traceback(self):
        from ldsc._logging import workflow_logging

        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "workflow.log"
            with self.assertRaisesRegex(RuntimeError, "boom"):
                with workflow_logging("unit", log_path, log_level="INFO"):
                    raise RuntimeError("boom")

            text = log_path.read_text(encoding="utf-8")

        self.assertIn("Failed", text)
        self.assertIn("Elapsed time:", text)
        self.assertIn("Traceback (most recent call last):", text)
        self.assertIn("RuntimeError: boom", text)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured && python -m pytest tests/test_logging_refactor.py::LoggingRefactorTest::test_workflow_logging_failed_footer_records_traceback -v`
Expected: FAIL — "Traceback (most recent call last):" not found in file.

- [ ] **Step 3: Add `import traceback` and write the traceback in `_write_footer`**

In `src/ldsc/_logging.py`, add `import traceback` to the stdlib import block (after `import time`):

```python
import threading
import time
import traceback
```

Rename the `__exit__` traceback parameter to avoid shadowing the module, and pass exception info to the footer. Replace the `__exit__` body's footer call (line ~183-184):

```python
    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> bool:
        if self._nested:
            return False

        status = "Failed" if exc_type is not None else "Finished"
        self._write_footer(status, exc_type, exc_value, exc_tb)
```

Replace `_write_footer` (lines 211-222) to accept and render the traceback before the closing border:

```python
    def _write_footer(
        self,
        status: str,
        exc_type: type[BaseException] | None = None,
        exc_value: BaseException | None = None,
        exc_tb: TracebackType | None = None,
    ) -> None:
        end_dt = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        elapsed = 0.0 if self._start_time is None else time.monotonic() - self._start_time
        minutes = int(elapsed // 60)
        seconds = int(round(elapsed - minutes * 60))
        if seconds == 60:
            minutes += 1
            seconds = 0
        self._write_line("")
        self._write_line(f"{status} {end_dt}")
        self._write_line(f"Elapsed time: {float(minutes):.1f}min:{seconds}s")
        if exc_type is not None:
            self._write_line("")
            formatted = "".join(traceback.format_exception(exc_type, exc_value, exc_tb))
            for line in formatted.rstrip("\n").split("\n"):
                self._write_line(line)
        self._write_line(_BORDER)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured && python -m pytest tests/test_logging_refactor.py -v`
Expected: PASS (all logging-refactor tests, including the renamed one).

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_logging.py tests/test_logging_refactor.py
git commit -m "feat(logging): capture run-aborting traceback in workflow log file"
```

---

## Task 2: CLI console handler + last-log-path helpers in `_logging`

**Files:**
- Modify: `src/ldsc/_logging.py` (module state ~38-41; `_WorkflowLoggingContext.__init__` 128-141; `__enter__` 143-172; `__exit__` 174-194; `__all__` 275-281)
- Test: `tests/test_logging_refactor.py`

- [ ] **Step 1: Write the failing tests for console helpers and level toggling**

Append to `tests/test_logging_refactor.py` (inside `LoggingRefactorTest`):

```python
    def test_cli_console_handler_install_remove_is_clean(self):
        from ldsc import _logging

        ldsc_logger = logging.getLogger("LDSC")
        original = list(ldsc_logger.handlers)
        handler = _logging.install_cli_console_handler()
        try:
            self.assertIn(handler, ldsc_logger.handlers)
            self.assertEqual(handler.level, logging.ERROR)
            self.assertIs(_logging.cli_console_handler(), handler)
        finally:
            _logging.remove_cli_console_handler()
        self.assertEqual(ldsc_logger.handlers, original)
        self.assertIsNone(_logging.cli_console_handler())

    def test_console_level_lowered_to_info_when_no_log_path(self):
        from ldsc import _logging
        from ldsc._logging import workflow_logging

        handler = _logging.install_cli_console_handler()
        try:
            self.assertEqual(handler.level, logging.ERROR)
            with workflow_logging("unit", None, log_level="INFO"):
                self.assertEqual(handler.level, logging.INFO)
            self.assertEqual(handler.level, logging.ERROR)
        finally:
            _logging.remove_cli_console_handler()

    def test_console_level_stays_error_when_log_path_present(self):
        from ldsc import _logging
        from ldsc._logging import workflow_logging

        handler = _logging.install_cli_console_handler()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                log_path = Path(tmpdir) / "workflow.log"
                with workflow_logging("unit", log_path, log_level="INFO"):
                    self.assertEqual(handler.level, logging.ERROR)
        finally:
            _logging.remove_cli_console_handler()

    def test_last_workflow_log_path_tracks_active_run(self):
        from ldsc import _logging
        from ldsc._logging import workflow_logging

        _logging.reset_workflow_log_path()
        self.assertIsNone(_logging.last_workflow_log_path())
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "workflow.log"
            with workflow_logging("unit", log_path, log_level="INFO"):
                pass
            self.assertEqual(_logging.last_workflow_log_path(), str(log_path))
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured && python -m pytest tests/test_logging_refactor.py -k "console or last_workflow" -v`
Expected: FAIL — `_logging` has no attribute `install_cli_console_handler`.

- [ ] **Step 3: Implement the console-handler helpers and path tracking**

In `src/ldsc/_logging.py`, extend the module state block (after line 41):

```python
_LDSC_LOGGER_NAME = "LDSC"
_CURRENT: threading.local = threading.local()
_BORDER = "=" * 51
_FLOAT_FORMAT = ".4g"
_CLI_CONSOLE_HANDLER: logging.Handler | None = None
```

Add these module-level functions (place them after `configure_package_logging`, before `workflow_logging`):

```python
def install_cli_console_handler(level: str = "ERROR") -> logging.Handler:
    """Attach a single stderr console handler to the ``LDSC`` logger for CLI runs.

    The handler is owned by the CLI boundary, not by individual workflows, so it
    outlives every :func:`workflow_logging` block and can still surface errors
    raised at the command boundary. It is a no-op-safe singleton: a second call
    returns the existing handler.
    """
    global _CLI_CONSOLE_HANDLER
    if _CLI_CONSOLE_HANDLER is not None:
        return _CLI_CONSOLE_HANDLER
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(_level_number(level))
    handler.setFormatter(logging.Formatter("%(message)s"))
    logging.getLogger(_LDSC_LOGGER_NAME).addHandler(handler)
    _CLI_CONSOLE_HANDLER = handler
    return handler


def remove_cli_console_handler() -> None:
    """Remove and close the CLI console handler if one is installed."""
    global _CLI_CONSOLE_HANDLER
    if _CLI_CONSOLE_HANDLER is None:
        return
    logging.getLogger(_LDSC_LOGGER_NAME).removeHandler(_CLI_CONSOLE_HANDLER)
    _CLI_CONSOLE_HANDLER.close()
    _CLI_CONSOLE_HANDLER = None


def cli_console_handler() -> logging.Handler | None:
    """Return the installed CLI console handler, or ``None``."""
    return _CLI_CONSOLE_HANDLER


def reset_workflow_log_path() -> None:
    """Clear the recorded most-recent workflow log path for this thread."""
    _CURRENT.last_log_path = None


def last_workflow_log_path() -> str | None:
    """Return the most recent workflow log path on this thread, or ``None``."""
    return getattr(_CURRENT, "last_log_path", None)
```

In `_WorkflowLoggingContext.__init__` (after line 141 `self._nested = False`), add console-toggle bookkeeping:

```python
        self._nested = False
        self._console_handler: logging.Handler | None = None
        self._saved_console_level: int | None = None
```

In `__enter__`, record the path and toggle the console level in the outermost (non-nested) branch. After `self._start_time = time.monotonic()` (line 151) add:

```python
        _CURRENT.last_log_path = str(self._log_path) if self._log_path is not None else None
        if self._log_path is None:
            handler = cli_console_handler()
            if handler is not None:
                self._console_handler = handler
                self._saved_console_level = handler.level
                handler.setLevel(self._level)
```

In `__exit__`, restore the console level during teardown. After the `if self._handler is not None:` block (lines 187-189) and before the level restore, add:

```python
        if self._console_handler is not None and self._saved_console_level is not None:
            self._console_handler.setLevel(self._saved_console_level)
```

Add the new public names to `__all__`:

```python
__all__ = [
    "cli_console_handler",
    "configure_package_logging",
    "format_float",
    "install_cli_console_handler",
    "last_workflow_log_path",
    "log_inputs",
    "log_outputs",
    "remove_cli_console_handler",
    "reset_workflow_log_path",
    "workflow_logging",
]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured && python -m pytest tests/test_logging_refactor.py -v`
Expected: PASS (all, including the four new tests).

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_logging.py tests/test_logging_refactor.py
git commit -m "feat(logging): add CLI console handler and log-path tracking helpers"
```

---

## Task 3: Route the CLI boundary through the console handler with concise errors

**Files:**
- Modify: `src/ldsc/cli.py` (imports 17-35; `run_cli` 164-186)
- Test: `tests/test_logging_refactor.py` (update `test_run_cli_logs_unexpected_errors_with_traceback`)

- [ ] **Step 1: Update the internal-error test to expect a concise, traceback-free console record**

Replace `test_run_cli_logs_unexpected_errors_with_traceback` (lines 153-166) with:

```python
    def test_run_cli_logs_unexpected_errors_concisely(self):
        from ldsc import cli

        with mock.patch.object(cli, "main", side_effect=RuntimeError("boom")):
            with self.assertLogs("LDSC.cli", level="ERROR") as caught:
                code = cli.run_cli(["ldscore"])

        self.assertEqual(code, 2)
        self.assertEqual(len(caught.records), 1)
        record = caught.records[0]
        self.assertEqual(record.levelname, "ERROR")
        self.assertIsNone(record.exc_info)
        self.assertIn("Internal error while running ldsc", record.getMessage())
        self.assertIn("boom", record.getMessage())
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured && python -m pytest tests/test_logging_refactor.py::LoggingRefactorTest::test_run_cli_logs_unexpected_errors_concisely -v`
Expected: FAIL — `record.exc_info` is currently set (the boundary still uses `LOGGER.exception`).

- [ ] **Step 3: Rewrite `run_cli` to own the console handler and emit concise errors**

In `src/ldsc/cli.py`, add the logging-helper import below the existing `from .errors import ...` (line 25):

```python
from . import annotation_builder, ldscore_calculator, ref_panel_builder
from ._logging import (
    install_cli_console_handler,
    last_workflow_log_path,
    remove_cli_console_handler,
    reset_workflow_log_path,
)
from .errors import LDSCError, LDSCUsageError, LDSCUserError
```

Replace `run_cli` (lines 164-186) with:

```python
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
```

Note: the full traceback is captured into the `.log` file by `workflow_logging.__exit__` (Task 1); the console deliberately shows only the concise line plus the logfile pointer.

- [ ] **Step 4: Run tests to verify they pass**

Run: `source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured && python -m pytest tests/test_logging_refactor.py -v`
Expected: PASS — including `test_run_cli_reports_user_errors_without_traceback` (unchanged) and `test_run_cli_logs_unexpected_errors_concisely`.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/cli.py tests/test_logging_refactor.py
git commit -m "feat(logging): route CLI console through LDSC handler with concise errors"
```

---

## Task 4: Stop the kernel `configure_logging` from installing a root handler

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` (`configure_logging` 286-291)
- Test: `tests/test_logging_refactor.py`

- [ ] **Step 1: Write the failing test that the kernel config adds no root handler**

Append to `tests/test_logging_refactor.py` (inside `LoggingRefactorTest`):

```python
    def test_kernel_configure_logging_adds_no_root_handler(self):
        from ldsc._kernel import ldscore

        root_logger = logging.getLogger()
        original = list(root_logger.handlers)
        ldscore.configure_logging("INFO")
        self.assertEqual(root_logger.handlers, original)
        self.assertEqual(logging.getLogger("LDSC").level, logging.INFO)
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured && python -m pytest tests/test_logging_refactor.py::LoggingRefactorTest::test_kernel_configure_logging_adds_no_root_handler -v`
Expected: FAIL — `basicConfig` adds a handler to the root logger.

- [ ] **Step 3: Replace `basicConfig` with an `LDSC` threshold set**

In `src/ldsc/_kernel/ldscore.py`, replace `configure_logging` (lines 286-291):

```python
def configure_logging(level: str) -> None:
    """Set the ``LDSC`` logger threshold for kernel-side execution.

    Console and file routing is owned by the public workflow boundary
    (``ldsc._logging`` / ``ldsc.cli``); the kernel only adjusts the level so
    direct kernel invocation does not install a root console handler.
    """
    logging.getLogger("LDSC").setLevel(getattr(logging, level.upper()))
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured && python -m pytest tests/test_logging_refactor.py::LoggingRefactorTest::test_kernel_configure_logging_adds_no_root_handler -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/ldscore.py tests/test_logging_refactor.py
git commit -m "fix(logging): kernel configure_logging sets level without root handler"
```

---

## Task 5: Full-suite verification

**Files:** none (verification only)

- [ ] **Step 1: Run the whole test suite**

Run: `source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured && python -m pytest -q`
Expected: PASS (no regressions; `test_liftover_logging.py` still passes because `LDSC.propagate` is unchanged).

- [ ] **Step 2: Confirm no root-handler `basicConfig` remains in the package**

Run: `cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured && grep -rn "basicConfig" src --include="*.py"`
Expected: no matches.

- [ ] **Step 3: Smoke-test the no-output console fallback path (optional manual check)**

Run a quick no-output regression (e.g. `ldsc h2 ...` without `--output-dir`) and confirm INFO progress appears on stderr; run a munge with `--output-dir` and confirm progress appears in `<output-dir>/diagnostics/*.log` and not on stderr.

---

## Self-Review

- **Spec coverage:** output-dir → file only (Task 2 keeps console at ERROR with file handler present); no-output → console INFO (Task 2 level toggle); errors → file full traceback (Task 1) + console concise (Task 3); internal error console pointer (Task 3 `_log_internal_error`); API → file only / no console (console handler only installed by `run_cli` when `argv is None`, never by the API path); kernel root-handler source removed (Task 4); `caplog` preserved (propagate untouched, verified Task 5).
- **Placeholder scan:** none — every code step shows full code.
- **Type consistency:** `install_cli_console_handler`/`remove_cli_console_handler`/`cli_console_handler`/`reset_workflow_log_path`/`last_workflow_log_path` names are used identically across `_logging.py`, `cli.py`, and tests. `_write_footer` signature `(status, exc_type, exc_value, exc_tb)` matches its single call site in `__exit__`.
