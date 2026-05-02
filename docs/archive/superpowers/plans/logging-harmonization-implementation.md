# Logging Harmonization — Implementation Plan
# `ldsc_py3_restructured`

**Date:** 2026-05-01
**Status:** Approved — ready for implementation
**Revision:** 2026-05-01b (post-review fixes — see Revision Notes at the end)
**Design doc:** `logging-harmonization-design.md` (same directory)

---

## Pre-Implementation Checklist

- [ ] Working directory: `ldsc_py3_restructured` (the `restructure` branch worktree)
- [ ] Conda env active: `ldsc3`
- [ ] Baseline tests pass: `pytest` (confirm green before starting)

---

## Task 1 — Create `src/ldsc/_logging.py`

**File to create:** `src/ldsc/_logging.py`

Implement the full central logging module:

```python
"""Centralized workflow-run logging infrastructure.

This module provides the `workflow_run` decorator and `WorkflowRunContext`
context manager used by every LDSC public workflow to emit a uniform log
file alongside scientific outputs.

Key design properties:

- **Single shared ancestor logger.** The FileHandler is attached to
  `logging.getLogger("LDSC")`. All module loggers (`LDSC.ldscore_calculator`,
  `LDSC.regression_runner`, ...) propagate up to it, so one handler captures
  every record produced during a run.

- **Nested-context guard.** The decorator can be applied to both an outer
  `run_*_from_args` function and an inner class method using the same
  `workflow_name`. Whichever decorated function is hit first opens the
  context; subsequent decorated calls in the same run detect the active
  context and just call through. There is exactly one log file per run.

- **Python-API safety.** When no `logging.basicConfig` has been called (Python
  API usage), the root logger defaults to WARNING and INFO records are
  swallowed. The context temporarily lifts the LDSC logger level to INFO and
  restores it on exit.
"""
from __future__ import annotations

import logging
import shlex
import sys
import threading
import time
from datetime import datetime
from functools import wraps
from pathlib import Path
from typing import Any


FLOAT_FMT = ".4g"
_BORDER = "=" * 51
_LDSC_LOGGER_NAME = "LDSC"
_CURRENT: threading.local = threading.local()


def format_float(value: float) -> str:
    """Consistent float formatting for all workflow log messages."""
    return format(value, FLOAT_FMT)


def _active_context() -> "WorkflowRunContext | None":
    """Return the currently-active WorkflowRunContext on this thread, or None."""
    return getattr(_CURRENT, "ctx", None)


class WorkflowRunContext:
    """Per-run logging context: file handler, header/footer, level lift.

    Lifecycle is owned by `workflow_run`; this should not be entered manually.
    Only the *outermost* decorated call enters this context for a given run.
    Nested decorated calls (e.g., `run_ldscore_from_args` →
    `LDScoreCalculator.run`) detect the active context via `_active_context()`
    and call through without re-entering.
    """

    def __init__(
        self,
        workflow_name: str,
        logger: logging.Logger,
        output_dir: str | Path | None,
    ) -> None:
        self._workflow_name = workflow_name
        self._logger = logger
        self._output_dir = Path(output_dir) if output_dir is not None else None
        self._handler: logging.FileHandler | None = None
        self._start_time: float | None = None
        self._saved_level: int | None = None

    def __enter__(self) -> "WorkflowRunContext":
        self._start_time = time.monotonic()

        # Lift the LDSC logger's level to at least INFO so Python API callers
        # (who may not have configured logging.basicConfig) still see workflow
        # records. The original level is restored on exit.
        ldsc_logger = logging.getLogger(_LDSC_LOGGER_NAME)
        self._saved_level = ldsc_logger.level
        if ldsc_logger.getEffectiveLevel() > logging.INFO:
            ldsc_logger.setLevel(logging.INFO)

        # Attach FileHandler when an output directory is available. Without an
        # output directory, the run logs only to the console handler set up by
        # cli.py (regression workflows allow --output-dir to be omitted).
        if self._output_dir is not None:
            self._output_dir.mkdir(parents=True, exist_ok=True)
            log_path = self._output_dir / f"{self._workflow_name}.log"
            self._handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
            self._handler.setFormatter(logging.Formatter("%(message)s"))
            ldsc_logger.addHandler(self._handler)

        _CURRENT.ctx = self

        start_dt = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self._logger.info(_BORDER)
        self._logger.info(f"LDSC {self._workflow_name} — Started {start_dt}")
        self._logger.info(_BORDER)
        self._logger.info("")
        self._logger.info("Call:")
        # shlex.join preserves quoting for argv entries containing whitespace,
        # producing a shell-reproducible command string.
        self._logger.info(f"  {shlex.join(sys.argv)}")
        self._logger.info("")
        return self

    def log_inputs(self, **kv: Any) -> None:
        """Log labeled input paths under an 'Inputs:' header."""
        self._logger.info("Inputs:")
        for label, value in kv.items():
            self._logger.info(f"  {label:<16} {value}")
        self._logger.info("")

    def log_outputs(self, **kv: Any) -> None:
        """Log labeled output paths under an 'Outputs:' header."""
        self._logger.info("Outputs:")
        for label, value in kv.items():
            self._logger.info(f"  {label:<16} {value}")
        self._logger.info("")

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:  # noqa: ANN001
        elapsed = time.monotonic() - (self._start_time or 0.0)
        end_dt = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        minutes = int(elapsed // 60)
        seconds = int(elapsed % 60)
        # Distinguish success from failure in the footer so log readers can
        # tell at a glance whether the run completed cleanly.
        status = "Failed" if exc_type is not None else "Finished"
        self._logger.info(f"{status} {end_dt}")
        self._logger.info(f"Elapsed: {minutes}:{seconds:02d}")
        self._logger.info(_BORDER)

        ldsc_logger = logging.getLogger(_LDSC_LOGGER_NAME)
        if self._handler is not None:
            ldsc_logger.removeHandler(self._handler)
            self._handler.close()
            self._handler = None
        if self._saved_level is not None:
            ldsc_logger.setLevel(self._saved_level)
            self._saved_level = None
        _CURRENT.ctx = None
        return False


# --- Module-level thread-local helpers ---

def log_inputs(**kv: Any) -> None:
    """Forward to the active WorkflowRunContext, if any. No-op otherwise."""
    ctx = _active_context()
    if ctx is not None:
        ctx.log_inputs(**kv)


def log_outputs(**kv: Any) -> None:
    """Forward to the active WorkflowRunContext, if any. No-op otherwise."""
    ctx = _active_context()
    if ctx is not None:
        ctx.log_outputs(**kv)


# --- Decorator ---

def workflow_run(workflow_name: str, *, get_output_dir=None, get_overwrite=None):
    """Decorator factory wrapping a workflow entry point with WorkflowRunContext.

    Behavior
    --------
    1. **Nested-context guard.** If a `WorkflowRunContext` is already active
       (i.e., this call is nested inside another decorated workflow function
       in the same call chain), the wrapper skips context setup, the
       overwrite check, and the header/footer entirely — it just calls the
       wrapped function. The outermost decorated call owns the context.
    2. **Top-level entry.** Otherwise, extracts `output_dir` / `overwrite`
       from the wrapped function's arguments via the supplied extractors (or
       defaults), enforces the output-bundle overwrite policy on the log
       file, and enters a new context.

    Output-bundle overwrite policy
    ------------------------------
    If `output_dir` is set, `overwrite` is False, and `<output_dir>/<workflow_name>.log`
    already exists, the wrapper raises `FileExistsError` *before* any file is
    written. This is the same gate scientific output paths use, applied to
    the log file as part of the "either all rewrite or none rewrite" bundle
    policy.

    Parameters
    ----------
    workflow_name : str
        Used in the log header and as the log filename stem
        (`<output_dir>/<workflow_name>.log`). When the same `workflow_name` is
        applied at multiple call-chain layers, the nested-context guard
        ensures only one log file is produced per run.
    get_output_dir, get_overwrite : callable or None
        Optional `(args_tuple, kwargs_dict) -> value` extractors. When
        omitted, defaults to `args[0].output_dir` and `args[0].overwrite`,
        which works for `run_*_from_args(args)` patterns. Class methods whose
        configuration lives deeper in the argument list (e.g., `output_config`
        as positional index 5) supply explicit lambdas.
    """
    def decorator(fn):
        @wraps(fn)
        def wrapper(*args, **kwargs):
            # Nested decorated call: outer context owns the run.
            if _active_context() is not None:
                return fn(*args, **kwargs)

            first = args[0] if args else None
            output_dir = (
                get_output_dir(args, kwargs)
                if get_output_dir is not None
                else getattr(first, "output_dir", None)
            )
            overwrite = (
                get_overwrite(args, kwargs)
                if get_overwrite is not None
                else getattr(first, "overwrite", False)
            )
            if output_dir:
                log_path = Path(output_dir) / f"{workflow_name}.log"
                if log_path.exists() and not overwrite:
                    raise FileExistsError(
                        f"Log file {str(log_path)!r} already exists. "
                        "Pass --overwrite to replace it."
                    )
            with WorkflowRunContext(
                workflow_name, logging.getLogger(_LDSC_LOGGER_NAME), output_dir
            ):
                return fn(*args, **kwargs)
        return wrapper
    return decorator
```

**Verification:** `python -c "from ldsc._logging import workflow_run, format_float, WorkflowRunContext; print('OK')"` — must not raise.

---

## Task 2 — Refactor `src/ldsc/_kernel/sumstats_munger.py`

This is the largest change. Work carefully — 35 call sites.

### 2a — Remove legacy infrastructure

Delete (or replace with nothing):

- Lines 59–71: `sec_to_str` function
- Lines 74–101: `Logger` class
- Line 36: `import time` — keep (still used for `time.time()` elsewhere? Check; remove if unused after changes)
- Line 47: `MASTHEAD` multiline string and `__version__` — remove (no longer used)

### 2b — Add float format lambda at module level (after imports)

```python
_fmt = lambda v: format(v, ".4g")
```

### 2c — Replace the kernel entry `munge_sumstats()` preamble

Around line 799–821, replace:

```python
START_TIME = time.time()
log = Logger(args.out + '.log')
try:
    if args.sumstats is None:
        raise ValueError(...)
    ...
    if p:
        defaults = vars(parser.parse_args(''))
        opts = vars(args)
        header = MASTHEAD
        header += "Call: \n"
        ...
        log.log(header)
```

With:

```python
try:
    if args.sumstats is None:
        raise ValueError(...)
    ...
    # (remove the entire 'if p:' MASTHEAD/Call block — now handled by WorkflowRunContext)
```

### 2d — Replace all 35 `log.*` call sites

Mechanical substitutions throughout the file:

| Old | New |
|---|---|
| `log.log(msg)` | `LOGGER.info(msg)` |
| `log.warning(msg)` | `LOGGER.warning(msg)` |
| `log.error(msg)` | `LOGGER.error(msg)` |

### 2e — Replace float formatting (3 sites near line 1010–1015)

```python
# Old:
log.log('Mean chi^2 = ' + str(round(mean_chisq, 3)))
log.log('Lambda GC = ' + str(round(CHISQ.median() / 0.4549, 3)))
log.log('Max chi^2 = ' + str(round(CHISQ.max(), 3)))

# New:
LOGGER.info(f"Mean chi^2 = {_fmt(mean_chisq)}")
LOGGER.info(f"Lambda GC = {_fmt(CHISQ.median() / 0.4549)}")
LOGGER.info(f"Max chi^2 = {_fmt(CHISQ.max())}")
```

### 2f — Remove the `finally` block cleanup

Near line 1023–1027:

```python
# Remove:
finally:
    log.log('\nConversion finished at {T}'.format(T=time.ctime()))
    log.log('Total time elapsed: {T}'.format(T=sec_to_str(round(time.time() - START_TIME, 2))))
    log.close()
```

Timing and the "finished" message are now handled by `WorkflowRunContext.__exit__()`.
Keep the `finally:` block only if there are other cleanup actions; otherwise remove it.

**Verification:** `python -c "from ldsc._kernel import sumstats_munger; print('OK')"` — must not raise.

---

## Task 3 — Refactor `src/ldsc/sumstats_munger.py`

### 3a — Add import

```python
from ._logging import workflow_run
import logging as _logging_module  # only if needed for setLevel
from . import _logging as _ldsc_logging
```

### 3b — Extract `run_munge_sumstats_from_args`

Create a new decorated function containing all the current logic from `main()`:

```python
@workflow_run("munge-sumstats")
def run_munge_sumstats_from_args(args: argparse.Namespace):
    """Run munge-sumstats from a parsed CLI namespace."""
    # Move all current main() body here (everything after parse_args):
    args.sumstats = resolve_scalar_path(args.raw_sumstats_file, label="raw sumstats")
    output_dir = ensure_output_directory(args.output_dir, label="output directory")
    args.out = str(output_dir / "sumstats")
    metadata_path = args.out + ".metadata.json"
    ensure_output_paths_available(
        [args.out + ".sumstats.gz", metadata_path],
        overwrite=getattr(args, "overwrite", False),
        label="munged output artifact",
    )
    # NOTE: "sumstats.log" removed from preflight — log file handled by decorator
    args.merge_alleles = None
    if getattr(args, "sumstats_snps_file", None):
        args.sumstats_snps = resolve_scalar_path(args.sumstats_snps_file, label="sumstats SNPs file")
    else:
        args.sumstats_snps = None
    config_snapshot = _resolve_main_global_config(args)
    sumstats_snps_label = "none" if args.sumstats_snps is None else str(args.sumstats_snps)
    _ldsc_logging.log_inputs(
        source=str(args.sumstats),
        output_dir=str(output_dir),
        snp_identifier=config_snapshot.snp_identifier,
        genome_build=str(config_snapshot.genome_build),
        sumstats_snps=sumstats_snps_label,
    )
    LOGGER.info(
        f"Munging summary statistics from '{args.sumstats}' into '{output_dir}' "
        f"with snp_identifier='{config_snapshot.snp_identifier}', "
        f"genome_build='{config_snapshot.genome_build}'."
    )
    data = kernel_munge.munge_sumstats(args, p=True)
    coordinate_metadata = dict(getattr(args, "_coordinate_metadata", {}))
    config_snapshot = _effective_sumstats_config(config_snapshot, coordinate_metadata)
    _write_sumstats_metadata(
        metadata_path,
        config_snapshot=config_snapshot,
        coordinate_metadata=coordinate_metadata,
        source_path=args.sumstats,
        sumstats_file=args.out + ".sumstats.gz",
    )
    LOGGER.info(f"Munged {len(data)} retained summary-statistics rows; wrote '{args.out}.sumstats.gz'.")
    _ldsc_logging.log_outputs(
        sumstats_gz=args.out + ".sumstats.gz",
        log=str(output_dir / "munge-sumstats.log"),
        metadata_json=metadata_path,
    )
    return data
```

### 3c — Slim `main()` to parse-and-delegate

```python
def main(argv: list[str] | None = None):
    """CLI entry point: parse args and delegate to run_munge_sumstats_from_args."""
    return run_munge_sumstats_from_args(build_parser().parse_args(argv))
```

### 3d — Decorate `SumstatsMunger.run()` (Python API path)

The actual class method is named `run()` (not `munge()`). It is decorated as
the Python API entry point for munging. The CLI does NOT call this method —
`run_munge_sumstats_from_args` calls the kernel directly — so there is no
nesting here.

The `munge_config` argument holds `output_dir` and `overwrite`. When called
with a single config (no explicit `munge_config`), `raw_sumstats_config` is
used as the munge config after normalization — so fall back to `a[1]`.

```python
from ._logging import workflow_run
from . import _logging as _ldsc_logging

@workflow_run(
    "munge-sumstats",
    get_output_dir=lambda a, kw: (
        getattr(kw.get("munge_config") or (a[2] if len(a) > 2 else None), "output_dir", None)
        or getattr(a[1] if len(a) > 1 else None, "output_dir", None)
    ),
    get_overwrite=lambda a, kw: bool(
        getattr(kw.get("munge_config") or (a[2] if len(a) > 2 else None), "overwrite", False)
        or getattr(a[1] if len(a) > 1 else None, "overwrite", False)
    ),
)
def run(
    self,
    raw_sumstats_config: MungeConfig,
    munge_config: MungeConfig | None = None,
    global_config: GlobalConfig | None = None,
) -> SumstatsTable:
    ...
```

Also update `MungeRunSummary.output_paths["log"]` to use `"munge-sumstats.log"` instead of `"sumstats.log"`.
Similarly update the preflight check in `SumstatsMunger.run()` body:
- Remove `fixed_output_stem + ".log"` from `ensure_output_paths_available` (log file is now handled by decorator)
- Update `_last_summary.output_paths["log"]` to point to `str(output_dir / "munge-sumstats.log")`

Add `log_inputs` / `log_outputs` calls inside `SumstatsMunger.run()` body as well (matching the CLI path).

**Verification:** `python -c "from ldsc.sumstats_munger import run_munge_sumstats_from_args; print('OK')"` — must not raise.

---

## Task 4 — Modify `src/ldsc/ldscore_calculator.py`

`run_ldscore_from_args` performs significant orchestration *before* calling
`LDScoreCalculator.run` (argument normalization, validation, regression-SNP
loading, annotation building, ref-panel loading, the
`"Starting LD-score workflow..."` INFO message). To capture all of these in
the log file, **decorate both** the outer function and the inner method with
the same `workflow_name="ldscore"`. The nested-context guard in
`workflow_run` ensures only one log file is opened per run.

### 4a — Add imports

```python
from ._logging import workflow_run
from . import _logging as _ldsc_logging
```

### 4b — Decorate `run_ldscore_from_args`

Uses the default extractor (`args[0].output_dir`, `args[0].overwrite`).

```python
@workflow_run("ldscore")
def run_ldscore_from_args(args: argparse.Namespace) -> LDScoreResult:
    ...
```

### 4c — Decorate `LDScoreCalculator.run`

Required to give direct Python API users (`LDScoreCalculator().run(...)`) a
log file when they don't go through `run_ldscore_from_args`. When called via
`run_ldscore_from_args`, the nested-context guard makes this decorator a
pass-through.

`output_config` is the 6th positional argument (index 5 after `self`), or keyword `output_config`.

```python
@workflow_run(
    "ldscore",
    get_output_dir=lambda a, kw: getattr(
        kw.get("output_config") or (a[5] if len(a) > 5 else None), "output_dir", None
    ),
    get_overwrite=lambda a, kw: getattr(
        kw.get("output_config") or (a[5] if len(a) > 5 else None), "overwrite", False
    ),
)
def run(self, annotation_bundle, ref_panel, ldscore_config, global_config,
        output_config=None, regression_snps=None, config_snapshot=None):
    ...
```

### 4d — Add `log_inputs` / `log_outputs`

Place `log_inputs` near the top of `run_ldscore_from_args` (after argument
normalization, before annotation building) so input paths are recorded
before any heavy I/O. Place `log_outputs` near the end of
`LDScoreCalculator.run` after the output writer has written its artifacts.

```python
# In run_ldscore_from_args, after _validate_run_args(...)
_ldsc_logging.log_inputs(
    ref_panel=str(normalized_args.ref_panel) if getattr(normalized_args, "ref_panel", None) else "none",
    baseline_annot=str(normalized_args.baseline_annot_sources or "none"),
    output_dir=str(normalized_args.output_dir),
)

# In LDScoreCalculator.run, near the bottom (after writing)
_ldsc_logging.log_outputs(
    ldscore_dir=str(output_config.output_dir) if output_config else "none",
)
```

(Adjust the input keys to match the actual normalized arg attribute names.)

**Verification:** `python -c "from ldsc.ldscore_calculator import LDScoreCalculator, run_ldscore_from_args; print('OK')"`.

---

## Task 5 — Modify `src/ldsc/regression_runner.py`

### 5a — Add imports

```python
from ._logging import workflow_run, format_float
from . import _logging as _ldsc_logging
```

### 5b — Add `--log-level` argument

In `_add_common_regression_arguments()`:

```python
parser.add_argument(
    "--log-level",
    default="INFO",
    choices=("DEBUG", "INFO", "WARNING", "ERROR"),
    help="Logging verbosity.",
)
```

### 5c — Wire `--log-level` in `_runner_from_args()`

```python
log_level = getattr(args, "log_level", "INFO")
logging.getLogger("LDSC").setLevel(getattr(logging, log_level))
```

### 5d — Decorate `run_h2_from_args`, `run_partitioned_h2_from_args`, `run_rg_from_args`

All three use `args.output_dir` and `args.overwrite` (or defaults) — so the default extractor works:

```python
@workflow_run("h2")
def run_h2_from_args(args): ...

@workflow_run("partitioned-h2")
def run_partitioned_h2_from_args(args): ...

@workflow_run("rg")
def run_rg_from_args(args): ...
```

### 5e — Add `log_inputs` / `log_outputs` in each function

`run_h2_from_args`:
```python
_ldsc_logging.log_inputs(sumstats=args.sumstats_file, ldscore_dir=args.ldscore_dir)
# ... workflow ...
_ldsc_logging.log_outputs(h2_tsv=str(Path(output_dir) / "h2.tsv") if output_dir else "none")
```

`run_partitioned_h2_from_args`:
```python
_ldsc_logging.log_inputs(sumstats=args.sumstats_file, ldscore_dir=args.ldscore_dir)
_ldsc_logging.log_outputs(output_dir=output_dir or "none")
```

`run_rg_from_args`:
```python
_ldsc_logging.log_inputs(
    sumstats_1=args.sumstats_1_file,
    sumstats_2=args.sumstats_2_file,
    ldscore_dir=args.ldscore_dir,
)
_ldsc_logging.log_outputs(rg_tsv=str(Path(output_dir) / "rg.tsv") if output_dir else "none")
```

### 5f — Add float formatting to key result log messages

Find the existing `LOGGER.info("Finished h2 regression with ...")` and similar lines. Where
float statistics (h2, SE, rg, p-value) are logged inline, wrap with `format_float()`:

```python
# Example — add after hsq is available in run_h2_from_args:
if hasattr(hsq, "tot") and hsq.tot is not None:
    LOGGER.info(f"h2 = {format_float(hsq.tot)} (SE: {format_float(hsq.tot_se)})")
```

**Verification:** `python -c "from ldsc.regression_runner import run_h2_from_args; print('OK')"`.

---

## Task 6 — Modify `src/ldsc/ref_panel_builder.py`

Same outer+inner decoration pattern as ldscore (Task 4): decorate **both**
`run_build_ref_panel_from_args` and `ReferencePanelBuilder.run` with
`workflow_name="build-ref-panel"`. The nested-context guard ensures one log
file per run, while still giving direct Python API users of the class method
a log file of their own.

### 6a — Add imports

```python
from ._logging import workflow_run
from . import _logging as _ldsc_logging
```

### 6b — Decorate `run_build_ref_panel_from_args`

Uses the default extractor (`args[0].output_dir`, `args[0].overwrite`).

```python
@workflow_run("build-ref-panel")
def run_build_ref_panel_from_args(args: argparse.Namespace) -> ReferencePanelBuildResult:
    ...
```

### 6c — Decorate `ReferencePanelBuilder.run`

`config` is the only non-self argument (positional index 1):

```python
@workflow_run(
    "build-ref-panel",
    get_output_dir=lambda a, kw: getattr(
        kw.get("config") or (a[1] if len(a) > 1 else None), "output_dir", None
    ),
    get_overwrite=lambda a, kw: getattr(
        kw.get("config") or (a[1] if len(a) > 1 else None), "overwrite", False
    ),
)
def run(self, config: ReferencePanelBuildConfig) -> ReferencePanelBuildResult:
    ...
```

### 6d — Add `log_inputs` / `log_outputs`

Place `log_inputs` near the top of `run_build_ref_panel_from_args` (after
argument resolution). Place `log_outputs` near the end of
`ReferencePanelBuilder.run` after the result is assembled.

```python
# In run_build_ref_panel_from_args
_ldsc_logging.log_inputs(
    plink_prefix=str(args.plink_prefix),
    output_dir=str(args.output_dir),
)

# In ReferencePanelBuilder.run, near the bottom
_ldsc_logging.log_outputs(output_dir=str(config.output_dir))
```

**Verification:** `python -c "from ldsc.ref_panel_builder import ReferencePanelBuilder, run_build_ref_panel_from_args; print('OK')"`.

---

## Task 7 — Full Test Suite

```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3
cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured
pytest
```

All tests must pass. If any fail, diagnose and fix before marking complete.

---

## Task 8 — Manual Smoke Tests

Run each workflow that has end-to-end test data available and verify:

1. Log file appears in output directory with correct name (`munge-sumstats.log`, `ldscore.log`, etc.)
2. Log file contains:
   - Border + `LDSC <name> — Started <datetime>` header
   - `Call:` section with `sys.argv`
   - `Inputs:` section
   - Workflow INFO messages
   - `Outputs:` section
   - `Finished` + `Elapsed: mm:ss` footer
3. Float values in log use `:.4g` format
4. Re-run without `--overwrite` raises `FileExistsError` on the log file
5. Re-run with `--overwrite` succeeds and overwrites all outputs including the log

---

## Verification Checklist

- [ ] Task 1: `_logging.py` imports cleanly
- [ ] Task 2: kernel imports cleanly; no `Logger` class, no `sec_to_str`
- [ ] Task 3: `run_munge_sumstats_from_args` exists and is exported; `main()` is one line
- [ ] Task 4: both `run_ldscore_from_args` and `LDScoreCalculator.run` decorated with `@workflow_run("ldscore")`; CLI run produces a single `ldscore.log` (no double header); direct `LDScoreCalculator().run(...)` produces a log too
- [ ] Task 5: `run_h2_from_args`, `run_partitioned_h2_from_args`, `run_rg_from_args` decorated; `--log-level` works
- [ ] Task 6: both `run_build_ref_panel_from_args` and `ReferencePanelBuilder.run` decorated with `@workflow_run("build-ref-panel")`; CLI run produces a single `build-ref-panel.log`
- [ ] Task 7: `pytest` green
- [ ] Task 8: Smoke test log files look correct

---

## Known Edge Cases

### Nested decorated calls (ldscore, build-ref-panel)
For `ldscore`, both `run_ldscore_from_args` and `LDScoreCalculator.run` are
decorated with `@workflow_run("ldscore")`. CLI/`run_ldscore(**kwargs)` paths
hit the outer decorator first; it enters the context, and the inner
decorator on `LDScoreCalculator.run` detects the active context (via
`_active_context()`) and calls through. Result: one `ldscore.log` per run,
one header/footer pair. Direct Python API calls to `LDScoreCalculator().run(...)`
hit only the inner decorator, which opens its own context. Same logic for
`build-ref-panel`. ✓

### No `output_dir` for regression
`run_h2_from_args`, `run_rg_from_args`, `run_partitioned_h2_from_args` have optional
`--output-dir`. When `None`, `output_dir` extracted by the decorator is `None`, so no
FileHandler is attached and no log file is written. Console logging continues normally. ✓

### `SumstatsMunger.run` with single config argument
When called as `munger.run(MungeConfig(raw_sumstats_file=..., output_dir=...))`, the
`munge_config` kwarg is `None`. The lambda falls back to `a[1]` (the `raw_sumstats_config`)
which holds `output_dir`. ✓

### Python API logger level
Python API users may not have called `logging.basicConfig`. Without it, the
root logger defaults to WARNING and INFO records would be silently dropped,
producing an empty log file. `WorkflowRunContext.__enter__` lifts the LDSC
logger level to INFO when `getEffectiveLevel()` is above INFO, and restores
it on exit. ✓

### Failed runs
When the wrapped function raises, `WorkflowRunContext.__exit__` writes
`Failed <datetime>` (instead of `Finished`) before re-raising the exception
(the context returns `False`, never suppressing). The handler is still
removed and the level still restored. ✓

### `sys.argv` quoting
The `Call:` line uses `shlex.join(sys.argv)`, which preserves quoting for
arguments containing whitespace and produces a shell-reproducible command
string. ✓

### Annotation subcommand
`ldsc annotate` / `annotation_builder.py` is entirely excluded from this task.
The `_kernel/annotation.py` is not modified. ✓

---

## Revision Notes (2026-05-01b)

Post-review fixes incorporated into this plan:

1. **Nested-context guard added to `workflow_run`.** Allows the same
   `workflow_name` at outer and inner call-chain layers; only the outermost
   call enters the context.
2. **Outer + inner decoration for `ldscore` and `build-ref-panel`.** Captures
   orchestration messages from `run_*_from_args` (validation, annotation
   loading, ref-panel loading) in the log file, while still giving direct
   class-method Python API users a log file of their own.
3. **LDSC logger level lift in `__enter__`.** Ensures Python API callers
   (without `logging.basicConfig`) see INFO records in the log file.
4. **`shlex.join(sys.argv)`** instead of `' '.join(sys.argv)` for the
   `Call:` line — preserves shell-safe quoting.
5. **`Failed` vs `Finished` footer** based on `exc_type` in `__exit__`.
6. **Method name corrected:** `SumstatsMunger.run()` (not `.munge()`) is the
   actual Python API method. All references updated.

### Out-of-scope notes
- *Output bundle atomicity:* If the log file is manually deleted while
  scientific outputs survive, the decorator's existence check passes and the
  log file is opened before `ensure_output_paths_available` raises on the
  surviving outputs. This leaves an orphan log file. Acceptable: this is a
  manually-induced inconsistent state.
- *Duplicate error logging through `cli.py`:* `cli.run_cli` catches
  user-facing exceptions and logs them. Kernel error sites also log via
  `LOGGER.error(...)` after this plan. Some errors may appear twice in the
  console. Fixing this requires changes to `cli.py` exception handling,
  which is out of scope here.
