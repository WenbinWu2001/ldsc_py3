# Logging Harmonization Design
# `ldsc_py3_restructured`

**Date:** 2026-05-01
**Status:** Approved — pending refresh after munge-sumstats layering refactor
**Revision:** 2026-05-02a (munge-sumstats now has workflow-owned `sumstats.log` and shared CLI/API orchestration)

---

## Problem

The package has 7 CLI subcommands with fragmented logging behavior:

- **`munge-sumstats`** — now writes its current `sumstats.log` from the public
  workflow layer via package logging. The old kernel `Logger` class, masthead,
  timing footer, and kernel `main()` have been removed. It still lacks the
  future shared header/timing/command format and still uses the old
  `sumstats.log` filename.
- **All other workflows** (`ldscore`, `h2`, `partitioned-h2`, `rg`, `build-ref-panel`) — use
  Python's `logging` module only. No log file, no timing, no command recording.
- Float formatting is inconsistent across the package (2–3 decimal places, `%.3f` in CSV,
  unformatted floats elsewhere).
- Logging levels are mostly `INFO`-only; `munge-sumstats` no longer logs
  re-raised kernel exceptions as duplicate errors, leaving user-facing errors
  to `cli.run_cli()`.
- No shared logging utilities exist.

---

## Goals

1. Every workflow writes a `.log` file alongside its output artifacts.
2. Every log file records the command (`sys.argv`) at the start.
3. Every log file records start time, end time, and elapsed duration (mm:ss).
4. Input and output paths are recorded in every log file.
5. Float formatting is unified across all log messages.
6. Logging levels follow modern conventions.

---

## Design Decisions (All Confirmed)

| Decision | Choice | Rationale |
|---|---|---|
| No `--output-dir` provided | Skip log file; console only | No surprise files; users who want a log must provide an output dir |
| Command format | Raw `sys.argv` as-is | Simple, accurate, reproducible |
| Float precision | `:.4g` (4 significant figures) | Handles both large stats and tiny p-values without switching formats |
| `annotate` subcommand | Deferred | Being refactored separately |
| `munge-sumstats` log filename | `munge-sumstats.log` (was `sumstats.log`) | Consistent with all other workflows |
| Python API argv | Write `sys.argv` as-is | Simple and honest; Python users see `python my_script.py` |
| Code structure intrusion | Decorator pattern; no wrapping `with` blocks | Minimizes changes to existing function bodies |
| Log file in output bundle | Yes — same `--overwrite` gate as scientific outputs | "Either all rewrite or none rewrite" |

---

## Logging Level Design

| Level | Used for |
|---|---|
| `DEBUG` | Per-chromosome progress, intermediate stats, detailed config values, algorithm parameters |
| `INFO` | Workflow start/end, key results (h2, rg, Lambda GC, chi^2 mean), elapsed timing, I/O paths |
| `WARNING` | Data quality issues (low SNP overlap, unusual statistics), deprecated usage |
| `ERROR` | User-facing errors surfaced before exceptions are raised |

**Already fixed:** `_kernel/sumstats_munger.py` no longer calls `LOGGER.debug(msg)` for
re-raised conversion failures. User-facing exceptions now propagate to
`cli.run_cli()`, and data-quality issues use ordinary warning/info records.

**Missing:** `h2`, `partitioned-h2`, `rg` subcommands currently have no `--log-level` argument.
Add it to match `ldscore`, `build-ref-panel`, and `annotate`.

---

## Log File Format

```
===================================================
LDSC ldscore — Started 2026-05-01 14:32:07
===================================================

Call:
  ldsc ldscore --ref-panel /data/1000G --output-dir /results/ldscore --genome-build hg38

Inputs:
  ref_panel      /data/1000G
  output_dir     /results/ldscore

[... workflow INFO log messages appear here ...]

Outputs:
  ldscore_dir    /results/ldscore

Finished 2026-05-01 14:38:43
Elapsed: 6:36
===================================================
```

The log file uses `%(message)s` format only (no `INFO:` prefix in file).
The console handler (set up by `cli.py` via `basicConfig`) continues to show `INFO: ...` prefixes.

---

## Central Module: `src/ldsc/_logging.py`

New module. Public surface:

### `format_float(value: float) -> str`
Uniform float formatter: `format(value, ".4g")`. Use everywhere log messages include float stats.

### `class WorkflowRunContext`
Context manager. On `__enter__`:
1. Creates `output_dir` (with `exist_ok=True`) so the FileHandler can open.
2. Saves the current level of `logging.getLogger("LDSC")` and forces it to
   `INFO` if `getEffectiveLevel()` is above `INFO`. This ensures that Python
   API callers, who may not have configured `logging.basicConfig`, still see
   workflow records (the root logger defaults to `WARNING`).
3. Attaches a `FileHandler` with `%(message)s` format to `logging.getLogger("LDSC")`.
4. Logs the border, `"LDSC <name> — Started <datetime>"`, and the call as
   `shlex.join(sys.argv)` (preserves quoting for paths with spaces).
5. Stores `self` in a `threading.local` for access by `log_inputs`/`log_outputs` helpers.

On `__exit__`:
1. Logs `"Finished <datetime>"` (success) or `"Failed <datetime>"` (when
   `exc_type` is not `None`) and `"Elapsed: mm:ss"`.
2. Removes and closes the FileHandler.
3. Restores the LDSC logger's saved level.
4. Clears the `threading.local` slot.
5. Returns `False` (never suppresses exceptions).

Why attach to `logging.getLogger("LDSC")` (the shared ancestor)? All module loggers
(`LDSC.ldscore_calculator`, `LDSC.regression_runner`, etc.) propagate records up to this
ancestor. One handler here captures the entire run without touching every module.

### `log_inputs(**kv)` / `log_outputs(**kv)` (module-level functions)
Thread-local wrappers that forward to the active `WorkflowRunContext`. No-ops when no
context is active (e.g., unit tests). Usage:

```python
from . import _logging as _ldsc_logging
_ldsc_logging.log_inputs(source="/path/to/input.gwas", output_dir="/results")
_ldsc_logging.log_outputs(sumstats_gz="/results/sumstats.sumstats.gz")
```

### `workflow_run(workflow_name, *, get_output_dir=None, get_overwrite=None)`
Decorator factory. Wraps a workflow entry function with `WorkflowRunContext`.

**Nested-context guard:** the wrapper checks `_active_context()` *first*. If a
`WorkflowRunContext` is already active (i.e., the wrapped function is being
called from inside another decorated workflow function), the wrapper skips
context setup, the overwrite check, and the header/footer entirely — it just
calls the function. This is what allows the *same* `workflow_name` to be
applied to both an outer `run_*_from_args` function and an inner class method
without producing two log files, two headers, or two overwrite checks. The
outermost decorated call owns the context for the run.

**Output-bundle enforcement** (only on top-level entry, before `__enter__()`):
```python
if output_dir:
    log_path = Path(output_dir) / f"{workflow_name}.log"
    if log_path.exists() and not overwrite:
        raise FileExistsError(f"Log file {str(log_path)!r} already exists. Pass --overwrite.")
```

Default `get_output_dir`: `getattr(args[0], "output_dir", None)` — works for all
`run_*_from_args(args)` functions.

Default `get_overwrite`: `getattr(args[0], "overwrite", False)` — same.

Custom lambdas for class methods (where output_dir is deeper in the signature):
```python
# ReferencePanelBuilder.run(self, config)
get_output_dir=lambda a, kw: getattr(kw.get("config") or (a[1] if len(a)>1 else None), "output_dir", None)
get_overwrite= lambda a, kw: getattr(kw.get("config") or (a[1] if len(a)>1 else None), "overwrite", False)

# LDScoreCalculator.run(self, annotation_bundle, ref_panel, ldscore_config, global_config, output_config=None, ...)
get_output_dir=lambda a, kw: getattr(kw.get("output_config") or (a[5] if len(a)>5 else None), "output_dir", None)
get_overwrite= lambda a, kw: getattr(kw.get("output_config") or (a[5] if len(a)>5 else None), "overwrite", False)

# SumstatsMunger.run(self, raw_sumstats_config, munge_config=None, global_config=None)
get_output_dir=lambda a, kw: (
    getattr(kw.get("munge_config") or (a[2] if len(a)>2 else None), "output_dir", None)
    or getattr(a[1] if len(a)>1 else None, "output_dir", None)
)
get_overwrite=lambda a, kw: (
    getattr(kw.get("munge_config") or (a[2] if len(a)>2 else None), "overwrite", False)
    or getattr(a[1] if len(a)>1 else None, "overwrite", False)
)
```

---

## Call-Chain Strategy (Outer + Inner Decoration with Nested-Context Guard)

All workflows share the same structural pattern:
`main(argv)` → `run_*_from_args(args)` → (optional class method) → work.

The `workflow_run` decorator's nested-context guard makes it safe to apply the
*same* `workflow_name` at multiple layers in the same call chain. Whichever
decorated function is hit first opens the context; subsequent decorated calls
in the same run detect the active context and call through without
re-entering. This gives us the "log captures the full run" property even when
significant orchestration (argument normalization, validation, annotation/ref-panel
loading) lives in the outer function, while still letting direct class-method
callers from the Python API get a log file of their own.

```
ldsc ldscore / run_ldscore(**kwargs)
  → run_ldscore_from_args()        [DECORATED @workflow_run("ldscore")] ← opens context
    → LDScoreCalculator.run()      [DECORATED @workflow_run("ldscore")] ← guarded; calls through

(direct Python API)
LDScoreCalculator().run(...)       [DECORATED @workflow_run("ldscore")] ← opens context

ldsc build-ref-panel / run_build_ref_panel(**kwargs)
  → run_build_ref_panel_from_args()[DECORATED @workflow_run("build-ref-panel")] ← opens context
    → ReferencePanelBuilder.run()  [DECORATED @workflow_run("build-ref-panel")] ← guarded; calls through

(direct Python API)
ReferencePanelBuilder().run(...)   [DECORATED @workflow_run("build-ref-panel")] ← opens context

ldsc munge-sumstats
  → run_munge_sumstats_from_args() [DECORATED @workflow_run("munge-sumstats")]
                                     (CLI path — maps args to MungeConfig, then calls SumstatsMunger.run)

SumstatsMunger.run()               [DECORATED @workflow_run("munge-sumstats")]
                                     (Python API path and shared CLI implementation)

ldsc h2 → run_h2_from_args()                          [DECORATED @workflow_run("h2")]
ldsc rg → run_rg_from_args()                          [DECORATED @workflow_run("rg")]
ldsc partitioned-h2 → run_partitioned_h2_from_args()  [DECORATED @workflow_run("partitioned-h2")]
```

**Why this is safe:** there is exactly one log file per run. The outermost
decorated call enters the context, and all inner decorated calls call through.
The inner decorator's `get_output_dir`/`get_overwrite` extractors are only
consulted when the inner method is the top-level entry point (direct Python
API usage of the class method). When invoked through the outer
`run_*_from_args`, the inner extractor is never called.

---

## Munge-Sumstats Structural Alignment

`sumstats_munger.main()` now follows the package pattern: it parses `argv` and
delegates to `run_munge_sumstats_from_args(args)`.

The remaining harmonization change is to decorate the parsed-args function and
`SumstatsMunger.run()` with the shared workflow logging context once that
context exists. `main()` is already:

```python
def main(argv=None):
    return run_munge_sumstats_from_args(build_parser().parse_args(argv))
```

This makes the decorator work correctly (it sees `args.output_dir` before the
function body runs) and aligns munge with the rest of the package.

---

## Kernel Cleanup (`_kernel/sumstats_munger.py`)

The legacy `Logger` class, `sec_to_str` function, `START_TIME`, `log =
Logger(...)`, `log.close()`, the MASTHEAD+Call header, and the "Conversion
finished" footer have already been removed. The `log.*()` call sites are now
normal `LOGGER.*()` calls.

Float formatting in the kernel uses a module-level lambda to avoid cross-layer imports:
```python
_fmt = lambda v: format(v, ".4g")
```

---

## Scope Summary

| File | Action |
|---|---|
| `src/ldsc/_logging.py` | **Create** |
| `src/ldsc/_kernel/sumstats_munger.py` | Legacy Logger infrastructure already removed; later pass only needs float/message harmonization if desired |
| `src/ldsc/sumstats_munger.py` | `run_munge_sumstats_from_args` and slim `main` already exist; later pass decorates it + `SumstatsMunger.run` and switches `sumstats.log` to the shared context |
| `src/ldsc/ldscore_calculator.py` | Decorate **both** `run_ldscore_from_args` and `LDScoreCalculator.run` with `@workflow_run("ldscore")` |
| `src/ldsc/regression_runner.py` | Decorate 3 entry functions; add `--log-level`; use `format_float` |
| `src/ldsc/ref_panel_builder.py` | Decorate **both** `run_build_ref_panel_from_args` and `ReferencePanelBuilder.run` with `@workflow_run("build-ref-panel")` |
| `annotation_builder.py`, `_kernel/annotation.py` | **Not modified** (annotate refactored separately) |
| `cli.py` | **Not modified** |

---

## Known Limitations (Out of Scope)

### Output bundle atomicity
The decorator's log-file existence check guards the common case (re-running an
existing workflow without `--overwrite`): the run aborts before any file is
touched. However, if the log file has been manually deleted while scientific
outputs remain, the decorator opens the log file before the function body's
`ensure_output_paths_available` check fires on the surviving outputs, and the
log file may be left behind from the failed run. This is an unusual,
manually-induced state; treating it as out of scope is acceptable.

### Duplicate error logging through `cli.py`
`cli.run_cli()` catches user-facing exceptions and logs them at `ERROR`.
After this plan, kernel error sites also call `LOGGER.error(...)` (the design
doc fix to the `LOGGER.debug()` suppression hack). Errors that propagate up to
`run_cli()` may therefore appear twice in the console output. A proper fix
belongs in a separate task that revisits CLI exception handling; this plan
does not modify `cli.py`.
