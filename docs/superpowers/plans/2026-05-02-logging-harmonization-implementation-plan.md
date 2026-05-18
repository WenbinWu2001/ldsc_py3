# Logging Harmonization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:executing-plans` or `superpowers:test-driven-development` to
> implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for
> tracking.

**Goal:** Harmonize per-run logging across public LDSC workflows with a small
shared context manager, while keeping workflow result contracts limited to data
artifacts.

**Architecture:** Add `src/ldsc/_logging.py` with an explicit
`workflow_logging(...)` context. Workflow entry points compute and preflight
scientific outputs plus the log path before entering the context. Direct class
APIs do not gain new log-file behavior.

**Tech Stack:** Python 3.11, standard `logging`, `unittest`

**Implementation status:** completed in the current logging harmonization
change set.

---

## File Map

| File | Change |
|---|---|
| `src/ldsc/_logging.py` | New shared workflow logging context and helpers |
| `src/ldsc/sumstats_munger.py` | Use shared context for `sumstats.log`; add `--log-level`; remove `"log"` from `MungeRunSummary.output_paths` |
| `src/ldsc/annotation_builder.py` | Include annotation in harmonization; write `annotate.log`; replace local `basicConfig(...)` |
| `src/ldsc/ldscore_calculator.py` | Preflight final outputs plus `ldscore.log`; wrap workflow execution |
| `src/ldsc/ref_panel_builder.py` | Preflight final outputs plus `build-ref-panel.log`; keep direct class API log-free |
| `src/ldsc/regression_runner.py` | Add `--log-level`; write `h2.log`, `partitioned-h2.log`, and `rg.log` when `output_dir` is present |
| `tests/test_logging_refactor.py` | Unit tests for shared logging behavior |
| Workflow tests | Assert expected log creation and that `output_paths` mappings do not include logs |
| `README.md`, `docs/current/*.md`, `tutorials/*` | Describe log filenames, preflight behavior, and the rule that result paths exclude logs |
| Workflow module docstrings | Document which APIs create logs and which direct class APIs stay log-file free |

---

## Task 1: Add Shared Logging Module

- [x] Add `src/ldsc/_logging.py`.
- [x] Implement `workflow_logging(workflow_name, log_path, log_level="INFO")`.
- [x] Attach a message-only `FileHandler` to `logging.getLogger("LDSC")`.
- [x] Write lifecycle audit lines directly to the file so start/end/failure
  footers appear even at `ERROR` level.
- [x] Restore logger level and remove the file handler on exit.
- [x] Add `log_inputs(**items)`, `log_outputs(**items)`, and
  `format_float(value)`.
- [x] Add cleanup if file setup fails during `__enter__`.

## Task 2: Test Shared Logging Behavior

- [x] Verify the LDSC ancestor handler captures child logger records.
- [x] Verify lifecycle audit lines appear at `ERROR` level while INFO module
  records are suppressed.
- [x] Verify escaping exceptions write a `Failed` footer without duplicating
  the exception message.
- [x] Verify levels and handlers are restored after context exit.

## Task 3: Update `munge-sumstats`

- [x] Add `--log-level`.
- [x] Keep writing `<output_dir>/sumstats.log`.
- [x] Replace the munge-specific kernel file handler with `workflow_logging`.
- [x] Preflight selected sumstats outputs, metadata, and `sumstats.log` before
  opening the log.
- [x] Preserve parquet-first output behavior.
- [x] Remove `"log"` from `MungeRunSummary.output_paths`.
- [x] Keep `sumstats.metadata.json["output_files"]` limited to curated sumstats
  data artifacts.

## Task 4: Update Other Workflows

- [x] `annotate`: write `annotate.log` when `output_dir` is present and use the
  shared LDSC logger configuration.
- [x] `ldscore`: preflight outputs plus `ldscore.log` and wrap calculation.
- [x] `build-ref-panel`: preflight outputs plus `build-ref-panel.log` and keep
  direct class API calls log-free.
- [x] `h2`: add `--log-level`, write `h2.log` only when `output_dir` is present.
- [x] `partitioned-h2`: add `--log-level`, write `partitioned-h2.log` only when
  `output_dir` is present.
- [x] `rg`: add `--log-level`, write `rg.log` only when `output_dir` is present.

## Task 5: Update Workflow Tests

- [x] `munge-sumstats` creates `sumstats.log`, but summary output paths exclude
  `log`.
- [x] `ldscore`, `annotate`, `build-ref-panel`, `h2`, `partitioned-h2`, and
  `rg` create their expected log files when `output_dir` is present.
- [x] Regression without `output_dir` creates no log file.
- [x] Non-munge result `output_paths` mappings do not gain log keys.

## Task 6: Documentation

- [x] Add `docs/superpowers/specs/2026-05-02-logging-harmonization-design.md`.
- [x] Add `docs/superpowers/plans/2026-05-02-logging-harmonization-implementation-plan.md`.
- [x] Leave archived logging docs unchanged.
- [x] Update source docstrings for workflow logging boundaries and result
  contracts.
- [x] Update active current docs and README.
- [x] Update tutorials so collision notes and expected outputs include logs
  without treating logs as scientific result paths.

---

## Verification Commands

Run the targeted suite:

```bash
/Users/wenbinwu/miniforge3/envs/ldsc-py3-Jerry/bin/python -m unittest \
  tests.test_logging_refactor \
  tests.test_sumstats_munger \
  tests.test_ldscore_workflow \
  tests.test_annotation \
  tests.test_ref_panel_builder \
  tests.test_regression_workflow \
  tests.test_package_layout
```

Run the full suite:

```bash
/Users/wenbinwu/miniforge3/envs/ldsc-py3-Jerry/bin/python -m unittest discover -s tests
```

Expected result: all tests pass. PyArrow may emit sandbox-related `sysctlbyname`
warnings on macOS; those warnings are not test failures.
