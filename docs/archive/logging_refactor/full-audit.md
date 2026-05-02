# Logging and Traceback Full Audit

Date: 2026-04-29

Scope: `src/ldsc` in `ldsc_py3_restructured`.

This audit follows the repository layer contract documented in
`docs/architecture.md` and `docs/layer-structure.md`: the CLI layer owns public
command dispatch and traceback routing; workflow modules own user-facing path,
schema, config, genome-build, and output-preflight behavior; `_kernel` modules
receive resolved primitive inputs and preserve legacy parsing or numerical
behavior; output modules own persisted artifact layout.

## 1. Current Logging Setup Summary

The refactor now has a package exception hierarchy in `src/ldsc/errors.py` and a
central CLI error boundary in `ldsc.cli.run_cli()`. `python -m ldsc` routes
through that boundary, `LDSCUserError` and common user-caused built-ins are
rendered as clean `ERROR: Error: ...` messages without traceback, and unexpected
internal failures are logged with `LOGGER.exception(...)` and exit code 2.

Operational logging now follows the module-level uppercase `LOGGER` convention.
Most workflow modules use stable logger names such as `LDSC.cli`,
`LDSC.ldscore_calculator`, `LDSC.sumstats_munger`,
`LDSC.regression_runner`, and `LDSC.ref_panel_builder`; the stale
`LDSC.new` logger has been removed. Public workflow entry points now log
start/completion context for LD-score calculation, summary-statistics munging,
regression, reference-panel building, config banners, genome-build inference,
and LD-score directory loading.

The summary-statistics munger layering debt has since been reduced:
`sumstats_munger.main()` now parses only, `run_munge_sumstats_from_args()` maps
CLI args into `MungeConfig`, and `SumstatsMunger.run()` owns `sumstats.log`,
metadata sidecars, path resolution, and result objects. The kernel no longer
owns a hand-rolled log writer; it emits normal package logger records while
retaining legacy-compatible QC and `.sumstats.gz` writing. The remaining
logging debt is workflow-level harmonization across subcommands: common log
file naming, headers, timing, command recording, and float formatting. One
standalone compatibility entry point in `_kernel.annotation` still raises
`SystemExit(main_bed_to_annot())` because it is a script-style boundary, not a
workflow-layer API.

## 2. Completed Since Previous Audit

| Previous ID | Status | Current implementation |
| --- | --- | --- |
| A-1, T-1, T-2 | Done | Added `LDSCError`, `LDSCUserError`, `LDSCUsageError`, `LDSCConfigError`, `LDSCInputError`, `LDSCDependencyError`, and `LDSCInternalError`; `run_cli()` centralizes user/internal error routing. |
| A-2 | Done | `ConfigMismatchError` now inherits from both `LDSCConfigError` and `ValueError`, preserving compatibility while allowing package-level user-error catches. |
| A-3, E-2 | Done | `ldsc annotate` argument validation raises `LDSCUsageError` instead of direct `SystemExit` inside dispatch. |
| A-4 | Done | Optional dependency failures in annotation, LD-score, and reference-panel parquet/liftover paths now raise `LDSCDependencyError`, which also subclasses `ImportError` for legacy compatibility. |
| L-1 | Done | `_kernel.ldscore` logger now uses `LDSC.ldscore`. |
| L-2 through L-10, L-12 through L-13 | Done | Standard logger calls in workflow modules and LD-score/reference-panel kernels now use f-strings. |
| M-1 through M-10 | Mostly done | CLI, LD-score, summary-statistics, regression, and LD-score directory loading now emit useful workflow start/completion logs. |
| P-1 | Done | Global config banners now use `LOGGER.info(...)`; suppression behavior is preserved. |
| P-3 through P-5 | Done | Legacy munger progress no longer writes progress dots/status directly to stdout. |
| P-6 | Done | IRWLS shape failure now reports expected and observed shapes in the exception. |
| E-9 | Done | `partitioned-h2` now tells users to rerun `ldsc ldscore` with query annotations and explicit baseline annotations. |
| T-3 | Done | Invalid manifest config provenance is logged at `DEBUG` with `exc_info=True` and warned with the original exception value. |
| T-4 | Superseded | Munger conversion failures now propagate to the CLI boundary without an extra kernel-owned error log; exception chains are preserved for callers. |
| T-5, T-6 | Done | Intentional parquet schema probe fallbacks are narrowed/documented as schema fallback logic. |

## 3. Residual Log Message Findings

Severity: High means the issue can hide failures or materially confuse users;
Medium means important operational context or consistency is missing; Low means
style, clarity, or maintainability.

| ID | File | Current pattern | Issue | Recommendation | Severity |
| --- | --- | --- | --- | --- | --- |
| L-1 | `src/ldsc/_kernel/ref_panel.py` | Broad fallback logs BIM-only metadata at `DEBUG` | The fallback is well logged with traceback, but the broad catch is still wider than the expected PLINK BED reader failures. | If the PLINK reader exposes stable exception types later, narrow this catch. | Low |
| L-2 | `src/ldsc/regression_runner.py`, `src/ldsc/sumstats_munger.py` | Provenance warnings use `warnings.warn(...)` plus selective debug logging | These are API warnings rather than operational logs, which is appropriate for Python callers but less visible in CLI logs. | Keep for API compatibility; optionally mirror to logger at `WARNING` only at CLI boundaries. | Low |
| L-3 | workflow modules | Per-workflow log-file setup is inconsistent | `munge-sumstats` still writes `sumstats.log` through a workflow-owned temporary handler, while other workflows do not yet share a common log-file context. | Implement the later logging harmonization pass with shared handlers, timing, command recording, and naming. | Medium |

## 4. Residual Error Message Findings

Track: U means user-facing; D means developer-facing.

| ID | File | Track | Current pattern | Issue | Recommendation | Severity |
| --- | --- | --- | --- | --- | --- | --- |
| E-1 | `src/ldsc/config.py` | U | `ValueError("Required path value is missing.")` | The helper does not name the field whose path is missing. | Add an optional `label` argument to `_normalize_required_path()` when a future config pass touches path validation. | Medium |
| E-2 | `src/ldsc/path_resolution.py` | U | Exact-one path resolution errors | Messages give token and match count but only some call sites explain accepted token forms. | For top-level CLI paths, include accepted forms: exact path, glob, or explicit `@` chromosome suite when supported. | Low |
| E-3 | `src/ldsc/regression_runner.py` | U | Empty regression merge | Now names identifier mode, sumstats source, LD-score row count, and the likely config mismatch. | No further action. | Done |
| E-4 | `src/ldsc/_kernel/sumstats_munger.py` | U | Missing signed statistic, required column, N, and allele errors | High-impact messages now include available columns or concrete fix hints for `--signed-sumstats`, sample-size flags, and `--no-alleles`. | Add targeted tests for these messages to prevent regression. | Low |
| E-5 | `src/ldsc/_kernel/ref_panel.py` | U | `ImportError("Chromosome discovery for parquet R2 requires explicit chromosomes or sidecar metadata...")` | This is not an optional dependency failure; `ImportError` is semantically misleading. | Convert to `LDSCInputError` or `ValueError` in a follow-up because the fix is to pass chromosomes/metadata, not install a package. | Medium |
| E-6 | `src/ldsc/_kernel/ref_panel.py` | U | `ImportError("ParquetR2RefPanel.load_metadata requires metadata sidecar files.")` | Same semantic issue as E-5. | Convert to `LDSCInputError` or `ValueError` in a follow-up and update tests if they assert `ImportError`. | Medium |

## 5. Exception Architecture Status

Current hierarchy:

```text
Exception
├── LDSCError
│   ├── LDSCUserError
│   │   ├── LDSCUsageError
│   │   ├── LDSCConfigError
│   │   │   └── ConfigMismatchError  (+ ValueError compatibility)
│   │   ├── LDSCInputError
│   │   └── LDSCDependencyError      (+ ImportError compatibility)
│   └── LDSCInternalError
├── ValueError
├── FileNotFoundError
├── FileExistsError
├── NotADirectoryError
└── SystemExit at true process boundaries
```

This hierarchy is now adequate for the current layered architecture. New
workflow-layer validation should prefer `LDSCUsageError`, `LDSCConfigError`, or
`LDSCInputError` when the error is meant to be shown cleanly at the CLI boundary.
Optional dependency failures should use `LDSCDependencyError`; it remains an
`ImportError` subclass so existing dependency-gated tests and callers can still
catch it as an import failure.

Do not introduce package exceptions named after Python built-ins such as
`SystemError`, `RuntimeError`, `ValueError`, `FileNotFoundError`,
`PermissionError`, or `ConnectionError`.

## 6. Traceback Enrichment Status

| ID | Boundary | Status | Notes |
| --- | --- | --- | --- |
| T-1 | `ldsc.__main__` | Done | `python -m ldsc` exits through `run_cli()`, not raw workflow propagation. |
| T-2 | `ldsc.cli.run_cli()` | Done | User errors produce clean messages; internal package and unexpected failures log traceback with `LOGGER.exception(...)`. |
| T-3 | Optional dependency paths | Done | Annotation, LD-score parquet/PLINK, reference-panel parquet, and liftover dependency failures now preserve exception chains with `from exc`. |
| T-4 | Manifest and metadata provenance | Done | Invalid LD-score provenance logs debug traceback and warns with the cause; invalid sumstats metadata raises from the original exception. |
| T-5 | Legacy munger conversion | Done | The `.log` artifact records the original conversion exception before re-raising, preserving the Python traceback chain. |
| T-6 | Intentional schema probes | Acceptable | Narrowed to `ValueError` and documented as fallback probes. |

## 7. Recommended Next Steps

### Priority 1: Keep covered

- [x] Keep `tests/test_logging_refactor.py` covering the error hierarchy and
      `run_cli()` user/internal split.
- [x] Keep dependency errors chained with `from exc`.
- [x] Keep raw tracebacks out of expected user-error paths.

### Priority 2: Follow-up refinements

- [ ] Convert the two semantic `ImportError` cases in `_kernel/ref_panel.py`
      to `LDSCInputError` or `ValueError`.
- [ ] Add focused tests for munger missing signed-statistic, required-column,
      sample-size, and allele messages.
- [ ] Add a `label` argument to `_normalize_required_path()` so missing config
      paths can name the field.
- [ ] Decide how much of the historical `sumstats.log` message text should
      remain stable before the shared logging context changes the file name,
      header, timing, and float formatting.

### Priority 3: Optional polish

- [ ] Mirror API-level `warnings.warn(...)` messages to the CLI logger only at
      command boundaries if CLI users are missing important provenance warnings.
- [ ] Narrow the broad PLINK BED fallback in `_kernel/ref_panel.py` if stable
      lower-level exception classes become available.

## 8. Verification Plan

Run these commands from the repository root after logging/error changes:

```bash
python -m unittest tests.test_logging_refactor tests.test_sumstats_munger tests.test_regression_workflow tests.test_annotation tests.test_ref_panel_builder tests.test_ldscore_workflow
python -m unittest discover -s tests -p 'test*.py'
python -m ldsc --help
```

Quick audit probes:

```bash
rg -n "LOGGER\.(info|warning|error|debug|exception)\([^f\n]*%|logger\.(info|warning|error|debug|exception)\([^f\n]*%" src/ldsc --glob '*.py'
rg -n "^\s*print\(|sys\.stdout\.write|sys\.stderr\.write" src/ldsc tests --glob '*.py'
rg -n "except Exception|except BaseException|except:\s*(#.*)?$" src/ldsc --glob '*.py'
```
