# Logging and Traceback Full Audit

Date: 2026-04-28

Scope: `src/ldsc` in `ldsc_py3_restructured`.

This audit follows the repository layer contract documented in `docs/architecture.md`
and `docs/layer-structure.md`: the CLI dispatch layer should expose one
user-facing command surface; workflow modules should own user-facing path,
schema, config, genome-build, and output-preflight behavior; `_kernel` modules
should receive resolved primitive inputs and perform legacy parsing or numerical
work; output modules should own persisted artifact layout.

## 1. Logging Setup Summary

Current setup: module-level `LOGGER = logging.getLogger(...)` exists in several
workflow and kernel modules, including `ldsc.ref_panel_builder`,
`ldsc.ldscore_calculator`, `ldsc.sumstats_munger`, `ldsc.column_inference`,
`ldsc.chromosome_inference`, `ldsc._kernel.annotation`,
`ldsc._kernel.identifiers`, `ldsc._kernel.ref_panel`, and
`ldsc._kernel.ldscore`. Logger naming is uppercase `LOGGER`, but logger names
mix stable names such as `LDSC.ref_panel_builder` with stale names such as
`LDSC.new`.

Logging configuration is not centralized. `logging.basicConfig(...)` is called
from `ReferencePanelBuilder._configure_logging()` and from standalone kernel
compatibility entry points in `_kernel.annotation` and `_kernel.ldscore`.
Most public CLI paths rely on argparse and raw exception propagation rather than
a single top-level logging and traceback boundary.

Key gap: the codebase does not distinguish user-facing errors from internal
failures at the CLI boundary. Most validation errors currently reach end users
as raw Python tracebacks, while operational progress sometimes goes to stdout
through `print()` or `sys.stdout.write()` rather than the package logger.

## 2. Log Message Findings

Severity: High means the issue can hide failures or materially confuse users;
Medium means important operational context or consistency is missing; Low means
style, clarity, or maintainability.

| ID | File | Line | Current message or pattern | Issue | Suggested fix | Severity |
| --- | --- | ---: | --- | --- | --- | --- |
| L-1 | `src/ldsc/_kernel/ldscore.py` | 235 | `LOGGER = logging.getLogger("LDSC.new")` | Stale logger name does not identify the module or layer. It makes filtered logs and test assertions harder to interpret. | Rename to `logging.getLogger("LDSC.ldscore")` or `logging.getLogger("LDSC._kernel.ldscore")`; update any tests/docs that refer to the old logger. | Medium |
| L-2 | `src/ldsc/ref_panel_builder.py` | 230 | `LOGGER.info("Loaded %d SNP restriction identifiers in %s mode from %s.", ...)` | Uses percent-style logger arguments, inconsistent with the requested f-string standard. | Build a message variable: `message = f"Loaded {len(restriction_values)} SNP restriction identifiers in {restriction_mode} mode from '{restriction_path}'."; LOGGER.info(message)`. | Medium |
| L-3 | `src/ldsc/ref_panel_builder.py` | 250 | `LOGGER.warning("No usable liftover chain was provided for %s -> %s; ...", ...)` | Percent-style logger arguments; long message would be clearer as a named message variable. | Convert to an f-string message that names source and target builds. | Medium |
| L-4 | `src/ldsc/ref_panel_builder.py` | 305 | `LOGGER.info("Building reference-panel artifacts for chromosome %s from %s.", chrom, prefix)` | Percent-style logger arguments; path is unquoted, making spaces or glob-resolved prefixes harder to read. | `LOGGER.info(f"Building reference-panel artifacts for chromosome {chrom} from '{prefix}'.")`. | Medium |
| L-5 | `src/ldsc/ref_panel_builder.py` | 325 | `LOGGER.info("Skipping chromosome %s because no SNPs remain after restriction.", chrom)` | Percent-style logger arguments; lacks retained/restricted count context. | Include counts when available, e.g. skipped chromosome, restriction mode, and total candidate SNPs. | Medium |
| L-6 | `src/ldsc/ref_panel_builder.py` | 447 | `LOGGER.info("Finished chromosome %s with %d retained SNPs.", chrom, len(metadata))` | Percent-style logger arguments; otherwise good completion context. | Convert to f-string: `LOGGER.info(f"Finished chromosome {chrom} with {len(metadata)} retained SNPs.")`. | Low |
| L-7 | `src/ldsc/ldscore_calculator.py` | 863 | `LOGGER.info("Resolved LD-score annotation genome build from %s.", sampled_path)` | Percent-style logger arguments; inferred build is logged by lower layer, but this message does not repeat the resolved build. | Include path and build in one message when possible. | Medium |
| L-8 | `src/ldsc/ldscore_calculator.py` | 882 | `LOGGER.info("Resolved LD-score reference-panel genome build from %s.", sampled_path)` | Percent-style logger arguments; message omits resolved build. | Include sampled path and resolved build. | Medium |
| L-9 | `src/ldsc/_kernel/identifiers.py` | 209 | `logger.warning("Dropped %d restriction rows with missing CHR or POS from %s.", dropped, path)` | Percent-style logger arguments. Message is operationally useful but should match project style. | Convert to f-string with quoted path. | Low |
| L-10 | `src/ldsc/_kernel/annotation.py` | 756 | `LOGGER.info("Resolved annotation genome build from %s.", sampled_path)` | Percent-style logger arguments; inferred build is not visible in this message. | Include sampled path and inferred build, or rely on a single lower-level inference message to avoid duplicates. | Medium |
| L-11 | `src/ldsc/_kernel/ldscore.py` | 1158 | `LOGGER.warning("Cannot apply --maf in %s because MAF metadata is unavailable.", context)` | Percent-style logger arguments; message is good but can include what behavior follows. | Convert to f-string and state that no MAF filter was applied for the context. | Medium |
| L-12 | `src/ldsc/_kernel/ldscore.py` | 1305 | Raw-schema parquet warning | Uses percent-style formatting and logs a severe performance warning from the kernel. The warning is important but should include the path and the remediation command context consistently. | Convert to f-string; mention that row-group pruning is disabled and recommend regenerating with `ldsc build-ref-panel`. | Medium |
| L-13 | `src/ldsc/_kernel/ldscore.py` | 1365 | Missing build metadata warning | Uses percent-style formatting and duplicates inference behavior. | Convert to a single f-string warning with path, inferred build, and regeneration hint. | Medium |
| L-14 | `src/ldsc/_kernel/ldscore.py` | 1950 | `LOGGER.warning("Skipping %s.l2.M_5_50 because MAF is unavailable.", prefix)` | Percent-style logger arguments; legacy artifact name is clear. | Convert to f-string. | Low |

### Missing Logging Coverage

These functions are long enough or operational enough that additional entry and
completion logging would help users diagnose runs without reading tracebacks.
Pure utilities and numerical kernels are intentionally excluded unless they
perform file I/O or public workflow orchestration.

| ID | Module / function | Lines | Reason flagged | Recommended entry log | Recommended exit log |
| --- | --- | ---: | --- | --- | --- |
| M-1 | `src/ldsc/cli.py::main` | 74 | The top-level public command dispatcher has no logging or error boundary. | Log selected command at DEBUG after parse, or INFO for actual workflow start through subcommands. | Log clean failure through centralized handler; do not print raw traceback for user errors. |
| M-2 | `src/ldsc/ldscore_calculator.py::run_ldscore_from_args` | 611 | Public workflow entry point resolves many inputs, builds configs, and writes outputs. | Inputs that matter: baseline/query token presence, reference backend, output directory, identifier mode, genome build. | Number of chromosomes, output path, baseline/query column counts. |
| M-3 | `src/ldsc/ldscore_calculator.py::LDScoreCalculator.run` | 229 | Main LD-score orchestration loop currently logs only indirect events. | Total chromosome count and annotation columns. | Number of completed/skipped chromosomes and output directory if written. |
| M-4 | `src/ldsc/ldscore_calculator.py::compute_chromosome` | 305 | Per-chromosome compute boundary delegates to expensive PLINK/parquet kernels. | Chromosome, backend, annotation row count, and LD window mode. | Retained row count and count-column summary. |
| M-5 | `src/ldsc/sumstats_munger.py::SumstatsMunger.run` | 231 | Public munging workflow resolves inputs and writes three fixed outputs but delegates most logging to legacy kernel stdout/log file. | Raw path, output directory, identifier mode, genome build, and overwrite state. | Input and retained row counts, output artifact paths. |
| M-6 | `src/ldsc/sumstats_munger.py::main` | 404 | CLI entry point performs path resolution and config inference, but only legacy kernel writes progress. | Raw path and output directory after argument normalization. | Retained row count and metadata sidecar path. |
| M-7 | `src/ldsc/regression_runner.py::run_h2_from_args` | 496 | Public regression CLI path loads artifacts and may write summary output without logger messages. | Sumstats path, LD-score directory, output directory. | Number of regression SNPs and output summary path if written. |
| M-8 | `src/ldsc/regression_runner.py::run_partitioned_h2_from_args` | 515 | Public regression CLI path loops over query annotations with no operational logs. | Sumstats path, LD-score directory, query annotation count. | Number of query summaries and output summary path if written. |
| M-9 | `src/ldsc/regression_runner.py::run_rg_from_args` | 535 | Public rg CLI path loads two sumstats artifacts and LD scores with no operational logs. | Trait paths/names and LD-score directory. | Number of merged SNPs and output summary path if written. |
| M-10 | `src/ldsc/regression_runner.py::load_ldscore_from_dir` | 606 | Public loader reads manifest and parquet artifacts but only emits warnings for provenance problems. | Directory path and manifest path at DEBUG or INFO. | Baseline row count, query row count, and provenance status. |

### `print()` and stdout Misuse

| ID | File | Line | Current call | Issue | Suggested replacement |
| --- | --- | ---: | --- | --- | --- |
| P-1 | `src/ldsc/config.py` | 192 | `print(f"{entrypoint} using {global_config!r}")` | Operational banner is written to stdout for Python and CLI workflows. This can pollute command output and scripted pipelines. | Replace with logger-backed banner. Keep suppression context for repeated nested calls. |
| P-2 | `src/ldsc/_kernel/sumstats_munger.py` | 80 | `print(msg, file=self.log_fh)` | File logging is intentional for legacy `.log` artifact; the class is not integrated with Python logging. | Keep file write but rename/adapt as `LegacyMungeLogWriter`; pair stdout mirroring with package logger or a verbosity flag. |
| P-3 | `src/ldsc/_kernel/sumstats_munger.py` | 82 | `print(msg)` | Legacy progress and diagnostics go directly to stdout. | Mirror through `LOGGER.info` or adapter so CLI boundary controls formatting and verbosity. |
| P-4 | `src/ldsc/_kernel/sumstats_munger.py` | 261 | `sys.stdout.write('.')` | Progress dots are not structured, cannot be filtered, and can interfere with shell pipelines. | Remove dots or emit DEBUG progress by chunk number. |
| P-5 | `src/ldsc/_kernel/sumstats_munger.py` | 325 | `sys.stdout.write(' done\n')` | Completes the progress-dot stream through stdout instead of a logger. | Replace with a structured completion log containing chunk and SNP counts. |
| P-6 | `src/ldsc/_kernel/_irwls.py` | 116 | `print('IRWLS update:', new_w.shape, w.shape)` | Debug print inside numerical kernel. It triggers only before raising but bypasses logging. | Include shapes in the raised `ValueError` and remove the print. |

## 3. Error Message Findings

Track: U means user-facing; D means developer-facing.

| ID | File | Line | Track | Current raise | Issue | Suggested replacement | Severity |
| --- | --- | ---: | --- | --- | --- | --- | --- |
| E-1 | `src/ldsc/cli.py` | 137 | U | `raise ValueError(f"Unsupported command: {args.command}")` | Should be impossible after argparse, but if reached it produces a raw traceback. | Raise `LDSCUsageError` or call `parser.error(...)` through the CLI boundary. | Medium |
| E-2 | `src/ldsc/cli.py` | 143 | U | `raise SystemExit("ldsc annotate requires --query-annot-bed-paths and --baseline-annot-paths.")` | Uses `SystemExit` directly inside dispatch. Message is clear, but it bypasses any structured error handler. | Raise a user/config error and let the CLI boundary render `Error: ...`. | Medium |
| E-3 | `src/ldsc/config.py` | 60 | U | `ValueError("Required path value is missing.")` | Generic message does not name the config field. | Include field/context at caller or change helper to accept `label`. | Medium |
| E-4 | `src/ldsc/path_resolution.py` | 123 | U | Exact-one path resolution error | Message gives count and token, but not examples of accepted token forms. | Include accepted forms where user-facing: exact path, glob, or explicit `@` suite when enabled. | Low |
| E-5 | `src/ldsc/path_resolution.py` | 179 | U | `FileNotFoundError(f"Could not resolve {label} token: {token}") from None` | `from None` intentionally suppresses lower-level token traceback, which is acceptable for user errors, but there is no domain exception. | Convert at workflow/CLI boundary to `LDSCInputError` or leave built-in but catch cleanly at boundary. | Medium |
| E-6 | `src/ldsc/ldscore_calculator.py` | 661 | U | Query requires baseline message | Strong actionable message; keep content. It should be routed cleanly without traceback at CLI. | Reclassify as `LDSCUsageError` or catch `ValueError` as user error at boundary. | Low |
| E-7 | `src/ldsc/ldscore_calculator.py` | 837 | U | Missing genome build for chr_pos | Actionable; mentions accepted fixes. | Keep content; route through clean CLI boundary. | Low |
| E-8 | `src/ldsc/regression_runner.py` | 139 | U | `No overlapping chr_pos SNPs remain after merging sumstats with LD scores.` | Good high-level cause, but lacks path/trait context. | Include sumstats source path or trait and LD-score directory when raised at workflow boundary. | Medium |
| E-9 | `src/ldsc/regression_runner.py` | 522 | U | `partitioned-h2 requires query annotations in --ldscore-dir.` | Actionable but can include how to create them. | Add hint: rerun `ldsc ldscore` with `--query-annot-paths` or `--query-annot-bed-paths` and explicit baseline annotations. | Medium |
| E-10 | `src/ldsc/regression_runner.py` | 640 | U | `LD-score directory does not exist or is not a directory: {root}` | Clear, user-facing. | Keep; route without traceback. | Low |
| E-11 | `src/ldsc/_kernel/sumstats_munger.py` | 382 | D | `Cannot determine N. This message indicates a bug...` | Internal bug message reaches users from a kernel path if pre-validation missed it. | Raise `LDSCInternalError` or add workflow pre-validation so this cannot be user-visible. | Medium |
| E-12 | `src/ldsc/_kernel/sumstats_munger.py` | 779 | U | `Could not find {C} column.` | Names missing internal canonical column but not the input header or accepted aliases. | Include available columns and accepted aliases for the missing field. | High |
| E-13 | `src/ldsc/_kernel/sumstats_munger.py` | 795 | U | `Could not determine N.` | Does not tell user which flags or columns can fix it. | Include required fixes: provide `--N`, `--N-cas/--N-con`, or an inferable sample-size column. | High |
| E-14 | `src/ldsc/_kernel/sumstats_munger.py` | 803 | U | `Could not find A1/A2 columns.` | Does not mention `--no-alleles` as a valid h2-only workaround. | Add hint: provide allele column flags or pass `--no-alleles` when alleles are not needed. | High |
| E-15 | `src/ldsc/_kernel/ldscore.py` | 1977 | U | `Specify exactly one reference-panel mode: parquet or PLINK.` | Clear but should name concrete flags. | Use: `Specify exactly one reference-panel mode: --r2-paths/--r2-table for parquet or --plink-path/--bfile for PLINK.` | Medium |
| E-16 | `src/ldsc/_kernel/ldscore.py` | 1188 | U | `Must specify exactly one of --ld-wind-snps, --ld-wind-kb, or --ld-wind-cm.` | Good user-facing validation. | Keep; route without traceback. | Low |

## 4. Exception Architecture Findings

Current hierarchy:

```text
Exception
├── ValueError
│   └── ConfigMismatchError
├── FileNotFoundError
├── FileExistsError
├── NotADirectoryError
├── ImportError
└── SystemExit
```

Evaluation: the hierarchy is insufficient for the current layered architecture.
The codebase has many user-caused validation and input errors, plus a smaller
set of internal or dependency failures, but callers cannot distinguish those
tracks by exception class. `ConfigMismatchError` is useful and semantic, but it
inherits from `ValueError` directly and does not establish a package-wide base.

Recommended hierarchy:

```python
class LDSCError(Exception):
    """Base class for package-specific failures."""


class LDSCUserError(LDSCError):
    """User input, CLI arguments, or configuration are invalid."""


class LDSCConfigError(LDSCUserError):
    """Configuration values or cross-artifact config provenance are invalid."""


class LDSCInputError(LDSCUserError):
    """Input paths, file schemas, or file contents are invalid."""


class LDSCDependencyError(LDSCUserError):
    """An optional dependency required for the requested workflow is missing."""


class LDSCInternalError(LDSCError):
    """Unexpected internal state; log with traceback at the boundary."""
```

Do not introduce names that shadow Python built-ins such as `SystemError`,
`RuntimeError`, `ValueError`, `FileNotFoundError`, `PermissionError`, or
`ConnectionError`.

| ID | Issue | Location | Description | Recommendation |
| --- | --- | --- | --- | --- |
| A-1 | No package base exception | `src/ldsc` | Most public workflow errors use built-ins directly. | Add `src/ldsc/errors.py`; export package exceptions from `ldsc.__init__`. |
| A-2 | Config mismatch not rooted in package base | `src/ldsc/config.py:44` | `ConfigMismatchError` inherits from `ValueError`, so boundary code cannot catch all package user errors by base class. | Make it inherit from `LDSCConfigError` while preserving `ValueError` compatibility if needed through multiple inheritance. |
| A-3 | Direct `SystemExit` in workflow dispatch | `src/ldsc/cli.py:143` | A specific subcommand validation exits directly instead of raising semantic user error. | Raise `LDSCUsageError` or `LDSCInputError`; central CLI boundary exits. |
| A-4 | Optional dependency failures use plain `ImportError` | `_kernel/annotation.py:769`, `_kernel/ldscore.py:282`, `_kernel/ref_panel_builder.py:601` | User sees dependency errors as Python import failures. | Wrap or raise `LDSCDependencyError` with install hint at workflow boundary. |
| A-5 | Intentional probe exceptions use `pass` | `column_inference.py:437`, `_kernel/ldscore.py:758`, `_kernel/ref_panel.py:307` | Silent catches are intentional schema probes but are indistinguishable from accidental swallowing. | Replace with named helper functions or comments such as `# Try build-specific POS aliases after generic POS fails.` |

## 5. Traceback Enrichment Findings

| ID | File | Line | Issue | Suggested fix |
| --- | --- | ---: | --- | --- |
| T-1 | `src/ldsc/__main__.py` | 6 | `python -m ldsc` calls `main()` directly, so all non-argparse errors become raw tracebacks. | Introduce `cli.run_cli(argv=None) -> int` with user/internal split; call it from `__main__`. |
| T-2 | `src/ldsc/cli.py` | 74 | Unified CLI dispatch has no `try/except` boundary around workflow calls. | Catch `LDSCUserError`, selected built-ins (`ValueError`, `FileNotFoundError`, `FileExistsError`, `NotADirectoryError`, `ImportError`) and render clean messages; catch `LDSCError` and unexpected `Exception` with `LOGGER.exception`. |
| T-3 | `src/ldsc/regression_runner.py` | 696 | `except Exception:` while rebuilding manifest config emits a warning and returns `None`; original exception type and value are discarded. | Catch `Exception as exc`; warning should include `exc`; optionally `LOGGER.debug(..., exc_info=True)` for developer trace. |
| T-4 | `src/ldsc/_kernel/sumstats_munger.py` | 889 | Broad catch logs a generic error banner and re-raises. Chain is preserved because it re-raises, but the log omits the exception message. | Log `ERROR converting summary statistics: {exc}` and re-raise. Remove unused `traceback` import. |
| T-5 | `src/ldsc/_kernel/ldscore.py` | 758 | Broad catch suppresses schema-probe failure. | Keep suppression only if intentional; narrow to `ValueError` or add comment explaining canonical-schema probe fallback. |
| T-6 | `src/ldsc/_kernel/ldscore.py` | 764 | Broad catch suppresses schema-probe failure. | Narrow to expected schema-resolution exceptions and document fallback to raw/unsupported schema. |

## 6. Recommended Changes

### Priority 1: High severity

- [ ] Add `src/ldsc/errors.py` with package exception classes that distinguish
      user/config/input/dependency failures from internal failures.
- [ ] Add a central CLI error boundary in `ldsc.cli` and route `python -m ldsc`
      through it from `ldsc.__main__`.
- [ ] Convert the most common user-facing validation errors in CLI and public
      workflow entry points to package user exceptions, or catch the existing
      built-ins cleanly during the first pass.
- [ ] Improve high-impact munger messages for missing signed statistic, sample
      size, and allele columns.

### Priority 2: Medium severity

- [ ] Replace `print_global_config_banner()` stdout output with logger-backed
      behavior while preserving `suppress_global_config_banner()`.
- [ ] Convert percent-style logger calls to f-strings in `ref_panel_builder.py`,
      `ldscore_calculator.py`, `_kernel/identifiers.py`, `_kernel/annotation.py`,
      and `_kernel/ldscore.py`.
- [ ] Rename `_kernel.ldscore` logger from `LDSC.new` to a stable module name.
- [ ] Add entry/completion logs to the public workflow entry points:
      `run_ldscore_from_args`, `SumstatsMunger.run`, regression `run_*_from_args`,
      and `load_ldscore_from_dir`.
- [ ] Replace operational `warnings.warn(...)` calls with logger warnings where
      the condition is runtime progress/status rather than Python API warning.

### Priority 3: Low severity and consistency

- [ ] Remove the IRWLS debug print and include shape details in the exception.
- [ ] Remove unused `traceback` import from `_kernel/sumstats_munger.py`.
- [ ] Document intentional schema-probe exception suppression with comments or
      helper names.
- [ ] Review log levels after conversion: expected milestones at INFO,
      recoverable data-quality drops at WARNING, internal traces at DEBUG.

## 7. Suggested Implementation Sequence

1. Add package exceptions and CLI boundary first. This reduces traceback noise
   immediately without changing numerical behavior.
2. Move global config banners from stdout to logging. This is low-risk but
   affects shell-script output, so update tests that assert stdout/stderr.
3. Convert logging style and logger names. This should be mostly mechanical.
4. Improve user-facing error messages in munger, LD-score, and regression
   workflow boundaries. Keep kernel behavior stable and avoid changing filters.
5. Clean up legacy stdout progress in the munger kernel behind the public
   workflow adapter.
6. Run targeted CLI tests and then the full test suite.

## 8. Verification Plan

Run these commands from the repository root after implementation:

```bash
python -m unittest discover -s tests -p 'test*.py' -v
python -m ldsc --help
python -m ldsc munge-sumstats --help
python -m ldsc ldscore --help
python -m ldsc h2 --help
```

Add targeted tests for:

- CLI user errors exit nonzero without raw tracebacks.
- Unexpected internal errors still log full traceback with `LOGGER.exception`.
- `print_global_config_banner()` no longer writes to stdout in normal workflows.
- Existing legacy munger `.log` artifact still contains conversion details.

## 9. Agent Execution Instructions

This section is written for an AI agent executing the cleanup plan above.

### Before You Start

- Read this full audit and the active design docs:
  `docs/architecture.md`, `docs/layer-structure.md`, and `AGENTS.md`.
- Do not change numerical behavior, filtering semantics, file formats, or public
  CLI flag names as part of this logging pass.
- Keep edits in the documented layers: CLI boundary in `ldsc.cli`, user-facing
  workflow context in public workflow modules, and only local message/logging
  improvements inside `_kernel`.

### Execution Order

1. Implement `src/ldsc/errors.py` and export exceptions.
2. Add the top-level CLI boundary and route `__main__.py` through it.
3. Move global config banners to logging.
4. Convert logger message style and stale logger names.
5. Improve high-impact user-facing error text.
6. Clean stdout/progress misuse in legacy kernels where tests allow it.
7. Run targeted and full verification.

### Per-Change Checklist

- For every `except` block you touch, preserve exception chains with `raise ...
  from exc` unless intentionally suppressing user-facing noise with `from None`.
- For user errors, include what was wrong, what was received, what was expected,
  and a concrete fix hint when practical.
- For developer/internal errors, include operation context and preserve the
  original exception.
- For log messages, prefer f-strings and include path/count/build/chromosome
  context when relevant.
- After each file edit, run the smallest targeted test that exercises it before
  the full suite.

Awaiting approval before applying code changes.
