# Logging Harmonization Design

**Date:** 2026-05-02
**Implementation status:** implemented in the current logging harmonization change set.
**Documentation status:** current docs, tutorials, and relevant workflow
docstrings updated with the final logging contract.
**Scope:** `src/ldsc/_logging.py`, `src/ldsc/sumstats_munger.py`,
`src/ldsc/annotation_builder.py`, `src/ldsc/ldscore_calculator.py`,
`src/ldsc/ref_panel_builder.py`, `src/ldsc/regression_runner.py`

---

## Problem

LDSC workflows had inconsistent logging behavior:

- `munge-sumstats` wrote `sumstats.log` through a private kernel-only handler.
- `annotate` and `build-ref-panel` used local `basicConfig(...)` helpers.
- `ldscore` and regression workflows logged to console but did not create
  per-run logs.
- Regression subcommands did not expose `--log-level` consistently.
- `MungeRunSummary.output_paths` exposed the log file path even though other
  workflow result contracts only expose scientific output artifacts.

The remaining problem was not layering. The public workflow modules now own
orchestration, so logging can be harmonized at those workflow boundaries without
decorating direct class APIs.

---

## Design Decisions

| Question | Decision |
|---|---|
| Abstraction shape | Add a small explicit `workflow_logging(...)` context manager in `src/ldsc/_logging.py`; do not introduce workflow decorators. |
| Handler target | Attach the per-run file handler to `logging.getLogger("LDSC")` so workflow and kernel child records are captured together. |
| File format | Message-only records. Lifecycle audit lines are file-only and are written directly by the context. |
| Invocation format | Write `Call:` followed by the executable and one shell-quoted option/value pair per continuation line where possible. |
| Output section spacing | Write a blank line before `Outputs:` so preceding metadata or QC summaries remain visually separate. |
| Elapsed-time format | Write footer duration as `Elapsed time: <minutes>.0min:<seconds>s`, for example `Elapsed time: 2.0min:12s`. |
| Log level | `--log-level` controls the LDSC logger threshold for module records in console and file. Lifecycle audit lines are always present in the file. |
| CLI boundary | `cli.run_cli()` remains the user-facing exception boundary. Workflows should not log an error and then re-raise only to duplicate CLI errors. |
| Python convenience wrappers | Wrappers that delegate through `run_*_from_args(...)` may create logs because they use the same workflow boundary. |
| Direct class APIs | Direct class methods such as `LDScoreCalculator.run(...)` and `ReferencePanelBuilder.run(...)` do not gain new per-run log file behavior. |
| Munge compatibility | `munge-sumstats` keeps `<output_dir>/sumstats.log`. |
| Output contracts | No workflow result `output_paths` mapping includes log files. `MungeRunSummary.output_paths` no longer has a `"log"` key. |
| Scientific formatting | Scientific output formatting is out of scope. `format_float(...)` is for log text only. |

---

## Log Filenames

| Workflow | Log path |
|---|---|
| `munge-sumstats` | `<output_dir>/sumstats.log` |
| `annotate` | `<output_dir>/annotate.log` |
| `ldscore` | `<output_dir>/ldscore.log` |
| `build-ref-panel` | `<output_dir>/build-ref-panel.log` |
| `h2` | `<output_dir>/h2.log` |
| `partitioned-h2` | `<output_dir>/partitioned-h2.log` |
| `rg` | `<output_dir>/rg.log` |

Regression workflows without `--output-dir` remain console-only and do not
create a log file.

---

## Shared Logging API

`src/ldsc/_logging.py` provides:

- `workflow_logging(workflow_name, log_path, log_level="INFO")`
- `log_inputs(**items)`
- `log_outputs(**items)`
- `format_float(value) -> str`
- `configure_package_logging(level="INFO")`

`workflow_logging(...)` saves the existing LDSC logger level, applies the
requested level, optionally attaches a message-only `FileHandler`, writes a
start header with a multi-line `Call:` block, and writes either a `Finished` or
`Failed` footer with elapsed time. It always restores the LDSC logger level and
removes its handler on exit.

Example header and footer:

```text
Call:
./munge_sumstats.py \
--genome-build auto \
--sumstats /path/raw.tsv \
--out /path/sumstats

Finished 2026-05-02 22:07:50
Elapsed time: 2.0min:12s
```

Nested contexts are treated as no-op wrappers around the active outer context.
This keeps helper code from accidentally replacing the active workflow log.

---

## Preflight And Failure Guarantees

Workflow modules must compute all known final output paths, including the log
path, before opening the log file.

Guarantees:

- If preflight fails, no new log file is created.
- If execution fails after preflight, the log file is kept.
- Escaping exceptions write a `Failed` footer and elapsed time.
- The workflow context does not write the exception message into the footer, so
  the CLI boundary can report the user-facing error once.

---

## Result Contract

`output_paths` means scientific or curated data outputs, not audit logs.

For `munge-sumstats`, `MungeRunSummary.output_paths` may contain:

- `sumstats_parquet`
- `sumstats_gz`
- `metadata_json`

It must not contain `log`.

`sumstats.metadata.json["output_files"]` remains limited to curated sumstats
data artifacts. It does not record `sumstats.log`.

The same rule applies outside munge: `LDScoreResult.output_paths` records
`manifest`, `baseline`, and optional `query`; `ReferencePanelBuildResult`
records emitted parquet and metadata artifacts; regression return values do not
gain log-path fields.

---

## Annotation Inclusion

`annotation_builder.py` is now a public workflow layer, so annotation is included
in harmonization. Its workflow writes `annotate.log` when `output_dir` is
present and uses the shared LDSC logging helper instead of a local
`basicConfig(...)` setup.

---

## Deferred Work

- Polish legacy message text that embeds severity words such as `WARNING:`
  inside warning log messages.
- Broaden structured input/output logging where it helps debugging.
- Revisit whether direct API methods should support opt-in log files.
