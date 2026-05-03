# Workflow Logging

Public workflow entry points share one logging policy:

- Per-run file handlers attach to the `LDSC` logger, so workflow and kernel
  child records are captured together.
- Known scientific outputs and the log path are preflighted before the log file
  is opened. If preflight fails, no new log file is created.
- If execution fails after preflight, the log is kept and ends with a `Failed`
  footer plus elapsed time.
- `--log-level` controls module records in console and file; lifecycle audit
  lines always appear in the file.
- Workflow result objects and `output_paths` mappings do not include log files.

## Log Layout

Workflow logs begin with a lifecycle banner, then a multi-line `Call:` block.
The executable is written first and each following option/value pair is written
on its own continuation line where possible:

```text
Call:
./munge_sumstats.py \
--genome-build auto \
--sumstats /path/raw.tsv \
--out /path/sumstats
```

`Inputs:` and `Outputs:` are audit sections, not machine-readable manifests.
`Outputs:` is separated from preceding workflow records by a blank line so
summary blocks such as `Metadata:` remain visually distinct.

The footer records final status and elapsed time with explicit units:

```text
Finished 2026-05-02 22:07:50
Elapsed time: 2.0min:12s
```

## Log Names

| Workflow | Log path |
| --- | --- |
| `munge-sumstats` | `<output_dir>/sumstats.log` |
| `annotate` | `<output_dir>/annotate.log` |
| `ldscore` | `<output_dir>/ldscore.log` |
| `build-ref-panel` | `<output_dir>/build-ref-panel.log` |
| `h2` | `<output_dir>/h2.log` |
| `partitioned-h2` | `<output_dir>/partitioned-h2.log` |
| `rg` | `<output_dir>/rg.log` |

Regression commands without `--output-dir` stay console-only and do not create
log files.

## API Boundary

Python convenience wrappers that delegate through parsed workflow functions may
create logs because they use the same output-directory contract as the CLI.
Direct computational class APIs remain data-oriented:

- `AnnotationBuilder.project_bed_annotations(...)` does not create
  `annotate.log` when called directly.
- `LDScoreCalculator.run(...)` does not create `ldscore.log`.
- `ReferencePanelBuilder.run(...)` does not create `build-ref-panel.log`.
- `RegressionRunner.estimate_*` methods do not create regression logs.

`SumstatsMunger.run(...)` remains a workflow-level API and keeps writing
`sumstats.log` for compatibility, but `MungeRunSummary.output_paths` excludes
the log path.

For the implementation rationale, see
`docs/superpowers/specs/2026-05-02-logging-harmonization-design.md`.
