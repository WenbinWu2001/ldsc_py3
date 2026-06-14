# Workflow Logging

Public workflow entry points share one logging policy:

- Per-run file handlers attach to the `LDSC` logger, so workflow and kernel
  child records are captured together.
- Known scientific outputs and the log path are preflighted before the log file
  is opened. For workflows with coherent output families, this preflight covers
  owned siblings that are not produced by the current run. If preflight fails,
  no new log file is created.
- If execution fails after preflight, the log is kept and ends with a `Failed`
  footer, the full Python traceback, and elapsed time, so failures are
  diagnosable from the log file alone.
- `--log-level` controls module-record verbosity. With an output directory, those
  records go to the workflow log file (see routing below); lifecycle audit lines
  always appear in the file. Supported levels are `DEBUG`, `INFO`, `WARNING`, and
  `ERROR`.
- Workflow result objects and `output_paths` mappings do not include log files.

## Console vs File Routing

The per-run `.log` file is the authoritative sink. Console output (stderr) is a
CLI-only concern, installed by `ldsc.cli.run_cli`; `stdout` is reserved for
deliverable command output (e.g. regression TSV tables), so log records never go
to `stdout`. The Python API never writes to the console.

| Context | Output dir | Module records (INFO/DEBUG) | Errors |
| --- | --- | --- | --- |
| CLI | provided | `.log` file only | full traceback to `.log`; concise line + logfile pointer to console |
| CLI | none | console (stderr) | console |
| Python API | provided | `.log` file only | full traceback to `.log`; exception propagates to caller |
| Python API | none | nowhere | exception propagates to caller |

Mechanism: `run_cli` installs one stderr `StreamHandler` on the `LDSC` logger at
`ERROR` level for the duration of the command. When a run has no log file (a
console-only quick run), `workflow_logging` temporarily lowers that handler to the
run's level so progress is visible; with a log file it stays at `ERROR`, so
ordinary records go to the file and only errors echo to the console. The `LDSC`
logger keeps `propagate = True` and the root logger is never given a handler, so
nothing is duplicated to the console and `caplog`-based tests keep working.

## Output-Family Preflight

Public workflow layers treat their fixed current-contract outputs plus workflow
log as one owned family. Without `--overwrite`, any existing owned artifact
rejects the run before the log is opened. With `--overwrite`, stale owned
siblings that the successful run did not produce are removed after the current
outputs are written. Legacy root log names such as `sumstats.log`,
`annotate.log`, `ldscore.log`, `build-ref-panel*.log`, or regression logs at the
output root are not part of the current owned family. For `build-ref-panel`,
concrete chromosome PLINK prefixes own only the matching chromosome-scoped log,
while `@` chromosome-suite prefixes own the full build-ref-panel log family.

This keeps an output directory from mixing artifacts from different
configurations, while preserving unrelated user files. Direct Python writer
APIs apply the same rule to their data artifacts and omit workflow logs unless
they are reached through a wrapper that creates one.

## Log Layout

Workflow logs begin with a lifecycle banner, then a multi-line `Call:` block.
The executable is written first and each following option/value pair is written
on its own continuation line where possible:

```text
Call:
ldsc munge-sumstats \
  --raw-sumstats-file /path/raw.tsv \
  --output-dir /path/sumstats \
  --source-genome-build auto \
  --output-genome-build hg38
```

`Inputs:` and `Outputs:` are audit sections, not machine-readable manifests.
`Outputs:` is separated from preceding workflow records by a blank line so
summary blocks such as `Metadata:` remain visually distinct.

Liftover drop summaries in ordinary logs are count-only. Row-level drop audit
records belong in the workflow's dropped-SNP sidecar, and example SNPs are
emitted only when the logger is set to `DEBUG`.

The footer records final status and elapsed time with explicit units:

```text
Finished 2026-05-02 22:07:50
Elapsed time: 2.0min:12s
```

The `Started` timestamp and elapsed timer are captured as one paired timepoint
when the workflow context is entered. The `Finished` or `Failed` timestamp and
elapsed timer are likewise captured as one paired timepoint when the context
exits. This keeps the wall-clock timepoints and `Elapsed time` footer aligned
with the actual interval covered by the workflow log, including setup before the
header is written and final work before the footer is written.

## Log Names

| Workflow | Log path |
| --- | --- |
| `munge-sumstats` | `<output_dir>/diagnostics/sumstats.log` |
| `annotate` | `<output_dir>/diagnostics/annotate.log` |
| `ldscore` | `<output_dir>/diagnostics/ldscore.log` |
| `build-ref-panel` | `<output_dir>/diagnostics/build-ref-panel.log`, or `<output_dir>/diagnostics/build-ref-panel.chr<chrom>.log` for concrete single-chromosome PLINK-prefix runs |
| `h2` | `<output_dir>/diagnostics/h2.log` |
| `partitioned-h2` | `<output_dir>/diagnostics/partitioned-h2.log` |
| `rg` | `<output_dir>/diagnostics/rg.log` |

Regression commands without `--output-dir` stay console-only and do not create
log files; their progress records print to the console (stderr) via the routing
described above.

## API Boundary

Python convenience wrappers that delegate through parsed workflow functions may
create logs because they use the same output-directory contract as the CLI. The
console handler is installed only by `run_cli`, so direct API use never emits
console output: records go to the workflow log file when one is created, and
otherwise nowhere, while exceptions propagate to the caller unchanged. Direct
computational class APIs remain data-oriented:

- `AnnotationBuilder.project_bed_annotations(...)` does not create
  `diagnostics/annotate.log` when called directly.
- `LDScoreCalculator.run(...)` does not create `diagnostics/ldscore.log`.
- `ReferencePanelBuilder.run(...)` does not create a build-ref-panel workflow log.
- `RegressionRunner.estimate_*` methods do not create regression logs.

`SumstatsMunger.run(...)` remains a workflow-level API and keeps writing
`diagnostics/sumstats.log`. `MungeRunSummary.output_paths` excludes the log
path but includes data artifacts such as the dropped-SNP audit sidecar.

For the implementation rationale, see
`docs/superpowers/specs/2026-05-02-logging-harmonization-design.md` and the
console/file routing change in
`docs/superpowers/plans/2026-06-04-logging-console-file-routing-plan.md`.
