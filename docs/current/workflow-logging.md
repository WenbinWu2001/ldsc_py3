# Workflow Logging

Last updated on: 2026-09-07

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
- LD-score BED/gene-list runs log the effective catalog projection build, one
  warning for every non-`ok` query, and a final query-status summary. Gene-list
  runs additionally log every rejected or zero-support row with role/source and
  every intentional MHC exclusion; the compressed audit is the machine-readable
  record of record.

## Console vs File Routing

The per-run `.log` file is the authoritative sink. Console output (stderr) is a
CLI-only concern, installed by `ldsc.cli.run_cli`; `stdout` is reserved for
explicit report output such as `munge-sumstats --infer-only`, so log records
never go to `stdout`. The Python API never writes to the console.

| Context | Output dir | Module records (INFO/DEBUG) | Errors |
| --- | --- | --- | --- |
| CLI | provided | `.log` file only | full traceback to `.log`; concise line + logfile pointer to console |
| CLI | none | not a public workflow state; every command requires `--output-dir` | argument error |
| Python API | provided | `.log` file only | full traceback to `.log`; exception propagates to caller |
| Python API | none | nowhere | exception propagates to caller |

Mechanism: `run_cli` installs one stderr `StreamHandler` on the `LDSC` logger at
`ERROR` level for the duration of the command. With a workflow log it stays at
`ERROR`, so ordinary records go to the file and only errors echo to the console.
Low-level Python calls that intentionally omit a log do not install a console
handler. The `LDSC`
logger keeps `propagate = True` and the root logger is never given a handler, so
nothing is duplicated to the console and `caplog`-based tests keep working.

Gene-list `ldscore` adds two deliberately bounded, direct stderr notices after a
successful CLI run: one when `resolved-only` omitted rows and one when Gate B
found zero-SNP support or skipped/warning query outcomes. These notices bypass
the ordinary logger threshold so a successful SLURM job cannot hide a
science-relevant subset/skip in a log that users may never open. They name the
diagnostic paths and cap affected outcomes at 10. Python entry points do not
emit these console notices; callers inspect the returned statuses and paths.

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
APIs apply the same rule to their data artifacts. Public materializing workflow
methods, including `AnnotationBuilder.run()`, create their canonical log.

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
| `build-gene-ldscore-index` | completed success: `<index_dir>/diagnostics/build-gene-ldscore-index.log`; running/failed: `<parent>/.<index-name>.build-state/build-gene-ldscore-index.log`; prior failed attempts move to hidden `history/` |
| `convert-ldsc2-ldscores` | `<output_dir>/diagnostics/convert-ldsc2-ldscores.log` |
| `h2` | `<output_dir>/diagnostics/h2.log` |
| `partitioned-h2` | `<output_dir>/diagnostics/partitioned-h2.log` |
| `quantile-h2` | `<output_dir>/diagnostics/quantile-h2.log` |
| `rg` | `<output_dir>/diagnostics/rg.log` |
| `query-r2` | `<output_dir>/diagnostics/query-r2.log` |

Every regression CLI command requires `--output-dir` and writes its command log
under that directory's `diagnostics/` tree.

LD-score logs list binary and quantitative fitted annotations. Partitioned-h2 logs repeat an actionable interpretation warning when quantitative annotations are present: legacy numerical proportion/enrichment summaries remain visible, but only coefficient-based fields retain their ordinary interpretation for those annotations. Quantile-h2 logs the selected fitted model, target, inherited common-MAF rule, common reference-SNP universe size, missing exclusions, verification level, and realized quantile bounds/counts. Row-addressable alignment issues are written to `snp_alignment_issues.tsv.gz` rather than expanded into the log.

## Exact gene-index build log

The exact gene-index builder uses the same lifecycle banner, `Call:`, `Inputs:`,
`Outputs:`, and `Finished`/`Failed` footer as the other artifact-building
commands. Its stable path is created before chromosome work so it can be
monitored live. Its concise INFO narrative records resolved configuration and input
counts, the baseline/PLINK identifier intersection and configured regression-row universe, start and
completion for each chromosome, catalog genes after exclusion, retained
and regression rows, atom/operator nonzeros, component bytes, and publication
state. The JSON sidecar is the machine-readable summary; the log renders the
same per-chromosome and aggregate measurements rather than recomputing them.

For this builder, `Finished chromosome N` is emitted only after the worker has closed every chromosome payload and atomically installed the shard in the hidden run-specific stage. It means the internal shard is durable and its in-memory record can be released; it does not report partial public publication. The public destination changes only after all chromosome and shared metadata are finalized and the complete stage reloads successfully.

For a successful run, the final lines identify `index_id`, staged reload
validation, complete atomic replacement, payload bytes, peak RSS, and the
validated index path. After the `Finished` footer closes the handler, the
successful log moves into the published index's `diagnostics/`. A failure keeps the `Failed` footer and traceback at the hidden live path without creating or modifying an index destination; a failed overwrite leaves the old scientific index loadable. A retained interrupted stage is never a resumable checkpoint and is discarded on retry. Once replacement
and reload validation complete, cleanup failures are warnings that name the
retained transaction directory and do not change the successful exit status.
The lifecycle log is the status authority, so no separate status JSON is
written.

## API Boundary

Python convenience wrappers that delegate through parsed workflow functions may
create logs because they use the same output-directory contract as the CLI. The
console handler is installed only by `run_cli`, so direct API use never emits
console output: records go to the workflow log file when one is created, and
otherwise nowhere, while exceptions propagate to the caller unchanged. Direct
computational class APIs remain data-oriented:

- `AnnotationBuilder.run(...)` and
  `AnnotationBuilder.project_bed_annotations(...)` create
  `diagnostics/annotate.log` when their optional output destination is set;
  output-free calls remain in memory.
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
