# Workflow Logging

Last updated on: 2026-09-14

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
- LD-score BED/gene-list runs log the effective catalog projection build, one warning for every non-`ok` query, and a final query-status summary. Gene-list runs additionally log every rejected or zero-support row with role/source. Intentional MHC exclusions use one INFO line per gene-list source and role, listing the count and physical line–gene pairs in input order; a submitted alias includes its canonical gene ID only when different. Audit chunks do not create extra summary lines. The compressed audit remains the complete machine-readable row record.

## Console vs File Routing

All command help pages reuse `LOG_LEVEL_HELP` from [`ldsc._logging`](../../src/ldsc/_logging.py). The default is `INFO`: progress and summaries. `DEBUG` adds troubleshooting details; `WARNING` records warnings and errors; `ERROR` records errors only. Each level includes more severe records. Run headers and completion status remain in the file at every level. This help wording describes the workflow log; it does not promise terminal progress at `INFO` or `DEBUG`.

The per-run `.log` file is the authoritative sink. Console output (stderr) is a
CLI-only concern, installed by `ldsc.cli.run_cli`; `stdout` is reserved for
explicit report output such as `munge-sumstats --infer-only` and the successful munging summary, so log records
never go to `stdout`. The Python API never writes to the console.

| Context | Output dir | Module records (INFO/DEBUG) | Errors |
| --- | --- | --- | --- |
| CLI | provided | `.log` file only | full traceback to `.log`; concise line + logfile pointer to console |
| CLI | derived | `plot` and `convert-h2-scale` write below their source result; other materializing commands require `--output-dir` | full traceback to the derived workflow log; concise line + logfile pointer to console |
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

After a successful `munge-sumstats` CLI run, stdout displays a `Munge-sumstats summary:` block with the selected restriction, liftover method, mapping/drop counts and reasons, and whole-run row accounting. The workflow writes the identical block to `diagnostics/sumstats.log` through `log_summary`, which bypasses the module-record threshold like lifecycle audit lines. This summary remains visible at every `--log-level`; Python API runs write it only to the log. See [liftover rules and count interpretation](munge-sumstats.md#liftover-rules).

## Output-Family Preflight

Public workflow layers treat their fixed current-contract outputs plus workflow
log as one owned family. Without `--overwrite`, any existing owned artifact
rejects the run before the log is opened. With `--overwrite`, stale owned
siblings that the successful run did not produce are removed after the current
outputs are written. Legacy root log names such as `sumstats.log`,
`annotate.log`, `ldscore.log`, `build-r2-panel*.log`, or regression logs at the
output root are not part of the current owned family. For `build-r2-panel`,
concrete chromosome PLINK prefixes own only the matching chromosome-scoped log,
while `@` chromosome-suite prefixes own the full build-r2-panel log family.

This keeps an output directory from mixing artifacts from different
configurations, while preserving unrelated user files. Direct Python writer
APIs apply the same rule to their data artifacts. Standalone `run_annotate()` owns the canonical annotation log; `AnnotationBuilder.run()` prepares private data without installing a workflow log.

The six directory writers in `ldsc.outputs` each expose `artifact_family()`, returning an `ArtifactFamily` with selected output paths and the complete owned scope. H2, partitioned-h2, rg, quantile-h2, query-r2, and LD-score workflow preflights consult their writer's declaration and add workflow-owned logs or audits. The LDSC2 converter and indexed LD-score CLI also consult the LD-score writer. Declarations only describe paths; the existing `path_resolution.preflight_output_artifact_family()` still performs collision checks.

Early checks discard the returned stale list because the final scientific and diagnostic outputs are not yet known. Each writer derives its final declaration from the result or write options, uses those same selected paths for output and metadata file entries, and removes stale owned siblings after its existing publication step. LD-score diagnostics-only writes use this same declaration; existing `diagnostics/dropped_snps/chr*_dropped.tsv.gz` reports are included in early collision checks and final reconciliation. There is no second workflow cleanup based on an earlier prediction.

Annotation, munging, and reference-panel build workflows retain their local ownership declarations and existing cleanup stages. Reference-panel builds keep their chromosome scope, and gene-index publication retains its dedicated transaction. Logs and workflow-only audit extensions do not enter scientific metadata through `ArtifactFamily`. Failure markers remain owned by the marker helper. Validation lives in `tests/test_artifact_declarations.py`, `tests/test_output.py`, `tests/test_derived_output_lifecycle.py`, `tests/test_failure_markers.py`, and the reference-panel/index workflow tests.

## Annotation preparation

Shared source preparation emits three INFO milestones in the workflow log: reading annotation inputs, checking SNP identities and preparing chromosome annotations, and successful completion. For example:

```text
Reading annotation inputs: baseline files=22, query files=22.
Checking SNP identities and preparing chromosome annotations.
Annotation preparation complete: chromosomes=22, retained SNPs=1,000,000, elapsed=30.00s.
```

The counts and time above are illustrative. The first line precedes source scanning; the second precedes global identity validation, retained-row selection, and shard writing. Completion appears only after usable annotation shards are ready, with the number of retained chromosomes, logical SNP rows after identity cleanup and any chromosome selection, and elapsed preparation time. Aligned baseline/query files describe the same logical rows and do not double the SNP count. For the gene-index builder these messages precede `Starting chromosome N`; preparation remains serial.

These are step boundaries, without per-chunk messages, periodic heartbeats, percentages, or an estimated finish time. A CM/MAF compatibility notice appears only on the first chunk of each annotation file read; later chunks do not repeat it, and later preparation calls still receive it. At `WARNING` or `ERROR`, INFO milestones and notices are suppressed by the existing log-level policy. Preparation failures retain the last reached milestone and the workflow's existing `Failed` footer and traceback, without a preparation-complete message. Preparation completion does not mean the enclosing workflow or output publication has finished.

Intentional gene exclusions use the following compact form:

```text
Genes intentionally excluded by region policy: role=focal source=pathway.txt count=2 line:gene=[31:ENSG00000196126, 33:ENSG00000196735]
```

Mechanisms: [`_annotation_sources.py`](../../src/ldsc/_annotation_sources.py), `prepare_annotation_sources()`; [`_annotation_parsing.py`](../../src/ldsc/_annotation_parsing.py), `normalize_annotation_chunk()`; [`query_annotations.py`](../../src/ldsc/query_annotations.py), `_log_gene_list_rejections()`.

## LD-score chromosome diagnostics

Direct query validation logs the validated baseline/reference chromosome sets, effective scope, and ordinary-glob selection caveat. The calculation log explicitly names chromosomes resolved and entering analysis. Indexed runs name their immutable validated scope. `diagnostics/chromosome_scope.json` persists this evidence; its `analysis_chromosomes` is empty when scope validation prevents analysis. The LD-score writer owns this file and `diagnostics/input_issues.tsv` along with existing query/gene diagnostics, so preflight, replacement, and diagnostics-only cleanup share one declaration. Scope is also carried on `LDScoreResult.chromosome_scope` and in root scientific metadata.

## Failed Overwrite Markers

Every public materializing CLI workflow and corresponding high-level Python
workflow wraps its existing action order with the marker helper in
`ldsc._logging`. If an attempt authorized with `--overwrite` or
`overwrite=True` fails, the helper writes a durable marker after the exception
escapes:

- ordinary commands: `<output_dir>/RUN_FAILED.txt`
- `plot`: `<result-dir>/plots/RUN_FAILED.txt`
- `convert-h2-scale`:
  `<h2-result-dir>/postprocessing/liability-scale/RUN_FAILED.txt`
- concrete chromosome `build-r2-panel`:
  `<output-dir>/RUN_FAILED.chr<chrom>.txt`

The marker records the failed command/API boundary, UTC timestamp, exception,
detailed log path or absence of one, and a conservative warning that the active
directory may contain incomplete or mixed artifacts. It is not a scientific
result and is not listed in metadata. Marker handling does not move, restore,
roll back, or quarantine files and does not change when a workflow writes its
ordinary artifacts or removes stale outputs. If the workflow log opened before
the failure, it remains with the usual `Failed` footer and traceback. A
successful materializing retry removes its applicable marker after normal
success. No-overwrite failures create no marker.

The marker helper uses `path_resolution.normalize_path_token`, matching output preflight and publication. CLI marker scope follows the last occurrence of a destination option, as argparse does. Marker-only directories do not collide with writer artifact families. Gene-index preflight, direct publication, and interrupted-publication recovery share recognition of owned diagnostics; this excludes symlinks, special files, and unrelated empty subdirectories. Existing recovery of a validated owned backup preserves the failure marker until the current build succeeds. See the [output-directory audit](../audits/2026-09-14-output-directory-retries.md).

CLI help and `munge-sumstats --infer-only` bypass the marker lifecycle entirely. Even with `--overwrite`, they create no output directories or files and never create, replace, or remove a failure marker. This also applies when inference or argument parsing fails.

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
| `build-r2-panel` | `<output_dir>/diagnostics/build-r2-panel.log`, or `<output_dir>/diagnostics/build-r2-panel.chr<chrom>.log` for concrete single-chromosome PLINK-prefix runs |
| `build-gene-ldscore-index` | completed success: `<index_dir>/diagnostics/build-gene-ldscore-index.log`; running/failed: `<parent>/.<index-name>.build-state/build-gene-ldscore-index.log`; prior failed attempts move to hidden `history/` |
| `convert-ldsc2-ldscores` | `<output_dir>/diagnostics/convert-ldsc2-ldscores.log` |
| `h2` | `<output_dir>/diagnostics/h2.log` |
| `partitioned-h2` | `<output_dir>/diagnostics/partitioned-h2.log` |
| `quantile-h2` | `<output_dir>/diagnostics/quantile-h2.log` |
| `rg` | `<output_dir>/diagnostics/rg.log` |
| `query-r2` | `<output_dir>/diagnostics/query-r2.log` |
| `plot` | `<result-dir>/plots/diagnostics/plot.log` |
| `convert-h2-scale` | `<h2-result-dir>/postprocessing/liability-scale/diagnostics/convert-h2-scale.log` |

Every regression CLI command requires `--output-dir` and writes its command log
under that directory's `diagnostics/` tree.

`convert-h2-scale` is post-processing rather than a regression fit and derives
its destination from an existing h2 result. `plot` does the same for any
supported result suite.

LD-score logs list binary and quantitative fitted annotations. Partitioned-h2 logs repeat an actionable interpretation warning when quantitative annotations are present: legacy numerical proportion/enrichment summaries remain visible, but only coefficient-based fields retain their ordinary interpretation for those annotations. Quantile-h2 logs the selected fitted model, target, inherited common-MAF rule, common reference-SNP universe size, missing exclusions, and realized quantile bounds/counts. Row-addressable alignment issues are written to `snp_alignment_issues.tsv.gz` rather than expanded into the log.

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

- `AnnotationBuilder.run(..., output_dir=...)` requires an output parent for private preparation but does not install a workflow log. `run_annotate(...)` owns standalone output preflight, writing, and `diagnostics/annotate.log`.
- `LDScoreCalculator.run(...)` does not create `diagnostics/ldscore.log`.
- `ReferencePanelBuilder.run(...)` does not create a build-r2-panel workflow log.
- `RegressionRunner.estimate_*` methods do not create regression logs.

`SumstatsMunger.run(...)` remains a workflow-level API and keeps writing
`diagnostics/sumstats.log`. `MungeRunSummary.output_paths` excludes the log
path but includes data artifacts such as the dropped-SNP audit sidecar.

For the implementation rationale, see
`docs/specs/2026-05-02-logging-harmonization-design.md` and the
console/file routing change in
`docs/plans/2026-06-04-logging-console-file-routing-plan.md`.

## Streaming memory-workflow diagnostics

Annotation, gene-resolution, direct chromosome-drop, and quantile alignment audits are replayed in bounded chunks into their complete output files. Results keep counts/statuses or persistent paths; they do not retain every chromosome's audit frame. Standalone annotate publishes query shards incrementally and detaches its returned bundle from construction scratch. Batch partitioned-h2 stages complete per-fit details privately and publishes them after final summary sorting, preserving overwrite and failure-marker contracts. See [the memory design](annotation-memory-design.md).
