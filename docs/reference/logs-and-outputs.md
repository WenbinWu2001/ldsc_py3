# Logs and Outputs

Last updated on: 2026-10-28

This page describes what a logged workflow writes besides its scientific results: logs, console messages, output preflight and failure markers.
The rules for writing log and error messages, and the logging API boundary, are in [logging design](../dev/design/logging.md).

## Read-only inspection exception

`ldsc infer-genome-build PATH` prints its report to stdout and actionable errors to stderr.
It never opens a workflow log, creates an output directory, or writes a failure marker, including on syntax errors.
It accepts no verbosity or logging options; errors must not suggest `--log-level DEBUG` or point to a previous workflow log.
An unknown build is a successful inspection (exit 0), not a logged workflow failure.
See the [inspection contract](snp-identity-and-genome-build.md#standalone-build-inspection) for the remaining exit codes.

## Workflow log lifecycle

Public workflow entry points share one logging policy:

- Per-run file handlers attach to the `LDSC` logger, so workflow and kernel
  child records are captured together.
- Known scientific outputs and the log path are preflighted before the log file
  is opened. For workflows with coherent output families, this preflight covers
  owned siblings that are not produced by the current run. If preflight fails,
  no new log file is created.
- If an exception escapes an opened workflow logging context, the log is kept with a `Failed` footer, the full Python traceback, and elapsed time.
  When the exception was raised in a calculation worker, the footer traceback ends with a `Calculation child traceback:` note that holds the worker's frames; the terminal still shows only the one-line error.
  Abrupt process termination can interrupt logging before a footer or traceback is written.
- `--log-level` controls module-record verbosity. With an output directory, those
  records go to the workflow log file (see routing below); lifecycle audit lines
  always appear in the file. Supported levels are `DEBUG`, `INFO`, `WARNING`, and
  `ERROR`.
- Workflow result objects and `output_paths` mappings do not include log files.
- INFO describes major phases, numbered batch start/completion, chromosome work, and aggregate outcomes. Successful query timings, source mappings, and full annotation lists belong at DEBUG.
- Gene exclusions, rejected rows, unsupported genes, and non-OK query outcomes use counts by role/reason and pointers to complete diagnostics. Exclusion counts describe submitted input rows, including repeats across lists; they are not counts of distinct genes. Failure messages retain the affected query/pair, stage, and traceback.
- A diagnostic pointer replaces detail only when the corresponding file is written.
  Warnings without an equivalent complete diagnostic, including raw summary-statistics QC warnings, remain self-contained.

Nested workflow logging contexts retain the active outer log rather than opening a second file or replacing its settings.
The owning context restores the previous logger level and removes its own handler on exit.

## Console vs File Routing

All command help pages reuse `LOG_LEVEL_HELP` from [`ldsc._runtime.workflow_logging`](../../src/ldsc/_runtime/workflow_logging.py).
The default is `INFO`: progress and summaries.
This help wording describes the workflow log; it does not promise terminal progress at `INFO` or `DEBUG`.

The per-run `.log` file is the authoritative workflow record.
The CLI installs its console handler on stderr; stdout is reserved for explicit reports such as `munge-sumstats --infer-only` and the successful munging summary.
Package-installed logging handlers do not send records to stdout.
Python APIs do not install a package console handler, but an application may configure its own handlers.

| Context | Output dir | Module records (INFO/DEBUG) | Errors |
| --- | --- | --- | --- |
| CLI, workflow log opened | provided | `.log` file | captured failure traceback to `.log`; concise terminal error, with a diagnostic/log pointer when supplied by the error message |
| CLI, workflow log opened | derived | `plot` and `convert-h2-scale` write below their source result; other materializing commands require `--output-dir` | same error routing as other logged CLI workflows |
| Python API, workflow log opened | provided | `.log` file; application handlers may also receive records | captured failure traceback to `.log`; exception propagates to caller |
| Python API, no workflow log | depends on API | no workflow file handler; application logging configuration controls other destinations | exception propagates to caller |

During a normal CLI invocation, `run_cli` installs one stderr `StreamHandler` on the `LDSC` logger at `ERROR` level for the duration of the command.
With a workflow log it stays at `ERROR`, so ordinary records go to the file and errors can reach the terminal.
CLI entry points that activate logging without a workflow file retain their existing verbosity fallback.
Errors before a workflow log opens cannot provide a traceback from that unopened file.
Current internal-error reports add the saved logfile pointer when available; expected user-error reports use the context and pointers supplied by their exception message.

The package does not install a root logger handler, and records from the `LDSC` logger may propagate to handlers installed by an application.
Claims of “file only” above describe the package's own routing, not suppression of application logging.
Calling `run_cli(argv=...)` programmatically does not install the CLI console handler.

After a successful CLI run, `ldscore` and `annotate` print deliberately bounded notices directly to stderr.
One appears when `--allow-unresolved-genes` (`resolved-only`) omitted gene-list rows, one when gene-list SNP-support checks found unsupported genes or skipped or warning queries, and one when BED query annotations were skipped: `WARNING: <k> of <n> BED query annotation(s) were skipped and are omitted from <annotation|LD-score> outputs: ...`.
These notices bypass the ordinary logger threshold so that a successful batch job cannot hide a science-relevant omission in a log that users may never open.
They name the diagnostic paths and cap the listed outcomes at 10.
In `annotate`, the SNP-support notice is headed `Gene-list SNP-support checks` and says skipped queries are omitted from annotation outputs.
In `ldscore` and `build-gene-ldscore-index`, a further notice appears when reference-panel SNPs not on the annotation SNP grid were excluded: `WARNING: <command> excluded <n> reference-panel SNPs on <k> chromosomes that are not on the annotation SNP grid; ...`.
It has no threshold, repeats the WARNING already written to the workflow log, and names `diagnostics/reference_snp_funnel.tsv`, or `diagnostics/build-gene-ldscore-index.json` for an index build ([reference SNPs on the annotation SNP grid](ldscore-snp-universes.md#reference-snps-on-the-annotation-snp-grid)).
Python entry points print no console notices; callers inspect the returned statuses and paths.

In `ldscore`, the successful gene-list notice distinguishes unsupported gene-list rows from queries skipped after LD-score computation.
It reports an unsupported-row count only when that count is positive, explains each SNP-support or zero-variance outcome in words, and states that skipped queries are omitted from LD-score outputs.
For example, a notice can report three queries skipped because their LD scores are constant on the regression SNPs (`zero_variance_ld_scores`), even when every evaluated gene has retained reference-SNP support.
SNP support is measured on the retained reference SNPs; LD-score variation is evaluated on the regression SNPs.
The reason codes, diagnostic schemas, query-selection policy, and console/file routing remain unchanged.

After a successful `munge-sumstats` CLI run, stdout displays a `Munge-sumstats summary:` block with the selected restriction, liftover method, mapping/drop counts and reasons, and whole-run row accounting. The workflow writes the identical block to `diagnostics/sumstats.log` through `log_summary`, which bypasses the module-record threshold like lifecycle audit lines. This summary remains visible at every `--log-level`; Python API runs write it only to the log. See [liftover rules and count interpretation](munge-sumstats.md#liftover-rules).

## Output-Family Preflight

Public workflow layers treat their fixed current-contract outputs plus workflow log as one owned family.
Without `--overwrite`, any existing owned artifact rejects the run before the log is opened.
With `--overwrite`, stale owned siblings that the successful run did not produce are removed after the current outputs are written.
Legacy root log names such as `sumstats.log`, `annotate.log`, `ldscore.log`, `build-r2-panel*.log`, or regression logs at the output root are not part of the current owned family.
For `build-r2-panel`, concrete chromosome PLINK prefixes own only the matching chromosome-scoped log, while `@` chromosome-suite prefixes own the full `build-r2-panel` log family.

This keeps an output directory from mixing artifacts from different
configurations, while preserving unrelated user files. Direct Python writer
APIs apply the same rule to their data artifacts. Standalone `run_annotate()` owns the canonical annotation log; `AnnotationBuilder.run()` prepares private data without installing a workflow log.

The directory writers in the domain `outputs` modules (`ldsc._annotation.outputs`, `ldsc._refpanel.outputs`, `ldsc._ldscore.outputs` and `ldsc._regression.outputs`) each expose `artifact_family()`, returning an `ArtifactFamily` with selected output paths and the complete owned scope.
`h2`, `partitioned-h2`, `rg`, `quantile-h2`, `query-r2`, and LD-score workflow preflights consult their writer's declaration and add workflow-owned logs or audits.
The legacy LDSC converter and indexed LD-score CLI also consult the LD-score writer.
Declarations only describe paths; the existing `_runtime.paths.preflight_output_artifact_family()` still performs collision checks.

Early checks discard the returned stale list because the final scientific and diagnostic outputs are not yet known. Each writer derives its final declaration from the result or write options, uses those same selected paths for output and metadata file entries, and removes stale owned siblings after its existing publication step. LD-score diagnostics-only writes use this same declaration; existing `diagnostics/dropped_snps/chr*_dropped.tsv.gz` reports are included in early collision checks and final reconciliation. There is no second workflow cleanup based on an earlier prediction.

Annotation, munging, and reference-panel build workflows retain their local ownership declarations and existing cleanup stages.
Reference-panel builds keep their chromosome scope, and gene LD-score index publication retains its dedicated transaction.
Logs and workflow-only audit extensions do not enter scientific metadata through `ArtifactFamily`.
Failure markers remain owned by the marker helper.
Validation lives in `tests/test_artifact_declarations.py`, `tests/test_output.py`, `tests/test_derived_output_lifecycle.py`, `tests/test_failure_markers.py`, and the reference-panel/index workflow tests.

## Staged input validation

Output ownership and input validity are separate gates.
Existing artifact-family checks authorize logs/diagnostics and prevent replacement without overwrite.
Input gates do not authorize scientific output replacement, remove stale science files, or imply that publication succeeded.
R² panel builds derive their chromosome-specific output family from validated PLINK contents; their early discovery repair table can therefore precede the full output-family check, as the previous PLINK repair table did.

[`_runtime.preflight.InputGate`](../../src/ldsc/_runtime/preflight.py) collects independent input errors across objects at a stage and retains `input_issues` on the exception.
Writing workflows persist the existing six-column `input_issues.tsv` repair table under diagnostics; gene LD-score index builds use the sibling `.<index>.build-state/` directory.
Existing PLINK and gene-specific audits remain available.
The table is a diagnostic, not a new scientific artifact or resumable checkpoint.

Input-preflight errors name the failed input, report the issue count and an actionable repair, and point to the written `diagnostics/input_issues.tsv`.
A missing raw summary-statistics path identifies `--raw-sumstats-file`, reports that path resolution found no matching file, and states that summary-statistics reading has not started.
The primary message shows that pathname once; the chained developer traceback retains the original exception and its full context.
The diagnostic table retains its six-column schema and preserves specialized nested diagnostics.
Other captured operation errors distinguish missing and unreadable inputs and record exception types and operating-system error numbers when available.
One captured failure becomes the error's direct cause; multiple captured failures remain available through an `ExceptionGroup` cause.
The declaration-validation phase reports `failed_inputs` as the number of distinct role/source/chromosome combinations with issues.
When a phase total is unknown, progress omits the denominator; any reported completed count uses `completed=N`.

Direct LD scores check declarations and annotation headers, resolve gene-list controls, and validate reference content before annotation staging.
PLINK companion mappings and reference scope evidence are reused; annotation headers supply the source scan widths.
Gene LD-score index and R² panel construction check baseline/auxiliary declarations and PLINK companions together.
Index consumers check every component path, small chromosome identity metadata, and Parquet header before opening the first numerical operator.
See [`test_staged_input_preflight.py`](../../tests/test_staged_input_preflight.py).

Validation remains staged: BIM/FAM rows and BED dimensions require reference inspection; annotation row alignment, duplicate identities, chromosome contents, and gene support require their content scans.
R² panel query-scope validation streams pair integrity before preparation and reuses that scope result.
Baseline-only R² panel runs check sidecar/reader compatibility early and retain pair-content checks in the consuming reader, avoiding an additional full pair scan.
Ordinary globs cannot reveal a chromosome that is absent consistently from all matches.
Scientific omission, ignored-option, overwrite, and failure-marker policies remain unchanged.

## Annotation preparation

Shared source preparation emits INFO milestones in the workflow log for reading annotation inputs, checking SNP identities and preparing chromosome annotations, and successful completion. A storage summary immediately before completion records binary/continuous column counts and actual value-file bytes. For example:

```text
Reading annotation inputs: baseline files=22, query files=22.
Checking SNP identities and preparing chromosome annotations.
Annotation storage: binary=50 packed-bit columns, continuous=3 float32 columns; value files=18,255,632 bytes.
Annotation preparation complete: chromosomes=22, retained SNPs=1,000,000, elapsed=30.00s.
```

The counts, bytes, and time above are illustrative.
The first line precedes the validation/classification pass; the second precedes global identity validation, retained-row selection, and a second input pass that writes final compact stores directly.
The storage line counts each supplied logical column once; bytes sum the completed NPY files across retained chromosomes, including headers but excluding Parquet metadata, identity/diagnostic scratch, and any subsequently generated queries.
These are file sizes, not peak RSS. Completion appears only after usable annotation shards are ready, with the number of retained chromosomes, logical SNP rows after identity cleanup and any chromosome selection, and elapsed preparation time.
Aligned baseline/query files describe the same logical rows and do not double the SNP count.
For the gene LD-score index builder these messages precede `Starting chromosome N`; preparation remains serial.

These milestones accompany structured phase records from `_runtime.progress.PhaseProgress`: start, complete, failed, and time-throttled progress at bounded chunk/object checkpoints after at least 30 seconds. INFO reports the phase, completed units, known totals, and elapsed time. Individual source/query names remain in structured records and appear in text at DEBUG. Annotation scans share one INFO phase across the source collection; individual scans have DEBUG phases. File totals count completed sources, while row totals in staging phases count scanned/staged rows, not unique genome-wide SNPs. A blocking call cannot emit an intermediate record until control returns.

One collection-level notice explains that supplied annotation CM/MAF fields are ignored because the reference panel is authoritative; DEBUG retains file-specific notices and inferred column mappings. BED source preparation has time-throttled collection progress. Chromosome milestones remain at INFO. Preparation completion means usable annotation shards are ready; computation and publication can still fail later. WARNING/ERROR levels suppress INFO progress through the existing policy.

Intentional gene exclusions use this form, after the audit is saved:

```text
Genes intentionally excluded by region policy: role=focal; sources=5607; input rows=18432. See diagnostics/gene_list_audit.tsv.gz.
```

The counts are illustrative. Rejection summaries point to the gene row audit or resolution summary, and query outcome summaries point to `diagnostics/query_annotation_status.tsv`. Live LD-score catalog defects also have a complete `diagnostics/gene_catalog_issues.tsv` table, using the existing catalog-issue schema shared with standalone `annotate`. This additive diagnostic is owned by the LD-score writer and linked from root metadata when present; scientific table schemas are unchanged.
Mechanisms: [`_annotation/sources.py`](../../src/ldsc/_annotation/sources.py), `prepare_annotation_sources()`; [`_annotation/parsing.py`](../../src/ldsc/_annotation/parsing.py), `normalize_annotation_chunk()`; [`_annotation/query_outcomes.py`](../../src/ldsc/_annotation/query_outcomes.py), `log_gene_list_rejections()`.

## LD-score chromosome diagnostics

Direct query validation logs the validated baseline/reference chromosome sets, effective scope, and ordinary-glob selection caveat. The calculation log explicitly names chromosomes resolved and entering analysis. Indexed runs name their immutable validated scope. `diagnostics/chromosome_scope.json` persists this evidence; its `analysis_chromosomes` is empty when scope validation prevents analysis. The LD-score writer owns this file and `diagnostics/input_issues.tsv` along with existing query/gene diagnostics, so preflight, replacement, and diagnostics-only cleanup share one declaration. Scope is also carried on `LDScoreResult.chromosome_scope` and in root scientific metadata.

## Failed Overwrite Markers

Every public materializing CLI workflow and corresponding high-level Python
workflow wraps its existing action order with the marker helper in
`ldsc._runtime.workflow_logging`. If an attempt authorized with `--overwrite` or
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

The marker helper uses `_runtime.paths.normalize_path_token`, matching output preflight and publication.
CLI marker scope follows the last occurrence of a destination option, as argparse does.
Marker-only directories do not collide with writer artifact families.
Gene LD-score index preflight, direct publication, and interrupted-publication recovery share recognition of owned diagnostics; this excludes symlinks, special files, and unrelated empty subdirectories.
Existing recovery of a validated owned backup preserves the failure marker until the current build succeeds.
See [`test_output_directory_retries.py`](../../tests/test_output_directory_retries.py).

CLI help and `munge-sumstats --infer-only` bypass the marker lifecycle entirely. Even with `--overwrite`, they create no output directories or files and never create, replace, or remove a failure marker. This also applies when inference or argument parsing fails.

## Log Layout

In the main workflow log, `INFO` records are plain message text.
The shared formatter prefixes other severities with `[LEVEL]`, for example `[DEBUG]`, `[WARNING]`, or `[ERROR]`.
Lifecycle audit sections bypass that formatter.
LD-score chromosome logs use timestamps and severity labels for detailed task records, as described under [LD-score chromosome logs](#ld-score-chromosome-logs).

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

Arguments are shell-quoted where necessary, and negative numeric option values remain paired with their option.
Do not abbreviate the recorded invocation when compacting long input or output summaries.

`Inputs:` and `Outputs:` are audit sections, not machine-readable manifests.
`Outputs:` is separated from preceding workflow records by a blank line so summary blocks such as `Metadata:` remain visually distinct.
Numbered LD-score query files are represented by their file count and a pointer to the complete `metadata.json.query_batches` manifest.
Diagnostic paths remain individually listed.
Large `annotate`/`rg` input groups use counts and provenance pointers; the exact `Call:` is preserved.

Liftover drop summaries in ordinary logs are count-only. Row-level drop audit
records belong in the workflow's dropped-SNP sidecar, and example SNPs are
emitted only when the logger is set to `DEBUG`.

The footer records final status and elapsed time with explicit units:

```text
Finished 2026-05-02 22:07:50
Elapsed time: 2.0min:12s
```

Runs shorter than one minute use seconds with two decimal places, for example `Elapsed time: 0.16s`.
Longer runs retain the minutes-and-seconds format shown above.

Counts should remain exact, and measurements should include units and enough precision to interpret the check or outcome.
The shared `format_float()` helper uses four significant digits when called; it is for log text and does not globally reformat every number in a message.
Logging refinements must not change scientific result-table precision or formats.

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
its destination from an existing `h2` result. `plot` does the same for any
supported result suite.

LD-score logs report binary and continuous annotation counts; names remain in `metadata.json.annotation_types` and DEBUG output.
`partitioned-h2` logs repeat an actionable interpretation warning when continuous annotations are present: legacy numerical proportion/enrichment summaries remain visible, but only coefficient-based fields retain their ordinary interpretation for those annotations.
`quantile-h2` logs the selected fitted model, target, inherited common-MAF rule, number of common reference SNPs, missing exclusions, and realized bin counts; per-bin bounds/counts remain in `quantile_h2.tsv` and DEBUG output.
Row-addressable alignment issues are written to `snp_alignment_issues.tsv.gz` rather than expanded into the log.

## Exact gene LD-score index build log

The exact gene LD-score index builder uses the same lifecycle banner, `Call:`, `Inputs:`, `Outputs:`, and `Finished`/`Failed` footer as the other artifact-building commands.
Its stable path is created before chromosome work so it can be monitored live.
Its concise INFO narrative records resolved configuration and input counts, the baseline/PLINK identifier intersection and configured regression SNPs, start and completion for each chromosome, catalog genes after exclusion, retained SNPs and regression SNPs, atom/operator nonzeros, component bytes, and publication state.
The JSON sidecar is the machine-readable summary; the log renders the same per-chromosome and aggregate measurements rather than recomputing them.

For this builder, `Finished chromosome N` is emitted only after the worker has closed every chromosome payload and atomically installed the shard in the hidden run-specific stage. It means the internal shard is durable and its in-memory record can be released; it does not report partial public publication. The public destination changes only after all chromosome and shared metadata are finalized and the complete stage reloads successfully.

For a successful run, the final lines identify `index_id`, staged reload
validation, complete atomic replacement, payload bytes, peak RSS, and the
validated index path. After the `Finished` footer closes the handler, the
successful log moves into the published index's `diagnostics/`. A failure keeps the `Failed` footer and traceback at the hidden live path without creating or modifying an index destination; a failed overwrite leaves the old scientific index loadable. A retained interrupted stage is never a resumable checkpoint and is discarded on retry. Once replacement
and reload validation complete, cleanup failures are warnings that name the
retained transaction directory and do not change the successful exit status.
The lifecycle log is the status authority, so no separate status JSON is
written.

## Batch and collection progress

Both LD-score backends report each batch start and completion.
Direct mode names LD-score calculation batches and logs the chromosome computation inside each batch.
Indexed mode logs chromosome computation for the active query batch, followed by that batch's numbered assembly phase.
The next batch's chromosome work starts only after the current batch has been streamed into private output and its raw fragments removed.
Batch completion describes completed private work, not a published or resumable result.
Each backend reports one whole-run query-status total after finalization.
Sources: `LDScoreCalculator.run()` in [_ldscore/calculator.py](../../src/ldsc/_ldscore/calculator.py) and `indexed_results()` in [_ldscore/gene_index/batches.py](../../src/ldsc/_ldscore/gene_index/batches.py).

Partitioned regression reports actual completed query outcomes, with status counts, through the same parent-owned progress stream for serial and parallel execution. Parallel collection advances when each result arrives, and timed waits can report unchanged completed counts while workers compute. Serial fits report when control returns. Numbered batch boundaries appear even for fast batches; per-query successful phases require DEBUG. For example, omitting elapsed times:

```text
Phase start: computation; query fitting; completed=0/3000
Phase start: batch; query fitting 1/3; completed=0/1000
Phase progress: computation; query fitting; completed=240/3000; success=240
Phase complete: batch; query fitting 1/3; completed=1000/1000
Phase start: batch; query fitting 2/3; completed=0/1000
...
Phase complete: computation; query fitting; completed=3000/3000; success=3000
```

A completed scan can contain query failures; the status counts, final diagnostic summary, and existing strict/continue policy determine publication.
Transport/process failures still terminate the phase with a failed event.
`rg` similarly reports per-trait heritability and pair-fitting totals, retaining pair-specific failure diagnostics.
Quantile and R² query workflows finish with result counts and table pointers.

Munging emits one immediate warning per invalid P, INFO, or FRQ category within a parsing run, then reports accumulated category counts. Counts can overlap and are separate from exclusive QC removal counts. Subsequent API runs receive their own warnings. The INFO accepted-range message describes the existing [0,2] check; no filter cutoff changed. Complete raw QC row details are not all present in the dropped-SNP sidecar, so the warnings and final QC summary retain their own counts and meaning. Source: `parse_dat()` and `_invalid_value_summary()` in [_sumstats/munge_kernel.py](../../src/ldsc/_sumstats/munge_kernel.py).

## LD-score chromosome logs

Direct and indexed writing workflows create `diagnostics/ldscore.log` and `diagnostics/chromosome-logs/chrN.log`, including for one worker.
The main log lists all chromosome log paths and overall progress.
Workers append detailed records and Python exception tracebacks to their own chromosome files across query batches within the run.
Numerical backend verification and independent Arrow CPU/I/O counts are recorded at the task boundary.
A handled chromosome failure retains logs, reports the chromosome and batch in the main log, cancels pending tasks where possible, and removes temporary results after normal worker shutdown.
A broken pool reports affected in-flight tasks and their paths without claiming a confirmed crash source.
Reruns follow `--overwrite`, resetting owned logs rather than appending an unrelated run.
`LDScoreCalculator.run(..., output_config=None)` retains its no-file-output contract.
See [parallel execution](parallel-execution.md) and [`_ldscore.diagnostics`](../../src/ldsc/_ldscore/diagnostics.py).
