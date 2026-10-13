# Workflow Logging

Last updated on: 2026-10-13

## Audience and authority

This is the authoritative, live contract for LDSC logging formats, message-writing rules, and diagnostic routing.
The writing rules guide new messages and the audit of existing messages; they do not imply that every current message already meets the standard.
The command-specific sections below describe current runtime behavior.

The default level, `INFO`, serves users who know what LDSC can do but need not understand its modules, architecture, or numerical implementation.
Messages should explain the task being performed, its outcome, and any action the user needs to take.
Familiar method terms such as SNP, LD score, genome build, and heritability are appropriate; explain specialized terms when their meaning is needed to understand a failure or remedy.

`DEBUG` adds evidence for advanced users and developers: source mappings, per-query timings, implementation stages, numerical checks, and worker/backend details.
Technical detail supplements the ordinary explanation rather than replacing it.
Errors and warnings must remain understandable at the default level.
Failure tracebacks belong in the saved logs even when `DEBUG` was not selected.

## Message-writing rules

- Describe user-visible operations and inputs.
  Prefer “checking the summary-statistics file” to an unexplained gate name or private function name.
  If an implementation term is necessary to interpret a diagnostic, explain its meaning first.
- Use concrete verbs and consistent terms.
  Distinguish a missing file, an unreadable file, invalid contents, incompatible inputs, and a failed calculation instead of calling all of them “invalid input.”
- Identify the relevant command, input role, filename, option, query, trait pair, or chromosome when it helps the user find the problem.
  Name options and diagnostic fields exactly so the user can locate them in help or the saved table.
  Avoid repeating the same pathname in the primary explanation.
- Give counts with their unit and scope.
  State whether a count refers to input rows, distinct genes, SNPs, queries, chromosomes, or failed inputs.
  Distinguish the reference SNPs contributing to LD scores from the SNPs used for regression when that distinction explains the result.
  Include expected and observed values, accepted formats, or units when they clarify a failed check.
- Keep routine progress concise and aggregate repeated events.
  Report meaningful phase, batch, and chromosome boundaries; place successful per-object details at `DEBUG`.
  Use counts and bounded examples for large collections only when complete detail is retained in a diagnostic or metadata file.
  Preserve the exact `Call:` block.
- Report the actual outcome.
  A completed validation, batch, or chromosome does not mean the full run succeeded or final results were written.
  Distinguish requested exclusions, skipped queries, failed queries, and a failed run.
  Do not imply that an unknown cause was established or that a partial result is complete.
- Let the shared formatter supply severity prefixes.
  Avoid embedding `WARNING:` or `[ERROR]` in message bodies that already use those logging levels.
  Keep the primary explanation short enough to scan, but do not sacrifice the failure, cause, or remedy to a fixed sentence limit.

### Levels and audiences

Keep levels in the order `DEBUG`, `INFO`, `WARNING`, `ERROR` in help, docstrings, and documentation.
Each level includes more severe records.

| Level | Intended content |
| --- | --- |
| `DEBUG` | Detailed troubleshooting evidence, including per-source/query records, module or function context, numerical diagnostics, and backend observations. |
| `INFO` | A clear run narrative: important settings, major phases, batch/chromosome progress, counts, and aggregate outcomes. |
| `WARNING` | A condition allowed by the workflow's existing policy that needs attention; explain its consequence and any required action. |
| `ERROR` | An operation that failed; identify the affected task/input, the observed problem, and the next useful action. |

Changing the log level changes detail, not scientific validation or the workflow's failure policy.
Run headers, final status, elapsed time, and captured failure tracebacks remain in an opened workflow log at every level.
The explicit CLI notices and munging summary described below also retain their current visibility.

## User-facing errors and warnings

A user-facing error must stand on its own in the terminal message or Python exception.
The reader should not need to open source code, enable `DEBUG`, or read a traceback to understand a known input or configuration failure.
Use the following content in a natural order; this is a writing standard, not a mandatory rigid template:

1. State what failed and identify the affected input, setting, or calculation.
2. State the observed problem or failed requirement, with useful expected/observed values.
3. Give the established cause when the check proves it.
   Otherwise label the best-supported explanation as a likely or possible cause and give a concrete way to check it.
4. Give the most useful actionable repair, and explain any consequence for the run or results that affects the user's next decision.
   Add a pointer to written diagnostics when they provide additional evidence or repairs.

A remedy should name the file, option, schema requirement, or compatibility check the user can act on.
“Check your input,” “invalid configuration,” and “rerun with DEBUG” alone are insufficient for a known user error.
For a path-resolution failure, suggest checking the filename and its visibility on the machine running LDSC; do not suggest column or delimiter changes before file reading has begun.
For an output collision, offer a fresh output directory, or `--overwrite` only when the user intends to replace previous outputs.
Do not recommend dropping inputs, relaxing filters, or enabling query continuation as a generic repair without explaining how that would change the requested analysis.

The following is an illustrative wording example, not a claim that every current error uses these exact lines.
The diagnostic pointer applies only after that table has been written:

```text
munge-sumstats could not find the raw summary-statistics file '/path/raw.tsv.gz'.
No file matched this path; summary-statistic reading has not started.
Check the filename and its visibility on the machine running LDSC, then set --raw-sumstats-file to an existing readable file.
Input issues and repairs: /path/results/diagnostics/input_issues.tsv.
```

For an unexpected internal failure, describe the interrupted user task and identify any known chromosome, batch, or query.
Do not attribute an unverified internal failure to bad user data.
Give the saved log location when available and ask for the command, package version, and relevant developer evidence if the user needs to report the problem.
If the cause is unknown, say so; do not invent a scientific repair.
An abnormal worker termination may leave no child traceback, so affected in-flight tasks are evidence of interrupted work rather than proof of the task that caused the termination.

Warnings explain what happened, whether processing continues, and how the condition affects the requested analysis or outputs.
An explicitly requested exclusion can be summarized with its reason, count, and saved audit.
A skip or omission must be explained in words; a status code such as `zero_variance_ld_scores` is not a sufficient explanation by itself.
The existing scientific-error and query-continuation policies decide whether processing may continue.
Unrequested input-resolution, filtering, ambiguity, or compatibility conditions that can silently change the analysis must fail before substantial computation; message cleanup must not turn them into successful-run warnings.

### Diagnostics and troubleshooting

Keep the common explanation and top remedy in the primary message.
A diagnostic pointer supplements that explanation and may replace exhaustive row or object lists only when the corresponding file was written and contains the relevant detail.
Name the actual path and tell the user what evidence or repairs to inspect there.
If failure occurred before a log or table was created, do not imply that it exists.
Warnings without complete row-level diagnostics, including raw summary-statistic QC warnings, remain self-contained.

For multiple input failures, summarize issue counts and affected input roles, retain enough representative context to identify the failed check, and point to the complete repair table.
Distinguish issue counts from counts of distinct failed inputs.
Bound the primary list when a complete table is available; do not replace specific failures with an unexplained “validation failed.”
The shared staged input gate shows at most 10 representative issues after writing the complete repair table.
When no complete table was written, it retains every issue in the primary error.

Use the existing [troubleshooting reference](../troubleshooting.md) for run-aborting errors with at least three distinct plausible causes.
Organize entries by the command the user runs, using `Common` for failures shared across commands and plain-language headings.
Rank causes by available evidence and pair each with a concrete check and remedy.
Simple one- or two-cause errors should be explained in the message without creating a separate reference entry.
Keep every in-code heading link synchronized with its target.
Human-readable messages do not require a new error-ID scheme; existing machine-readable reason codes remain in their diagnostics.

### Developer evidence and error boundaries

Preserve the original exception chain when a translation adds useful user context.
Keep exception types, operating-system error numbers, module/function details, numerical evidence, and tracebacks in developer logs or structured diagnostics where they help establish the cause.
The ordinary error should explain the failure in user terms without requiring the reader to interpret those details.
Full tracebacks for captured workflow or task failures do not depend on selecting `DEBUG`.

Reusable APIs raise to their caller.
The CLI boundary owns the final terminal error report; avoid logging the same error at every layer and then reporting it again.
Task-specific tracebacks and query-failure diagnostics may still be recorded at their existing boundaries.
A single concise terminal report and richer file evidence serve different purposes.
Cleanup, wrapping, or continuation must not hide the original failure or misreport success.

## Workflow log lifecycle

Public workflow entry points share one logging policy:

- Per-run file handlers attach to the `LDSC` logger, so workflow and kernel
  child records are captured together.
- Known scientific outputs and the log path are preflighted before the log file
  is opened. For workflows with coherent output families, this preflight covers
  owned siblings that are not produced by the current run. If preflight fails,
  no new log file is created.
- If an exception escapes an opened workflow logging context, the log is kept with a `Failed` footer, the full Python traceback, and elapsed time.
  Abrupt process termination can interrupt logging before a footer or traceback is written.
- `--log-level` controls module-record verbosity. With an output directory, those
  records go to the workflow log file (see routing below); lifecycle audit lines
  always appear in the file. Supported levels are `DEBUG`, `INFO`, `WARNING`, and
  `ERROR`.
- Workflow result objects and `output_paths` mappings do not include log files.
- INFO describes major phases, numbered batch start/completion, chromosome work, and aggregate outcomes. Successful query timings, source mappings, and full annotation lists belong at DEBUG.
- Gene exclusions, rejected rows, unsupported genes, and non-OK query outcomes use counts by role/reason and pointers to complete diagnostics. Exclusion counts describe submitted input rows, including repeats across lists; they are not counts of distinct genes. Failure messages retain the affected query/pair, stage, and traceback.
- A diagnostic pointer replaces detail only when the corresponding file is written. Warnings without an equivalent complete diagnostic, including raw summary-statistic QC warnings, remain self-contained.

Nested workflow logging contexts retain the active outer log rather than opening a second file or replacing its settings.
The owning context restores the previous logger level and removes its own handler on exit.

## Console vs File Routing

All command help pages reuse `LOG_LEVEL_HELP` from [`ldsc._logging`](../../src/ldsc/_logging.py).
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

The package does not install a root logger handler, and LDSC records may propagate to handlers installed by an application.
Claims of “file only” above describe the package's own routing, not suppression of application logging.
Calling `run_cli(argv=...)` programmatically does not install the CLI console handler.

Gene-list `ldscore` adds two deliberately bounded, direct stderr notices after a
successful CLI run: one when `--allow-unresolved-genes` selected `resolved-only` and omitted rows, and one when Gate B
found zero-SNP support or skipped/warning query outcomes. These notices bypass
the ordinary logger threshold so a successful SLURM job cannot hide a
science-relevant subset/skip in a log that users may never open. They name the
diagnostic paths and cap affected outcomes at 10. Python entry points do not
emit these console notices; callers inspect the returned statuses and paths.

The successful gene-list notice distinguishes unsupported gene-list rows from queries skipped after LD-score computation.
It reports an unsupported-row count only when that count is positive, explains each SNP-support or zero-variance outcome in words, and states that skipped queries are omitted from LD-score outputs.
For example, a notice can report three queries skipped because their LD scores are constant on regression SNP rows (`zero_variance_ld_scores`), even when every evaluated gene has retained reference-SNP support.
SNP support is measured on the retained reference-SNP universe; LD-score variation is evaluated on regression SNP rows.
The reason codes, diagnostic schemas, query-selection policy, and console/file routing remain unchanged.

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

## Staged input validation

Output ownership and input validity are separate gates. Existing artifact-family checks authorize logs/diagnostics and prevent replacement without overwrite. Input gates do not authorize scientific output replacement, remove stale science files, or imply that publication succeeded. R2 builds derive their chromosome-specific output family from validated PLINK contents; their early discovery repair table can therefore precede the full output-family check, as the previous PLINK repair table did.

[`_input_preflight.InputGate`](../../src/ldsc/_input_preflight.py) collects independent input errors across objects at a stage and retains `input_issues` on the exception. Writing workflows persist the existing six-column `input_issues.tsv` repair table under diagnostics; gene-index builds use the sibling `.<index>.build-state/` directory. Existing PLINK and gene-specific audits remain available. The table is a diagnostic, not a new scientific artifact or resumable checkpoint.

Input-preflight errors name the failed input, report the issue count and an actionable repair, and point to the written `diagnostics/input_issues.tsv`.
A missing raw-summary-statistic path identifies `--raw-sumstats-file`, reports that path resolution found no matching file, and states that summary-statistic reading has not started.
The primary message shows that pathname once; the chained developer traceback retains the original exception and its full context.
The diagnostic table retains its six-column schema and preserves specialized nested diagnostics.
Other captured operation errors distinguish missing and unreadable inputs and record exception types and operating-system error numbers when available.
One captured failure becomes the error's direct cause; multiple captured failures remain available through an `ExceptionGroup` cause.
The declaration-validation phase reports `failed_inputs` as the number of distinct role/source/chromosome combinations with issues.
When a phase total is unknown, progress omits the denominator; any reported completed count uses `completed=N`.

Direct LD scores check declarations and annotation headers, resolve gene-list controls, and validate reference content before annotation staging. PLINK companion mappings and reference scope evidence are reused; annotation headers supply the source scan widths. Gene-index and R2 construction check baseline/auxiliary declarations and PLINK companions together. Index consumers check every component path, small chromosome identity metadata, and Parquet header before opening the first numerical operator. See the [verified call-order audit](../audits/2026-09-14-staged-input-preflight.md).

Validation remains staged: BIM/FAM rows and BED dimensions require reference inspection; annotation row alignment, duplicate identities, chromosome contents, and gene support require their content scans. R2 query-scope validation streams pair integrity before preparation and reuses that scope result. Baseline-only R2 runs check sidecar/reader compatibility early and retain pair-content checks in the consuming reader, avoiding an additional full pair scan. Ordinary globs cannot reveal a chromosome that is absent consistently from all matches. Scientific omission, ignored-option, overwrite, and failure-marker policies remain unchanged.

## Annotation preparation

Shared source preparation emits INFO milestones in the workflow log for reading annotation inputs, checking SNP identities and preparing chromosome annotations, and successful completion. A storage summary immediately before completion records binary/continuous column counts and actual value-file bytes. For example:

```text
Reading annotation inputs: baseline files=22, query files=22.
Checking SNP identities and preparing chromosome annotations.
Annotation storage: binary=50 packed-bit columns, continuous=3 float32 columns; value files=18,255,632 bytes.
Annotation preparation complete: chromosomes=22, retained SNPs=1,000,000, elapsed=30.00s.
```

The counts, bytes, and time above are illustrative. The first line precedes the validation/classification pass; the second precedes global identity validation, retained-row selection, and a second input pass that writes final compact stores directly. The storage line counts each supplied logical column once; bytes sum the completed NPY files across retained chromosomes, including headers but excluding Parquet metadata, identity/diagnostic scratch, and any subsequently generated queries. These are file sizes, not peak RSS. Completion appears only after usable annotation shards are ready, with the number of retained chromosomes, logical SNP rows after identity cleanup and any chromosome selection, and elapsed preparation time. Aligned baseline/query files describe the same logical rows and do not double the SNP count. For the gene-index builder these messages precede `Starting chromosome N`; preparation remains serial.

These milestones accompany structured phase records from `_progress.PhaseProgress`: start, complete, failed, and time-throttled progress at bounded chunk/object checkpoints after at least 30 seconds. INFO reports the phase, completed units, known totals, and elapsed time. Individual source/query names remain in structured records and appear in text at DEBUG. Annotation scans share one INFO phase across the source collection; individual scans have DEBUG phases. File totals count completed sources, while row totals in staging phases count scanned/staged rows, not unique genome-wide SNPs. A blocking call cannot emit an intermediate record until control returns.

One collection-level notice explains that supplied annotation CM/MAF fields are ignored because the reference panel is authoritative; DEBUG retains file-specific notices and inferred column mappings. BED source preparation has time-throttled collection progress. Chromosome milestones remain at INFO. Preparation completion means usable annotation shards are ready; computation and publication can still fail later. WARNING/ERROR levels suppress INFO progress through the existing policy.

Intentional gene exclusions use this form, after the audit is saved:

```text
Genes intentionally excluded by region policy: role=focal; sources=5607; input rows=18432. See diagnostics/gene_list_audit.tsv.gz.
```

The counts are illustrative. Rejection summaries point to the gene row audit or resolution summary, and query outcome summaries point to `diagnostics/query_annotation_status.tsv`. Live LD-score catalog defects also have a complete `diagnostics/gene_catalog_issues.tsv` table, using the existing catalog-issue schema shared with standalone `annotate`. This additive diagnostic is owned by the LD-score writer and linked from root metadata when present; scientific table schemas are unchanged.
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
`Outputs:` is separated from preceding workflow records by a blank line so summary blocks such as `Metadata:` remain visually distinct. Numbered LD-score query files are represented by their file count and a pointer to the complete `metadata.json.query_batches` manifest. Diagnostic paths remain individually listed. Large annotate/rg input groups use counts and provenance pointers; the exact `Call:` is preserved.

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
its destination from an existing h2 result. `plot` does the same for any
supported result suite.

LD-score logs report binary and quantitative annotation counts; names remain in `metadata.json.annotation_types` and DEBUG output. Partitioned-h2 logs repeat an actionable interpretation warning when quantitative annotations are present: legacy numerical proportion/enrichment summaries remain visible, but only coefficient-based fields retain their ordinary interpretation for those annotations. Quantile-h2 logs the selected fitted model, target, inherited common-MAF rule, common reference-SNP universe size, missing exclusions, and realized bin counts; per-bin bounds/counts remain in `quantile_h2.tsv` and DEBUG output. Row-addressable alignment issues are written to `snp_alignment_issues.tsv.gz` rather than expanded into the log.

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
create logs because they use the same output-directory contract as the CLI.
The package console handler is installed only by a normal CLI invocation of `run_cli`.
Python APIs retain application logging configuration and propagate exceptions to their caller.
Log creation depends on the API's output contract:

- `AnnotationBuilder.run(..., output_dir=...)` requires an output parent for private preparation but does not install a workflow log. `run_annotate(...)` owns standalone output preflight, writing, and `diagnostics/annotate.log`.
- `LDScoreCalculator.run(..., output_config=...)` creates `diagnostics/ldscore.log` and chromosome logs for writing runs.
  `LDScoreCalculator.run(..., output_config=None)` creates no files.
- `ReferencePanelBuilder.run(...)` does not create a build-r2-panel workflow log.
- `RegressionRunner.estimate_*` methods do not create regression logs.

`SumstatsMunger.run(...)` remains a workflow-level API and keeps writing
`diagnostics/sumstats.log`. `MungeRunSummary.output_paths` excludes the log
path but includes data artifacts such as the dropped-SNP audit sidecar.

For historical implementation rationale, see the [logging harmonization design](../specs/2026-05-02-logging-harmonization-design.md) and the [console/file routing plan](../plans/2026-06-04-logging-console-file-routing-plan.md).
The current behavior above supersedes their earlier assumptions.

## Streaming memory-workflow diagnostics

Annotation, gene-resolution, direct chromosome-drop, and quantile alignment audits are replayed in bounded chunks into their complete output files. Results keep counts/statuses or persistent paths; they do not retain every chromosome's audit frame. Standalone annotate publishes query shards incrementally and detaches its returned bundle from construction scratch. Batch partitioned-h2 stages complete per-fit details privately and publishes them after final summary sorting, preserving overwrite and failure-marker contracts. Explicit `--continue-on-query-error` permits publication of successful query fits after recording failures in `diagnostics/query_status.tsv`. Each query exception is logged with its name, stage, and traceback; singular jackknife errors also log block indices, retained-row intervals, available panel genomic spans, normal-matrix ranks, and before/after-deletion support counts. A final warning summarizes successful and failed counts. Strict scans publish diagnostics alone when any query fails. See [query failure handling](partitioned-h2-results.md#query-failures-and-continuation) and [the memory design](annotation-memory-design.md).

### Partitioned query worker diagnostics

The INFO execution summary reports the requested concurrency, detected available CPU count from `joblib.cpu_count()`, task count, batch cap, and actual worker count.
Every worker targets 1 BLAS thread, including serial execution; worker records describe observed backends and verification status.
Regression forwards its per-model log records to the parent workflow log as before.
Existing phase timings, deterministic status order, and final sorting remain unchanged.
Summed worker timings can exceed elapsed wall time.

Ordinary query-model exceptions include the query name and stage. Process death, transport failures, worker initialization failures, interrupts, and output errors remain fatal under both strict and continue policies and can stop the scan before the final status table is written. The coordinator reports affected in-flight query identities and reaps workers before removing scratch. Sources: `RegressionRunner.estimate_partitioned_h2_batch()` in [regression_runner.py](../../src/ldsc/regression_runner.py), `_WorkerLog` and `_QueryWorkers` in [_partitioned_h2_parallel.py](../../src/ldsc/_partitioned_h2_parallel.py), and [query-worker policy](regression-configuration.md#43-query-workers-and-memory).


## Batch and collection progress

Both LD-score backends report each batch start and completion. Direct mode names LD-score calculation batches and logs the chromosome computation inside each batch. Indexed mode logs chromosome computation first, then numbers genome-wide assembly batches; one chromosome computes all of its query batches before assembly. Batch completion describes completed private work, not a published or resumable result. Each backend reports one whole-run query-status total after finalization. Sources: `LDScoreCalculator.run()` in [ldscore_calculator.py](../../src/ldsc/ldscore_calculator.py) and `indexed_results()` in [_indexed_ldscore_batches.py](../../src/ldsc/_indexed_ldscore_batches.py).

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

A completed scan can contain query failures; the status counts, final diagnostic summary, and existing strict/continue policy determine publication. Transport/process failures still terminate the phase with a failed event. `rg` similarly reports per-trait h2 and pair-fitting totals, retaining pair-specific failure diagnostics. Quantile and R² query workflows finish with result counts and table pointers.

Munging emits one immediate warning per invalid P, INFO, or FRQ category within a parsing run, then reports accumulated category counts. Counts can overlap and are separate from exclusive QC removal counts. Subsequent API runs receive their own warnings. The INFO accepted-range message describes the existing [0,2] check; no filter cutoff changed. Complete raw QC row details are not all present in the dropped-SNP sidecar, so the warnings and final QC summary retain their own counts and meaning. Source: `parse_dat()` and `_invalid_value_summary()` in [_kernel/sumstats_munger.py](../../src/ldsc/_kernel/sumstats_munger.py).

## LD-score chromosome logs

Direct and indexed writing workflows create `diagnostics/ldscore.log` and `diagnostics/chromosome-logs/chrN.log`, including for one worker.
The main log lists all chromosome log paths and overall progress.
Children append detailed records and Python exception tracebacks to their own chromosome files across query batches within the run.
Numerical backend verification and independent Arrow CPU/I/O counts are recorded at the task boundary.
A handled chromosome failure retains logs, reports the chromosome and batch in the main log, cancels pending tasks where possible, and removes temporary results after normal worker shutdown.
A broken pool reports affected in-flight tasks and their paths without claiming a confirmed crash source.
Reruns follow `--overwrite`, resetting owned logs rather than appending an unrelated run.
`LDScoreCalculator.run(..., output_config=None)` retains its no-file-output contract.
See [parallel execution](parallel-execution.md) and [`_ldscore_diagnostics`](../../src/ldsc/_ldscore_diagnostics.py).

## Maintaining this contract

Use this page for current rules rather than copying dated implementation instructions.
The [May harmonization design](../specs/2026-05-02-logging-harmonization-design.md) remains useful for the shared workflow context, exact invocation, input/output sections, and separation of logs from scientific outputs.
Its message-only formatter, console-verbosity rules, regression calls without output directories, direct LD-score API log exclusion, and failure footer without exception details are superseded.
The current command is `build-r2-panel`, and short-run elapsed times use seconds as specified above.

The [June refinement plan](../plans/2026-06-04-logging-refinement-plan.md) contributes self-contained errors, command-oriented troubleshooting, original exception context, and synchronized reference links.
Its source counts, line numbers, command map, companion-document paths, and batch/commit instructions describe that historical implementation effort.
They are not current message-writing requirements.
Mandatory conversion of every logging call to f-strings, a fixed two-sentence cap, and logging every otherwise silent public I/O function are not standing rules.
Use the existing interpolation convention and log meaningful operational events.
Expected input/configuration failures should use the appropriate existing [user-error subclass](../../src/ldsc/errors.py); do not classify an unexpected internal exception as a user mistake solely because its type is `ValueError`.

An audit must check each message at its actual command/API boundary, including the default-level output, available diagnostics, and failure stage.
Separate confirmed failure evidence from suggested causes and wording.
Record any mismatch between the writing standard and current messages as work to be implemented, rather than claiming documentation changes repaired runtime output.
Message edits must preserve scientific meaning, counts, result formats, and established continuation policies.
Changing exception classes or chaining, exit status, warning visibility, handler routing, or control flow requires explicit behavioral scope and appropriate validation.
When runtime changes are authorized, test the meaningful contract: exit status, error context, preserved causes, diagnostic availability, and bounded output.
Avoid tests of logging-library internals or exact incidental prose unless the rendered format itself is a contract.
