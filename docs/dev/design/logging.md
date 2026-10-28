# Logging Design

Last updated on: 2026-10-27

This page holds the rules for writing log, warning and error messages, and the logging API boundary.
The runtime behavior users see (log files, console routing, output preflight and failure markers) is in [logs and outputs](../../reference/logs-and-outputs.md).

## Audience and authority

The writing rules guide new messages and the audit of existing messages; they do not imply that every current message already meets the standard.

The default level, `INFO`, serves users who know what LDSC3 can do but need not understand its modules, architecture, or numerical implementation.
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
For a path-resolution failure, suggest checking the filename and its visibility on the machine running LDSC3; do not suggest column or delimiter changes before file reading has begun.
For an output collision, offer a fresh output directory, or `--overwrite` only when the user intends to replace previous outputs.
Do not recommend dropping inputs, relaxing filters, or enabling query continuation as a generic repair without explaining how that would change the requested analysis.

The following is an illustrative wording example, not a claim that every current error uses these exact lines.
The diagnostic pointer applies only after that table has been written:

```text
munge-sumstats could not find the raw summary-statistics file '/path/raw.tsv.gz'.
No file matched this path; summary-statistic reading has not started.
Check the filename and its visibility on the machine running LDSC3, then set --raw-sumstats-file to an existing readable file.
Input issues and repairs: /path/results/diagnostics/input_issues.tsv.
```

For an unexpected internal failure, describe the interrupted user task and identify any known chromosome, batch, or query.
Do not attribute an unverified internal failure to bad user data.
Give the saved log location when available and ask for the command, package version, and relevant developer evidence if the user needs to report the problem.
If the cause is unknown, say so; do not invent a scientific repair.
An abnormal worker termination may leave no worker traceback, so affected in-flight tasks are evidence of interrupted work rather than proof of the task that caused the termination.

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

Use the existing [troubleshooting reference](../../user-guide/troubleshooting.md) for run-aborting errors with at least three distinct plausible causes.
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

## Validation timing and completeness

Detect run-stopping conditions before substantial computation.
When a check needs only cheap evidence (paths, headers, BIM tables, artifact metadata, or identifier sets), run it for every input before expensive work such as genotype reading, R² computation, LD-score accumulation, or regression fitting.
A check that needs expensive evidence, such as MAF computed from genotypes, fails as soon as that evidence exists rather than after unrelated work continues.

Collect every safely discoverable issue before aborting.
At each validation stage, check all independent inputs (files, chromosomes, queries, and focal/control lists) and raise one error that summarizes every issue found, instead of stopping at the first.
Use `_runtime.preflight.InputGate`, write the complete repair table to `diagnostics/input_issues.tsv`, and show at most 10 representative issues in the primary message.
An unreadable structural contract may stop its own stage immediately; issues in other independent inputs at that stage are still collected.

## API Boundary

Python convenience wrappers that delegate through parsed workflow functions may
create logs because they use the same output-directory contract as the CLI.
The package console handler is installed only by a normal CLI invocation of `run_cli`.
Python APIs retain application logging configuration and propagate exceptions to their caller.
Log creation depends on the API's output contract:

- `AnnotationBuilder.run(..., output_dir=...)` requires an output parent for private preparation but does not install a workflow log. `run_annotate(...)` owns standalone output preflight, writing, and `diagnostics/annotate.log`.
- `LDScoreCalculator.run(..., output_config=...)` creates `diagnostics/ldscore.log` and chromosome logs for writing runs.
  `LDScoreCalculator.run(..., output_config=None)` creates no files.
- `ReferencePanelBuilder.run(...)` does not create a `build-r2-panel` workflow log.
- `RegressionRunner.estimate_*` methods do not create regression logs.

`SumstatsMunger.run(...)` remains a workflow-level API and keeps writing
`diagnostics/sumstats.log`. `MungeRunSummary.output_paths` excludes the log
path but includes data artifacts such as the dropped-SNP audit sidecar.

The May 2026 logging harmonization design and the June console/file routing plan that preceded this contract remain in git history.
The current behavior above supersedes their earlier assumptions.

### Standalone annotation row counts and downstream inputs

After successful `annotate`, the log reports `Baseline SNP rows: input=N, output=M, dropped=D` for annotation-grid input, or `BIM grid SNP rows: ...` for BIM input.
Counts include zero drops and describe logical grid rows, not annotation columns or the regression SNP subset.
`diagnostics/metadata.json` stores the totals and per-chromosome counts in `snp_rows`; the full row audit remains under `diagnostics/dropped_snps/`.

For annotation-grid input, the log lists output paths and states:

```text
For downstream ldscore, use the exported baseline files with --baseline-annot-sources and the exported query files with --query-annot-sources; do not use the original baseline inputs. Exported paths are listed above.
```

Both exports share the cleaned grid, preventing an alignment failure caused by pairing cleaned queries with unmodified input baselines.
BIM mode emits query files only and has no baseline-value output.

## Streaming memory-workflow diagnostics

Annotation, gene-resolution, direct chromosome-drop, and quantile alignment audits are replayed in bounded chunks into their complete output files.
Results keep counts/statuses or persistent paths; they do not retain every chromosome's audit frame.
Standalone annotate publishes query shards incrementally and detaches its returned bundle from construction scratch.
Batch `partitioned-h2` stages complete per-fit details privately and publishes them after final summary sorting, preserving overwrite and failure-marker contracts.
Explicit `--continue-on-query-error` permits publication of successful query fits after recording failures in `diagnostics/query_status.tsv`.
Each query exception is logged with its name, stage, and traceback; singular jackknife errors also log block indices, retained-row intervals, available panel genomic spans, normal-matrix ranks, and before/after-deletion support counts.
A final warning summarizes successful and failed counts.
Strict scans publish diagnostics alone when any query fails.
See [query failure handling](../../reference/results/partitioned-h2.md#query-failures-and-continuation) and [the memory design](memory-and-cpu.md).

### Partitioned query worker diagnostics

The INFO execution summary reports the requested concurrency, detected available CPU count from `joblib.cpu_count()`, task count, batch cap, and actual worker count.
Every worker targets 1 BLAS thread, including serial execution; worker records describe observed backends and verification status.
Regression forwards its per-model log records to the parent workflow log as before.
Existing phase timings, deterministic status order, and final sorting remain unchanged.
Summed worker timings can exceed elapsed wall time.

Ordinary query-model exceptions include the query name and stage. Process death, transport failures, worker initialization failures, interrupts, and output errors remain fatal under both strict and continue policies and can stop the scan before the final status table is written. The coordinator reports affected in-flight query identities and reaps workers before removing scratch. Sources: `RegressionRunner.estimate_partitioned_h2_batch()` in [regression_runner.py](../../../src/ldsc/regression_runner.py), `_WorkerLog` and `_QueryWorkers` in [_partitioned_h2_parallel.py](../../../src/ldsc/_partitioned_h2_parallel.py), and [query-worker policy](../../reference/regression-options.md#43-query-workers-and-memory).

## Maintaining this contract

Use this page for current rules rather than copying dated implementation instructions.
The May 2026 harmonization design introduced the shared workflow context, exact invocation, input/output sections, and separation of logs from scientific outputs.
Its message-only formatter, console-verbosity rules, regression calls without output directories, direct LD-score API log exclusion, and failure footer without exception details are superseded.
The current command is `build-r2-panel`, and short-run elapsed times use seconds as specified above.

The June 2026 refinement work contributed self-contained errors, command-oriented troubleshooting, original exception context, and synchronized reference links.
Its source counts, line numbers, command map, companion-document paths, and batch/commit instructions describe that historical implementation effort.
They are not current message-writing requirements.
Mandatory conversion of every logging call to f-strings, a fixed two-sentence cap, and logging every otherwise silent public I/O function are not standing rules.
Use the existing interpolation convention and log meaningful operational events.
Expected input/configuration failures should use the appropriate existing [user-error subclass](../../../src/ldsc/errors.py); do not classify an unexpected internal exception as a user mistake solely because its type is `ValueError`.

An audit must check each message at its actual command/API boundary, including the default-level output, available diagnostics, and failure stage.
Separate confirmed failure evidence from suggested causes and wording.
Record any mismatch between the writing standard and current messages as work to be implemented, rather than claiming documentation changes repaired runtime output.
Message edits must preserve scientific meaning, counts, result formats, and established continuation policies.
Changing exception classes or chaining, exit status, warning visibility, handler routing, or control flow requires explicit behavioral scope and appropriate validation.
When runtime changes are authorized, test the meaningful contract: exit status, error context, preserved causes, diagnostic availability, and bounded output.
Avoid tests of logging-library internals or exact incidental prose unless the rendered format itself is a contract.
