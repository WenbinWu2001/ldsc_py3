# Logging audit for large annotation batches

Last updated on: 2026-09-15

## Executive summary

Keep default `INFO` logs focused on the run: input counts, scientific configuration, major stages, periodic completed-work counts, final outcomes, and output locations. Move successful per-query timings and exhaustive name lists to `DEBUG`; keep complete gene/query details in the existing diagnostic tables. Removing query-name prefixes alone would leave the underlying repetition intact.

The largest observed contributors are different across workflows:

- `partitioned-h2`: four start/completion pairs for each successful query, forwarded from workers into the main log.
- Direct `ldscore`: thousands of intentional-exclusion summaries, enormous annotation-classification lines, and repeated configuration messages for successive batches.
- Indexed `ldscore`: the same exclusion summaries, followed by batch-local status summaries that never present a clearly labeled whole-run total.

The audit below records the original evidence and recommendations. The user subsequently approved implementation across all affected commands, including chromosome progress and explicit batch boundaries; see the implementation status below. The local repository is `ldsc_py3_restructured` on branch `restructure`. Existing changes to `src/ldsc/_logging.py` and `tests/test_logging_refactor.py` concern negative command-line argument formatting and were preserved.

## Implementation status

Implemented the approved policy across direct/indexed LD scores, partitioned-h2, rg, standalone annotate, shared gene-index annotation preparation, munging, quantile-h2, and query-r2. The other inspected commands had bounded output and did not need artificial edits. Successful per-query phases and source mappings are DEBUG; INFO contains major stages, numbered batch start/completion, chromosome work, time-throttled completed counts, and final outcomes. Regression parent progress advances on actual worker returns, preserving returned order and strict/continue decisions. Timed waits can report unchanged completed counts; no separate active-worker metric is added.

Exclusion/rejection/support summaries use counts and existing diagnostic pointers. Large annotation/trait/output lists use counts and metadata references, preserving the exact invocation. Quantitative-annotation interpretation warnings retain their scientific guidance. Munging warnings remain self-contained because the dropped-SNP audit does not cover every raw QC row: one warning per invalid category is followed by cumulative counts, with no filter change. INFO range wording now matches the existing [0,2] predicate.

One additional audit gap was confirmed during implementation: live LD-score catalog defects existed only as repeated log rows. The LD-score writer now owns and writes `diagnostics/gene_catalog_issues.tsv` when defects are present, including diagnostic-only gene gates, and links it from root metadata. It reuses the standalone annotate catalog-issue schema. This is an additive diagnostic; scientific schemas and numerical behavior are preserved.

Behavior tests exercise real direct/indexed LD scores, inline/spawned regression, many annotation/BED sources, rg totals, quantile results, repeated munging runs at different chunk sizes, complete catalog diagnostics, and compact output listings. Final verification results are recorded at the end of this document. The examples and proposed details later in this audit describe the original recommendation; the current [logging contract](../current/workflow-logging.md) describes implemented messages.

## Evidence and coverage

### Longleaf snapshots

Only the three authorized logs were read. Paths were resolved and checked against the permitted roots and the `data/` exclusion before reading. No remote package code, input data, scheduler jobs, or results were changed. Paths mentioned inside log messages were not followed.

All log paths below are relative to:

```text
/users/w/e/wenbinwu/Sullivan/scz2026_analysis/results/h2-enrichment-all-gene-sets/ldsc-seg_on_server/new/
```

| Log | Snapshot in UTC | Observed volume | Dominant messages |
| --- | --- | --- | --- |
| `ldscores/diagnostics/ldscore.log` | 2026-09-15 14:55:30 | 5,766 lines; 1,978,341 bytes; no completion footer | 5,607 intentional-exclusion lines; 14 annotation-classification lines containing 735,470 characters; the longest classification line has 68,090 characters |
| `ldscores-indexed-mode/diagnostics/ldscore.log` | 2026-09-15 14:55:30 | 5,706 lines; 1,238,520 bytes; completed footer | 5,607 intentional-exclusion lines; 21 non-OK query notices; 19 batch-local status summaries |
| `regression-indexed-mode/diagnostics/partitioned-h2.log` | 2026-09-15 14:53:35 | 3,707 lines; 373,599 bytes; no completion footer | 3,629 phase-start/timing lines, or 97.9% of all lines; 3,271 lines carry a query prefix |

The exclusion messages alone account for 98.3% of lines in the completed indexed LD-score log. Summing its 19 status summaries gives 18,809 submitted query outcomes: 18,788 `ok`, zero `warning`, and 21 `skipped`. These are different quantities from the 18,788 queries entering regression.

The unfinished regression log grew during inspection. Counts above refer to the stated snapshot, not the eventual run total. The remote LD-score logs also lack some progress messages present in the current local source. Findings distinguish observed output from patterns established by inspecting the local code; an exact remote code revision was not verified. The pasted excerpt independently demonstrates both query-prefixed and generic repeated estimator timings.

### Local coverage

Inspected all 13 CLI workflow surfaces and their relevant shared logging, progress, input preparation, worker coordination, and diagnostic-writing paths. Searched package logging calls, direct printing, and Python warnings; traced the high-volume call sites through their loops. Inspected the current logging contract, existing troubleshooting sections, and relevant message/diagnostic tests.

## Findings

### F1 — High: per-query regression timers overwhelm batch progress

**Evidence.** `_fit_partitioned_query()` surrounds model preparation, estimation, summarization, and staging with `_log_phase_timing()`, which emits both boundaries at `INFO`. That is eight routine records per successful query: 150,304 records if all 18,788 requested models succeed, before batch and other messages. `_WorkerLog` buffers each model's messages; `_QueryWorkers.map()` replays them with a query prefix when that model returns. The messages arrive in bursts and do not report overall completion. See [regression_runner.py](../../src/ldsc/regression_runner.py), `_log_phase_timing()`, `_fit_partitioned_query()`, and `estimate_partitioned_h2_batch()`; [_partitioned_h2_parallel.py](../../src/ldsc/_partitioned_h2_parallel.py), `_WorkerLog`, `_query_worker()`, and `_QueryWorkers.map()`.

**Impact.** Users cannot quickly distinguish preparation, fitting, and publication, or tell how much of the scan is complete. Query names and millisecond staging timings consume nearly the entire log.

**Suggested direction.** Retain shared preparation, overall query fitting, and final writing at `INFO`. Put individual successful model phases at `DEBUG`. Emit aggregate progress from the parent when outcomes arrive, with attempted/total, successful, unestimable, failed, active workers, and elapsed wall time. Inline execution should use the same counters. Update the parallel counter inside `map()` as each outcome arrives; updating only after it returns would leave a whole loading batch without progress. Keep query identity, stage, exception, and jackknife evidence for failures.

**Behavioral risk.** Changing timing levels is a logging-contract change. Progress coordination must preserve task order in returned results, deterministic publication, strict/continue error policy, worker cleanup, and numerical work. Do not classify every worker `INFO` record as disposable: some records describe scientific filtering. **Confidence: High. Audience: users and operators.**

### F2 — High: gene-list auditing duplicates thousands of detail records into the log

**Evidence.** `_log_gene_list_rejections()` prints one INFO line per source/role with every intentionally excluded line/gene pair. It also prints every rejected row at WARNING. `_log_gene_list_snp_support()` prints every unsupported row, and `_log_query_annotation_statuses()` prints every non-OK query. The supplied logs each contain 5,607 intentional-exclusion records. See [query_annotations.py](../../src/ldsc/query_annotations.py), those four helpers. The existing per-source format is explicitly documented in [Workflow Logging, Annotation preparation](../current/workflow-logging.md#annotation-preparation) and asserted by [test_gene_resolution_storage.py](../../tests/test_gene_resolution_storage.py), `test_exclusion_log_summarizes_each_source_across_audit_chunks()`.

**Impact.** A rule intended to be concise per source becomes the dominant output across thousands of sources. Row-level audit replay also does work purely to build log messages.

**Suggested direction.** Report intentional exclusions once per stage and role: policy, affected source count, and excluded input-row count. Label these as input rows, since the same gene can occur in multiple lists. For rejected, unsupported, and skipped outcomes, keep visible summaries by reason with totals and a diagnostic pointer; include at most a small fixed number of examples, such as 10 across a stage. Surface the first scientific issue promptly and final totals when the gate or stage finishes. Preserve every source, gene, physical line, reason, and query status in `gene_list_audit.tsv.gz`, `gene_list_resolution_summary.tsv`, and `query_annotation_status.tsv`; optional detailed logging belongs at `DEBUG`.

Bound the source breakdown in fatal gate messages as well: `_gene_list_gate_a_message()` currently caps some row examples but can still concatenate every source, and `_all_query_annotations_skipped_message()` lists every query. A concise exception should retain counts, representative identities, remediation, and the complete diagnostic location.

**Behavioral risk.** This deliberately revises the current requirement to print every rejected/unsupported/non-OK item. It must preserve exhaustive validation and audit contents, warning severity, scientific omission notices, and abort/continue decisions. Consolidate diagnostics only after confirming that the referenced audit has been written; do not replace an actionable issue with an unavailable file pointer. **Confidence: High. Audience: users.**

### F3 — High: LD-score batch internals present themselves as repeated whole runs

**Evidence.** `LDScoreCalculator.run()` invokes `_run_batch()` repeatedly; `_run_batch()` repeats the global configuration, reference CM/MAF source, start/completion messages, and complete binary/quantitative annotation lists. The direct snapshot has 15 configuration/start blocks and 14 completed classification blocks. The indexed generator invokes `_log_query_annotation_statuses()` separately for each output batch. See [ldscore_calculator.py](../../src/ldsc/ldscore_calculator.py), `LDScoreCalculator.run()` and `_run_batch()`; [_indexed_ldscore_batches.py](../../src/ldsc/_indexed_ldscore_batches.py), `indexed_results()`.

**Impact.** Repeated setup messages obscure the actual batch ordinal and total. An unlabeled `ok=1000` is neither a global total nor a useful indication of remaining work. Large annotation lists can overwhelm a log even when line count is modest.

**Suggested direction.** Report resolved configuration and reference provenance once per workflow. Log annotation counts by role and type, keeping full names in existing metadata and DEBUG output. Identify batch number, total batches, query-column count, and completed chromosome work in a bounded aggregate progress stream. Emit one final global query-status summary after all batches are finalized; do not repeatedly count shared baseline columns or early skipped inputs.

Direct mode computes chromosome work within each query batch. Indexed mode computes all query batches for a chromosome before assembling genome-wide output batches. Reflect that distinction in progress labels: direct mode can report current batch plus completed chromosomes; indexed mode should first report chromosome computation and then batch assembly. Neither should imply scientific outputs are published while they are still private fragments. See `_run_chromosomes()` and `_write_chromosomes()` in the sources above.

**Behavioral risk.** Preserve batching, reference traversal, output ownership, status ordering, and scientific counts. This calls for logging at the outer workflow boundary, not changing the calculation schedule. **Confidence: High. Audience: users and operators.**

### F4 — Medium: existing time throttling does not bound per-source events

**Evidence.** [PhaseProgress](../../src/ldsc/_progress.py) throttles ordinary `advance()` calls, but emits every start/completion and every `force=True` call. [_annotation_sources.py](../../src/ldsc/_annotation_sources.py), `_scan()`, opens a phase and forces a record for each annotation file. [_annotation_queries.py](../../src/ldsc/_annotation_queries.py), `prepare_bed_queries()`, forces a record for every BED source. Alias inference logs each distinct source/field mapping in [column_inference.py](../../src/ldsc/column_inference.py), `_log_if_inferred()`.

**Impact.** Many tiny inputs still generate thousands of INFO lines, despite the nominal 30-second interval. This affects `annotate`, direct `ldscore`, and shared annotation preparation in gene-index construction.

**Suggested direction.** Place the INFO phase around the whole source collection; increment completed file/row counts at bounded checkpoints. Keep per-file scope and detailed header mappings at DEBUG. Retain meaningful chromosome milestones where the number of chromosomes is bounded, and do not use `force=True` for every member of an arbitrarily large source list. Aggregate repeated compatibility notices per preparation run while retaining their scientific meaning and source-specific evidence.

**Validation performed.** A standalone harness loaded the actual `_progress.py`, fixed its clock at zero, and counted INFO records: 10,000 ordinary advances in one phase produced 2 records; 10,000 separate phases produced 20,000; 10,000 forced advances in one phase produced 10,002. These checks demonstrate event-volume behavior, not end-to-end workflow performance.

**Behavioral risk.** Preserve failure boundaries and structured event fields. A time interval is not a guarantee of a heartbeat during a blocking library call. **Confidence: High. Audience: users and operators.**

### F5 — Medium: batch `rg` repeats trait-level messages inside the pair loop

**Evidence.** `estimate_rg_pairs()` fits per-trait h2 and then every selected pair. `build_rg_dataset()` calls `_regression_genome_build()` for both traits in every pair, producing repeated provenance lines for coordinate-based inputs. All-pairs execution scales to `T*(T-1)/2` pairs; with 100 traits, those two provenance messages alone can produce 9,900 lines. There is no aggregate progress counter around pair fitting. The CLI already suppresses repeated global configuration banners, so those banners are not an additional CLI finding. See [regression_runner.py](../../src/ldsc/regression_runner.py), `estimate_rg_pairs()`, `build_rg_dataset()`, `_regression_genome_build()`, and `run_rg_from_args()`.

**Suggested direction.** Report trait and pair counts once, validate scientific compatibility as currently required, and report trait provenance once per run. Add aggregate per-trait h2 and pair-fitting progress. Keep pair-specific failures and actual pair-specific filtering outcomes identifiable. The full pair table and metadata retain details; do not hide a changed SNP population by treating all pair messages as identical.

**Behavioral risk.** Preserve per-pair checks, recoverable-pair failure behavior, allele handling, and metadata. Deduplicate presentation rather than removing validation. **Confidence: High. Audience: users.**

### F6 — Medium: unbounded lists also appear outside the query loop

**Evidence.** The quantitative-annotation interpretation warning joins every annotation name; `rg` logs every trait's prevalence and identity provenance in joined strings; `annotate` passes potentially large lists directly into `log_inputs()`. `_WorkflowLoggingContext._write_items()` prints arbitrary values verbatim, and input/output audit sections bypass `--log-level`. See [regression_runner.py](../../src/ldsc/regression_runner.py), `_log_quantitative_annotation_interpretation()`, `_log_effective_regression_identity()`, and `run_rg_from_args()`; [annotate_workflow.py](../../src/ldsc/annotate_workflow.py), `run_annotate()`; [_logging.py](../../src/ldsc/_logging.py), `_write_items()`.

**Suggested direction.** Use counts and bounded examples for large input groups and output families. Keep interpretation guidance visible even when its annotation-name list is abbreviated. Link to the applicable existing metadata, status table, or manifest for full membership. Preserve the exact `Call:` block for reproducibility; do not indiscriminately truncate every string in the formatter.

**Behavioral risk.** Compact only collections with a complete retained source of detail. Exact invocation, scientific provenance, and failure messages need deliberate handling. **Confidence: High. Audience: users.**

### F7 — Medium: munging warnings can repeat for every data chunk

**Evidence.** `parse_dat()` calls `filter_pvals()`, `filter_info()`, and `filter_frq()` per chunk; each warns whenever that chunk has invalid values. Routine drop totals are already aggregated after parsing. See [_kernel/sumstats_munger.py](../../src/ldsc/_kernel/sumstats_munger.py), those functions.

**Suggested direction.** Keep an immediate warning when an invalid-value category is first encountered; aggregate subsequent counts and emit a final category summary. Keep thresholding and filter results unchanged. Avoid a process-global suppression set: successive API runs need their own warnings.

An adjacent wording defect should be tracked separately: `filter_info()` tests values greater than `2.0` but warns about values outside `[0,1.5]`. Correcting that wording must follow the intended existing threshold; a logging cleanup must not silently choose a new numerical cutoff.

**Behavioral risk.** Warning occurrence/count changes are observable to callers and tests. Preserve the warning's scientific significance and distinguish invalid-value counts from rows removed by a user-requested threshold. **Confidence: High. Audience: users and developers.**

## Command-by-command disposition

| Command | Audit result and recommended scope |
| --- | --- |
| `ldscore`, direct | Highest priority: F2–F4 and compact annotation/output lists; parent-owned progress across batches and chromosomes. |
| `ldscore`, indexed | Highest priority: F2–F3; add computation progress before the batch assembly summaries. |
| `partitioned-h2` | Highest priority: F1 and F6; retain the useful regime/interpretation banner and exact failure diagnostics. |
| `annotate` | F4 and F6: shared annotation scans, forced BED-source updates, and potentially large input lists. The current standalone workflow does not call the LD-score per-row rejection logger. |
| `rg` | F5–F6: pair progress, repeated per-trait provenance, and large trait lists. |
| `munge-sumstats` | F7: repeated invalid-value warnings. Existing whole-run QC/drop summaries and the CLI completion summary are useful. |
| `build-gene-ldscore-index` | Shared per-file preparation can inherit F4. Existing chromosome start/completion and build metrics are useful; no per-gene success-log flood was found in the inspected build loop. |
| `build-r2-panel` | Mostly bounded chromosome-level messages. Preserve reference build/window/liftover decisions and counts; consolidate overlapping generic phase and chromosome messages only where they duplicate one another. |
| `h2` | Single fit: no analogous query-scale repetition found. Preserve effective identity, retained SNP count, and fit summary; shared collection formatting still applies if many baseline columns are reported. |
| `quantile-h2` | One INFO row per realized quantile. Useful at ordinary sizes; summarize counts and move individual bin details to DEBUG for large requested quantile counts. Full bounds/counts already belong in the result table. See [quantile_h2.py](../../src/ldsc/quantile_h2.py), `run_quantile_h2_from_args()`. |
| `query-r2` | No per-pair INFO loop found in the CLI path. Existing pair result statuses and metadata provide details; a compact final pair-status count would improve visibility. See [r2_query.py](../../src/ldsc/r2_query.py), `_write_query_r2_directory()` and `_query_r2_metadata()`. |
| `convert-ldsc2-ldscores` | No analogous routine per-SNP flood found. Preserve conversion issues in the existing issue table. See [legacy_ldscore_converter.py](../../src/ldsc/legacy_ldscore_converter.py), `LegacyLDScoreConverter.convert()`. |
| `convert-h2-scale` | Fixed input/output logging; no per-prevalence-point log loop found. See [h2_scale.py](../../src/ldsc/h2_scale.py), `convert_h2_scale()`. |
| `plot` | Fixed artifact input/output logging; no per-plotted-row log flood found. See [plotting/__init__.py](../../src/ldsc/plotting/__init__.py), `plot_result()`. |

## Proposed default log policy

1. **`INFO`: run-level narrative.** Report input counts, resolved scientific choices, each major stage's start/end, and periodic progress. Successful queries do not need individual log records or names.
2. **Progress: at most one periodic INFO update per active workflow stage every 30 seconds, plus stage start/end.** Reuse the existing monotonic progress mechanism at collection boundaries. Include completed/total with explicit units and elapsed wall time. Do not reset the cadence for every source, model, or small batch. Parallel coordinators can use timed waits to report unchanged completed counts with active workers; serial native calls can report only when control returns. Never increment completion just to show activity.
3. **`DEBUG`: troubleshooting detail.** Retain per-query preparation/estimation/staging times, per-source mappings, and exhaustive annotation names when requested. Keep context on those detailed records. Guard expensive formatting so INFO runs do not construct giant DEBUG strings.
4. **Warnings: concise, visible scientific consequences.** Summarize known repeat categories by stage/reason/role with exact counts and a bounded number of examples. Preserve the existing direct stderr notices for gene-list omissions/skips. Keep failures specific and retain full traceback/jackknife evidence at the current diagnostic boundary; the status tables do not replace that evidence.
5. **Complete diagnostics: existing tables and metadata.** Use the existing gene audits, query statuses, per-query manifest, regression tables, and annotation-type metadata. A new per-query timing file is optional and unnecessary for the initial cleanup.
6. **One owner per progress stream.** The parent/coordinator owns INFO progress; workers return outcomes and diagnostic evidence. Keep existing file/console routing and supported log-level flags. Reuse the current mechanisms instead of adding a generic text filter that guesses which lines to suppress.

This separation follows Python's documented distinction between operational INFO and diagnostic DEBUG, and its guidance to centralize file logging from multiple processes. The proposed 30-second cadence and message contents are recommendations for this package, not requirements of Python. Sources: [Python Logging HOWTO](https://docs.python.org/3/howto/logging.html#when-to-use-logging), [Python Logging Cookbook, logging to a single file from multiple processes](https://docs.python.org/3/howto/logging-cookbook.html#logging-to-a-single-file-from-multiple-processes).

### Illustrative regression log

The following is a proposed format, not a reported run outcome. Counts and timings other than the example query total are illustrative.

```text
Inputs: queries=18,788; baseline annotations=53; workers=8.
Phase start: shared regression preparation.
Phase complete: shared regression preparation; elapsed=18.4s.
Phase start: query fitting; total=18,788.
Query fitting: completed=420/18,788; successful=418; unestimable=2; failed=0; active=8; elapsed=10m30s.
...
Phase complete: query fitting; completed=18,788/18,788; successful=18,786; unestimable=2; failed=0; elapsed=4h12m.
WARNING: 2 queries were unestimable. See diagnostics/query_status.tsv and the failure details in this log.
Phase start: writing results.
Phase complete: writing results; elapsed=12.1s.
Outputs: partitioned_h2.tsv; per-query directory and manifest; diagnostics/metadata.json.
Finished ...
```

Individual error records would still identify the affected queries and stages when those errors occur. A fitting-stage completion is not a whole-run success; writing can still fail, and the existing workflow footer records the final outcome.

## Healthy patterns and documentation

Preserve the workflow header, exact invocation, final status/elapsed time, failure traceback, single authoritative log file, API exception propagation, and existing CLI routing. The shared liftover logger already provides a useful model: counts at INFO and example rows at DEBUG. See [_logging.py](../../src/ldsc/_logging.py), `_WorkflowLoggingContext`; [_kernel/liftover.py](../../src/ldsc/_kernel/liftover.py), `log_liftover_drop_report()`; [test_liftover_logging.py](../../tests/test_liftover_logging.py).

The established [troubleshooting document](../troubleshooting.md) already has command-specific validation, query-failure, and worker-failure guidance. No new troubleshooting surface is needed. A subsequent implementation should update [Workflow Logging](../current/workflow-logging.md), especially annotation preparation and partitioned-query diagnostics, because their current detailed INFO requirements are explicit. Update affected help/tutorial wording only where it promises the old detail level.

## Remediation boundary and verification

Recommended first scope: F1–F4 plus the unbounded annotation-name warning in F6. Then apply the same policy to `rg` and munging warnings. These were the initial audit recommendations; the implementation status above records the subsequent authorized work.

Validation performed here: read-only log counts and message-category measurements; local call-site and loop inspection; review of existing diagnostic contracts/tests; and the three deterministic progress-volume checks described in F4. No end-to-end workflow or full test suite was run during the original documentation-only audit; implementation verification is separate below. Findings of “no analogous flood found” refer to inspected paths, not an exhaustive runtime proof for every input combination.

Acceptance checks for implementation should cover:

- INFO volume remains bounded for many successful queries/files, both inline and in parallel, and the parent reports actual completion before a whole query batch finishes.
- DEBUG retains query timing and source detail; WARNING/ERROR levels still preserve their intended messages and lifecycle footer.
- Every scientific omission/failure remains visible, with exact totals and complete existing diagnostics. Interrupted/fatal runs do not claim successful fitting or publication.
- Batch sizes and worker counts do not change global totals; shared baseline annotations and early skipped inputs are counted once.
- Serial/parallel numerical outputs, artifact schemas, strict/continue policy, warning categories, and output cleanup remain unchanged.
- Existing message assertions are deliberately revised where the policy changes. Relevant coverage includes [test_logging_refactor.py](../../tests/test_logging_refactor.py), [test_staged_input_preflight.py](../../tests/test_staged_input_preflight.py), [test_gene_resolution_storage.py](../../tests/test_gene_resolution_storage.py), [test_gene_query_integration.py](../../tests/test_gene_query_integration.py), [test_annotation_preparation.py](../../tests/test_annotation_preparation.py), [test_partitioned_query_failures.py](../../tests/test_partitioned_query_failures.py), and [test_ldscore_parallelism.py](../../tests/test_ldscore_parallelism.py).


## Final implementation verification

- Python 3.13.13 in the existing `ldsc3-dev` environment: full pytest suite **1,882 passed, 1 skipped, 133 subtests passed**; 312 warnings. The skip is the missing-PyArrow dependency test because PyArrow is installed.
- Sequential unittest compatibility run: **1,001 tests, OK, 1 skipped**. It ran after pytest to avoid concurrent tests interfering with shared temporary BED resources.
- Both real CLI entry points, `ldsc --help` and `python -m ldsc --help`, succeeded and listed all 13 subcommands. Python compilation, `git diff --check`, and local documentation link targets passed.
- Behavior coverage includes INFO/DEBUG and inline/spawned regression batches; direct/indexed LD scores with two execution batches; collection-level progress for 40 annotation and 40 BED sources; row-identical munging results at different chunk sizes; rg and quantile totals; complete catalog issues on success and Gate A failure; removal of a stale catalog diagnostic after successful overwrite; and bounded trait/output lists. Existing numerical, scientific filtering, query-failure, artifact, and worker-cleanup checks pass.
- During validation, a new startup message initially accessed the annotation bundle before the existing missing-output argument check. The existing regression test exposed this; the check now precedes that access. INFO computation boundaries remain visible at the outer calculator, with inner batch details at DEBUG.
- Gene-index catalog failures now rely on their existing complete repair table and aggregate error instead of replaying every defect into warnings. No-analogue findings for the remaining conversion/plotting/reference workflows required no runtime change.
- No production Longleaf calculation was rerun and no remote logs or package files were changed. The local branch remains `restructure`. The user approved committing the completed cleanup together with the earlier negative-argument formatting changes in `_logging.py` and `test_logging_refactor.py`.
