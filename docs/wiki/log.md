# Wiki Update Log

Last updated on: 2026-09-15

## 2026-09-15 update | query regression failures

Updated the partitioned-h2 and gene-list functionality pages, guided tutorial, LDSC-SEG PC-gene walkthrough, and continuous-annotation guide to explain strict default behavior and the single `--continue-on-query-error` flag. Documented the attempted-query status table, successful-only scientific artifacts, detailed failure logs, fatal shared failures, and unchanged jackknife method. Aligned the Python docstrings, current reference docs, Markdown tutorials, and two notebook tutorials.

Source: `RegressionRunner.estimate_partitioned_h2_batch` in [regression_runner.py](../../src/ldsc/regression_runner.py), `PartitionedH2DirectoryWriter` in [outputs.py](../../src/ldsc/outputs.py), and [query failures and continuation](../current/partitioned-h2-results.md#query-failures-and-continuation).

Validation: both edited notebooks pass JSON and code-cell syntax checks. The synthetic cell-specific notebook executes top to bottom, including its API/CLI numerical comparison. The real-input partitioned notebook requires user-supplied inputs and was syntax-checked only. All 293 local file links across changed Markdown documents resolve; edited Python modules parse, and `git diff --check` passes.

## 2026-09-14 update | annotation preparation and concise logs

Reason: bring annotation, direct scoring, gene-list, and index-construction guidance up to date with group-aware chunks, metadata-only global identity cleanup, and the completed decision to retain serial preparation. Explain the three preparation milestones, one CM/MAF notice per file read, and one intentional-exclusion summary per gene set and role. Detailed row audits and workflow failure/publication semantics are unchanged. The partitioned and cell-specific tutorials and notebook prose use the same guidance.

Source: `prepare_annotation_sources()` in [_annotation_sources.py](../../src/ldsc/_annotation_sources.py), `_log_gene_list_rejections()` in [query_annotations.py](../../src/ldsc/query_annotations.py), [workflow logging](../current/workflow-logging.md#annotation-preparation), and the [closed parallelism decision](../audits/annotation-memory/preparation-parallelism.md#decision). Recent implementation commits reviewed: `cf1dfac`, `ce2bce9`, `b91acbd`, and `6c90636`.

Validation: the full pytest suite passed 1,664 tests with one skip and 132 passing subtests; unittest compatibility ran 1,000 tests successfully with one skip. Updated Python files parse, public signatures and docstrings were inspected, and CLI/root and index-build help succeed. All 233 local links and anchors across the 16 changed documentation files resolve. The two notebook changes affect Markdown only; code cells, outputs, and execution metadata are unchanged. No documentation builder is configured, and notebooks were not re-executed for prose-only changes. Work stayed local and did not modify HPC jobs.

## 2026-09-14 update | sequential query batches

Reason: align scoring, gene-list, index-reuse, guided-analysis, and continuous-annotation guidance with sequential output batching. Document numbered genome-wide query files, the mandatory ordered manifest, explicit uncached Python reads, the small zero-write Python route, and indexed chromosome workers capped at chromosome count. Correct the earlier quantitative-count limitation: normalized float32 annotations now use float64 count accumulation with existing tolerances. Batch resumption remains a possible future feature.

Source: `LDScoreCalculator.run`, `LDScoreSource.read_queries`, `indexed_results` in [_indexed_ldscore_batches.py](../../src/ldsc/_indexed_ldscore_batches.py), `annotation_statistics` in [_kernel/overlap.py](../../src/ldsc/_kernel/overlap.py), and the [numerical verification report](../audits/annotation-memory/sequential-query-batches.md). Historical entries below describe their original verification state.

Validation: current CLI parsing succeeds for all 40 complete Bash examples across the updated wiki/tutorial pages; one abbreviated example is intentionally skipped. Local links and heading anchors resolve. Both updated tutorial notebooks execute successfully with their built-in or substituted small fixtures, including cross-file reads and an API/CLI numerical comparison. See [documentation verification](../audits/annotation-memory/sequential-query-batches.md#documentation-verification).

## 2026-09-13 update | wiki consistency corrections

Reason: fix the guided tutorial's MDD-to-ADHD variable leak; correct identity defaults, baseline requirements, path expansion, private annotation disk use, index destinations, genetic-coordinate handling, legacy compatibility, and removed-flag wording. Replace branch-specific wiki URLs with relative links, replace unsupported resource estimates with benchmark context, complete empty navigation pages, and add baseline-only and saved-h2 conversion examples. Existing gene-list edits were preserved.

Validation: all 36 complete Bash-fenced LDSC examples parse with the current CLI; all local wiki file links and heading anchors resolve; all pages contain content and update dates. The actual guided-tutorial Bash blocks, with scientific commands replaced by an argument recorder, select MDD after the three-trait preparation sequence and when Analysis 3 starts with stale ADHD variables. Both index variants resolve matching build/monitor/Python paths. The Parquet preview reads 1,000 rows from one selected column in a synthetic file. `git diff --check` passed. No production datasets or HPC jobs were used.

Source: [audit and remediation record](../audits/2026-09-13-wiki-consistency.md), current CLI parsers, `build_query_shards`, `_resolve_genetic_map`, `normalize_annotation_chunk`, shared output writers, and the current path and result-schema documentation. The separate quantitative-count limitation remains open and is now linked from the continuous-annotation page.

## 2026-09-11 update | build-r2-panel, query-r2, examine-parquet-and-npz

Reason: replace the build/query placeholders with current usage and artifact contracts, rename the build page, explain SIGN_R versus sign_r and required sidecars, and correct inspection examples for quantized R².

Source: `src/ldsc/ref_panel_builder.py`, `src/ldsc/_kernel/ref_panel_builder.py`, `src/ldsc/r2_query.py`, `src/ldsc/outputs.py`, and the current Parquet and query format documentation.

## 2026-09-11 update | munge-sumstats, guided-tutorial

Reason: remove obsolete HM3 switches, explain default restriction and automatic mapping with explicit chain precedence, document stdout/log count summaries, require an output directory for inference previews, and use trait-named data files in downstream commands.

Source: `src/ldsc/sumstats_munger.py` (`build_parser`, `run_munge_sumstats_from_args`, `_render_munge_summary`), `src/ldsc/config.py` (`MungeConfig`), and [the current munging guide](../current/munge-sumstats.md#liftover-rules).

## 2026-09-14 update | guided-tutorial, LDSC-SEG-PC-genes

Reason: record the completed 1,000-query chromosome-22 comparison at generation widths 1,000 and 100. Explain the measured scratch/runtime trade-off and essentially unchanged peak RAM, and distinguish these pilots from full-genome projections and the gene-list tutorial's different padding.

Source: [batch benchmark audit](../audits/annotation-memory/sequential-query-batches.md#resource-measurements-and-limits) and its [portable measured/projection evidence](../audits/annotation-memory/sequential-query-batch-benchmarks.json), measured at package revision `cf1dfac`.

## 2026-09-15 update | gene-list resolution flag

Reason: align gene-list LD-score, standalone annotation, and continuous-annotation guidance with the value-free `--allow-unresolved-genes` flag, strict omission default, and rejection of the retired option without a compatibility alias. Explain unchanged Python policy strings and metadata, and reproduction of the original resolution policy for quantile reconstruction. Update the corresponding docstrings, current references, and LD-score tutorial.

Source: `build_parser()` and `run_ldscore()` in [ldscore_calculator.py](../../src/ldsc/ldscore_calculator.py), `add_annotate_arguments()` in [annotation_builder.py](../../src/ldsc/annotation_builder.py), `run_annotate()` in [annotate_workflow.py](../../src/ldsc/annotate_workflow.py), and `add_quantile_h2_arguments()` in [quantile_h2.py](../../src/ldsc/quantile_h2.py).


## 2026-09-15 update | parallel partitioned query workers

Reason: document whole-query process workers and the shared `--threads` parser, validator, and resolver across the partitioned, gene-set, guided, LDSC-SEG, and continuous-annotation pages. Distinguish query loading width from worker count and LD-score generation files from regression scratch. Explain CPU affinity/negative offsets, work caps, positive scheduler requests, one native thread per parallel worker, memory ownership, deterministic ordering, and fatal process failures. Update API/output docstrings, argument/metadata inventories, data flow, logging, feature navigation, Markdown tutorials, and both notebook tutorials.

Source: `RegressionRunner.estimate_partitioned_h2_batch()` in [regression_runner.py](../../src/ldsc/regression_runner.py), the shared resolver in [_parallelism.py](../../src/ldsc/_parallelism.py), [_partitioned_h2_parallel.py](../../src/ldsc/_partitioned_h2_parallel.py), and [query-worker policy](../current/regression-configuration.md#43-query-workers-and-memory). Local runtime and memory guidance links to the existing [benchmark report](../audits/2026-09-15_partitioned-query-workers.md).

Validation: the synthetic cell-specific notebook executes top to bottom, comparing the inline API summary and query-status ordering against a real two-worker CLI run and checking requested/effective/native-thread metadata. Both notebooks pass JSON and code-cell syntax checks; the real-input partitioned notebook remains an input template and was not executed without its user-supplied data. Complete CLI examples parse, local file links and heading anchors resolve, update dates and `git diff --check` pass. An AST comparison confirms that the three source modules edited in this documentation pass changed only in docstrings. The previously completed full pytest/unittest and 1/2/4-worker benchmark results remain recorded in the implementation report.
