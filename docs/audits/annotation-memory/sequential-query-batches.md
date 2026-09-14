# Sequential query batch verification

Last updated on: 2026-09-14

The [confirmed September 14 design](../../current/annotation-memory-decisions.md#confirmed-batch-integration-2026-09-14) is implemented locally on `restructure`. Direct calculation now prepares, computes, privately writes, and releases each query batch. Indexed generation retains one operator per chromosome worker, writes float64 batch fragments, and assembles them in deterministic chromosome order. Both writing routes return `LDScoreSource` without retained query LD tables. Single-batch prepared Python inputs can return a materialized result and complete diagnostics without filesystem writes.

## Numerical evidence

No existing numerical tolerance was loosened. Annotation normalization and final Parquet LD values remain float32; LD accumulation and overlap products remain float64. Annotation sums now also accumulate in float64, correcting the previously recorded batch/layout-dependent quantitative count error.

| Check | Evidence |
| --- | --- |
| Quantitative counts | Independent `math.fsum` of normalized float32 annotations, including signed cancellation and high-offset values; 6,000 rows, 1,001 columns, C/F layouts, query widths 1 and 1,000; `rtol=1e-6`, `atol=1e-8` |
| Direct LD calculation | Two-chromosome R2 fixture with binary and signed quantitative annotations; widths 1, 2, and 10, with one/two workers; LD values, counts, baseline, and overlaps agree with a materialized single batch at existing tolerances |
| Indexed LD calculation | Two chromosome operators, three overlapping gene sets, widths 1, 2, and 1,000; known union scores and all/common counts/overlaps preserved; worker requests 1, 2, -1, -2, and 10 agree, with requests capped at two chromosomes |
| Direct/indexed parity | Existing real computation integration tests retain their original control, focal-score, count, overlap, and diagnostic comparisons |
| Downstream regression | Three models loaded from single versus multiple query files; generation widths 1/2 and regression width 2; summaries, full category tables, and coefficient jackknife delete values agree at `rtol=1e-10`, `atol=1e-12` |

The main new checks are in [test_ldscore_query_batches.py](../../../tests/test_ldscore_query_batches.py), [test_ldscore_streaming.py](../../../tests/test_ldscore_streaming.py), [test_gene_index_streaming.py](../../../tests/test_gene_index_streaming.py), and [test_regression_streaming.py](../../../tests/test_regression_streaming.py).

## Ownership, files, and compatibility

Tests verify zero writes for prepared single-batch calculation, retained in-memory diagnostics after input closure, global duplicate cleanup, and rejection of incomplete allele metadata. Lifetime checks assert that completed query tables, chromosome operators, prepared BED batches, and prebuilt-query diagnostic scratch are released at their intended boundaries.

Publication checks cover private query files before complete success, producer closure and scratch cleanup on write failure, selected reads spanning files in requested order, stale numbered/single-file cleanup after overwrite, preservation of unrelated files, and rejection of directories missing the current `query_batches` manifest. Complete source, schema, allele, SNP-universe, count, overlap, regression, and legacy-conversion tests remain active. Existing old-format test fixtures were migrated to the current directory contract.

Full pytest verification: **1,634 passed, 1 expected skip, 132 subtests passed**. The unittest compatibility command also completed successfully: **1,000 tests run, one expected skip**. `git diff --check` passed. Both the installed `ldsc --help` entry point and `python -m ldsc ldscore --help` succeeded. The latter documents indexed/direct chromosome workers, chromosome-count caps, and query output batching.

## Documentation verification

Public calculator, configuration, source-reader, indexed-generation, and batch-writer docstrings now describe the actual return types, query-file layout, worker controls, restrictions, diagnostics, and zero-write boundary. README, current design/schema/workflow documents, the affected wiki pages, and Markdown tutorials use the same contract. The SNP-universe reference explicitly supersedes the old one-traversal-across-all-queries rule while retaining contributor and output-row definitions.

Both updated notebooks executed their code cells in order in fresh temporary directories. The [cell-specific notebook](../../../tutorials/cell-specific-ldsc.ipynb) used its built-in 20-SNP fixture, wrote two query files, checked cross-file requested-column order and exact stored values, and compared API/CLI regression summaries at `rtol=1e-10`, `atol=1e-12`. The [partitioned notebook](../../../tutorials/partitioned-ldsc.ipynb) remains a user-input template; verification substituted an 80-SNP canonical R²/annotation/sumstats fixture and generation width 1 into its parameter cells, then ran the workflow cells unchanged. It produced two manifest-declared query files, preserved all expected rows, fitted both queries, and cleaned private staging. Execution used ordinary Python code-cell evaluation because the development environment has no notebook execution library; notebook JSON and all code cells were also validated.

Forty complete Bash LDSC examples across the updated wiki/tutorial pages parsed with the current CLI; one explicitly abbreviated example was skipped. Local Markdown/notebook links and heading anchors resolve, all edited Markdown has the current update date, changed Python files compile, and `git diff --check` passes. These documentation checks add no new resource-performance evidence.

## Resource measurements and limits

Both approved local chromosome-22 benchmarks completed at committed package revision `cf1dfac052728df6ef4d919a09ff814f41e95b77`. The [benchmark runner](../../../benchmarks/sequential_query_batches.py) replayed the same preserved CLI command into fresh output directories, sampling simultaneous process-tree RSS and private/persistent logical disk bytes separately. Each run used 1,000 chr22-supported gene lists, 53 binary baselines, 489 individuals, 141,123 reference SNPs, 17,380 output SNPs, a 1-cM window, and 35-kb gene padding. The requested 22 workers were capped at one chromosome; BLAS pools used one thread.

| Generation width | Query files | Wall time | Peak sampled process-tree RSS | Peak private logical disk | Persistent output |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1,000 | 1 | 163.407 s | 1.166 GiB | 0.568 GiB | 17.21 MiB |
| 100 | 10 | 310.996 s | 1.164 GiB | 0.101 GiB | 19.46 MiB |

Width 100 used **82.2% less peak private scratch** and took **1.903 times as long**. Peak RAM was effectively unchanged: sampled RSS differed by about 0.1%, while the independent single-process high-water marks differed in the opposite direction (1.169 versus 1.173 GiB). Both runs cleaned all private scratch. This small-chromosome pilot demonstrates the scratch/runtime trade-off; the larger-workload memory benefit remains a projection supported by ownership tests.

Every one of the **17,380,000 query LD values per run** matched the preserved pre-change width-1,000 result exactly. Baseline tables, row identities/order, all/common counts, and overlaps also matched exactly. Requested-order reads across numbered files passed, and loaded sources retained no query or chromosome-result tables. The [portable evidence snapshot](sequential-query-batch-benchmarks.json) records measured bytes, validation, dimensions, and conditional projections. The full local report and reproduction artifacts remain in the workspace at `runtime_estimates/2026-09-14_gene_sets_chr22/REPORT.md`, outside this package repository; original measurements and unsuccessful attempts are preserved there.

For **18,000 eligible queries, 53 baselines, chromosomes 1–22, and 22 chromosome workers**, the revised width-1,000 planning estimate is **8–16 hours and 32–48 GiB peak RAM**, with **64 GiB RAM and 150–200 GiB free output space** provisioned. Width 100 has a **16–32-hour** runtime projection. Runtime calibrates the earlier repeated-pipeline phase model using these new pilot durations; memory uses active-batch dimensions and the earlier width-1,000 aggregation copy probe. These are conditional planning ranges, not measured bounds or confidence intervals. The previous 400–600 GiB LD-score RAM estimate assumed all-query aggregation and is superseded. The earlier regression estimate remains **10–20 hours per trait and 32–48 GiB RAM**; no new regression timing was run.

Each width was measured once, in order 1,000 then 100, on the same 16-GiB Apple M1 Pro. Cache state, thermals, and host load were not controlled. The new RSS/disk sampler was absent from the older timing run, so the new-versus-old runtime ratio does not isolate implementation overhead. Full-genome list density, worker overhead, storage behavior, output filtering, and compression remain unmeasured. The initial new launch failed missing-input preflight because its working directory did not resolve the saved relative resource paths; it is excluded from all projections. No package code changed between successful runs.

At fixed query width and worker count, the implementation bounds live query LD values by the active batch, while compact source metadata, counts, and baseline-by-query overlaps still grow with total query count. Direct batching repeats reference preparation and numerical traversal. Increasing chromosome workers increases the number of simultaneously resident chromosome states. The full 18,000-query, 22-chromosome workload has not been measured here. Resume from the last completed batch remains a documented future feature.
