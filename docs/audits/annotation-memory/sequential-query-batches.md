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

No new measured RAM or runtime improvement is claimed. The prepared [benchmark runner](../../../benchmarks/sequential_query_batches.py) replays a preserved CLI command into a fresh output directory, sampling simultaneous process-tree RSS and private/persistent logical disk bytes separately. The local chromosome-22, 1,000-query comparison at generation widths 1,000 and 100 remains pending explicit execution approval. Earlier resource reports predate this implementation.

At fixed query width and worker count, the implementation bounds live query LD values by the active batch, while compact source metadata, counts, and baseline-by-query overlaps still grow with total query count. Direct batching repeats reference preparation and numerical traversal. Increasing chromosome workers increases the number of simultaneously resident chromosome states. The full 18,000-query, 22-chromosome workload has not been measured here. Resume from the last completed batch remains a documented future feature.
