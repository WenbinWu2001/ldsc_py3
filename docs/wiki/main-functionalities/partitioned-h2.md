# Partitioned heritability

Last updated on: 2026-10-26

`ldsc partitioned-h2` reads a canonical LD-score directory and tests how its
annotations contribute to SNP heritability.

- With baseline annotations only, it reports a functional partitioning model.
- With focal query annotations, it fits one baseline-plus-one-query model per
  focal annotation.
- Gene-list analyses add no gene control by default. Supplying
  `--control-gene-list-file` during LD-score calculation adds a fixed
  `gene_control` baseline annotation from that one-column file.

For a complete gene-set workflow—from an explicit complete index through
LD-score calculation and regression—follow [Calculate LD scores for gene lists
with an index](ldscore-from-gene-list.md). For the general end-to-end workflow and output
interpretation, see [the guided tutorial](../guided-tutorial.md) and the
[partitioned LDSC technical reference](../../current/partitioned-ldsc-workflow.md).


LD scoring accepts prebuilt annotation columns, BED intervals, or gene lists as mutually exclusive query routes. BED files can be prepared by users or obtained from an external resource in the required format, and support optional padding. Gene lists require a coordinate catalog for direct projection. See [annotation preparation](../utility-functionalities/annotate.md).

## Baseline-only functional partitioning

Coordinate-family summary statistics and LD-score artifacts must agree on their genome build before fitting.
This command has no build override; supply matching upstream artifacts if the grouped build check fails.
See the [genome-build policy](../../current/snp-identifier-genome-build-defaults.md#genome-build-policy), including the separate legacy LDSC2 text projection path.

To fit all supplied baseline categories jointly, use a canonical LD-score directory built with baseline annotations and no query columns:

```bash
ldsc partitioned-h2 \
  --ldscore-dir results/baseline_ldscores \
  --sumstats-file results/trait/sumstats.parquet \
  --output-dir results/functional_partitioning
```

The root `partitioned_h2.tsv` reports one row per baseline category. For binary categories, interpret `enrichment` with its two-sided `enrichment_p`. For focal query models, prefer the conditional `coefficient` and one-sided `coefficient_p`; each query is fitted separately with the baseline. Quantitative annotations require [continuous-annotation interpretation](../continuous-annotation-partitioned-ldsc.md). The [result schema](../../current/partitioned-h2-results.md) explains the saved columns, complete-model diagnostics, and distinctions between category and total heritability. Source: `run_partitioned_h2_from_args` in [regression_runner.py](../../../src/ldsc/regression_runner.py).

## Testing enrichment for a large batch of pathways

The memory optimization is designed for testing many pathways at once—for example, 1,000 pathways in one run, with each pathway tested separately against the same baseline categories. Put one pathway in each gene-list or BED source, or one pathway in each prebuilt query annotation column.

`ldscore --query-batch-size 100` bounds the active query calculation batch and defaults to 100. Completed batches are written and their LD tables released. Multiple batches produce numbered query Parquet files, each covering the full regression SNP universe; root `metadata.json.query_batches` records their columns and chromosome row groups. Direct and indexed `--workers` control chromosome workers, default to 1, and are capped at the chromosome count. Larger worker counts need more memory. See [LD-score memory controls](ldscore.md#memory-for-many-pathways).

After LD scoring, fit the pathways with:

```bash
ldsc partitioned-h2 \
  --ldscore-dir results/pathway_ldscores \
  --sumstats-file results/trait/sumstats.parquet \
  --query-batch-size 100 \
  --output-dir results/pathway_enrichment
```

The command prepares shared trait/baseline alignment once and reads query columns in batches.
Regression batch width is independent of generation batch width; one read can span several saved query files.
Every pathway still fits its complete retained genome-wide SNP set with model-specific filtering, weights, and jackknife calculations.
A smaller batch such as `--query-batch-size 100` reduces loaded query values and private mapping space, but also limits how many query workers can be active.
It does not reduce the SNP rows or baseline categories in each active fit.
Loading and mapping overhead can change with batch width; the best value depends on the workload.
Sources: `RegressionRunner.estimate_partitioned_h2_batch` in [regression_runner.py](../../../src/ldsc/regression_runner.py) and `LDScoreSource.read_queries` in [ldscore_source.py](../../../src/ldsc/ldscore_source.py).

Current regression automatically borrows aligned baseline/query column views, applies row filtering before materializing the numerical design, and reuses owned numerical buffers where dtype promotion permits.
Shared worker inputs remain read-only.
These changes reduce copying for regression over both direct and indexed LD-score inputs, without changing float64 calculations, SNP/query ordering, counts, overlaps, or error policies.
No new option is needed; batch and worker defaults are unchanged.
See the [regression memory contract](../../current/annotation-memory-design.md) and `_dataset_from_query_values` / `_fit_h2_dataset` in [regression_runner.py](../../../src/ldsc/regression_runner.py).

The aggregate `partitioned_h2.tsv` keeps the requested final ordering.
Per-query category tables, coefficient delete values, and metadata are written as fits finish, staged privately until successful sorted publication under `diagnostics/query_annotations/`.
Publication moves disposable staged detail files into place and stores run-wide provenance once, reducing repeated writing and serialization.
Python `estimate_partitioned_h2_batch(..., output_dir=...)` returns the aggregate summary and persistent `per_query_artifacts` paths, rather than every detailed table in memory.
Republishing those persistent artifacts copies them and preserves the original results.
Source: `PartitionedH2DirectoryWriter` in [outputs.py](../../../src/ldsc/outputs.py).

The saved root summary appends `query_result_dir`, the exact per-query folder relative to the output directory; `category` retains the original gene-set name. Use this path to find detailed results after re-sorting the table. Baseline-only summaries use `NaN` because no per-query folders exist. See [the output contract](../../current/partitioned-h2-results.md#partitioned_h2tsv).

New temporary files stay below the selected output directory and are cleaned on handled completion or failure. Existing overwrite/failure-marker behavior is preserved. Released allocations can remain reserved by Python/NumPy, so RSS need not fall after every fit. See the [memory design](../../current/annotation-memory-design.md) for ownership, indexing, and exact-quantile details.

## Keep fitted results together

Copy the complete `partitioned-h2` output directory when moving or sharing fitted queries.
Run-wide provenance is stored in `diagnostics/metadata.json`; each query's `metadata.json` contains fit-specific fields and `"global_metadata": "../../metadata.json"`, resolved relative to that query folder.
A query folder alone is therefore incomplete for downstream model loading.
The `files` entries in query metadata remain relative to the aggregate result root.

For `quantile-h2`, still select one successful `diagnostics/query_annotations/<query>/` folder inside the complete result tree.
The current `load_fitted_partitioned_model()` reader resolves shared provenance automatically, with query-specific values taking precedence.
If the shared metadata is missing, restore the enclosing `diagnostics/metadata.json` from the same fitted run and preserve its directory layout.
See [shared-metadata troubleshooting](../../troubleshooting.md#quantile-h2-shared-fitted-model-metadata-is-missing), the [result contract](../../current/partitioned-h2-results.md), and the reader in [quantile_h2.py](../../../src/ldsc/quantile_h2.py).

## Memory and runtime measurements

Longleaf pilots measured the two regression changes separately before integration in `7a478e6`.
The following values are phase medians, comparing the earlier implementation with the optimized implementation:

| Phase and workload | Wall time, s | CPU time, s | Measured resource change |
| --- | ---: | ---: | --- |
| Model fitting: 60 queries, 20 workers | 60.83 → 31.74 | 898.84 → 401.08 | Peak anonymous memory: 37.25 → 27.13 GiB |
| Publication: 100 completed fits | 5.456 → 1.917 | 3.355 → 0.650 | Logical writes: 139.856 → 1.667 MiB; anonymous memory roughly unchanged |

The model comparison used two baseline and two optimized runs on the same node and BLAS backend, retaining 1,127,525 aligned SNP rows and 53 baseline columns per fit.
It used direct LD-score inputs; indexed-input regression uses the same fitting code.
The publication comparison used four baseline and two combined-optimization runs with the full 18,841-entry annotation dictionary.
Scientific results agreed exactly within each tested BLAS backend; cross-backend differences were also present between unchanged baselines.
Both changes reduced measured phase wall time and CPU time, but these pilots do not establish a full production-job speedup or a smaller safe RAM allocation.
In particular, reduced publication writes do not imply a proportional reduction in anonymous memory or charged file cache.
Source: the workspace [pilot report](../../../../diagnostics/20261014_memory_pilots/report.md), including workload boundaries, repetitions, cache measurements, and numerical checks.
The report link requires the enclosing research workspace; the table above summarizes its findings for standalone package readers.

## Fit queries concurrently

Use `--workers` when several independent queries are ready to fit. For example, with four allocated CPUs:

```bash
ldsc partitioned-h2 \
  --ldscore-dir results/pathway_ldscores \
  --sumstats-file results/trait/sumstats.parquet \
  --workers 4 \
  --query-batch-size 100 \
  --output-dir results/pathway_parallel
```

This loads at most 100 query columns at once and fits up to four complete baseline-plus-query models concurrently. Each process performs its own filtering, weights, regression, jackknife, summaries, and private output staging. SNP ordering and jackknife blocks follow the serial calculation; a model is never divided into chromosome fits.

| Setting | Rule |
| --- | --- |
| `--workers 1` | Default: one managed child, with 1 BLAS thread. |
| Positive `N` | Reject above detected CPU count, then cap by available models and batch width. |
| Negative `-k` | `max(1, joblib.cpu_count() + 1 - k)`, then apply the work cap. |
| CPU detection | Detected available CPUs, not idle CPUs or the requested scheduler allocation. |
| One query or batch width 1 | One managed child; no concurrent models. |

Every worker uses 1 BLAS thread, including serial and baseline-only calculation.
OpenMP and Arrow CPU/I/O limits are separately one; unsupported numerical limits are recorded and computation continues.
Zero, Boolean, and non-integer requests are invalid.
See [the execution policy](../../current/parallel-execution.md) for allocation limits, process startup, and diagnostics.

Workers borrow read-only maps of the shared baseline/trait arrays and current query batch from private output scratch.
More workers still need more RAM for private filtered model matrices and estimators, and the maps need disk space.
Query completion order does not change the requested status order or the stable summary/manifest order.
The INFO log reports requested/effective workers and phase timings; run-wide worker counts are stored in root metadata and available through each successful query's shared metadata reference.
A historical local 1,000-query chromosome-22 benchmark, before the current copy and publication improvements, measured 74.17, 47.63, and 30.47 seconds with 1, 2, and 4 workers, respectively, including preparation and publication.
Its sampled peak process-tree RSS was 472, 1,048, and 1,507 MiB; these are not full-genome RAM estimates.

The Python option is `runner.estimate_partitioned_h2_batch(..., workers=4, query_batch_size=100)`. Put parallel calls in an importable script under `if __name__ == "__main__":`; from a notebook, use the CLI in a subprocess for parallel fitting. Both error policies below apply to ordinary query-model exceptions. Worker death, transport errors, output errors, and interrupts abort under either policy, identify affected queries, and trigger child-process cleanup. Sources: `RegressionRunner.estimate_partitioned_h2_batch()` in [regression_runner.py](../../../src/ldsc/regression_runner.py), the shared resolver in [_parallelism.py](../../../src/ldsc/_parallelism.py), and process ownership in [_partitioned_h2_parallel.py](../../../src/ldsc/_partitioned_h2_parallel.py).

## Choose what happens when a query fails

One flag controls this behavior:

| Regression command | Query-error behavior |
| --- | --- |
| Without `--continue-on-query-error` (default) | Collect per-query errors and fail the run without publishing new scientific results if any query fails. |
| With `--continue-on-query-error` | Skip and mark failed queries, finish the scan, and publish successful fits. |

For a scan where successful fits should survive query errors:

```bash
ldsc partitioned-h2 \
  --ldscore-dir results/pathway_ldscores \
  --sumstats-file results/trait/sumstats.parquet \
  --continue-on-query-error \
  --output-dir results/pathway_enrichment
```

Each completed query scan writes `diagnostics/query_status.tsv` in requested order. Its columns are `query_annotation`, `status`, `stage`, `error_type`, and `error_message`. Status is `success`, `unestimable` for singular jackknife deletions, or `failed` for other exceptions. Failed queries have no scientific summary row or result folder; their absence is not a null association. Keep the full status table when accounting for the requested testing family.

The CLI log, `diagnostics/partitioned-h2.log`, records query names, stages, exceptions, and tracebacks. Singular delete-block failures also report block numbers, retained-row intervals, available genomic spans, normal-matrix ranks, and support counts. A query can fit on all SNPs but become unidentifiable when its only supported block is deleted. The flag does not change filtering, weights, blocks, or the solver.

Shared input loading, output failures, and interrupts remain fatal. A scan with no successful fits also fails, retaining its query diagnostics. Baseline-only regressions have no independent query to skip. This flag is separate from `ldscore --allow-unresolved-genes`, which controls gene identifier resolution before regression.

For Python, use `runner.estimate_partitioned_h2_batch(..., continue_on_query_error=True)` and inspect `result.query_status`. The method logs through the configured LDSC logger; the CLI owns the log-file handler. Sources: `RegressionRunner.estimate_partitioned_h2_batch()` in [regression_runner.py](../../../src/ldsc/regression_runner.py) and the [result contract](../../current/partitioned-h2-results.md#query-failures-and-continuation).
