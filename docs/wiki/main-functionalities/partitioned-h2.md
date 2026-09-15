# Partitioned heritability

Last updated on: 2026-09-15

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

`ldscore --query-batch-size 1000` bounds the active query calculation batch and defaults to 1000. Completed batches are written and their LD tables released. Multiple batches produce numbered query Parquet files, each covering the full regression SNP universe; root `metadata.json.query_batches` records their columns and chromosome row groups. Direct and indexed `--threads` control chromosome workers, default to 1, and are capped at the chromosome count. Larger worker counts need more memory. See [LD-score memory controls](ldscore.md#memory-for-many-pathways).

After LD scoring, fit the pathways with:

```bash
ldsc partitioned-h2 \
  --ldscore-dir results/pathway_ldscores \
  --sumstats-file results/trait/sumstats.parquet \
  --query-batch-size 1000 \
  --output-dir results/pathway_enrichment
```

The command prepares shared trait/baseline alignment once and reads query columns in batches. Regression batch width is independent of generation batch width; one read can span several saved query files. Every pathway still fits its complete retained genome-wide SNP set with model-specific filtering, weights, and jackknife calculations. A smaller batch such as `--query-batch-size 100` reduces active query memory without changing the separate models. It can increase I/O or reduce multiplication throughput; the best value depends on the workload. Sources: `RegressionRunner.estimate_partitioned_h2_batch` in [regression_runner.py](../../../src/ldsc/regression_runner.py) and `LDScoreSource.read_queries` in [ldscore_source.py](../../../src/ldsc/ldscore_source.py).

The aggregate `partitioned_h2.tsv` keeps the requested final ordering. Per-query category tables, coefficient delete values, and metadata are written as fits finish, staged privately until successful sorted publication under `diagnostics/query_annotations/`. Python `estimate_partitioned_h2_batch(..., output_dir=...)` returns the aggregate summary and persistent `per_query_artifacts` paths, rather than every detailed table in memory.

New temporary files stay below the selected output directory and are cleaned on handled completion or failure. Existing overwrite/failure-marker behavior is preserved. Released allocations can remain reserved by Python/NumPy, so RSS need not fall after every chromosome. See the [memory design](../../current/annotation-memory-design.md) for ownership, indexing, and exact-quantile details.

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

Shared input loading, output failures, and interrupts remain fatal. A scan with no successful fits also fails, retaining its query diagnostics. Baseline-only regressions have no independent query to skip. This flag is separate from `ldscore --gene-list-resolution-policy`, which controls gene identifier resolution before regression.

For Python, use `runner.estimate_partitioned_h2_batch(..., continue_on_query_error=True)` and inspect `result.query_status`. The method logs through the configured LDSC logger; the CLI owns the log-file handler. Sources: `RegressionRunner.estimate_partitioned_h2_batch()` in [regression_runner.py](../../../src/ldsc/regression_runner.py) and the [result contract](../../current/partitioned-h2-results.md#query-failures-and-continuation).
