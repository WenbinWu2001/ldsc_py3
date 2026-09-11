# Partitioned heritability

Last updated on: 2026-09-10

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

## Testing enrichment for a large batch of pathways

The memory optimization is designed for testing many pathways at once—for example, 1,000 pathways in one run, with each pathway tested separately against the same baseline categories. Put one pathway in each gene-list or BED source, or one pathway in each prebuilt query annotation column.

`ldscore --query-batch-size 1000` bounds the active query projection batch; this is the default. `--threads 1` processes chromosomes sequentially and releases their working arrays before the next chromosome. Larger worker counts process chromosomes concurrently and need more memory. The final HM3 LD tables remain aggregate Parquet files and may stay in memory; they are never split into public chromosome LD tables.

After LD scoring, fit the pathways with:

```bash
ldsc partitioned-h2 \
  --ldscore-dir results/pathway_ldscores \
  --sumstats-file results/trait/sumstats.parquet \
  --query-batch-size 1000 \
  --output-dir results/pathway_enrichment
```

The command prepares shared trait/baseline alignment once and reads query columns in batches. Every pathway still fits its complete retained genome-wide SNP set with model-specific filtering, weights, and jackknife calculations. A smaller batch such as `--query-batch-size 100` reduces active query memory without changing the separate models. It can increase I/O or reduce multiplication throughput; the best value depends on the workload.

The aggregate `partitioned_h2.tsv` keeps the requested final ordering. Per-query category tables, coefficient delete values, and metadata are written as fits finish, staged privately until successful sorted publication under `diagnostics/query_annotations/`. Python `estimate_partitioned_h2_batch(..., output_dir=...)` returns the aggregate summary and persistent `per_query_artifacts` paths, rather than every detailed table in memory.

New temporary files stay below the selected output directory and are cleaned on handled completion or failure. Existing overwrite/failure-marker behavior is preserved. Released allocations can remain reserved by Python/NumPy, so RSS need not fall after every chromosome. See the [memory design](../../current/annotation-memory-design.md) for ownership, indexing, and exact-quantile details.
