# Calculate LD scores

Last updated on: 2026-09-13

`ldsc ldscore` computes a reusable LD-score directory from PLINK genotypes or a precomputed R² panel. Regression consumes that directory. Ordinary unpartitioned scoring uses a synthetic all-ones `base` annotation:

```bash
ldsc ldscore \
  --plink-prefix /data/plink/panel_chr22 \
  --snp-identifier chr_pos \
  --genome-build hg19 \
  --ld-wind-kb 1000 \
  --output-dir results/unpartitioned_ldscores
```

This illustrative command selects one concrete chromosome; use the matching complete PLINK suite for a genome-wide analysis. For R² input, replace `--plink-prefix` with `--r2-dir /data/r2/hg19` and keep each Parquet with its required matching sidecar. See [building R² panels](build-r2-panel.md).

The output root contains `metadata.json` and `ldscore.baseline.parquet`. Query runs additionally write `ldscore.query.parquet` and `ldscore.overlap.parquet`; `diagnostics/` contains the log and applicable audits. Default regression rows are packaged HapMap3 candidates after MHC-and-centromere exclusion, while reference contributors and annotation counts use their own retained reference universe.

Direct query scoring requires explicit `--baseline-annot-sources` plus one query route: `--query-annot-sources`, `--query-annot-bed-sources`, or `--query-annot-gene-list-sources`. See the [guided tutorial](../guided-tutorial.md#analysis-3-partition-heritability-with-functional-annotations-known-as-cell-type-specific-regression), [gene-list workflow](ldscore-from-gene-list.md), and [standalone annotation](../utility-functionalities/annotate.md). An exact gene index supplies its own configuration and baseline; omit live overrides in indexed mode.

## Effective CM Coordinates and Metadata Export

For PLINK input with `--ld-wind-cm`, LDSC interpolates CM from a matching explicit genetic map when supplied; otherwise BIM CM must be informative. `--genetic-map-hg19-sources` and `--genetic-map-hg38-sources` accept comma-separated exact paths, without `*` or `@` expansion. With R² input, CM comes from the authoritative metadata sidecar; map flags do not override it.

For PLINK input, `--export-ref-metadata` writes `ref_metadata/chrN_meta.tsv.gz` with the effective CM and reference metadata without changing the BIM. It is not an indexed-mode option. Sources: `build_parser` in [ldscore_calculator.py](../../../src/ldsc/ldscore_calculator.py) and `_resolve_genetic_map` in [_kernel/ref_panel.py](../../../src/ldsc/_kernel/ref_panel.py).

## Memory for many pathways

`ldscore --threads` controls chromosome worker processes. `build-gene-ldscore-index --threads` controls a chromosome thread pool. Both default to `1` (sequential); the name is shared for consistency.

This optimization targets workloads such as 1,000 pathways in one run, each later tested separately against shared baseline categories. Direct and indexed `ldscore --query-batch-size 1000` use a positive maximum query batch size; 1000 is the default. Smaller batches reduce query workspace. Direct `--threads 1` releases chromosome working data before advancing; larger counts use bounded worker processes. Freed memory can remain reserved by the allocator, so RSS may not drop immediately. Final HM3 baseline/query LD tables remain aggregate Parquet files and may stay materialized. Temporary annotation data stay below the output directory. See the [memory design](../../current/annotation-memory-design.md) and [batch regression guide](partitioned-h2.md#testing-enrichment-for-a-large-batch-of-pathways).
