# Calculate LD scores

Last updated on: 2026-09-14

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

The output root contains `metadata.json` and `ldscore.baseline.parquet`. One query batch writes `ldscore.query.parquet`; multiple batches write `ldscore.query.batch00001.parquet`, `ldscore.query.batch00002.parquet`, and so on. Query runs also write `ldscore.overlap.parquet`; `diagnostics/` contains the log and applicable audits. Every query file covers the same genome-wide regression rows, with chromosome row groups. The ordered `query_batches` manifest in root metadata records each file, its columns, and its row groups. Current readers require this manifest; regenerate older directories. Default regression rows are packaged HapMap3 candidates after MHC-and-centromere exclusion, while reference contributors and annotation counts use their own retained reference universe.

Direct query scoring requires explicit `--baseline-annot-sources` plus one query route: `--query-annot-sources`, `--query-annot-bed-sources`, or `--query-annot-gene-list-sources`. See the [guided tutorial](../guided-tutorial.md#analysis-3-partition-heritability-with-functional-annotations-known-as-cell-type-specific-regression), [gene-list workflow](ldscore-from-gene-list.md), and [standalone annotation](../utility-functionalities/annotate.md). An exact gene index supplies its own configuration and baseline; omit live overrides in indexed mode.

## Effective CM Coordinates and Metadata Export

For PLINK input with `--ld-wind-cm`, LDSC interpolates CM from a matching explicit genetic map when supplied; otherwise BIM CM must be informative. `--genetic-map-hg19-sources` and `--genetic-map-hg38-sources` accept comma-separated exact paths, without `*` or `@` expansion. With R² input, CM comes from the authoritative metadata sidecar; map flags do not override it.

For PLINK input, `--export-ref-metadata` writes `ref_metadata/chrN_meta.tsv.gz` with the effective CM and reference metadata without changing the BIM. It is not an indexed-mode option. Sources: `build_parser` in [ldscore_calculator.py](../../../src/ldsc/ldscore_calculator.py) and `_resolve_genetic_map` in [_kernel/ref_panel.py](../../../src/ldsc/_kernel/ref_panel.py).

## Memory for many pathways

Direct and indexed `ldscore --threads` control chromosome worker processes. The default is `1`; positive values request that many workers, `-1` uses available cores, and `-2` leaves one core free. The effective count is capped at the chromosome count, and `0` is invalid. Each concurrent chromosome increases memory use. `build-gene-ldscore-index --threads` controls its existing chromosome thread pool for index construction.

This optimization targets workloads such as 1,000 pathways, each later tested separately against shared baseline categories. `--query-batch-size 1000` sets the maximum active query batch and defaults to 1000. Direct calculation prepares, computes, privately writes, and releases one batch before starting the next; smaller batches repeat reference/genotype work. Indexed calculation keeps one chromosome operator per worker through that chromosome's batches, writing and releasing each batch before advancing. The coordinator assembles one genome-wide batch at a time. Baseline values remain resident, but completed query LD tables do not accumulate in the returned result.

Python writing workflows return `LDScoreSource`. Use `source.read_queries(["pathway_A", "pathway_B"])` to load named columns, even across files; the call preserves requested column order, does not cache values, and imposes no read-width limit. Callers own the RAM for these explicit reads. For a small calculation that must write nothing, use prepared in-memory annotations with `LDScoreCalculator.run(..., output_config=None)` and keep all queries in one batch. See the [Python example](../../../tutorials/ld-score-calculation.md#small-python-calculations-without-writes).

Private scratch stays below the output directory and is cleaned on handled success or failure. Metadata is published after all batches succeed. Resuming from a completed batch is a possible future feature and is not supported. Freed memory can remain reserved by the allocator, so RSS need not drop immediately. See the [memory design](../../current/annotation-memory-design.md), [numerical verification](../../audits/annotation-memory/sequential-query-batches.md), and [batch regression guide](partitioned-h2.md#testing-enrichment-for-a-large-batch-of-pathways). Sources: `LDScoreCalculator.run`, `LDScoreSource.read_queries`, and `indexed_results` in [_indexed_ldscore_batches.py](../../../src/ldsc/_indexed_ldscore_batches.py).
