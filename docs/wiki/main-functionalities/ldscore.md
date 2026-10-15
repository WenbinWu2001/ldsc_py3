# Calculate LD scores

Last updated on: 2026-10-15

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

Baseline inputs are prebuilt annotation tables. To prepare binary baseline categories from BED files, follow [building baseline annotations with `annotate`](../utility-functionalities/annotate.md#build-baseline-annotations-from-bed-files), including target SNP-grid creation, category-specific interval boundaries, and the all-ones `base` column.

## Effective CM Coordinates and Metadata Export

For PLINK input with `--ld-wind-cm`, LDSC interpolates CM from a matching explicit genetic map when supplied; otherwise BIM CM must be informative. `--genetic-map-hg19-sources` and `--genetic-map-hg38-sources` accept comma-separated exact paths, without `*` or `@` expansion. With R² input, CM comes from the authoritative metadata sidecar; map flags do not override it.

For PLINK input, `--export-ref-metadata` writes `ref_metadata/chrN_meta.tsv.gz` with the effective CM and reference metadata without changing the BIM. It is not an indexed-mode option. Sources: `build_parser` in [ldscore_calculator.py](../../../src/ldsc/ldscore_calculator.py) and `_resolve_genetic_map` in [_kernel/ref_panel.py](../../../src/ldsc/_kernel/ref_panel.py).

## Choosing workers for indexed scoring

For indexed gene-list scoring on Longleaf with a similar workload, start with **`--workers 8` and 8 allocated CPUs**, retaining **128 GiB RAM** and **one native numerical thread per worker**.
Use Slurm `--cpus-per-task=8 --mem=128G`; workers are managed processes, not a synonym for CPU threads.
This is a provisional starting point from the historical benchmark below, which predates the bounded indexed-batch implementation in `659738a`.
It does not establish the optimum worker count or RAM requirement for the current implementation.
The package default remains one worker; this recommendation does not cover direct LD scoring, index construction, or regression.

A matched Longleaf comparison at revision `e4cefba` used the same first 1,000 query files, all 22 chromosomes, 1,134,471 SNP rows, 53 baseline columns, saved index, and `--query-batch-size 1000`.
Every run requested 128 GiB and used one native numerical thread per worker.

| Workers / allocated CPUs | Total allocation elapsed | Chromosome interval | Reserved CPU-hours |
| --- | ---: | ---: | ---: |
| 22 | 4:47 | 6.00 s | 1.7539 |
| 8 | 3:00 | 25.17 s | 0.4000 |
| 4 | 5:28 | 40.03 s | 0.3644 |
| 1 | 8:05 | 140.13 s | 0.1347 |

Eight workers had the shortest observed total elapsed time and reserved 77.2% fewer CPU-hours than 22.
Four workers are a reasonable alternative when CPU reservation matters more than the observed turnaround advantage; one worker minimized reserved CPU-hours but finished slowest.
At the benchmark's verified Longleaf billing settings, 128 GiB contributed 32 billing units for every allocation, so fewer CPUs did not lower the hourly billing rate.

Twenty-two workers computed chromosomes fastest, while startup, input preparation, result assembly, and output writing dominated total time.
Increasing workers does not automatically accelerate those serial stages.
These were single runs on shared nodes in fixed order, with uncontrolled cache, filesystem, and load effects; the faster 8-worker total does not establish that worker contention slowed the 22-worker run.
All checked outputs matched exactly across the four settings, including SNP/query order, LD-score values, counts, overlaps, and query statuses.

The one-batch comparison does not establish a full-workload optimum or a lower RAM requirement.
With that earlier implementation, the separate full 22-worker run completed in 2:03:14, whereas the full 8-worker run timed out at 4:00:15 during serial assembly; its completed full-workload runtime remains unknown.
Those runs staged approximately 159 GiB of raw score payload before assembly.
The current implementation bounds raw fragments to one active query batch, as described below, so that historical staging volume is no longer its scratch model.
Keep the existing 128 GiB allocation until the current full workload is measured; neither these historical runs nor the bounded design alone establishes a smaller safe allocation.
Do not extrapolate one-batch time by multiplying by the number of batches: repeated operator reads, final output volume, and cache pressure affect the full run.

Source: the workspace [matched benchmark report](../../../../diagnostics/20261014_indexed_workers/fourway_report.md), with [measurements and commands](../../../../diagnostics/20261014_indexed_workers/fourway_details.md).
Benchmark jobs were 4933958 (22), 4933979 (8), 4930922 (4), and 4930923 (1); exact validator 4934033 passed with zero numeric differences.
These evidence links require the enclosing research workspace; the table and limitations above summarize the results for standalone package readers.

## Memory for many pathways

Live annotation inputs are prepared serially before chromosome computation. Whole-genome and chromosome-sharded inputs use bounded chunks sized for aligned source groups; duplicate cleanup remains global, including rsIDs repeated across chromosomes. At INFO, `diagnostics/ldscore.log` marks input reading, SNP identity checks and chromosome preparation, and completion with retained counts and elapsed time. See [preparation logging](../../current/workflow-logging.md#annotation-preparation). The preparation milestone does not report completed LD-score calculation.

`--workers` / `workers` defaults to one managed child process, including on HPC.
Every worker uses **1 BLAS thread**, in serial and parallel execution.
`joblib.cpu_count()` supplies the detected available CPU count; it does not measure idle CPUs or read the requested scheduler allocation.
Positive requests above that count fail before capping by available tasks; negative `-k` uses `max(1, detected_count + 1 - k)`.
LD scoring and index construction cap by chromosomes; regression caps by complete query-plus-baseline models and query batch size.
Zero, Boolean, and non-integer Python requests are invalid.
See [the execution policy](../../current/parallel-execution.md) for native-limit verification, memory costs, process startup, and failure handling.

This optimization targets workloads such as 1,000 pathways, each later tested separately against shared baseline categories.
`--query-batch-size 100` sets the maximum active query batch and defaults to 100.
Direct calculation prepares, computes, privately writes, and releases one batch before starting the next; smaller batches repeat reference/genotype work.
Direct workers now write their existing post-kernel float32 query scores to private chromosome fragments and return baseline data with scientific and diagnostic summaries.
Indexed calculation also admits only one query batch at a time; each chromosome task reloads its operator for that batch and writes a lossless float64 Parquet fragment.
Both modes use a persistent worker pool and determine global query filtering from chromosome summaries at their existing precision boundaries.
The shared writer consumes fragments in chromosome order, streams them into the staged public output, and removes that batch's fragments before submitting the next batch.
This avoids accumulating raw fragments for all queries or concatenating a full genome-wide query batch in memory.
Numerical accumulation, counts, and overlaps remain float64, and public LD-score values remain float32.
Streaming preserves both the broader reference contributor universe and the separately selected regression/output rows.

Raw-fragment scratch therefore scales with batch width, while staged final output still grows with the total number of queries.
Smaller batches can reduce scratch and active query arrays but add repeated direct reference work or indexed operator reads; reducing workers alone does not reduce the number of query batches.
Baseline values remain resident, and every active worker needs its reference/operator state and numerical workspace.
No new CLI option is required.

Python writing workflows return `LDScoreSource`. Use `source.read_queries(["pathway_A", "pathway_B"])` to load named columns, even across files; the call preserves requested column order, does not cache values, and imposes no read-width limit. Callers own the RAM for these explicit reads. For a small calculation that must write nothing, use prepared in-memory annotations with `LDScoreCalculator.run(..., output_config=None)` and keep all queries in one batch. See the [Python example](../../../tutorials/ld-score-calculation.md#small-python-calculations-without-writes).

Private scratch stays below the output directory and is cleaned on handled success or failure, after worker shutdown when a task fails.
Metadata is published after all batches succeed.
Resuming from a completed batch is a possible future feature and is not supported.
Freed memory can remain reserved by the allocator, so RSS need not drop immediately.
See the [memory design](../../current/annotation-memory-design.md), [direct integration validation](../../audits/2026-10-15_direct-streaming.md), and [batch regression guide](partitioned-h2.md#testing-enrichment-for-a-large-batch-of-pathways).
Sources: `LDScoreCalculator.run`, `LDScoreSource.read_queries`, [direct fragment staging](../../../src/ldsc/_direct_ldscore_batches.py), and `indexed_results` in [_indexed_ldscore_batches.py](../../../src/ldsc/_indexed_ldscore_batches.py).

## Direct assembly measurements

The October 15 Longleaf replays compared materialized direct assembly with chromosome-fragment assembly at frozen revision `03aa8ca`.
Each used precomputed scores for 22 chromosomes, 22 persistent workers, 53 baseline columns, 1,134,471 output SNP rows, and 1,000 queries per batch.
They measured assembly and writing, without real genotype/reference preparation or LD calculation.
Both variants used the same output rows; the improvement did not come from dropping SNPs.

| Replay workload | Median whole-job cgroup peak, materialized → streamed | Median wall seconds | Median combined CPU seconds |
| --- | ---: | ---: | ---: |
| One batch | 42.15 → 22.74 GiB | 155.32 → 102.25 | 210.35 → 168.81 |
| Three batches, original node | 43.84 → 25.99 GiB | 340.69 → 178.92 | 472.50 → 340.23 |
| Three batches, independent node | 43.62 → 25.71 GiB | 399.80 → 172.22 | 503.53 → 265.83 |

All 28 valid replays passed exact within-backend comparisons of saved values, schemas, SNP/query order, row groups, filtering statuses, counts, overlaps, and scientific metadata.
System-time outliers were retained and the two node groups were analyzed separately.
No comparison group showed an observed wall/CPU increase above the 5% screening threshold; these measurements do not guarantee a production speedup.

Active fragments totaled about 257 MiB across all 22 files, not per file.
Sampled workload-file peaks were approximately equal between variants: 0.644 GiB for one batch and 1.115 GiB for three batches, with 1.259 GiB allocated blocks for three batches.
These figures exclude fixed inputs and measurement logs and are not complete production scratch peaks; compression depends on the queries.

The integrated implementation has separate local tests using real small PLINK and parquet-R2 calculations, including metadata parity and failure cleanup.
Production memory savings and a smaller safe Slurm allocation remain unproven because the replay omitted reference preparation and kernel computation.
Source: the workspace [full replay report](../../../../diagnostics/20261015_direct_streaming_scale/report.md) and [lifecycle record](../../../../diagnostics/20261015_direct_streaming_scale/lifecycle.md), plus the [local integration audit](../../audits/2026-10-15_direct-streaming.md).
The replay links require the enclosing research workspace; the table and limitations above summarize the evidence for standalone package readers.
