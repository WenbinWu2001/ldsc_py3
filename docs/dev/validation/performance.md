# Performance Measurements

Last updated on: 2026-10-28

This page records the measurements behind the package's resource guidance, each with its date, revision and limits.
The raw reports and probe scripts are kept outside the repository; the summaries below are the in-repository record.
None of these measurements is a performance guarantee for other data or systems.

## Parallel execution

The runtime contract these measurements support is in [parallel execution](../../reference/parallel-execution.md).

### Gene LD-score index transport probe

A local probe on 2026-10-12, at the revision that moved gene LD-score index construction from threads to spawned processes (`8f90165`), checked that migration.
It used two synthetic chromosomes, each with four PLINK SNPs plus 25,000 annotation rows and 64 continuous columns, and two workers.

| Execution | Wall time | Sampled peak summed RSS |
| --- | --- | --- |
| Threads | 2.01 s | 329.3 MiB |
| Spawned processes | 2.29 s | 638.2 MiB |

Each task sent about 5.1 KB of storage descriptors and gene-coordinate-catalog metadata, not baseline arrays.
Both runs produced the same index identity and gene-support counts, and twenty exact comparisons of baseline tables, sparse operators, counts, overlaps and atom statistics passed.
Each setting ran once, and summed RSS counts shared pages more than once, so the numbers show the direction of the memory cost rather than a stable overhead estimate.
Full-genome memory needs were not measured.

### Persistent worker cost

A local pilot on 2026-10-12 measured the cost of the worker that runs even when `--workers` is one.
Worker startup took about 0.8 to 1.0 seconds, which was accepted as the cost of the design.
The runtime-only control added about 158 MB of sampled RSS and 108 MB of footprint accounting; these are two different memory measures, and further input duplication is possible.

### BLAS threads

A Longleaf pilot on 2026-10-12 found that single-threaded regression performs well and that higher BLAS thread limits gave no convincing benefit for indexed LD scoring on its compact workload.
The direct LD-score serial comparison did not identify a fastest setting; one BLAS thread per worker was chosen for simplicity and resource efficiency.
A workload-specific four-worker result from the same pilot does not change the default of one worker.

### Limits

These pilots do not establish LD-score scaling across chromosomes, genome-wide memory needs, or representative input and output performance on HPC storage.
The initial Arrow CPU and I/O thread defaults have not been validated on representative HPC storage.

## Partitioned heritability regression

Longleaf pilots measured two regression changes separately before they were integrated at revision `7a478e6`.
The values are phase medians, comparing the earlier implementation with the optimized one:

| Phase and workload | Wall time, s | CPU time, s | Measured resource change |
| --- | ---: | ---: | --- |
| Model fitting: 60 queries, 20 workers | 60.83 → 31.74 | 898.84 → 401.08 | Peak anonymous memory: 37.25 → 27.13 GiB |
| Publication: 100 completed fits | 5.456 → 1.917 | 3.355 → 0.650 | Logical writes: 139.856 → 1.667 MiB; anonymous memory roughly unchanged |

The model comparison used two earlier and two optimized runs on the same node and BLAS backend, with 1,127,525 aligned SNP rows and 53 baseline columns per fit, from direct LD-score inputs; indexed inputs use the same fitting code.
The publication comparison used four earlier and two optimized runs with the full 18,841-entry annotation dictionary.
Scientific results agreed exactly within each tested BLAS backend; differences across backends were also present between unchanged runs.
Both changes reduced measured phase wall time and CPU time, but these pilots do not establish a full production-job speedup or a smaller safe memory request, and fewer publication writes do not imply proportionally less anonymous memory or file cache.
The raw pilot report, with workload boundaries, repetitions, cache measurements and numerical checks, is kept outside the repository.

An earlier local benchmark on chromosome 22 with 1,000 queries, before the current copy and publication improvements, measured 74.17, 47.63 and 30.47 seconds with 1, 2 and 4 workers, including preparation and publication.
Its sampled peak process-tree RSS was 472, 1,048 and 1,507 MiB; these are not full-genome memory estimates.

## Indexed gene-list LD scoring

A matched Longleaf comparison at revision `e4cefba` used the first 1,000 query files, all 22 chromosomes, 1,134,471 SNP rows, 53 baseline columns, a saved gene LD-score index and `--query-batch-size 1000`.
Every run requested 128 GiB and used one native numerical thread per worker.
It predates the bounded indexed-batch implementation (`659738a`), so it does not establish the optimum worker count or memory request for the current implementation.

| Workers and allocated CPUs | Total elapsed | Chromosome interval | Reserved CPU-hours |
| --- | ---: | ---: | ---: |
| 22 | 4:47 | 6.00 s | 1.7539 |
| 8 | 3:00 | 25.17 s | 0.4000 |
| 4 | 5:28 | 40.03 s | 0.3644 |
| 1 | 8:05 | 140.13 s | 0.1347 |

Eight workers had the shortest total elapsed time and reserved 77.2% fewer CPU-hours than 22; one worker reserved the fewest CPU-hours but finished slowest.
At the benchmark's Longleaf billing settings, 128 GiB contributed 32 billing units for every allocation, so fewer CPUs did not lower the hourly rate.
Twenty-two workers computed chromosomes fastest, while start-up, input preparation, result assembly and output writing dominated total time; more workers do not speed up those serial stages.
These were single runs on shared nodes in fixed order, with uncontrolled cache, file-system and load effects.
All checked outputs matched exactly across the four settings, including SNP and query order, LD-score values, counts, overlaps and query statuses.

With that earlier implementation, a full 22-worker run finished in 2:03:14, and a full 8-worker run timed out at 4:00:15 during serial assembly, after staging about 159 GiB of raw scores.
The current implementation keeps raw fragments to one active query batch, so that staging volume no longer applies, but neither these runs nor the bounded design establishes a memory request below 128 GiB.
One-batch time should not be multiplied by the number of batches: repeated operator reads, output volume and cache pressure affect a full run.
The raw benchmark report, commands and validator output are kept outside the repository.

## Direct LD-score assembly

Longleaf replays on 2026-10-15 at revision `03aa8ca` compared materialized direct assembly with chromosome-fragment assembly.
Each used precomputed scores for 22 chromosomes, 22 persistent workers, 53 baseline columns, 1,134,471 output SNP rows and 1,000 queries per batch, and measured assembly and writing only, without reference preparation or LD-score computation.

| Replay workload | Median whole-job cgroup peak, materialized → streamed | Median wall seconds | Median combined CPU seconds |
| --- | ---: | ---: | ---: |
| One batch | 42.15 → 22.74 GiB | 155.32 → 102.25 | 210.35 → 168.81 |
| Three batches, original node | 43.84 → 25.99 GiB | 340.69 → 178.92 | 472.50 → 340.23 |
| Three batches, independent node | 43.62 → 25.71 GiB | 399.80 → 172.22 | 503.53 → 265.83 |

All 28 valid replays passed exact comparisons of saved values, schemas, SNP and query order, row groups, filtering statuses, counts, overlaps and scientific metadata within each backend, and both variants wrote the same output rows.
No comparison group showed a wall or CPU increase above the 5% screening threshold; this does not guarantee a production speedup.
Active fragments totaled about 257 MiB across all 22 files, and sampled workload-file peaks were about 0.644 GiB for one batch and 1.115 GiB for three batches (1.259 GiB allocated blocks), excluding fixed inputs and logs.
Because the replays omitted reference preparation and computation, they do not establish production memory savings or a smaller safe memory request.
The integrated implementation is also covered by local tests with small PLINK and R² matrix calculations, including metadata parity and failure cleanup ([`test_direct_batch_staging.py`](../../../tests/ldscore/test_direct_batch_staging.py)); the raw replay report is kept outside the repository.

