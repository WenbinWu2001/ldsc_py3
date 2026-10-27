# Performance Measurements

Last updated on: 2026-10-27

This page records the measurements behind the package's resource guidance, each with its date, revision and limits.
The raw reports and probe scripts are kept outside the repository; the summaries below are the in-repository record.
None of these measurements is a performance guarantee for other data or systems.

## Parallel execution

The runtime contract these measurements support is in [parallel execution](../../reference/parallel-execution.md).

### Gene-index transport probe

A local probe on 2026-10-12, at the revision that moved gene-index construction from threads to spawned processes (`8f90165`), checked that migration.
It used two synthetic chromosomes, each with four PLINK SNPs plus 25,000 annotation rows and 64 continuous columns, and two workers.

| Execution | Wall time | Sampled peak summed RSS |
| --- | --- | --- |
| Threads | 2.01 s | 329.3 MiB |
| Spawned processes | 2.29 s | 638.2 MiB |

Each task sent about 5.1 KB of storage descriptors and catalog metadata, not baseline arrays.
Both runs produced the same index identity and gene-support counts, and twenty exact comparisons of baseline tables, sparse operators, counts, overlaps and atom statistics passed.
Each setting ran once, and summed RSS counts shared pages more than once, so the numbers show the direction of the memory cost rather than a stable overhead estimate.
Full-genome memory needs were not measured.

### Persistent worker cost

A local pilot on 2026-10-12 measured the managed child process that runs even with one worker.
Child startup took about 0.8 to 1.0 seconds, which was accepted as the cost of the design.
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

