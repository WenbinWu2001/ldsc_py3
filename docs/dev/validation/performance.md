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
