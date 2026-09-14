# Annotation preparation parallelism evaluation

Last updated on: 2026-09-14

The four serial preparation optimizations are committed as `ce2bce9` (`perf(annotation): streamline preparation`). This follow-up evaluates concurrency before integrating it into the package. Production preparation remains serial; the only executable additions are the [benchmark adapter](../../../benchmarks/annotation_preparation_parallel.py) and its [contract checks](../../../benchmarks/annotation_preparation_parallel_check.py). All work ran locally, without HPC access or job changes.

## Decision

**Process workers are worthwhile for the measured large, narrow chromosome shards.** Four processes reduced median preparation time from **22.77 to 17.01 seconds (25.3%)**, while the post-identity phase fell from **9.15 to 3.57 seconds (61.0%)**. Median sampled process-tree peak RSS increased from **399 to 1,075 MiB**. This supports a focused production follow-up for selecting rows and writing shards after the shared global index is committed.

Four processes offer a useful local balance. Eight reduced total time to 16.18 seconds, another 0.83 seconds, but increased RSS to 1,748 MiB. Two processes used 682 MiB and still saved 20.1% of total time. These are measured trade-offs, not a proposed hard-coded worker limit.

**Skip numeric-only parallelism and thread-based selection.** Numeric writing represented only 5.3% of median serial time on the real inputs and 10.2% on the wide inputs. It did not dominate either workload. Real-data selection alone accounted for 6.63 seconds, so useful concurrency must include identity selection. Thread workers increased the real post-identity phase from 9.15 seconds to 10.03, 16.34, and 53.76 seconds at 2, 4, and 8 workers.

**The wide synthetic workload does not justify concurrency.** Process-worker medians were 2.4–6.5% slower overall, with higher memory and scratch use. Input scanning dominated. Do not treat annotation width, chromosome count, or the availability of more CPUs as evidence that preparation parallelism will help.

## Measured scope and design

The serial comparator is the optimized [`prepare_annotation_sources()`](../../../src/ldsc/_annotation_sources.py), not the earlier pickle-heavy implementation. The initial characterization run took 25.03 seconds: 9.85 in scanning, 5.26 in `DiskIdentityIndex.add()`, 8.07 in selection plus metadata writing, and 1.45 in numeric writing. Of selection time, 6.66 seconds were inside `DiskIdentityIndex.select()`. This justified testing the combined post-identity phase instead of starting with numeric-only workers.

The benchmark temporarily substitutes the coordinator calls to [`_select_rows()` and `_write_values()`](../../../src/ldsc/_annotation_sources.py). Workers reuse those real serial functions and the existing [`DiskIdentityIndex.select()`](../../../src/ldsc/_annotation_identity.py) policy. It does not duplicate the annotation pipeline or create independent chromosome identity indexes.

The tested design has these boundaries:

- Scan, validate, discover layouts, align logical rows, and build the global identity index serially. Explicitly commit its last transaction before launching readers.
- Submit one coarse task per discovered chromosome group. Each task receives source descriptors and opens its own read-only SQLite connection, with the existing 8 MiB page-cache limit. It selects rows and writes metadata, numeric values, and private diagnostics under its own writable directory.
- Bound pending tasks by the effective worker count. Reuse the existing worker-count convention from [`_resolve_worker_count()`](../../../src/ldsc/ldscore_calculator.py), and use spawned processes as in existing LD-score computation. Benchmark Arrow CPU and I/O pools are each limited to one thread.
- Replay private diagnostic spools in group order through `IdentityDropSpool`, preserving global reason order and row order. Completion order does not control scientific output order.
- Keep the shared database until all workers finish. Each worker releases only its own source staging after reading it. On an observed failure, stop submission, cancel pending work, and wait for active workers before the parent workspace closes.
- Use the original serial calls for one worker or a single group, including whole-genome column sources. No repartitioning, nested pool, new public configuration, or publication path is added.

The gene-index builder currently calls preparation before creating its chromosome computation thread pool in [`gene_ldscore_index.py`](../../../src/ldsc/gene_ldscore_index.py). The benchmark does **not** connect preparation to that command's `--threads` setting. Any production follow-up should preserve this phase separation and default serial path, reuse its worker budget, and validate through the public builder before release.

## Inputs and measurement

Hardware: Apple M1 Pro, eight logical CPUs, 16 GiB RAM, macOS 15.3.2. Software: Python 3.13.13, NumPy 2.4.5, pandas 2.3.3, PyArrow 23.0.1. Exact versions, input paths, row counts, raw phase timings, per-group worker durations, and all 31 runs are saved in [the machine-readable results](preparation-parallelism.json).

| Workload | Logical rows | Chromosomes | Baseline/query columns | Compressed input | Identity drops |
| --- | ---: | ---: | ---: | ---: | ---: |
| Real `baseline_v1.2` | 1,892,184 | Complete shards 15–22 | 53 / 0 | 30.79 MiB | 0 |
| Synthetic wide shards | 48,000 | 1–8; 6,000 rows each | 8 / 1,000 | 4.23 MiB | 8 |

The real files are the unmodified local `data/ld_refs/baseline_v1.2/baseline.<chrom>.annot.gz` sources. The wide data reuse the [existing deterministic generator](../../../benchmarks/annotation_preparation.py): annotation value at zero-based row `i`, within-source column `j` is `((7*i+j)%41-20)/8`. The first rsID is shared by all chromosomes. Aligned baseline/query records describe one logical row, so eight records are dropped globally, not sixteen. Input generation is outside measurement.

Each run uses a fresh Python process and fresh owned scratch; runs do not overlap. The timer encloses preparation, including worker startup, shutdown, and diagnostic collection. Reference snapshots and comparisons occur afterward. Process-tree RSS and private scratch are sampled every 50 ms. RSS is a simultaneous sum over the benchmark process and descendants, including pool helpers and the sampling subprocess; shared pages can be counted more than once. Scratch bytes are logical file sizes, including prepared shards, not physical disk blocks or transferred bytes. Snapshots and original inputs are outside scratch measurements. Samples can miss short peaks.

Tables report medians of three fresh runs for each serial/process configuration. The initial real characterization run is reported separately; matched real serial medians use runs 2–4. Thread configurations were screened once each and were not repeated after showing no worthwhile gain. Filesystem caches were not flushed; these are local cached-input measurements, not cold-cache tests. Ordinary desktop timing variation remains, including a 9.07-second wide serial run whose scan alone took 7.29 seconds. All observations are retained in the JSON.

## Results

`Post-identity` means selection, retained metadata writing, numeric writing, and, for pool variants, startup/shutdown and diagnostic merging. Worker durations overlap and must not be summed as wall time. `Add` times the actual identity `add()` calls; metadata reblocking and other coordinator work outside the named calls remain in total runtime. Phase medians need not sum to the median total.

### Real chromosome-sharded annotations

| Execution | Runs | Total seconds (range) | Scan / add seconds | Post-identity seconds | Peak RSS MiB | Peak scratch MiB | Peak scratch files |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Serial | 3 | 22.77 (21.96–24.42) | 8.57 / 4.74 | 9.15 | 399 | 553 | 33 |
| 2 processes | 3 | 18.19 (17.79–18.25) | 8.47 / 4.65 | 4.66 | 682 | 560 | 22 |
| 4 processes | 3 | 17.01 (16.00–18.65) | 8.34 / 4.65 | 3.57 | 1,075 | 681 | 28 |
| 8 processes | 3 | 16.18 (15.81–16.75) | 8.17 / 4.71 | 2.99 | 1,748 | 669 | 35 |
| 2 threads | 1 | 24.45 | 9.09 / 4.97 | 10.03 | 400 | 618 | 23 |
| 4 threads | 1 | 29.06 | 7.89 / 4.52 | 16.34 | 494 | 614 | 26 |
| 8 threads | 1 | 67.18 | 8.45 / 4.62 | 53.76 | 501 | 601 | 34 |

Four processes use 2.69 times the sampled serial RSS and 23.2% more peak logical scratch. Fewer concurrent files at some counts result from selecting, writing, and releasing each group together; that does not imply fewer bytes or I/O operations. The shared database remains live during parallel numeric writing. File and byte peaks need not occur together or increase monotonically with worker count.

### Wide annotations

| Execution | Runs | Total seconds (range) | Scan / add seconds | Post-identity seconds | Peak RSS MiB | Peak scratch MiB | Peak scratch files |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Serial | 3 | 6.15 (6.11–9.07) | 4.90 / 0.07 | 0.99 | 415 | 212 | 50 |
| 2 processes | 3 | 6.55 (6.37–6.65) | 5.04 / 0.06 | 1.07 | 694 | 235 | 41 |
| 4 processes | 3 | 6.30 (6.15–6.64) | 4.89 / 0.08 | 0.94 | 881 | 281 | 49 |
| 8 processes | 3 | 6.32 (6.08–6.55) | 4.85 / 0.08 | 1.04 | 1,287 | 328 | 63 |
| 2 threads | 1 | 6.22 | 4.94 / 0.06 | 0.81 | 410 | 235 | 41 |
| 4 threads | 1 | 6.45 | 4.89 / 0.06 | 1.12 | 438 | 281 | 49 |
| 8 threads | 1 | 6.96 | 4.73 / 0.06 | 1.79 | 402 | 374 | 65 |

The small total differences do not establish a reliable wide-input speedup. In particular, four processes more than doubled sampled RSS without improving median total runtime. This synthetic case is smaller in row count and highly compressible; it isolates a wide-column regime rather than representing every large pathway dataset.

## Correctness and cleanup

All 29 comparison runs matched their serial reference for every numeric value, float32 dtype, metadata field/dtype, row and column order, chromosome scope/order, and complete decompressed diagnostic output. The other two runs created those serial references. All 31 runs removed their owned scratch.

The benchmark adapter passed 14 independent known-value and failure checks across threads/processes and 1/2/4/8 requested workers. They cover global cross-chromosome duplicates, allele conflicts and invalid alleles, an allele-free shard, aligned baseline/query observations, unsorted retained rows, chromosome restriction after global selection, a chromosome with no retained rows, deterministic diagnostics, chromosome-count worker caps, one-worker execution, whole-genome serial fallback, and SQLite write rejection on read-only connections. Injected failures propagated and all started workers finished while the shared database still existed; only then did workspace cleanup finish.

The unchanged production preparation/storage suites passed **33 tests in 7.30 seconds**. The full production suite and unittest compatibility results for `ce2bce9` are in the [serial optimization audit](annotation-preparation.md); those broad suites were not rerun for this benchmark-only follow-up. The new scripts also compile, and their command help succeeds. No production concurrency or gene-index publication change is claimed as validated by this experiment.

## Reproduction

Run from the repository root in the editable `ldsc3-dev` environment. Process-tree RSS sampling requires permission to run `ps`. Use a fresh output root; scripts refuse existing per-run destinations.

```bash
annotation_python=/Users/wenbinwu/miniforge3/envs/ldsc3-dev/bin/python
annotation_data=/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/data/ld_refs/baseline_v1.2
annotation_runs=$(mktemp -d /private/tmp/ldsc-preparation-parallel-XXXXXXXX)

"$annotation_python" benchmarks/annotation_preparation_parallel.py \
  --baseline "$annotation_data"/baseline.{15..22}.annot.gz \
  --output "$annotation_runs/serial" --save-reference

for annotation_workers in 2 4 8; do
  "$annotation_python" benchmarks/annotation_preparation_parallel.py \
    --baseline "$annotation_data"/baseline.{15..22}.annot.gz \
    --executor process --workers "$annotation_workers" \
    --output "$annotation_runs/process-$annotation_workers" \
    --reference "$annotation_runs/serial/reference"
done

"$annotation_python" benchmarks/annotation_preparation_parallel_check.py "$annotation_runs/checks"
"$annotation_python" -m pytest -q tests/test_annotation_preparation.py tests/test_annotation_storage.py --tb=short
```

Use `--executor thread` for thread screens and unique output names for repeats. Generate the wide inputs in a separate process with `benchmarks` on `sys.path`, set `annotation_preparation.CASES = {"wide-eight": (8, 6000, 8, 1000, True)}`, and call its `generate(Path(fresh_directory))`. Pass the ordered baseline/query paths from the resulting `wide-eight/inputs.json` through `--baseline` and `--query`. Save a serial reference for that dataset before testing worker counts. The exact locally used paths and each measured configuration are retained in the JSON.

## Remaining limits

Scanning and global identity construction remain serial and account for most runtime after the successful process optimization. Numeric-only workers would leave the large real-data identity-selection cost untouched. Additional input-scan parallelism is outside this evaluation.

These preparation timings do not measure an entire gene-index build or LD-score run, and a 25.3% preparation saving is not a 25.3% whole-workflow saving. Only eight real chromosome shards and one synthetic width/row regime were measured. No HPC runtime, filesystem throughput, or memory result is inferred from them. SQLite statement execution does not imply one physical disk read per statement; local OS and SQLite caches matter. More readers and writers can contend on a shared HPC filesystem, so additional workers must not be assumed to improve its throughput.
