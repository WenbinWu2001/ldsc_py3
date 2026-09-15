# Parallel partitioned-h2 query workers

Last updated on: 2026-09-15

## Scope and implementation

Implemented locally on `restructure`, starting from clean revision `63ac296`. The existing `cf0a4ac` query-failure handling is retained. No HPC package copies, monitoring, or production jobs were changed. The user confirmed the worker, memory, and numerical-thread design before implementation.

`partitioned-h2 --threads N` and `RegressionRunner.estimate_partitioned_h2_batch(..., threads=N)` default to one inline worker. Parsing, validation, and resolution share the exact helpers used by `ldscore` and `build-gene-ldscore-index`: positive N requests N workers; negative offsets use process CPU affinity with a machine CPU-count fallback. All requests are capped by query count and loading-batch width. Positive N must fit the caller's allocation; `SLURM_CPUS_PER_TASK` is not read separately. Baseline-only runs remain inline. Each spawned worker calls the same `_fit_partitioned_query()` routine used by inline execution: complete model preparation, filtering, weighting, regression, jackknife, summaries, and private output staging. No chromosome-wise fits or new block boundaries are introduced.

[`_partitioned_h2_parallel.py`](../../src/ldsc/_partitioned_h2_parallel.py) owns read-only, dtype-preserving numeric maps and a bounded set of spawned processes. Shared arrays are written once; tasks carry one query-column descriptor and one unique ordinal staging path. Labels and count/overlap metadata are loaded once per worker. Nonnumeric/extension query columns use separate private snapshots so their original per-query preparation failures are preserved. A worker receives its next task only after returning the previous compact outcome. The parent restores input order before the existing stable sorting/publication path.

File-backed maps were chosen over named shared memory because they fit the existing output-owned scratch lifecycle, avoid a separate shared-memory quota/unlink owner, and permit normal filesystem cleanup after workers exit. The cost is staging I/O and disk space. A single workspace owns baseline maps and staged fits; each loading batch owns its query maps. Public `multiprocessing.Process`, pipe, and sentinel APIs provide explicit worker-death detection and termination on the supported Python 3.11–3.13 range. Workers are reaped before mapped files are removed on handled failure or interruption. SIGKILL of the parent cannot guarantee cleanup; there is no resume mechanism.

`threadpoolctl>=3.6,<4` controls supported BLAS/OpenMP libraries after imports. Apple Accelerate needs its separate macOS 15 `BLASSetThreading` API; older macOS requires `VECLIB_MAXIMUM_THREADS=1` at launch for parallel workers. Worker Arrow pools are also limited to one. Native limits are changed only in children; inline execution keeps caller settings. See the [complete runtime policy](../current/regression-configuration.md#43-query-workers-and-memory).

Strict runs still attempt all query models and publish no new scientific outputs if any fail. Continuation publishes successful models and retains all attempted statuses. Model exceptions retain stage, type, message, traceback, and available jackknife diagnostics. Output failures, process death, transport failures, and interrupts remain fatal. Structured `JackknifeIdentifiabilityError` pickling preserves its original type, failures, separators, and notes. The kernel's scientific calculations and canonical table schemas are unchanged.

## Profile before implementation

The original inline workflow was profiled with 1,000 queries using `cProfile`. It produced 999 successful models and one full-data collinearity failure. Total profiled process time was 125.46 seconds. A separate full test suite was running during this diagnostic profile, and profiling itself adds overhead; these values identify expensive phases rather than establish a clean speedup baseline.

| Original phase | Cumulative seconds |
| --- | ---: |
| Focal/full summaries | 48.86 |
| Estimator wrapper | 41.89 |
| Per-model preparation | 12.53 |
| Private fit staging | 10.56 |
| Final query-tree publication | 5.86 |
| Query-batch loading | 0.83 |
| Shared preparation | 0.01 |

Summaries, preparation, and private writing were substantial, so the worker boundary includes all of them. The final aggregate publication stays in the parent. No estimator-only speedup is presented as an end-to-end result.

## End-to-end benchmark

Workload: existing local chromosome-22 LD-score artifacts, **17,380 SNP rows, 53 baseline annotations, 1,000 queries**, loading width 100, and 200 requested jackknife blocks. Summary statistics are the existing synthetic fixture (seed 20260914, N=100000). These are computational measurements with no biological interpretation. Host: macOS 15.3.2 arm64, Python 3.13.13, NumPy 2.4.5 linked to Accelerate.

[`benchmarks/partitioned_h2_workers.py`](../../benchmarks/partitioned_h2_workers.py) launches a fresh real CLI process per configuration and measures through process exit, including imports, input loading/projection, shared preparation, maps, model work, summaries, private staging, final publication, and cleanup. It samples simultaneous parent-plus-descendant RSS using POSIX `ps` about every 0.1 seconds. Shared mapped pages may be counted repeatedly; sampling can miss brief peaks. The sampler is excluded. Numerical environment limits are set to one for every launch, including the inline reference. No other test suite ran during these measurements.

Two initial rounds ran in orders 1/2/4 and 4/2/1. After verifying Apple's explicit runtime limit, the final implementation was measured again in order 1/2/4:

| Workers | Final elapsed seconds | Speedup vs inline | Peak sampled process-tree MiB |
| ---: | ---: | ---: | ---: |
| 1 | 74.17 | 1.00× | 472 |
| 2 | 47.63 | 1.56× | 1,048 |
| 4 | 30.47 | 2.43× | 1,507 |

The earlier two rounds measured 74.23–76.35 s at one worker, 47.95–50.25 s at two, and 29.48–30.23 s at four. Peak RSS across all nine runs reached 472 MiB, 1,048 MiB, and 1,627 MiB respectively. CPU work per model remains similar; more workers increase simultaneous private model memory. Mapping adds roughly 0.1 seconds over the ten loading batches in the initial parallel rounds. Final publication remains a serial cost of roughly 3.3–3.6 seconds.

Every run produced the same 999 successes and one collinearity failure, with all private scratch removed. Comparisons cover every compact and full scientific table, including coefficients, standard errors and p-values; every coefficient-delete Parquet array; exact status and manifest ordering; and retained-column, SNP-count, block-count, filtering, and identity metadata. Numeric tolerance is `rtol=1e-10`, `atol=1e-12`. The original preimplementation inline output was also compared against the final inline output for all 999 successful models.

These measurements establish useful local concurrency on this bounded workload. They do not establish full-genome/18,000-query RAM requirements or Longleaf runtime. Whole-genome fits have larger private matrices, and network-filesystem map/output I/O, CPU allocation, cache state, and host load can change scaling. Do not extrapolate the measured memory as a production allocation.

Reproduce with an installed development environment and compatible local inputs:

```bash
python benchmarks/partitioned_h2_workers.py \
  --ldscore-dir /path/to/canonical/ldscores \
  --sumstats-file /path/to/compatible/synthetic.sumstats.gz \
  --query-batch-size 100 \
  --output-dir /path/to/new/benchmark-directory \
  --repeats 2
```

Raw measurements, profile totals, and comparison settings are in the [JSON evidence](2026-09-15_partitioned-query-workers.json). Input paths in that portable evidence use descriptive labels; host-specific locations remain in the local benchmark records. Detailed logs/results remain in the workspace's `runtime_estimates/2026-09-15_partitioned_workers/` directory, separate from source inputs and earlier session artifacts.

## Verification

Test-first checks observed the missing `threads` API, structured-exception unpickling failure, nonnumeric-query transport failure, and missing Accelerate runtime limit before their implementations passed. Focused coverage also exercises cross-file query reads, shuffled trait rows, both error policies, all-query failure, unchanged baseline behavior, strict overwrite preservation, actual worker crashes/exceptions, parent/worker interrupts, worker staging and parent publication failures, immutable dtype-preserving maps, effective-one pool avoidance, shared positive/negative worker resolution, and child/scratch cleanup.

Initial implementation verification in `ldsc3-dev` (Python 3.13.13), before the worker-policy alignment below:

- Full `python -m pytest -q`: **1,869 passed, 131 subtests passed, 1 skipped**, 311 warnings, 153.53 seconds. The original full suite passed 1,827 tests with the same skip and warning count before implementation.
- `python -m unittest discover -s tests -p 'test*.py' -v`: **1,000 tests, OK, 1 skipped**, 53.17 seconds. This check ran after pytest so their temporary files could not interfere.
- Expanded regression/kernel/parallelism selection: **109 passed**. Final native-runtime/failure selection: **32 passed**; the full suite also verifies workers override inherited seven-thread environment requests.
- `ldsc --help`, `python -m ldsc --help`, and `python -m ldsc partitioned-h2 --help` succeeded; the latter documents the default, negative counts, affinity discovery for negative counts, query/batch work caps, and native-thread policy.
- Python compilation, `git diff --check`, changed Markdown dates, and local links in the new audit/runtime-policy documents passed.

### Shared worker-policy alignment

The user subsequently requested exact parsing and resolution consistency with `ldscore` and `build-gene-ldscore-index`. Removed the query-specific CPU/SLURM resolver and reused `_resolve_worker_count(threads, min(query_count, query_batch_size))`. The new public-API test first failed because four requested workers and three queries resolved to one under restricted affinity/SLURM; it now resolves to three, matching the existing positive-count convention. Negative values retain affinity discovery and CPU offsets. Baseline-only and effective-one execution remain inline, and native worker limits are unchanged.

The shared workflow tests now exercise partitioned fitting alongside direct/indexed LD scoring and index construction under identical CPU affinity and a conflicting `SLURM_CPUS_PER_TASK`. CLI checks compare all three commands for defaults, positive/negative integers, and rejected values. Focused checks passed **81 tests** in 23.14 seconds. The full pytest suite then passed **1,866 tests and 131 subtests**, with one skip and 311 warnings, in 174.25 seconds. Redundant query-specific resolver cases were replaced by shared workflow and public-API checks. The subsequent unittest compatibility check passed **1,000 tests, one skipped**, in 50.76 seconds. Updated CLI help, Python compilation, Markdown dates/links, and `git diff --check` also passed. The earlier benchmark requests still resolve to 1/2/4 workers on the benchmark host, so these measurements remain applicable; no benchmark rerun was needed for this resolution change.
