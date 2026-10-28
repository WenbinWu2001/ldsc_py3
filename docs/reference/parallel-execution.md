# Parallel execution

Last updated on: 2026-10-28

Use `workers` in Python or `--workers` in the CLI to choose task concurrency.
The default is one worker on both laptops and HPC systems.
We use **1 BLAS thread for each worker**, in serial and parallel calculation.
The retired `threads` / `--threads` spelling is not accepted.

## Task and CPU budgets

| Workflow | One task | Concurrency limit |
| --- | --- | --- |
| Direct LD scoring | One chromosome in the active query batch | Requested chromosomes |
| Indexed LD scoring | One chromosome in the active query batch | Index chromosomes |
| Gene LD-score index construction | One chromosome | Requested build chromosomes |
| Batch regression | One complete query-plus-baseline model across all chromosomes | Query count and `query_batch_size` |

The CPU budget is the **detected available CPU count**, returned directly by [`joblib.cpu_count()`](https://joblib.readthedocs.io/en/stable/generated/joblib.cpu_count.html).
It counts logical CPUs, not physical cores.
Supported affinity and quota restrictions can affect this count.
It does not measure idle CPUs or read how many CPUs a user requested from Slurm or another scheduler.
An allocation not exposed through supported process restrictions cannot be inferred; choose `--workers` consistently with that allocation.
The package does not detect Slurm, reserve CPUs, or implement CPU binding.

| Request | Behavior before the task cap |
| --- | --- |
| Omitted or `1` | One managed worker process |
| Positive integer `N` | At most `N` workers; error if `N` exceeds the detected CPU budget |
| `-1` | The detected CPU budget |
| `-2` | One fewer than that budget, with a minimum of one |
| Negative integer `-k` | `max(1, detected_budget + 1 - k)` |
| `0`, Boolean, or non-integer Python value | Error |

Validate an explicit positive request against the CPU budget **before** capping by work.
With eight detected CPUs and two chromosomes, `--workers 4` resolves to two workers.
With two detected CPUs and one chromosome, `--workers 4` is an error.
The same resolver applies to regression, with its query-batch cap.
It reports the request, detected count and source, task unit/count, batch cap, and resolved worker count.
See [`_runtime/parallelism.py`](../../src/ldsc/_runtime/parallelism.py).

## Processes, numerical libraries, and memory

All task workflows use reusable workers, including one persistent worker for serial execution.
A ready worker receives the next task; chromosomes and models do not permanently own workers.
Direct and indexed LD keep their pools across query batches, and regression keeps its pool across query-loading batches.
There are no nested chromosome/model pools.
The parent coordinates inputs, collection, aggregation, and publication.

Direct scoring, indexed scoring, and index construction share the private `run_chromosome_tasks()` function in [`_runtime/parallelism.py`](../../src/ldsc/_runtime/parallelism.py).
It bounds outstanding futures by the resolved worker count, tracks their chromosomes, and collects each ready group before submitting replacements, using the same path for one or multiple workers.
On an observed exception it attempts cancellation and passes the original error, current chromosome, and outstanding chromosomes to the workflow's failure callback.
The workflow retains task preparation, result processing, progress and diagnostics, executor ownership, and all temporary-file and publication responsibilities.
Direct and indexed LD explicitly wait for their persistent executors to shut down on failure before the current query-batch workspace closes; an outer pool context alone would close too late.
Index construction exits its executor context before its enclosing workspace is cleaned up.
Regression retains its separate process communication and memory-mapping implementation.

BLAS and numerical OpenMP each target one thread throughout every worker's lifetime.
Arrow's CPU and I/O pools are controlled independently and each use one thread per worker.
Higher inherited numerical thread settings do not override this policy inside workers.
These settings do not modify the notebook or Python caller's global numerical or Arrow configuration.
Standalone numerical entry points without a concurrency option use a controlled serial worker and reuse the current worker when called from a worker.
This includes standalone regression fits and summaries, R² panel construction, quantile calculations, and legacy conversion's annotation overlap products.
See [`_runtime/worker_runtime.py`](../../src/ldsc/_runtime/worker_runtime.py) and [`_runtime/serial_calculation.py`](../../src/ldsc/_runtime/serial_calculation.py).

Runtime inspection checks supported loaded backends after importing NumPy and SciPy; these can use separate BLAS libraries.
Apple Accelerate uses its supported single-threaded mode, rather than an invented exact automatic thread count.
If enforcement or verification is unavailable, computation continues with an explicit unverified/unmet diagnostic and any observed settings.
Environment variables alone are not reported as proof of enforcement.
The implementation follows the runtime inspection and limiting APIs documented by [threadpoolctl](https://github.com/joblib/threadpoolctl#python-runtime-programmatic-introspection), and Arrow's separate [CPU](https://arrow.apache.org/docs/python/generated/pyarrow.set_cpu_count.html) and [I/O](https://arrow.apache.org/docs/python/generated/pyarrow.set_io_thread_count.html) setters.

A worker process, its main software thread, native helper threads, logical CPUs, and physical cores are different concepts.
Worker count is not the total operating-system thread count, a reservation of logical CPUs, or a guarantee of sustained CPU utilization.
Separate native and Arrow pools do not establish a whole-job CPU ceiling.

Regression shares suitable numerical inputs through read-only maps and stages each model separately.
LD and index-building tasks use chromosome-specific staged inputs and outputs; index construction returns compact component metadata rather than retaining completed operators in the parent.
Direct writing workers stage post-kernel float32 query scores and return baseline data, scientific statistics, variation summaries, and diagnostic descriptors.
Indexed scoring reloads operators for each batch, stages lossless float64 Parquet fragments, and returns only compact counts, overlaps, and variation summaries.
Both scoring workflows keep their worker pool alive across batches and final publication.
The shared writer consumes one chromosome fragment at a time in the producer's chromosome order and removes the active batch before submitting another.
Raw query scratch therefore scales with batch width rather than the total number of queries; privately staged final output files still grow with the run's total output.
The direct no-output API retains its materialized result path.
More workers require more private model/chromosome memory.
The package does not select workers from an automatic memory estimate.
The [focused index transport probe](../dev/validation/performance.md#gene-ld-score-index-transport-probe) supports the process migration but does not establish full-genome memory needs or HPC scaling.

Python scripts must put workflow calls under `if __name__ == "__main__":` for both serial and parallel execution because workers use `spawn`.
The CLI handles process startup; notebook examples also show a CLI subprocess pattern for environments that cannot launch importable Python workers.

## LD logs and failure cleanup

Both direct and indexed writing workflows use this layout, at every worker count:

```text
<output-dir>/diagnostics/
├── ldscore.log
└── chromosome-logs/
    ├── chr1.log
    ├── chr2.log
    └── ...
```

The main log records overall progress and each chromosome's diagnostic path.
Chromosome logs append across batches in the same run and record detailed execution and full Python exception tracebacks.
A subsequent run follows `--overwrite`: refuse existing owned artifacts by default, or reset owned logs for an explicitly authorized replacement.
An abnormal worker termination may produce no traceback; the main log lists affected in-flight tasks without guessing which chromosome caused the crash.
Low-level `LDScoreCalculator.run(..., output_config=None)` remains a no-file-output API for already-prepared single-batch inputs.

Any required chromosome failure fails the run.
Stop submitting tasks, attempt to cancel pending tasks, let running workers finish through normal executor shutdown, and then remove temporary outputs, including completed chromosome shards.
Cancellation is best effort: `Future.cancel()` cannot stop work already marked running, which can include dispatched work, following [Python's executor contract](https://docs.python.org/3.13/library/concurrent.futures.html#concurrent.futures.Future.cancel).
Logs and existing error diagnostics remain available.
There is no partial-success, checkpoint, resume, or automatic retry workflow for chromosomes.
This failure handling does not change intentional skip and eligibility rules, such as skipped query annotations.
All required chromosomes must succeed before cross-chromosome aggregation, and all batches must succeed before final scientific publication.

Optional direct `--export-ref-metadata` files are staged with temporary results and published only after successful computation.
Handled computation failures remove these files with the run's scratch.
Final file publication retains its existing behavior, without additional transactions or rollback; the gene LD-score index builder retains its existing atomic publication transaction.
Interruption handling is unchanged, and forced termination can leave temporary files.
Regression's separate explicit `continue_on_query_error` policy is unchanged.
See [`_ldscore/diagnostics.py`](../../src/ldsc/_ldscore/diagnostics.py), [`LDScoreCalculator.run`](../../src/ldsc/_ldscore/calculator.py), and [`indexed_results`](../../src/ldsc/_ldscore/gene_index/batches.py).

Representative chromosome scaling, full-genome memory, HPC Arrow I/O, and repeated serial work remain benchmarking topics.
The uniform BLAS policy is a confirmed resource policy, not a claim of universally optimal performance.
Tested behavior is in [`test_worker_policy.py`](../../tests/runtime/test_worker_policy.py), [`test_execution_policy.py`](../../tests/runtime/test_execution_policy.py), and [`test_ldscore_parallelism.py`](../../tests/ldscore/test_ldscore_parallelism.py).
