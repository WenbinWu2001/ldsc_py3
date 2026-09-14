# Thread-option consistency audit

Last updated on: 2026-09-14

Audited revision `6454dd8` on local `restructure`. The findings below describe that revision. The correction is now implemented in [`_parallelism.py`](../../src/ldsc/_parallelism.py): shared strict validation and affinity-aware resolution, followed by inline execution for an effective count of one. [`test_worker_policy.py`](../../tests/test_worker_policy.py) covers all three workflows with real computation and executor-construction checks. No HPC access was used.

## Command inventory

The constructed CLI parser exposes `--threads` on exactly two commands: `ldscore` and `build-gene-ldscore-index`. Neither `build-r2-panel` nor any regression command exposes this flag.

| Behavior | Direct `ldscore` | Indexed `ldscore` | `build-gene-ldscore-index` |
| --- | --- | --- | --- |
| Default | 1 | 1 | 1 |
| CLI value type | Integer | Integer | Integer |
| Zero | Rejected | Rejected | Rejected |
| Positive N | Requested workers, capped by chromosomes | Same | Same |
| Negative -k | Available CPUs + 1 - k, minimum 1, capped by chromosomes | Same shared resolver | Same arithmetic, but machine-wide CPU count |
| CPU discovery for negative values | CPU affinity, then `os.cpu_count()` fallback | Same shared resolver | `os.cpu_count()` |
| Parallel execution | Spawned processes | Spawned processes | Threads |
| Effective count of 1 | Inline | Inline | Inline only for requested 1 or one chromosome; negative values resolving to 1 still construct a pool |

Sources: [`LDScoreCalculator._run_batch`, `_run_chromosomes`, `_available_cpu_count`, and `_resolve_worker_count`](../../src/ldsc/ldscore_calculator.py); [`_write_chromosomes`](../../src/ldsc/_indexed_ldscore_batches.py); [`_run_gene_ldscore_index_build`](../../src/ldsc/gene_ldscore_index.py), around the `ThreadPoolExecutor` construction.

## Confirmed inconsistencies

1. **Automatic CPU budgeting differs.** Index building ignores CPU affinity for every negative value. It can exceed the process's available CPU set when `ldscore` would stay within it. This is the primary command-level defect.
2. **Effective sequential execution differs.** The index builder tests the requested value before resolving negative values. A negative value that resolves to one still constructs a one-worker thread pool when multiple chromosomes exist. This is unnecessary executor overhead, not a demonstrated numerical defect.
3. **Python API validation differs.** `LDScoreConfig` and `GeneLDScoreIndexBuildConfig` reject zero but accept `True`, `1.5`, and `'2'`. Public `run_indexed_ldscore` explicitly rejects booleans and non-integers before I/O. CLI parsing already enforces integer values, so this does not affect ordinary CLI syntax; it affects direct Python callers and when their errors appear.

The process-versus-thread distinction is explicitly documented in command help and is intentional. The flag represents concurrent chromosome workers, not an exact count of all operating-system threads.

## Shared limitations

- Explicit positive N can exceed CPU affinity in both commands; only the chromosome cap applies. This matches the current explicit-count contract.
- `ldscore` respects scheduler allocation only when that allocation is reflected in CPU affinity. Its resolver does not read `SLURM_CPUS_PER_TASK` or cgroup CPU quota files. A claim that it always respects every scheduler/container CPU limit would be too strong.
- Spawned LD-score workers call `_init_worker`, which sets default `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, and `MKL_NUM_THREADS` values to 1 while preserving user-specified values. Inline execution and the index builder do not call that initializer. These are environment defaults, not a verified runtime limit on already initialized native libraries. Native BLAS pool sizes were not measured because the local environment lacks `threadpoolctl`.

## Verification

Ran five real miniature index builds with two chromosomes, mocked machine CPU count 8, and mocked affinity containing one CPU. Intercepted executor construction while still using the real thread executor and build calculation; all five saved indexes retained expected gene support `[4, 4]`.

| Requested value | LD-score resolved workers | Index-builder execution |
| --- | --- | --- |
| 1 | 1 | Inline |
| 4 | 2 | Pool with 2 workers |
| -1 | 1 | Pool with 2 workers |
| -2 | 1 | Pool with 2 workers |
| -20 | 1 | Pool with 1 worker |

The fixture and argument setup are `write_inputs` and `index_args` in [`test_plink_workflow_resolution.py`](../../tests/test_plink_workflow_resolution.py). CPU conditions were supplied with `unittest.mock.patch`; the observed index pool sizes came from the actual builder, not a reimplementation of its arithmetic. Indexed LD-score routing to the shared resolver was verified in `_indexed_ldscore_batches._write_chromosomes`.

`python -m pytest tests/test_ldscore_parallelism.py -q`: **24 passed**. Existing tests cover defaults, zero rejection, negative/positive resolution, affinity, worker initialization, and numerical parallel-versus-sequential equivalence. Python config constructor probes confirmed the accepted invalid types described above.

## Recommended correction

Centralize nonzero-integer validation, affinity-aware CPU discovery, and worker-count resolution in a shared helper used by both commands and both LD-score modes. Resolve first, then choose inline execution when the effective count is one. Preserve the documented process/thread backends, positive-count behavior, chromosome cap, and negative-value arithmetic. Add cross-workflow tests under restricted affinity, fallback CPU discovery, and negative values that reach the minimum.

Treat native-library thread control as a separate explicit policy; do not claim the chromosome-worker option caps every native runtime thread. For scheduler allocations not reflected in CPU affinity, use an explicit allocation-derived count such as `--threads "${SLURM_CPUS_PER_TASK:-1}"` in SLURM scripts.
