# Parallel execution implementation validation

Last updated on: 2026-10-12

This validates the local `restructure` working tree against the [confirmed execution policy](../../current/parallel-execution.md).
Existing unrelated working-tree edits were retained.
The test environment was macOS arm64, Python 3.13.13, NumPy 2.4.5, SciPy 1.17.1, pandas 2.3.3, PyArrow 23.0.1, threadpoolctl 3.6.0, and joblib 1.6.0.
The loaded NumPy and SciPy numerical backend was Apple Accelerate.

## Results

| Check | Result |
| --- | --- |
| Full `python -m pytest -q --tb=short` | 1,914 passed; 133 subtests passed; one conditional skip |
| Separate `python -m unittest discover -s tests -p 'test*.py' -v` | 1,001 tests; passed with one skip |
| Cell-specific tutorial | All five notebook code cells executed in order in a guarded Python script, including the CLI subprocess example |
| Changed notebook/Python source | Parsed successfully |
| LD-score, gene-index, and batch-regression CLI help | `--workers` and the one-BLAS-thread policy present; retired `--threads` absent |
| Changed Markdown local links | All targets resolved |
| `git diff --check` | Passed |

The notebook check executed the code without a Jupyter kernel and wrote its generated results into temporary storage.
Tutorials requiring user-supplied scientific paths were checked for syntax and interface consistency, without executing those external workloads.
The two test runners ran sequentially because shared pybedtools cleanup can interfere across concurrent test processes.

## Behavior covered

- CPU detection delegates to `joblib.cpu_count()`; explicit positive requests are checked before work caps, and scheduler environment variables do not select workers.
- Serial and parallel workers use the same numerical target, with runtime inspection inside real children and unchanged parent numerical/Arrow settings.
- Direct LD reuses workers across query batches, and regression retains read-only numerical maps and complete-model tasks.
- Scientific reference and equivalence checks preserve LD scores, annotation counts/overlaps, and regression results.
- Direct and indexed chromosome logs persist on failure, with batch context and complete Python tracebacks.
- Forced child termination reports affected in-flight tasks and an unconfirmed crash source; workers are reaped before scratch cleanup.
- Pending direct futures are cancelled, later tasks are not submitted, and reference metadata remains private until successful computation.
- A real index-builder child failure removes an already-completed chromosome shard and leaves the build log available.
- Low-level prepared LD calculation without an output configuration remains file-free.
- Output collisions, authorized overwrites, and regression's strict/continue policies retain their existing semantics.

The principal policy tests are [test_execution_policy.py](../../../tests/test_execution_policy.py), [test_worker_policy.py](../../../tests/test_worker_policy.py), and [test_partitioned_query_workers.py](../../../tests/test_partitioned_query_workers.py).
Mock-based orchestration tests use explicit local calculation seams; resource, transport, and process-failure tests execute real spawned children.

## Empirical limits

The [focused index transport probe](index-transport.md) passed twenty exact scientific comparisons and measured the additional process RSS before migration.
Its small synthetic chromosomes do not establish genome-wide memory or useful LD chromosome scaling.
This implementation validation did not run on Longleaf, retest Linux numerical backends, or establish representative HPC Arrow I/O performance.
Those performance questions and repeated-work optimizations remain later benchmarking work.

## Suggested commit

```text
refactor(runtime)!: standardize compute workers

Use managed children and one BLAS thread for every calculation.
Share joblib CPU budgets and retain chromosome failure diagnostics.

BREAKING CHANGE: replace threads/--threads with workers/--workers.
```
