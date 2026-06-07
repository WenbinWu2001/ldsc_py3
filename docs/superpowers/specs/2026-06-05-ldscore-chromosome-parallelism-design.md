# ldscore cross-chromosome parallelism — design

Date: 2026-06-05
Status: approved (design); implementation plan to follow.

## 1. Problem and goal

LD-score calculation is the most time-consuming step in the pipeline. The
chromosome loop in `LDScoreCalculator.run` is sequential
(`src/ldsc/ldscore_calculator.py:348`): each chromosome's LD scores are computed
one after another, leaving all but one CPU core idle on a multi-core machine.

Per-chromosome compute is already fully independent — `compute_chromosome`
(`ldscore_calculator.py:388`) takes a per-chromosome annotation slice plus
configs and returns a self-contained `ChromLDScoreResult`, using no mutable
instance state. The chromosomes are therefore an embarrassingly parallel
workload.

**Goal.** Run the per-chromosome computations concurrently on one machine via a
process pool, then aggregate exactly as today. On an 8-core machine this is
expected to cut the LD-score wall-clock by ~4–6× (capped near ~10–12× by the
largest chromosome, chr1 ≈ 8% of the genome).

### Non-goals

- **No change to the public output contract.** The output stays a single
  canonical result directory (`metadata.json`, `ldscore.baseline.parquet`,
  optional `ldscore.query.parquet`) written by `LDScoreDirectoryWriter`. This is
  Approach A from brainstorming; per-chromosome public output files (Approach C)
  are explicitly rejected.
- **No change to LD-score math or numerical output.** Results must be identical
  to the sequential path (same row order, same values) for any `threads`.
- **No within-chromosome parallelism**, no thread pool, and no change to the
  per-window kernel (strategies #1/#2/#5 from the speedup analysis are separate
  later efforts).
- **No cluster / job-array support** (Approach B). Single-machine only.

## 2. Background: current run flow and the aggregation seam

`LDScoreCalculator.run` (`ldscore_calculator.py:294`) does:

1. Validate config, derive `chromosomes = _chromosomes_from_bundle(bundle)`.
2. **Sequential loop** (`:348`): for each chromosome, slice the bundle
   (`_slice_annotation_bundle`), call `compute_chromosome`, append the
   `ChromLDScoreResult` to an ordered list. A per-chromosome
   `LDSCInputError`/`ValueError` is routed through
   `_warn_and_skip_empty_intersection`: empty-intersection chromosomes are
   warned-and-skipped, all other errors re-raised.
3. If no results survive, raise (`:364`).
4. `_aggregate_chromosome_results(...)` (`:506`) concatenates per-chromosome
   tables into the aggregate `LDScoreResult`.
5. If `output_config` is set, `self.output_writer.write(...)` writes the
   canonical directory; `_write_chromosome_aligned_parquet` (`outputs.py:62`)
   already lays the single parquet out as **one row group per chromosome**.

Only step 2 changes. Steps 3–5 are untouched; the single-file,
chromosome-row-grouped output already preserves per-chromosome locality, so there
is nothing to gain from per-chromosome files on one machine.

### What a worker needs (cross-process serialization)

The macOS default multiprocessing start method is **spawn**, so workers do not
inherit parent memory: the worker callable must be a module-level function and
every argument must be serializable by `multiprocessing`. Confirmed
serializable inputs:

| Input | Type | Notes |
|---|---|---|
| `chrom` | `str` | — |
| per-chrom annotation slice | `AnnotationBundle` of DataFrames | from `_slice_annotation_bundle` |
| `RefPanelConfig` (`ref_panel.spec`) | frozen dataclass | worker rebuilds its own `RefPanel` via `RefPanelLoader`; **the live `RefPanel` object with open pyarrow handles is not sent** |
| `LDScoreConfig`, `GlobalConfig` | frozen dataclasses | — |
| `regression_snps` | `set[str]` / `RestrictionIdentityKeys` / `None` | sent once via pool initializer, not per task |
| `ChromLDScoreResult` (return) | DataFrames | already aggregated today; same peak memory |

Workers rebuild the reference panel locally (pass `RefPanelConfig`, not the live
adapter). Each worker handles a disjoint subset of chromosomes, so rebuilding the
per-chromosome metadata cache in each worker is correct and not wasteful.

## 3. Design

### 3.1 New knob: `LDScoreConfig.threads`

Add one integer field to `LDScoreConfig` (`config.py:431`), following the
scikit-learn/joblib `n_jobs` convention — the industrial standard for a compute
library (one integer, default sequential, opt-in parallelism):

```
threads: int = 1
```

Semantics (resolved by `_resolve_worker_count(threads, n_chromosomes)`):

- `1` (default) → **sequential in-process path**, no pool created. The default
  is sequential/opt-in, matching scikit-learn, BWA, STAR, bowtie2, samtools.
- `N > 1` → pool of `N` workers, capped at the chromosome count.
- `-1` → all available cores; `-2` → all but one; `-k` → `n_cpus + 1 - k`
  (joblib semantics), capped at the chromosome count.
- `0` → rejected by `__post_init__` (ambiguous).

**Core counts respect CPU affinity**, not the raw machine size:
`_available_cpu_count()` prefers `os.sched_getaffinity(0)` (Linux; honors
SLURM/cgroup/cpuset/Docker allocations) and falls back to `os.cpu_count()` on
platforms without affinity support. This prevents oversubscribing shared HPC
nodes — the key correctness rule, more important than the default value.

The CLI exposes a single integer `--threads` (default `1`), wired through
`run_ldscore_from_args` alongside the existing `--snp-batch-size`. No boolean
on/off switch is added: `1` already means "no parallelism", so a separate flag
would be redundant.

### 3.2 Module-level worker

A new module-level function (spawn-safe, serializable inputs) performs one
chromosome's work end-to-end:

```
def _compute_one_chromosome(chrom, chrom_bundle, ref_panel_spec,
                            ldscore_config, global_config) -> _ChromOutcome
```

It rebuilds `ref_panel = RefPanelLoader(global_config).load(ref_panel_spec)`,
runs the same body as today's `compute_chromosome`, and returns a tagged outcome
rather than raising across the process boundary:

- success → `_ChromOutcome(chrom, result=ChromLDScoreResult, skipped=False)`
- empty-intersection (the error class `_warn_and_skip_empty_intersection`
  currently swallows) → `_ChromOutcome(chrom, result=None, skipped=True,
  skip_message=...)`
- any other exception → re-raised in the worker; the pool surfaces it to the
  parent, which fails the whole run (matching today's "re-raise" branch).

`regression_snps` is injected into each worker once through a pool
**initializer** that stores it in a module-level global, avoiding repeated
serialization of a potentially large SNP set per task.

### 3.3 Parent orchestration (replaces the loop body)

`run` builds the ordered chromosome list, resolves `worker_count` from
`threads` and the chromosome count, then branches on the **resolved**
`worker_count` (not the raw `threads`, so `threads=-1` resolving to `1` on a
single-core or single-chromosome run still takes the inline path):

- **`worker_count == 1`** → call `_compute_one_chromosome` inline in a `for` loop
  (no pool). Preserves the current path exactly for the default.
- **`worker_count > 1`** → create a `concurrent.futures.ProcessPoolExecutor`
  (spawn context) with the resolved worker count and the regression-keys
  initializer. Submit one task per chromosome. **Collect results into a dict
  keyed by chromosome, then re-order by the original chromosome list** before
  aggregation, so output order is deterministic regardless of completion order.

Skip handling mirrors today: a `skipped` outcome logs the same warning and is
dropped; the "no surviving results" error (`:364`) is unchanged. Aggregation
(`_aggregate_chromosome_results`) and writing run exactly as before on the
ordered surviving results.

### 3.4 Logging and start method

- Use an explicit **spawn** context (`multiprocessing.get_context("spawn")`) so
  behavior is identical across macOS and Linux and never relies on fork
  inheritance.
- The pool initializer reconfigures logging in each worker (level from
  `GlobalConfig`) so worker progress lines (`"Computing chromosome ..."`,
  `"Finished chromosome ..."`) still appear. Interleaving across workers is
  expected and acceptable; each line is chromosome-tagged.
- BLAS thread oversubscription guard: with `W` worker processes each running
  multi-threaded BLAS `np.dot`, the machine can oversubscribe `W × BLAS_threads`
  cores. The initializer sets `OMP_NUM_THREADS`/`OPENBLAS_NUM_THREADS` to a small
  value (default 1) unless the user has set them, so process-level parallelism is
  not fighting thread-level BLAS parallelism. This guard is documented and
  overridable.

## 4. Data flow

```
run()
 ├─ chromosomes = ordered list
 ├─ resolve worker_count(threads, n_chromosomes)   # affinity-aware
 ├─ if worker_count == 1:
 │     for chrom: outcome = _compute_one_chromosome(...)   # inline, unchanged path
 └─ else:
       with ProcessPoolExecutor(spawn, worker_count, initializer=set regression_snps + logging + BLAS threads):
          futures = { submit(_compute_one_chromosome, chrom, slice, spec, cfgs): chrom }
          outcomes = gather()                              # any non-skip exception → propagate
 ├─ surviving = [outcomes[c].result for c in chromosomes if not skipped]  # deterministic order
 ├─ if not surviving: raise (unchanged)
 ├─ result = _aggregate_chromosome_results(surviving, ...)  # unchanged
 └─ output_writer.write(result, output_config)              # unchanged: one row group / chrom
```

## 5. Error handling

| Situation | Today (sequential) | With pool |
|---|---|---|
| Empty intersection on a chromosome | warn + skip | worker returns `skipped` outcome → parent warns + skips (same message) |
| Other `LDSCInputError`/`ValueError` | re-raise, abort run | worker re-raises → propagates to parent → abort run |
| Unexpected worker crash (segfault/OOM) | n/a | `ProcessPoolExecutor` raises `BrokenProcessPool`; parent wraps it in an `LDSCInternalError` advising a lower `--threads` |
| `threads == 0` | n/a | `LDScoreConfig.__post_init__` raises `LDSCConfigError` |

## 6. Testing

Correctness is defined as **identical output to the sequential path**.

1. **Equivalence test (core):** on a small fixture panel + annotation, compute
   with `threads=1` and `threads=3` and assert the written
   `ldscore.baseline.parquet`, `ldscore.query.parquet`, and the LD-score columns
   of `metadata.json` are identical (same row order, same float values). This is
   the gate. A second equivalence test uses `threads=-1` (all cores) to exercise
   the affinity-resolved pool path.
2. **Determinism:** `threads=3` produces byte-identical baseline parquet across
   two runs (chromosome ordering is stable).
3. **Skip path:** a fixture where one chromosome has empty annotation
   intersection still completes and skips exactly that chromosome under
   `threads>1`, matching sequential.
4. **Config validation:** `LDScoreConfig(threads=0)` raises `LDSCConfigError`;
   negative and positive values are accepted.
5. **Resolution:** `_resolve_worker_count` honors the joblib convention
   (`-1`→all, `-2`→all-but-one via a mocked `_available_cpu_count`), caps at the
   chromosome count, and `_available_cpu_count` prefers `os.sched_getaffinity`.
6. **Dispatch:** `threads=1` constructs no `ProcessPoolExecutor`; `threads>1`
   constructs exactly one (spy assertions).
7. **CLI plumbing:** `--threads 2` / `--threads -1` reach `LDScoreConfig.threads`;
   bare default is `1`.

Because the default is sequential, in-process orchestration tests that **mock**
`compute_chromosome`/`compute_chrom_from_parquet` or stub `ref_panel` work
unchanged (`threads=1`). Skip-under-pool and equivalence under a real pool are
covered by dedicated real-panel fixture tests in `tests/test_ldscore_parallelism.py`.

## 7. Risks and mitigations

- **BLAS oversubscription** (most likely performance footgun) → initializer pins
  BLAS threads to 1 by default; documented override.
- **Spawn re-import cost / serializability** → worker is module-level, inputs are
  dataclasses + DataFrames; an equivalence test exercises the real cross-process
  transfer path.
- **Peak memory** → unchanged from today (the sequential path already holds all
  surviving `ChromLDScoreResult`s before aggregating). Worker-local peak is one
  chromosome, same as today.
- **Non-determinism from completion order** → results re-ordered by the original
  chromosome list before aggregation.
