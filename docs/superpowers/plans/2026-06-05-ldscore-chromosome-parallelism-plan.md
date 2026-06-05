# ldscore cross-chromosome parallelism — Implementation Plan

> **Superseded note (final API):** this plan was executed, but the public knob
> was simplified post-implementation to a single integer `LDScoreConfig.threads`
> / `--threads` flag (joblib `n_jobs` convention: `1`=sequential default, `-1`=all
> cores, `-2`=all but one), affinity-aware via `_available_cpu_count`. The earlier
> `num_workers`/`parallel` two-knob design below is historical. See the design doc
> §3.1 for the final semantics.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compute per-chromosome LD scores concurrently on one machine via a process pool, controlled by a new `num_workers` knob, with output bit-identical to the sequential path.

**Architecture:** Add `LDScoreConfig.num_workers`. Extract the per-chromosome body into a module-level, spawn-safe worker `_compute_one_chromosome` returning a tagged `_ChromOutcome`. `LDScoreCalculator.run` resolves a worker count and either runs the worker inline (count == 1, exact current path) or fans out over a spawn `ProcessPoolExecutor`, then re-orders results by chromosome and aggregates/writes unchanged.

**Tech Stack:** Python `concurrent.futures.ProcessPoolExecutor`, `multiprocessing` spawn context, frozen dataclasses, pytest.

Design: `docs/superpowers/specs/2026-06-05-ldscore-chromosome-parallelism-design.md`

---

## File Structure

- `src/ldsc/config.py` — add `num_workers` field + validation to `LDScoreConfig` (frozen dataclass at `config.py:431`).
- `src/ldsc/ldscore_calculator.py` — module-level worker `_compute_one_chromosome`, `_ChromOutcome`, `_resolve_worker_count`, `_is_empty_intersection`, pool initializer `_init_worker`, rewritten `run` loop body (`:347-363`), CLI flag (`:810`), `_ldscore_config_from_args` (`:1378`).
- `tests/test_ldscore_parallelism.py` — new test module for all parallelism behavior.

No change to `outputs.py`, `_aggregate_chromosome_results`, or the regression layer.

---

## Task 1: Add `num_workers` to `LDScoreConfig`

**Files:**
- Modify: `src/ldsc/config.py:431` (dataclass fields + `__post_init__` near `:489`, docstring near `:471`)
- Test: `tests/test_ldscore_parallelism.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_ldscore_parallelism.py
import pytest
from ldsc.config import LDScoreConfig
from ldsc.errors import LDSCConfigError


def _cfg(**kw):
    base = dict(ld_wind_cm=1.0)
    base.update(kw)
    return LDScoreConfig(**base)


def test_num_workers_defaults_to_one():
    assert _cfg().num_workers == 1


def test_num_workers_negative_rejected():
    with pytest.raises(LDSCConfigError):
        _cfg(num_workers=-1)


def test_num_workers_zero_allowed_as_auto_sentinel():
    assert _cfg(num_workers=0).num_workers == 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_ldscore_parallelism.py -k num_workers -v`
Expected: FAIL — `TypeError: __init__() got an unexpected keyword argument 'num_workers'`.

- [ ] **Step 3: Add the field, validation, and docstring entry**

In `LDScoreConfig` fields block (after `whole_chromosome_ok: bool = False`):

```python
    num_workers: int = 1
```

In `__post_init__`, after the `snp_batch_size` check:

```python
        if self.num_workers < 0:
            raise LDSCConfigError(_positive_number_message("LDScoreConfig", "num_workers", self.num_workers))
```

In the class docstring, after the `whole_chromosome_ok` entry:

```
    num_workers : int, optional
        Number of worker processes for cross-chromosome parallelism on one
        machine. ``1`` (default) runs sequentially in-process. ``0`` means
        auto: ``min(os.cpu_count(), n_chromosomes)``. Values ``>1`` cap at the
        chromosome count. Output is identical regardless of this value.
```

Note: `_positive_number_message` rejects `<= 0`, but `0` is the valid auto sentinel, so use the explicit `< 0` check shown above (do not reuse the `<= 0` snp_batch_size pattern verbatim).

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_ldscore_parallelism.py -k num_workers -v`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/config.py tests/test_ldscore_parallelism.py
git commit -m "feat(ldscore): add num_workers field to LDScoreConfig"
```

---

## Task 2: Worker-count resolution helper

**Files:**
- Modify: `src/ldsc/ldscore_calculator.py` (add helper near the other module-level helpers, e.g. above `_chromosomes_from_bundle` at `:1432`)
- Test: `tests/test_ldscore_parallelism.py`

- [ ] **Step 1: Write the failing test**

```python
def test_resolve_worker_count():
    from ldsc.ldscore_calculator import _resolve_worker_count
    # explicit 1 stays sequential
    assert _resolve_worker_count(1, n_chromosomes=22) == 1
    # explicit >1 caps at chromosome count
    assert _resolve_worker_count(8, n_chromosomes=3) == 3
    assert _resolve_worker_count(4, n_chromosomes=22) == 4
    # auto (0) is positive and capped by chromosome count
    auto = _resolve_worker_count(0, n_chromosomes=2)
    assert 1 <= auto <= 2
    # single chromosome never spawns a pool
    assert _resolve_worker_count(8, n_chromosomes=1) == 1
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_ldscore_parallelism.py -k resolve_worker_count -v`
Expected: FAIL — `ImportError: cannot import name '_resolve_worker_count'`.

- [ ] **Step 3: Implement the helper**

Add near the top-of-module imports (verify `os` is imported; add `import os` if absent):

```python
def _resolve_worker_count(num_workers: int, n_chromosomes: int) -> int:
    """Resolve the configured ``num_workers`` against the chromosome count.

    ``0`` means auto (``min(os.cpu_count(), n_chromosomes)``). Any value is
    capped at ``n_chromosomes`` and floored at ``1`` so a single chromosome
    never spawns a pool.
    """
    if n_chromosomes <= 0:
        return 1
    if num_workers == 0:
        num_workers = os.cpu_count() or 1
    return max(1, min(num_workers, n_chromosomes))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_ldscore_parallelism.py -k resolve_worker_count -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/ldscore_calculator.py tests/test_ldscore_parallelism.py
git commit -m "feat(ldscore): add worker-count resolution helper"
```

---

## Task 3: Module-level chromosome worker and outcome type

**Files:**
- Modify: `src/ldsc/ldscore_calculator.py` (add `_ChromOutcome`, `_is_empty_intersection`, `_compute_one_chromosome`; refactor `_warn_and_skip_empty_intersection` at `:1571`)
- Test: `tests/test_ldscore_parallelism.py`

- [ ] **Step 1: Write the failing test**

This uses the existing LD-score parquet fixtures. Locate them first with `grep -rn "r2_dir\|resolve_r2_paths\|build-ref-panel" tests/ | head`; reuse whatever small parquet panel + annotation an existing parquet LD-score test uses (e.g. `tests/test_ldscore_calculator*.py`). The test builds the same inputs that test sets up, then drives one chromosome through the worker.

```python
def test_compute_one_chromosome_returns_success_outcome(small_parquet_ldscore_inputs):
    # small_parquet_ldscore_inputs: a fixture/helper returning
    #   (chrom, chrom_bundle, ref_panel_spec, ldscore_config, global_config, regression_snps)
    from ldsc.ldscore_calculator import _compute_one_chromosome
    chrom, chrom_bundle, spec, ldcfg, gcfg, reg = small_parquet_ldscore_inputs
    outcome = _compute_one_chromosome(chrom, chrom_bundle, spec, ldcfg, gcfg, regression_snps=reg)
    assert outcome.chrom == chrom
    assert outcome.skipped is False
    assert outcome.result is not None
    assert len(outcome.result.baseline_table) > 0
```

If no reusable fixture exists, build inputs inline in the test by copying the smallest existing parquet LD-score test's setup (panel dir + annotation), then derive `ref_panel.spec`, `_ldscore_config_from_args`/`LDScoreConfig`, `GlobalConfig`, and `regression_snps=None`.

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_ldscore_parallelism.py -k compute_one_chromosome -v`
Expected: FAIL — `ImportError: cannot import name '_compute_one_chromosome'`.

- [ ] **Step 3: Implement outcome type, predicate, and worker**

Add a serializable sentinel + worker state near the module helpers (`dataclass`
is already imported at the top of the file via `from dataclasses import dataclass, field`):

```python
_WORKER_UNSET = "__use_worker_global__"
_WORKER_STATE: dict[str, Any] = {}


@dataclass(frozen=True)
class _ChromOutcome:
    """Tagged result of one chromosome's LD-score computation in a worker."""

    chrom: str
    result: ChromLDScoreResult | None
    skipped: bool
    skip_message: str | None = None
```

Refactor `_warn_and_skip_empty_intersection` into a pure predicate plus a thin warner so the worker can classify without emitting a warning across the process boundary:

```python
def _is_empty_intersection(error: Exception, chrom: str) -> bool:
    """Return whether ``error`` is the recoverable empty-intersection case."""
    message = str(error)
    old_shape = message.startswith("No retained annotation SNPs remain on chromosome ")
    new_shape = message.startswith("ldscore retained no annotation SNPs on chromosome ")
    if not (old_shape or new_shape):
        return False
    if old_shape and not any(suffix in message for suffix in (" after parquet intersection.", " after PLINK intersection.")):
        return False
    if new_shape and "after intersecting with the " not in message:
        return False
    return True
```

Add the worker (rebuilds the ref panel locally; reuses the existing instance method `compute_chromosome`, which holds no mutable state):

```python
def _compute_one_chromosome(
    chrom: str,
    chrom_bundle,
    ref_panel_spec,
    ldscore_config: LDScoreConfig,
    global_config: GlobalConfig,
    regression_snps=_WORKER_UNSET,
) -> _ChromOutcome:
    """Compute one chromosome end-to-end. Returns a tagged outcome.

    ``regression_snps`` defaults to the pool-initializer-provided global so a
    large SNP set is transferred once per worker, not once per task. The inline
    (single-worker) caller passes it explicitly.
    """
    from ._kernel.ref_panel import RefPanelLoader

    if regression_snps == _WORKER_UNSET:
        regression_snps = _WORKER_STATE.get("regression_snps")
    ref_panel = RefPanelLoader(global_config).load(ref_panel_spec)
    calculator = LDScoreCalculator()
    try:
        result = calculator.compute_chromosome(
            chrom=chrom,
            annotation_bundle=chrom_bundle,
            ref_panel=ref_panel,
            ldscore_config=ldscore_config,
            global_config=global_config,
            regression_snps=regression_snps,
        )
    except (ValueError, LDSCInputError) as exc:
        if _is_empty_intersection(exc, chrom):
            return _ChromOutcome(chrom=chrom, result=None, skipped=True, skip_message=str(exc))
        raise
    return _ChromOutcome(chrom=chrom, result=result, skipped=False)
```

Note: `regression_snps == _WORKER_UNSET` uses `==` not `is`, because spawn re-imports the module so the sentinel object identity differs across processes; comparing by string value is spawn-safe. `RestrictionIdentityKeys`/`set` compared against a string returns `False`, which is correct.

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_ldscore_parallelism.py -k compute_one_chromosome -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/ldscore_calculator.py tests/test_ldscore_parallelism.py
git commit -m "feat(ldscore): add module-level per-chromosome worker"
```

---

## Task 4: Pool initializer (regression keys, logging, BLAS threads)

**Files:**
- Modify: `src/ldsc/ldscore_calculator.py` (add `_init_worker`)
- Test: `tests/test_ldscore_parallelism.py`

- [ ] **Step 1: Write the failing test**

```python
def test_init_worker_sets_state_and_blas_env(monkeypatch):
    from ldsc.ldscore_calculator import _init_worker, _WORKER_STATE
    monkeypatch.delenv("OMP_NUM_THREADS", raising=False)
    monkeypatch.delenv("OPENBLAS_NUM_THREADS", raising=False)
    _init_worker(regression_snps={"rs1", "rs2"}, log_level="WARNING")
    assert _WORKER_STATE["regression_snps"] == {"rs1", "rs2"}
    import os
    assert os.environ["OMP_NUM_THREADS"] == "1"
    assert os.environ["OPENBLAS_NUM_THREADS"] == "1"


def test_init_worker_respects_user_blas_env(monkeypatch):
    from ldsc.ldscore_calculator import _init_worker
    monkeypatch.setenv("OMP_NUM_THREADS", "4")
    _init_worker(regression_snps=None, log_level="WARNING")
    import os
    assert os.environ["OMP_NUM_THREADS"] == "4"  # not overridden
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_ldscore_parallelism.py -k init_worker -v`
Expected: FAIL — `ImportError: cannot import name '_init_worker'`.

- [ ] **Step 3: Implement the initializer**

```python
def _init_worker(regression_snps, log_level: str) -> None:
    """Initialize a pool worker: shared regression keys, logging, BLAS threads.

    Pins BLAS thread counts to 1 unless the user already set them, so ``W``
    worker processes do not oversubscribe ``W x BLAS_threads`` cores.
    """
    import os

    from ._logging import configure_package_logging

    _WORKER_STATE["regression_snps"] = regression_snps
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(var, "1")
    configure_package_logging(log_level)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_ldscore_parallelism.py -k init_worker -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/ldscore_calculator.py tests/test_ldscore_parallelism.py
git commit -m "feat(ldscore): add pool worker initializer with BLAS guard"
```

---

## Task 5: Rewrite `run` loop body to branch sequential vs pool

**Files:**
- Modify: `src/ldsc/ldscore_calculator.py:347-363` (the loop), imports at top of file
- Test: `tests/test_ldscore_parallelism.py`

- [ ] **Step 1: Write the failing equivalence + determinism tests**

```python
def test_parallel_output_matches_sequential(small_parquet_ldscore_run, tmp_path):
    # small_parquet_ldscore_run(num_workers, out_dir) -> path to written result dir.
    # Build it from the same inputs an existing parquet LD-score E2E test uses,
    # writing a canonical directory via LDScoreOutputConfig.
    import pandas as pd
    seq = small_parquet_ldscore_run(num_workers=1, out_dir=tmp_path / "seq")
    par = small_parquet_ldscore_run(num_workers=3, out_dir=tmp_path / "par")

    seq_base = pd.read_parquet(seq / "ldscore.baseline.parquet")
    par_base = pd.read_parquet(par / "ldscore.baseline.parquet")
    pd.testing.assert_frame_equal(seq_base, par_base)


def test_parallel_run_is_deterministic(small_parquet_ldscore_run, tmp_path):
    import pandas as pd
    a = small_parquet_ldscore_run(num_workers=3, out_dir=tmp_path / "a")
    b = small_parquet_ldscore_run(num_workers=3, out_dir=tmp_path / "b")
    pd.testing.assert_frame_equal(
        pd.read_parquet(a / "ldscore.baseline.parquet"),
        pd.read_parquet(b / "ldscore.baseline.parquet"),
    )
```

The `small_parquet_ldscore_run` helper constructs `LDScoreCalculator().run(...)` (or calls `run_ldscore(...)`) with `LDScoreConfig(..., num_workers=N)` and an `LDScoreOutputConfig(output_dir=out_dir)`, returning `out_dir`. Reuse the smallest existing multi-chromosome parquet fixture (need >= 2 chromosomes for the pool path to engage).

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_ldscore_parallelism.py -k "matches_sequential or deterministic" -v`
Expected: FAIL — `num_workers` is accepted but ignored (pool path not wired), so either both runs are sequential (test passes trivially — guard against this by also asserting the pool path executes) or `run` has no branch yet. To make the failure meaningful, first confirm the branch does not exist; the test will pass only after Step 3 wires it. If it passes before Step 3, add `assert hasattr(LDScoreCalculator, "_run_chromosomes")` to force red.

- [ ] **Step 3: Implement the branch**

Ensure these imports exist at the top of `ldscore_calculator.py`. `warnings` is
already imported (line 30); `os` is added in Task 2. Add the two still missing:

```python
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
```

Replace the loop body at `:347-363` with:

```python
        worker_count = _resolve_worker_count(ldscore_config.num_workers, len(chromosomes))
        outcomes = self._run_chromosomes(
            chromosomes=chromosomes,
            annotation_bundle=annotation_bundle,
            ref_panel=ref_panel,
            ldscore_config=ldscore_config,
            global_config=global_config,
            regression_snps=regression_snps,
            worker_count=worker_count,
        )
        chromosome_results: list[ChromLDScoreResult] = []
        for chrom in chromosomes:
            outcome = outcomes[chrom]
            if outcome.skipped:
                warnings.warn(f"Skipping chromosome {chrom}: {outcome.skip_message}", UserWarning, stacklevel=2)
                continue
            chromosome_results.append(outcome.result)
```

Add the dispatch method on `LDScoreCalculator` (e.g. right after `run`):

```python
    def _run_chromosomes(
        self,
        chromosomes,
        annotation_bundle,
        ref_panel,
        ldscore_config,
        global_config,
        regression_snps,
        worker_count: int,
    ) -> dict[str, _ChromOutcome]:
        """Compute every chromosome, sequentially or via a spawn process pool.

        Returns outcomes keyed by chromosome; callers re-order by the original
        chromosome list for deterministic aggregation.
        """
        ref_panel_spec = ref_panel.spec
        if worker_count == 1:
            outcomes: dict[str, _ChromOutcome] = {}
            for chrom in chromosomes:
                chrom_bundle = _slice_annotation_bundle(annotation_bundle, chrom)
                outcomes[chrom] = _compute_one_chromosome(
                    chrom, chrom_bundle, ref_panel_spec, ldscore_config, global_config,
                    regression_snps=regression_snps,
                )
            return outcomes

        ctx = mp.get_context("spawn")
        outcomes = {}
        with ProcessPoolExecutor(
            max_workers=worker_count,
            mp_context=ctx,
            initializer=_init_worker,
            initargs=(regression_snps, global_config.log_level),
        ) as pool:
            futures = {
                pool.submit(
                    _compute_one_chromosome,
                    chrom,
                    _slice_annotation_bundle(annotation_bundle, chrom),
                    ref_panel_spec,
                    ldscore_config,
                    global_config,
                ): chrom
                for chrom in chromosomes
            }
            try:
                for future in futures:
                    outcome = future.result()
                    outcomes[outcome.chrom] = outcome
            except mp.ProcessError as exc:  # BrokenProcessPool subclasses this
                raise LDSCInternalError(
                    "ldscore parallel chromosome computation aborted: a worker process "
                    "crashed (most likely out-of-memory or a native library fault). "
                    "Re-run with a lower `--num-workers`, or `--num-workers 1` to isolate "
                    "the failing chromosome."
                ) from exc
        return outcomes
```

Note: iterating `for future in futures` (submission order) and keying outcomes by `outcome.chrom` makes ordering independent of completion order; the parent loop re-reads them in `chromosomes` order. The single-worker branch passes `regression_snps` explicitly (no initializer); the pool branch relies on `_init_worker` setting the worker global.

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_ldscore_parallelism.py -k "matches_sequential or deterministic" -v`
Expected: PASS.

- [ ] **Step 5: Run the full existing LD-score suite (regression guard)**

Run: `pytest tests/ -k ldscore -v`
Expected: PASS — default `num_workers=1` path unchanged.

- [ ] **Step 6: Commit**

```bash
git add src/ldsc/ldscore_calculator.py tests/test_ldscore_parallelism.py
git commit -m "feat(ldscore): parallelize chromosome computation across workers"
```

---

## Task 6: Skip-path equivalence under the pool

**Files:**
- Test: `tests/test_ldscore_parallelism.py`

- [ ] **Step 1: Write the failing test**

```python
def test_skip_path_matches_under_pool(parquet_inputs_with_one_empty_chrom, tmp_path):
    # Fixture: a multi-chromosome input where exactly one chromosome's annotation
    # has no intersection with the panel (e.g. annotation rows removed for chrom 'X').
    import pandas as pd
    seq = parquet_inputs_with_one_empty_chrom(num_workers=1, out_dir=tmp_path / "seq")
    par = parquet_inputs_with_one_empty_chrom(num_workers=3, out_dir=tmp_path / "par")
    pd.testing.assert_frame_equal(
        pd.read_parquet(seq / "ldscore.baseline.parquet"),
        pd.read_parquet(par / "ldscore.baseline.parquet"),
    )
```

Build the empty-intersection chromosome by dropping all annotation SNPs for one chromosome from the small fixture so `_is_empty_intersection` triggers, while >=2 chromosomes survive.

- [ ] **Step 2: Run test to verify it fails or errors**

Run: `pytest tests/test_ldscore_parallelism.py -k skip_path -v`
Expected: PASS if Task 5 is correct; if it FAILs (e.g. a worker raises instead of returning a skip outcome, or warning side effects differ), fix `_compute_one_chromosome`/`_run_chromosomes` until both paths produce identical surviving output. The test is the gate that the skip classification works across the process boundary.

- [ ] **Step 3: Commit**

```bash
git add tests/test_ldscore_parallelism.py
git commit -m "test(ldscore): verify empty-intersection skip parity under pool"
```

---

## Task 7: CLI flag `--num-workers`

**Files:**
- Modify: `src/ldsc/ldscore_calculator.py:810` (argparse), `:1378` (`_ldscore_config_from_args`)
- Test: `tests/test_ldscore_parallelism.py`

- [ ] **Step 1: Write the failing test**

```python
def test_cli_num_workers_reaches_config():
    from ldsc.ldscore_calculator import build_parser, _ldscore_config_from_args
    parser = build_parser()  # use the actual parser factory name found below
    args = parser.parse_args([
        "--output-dir", "x", "--ld-wind-cm", "1", "--num-workers", "2",
    ])
    cfg = _ldscore_config_from_args(args)
    assert cfg.num_workers == 2


def test_cli_num_workers_defaults_to_one():
    from ldsc.ldscore_calculator import build_parser, _ldscore_config_from_args
    parser = build_parser()
    args = parser.parse_args(["--output-dir", "x", "--ld-wind-cm", "1"])
    assert _ldscore_config_from_args(args).num_workers == 1
```

Before writing, run `grep -n "def .*parser\|add_argument(\"--snp-batch-size\"" src/ldsc/ldscore_calculator.py` to find the exact parser factory function name and use it in the test instead of `build_parser` if it differs.

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_ldscore_parallelism.py -k cli_num_workers -v`
Expected: FAIL — `--num-workers` is an unrecognized argument.

- [ ] **Step 3: Add the flag and plumb it**

Next to the `--snp-batch-size` argument (`:810`):

```python
    parser.add_argument("--num-workers", default=1, type=int, help="Worker processes for cross-chromosome parallelism. 1=sequential, 0=auto (all cores).")
```

In `_ldscore_config_from_args` (`:1378`), add to the `LDScoreConfig(...)` call:

```python
        num_workers=getattr(args, "num_workers", 1),
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_ldscore_parallelism.py -k cli_num_workers -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/ldscore_calculator.py tests/test_ldscore_parallelism.py
git commit -m "feat(ldscore): add --num-workers CLI flag"
```

---

## Task 8: Documentation and full-suite verification

**Files:**
- Modify: `docs/current/` data-flow/architecture doc that describes the LD-score run loop; `design_map.md`; `tutorials/` LD-score tutorial if it lists CLI flags; `lessons.md` only if a non-obvious bug was hit during implementation.

- [ ] **Step 1: Update docstrings**

Confirm `LDScoreCalculator.run` docstring and the new `_run_chromosomes`, `_compute_one_chromosome`, `_init_worker`, `_resolve_worker_count` all have complete docstrings (use `my-skills:fun-doc` if any are thin).

- [ ] **Step 2: Update design + tutorial docs**

Add a short "Parallelism" note to the LD-score section of the relevant `docs/current/` doc: single-machine process pool over chromosomes, `--num-workers`, output unchanged, BLAS-thread guard. Update `design_map.md` to map the design/plan to the new functions. If the LD-score tutorial enumerates flags, add `--num-workers` with a one-line example.

- [ ] **Step 3: Run the entire test suite**

Run: `pytest`
Expected: PASS (all green, including the new `tests/test_ldscore_parallelism.py` and the unchanged LD-score suite).

- [ ] **Step 4: Manual smoke check (real multiprocessing)**

Run an actual 2+-chromosome parquet LD-score job twice and diff:

```bash
ldsc ldscore --r2-dir <panel> --annot <annot> --ld-wind-cm 1 --output-dir /tmp/seq --num-workers 1
ldsc ldscore --r2-dir <panel> --annot <annot> --ld-wind-cm 1 --output-dir /tmp/par --num-workers 4
python -c "import pandas as pd; pd.testing.assert_frame_equal(pd.read_parquet('/tmp/seq/ldscore.baseline.parquet'), pd.read_parquet('/tmp/par/ldscore.baseline.parquet')); print('identical')"
```

Expected: prints `identical`. Note wall-clock difference between the two runs.

- [ ] **Step 5: Commit**

```bash
git add docs/ design_map.md tutorials/ 2>/dev/null; git commit -m "docs(ldscore): document cross-chromosome parallelism"
```

---

## Self-Review notes (for the implementer)

- **Spec coverage:** `num_workers` field + validation (T1), auto/cap resolution (T2), module-level worker + skip classification (T3), initializer with BLAS guard + shared regression keys (T4), sequential/pool branch with deterministic ordering + crash handling (T5), skip parity (T6), CLI (T7), docs + real-multiprocessing smoke (T8). All §3–§6 design items are covered.
- **Fixture dependency:** Tasks 3/5/6 depend on a small, multi-chromosome parquet LD-score fixture. Before Task 3, run `grep -rn "ldscore.baseline.parquet\|r2_dir=" tests/ | head` to find the nearest existing one and adapt it; do not invent a new panel format.
- **Spawn-safety:** the worker, outcome type, initializer, and sentinel comparison (`==`, not `is`) are all spawn-safe; the equivalence test in Task 5 exercises the real pool, not a mock.
