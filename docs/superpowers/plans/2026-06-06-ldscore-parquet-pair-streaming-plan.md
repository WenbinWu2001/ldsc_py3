# LD-score Parquet Pair-Streaming Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the parquet LD-score backend's sliding-block GEMM imitation with a direct pair-streaming accumulation (`cor_sum = R @ annot` by scattering each stored pair once), eliminating the `O(b³/c)` query/scatter blowup.

**Architecture:** The parquet panel stores each within-window R² pair once as `(i<j, r2)`. Instead of tiling R into blocks and rebuilding each tile, stream every stored pair and accumulate `cor_sum = annot + U·annot + Uᵀ·annot` (U = strict-upper-triangular pair matrix), with an `i >= block_left[j]` window filter. Parquet-only; the PLINK genotype backend is untouched. Not bit-identical (~5e-4 vs the block path), guarded by exact small-fixture unit tests plus the existing parity suites at tolerance.

**Status:** Tasks 1–5 (the streaming structure with a `np.add.at` scatter) are **complete** (commits `d84d4cd`…`2fe1dcc`). A 2026-06-06 re-profile showed `np.add.at` is a ~10× regression on chr22 (89% of time in unbuffered scatter); **Task 6** replaces it with the chunked float64 `scipy.sparse` CSR SpMM (benchmarked ~75× faster, exact). Tasks 1–5 below are kept for history; Task 6 is the active work.

**Tech Stack:** Python, numpy, `scipy.sparse` (chunked float64 CSR SpMM — Task 6), pyarrow parquet, pytest/unittest.

**Design source of truth:** `docs/superpowers/specs/2026-06-06-ldscore-parquet-pair-streaming-design.md` and `docs/current/parquet-r2-format-and-read-pipeline.md`.

**Environment for every command:**
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev
```

---

## File Structure

- Modify: `src/ldsc/_kernel/ldscore.py`
  - Add `SortedR2BlockReader.iter_all_pairs()` (decode every row group once).
  - Add module-level `_accumulate_pair_contributions(...)` and `ld_score_streaming_from_r2_reader(...)`.
  - Switch `compute_chrom_from_parquet` to call the streaming driver.
  - Delete (Task 5): `ld_score_var_blocks_from_r2_reader`, `cross_block_matrix`, `within_block_matrix`, `_query_union_arrays`, `_RowGroupLRUCache`, `configure_auto_row_group_cache`, `_sliding_query_index_windows`, `_row_group_indices_for_index_window`, `_get_decoded_row_group`, `log_row_group_cache_summary`, `_rg_bounds` construction.
- Modify: `tests/test_plink_io.py` (add streaming/`iter_all_pairs` tests; remove obsolete query/matrix/cache tests in Task 5).
- Already revised in the design phase (verify in Task 5): `docs/current/parquet-r2-format-and-read-pipeline.md` (§3.2–§3.6, §6) and `docs/current/ld-window-parquet-r2-sidecar-behavior.md` (read-side/Memory) now describe streaming.
- Modify: `design_map.md` (Task 5).

---

## Task 1: `iter_all_pairs()` on the reader

Streams every stored pair once by decoding each row group exactly once (no LRU cache, no window pruning).

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` (add method to `SortedR2BlockReader`, near `_decode_index_row_group`)
- Test: `tests/test_plink_io.py` (`IndexParquetRuntimeTest`)

- [ ] **Step 1: Write the failing test**

Add to `IndexParquetRuntimeTest`:

```python
    def test_iter_all_pairs_yields_every_stored_pair_once(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            snps = ["rs1", "rs2", "rs3", "rs4"]
            bps = [100, 120, 140, 160]
            panel = self._panel_meta(snps, bps)
            pairs = [
                {"i": 0, "j": 1, "R2": 0.4, "sign": "+"},
                {"i": 0, "j": 2, "R2": 0.2, "sign": "+"},
                {"i": 1, "j": 2, "R2": 0.6, "sign": "+"},
                {"i": 2, "j": 3, "R2": 0.5, "sign": "+"},
            ]
            r2 = self._write_index_panel(tmpdir, panel, pairs, row_group_size=2)
            reader = self._reader(r2, self._reader_meta(snps, bps))
            seen = {}
            for i, j, r2v in reader.iter_all_pairs():
                self.assertEqual(i.dtype, np.dtype("int32"))
                self.assertEqual(r2v.dtype, np.dtype("float32"))
                for a, b, v in zip(i.tolist(), j.tolist(), r2v.tolist()):
                    seen[(a, b)] = v
            self.assertEqual(sorted(seen), [(0, 1), (0, 2), (1, 2), (2, 3)])
            # int16-quantized R2: tolerance is the half-step
            np.testing.assert_allclose(
                [seen[k] for k in [(0, 1), (0, 2), (1, 2), (2, 3)]],
                [0.4, 0.2, 0.6, 0.5], rtol=0, atol=2e-5,
            )
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```bash
pytest tests/test_plink_io.py::IndexParquetRuntimeTest::test_iter_all_pairs_yields_every_stored_pair_once -v
```
Expected: FAIL — `AttributeError: ... has no attribute 'iter_all_pairs'`.

- [ ] **Step 3: Add the method**

In `SortedR2BlockReader`, immediately after `_decode_index_row_group`, add:

```python
    def iter_all_pairs(self):
        """Yield ``(i, j, r2)`` arrays for every row group, decoding each once.

        Streams the whole chromosome's stored pairs in IDX_1 order with no
        window pruning and no row-group cache: each row group is decoded a
        single time, remapped to retained indices, and dropped endpoints
        removed. Empty groups (all endpoints outside the analysis universe) are
        skipped.
        """
        if self._pf is None:
            raise LDSCInternalError(
                "LD-score parquet reader failed in SortedR2BlockReader.iter_all_pairs(): "
                "the parquet file handle is not initialized. Re-run with DEBUG logging "
                "and report the traceback."
            )
        for rg_index in range(self._pf.metadata.num_row_groups):
            group = self._decode_index_row_group(rg_index)
            if group.i.size == 0:
                continue
            yield group.i, group.j, group.r2
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```bash
pytest tests/test_plink_io.py::IndexParquetRuntimeTest::test_iter_all_pairs_yields_every_stored_pair_once -v
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/ldscore.py tests/test_plink_io.py
git commit -m "feat(ldscore): add iter_all_pairs streaming over row groups"
```

---

## Task 2: Streaming accumulator core (pure function)

The numeric heart, isolated for exact unit testing: diagonal seed, both orientations, duplicate-left-index correctness (`np.add.at`), and the `i >= block_left[j]` window filter.

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` (module-level function, near `ld_score_var_blocks_from_r2_reader`)
- Test: `tests/test_plink_io.py` (`IndexParquetRuntimeTest`)

- [ ] **Step 1: Write the failing tests**

Add to `IndexParquetRuntimeTest`:

```python
    def test_accumulate_full_window_matches_R_at_annot(self):
        annot = np.ones((4, 1), dtype=np.float32)
        cor_sum = annot.astype(np.float64).copy()  # diagonal R2 = 1
        block_left = np.zeros(4, dtype=np.int64)
        i = np.array([0, 0, 1, 2], dtype=np.int32)
        j = np.array([1, 2, 2, 3], dtype=np.int32)
        r2 = np.array([0.4, 0.2, 0.6, 0.5], dtype=np.float32)
        ld._accumulate_pair_contributions(cor_sum, i, j, r2, annot, block_left)
        # Row sums of symmetric R (unit diagonal):
        np.testing.assert_allclose(cor_sum[:, 0], [1.6, 2.0, 2.3, 1.5], atol=1e-6)

    def test_accumulate_handles_duplicate_left_index(self):
        # SNP 0 pairs with SNP 1 and SNP 2 — both must accumulate into cor_sum[0].
        annot = np.array([[1.0], [2.0], [3.0], [0.0]], dtype=np.float32)
        cor_sum = np.zeros((4, 1), dtype=np.float64)
        block_left = np.zeros(4, dtype=np.int64)
        i = np.array([0, 0], dtype=np.int32)
        j = np.array([1, 2], dtype=np.int32)
        r2 = np.array([0.5, 0.5], dtype=np.float32)
        ld._accumulate_pair_contributions(cor_sum, i, j, r2, annot, block_left)
        np.testing.assert_allclose(cor_sum[0, 0], 0.5 * 2.0 + 0.5 * 3.0)  # 2.5
        np.testing.assert_allclose(cor_sum[1, 0], 0.5 * 1.0)
        np.testing.assert_allclose(cor_sum[2, 0], 0.5 * 1.0)

    def test_accumulate_window_filter_drops_out_of_window_pairs(self):
        annot = np.ones((4, 1), dtype=np.float32)
        cor_sum = np.zeros((4, 1), dtype=np.float64)
        block_left = np.array([0, 0, 1, 2], dtype=np.int64)  # SNP2 window starts at 1
        i = np.array([0, 0, 1, 2], dtype=np.int32)
        j = np.array([1, 2, 2, 3], dtype=np.int32)
        r2 = np.array([0.4, 0.2, 0.6, 0.5], dtype=np.float32)
        # (0,2) dropped: 0 < block_left[2]=1. Others kept.
        ld._accumulate_pair_contributions(cor_sum, i, j, r2, annot, block_left)
        np.testing.assert_allclose(cor_sum[:, 0], [0.4, 1.0, 1.1, 0.5], atol=1e-6)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run:
```bash
pytest tests/test_plink_io.py::IndexParquetRuntimeTest::test_accumulate_full_window_matches_R_at_annot -v
```
Expected: FAIL — `AttributeError: module 'ldsc._kernel.ldscore' has no attribute '_accumulate_pair_contributions'`.

- [ ] **Step 3: Add the function**

In `src/ldsc/_kernel/ldscore.py`, just above `ld_score_var_blocks_from_r2_reader`, add:

```python
def _accumulate_pair_contributions(
    cor_sum: np.ndarray,
    i: np.ndarray,
    j: np.ndarray,
    r2: np.ndarray,
    annot: np.ndarray,
    block_left: np.ndarray,
) -> None:
    """Scatter-add one chunk of stored pairs ``(i<j, r2)`` into ``cor_sum``.

    Implements ``cor_sum += R @ annot`` for the off-diagonal pairs of one chunk:
    each pair contributes ``r2·annot[j]`` to row ``i`` and ``r2·annot[i]`` to row
    ``j`` (R is symmetric). Pairs outside the ldscore window are dropped via
    ``i >= block_left[j]``. ``np.add.at`` is required because ``i`` (and ``j``)
    contain repeated indices within a chunk; a buffered ``+=`` would collapse them.
    """
    keep = i >= block_left[j]
    if not keep.all():
        i, j, r2 = i[keep], j[keep], r2[keep]
    if i.size == 0:
        return
    np.add.at(cor_sum, i, r2[:, None] * annot[j])
    np.add.at(cor_sum, j, r2[:, None] * annot[i])
```

- [ ] **Step 4: Run the tests to verify they pass**

Run:
```bash
pytest tests/test_plink_io.py -k "accumulate" -v
```
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/ldscore.py tests/test_plink_io.py
git commit -m "feat(ldscore): add pair-streaming accumulation core"
```

---

## Task 3: Streaming driver + switch the caller

Add the chromosome-level driver, point `compute_chrom_from_parquet` at it, and prove it matches the existing block path within tolerance (both paths still present).

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` (add `ld_score_streaming_from_r2_reader`; edit `compute_chrom_from_parquet`)
- Test: `tests/test_plink_io.py` (`IndexParquetRuntimeTest`)

- [ ] **Step 1: Write the failing parity test**

Add to `IndexParquetRuntimeTest` (compares the new streaming driver to the existing block driver on a random fixture, within tolerance):

```python
    def test_streaming_matches_block_path(self):
        rng = np.random.default_rng(0)
        n = 12
        snps = [f"rs{k+1}" for k in range(n)]
        bps = [100 + 10 * k for k in range(n)]
        panel = self._panel_meta(snps, bps)
        # Dense banded pairs within a width-4 window.
        pairs = []
        for a in range(n):
            for b in range(a + 1, min(a + 4, n)):
                pairs.append({"i": a, "j": b, "R2": float(rng.uniform(0.05, 0.95)), "sign": "+"})
        with tempfile.TemporaryDirectory() as tmpdir:
            r2 = self._write_index_panel(tmpdir, panel, pairs, row_group_size=3)
            reader = self._reader(r2, self._reader_meta(snps, bps))
            block_left = np.array([max(0, k - 3) for k in range(n)], dtype=np.int64)
            annot = rng.uniform(0, 1, size=(n, 3)).astype(np.float32)
            block_out = ld.ld_score_var_blocks_from_r2_reader(
                block_left=block_left, snp_batch_size=2, annot=annot, block_reader=reader,
            )
            stream_out = ld.ld_score_streaming_from_r2_reader(
                block_left=block_left, annot=annot, block_reader=reader,
            )
            np.testing.assert_allclose(stream_out, block_out, rtol=1e-5, atol=1e-5)
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```bash
pytest tests/test_plink_io.py::IndexParquetRuntimeTest::test_streaming_matches_block_path -v
```
Expected: FAIL — `AttributeError: ... has no attribute 'ld_score_streaming_from_r2_reader'`.

- [ ] **Step 3: Add the streaming driver**

In `src/ldsc/_kernel/ldscore.py`, just below `_accumulate_pair_contributions`, add:

```python
def ld_score_streaming_from_r2_reader(
    block_left: np.ndarray,
    annot: np.ndarray,
    block_reader: SortedR2BlockReader,
) -> np.ndarray:
    """Compute ``cor_sum = R @ annot`` by streaming the parquet's stored R2 pairs.

    The diagonal (R2=1) seeds ``cor_sum`` from ``annot``; each stored within-window
    pair is then scattered into both endpoints via ``_accumulate_pair_contributions``,
    restricted to the ldscore window by ``block_left``. Replaces the sliding-block
    GEMM imitation: each stored pair is touched exactly once (O(nnz · n_a)).
    """
    m, _n_a = annot.shape
    block_left = np.asarray(block_left, dtype=np.int64)
    cor_sum = annot.astype(np.float64, copy=True)  # diagonal R2 = 1 for every SNP
    for i, j, r2 in block_reader.iter_all_pairs():
        _accumulate_pair_contributions(cor_sum, i, j, r2, annot, block_left)
    return np.asarray(cor_sum, dtype=np.float32)
```

- [ ] **Step 4: Run the parity test to verify it passes**

Run:
```bash
pytest tests/test_plink_io.py::IndexParquetRuntimeTest::test_streaming_matches_block_path -v
```
Expected: PASS. (If it fails on tolerance, inspect the diff — a true mismatch indicates a window-filter or orientation bug, not float drift.)

- [ ] **Step 5: Switch `compute_chrom_from_parquet` to the streaming driver**

In `compute_chrom_from_parquet`, replace the call:

```python
    combined_scores = ld_score_var_blocks_from_r2_reader(
        block_left=block_left,
        snp_batch_size=args.snp_batch_size,
        annot=combined_annot,
        block_reader=block_reader,
    )
```

with:

```python
    combined_scores = ld_score_streaming_from_r2_reader(
        block_left=block_left,
        annot=combined_annot,
        block_reader=block_reader,
    )
```

- [ ] **Step 6: Run the workflow + reader suites**

Run:
```bash
pytest tests/test_plink_io.py tests/test_ldscore_workflow.py -v 2>&1 | tail -20
```
Expected: PASS. (`test_streaming_matches_block_path` green; workflow LD-score tests green within their tolerances.)

- [ ] **Step 7: Commit**

```bash
git add src/ldsc/_kernel/ldscore.py tests/test_plink_io.py
git commit -m "perf(ldscore): stream R2 pairs for parquet LD scores"
```

---

## Task 4: Full parity verification (within tolerance)

Confirm end-to-end agreement with the genotype path and across identifier modes, and record/adjust tolerances if any parity test was exact-equality.

**Files:** none unless a tolerance needs widening (then `tests/...`).

- [ ] **Step 1: Run the parity-critical suites**

Run:
```bash
pytest tests/test_ldscore_workflow.py tests/test_ldscore_parallelism.py tests/test_ref_panel.py tests/test_ref_panel_builder.py -v 2>&1 | tail -20
```
Expected: PASS. The build→read parity (`test_hm3_chr22_subset_runs_direct_and_parquet_ldscore_paths`) and cross-mode parity gates exercise the streaming path end-to-end.

- [ ] **Step 2: If any parity test fails only on tolerance**

Inspect the max abs/rel diff. If it is ≤ ~1e-4 absolute on LD scores (consistent with float accumulation order, not a logic error), widen that test's tolerance and add an inline comment: `# pair-streaming changes accumulation order; drift is float-rounding (~1e-6), not algorithmic`. If the diff is larger, STOP — it indicates a window-filter/orientation bug; debug before proceeding. Re-run Step 1.

- [ ] **Step 3: Run the full suite**

Run:
```bash
pytest 2>&1 | tail -8
```
Expected: PASS (the obsolete block/cache tests still pass here because the old code is still present — they are removed in Task 5).

- [ ] **Step 4: Commit (only if a tolerance was adjusted)**

```bash
git add tests/
git commit -m "test(ldscore): tolerance for pair-streaming accumulation order"
```

---

## Task 5: Remove dead block/query/cache machinery + obsolete tests + docs

The streaming path makes the sliding-block query/scatter/cache code unreachable. Delete it in one cleanup commit, replace the obsolete tests, and update docs.

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` (delete members below)
- Modify: `tests/test_plink_io.py` (remove obsolete tests)
- Modify: `docs/current/parquet-r2-format-and-read-pipeline.md`, `design_map.md`

- [ ] **Step 1: Delete the now-dead production members**

In `src/ldsc/_kernel/ldscore.py`, delete in full:
- function `ld_score_var_blocks_from_r2_reader`
- class `_RowGroupLRUCache`
- `SortedR2BlockReader` methods: `_query_union_arrays`, `cross_block_matrix`, `within_block_matrix`, `configure_auto_row_group_cache`, `log_row_group_cache_summary`, `_get_decoded_row_group`, `_row_group_indices_for_index_window`, and the static `_sliding_query_index_windows`
- in `_init_index_path`: the `self._rg_bounds` construction loop and the `self._row_group_cache` initialization in `__init__` (and the `self._rg_bounds = []` / `self._row_group_cache = None` initializers)

Keep: `_DecodedR2RowGroup`, `_decode_index_row_group`, `iter_all_pairs`, `_resolve_r2_scale`, `_transform_r2`, `build_index_remap`, the coarse-row-group startup warning in `_init_index_path`.

- [ ] **Step 2: Remove obsolete tests**

In `tests/test_plink_io.py` (`IndexParquetRuntimeTest`), delete these methods (they exercise deleted code):
`test_index_within_block_matrix_correctness`, `test_index_cross_block_matrix_correctness`, `test_index_query_reads_only_overlapping_row_groups`, `test_query_union_arrays_returns_numpy_window_slice`, `test_query_union_arrays_empty_window_returns_empty_arrays`, `test_auto_row_group_cache_capacity_uses_adjacent_query_union`, `test_auto_row_group_cache_capacity_uses_actual_cm_block_left`, `test_cached_ld_score_matches_expected_and_preserves_absent_pairs`, `test_tiny_row_group_cache_preserves_results_but_reads_more`, `test_overlapping_windows_reuse_cached_decoded_row_groups`, and `test_streaming_matches_block_path` (its purpose is now served by the kept parity suites + the Task 2 exact unit tests).

Keep: `test_index_init_builds_rg_bounds` only if `_rg_bounds` is retained; if `_rg_bounds` was removed in Step 1, delete `test_index_init_builds_rg_bounds` too. Keep `test_index_init_raises_on_build_mismatch`, `test_index_init_warns_on_coarse_row_groups`, `test_index_init_hard_fails_on_multiple_paths`, `test_build_index_remap_*`, `test_iter_all_pairs_*`, `test_accumulate_*`.

- [ ] **Step 3: Verify nothing references the deleted members**

Run:
```bash
grep -rn "ld_score_var_blocks_from_r2_reader\|_RowGroupLRUCache\|_query_union_arrays\|cross_block_matrix\|within_block_matrix\|configure_auto_row_group_cache\|_row_group_indices_for_index_window\|_sliding_query_index_windows\|_get_decoded_row_group\|log_row_group_cache_summary" src/ldsc tests
echo "exit: $? (1 = no matches = good)"
```
Expected: `exit: 1` (no matches).

- [ ] **Step 4: Resolve `--snp-batch-size` (now a parquet-read no-op)**

Streaming does not use `snp_batch_size`; the parquet LD-score read no longer
consumes it (`build-ref-panel` still does). Decide and apply:
- **Default (recommended):** keep the CLI flag accepted, but document that it has
  no effect on parquet LD-score reads. Add a one-line note where `--snp-batch-size`
  is parsed in the `ldscore` workflow parser, and update the help text.

Confirm no remaining ldscore-read use:
```bash
grep -rn "snp_batch_size" src/ldsc/_kernel/ldscore.py src/ldsc/ldscore_calculator.py
```
Expected: only the CLI/config plumbing remains (no use inside the streaming read path).

- [ ] **Step 5: Sweep ALL design docs for the removed mechanism**

The streaming change touches more than the two docs revised in the design phase.
For each file below, replace any "decoded row-group cache", "sliding/window read
path", "dense block matrix", or "`snp_batch_size` sizes the cache" statement with
the streaming model (`iter_all_pairs` → `cor_sum = R · annot` scatter-add;
read-side RSS bounded by `cor_sum`, flat in `--ld-wind-*`; `snp_batch_size` does
not affect parquet reads):

- `docs/current/parquet-r2-format-and-read-pipeline.md` — verify (rewritten in design phase).
- `docs/current/ldscore-parquet-accumulation.md` — verify the math spec matches the landed function/symbol names.
- `docs/current/ld-window-parquet-r2-sidecar-behavior.md` — verify (rewritten in design phase).
- `docs/current/architecture.md` (line ~115: "sizes a decoded row-group cache once per chromosome from … `snp_batch_size`; cache hits only avoid parquet rereads").
- `docs/current/layer-structure.md` (line ~21: "the index-format parquet decoded row-group cache").
- `docs/current/io-argument-inventory.md` (line ~169: "`--snp-batch-size` … sizes its decoded row-group cache automatically").
- `docs/current/inference-genome-build.md` (line ~210: confirm the referenced reader method name, e.g. `_init_index_path`, is current — fix if stale).
- `design_map.md` (line ~91: "`SortedR2BlockReader` decode/window read path" → streaming read path / `iter_all_pairs` + `ld_score_streaming_from_r2_reader`).

Verification gate — no live references to the deleted mechanism remain in any design doc:
```bash
grep -rn "decoded row-group cache\|row group cache\|cross_block_matrix\|within_block_matrix\|ld_score_var_blocks\|_RowGroupLRUCache\|configure_auto_row_group_cache\|sliding-block\|window read path" docs/current design_map.md
echo "exit: $? (1 = no matches = good)"
```
Expected: `exit: 1`.

- [ ] **Step 6: Run the full suite**

Run:
```bash
pytest 2>&1 | tail -8
```
Expected: PASS (no references to deleted code; streaming path covered by kept tests).

- [ ] **Step 7: Commit**

```bash
git add src/ldsc/_kernel/ldscore.py src/ldsc/ldscore_calculator.py tests/test_plink_io.py docs/current design_map.md
git commit -m "refactor(ldscore): remove block/query/cache machinery for streaming"
```

---

## Self-Review Notes

- **Spec coverage:** §2 streaming idea → Tasks 1-3; §3 window filter → Task 2 filter test; §4 numerics (exact small + parity at tolerance) → Tasks 2/4; §5 memory (no cache) → Task 5 deletion; §7 parquet-only → no PLINK file touched; §8 affected modules → Tasks 3/5.
- **Correctness landmine:** repeated left indices within a chunk must accumulate (not overwrite) — the Task 2 `np.add.at` and the Task 6 CSR SpMM both sum shared rows; locked by `test_accumulate_handles_duplicate_left_index`. (Task 6 adds a second landmine: the SpMM must be **float64**, else ~3e-3 precision loss.)
- **Not bit-identical by design:** Task 3 parity test and Task 4 use tolerance (`rtol/atol ~1e-5`); a *large* diff means a logic bug (window filter / orientation), not drift — Task 4 Step 2 makes this distinction explicit.
- **Design docs kept current:** Task 5 Step 5 sweeps *all* `docs/current` files + `design_map.md` (not just the two revised in the design phase) with a grep gate that must return no live references to the removed cache/block/query mechanism (`architecture.md`, `layer-structure.md`, `io-argument-inventory.md`, `inference-genome-build.md` are explicitly listed). This enforces "all design docs up to date after the code lands."
- **`--snp-batch-size` becomes a parquet-read no-op** (Task 5 Step 4): the flag is kept but documented as not affecting parquet LD-score reads (still used by `build-ref-panel`). A later cleanup could deprecate it if desired.
- **Type consistency:** `iter_all_pairs` yields `(int32, int32, float32)`; `_accumulate_pair_contributions(cor_sum:float64, i, j, r2, annot:float32, block_left:int64)`; `ld_score_streaming_from_r2_reader(block_left, annot, block_reader) -> float32` (note: drops the `snp_batch_size` parameter the old driver had — the single caller is updated in Task 3 Step 5).
- **Ordering:** the old block driver is kept through Task 4 so the transitional parity test can compare against it; it is deleted only in Task 5.

---

## Task 6: Replace the `np.add.at` scatter with chunked float64 CSR SpMM

The streaming structure (Tasks 1–5) is correct but the `np.add.at` scatter is ~10× slower than the block path on chr22 (89% of time in numpy's unbuffered `ufunc.at`). Replace it with a chunked `scipy.sparse` CSR SpMM in **float64** — benchmarked ~75× faster than `np.add.at` and exact (a float32 SpMM accumulates row sums in float32 and loses ~3e-3; see `lessons.md`). `scipy` is already a kernel dependency.

**Current line numbers** (`src/ldsc/_kernel/ldscore.py`): `iter_all_pairs` L1273, `_accumulate_pair_contributions` L1295 (with `np.add.at` at L1316–1317), `ld_score_streaming_from_r2_reader` L1320, `compute_chrom_from_parquet` L1380. `scipy` is imported in sibling kernels (`regression.py`, `_jackknife.py`) but **not** in `ldscore.py` yet.

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` (`_accumulate_pair_contributions` body, `ld_score_streaming_from_r2_reader`, add `_CSR_CHUNK_PAIRS` constant)
- Test: `tests/test_plink_io.py` (`IndexParquetRuntimeTest`)

- [ ] **Step 1: Write the failing chunking-equivalence test**

`ld_score_streaming_from_r2_reader` will gain a `chunk_pairs` parameter (default `_CSR_CHUNK_PAIRS`). Different chunk sizes must agree to float rounding (chunking is exact by linearity). Add to `IndexParquetRuntimeTest`:

```python
    def test_streaming_chunk_size_invariant(self):
        rng = np.random.default_rng(1)
        n = 12
        snps = [f"rs{k+1}" for k in range(n)]
        bps = [100 + 10 * k for k in range(n)]
        panel = self._panel_meta(snps, bps)
        pairs = []
        for a in range(n):
            for b in range(a + 1, min(a + 4, n)):
                pairs.append({"i": a, "j": b, "R2": float(rng.uniform(0.05, 0.95)), "sign": "+"})
        with tempfile.TemporaryDirectory() as tmpdir:
            r2 = self._write_index_panel(tmpdir, panel, pairs, row_group_size=3)
            reader = self._reader(r2, self._reader_meta(snps, bps))
            block_left = np.array([max(0, k - 3) for k in range(n)], dtype=np.int64)
            annot = rng.uniform(0, 1, size=(n, 3)).astype(np.float32)
            tiny = ld.ld_score_streaming_from_r2_reader(block_left, annot, reader, chunk_pairs=2)
            whole = ld.ld_score_streaming_from_r2_reader(block_left, annot, reader, chunk_pairs=10**9)
            np.testing.assert_allclose(tiny, whole, rtol=0, atol=1e-6)
```

- [ ] **Step 2: Run it to verify it fails**

Run:
```bash
pytest tests/test_plink_io.py::IndexParquetRuntimeTest::test_streaming_chunk_size_invariant -v
```
Expected: FAIL — `ld_score_streaming_from_r2_reader()` has no `chunk_pairs` keyword.

- [ ] **Step 3: Implement the chunked float64 CSR scatter**

In `src/ldsc/_kernel/ldscore.py`, add the chunk constant near the top of the module (with the other module constants):

```python
_CSR_CHUNK_PAIRS = 16_000_000  # K: pairs buffered before one CSR SpMM step (~0.45 GiB/chunk)
```

Replace `_accumulate_pair_contributions` (L1295–1317) with a float64 CSR builder (same signature; `m` is derived from `cor_sum`):

```python
def _accumulate_pair_contributions(
    cor_sum: np.ndarray,
    i: np.ndarray,
    j: np.ndarray,
    r2: np.ndarray,
    annot: np.ndarray,
    block_left: np.ndarray,
) -> None:
    """Add one chunk of stored pairs ``(i<j, r2)`` to ``cor_sum`` via a float64 CSR SpMM.

    Builds the strict-upper-triangular sparse matrix ``U`` of the window-filtered
    chunk and accumulates ``cor_sum += U @ annot + U.T @ annot`` (R is symmetric).
    The SpMM is float64 (``U.data`` is float64 and the driver passes a float64
    ``annot``): scipy accumulates each row sum in the operand dtype, so a float32
    SpMM would lose ~3e-3 (see lessons.md). Window filter: ``i >= block_left[j]``.
    """
    from scipy import sparse

    keep = i >= block_left[j]
    if not keep.all():
        i, j, r2 = i[keep], j[keep], r2[keep]
    if i.size == 0:
        return
    m = cor_sum.shape[0]
    u = sparse.csr_matrix((r2.astype(np.float64, copy=False), (i, j)), shape=(m, m))
    cor_sum += u @ annot
    cor_sum += u.T @ annot
```

Replace `ld_score_streaming_from_r2_reader` (L1320–1336) with a chunk-buffering driver:

```python
def ld_score_streaming_from_r2_reader(
    block_left: np.ndarray,
    annot: np.ndarray,
    block_reader: SortedR2BlockReader,
    chunk_pairs: int = _CSR_CHUNK_PAIRS,
) -> np.ndarray:
    """Compute ``cor_sum = R @ annot`` by streaming the parquet's stored R2 pairs.

    The diagonal (R2=1) seeds ``cor_sum`` from ``annot``; stored within-window pairs
    are buffered into chunks of ~``chunk_pairs`` and accumulated with a float64 CSR
    SpMM (`_accumulate_pair_contributions`), restricted to the ldscore window by
    ``block_left``. Each stored pair is touched exactly once (O(nnz * n_a)).
    """
    block_left = np.asarray(block_left, dtype=np.int64)
    annot64 = annot.astype(np.float64, copy=False)
    cor_sum = annot64.copy()  # diagonal R2 = 1 for every SNP
    buf_i: list[np.ndarray] = []
    buf_j: list[np.ndarray] = []
    buf_r2: list[np.ndarray] = []
    buffered = 0

    def flush() -> None:
        nonlocal buffered
        if buffered == 0:
            return
        _accumulate_pair_contributions(
            cor_sum,
            np.concatenate(buf_i),
            np.concatenate(buf_j),
            np.concatenate(buf_r2),
            annot64,
            block_left,
        )
        buf_i.clear(); buf_j.clear(); buf_r2.clear()
        buffered = 0

    for i, j, r2 in block_reader.iter_all_pairs():
        buf_i.append(i); buf_j.append(j); buf_r2.append(r2)
        buffered += int(i.size)
        if buffered >= chunk_pairs:
            flush()
    flush()
    return np.asarray(cor_sum, dtype=np.float32)
```

- [ ] **Step 4: Run the accumulator + chunking + reader tests**

Run:
```bash
pytest tests/test_plink_io.py -k "accumulate or streaming or iter_all_pairs" -v
```
Expected: PASS. The three existing `test_accumulate_*` tests pass unchanged (same signature; `m = cor_sum.shape[0]`; the n=4 fixtures are exact at that scale), and `test_streaming_chunk_size_invariant` passes.

- [ ] **Step 5: Run the parity-critical + full suites**

Run:
```bash
pytest tests/test_ldscore_workflow.py tests/test_ldscore_parallelism.py tests/test_plink_io.py 2>&1 | tail -6
pytest 2>&1 | tail -6
```
Expected: PASS, including the build→read parity and cross-mode oracle (`atol≈5e-5`) — float64 CSR is *more* precise than the prior `np.add.at`, so these should not regress. If the oracle test tightens below tolerance, that is expected improvement, not a failure.

- [ ] **Step 6: Commit**

```bash
git add src/ldsc/_kernel/ldscore.py tests/test_plink_io.py
git commit -m "perf(ldscore): chunked float64 CSR scatter for parquet LD scores"
```

- [x] **Step 7: Re-profile on the server** — DONE (`7568f72`, `--threads 1`).

Measured: chr22 **45 s** (block parity, 10.7× vs `np.add.at`); chr6 **2 min 36 s**
(42.75× vs `np.add.at`, vs block DNF); `ufunc.at` 0%, SpMM+CSR build 5.5% — now
setup-bound (annotation parsing ≈25%); chr6 peak RSS 3.35 GiB; chr22 max |Δ| vs
block 4.88e-4 (below the ~2e-3 quantization floor). The scatter optimization is
complete; the remaining lever is annotation/identifier parsing setup (~11 s/chrom).

### Task 6 notes
- **Signature preserved:** `_accumulate_pair_contributions` keeps its 6-arg signature, so the Task 2 exact unit tests are unchanged; only the body (np.add.at → CSR) and the driver (chunk buffering) change.
- **float64 is mandatory, not optional:** `U.data` must be float64 and the driver passes `annot64`. A float32 SpMM loses ~3e-3 (lessons.md). This is the one easy way to get a *correct-looking but wrong* result.
- **`chunk_pairs` is a test seam, not a CLI flag:** the parameter exists so `test_streaming_chunk_size_invariant` can force many small chunks; production uses the `_CSR_CHUNK_PAIRS=16M` default. Do not expose it on the CLI.
- **Docs already updated** (commits `5641286`, `db0712b`): math spec §5/§6.5, format doc §3.4/§3.6, design spec §5/§7 describe chunked float64 CSR with `K=16M`. No further doc sweep needed for Task 6 beyond verifying these still match.
