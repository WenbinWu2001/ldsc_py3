# LD-score R2 Query Optimization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the per-window pandas deduplication and DataFrame round-trip from the parquet-R2 LD-score read path, replacing them with an up-front injective-remap check and a numpy `searchsorted` window slice.

**Architecture:** The `SortedR2BlockReader` hot loop currently, for every sliding window, concatenates whole decoded row groups, masks them, builds a `DataFrame`, and runs a `groupby`-based dedup. We move the uniqueness guarantee to a one-time check in `build_index_remap` (reject any identifier-mode collapse, which also fixes a latent `retained_build_idx` last-write-wins bug), then make the window query return plain numpy arrays sliced by binary search on the sorted `i` array. The decoded-row-group LRU cache, row-group pruning, and the dense-matrix scatter are unchanged. All changes live inside per-worker reader code, so they compose cleanly with the existing cross-chromosome process pool.

**Tech Stack:** Python, numpy, pandas (only at metadata edges now), pyarrow parquet, pytest/unittest.

**Design source of truth:** `docs/current/parquet-r2-format-and-read-pipeline.md` (§2.2 pair uniqueness, §3.1 injective remap, §3.2–§3.4 numpy/searchsorted read mechanism). This plan implements that spec.

**Environment for every command below:**
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev
```

---

## File Structure

- Modify: `src/ldsc/_kernel/ldscore.py`
  - `build_index_remap` (~L386–419): add injective-remap rejection + bug fix.
  - `SortedR2BlockReader._query_union_rows_canonical_by_index` (~L1446–1468): replace with `_query_union_arrays` (numpy + searchsorted).
  - `SortedR2BlockReader._deduplicate_pairs` (~L1470–1491): delete.
  - `SortedR2BlockReader._empty_pair_rows` (~L1275–1284): delete.
  - `SortedR2BlockReader.cross_block_matrix` (~L1493–1527): consume numpy arrays.
  - `SortedR2BlockReader.within_block_matrix` (~L1529–1555): consume numpy arrays.
- Test: `tests/test_plink_io.py` (class `IndexParquetRuntimeTest`, helpers `_panel_meta`, `_reader_meta`, `_write_index_panel`, `_reader`).
- Docs already updated: `docs/current/parquet-r2-format-and-read-pipeline.md`.

---

## Task 1: Injective-remap rejection in `build_index_remap`

Adds the up-front, mode-aware collapse check that lets the per-window dedup be removed later, and fixes the latent last-write-wins assignment in `retained_build_idx`.

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py:386-419`
- Test: `tests/test_plink_io.py` (class `IndexParquetRuntimeTest`)

- [ ] **Step 1: Write the failing tests**

Add these two methods to `IndexParquetRuntimeTest` in `tests/test_plink_io.py`:

```python
    def test_build_index_remap_rejects_collapsing_sidecar(self):
        # Two distinct panel rows share an rsID; under rsid mode they collapse
        # onto the single retained "rsX", which must be a hard error.
        full = pd.DataFrame({"SNP": ["rsX", "rsY", "rsX"]})
        retained = pd.DataFrame({"SNP": ["rsX", "rsY"]})
        with self.assertRaises(LDSCInputError) as ctx:
            ld.build_index_remap(full, retained, "rsid")
        self.assertIn("collapse", str(ctx.exception).lower())

    def test_build_index_remap_injective_builds_retained_build_idx(self):
        # rsB is absent from the retained universe (maps to -1); rsA, rsC are kept.
        full = pd.DataFrame({"SNP": ["rsA", "rsB", "rsC"]})
        retained = pd.DataFrame({"SNP": ["rsA", "rsC"]})
        remap, retained_build_idx = ld.build_index_remap(full, retained, "rsid")
        np.testing.assert_array_equal(remap, [0, -1, 1])
        np.testing.assert_array_equal(retained_build_idx, [0, 2])
```

- [ ] **Step 2: Run the tests to verify they fail**

Run:
```bash
pytest tests/test_plink_io.py::IndexParquetRuntimeTest::test_build_index_remap_rejects_collapsing_sidecar -v
```
Expected: FAIL — no exception raised (current code silently overwrites `retained_build_idx`).

- [ ] **Step 3: Add the injectivity check**

In `src/ldsc/_kernel/ldscore.py`, replace the tail of `build_index_remap` (from `remap = retained_index.get_indexer(...)` through `return remap, retained_build_idx`) with:

```python
    remap = retained_index.get_indexer(full_keys).astype(np.int32, copy=False)

    valid = remap >= 0
    valid_targets = remap[valid]
    if np.unique(valid_targets).size != valid_targets.size:
        raise LDSCInputError(
            "ldscore could not align the reference panel to its sidecar: distinct "
            f"panel SNPs collapse onto the same retained SNP under mode '{mode}'. "
            "Most likely the panel sidecar has duplicate SNP identities at this "
            "identifier resolution. Use an allele-aware `--snp-identifier` mode, or "
            "rebuild the reference panel after removing duplicate SNP identities. "
            f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
        )

    m = len(retained_metadata)
    retained_build_idx = np.empty(m, dtype=np.int64)
    retained_build_idx[valid_targets] = np.nonzero(valid)[0]
    return remap, retained_build_idx
```

Then update the `build_index_remap` docstring's final sentence to add:

```
    The remap must be injective over retained indices: two distinct panel rows
    sharing an identity key under ``identifier_mode`` (a collapse) is a hard
    error, which also keeps ``retained_build_idx`` a well-defined permutation.
```

- [ ] **Step 4: Run the tests to verify they pass**

Run:
```bash
pytest tests/test_plink_io.py::IndexParquetRuntimeTest::test_build_index_remap_rejects_collapsing_sidecar tests/test_plink_io.py::IndexParquetRuntimeTest::test_build_index_remap_injective_builds_retained_build_idx -v
```
Expected: PASS (2 passed).

- [ ] **Step 5: Run the existing reader suite to confirm no regression**

Run:
```bash
pytest tests/test_plink_io.py -v
```
Expected: PASS (all existing `IndexParquetRuntimeTest` tests still green).

- [ ] **Step 6: Commit**

```bash
git add src/ldsc/_kernel/ldscore.py tests/test_plink_io.py
git commit -m "fix(ldscore): reject identifier collapse in build_index_remap"
```

---

## Task 2: Numpy `searchsorted` window query (`_query_union_arrays`)

Introduces the replacement query method alongside the old one. The old `_query_union_rows_canonical_by_index` and `_deduplicate_pairs` stay until Task 3 removes them, so this task is independently testable.

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` (add method to `SortedR2BlockReader`, near `_query_union_rows_canonical_by_index` ~L1468)
- Test: `tests/test_plink_io.py` (class `IndexParquetRuntimeTest`)

- [ ] **Step 1: Write the failing tests**

Add to `IndexParquetRuntimeTest` in `tests/test_plink_io.py`:

```python
    def test_query_union_arrays_returns_numpy_window_slice(self):
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
            i, j, values = reader._query_union_arrays(1, 3)  # window {1, 2}
            self.assertEqual(i.dtype, np.dtype("int32"))
            self.assertEqual(j.dtype, np.dtype("int32"))
            self.assertEqual(values.dtype, np.dtype("float32"))
            # Only the (1, 2) pair has both endpoints inside [1, 3).
            self.assertEqual(sorted(zip(i.tolist(), j.tolist())), [(1, 2)])
            np.testing.assert_allclose(values, [0.6], rtol=0, atol=2e-5)

    def test_query_union_arrays_empty_window_returns_empty_arrays(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            snps = ["rs1", "rs2", "rs3"]
            bps = [100, 120, 140]
            panel = self._panel_meta(snps, bps)
            pairs = [{"i": 0, "j": 1, "R2": 0.4, "sign": "+"}]
            r2 = self._write_index_panel(tmpdir, panel, pairs, row_group_size=2)
            reader = self._reader(r2, self._reader_meta(snps, bps))
            i, j, values = reader._query_union_arrays(2, 2)  # empty window
            self.assertEqual(i.size, 0)
            self.assertEqual(j.size, 0)
            self.assertEqual(values.size, 0)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run:
```bash
pytest tests/test_plink_io.py::IndexParquetRuntimeTest::test_query_union_arrays_returns_numpy_window_slice -v
```
Expected: FAIL with `AttributeError: ... has no attribute '_query_union_arrays'`.

- [ ] **Step 3: Add the `_query_union_arrays` method**

In `src/ldsc/_kernel/ldscore.py`, add this method to `SortedR2BlockReader` immediately after `_query_union_rows_canonical_by_index` (do not delete the old method yet):

```python
    def _query_union_arrays(self, start: int, stop: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(i, j, r2)`` numpy arrays for retained-index window ``[start, stop)``.

        Each selected decoded row group has ``i`` sorted ascending (the remap
        preserves build order), so the window's rows are located by binary
        search on ``i`` and only the slice with ``j`` in ``[start, stop)`` is
        kept. Canonical artifacts are duplicate-free (parquet pair uniqueness +
        the injective remap from ``build_index_remap``), so no deduplication is
        performed.
        """
        start = max(0, int(start))
        stop = min(int(stop), self.m)
        if stop <= start:
            return (
                np.empty(0, dtype=np.int32),
                np.empty(0, dtype=np.int32),
                np.empty(0, dtype=np.float32),
            )

        parts_i: list[np.ndarray] = []
        parts_j: list[np.ndarray] = []
        parts_r2: list[np.ndarray] = []
        for index in self._row_group_indices_for_index_window(start, stop):
            group = self._get_decoded_row_group(index)
            if group.i.size == 0:
                continue
            lo = int(np.searchsorted(group.i, start, side="left"))
            hi = int(np.searchsorted(group.i, stop, side="left"))
            if hi <= lo:
                continue
            gj = group.j[lo:hi]
            keep = (gj >= start) & (gj < stop)
            if not keep.any():
                continue
            parts_i.append(group.i[lo:hi][keep])
            parts_j.append(gj[keep])
            parts_r2.append(group.r2[lo:hi][keep])

        if not parts_i:
            return (
                np.empty(0, dtype=np.int32),
                np.empty(0, dtype=np.int32),
                np.empty(0, dtype=np.float32),
            )
        return (
            np.concatenate(parts_i),
            np.concatenate(parts_j),
            np.concatenate(parts_r2),
        )
```

- [ ] **Step 4: Run the tests to verify they pass**

Run:
```bash
pytest tests/test_plink_io.py::IndexParquetRuntimeTest::test_query_union_arrays_returns_numpy_window_slice tests/test_plink_io.py::IndexParquetRuntimeTest::test_query_union_arrays_empty_window_returns_empty_arrays -v
```
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/ldscore.py tests/test_plink_io.py
git commit -m "feat(ldscore): add numpy searchsorted R2 window query"
```

---

## Task 3: Rewire matrix builders to numpy; delete dedup and old query

Switches `cross_block_matrix` / `within_block_matrix` to `_query_union_arrays`, then removes `_query_union_rows_canonical_by_index`, `_deduplicate_pairs`, and `_empty_pair_rows`. Updates the one existing test that referenced the old query method.

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` (`cross_block_matrix`, `within_block_matrix`; delete three methods)
- Test: `tests/test_plink_io.py` (update `test_index_query_reads_only_overlapping_row_groups`; add `cross_block_matrix` correctness test)

- [ ] **Step 1: Update the existing query-pruning test and add a cross-block test**

In `tests/test_plink_io.py`, in `test_index_query_reads_only_overlapping_row_groups`, replace:

```python
            rows = reader._query_union_rows_canonical_by_index(0, 3)
            self.assertEqual(sorted(set(rg_idxs_used)), [0])
            self.assertEqual(len(rows), 2)
```

with:

```python
            i, j, r2 = reader._query_union_arrays(0, 3)
            self.assertEqual(sorted(set(rg_idxs_used)), [0])
            self.assertEqual(i.size, 2)
```

Then add this new correctness test to `IndexParquetRuntimeTest`:

```python
    def test_index_cross_block_matrix_correctness(self):
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
            # Block A = rows [0, 2), Block B = rows [2, 4): cross R2 of A×B.
            matrix = reader.cross_block_matrix(l_A=0, b=2, l_B=2, c=2)
            expected = np.array([[0.2, 0.0], [0.6, 0.0]], dtype=np.float32)
            self.assertEqual(matrix.dtype, np.dtype("float32"))
            np.testing.assert_allclose(matrix, expected, rtol=0, atol=2e-5)
```

- [ ] **Step 2: Run the new/updated tests to verify they fail**

Run:
```bash
pytest tests/test_plink_io.py::IndexParquetRuntimeTest::test_index_query_reads_only_overlapping_row_groups tests/test_plink_io.py::IndexParquetRuntimeTest::test_index_cross_block_matrix_correctness -v
```
Expected: FAIL — `test_index_query_reads_only_overlapping_row_groups` fails on the renamed method; `test_index_cross_block_matrix_correctness` passes already (cross_block_matrix still works via the old path). The failing one drives Step 3.

- [ ] **Step 3: Rewire `cross_block_matrix` to numpy arrays**

In `src/ldsc/_kernel/ldscore.py`, replace the body of `cross_block_matrix` (from `query_rows = ...` through the end of the pair-scatter block, i.e. the lines that call `_query_union_rows_canonical_by_index`, `_deduplicate_pairs`, and read `rows["i"]/["j"]/["R2"]`) so the method reads:

```python
    def cross_block_matrix(self, l_A: int, b: int, l_B: int, c: int) -> np.ndarray:
        """Build the dense cross-block R2 matrix used by the parquet backend."""
        if b <= 0 or c <= 0:
            return np.zeros((max(b, 0), max(c, 0)), dtype=np.float32)

        a_start, a_stop = l_A, l_A + b
        b_start, b_stop = l_B, l_B + c
        i, j, values = self._query_union_arrays(min(a_start, b_start), max(a_stop, b_stop))

        matrix = np.zeros((b, c), dtype=np.float32)
        if i.size > 0:
            left_i = (a_start <= i) & (i < a_stop) & (b_start <= j) & (j < b_stop)
            if left_i.any():
                matrix[i[left_i] - a_start, j[left_i] - b_start] = values[left_i]

            left_j = (a_start <= j) & (j < a_stop) & (b_start <= i) & (i < b_stop)
            if left_j.any():
                matrix[j[left_j] - a_start, i[left_j] - b_start] = values[left_j]

        overlap_start = max(a_start, b_start)
        overlap_stop = min(a_stop, b_stop)
        for idx in range(overlap_start, overlap_stop):
            matrix[idx - a_start, idx - b_start] = 1.0

        return matrix
```

- [ ] **Step 4: Rewire `within_block_matrix` to numpy arrays**

Replace `within_block_matrix` so it reads:

```python
    def within_block_matrix(self, l_B: int, c: int) -> np.ndarray:
        """Build the dense within-block R2 matrix used by the parquet backend."""
        if c <= 0:
            return np.zeros((0, 0), dtype=np.float32)

        b_start, b_stop = l_B, l_B + c
        # _query_union_arrays already restricts both endpoints to [b_start, b_stop).
        i, j, values = self._query_union_arrays(b_start, b_stop)

        matrix = np.zeros((c, c), dtype=np.float32)
        if i.size > 0:
            ii = i - b_start
            jj = j - b_start
            matrix[ii, jj] = values
            matrix[jj, ii] = values

        np.fill_diagonal(matrix, 1.0)
        return matrix
```

- [ ] **Step 5: Delete the dead methods**

In `src/ldsc/_kernel/ldscore.py`, delete these three methods in full:
- `_query_union_rows_canonical_by_index` (the old DataFrame-returning query)
- `_deduplicate_pairs`
- `_empty_pair_rows`

- [ ] **Step 6: Confirm no remaining references to the deleted methods**

Run:
```bash
grep -rn "_query_union_rows_canonical_by_index\|_deduplicate_pairs\|_empty_pair_rows" src/ldsc tests
```
Expected: no output (all references removed).

- [ ] **Step 7: Run the reader suite to verify it passes**

Run:
```bash
pytest tests/test_plink_io.py -v
```
Expected: PASS (all `IndexParquetRuntimeTest` tests green, including the updated pruning test and the new cross-block test).

- [ ] **Step 8: Commit**

```bash
git add src/ldsc/_kernel/ldscore.py tests/test_plink_io.py
git commit -m "perf(ldscore): drop per-window dedup; numpy searchsorted query"
```

---

## Task 4: Full regression and parity verification

Confirms LD scores are unchanged end-to-end (the optimization must be numerically identical for canonical artifacts) and that parallel runs still match sequential.

**Files:** none (verification only).

- [ ] **Step 1: Run the LD-score workflow and parallelism suites**

Run:
```bash
pytest tests/test_ldscore_workflow.py tests/test_ldscore_parallelism.py tests/test_ref_panel.py tests/test_ref_panel_builder.py -v
```
Expected: PASS. These include the cross-mode parity gate (one parquet, four identifier modes → identical LD scores) and the build→read parity test, which together prove the numpy path reproduces the previous results.

- [ ] **Step 2: Run the full suite**

Run:
```bash
pytest
```
Expected: PASS (no new failures vs. the pre-change baseline).

- [ ] **Step 3: Update `design_map.md` if it references the removed methods**

Run:
```bash
grep -n "_deduplicate_pairs\|_query_union_rows_canonical_by_index\|build_index_remap" design_map.md
```
If any line references the deleted methods, update it to point to `_query_union_arrays` / the injective `build_index_remap`. If there is no output, skip.

- [ ] **Step 4: Commit any doc/map updates**

```bash
git add -A
git commit -m "docs(ldscore): sync design_map with R2 query optimization"
```
(If Step 3 produced no changes, skip this commit.)

---

## Self-Review Notes

- **Spec coverage:** §2.2 pair uniqueness (assumed, not scanned) → Tasks 3 deletes the scan; §3.1 injective remap → Task 1; §3.2 sorted-`i` invariant → relied on in Task 2; §3.3 searchsorted slice + numpy return → Task 2; §3.4 numpy matrix build, no dedup → Task 3; cross-mode/build parity → Task 4.
- **Behavior change:** a collapsing *full sidecar* (duplicate identities where only one copy is retained) now errors up front in `build_index_remap` instead of (sometimes) erroring later inside `_deduplicate_pairs`. For canonical artifacts this path never triggers. The existing workflow test `test_package_parquet_metadata_duplicates_remain_invariant_failures` is unaffected — it fails earlier at metadata validation (`"non-unique SNP identifiers"`).
- **Type consistency:** `_query_union_arrays` returns `(int32, int32, float32)`; `cross_block_matrix`/`within_block_matrix` consume those directly; `build_index_remap` returns `(remap: int32, retained_build_idx: int64)` unchanged in signature.
- **No external callers:** `_deduplicate_pairs`, `_empty_pair_rows`, and `_query_union_rows_canonical_by_index` are private and used only within `cross_block_matrix`/`within_block_matrix` and one test (updated in Task 3).
