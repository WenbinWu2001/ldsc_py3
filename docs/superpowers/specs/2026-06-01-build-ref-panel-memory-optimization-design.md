# build-ref-panel Memory Optimization: Six Refactors

**Date:** 2026-06-01
**Implementation status:** designed; not yet implemented on `restructure`.
**Scope:** `src/ldsc/_kernel/ref_panel_builder.py`, `src/ldsc/ref_panel_builder.py`,
`src/ldsc/_kernel/ldscore.py`, `docs/current/parquet-r2-format-and-read-pipeline.md`
**Companion execution plan:** approved plan at the session plan file (six-task roadmap).

---

## Problem

`build-ref-panel` reads PLINK genotypes and writes pairwise R² parquet tables.
With 1000 Genomes 30× deep-sequencing input (~30–100× more SNPs than HapMap3),
the builder hits out-of-memory before completion and the output parquets become
intractably large.

The dominant memory cost is the **pair buffer**: `yield_pairwise_r2_rows`
accumulates emitted pairs in `pending: dict[int, list[dict]]`
([ref_panel_builder.py:323-341](../../../src/ldsc/_kernel/ref_panel_builder.py)),
one Python `dict` per retained SNP pair (`{"i","j","R2","sign"}`). For a wide LD
window over tens of millions of SNPs this is billions of dict objects, each with
~hundreds of bytes of interpreter overhead. Secondary costs: float64 genotype
matrices in the sliding block, the `pd.DataFrame → pa.Table.from_pandas` write
path, and the full BED bitarray retained after SNP filtering.

This design specifies six refactors (Tier 1 & 2). Streaming BED reads (#7) are
deferred.

---

## Design Decisions (confirmed with user)

| # | Refactor | Decision |
|---|---|---|
| 1 | `--min-r2` threshold | **Opt-in, default disabled.** `min_r2 <= 0` → no filtering (exact backward-compat, negative unbiased R² still written). `min_r2 > 0` → emit only pairs whose **unbiased** R² `>= min_r2`. Record `ldsc:min_r2` schema metadata. |
| 2 | Compact pair buffers | Replace dict-of-dicts with per-left-index **numpy column arrays** (`i:int64`, `j:int64`, `R2:float32`, `sign:int8`). |
| 3 | Vectorized extraction | Replace per-element Python loops in the emit functions with numpy mask/gather operations. |
| 4 | float32 genotypes | **Builder-scoped only**, via a `dtype` parameter on `nextSNPs` (default `float64`). Legacy in-PLINK LD-score path stays bit-for-bit. |
| 5 | Direct PyArrow writes | Build `pa.Table` directly from numpy column buffers, bypassing the `pd.DataFrame` intermediary. |
| 6 | Early bitarray release | `del` the source bitarray inside `__filter_snps_maf__` once the filtered copy exists. |
| — | `sign` field | **KEEP.** Never reaches disk but is pinned by a unit test. Store compactly as `int8` internally; reconstruct `"+"`/`"-"` only at the yielded-dict boundary. |

---

## Public contracts (signatures)

```python
# _kernel/ref_panel_builder.py

def yield_pairwise_r2_rows(
    *,
    block_left: np.ndarray,
    snp_batch_size: int,
    standardized_snp_getter,
    m: int,
    n: int,
    min_r2: float = 0.0,            # NEW
) -> Iterator[dict[str, float | int | str]]: ...

def iter_pairwise_r2_rows(
    *, block_left, snp_batch_size, standardized_snp_getter, m, n,
    min_r2: float = 0.0,            # NEW
) -> list[dict[str, float | int | str]]: ...

def write_r2_parquet(
    *,
    pair_rows,
    reference_snp_table,
    path,
    genome_build,
    n_samples,
    snp_identifier,
    batch_size: int = 100_000,
    row_group_size: int = 50_000,
    min_r2: float = 0.0,            # NEW — recorded as metadata only
) -> str: ...
```

```python
# _kernel/ldscore.py — PlinkBEDFile
def nextSNPs(self, b, minorRef=None, dtype=np.float64): ...   # dtype NEW
```

**Yielded row dict is unchanged**: `{"i": int, "j": int, "R2": float, "sign": "+"|"-"}`.
This preserves the contract that `iter_pairwise_r2_rows` and its pinning test
([test_ref_panel_builder.py:253](../../../tests/test_ref_panel_builder.py)) assert.

---

## #1 — R² minimum-threshold filter

**Semantics.** A module-level helper decides filtering:

```python
def _r2_filter_active(min_r2: float) -> bool:
    return min_r2 is not None and min_r2 > 0.0
```

- `min_r2 <= 0` (default `0.0`): **no filtering**. Every pair in the LD window is
  emitted, including pairs with negative unbiased R². Byte-for-byte identical to
  current output.
- `min_r2 > 0`: emit a pair only when its **unbiased** R²
  (`_unbiased_r2_from_correlation(corr, n)`) is `>= min_r2`.

**Where the guard lives.** Inside the vectorized emit functions (#3): after
computing the unbiased-R² array for a batch, apply a boolean mask
`keep = r2 >= min_r2` when `_r2_filter_active(min_r2)`. This keeps the filter on
the array path (no per-row Python branching).

**Metadata.** In `write_r2_parquet`, add to `pa_meta` alongside the existing keys
([ref_panel_builder.py:668-673](../../../src/ldsc/_kernel/ref_panel_builder.py)):

```python
b"ldsc:min_r2": str(min_r2).encode("utf-8"),   # NEW
```

So a reader can tell whether a panel is complete (`min_r2 <= 0`) or thresholded.

**Plumbing.** `build_parser` ([src/ldsc/ref_panel_builder.py:1542+](../../../src/ldsc/ref_panel_builder.py))
gains `--min-r2` (`type=float, default=0.0`). The build config carries it to
`_build_chromosome`, which passes it to **both** `yield_pairwise_r2_rows(...)`
(actual filtering) and `write_r2_parquet(...)` (metadata).

**Scientific note.** Thresholding biases downstream LD scores **upward**, because
the parquet read path treats any missing pair as exactly `R²=0.0` while the
unbiased estimator can be negative or small-positive. This is why the default
disables filtering and the threshold is recorded in metadata — completeness is
guaranteed only when `ldsc:min_r2 <= 0`.

---

## #2 — Compact numpy pair buffers

**Replace** `pending: dict[int, list[dict]]` and the helpers `_stash_pair_rows` /
`_pop_pair_rows_before` ([ref_panel_builder.py:323-341](../../../src/ldsc/_kernel/ref_panel_builder.py)).

**New internal buffer** keyed by left index `i`, storing parallel numpy arrays:

```python
@dataclass
class _PendingPairs:
    j:    list[np.ndarray]   # int64 chunks
    r2:   list[np.ndarray]   # float32 chunks
    sign: list[np.ndarray]   # int8 chunks (+1 / -1)
```

A `dict[int, _PendingPairs]` maps `i → its pending right-endpoints`. Stash appends
array chunks (from the vectorized emit, #3). Flush (`_pop_pair_rows_before`) for
each finalized `i` concatenates its chunks, sorts by `j`, and yields one dict per
element:

```python
sign_str = "+" if s == 1 else "-"
yield {"i": i, "j": int(j_val), "R2": float(r2_val), "sign": sign_str}
```

**Memory effect.** A Python `dict` row (~232 B + str/key overhead) collapses to
`8 + 4 + 1 = 13` bytes of array storage per pair, plus amortized list overhead.
The `sign` string is materialized only transiently at yield.

---

## #3 — Vectorized pair extraction

**Rewrite** `_emit_cross_block_pairs` ([ref_panel_builder.py:274](../../../src/ldsc/_kernel/ref_panel_builder.py))
and `_emit_within_block_pairs` (298) to return numpy column arrays instead of
`list[dict]`.

**Cross-block** (carry-over `A` columns × batch `B` columns):

```python
def _emit_cross_block_pairs(correlation_matrix, a_indices, b_indices,
                            block_left, n_samples, min_r2):
    # valid[local_i, local_j] = (a_indices[i] >= block_left[b_indices[j]])
    #                           & (a_indices[i] < b_indices[j])
    bl   = block_left[b_indices]                      # (nb,)
    valid = (a_indices[:, None] >= bl[None, :]) & (a_indices[:, None] < b_indices[None, :])
    li, lj = np.nonzero(valid)
    corr = correlation_matrix[li, lj]
    r2   = _unbiased_r2_array(corr, n_samples)        # vectorized
    sign = np.where(corr >= 0, np.int8(1), np.int8(-1))
    i    = a_indices[li].astype(np.int64)
    j    = b_indices[lj].astype(np.int64)
    if min_r2 is not None and min_r2 > 0.0:
        keep = r2 >= min_r2
        i, j, r2, sign = i[keep], j[keep], r2[keep], sign[keep]
    return i, j, r2.astype(np.float32), sign
```

**Within-block** (upper-triangular over batch columns, subject to `block_left`):

```python
def _emit_within_block_pairs(correlation_matrix, b_indices, block_left,
                             n_samples, min_r2):
    nb = len(b_indices)
    li, lj = np.triu_indices(nb, k=1)                 # local i < local j
    gi, gj = b_indices[li], b_indices[lj]
    valid  = gi >= block_left[gj]
    li, lj, gi, gj = li[valid], lj[valid], gi[valid], gj[valid]
    corr = correlation_matrix[li, lj]
    r2   = _unbiased_r2_array(corr, n_samples)
    sign = np.where(corr >= 0, np.int8(1), np.int8(-1))
    if min_r2 is not None and min_r2 > 0.0:
        keep = r2 >= min_r2
        gi, gj, r2, sign = gi[keep], gj[keep], r2[keep], sign[keep]
    return gi.astype(np.int64), gj.astype(np.int64), r2.astype(np.float32), sign
```

**New vectorized helper** (alongside the scalar `_unbiased_r2_from_correlation`,
which stays for callers that pass scalars):

```python
def _unbiased_r2_array(correlation: np.ndarray, n_samples: int) -> np.ndarray:
    sq = correlation * correlation
    denom = n_samples - 2 if n_samples > 2 else n_samples
    return sq - (1.0 - sq) / denom
```

**Invariants preserved.** Output always has `i < j` (cross-block: `a_indices < b_indices`
by construction; within-block: `triu_indices(k=1)`). The per-left-index flush in #2
preserves non-decreasing `POS_1` ordering required by the read path.

---

## #4 — float32 carry-over matrix (builder-scoped)

**`nextSNPs` gains a `dtype` parameter** ([ldscore.py:611](../../../src/ldsc/_kernel/ldscore.py)),
defaulting to `np.float64`:

```python
def nextSNPs(self, b, minorRef=None, dtype=np.float64):
    ...
    X = np.array(list(slice.decode(self._bedcode)), dtype=dtype).reshape((b, nru)).T
    X = X[0:n, :]
    Y = np.zeros(X.shape, dtype=dtype)
    ...
```

All intermediate arithmetic (`mean`, `std`, standardization) follows `X`'s dtype.

**Builder passes float32.** At the `_build_chromosome` call site
([src/ldsc/ref_panel_builder.py:858-871](../../../src/ldsc/ref_panel_builder.py)),
wrap the getter so the sliding block requests float32:

```python
standardized_snp_getter=lambda b: geno.nextSNPs(b, dtype=np.float32)
```

Then carry-over `A` and batch `B` inside `yield_pairwise_r2_rows` are float32 from
birth (lowest peak memory).

**Legacy path untouched.** `ldScoreVarBlocks` / `__corSumVarBlocks__`
([ldscore.py:434-517](../../../src/ldsc/_kernel/ldscore.py)) call `nextSNPs`
without `dtype`, so they remain float64. Their `rtol=1e-7` equivalence tests
([test_ldscore_workflow.py:2717-2718](../../../tests/test_ldscore_workflow.py))
must stay green — this is the proof that #4 is correctly scoped.

**Precision.** Builder R² shifts ~1e-5–1e-6 relative vs float64 (variance/dot-product
in single precision). Acceptable for R² sums; validated with a tolerant
equivalence test (`rtol≈1e-4`), not a 1e-7 assert.

---

## #5 — Direct PyArrow record-batch writes

**`write_r2_parquet`** ([ref_panel_builder.py:629](../../../src/ldsc/_kernel/ref_panel_builder.py))
builds each batch as a `pa.Table` directly from numpy columns, replacing
`build_standard_r2_table(...) → pa.Table.from_pandas`.

Per batch: gather endpoint metadata via the i/j index arrays into the canonical
10-column layout (`_STANDARD_R2_COLUMNS`), producing pyarrow arrays:

```python
left  = reference_snp_table.iloc[i_arr]
right = reference_snp_table.iloc[j_arr]
table = pa.table({
    "CHR":   pa.array(left[chr_col].astype(str)),
    "POS_1": pa.array(left[pos_col].to_numpy(np.int64)),
    "POS_2": pa.array(right[pos_col].to_numpy(np.int64)),
    "SNP_1": pa.array(left["rsID"].astype(str)),
    "SNP_2": pa.array(right["rsID"].astype(str)),
    "A1_1":  pa.array(left["A1"].astype(str)),
    "A2_1":  pa.array(left["A2"].astype(str)),
    "A1_2":  pa.array(right["A1"].astype(str)),
    "A2_2":  pa.array(right["A2"].astype(str)),
    "R2":    pa.array(r2_arr, type=pa.float32()),
}).select(_STANDARD_R2_COLUMNS)
```

**Schema/metadata unchanged.** Same column order, dtypes, and all `ldsc:*` keys
(now including `ldsc:min_r2`). `build_standard_r2_table` ([534](../../../src/ldsc/_kernel/ref_panel_builder.py))
is retained for the empty-batch path and existing schema tests, or thinned to
share the column-gather logic.

**Sort invariant preserved** ([ref_panel_builder.py:685-693](../../../src/ldsc/_kernel/ref_panel_builder.py)):
the non-decreasing `POS_1` check is applied per batch, vectorized
(`np.all(POS_1[1:] >= POS_1[:-1])` plus cross-batch boundary check against
`prev_pos1`).

---

## #6 — Early release of original BED bitarray

In `PlinkBEDFile.__filter_snps_maf__` ([ldscore.py:584-609](../../../src/ldsc/_kernel/ldscore.py)),
once the filtered bitarray `y` is fully built, drop the reference to the source
`geno` argument before returning so the original is freed promptly (otherwise the
pre-filter and post-filter bitarrays coexist at peak):

```python
# after y is fully populated, before return
del geno
return (y, m_poly, n, kept_snps, freq)
```

**Confirm** no attribute retains the pre-filter bitarray after the `__init__`
rebind at [ldscore.py:411](../../../src/ldsc/_kernel/ldscore.py)
(`self.geno` is reassigned to the filtered result).

---

## Testing Strategy

### Preserve existing pins — `tests/test_ref_panel_builder.py`
- `test_iter_pairwise_r2_rows_emits_one_unordered_pair_per_window` (253): pair set,
  `i<j`, R² `assertAlmostEqual`, **`sign`** (287) all still pass.
- `test_yield_pairwise_r2_rows_emits_non_decreasing_left_index` (289).
- `test_build_standard_r2_table_uses_exact_schema`, sort-invariant, schema-metadata,
  n_samples/r2_bias tests — adjust only for the new `min_r2` param and the new
  `ldsc:min_r2` key.

### Preserve legacy equivalence — `tests/test_ldscore_workflow.py`
- In-PLINK LD-score `assert_allclose` at `rtol=1e-7` (2717-2718) **must stay green**,
  proving #4 is builder-scoped.

### New tests
1. **#1 backward-compat:** with `min_r2=0.0`, a fixture producing a negative
   unbiased R² pair still emits that pair (assert the negative-R² pair is present).
2. **#1 filtering + metadata:** with `min_r2=0.2`, sub-threshold pairs are dropped
   and `pyarrow.parquet.read_schema(path).metadata[b"ldsc:min_r2"] == b"0.2"`.
3. **#3 equivalence:** vectorized `_emit_cross_block_pairs` / `_emit_within_block_pairs`
   match the previous loop output (pairs, R², sign) on a small deterministic fixture.
4. **#4 precision oracle:** builder R² output (float32) matches a float64 reference
   (same panel via the in-PLINK `cor_sum` path) within `rtol≈1e-4`.
5. **#2 flush ordering:** assert yielded rows preserve non-decreasing `POS_1` and
   `j`-sorted order within each `i`.

### Full suite + manual
- `pytest` green; `ldsc build-ref-panel --help` shows `--min-r2`.
- Build a small panel; confirm parquet schema unchanged, row counts and `ldsc:*`
  metadata correct, and peak RSS reduced (`/usr/bin/time -l`).

---

## Risks & Mitigations

| Risk | Mitigation |
|---|---|
| `min_r2` default accidentally drops negative-R² pairs | Gate on `min_r2 > 0`; regression test #1 above. |
| Vectorized emit changes pair ordering | `i<j` by construction; per-`i` flush sorts by `j`; test #5. |
| float32 leaks into legacy path | `dtype` defaults to float64; only the builder lambda passes float32; legacy 1e-7 tests guard it. |
| `sign` accidentally dropped or leaked to parquet | Keep `sign` in yielded dict (test 287); never add to `_STANDARD_R2_COLUMNS`. |
| Direct-arrow path diverges from canonical schema | Reuse `_STANDARD_R2_COLUMNS` + `.select(...)`; schema-metadata test guards it. |
