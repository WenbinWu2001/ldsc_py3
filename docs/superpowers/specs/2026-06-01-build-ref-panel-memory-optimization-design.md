# build-ref-panel Memory & Speed Optimization

**Date:** 2026-06-01
**Implementation status:** **#1–#10 implemented on `restructure`** (six memory
refactors #1–#6 plus per-build BED release, zstd-L9 compression, vectorized decode
lookup; full suite green, 870 passed, 1 skipped). **#11 (int16 R² quantization +
BYTE_STREAM_SPLIT) planned** on `feat/int16-r2-quantization` — see the §`#11` section
and the companion plan `docs/superpowers/plans/2026-06-02-int16-r2-quantization-plan.md`.
**Scope:** `src/ldsc/_kernel/ref_panel_builder.py`, `src/ldsc/ref_panel_builder.py`,
`src/ldsc/_kernel/ldscore.py`, `docs/current/parquet-r2-format-and-read-pipeline.md`
**Companion execution plan:** approved plan at the session plan file (six-task roadmap).

---

## Implementation status by refactor

| # | Refactor | Commit | Status |
|---|---|---|---|
| 1 | `--min-r2` unbiased-R² floor (opt-in, default off) | `b80c3a8` | ✅ implemented |
| 2 | Compact numpy pair buffers | `3a636bc` | ✅ implemented |
| 3 | Vectorized pair extraction | `3a636bc` | ✅ implemented |
| 4 | float32 genotypes (builder-scoped) | `338837f` | ✅ implemented |
| 5 | Direct PyArrow record-batch writes | `d4b6194` | ✅ implemented |
| 6 | Early pre-filter BED bitarray release (MAF filter) | `1dc9ad7` | ✅ implemented |
| 8 | Per-build BED release across emitted genome builds | `1bf1452` | ✅ implemented |
| 9 | zstd compression (raised to **level 9**) | `1bf1452`, `45f0ec0` | ✅ implemented |
| 10 | Vectorized canonical decode endpoint-index lookup (read path) | `b112013` | ✅ implemented |
| 11 | int16 R² quantization + BYTE_STREAM_SPLIT (lossy, ~63% smaller) | — | 🔜 planned (`feat/int16-r2-quantization`) |
| 7 | Streaming BED reads | — | ⏸ deferred |

Refactors #1–#6 are detailed below; the session additions #8–#10 follow in
["Additional optimizations"](#additional-optimizations).

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
| 11 | int16 R² quantization | **Symmetric int16, scale 32767, lossy.** Store `R2` as `int16` (`round(R2·32767)`, clipped `[-32767, 32767]`) with `BYTE_STREAM_SPLIT`; reader dequantizes (`q/32767`) auto-detecting int16 columns. Clip unbiased R² at `1.0` upstream so the endpoint maps to `32767` and decodes to **exactly `1.0`**. Negatives (~15%) preserved. Float32 panels still read via auto-detect. |

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

## Additional optimizations

These were implemented in the same effort but sit outside the original six-task
roadmap. #8 and #9 are build-path; #10 is a read-path companion that pays off the
storage gains from #9 without slowing LD-score calculation.

### #8 — Per-build BED release across emitted genome builds

**Problem.** When a panel emits more than one genome build (e.g. source hg38 plus
an hg38→hg19 liftover), the per-chromosome build loop in
[`src/ldsc/ref_panel_builder.py`](../../../src/ldsc/ref_panel_builder.py)
constructs a fresh `PlinkBEDFile` for each emitted build. Without an explicit
release, the previous build's decoded genotype bitarray stays alive while the
next build's bitarray is decoded, so two full bitarrays briefly coexist at peak.

**Fix.** At the end of each iteration of `for build in _emitted_genome_builds(config)`,
drop the genotype reference before the next build loads:

```python
del geno
gc.collect()
```

This is distinct from #6: #6 frees the *pre-filter* bitarray *within* one
`PlinkBEDFile`; #8 frees an *entire* `PlinkBEDFile` *between* builds.

**Test.** `test_each_genome_build_releases_prior_bed_before_loading_next`
([tests/test_ref_panel_builder.py](../../../tests/test_ref_panel_builder.py))
patches `PlinkBEDFile` with a tracking wrapper that records, via `weakref` +
`gc.collect()`, how many prior instances are still alive at each construction;
it asserts `[0, 0]` (no prior BED alive when the next build starts).

### #9 — zstd compression at level 9

**Decision.** R² parquet column chunks are compressed with **zstd level 9** when
the codec is available in the local `pyarrow` build, falling back to `snappy`
otherwise (snappy has no level knob). No new CLI flag — this is a silent default.

```python
# write_r2_parquet, writer construction
if pa.Codec.is_available("zstd"):
    writer = pq.ParquetWriter(str(path), schema, compression="zstd", compression_level=9)
else:
    writer = pq.ParquetWriter(str(path), schema, compression="snappy")
```

**Rationale.** zstd yields ~1.5–2× smaller artifacts than the prior snappy default
with comparable decode speed. Level 9 trades a slower write (~2.4× vs level 1) for
~3.5% smaller files than level 1, with **identical read speed and memory** (zstd
decode is level-independent). The parquet footer records the codec but **not** the
level, so the reader needs no change and the level is observable only at write time.

**Tests.** `test_write_r2_parquet_uses_zstd_compression` reads the footer and
asserts every column chunk reports `ZSTD`; `test_write_r2_parquet_uses_zstd_level_9`
spies on the `pq.ParquetWriter` constructor (the only way to observe the
write-time level) and asserts `compression="zstd"`, `compression_level=9`.

### #10 — Vectorized canonical decode endpoint-index lookup (read path)

**Problem.** `SortedR2BlockReader._decode_canonical_row_group`
([ldscore.py](../../../src/ldsc/_kernel/ldscore.py)) resolved each row endpoint to
its retained-SNP index with a per-row `index_map.get(key, -1)` Python loop
(`np.fromiter(...)`). This interpreted loop dominated row-group decode time — the
codec (zstd) decompression is comparatively negligible thanks to the decoded
row-group LRU cache — so it is the relevant lever for read speed once #9 is in.

**Fix.** Two pure module-level helpers replace the loop:

```python
def _vectorized_pos_lookup(sorted_pos, query):
    # chr_pos base mode: binary search + equality check, -1 for non-matches.
    query = np.asarray(query, dtype=np.int64)
    if sorted_pos.size == 0:
        return np.full(query.shape, -1, dtype=np.int64)
    idx = np.searchsorted(sorted_pos, query)
    idx_clipped = np.clip(idx, 0, sorted_pos.size - 1)
    matched = sorted_pos[idx_clipped] == query
    return np.where(matched, idx, -1).astype(np.int64, copy=False)

def _vectorized_key_lookup(key_index, key_values, query):
    # rsID / allele-aware modes: hashtable join over the index_map keys.
    loc = key_index.get_indexer(np.asarray(query))
    result = np.full(loc.shape, -1, dtype=np.int64)
    found = loc >= 0
    result[found] = key_values[loc[found]]
    return result
```

`chr_pos` base mode reuses the already-sorted `self.pos` (zero extra memory). The
string-keyed modes build `self._index_keys` (a `pd.Index` over the `index_map`
keys) and `self._index_values` (a parallel `int64` array, so values stay correct
even when NaN keys were excluded from the map) **once per chromosome**, reused
across every row-group decode.

**Correctness rests on existing invariants.** `self.pos` is ascending and its
array index equals the SNP index — the same assumption the row-group windowing
already uses (`np.searchsorted(self.pos, …)`). Within a chromosome, retained
identifiers are unique (`validate_retained_identifier_uniqueness`), so each key
maps to exactly one SNP. Per-chromosome scoping means positions on different
chromosomes never share a lookup structure.

**Tests.** `VectorizedIndexLookupTest`
([tests/test_ldscore_workflow.py](../../../tests/test_ldscore_workflow.py)) pins
both helpers bit-for-bit to the dict loop across all four identity modes,
including not-found (`-1`), empty-query/empty-panel, and excluded-NaN-key cases.
The 100 existing workflow tests decode real parquet row groups through the new
path and assert LD-score values against references, guarding the wiring.

**Measured.** ~1.3× faster (chr_pos int keys) and ~2.2× faster (rsID string keys)
on the lookup itself, with negligible added peak memory and no change to decoded
values.

### #11 — int16 R² quantization + BYTE_STREAM_SPLIT

**Problem.** The `R2` float32 column dominates the index-format parquet — ~77% of
the file in the production panel. It is high-entropy (unbiased correction noise in
the mantissa LSBs) so generic compression barely touches it. It is the only large
lever left after #9/#10.

**Supplementary evidence (the data that cleared the on-hold).** A genome-wide
evaluation (all 22 autosomes, 473.5M pairs; notebook
`2026-06-01_int16_r2_quantization.ipynb`) measured the *downstream LD-score* impact
of symmetric int16 quantization, not just raw R² error:

| Metric | Value |
|---|---|
| Encoding | `stored = clip(round(R²·32767), −32767, 32767)`, `decoded = stored/32767` |
| Saturation (clipped values) | **0** across every chromosome (R² range `[−0.000312, 1.00001]`) |
| Per-pair max abs error | `1.5e-5` (= ½ step; pure rounding, no fat tail) |
| Per-pair mean abs error | `8e-6` (≈ step/4, textbook uniform-rounding fingerprint) |
| LD-score mean abs delta | `2.1e-4` |
| LD-score max abs delta (genome-wide) | `2.1e-3` (chr3) |

The per-pair ceiling (1.5e-5) is ~20× below the unbiased-correction noise floor
(`1/N ≈ 3.1e-4` at `N=3202`), so the quantization adds noise an order of magnitude
smaller than the statistical noise already in the estimate. Zero-mean independent
per-pair errors largely cancel when summed into a SNP's LD score (√N, not N
growth), collapsing the 1.5e-5 per-pair ceiling to a ~2e-3 LD-score ceiling —
negligible against GWAS χ² sampling variance. Per-chromosome QC is flat (no outlier
in error, saturation, or LD delta). **Verdict: effectively lossless for LD-score
regression.**

**Storage (genome-wide aggregate from the same notebook):**

| Variant | Size (MiB) | Reduction vs source float32 |
|---|---|---|
| source float32 | 2710.9 | — |
| float32 + BSS (lossless) | 1856.1 | 31.5% |
| int16 symmetric | 1076.8 | 60.3% |
| **int16 symmetric + BSS** | **1014.2** | **62.6%** ← chosen |

BSS adds only ~2.3 pts on top of int16 (vs 31.5 pts on float32) because int16 has
only two bytes and the high byte is already low-entropy (R² clusters near 0), but
it is lossless and free, so it is taken.

**Decision.** **int16 symmetric + BYTE_STREAM_SPLIT** (option b). Rejected: lossless
BSS-only (leaves ~840 MiB for precision below the noise floor); int8 (max_err 2e-3
is visible at the small-R² mass where ~59% of pairs live).

**Endpoint invariant (user-required).** The maximum R² value must round-trip to
**exactly `1.0`**. The symmetric scale guarantees this: clip at `1.0` upstream →
`round(1.0·32767) = 32767` → `32767/32767 = 1.0` exactly (32767 < 2²⁴, exactly
representable in float32, and `x/x = 1.0` is exact). Using `32767` (not `32768`) as
both scale and clip bound is what makes the endpoint exact; the `−32767` lower clip
keeps the mapping symmetric so any negative also round-trips by the same rule.

**Writer.**

```python
# module constants
R2_QUANT_SCALE = 32767            # int16 max; symmetric, endpoint-exact
R2_ENCODING = "int16_symmetric"

# _unbiased_r2_array: upper-clip only, keep negatives
sq = correlation * correlation
r2 = sq - (1.0 - sq) / denom
return np.minimum(r2, 1.0)

# _standard_r2_index_table: quantize instead of float32 R2
q = np.clip(np.rint(r2 * R2_QUANT_SCALE), -R2_QUANT_SCALE, R2_QUANT_SCALE).astype(np.int16)

# write_r2_parquet: schema + encoding + metadata
("R2", pa.int16())
column_encoding={"IDX_2": "DELTA_BINARY_PACKED", "R2": "BYTE_STREAM_SPLIT"}
use_dictionary=["IDX_1"]   # unchanged; keeps R2 off the dictionary path
b"ldsc:r2_encoding": b"int16_symmetric"
b"ldsc:r2_scale": b"32767"
```

`schema_version` stays `"1"` (the new keys are additive provenance, not a contract
break). The upstream clip means the quantizer's `[-32767, 32767]` clip never fires
on real data — it is defensive.

**Reader (auto-detect by column type).** At open, read `ldsc:r2_encoding` /
`ldsc:r2_scale` from `self._pf.schema_arrow.metadata`. Set `self._r2_scale` to the
scale when the on-disk `R2` column is an **integer** type, else `None`. In
`_decode_index_row_group`, when `self._r2_scale` is set, dequantize
`r2 = col.astype(float32) / scale` **before** `_transform_r2`. If the column is
int16 but metadata is missing, fall back to `32767` with a warning. float32 columns
take the existing path unchanged — old panels keep working.

`_transform_r2` raw branch also gets the `np.minimum(..., 1.0)` clip so external raw
float32 panels share the `R² ≤ 1` invariant. For our int16 panels `r2_bias` resolves
to `unbiased`, so `_transform_r2` is a no-op after dequant — order is safe regardless.

**Pin.** `environment.yml`: `pyarrow>=14` (int16 BYTE_STREAM_SPLIT support; verified
on 23.0.1). Below 14, BSS-on-int16 is unavailable.

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
6. **#8 per-build release:** `test_each_genome_build_releases_prior_bed_before_loading_next`
   asserts no prior `PlinkBEDFile` is alive when the next build starts (`[0, 0]`).
7. **#9 compression:** `test_write_r2_parquet_uses_zstd_compression` (footer codec ==
   `ZSTD`) and `test_write_r2_parquet_uses_zstd_level_9` (writer spy: level 9).
8. **#10 decode lookup:** `VectorizedIndexLookupTest` pins the vectorized lookups
   bit-for-bit to the dict loop across all four identity modes (incl. `-1` and
   NaN-key cases).
9. **#11 endpoint exactness (user-required):** an R² of `1.0` (and a corr-roundoff
   value `>1.0`) stores as int16 `32767` and the reader dequantizes to **exactly
   `1.0`** (`==`, not `assertAlmostEqual`).
10. **#11 schema + metadata:** written parquet has `R2` dtype `int16`,
    `ldsc:r2_encoding == b"int16_symmetric"`, `ldsc:r2_scale == b"32767"`, and the
    `R2` column chunk reports `BYTE_STREAM_SPLIT` encoding in the footer.
11. **#11 round-trip precision:** a known float32 R² array written then read back
    matches within per-pair `atol=2e-5`; a negative unbiased R² survives with sign.
12. **#11 backward-compat:** a hand-built **float32** `R2` parquet reads unchanged
    (no dequant; `self._r2_scale is None`).
13. **#11 clip:** `_unbiased_r2_array` upper-clips at `1.0` (a perfect-LD pair gives
    exactly `1.0`, never `>1`); negatives are not floored. `_transform_r2` raw branch
    applies the same upper clip.
14. **#11 parity (relaxed):** `ReferencePanelBuilderParityTest` /
    `IndexCrossModeParityTest` on the quantized path assert per-pair R² `atol=2e-5`
    and LD-score `atol=3e-3` (no longer bit-identical).

### Full suite + manual
- `pytest` green — **870 passed, 1 skipped** on `restructure`; `ldsc build-ref-panel
  --help` shows `--min-r2`.
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
| #11 int16 saturates on a smaller-N panel (R² pushed `>1`) | Upstream `min(R²,1.0)` clip + defensive `[-32767,32767]` quantizer clip; endpoint test pins `1.0→32767→1.0`. |
| #11 endpoint drifts off exactly `1.0` | Scale **and** clip bound are both `32767` (`x/x=1.0` exact, 32767 < 2²⁴); dedicated `==` test (#9 above). |
| #11 reader mis-reads old float32 panels as int16 | Dequant gated on integer column dtype; float32 → `_r2_scale=None`; backward-compat test (#12). |
| #11 BSS-on-int16 unsupported in older pyarrow | `environment.yml` pins `pyarrow>=14`; writer already warns + falls back on missing zstd. |
| #11 lossy values break a downstream consumer that reads R² directly | Per-pair error bounded `≤1.5e-5` (20× below correction floor); documented in `parquet-r2-format-and-read-pipeline.md` §2.3. |
