# LD-score Parquet Pair-Streaming — Design

**Status:** approved design (brainstorming resolved 2026-06-06)
**Scope:** parquet (`parquet_r2`) LD-score backend only. PLINK backend unchanged.
**Companion plan:** `docs/superpowers/plans/2026-06-06-ldscore-parquet-pair-streaming-plan.md`

---

## 0. Relationship to the current code

The parquet LD-score backend **currently** uses the sliding-block algorithm:
`ld_score_var_blocks_from_r2_reader` drives `cross_block_matrix` /
`within_block_matrix`, which query windows via `_query_union_arrays` over a
decoded-row-group LRU cache (`_RowGroupLRUCache`, `configure_auto_row_group_cache`)
and multiply dense `(b×c)` strips against `annot`. The recent dedup-removal and
`searchsorted` work happened *inside* this structure; the blocking and the cache
are both still present. This design **replaces that whole structure** for the
parquet backend with a single streaming pass — the block driver, the dense matrix
builders, the windowed query, and the LRU cache are all removed. The PLINK
(genotype) backend is a physically separate kernel and is untouched.

## 1. Problem

Profiling the parquet LD-score path (complete chr22 cProfile; censored chr6) showed
the cost is **not** GEMM and **not** decode:

| cProfile self-time (chr22) | seconds |
|---|---|
| `_query_union_arrays` | 6.42 |
| `cross_block_matrix` (scatter) | 4.59 |
| `read_row_group` | 0.80 |
| `_decode_index_row_group` | 0.41 |
| GEMM / `np.dot` | negligible |

chr6 (MHC, LD windows up to 21,465 SNPs) did not finish in 38 minutes.

**Root cause.** The parquet backend *imitates* the genotype backend's sliding-block
GEMM algorithm (`ld_score_var_blocks_from_r2_reader` — "Mirror the old LDSC
sliding-block accumulation"). For each block it re-queries and re-scatters an
`O(b)`-wide band once per `c`-wide batch → `O(b³/c)` per block in window width `b`.
The redundancy is inherent to the block tiling, which only existed because the
genotype backend must *recompute* R² per window. The parquet backend already has R²
as an explicit pair list, so the tiling is pure inherited cost.

## 2. Key idea

The LD-score accumulation is exactly:

```
cor_sum = R @ annot
```

where `R` is the `m×m` symmetric within-window R² matrix with unit diagonal. The
parquet stores every within-window pair exactly once as `(i, j, r2)` with `i < j`
(retained-matrix indices). Compute `R @ annot` by **streaming each stored pair once**:

```
cor_sum[k]  = annot[k]                     # diagonal R²(k,k) = 1
cor_sum[i] += r2 · annot[j]                # each stored pair, both orientations
cor_sum[j] += r2 · annot[i]
```

This is `O(nnz · n_a)` — the theoretical floor (each stored pair contributes once per
annotation). No blocks, no GEMM, no sliding-window queries, no row-group LRU cache.

## 3. Window semantics

The ldscore LD window (`--ld-wind-*`) may be **smaller** than the panel's build
window, so the stored pairs are a superset. A stored pair `(i, j)` with `i < j` is
inside the ldscore window iff `i >= block_left[j]` (the same mutual-membership rule
the block traversal enforced via `block_left`). The streaming accumulator applies
this filter per chunk. When the ldscore window equals the build window the filter is
a no-op.

## 4. Numerics

The accumulation order differs from the blocked GEMM, so results are **not
bit-identical** — expect ≈1e-6 drift (float32 products accumulated into a float64
`cor_sum`). This is well within the LD-score tolerances already used by the parity
tests (`atol≈2e-5` per pair; LD-score deltas ≤~2e-3 from int16 quantization alone).

Guardrails:
- **Exact** hand-computed unit tests on tiny fixtures (no float drift at small scale)
  for the diagonal, both orientations, duplicate left-index correctness
  (`np.add.at`), and the window filter.
- The existing **build→read parity** and **cross-mode parity** suites, run within
  tolerance, prove end-to-end agreement with the genotype path and across identifier
  modes.

## 5. Memory

`O(m · n_a)` for `cor_sum` plus one decoded row group at a time while streaming —
**lower** than the prior sliding-band decoded cache (~1.15 GiB on chr6). Streaming
reads each row group exactly once.

## 6. Performance expectation

- chr22: the ~11 s of query+scatter collapse to a single `O(nnz·n_a)` pass; the run
  becomes decode/IO-bound.
- chr6: removes the `O(b³/c)` blowup entirely; the MHC stops being pathological.
  Target: minutes, not 38+ minutes.

## 7. Scope and non-goals

**In scope (parquet backend only):**
- Replace `ld_score_var_blocks_from_r2_reader` with a streaming accumulator.
- Add `SortedR2BlockReader.iter_all_pairs()` (decode every row group once).
- Remove the now-dead sliding-window query + dense matrix builders + LRU cache.

**Out of scope / unchanged:**
- PLINK backend (`compute_chrom_from_plink` → genotype kernel) — physically separate
  code, never touched.
- Panel format, the index↔sidecar binding, regression, output writer.
- The injective-remap collapse check and `iter`-time decode/dequant (reused as-is).

**Deferred:** replacing `np.add.at` with a `scipy.sparse` CSR `R @ annot` is a later
optional optimization; the initial implementation uses `np.add.at`, which is simple,
obviously correct, and bounded-memory. (`scipy` is already a dependency.)

## 8. Affected modules

| Module | Change |
|---|---|
| `_kernel/ldscore.py` | add `iter_all_pairs()` + streaming accumulator; switch `compute_chrom_from_parquet`; delete block/query/LRU-cache machinery |
| `docs/current/parquet-r2-format-and-read-pipeline.md` | §3.2–§3.6 + §6: streaming replaces windowed queries, dense matrices, and the decoded-row-group cache |
| `docs/current/ldscore-parquet-accumulation.md` | math spec of the streaming accumulation (`C = R·A`); notation glossary |
| `docs/current/ld-window-parquet-r2-sidecar-behavior.md` | read-side/Memory: read-side RSS flat in the LD window |
| `docs/current/architecture.md`, `layer-structure.md`, `io-argument-inventory.md`, `inference-genome-build.md` | drop "decoded row-group cache" / `--snp-batch-size`-sizes-cache statements (Task 5 sweep) |
| `tests/test_plink_io.py` | replace obsolete query/matrix/cache tests with streaming + `iter_all_pairs` tests |
| `design_map.md` | update kernel ↔ design mapping |
