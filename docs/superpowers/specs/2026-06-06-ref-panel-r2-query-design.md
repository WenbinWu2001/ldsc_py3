# Reference-Panel R² Pair Query and R²→r Conversion

**Date:** 2026-06-06
**Status:** Approved design (pre-implementation)
**Topic:** `ref-panel-r2-query`
**Companion plan:** `docs/superpowers/plans/2026-06-06-ref-panel-r2-query-plan.md`

---

## 1. Motivation

`build-ref-panel` writes pairwise-R² reference panels in the canonical **index
format** (`docs/current/parquet-r2-format-and-read-pipeline.md`): per chromosome,
a `chrN_r2.parquet` of `(IDX_1, IDX_2, R2, SIGN)` rows paired with an
authoritative `chrN_meta.tsv.gz` sidecar that defines the index space.

Today the only consumer of that artifact is `SortedR2BlockReader`, which does a
**full streaming sweep** (`iter_all_pairs`) to accumulate LD scores. There is no
way to ask *"what is the R² for this specific list of SNP pairs?"* — a common
user need (e.g. checking LD between candidate variants, building a custom LD
matrix for a locus, sanity-checking the panel).

This feature adds two user-facing capabilities:

1. **Pair query** — a vectorized API + CLI command that, given a list of SNP
   pairs, returns the stored (adjusted/unbiased) R² for each pair, plus a
   harmonized correlation sign and a per-row status.
2. **R²→r conversion** — a standalone function that converts adjusted R² back to
   a signed Pearson correlation `r`, using the reference sample size recorded in
   the parquet metadata.

## 2. Scope

### In scope
- A reusable `R2Panel` handle that opens a panel once and serves many queries.
- A one-shot functional wrapper `query_r2(...)`.
- A pure, vectorized converter `unbiased_r2_to_pearson_r(...)`.
- Structured-column pair input (CHR/POS/A1/A2 and/or SNP per endpoint), keyed
  through the package's existing identity system in the panel's mode.
- Allele-order/strand-agnostic **matching** (inherited from the identity system)
  and **sign harmonization to the query's allele orientation** in allele-aware
  modes (base modes return `NA` sign; see §9).
- Auto-selected lookup strategy (random-access row-group pruning vs. streaming).
- A `query-r2` CLI subcommand.
- Validation, error handling, logging, and tests.

### Out of scope (deferrable without rework — see §13)
- Other query shapes (SNP-set → pairwise matrix; focal-SNP → all partners).
- Single pre-built key-string input (only structured columns are supported).
- Reading legacy 10-column / external raw-R² panels for *query* (the binding and
  schema contracts here target package-written index panels; raw-R² conversion
  math is still supported by reading `ldsc:r2_bias`).
- Returning raw (biased) R² as an output column (the converter recovers it
  internally; not surfaced).

## 3. Design principles & rationale

1. **Reuse, don't reinvent.** The parquet adapter already has the helpers this
   tool needs: directory→path resolution (`_resolve_r2_build_dir`,
   `_r2_dir_r2_paths`, `_r2_dir_metadata_paths`), schema-metadata parsing
   (`_read_r2_schema_meta`, `_read_identity_schema_meta`), binding hash
   (`sidecar_identity_sha256`), and identity keys (`build_snp_id_series` /
   `effective_merge_key_series`). The query tool composes these rather than
   duplicating them.

2. **Random access is a different access pattern from the LD-score sweep.** The
   streaming reader is built around the *retained-metadata remap* for whole-file
   accumulation. Pair query is point lookup. Keeping them in separate modules
   preserves clean boundaries (see §4); the streaming reader is untouched.

3. **Matching is governed by the package identity system, not new logic.** In
   allele-aware modes the identity key already normalizes allele *order* and
   *strand* (and rejects strand-ambiguous SNPs). In base modes alleles are not
   part of identity. So allele-order-agnostic matching is automatic.

4. **R² magnitude is orientation-independent; only the sign is.** Swapping A1/A2
   maps dosage `x → 2−x`, negating `r` but leaving `r²` unchanged. So allele
   orientation affects only the returned sign / signed `r`, never `r2`.

## 4. Module layout & public surface

```
src/ldsc/_kernel/r2_query.py   # lookup kernel (compute, no user-facing API)
src/ldsc/r2_query.py           # public: R2Panel + query_r2 + unbiased_r2_to_pearson_r
src/ldsc/cli.py                # + `query-r2` subcommand dispatch
src/ldsc/__init__.py           # export R2Panel, query_r2, unbiased_r2_to_pearson_r
```

Public objects (re-exported from `ldsc`):

- `R2Panel` — reusable panel handle (class).
- `query_r2(pairs, panel_dir=..., ...)` — one-shot wrapper (opens a handle,
  queries, discards).
- `unbiased_r2_to_pearson_r(r2_adj, n, sign=None)` — pure converter.

## 5. `R2Panel` handle

```python
@classmethod
def open(
    cls,
    panel_dir: str | Path | None = None,
    *,
    meta_path: str | Path | None = None,
    parquet_path: str | Path | None = None,
    snp_identifier: str | None = None,   # default: ldsc:snp_identifier from parquet schema
    genome_build: str | None = None,     # for hg19/hg38 sub-dir resolution
    validate_binding: bool = True,
) -> "R2Panel": ...
```

### 5.1 Two input modes (mutually exclusive)
- **`panel_dir`** — a `build-ref-panel` output directory (optionally with
  `hg19/`/`hg38/` sub-dirs, resolved by `_resolve_r2_build_dir` using
  `genome_build`). The handle discovers every `chr{N}_r2.parquet` and its
  `chr{N}_meta.tsv.gz` sidecar via the existing helpers → `{chrom: (meta, parquet)}`.
- **`meta_path` + `parquet_path`** — one explicit chromosome pair (chromosome
  taken from the sidecar). Both required together.

Supplying neither, or `panel_dir` together with `meta_path`/`parquet_path`, is a
`LDSCUsageError`.

### 5.2 Lazy per-chromosome loading
The handle resolves the path map at `open()` but loads chromosome state **on
first touch** (queries usually hit one or a few chromosomes; eagerly reading 22
sidecars + hashing them would be wasteful). Per chromosome, first touch:
1. Reads the sidecar (`CHR, POS, SNP, A1, A2, CM, MAF`), preserving build/row order.
2. Parses parquet schema metadata: `n_samples`, `r2_bias`, `r2_scale`,
   `r2_encoding`, `n_snps`, `sorted_by_build`, `sidecar_identity_sha256`,
   `snp_identifier` (provenance).
3. If `validate_binding`: recomputes `sidecar_identity_sha256` and checks it +
   `len(sidecar) == n_snps`; mismatch is a hard `LDSCInputError` (reused message).
4. Builds the identity-key → row-index map under the active `snp_identifier`
   (`build_snp_id_series`). Rejects key collapse (two sidecar rows sharing a key
   under the active mode) with the same injectivity hard-error the reader uses.
5. Opens the `pyarrow.parquet.ParquetFile` and caches its per-row-group `IDX_1`
   min/max footer statistics (for pruning, §8).

All of the above is cached for the chromosome's lifetime on the handle.

### 5.3 Attributes / accessors
- `panel.chromosomes -> list[str]`
- `panel.snp_identifier -> str`
- `panel.n_samples -> int | None`  (consistent across chromosomes; first
  chromosome wins, with a warning if a later chromosome disagrees)
- `panel.query_pairs(pairs, *, with_r=False, strategy="auto", strategy_threshold=50_000) -> pd.DataFrame`

## 6. Pair input format

`pairs` is a `pandas.DataFrame` (the CLI reads a TSV/CSV into one). Each row is
one queried pair; the two endpoints are columns suffixed **`_1`** and **`_2`**:

| Endpoint 1 | Endpoint 2 | Required when |
|---|---|---|
| `CHR_1, POS_1` | `CHR_2, POS_2` | `chr_pos*` modes |
| `A1_1, A2_1` | `A1_2, A2_2` | allele-aware modes only (matching + sign) |
| `SNP_1` | `SNP_2` | `rsid*` modes |

Resolution splits the frame into two endpoint tables, renames `*_1`/`*_2` to the
canonical `CHR/POS/SNP/A1/A2`, and calls `build_snp_id_series` on each in the
panel's mode. Column-name resolution reuses the package's inference aliases where
practical; the canonical suffixed names above are the documented contract.

**Base modes ignore alleles entirely.** In `rsid` / `chr_pos` modes the tool
consumes only the matching columns (`SNP_*` or `CHR_*/POS_*`); any `A1_*`/`A2_*`
columns present in the input are **dropped before resolution and never read**.
A base-mode query therefore behaves identically whether or not allele columns
are supplied — same matching, same `r2`, and `sign = NA` either way (§9). Allele
columns are read only in allele-aware modes.

## 7. Identity resolution, chromosome routing & status

For each pair, each endpoint resolves to `(chrom, row_index)`:
- **Routing.** If a `CHR` column is present (always in `chr_pos*` modes), route
  directly to that chromosome's key map. For `rsid*` panels lacking a `CHR`
  column, fall back to a cross-chromosome `rsid → (chrom, row)` map built lazily
  by merging per-chromosome maps.
- **Outcomes** (drives `r2` and `status`):

| Case | `r2` | `status` |
|---|---|---|
| both endpoints → **same row** (diagonal) | `1.0` | `""` |
| both → same chrom, different rows, **pair stored** | dequantized R² | `""` |
| both → same chrom, different rows, **not stored** (out of window / below `min_r2`) | `NaN` | `absent` |
| both resolve, **different chromosomes** | `NaN` | `cross_chromosome` |
| either endpoint **unresolved** (key not in sidecar, incl. invalid/ambiguous query alleles in allele-aware mode) | `NaN` | `not_in_panel` |

`status` is **always present** and is empty (`""`) for every row with a valid
numeric `r2` (including the diagonal). It carries the cause only for `NaN` rows.
The status vocabulary is the closed set `{not_in_panel, cross_chromosome, absent}`.

The `sign` column is independent of `status` and is governed solely by the
identifier mode (§9): determinate `±1` in allele-aware modes, always `NA` in base
modes. `status` never reports a sign-related cause.

## 8. Lookup kernel (`_kernel/r2_query.py`)

Same-chromosome, non-diagonal pairs are grouped by chromosome. Each group is
resolved against that chromosome's parquet using one shared primitive and a
strategy switch.

### 8.1 Shared int64-key match
Canonicalize each query pair to `(i, j) = (min, max)` of the two row indices.
Encode as `key = i * n_snps + j` (fits `int64`, since `n_snps < 2³¹` ⇒
`key < 2⁶²`). Build the sorted unique array of query keys once. For each decoded
row group, form the same key from its `IDX_1, IDX_2`, locate matches with
`np.searchsorted` (+ equality check), and scatter the dequantized
`r2 = R2_int16 / r2_scale` (`float32`) and `SIGN` back to the matching query
rows. A legacy `float32` R² column (no `r2_scale`) is read unscaled.

### 8.2 Two strategies over the primitive
- **random-access (small):** prune row groups using cached `IDX_1` footer
  min/max against the set of needed left indices `{i}`; run the match only on
  surviving groups. Low I/O for small/interactive queries.
- **streaming (large):** scan every row group once in `IDX_1` order; run the
  match on each. Bounded memory (one row group at a time), no pruning overhead.

### 8.3 Auto-select
`strategy="auto"` (default): **random-access** when the group's
`n_query_pairs ≤ strategy_threshold` (default **50 000**), else **streaming**.
Overridable via `strategy ∈ {"auto","random","stream"}` and `strategy_threshold`.
A parity test asserts `random` and `stream` produce identical results.

## 9. Sign (allele-aware modes only)

**Sign is produced only in allele-aware modes** (`rsid_allele_aware`,
`chr_pos_allele_aware`). In base modes (`rsid`, `chr_pos`) **no allele
information is consulted at all** — not for matching, not for orientation — so
`sign` is **always `NA`** (direction unknown), while `r2` is returned normally by
position/rsID matching. Any A1/A2 columns supplied in a base-mode query are
ignored entirely (so the row-B "no alleles" and row-C "mismatched alleles" cases
from brainstorming behave identically).

In allele-aware modes the returned sign is harmonized to the **query's** allele
coding. The stored `SIGN` is the sign of Pearson `r` in the panel's A1/A2
orientation; we re-orient it to the query's:

- For each endpoint, classify query `(A1,A2)` against panel `(A1,A2)`, allowing
  strand complement:
  - **aligned** iff query `A1 ∈ {panel A1, complement(panel A1)}`,
  - **swapped** iff query `A1 ∈ {panel A2, complement(panel A2)}`.
  Exactly one holds (strand-ambiguous pairs are already excluded from matching).
- Pair sign multiplier is `(−1)^(number of swapped endpoints)` — parity over the
  two endpoints (two swaps cancel; a strand flip *without* order swap does not
  flip).

Worked cases (allele-aware mode; panel SNP1 `A/G`, SNP2 `C/T`, stored `r = +0.8`):

| query SNP1 | query SNP2 | swapped count | harmonized r |
|---|---|---|---|
| `A/G` | `C/T` | 0 | `+0.8` |
| `G/A` | `C/T` | 1 | `−0.8` |
| `G/A` | `T/C` | 2 | `+0.8` |
| `T/C` (strand flip) | `C/T` | 0 | `+0.8` |
| `C/T` (flip+swap) | `C/T` | 1 | `−0.8` |

**No `NA` sign on a valid `r2` in allele-aware modes.** Query alleles are
mandatory there and a matched pair's alleles correspond to the panel by
construction, so every numeric-`r2` pair gets a determinate `±1`. Therefore an
`NA` sign on a valid `r2` has exactly one meaning — **base mode** — which is
fully documented behavior and needs no extra status code or diagnostic column.

## 10. Output schema

`query_pairs` / `query_r2` return a `DataFrame` that **echoes all input columns**
in original row order, plus:

| Column | Dtype | Meaning |
|---|---|---|
| `r2` | `float32` | adjusted (unbiased) R² as stored, dequantized; `NaN` per §7 |
| `sign` | `Int8` | allele-aware modes: harmonized `+1`/`−1` (`NA` only where `r2` is `NaN`). Base modes: always `NA` (§9) |
| `status` | `string` | `""` for valid `r2`, else the cause (`not_in_panel` / `cross_chromosome` / `absent`) |
| `r` | `float64` | **only when `with_r=True`** — signed Pearson `r` (§11); `NaN` where `r2` is `NaN` or `sign` is `NA` |

## 11. R²→r conversion

```python
def unbiased_r2_to_pearson_r(r2_adj, n, sign=None):
    """Convert adjusted (unbiased) R² to signed Pearson r.

    Inverts the unbiased correction to recover raw R², then takes the root and
    applies the sign. Vectorized; ``n`` may be scalar or array-like.
    """
    r2_raw = (r2_adj * (n - 2) + 1) / (n - 1)   # invert r2_adj = r2_raw - (1 - r2_raw)/(n - 2)
    r2_raw = np.clip(r2_raw, 0.0, 1.0)          # guard rounding / tiny negatives
    r = np.sqrt(r2_raw)
    return r * sign if sign is not None else r
```

- Pure NumPy; no I/O. `sign=None` returns magnitude `|r|`.
- The handle path `panel.query_pairs(pairs, with_r=True)` (and CLI `--with-r`)
  feeds `panel.n_samples` and the harmonized `sign` column to add the `r` column.
- If `n_samples` is absent (legacy panel without `ldsc:n_samples`) and `with_r`
  is requested, raise an actionable `LDSCInputError`.
- In a **base mode** the `sign` column is all `NA` (§9), so `with_r` yields an
  all-`NaN` `r` column. This is logged as a warning pointing users to an
  allele-aware panel/mode for signed `r`; it is not an error (the magnitude is in
  `r2`).

Inversion derivation: from `r2_adj = r2_raw − (1 − r2_raw)/(n − 2)`,
`r2_raw = (r2_adj·(n − 2) + 1)/(n − 1)`.

## 12. CLI — `ldsc query-r2`

```
ldsc query-r2 (--panel-dir DIR | --meta M --parquet P)
              --pairs FILE [--out FILE]
              [--snp-identifier MODE] [--genome-build {hg19,hg38}]
              [--with-r]
              [--strategy {auto,random,stream}] [--strategy-threshold N]
```

- `--pairs` is a TSV/CSV with the `_1`/`_2` endpoint columns (`-` = stdin).
- Output is a TSV echoing the input columns + `r2` + `sign` + `status`
  (+ `r` when `--with-r`), to `--out` or stdout.
- Dispatch is added to `cli.py` alongside the existing subcommands; argument
  parsing lives in `r2_query.py` (mirroring how other modules own their parsers).

## 13. Error handling

Actionable errors in the package style (cause + most-likely + fix):
- neither / both input modes supplied (`LDSCUsageError`);
- binding hash or `n_snps` length mismatch (hard `LDSCInputError`, reused);
- key collapse under the active mode (hard error, reused injectivity check);
- allele-aware mode but pairs lack `A1_*`/`A2_*` (`LDSCInputError` from
  `build_snp_id_series`);
- `--with-r` / `with_r=True` without `n_samples` (`LDSCInputError`);
- `panel_dir` with no `chr*_r2.parquet` (reused directory error).

## 14. Testing plan (TDD)

- **Converter:** round-trip `raw → adj → raw`; sign application; negative-adj
  input still yields real `r`; `sign=None` returns magnitude; scalar & array `n`.
- **Key resolution:** each of the four identifier modes resolves endpoints to the
  expected rows; allele-order/strand variants collapse to the same row.
- **Sign:** the five §9 allele-aware cases; base mode → `sign` always `NA`
  (alleles ignored); `with_r` in base mode → all-`NaN` `r` + warning.
- **Integration:** build a tiny panel (reuse builder fixtures), query known
  pairs; assert diagonal `= 1.0`, stored pairs match within quantization
  tolerance (`|err| ≤ 1.5e-5`), `absent`/`cross_chromosome`/`not_in_panel` →
  `NaN` with the right `status`.
- **Strategy parity:** `random` vs `stream` identical; `auto` picks the expected
  branch around the threshold.
- **CLI smoke:** end-to-end TSV in → TSV out, with and without `--with-r`, and
  both input modes.

## 15. File-by-file change list

- **new** `src/ldsc/_kernel/r2_query.py` — int64-key match primitive, the two
  strategies + auto-select, row-group pruning helpers.
- **new** `src/ldsc/r2_query.py` — `R2Panel`, `query_r2`, the converter, the
  argument parser + `run_query_r2_from_args` / `main`.
- **edit** `src/ldsc/cli.py` — register and dispatch the `query-r2` subcommand.
- **edit** `src/ldsc/__init__.py` — export `R2Panel`, `query_r2`,
  `unbiased_r2_to_pearson_r`.
- **new** `tests/test_r2_query.py` (+ converter unit tests) — §14.
- **edit/new docs** — a `docs/current/` page for the query tool; link from
  README / tutorials as appropriate.

## 16. Reused helpers (no duplication)

`_resolve_r2_build_dir`, `_r2_dir_r2_paths`, `_r2_dir_metadata_paths`,
`_read_r2_schema_meta`, `_read_identity_schema_meta`,
`_r2_path_has_ldsc_package_schema`, `sidecar_identity_sha256`,
`build_snp_id_series` / `effective_merge_key_series`, `normalize_chromosome`,
and the sidecar reader used by `ParquetR2RefPanel.load_metadata`.

## 17. Open questions / future work

- Optional pre-built key-string input (Option 1 from brainstorming) could be
  added later behind the same API without breaking structured-column callers.
- Additional query shapes (SNP-set matrix, focal-SNP neighbors) can build on the
  same kernel.
