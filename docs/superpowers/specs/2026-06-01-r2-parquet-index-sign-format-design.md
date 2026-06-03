# R2 Parquet Index+Sign Format — Design

Status: approved design (pre-implementation)
Date: 2026-06-01
Scope owner: build-ref-panel / ldscore parquet backend

## 1. Motivation

The package's pairwise-R² reference panels are stored as per-chromosome parquet
files. The current canonical schema carries ten columns — `CHR`, `POS_1`,
`POS_2`, `SNP_1`, `SNP_2`, `A1_1`, `A2_1`, `A1_2`, `A2_2`, `R2` — even though the
reader keys on only a subset depending on the SNP identifier mode. Measured on a
production 1000G 30× chr22 panel (`chr_pos` mode, 521 M pairs, 3.87 GB):

| Column group | Compressed | Share |
|---|---|---|
| `R2` (float32) | 2386 MB | 61.9% |
| `SNP_2` (rsID) | 684 MB | 17.7% |
| `POS_2` (int64) | 596 MB | 15.5% |
| `A1_2` + `A2_2` | 175 MB | 4.5% |
| everything else | small (RLE runs) | <0.5% |

The endpoint identity columns (`SNP_*`, `A*_*`, and most of `POS_*`) are
redundant: the builder *computes* each pair as a pair of integer indices into
the retained reference-SNP table, then **expands** those indices into
POS/rsID/allele strings before writing
([`_standard_r2_arrow_table`](../../../src/ldsc/_kernel/ref_panel_builder.py)).
The reader does the inverse — it reads those strings and **collapses** them back
to indices during decode
([`_decode_canonical_row_group`](../../../src/ldsc/_kernel/ldscore.py)). The
expansion and the collapse are inverse operations wrapped around a value the
builder already holds.

This design removes both. The parquet stores the indices directly, plus a `SIGN`
bit. The per-SNP metadata sidecar — already the authoritative per-SNP table for
the panel — defines the index space.

### Goals

1. Shrink the identity overhead by storing endpoint **indices** instead of
   expanded identifiers.
2. Persist the correlation **sign** (currently computed then discarded), useful
   to downstream consumers even though LD-score regression ignores it.
3. Make the parquet **identifier-mode-agnostic**: one file serves `rsid`,
   `chr_pos`, `rsid_allele_aware`, and `chr_pos_allele_aware`.
4. Simplify and speed up the read hot path (per-pair decode becomes a pure
   integer gather; mode-specific matching happens once per chromosome).

### Non-goals (explicitly deferred / out of scope)

- **External R² ingest.** This package reads and writes only its own canonical
  index format. Foreign formats (TOP-LD pairwise `tsv.gz`, SBayesR/GCTB binary
  LD matrices, PLINK `.ld`, LDstore `.bcor`, …) are not supported. A future
  converter (foreign pairwise table → canonical index parquet + sidecar) is a
  separate spec. The current inline raw-external read path is removed.
- **R² value-level size reduction.** `BYTE_STREAM_SPLIT` encoding and R²
  quantization are out of scope; the `R2` column is written exactly as today.
  They compose with this format and can be added later.
- **Panel migration.** Clean break: existing ten-column panels are not read or
  migrated. Users regenerate panels with `build-ref-panel`.

## 2. The index format

### 2.1 Parquet columns (exactly four)

| Column | Type | Meaning |
|---|---|---|
| `IDX_1` | `int32` | left endpoint: row number in the panel sidecar (build/position order) |
| `IDX_2` | `int32` | right endpoint: row number in the panel sidecar |
| `R2` | `float32` | unbiased R² (unchanged semantics; `min_r2` floor still honored upstream) |
| `SIGN` | `bool` (bit-packed) | `True` ⇔ Pearson r ≥ 0 in the sidecar's A1/A2 orientation |

`CHR`, `POS_*`, `SNP_*`, and `A*_*` are gone. One chromosome per file; the
chromosome is defined by the paired sidecar.

`SIGN` encodes the sign of the Pearson correlation r between the two SNPs'
allele dosages, in the panel's A1/A2 orientation — *not* the sign of the stored
unbiased R² (which can be negative even when r > 0). It is the existing
`corr >= 0` computation in
[`_finalize_pair_columns`](../../../src/ldsc/_kernel/ref_panel_builder.py),
preserved instead of discarded. The sign is interpretable only together with the
sidecar's A1/A2 alleles. Arrow `bool` is bit-packed (1 bit/value); sign is
locally correlated within LD blocks, so it compresses well, but even
uncorrelated it costs ≈ 1 bit/pair.

### 2.2 Invariants

- `IDX_1 < IDX_2` for every row (canonical orientation).
- Rows sorted by non-decreasing `IDX_1`. Because sidecar row order is
  position-sorted, `IDX_1` order ≡ `POS_1` order, so row-group pruning works on
  `IDX_1` footer statistics exactly as it did on `POS_1`.
- One chromosome per file.
- `int32` is sufficient: asserted at write time (`n_snps < 2³¹`; per-chromosome
  SNP counts are far below this).
- All pairs lie within the configured LD window at write time.

The writer asserts the sort invariant per batch on `IDX_1` (replacing the
current `POS_1` assertion), because row-group pruning depends on monotonic footer
statistics.

### 2.3 Parquet schema metadata

Retained: `ldsc:artifact_type` (`ref_panel_r2`), `ldsc:genome_build`,
`ldsc:sorted_by_build`, `ldsc:n_samples`, `ldsc:r2_bias`, `ldsc:min_r2`,
`ldsc:row_group_size`.

Changed / added:

| Key | Value | Purpose |
|---|---|---|
| `ldsc:schema_version` | `2` | bumped; marks the index layout |
| `ldsc:n_snps` | string(int) | full sidecar length = size of the index space |
| `ldsc:sidecar_identity_sha256` | hex string | binding hash (see §4) |
| `ldsc:snp_identifier` | mode string | **provenance only** — no longer gates how the file is read |

`ldsc:snp_identifier` records the mode the panel was built for, but the index
parquet is mode-agnostic at read time (§3.3). All values remain UTF-8 byte
strings.

## 3. Read pipeline

### 3.1 The sidecar is mandatory and authoritative

The sidecar (`chrN_meta.tsv.gz`, columns `CHR, POS, SNP, A1, A2, CM, MAF`) is
unchanged on disk except its `# ldsc:*` header bumps `schema_version` to `2`.
There is **no `IDX` column** — the index binding is implicit row order.

Because indices are meaningless without the exact sidecar, the
synthesize-metadata-from-parquet-endpoints fallback
([`ref_panel.py` `load_metadata`](../../../src/ldsc/_kernel/ref_panel.py)) is
**removed** for the index format. A missing sidecar is a hard error.

### 3.2 Reader changes (Option A: reader loads the full sidecar)

The reader receives the **retained** metadata (the matrix universe and order)
exactly as today, and builds `index_map: identity → retained_idx` from it
unchanged. New steps, all additive and confined to the reader:

1. **Load the full panel sidecar** (`chrN_meta.tsv.gz`, all `n_snps` rows in
   build order). This is the index space the parquet `IDX` values reference.
   *Note: the LD-score path already reads this sidecar inside
   `ref_panel.load_metadata`, but does not retain the pre-restriction frame, so
   Option A re-reads it once per chromosome — a ~10 MB-gz read, sub-second,
   negligible against a multi-GB parquet. Keeping the reader self-contained is
   the explicit trade for not threading state through the pipeline.*
2. **Validate the binding** (§4): recompute the identity hash from the full
   sidecar; hard-error on mismatch, on `len(sidecar) != n_snps`, or on any
   observed `IDX ≥ n_snps`.
3. **Build the remap** `remap: build_idx → retained_idx` (int32 array of length
   `n_snps`, default −1): for each full-sidecar row `b`, set
   `remap[b] = index_map.get(identity_key(full_sidecar_row_b, mode), -1)`, using
   the same mode-dependent `identity_key` that built `index_map`. This is the
   single point where mode-specific matching happens — once per chromosome.
4. **Build `self.build_idx`** (build index of each retained SNP, monotonic
   ascending) for row-group pruning: pruning maps a retained-index window
   `[start, stop)` to the build-index range
   `[build_idx[start], build_idx[stop-1]]` and selects row groups whose `IDX_1`
   footer bounds overlap it.

Per-pair decode collapses to a pure gather:

```
i = remap[IDX_1]; j = remap[IDX_2]; keep = (i >= 0) & (j >= 0)
```

All per-mode branching, `to_pandas`, and per-row key/pos lookups are removed from
the decode path. R² bias correction (`_transform_r2`) still applies to the `R2`
column; `SIGN` is read through unchanged.

### 3.3 Mode-agnosticism

The parquet carries no identity, only indices into the sidecar. The identifier
mode affects only how `index_map`/`remap` keys are formed at reader init
(`rsid` → `SNP`; `chr_pos` → `POS`; allele-aware → effective `CHR:POS:A1:A2`
key). The same index parquet therefore serves all four modes; the analysis mode
comes from the run config, not the file. The retained set still depends on mode
(allele-aware matching can retain/drop differently), but that is decided entirely
in the sidecar-keyed `index_map`/`remap`, never in the parquet.

## 4. Index↔sidecar binding (hash)

The binding breaks if the panel's SNP **set or order** changes, so the hash
covers exactly that.

- **Basis:** newline-joined `f"{CHR}:{POS}:{A1}:{A2}"` over sidecar rows in order,
  computed **in-memory** at build and at load (never written as a sidecar
  column). `CHR:POS:A1:A2` is the stable physical identity and disambiguates
  multi-allelic sites; the rsID is deliberately excluded (it is the most
  volatile identifier across dbSNP builds and changing it does not change the
  index binding). CM/MAF are excluded so a legitimate MAF/CM recompute on the
  same SNP set does not trigger a false mismatch.
- **Algorithm:** SHA-256 (`hashlib`, one pass; sub-second for ~1 M SNPs).
- **Storage:** hex digest in `ldsc:sidecar_identity_sha256`, length in
  `ldsc:n_snps` (parquet metadata only).
- **Validation (load):** recompute from the full sidecar; hard-error on any of:
  digest mismatch (wrong or reordered sidecar), `len(sidecar) != n_snps`
  (truncated/extended), or `IDX ≥ n_snps` (range violation).

## 5. Write pipeline

The builder already yields each pair as `{"i", "j", "R2", "sign"}` from
[`yield_pairwise_r2_rows`](../../../src/ldsc/_kernel/ref_panel_builder.py); today
`"sign"` is computed and dropped, and `i`/`j` are expanded to identifier strings.

Changes:

- The writer emits `IDX_1 = i`, `IDX_2 = j`, `R2`, `SIGN = sign` directly,
  dropping the `reference_snp_table` expansion entirely. `sign` is plumbed
  through to the writer instead of discarded.
- Schema becomes the four columns of §2.1; metadata gains `ldsc:n_snps` and
  `ldsc:sidecar_identity_sha256`, and `ldsc:schema_version` becomes `2`.
- The sort assertion validates `IDX_1` monotonicity.
- The sidecar writer is unchanged except the `# ldsc:*` header `schema_version`
  bump. The builder computes the identity hash from the same retained metadata
  that produces the sidecar, guaranteeing the parquet's recorded hash matches the
  shipped sidecar.

The builder's pair index `i` already equals row `i` of both `reference_snp_table`
and the sidecar (all derived from the same retained `metadata` in the same
position-sorted order — see
[`ref_panel_builder.py`](../../../src/ldsc/ref_panel_builder.py) where
`reference_snp_table`, the pair stream, and `runtime_metadata` share one
`metadata`). This is the load-bearing build-side invariant; the spec relies on
it and the implementation must assert it.

## 6. Correctness invariant: preserve the current pipeline order

The retained matrix universe and its row order must be **bit-identical** to the
current implementation. The refactor changes only the reader's internal endpoint
resolution; it must not reorder or alter any filtering, restriction, matching, or
merging step. The current order (LD-score-from-parquet path) is:

1. `ref_panel.load_metadata(chrom)`
   ([`ref_panel.py`](../../../src/ldsc/_kernel/ref_panel.py)):
   read sidecar(s) → concat → `_apply_snp_restriction`
   (`ref_panel_snps_file` / HM3) → `_apply_maf_filter` → external identity
   cleanup (external inputs only) → `_validate_metadata` (reset index, validate
   uniqueness) → A′.
2. `_align_annotation_bundle_to_ref_panel(annotation_bundle, ref_panel, chrom)`
   ([`ldscore_calculator.py`](../../../src/ldsc/ldscore_calculator.py)):
   intersect annotation bundle B with A′ → retained universe B ∩ A′.
3. `compute_chrom_from_parquet`
   ([`ldscore.py`](../../../src/ldsc/_kernel/ldscore.py)):
   `merge_frequency_metadata` → `apply_maf_filter` →
   `validate_retained_identifier_uniqueness` → `build_window_coordinates` /
   `get_block_lefts` / `check_whole_chromosome_window` →
   `regression_mask_from_keys` → construct `SortedR2BlockReader(metadata=…)` →
   `ld_score_var_blocks_from_r2_reader`.

The matrix order is the row order of the final `metadata` handed to the reader at
step 3. **None of steps 1–3 change.** The reader still builds `index_map` from
that same retained metadata.

**Parity equivalence.** For every build SNP `b`, the new decode yields
`remap[b] = index_map.get(identity_key(full_sidecar_row_b))`. In the current
implementation, decode reads endpoint `b`'s identity from the parquet's per-pair
columns — which were exactly the expansion of build-sidecar row `b` — and looks
it up in the same `index_map`. The identities are therefore the same value, so
the resulting `(i, j)` and the constructed dense matrices are identical. Storing
the index and remapping once per chromosome reproduces, rather than reinterprets,
the existing per-pair mapping.

## 7. Removed surface (clean break)

- The ten-column canonical writer and reader path.
- The inline raw external-input read path (`_query_union_rows_raw` and the raw
  `Dataset` branch in `SortedR2BlockReader`).
- The synthesize-metadata-from-parquet fallback (incompatible with indices).
- Tier-3 build inference from `POS_1` (no positions in the parquet): the index
  format always carries `ldsc:sorted_by_build`; its absence is a hard error.

`_transform_r2` raw-bias correction is retained (cheap; dormant for unbiased
package panels, available for a future raw-R² index panel).

Encountering a non-index schema raises a clear, actionable error directing the
user to regenerate with `build-ref-panel`.

## 8. Error handling

Hard errors with actionable messages for: missing sidecar; identity-hash
mismatch; `n_snps` length mismatch; `IDX` out of `[0, n_snps)`; missing
`ldsc:sorted_by_build`; build mismatch versus the analysis build; and any
non-index / legacy / foreign schema.

## 9. Testing strategy

Test-driven, parity-first.

- **Primary gate — build→read LD-score parity.** On a real fixture, regenerate a
  panel in the index format and assert LD scores are identical to the current
  implementation's output, across **all four identifier modes** (demonstrating
  mode-agnosticism from one file).
- **Writer.** Emits exactly the four columns with correct dtypes; metadata
  carries `schema_version=2`, `n_snps`, and a hash matching the sidecar; `IDX_1`
  sort assertion fires on out-of-order input; `int32` range assertion.
- **Reader remap.** Endpoints absent from the retained set map to −1 and are
  dropped; retained endpoints map to the correct matrix indices; row-group
  pruning on `IDX_1` selects the correct groups.
- **Binding.** Hash mismatch, `n_snps` mismatch, and `IDX` range violation each
  raise; a correct pair validates silently.
- **SIGN.** Round-trips bit-exactly; `True ⇔ r ≥ 0` matches the builder's
  `corr >= 0`.
- **Removed paths.** Feeding a legacy/foreign schema raises the regenerate
  error; a missing sidecar raises.
- Old ten-column fixtures are removed (clean break). The
  `tests/test_ldscore_workflow.py` integration suite remains the parity anchor,
  re-baselined on regenerated index-format fixtures.

## 10. Affected modules (summary)

| Module | Change |
|---|---|
| `_kernel/ref_panel_builder.py` | new 4-col schema + metadata; plumb `sign`; drop identifier expansion; `IDX_1` sort assert; compute identity hash |
| `ref_panel_builder.py` (workflow) | wire hash/`n_snps`; assert pair-index ≡ sidecar-row invariant |
| `_kernel/ldscore.py` (`SortedR2BlockReader`) | load full sidecar; validate binding; build `remap` + `build_idx`; gather-based decode; `IDX_1` pruning; remove raw + legacy-canonical paths |
| `_kernel/ref_panel.py` | remove synth-from-parquet fallback; sidecar mandatory |
| `docs/current/parquet-r2-format-and-read-pipeline.md` | rewrite format/read sections for the index layout |
| `tests/` | new fixtures; parity gate; binding/error tests; remove old-format fixtures |

## 11. Open follow-ups (separate specs)

- External→index converter (restores foreign pairwise ingest as an explicit
  offline step).
- R² column encoding (`BYTE_STREAM_SPLIT`) and quantization for further size
  reduction; both compose with this format.
