# build-ref-panel genotype-reader memory refactor — design

Date: 2026-06-02
Status: approved (design); implementation plan to follow.

## 1. Problem and goal

`build-ref-panel` peaks at ~680 MiB on a restricted chr22 build and is projected
to ~940 MB+ on an unrestricted chr1 build. A memory profile attributed the peak
to one event: the **whole-chromosome PLINK `.bed` payload loaded into a single
in-RAM bitarray before any SNP/MAF filtering** (`PlinkBEDFile.__read__`). The
pairwise R2 computation itself never exceeds ~370 MiB; the read is the cost.

**Goal.** Stop materializing the whole-chromosome genotype payload, so that:

- Restricted builds (HM3 / `--ref-panel-snps-file`) read only the kept SNP
  blocks (peak ~680 → ~370 MiB).
- Unrestricted builds (the default, keep-all) stream the genotypes so only a
  sliding window is ever resident (chr1 ~940 MB → ~one window + workflow floor).

This optimization is limited to PLINK genotype payload reading. The SNP
restriction file itself is not streamed by chromosome: `--ref-panel-snps-file`
or HM3 is resolved up front into an in-memory identity key set, and each
chromosome's PLINK metadata is matched against that set. Very large custom
restriction files therefore still contribute O(number of restriction SNPs)
memory before selective genotype reading begins.

The same reader is shared with the legacy LD-score-with-PLINK path
(`compute_chrom_from_plink`) and the `PlinkRefPanel` API, so the change must be
made at the shared reader **without altering LD-score output and without new
user-facing flags**.

### Non-goals

- Output-artifact encoding (int16 quantize / `BYTE_STREAM_SPLIT` /
  `DELTA_BINARY_PACKED` / zstd) — untouched.
- Genetic-map loading, `.bim` parsing, columnar pair emission — separate efforts.
- Changing the LD-score math, the CLI surface, or splitting the LD-score math
  out of the genotype-array class (explicitly deferred).

## 2. Background: the shared reader and its consumers

`PlinkBEDFile` (subclass of `__GenotypeArrayInMemory__`) in
`src/ldsc/_kernel/ldscore.py` is the single PLINK genotype reader. Four
consumers depend on it:

| Consumer | Access | Passes | Precision | Filters |
|---|---|---|---|---|
| `build-ref-panel` (`ref_panel_builder.py:802`) | 1 forward pass | float32, no `minorRef` | keep_snps + optional keep_indivs + mafMin | target |
| legacy LD-score-PLINK (`compute_chrom_from_plink`, CLI-wired) | **2 passes** (rewind `_currentSNP=0`) | float64 (must stay bit-for-bit) | keep_snps + keep_indivs + mafMin | must not regress |
| `PlinkRefPanel.build_reader` (`_kernel/ref_panel.py:340`) | reader returned | — | keep_snps/keep_indivs/maf | public API, tested |
| `PlinkRefPanel._load_genotype_metadata` (`_kernel/ref_panel.py:359`) | 0 (MAF metadata only) | — | keep_snps + maf | reads whole BED for `bed.df` only |

Both iterating consumers walk SNPs **strictly forward with a carry-over window**
(`__corSumVarBlocks__` and `yield_pairwise_r2_rows` share this structure). The
LD-score path additionally **rewinds and replays** the genotypes a second time.

### Frozen post-`__init__` contract

Any reader change must preserve these attributes and semantics exactly, since
all consumers depend on them (`ldscore.py:388-397`):

- `geno`, `m`, `n`, `nru`
- `_currentSNP` cursor, with **rewind-to-0-and-replay** support
- `kept_snps`, `freq`, `maf`, `sqrtpq`
- `df` subset to `kept_snps` with an appended `MAF` column; `colnames`
- `nextSNPs(b, minorRef=None, dtype=...)` returning mean-imputed, z-scored
  columns — float64 default (LD-score), float32 for ref-panel

The `minorRef` branch (the only in-iteration reader of `self.freq`,
`ldscore.py:628`) is **dead**: no consumer passes `minorRef`. `freq` must still
be computed for `df`/`MAF`/counts, but need not be hot during iteration.

## 3. Approach

### 3.1 Module extraction (pure move, no logic change)

Move `__GenotypeArrayInMemory__`, `PlinkBEDFile`, and the `bitarray`-absent
fallback verbatim from `ldscore.py` into a new `src/ldsc/_kernel/plink_bed.py`.
The in-class LD-score methods (`ldScoreVarBlocks`, `__corSumVarBlocks__`,
`__l2_unbiased__`, `ldScoreBlockJackknife`) travel with the class — separating
them is deferred.

`ldscore.py` re-imports the classes (`from .plink_bed import
__GenotypeArrayInMemory__, PlinkBEDFile`). Because `get_legacy_ld_module()`
returns `sys.modules[__name__]` (the `ldscore` module), the re-export keeps
`legacy_ld.PlinkBEDFile` resolving and leaves `compute_chrom_from_plink`
unchanged. Coordinate helpers (`getBlockLefts`, `get_block_lefts`,
`block_left_to_right`) stay in `ldscore.py` — they are not part of the reader.

This is its own commit, with the full suite green and zero behavior change.

### 3.2 O1 + O3: lazy, selective read (in-RAM mode — the default)

- `__read__`: open the file, validate the header, store `nru` — **do not load
  the payload**. The reader holds the open handle for its lifetime.
- `__filter_snps_maf__`: for each candidate SNP, seek to its source offset
  (`3 + j * nru / 4` bytes) and read only that SNP's bytes; when `keep_indivs`
  is set, **subset individuals while packing** (this is O3 — the second full
  bitarray never forms); compute `freq`/MAF on the kept individuals using the
  exact current arithmetic and apply the drop; append survivors to the retained
  `self.geno`.
- `nextSNPs` is **unchanged**: it reads from the now-small retained `self.geno`.

**Ordering invariant (bit-for-bit safety).** MAF is computed on kept-individual
genotypes, preserving today's "filter individuals, then compute MAF" semantics,
so the retained genotypes are byte-identical to the current implementation.

For restricted `build-ref-panel`, `compute_chrom_from_plink`, and the MAF-only
path, the whole-chromosome bitarray is never materialized — they all pass
`keep_snps`, so selective read applies automatically with no caller change.

### 3.3 O2: streaming mode (opt-in; `build-ref-panel` unrestricted only)

`PlinkBEDFile.__init__` gains `streaming: bool = False`. Only `build-ref-panel`
ever sets it `True`, and only when no SNP restriction was supplied.

- **Trigger (call site).** In `ReferencePanelBuilder`, a build is unrestricted
  iff no restriction resolved: `config.ref_panel_snps_file is None and not
  config.use_hm3_snps` (equivalently `restriction_values is None` on the
  resolved build state, `ref_panel_builder.py:568`). The reader receives only
  the boolean; the intent logic lives in the workflow.
- **Pass 1 (no retention).** Stream the file SNP-by-SNP and **bit-count**
  alleles (no float decode) with the exact current freq formula; apply the
  monomorphic/MAF drop; build `freq`/`maf`/`sqrtpq`/`kept_snps`/`df`+`MAF`.
  `self.geno` is never packed.
- **`nextSNPs` (streaming branch).** Seek to the source offset of the kept SNP
  at `_currentSNP`, decode `b` SNPs, advance. Identical standardization and
  dtype handling — only the byte source differs. A small `kept → source` index
  array (one int per kept SNP) drives the seeks. `_currentSNP = 0` re-seeks to
  the first kept SNP; rewind is preserved for contract safety even though
  `build-ref-panel` is single-pass.

The LD-score path never opts in, so its in-RAM rewind-replay stays free — no
second disk pass is imposed on it.

### 3.4 Single-class shape

`_streaming` gates exactly two spots: whether `__filter_snps_maf__` packs
`self.geno`, and which source `nextSNPs` reads. Decode/standardization, the
cursor, and the whole post-`__init__` contract are one shared code path. If the
branch ever threatens to sprawl, fall back to a shared decode helper; the target
is two small `if self._streaming` sites.

## 4. Data flow

**Restricted / in-RAM (default):** open handle → per-SNP selective read +
individual subset + MAF drop → pack retained `self.geno` → `nextSNPs` reads
`self.geno` → (LD-score: rewind replays from `self.geno`, free).

**Unrestricted / streaming (`build-ref-panel` only):** open handle → pass-1
bit-count metadata (no pack) → `nextSNPs` seeks the file per batch → single
pass, peak genotype RAM ≈ one window.

## 5. Error handling and edge cases

- **File handle lifetime:** held open for the reader's life; closed via
  context-manager / `__del__` / explicit close. Per-chromosome
  construct-use-discard is unchanged.
- **`bitarray` absent:** the dependency-gated fallback class is preserved.
- **Empty kept set (`m <= 0`):** raises as today.
- **Out-of-range batch request:** raises as today.
- **Monomorphic/all-missing drop:** identical in both modes — pass-1 bit-count
  reuses the current freq arithmetic, so drop decisions match exactly.
- **Seek order / I/O locality:** `keep_snps` arrives sorted by genome-build
  position, which may differ from source `.bed` order under liftover; seeks
  remain correct but may be non-sequential. Locality only, not correctness.

## 6. Testing (TDD guardrails)

Fixtures already present: `tests/fixtures/plink/plink.{bed,bim,fam}` and
`tests/fixtures/minimal_external_resources/plink/hm3_chr22_subset.{bed,bim,fam}`.

- **Pure-move gate:** full suite green after §3.1 with no behavior change.
- **Golden bit-for-bit:** snapshot current-`main` outputs from both fixtures —
  LD-score `.l2.ldscore` / `M` / `M_5_50` (via `compute_chrom_from_plink`) and
  `build-ref-panel` parquet (values + schema) — and assert the refactor
  reproduces them.
- **Mode equivalence:** streaming vs in-RAM produce identical `build-ref-panel`
  parquet for the same unrestricted fixture.
- **Filter matrix:** `keep_snps {None-equivalent full, subset} × keep_indivs
  {None, subset} × mafMin {0, 0.05}`.
- **MAF-on-kept-individuals:** explicit test that MAF uses kept individuals.

Per `CLAUDE.md`, numeric tests assert against reference output, not just
"runs without error."

## 7. Implementation sequence (commits)

1. **Pure move** → `_kernel/plink_bed.py` with re-export; suite green.
2. **O1 + O3** — lazy selective read + fused individual filter (in-RAM mode);
   golden + matrix tests green. Benefits restricted `build-ref-panel`,
   `compute_chrom_from_plink`, and the MAF-only path.
3. **O2** — streaming mode + `build-ref-panel` unrestricted trigger; mode
   equivalence tests green.
4. **Targeted cleanup** — remove the dead `minorRef` / `freq`-in-`nextSNPs`
   branch.

## 8. Risks and mitigations

- **LD-score regression.** Mitigated by O1/O3 not touching `nextSNPs` decode and
  by the golden bit-for-bit gate. O2 is gated away from the LD-score path.
- **Streaming bit-for-bit vs in-RAM.** Mitigated by the mode-equivalence test
  and by sharing the decode/standardize code path between modes.
- **Handle leaks.** Mitigated by context-manager/`__del__` and per-chromosome
  lifetime; covered by existing reader teardown paths.
- **Scope creep into the LD-score-math split.** Explicitly out of scope; the
  move is the reader class only.
