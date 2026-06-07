# Allele Orientation Canonicalization (A1/A2/Frequency) — Design

**Status:** approved design (brainstorming resolved 2026-06-07)
**Scope:** package-wide allele (`A1`/`A2`) and frequency (`MAF`/`FRQ`) semantics for
every non-sumstats artifact; sumstats are an explicit carve-out. Touches
`build-ref-panel`, LD-score artifacts, annotation/frequency sidecars, `query-r2`,
restriction/keep matching, and the documentation of `munge-sumstats`.
**Companion plan:** `docs/superpowers/plans/2026-06-07-allele-orientation-canonicalization-plan.md`

---

## 0. Relationship to the current code

Today the package has **no global allele-orientation convention**:

- `build-ref-panel` copies `A1`/`A2` verbatim from the PLINK `.bim` and folds
  frequency to a minor-allele value: `self.maf = np.minimum(self.freq, 1 - self.freq)`
  ([`_kernel/plink_bed.py:67`](../../../src/ldsc/_kernel/plink_bed.py)), where
  `self.freq` is the frequency of the `.bim` **second** allele (`A2`). The meta
  sidecar therefore carries identity-only `A1`/`A2` and a folded `MAF` that is
  **not tied to `A1`**.
- The R² parquet stores `SIGN = (Pearson r >= 0)` **in the panel's `A1`/`A2`
  orientation** ([`_kernel/ref_panel_builder.py:371,631,718`](../../../src/ldsc/_kernel/ref_panel_builder.py)).
  Because orientation is whatever PLINK happened to emit, the sign of `query-r2
  --with-r` has no stable, documented meaning.
- The per-read `minorRef` re-orientation that legacy LDSC carried on
  `nextSNPs` was production-dead and was removed in `ce341b6`. `nextSNPs` now
  standardizes raw genotypes (coding the `.bim` `A2` allele) with no flip.
- `munge-sumstats` already treats `A1` as the **effect allele** and folds `FRQ`
  *only* for the `--maf-min` mask, never overwriting the column or reorienting
  alleles ([`_kernel/sumstats_munger.py:120,252`](../../../src/ldsc/_kernel/sumstats_munger.py)).

This design introduces **one** convention — `A1` is the minor allele everywhere
except sumstats — applied by canonicalizing artifacts at write time and
defensively normalizing externally-supplied artifacts at read time.

Empirical motivation (chr22 of `resources/example_1kg_30x`, 438,232 SNPs, 3,202
samples): the `.bim` `A2` allele has frequency `> 0.5` for **94.74%** of SNPs and
`<= 0.5` for **5.26%** (5 SNPs exactly at 0.5). So PLINK's usual `A1`=minor /
`A2`=major holds for ~95% but **not universally** — the convention must be
*enforced*, never *assumed*.

## 1. The unifying rule

> **The frequency column always means "frequency of `A1`."**
> The two regimes differ only in what `A1` *is* and what that implies:
>
> - **Outside sumstats:** `A1 ≡ minor allele` ⇒ `MAF = freq(A1)` and is
>   **guaranteed `<= 0.5`**.
> - **In sumstats:** `A1 ≡ effect allele` ⇒ `FRQ = freq(A1)`, **unconstrained**
>   (may be `> 0.5`), never folded, never reoriented.

The PLINK association case (`A1` = minor = effect) is the special point where the
two regimes coincide, which is why PLINK output reads cleanly under either
interpretation. The package keeps them distinct so a GWAS whose effect allele is
*major* still munges correctly (effect sign preserved) instead of being silently
flipped.

### Canonical rule (non-sumstats)

- **`A1` is the minor allele; `A2` is the major allele.**
- Formally the enforced invariant is **`freq(A1) <= 0.5`** (`<=`, not `<`).
- **Tie-break:** when `freq(A2) == 0.5` exactly, keep the PLINK `.bim` order (no
  flip). Both alleles have frequency 0.5, so `freq(A1) = 0.5 <= 0.5` still holds.

## 2. Definitions of `A1`, `A2`, and frequency across artifacts

| Artifact | `A1` | `A2` | Frequency column | Frequency meaning | Always `<= 0.5`? |
|---|---|---|---|---|---|
| **PLINK `.bim` / `.bed`** (input) | `.bim` allele 1 (PLINK *usually* sets this to the minor/counted allele, but not guaranteed) | `.bim` allele 2 (usually major) | — (computed from `.bed`; package's internal `freq` is the `A2` frequency) | n/a (identity only) | n/a |
| **LD reference panel** — meta sidecar `chr*_meta.tsv.gz` (`CHR POS SNP CM MAF A1 A2`) | **minor allele** (canonical) | **major allele** (canonical) | `MAF` | `freq(A1)` = minor-allele frequency | **Yes** (`<= 0.5` by construction) |
| **LD reference panel** — R² parquet `chr*_r2.parquet` (`IDX_1 IDX_2 R2 SIGN`) | (alleles via the meta sidecar rows the indices reference) | — | — | — | n/a |
| **LD-score artifacts** (`baseline.parquet` / `query.parquet`, allele-aware mode) | **minor allele** (inherited from canonical ref panel) | **major allele** (inherited) | `MAF` | `freq(A1)` = minor-allele frequency | **Yes** |
| **Annotation / frequency sidecars** (optional `A1`/`A2`, optional freq) | **minor allele** (canonical on write; normalized on read) | **major allele** | `MAF`/freq column if present | `freq(A1)` = minor-allele frequency | **Yes** |
| **Restriction / keep files** (optional `A1`/`A2`) | matched as an **unordered, strand-aware allele set** — order is not interpreted | same | — | — | n/a |
| **Curated sumstats** (`.sumstats.gz`: `SNP CHR POS A1 A2 Z N`) | **effect allele** (signed-statistic reference) | counterpart allele | `FRQ` only if `--keep-maf`; dropped by default | `freq(A1)` = effect-allele frequency | **No** (unconstrained) |
| **Raw GWAS sumstats** (input to munge) | aliases `EA`/`EFFECT_ALLELE` → `A1` | `NEA`/`NON_EFFECT_ALLELE` → `A2`; `REF`/`ALT` are **not** global aliases | `FRQ`/`MAF` alias | `freq(A1)` (effect-allele frequency in the source) | source-dependent |

Notes:

- In **all** rows the frequency column, when present, is the frequency of `A1`.
  The label `MAF` (reference panel / LD-score / annotation) advertises the extra
  guarantee that `A1` is the minor allele, so `MAF <= 0.5`. The label `FRQ`
  (sumstats) advertises no such guarantee.
- Reference-panel / LD-score / annotation `A1` carries **no effect or reference
  semantics** — it is the minor allele for identity and frequency, nothing more.
- Sumstats `A1` carries **no minor/major or reference semantics** — it is the
  effect allele for statistical direction, nothing more.

### 2.1 How a sidecar signals that its frequency is *oriented*

The defensive-read helper (§7) can only normalize a table whose frequency is the
**oriented frequency of `A1`** (which may exceed 0.5). A **folded** `min(p, 1-p)`
value cannot be normalized because folding already discarded the orientation. A
sidecar declares which kind it provides **by column name**:

| Column name (and aliases) | Meaning | Range | Orientation source? |
|---|---|---|---|
| `FRQ` (`FREQ`, `EAF`, `A1_FRQ`, `FRQ_A1`, `FREQ_A1`) | **Oriented** frequency of the allele in the `A1` column | `[0, 1]` — values `> 0.5` allowed | **Yes** — helper swaps `A1`/`A2` and sets `freq := 1 - freq` for rows with `freq(A1) > 0.5` |
| `MAF` | **Folded** minor-allele frequency (orientation-free) | `[0, 0.5]` | **No** — identity-only; cannot drive `A1`/`A2` normalization |

Requirements for an oriented (`FRQ`-family) frequency sidecar:

- It carries `A1` and `A2` columns alongside the frequency, so the helper can
  swap both alleles and the value together.
- The frequency is the frequency of the allele named in `A1` (not of a reference
  or alternate allele); `REF`/`ALT` frequencies are a different concept and are
  not accepted as `A1` frequencies.
- A sidecar that provides only a `MAF` column (no oriented `FRQ`) is read as
  identity-only: its alleles still participate in allele-aware matching, but its
  frequency is not used to re-orient `A1`/`A2`.

**Registry implication.** Today `MAF`, `FRQ`, `FREQ`, `FREQUENCY` are mutual
aliases for a single `MAF` column ([`column_inference.py:112`](../../../src/ldsc/column_inference.py)).
For the defensive-read path the `FRQ`-family (oriented) names must be separated
from `MAF` (folded) so the helper can tell them apart. Inside sumstats the
existing conflation is fine and stays (everything there is the oriented `A1`
frequency regardless of source name); the split applies to the non-sumstats
sidecar reader.

## 3. Is the frequency always the MAF (`<= 0.5`)?

- **Reference panel, LD-score artifacts, annotation/frequency sidecars:** **Yes.**
  By the canonical rule `A1` is the minor allele, so the reported `MAF = freq(A1)`
  is `<= 0.5` by construction. (Numerically this equals today's folded
  `min(freq, 1 - freq)`; the change is the *guarantee that the value belongs to
  `A1`*, plus the ~5% of rows whose `A1`/`A2` labels are swapped to make it so.)
- **Sumstats (`FRQ`):** **No.** `FRQ = freq(A1)` where `A1` is the effect allele.
  It is *not* folded and may exceed 0.5. The `--maf-min` filter folds a *local
  copy* for the threshold test only and never mutates the stored column or the
  allele orientation.

## 4. Sign convention (signed r in the reference panel)

### Definition

- **`A1` is the minor allele** (canonical reference-panel orientation).
- **Signed r is the Pearson correlation between the two SNPs' `A1` (minor) allele
  dosages.** `SIGN` stores `True` when this r `>= 0`.
- **Positive sign** → the minor alleles of the two SNPs co-occur on haplotypes
  more than chance → **positive LD between minor alleles**.
- **Negative sign** → a SNP's minor allele tends to ride with the *other* SNP's
  major allele.

### Why the sign is well-defined and implementation-robust

For a biallelic SNP, `A1-dosage = 2 - A2-dosage`, and
`corr(2 - x, 2 - y) = corr(x, y)`. So orienting **both** SNPs of a pair by the
minor allele and orienting **both** by the major allele give the **identical**
sign. The sign only flips if the two endpoints are oriented *differently* — which
is exactly the inconsistency canonicalization eliminates.

Consequence for the implementation: it may take the cheap path of negating only
the ~5% of standardized columns whose `.bim` `A2` is the minor allele (so the
*coded* genotype lands on the canonical `A2`/major side for every SNP), and the
resulting `SIGN` is **exactly** the `A1`-referenced (minor-allele) sign. The
allele the implementation happens to *code* (major) and the allele the contract
*documents* (`A1`/minor) differ, but yield the same number.

### Modes without a sign

Base / allele-blind identity modes carry no `A1`/`A2` and therefore no sign;
`query-r2 --with-r` already errors in those modes
([`r2_query.py:391`](../../../src/ldsc/r2_query.py)) and continues to.

## 5. Build-stage canonicalization (Approach A)

The flip is applied in `build-ref-panel`, strictly **between `nextSNPs` and the
correlation / `SIGN` / metadata writes**. `nextSNPs` stays orientation-free.

```
geno.nextSNPs(b)               # raw, A2-coded standardized columns (orientation-FREE)
   │
   ▼  per-SNP sign flip:   flip = (A2_freq < 0.5)        # strict <, ties keep order
   ▼     negate standardized column for flipped SNPs
   ▼     swap A1 <-> A2 in the meta sidecar for the same SNPs
   │
   ▼
np.dot(A.T, B / n)  →  corr  →  SIGN          # CANONICAL by construction
meta sidecar: A1=minor, A2=major, MAF=freq(A1) (<= 0.5)   # CANONICAL
```

### Mechanics

1. **Compute the flip mask** from the reader's per-SNP `A2` frequency
   (`geno.freq`): `flip = freq < 0.5`. `freq == 0.5` keeps PLINK order.
2. **Orient the standardized genotypes.** The windowed correlation pulls blocks
   `A` (left) and `B` (right) via `standardized_snp_getter`
   ([`ref_panel_builder.py:1003`](../../../src/ldsc/ref_panel_builder.py),
   consumed at [`_kernel/ref_panel_builder.py:473,511,524`](../../../src/ldsc/_kernel/ref_panel_builder.py)).
   Wrap the getter so each returned column `k` is multiplied by `sign[k]`
   (`-1` if flipped, else `+1`), tracking the sequential SNP position. Because
   `corr_canonical(i,j) = sign[i]·sign[j]·corr_raw(i,j)`, applying the per-column
   sign before the dot-products yields canonical `corr` and hence canonical
   `SIGN`, with **no change to R²** (`sign²` cancels).
3. **Swap labels in metadata.** Where `flip` is true, swap the `A1` and `A2`
   strings during meta assembly ([`ref_panel_builder.py:800`](../../../src/ldsc/ref_panel_builder.py)).
   `MAF` already equals `min(freq, 1 - freq)`; after the swap it is, by
   construction, `freq(A1)`.

### What this does and does not change

- **Unchanged (sign-invariant):** R², LD scores, partitioned LD scores, `h2`,
  `rg`. These never depend on `A1`/`A2` orientation. The regression hot path
  cannot regress from this work.
- **Changed (observable):** meta `A1`/`A2` order for the ~5% flipped SNPs; the
  guarantee `MAF = freq(A1) <= 0.5`; `SIGN` (now canonical and meaningful);
  `query-r2 --with-r` output sign (now stable and documented).

### `nextSNPs` and other callers

`nextSNPs` is left orientation-free. Its only other caller is the **r²-based**
in-PLINK LD-score path (`ldScoreVarBlocks` / `ldScoreBlockJackknife`,
[`_kernel/plink_bed.py:87,93`](../../../src/ldsc/_kernel/plink_bed.py)), which is
sign-invariant and neither needs nor is harmed by the absence of orientation.
Orientation is a build-policy decision, not a read primitive; pushing it into
`nextSNPs` would re-introduce the hidden coupling that made `minorRef` dead
weight.

## 6. Downstream artifacts

### LD-score artifacts (write)

In allele-aware mode `baseline.parquet` / `query.parquet` carry `A1`/`A2` and
`MAF` sourced from `ref_panel.load_metadata` (`MAF` rides in
`_LDSCORE_SUFFIX_COLUMNS`, [`ldscore_calculator.py:77`](../../../src/ldsc/ldscore_calculator.py)).
A canonical ref panel makes these **canonical for free**. The LD-score writer
adds a cheap assertion of the invariant (`freq(A1) = MAF <= 0.5`) rather than its
own flip logic.

### Annotation / frequency sidecars

- **Write:** sidecars produced by the package that carry `A1`/`A2` (and a
  frequency) are written canonical.
- **Read (defensive):** an externally-supplied sidecar that carries an *oriented*
  `A1` frequency is normalized to canonical via the shared helper (below). A
  sidecar that carries no oriented frequency cannot be re-oriented (folding
  already destroyed the information); such inputs are treated as identity-only
  and must not be used as a frequency source.

### `query-r2`

- `--with-r` is documented as returning **signed r = correlation of `A1` (minor)
  allele dosages**, with `A1` guaranteed minor by the canonical panel
  ([`r2_query.py:376`](../../../src/ldsc/r2_query.py)). No code change beyond
  documentation if it reads a canonical panel; for robustness it relies on the
  same invariant the panel guarantees.

### Restriction / keep files and regression alignment

Matching uses the **unordered, strand-aware allele set**
([`_row_alignment.py:58`](../../../src/ldsc/_row_alignment.py)); `A1`/`A2` order
is never interpreted for identity. Canonicalization is therefore **orientation-
neutral** for matching — a canonical panel and a sumstats file with `A1` = effect
still match correctly, with any order disagreement resolved as a sign flip at
alignment time, not a match failure. No change required; documented as such.

## 7. Defensive-read helper

A single shared function canonicalizes any table that carries `A1`, `A2`, and an
**oriented `A1` frequency**:

```
canonicalize_alleles(df, a1_col="A1", a2_col="A2", freq_col="MAF") -> df
    # freq_col holds the oriented frequency of A1 (caller names the column).
    # For rows where freq(A1) > 0.5: swap A1<->A2 and set freq := 1 - freq.
    # Idempotent: a canonical table is returned unchanged.
    # freq == 0.5 rows are left as-is (tie-break: keep order).
```

- Lives alongside the allele-aware identity utilities in
  [`_kernel/snp_identity.py`](../../../src/ldsc/_kernel/snp_identity.py).
- Used on **externally-supplied** artifacts only; package-produced artifacts are
  already canonical at write, so the helper is a no-op verification for them.
- It **cannot** repair a table that lacks an oriented frequency (e.g. a legacy
  folded `MAF` with arbitrary `A1`/`A2`) — the orientation is unrecoverable. Such
  tables are out of scope for repair (see §8). The oriented vs folded distinction
  is signaled by column name per §2.1: only an `FRQ`-family column drives
  normalization; a `MAF` column is identity-only.

## 8. Legacy artifacts

All non-sumstats artifacts the package previously produced (folded `MAF`,
arbitrary PLINK `A1`/`A2`, `SIGN` in PLINK orientation) will be **regenerated**
with the updated code. There is therefore **no migration machinery, no version
gate, and no compatibility shim** for package-produced artifacts — they are
simply rebuilt canonical. The defensive-read helper exists for *externally*
supplied inputs that carry an oriented frequency, not for repairing legacy
package output.

## 9. Sumstats carve-out (no behavioral change)

`munge-sumstats` already implements the sumstats regime correctly:

- `A1` = effect allele; never reoriented to minor
  ([`sumstats_munger.py:120,373`](../../../src/ldsc/_kernel/sumstats_munger.py)).
- `FRQ` = `freq(A1)` (effect-allele frequency); folded only as a local mask for
  `--maf-min` ([`:252`](../../../src/ldsc/_kernel/sumstats_munger.py)), never
  overwritten; retained as the oriented value under `--keep-maf`, dropped by
  default.
- `A1`/`A2` matching uses the unordered strand-aware set.

The only deliverable here is **documentation**: an explicit note that the
minor-allele canonical rule **must not** be applied to sumstats, so a future
"consistency" pass does not mistakenly fold or reorient `A1`/`A2`/`FRQ`.

## 10. Documentation deliverables

- `docs/current/column-schema.md`: define `A1`/`A2`/`FRQ`/`MAF` per artifact per
  the table in §2; state the unifying rule (§1); state the `MAF <= 0.5` guarantee
  and its sumstats exception (§3); state the sign convention (§4); document the
  oriented-vs-folded naming convention for external frequency sidecars (§2.1).
- `docs/current/data-flow.md` / `architecture.md`: note the build-stage flip
  point and that `nextSNPs` is orientation-free.
- `query-r2` CLI help / docstring: document signed r as the `A1` (minor) allele
  correlation.
- `docs/troubleshooting.md`: if a "non-canonical external frequency sidecar"
  error is added, give it a section.

## 11. Test plan

Tests precede implementation (red → green). Numerical tests assert against known
values, not merely "runs without error."

1. **Flip mask correctness.** Synthetic panel with known `A2` frequencies above,
   below, and exactly at 0.5 → `flip = freq < 0.5`; the 0.5 case keeps order.
2. **Canonical meta.** After build, every row has `MAF = freq(A1) <= 0.5` and
   `A1`/`A2` swapped exactly on the flipped rows; the unordered `{A1, A2}` set is
   unchanged for every SNP.
3. **R² invariance.** R² (and resulting LD scores) are **bit-for-bit identical**
   to a pre-change build on the same input (sign-invariance check).
4. **SIGN canonical.** For a small hand-computable panel, `SIGN` equals the sign
   of the correlation of `A1` (minor) allele dosages; verify a pair that includes
   one flipped and one non-flipped SNP gets the correct sign.
5. **Signed-r equivalence.** `query-r2 --with-r` returns the same signed r whether
   the underlying coding negates the minority or majority of columns (minor- vs
   major-coded equivalence).
6. **Defensive read.** `canonicalize_alleles` flips rows with `freq(A1) > 0.5`,
   leaves `== 0.5` and already-canonical rows unchanged, and is idempotent.
7. **Regression neutrality.** `h2`/`rg` on a canonical panel match results on a
   pre-change panel within numerical tolerance (matching is orientation-neutral).
8. **Sumstats untouched.** Munge output `A1`/`A2`/`FRQ` are unchanged by this
   work; an effect-allele-is-major fixture keeps `FRQ > 0.5` and the original
   `A1`.

## 12. Risks and mitigations

- **Per-column sign bookkeeping in the windowed builder.** The sign vector must
  index SNPs consistently as they appear in both `A` (left) and `B` (right)
  blocks. Mitigation: apply the sign at the single `standardized_snp_getter`
  boundary so it is impossible for a SNP to be oriented one way as a left index
  and another as a right index; test #4 covers a mixed pair.
- **Streaming build path.** `build-ref-panel` has a streaming reader mode; the
  flip must apply identically in streaming and in-RAM modes. Mitigation: the flip
  is a property of the standardized-column getter, which both modes share.
- **Accidental sumstats folding.** Mitigation: §9 documentation + test #8.
