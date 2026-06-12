# CM/MAF Source-of-Truth Across LD-Score Backends Design

**Date:** 2026-06-11
**Implementation status:** not yet implemented (design approved; plan to follow).
**Documentation status:** to be updated on implementation — `docs/current/`,
`tutorials/`, and `docs/troubleshooting.md`.
**Scope:** `src/ldsc/annotation_builder.py`, `src/ldsc/_kernel/ldscore.py`,
`src/ldsc/_kernel/ref_panel.py`, `src/ldsc/ref_panel_builder.py`,
`src/ldsc/ldscore_calculator.py`, `src/ldsc/cli.py`, `src/ldsc/config.py`,
`src/ldsc/outputs.py`, `docs/troubleshooting.md`, and the `tests/` LD-score suite.

---

## Problem

`ldsc ldscore` has two reference-panel backends: a parquet-R2 backend
(`compute_chrom_from_parquet`) and a PLINK backend (`compute_chrom_from_plink`).
They currently draw the genetic-map coordinate `CM` and the minor-allele
frequency `MAF` from **different sources**, so the same annotation inputs can
produce different LD windows and different common-SNP counts depending only on
the backend chosen.

- **`CM` (used by `--ld-wind-cm`).**
  - Parquet: annotation `CM` is authoritative; the panel/frequency sidecar only
    backfills missing values.
  - PLINK: the window uses `geno_meta["CM"]`, sourced solely from the `.bim`
    third column. Standard PLINK `.bim` files frequently carry an all-zero or
    otherwise uninformative `CM` column. With all-zero `CM` and any positive
    `--ld-wind-cm`, every SNP collapses to the same coordinate, producing a
    whole-chromosome window — which today either trips the generic
    whole-chromosome guard (a misleading error) or, under `--yes-really`,
    silently computes whole-chromosome LD.

- **`MAF` (used by `--maf-min` and `--common-maf-min`).**
  - Parquet: annotation `MAF` is authoritative; the sidecar backfills missing
    values.
  - PLINK: genotype-derived `MAF` is used; annotation `MAF` is ignored.

This makes backend choice change biological/statistical semantics. It also
diverges from legacy LDSC (ldsc2), where the annotation file's `CM`/`MAF` are
dropped and both quantities come exclusively from the reference panel
(`.bim` `CM` for windows, genotype frequencies for `MAF`).

## Goals

1. One consistent, well-defined source of truth for `CM` and `MAF` shared by
   both backends, so backend choice does not change LD windows or counts.
2. Make annotation files population-agnostic: annotation membership is a
   property of a locus, not of a population.
3. Fail loudly and helpfully when `--ld-wind-cm` is requested but no usable
   genetic-map coordinate is available, with concrete remedies — never silently
   compute a meaningless cM window.
4. Preserve correctness of `M`, `M_5_50`, and the partitioned-h2 overlap
   matrices, computed over the retained reference universe with reference-panel
   `MAF`.

## Non-goals

- No change to the R2/LD-score numerical kernels, parquet R2 storage, or the
  public LD-score output layout (`metadata.json`, `baseline.parquet`,
  `query.parquet`) beyond an additive provenance record and an opt-in sidecar.
- No change to `--ld-wind-kb` / `--ld-wind-snps` semantics (they never use `CM`).
- No change to regression (`h2`, `partitioned-h2`, `rg`) consumption of the
  aggregated LD-score directory.

## Core principle

> **The reference panel is the sole source of truth for `CM` and `MAF`.
> Annotation files are population-agnostic and contribute only SNP membership
> and identity keys.**

`CM` (recombination distance) and `MAF` (allele frequency) are population-
specific; they belong to the reference panel, not the annotation. This restores
legacy ldsc2 semantics and makes the two backends agree by construction.

---

## Specification

### 1. Annotation file contract

Annotation files carry **only**:

- **Required:** `CHR`, `POS`.
- **Identity-dependent:** `SNP` is required only in the `rsid` /
  `rsid_allele_aware` identifier modes; `A1` / `A2` are required only in the
  allele-aware modes (`rsid_allele_aware`, `chr_pos_allele_aware`).
- **Annotation value columns:** one or more numeric columns.

`CM` is **no longer a required column**. If `CM` or `MAF` columns are present in
an annotation file they are **ignored** (exactly as legacy ldsc2's
`annot_parser` drops `CHR/BP/SNP/CM`), with a one-time informational log noting
that the reference panel is authoritative for `CM`/`MAF`. Column membership is
resolved by header name, never by positional index. Thinned annotations (value
columns without `CHR`/`POS`/identity columns) remain rejected.

**Why:** the same annotation file can be reused across populations and across
reference panels; `CM`/`MAF` are supplied by whichever panel the run uses.

### 2. Reference panel as the `CM`/`MAF` authority

Both backends source `CM` and `MAF` from the reference panel. Because of the
backends' different data timing, `CM` is resolved by a single shared helper
while `MAF` stays at each backend's natural source — both still "from the
reference panel":

| Quantity | Parquet backend | PLINK backend |
|---|---|---|
| `CM` | panel metadata sidecar (`chrN_meta.tsv.gz`) | `.bim` `CM`, or interpolated from a supplied genetic map (see §4) |
| `MAF` | panel metadata sidecar (folded) | genotype-derived, computed during the BED read (folded) |

The annotation's `CM`/`MAF` are never consulted in either backend.

### 3. `CM` resolution and the `--ld-wind-cm` guard

`CM` is only consulted under `--ld-wind-cm`. A single shared resolver produces
the authoritative per-SNP `CM` for the retained reference universe; the existing
window builder then consumes it unchanged.

- **Parquet.** Use sidecar `CM`. If **any** retained SNP has a missing (`NaN`)
  `CM`, raise a hard error — **never silently drop SNPs** (dropping SNPs from
  the LD reference universe would alter `M`, `M_5_50`, and neighbors' LD
  scores). A properly built panel has no partial-`NaN` `CM`: build-time
  interpolation clamps out-of-range positions to endpoints and interpolates
  across internal gaps, so `NaN` indicates either a panel built without a
  genetic map (entire column `NA`, intended for kb/snps use) or a corrupted
  sidecar — both warrant a clear failure.

- **PLINK.** Resolve `CM` as in §4. If the resolved `.bim` `CM` is **unusable**,
  raise a hard error (see §4).

- **Unusable `CM` definition.** Per chromosome, `CM` is unusable when it has
  **fewer than two distinct finite values** (covers all-zero, all-identical,
  and all-`NaN`). This is the true failure mode: such a column cannot order SNPs
  and so cannot define any window.

- **`--yes-really` interaction.** The unusable-`CM` check is **unconditional**
  and runs **before** the whole-chromosome guard. `--yes-really` continues to
  authorize an expensive whole-chromosome window only when `CM` is *valid*
  (informative and ordered); it does **not** authorize computing a cM window
  from a degenerate `CM` column. Rationale: `--yes-really` means "I accept the
  cost of whole-chromosome LD," not "I accept fabricated coordinates."

### 4. PLINK genetic-map remedy

When `--ld-wind-cm` is used with a PLINK panel, the user may supply a genetic map
so cM windows work even when the `.bim` `CM` is uninformative. This reuses the
`build-ref-panel` machinery (`--genetic-map-hg19-sources` /
`--genetic-map-hg38-sources`, the `load_genetic_map_group` loader, and the
`interpolate_genetic_map_cm` interpolator).

- **New `ldscore` options:** `--genetic-map-hg19-sources`,
  `--genetic-map-hg38-sources` (same names and semantics as `build-ref-panel`).

- **Precedence — an explicit genetic map always wins.** If a genetic map is
  supplied, `CM` is interpolated at the `.bim` positions and used for **all**
  chromosomes, ignoring the `.bim` `CM` entirely. This yields one uniform `CM`
  source per run (reproducible), regardless of `.bim` `CM` quality. If no map is
  supplied: use `.bim` `CM` when usable; otherwise error.

- **Genome-build selection.** The genetic map must match the build of the
  `.bim` positions. The build is the resolved genome build of the run:
  `--genome-build` when given, otherwise inferred (`genome_build="auto"`, valid
  only in `chr_pos`-family identifier modes). The matching
  `--genetic-map-<build>-sources` must be present. In `rsid`-family modes (where
  `auto` inference is unavailable), the user must pass `--genome-build` to use a
  genetic map.

- **Interpolation behavior.** `interpolate_genetic_map_cm` uses linear
  interpolation with endpoint clamping (no `NaN` within a covered chromosome),
  and errors if a chromosome has no map coverage at all.

- **Error when no usable `CM` and no map** (`--ld-wind-cm`, PLINK, unusable
  `.bim` `CM`): a single actionable error offering three remedies —
  (a) provide a `.bim` whose `CM` carries real genetic-map positions;
  (b) add `--genetic-map-hg38-sources` / `--genetic-map-hg19-sources` for the
  panel's build; or (c) use `--ld-wind-kb` / `--ld-wind-snps`. A dedicated
  section is added to `docs/troubleshooting.md` and linked from the message.

### 5. `MAF`: always required, inclusive thresholds

- **Always required.** A reference panel must provide `MAF`; its absence is a
  hard error. For PLINK this is automatic (MAF is intrinsic to genotypes); for
  parquet the sidecar-`MAF` requirement becomes **unconditional** (no longer
  gated on a "strict metadata" flag). `build-ref-panel` always emits `MAF` in
  the sidecar. Rationale: `M_5_50` and `--maf-min` are meaningless without
  `MAF`, and the partitioned-h2 common-overlap correction depends on it.

- **Inclusive `>=` everywhere.** Both `--maf-min` (the SNP filter) and
  `--common-maf-min` (the common-count / overlap threshold) use inclusive `>=`
  comparisons. This is already the implemented behavior across both backends,
  the reference-panel filter, and munge-sumstats; it is locked here as an
  invariant with boundary tests. (Note: this is a deliberate, documented
  departure from legacy ldsc2's strict `MAF > 0.05` for common counts.)

- **Distinction preserved.** `--maf-min` removes SNPs *before* LD-score
  computation (changes which SNPs contribute to LD). `--common-maf-min` only
  partitions the *retained* SNPs into the common-count `M_5_50` and the common
  overlap; it does not drop SNPs.

- **`--maf-min` applies identically in both backends.** Today `--maf-min` is
  applied as a reference-panel filter only for parquet (in
  `ref_panel.load_metadata` → `_apply_maf_filter`); the PLINK path passes
  `maf_min=None` to the kernel, so `--maf-min` is silently *not* applied for
  PLINK. This design makes `--maf-min` a reference-panel filter for **both**
  backends so the retained reference universe is identical given the same
  `--maf-min`. (The kernel's now-dead `apply_maf_filter` path is removed or made
  consistent as part of the cleanup.)

### 6. Counts and overlap consistency

`M`, `M_5_50`, the total overlap matrix, and the common overlap matrix are
computed over the **retained reference universe**, identically in both backends.
The retained reference universe is:

> reference-panel SNPs ∩ (`ref-panel-snps` / curated-HM3 keep-restriction),
> minus excluded regions, MAF-filtered by `--maf-min`, intersected with the
> annotation SNP set.

The common subset uses reference-panel `MAF` with `MAF >= common_maf_min`. This
universe and these counts are already implemented; the design's contribution is
to make the `MAF` feeding the common split come from the reference panel in both
backends, so `M_5_50` and the overlap matrices agree across backends for
identical inputs. Locked as an invariant with a cross-backend regression test.

### 7. Provenance and optional metadata export

Because PLINK `CM`/`MAF` are otherwise transient (recomputed each run, not
persisted per-SNP in the output), the design adds:

- **Run provenance (always, additive).** the LD-score root `metadata.json` and the log record the
  `CM` source (`.bim` vs genetic-map file(s) + resolved genome build) and the
  `MAF` source (genotype-derived for PLINK; sidecar for parquet).

- **Opt-in metadata export.** A new `ldscore` flag `--export-ref-metadata`
  (**default `false`**, PLINK-scoped) writes a per-chromosome
  `chrN_meta.tsv.gz` sidecar in the **same schema as the parquet panel sidecar**
  — gzip TSV with columns `CHR POS SNP A1 A2 CM MAF` (folded `MAF`,
  validated/interpolated `CM`). This unifies the metadata format across backends
  and lets an exported sidecar seed a future `build-ref-panel` or a QC diff.
  Parquet runs already ship this sidecar, so the flag is unnecessary there.

For durable, reusable `CM`/`MAF`, the recommended path remains
`ldsc build-ref-panel` (a parquet panel persists `CM`/`MAF` + R2 as canonical
artifacts); PLINK mode is "bring-your-own-panel, compute on the fly."

---

## Error and UX summary

| Situation | Behavior |
|---|---|
| Annotation file contains `CM`/`MAF` columns | Ignored; one-time info log. |
| `rsid`-mode annotation missing `SNP` | Clear error (identity key absent). |
| `chr_pos`-mode annotation missing `SNP` | Works (SNP optional). |
| `--ld-wind-cm`, parquet sidecar has any `NaN` `CM` | Hard error; rebuild panel with a genetic map, or use kb/snps. No silent drop. |
| `--ld-wind-cm`, PLINK, unusable `.bim` `CM`, no map | Hard error with 3 remedies (real-`CM` `.bim` / genetic map / kb-snps). Fires even under `--yes-really`. |
| `--ld-wind-cm`, PLINK, genetic map supplied | Interpolate at `.bim` positions for all chromosomes (map wins). |
| Genetic map supplied but build unresolved (`rsid` mode, no `--genome-build`) | Error asking for explicit `--genome-build`. |
| Reference panel lacks `MAF` | Hard error (always required). |
| SNP with `MAF == maf_min` / `MAF == common_maf_min` | Kept / counted (inclusive `>=`). |
| `--maf-min` with a PLINK panel | Drops `MAF < maf_min` SNPs (same as parquet); previously silently ignored. |
| `--genetic-map-*-sources` with a parquet panel | Warning; flags ignored. CM comes from the sidecar (authoritative); genetic-map flags apply only to PLINK. |
| PLINK chromosome with no SNPs left after genotype filtering | Distinct "retained no PLINK reference SNPs" error (all monomorphic / `--maf-min` too high / no overlap) — not the missing-`MAF` error. |

## Backward compatibility and migration

- **Annotation files with a `CM` column keep loading** — the column is ignored,
  not rejected. Files that previously relied on annotation `CM` to drive parquet
  cM windows must now use a reference panel that provides `CM` (its sidecar, or
  a genetic map for PLINK). This is the intended correction, not a regression:
  population-specific `CM` should come from the panel.
- **Parquet panels built without a genetic map** (sidecar `CM = NA`) can no
  longer be used with `--ld-wind-cm`; they continue to work with
  `--ld-wind-kb` / `--ld-wind-snps`. Rebuild with a genetic map for cM windows.
- **Parquet panels built without `MAF`** (if any exist) now hard-error; rebuild
  with `MAF`. New panels always include `MAF`.
- **PLINK runs are unchanged for `MAF`** and gain the genetic-map remedy plus
  the unusable-`CM` guard for `--ld-wind-cm`.
- **`annotate` output keeps the `CM` column.** This is a *read*-side change for
  `ldscore`; the `annotate` *write* format is unchanged. Legacy ldsc2 reads
  full-annot files **positionally** — `np.array(annot.df.iloc[:, 4:])` assuming
  columns `CHR BP SNP CM` precede the annotations (ldsc2 `ldsc.py:153`). To stay
  consumable by legacy ldsc2, `annotate` must keep emitting
  `CHR BP SNP CM <annotations…>`. The `CM` value is a placeholder (legacy skips
  it positionally; our reader ignores it by name), so `annotate` needs no genetic
  map. Both annotation readers in this package (`AnnotationBuilder.parse_annotation_file`
  and the kernel `parse_annotation_file`) select columns **by name**, so they
  ignore the placeholder `CM` regardless.

## Alignment with legacy ldsc2

This design restores legacy semantics rather than inventing new ones:

- ldsc2's `annot_parser` drops `SNP/CHR/BP/CM` from annotation files; only
  annotation value columns are used.
- ldsc2 cM windows use the `.bim` `CM` (`array_snps.df['CM']`), and `M_5_50`
  uses genotype `MAF` (`geno_array.maf`).

The parquet backend's "annotation-first `CM`/`MAF`" behavior was introduced in
the Python-3 port; this design removes it so both backends match legacy and each
other.

## Testing strategy (TDD; tests written first)

Numerical/behavioral tests against constructed fixtures:

1. **Cross-backend window consistency.** Identical annotation + identical
   reference `CM` → parquet and PLINK produce identical block-left arrays /
   cM windows.
2. **Cross-backend count consistency.** Identical inputs → identical `M`,
   `M_5_50`, total overlap, and common overlap across backends.
3. **All-zero `.bim` `CM` + genetic map.** PLINK uses interpolated `CM` and
   matches the parquet result built from the same map.
4. **All-zero `.bim` `CM`, no map.** Clear unusable-`CM` error, **including
   under `--yes-really`**; whole-chromosome guard remains bypassable for valid
   `CM`.
5. **Annotation `CM` ignored.** An annotation file with a populated `CM` column
   yields a result identical to the same file without `CM`.
6. **Parquet `NaN` `CM` under `--ld-wind-cm`.** Hard error; no silent drop.
7. **Missing reference `MAF`.** Hard error (parquet sidecar without `MAF`).
8. **Inclusive thresholds.** A SNP at exactly `maf_min` is retained; a SNP at
   exactly `common_maf_min` is counted as common.
9. **Identity-mode required columns.** `rsid`-mode annotation without `SNP`
   errors; `chr_pos`-mode annotation without `SNP` works.
10. **Provenance + export.** the root `metadata.json` records `CM`/`MAF` sources;
    `--export-ref-metadata` writes a sidecar matching the parquet schema.
11. **Cross-backend `--maf-min`.** `--maf-min t` drops the same SNPs (and yields
    the same `M`/`M_5_50`) in parquet and PLINK for identical inputs.

## Out of scope / future work

- A genetic-map interpolation path for parquet beyond what `build-ref-panel`
  already provides.
- Persisting per-SNP `CM`/`MAF` inside `baseline.parquet` (kept lean).
- Unifying `--common-maf-min` with legacy strict `>` (explicitly rejected:
  inclusive `>=` is the project convention).

## Resolved decisions (brainstorm record)

1. Annotation `CM`/`MAF`: **ignored entirely** (legacy-style), reference panel
   authoritative.
2. CM/MAF sourcing architecture: **shared `CM` resolver + per-backend `MAF`**.
3. PLINK cM remedy: **error + genetic-map interpolation** via
   `--genetic-map-{hg19,hg38}-sources`, with the three-remedy error message.
4. CM precedence (PLINK): **explicit genetic map always wins**.
5. Unusable `CM`: **fewer than two distinct finite values per chromosome**;
   check is **unconditional** (fires under `--yes-really`), runs before the
   whole-chromosome guard.
6. `MAF`: **always required** (hard error if absent); `--maf-min` and
   `--common-maf-min` both inclusive `>=`.
7. Parquet `NaN` `CM` under `--ld-wind-cm`: **reject** (hard error), never drop.
8. Counts/overlap over the **retained reference universe** with reference-panel
   `MAF` (already implemented; locked as invariant).
9. Persistence: **manifest provenance + opt-in `--export-ref-metadata`**
   (default `false`), sidecar in the unified `chrN_meta.tsv.gz` schema.
10. `--maf-min` is applied as a **reference-panel filter in both backends**
    (fixes the PLINK gap where `--maf-min` was silently ignored), with a
    cross-backend test asserting identical retained universes.
11. **Both** annotation readers adopt the new contract: the public
    `AnnotationBuilder.parse_annotation_file` and the kernel
    `parse_annotation_file` (`_kernel/ldscore.py:636`). The kernel reader is also
    checked for reachability; if dead, it is flagged for separate removal.
12. **`annotate` keeps writing `CM`** (placeholder) to preserve legacy ldsc2
    full-annot positional compatibility (`iloc[:, 4:]`). `annotate`'s write
    format is otherwise unchanged and out of scope here.
13. The **`CM`-not-used-downstream caveat** is documented both in a log message
    (when `annotate` writes the placeholder column) and in `docs/current/`.
14. **`--genetic-map-*-sources` is ignored in parquet mode** with a warning
    (sidecar `CM` is authoritative).
15. A PLINK chromosome left with **no SNPs after genotype filtering** raises a
    distinct empty-universe error, separate from the missing-`MAF` error.
