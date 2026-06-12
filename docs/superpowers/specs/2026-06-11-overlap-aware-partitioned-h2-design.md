# Overlap-Aware Partitioned Heritability: Design Spec

Status: proposed (2026-06-11)
Scope: `ldscore` (overlap-matrix computation + storage), `partitioned-h2`
regression (overlap-aware category summaries), docs, tests.
Non-goals: changing `h2`/`rg`; changing the LD-score numerical kernel; joint
co-fitting of multiple query annotations in one model.

## 1. Problem

Legacy LDSC2 `--overlap-annot` does more than fit a multiple regression. After
fitting, it re-reads the `.annot` files, builds the annotation overlap matrix
`O = AᵀA`, and uses it to compute overlap-aware category summaries
(`Prop._SNPs`, `Prop._h2` ± SE, `Enrichment` ± SE, `Enrichment_p`). LDSC3's
`partitioned-h2` carries forward only the regression-ready LD scores and the
marginal count vectors; it never carries the annotation matrix `A`, so it cannot
form `O`. As a result LDSC3 currently reports **disjoint-style** category
summaries that are correct only when annotations do not overlap — which is
essentially never, because a baseline model almost always contains an all-ones
`base` annotation that overlaps every category.

What LDSC3 reports correctly today is the **conditional** coefficient family
(`Coefficient = tau_c`, its SE, its p-value), because those come straight from the
fitted regression and need no overlap matrix.

The only missing ingredient is the overlap matrix `O`. Every other input to the
legacy formula (`prop`, `prop_cov`, `coef`, `coef_cov`, `n_blocks`, `cat`,
`cat_se`) already exists on the fitted `Hsq`, and the legacy formula itself is
already ported verbatim as `Hsq._overlap_output`
(`src/ldsc/_kernel/regression.py:471`).

## 2. Background: one model, two analyses

Both analyses share the LDSC model

```
E[χ²_j] = N · Σ_c tau_c · ℓ(j,c) + N·a + 1
```

where `ℓ(j,c) = Σ_k A[k,c]·r²_{jk}` is SNP `j`'s LD score for annotation `c`,
`tau_c` is the per-SNP heritability contribution of annotation `c`, and `a` is a
confounding term. They differ in how annotations enter the regression and which
quantity is the headline.

| | Functional-category partitioned h² | Cell-type-specific |
|---|---|---|
| Reference | Finucane et al. 2015 | Finucane et al. 2018 |
| Legacy entry point | `--overlap-annot` (`ldscore/sumstats.py:estimate_h2`) | `--h2-cts` (`ldscore/sumstats.py:cell_type_specific`) |
| Regression structure | one **joint** fit of all categories | one fit per cell type: **baseline + one** annotation |
| `tau_c` conditional on | all other categories | baseline only |
| Headline | **Enrichment** (+ `Enrichment_p`) | **Coefficient `tau_c`** (+ p) |
| Headline p-value | two-sided t, df = `n_blocks` (`2·t.sf(|·|, n_blocks)`) | one-sided normal `norm.sf(tau/se)`, H₁: `tau_c > 0` |

The two regressions exist for a practical reason, not a mathematical one:
cell-type annotations are too collinear across tissues to co-fit, so they are
tested one at a time against a fixed baseline. The estimator `(tau, cov(tau))` is
identical in both; only the design matrix and the chosen summary differ.

### 2.1 Mapping to LDSC3

LDSC3 `partitioned-h2` always fits `baseline + one query` at a time
(`regression_runner.py:estimate_partitioned_h2_batch`). This is mechanically the
cell-type-specific model. Therefore LDSC3 detects the analysis **by input
shape**, with no user flag:

- **Functional regime** — the LD-score directory has **no query columns**
  (the user supplied only baseline annotations). `partitioned-h2` performs the
  single joint fit of all baseline annotations and reports one overlap-aware row
  per baseline category. This reproduces Finucane-2015 / legacy `--overlap-annot`.
- **Cell-type regime** — the LD-score directory **has query columns**.
  `partitioned-h2` fits `baseline + one query` per query and reports one row per
  query (full baseline-plus-query tables available on request).

Both regimes emit the **same column schema** and call the **same** ported
`_overlap_output`; the regime only selects which model(s) to fit and which rows
to surface. The classic Finucane-2015 joint enrichment of baseline categories is
always recoverable from the baseline rows of a cell-type-regime full table,
because every per-query model fits the entire baseline jointly.

This is a deliberate behavior change: `partitioned-h2` previously **rejected**
baseline-only directories (`PARTITIONED_H2_REQUIRES_QUERY_ANNOTATIONS_MESSAGE`,
`regression_runner.py:140`). It now accepts them as the functional regime.

## 3. Overlap-aware math

Let a fitted model have `K` annotations over a SNP universe of size `M_tot`
(the number of reference SNPs; for the common universe, the number of common
reference SNPs). Define:

- `A ∈ ℝ^{M_tot × K}` — annotation matrix (binary or continuous), `A[s,c]`.
- `O = AᵀA ∈ ℝ^{K×K}` — overlap matrix, `O[i,j] = Σ_s A[s,i]·A[s,j]`.
  For binary annotations `O[i,j]` is the number of SNPs in both `i` and `j`, and
  `O[i,i] = M_i`. For continuous annotations it is a co-sum (see §6).
- `M ∈ ℝ^K` — marginal counts, `M[c] = Σ_s A[s,c]` (the count vector LDSC3
  already stores). Note `diag(O)[c] = Σ_s A[s,c]²`, which equals `M[c]` only for
  binary annotations.
- `tau = coef`, `cov(tau) = coef_cov`, `prop`, `prop_cov`, `cat = M⊙tau`,
  `cat_se`, `n_blocks` — all taken from the fitted `Hsq`.

`prop[c] = cat[c]/Σ_k cat[k] = M[c]·tau[c] / h²_tot`, where
`h²_tot = Σ_k cat[k] = hsq.tot` is the model's total heritability.

The overlap-aware columns below — `Prop._SNPs`, `Prop._h2` ± SE, `Enrichment`
± SE, `Enrichment_p`, `Coefficient`, `Coefficient_std_error`, and the z-score —
are produced verbatim by the ported `Hsq._overlap_output`
(`src/ldsc/_kernel/regression.py:471`, mirroring legacy `regressions.py:394`).
The runner then augments that table with columns the ported function does not
emit: the one-sided `Coefficient_p`, the conditional `Category_h2` /
`Category_h2_std_error` (from `hsq.cat`/`hsq.cat_se`), and `overlap_aware`.
Legacy names the z-score `Coefficient_z-score`; LDSC3 renames it `Coefficient_z`.

**Proportion of SNPs.** `Prop._SNPs[i] = M[i] / M_tot`.
For an all-ones `base` category, `M[base] = M_tot`, so `Prop._SNPs[base] = 1`.

**Proportion of h² (overlap-aware, marginal).** With
`P[i,j] = O[i,j]/M[j]`,

```
Prop._h2[i] = Σ_j P[i,j]·prop[j] = (Σ_j O[i,j]·tau[j]) / h²_tot
            = (Σ_s A[s,i]·h²_s) / h²_tot,   where h²_s = Σ_j A[s,j]·tau[j].
```

This is the fraction of total h² attributable to **SNPs that fall in category
`i`**, folding in every category those SNPs also belong to. It is *marginal*,
not a clean partition; values across categories do not sum to 1 under overlap.
Standard error by the delta method:
`Cov(Prop._h2) = P · prop_cov · Pᵀ`,
`Prop._h2_std_error[i] = sqrt(max(0, [P·prop_cov·Pᵀ]_{ii}))`.

**Enrichment.** `Enrichment[i] = Prop._h2[i] / Prop._SNPs[i]`,
`Enrichment_std_error[i] = Prop._h2_std_error[i] / Prop._SNPs[i]`.

**Enrichment p-value (coefficient contrast).** With

```
D[i,:] = O[i,:]/M[i] − (M − O[i,:])/(M_tot − M[i])      if M_tot ≠ M[i]
D[i,:] = 0                                              if M_tot =  M[i]
```

`diff_est = D·tau`, `diff_cov = D·cov(tau)·Dᵀ`, `diff_se = sqrt(diag(diff_cov))`,
`Enrichment_p[i] = 2·t.sf(|diff_est[i]/diff_se[i]|, n_blocks)`
(NA if `diff_se[i] = 0`). The contrast `diff_est[i]` is the mean per-SNP h²
**inside** category `i` minus the mean per-SNP h² **outside** it. The
`M_tot = M[i]` guard zeroes out a category that contains every SNP (e.g. `base`),
which therefore has no "outside" and reports `Enrichment_p = NA`.

**Coefficient (conditional).** `Coefficient[i] = tau[i]`,
`Coefficient_std_error[i] = coef_se[i]`,
`Coefficient_z[i] = tau[i]/coef_se[i]`,
`Coefficient_p[i] = norm.sf(Coefficient_z[i])` — **one-sided**, H₁: `tau_i > 0`,
matching legacy `cell_type_specific` (`sumstats.py:302`). The sign is visible
from `Coefficient` itself.

**Category h² (conditional).** `Category_h2[i] = cat[i] = M[i]·tau[i]`,
`Category_h2_std_error[i] = cat_se[i]`. This is the category's own
coefficient-weighted contribution and can be negative under overlap. It is a
different quantity from `Prop._h2`: conditional vs marginal.

**Overlap indicator.** `overlap_aware = (max_{i≠j} O[i,j] > 0)` over the
retained model columns — i.e. some SNP belongs to ≥2 retained categories. For a
baseline that includes `base`, this is effectively always true. It is a
per-fitted-model property: in the functional regime one value labels the whole
table; in the cell-type regime each query is its own model, so it is reported
per row and may differ across queries. The flag is informational: the
overlap-aware formula is the correct generalization either way (the no-overlap
case is just `O` diagonal), so there is **no formula branching** — we always
compute overlap-aware and report the flag.

### 3.1 Why marginal and conditional are both reported

For a query annotation the two views answer different questions:

- `Enrichment` (marginal, baseline-confounded): "do this annotation's SNPs carry
  disproportionate heritability?" If the annotation's SNPs sit in
  generally-enriched regions (coding, conserved), it looks enriched for that
  reason alone.
- `Coefficient tau` (conditional on baseline, de-confounded): "does this
  annotation add heritability **beyond** the baseline?" This is the
  cell-type-specific headline.

The full table carries both `Category_h2` (conditional) and `Prop._h2`
(marginal) so neither interpretation is lost.

## 4. The overlap matrix: storage and block reuse

### 4.1 What must be stored

LDSC3 never co-fits two query annotations, so cross-query overlaps `O[q,q']`
(`q ≠ q'`) are never needed. The model overlap matrix for query `q` is

```
        baseline      query q
      ┌─────────────┬──────────┐
 bl   │   O_BB      │  O_Bq    │
      ├─────────────┼──────────┤
 q    │   O_qB=O_Bqᵀ│  O_qq    │
      └─────────────┴──────────┘
```

so the only data required across the whole batch is:

- `O_BB` — baseline×baseline block (`B×B`), reused across all queries and the
  sole block needed by the functional regime.
- `O_Bq` — baseline×query block (`B×Q`), one column per query.
- `O_qq` — each query's self-overlap `Σ_s A[s,q]²` (`Q`-vector). Equals the query
  count for binary annotations but **not** for continuous ones, so it is stored
  explicitly.

Storage is `O(B·(B+Q))`, linear in the number of query annotations, not
quadratic. For `B≈100`, `Q≈10⁴` this is ~1M small rows — trivial.

### 4.2 Computation (at `ldscore` time)

The annotation matrix `A` is in memory next to `compute_counts`
(`_kernel/ldscore.py:1637` parquet, `:1780` PLINK), columns ordered
`baseline + query` (write `A_B` for the baseline columns of `A`). Per chromosome
we accumulate, alongside `M`/`M_5_50`:

```
G_all      += A_Bᵀ · A                         # (B × (B+Q)) baseline-rows block
G_common   += (A_B[c])ᵀ · A[c]                 # c = common mask
dq_all[q]  += Σ_s A[s,q]²                       # query self-overlaps
dq_common  += Σ_{s∈c} A[s,q]²
n_all      += rows;   n_common += c.sum()       # M_tot per universe
```

where the common mask `c` is `MAF >= common_maf_min` (inclusive; see §5). These are
summed across chromosomes in `_aggregate_chromosome_results` exactly like the
count vectors. The cost is one extra matrix product on data already resident in
memory, negligible against the LD-score (R²) computation.

### 4.3 Persisted form (storage decision: parquet sidecar)

A dedicated sidecar `ldscore.overlap.parquet`, long/tidy, name-keyed:

| column | meaning |
|---|---|
| `row_annotation` | a baseline annotation, or a query annotation (diagonal rows) |
| `col_annotation` | any baseline/query annotation |
| `overlap_all_snps` | `O[row,col]` over the all-SNP universe |
| `overlap_common_snps` | `O[row,col]` over the common-SNP universe; null if MAF unavailable |

Rows: `(baseline_i, c)` for every baseline `i` and every column `c`
(the baseline-rows block `O_BB | O_Bq`), plus `(query_q, query_q)` for each query
(the self-overlaps). `metadata.json` gains:

```json
"files":  { "overlap": "ldscore.overlap.parquet", ... },
"overlap_config": {
  "total_all_reference_snps": 9997231,
  "total_common_reference_snps": 5821004,
  "common_maf_min": 0.05,
  "common_maf_operator": ">=",
  "stored_block": "baseline_rows_plus_query_diagonal"
}
```

Rationale: this mirrors the existing parquet-data / JSON-metadata split, aligns
by **name** (robust to zero-variance column dropping and reordering), scales to
thousands of query annotations, and keeps `metadata.json` lean and readable.
Embedding the arrays in `metadata.json` was rejected because `B·(B+Q)` floats
bloat the JSON to tens of MB for cell-type panels; per-column-record or `.npz`
encodings were rejected as ad-hoc or position-aligned (fragile).

### 4.4 Reconstruction (at regression time)

For a fitted model with retained columns `R ⊆ baseline ∪ {q}` (after
zero-variance dropping), the model overlap matrix `O_R` is assembled by name:

- baseline–baseline entries from `O_BB`,
- baseline–query entries from `O_Bq` column `q`,
- query–baseline by symmetry,
- query–query from the stored self-overlap `O_qq`.

The count universe is chosen consistently: `common` columns of the sidecar with
`common` counts (default), `all` columns with `all` counts (or when common is
unavailable), matching `dataset.count_key_used_for_regression`. `M_tot` is the
matching `total_*_reference_snps` scalar. The assembled `O_R`, `M`, and `M_tot`
go straight into `hsq._overlap_output`.

## 5. SNP universes and the common-MAF filter

Counts and the overlap matrix are computed over the reference SNP universe
(`ld_reference_snps`), the same universe legacy uses (it reads `.M`/`.M_5_50` and
re-reads `.annot`, never the regression subset). Two universes:

- **all** — every retained reference SNP. `M_tot = total_all_reference_snps`.
- **common** — `MAF >= common_maf_min` (default 0.05).
  `M_tot = total_common_reference_snps`.

The common filter uses **inclusive `MAF >= common_maf_min`** (`_kernel/ldscore.py`),
following the project convention that `-min` thresholds are inclusive (`>=`) and
`-max` thresholds inclusive (`<=`); the exact boundary `MAF == 0.05` is therefore
included. The `count_config.common_reference_snp_maf_operator` string is `">="`.
The overlap common mask reuses this same mask, so its diagonal stays equal to the
common-count vector `M_5_50`. (An earlier draft of this spec used strict `>` to
mirror legacy's `0.05 < FRQ`; that was reconsidered in favor of the inclusive
convention — see `docs/superpowers/specs/2026-06-11-cm-maf-source-of-truth-design.md`.)

Regression continues to default to common counts when available and fall back to
all counts otherwise (`--count-kind`), and the overlap universe follows that
choice.

## 6. Continuous annotations

LDSC annotations come in three common shapes, all of which the program treats
identically as real-valued columns:

- **Binary functional categories** — coding, UTR, promoter, enhancer, conserved,
  histone marks, DHS, etc. Values in `{0,1}`.
- **Binned continuous features** — e.g. **MAF** and several LD-related features
  are split into quantile **bins**, each bin a separate 0/1 indicator column
  (e.g. `MAFbin_1 … MAFbin_10`). The binning is performed **upstream** when the
  `.annot` file is constructed, not by this program.
- **Raw continuous features** (baselineLD; Gazal et al. 2017) — recombination
  rate, nucleotide diversity, background-selection statistic, CpG content, GERP,
  predicted allele age, level of LD. These are real-valued, typically
  standardized upstream.

The program performs **no binning, binarization, or standardization**. Every
annotation column is read as `float32` (`_kernel/ldscore.py:706`) and
`M[c] = Σ_s A[s,c]`, `O[i,j] = Σ_s A[s,i]·A[s,j]` are computed verbatim. The
overlap-aware formula is pure linear algebra and is therefore unchanged. The only
shift is interpretive: for a continuous annotation `M[c]` is the summed
annotation weight rather than a SNP count, `O[i,i] = Σ_s A[s,c]² ≠ M[c]`, and
`Prop._SNPs` reads as "proportion of summed annotation weight." This matches
legacy exactly (`annot_parser`/`annot`, `parse.py:125,180`). Because continuous
columns have generically nonzero off-diagonal `O`, they always report
`overlap_aware = true`, which is harmless.

## 7. Output schema

Both regimes write **one** schema to `partitioned_h2.tsv`; they differ only in
which rows appear and the default sort, never in columns:

```
Category, Prop._SNPs,
Category_h2, Category_h2_std_error,          # conditional M_c·tau_c
Prop._h2, Prop._h2_std_error,                # marginal, overlap-aware
Enrichment, Enrichment_std_error, Enrichment_p,
Coefficient, Coefficient_std_error, Coefficient_z, Coefficient_p,
overlap_aware
```

A single schema keeps the artifact self-contained (the headline metric's SE is
always present, so no `--write-per-query-results` is needed just to read a
standard error) and keeps both regimes — and the writer code — identical apart
from rows and sort.

Row meaning, files, and default ordering by regime:

- Functional (no query): one joint fit; rows are **all baseline categories**,
  written in baseline-annotation order (`--summary-sort-by` resolves to
  `category`, matching legacy `.results`). Headline column is `Enrichment`
  (+ two-sided `Enrichment_p`).
- Cell-type (query present): one row **per query**, sorted by `Coefficient_p`
  ascending by default (most-significant first), matching legacy
  `cell_type_specific` (`sumstats.py:310`). Headline column is `Coefficient`
  (+ one-sided `Coefficient_p`). `--write-per-query-results` additionally emits,
  per query under `diagnostics/query_annotations/<q>/`, the full
  baseline-plus-query table (all `B+1` model rows), whose baseline rows are the
  joint functional enrichments conditional on that query.

`--summary-sort-by` defaults to `auto`, which resolves to `category` in the
functional regime and `coefficient-p` in the cell-type regime; any explicit
choice is honored in both.

`Enrichment_p` is two-sided t (df = `n_blocks`) in both regimes; `Coefficient_p`
is one-sided normal in both.

### 7.1 Headline signaling and interpretation

The TSV is never decorated (stable, scriptable columns); the default row sort is
its only in-file emphasis. The regime and how to read it are surfaced two ways.

**Log banner** (human-facing), one block per run, naming the regime and the
column to focus on:

- Functional: "Functional-category regime: joint fit of B baseline categories.
  Focus on `Enrichment` and `Enrichment_p`: `Enrichment > 1` means the category's
  SNPs explain a larger share of heritability than their share of SNPs (`< 1`
  means a smaller share); a small `Enrichment_p` indicates the enrichment is
  significantly different from 1, i.e. significantly larger or smaller."
- Cell-type: "Cell-type-specific regime: baseline + one query per model
  (N queries). Focus on `Coefficient` together with the one-sided `Coefficient_p`
  (a one-sided test of whether `Coefficient` > 0): a **positive** `Coefficient`
  with a small `Coefficient_p` means
  the query annotation contributes heritability **beyond** the baseline
  annotations. The `Enrichment` column is confounded by the query annotation's
  overlap with the baseline annotations, so use `Coefficient_p` to judge whether
  the additional contribution is significant."

**`diagnostics/metadata.json`** (machine-readable, self-describing):
`analysis_type` (`functional_category` | `cell_type_specific`),
`headline_metric` (`enrichment` | `coefficient`),
`enrichment_p_test` (`two_sided_t`), `coefficient_p_test` (`one_sided_greater`).

## 8. Collinearity hard error

Restore legacy's condition-number guard (`_check_ld_condnum`, `sumstats.py:175`),
which LDSC3 dropped, as a **hard error** (legacy's behavior). At regression time,
after the retained LD-score design matrix `X` is assembled, if `cond(X) > 1e5`
(legacy's threshold on the LD-score design matrix itself) raise `LDSCInputError`
naming the most collinear annotation pair — derived from the overlap matrix as
`O[i,j]/sqrt(O[i,i]·O[j,j])` — with a suggested action (drop one of the pair, or
coarsen the annotation). Aborting (rather than warning) keeps every emitted
partition trustworthy: collinear annotations cannot be separated, so the run
stops instead of asking a non-statistician user to adjudicate a warning. Checked
once for the functional joint fit; checked per model in the cell-type regime,
where the first collinear model aborts the batch. Location: the regression
workflow layer (`regression_runner.py`), where `X` and the overlap matrix are
both available.

## 9. Backward compatibility

New `ldscore` runs always write `ldscore.overlap.parquet`. `partitioned-h2`
against a directory lacking it fails fast with a regenerate message and a
`docs/troubleshooting.md` pointer; `h2` and `rg` are unaffected (they never need
the overlap matrix). No silent fallback to disjoint numbers.

## 10. Testing strategy

- **Formula parity.** Build a small hand-specified `A`, form `O = AᵀA`, construct
  a fitted `Hsq` with known `prop/prop_cov/coef/coef_cov`, and assert the
  augmented table equals values from the reference formula (the legacy
  `_overlap_output` expressions, transcribed into the test as a pure-NumPy
  oracle since legacy is Python 2). Covers `Prop._SNPs`, `Prop._h2` ± SE,
  `Enrichment` ± SE, two-sided `Enrichment_p`, one-sided `Coefficient_p`,
  conditional `Category_h2`.
- **Hand-computed ground truth.** A 6-SNP, 3-annotation overlapping example with
  every quantity computed by hand, including the `base` `Enrichment_p = NA` guard.
- **Block-assembly equivalence.** Assert that `O_R` reconstructed from the stored
  baseline-rows block + query diagonal equals a direct `AᵀA` on the full model.
- **Non-overlap reduction.** With disjoint annotations (plus `base`), the
  overlap-aware numbers match the closed-form disjoint values and
  `overlap_aware` is reported per the off-diagonal rule.
- **Continuous annotations.** A continuous column reproduces the legacy co-sum
  behavior (`O[i,i] ≠ M[i]`).
- **Strict-MAF boundary.** A SNP at exactly `MAF == 0.05` is excluded from common
  counts and the common overlap universe.
- **Regimes & I/O.** Functional (no-query) end-to-end produces a baseline-category
  table; the cell-type primary and per-query model tables carry the new columns; an
  old directory without the sidecar raises the regenerate error.

## 11. Touch points (summary; details in the plan)

- `_kernel/ldscore.py` — compute overlap blocks beside `compute_counts`; keep the
  inclusive `MAF >= common_maf_min` common mask; thread blocks through `ChromComputationResult`.
- `ldscore_calculator.py` — carry/sum overlap blocks in `_LegacyChromResult` →
  `ChromLDScoreResult` → `LDScoreResult`; inclusive `>=` operator in `count_config`.
- `outputs.py` — write `ldscore.overlap.parquet`; add `files.overlap` +
  `overlap_config` to `build_metadata`; `PartitionedH2DirectoryWriter` writes the
  single schema as the primary `partitioned_h2.tsv` in both regimes, with
  `--write-per-query-results` adding the per-query model-detail folders.
- `regression_runner.py` — load overlap sidecar in `load_ldscore_from_dir`;
  assemble `O_R` aligned to retained columns/universe; add the two regimes;
  replace disjoint logic in `summarize_partitioned_h2` with `_overlap_output` +
  augmentation; one-sided `Coefficient_p`; collinearity hard error; remove the
  "requires query annotations" rejection; change `--summary-sort-by` default to
  `auto` (regime-aware: `category` functional, `coefficient-p` cell-type); emit
  the regime log banner with interpretation guidance and the `analysis_type` /
  `headline_metric` / p-test fields in `diagnostics/metadata.json`.
- `docs/current/partitioned-ldsc-workflow.md`, `docs/troubleshooting.md`,
  `design_map.md` — update.
- `tests/` — per §10.
