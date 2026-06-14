# Regression Configuration Reference (Advanced)

This document describes the tunable configuration of the regression module
(`ldsc h2`, `ldsc rg`, `ldsc partitioned-h2`). **You do not need any of this for
a standard run.** The defaults reproduce the original LDSC behavior, and the
common workflow is simply to point each command at a munged `.sumstats` file and
a canonical LD-score directory. The parameters below are opt-in controls for
advanced users who need to constrain the intercept, exclude outliers, tune the
jackknife, or change the per-SNP normalization.

For the plain argument list (paths, output flags) see
[`io-argument-inventory.md`](io-argument-inventory.md). This document focuses on
the *statistical* configuration and what each knob does to the estimator.

Source of truth in code: `RegressionConfig` in `src/ldsc/config.py`,
`_runner_from_args` / `estimate_h2` / `_fit_rg_dataset` /
`estimate_partitioned_h2_batch` in `src/ldsc/regression_runner.py`, and the
kernel in `src/ldsc/_kernel/regression.py`.

---

## 1. Math backbone

LDSC fits a **weighted linear regression** of a per-SNP association statistic on
LD scores, then turns the slope and intercept into heritability / genetic
correlation, and uses a **block jackknife** for standard errors.

### 1.1 Single-trait heritability (`h2`)

Under a polygenic model (Bulik-Sullivan et al., *Nat Genet* 2015):

```
E[χ²_j]  =  N · (h² / M) · ℓ_j  +  intercept
```

| Symbol | Meaning | Source |
|---|---|---|
| `χ²_j` | association chi-square of SNP *j* (`= Z_j²`) | sumstats `Z` |
| `ℓ_j` | LD score of SNP *j* (`Σ_k r²_jk`) | LD-score directory |
| `N` | per-SNP sample size | sumstats `N` |
| `M` | number of SNPs used for normalization (the heritability denominator) | LD-score count metadata; see `--count-kind` |
| `h²` | SNP-heritability — recovered from the **slope** | estimated |
| `intercept` | `1 + N·a`, where `a` captures confounding (population stratification, cryptic relatedness). Equals 1 under no confounding | estimated, or fixed |

The regression is **weighted** to correct two forms of heteroscedasticity:
(1) the variance of `χ²` grows with its mean, and (2) SNPs in LD are not
independent. The weights use a separate LD score over the regression-SNP set
(`regression_ld_scores`, the historical `w_ld`). These weights are computed
internally and are **not** user-tunable.

**Two-step estimator.** When the polygenic signal is strong, the slope and
intercept become correlated and the intercept is biased upward. The two-step
estimator first estimates the intercept using only low-`χ²` SNPs
(`χ² ≤ twostep`, default `30`), fixes it, then estimates the slope from all
retained SNPs. This is the default for single-annotation `h2`.

**Outlier cap (`chisq_max`).** A handful of very large `χ²` SNPs (major-effect
loci, the MHC) violate the infinitesimal/weighting assumptions and can leverage
the fit. Excluding `χ² > cap` mitigates this. See each command below for the
default policy.

### 1.2 Genetic correlation (`rg`)

For a pair of traits (Bulik-Sullivan et al., *Nat Genet* 2015):

```
E[Z_1j · Z_2j]  =  (√(N₁N₂) · ρ_g / M) · ℓ_j  +  intercept_gencov
```

| Symbol | Meaning |
|---|---|
| `Z_1j, Z_2j` | Z-scores of the two traits at SNP *j* |
| `ρ_g` | genetic covariance — recovered from the **slope** |
| `intercept_gencov` | sample-overlap term `≈ N_s·ρ / √(N₁N₂)`; equals 0 when the two GWAS share no samples |

The reported correlation is `rg = ρ_g / √(h²₁ · h²₂)`, so an `rg` run also fits a
single-trait `h2` regression per trait (each following Section 1.1). Only
baseline LD scores are used for `rg`.

### 1.3 Partitioned heritability (`partitioned-h2`)

For a set of (possibly overlapping) annotation categories *C*
(Finucane et al., *Nat Genet* 2015):

```
E[χ²_j]  =  N · Σ_C τ_C · ℓ(j, C)  +  intercept
```

| Symbol | Meaning |
|---|---|
| `ℓ(j, C)` | LD score of SNP *j* to category *C* (`Σ_{k∈C} r²_jk`) |
| `τ_C` | per-SNP heritability contribution of category *C* — the **coefficient** |
| `enrichment` | `(proportion of h² in C) / (proportion of SNPs in C)`, overlap-aware |

Multi-annotation models always use the "old" (legacy) regression weights, and
the two-step estimator is **not applicable** (the kernel rejects it for
`n_annot > 1`). Outlier control is therefore done with the `chisq_max` cap
instead of the two-step estimator.

### 1.4 Parameter glossary

| Parameter (CLI / config) | Backbone role | How to choose (empirical) |
|---|---|---|
| `--n-blocks` / `n_blocks` | Number of block-jackknife blocks for standard errors | Default `200` is the published setting and is almost always correct. Only worth revisiting with unusually small SNP sets (blocks become tiny) — keep the default otherwise. |
| `--count-kind` / `use_common_counts` | Selects `M`, the heritability denominator: `common` = common-variant count (`.M_5_50`, MAF ≥ 5%), `all` = full reference count (`.M`) | Use `common` (default) to match standard LDSC heritability reported in the literature. Use `all` only if your LD scores and target heritability are defined over the full SNP set. |
| `--no-intercept` / `use_intercept` | Constrains the intercept instead of estimating it (h2 → 1, gencov → 0) | Constrain only when you are confident there is no confounding/overlap (e.g. in-sample LD, no stratification); this increases power but **biases the estimate if confounding is present**. Default (estimate) is safer and lets the intercept flag confounding (intercept ≫ 1). |
| `--intercept-h2` / `intercept_h2` | Fixes the h2 intercept to a specific value | Use to pin the intercept to a value from a prior unconstrained run, or to a known stratification level. Leave unset to estimate. |
| `--intercept-gencov` / `intercept_gencov` (`rg` only) | Fixes the genetic-covariance intercept (sample-overlap term) | Set to `0` when the two GWAS share no samples to recover power; leave unset to estimate the overlap. |
| `--two-step-cutoff` / `two_step_cutoff` | Step-1 `χ²` ceiling for two-step intercept estimation | Default `30` is the published value. Lower it if moderate-`χ²` SNPs still carry polygenic signal that inflates the intercept. Applies to single-annotation `h2` and single-annotation `rg` only. |
| `--chisq-max` / `chisq_max` | Excludes SNPs with `χ²` above the cap (for `rg`, on the product `Z₁²·Z₂²`) | Set explicitly to be more aggressive about large-effect loci (e.g. when the MHC or a single locus dominates). See per-command defaults below. |
| `--samp-prev` / `samp_prev`, `--pop-prev` / `pop_prev` | Observed→liability-scale conversion inputs for binary (case-control) traits. `samp_prev` = case fraction in the GWAS sample (`P`); `pop_prev` = disease prevalence in the population (`K`). Output tables report both `*_obs` and `*_liab` columns plus the applied prevalences; the `rg` ratio is scale-invariant and unchanged. | Required only for binary traits; supply both or neither. Each is a probability in the open interval `(0, 1)`, or `nan` for a quantitative trait. For `h2` / `partitioned-h2` these are scalars; for `rg` they are comma-separated per-trait lists, or an `rg` `--prevalence-manifest` TSV (see below). Default (no prevalence) reports observed scale with `*_liab` columns `NaN`. See `tutorials/liability-scale-conversion.md` for worked examples. |
| `--prevalence-manifest` (`rg` only) | TSV of per-trait prevalences looked up by exact munged trait name, as an alternative to the positional `--samp-prev`/`--pop-prev` lists (mutually exclusive). | Use for multi-trait/glob runs where positional alignment is error-prone. Columns `trait_name`, `samp_prev`, `pop_prev` (whitespace/tab-delimited; `#` comment lines ignored). May contain extra traits (a standing repository); every resolved trait must have a row. Duplicate resolved munged names abort the run (use the positional list, or re-munge with distinct names). |
| `--allow-identity-downgrade` / `allow_identity_downgrade` | SNP-identity compatibility override (not a statistical knob) | Mechanical/compatibility only; see `io-argument-inventory.md`. |

> **Boundary convention.** Every `χ²` threshold in this module is **inclusive**
> (`χ² ≤ cutoff` / `χ² ≤ cap` / `Z₁²·Z₂² ≤ chisq_max²`), a deliberate deviation
> from legacy LDSC's strict `<`. This follows the project rule that `-max`
> arguments retain the boundary value. All other behavior matches legacy ldsc2.

---

## 2. `h2` — single-trait heritability

Fits one regression per trait over the baseline LD scores. The number of
annotations equals the number of retained baseline LD-score columns
(`n_annot`), which drives the default policy below.

| Aspect | Default behavior |
|---|---|
| Intercept | Estimated freely. |
| Two-step estimator | **On** when single-annotation (`n_annot == 1`), `--two-step-cutoff` unset, and the intercept is free → `twostep = 30`. Disabled if the intercept is fixed (`--intercept-h2` or `--no-intercept`). |
| `chisq_max` | **No cap** when single-annotation (the two-step estimator handles outliers). When multi-annotation (a baseline directory with several annotations), the partitioned default cap applies — see Section 4. An explicit `--chisq-max` always applies and is reported. |
| Counts | `M_5_50` common counts (`--count-kind common`). |
| Jackknife | `n_blocks = 200`. |

Reported diagnostics (`metadata.json`) include the **post-filter** SNP count
(`n_snps`) and the `effective_chisq_max` actually applied (`null` when uncapped).

*Mechanism:* `RegressionRunner.estimate_h2` in `regression_runner.py`.

---

## 3. `rg` — genetic correlation

Fits a single-trait `h2` regression per trait (Section 2 rules apply to those),
then a cross-trait genetic-covariance regression per pair over the baseline LD
scores.

| Aspect | Default behavior |
|---|---|
| h2 intercepts | Estimated freely; `--intercept-h2` broadcasts a fixed value to every pair. |
| gencov intercept | Estimated freely; `--intercept-gencov` fixes the sample-overlap term. |
| Two-step estimator | **On** for the rg fit when single-annotation, `--two-step-cutoff` unset, and the h2 intercept is free → `twostep = 30` (matches legacy `estimate_rg`). |
| `chisq_max` | **No default cap** (matches legacy ldsc2 `_rg`). The filter is opt-in: when `--chisq-max` is given, SNPs are dropped on the **product** `Z₁²·Z₂² ≤ chisq_max²`. |
| Counts / Jackknife | `M_5_50` common counts; `n_blocks = 200`. |

Per-pair metadata reports the post-filter `n_snps_used` / `n_blocks_used` and an
`effective_chisq_max` field. Because `rg` has no default cap, that field simply
**echoes** `config.chisq_max` (or `null`); there is no hidden default to reveal.

*Mechanism:* `RegressionRunner._fit_rg_dataset` and `estimate_rg_pairs` in
`regression_runner.py`.

---

## 4. `partitioned-h2`

Always multi-annotation by construction (overlap-aware). The two-step estimator
does not apply; outlier control is the `chisq_max` cap. Two usage regimes are
auto-detected from the LD-score directory (no user flag):

### 4.1 Functional-category regime (no query columns)

One joint fit of all baseline annotations; one output row per baseline category,
headline metric `enrichment`. Reproduces Finucane-2015 / legacy `--overlap-annot`.

### 4.2 Cell-type-specific regime (query columns present)

For each query annotation, a `baseline + one query` model is fit; one output row
per query, headline metric `coefficient` (the conditional `τ`).

### Defaults (both regimes)

| Aspect | Default behavior |
|---|---|
| Intercept | Estimated freely (`--no-intercept` / `--intercept-h2` available). |
| Two-step estimator | **Off** — not applicable to multi-annotation models (the kernel rejects it). |
| `chisq_max` | **Default cap `max(0.001 · N.max(), 80)`** when `--chisq-max` is unset, to keep extreme-`χ²` SNPs from dominating the regression. The cap is inclusive (`χ² ≤ cap`); dropped SNPs are logged as `Removed N SNPs with chi^2 > C (...)`. An explicit `--chisq-max` overrides the default. |
| Counts | `M_5_50` common counts. `N.max()` for the cap is taken **before** filtering. |
| Jackknife | `n_blocks = 200`. |

Reported diagnostics: the cell-type regime records per-query `n_snps` +
`effective_chisq_max`; the functional regime records the same pair in the
aggregate `metadata.json`. In both, `effective_chisq_max` reveals the silent
default cap when `--chisq-max` was not passed.

> The policy above assumes the genuine multi-annotation case (`n_annot > 1`),
> which holds for any real partitioning. A degenerate functional run over a
> single retained annotation collapses to the single-annotation `h2` policy
> (two-step `30`, no default cap), because the default policy keys off the
> retained annotation count, not the command name. In that case the run logs a
> `WARNING` ("...functional regime retained a single annotation...degenerate...")
> and you should supply multiple baseline annotations for a meaningful
> partitioned analysis. The cell-type regime always fits `baseline + one query`
> (`n_annot >= 2`) and so never degenerates.

*Mechanism:* `RegressionRunner.estimate_partitioned_h2_batch` in
`regression_runner.py`.

---

## 5. Internal constants (not user-exposed)

These are fixed in code and relevant to the regression configuration, but are
**not** CLI flags. They are listed so advanced users understand the full
estimator state.

| Constant | Value | Where | Role |
|---|---|---|---|
| Default two-step cutoff | `30` | `estimate_h2`, `_fit_rg_dataset` | Step-1 `χ²` ceiling when two-step is enabled and `--two-step-cutoff` is unset. |
| Default partitioned cap | `max(0.001 · N.max(), 80)` | `_resolve_default_chisq_max` | Outlier cap for multi-annotation `h2` when `--chisq-max` is unset. |
| Constrained h2 intercept | `1` | `Hsq.__null_intercept__` | Value the h2 intercept is fixed to under `--no-intercept`. |
| Constrained gencov intercept | `0` | `RG`/`Gencov.__null_intercept__` | Value the gencov intercept is fixed to under `--no-intercept`. |
| `old_weights` | `True` iff `n_annot > 1` | `estimate_h2` | Forces legacy regression weights for partitioned models (newer weights unimplemented for `n_annot > 1`). |
| Two-step mask | `χ² ≤ cutoff` (h2); `Z₁² ≤ cutoff and Z₂² ≤ cutoff` (rg) | `_kernel/regression.py` | Step-1 SNP selection (inclusive; legacy used `<`). |
| `λ_GC` divisor | `0.4549` | `_kernel/regression.py` | `median(χ²)/0.4549`, the genomic-control inflation summary. |
| Regression weights (`w_ld`) | from `regression_ld_scores` column | LD-score directory | Heteroscedasticity/LD weights; computed internally, not tunable. |

---

## 6. Summary of default policies

| Command | Default two-step | Default `chisq_max` | Default intercept | Counts |
|---|---|---|---|---|
| `h2` (single-annotation) | `30` (free intercept) | none | estimated | `M_5_50` |
| `h2` (multi-annotation baseline) | none | `max(0.001·N.max(), 80)` | estimated | `M_5_50` |
| `rg` | `30` (free h2 intercept) | none (opt-in, product form) | estimated (h2 + gencov) | `M_5_50` |
| `partitioned-h2` (functional) | none | `max(0.001·N.max(), 80)` | estimated | `M_5_50` |
| `partitioned-h2` (cell-type) | none | `max(0.001·N.max(), 80)` | estimated | `M_5_50` |

All defaults reproduce legacy ldsc2 exactly, except that every `χ²` threshold is
applied inclusively (`≤`) rather than with legacy's strict `<`.
