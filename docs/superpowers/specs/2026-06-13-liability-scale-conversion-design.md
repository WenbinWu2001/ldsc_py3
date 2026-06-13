# Observed-to-Liability Scale Conversion -- Design

**Status:** approved (brainstormed 2026-06-12/13)
**Modules affected:** `h2`, `partitioned-h2`, `rg`
**Legacy reference:** `ldsc_py2_Bulik/ldscore/regressions.py` (`h2_obs_to_liab`,
`gencov_obs_to_liab`, `Hsq.summary`, `Gencov.summary`) and
`ldscore/sumstats.py` (`estimate_h2`, `estimate_rg`, `_get_rg_table`,
`_print_gencor`).

---

## 1. Goal

Report binary-trait heritability and genetic covariance on the **liability
scale** (comparable across studies with different prevalences), wiring the
already-dead `RegressionConfig.samp_prev` / `pop_prev` fields to real CLI flags,
parsing, conversion, and reporting. The numbers must match legacy ldsc2.

Unlike legacy (which *replaces* observed with liability in a text log), the
restructured tabular output reports **both scales side by side**, plus the
prevalences actually applied.

## 2. The conversion (grounded in legacy)

For a binary trait with sample prevalence `P` and population prevalence `K`:

```
c(P, K) = K^2 * (1 - K)^2 / ( P * (1 - P) * phi(t)^2 ),   t = Phi^-1(1 - K)
```

where `phi` is the standard normal pdf and `Phi^-1` the inverse CDF
(`scipy.stats.norm.isf(K)`). Source: `regressions.py:107` `h2_obs_to_liab`.

- **Heritability** of trait *t*: `h2_liab = c(P_t, K_t) * h2_obs`; the SE scales
  by the same factor.
- **Genetic covariance** of pair *(i, j)*:
  `gencov_liab = gencov_obs * sqrt(c_i) * sqrt(c_j)` (`regressions.py:75`).
- **rg ratio** `gencov / sqrt(h2_i * h2_j)` is **scale-invariant** -- the
  `sqrt(c)` factors cancel, so `rg`, `rg_se`, `z`, `p` are never converted.
- `c = 1` (no-op) when `P` and `K` are both NaN/unset -- preserves observed
  scale.

Numeric anchor (regression test): `c(0.5, 0.01) ~= 0.5519`.

### What converts vs. what does not

| Quantity | Converts? | Why |
| --- | --- | --- |
| total h2, per-trait h2, per-category h2 (absolute) | yes (`x c`) | liability-variance units |
| genetic covariance (absolute) | yes (`x sqrt(c_i) sqrt(c_j)`) | liability-variance units |
| rg, rg_se, z, p | no | scale-invariant ratio |
| LDSC intercept, intercept_se | no | chi-square scale; legacy never scales it |
| ratio, lambda_gc, mean_chisq | no | scale-invariant |
| Prop._SNPs, Prop._h2, Enrichment (+ SE, p) | no | dimensionless ratios |
| Coefficient (tau) + se/z/p | no | legacy `_overlap_output` leaves tau observed |

## 3. Vectorizable conversion primitive

The factor `c(P, K)` is built from `norm.isf`/`norm.pdf`, which are already
numpy-vectorized. The only blocker to vectorization is scalar control flow.

**Design:** extract a standalone primitive
`liability_conversion_factor(samp_prev, pop_prev)` in
`src/ldsc/_kernel/regression.py` that:

- Accepts scalars **or** array-likes; broadcasts; scalar in -> scalar out.
- Returns `1.0` for entries where both `P` and `K` are NaN (observed scale),
  via an explicit mask (NaN would otherwise propagate through `norm.isf`).
- Validates element-wise: any non-NaN entry outside `(0, 1)` raises; a per-entry
  `P`/`K` NaN-mask that disagrees (one NaN, the other finite) raises. Error
  messages list the offending indices/values.

`h2_obs_to_liab` and `gencov_obs_to_liab` are rewritten to delegate to this
primitive, preserving their public signatures and legacy semantics.

This enables a future (out-of-scope) sensitivity module to call
`liability_conversion_factor(P_scalar, K_array)` and obtain the full
`h2_liab(K)` curve and its SE band in one broadcast. The regression CLI stays
scalar; nothing in this feature passes arrays.

**Per-trait consistency rule:** a trait's `(P, K)` must be **both finite in
`(0, 1)`** (binary case-control) or **both NaN** (quantitative). One supplied
without the other is underspecified and is a hard error.

## 4. CLI

Valid range documented in every `help=` and docstring: *"a probability in the
open interval (0, 1), or `nan` for a quantitative trait"*.

### 4.1 `h2`, `partitioned-h2` -- scalar

- `--samp-prev FLOAT` (P), `--pop-prev FLOAT` (K). `type=float` accepts `nan`.
- Both-or-neither. If both omitted -> observed scale only.
- Validated `(0, 1)`-or-NaN at parse time, before inputs load.

### 4.2 `rg` -- two mutually exclusive schemes

`rg` traits are resolved from `--sumstats-sources` (a `nargs="+"` list that may
contain globs), then disambiguated. Prevalence is supplied by exactly one of:

**(a) Positional comma-list** (legacy-compatible):

```
--samp-prev 0.5,0.5,nan  --pop-prev 0.01,0.01,nan
```

- Single comma-separated string each (`type=str`, parsed by us -- *not*
  `nargs="+"`, to avoid argparse ambiguity next to `--sumstats-sources`).
- `nan` / `na` / `none` / empty (`0.5,,0.3`) -> quantitative (None),
  case-insensitive.
- Length must **exactly** equal the resolved trait count, aligned to resolved
  order. For pair `(i, j)`: `P = [prev[i], prev[j]]`.

**(b) Prevalence manifest** (`--prevalence-manifest FILE.tsv`):

- Plain TSV; **flexible delimiter** (tab or runs of whitespace); lines whose
  first non-space char is `#` are ignored; a header row naming
  `trait_name`, `samp_prev`, `pop_prev` (order-independent).
- Lookup is **by munged trait name, exact string match** (never prefix; so
  `scz2022` and `scz2014` never collide).
- **Superset semantics:** the manifest *may* contain extra traits not in this
  run (a standing prevalence repository) -- extras are ignored. Every resolved
  trait **must** have a matching row, else a hard error listing the missing
  name(s).
- `nan` in a cell -> quantitative.

**Mutual exclusivity:** supplying both `--prevalence-manifest` and
`--samp-prev`/`--pop-prev` is a hard error. Supplying only one of
`--samp-prev`/`--pop-prev` (positional) is a hard error.

**Duplicate-name guard (manifest only):** if the resolved traits contain
duplicate munged names, the by-name manifest lookup is ambiguous. Raise a hard
error printing the duplicated names and their file paths, instructing the user
to (1) use the positional comma-list instead, or (2) re-munge with distinct
trait names (the `trait_name` metadata field of the munged sumstats). No silent
disambiguation when it would affect prevalence matching. (The positional scheme
is unaffected -- it aligns by position.)

### 4.3 Logging and provenance

Whichever scheme is used, after resolution the effective per-trait
`(samp_prev, pop_prev)` is echoed into the run log (input echo) and recorded in
`diagnostics/metadata.json`, so the log is self-contained regardless of source.

## 5. Output schema (both scales)

Convention: existing absolute columns gain an explicit `_obs` suffix; a parallel
`_liab` column is added. Liability columns hold `NaN` when no prevalence was
supplied. No separate `scale` column (the populated `_liab` + prevalence columns
are self-describing). Prevalences appear as columns in the per-trait/per-category
tables and in `rg_full.tsv`, but **not** in the concise `rg.tsv` (which carries
only scale-invariant columns).

### 5.1 `h2.tsv` and `h2_per_trait.tsv`

`summarize_total_h2` columns become:

```
trait_name, n_snps,
total_h2_obs, total_h2_obs_se, total_h2_liab, total_h2_liab_se,
samp_prev, pop_prev,
intercept, intercept_se, mean_chisq, lambda_gc, ratio, ratio_se
```

(`total_h2`/`total_h2_se` renamed to `_obs`; `_liab` + `samp_prev` + `pop_prev`
added. Intercept/ratio/lambda/mean_chisq unchanged.)

### 5.2 `partitioned_h2.tsv` (`PARTITIONED_H2_COLUMNS`)

```
Category, Prop._SNPs,
Category_h2_obs, Category_h2_obs_std_error,
Category_h2_liab, Category_h2_liab_std_error,
samp_prev, pop_prev,
Prop._h2, Prop._h2_std_error,
Enrichment, Enrichment_std_error, Enrichment_p,
Coefficient, Coefficient_std_error, Coefficient_z, Coefficient_p,
overlap_aware
```

(`Category_h2`/`Category_h2_std_error` renamed to `_obs`; `_liab` +
`samp_prev`/`pop_prev` added. Proportions, enrichment, coefficients unchanged.)

### 5.3 `rg_full.tsv` (`RG_FULL_COLUMNS`)

`h2_1`/`h2_2`/`gencov` (+ their `_se`) renamed to `_obs`; `_liab` variants and
the four prevalence columns added:

```
... n_snps_used, rg, rg_se, z, p,
h2_1_obs, h2_1_obs_se, h2_1_liab, h2_1_liab_se,
h2_2_obs, h2_2_obs_se, h2_2_liab, h2_2_liab_se,
gencov_obs, gencov_obs_se, gencov_liab, gencov_liab_se,
samp_prev_1, pop_prev_1, samp_prev_2, pop_prev_2,
intercept_h2_1, ... (unchanged) ..., pair_kind, status, error
```

### 5.4 `rg.tsv` (concise) -- unchanged

`trait_1, trait_2, n_snps_used, rg, rg_se, p, note`. No prevalence columns; rg is
scale-invariant so no information is lost.

## 6. Wiring summary

- **h2 / partitioned-h2:** the scalar `(P, K)` is parsed in `_runner_from_args`
  into `RegressionConfig.samp_prev` / `pop_prev` (floats / NaN). The summary
  functions read them.
- **rg:** the per-trait list is resolved in `run_rg_from_args` *after* sumstats
  and trait names are known (it needs them for manifest lookup + length checks),
  then passed explicitly to `RegressionRunner.estimate_rg_pairs(...,
  prevalences=...)`. It is not stored on the shared `RegressionConfig` (which is
  per-run, not per-trait).

## 7. Validation summary (all fail-fast, before LD scores load)

- Each supplied prevalence is in `(0, 1)` or NaN.
- Per trait, `P` and `K` are both finite or both NaN.
- `h2`/`partitioned`: both-or-neither scalar flags.
- `rg` positional: both lists present, equal length, length == resolved trait
  count.
- `rg` manifest: every resolved trait present (superset allowed), no duplicate
  resolved names.
- `rg`: manifest and positional are mutually exclusive.

## 8. Default behavior (no prevalence)

Numbers are identical to today's observed-scale output. The only changes to a
default run are the column **renames** (`total_h2` -> `total_h2_obs`, etc.) and
the added `_liab` (NaN) + prevalence (NaN) columns. This is an intentional,
accepted schema change for the pre-release restructure.

## 9. Out of scope

The K-sweep sensitivity-analysis output and plotting module. This design only
guarantees `liability_conversion_factor` is vectorizable so that module can be
built later without touching the kernel.
