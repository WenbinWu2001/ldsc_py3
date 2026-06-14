# Liability-Scale Heritability and Genetic Covariance

LDSC estimates heritability and genetic covariance on the **observed scale** by
default. For binary (case-control) traits, the observed-scale value depends on
the case fraction in your sample, so it is not comparable across studies. The
**liability scale** removes that dependence and is the standard scale to report.

You enable the conversion by giving LDSC two prevalences per binary trait:

- `--samp-prev` (P): the **case fraction in the GWAS sample** (e.g. 0.5 for a
  balanced case-control study).
- `--pop-prev` (K): the **disease prevalence in the population** (e.g. 0.01 for
  schizophrenia).

Each value is a probability in the open interval **(0, 1)**, or `nan` for a
**quantitative trait** (which has no liability scale). Supply both or neither.

The genetic correlation (`rg`) is scale-invariant -- it is identical on either
scale -- so conversion changes only heritability and genetic covariance. Output
tables always report both scales side by side (`*_obs` and `*_liab` columns) plus
the prevalences applied; `*_liab` columns are `NaN` when no prevalence is given.

## h2: one binary trait

```bash
ldsc h2 \
  --sumstats-file scz.sumstats.gz \
  --ldscore-dir eur_ldscores/ \
  --output-dir scz_h2/ \
  --samp-prev 0.5 \
  --pop-prev 0.01
```

`h2.tsv` reports `total_h2_obs` and `total_h2_liab` (= observed x conversion
factor), with `samp_prev` / `pop_prev` columns. Omit both flags for observed
scale only (`total_h2_liab` is then `NaN`).

## partitioned-h2: one binary trait

```bash
ldsc partitioned-h2 \
  --sumstats-file scz.sumstats.gz \
  --ldscore-dir baseline_ldscores/ \
  --output-dir scz_part/ \
  --samp-prev 0.5 \
  --pop-prev 0.01
```

`partitioned_h2.tsv` reports `category_h2_obs` and `category_h2_liab` per
category, plus the whole-model `total_h2_obs`/`total_h2_liab`. Proportions,
enrichment, and coefficients are scale-invariant.

## rg: two binary traits (ordered list)

Give a comma-separated value per trait, in the same order as
`--sumstats-sources`. The population prevalence of schizophrenia and bipolar
disorder are around 1-2% and the sample prevalence about 50% here:

```bash
ldsc rg \
  --sumstats-sources scz.sumstats.gz bip.sumstats.gz \
  --ldscore-dir eur_ldscores/ \
  --output-dir scz_bip_rg/ \
  --samp-prev 0.5,0.5 \
  --pop-prev 0.01,0.02
```

`rg_full.tsv` reports `h2_1_liab`, `h2_2_liab`, and `gencov_liab`; `rg` is
unchanged. The run log echoes the prevalence applied to each trait.

## rg: one binary trait and one quantitative trait

Mark the quantitative trait with `nan` in both lists (its heritability has no
liability scale, but the genetic covariance is still partially rescaled):

```bash
ldsc rg \
  --sumstats-sources scz.sumstats.gz height.sumstats.gz \
  --ldscore-dir eur_ldscores/ \
  --samp-prev 0.5,nan \
  --pop-prev 0.01,nan
```

## rg: many traits via a prevalence manifest

For multi-trait runs (especially with glob inputs where the resolved order is
not obvious), provide a manifest looked up by **trait name** instead of position:

`prevalences.tsv`:

```text
# standing prevalence repository; nan = quantitative trait
trait_name   samp_prev   pop_prev
scz2022      0.43        0.01
bip2021      0.41        0.02
height       nan         nan
```

```bash
ldsc rg \
  --sumstats-sources 'sumstats/*.sumstats.gz' \
  --ldscore-dir eur_ldscores/ \
  --output-dir panel_rg/ \
  --prevalence-manifest prevalences.tsv
```

Manifest notes:

- Whitespace- or tab-delimited; header columns `trait_name`, `samp_prev`,
  `pop_prev` in any order; lines starting with `#` are ignored.
- Matched by **exact** munged trait name, so `scz2022` and `scz2014` never
  collide. The munged trait name is the `trait_name` recorded when you ran
  `ldsc munge-sumstats`.
- The manifest may list **more** traits than this run uses; extras are ignored.
  Every resolved trait must have a row, or the run aborts listing the missing
  name(s).
- `--prevalence-manifest` and `--samp-prev`/`--pop-prev` are mutually exclusive.

If two inputs share the same munged trait name, the manifest cannot tell them
apart; LDSC aborts and asks you to use the ordered list, or re-munge with
distinct trait names.
