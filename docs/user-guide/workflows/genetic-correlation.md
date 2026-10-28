# Estimate genetic correlation

Last updated on: 2026-10-28

`ldsc rg` estimates the genetic correlation between traits with cross-trait LD score regression (LDSC; Bulik-Sullivan et al. 2015).
It fits every pair of the traits you supply, or each trait against one anchor trait, and also fits each trait's SNP heritability on its own.
Use it to ask which traits share genetic effects, and in which direction.

## What you need

- **Munged summary statistics for two or more traits**, written by [`ldsc munge-sumstats`](munge-sumstats.md) as `<trait>.parquet`.
  Legacy `.sumstats` or `.sumstats.gz` files with `SNP`, `A1`, `A2`, `Z` and `N` are also accepted and can be mixed with Parquet files; see [legacy summary statistics](../../reference/formats/legacy-sumstats.md).
  Every trait needs `A1` and `A2`: `rg` aligns the second trait's signed statistics to the first trait's alleles and drops strand-ambiguous or mismatched SNPs.
  Munged summary statistics written in a base mode from raw summary statistics without allele columns lack them, and `rg` stops before fitting and names those inputs; munge them again with allele columns.
- **One unpartitioned LD-score directory**, usually the one you use for [SNP heritability](heritability.md).
  `rg` uses only the baseline LD scores, so query LD scores in the directory are ignored.

All traits and the LD scores must share the SNP identifier mode and, in coordinate-family modes, the genome build.
`rg` checks the genome build of every trait and the LD scores together before fitting anything, and reports all conflicting or unresolved inputs in one error.
It has no option to override the build; supply matching inputs, as described in the [genome-build policy](../../reference/snp-identity-and-genome-build.md#genome-build-policy).
Legacy text files are matched to the LD scores by rsID instead.

Give each trait a distinct `--trait-name` when you munge it; that label names the trait in every output.

## Estimate all pairs

Munge each trait, then pass all of them to `rg`:

```bash
ldsc munge-sumstats \
  --raw-sumstats-file data/trait_1.tsv.gz \
  --trait-name trait_1 \
  --output-genome-build hg38 \
  --output-dir results/traits/trait_1

ldsc munge-sumstats \
  --raw-sumstats-file data/trait_2.tsv.gz \
  --trait-name trait_2 \
  --output-genome-build hg38 \
  --output-dir results/traits/trait_2

ldsc rg \
  --sumstats-sources results/traits/trait_1/trait_1.parquet results/traits/trait_2/trait_2.parquet \
  --ldscore-dir results/unpartitioned_ldscores \
  --output-dir results/trait_1_trait_2_rg
```

For many traits, pass a quoted glob; all unordered pairs are fitted in input order:

```bash
ldsc rg \
  --sumstats-sources "results/traits/*/*.parquet" \
  --ldscore-dir results/unpartitioned_ldscores \
  --output-dir results/panel_rg

ldsc plot --result-dir results/panel_rg
```

## Compare every trait with one anchor

```bash
ldsc rg \
  --sumstats-sources "results/traits/*/*.parquet" \
  --anchor-trait trait_1 \
  --ldscore-dir results/unpartitioned_ldscores \
  --output-dir results/trait_1_anchor_rg

ldsc plot --result-dir results/trait_1_anchor_rg
```

`--anchor-trait` takes a trait label or the path of one input; a label is matched first, then a path.
Only the pairs that include the anchor are fitted.

## Read the outputs

```text
results/panel_rg/
├── rg.tsv
├── rg_full.tsv
├── h2_per_trait.tsv
└── diagnostics/
    ├── metadata.json
    └── rg.log
```

`rg.tsv` has one row per fitted pair: `trait_1`, `trait_2`, `n_snps_used`, `rg` (the genetic correlation), `rg_se` (its block-jackknife standard error), `p` and `note`.
The p-values are nominal two-sided values; apply any multiple-testing correction yourself.
A pair whose fit fails gets `NaN` estimates and a note, and the remaining pairs are still fitted; the reason is in the `error` column of `rg_full.tsv` and in `diagnostics/rg.log`.

`rg_full.tsv` adds, for each pair, `z`, the two traits' SNP heritability (`h2_1_obs`, `h2_2_obs`) and genetic covariance (`gencov_obs`) with standard errors, their liability-scale versions, the prevalences used, the intercepts, `ratio_1`, `ratio_2`, `lambda_gc_*`, `mean_chisq_*`, `pair_kind` (`all_pairs` or `anchor`) and `status`.
The heritabilities in `rg_full.tsv` come from each pair's shared SNPs.

`h2_per_trait.tsv` has one row per trait from a separate single-trait fit on that trait's own SNPs, with the same columns as [`h2.tsv`](heritability.md#read-the-outputs).
These estimates can differ from the per-pair heritabilities in `rg_full.tsv`.

`diagnostics/metadata.json` records the trait labels, the number of pairs and the pair kind; see [artifact metadata](../../reference/formats/artifact-metadata.md#rg).
Add `--write-per-pair-detail` to also write `diagnostics/pairs/`, with a `manifest.tsv` and one folder per pair holding its `rg_full.tsv` and a `metadata.json` that records the two source files, SNP and block counts and the intercept settings.
An existing result is not replaced unless you add `--overwrite`; an overwrite without `--write-per-pair-detail` removes an old `diagnostics/pairs/` folder after the new tables are written.
See [logs and outputs](../../reference/logs-and-outputs.md).

## Plot the result

`ldsc plot --result-dir` reads the saved result without refitting and writes one figure below `plots/`:

- All pairs: `plots/rg_heatmap.png`, a lower-triangle heatmap of genetic correlations with their standard errors; the diagonal shows each trait's observed-scale SNP heritability and standard error from `h2_per_trait.tsv`.
- Anchor run: `plots/rg_anchor_forest.png`, each partner's genetic correlation with the anchor and a one-standard-error bar, with partner heritabilities in a text column and the anchor's in the subtitle.

The heritability annotations always use the observed scale, matched by `trait_name`; a missing or unusable value is shown as `failed` and does not stop the plot.
The plot also writes `plots/diagnostics/metadata.json` and `plot.log`; add `--overwrite` to redraw an existing figure.
See [plotting](plotting.md) for how to read both figures, and [troubleshooting](../troubleshooting.md#plot-rg-heritability-source-is-invalid) if plotting stops with an error.

## Inputs that share a trait label

Each input is labeled by the trait name stored at munging, or by its file name if it has none.
When several inputs share a label, `rg` relabels each of them as `label@<suffix>`, where the suffix is the shortest trailing part of its path that tells that group apart.
Parent directories are tried first, so `pgc/v3/SCZ.parquet` and `ukb/v3/SCZ.parquet` become `SCZ@pgc/v3` and `SCZ@ukb/v3`.
Inputs in the same directory get a suffix ending in the file name, such as `SCZ@SCZ_a.parquet` and `SCZ@SCZ_b.parquet`.
One suffix length is used for the whole group, so the labels do not depend on input order; the suffix comes from the path as you typed it, relative or absolute.
Python `RegressionRunner.estimate_rg_pairs()` applies the same rule to each table's `source_path`.

The new labels appear in all output tables and figures and in `diagnostics/metadata.json`; the log reports how many labels were changed.
To use such an input as the anchor, pass its new label or its path to `--anchor-trait`.
`--prevalence-manifest` matches the original munged names, so it cannot be used with inputs that share a name; use the ordered `--samp-prev` and `--pop-prev` lists instead (see [liability-scale conversion](liability-scale.md)).
If a trait name contains `@` and equals one of the generated labels, `rg` stops and asks you to munge one input again with a distinct `--trait-name`.

## Options you may change

| Option | Default and purpose |
| --- | --- |
| `--anchor-trait` | Fit only the pairs that include this trait; default: all pairs |
| `--count-kind` | `common` (default) or `all` SNP counts, as for [SNP heritability](heritability.md#options-you-may-change) |
| `--intercept-h2`, `--intercept-gencov` | Intercepts are estimated by default; `--intercept-h2 1 --intercept-gencov 0` requests the standard fixed intercepts, replacing the legacy LDSC `--no-intercept` |
| `--two-step-cutoff` | 30 by default for single-annotation fits with free intercepts; cannot be combined with fixed intercepts |
| `--chisq-max` | No cap by default; when set, keeps SNPs with $Z_1^2 Z_2^2 \le$ value$^2$ in each pair |
| `--num-blocks` | 200 block-jackknife blocks, capped at the number of SNPs |
| `--samp-prev`, `--pop-prev`, `--prevalence-manifest` | Liability-scale heritability and genetic covariance; see [liability-scale conversion](liability-scale.md) |
| `--allow-identity-downgrade` | Lets allele-aware and base-mode inputs of the same family match without alleles |
| `--write-per-pair-detail` | Also write `diagnostics/pairs/` |

The genetic-covariance intercept absorbs sample overlap between the two studies; it is near 0 when they share no samples.
Fix intercepts only when that assumption is justified for your analysis.
With a single annotation, fixing only the covariance intercept fails, because the automatic two-step estimator conflicts with it; fix the heritability intercept too.
See [regression options](../../reference/regression-options.md#3-rg--genetic-correlation), [fixed intercepts](../../reference/regression-options.md#15-requesting-fixed-intercepts) and the [command-line reference](../../reference/cli.md#ldsc-rg).

Liability-scale conversion changes the heritability and genetic-covariance columns only; the genetic correlation is the same on both scales.

For the Python API, see [Python usage](../python-usage.md#genetic-correlation).
