# Estimate SNP heritability

Last updated on: 2026-10-28

`ldsc h2` estimates the SNP heritability of one trait with LD score regression (LDSC; Bulik-Sullivan et al. 2015): it regresses the trait's $\chi^2$ statistics on LD scores and turns the slope into SNP heritability.
It writes the estimate with its block-jackknife standard error, the regression intercept, and the binned data behind a diagnostic plot.
To divide SNP heritability among annotations, use [partitioned heritability](partitioned-heritability.md) instead.

## What you need

- **Munged summary statistics** for the trait, written by [`ldsc munge-sumstats`](munge-sumstats.md): `<trait>.parquet` when munged with `--trait-name`, otherwise `sumstats.parquet`.
  Legacy `.sumstats` or `.sumstats.gz` files are also accepted; they are matched to the LD scores by rsID, as described in [legacy summary statistics](../../reference/formats/legacy-sumstats.md).
  Parquet files without LDSC3 metadata are rejected and must be munged again.
- **An unpartitioned LD-score directory**, written by [`ldsc ldscore`](ld-scores.md) without annotations, or converted from legacy files with `convert-legacy-ldscores` (see [migrating from legacy LDSC](migrating-from-legacy-ldsc.md)).
  Without annotations, `ldscore` computes LD scores for one all-ones annotation named `base`.

Both inputs must use the same SNP identifier mode and, for coordinate-family modes, the same genome build.
`h2` has no option to choose either: it reads them from the two files and stops if they disagree (`--allow-identity-downgrade` is the one exception, described under [options](#options-you-may-change)).
Choose them when you munge the summary statistics and compute the LD scores; see [SNP identity and genome build](../concepts/snp-identity-and-genome-build.md).
The regression SNPs (the HM3 map by default) are fixed when the LD scores are computed; see [reference and regression SNPs](../../reference/ldscore-snp-universes.md).

## Run the analysis

From raw summary statistics and an [R² matrix](r2-matrices.md), the full analysis is three commands and a plot:

```bash
ldsc munge-sumstats \
  --raw-sumstats-file "data/trait*.tsv.gz" \
  --trait-name trait \
  --snp-identifier chr_pos_allele_aware \
  --source-genome-build auto \
  --output-genome-build hg38 \
  --output-dir results/trait_sumstats

ldsc ldscore \
  --r2-matrix-dir /path/to/resources/r2_matrix/hg38 \
  --snp-identifier chr_pos_allele_aware \
  --ld-wind-cm 1.0 \
  --output-dir results/unpartitioned_ldscores

ldsc h2 \
  --sumstats-file results/trait_sumstats/trait.parquet \
  --ldscore-dir results/unpartitioned_ldscores \
  --output-dir results/trait_h2

ldsc plot --result-dir results/trait_h2
```

Run `munge-sumstats --infer-only` first to check how the raw file is read, and add repair flags such as `--N-col NEFF` only when it reports they are needed; see [munge summary statistics](munge-sumstats.md).
`ldscore` with an R² matrix needs the R² matrix's SNP identifier mode and a window no wider than the R² matrix's; `--workers` can split chromosomes over several workers (see [parallel execution](../../reference/parallel-execution.md)).
The same LD-score directory serves every trait analysed with that reference panel.

## Read the outputs

```text
results/trait_h2/
├── h2.tsv
└── diagnostics/
    ├── ld_score_regression_bins.tsv
    ├── metadata.json
    └── h2.log
```

`h2.tsv` has one row:

| Column | Meaning |
| --- | --- |
| `trait_name` | Trait label: `--trait-name`, else the label stored at munging, else the file name |
| `n_snps` | SNPs in the fit, after the $\chi^2$ filter |
| `total_h2_obs`, `total_h2_obs_se` | SNP heritability on the observed scale and its block-jackknife standard error: the headline result |
| `total_h2_liab`, `total_h2_liab_se` | SNP heritability on the liability scale; `NaN` unless `--samp-prev` and `--pop-prev` are given |
| `samp_prev`, `pop_prev` | Sample and population prevalence used for the liability scale, or `NaN` |
| `intercept`, `intercept_se` | Regression intercept; 1 means no confounding, and larger values indicate confounding such as population stratification; the standard error is `NaN` when the intercept is fixed |
| `mean_chisq`, `lambda_gc` | Mean $\chi^2$ and genomic-control $\lambda$ of the fitted SNPs |
| `ratio`, `ratio_se` | $(\text{intercept} - 1)/(\text{mean } \chi^2 - 1)$, the share of the inflation in mean $\chi^2$ that the intercept attributes to confounding; `NaN` when the intercept is fixed or mean $\chi^2 \le 1$ |

SNP heritability refers to the SNPs counted in the denominator: common-SNP counts by default (see [options](#options-you-may-change)).

`diagnostics/ld_score_regression_bins.tsv` summarizes the exact SNPs of the final fit in LD-score bins.
`ldsc plot --result-dir results/trait_h2` draws it as `plots/ld_score_regression.png`, a binned LD score regression diagnostic, without refitting; see [plotting](plotting.md).
`diagnostics/metadata.json` records provenance such as the SNP identifier mode used, the count vector (`count_key_used_for_regression`), the $\chi^2$ cap applied and the number of block-jackknife blocks used; its fields are listed in [artifact metadata](../../reference/formats/artifact-metadata.md#h2).

An existing result is not replaced unless you add `--overwrite`.
A successful overwrite also removes the result's old `plots/` and `postprocessing/` folders, so no figure or conversion describes an earlier fit; see [output families](../../reference/input-paths.md#coherent-output-artifact-families) and [logs and outputs](../../reference/logs-and-outputs.md).

## Liability scale for binary traits

For a case-control trait, report SNP heritability on the liability scale.
Either give `--samp-prev` (case fraction in the GWAS sample) and `--pop-prev` (population prevalence) to `h2`, or convert a saved result without refitting:

```bash
ldsc convert-h2-scale \
  --h2-result-dir results/trait_h2 \
  --samp-prev 0.5 \
  --pop-prev 0.01
```

The prevalences above are illustrative.
The conversion writes `h2_scale_conversion.tsv` under `results/trait_h2/postprocessing/liability-scale/`; there is no option to write elsewhere.
Replace `--pop-prev` with `--pop-prev-range MIN MAX` to tabulate and plot the estimate over an inclusive grid of population prevalences (`--num-points`, default 201).
See [liability-scale conversion](liability-scale.md).

## Options you may change

The defaults follow the published LDSC settings; change them only for a specific reason.

| Option | Default and purpose |
| --- | --- |
| `--count-kind` | `common` (default) scales by common-SNP counts and `all` by all-SNP counts; the requested counts must be in the LD-score directory, and the other kind is never substituted |
| `--intercept-h2` | Estimated by default; `--intercept-h2 1` fixes the intercept at the standard value (this replaces the legacy LDSC `--no-intercept`) |
| `--two-step-cutoff` | 30 by default for a single-annotation model with a free intercept; cannot be combined with `--intercept-h2` |
| `--chisq-max` | No cap for a single-annotation model; with several annotations, $\max(0.001 \cdot \max N, 80)$ |
| `--num-blocks` | 200 block-jackknife blocks, capped at the number of SNPs |
| `--allow-identity-downgrade` | Lets an allele-aware input meet a base-mode input of the same family by matching without alleles |
| `--trait-name` | Overrides the trait label in the output |
| `--samp-prev`, `--pop-prev` | Liability-scale columns; give both or neither |

`h2` fits every baseline LD-score column in the directory; with more than one, it fits a multi-annotation model, which uses the $\chi^2$ cap and no two-step estimator.
What each setting does to the estimator is explained in [regression options](../../reference/regression-options.md#2-h2--single-trait-heritability); every option is listed in the [command-line reference](../../reference/cli.md#ldsc-h2).

## Caveats

- Summary statistics and LD scores built on different genome builds or SNP identifier modes are rejected before fitting; rebuild the mismatched input rather than forcing a match.
- If no SNPs remain after the two inputs are merged, the SNP identifier mode, genome build or alleles usually differ; see [troubleshooting](../troubleshooting.md#regression-no-snps-remain-after-merging-regression-inputs).
- The diagnostic plot needs `diagnostics/ld_score_regression_bins.tsv`; older results without it must be refitted.

For the Python API, see [Python usage](../python-usage.md#snp-heritability).
