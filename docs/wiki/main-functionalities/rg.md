# Cross-trait genetic correlation

Last updated on: 2026-09-11

`ldsc rg` estimates genetic correlation for all unordered pairs of supplied traits, or for each trait against one anchor. It writes both correlation results and separately fitted single-trait heritabilities. The [cross-trait tutorial](../../../tutorials/cross-trait-genetic-correlation.md) covers input preparation and the Python API.

## Estimate and plot all trait pairs

```bash
ldsc rg \
  --sumstats-sources "tutorial_outputs/traits/*.parquet" \
  --ldscore-dir tutorial_outputs/baseline_ldscores \
  --output-dir tutorial_outputs/panel_rg

ldsc plot --result-dir tutorial_outputs/panel_rg
```

The rg result root contains `rg.tsv`, `rg_full.tsv`, `h2_per_trait.tsv`, and `diagnostics/metadata.json` plus `rg.log`. Plotting reads the metadata-declared sources and writes `plots/rg_heatmap.png`; it does not refit the data. Add `--overwrite` to the plot command to replace an existing figure.

The lower triangle shows \(r_g\) with jackknife SE underneath in parentheses. Light-gray diagonal cells show each trait's observed-scale \(h^2\) with its jackknife SE in the same two-line format. The upper triangle is empty. Correlation color indicates sign and magnitude; diagonal heritabilities have a neutral background.

## Estimate and plot anchor comparisons

```bash
ldsc rg \
  --sumstats-sources "tutorial_outputs/traits/*.parquet" \
  --anchor-trait trait_1 \
  --ldscore-dir tutorial_outputs/baseline_ldscores \
  --output-dir tutorial_outputs/trait_1_anchor_rg

ldsc plot --result-dir tutorial_outputs/trait_1_anchor_rg
```

Use a recovered trait label or resolved source path for `--anchor-trait`. The resulting `plots/rg_anchor_forest.png` shows correlation points with ± one jackknife SE. An aligned `Observed h² (SE)` column sits between partner names and the correlation axis, and the subtitle gives the anchor's observed heritability and SE. The heritability column is text and has no numerical axis.

## Interpret heritability annotations

Both plots match `trait_name` to the saved single-trait `total_h2_obs` and `total_h2_obs_se` columns in `h2_per_trait.tsv`. These estimates use each trait's eligible SNPs; they can differ from heritabilities estimated within rg pair fits, which use each pair's retained SNPs. Liability-scale fields are not used for these annotations. The correlation ratio itself is unchanged by observed-to-liability conversion.

Every h2 annotation uses two decimals. Finite estimates outside 0–1 are displayed without clipping, and an SE of zero is valid. `failed` means no usable estimate–SE pair is available: the source or row is missing, an estimate/SE is nonnumeric or nonfinite, or an SE is negative. Missing h2 data does not prevent plotting the correlations. Malformed tables, missing required columns, duplicate trait rows, and unsafe declared paths stop plotting with an actionable error.

See the [plotting manual](../../../tutorials/plotting-results.md#genetic-correlation-heatmap), [developer annotation contract](../../current/plotting-module.md#heritability-annotations-in-rg-plots), and [troubleshooting guide](../../troubleshooting.md#plot-rg-heritability-source-is-invalid). The source implementations are `RegressionRunner.estimate_rg_pairs()` in [regression_runner.py](../../../src/ldsc/regression_runner.py), `plot_result()` in [plotting](../../../src/ldsc/plotting/__init__.py), and the rg builders in [_builders.py](../../../src/ldsc/plotting/_builders.py).
