# Plotting LDSC Results

Last updated on: 2026-09-08

The plotting command turns a completed LDSC result directory into one concise exploratory figure. It chooses the plot from the result metadata, so you provide the result root rather than a table or plot type.

## Installation

Matplotlib is included in the default LDSC installation; no plotting extra is needed:

```bash
python -m pip install ldsc
```

Figures use a noninteractive backend, so the same command works in a terminal, notebook, or SLURM job without a display server. Plotting remains an explicit post-processing step and never runs automatically with regression.

## Create a plot

```bash
ldsc plot --result-dir tutorial_outputs/trait_h2
```

The figure is saved below the input result:

```text
tutorial_outputs/trait_h2/plots/
  ld_score_regression.png
  diagnostics/
    metadata.json
    plot.log
```

If the fixed plot already exists, the command stops before replacing it. Rerun with `--overwrite` only when you intend to regenerate that plot:

```bash
ldsc plot --result-dir tutorial_outputs/trait_h2 --overwrite
```

The plot command never modifies `h2.tsv`, `rg.tsv`, `partitioned_h2.tsv`, `quantile_h2.tsv`, or another source result. Conversely, overwriting a core result removes its default nested plots so an old figure cannot silently remain beside new numerical results.

## Which result directory produces which plot?

| Analysis you ran | Pass this directory to `--result-dir` | Figure | Scientific question |
| --- | --- | --- | --- |
| `ldsc h2` | Root directory containing `h2.tsv` | Binned LD Score regression diagnostic | Do observed mean chi-square statistics increase with LD Score in the pattern fitted by LDSC? |
| `ldsc rg`, all pairs | Root directory containing `rg.tsv` | Lower-triangular genetic-correlation heatmap | Which trait pairs share positive or negative genetic effects, and how uncertain is each estimate? |
| `ldsc rg --anchor-trait ...` | Root directory containing `rg.tsv` | Anchor-trait forest plot | How strongly is each partner trait genetically correlated with one focal trait? |
| Baseline-only functional `ldsc partitioned-h2` | Root directory containing `partitioned_h2.tsv` | Functional heritability-enrichment bars | Which functional categories carry more or less heritability than expected from their SNP share? |
| Query/cell-type `ldsc partitioned-h2` | Aggregate root containing the query summary `partitioned_h2.tsv` | Query-annotation evidence scatter plot | Which separately fitted queries show stronger nominal evidence after conditioning on the same baseline annotations? |
| `ldsc quantile-h2` | Root directory containing `quantile_h2.tsv` | Continuous-annotation quantile-enrichment bars | Does heritability enrichment change across increasing values of the target annotation? |

Do not pass a cell-type result's internal `diagnostics/query_annotations/<query>/` directory. Each such folder is one baseline-plus-query fit. The plot requires the aggregate root that summarizes every query.

## Reading the plots

### Binned LD Score regression diagnostic

Each point is one rank bin of SNPs. The horizontal coordinate is mean LD Score and the vertical coordinate is mean chi-square statistic. The line connects the LDSC fitted expectation saved during the original h2 fit; it is not a new trend fit. Darker point colors indicate greater relative regression weight.

Use this figure to inspect whether the binned observed pattern broadly follows the fitted LDSC relationship and whether deviations concentrate in particular LD Score ranges. It is a diagnostic, not an additional h2 estimate. Current h2 runs always save the exact bin table needed by the plot. Older runs without `diagnostics/ld_score_regression_bins.tsv` must be rerun; the plot command will not reconstruct or approximate the bins.

### Functional heritability enrichment

Each dark-grey horizontal bar is an annotation's heritability enrichment. The error bar is one block-jackknife standard error. The dashed grey line at enrichment 1 is the null expectation that an annotation's share of heritability equals its share of SNPs.

Values above 1 indicate enrichment; values below 1 indicate depletion. Overlapping annotations and correlated LD Scores mean bars are not independent. The plot intentionally does not display functional coefficients or p-values.

### Query/cell-type annotation evidence

Each point is nominal one-sided \(-\log_{10}(P)\) evidence for one query annotation's positive conditional coefficient. Each query comes from a separate fit that includes the same baseline annotations. A point farther to the right has a smaller nominal p-value.

The plot is an exploratory summary, not a ranking of comparable effect sizes. It shows no multiple-testing correction, threshold, enrichment, or raw coefficient. Apply any study-specific multiplicity procedure outside this plotting utility.

### Continuous-annotation quantile enrichment

Each horizontal bar is the heritability enrichment for one realized quantile of the target annotation. Quantiles are ordered from low to high annotation values, and bar colors become darker at higher quantiles. Bracketed values in the row labels show the realized annotation-value bounds. Error bars are one block-jackknife standard error; the dashed grey line marks enrichment 1.

Look for increasing, decreasing, or non-monotonic patterns across the ordered quantiles. Interpret each bar within the fitted joint partitioned model rather than as a marginal association.

### Genetic-correlation heatmap

Each lower-triangle cell shows genetic correlation \(r_g\), with its block-jackknife SE in parentheses. Blue indicates negative correlation and red indicates positive correlation. The upper triangle is intentionally empty because it would repeat the same trait pairs, and no color bar is shown because every available value is printed in its cell. A grey cell labeled `failed` represents a pair without an available estimate.

The figure displays estimates and uncertainty only; it does not add p-values, symbols, or multiplicity correction.

### Anchor-trait genetic-correlation plot

Each point is \(r_g\) between the anchor trait and one partner trait, and each horizontal error bar spans one block-jackknife SE. The dashed grey line marks zero genetic correlation. This view is more readable than a sparse matrix when every comparison shares one anchor.

## Explore liability-scale h2 over population prevalence

If a case-control trait's population prevalence \(K\) is uncertain, convert an existing observed-scale h2 result over a fine grid without rerunning regression:

```bash
ldsc convert-h2-scale \
  --h2-result-dir tutorial_outputs/trait_h2 \
  --samp-prev 0.50 \
  --pop-prev-range 0.01 0.20
```

The default 201-point table and figure are written to:

```text
tutorial_outputs/trait_h2/postprocessing/liability-scale/
  h2_scale_conversion.tsv
  h2_prevalence_sensitivity.png
  diagnostics/
    metadata.json
    convert-h2-scale.log
```

The curve shows liability-scale h2 across assumed population prevalence, with a ribbon of one block-jackknife SE. The prevalence values are assumptions, not estimates, and their uncertainty is not propagated. The plot therefore describes sensitivity to \(K\); it does not correct or measure uncertainty in \(K\).

For one known population prevalence, exact mode does not import Matplotlib and writes only a one-row table:

```bash
ldsc convert-h2-scale \
  --h2-result-dir tutorial_outputs/trait_h2 \
  --samp-prev 0.50 \
  --pop-prev 0.10
```

This post-processing route complements the exact shortcut on regression commands. Use `--samp-prev` and `--pop-prev` directly with `h2` or `partitioned-h2` when you already know \(K\) at fit time. For `rg`, supply one prevalence pair per input trait using its existing ordered lists or prevalence manifest. Changing \(K\) later does not require refitting; run `convert-h2-scale` on the saved h2 result.

An exact conversion with `--overwrite` removes a stale sensitivity figure from an earlier range conversion, because exact mode has no sensitivity plot.

## Python use

Python callers receive paths plus live Matplotlib objects for limited notebook customization:

```python
from ldsc import convert_h2_scale, plot_result

plot = plot_result("tutorial_outputs/trait_h2")
plot.axes.set_title("Trait A: LD Score regression diagnostic")
plot.figure.savefig(plot.path, dpi=300, bbox_inches="tight")

conversion = convert_h2_scale(
    "tutorial_outputs/trait_h2",
    samp_prev=0.50,
    pop_prev_range=(0.01, 0.20),
)
```

The CLI always uses the fixed nested destination. Python callers may pass `output_dir=` to either function, but that location is unmanaged: a later overwrite of the source result cannot find or remove it.

## Failed overwrites

If an authorized overwrite fails, inspect `RUN_FAILED.txt` in the affected output directory. It names the failed attempt, the detailed workflow log when one opened, and warns when the directory may contain incomplete or mixed artifacts. The package does not roll back a failed overwrite; files written before the failure remain, and untouched files remain as they were. A successful retry removes the marker.
