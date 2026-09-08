# Plotting and Post-processing Integration Specification

Last updated on: 2026-09-08

Status: approved; design closed; implemented and validated

## Problem and goal

LDSC3 produces structured numerical results but does not provide a small, coherent set of plots for scientists who want to explore and interpret those results. Plot generation remains an explicit post-processing action, while Matplotlib is installed as a required default dependency.

The successful implementation adds a lightweight Matplotlib-only plotting layer, a post-fit observed-to-liability-scale conversion workflow, and the upstream diagnostic data needed for a faithful binned LD Score regression plot. Every figure is derived from a canonical saved result without refitting LDSC or reconstructing scientific quantities from incomplete summaries.

This specification extends the approved observed-to-liability conversion design in [`2026-06-13-liability-scale-conversion-design.md`](2026-06-13-liability-scale-conversion-design.md). It supersedes only that document's statement that post-fit prevalence sensitivity is out of scope; its formulas, estimands, and existing regression-command prevalence semantics remain authoritative.

The detached [`ldsc_plotting_sandbox`](../../../ldsc_plotting_sandbox/README.md) is visual and behavioral reference material. Its standalone package structure, per-analysis plotting commands, required conversion destination, and optional plot destinations are prototypes and are not the integrated public interface.

## Observable behavior

### Plotting installation and execution

- A default LDSC3 installation requires Matplotlib; there is no plotting extra or opt-out mode. Failure to install a compatible Matplotlib fails normal LDSC installation.
- Plotting is invoked explicitly through `ldsc plot --result-dir RESULT_DIR`; regression workflows never generate figures automatically.
- The plot command identifies the supported analysis from canonical `diagnostics/metadata.json`; users do not select a plot type or provide an internal table path.
- The command writes one scientifically appropriate PNG below `<result-dir>/plots/`, plus plot metadata and a plot log.
- Matplotlib remains lazily imported. Ordinary imports and numerical commands do not initialize it; a damaged environment missing the required dependency fails a figure-producing path with actionable repair guidance.

### Supported result suites

| Canonical input directory | Required recorded regime | Plot | Output filename |
| --- | --- | --- | --- |
| Root `h2` result | `artifact_type=h2_result` | Binned LD Score regression diagnostic | `ld_score_regression.png` |
| Root `rg` result | `pair_kind=all_pairs` | Lower-triangular genetic-correlation heatmap | `rg_heatmap.png` |
| Root `rg` result | `pair_kind=anchor` | Anchor-trait genetic-correlation forest plot | `rg_anchor_forest.png` |
| Root `partitioned-h2` result | `analysis_type=functional_category` | Functional partitioned-heritability enrichment | `functional_h2_enrichment.png` |
| Root `partitioned-h2` result | `analysis_type=cell_type_specific` | Suite-level query-annotation nominal p-value summary | `cell_type_query_pvalues.png` |
| Root `quantile-h2` result | `artifact_type=quantile_h2_result` | Continuous-annotation quantile enrichment | `continuous_annotation_quantile_enrichment.png` |

The aggregate cell-type-specific `partitioned-h2` root is the plotting input. An individual `diagnostics/query_annotations/<query>/` directory is rejected because it represents one baseline-plus-query fit rather than the suite-level scientific question.

The command rejects raw LDSC2 output, loose TSV files, manually assembled result folders, unsupported artifact types or analysis regimes, and incomplete or contradictory canonical artifacts. It does not infer the analysis from filenames or table shape.

### Liability-scale post-processing

`ldsc convert-h2-scale --h2-result-dir H2_RESULT --samp-prev P` requires exactly one of `--pop-prev K` or `--pop-prev-range MIN MAX`.

- Exact mode writes a one-row conversion table and no plot.
- Range mode writes an inclusive, linearly spaced conversion table and `h2_prevalence_sensitivity.png`. `--num-points` defaults to 201.
- The CLI has no output-directory option. Its destination is `<h2-result-dir>/postprocessing/liability-scale/`.
- Exact mode does not import Matplotlib. Range mode uses the default dependency and fails preflight without writing a partial conversion family if the environment is damaged and Matplotlib cannot be imported.
- Every conversion starts from `total_h2_obs` and `total_h2_obs_se`; populated liability-scale fields are never conversion inputs.

The existing one-command exact-\(K\) shortcut remains unchanged: `h2` and `partitioned-h2` accept one scalar `--samp-prev`/`--pop-prev` pair, while `rg` accepts one pair per input trait through its existing positional or manifest mapping.

## Public interfaces

### Command line

The integrated command surface is:

```text
ldsc plot --result-dir RESULT_DIR [--overwrite] [--log-level LEVEL]
ldsc convert-h2-scale --h2-result-dir H2_RESULT --samp-prev P (--pop-prev K | --pop-prev-range MIN MAX) [--num-points N] [--overwrite] [--log-level LEVEL]
```

Neither command exposes `--output-dir`.

### Python

`plot_result(result_dir, output_dir=None, overwrite=False)` is the sole public plotting dispatcher. It returns an immutable `PlotArtifact` containing the plot kind, saved path, and live Matplotlib figure and axes objects.

`convert_h2_scale(..., output_dir=None, overwrite=False)` is the public post-fit conversion entry point. It returns an immutable artifact containing the table, metadata, log, and optional plot paths.

Both APIs use the CLI destination when `output_dir` is omitted. An explicit Python `output_dir=` is an advanced unmanaged destination: a later core overwrite cannot discover, validate, or remove artifacts written there.

Plot-specific builders remain internal. Core modules and the root package import surface do not import Matplotlib merely to expose non-plotting functionality.

## Source-artifact contracts

Plot dispatch requires a canonical `diagnostics/metadata.json`. The stable plotting subset is deliberately small:

- Every supported source requires `artifact_type` and the applicable source path under `files`.
- `h2` requires `files.ld_score_regression_bins`.
- `rg` requires `pair_kind` and its trait identifiers.
- `partitioned-h2` requires `analysis_type` and `headline_metric`.
- `quantile-h2` requires `target_annotation`.

The declared table must exist and contain the columns required by the selected plot. Plotting does not require or inspect `schema_version`.

There is no compatibility fallback for an older `h2` result lacking the new bin diagnostic. The plot command fails with a nonzero error that names `diagnostics/ld_score_regression_bins.tsv` and instructs the user to rerun `h2`. It never reloads original inputs or reconstructs bins from `h2.tsv`.

## Binned h2 regression diagnostic

Every successful unpartitioned `h2` run writes `diagnostics/ld_score_regression_bins.tsv` without requiring Matplotlib. The h2 metadata records it as `files.ld_score_regression_bins`.

The diagnostic uses exactly the post-chi-square-filter SNP population entering the final slope fit. Under two-step estimation, this means all SNPs used by the second-step slope fit, not only the first-step intercept subset.

SNPs are ordered by the single retained unpartitioned LD Score and divided deterministically into at most 50 nonempty, nearly equal-count rank bins. Analyses with at least 50 fitted SNPs have 50 bins; smaller analyses have one bin per SNP.

The TSV contains one row per bin:

| Column | Meaning |
| --- | --- |
| `bin` | One-based ascending LD Score bin |
| `n_snps` | Number of fitted SNPs in the bin |
| `ld_score_min`, `ld_score_max`, `mean_ld_score` | LD Score range and mean |
| `mean_chi_square`, `sd_chi_square` | Observed chi-square summary |
| `mean_sample_size` | Mean per-SNP sample size |
| `mean_fitted_chi_square` | Mean of per-SNP fitted LDSC expectations |
| `mean_regression_weight` | Mean final-model heritability regression weight |

Per-SNP fitted expectations use the final fitted intercept and slope with that SNP's exact sample size and LD Score. Diagnostic weights use the existing heritability weight function evaluated with the final fitted model and the exact fitted LD Score, regression-weight LD Score, sample size, and reference-SNP count. Computing and aggregating these diagnostics is not a second regression.

The plot places one point at saved mean LD Score and saved mean chi-square, colors points by the saved mean final-model weight, and connects the saved fitted mean chi-square values across bins. The plotting layer may normalize positive weights for the displayed relative-weight scale but may not derive fitted expectations or weights from the scalar h2 summary.

## Visual contracts

All command-generated figures are fixed 300-dpi PNGs produced with Matplotlib's noninteractive backend. Labels use plain scientific terminology supplemented by standard LDSC notation where it improves correspondence with result fields and papers.

- Binned regression: mean chi-square versus mean LD Score; fitted expectation line; points colored by relative regression weight with a compact weight scale; no p-values and no refitted trend line.
- Functional enrichment: dark-grey horizontal bars, one block-jackknife SE, and a medium-grey dashed null line at enrichment 1.
- Continuous-annotation quantiles: horizontal bars with an ordered gradient that darkens toward larger quantiles, one block-jackknife SE, and a medium-grey dashed null line at enrichment 1.
- Cell-type/query suite: horizontal scatter plot of nominal one-sided `-log10(P)`, one query per row. Each query originates from a separate fit conditional on the baseline annotations. There is no correction, significance marker, or threshold.
- All-pairs genetic correlation: one lower-triangular matrix; each available cell prints \(r_g\) and its SE; no color bar and no p-value symbols.
- Anchor genetic correlation: forest plot of \(r_g\) with one SE and a null line at zero.
- Liability sensitivity: liability-scale h2 curve with a one-SE block-jackknife ribbon; population prevalence on the horizontal axis; legend in the lower-right corner.

Null lines are drawn behind bars. Failed or unavailable estimates remain explicit rather than being silently removed.

## Derived outputs and metadata

A default plot family is:

```text
<result-dir>/plots/
  <fixed-figure-name>.png
  diagnostics/
    metadata.json
    plot.log
```

Plot metadata uses `artifact_type="plot_result"` and records the plot kind, source result location and artifact type, selected source table, produced files, and uncertainty convention.

The default conversion family is:

```text
<h2-result-dir>/postprocessing/liability-scale/
  h2_scale_conversion.tsv
  h2_prevalence_sensitivity.png  # range mode only
  diagnostics/
    metadata.json
    convert-h2-scale.log
```

Conversion metadata uses `artifact_type="h2_scale_conversion_result"` and records the source h2 result, mode, prevalence assumptions, grid, produced files, and uncertainty convention. Derived metadata does not hash inputs or copy unrestricted source metadata.

## Collision, overwrite, and failure behavior

The nested `plots/` and `postprocessing/` roots are reserved package-owned derived outputs. Users must not store irreplaceable files there.

Without `--overwrite`, plot or conversion preflight reports every conflicting path in the requested fixed artifact family before writing. An existing parent directory is not itself a conflict.

If a reserved derived subtree exists while its canonical core result files are absent, a new core run without `--overwrite` treats that orphaned subtree as a collision. Core `--overwrite` authorizes its removal through the command's ordinary stale-output cleanup.

With `--overwrite`:

- `ldsc plot` replaces only the selected fixed PNG, plot metadata, and plot log. It never alters source numerical artifacts.
- A range conversion replaces its table, metadata, log, and sensitivity plot.
- An exact conversion replaces its table, metadata, and log and removes a stale sensitivity plot from a prior range conversion.
- A successful `h2 --overwrite` replaces its bin diagnostic and removes default `plots/` and `postprocessing/` after publishing the new h2 artifacts.
- Successful `partitioned-h2`, `rg`, and `quantile-h2` overwrites remove default `plots/` after publishing their new core artifacts.

No command gains quarantine, rollback, whole-family staging, or reordered publication solely for this integration. Existing command-specific action order and failure behavior remain unchanged. Files already published by a failed attempt remain published; untouched files remain as they were.

If a core overwrite fails before its existing stale-output cleanup, nested plots or conversions remain on disk and may describe the earlier result rather than the failed new attempt. The failure marker must report that possibility plainly.

Every failed authorized overwrite in a public materializing CLI workflow or corresponding high-level Python workflow writes a durable failure marker in the applicable output scope. A no-overwrite collision that authorized no mutation does not create one.

- Ordinary workflows use `<output-dir>/RUN_FAILED.txt`.
- `ldsc plot` uses `<result-dir>/plots/RUN_FAILED.txt`.
- `convert-h2-scale` uses `<h2-result-dir>/postprocessing/liability-scale/RUN_FAILED.txt`.
- Python output overrides place the marker in the overridden root.
- Concrete chromosome-scoped `build-ref-panel` attempts use `RUN_FAILED.chr<chrom>.txt`; a full-suite attempt uses `RUN_FAILED.txt`.

The plain-text marker records the failed command, timestamp, detailed log path or the absence of a log, a conservative warning that results may be incomplete or mixed, and the recommended next action. It is not a scientific artifact, does not enter metadata, and does not block a retry. A successful materializing retry removes its applicable marker only after the command's existing success path completes.

If a failed overwrite opened its ordinary workflow log, that log remains at its canonical path with the existing failure footer and traceback. The marker does not promise restoration or rollback.

## Scientific and operational constraints

- Plots must faithfully represent package-written numerical results and must not substitute a new estimator, significance procedure, or input-alignment path.
- Error bars and ribbons represent one block-jackknife SE with inputs and prevalence assumptions held fixed. They are not confidence intervals unless explicitly labeled as such.
- Population prevalence is treated as fixed. The conversion does not propagate uncertainty in \(K\).
- Functional enrichment is the only plotted headline metric for ordinary functional partitioning.
- Query-annotation p-values are one-sided conditional coefficient tests from separate baseline-plus-query fits. Their summary plot does not imply direct effect-size comparability among queries.
- The genetic-correlation heatmap and forest plot report the scale-invariant \(r_g\); liability conversion does not alter it.
- Headless and SLURM execution must not require a display server.
- Plotting failures must not affect core result production because plotting is never automatic.

## Validation strategy

### Dependency and import isolation

- Inspect built distribution metadata and verify Matplotlib is a required dependency and no `plot` extra is published.
- Verify imports, core help, and exact-\(K\) conversion do not import Matplotlib eagerly.
- Exercise all supported dispatch paths under a headless backend.
- Simulate a damaged environment missing Matplotlib and verify that `ldsc plot` and range conversion fail before output with repair guidance.

### h2 diagnostic numerics

- Use deterministic fitted arrays to verify rank-bin membership, counts, boundaries, means, and chi-square SDs.
- Independently evaluate per-SNP fitted expectations and existing `Hsq.weights` values, then verify their bin aggregates.
- Cover free-intercept, constrained-intercept, default two-step, explicit chi-square filtering, variable sample size, tied LD Scores, fewer than 50 SNPs, and ordinary 50-bin output.
- Verify that every row contributing to the final slope appears in exactly one bin and the summed `n_snps` equals the h2 metadata count.

### Dispatch and artifacts

- Test every metadata regime against its fixed filename and required table schema.
- Verify hard failures for missing metadata fields, missing files, missing required columns, unsupported regimes, per-query directories, and old h2 results without the bin diagnostic.
- Verify plot and conversion metadata, logs, returned Python artifacts, default destinations, unmanaged Python overrides, and 300-dpi PNG creation.

### Visual semantics

- Use deterministic example tables and assert plotted values, ordering, orientation, labels, reference-line positions and z-order, error-bar magnitudes, quantile color ordering, heatmap triangle masking, absence of the heatmap color bar, and explicit unavailable cells.
- Retain example-image review as a focused visual regression check; do not use pixel-perfect snapshots as the sole numerical test.

### Conversion

- Reuse the vectorized kernel conversion primitive and its numerical anchors rather than duplicating the formula.
- Verify exact and inclusive range grids, observed-source-only conversion, scalar/array parity, SE scaling, fixed-prevalence semantics, exact-mode stale-plot removal, and range-mode dependency preflight.

### Lifecycle

- Cover no-overwrite collision preflight, command-owned replacement boundaries, core removal of default derived outputs, failed authorized overwrite markers, retained failed logs, absence of rollback, and marker cleanup after a successful retry.
- Run the focused tests, the complete pytest suite, and the standard-library unittest compatibility suite before declaring implementation complete.

## Documentation deliverables

Implementation must add:

- `docs/current/plotting-module.md`, a developer-facing description of module boundaries, optional dependency isolation, dispatch metadata, h2 bin construction, compatible result directories, outputs, lifecycle, and failure behavior.
- `tutorials/plotting-results.md`, a scientist-facing manual covering installation, commands, plot-to-analysis mapping, required result directories, and appropriate interpretation.

The README and relevant analysis tutorials must link to the user manual. Public Python objects and modules require scientific Python docstrings consistent with their implemented contracts.

## Out of scope

- Automatic or mandatory plot generation by numerical commands.
- A separate plotting distribution or one public command per plot type.
- A CLI output-directory override.
- A munged-summary-statistics Manhattan plot.
- Plotting legacy LDSC2 output or adding compatibility for old LDSC3 h2 runs.
- Plotting the legacy textual LD-score construction summary.
- Raw coefficient or p-value plots for ordinary functional partitioning.
- Coefficient, p-value, or enrichment plots for an individual baseline-plus-query result directory.
- Multiple-testing correction, corrected significance symbols, or significance thresholds.
- Publication-oriented themes, exhaustive plot customization, interactive plotting, or an exhaustive figure suite.
- Quarantine, rollback, or a new transactional publication framework.

## Risks and open questions

There are no unresolved product, scientific, or architectural decisions. Implementation must guard three result-affecting risks: retaining the exact final-slope SNP population under two-step h2, evaluating diagnostic fitted values and weights from the final model rather than an intermediate IRWLS state, and preventing eager Matplotlib imports from crossing into the core dependency path. These are validation obligations, not open interface choices.
