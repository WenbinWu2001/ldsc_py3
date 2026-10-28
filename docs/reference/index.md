# Reference

Last updated on: 2026-10-27

This section is the exact reference for LDSC3: commands, the Python API, file formats, result columns and methods.
Task guides are in the [user guide](../user-guide/index.md), and the package internals are in the [developer guide](../dev/index.md).

## Workflow map

| Workflow | Command | Python entry points | Main inputs | Main outputs |
| --- | --- | --- | --- | --- |
| Inspect a genome build | `ldsc infer-genome-build` | none | One file or artifact directory | A report on stdout; no files ([genome build](snp-identity-and-genome-build.md#standalone-build-inspection)) |
| Build query annotations | `ldsc annotate` | `run_annotate()`, `AnnotationBuilder` | A BIM or annotation grid; BED files or gene lists | `query.<chrom>.annot.gz`, aligned `baseline.<chrom>.annot.gz` for an annotation grid ([`annotate`](annotate.md)) |
| Build an R² panel | `ldsc build-r2-panel` | `run_build_ref_panel()`, `ReferencePanelBuilder` | PLINK files; optional chain file, genetic maps, keep and restriction files | One Parquet file and one sidecar per chromosome ([R² panels](formats/r2-panels.md)) |
| Query R² for SNP pairs | `ldsc query-r2` | `query_r2()`, `R2Panel` | An R² panel directory and a SNP pair table | `query_r2.tsv` ([`query-r2`](query-r2.md)) |
| Build a gene LD-score index | `ldsc build-gene-ldscore-index` | `build_gene_ldscore_index()`, `load_gene_ldscore_index()`, `GeneLDScoreIndexBuildConfig` | Baseline annotations, PLINK files, a gene coordinate catalog | A complete index directory with diagnostics ([gene LD-score index](gene-ldscore-index.md)) |
| Compute LD scores | `ldsc ldscore` | `run_ldscore()`, `LDScoreCalculator`, `load_ldscore_from_dir()` returning `LDScoreSource` | A reference panel; optional baseline annotations; query annotations, BED files, gene lists or a gene LD-score index | An LD-score directory ([artifact metadata](formats/artifact-metadata.md#ldscore)) |
| Convert legacy LD scores | `ldsc convert-ldsc2-ldscores` | `convert_ldsc2_ldscores()`, `LegacyLDScoreConverter` | Complete chromosome 1 to 22 legacy reference and weight suites | An LD-score directory and `diagnostics/conversion_issues.tsv.gz` ([legacy LD scores](formats/legacy-ldscores.md)) |
| Munge summary statistics | `ldsc munge-sumstats` | `SumstatsMunger`, `infer_raw_sumstats()`, `load_sumstats()` | Raw summary statistics | Munged summary statistics as Parquet, optionally also `.sumstats.gz` ([`munge-sumstats`](munge-sumstats.md)) |
| Estimate SNP heritability | `ldsc h2` | `run_h2()`, `RegressionRunner.estimate_h2()` | Munged summary statistics and an LD-score directory | `h2.tsv` ([regression options](regression-options.md)) |
| Estimate partitioned heritability | `ldsc partitioned-h2` | `run_partitioned_h2()`, `RegressionRunner.estimate_partitioned_h2()`, `RegressionRunner.estimate_partitioned_h2_batch()` | Munged summary statistics and an LD-score directory with its overlap file | `partitioned_h2.tsv` and `partitioned_h2_full.tsv` ([result columns](results/partitioned-h2.md)) |
| Summarize heritability across annotation quantiles | `ldsc quantile-h2` | `run_quantile_h2()`, `compute_quantile_h2()`, `compute_standardized_coefficients()`, `load_fitted_partitioned_model()` | A fitted `partitioned-h2` result, its annotation sources and a target annotation | `quantile_h2.tsv` and `standardized_coefficients.tsv` ([`quantile-h2`](results/quantile-h2.md)) |
| Estimate genetic correlation | `ldsc rg` | `run_rg()`, `RegressionRunner.estimate_rg()`, `RegressionRunner.estimate_rg_pairs()` | Two or more munged summary statistics files and an LD-score directory | `rg.tsv`, `rg_full.tsv` and `h2_per_trait.tsv` ([artifact metadata](formats/artifact-metadata.md#rg)) |
| Plot a result | `ldsc plot` | `plot_result()` | A result directory from `h2`, `partitioned-h2`, `quantile-h2` or `rg` | One PNG under `<result-dir>/plots/` |
| Convert to the liability scale | `ldsc convert-h2-scale` | `convert_h2_scale()` | An `h2` result directory and prevalences | A table under `<h2-result-dir>/postprocessing/liability-scale/`, optionally a prevalence-sensitivity PNG |

## Pages in this section

| Page | Content |
| --- | --- |
| [Command-line reference](cli.md) | Every subcommand and option. |
| [Python API reference](python-api.md) | Every name exported by `ldsc`. |
| [Input paths](input-paths.md) | Exact paths, globs and chromosome suites. |
| [SNP identity and genome build](snp-identity-and-genome-build.md) | SNP identifier modes, genome-build policy and inference, liftover. |
| [LD-score SNP sets](ldscore-snp-universes.md) | Reference SNPs, regression SNPs and counts. |
| [`annotate`](annotate.md), [`munge-sumstats`](munge-sumstats.md), [gene LD-score index](gene-ldscore-index.md), [`query-r2`](query-r2.md) | Command contracts. |
| [Regression options](regression-options.md) | Estimator options and defaults for `h2`, `partitioned-h2` and `rg`. |
| [Parallel execution](parallel-execution.md) | Workers and CPU budgets. |
| [Logs and outputs](logs-and-outputs.md) | Logs, console notices, output preflight and failure markers. |
| [legacy LDSC flag map](legacy-flag-map.md) | Legacy LDSC options and their LDSC3 equivalents. |
| [Formats](formats/) | Column schema, artifact metadata, R² panels, BED and gene-list inputs, region presets, legacy files. |
| [Results](results/) | `partitioned-h2`, `quantile-h2` and gene-list diagnostic columns. |
| [Methods](methods/) | LD-score computation and the gene LD-score index mathematics. |
