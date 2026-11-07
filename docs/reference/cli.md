# Command-line reference

Last updated on: 2026-11-07

LDSC3 provides one command, `ldsc`, with one subcommand per workflow; `python -m ldsc` is equivalent.
The documentation build generates the options of each subcommand below from the parser that `ldsc.cli._command_parser()` returns, so they match `ldsc <command> --help` for the documented version.
Options appear in help order, where the help groups them into numbered steps followed by optional settings, computing resources and logging.
Task guides are in the [user guide](../user-guide/index.md).

## Conventions

- Options must be spelled in full; abbreviated options are rejected.
- Commands that write results require `--output-dir`; `plot` and `convert-h2-scale` instead write to fixed folders inside the result directory they read (`<result-dir>/plots/` and `<h2-result-dir>/postprocessing/liability-scale/`).
- `--overwrite` replaces only the files the command owns; see [logs and outputs](logs-and-outputs.md).
- `--log-level` sets the detail of the workflow log, not of the console.
- Path options accept the path tokens described in [input paths](input-paths.md).
- `--snp-identifier` and `--genome-build` follow [SNP identity and genome build](snp-identity-and-genome-build.md).

## Subcommands

| Command | Purpose |
| --- | --- |
| [`infer-genome-build`](#ldsc-infer-genome-build) | Inspect an artifact's genome build without writing files. |
| [`annotate`](#ldsc-annotate) | Find SNPs covered by BED intervals or gene boundaries. |
| [`ldscore`](#ldsc-ldscore) | Compute LD scores. |
| [`build-r2-matrix`](#ldsc-build-r2-matrix) | Build pairwise R2 matrices from PLINK genotypes. |
| [`build-gene-ldscore-index`](#ldsc-build-gene-ldscore-index) | Build a reusable gene LD-score index. |
| [`convert-legacy-ldscores`](#ldsc-convert-legacy-ldscores) | Convert selected legacy LD-score suites to LDSC3 format. |
| [`munge-sumstats`](#ldsc-munge-sumstats) | Munge GWAS summary statistics. |
| [`h2`](#ldsc-h2) | Estimate heritability from munged summary statistics and LD scores. |
| [`partitioned-h2`](#ldsc-partitioned-h2) | Estimate partitioned heritability with independent query models. |
| [`quantile-h2`](#ldsc-quantile-h2) | Summarize fitted heritability across annotation quantiles. |
| [`rg`](#ldsc-rg) | Estimate genetic correlation. |
| [`query-r2-matrix`](#ldsc-query-r2-matrix) | Look up R2 for SNP pairs in an R2 matrix. |
| [`convert-h2-scale`](#ldsc-convert-h2-scale) | Convert a saved observed-scale SNP heritability estimate to liability scale. |
| [`plot`](#ldsc-plot) | Plot a saved LDSC3 result. |

## `ldsc infer-genome-build`

Inspect an artifact's genome build without writing files.

```{argparse}
:module: ldsc_docs
:func: cli_infer_genome_build
:prog: ldsc infer-genome-build
```

## `ldsc annotate`

Find SNPs covered by BED intervals or gene boundaries.

```{argparse}
:module: ldsc_docs
:func: cli_annotate
:prog: ldsc annotate
```

## `ldsc ldscore`

Compute LD scores.

```{argparse}
:module: ldsc_docs
:func: cli_ldscore
:prog: ldsc ldscore
```

## `ldsc build-r2-matrix`

Build pairwise R2 matrices from PLINK genotypes.

```{argparse}
:module: ldsc_docs
:func: cli_build_r2_matrix
:prog: ldsc build-r2-matrix
```

## `ldsc build-gene-ldscore-index`

Build a reusable gene LD-score index.

```{argparse}
:module: ldsc_docs
:func: cli_build_gene_ldscore_index
:prog: ldsc build-gene-ldscore-index
```

## `ldsc convert-legacy-ldscores`

Convert selected legacy LD-score suites to LDSC3 format.

```{argparse}
:module: ldsc_docs
:func: cli_convert_legacy_ldscores
:prog: ldsc convert-legacy-ldscores
```

## `ldsc munge-sumstats`

Munge GWAS summary statistics.

```{argparse}
:module: ldsc_docs
:func: cli_munge_sumstats
:prog: ldsc munge-sumstats
```

## `ldsc h2`

Estimate heritability from munged summary statistics and LD scores.

```{argparse}
:module: ldsc_docs
:func: cli_h2
:prog: ldsc h2
```

## `ldsc partitioned-h2`

Estimate partitioned heritability with independent query models.

```{argparse}
:module: ldsc_docs
:func: cli_partitioned_h2
:prog: ldsc partitioned-h2
```

## `ldsc quantile-h2`

Summarize fitted heritability across annotation quantiles.

```{argparse}
:module: ldsc_docs
:func: cli_quantile_h2
:prog: ldsc quantile-h2
```

## `ldsc rg`

Estimate genetic correlation.

```{argparse}
:module: ldsc_docs
:func: cli_rg
:prog: ldsc rg
```

## `ldsc query-r2-matrix`

Look up R2 for SNP pairs in an R2 matrix.

```{argparse}
:module: ldsc_docs
:func: cli_query_r2_matrix
:prog: ldsc query-r2-matrix
```

## `ldsc convert-h2-scale`

Convert a saved observed-scale SNP heritability estimate to liability scale.

```{argparse}
:module: ldsc_docs
:func: cli_convert_h2_scale
:prog: ldsc convert-h2-scale
```

## `ldsc plot`

Plot a saved LDSC3 result.

```{argparse}
:module: ldsc_docs
:func: cli_plot
:prog: ldsc plot
```
