# User Guide

Last updated on: 2026-10-28

This guide is for scientists running analyses with LDSC3.
It explains what each analysis does, which inputs it needs, how to run it from the command line, and how to read the results.
Exact option lists, file formats and result columns are in the [reference](../reference/index.md).

## Start here

1. [Installation](installation.md): set up the package and its conda environment.
2. [Guided tutorial](guided-tutorial.md): one end-to-end analysis, from raw summary statistics to SNP heritability, genetic correlation and partitioned heritability.
3. [SNP identity and genome build](concepts/snp-identity-and-genome-build.md) and [inputs and outputs](concepts/inputs-and-outputs.md): the two conventions every command shares.
4. [Reference data](reference-data.md): which reference panels, baseline annotations and maps you need.

## Workflows

| Task | Page |
| --- | --- |
| Prepare GWAS summary statistics | [Munge summary statistics](workflows/munge-sumstats.md) |
| Build or query an LD reference | [R² matrices](workflows/r2-matrices.md) |
| Compute LD scores | [LD scores](workflows/ld-scores.md), [gene-list LD scores](workflows/gene-list-ld-scores.md) |
| Make annotation files from BED files or gene lists | [Annotations](workflows/annotations.md) |
| Estimate SNP heritability | [SNP heritability](workflows/heritability.md) |
| Estimate genetic correlation | [Genetic correlation](workflows/genetic-correlation.md) |
| Partition heritability across annotations | [Partitioned heritability](workflows/partitioned-heritability.md) |
| Test cell types or gene sets for enrichment | [Cell-type enrichment](workflows/cell-type-enrichment.md) |
| Analyze continuous annotations by quantile | [Continuous annotations](workflows/continuous-annotations.md) |
| Report binary traits on the liability scale | [Liability scale](workflows/liability-scale.md) |
| Plot a saved result | [Plotting](workflows/plotting.md) |
| Inspect Parquet and NPZ files | [Inspect artifacts](workflows/inspect-artifacts.md) |
| Move from legacy LDSC | [Migrating from legacy LDSC](workflows/migrating-from-legacy-ldsc.md) |

All Python API usage is collected in [Python usage](python-usage.md); the workflow pages show the command line.

## Concepts

- [SNP identity and genome build](concepts/snp-identity-and-genome-build.md)
- [Inputs and outputs](concepts/inputs-and-outputs.md)
- [SNP heritability from LDSC versus partitioned heritability](concepts/ldsc-vs-partitioned-h2.md)
- [Glossary](../dev/DOCS_VOCABULARY.md): the meaning of each term used in these pages

## When something goes wrong

Error messages point to the matching section of [troubleshooting](troubleshooting.md), which lists likely causes, how to check each, and remedies.
