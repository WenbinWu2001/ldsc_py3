# LDSC3 Beta Version 2026-08-09

Last updated: 2026-08-10

This branch is a beta release of the refactored LDSC package. The codebase and documentation are still messy and under development.

## Installation

Clone the beta branch:

```bash
git clone --branch ldsc3-beta https://github.com/WenbinWu2001/ldsc_py3.git ldsc3-beta
cd ldsc3-beta
```

Create a environment `ldsc3-beta` and install the package:

```bash
conda env create -f environment.yml -n ldsc3-beta
conda activate ldsc3-beta
python -m pip install -e ".[dev]"
```

Verify the installation. This should display the available commands and options:

```bash
ldsc --help
```

Activate the environment before later runs:

```bash
conda activate ldsc3-beta
```

To update the beta installation:

```bash
cd ldsc3-beta
git pull origin ldsc3-beta
conda activate ldsc3-beta
python -m pip install -e ".[dev]"
```

LDSC3 supports Python 3.11–3.13. The Conda environment includes the optional
dependencies needed by all workflows.

## Use

Main commands:

- `munge-sumstats`: prepare summary statistics
- `ldscore`: calculate LD scores
- `h2`: estimate SNP heritability
- `partitioned-h2`: estimate partitioned heritability
- `rg`: estimate genetic correlation

A complete list of commands is available with `ldsc --help`. To view help for a specific command, run `ldsc <command> --help`.


## Tutorials

Choose the tutorial that matches your analysis:

- **Guided analysis tutorial:** Follow the general LDSC3 workflow, from preparing summary statistics to heritability, genetic-correlation, and partitioned-heritability analyses. [Read on the GitHub Wiki](https://github.com/WenbinWu2001/ldsc_py3/wiki/Guided-tutorial) or [view the Markdown source](docs/wiki/guided-tutorial.md).
- **LDSC-SEG for protein-coding gene lists:** Test whether a protein-coding gene list is enriched for trait heritability using direct or indexed LD-score calculation. [Read on the GitHub Wiki](https://github.com/WenbinWu2001/ldsc_py3/wiki/LDSC%E2%80%90SEG-for-Protein%E2%80%90Coding-Gene-Lists) or [view the Markdown source](docs/wiki/LDSC-SEG-PC-genes.md).

## Beta feedback

When reporting a problem, include the command you ran, the full error message,
the generated log file, and the `diagnostics/` folder. To run the test suite:

```bash
pytest
```
