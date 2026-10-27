# ldsc

Last updated on: 2026-10-27

LDSC3 is a Python 3 implementation of LD score regression (LDSC).
It estimates SNP heritability, genetic correlation and partitioned heritability from GWAS summary statistics, computes LD scores from PLINK files or package-built R² panels, and builds annotations from BED files or gene lists.
It provides one command, `ldsc`, and a Python package of the same name (distribution `ldsc`, version 2.0b0).
LDSC2, the original Python 2 software, is described with its differences in [migrating from LDSC2](docs/user-guide/workflows/migrating-from-ldsc2.md).

## Install

```bash
git clone https://github.com/WenbinWu2001/ldsc_py3.git ldsc3
cd ldsc3
conda env create -f environment.yml
conda activate ldsc3
python -m pip install -e ".[dev]"
ldsc --help
```

LDSC3 supports Python 3.11 to 3.13.
The development branch, updates and optional extras are covered in [installation](docs/user-guide/installation.md).

## Quick start

Munge raw summary statistics, compute LD scores from a PLINK reference panel, and estimate SNP heritability:

```bash
ldsc munge-sumstats --raw-sumstats-file /path/to/trait.tsv.gz --output-genome-build hg19 --output-dir results/munged
ldsc ldscore --plink-prefix /path/to/1000G_EUR_hg19/1000G.EUR.@ --snp-identifier chr_pos_allele_aware --genome-build hg19 --ld-wind-cm 1 --output-dir results/ldscores
ldsc h2 --sumstats-file results/munged/sumstats.parquet --ldscore-dir results/ldscores --output-dir results/h2
```

The [guided tutorial](docs/user-guide/guided-tutorial.md) walks through a complete analysis, and `ldsc <command> --help` lists every option.

## Documentation

| Section | Content |
| --- | --- |
| [User guide](docs/user-guide/index.md) | Installation, guided tutorial, task guides, Python usage, troubleshooting. |
| [Reference](docs/reference/index.md) | Every command and option, the Python API, file formats, result columns, methods. |
| [Developer guide](docs/dev/index.md) | Code structure, architecture, design guides, writing rules, validation, release. |

Start at the [documentation home](docs/index.md); terms are defined in the [glossary](docs/dev/DOCS_VOCABULARY.md).

## Authors and maintainer

Authors: **Wenbin Wu, Anthony Abrantes, Brendan Bulik-Sullivan, and Hilary Finucane**.

Maintainer: **Wenbin Wu** ([drkwu999@gmail.com](mailto:drkwu999@gmail.com)).
Report package issues in the [issue tracker](https://github.com/WenbinWu2001/ldsc_py3/issues).

Brendan Bulik-Sullivan and Hilary Finucane developed upstream LDSC; Anthony Abrantes contributed the Python 3 port and performance work; Wenbin Wu develops and maintains this refactored package.
See [NOTICE](NOTICE) for preserved upstream attribution and dated modifications.

## License

This project continues to use the **GNU General Public License version 3 (GPLv3)**, consistent with upstream LDSC.
The complete license is in [LICENSE](LICENSE), which is included in source distributions and wheels together with [NOTICE](NOTICE).
Upstream copyright notices remain applicable to derived code; current maintainer information does not replace them.

GPLv3 permits use, modification, and commercial use.
Distributing covered derivative software requires compliance with GPL terms, including applicable licensing and corresponding-source obligations.
Merely using the package for analysis does not automatically place papers, input data, or ordinary analysis results under GPL.
See [LICENSE](LICENSE), sections 0, 2, and 4–6, for the governing terms.

The upstream declarations inspected specify GPLv3 but do not explicitly distinguish `GPL-3.0-only` from `GPL-3.0-or-later`.
We retain `GPLv3` in package metadata without assuming either SPDX identifier; the evidence is recorded in [NOTICE](NOTICE).
Bundled datasets are assessed independently in [the resource attribution record](src/ldsc/data/ATTRIBUTION.txt), including unresolved source and reuse-term questions.

## Citation

If you use this software, cite **Wenbin Wu, Anthony Abrantes, Brendan Bulik-Sullivan, and Hilary Finucane. ldsc (version 2.0b0)**, with [the repository URL](https://github.com/WenbinWu2001/ldsc_py3) and the exact version or commit used.
[CITATION.cff](CITATION.cff) provides the software citation metadata.
The installed distribution, import, and command are named `ldsc`, which is also the citation title.

Also cite the scientific methods used in your analysis:

| Analysis | Method reference |
| --- | --- |
| LDSC, SNP heritability, or regression intercept | [Bulik-Sullivan et al. (2015), basic LD Score regression](https://doi.org/10.1038/ng.3211) |
| Genetic correlation | [Bulik-Sullivan et al. (2015), genetic correlations](https://doi.org/10.1038/ng.3406) |
| Partitioned heritability | [Finucane et al. (2015), functional annotations](https://doi.org/10.1038/ng.3404) |
| Continuous annotations | [Gazal et al. (2017), LD-dependent architecture](https://doi.org/10.1038/ng.3954) |
| Cell-type-specific expression analysis | [Finucane et al. (2018), specifically expressed genes](https://doi.org/10.1038/s41588-018-0081-4) |

These academic credit instructions are separate from the software license.
Cite source datasets and any additional methods actually used as appropriate.
