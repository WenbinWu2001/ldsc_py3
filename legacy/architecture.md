# Architecture

`ldsc_py3_Jerry` is a Python 3 command-line repository for LD Score workflows: computing LD Scores from PLINK genotype data, munging GWAS summary statistics into LDSC format, and running LD Score regressions for heritability, partitioned heritability, cell-type-specific analysis, and genetic correlation.

The goal of this document is to help readers quickly locate the parts of the current codebase that matter for changes, while ignoring newer experimental additions that are intentionally out of scope for this roadmap.

---

## Bird's Eye View

At a high level, the repository reads either genotype-side reference inputs (`.bed/.bim/.fam`, optional `.annot`) or summary-statistics-side inputs (`.sumstats`, `.l2.ldscore`, `.M`), normalises them into aligned SNP-level matrices, and then applies weighted linear regression with block jackknife standard errors. The main outputs are gzipped LD Score tables, `.M` count files, munged `.sumstats` files, covariance and delete-value diagnostics, and regression summary tables.

At a high level, the library implements the following features:

- **LD Score estimation**: reads PLINK genotypes, optional annotations, and window definitions, then computes per-SNP LD Scores and `.M` summaries. Entry point: `ldscore(args)`
  CLI: `./ldsc.py --bfile ... --l2 --out ...`
- **Single-trait heritability and partitioned heritability**: loads munged summary statistics with reference LD Scores and fits `Hsq`. Entry point: `ldscore.sumstats.estimate_h2()`
  CLI: `./ldsc.py --h2 ... --ref-ld ... --w-ld ... --out ...`
- **Cross-trait genetic correlation**: aligns two traits to the same SNP and allele convention, then fits two `Hsq` models plus one `Gencov` model through `RG`. Entry point: `ldscore.sumstats.estimate_rg()`
  CLI: `./ldsc.py --rg trait1.sumstats.gz,trait2.sumstats.gz --ref-ld ... --w-ld ... --out ...`
- **Cell-type-specific regression**: augments a baseline model with one cell-type-specific LD Score set at a time and writes coefficient-level results. Entry point: `ldscore.sumstats.cell_type_specific()`
  CLI: `./ldsc.py --h2-cts ... --ref-ld-chr ... --ref-ld-chr-cts ... --w-ld ... --out ...`
- **Summary-statistics munging**: normalises heterogeneous GWAS columns, applies QC filters, and writes LDSC-ready `.sumstats` files. Entry point: `munge_sumstats(args, p=True)`
  CLI: `./munge_sumstats.py --sumstats ... --out ...`
- **Annotation generation**: converts a BED file or a gene set into a binary `.annot` file aligned to a BIM file. Entry point: `make_annot_files()`
  CLI: `./make_annot.py --bed-file ... --bimfile ... --annot-file ...`

The repository is organized into 4 layers:

- CLI entry points: `ldsc.py`, `munge_sumstats.py`, `make_annot.py`
- Parsing and workflow assembly: `ldscore.parse`, `ldscore.sumstats`
- Numerical core: `ldscore.ldscore`, `ldscore.regressions`, `ldscore.irwls`, `ldscore.jackknife`
- Packaging and environment metadata: `setup.py`, `requirements.txt`, `environment.yml`

The main architectural decision is still to keep the regression modes on one shared numerical core. `Hsq`, `Gencov`, and `RG` all rely on the same regression, weighting, and jackknife machinery instead of living in separate analysis stacks.

This repository is still CLI-first rather than package-first. `ldscore/__init__.py` is empty, `setup.py` only installs the `ldscore` package plus the `ldsc.py` and `munge_sumstats.py` scripts, and most user-facing behavior is defined by command-line arguments and on-disk file formats rather than a stable import API.

---

## Code Map

This section describes the main modules and the role each one plays.

**Architecture Invariant** notes call out rules the current code relies on.

**API Boundaries** mark the parts users invoke directly.

### `ldsc.py`

Main command-line driver for LD Score estimation and all regression workflows. It defines the shared masthead, helper formatting functions, a legacy `Logger` class that `munge_sumstats.py` still imports, the `ldscore(args)` workflow for genotype-side computation, and the global CLI parser used in `__main__`. The bottom of the file now configures a standard-library `logging` logger named `LDSC`, and `ldscore.sumstats` logs through its child logger `LDSC.sumstats`, so running `ldsc.py` sets up logging for both layers. **API Boundary:** primary user-facing interface. **Architecture Invariant:** this file assembles workflows and writes outputs, but the regression math itself lives in `ldscore.regressions`, `ldscore.irwls`, and `ldscore.jackknife`.

### `munge_sumstats.py`

Preprocessing script for turning heterogeneous GWAS summary statistics into LDSC-ready `.sumstats` files. The key pieces are the header normalization maps (`default_cnames`, `describe_cname`), chunked parsing and filtering in `parse_dat()`, sample-size handling in `process_n()`, and allele alignment through `allele_merge()`. This file is not fully standalone: it imports `MASTHEAD`, `Logger`, and `sec_to_str` from `ldsc.py`, and it imports allele-validity rules from `ldscore.sumstats`. **API Boundary:** user-facing CLI for summary-statistics preparation. **Architecture Invariant:** downstream regression code expects the output schema from this script or an equivalent preprocessing step, especially `SNP`, `Z`, `N`, and optionally `A1`/`A2`.

### `make_annot.py`

Small standalone script for creating binary annotation files from either a BED file or a gene set plus gene-coordinate table. The two main functions are `gene_set_to_bed()` and `make_annot_files()`, both built around `pybedtools` and Pandas. This script prepares inputs for LD Score estimation but is not part of the regression stack. **API Boundary:** direct CLI entry point for annotation generation, although it is not listed in `setup.py`'s installed `scripts`. **Architecture Invariant:** it writes `.annot` files only; it does not call the regression or LD Score core.

### `ldscore/parse.py`

Low-level readers for LDSC-defined files and PLINK side files. Important entry points are `sumstats()`, `ldscore()`, `M()`, `annot()`, chromosome-prefix helpers such as `sub_chr()` and `get_present_chrs()`, and the `__ID_List_Factory__` used to construct file readers such as `PlinkBIMFile`, `PlinkFAMFile`, `FilterFile`, `AnnotFile`, and `ThinAnnotFile`. **Architecture Invariant:** this module owns file parsing, compression detection, and chromosome-split filename conventions; it should stay free of regression logic.

### `ldscore/sumstats.py`

Workflow assembly layer for regression analyses. It reads reference LD Scores, regression weights, `.M` files, optional annotation overlap matrices, and munged summary statistics; merges them on SNPs; checks allele compatibility for `rg`; and dispatches to `reg.Hsq`, `reg.Gencov`, or `reg.RG`. The public workflow functions are `estimate_h2()`, `estimate_rg()`, and `cell_type_specific()`. **Architecture Invariant:** this module is responsible for file-to-matrix preparation and workflow control, while the actual estimators live in `ldscore.regressions`.

### `ldscore/regressions.py`

Shared estimator layer for heritability, genetic covariance, and genetic correlation. `LD_Score_Regression` implements the common setup, coefficient scaling, intercept handling, category summaries, and jackknife-derived totals; `Hsq` adds heritability-specific weighting and summary statistics; `Gencov` adds cross-trait covariance weighting; and `RG` combines two `Hsq` fits with one `Gencov` fit to report genetic correlation. This is the main place to change model behavior, coefficient summaries, or weighting formulas. **Architecture Invariant:** inputs here are already in-memory numeric arrays with LDSC shape conventions; this module does not read files.

### `ldscore/irwls.py`

Iterative re-weighted least squares implementation used by the regression classes. `IRWLS.irwls()` alternates between weighted least-squares fits and weight updates, then hands the weighted system to the jackknife classes in `ldscore.jackknife`. The module is small, but changes here affect both single-trait and cross-trait regression workflows. **Architecture Invariant:** it operates only on matrices, weights, and callbacks; it has no knowledge of LDSC file formats.

### `ldscore/jackknife.py`

Block jackknife utilities for linear-regression estimates and ratio estimates. `Jackknife` defines separator handling and pseudovalue conversion, `LstsqJackknifeFast` is the default fast path, `LstsqJackknifeSlow` is a direct delete-block implementation used mainly for validation or special cases, and `RatioJackknife` supports derived quantities such as enrichment proportions and genetic correlation ratios. **Architecture Invariant:** everything here assumes 2D NumPy arrays and treats the first dimension as observations or blocks.

### `ldscore/ldscore.py`

Genotype-side LD Score computation code. `PlinkBEDFile` reads SNP-major PLINK `.bed` files into a bitarray-backed in-memory genotype representation, applies individual and SNP filtering, derives allele frequencies, and exposes methods used by the correlation-summing routines. `getBlockLefts()` defines sliding windows, and `__GenotypeArrayInMemory__.ldScoreVarBlocks()` drives the chunked correlation accumulation used to produce LD Scores. **Architecture Invariant:** this module is performance-sensitive and focused on genotype arrays plus annotation matrices; it does not know about summary-statistics regression workflows.

### `ldscore/__init__.py`

Currently empty. The package therefore does not define a curated Python import surface, and internal modules are imported directly. **Architecture Invariant:** adding exports here would create a new public interface that the current repository does not promise.

### `setup.py`, `requirements.txt`, and `environment.yml`

Packaging and environment metadata for the Python 3 port. `setup.py` declares the `ldscore` package and installs `ldsc.py` plus `munge_sumstats.py` as scripts, while `requirements.txt` and `environment.yml` pin the NumPy, Pandas, SciPy, and bitarray versions expected by this code. These files do not change the algorithms, but they do define the supported runtime assumptions for contributors.

---

## Cross-Cutting Concerns

### Array conventions

The numerical core consistently expects 2D NumPy arrays. Responses, weights, and sample sizes use shape `(n_snp, 1)`, while predictor matrices use `(n_snp, n_annot)`. `ldscore.jackknife`, `ldscore.irwls`, and `ldscore.regressions` all enforce these conventions, so many seemingly small reshaping mistakes surface as shape errors.

### File-format contracts

The stable interface of this repository is mostly file-based. `munge_sumstats.py` writes `.sumstats`; `ldsc.py --l2` writes `.l2.ldscore.gz`, `.l2.M`, and `.l2.M_5_50`; `make_annot.py` writes `.annot`; and `ldscore.parse` is built around those naming and column conventions. Changes to filenames, compression handling, or column layouts need coordinated updates across both the writers and the parsers.

### Logging and error handling

There is no custom exception hierarchy. Most failures surface as `ValueError`, `IOError`, or library exceptions from NumPy, Pandas, SciPy, or `pybedtools`. Logging is split across two styles: `ldsc.py` now configures standard-library logging for the main regression workflows, while `munge_sumstats.py` still uses the older file-and-stdout `Logger` class imported from `ldsc.py`.

### Performance-sensitive code

The most performance-sensitive path is genotype-side LD Score computation in `ldscore/ldscore.py`, especially the bitarray-backed genotype decoding and the chunked correlation sums inside `ldScoreVarBlocks()`. Regression fitting and jackknife work in `ldscore/regressions.py`, `ldscore/irwls.py`, and `ldscore/jackknife.py` are also runtime-critical on large datasets, so changes there should preserve vectorized NumPy-style operations where possible.

### Testing status in this scoped folder

Within the scoped `ldsc_py3_Jerry` folder, there is no in-repository `tests/` or `test/` tree to describe. That means the current architecture is defined by the scripts, package modules, and packaging metadata present here, not by an adjacent automated test layout. If you add or restructure modules, there is no local test map in this folder that will automatically show you what else needs to move.
