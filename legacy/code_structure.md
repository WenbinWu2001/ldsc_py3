# Code Structure

## Scope

This document is a coding-oriented map of the current repository, centered on the code under `ldsc_py3_Jerry/` and the newer additions that now live beside the legacy LDSC core.

Goals:

- explain what each feature does
- show the real entry points, functions, classes, and data structures
- trace the input -> transform -> output flow for each workflow
- expose module boundaries and reuse points before restructuring

Non-goals:

- statistical interpretation of LDSC output
- design proposals for the refactor

---

## Quick Repository Map

```text
ldsc_py3_Jerry/
├── ldsc.py                     # legacy CLI: LD scores + h2/rg/h2-cts workflows
├── munge_sumstats.py           # legacy CLI: GWAS summary-statistics normalization
├── make_annot.py               # legacy CLI: BED/gene-set -> single-column .annot
├── ldsc_new.py                 # thin CLI wrapper for the new LD-score workflow
├── ldscore/
│   ├── __init__.py             # empty; no curated import API
│   ├── parse.py                # file readers and ID-list helpers
│   ├── ldscore.py              # legacy PLINK-backed LD-score computation kernel
│   ├── regressions.py          # LDSC estimators: Hsq, Gencov, RG
│   ├── irwls.py                # iterative re-weighted least squares
│   ├── jackknife.py            # fast/slow block jackknife and ratio jackknife
│   ├── sumstats.py             # regression workflow assembly layer
│   └── ldscore_new.py          # new LD-score engine for PLINK or parquet R2 input
├── utils/
│   └── run_bed_to_annot.py     # package-style BED -> .annot converter
├── setup.py                    # package metadata; installs ldscore + legacy scripts
├── requirements.txt            # pip dependency pins
├── environment.yml             # conda environment
├── README.md                   # project notes/status
├── architecture.md             # architecture-level roadmap
└── plans/plan_ldscore.md       # planning note, not runtime code
```

---

## Feature Inventory

| Feature | Main user entry point | Main implementation path | Main inputs | Main outputs |
| --- | --- | --- | --- | --- |
| Legacy LD-score estimation from PLINK | `python ldsc.py --bfile ... --l2 ...` | `ldsc.py -> ldscore(args) -> ldscore.ldscore.PlinkBEDFile -> ldScoreVarBlocks()` | PLINK `.bed/.bim/.fam`, optional `.annot`, optional SNP/individual filters | `.l2.ldscore.gz`, `.l2.M`, `.l2.M_5_50`, optional `.annot.gz`, log |
| Legacy heritability / partitioned h2 | `python ldsc.py --h2 ...` | `ldsc.py -> ldscore.sumstats.estimate_h2() -> reg.Hsq` | `.sumstats`, reference `.l2.ldscore`, `.M`, weight LD scores | text/log summaries, optional `.results`, `.cov`, `.delete`, `.part_delete` |
| Legacy genetic correlation | `python ldsc.py --rg trait1,trait2,...` | `ldsc.py -> ldscore.sumstats.estimate_rg() -> reg.RG` | multiple `.sumstats`, reference LD scores, weight LD scores | logged summaries, optional covariance/delete-value files |
| Legacy cell-type-specific regression | `python ldsc.py --h2-cts ...` | `ldsc.py -> ldscore.sumstats.cell_type_specific() -> reg.Hsq` | one `.sumstats`, baseline LD scores, cell-type LD-score prefixes | `.cell_type_results.txt` |
| Summary-statistics munging | `python munge_sumstats.py --sumstats ...` | `munge_sumstats.py -> munge_sumstats()` | raw GWAS summary stats, header hints, optional allele merge file | `.sumstats.gz`, `.log` |
| Legacy single-file annotation generation | `python make_annot.py ...` | `make_annot.py -> gene_set_to_bed() / make_annot_files()` | BED or gene set + gene coordinates, BIM file | one `.annot` or `.annot.gz` |
| New LD-score workflow from SNP-level annotations + PLINK | `python ldsc_new.py ... --bfile ...` | `ldsc_new.py -> ldscore.ldscore_new.main() -> run_ldscore_from_args() -> compute_chrom_from_plink()` | SNP-level baseline/query `.annot` files, PLINK reference panel, optional freq metadata | aggregated or per-chromosome `.l2.ldscore.gz`, `.w.l2.ldscore.gz`, `.l2.M`, optional `.l2.M_5_50`, `.annotation_groups.tsv` |
| New LD-score workflow from SNP-level annotations + parquet R2 | `python ldsc_new.py ... --r2-table...` | `ldscore.ldscore_new.run_ldscore_from_args() -> compute_chrom_from_parquet()` | SNP-level annotations, sorted parquet R2 tables, optional freq metadata | same output family as above |
| R2 table normalization | Python call only | `ldscore.ldscore_new.convert_r2_table_to_sorted_parquet()` | common tabular R2 table | normalized sorted parquet |
| New BED -> `.annot.gz` conversion | `python ldsc_py3_Jerry/utils/run_bed_to_annot.py ...` | `run_bed_to_annot() -> process_baseline_file()` | BED files, baseline `.annot` templates, optional restriction file | batch or per-BED chromosome `.annot.gz` files |

---

## Layer Structure

| Layer | Files | Responsibility | Notes |
| --- | --- | --- | --- |
| User-facing CLI layer | `ldsc.py`, `munge_sumstats.py`, `make_annot.py`, `ldsc_new.py`, `utils/run_bed_to_annot.py` | parse arguments, resolve file choices, log progress, write outputs | legacy code is CLI-first; new code is both CLI- and package-friendly |
| File and metadata parsing layer | `ldscore/parse.py`, parts of `munge_sumstats.py`, parts of `ldscore/ldscore_new.py`, parts of `utils/run_bed_to_annot.py` | read LDSC formats, PLINK side files, tabular metadata, BED/restriction files | duplicated parsing conventions exist in old and new code |
| Workflow assembly layer | `ldscore/sumstats.py`, `ldscore/ldscore_new.py`, `ldsc.py` | align inputs, build matrices, choose backend, dispatch compute | this is where most orchestration complexity lives |
| Numerical kernel layer | `ldscore/ldscore.py`, `ldscore/regressions.py`, `ldscore/irwls.py`, `ldscore/jackknife.py` | LD-score accumulation, WLS/IRWLS fitting, jackknife SEs | mostly in-memory NumPy/bitarray logic |
| Packaging/environment layer | `setup.py`, `requirements.txt`, `environment.yml` | runtime assumptions and install behavior | only legacy scripts are listed in `setup.py` |

---

## Entry Points

| Entry point | Type | Calls into | Purpose |
| --- | --- | --- | --- |
| `ldsc.py` | CLI script | `ldscore(args)`, `sumstats.estimate_h2()`, `sumstats.estimate_rg()`, `sumstats.cell_type_specific()` | legacy top-level LDSC driver |
| `munge_sumstats.py` | CLI script and function | `munge_sumstats(args, p=True)` | normalize GWAS summary stats |
| `make_annot.py` | CLI script | `gene_set_to_bed()`, `make_annot_files()` | simple annotation-file generator |
| `ldsc_new.py` | CLI script | `ldscore.ldscore_new.main()` | new LD-score CLI wrapper |
| `ldscore.ldscore_new.run_ldscore()` | Python API | `run_ldscore_from_args()` | package-style new LD-score API |
| `ldscore.ldscore_new.convert_r2_table_to_sorted_parquet()` | Python API | internal parquet normalizers | prepare parquet R2 input |
| `utils.run_bed_to_annot.run_bed_to_annot()` | Python API | `process_baseline_file()` | package-style BED -> `.annot.gz` converter |

---

## Key Data Structures

| Name | Defined in | Kind | Main fields / attributes | Used by |
| --- | --- | --- | --- | --- |
| `Logger` | `ldsc.py` | class | `file_handle`, `log_fh`, `log()` | legacy `munge_sumstats.py`; old-style stdout+file logging |
| `IDContainer` family from `__ID_List_Factory__` | `ldscore/parse.py` | generated class | `df`, `IDList`, `n`, `loj()` | `PlinkBIMFile`, `PlinkFAMFile`, `FilterFile`, `AnnotFile`, `ThinAnnotFile` |
| `__GenotypeArrayInMemory__` | `ldscore/ldscore.py` | base class | `geno`, `df`, `m`, `n`, `freq`, `maf`, `kept_snps`, `sqrtpq` | base abstraction for genotype-backed LD-score computations |
| `PlinkBEDFile` | `ldscore/ldscore.py` | class | inherits genotype storage/filtering; PLINK `.bed` reader | legacy LD-score kernel and new PLINK backend |
| `LD_Score_Regression` | `ldscore/regressions.py` | base estimator | `coef`, `cat`, `tot`, `prop`, `enrichment`, `intercept`, jackknife-derived arrays | parent for `Hsq`, `Gencov` |
| `Hsq` | `ldscore/regressions.py` | estimator class | adds `mean_chisq`, `lambda_gc`, `ratio` | single-trait heritability and partitioned h2 |
| `Gencov` | `ldscore/regressions.py` | estimator class | cross-trait covariance estimate + summaries | component used by `RG` |
| `RG` | `ldscore/regressions.py` | wrapper/result class | `hsq1`, `hsq2`, `gencov`, `rg_ratio`, `rg_se`, `p`, `z` | genetic correlation workflow |
| `IRWLS` | `ldscore/irwls.py` | fitter class | `est`, `jknife_est`, `jknife_se`, `delete_values`, `separators` | used by regression estimators |
| `Jackknife` / `LstsqJackknifeFast` / `LstsqJackknifeSlow` / `RatioJackknife` | `ldscore/jackknife.py` | estimator/helper classes | separators, delete values, pseudovalues, SE/covariance | regression uncertainty estimation |
| `AnnotationBundle` | `ldscore/ldscore_new.py` | dataclass | `metadata`, `annotations`, `baseline_columns`, `query_columns` | new LD-score input assembly |
| `ChromComputationResult` | `ldscore/ldscore_new.py` | dataclass | `chrom`, `metadata`, `ld_scores`, `w_ld`, `M`, `M_5_50`, column-group fields | per-chromosome new LD-score output |
| `SortedR2BlockReader` | `ldscore/ldscore_new.py` | class | parquet dataset handle, identifier index map, cached query window | parquet-backed block-matrix provider |
| `BaselineRow` | `utils/run_bed_to_annot.py` | dataclass | `chrom`, `bp`, `snp`, `cm` | baseline template row representation |
| `RestrictResource` | `utils/run_bed_to_annot.py` | dataclass | `mode`, `bed_path`, `snp_ids` | optional restriction resource for BED -> annot |

---

## Module-by-Module Code Map

### `ldsc_py3_Jerry/ldsc.py`

| Item | Kind | Role |
| --- | --- | --- |
| `sec_to_str()` | helper | convert elapsed seconds into `d:h:m:s` text |
| `_remove_dtype()` | helper | strip pandas dtype trailers from log strings |
| `_resolve_log_path()` | helper | avoid appending to an existing log file |
| `Logger` | class | legacy file+stdout logger, still imported by `munge_sumstats.py` |
| `__filter__()` | helper | read include-list files and intersect with ID containers |
| `annot_sort_key()` | helper | sort binned continuous-annotation column labels |
| `ldscore(args)` | main function | legacy LD-score estimation workflow from PLINK input |
| global `parser` | CLI parser | defines all legacy CLI flags for LD score and regression workflows |

Main dependencies:

- imports `ldscore.ldscore` as `ld`
- imports `ldscore.parse` as `ps`
- imports `ldscore.sumstats` as `sumstats`
- imports `ldscore.regressions` as `reg`

Important coding notes:

- `__main__` configures a standard `logging` logger named `LDSC`
- `ldscore(args)` expects a module-global `logger`; it is not a pure function
- output compression still shells out to `gzip` via `subprocess.call(['gzip', '-f', ...])`
- this file mixes old and new logging styles

### `ldsc_py3_Jerry/munge_sumstats.py`

| Item group | Main names | Role |
| --- | --- | --- |
| header normalization | `default_cnames`, `describe_cname`, `get_cname_map()`, `clean_header()` | map many source column names into LDSC’s expected schema |
| file/compression handling | `read_header()`, `get_compression()` | detect input format |
| QC filters | `filter_pvals()`, `filter_info()`, `filter_frq()`, `filter_alleles()` | row-level variant filtering |
| chunked parse path | `parse_dat()` | read large summary-statistics files in chunks and apply filters |
| sample-size handling | `process_n()` | infer or standardize `N` |
| validation helpers | `check_median()`, `p_to_z()` | sanity checks and z-score conversion |
| column-flag parsing | `parse_flag_cnames()` | combine CLI-provided header hints with default mapping |
| allele merge | `allele_merge()` | align to merge-alleles reference and null out mismatches |
| top-level workflow | `munge_sumstats(args, p=True)` | end-to-end munging pipeline |

Main dependencies:

- imports `MASTHEAD`, `Logger`, `sec_to_str` from `ldsc.py`
- imports `sumstats` from `ldscore` for allele-validity tables

Important coding notes:

- this module is the most coupled script in the repo: it depends on both `ldsc.py` and `ldscore.sumstats`
- output schema drives downstream `ldscore.parse.sumstats()`
- many behaviors are controlled by header inference rather than explicit schemas

### `ldsc_py3_Jerry/make_annot.py`

| Function | Role | Inputs | Outputs |
| --- | --- | --- | --- |
| `gene_set_to_bed(args)` | convert gene-set + coordinates into merged intervals | gene names, gene coordinate table, window size | `BedTool` interval set |
| `make_annot_files(args, bed_for_annot)` | intersect BIM SNP positions with annotation intervals | BIM file + BED-like intervals | one binary `.annot` column |

Important coding notes:

- requires `pybedtools`
- produces a single annotation column named `ANNOT`
- much narrower than `utils/run_bed_to_annot.py`

### `ldsc_py3_Jerry/ldscore/parse.py`

Core responsibilities:

- read LDSC-defined files: `.sumstats`, `.l2.ldscore`, `.M`, `.annot`, `.frq`
- handle chromosome-split naming conventions with `@`
- construct lightweight file-reader classes for PLINK and filter files

Most important functions/classes:

| Name | Role |
| --- | --- |
| `sub_chr()` | replace or append chromosome token |
| `which_compression()` / `get_compression()` | detect gzip/bz2/plain |
| `sumstats()` | read munged LDSC sumstats |
| `ldscore_fromlist()` | concatenate multiple LD-score files horizontally |
| `ldscore()` | read one or many `.l2.ldscore` files |
| `M()` / `M_fromlist()` | read annotation counts |
| `annot()` | compute annotation overlap matrix |
| `__ID_List_Factory__()` | create readers such as `PlinkBIMFile` and `FilterFile` |

### `ldsc_py3_Jerry/ldscore/sumstats.py`

This is the legacy regression orchestration layer.

Key responsibilities:

- read reference LD scores and regression weights
- read and merge munged summary statistics
- align alleles for `rg`
- validate LD-score matrix quality
- call `reg.Hsq`, `reg.Gencov`, `reg.RG`
- write covariance/delete-value outputs

Most important public workflows:

| Function | Purpose |
| --- | --- |
| `estimate_h2(args)` | h2 / partitioned h2 |
| `estimate_rg(args)` | cross-trait genetic correlation |
| `cell_type_specific(args)` | one cell-type model at a time |

Important helper clusters:

| Cluster | Main names |
| --- | --- |
| input readers | `_read_ref_ld()`, `_read_annot()`, `_read_M()`, `_read_w_ld()`, `_read_sumstats()` |
| merges and validation | `_merge_and_log()`, `_check_ld_condnum()`, `_check_variance()`, `_warn_length()` |
| rg-specific alignment | `_read_other_sumstats()`, `_merge_sumstats_sumstats()`, `_filter_alleles()`, `_align_alleles()`, `_rg()` |
| output writers | `_print_cov()`, `_print_delete_values()`, `_print_part_delete_values()`, `_print_rg_cov()`, `_print_rg_delete_values()` |

### `ldsc_py3_Jerry/ldscore/regressions.py`

This is the statistical core of the legacy regression suite.

Important helpers:

| Helper | Role |
| --- | --- |
| `append_intercept()` / `remove_intercept()` | manage design matrices |
| `p_z_norm()` | convert estimate + SE to p-value and z-score |
| `h2_obs_to_liab()` / `gencov_obs_to_liab()` | liability-scale conversion |
| `update_separators()` | remap jackknife separators after masking |

Main classes:

| Class | Purpose |
| --- | --- |
| `LD_Score_Regression` | shared estimator state, coefficient/category/total summaries |
| `Hsq` | single-trait heritability and partitioned heritability |
| `Gencov` | cross-trait covariance |
| `RG` | wraps `Hsq + Hsq + Gencov` into one genetic-correlation result |

### `ldsc_py3_Jerry/ldscore/irwls.py`

Single class:

| Class | Key methods | Purpose |
| --- | --- | --- |
| `IRWLS` | `irwls()`, `wls()`, `_weight()` | iterative re-weighted least squares with jackknife integration |

### `ldsc_py3_Jerry/ldscore/jackknife.py`

Main classes/functions:

| Name | Purpose |
| --- | --- |
| `_check_shape()` / `_check_shape_block()` | enforce 2D array conventions |
| `Jackknife` | base separator/pseudovalue logic |
| `LstsqJackknifeSlow` | explicit delete-block regression jackknife |
| `LstsqJackknifeFast` | fast block jackknife from block summaries |
| `RatioJackknife` | ratio estimates such as enrichment proportions and rg |

### `ldsc_py3_Jerry/ldscore/ldscore.py`

This is the legacy genotype-backed LD-score kernel.

| Name | Purpose |
| --- | --- |
| `getBlockLefts()` | derive left boundary for each SNP window |
| `block_left_to_right()` | inverse-style helper for block bounds |
| `__GenotypeArrayInMemory__` | shared genotype-array container and LD-score accumulation methods |
| `PlinkBEDFile` | PLINK `.bed` reader plus filtering and allele-frequency derivation |

Important coding notes:

- uses `bitarray` for compact genotype storage
- `ldScoreVarBlocks()` is the core path used by both legacy and new PLINK-backed LD-score workflows

### `ldsc_py3_Jerry/ldscore/ldscore_new.py`

This is the new LD-score estimation engine. It is package-friendly and supports two reference-panel backends.

Main feature groups:

| Group | Main names | Purpose |
| --- | --- | --- |
| output/result types | `AnnotationBundle`, `ChromComputationResult` | standardize per-chromosome state |
| general helpers | `split_arg_list()`, `normalize_chromosome()`, `identifier_keys()`, `sort_frame_by_genomic_position()` | shared normalization |
| input resolution | `resolve_annotation_files()`, `resolve_chr_files()`, `resolve_parquet_files()`, `resolve_bfile_prefix()`, `resolve_frequency_files()` | file discovery |
| R2 normalization | `read_common_tabular_r2()`, `validate_r2_source_columns()`, `canonicalize_r2_pairs()`, `deduplicate_normalized_r2_pairs()`, `convert_r2_table_to_sorted_parquet()` | prepare parquet R2 data |
| annotation assembly | `parse_annotation_file()`, `combine_annotation_groups()` | build aligned baseline/query annotation matrices |
| metadata merge | `read_identifier_list()`, `load_regression_keys()`, `parse_frequency_metadata()`, `merge_frequency_metadata()`, `apply_maf_filter()` | bring in CM/MAF/regression masks |
| windowing | `chromosome_set_from_annotation_inputs()`, `build_window_coordinates()`, `check_whole_chromosome_window()` | per-chromosome traversal setup |
| parquet backend | `SortedR2BlockReader`, `ld_score_var_blocks_from_r2_reader()`, `compute_chrom_from_parquet()` | LD-score computation from sorted parquet R2 |
| PLINK backend | `compute_chrom_from_plink()` | reuse legacy PLINK kernel |
| output writing | `result_to_dataframe()`, `weight_result_to_dataframe()`, `write_ldscore_file()`, `write_counts()`, `write_annotation_groups()`, `aggregate_results()`, `emit_outputs()` | LDSC-compatible outputs |
| top-level workflow | `validate_args()`, `build_parser()`, `run_ldscore_from_args()`, `run_ldscore()`, `main()` | end-to-end execution |

### `ldsc_py3_Jerry/utils/run_bed_to_annot.py`

This is the newer BED -> `.annot.gz` converter. It is more reusable and more configurable than `make_annot.py`.

Main feature groups:

| Group | Main names | Purpose |
| --- | --- | --- |
| argument/config handling | `parse_args()`, `configure_logging()` | CLI setup |
| data types | `BaselineRow`, `RestrictResource` | normalized row/restriction state |
| file parsing helpers | `detect_delimiter()`, `iter_text_lines()`, `infer_column_index()`, `expand_bed_paths()` | robust text input handling |
| baseline-template loading | `list_baseline_annots()`, `read_baseline_annot()` | baseline SNP universe |
| BED normalization | `write_baseline_bed()`, `write_normalized_bed()` | convert to bedtools-friendly inputs |
| restriction logic | `build_restrict_resource()`, `load_restrict_snp_ids()`, `write_restrict_table_as_bed()`, `build_restrict_mask()` | optional SNP restriction |
| overlap computation | `compute_bed_overlap_mask()`, `validate_and_convert_intersection()`, `combine_masks()` | BED-to-baseline overlap masks |
| output writing | `write_annot_file()`, `baseline_output_name()` | emit `.annot.gz` |
| top-level workflow | `process_baseline_file()`, `run_bed_to_annot()`, `main()` | end-to-end conversion |

---

## End-to-End Workflows

## 1. Legacy LD-Score Estimation From PLINK (`ldsc.py --bfile --l2`)

### Input contract

| Input | Required | Meaning |
| --- | --- | --- |
| `--bfile` | yes | PLINK prefix; expects `.bed/.bim/.fam` |
| one of `--ld-wind-snps`, `--ld-wind-kb`, `--ld-wind-cm` | yes | LD window definition |
| `--annot` or `--extract` or `--cts-bin` | optional | annotation or SNP subset strategy |
| `--maf` | optional | MAF filter |
| `--keep` / `--extract` | optional | individual/SNP filters |
| `--print-snps` | optional | output subset only |
| `--out` | yes | output prefix |

### Main pipeline

1. Parse CLI and configure `logger` in `__main__`.
2. In `ldscore(args)`, resolve PLINK side files and create reader classes from `ldscore.parse`.
3. Load SNP table from `.bim` and individual table from `.fam`.
4. If provided, load annotations from `.annot`, `ThinAnnotFile`, or build binned continuous annotations from `--cts-bin`.
5. Instantiate `ldscore.ldscore.PlinkBEDFile`, which:
   - reads the `.bed` file into a `bitarray`
   - filters individuals
   - filters SNPs / MAF
   - computes allele frequencies and per-SNP metadata
6. Align annotation matrix to the retained SNP order.
7. Build `block_left` from SNP index, BP, or CM coordinates.
8. Call `geno_array.ldScoreVarBlocks(block_left, args.chunk_size, annot=...)`.
9. Build output DataFrame with metadata columns plus LD-score columns.
10. Optionally subset to `--print-snps`.
11. Write:
   - `<out>.l2.ldscore.gz`
   - `<out>.l2.M`
   - `<out>.l2.M_5_50`
   - optional `<out>.annot.gz` for `--cts-bin`
12. Log summary statistics and correlation matrices.

### Output contract

| Output | Meaning |
| --- | --- |
| `<out>.l2.ldscore.gz` | reference LD scores |
| `<out>.l2.M` | per-annotation SNP counts |
| `<out>.l2.M_5_50` | per-annotation common-SNP counts (`MAF > 0.05`) |
| `<out>.annot.gz` | optional binned annotation matrix from `--cts-bin` |
| log file | run metadata and summaries |

---

## 2. Summary-Statistics Munging (`munge_sumstats.py`)

### Input contract

| Input | Required | Meaning |
| --- | --- | --- |
| `--sumstats` | yes | raw summary-statistics file |
| `--out` | yes | output prefix |
| header hints such as `--snp`, `--a1`, `--p`, `--signed-sumstats` | optional | override column inference |
| `--merge-alleles` | optional | enforce allele alignment to a reference list |
| `--N`, `--N-cas`, `--N-con`, `--N-col`, `--N-cas-col`, `--N-con-col` | optional | sample-size hints |
| `--info-min`, `--maf-min`, `--n-min` | optional | QC thresholds |

### Main pipeline

1. Validate required CLI flags in `munge_sumstats(args)`.
2. Write run header via legacy `Logger`.
3. Read first line using `read_header()`.
4. Build source-column -> LDSC-column translation with `parse_flag_cnames()` and `get_cname_map()`.
5. Stream the input file in chunks through `pandas.read_csv(..., chunksize=...)`.
6. For each chunk, `parse_dat()`:
   - renames columns
   - checks numeric columns
   - optionally intersects with `--merge-alleles`
   - filters on INFO, FRQ, p-value, allele validity
   - drops missing data
7. Concatenate kept rows.
8. Call `process_n()` to infer or normalize `N`.
9. Convert the chosen signed statistic into standardized `Z`.
10. Optionally call `allele_merge()` against the merge-alleles reference.
11. Write standardized `.sumstats.gz`.

### Output schema

Primary downstream columns:

- `SNP`
- `Z`
- `N`
- `A1`, `A2` when alleles are retained
- optional extra retained columns such as `FRQ`

---

## 3. Legacy Regression Workflows (`ldsc.py --h2/--rg/--h2-cts`)

### Shared regression input family

| Input family | Used by | Reader path |
| --- | --- | --- |
| munged summary statistics | h2, rg, h2-cts | `ldscore.parse.sumstats()` via `ldscore.sumstats._read_sumstats()` |
| reference LD scores | h2, rg, h2-cts | `ldscore.parse.ldscore_fromlist()` |
| annotation counts `.M` | h2, rg, h2-cts | `ldscore.parse.M_fromlist()` |
| weight LD scores | h2, rg, h2-cts | `ldscore.parse.ldscore_fromlist()` |
| annotation overlap matrix | overlap-annot path only | `ldscore.parse.annot()` |

### 3A. Heritability / Partitioned Heritability

Pipeline:

1. `ldsc.py` validates flag combinations and calls `sumstats.estimate_h2(args)`.
2. `estimate_h2()` deep-copies args and normalizes prevalence/intercept options.
3. `_read_ld_sumstats()` reads and merges:
   - `.sumstats`
   - reference LD scores
   - `.M`
   - regression-weight LD scores
4. `_check_variance()` removes zero-variance LD-score columns.
5. Optionally filter high chi-square SNPs.
6. Build 2D arrays for `chisq`, `ref_ld`, weights, and `N`.
7. Instantiate `reg.Hsq(...)`.
8. Optionally write covariance and delete-value files.
9. If `--overlap-annot`, load overlap matrix and write `<out>.results`.
10. Log summary string from `hsqhat.summary(...)`.

### 3B. Genetic Correlation

Pipeline:

1. `ldsc.py` calls `sumstats.estimate_rg(args)`.
2. Parse `--rg` into one anchor trait plus one or more target traits.
3. Read anchor trait and shared LD/weight inputs through `_read_ld_sumstats()`.
4. For each additional trait:
   - `_read_other_sumstats()` reads and merges the second trait
   - `_merge_sumstats_sumstats()` renames columns to `Z1/Z2`, `N1/N2`, `A1/A1x`, `A2/A2x`
   - `_filter_alleles()` / `_align_alleles()` ensure compatible allele orientation
   - `_rg()` instantiates `reg.RG(...)`
5. Log per-pair summaries and append table rows.
6. Optionally write covariance/delete-value files for `hsq1`, `hsq2`, and `gencov`.

### 3C. Cell-Type-Specific Analysis

Pipeline:

1. `ldsc.py` calls `sumstats.cell_type_specific(args)`.
2. Read baseline LD-score inputs through `_read_ld_sumstats()`.
3. For each line in `--ref-ld-chr-cts`:
   - load the cell-type LD scores
   - merge them onto the retained SNP set
   - prepend the cell-type predictors to the baseline predictor matrix
   - read corresponding `.M` counts
   - run `reg.Hsq(...)`
4. Collect coefficient, SE, and p-value.
5. Write `<out>.cell_type_results.txt`.

---

## 4. New LD-Score Workflow (`ldsc_new.py` / `ldscore.ldscore_new`)

### Conceptual differences vs legacy `ldsc.py`

| Legacy path | New path |
| --- | --- |
| starts from genotype reference + optional `.annot` | starts from SNP-level baseline/query annotation files |
| only PLINK reference-panel backend | PLINK or sorted parquet R2 backend |
| output only reference LD scores | outputs both reference LD scores and regression-weight LD scores |
| script-oriented | script + Python API |

### Input contract

| Input group | Main flags |
| --- | --- |
| output | `--out` |
| annotation files | `--query-annot`, `--query-annot-chr`, `--baseline-annot`, `--baseline-annot-chr` |
| reference panel mode | exactly one of `--bfile/--bfile-chr` or `--r2-table/--r2-table-chr` |
| identifier mode | `--snp-identifier {chr_pos,rsID}` |
| parquet-only metadata | `--genome-build`, `--r2-bias-mode`, optional `--r2-sample-size` |
| regression-SNP selection | `--regression-snps` |
| optional MAF/CM metadata | `--frqfile`, `--frqfile-chr` |
| windowing | `--ld-wind-snps`, `--ld-wind-kb`, `--ld-wind-cm` |
| output form | `--per-chr-output` |

### Top-level pipeline

1. `run_ldscore_from_args(args)` validates flags with `validate_args()`.
2. `chromosome_set_from_annotation_inputs(args)` discovers chromosomes.
3. `load_regression_keys(args)` reads optional regression-SNP set.
4. For each chromosome:
   - resolve baseline and query annotation files
   - `combine_annotation_groups(...)` reads all files and requires identical SNP rows
   - choose backend:
     - `compute_chrom_from_parquet(...)`, or
     - `compute_chrom_from_plink(...)`
5. `emit_outputs(results, args)` writes aggregated or per-chromosome outputs.

### 4A. New PLINK backend

Pipeline inside `compute_chrom_from_plink()`:

1. Resolve PLINK prefix with `resolve_bfile_prefix()`.
2. Read `.bim` and `.fam` with legacy parse helpers.
3. Normalize panel metadata and build identifier keys.
4. Merge optional frequency metadata into annotation metadata.
5. Apply MAF filter if possible.
6. Intersect annotation SNPs with panel SNPs.
7. Create `keep_snps` index list in panel order.
8. Instantiate legacy `PlinkBEDFile`.
9. Reindex annotation matrix to the retained genotype order.
10. Build LD window coordinates and `block_left`.
11. Compute:
   - partitioned reference LD scores using the annotation matrix
   - one-column `w_ld` using `regression_mask_from_keys(...)`
12. Compute `M` and `M_5_50`.
13. Return `ChromComputationResult`.

### 4B. New parquet backend

Pipeline inside `compute_chrom_from_parquet()`:

1. Merge optional CM/MAF metadata from frequency files.
2. Apply MAF filter when available.
3. `read_sorted_r2_presence()` loads present SNP IDs or positions from parquet.
4. `filter_reference_to_present_r2()` removes annotation rows absent from the parquet reference.
5. Validate identifier uniqueness for the retained SNP set.
6. Build LD window coordinates and `block_left`.
7. Build `regression_mask` and append it as an extra annotation column.
8. Create `SortedR2BlockReader(...)`.
9. Call `ld_score_var_blocks_from_r2_reader(...)`, which mirrors the legacy sliding-block accumulation logic but queries dense block-local matrices from parquet.
10. Split the combined result into:
   - reference LD-score columns
   - single-column `w_ld`
11. Compute `M` and optional `M_5_50`.
12. Return `ChromComputationResult`.

### Output contract

| Output | Meaning |
| --- | --- |
| `<out>.l2.ldscore.gz` | aggregated reference LD-score file |
| `<out>.w.l2.ldscore.gz` | aggregated regression-weight LD-score file |
| `<out>.l2.M` | annotation counts |
| `<out>.l2.M_5_50` | common-SNP counts when MAF is available |
| `<out>.annotation_groups.tsv` | annotation -> group mapping (`baseline` or `query`) |
| `<out>.<chrom>.*` | same files in per-chromosome mode |

---

## 5. BED -> `.annot.gz` Workflows

## 5A. Legacy `make_annot.py`

Use when:

- only one annotation column is needed
- inputs are either one BED file or one gene-set file

Pipeline:

1. Read BED directly or convert gene-set + coordinates into BED via `gene_set_to_bed()`.
2. Read BIM SNP positions.
3. Convert BIM rows into 1bp intervals.
4. Intersect with annotation intervals using `pybedtools`.
5. Build one binary `ANNOT` column.
6. Write `.annot` or `.annot.gz`.

## 5B. New `utils/run_bed_to_annot.py`

Use when:

- multiple BED files should become multiple annotation columns
- you want to preserve the SNP universe and order from baseline `.annot` templates
- you need optional SNP restriction logic
- you want package-style reuse instead of a thin script

Pipeline:

1. Resolve BED inputs with `expand_bed_paths()`.
2. Resolve chromosome-specific baseline `.annot` templates with `list_baseline_annots()`.
3. Normalize BED chromosome labels with `write_normalized_bed()`.
4. Build optional restriction resource:
   - BED-like restriction for `chr_pos`
   - SNP-ID set for `rsID`
5. For each baseline template:
   - read baseline rows with `read_baseline_annot()`
   - write baseline SNPs as 1bp BED intervals
   - compute reusable restriction mask
   - compute one overlap mask per BED file
   - combine overlap and restriction masks
   - write batch or per-BED `.annot.gz`

Outputs:

- batch mode: one chromosome output with one column per BED file
- non-batch mode: one subdirectory per BED file, each with one column

---

## Reused Helpers and Shared Logic

| Reused element | Defined in | Reused by | Why it matters for refactor |
| --- | --- | --- | --- |
| `legacy_parse.sub_chr()` | `ldscore/parse.py` | `ldsc.py`, `ldscore/sumstats.py`, `ldscore/ldscore_new.py` | chromosome-split filename convention is shared across old and new code |
| PLINK readers from `__ID_List_Factory__` | `ldscore/parse.py` | `ldsc.py`, `ldscore/ldscore_new.py` | old file-reader abstractions still underpin new PLINK backend |
| `PlinkBEDFile` + `getBlockLefts()` | `ldscore/ldscore.py` | `ldsc.py`, `ldscore/ldscore_new.py` | central LD-score kernel for all PLINK-backed computation |
| regression estimators (`Hsq`, `Gencov`, `RG`) | `ldscore/regressions.py` | `ldscore/sumstats.py` | all old regression workflows share one math core |
| allele compatibility tables (`VALID_SNPS`, `MATCH_ALLELES`, `FLIP_ALLELES`) | `ldscore/sumstats.py` | `munge_sumstats.py`, `ldscore/sumstats.py` | one shared allele-logic source, but located in workflow module rather than dedicated utility module |
| legacy `Logger` | `ldsc.py` | `munge_sumstats.py` | cross-script coupling that makes logging harder to modularize |

---

## Important Input and Output Shapes

| Location | Convention |
| --- | --- |
| `ldscore/regressions.py` | predictors `x` have shape `(n_snp, n_annot)` |
| `ldscore/regressions.py` | responses `y`, weights `w`, sample sizes `N` have shape `(n_snp, 1)` |
| `ldscore/jackknife.py` | all data are 2D arrays; 1D arrays are treated as invalid |
| `ldscore/ldscore.py` | annotation matrices passed to `ldScoreVarBlocks()` have shape `(m_snp, n_annot)` |
| `ldscore/ldscore_new.py` | `ld_scores` in `ChromComputationResult` has shape `(m_snp, n_annotation_columns)` |
| `ldscore/ldscore_new.py` | `w_ld` in `ChromComputationResult` has shape `(m_snp, 1)` |

---

## Cross-Cutting Coding Constraints

| Concern | Current rule |
| --- | --- |
| Public import surface | there is effectively no curated package API because `ldscore/__init__.py` is empty |
| Logging | split across standard `logging` (`ldsc.py`, `ldscore/sumstats.py`, new modules) and legacy custom `Logger` (`munge_sumstats.py`) |
| File contracts | old code assumes LDSC-compatible filenames and column names; many functions hard-code suffix conventions |
| Compression | old parse layer supports plain/gzip/bz2 for several inputs; output paths are less consistent across modules |
| Type annotations | legacy core is largely untyped; new modules use modern type hints and dataclasses |
| Performance-sensitive code | `ldscore/ldscore.py`, `ldscore/jackknife.py`, and parquet block assembly in `ldscore/ldscore_new.py` |
| Dependency split | `setup.py` installs only the legacy scripts; newer utilities are present in the repo but not wired into packaging |

---

## Refactor-Relevant Tension Points

These are not recommendations yet; they are the places where the current structure is most uneven.

| Tension point | Evidence in code |
| --- | --- |
| Two LD-score entry stacks | legacy `ldsc.py` and new `ldscore_new.py` solve overlapping LD-score problems with different interfaces |
| Old/new logging mismatch | `munge_sumstats.py` imports `Logger` from `ldsc.py`, while the rest of the code increasingly uses `logging` |
| Parsing logic is duplicated | legacy parse helpers live in `ldscore/parse.py`, but new modules also implement table parsing, delimiter detection, alias resolution, and metadata normalization |
| Old core + new wrappers are mixed | new PLINK backend in `ldscore_new.py` still depends heavily on legacy parse/LD-score code |
| File-format knowledge is spread out | `.annot`, `.sumstats`, `.l2.ldscore`, BED, parquet R2, and frequency metadata handling are split across several modules |
| Packaging does not match code reality | `setup.py` exposes only part of the runnable codebase |

---

## Minimal “Where Do I Change X?” Guide

| If you need to change... | Start here |
| --- | --- |
| legacy LD-score calculation from PLINK | `ldsc.py`, `ldscore/ldscore.py` |
| summary-statistics parsing or QC | `munge_sumstats.py` |
| how h2/rg workflows assemble files and merges | `ldscore/sumstats.py` |
| regression math, weights, intercept handling | `ldscore/regressions.py` |
| jackknife behavior or delete values | `ldscore/jackknife.py` |
| IRWLS fitting loop | `ldscore/irwls.py` |
| new LD-score CLI/API behavior | `ldscore/ldscore_new.py` |
| parquet R2 normalization or block querying | `ldscore/ldscore_new.py` |
| BED -> `.annot.gz` workflow | `utils/run_bed_to_annot.py` |
| package/install exposure | `setup.py`, possibly new package entry points |

