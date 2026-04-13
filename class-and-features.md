# Class And Feature Design

This document defines the target class model and feature surface for the refactored `ldsc_py3_Jerry` codebase.

It is based on the current implementation in:
- `architecture.md`
- `code_structure.md`
- the legacy LDSC modules under `ldscore/`
- the newer LD-score pipeline in `ldscore/ldscore_new.py`
- the newer BED-to-annotation utility in `utils/run_bed_to_annot.py`

The goal is a codebase with:
- clear layer boundaries
- practical class responsibilities
- minimal confusion between data containers and computation logic
- simple user-facing arguments and APIs
- direct traceability from the old implementation to the refactored one

---

## 1. Recommended Design Style

Use a **hybrid model**.

### 1.1 What this means

- Use **data classes** for stable domain state, validated inputs, and results.
- Use **service classes** or **workflow modules** for expensive computation, file conversion, alignment, and orchestration.
- Keep the public API small and explicit.
- Keep low-level helpers as internal module functions, not as public methods on domain objects.

### 1.2 What not to do

- Do **not** make every concept a rich class with many methods. That usually creates large “god objects”.
- Do **not** make everything a bare dataclass either. That pushes important logic into disconnected helper files.
- Do **not** use a mutable global singleton config object shared across the whole codebase.

### 1.3 Practical rule

Use this rule for each new class:

| If the object mainly represents... | Use | Typical contents |
| --- | --- | --- |
| validated input specification | dataclass | paths, options, identifiers, thresholds |
| aligned in-memory data | dataclass | metadata table, annotation matrix, retained SNP universe |
| final output/result | dataclass | file-ready tables, counts, summary values |
| computation workflow | service class or workflow function | `run`, `load`, `compute`, `write`, `validate` |
| low-level stateless utilities | internal module functio | parsing helpers, merge helpers, formatting helpers |

---

## 2. Target Layer Structure

This should become the main structural rule for the refactored package.

| Layer | Responsibility | What belongs here | What should not live here |
| --- | --- | --- | --- |
| Public API layer | the small set of functions/classes users call directly | CLI entry points, top-level `run_*` functions, public config dataclasses | parsing details, internal matrix logic |
| Workflow/service layer | orchestration of end-to-end features | annotation builder, ref-panel loader, LD-score calculator, regression runner, sumstats munger | heavy numerical kernels, global mutable state |
| Domain model layer | validated state objects and result objects | `AnnotationBundle`, `RefPanelSpec`, `LDScoreResult`, `SumstatsTable` | unrelated helper logic |
| Adapter / IO layer | file-format specific readers/writers | PLINK reader, parquet-R2 reader, `.annot` parser, `.sumstats` parser | high-level feature decisions |
| Numerical kernel layer | performance-sensitive math | LD-score accumulation, WLS/IRWLS, jackknife, regression formulas | file reading, path handling |

### 2.1 CLI-first design principle

The design should favor **CLI usage first, Python API second**.

That means:
- every major feature should have one clear CLI entry point
- CLI arguments should map cleanly onto validated config dataclasses
- internal service classes should support the CLI, not compete with it
- Python users should still have access to the same workflows through thin wrapper functions

Recommended pattern:
- CLI parses arguments
- CLI builds config dataclasses and source/spec objects
- workflow service runs
- result object is returned and optionally written

This matches the current reality of the codebase much better than a package-first design.
### 2.2 Recommendation on exposed methods

Answer to your question: **No, not all methods should be exposed to users.**

Use three visibility levels:

| Level | Naming | Intended user |
| --- | --- | --- |
| Public | no underscore | end users and downstream developers |
| Protected/internal | single leading underscore | internal package use |
| Private implementation detail | local helper functions or nested helpers | only the module that defines them |

### 2.3 Recommendation on helper functions

Answer to your question: **store helper functions by layer, not by convenience.**

| Helper type | Best location |
| --- | --- |
| file parsing / path resolution | adapter or IO modules |
| dataframe alignment / SNP identifier normalization | domain-adjacent workflow modules |
| numerical formulas / matrix helpers | numerical kernel modules |
| formatting / logging / CLI glue | public API or CLI modules |

Do **not** place helper functions inside domain dataclasses unless they directly operate on that object’s own validated state and do not pull in unrelated dependencies.

---

## 3. Canonical SNP Identifier Rules

This must be decided early because it affects annotations, reference panels, restriction lists, LD-score outputs, and regression alignment.

### 3.1 Recommendation

Support two user-facing identifier modes only:
- `rsID`
- `chr_pos`

Internally, create a **canonical SNP identifier** column for every workflow.

### 3.2 Required columns by mode

| Identifier mode | Required columns | Matching rule | Uniqueness rule |
| --- | --- | --- | --- |
| `rsID` | `SNP` | exact string match on SNP ID | retained rows must have unique `SNP` |
| `chr_pos` | `CHR`, `BP` | exact chromosome + 1-based base-pair position | retained rows must have unique `(CHR, BP)` within chromosome |

### 3.3 Unified metadata rule

For SNP-level tables used by LD-score workflows, require the following metadata columns whenever possible:
- `CHR`
- `BP`
- `SNP`
- `CM`
- optional `MAF`

This is the most practical choice because:
- the current code already assumes these columns in multiple places
- `CM` is needed for `--ld-wind-cm`
- `SNP` is still useful even when `chr_pos` is the canonical identifier
- keeping both `SNP` and `(CHR, BP)` makes cross-format debugging much easier

### 3.4 Restriction-list rules

Your note about a unified format is correct. The practical rule should be:

| Restriction input type | Allowed formats | Interpretation |
| --- | --- | --- |
| SNP list in `rsID` mode | one-column SNP list or table with inferable SNP column | matched to `SNP` |
| SNP list in `chr_pos` mode | table with inferable `CHR` and `BP` columns, or BED-like file | matched to `(CHR, BP)` |
| BED interval restriction | allowed only in workflows that explicitly support interval projection | should not be silently treated as SNP IDs |

### 3.5 Internal recommendation

Always derive one internal SNP identifier column, for example:
- `snp_id = SNP` for `rsID`
- `snp_id = f"{CHR}:{BP}"` for `chr_pos`

All internal merges should use that canonical SNP identifier after validation.

Recommendation on naming:
- avoid `key` in public class fields and public method names
- prefer `snp_id`, `reference_snps`, `regression_snps`, or `snp_identifiers`
- if the identifier mode is `chr_pos`, those SNP identifiers are still the retained SNP labels for the workflow

**Mapping from current implementation**
- `ldscore.ldscore_new.identifier_keys()`
- `ldscore.ldscore_new.normalize_chromosome()`
- `ldscore.ldscore_new.validate_retained_identifier_uniqueness()`
- `utils.run_bed_to_annot.build_restrict_resource()`
- `utils.run_bed_to_annot.build_restrict_mask()`

---

## 4. Configuration Model

Your draft has a `config` class. The idea is right, but the design should be adjusted.

### 4.1 Recommendation

Do **not** use one mutable global configuration object.

Instead use:
- `CommonConfig`: shared options used across workflows
- one workflow-specific config per feature

Suggested config classes:
- `CommonConfig`
- `AnnotationBuildConfig`
- `RefPanelConfig`
- `LDScoreConfig`
- `MungeConfig`
- `RegressionConfig`

### 4.2 Why this is better

A single global config becomes hard to reason about because:
- many options are irrelevant for many workflows
- validation logic becomes tangled
- function signatures become less clear, not more clear

### 4.3 Recommended shared config

#### `CommonConfig` (dataclass)

| Field | Type | Purpose |
| --- | --- | --- |
| `snp_identifier` | `Literal["rsID", "chr_pos"]` | matching mode |
| `genome_build` | `Literal["hg19", "hg38"] \| None` | required for build-specific files when needed |
| `restrict_snps_path` | `str \| None` | optional global SNP restriction resource |
| `log_level` | `str` | logging verbosity |
| `fail_on_missing_metadata` | `bool` | strictness on missing `CM`, `MAF`, etc. |

#### `LDScoreConfig` (dataclass)

| Field | Type | Purpose |
| --- | --- | --- |
| `ld_wind_snps` | `int \| None` | SNP-count window |
| `ld_wind_kb` | `float \| None` | base-pair window |
| `ld_wind_cm` | `float \| None` | genetic-distance window |
| `maf_min` | `float \| None` | SNP filter threshold |
| `chunk_size` | `int` | chunking / block processing size |
| `compute_m5_50` | `bool` | whether `.M_5_50` should be emitted when possible |
| `whole_chromosome_ok` | `bool` | equivalent of `--yes-really` behavior |

#### `RegressionConfig` (dataclass)

| Field | Type | Purpose |
| --- | --- | --- |
| `n_blocks` | `int` | jackknife block count |
| `use_intercept` | `bool` | whether intercept is estimated |
| `intercept_h2` | `float \| list[float] \| None` | constrained intercept for h2 |
| `intercept_gencov` | `float \| list[float] \| None` | constrained intercept for covariance |
| `two_step_cutoff` | `float \| None` | two-step estimator cutoff |
| `chisq_max` | `float \| None` | high-chi-square filter |
| `samp_prev` | `float \| list[float] \| None` | liability conversion |
| `pop_prev` | `float \| list[float] \| None` | liability conversion |

**Mapping from current implementation**
- `ldscore.ldscore_new.build_parser()`
- legacy parser in `ldsc.py`
- parser in `munge_sumstats.py`
- parser in `utils/run_bed_to_annot.py`

---

## 5. Target Main Classes

## 5.1 `AnnotationBundle`

### Role

`AnnotationBundle` should represent **aligned SNP-level annotations already projected onto a SNP universe**.

It should **not** also be the class that stores raw BED files or gene-list files.

That is the biggest correction to the current draft.

### Recommendation

Keep `AnnotationBundle` as a dataclass.

### Suggested fields

| Field | Type | Meaning |
| --- | --- | --- |
| `metadata` | `pd.DataFrame` | SNP metadata rows for the retained annotation universe; should include `CHR/BP/SNP/CM` and optional `MAF` |
| `baseline_annotations` | `pd.DataFrame` | aligned baseline annotation columns |
| `query_annotations` | `pd.DataFrame` | aligned query annotation columns |
| `baseline_columns` | `list[str]` | ordered baseline column names |
| `query_columns` | `list[str]` | ordered query column names |
| `identifier_mode` | `str` | `rsID` or `chr_pos` |
| `chromosomes` | `list[str]` | chromosomes represented |
| `source_summary` | `dict[str, object]` | bookkeeping for traceability |

### Suggested methods

Keep methods light:
- `validate()`
- `reference_snps()`
- `has_full_baseline_cover()`
- `annotation_matrix(include_query: bool = True)`
- `summary()`

### Methods that should **not** live here

Do **not** put these on `AnnotationBundle`:
- `run_bed_to_annot`
- file discovery
- BED parsing
- path resolution
- bulk writing to disk

Those belong in a separate service.

### Supporting service

Create `AnnotationBuilder` or `AnnotationProjector`.

Suggested responsibilities:
- load baseline `.annot` files
- project BED inputs to SNP-level annotations
- combine baseline + query annotations
- enforce identical row order / identifier alignment
- apply optional restriction set

### Current-to-target mapping

| Target responsibility | Current implementation |
| --- | --- |
| aligned SNP-level annotation bundle | `ldscore.ldscore_new.AnnotationBundle` |
| annotation file parsing | `ldscore.ldscore_new.parse_annotation_file()` |
| combine baseline and query groups | `ldscore.ldscore_new.combine_annotation_groups()` |
| BED -> SNP-level annotation projection | `utils.run_bed_to_annot.run_bed_to_annot()` and `process_baseline_file()` |
| legacy simple single-column projection | `make_annot.py` |

### Design decision

For future cell-type-specific workflows, this class can be reused. That part of your draft is correct.

---

## 5.2 `AnnotationSourceSpec`

This class is missing from the current draft but will make the design much cleaner.

### Role

Represents raw annotation inputs before they become an `AnnotationBundle`.

### Recommendation

Create it as a small dataclass.

### Suggested fields

| Field | Type |
| --- | --- |
| `baseline_annot_paths` | `list[str]` |
| `query_annot_paths` | `list[str]` |
| `bed_paths` | `list[str]` |
| `gene_set_paths` | `list[str]` |
| `allow_missing_query` | `bool` |

### Why it helps

It cleanly separates:
- raw inputs
- aligned domain data
- workflows that create aligned data

### Current-to-target mapping

| Target responsibility | Current implementation |
| --- | --- |
| resolve annotation file inputs | `ldscore.ldscore_new.resolve_annotation_files()` |
| resolve chromosome-split files | `ldscore.ldscore_new.resolve_chr_files()` |
| resolve BED paths | `utils.run_bed_to_annot.expand_bed_paths()` |

---

## 5.3 `RefPanelSpec`

### Role

Describes the reference panel source and loading rules, without eagerly loading the data.

### Recommendation

Create `RefPanelSpec` as a dataclass.

### Suggested fields

| Field | Type | Meaning |
| --- | --- | --- |
| `backend` | `Literal["plink", "parquet_r2"]` | which backend is used |
| `bfile_prefix` | `str \| None` | PLINK prefix |
| `r2_table_paths` | `list[str]` | normalized sorted parquet files |
| `sample_size` | `int \| None` | needed for raw-to-unbiased R2 correction |
| `maf_metadata_paths` | `list[str]` | optional frequency metadata |
| `chromosomes` | `list[str] \| None` | explicit chromosome scope |
| `identifier_mode` | `str` | `rsID` or `chr_pos` |
| `genome_build` | `str \| None` | build-specific interpretation |

### Current-to-target mapping

| Target responsibility | Current implementation |
| --- | --- |
| resolve PLINK prefix | `ldscore.ldscore_new.resolve_bfile_prefix()` |
| resolve parquet files | `ldscore.ldscore_new.resolve_parquet_files()` |
| resolve freq files | `ldscore.ldscore_new.resolve_frequency_files()` |
| parser validation | `ldscore.ldscore_new.validate_args()` |

---

## 5.4 `RefPanel`

### Role

Represents a validated reference panel interface after the source has been resolved.

### Important design decision

Answer to your question: **store the path/spec and lightweight metadata, not the full raw data object, as the long-lived class state.**

Reason:
- PLINK and parquet R2 backends can be very large
- eager loading will waste memory
- many workflows only need chromosome-wise or block-wise access

### Recommendation

Use an abstract base class or protocol-like interface.

### Suggested interface

| Method | Purpose |
| --- | --- |
| `available_chromosomes()` | list chromosomes with usable reference data |
| `load_metadata(chrom)` | return SNP metadata for one chromosome |
| `filter_to_keys(chrom, keys)` | reduce to requested SNP universe |
| `build_reader(chrom)` | create backend-specific reader used by LD-score calculation |
| `summary()` | backend and coverage summary |

### Suggested concrete implementations

#### `PlinkRefPanel`
- stores PLINK prefixes and optional include/keep filters
- lazily constructs `PlinkBEDFile` readers

#### `ParquetR2RefPanel`
- stores parquet paths and optional frequency metadata
- lazily constructs `SortedR2BlockReader`

### What should not be stored permanently

Avoid storing the whole raw matrix / all pairwise R2 rows as a long-lived attribute on `RefPanel`.

Store instead:
- source paths
- backend type
- optional metadata caches
- backend-specific lightweight readers created per chromosome

### Current-to-target mapping

| Target responsibility | Current implementation |
| --- | --- |
| PLINK genotype-backed panel | `ldscore.ldscore.PlinkBEDFile` and parser helpers in `ldscore.parse` |
| parquet-R2 block access | `ldscore.ldscore_new.SortedR2BlockReader` |
| frequency merge for parquet | `ldscore.ldscore_new.parse_frequency_metadata()` and `merge_frequency_metadata()` |
| R2 preprocessing | `ldscore.ldscore_new.convert_r2_table_to_sorted_parquet()` |

---

## 5.5 Why Not `LDScoreRequest`

### Recommendation

Do **not** introduce `LDScoreRequest` as a first-class target object unless future complexity clearly demands it.

### Why I changed this

Your concern is valid: in this codebase, `LDScoreRequest` would mostly be a bundle of arguments passed once into `LDScoreCalculator`.

That usually adds one more layer to read without adding much real structure.

For this project, the simpler design is:
- keep `AnnotationBundle`
- keep `RefPanel`
- keep `LDScoreConfig` and `CommonConfig`
- pass them directly into `LDScoreCalculator.run(...)`
- keep output path information in a small `OutputSpec` dataclass or as explicit function arguments

### Recommended replacement

Use a small `OutputSpec` dataclass instead.

Suggested fields:

| Field | Type |
| --- | --- |
| `out_prefix` | `str` |
| `write_per_chrom` | `bool` |
| `write_annotation_manifest` | `bool` |

### Recommended calculator signature

```python
result = LDScoreCalculator().run(
    annotation_bundle=bundle,
    ref_panel=ref_panel,
    ldscore_config=ldscore_config,
    common_config=common_config,
    output_spec=output_spec,
    regression_snps=regression_snps,
)
```

### When `LDScoreRequest` would become justified

Only introduce it later if one of these becomes important:
- job queuing
- caching and replay
- distributed execution
- a large number of optional run-time switches that are awkward in a direct signature

### Current-to-target mapping

| Target responsibility | Current implementation |
| --- | --- |
| assemble run inputs | `ldscore.ldscore_new.run_ldscore_from_args()` |
| regression SNP loading | `ldscore.ldscore_new.load_regression_keys()` |

---

## 5.6 `LDScoreResult`

### Role

Represents the outputs of LD-score estimation.

### Correction to the current draft

Your current `LDScore` class mixes:
- input configuration
- retained SNP universe
- computation logic
- output files/results

That is too much for one class.

### Recommendation

Replace it with:
- direct validated inputs (`AnnotationBundle`, `RefPanel`, configs, optional `OutputSpec`)
- `LDScoreCalculator` for computation
- `LDScoreResult` for outputs

### Suggested fields

| Field | Type | Meaning |
| --- | --- | --- |
| `metadata` | `pd.DataFrame` | retained reference SNP table |
| `ld_scores` | `pd.DataFrame` | reference LD-score columns |
| `w_ld` | `pd.DataFrame` | regression-weight LD scores |
| `M` | `np.ndarray` | annotation counts |
| `M_5_50` | `np.ndarray \| None` | common-SNP counts |
| `baseline_columns` | `list[str]` | baseline annotation columns |
| `query_columns` | `list[str]` | query annotation columns |
| `reference_snps` | `set[str]` | retained reference SNP identifiers under the selected identifier mode |
| `regression_snps` | `set[str]` | retained regression SNP identifiers; must be a subset of `reference_snps` |
| `per_chrom_results` | `list[ChromLDScoreResult] \| None` | optional chromosome-level outputs |
| `output_paths` | `dict[str, str]` | written files |

### Suggested methods

Light methods only:
- `validate()`
- `to_ldscore_table()`
- `to_weight_table()`
- `summary()`

### Answer to your question about output attributes

**Yes, the final LD-score outputs should exist as attributes, but on a result object, not on the calculator itself.**

That gives you:
- easy reuse for regression workflows
- easy testing
- clear separation between “computation engine” and “result data”

### Current-to-target mapping

| Target responsibility | Current implementation |
| --- | --- |
| per-chromosome result object | `ldscore.ldscore_new.ChromComputationResult` |
| aggregate result assembly | `ldscore.ldscore_new.aggregate_results()` |
| output writing | `ldscore.ldscore_new.emit_outputs()` |
| legacy output writing | output block inside `ldsc.py:ldscore(args)` |

---

## 5.7 `LDScoreCalculator`

### Role

Service class that computes LD scores from validated inputs.

### Suggested public methods

| Method | Purpose |
| --- | --- |
| `run(annotation_bundle, ref_panel, ldscore_config, common_config, output_spec=None, regression_snps=None)` | main public entry point |
| `compute_chromosome(chrom, annotation_bundle, ref_panel, ldscore_config, common_config, regression_snps=None)` | backend-independent chromosome step |
| `write_outputs(result, out_prefix)` | emit LDSC-compatible files |

### Suggested internal methods

| Method | Purpose |
| --- | --- |
| `_compute_from_plink(...)` | PLINK backend |
| `_compute_from_parquet(...)` | parquet backend |
| `_build_reference_universe(...)` | intersect annotation SNPs with ref-panel presence |
| `_build_regression_subset(...)` | enforce regression SNP subset rule |
| `_compute_counts(...)` | produce `.M` and `.M_5_50` |

### Reference-SNP universe rule

Your statement is correct, but I recommend expressing it more strictly:

1. Start from the SNP universe defined by aligned annotation rows.
2. Intersect with SNPs actually present in the reference panel.
3. Apply optional global restriction set.
4. Define regression SNPs as a subset of the retained reference SNP universe.

This rule should be enforced by validation, not left as an informal convention.

### Current-to-target mapping

| Target responsibility | Current implementation |
| --- | --- |
| top-level new LD-score workflow | `ldscore.ldscore_new.run_ldscore_from_args()` |
| parquet chromosome compute | `ldscore.ldscore_new.compute_chrom_from_parquet()` |
| PLINK chromosome compute | `ldscore.ldscore_new.compute_chrom_from_plink()` |
| count computation | `ldscore.ldscore_new.compute_counts()` |
| whole-chromosome window check | `ldscore.ldscore_new.check_whole_chromosome_window()` |

---

## 5.8 `SumstatsTable`

Your draft says “skip for now”, but for a high-level design document this class should be defined at least minimally.

### Role

Represents validated, LDSC-ready summary statistics after loading or munging.

### Recommendation

Use a dataclass.

### Suggested fields

| Field | Type |
| --- | --- |
| `data` | `pd.DataFrame` |
| `identifier_mode` | `str` |
| `has_alleles` | `bool` |
| `source_path` | `str \| None` |
| `trait_name` | `str \| None` |

### Required columns

At minimum:
- `SNP` or canonical SNP identifier representation
- `Z`
- `N`
- optional `A1`, `A2`
- optional additional metadata columns retained for debugging

### Suggested methods

- `validate()`
- `snp_identifiers()`
- `subset_to(keys)`
- `align_to_metadata(metadata)`
- `summary()`

### Current-to-target mapping

| Target responsibility | Current implementation |
| --- | --- |
| read munged sumstats | `ldscore.parse.sumstats()` |
| legacy merge/read helpers | `ldscore.sumstats._read_sumstats()` |
| allele handling constants | `ldscore.sumstats.MATCH_ALLELES`, `FLIP_ALLELES`, `VALID_SNPS` |

---

## 5.9 `SumstatsMunger`

### Role

Service class for converting raw GWAS summary statistics into `SumstatsTable`.

### Suggested public methods

| Method | Purpose |
| --- | --- |
| `run(config)` | full munging workflow |
| `write_output(sumstats, out_prefix)` | write LDSC-ready `.sumstats.gz` |

### Suggested internal methods

| Method | Purpose |
| --- | --- |
| `_infer_columns(...)` | header inference |
| `_filter_chunks(...)` | chunked QC logic |
| `_process_n(...)` | sample-size handling |
| `_restore_signed_z(...)` | sign logic |
| `_merge_alleles(...)` | allele harmonization |

### Current-to-target mapping

| Target responsibility | Current implementation |
| --- | --- |
| full munging workflow | `munge_sumstats.munge_sumstats()` |
| header inference | `read_header()`, `get_cname_map()`, `parse_flag_cnames()` |
| row QC | `filter_pvals()`, `filter_info()`, `filter_frq()`, `filter_alleles()` |
| chunk processing | `parse_dat()` |
| sample-size handling | `process_n()` |
| allele merge | `allele_merge()` |

---

## 5.10 `RegressionDataset`

This is another class worth adding even though it does not exist in the current draft.

### Role

Represents the merged, regression-ready dataset built from `SumstatsTable` and `LDScoreResult`.

### Suggested fields

| Field | Type |
| --- | --- |
| `merged` | `pd.DataFrame` |
| `ref_ld_columns` | `list[str]` |
| `weight_column` | `str` |
| `M` | `np.ndarray` |
| `novar_mask` | `np.ndarray \| None` |
| `trait_names` | `list[str]` |

### Why it is useful

It gives a clean handoff between:
- IO and alignment logic
- regression kernels

### Current-to-target mapping

| Target responsibility | Current implementation |
| --- | --- |
| merged h2 input bundle | `ldscore.sumstats._read_ld_sumstats()` |
| rg merged table | `_read_other_sumstats()` and `_merge_sumstats_sumstats()` |
| variance / condition checks | `_check_variance()`, `_check_ld_condnum()` |

---

## 5.11 `RegressionRunner`

### Role

Service class wrapping the regression suite.

### Suggested public methods

| Method | Purpose |
| --- | --- |
| `estimate_h2(dataset, config)` | single-trait h2 |
| `estimate_partitioned_h2(dataset, annotation_bundle, config)` | partitioned h2 |
| `estimate_rg(dataset1, dataset2, config)` | genetic correlation |
| `estimate_cell_type_specific(...)` | future extension |

### Important design recommendation

Keep the existing numerical estimator classes (`Hsq`, `Gencov`, `RG`) in the numerical layer. Do **not** replace them with high-level workflow objects.

Instead, `RegressionRunner` should:
- build `RegressionDataset`
- decide which estimator to call
- convert outputs into result objects or tables
- write or return summaries

### Current-to-target mapping

| Target responsibility | Current implementation |
| --- | --- |
| h2 workflow | `ldscore.sumstats.estimate_h2()` |
| rg workflow | `ldscore.sumstats.estimate_rg()` |
| cts workflow | `ldscore.sumstats.cell_type_specific()` |
| estimators | `ldscore.regressions.Hsq`, `Gencov`, `RG` |

---

## 6. Feature Design

## 6.1 Feature: Build SNP-Level Annotations

### Public feature goal

Convert user-supplied annotation sources into aligned SNP-level annotation tables suitable for LD-score estimation.

### Inputs

| Input | Required | Notes |
| --- | --- | --- |
| baseline SNP-level `.annot` files | yes | define initial annotation row universe |
| query SNP-level `.annot` files | optional | same row universe required |
| BED files | optional | should be converted via `AnnotationBuilder`, not stored in `AnnotationBundle` |
| gene-set files | optional | only if you choose to support them explicitly |
| optional restriction set | optional | applied after row-universe validation |

### Outputs

- `AnnotationBundle`
- optional written `.annot.gz` files for reuse

### Recommended public API

```python
bundle = AnnotationBuilder(common_config, annotation_build_config).run(source_spec)
```

### Mapping from current implementation

- `ldscore.ldscore_new.parse_annotation_file()`
- `ldscore.ldscore_new.combine_annotation_groups()`
- `utils.run_bed_to_annot.run_bed_to_annot()`
- `make_annot.py`

---

## 6.2 Feature: Prepare Reference Panel

### Public feature goal

Load and validate a reference panel from either PLINK or normalized sorted parquet R2 data.

### Inputs

| Input | Required | Notes |
| --- | --- | --- |
| PLINK prefix or parquet-R2 files | yes | exactly one backend |
| identifier mode | yes | must match annotation workflow |
| genome build | conditional | needed when build-dependent files are used |
| optional MAF/frequency metadata | optional | required for `.M_5_50` in parquet mode |

### Outputs

- `RefPanel`
- optional normalized parquet R2 file if preprocessing is requested

### Recommended public API

```python
ref_panel = RefPanelLoader(common_config, ref_panel_config).load(ref_panel_spec)
```

### Mapping from current implementation

- `ldscore.ldscore_new.convert_r2_table_to_sorted_parquet()`
- `ldscore.ldscore_new.SortedR2BlockReader`
- `ldscore.ldscore.PlinkBEDFile`
- `ldscore.parse.PlinkBIMFile`, `PlinkFAMFile`

---

## 6.3 Feature: LD Score Estimation

### Public feature goal

Compute LDSC-compatible reference LD scores and regression-weight LD scores from a validated annotation bundle and reference panel.

### Inputs

| Input | Required | Notes |
| --- | --- | --- |
| `AnnotationBundle` | yes | SNP-level aligned annotations |
| `RefPanel` | yes | validated backend-specific ref panel |
| regression SNP list | optional | if absent, use retained reference SNP universe |
| LD window definition | yes | exactly one of SNP/kb/cM |

### Outputs

- `LDScoreResult`
- on disk: `.l2.ldscore.gz`, `.w.l2.ldscore.gz`, `.l2.M`, optional `.l2.M_5_50`, annotation manifest

### Important design rules

1. The retained reference SNP universe is the intersection of:
   - annotation rows
   - SNPs present in the reference panel
   - optional global restriction

2. Regression SNPs must be a subset of retained reference SNPs.

3. The result object should carry the final retained universes and the file-ready outputs.

### Recommended public API

```python
result = LDScoreCalculator().run(
    annotation_bundle=bundle,
    ref_panel=ref_panel,
    ldscore_config=ldscore_config,
    common_config=common_config,
    output_spec=output_spec,
    regression_snps=regression_snps,
)
```

### Mapping from current implementation

- `ldscore.ldscore_new.run_ldscore_from_args()`
- `compute_chrom_from_plink()`
- `compute_chrom_from_parquet()`
- legacy `ldsc.py:ldscore(args)`

---

## 6.4 Feature: Summary-Statistics Munging

### Public feature goal

Normalize raw GWAS summary statistics into a simple, validated LDSC-ready table.

### Inputs

- raw summary-statistics file
- optional column-name hints
- optional merge-alleles file
- QC thresholds

### Outputs

- `SumstatsTable`
- `.sumstats.gz`
- log/summary report

### Recommended public API

```python
sumstats = SumstatsMunger(common_config, munge_config).run(raw_source)
```

### Mapping from current implementation

- `munge_sumstats.py`
- `ldscore.parse.sumstats()`

---

## 6.5 Feature: Regression Suite

### Public feature goal

Run heritability, partitioned heritability, and cross-trait genetic correlation from validated LD-score outputs and validated summary statistics.

### Inputs

| Feature | Required inputs |
| --- | --- |
| heritability | one `SumstatsTable`, one `LDScoreResult` |
| partitioned heritability | one `SumstatsTable`, one `LDScoreResult`, annotation grouping metadata |
| genetic correlation | two `SumstatsTable` objects, one compatible `LDScoreResult` |

### Outputs

| Feature | Output |
| --- | --- |
| heritability | `HeritabilityResult` or equivalent summary object |
| partitioned heritability | partitioned result table + summary object |
| genetic correlation | `GeneticCorrelationResult` or equivalent summary object |

### Important practical recommendation

Keep one shared regression backbone.

Your note is correct:
- heritability and partitioned heritability should share the same base workflow
- partitioned h2 should loop over query annotation sets against a fixed baseline model when needed
- cell-type-specific regression should be treated as a later specialization of the same pattern

### Important rule on “all SNPs” baseline

Your note is also correct, but I recommend making it an explicit validation rule:

If the annotation design used for regression does not span the full SNP universe, require an explicit all-ones baseline category and raise a clear error if it is missing.

### Recommended public API

```python
runner = RegressionRunner(common_config, regression_config)
h2 = runner.estimate_h2(sumstats, ldscore_result)
part = runner.estimate_partitioned_h2(sumstats, ldscore_result, annotation_bundle)
rg = runner.estimate_rg(sumstats1, sumstats2, ldscore_result)
```

### Mapping from current implementation

- `ldscore.sumstats.estimate_h2()`
- `ldscore.sumstats.estimate_rg()`
- `ldscore.sumstats.cell_type_specific()`
- `ldscore.regressions.Hsq`
- `ldscore.regressions.Gencov`
- `ldscore.regressions.RG`

---

## 7. Public API Recommendation

If the goal is easiest usage for users, the public interface should be **CLI-first and small**.

### 7.1 Primary public interface: CLI

For this project, the main public interface should be command-line commands, not direct class construction.

Recommended primary commands:
- `ldsc annotate ...`
- `ldsc ldscore ...`
- `ldsc munge-sumstats ...`
- `ldsc regress h2 ...`
- `ldsc regress partitioned-h2 ...`
- `ldsc regress rg ...`

Each command should:
- expose concrete argument names using `snp`, `annotation`, `reference-panel`, `regression-snps`, `window`, `output`
- build the corresponding config dataclasses internally
- call the workflow services
- write user-readable output files and logs

### 7.2 Secondary public interface: Python API

| Public item | Purpose |
| --- | --- |
| `AnnotationBuilder.run(...)` | build `AnnotationBundle` |
| `RefPanelLoader.load(...)` | create `RefPanel` |
| `LDScoreCalculator.run(...)` | compute LD scores |
| `SumstatsMunger.run(...)` | munge summary statistics |
| `RegressionRunner.estimate_*` | run regressions |
| workflow-specific config dataclasses | explicit parameter containers |

### 7.3 Recommended internal-only classes/functions

| Internal item | Why internal |
| --- | --- |
| parser helpers | too low-level for users |
| identifier normalization helpers | internal merge mechanics |
| backend-specific block readers | advanced implementation detail |
| direct `compute_chrom_*` functions | should sit behind calculator APIs |
| file emitters | should be called through result/calculator APIs |

---

## 8. Suggested Implementation Order

This is included because the document will be used for refactoring planning.

1. Freeze identifier rules and config model.
2. Introduce shared domain dataclasses:
   - `CommonConfig`
   - `RefPanelSpec`
   - `AnnotationSourceSpec`
   - `AnnotationBundle`
   - `OutputSpec`
   - `LDScoreResult`
   - `SumstatsTable`
3. Wrap current new LD-score code behind `RefPanelLoader` and `LDScoreCalculator`.
4. Separate BED projection logic into `AnnotationBuilder`.
5. Introduce `RegressionDataset` and `RegressionRunner` without changing `Hsq/Gencov/RG` internals at first.
6. Only after behavior is stable, simplify or retire overlapping legacy entry points.

---

## 9. Final Recommendations And Corrections To Current Draft

### 9.1 Main corrections

- `AnnotationBundle` should hold **aligned SNP-level annotations**, not raw BED files.
- `RefPanel` should hold **paths/specs and lightweight metadata**, not large raw matrices as persistent state.
- `LDScore` should be split into **request + calculator + result**, not one mixed class.
- `config` should be replaced by **explicit immutable workflow config dataclasses**, not one global mutable object.
- `Sumstats` should not be fully skipped in the design; at least define the domain object and workflow boundary now.

### 9.2 Most practical architecture choice

If you want the easiest codebase to maintain and the easiest user experience:
- keep domain objects simple and validated
- keep computation in a small number of workflow services
- keep the statistical kernel mostly unchanged at first
- make user-facing arguments live in workflow config dataclasses and narrow CLI wrappers

That gives the cleanest bridge from the current implementation to the refactored package.
