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
- `rsid`
- `chr_pos`

Internally, create a **canonical SNP identifier** column for every workflow.

Command-line and file-reading recommendation:
- for the user-facing config value, assume the command-line input is simply `rsid`
- when reading files, be flexible about capitalization and common SNP-column aliases
- reuse the current style of column-name inference so headers such as `rsid`, `rsID`, `RSID`, `snpid`, `SNPID`, `snp_id`, `SNP`, and similar variants can all be recognized as the SNP identifier column

### 3.2 Required columns by mode

| Identifier mode | Required columns | Matching rule | Uniqueness rule |
| --- | --- | --- | --- |
| `rsid` | `SNP` | exact string match on SNP ID | retained rows must have unique `SNP` |
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
| SNP list in `rsid` mode | one-column SNP list or table with inferable SNP column | matched to `SNP` |
| SNP list in `chr_pos` mode | table with inferable `CHR` and `BP` columns, or BED-like file | matched to `(CHR, BP)` |
| BED interval restriction | allowed only in workflows that explicitly support interval projection | should not be silently treated as SNP IDs |

### 3.5 Internal recommendation

Always derive one internal SNP identifier column, for example:
- `snp_id = SNP` for `rsid`
- `snp_id = f"{CHR}:{BP}"` for `chr_pos`

All internal merges should use that canonical SNP identifier after validation.

Recommendation on naming:
- avoid `key` in public class fields and public method names
- prefer `snp_id`, `reference_snps`, `regression_snps`, or `snp_identifiers`
- if the identifier mode is `chr_pos`, those SNP identifiers are still the retained SNP labels for the workflow

**Mapping from current implementation**
- the current implementation already has useful column-name inference logic for SNP-like headers; the refactor should keep that flexibility while standardizing the stored config value as `rsid`
- `ldscore.ldscore_new.identifier_keys()`
- `ldscore.ldscore_new.normalize_chromosome()`
- `ldscore.ldscore_new.validate_retained_identifier_uniqueness()`
- `munge_sumstats.default_cnames`
- `munge_sumstats.clean_header()`
- `munge_sumstats.get_cname_map()`
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
| `snp_identifier` | `Literal["rsid", "chr_pos"]` | global SNP identifier mode used across annotations, reference panel loading, sumstats alignment, and output labeling |
| `genome_build` | `Literal["hg19", "hg38"] \| None` | required for build-specific files when needed |
| `global_snp_restriction_path` | `str \| None` | optional global SNP restriction resource applied to annotation SNPs, reference SNPs, and regression SNPs |
| `log_level` | `str` | logging verbosity |
| `fail_on_missing_metadata` | `bool` | strictness on missing `CM`, `MAF`, etc. |

`global_snp_restriction_path` should be interpreted as a **global constraint on SNP universe**, not as a local filter used by only one feature.

That means the following sets are all intersected with it:
- baseline annotation SNPs
- retained reference SNPs used to sum R2 for LD-score computation
- regression SNPs used for regression-weight LD scores

#### `LDScoreConfig` (dataclass)

| Field | Type | Purpose |
| --- | --- | --- |
| `ld_wind_snps` | `int \| None` | SNP-count window |
| `ld_wind_kb` | `float \| None` | base-pair window |
| `ld_wind_cm` | `float \| None` | genetic-distance window |
| `maf_min` | `float \| None` | SNP filter threshold |
| `chunk_size` | `int` | chunking / block processing size |
| `compute_m5_50` | `bool` | whether to compute and emit the extra common-SNP count artifact analogous to `.M_5_50` when MAF is available |
| `whole_chromosome_ok` | `bool` | allow LD-score computation to proceed when the chosen window effectively spans a whole chromosome; this is a safety override, not a numerical option |

`whole_chromosome_ok` is only a guardrail. It does **not** change the LD-score formula. It only controls whether the run is allowed to continue when the chosen window is so large that nearly the whole chromosome would be traversed.

`compute_m5_50` should affect **counts only**, not the main LD-score computation. The main LD scores and regression weights should still be computed on the retained SNP universe after the normal filters. The extra common-SNP artifact should record how many retained reference SNPs per annotation also satisfy `MAF > 0.05` when MAF is available.

#### `RegressionConfig` (dataclass)

| Field | Type | Purpose |
| --- | --- | --- |
| `n_blocks` | `int` | jackknife block count |
| `use_m_5_50` | `bool` | whether regression should use the common-SNP count vector analogous to legacy `.M_5_50`; this should default to `True` to preserve the original LDSC behavior |
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

### 4.4 Avoid repeating global config inside domain classes

Since `snp_identifier` and `genome_build` are global workflow settings, do **not** repeat them as normal stored fields on every class unless there is a concrete reason to snapshot them for provenance.

Recommendation:
- `AnnotationBundle`, `RefPanelSpec`, and `SumstatsTable` should assume they are validated under one `CommonConfig`
- if provenance is needed, store a small immutable `config_snapshot` or `provenance` record at write time rather than duplicating the same field across many classes

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
- apply the global SNP restriction from `CommonConfig`

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

Answer to your question: **store the path/spec and lightweight metadata access strategy, not the full raw data object, as the long-lived class state.**

Reason:
- PLINK and parquet R2 backends can be very large
- eager loading will waste memory
- many workflows only need chromosome-wise or block-wise access
- even metadata tables can be large enough that always embedding them in the long-lived object is wasteful

### Recommendation

Use an abstract base class or protocol-like interface.

### Suggested interface

| Method | Purpose |
| --- | --- |
| `available_chromosomes()` | list chromosomes with usable reference data |
| `load_metadata(chrom)` | return SNP metadata for one chromosome |
| `filter_to_snps(chrom, snps)` | reduce to requested SNP universe |
| `build_reader(chrom)` | create backend-specific reader used by LD-score calculation |
| `summary()` | backend and coverage summary |

### Suggested concrete implementations

#### `PlinkRefPanel`
- stores PLINK prefixes and optional include/keep filters
- lazily constructs `PlinkBEDFile` readers
- may cache chromosome-level SNP metadata DataFrames after first load if repeated access is common

#### `ParquetR2RefPanel`
- stores parquet paths and optional frequency metadata
- lazily constructs `SortedR2BlockReader`
- may cache compact metadata indexes or memory-light DataFrames per chromosome

### What should not be stored permanently

Avoid storing the whole raw matrix / all pairwise R2 rows as a long-lived attribute on `RefPanel`.

Store instead:
- source paths
- backend type
- optional chromosome-level metadata caches
- backend-specific lightweight readers created per chromosome

Practical suggestion on metadata:
- do not eagerly store one giant metadata table on `RefPanel`
- expose `load_metadata(chrom)` and allow lazy caching per chromosome
- if the workflow later needs fully materialized retained metadata, store that on `LDScoreResult`, not on `RefPanel`

Recommended backend-specific metadata strategy:

| Backend | Preferred metadata strategy | Why |
| --- | --- | --- |
| `PlinkRefPanel` | read `.bim` metadata per chromosome and optionally cache the chromosome-level DataFrame in memory | `.bim` access is simple and chromosome filtering is cheap |
| `ParquetR2RefPanel` | keep path-based sidecar metadata resources per chromosome; load only requested columns and optionally cache a compact index/DataFrame | parquet metadata and frequency sidecars can be large, so path-first access scales better |

Practical rule:
- `RefPanel` owns metadata access
- `LDScoreResult` owns fully materialized retained metadata actually used in a run
- `OutputManager` decides whether that retained metadata is written to disk

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

`OutputSpec` is the part of the old “request object” idea that is actually useful: it says what artifacts to write, where to write them, and how to organize them. It should not carry runtime data or numerical algorithm settings.

Suggested fields:

| Field | Type |
| --- | --- |
| `out_prefix` | `str` |
| `output_dir` | `str \| None` |
| `artifact_layout` | `Literal["flat", "by_chrom", "run_dir"]` |
| `write_ldscore` | `bool` |
| `write_w_ld` | `bool` |
| `write_counts` | `bool` |
| `write_per_chrom` | `bool` |
| `aggregate_across_chromosomes` | `bool` |
| `write_annotation_manifest` | `bool` |
| `compression` | `Literal["gzip", "none"]` |
| `overwrite` | `bool` |
| `log_path` | `str \| None` |
| `write_summary_json` | `bool` |
| `write_summary_tsv` | `bool` |
| `write_run_metadata` | `bool` |
| `enabled_artifacts` | `list[str] \| None` |

Design note:
- always include `CM` and `MAF` in exported tables when available
- if missing, write `NA`
- do not make `include_cm` and `include_maf` separate user options unless there is a real downstream need; default clean usage matters more than option count
- `enabled_artifacts=None` should mean “use the clean built-in default artifact set”
- keep artifact-specific plotting or summary options out of `OutputSpec`; if needed later, add a separate advanced artifact config rather than bloating the core output config

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
| `reference_metadata` | `pd.DataFrame` | fully materialized retained reference SNP metadata aggregated across chromosomes in final genomic order and aligned row-for-row with `ld_scores` |
| `ld_scores` | `pd.DataFrame` | reference LD-score columns aggregated across chromosomes and aligned row-for-row with `reference_metadata` |
| `regression_metadata` | `pd.DataFrame` | retained regression SNP metadata aggregated across chromosomes in final genomic order and aligned row-for-row with `w_ld` |
| `w_ld` | `pd.DataFrame` | regression-weight LD-score columns aggregated across chromosomes and aligned row-for-row with `regression_metadata` |
| `snp_count_totals` | `dict[str, np.ndarray]` | named SNP-count arrays obtained by summing the chromosome-specific count arrays across chromosomes |
| `baseline_columns` | `list[str]` | baseline annotation columns |
| `query_columns` | `list[str]` | query annotation columns |
| `reference_snps` | `set[str]` | retained reference SNP identifiers under the selected identifier mode |
| `regression_snps` | `set[str]` | retained regression SNP identifiers; must be a subset of `reference_snps` |
| `chromosome_results` | `list[ChromLDScoreResult]` | chromosome-level LD-score results produced during the run and then aggregated into the combined result used by regression |
| `output_paths` | `dict[str, str]` | written files |

### Internal helper result: `ChromLDScoreResult`

Define one small chromosome-level result object and treat it as a first-class intermediate in the LD-score pipeline.

Suggested fields:
- `chrom`
- `reference_snps`
- `regression_snps`
- `reference_metadata`
- `ld_scores`
- `regression_metadata`
- `w_ld`
- `snp_count_totals`

Practical meaning:
- it is the chromosome-local version of `LDScoreResult`
- it is the natural unit of LD-score computation because the current pipeline operates chromosome by chromosome for efficiency
- it exists so `OutputManager` can either write one combined artifact set or one artifact set per chromosome
- the final combined `LDScoreResult` is built by concatenating and summing these chromosome-level results in genomic order
- it should not be the main user-facing result type

### Suggested methods

Light methods only:
- `validate()`
- `to_ldscore_table()`
- `to_weight_table()`
- `summary()`

Alignment invariants:
- `reference_metadata` and `ld_scores` must have the same number of rows and the same SNP order
- `regression_metadata` and `w_ld` must have the same number of rows and the same SNP order
- `regression_snps` must correspond exactly to the SNPs represented by `regression_metadata`
- if `w_ld` is derived from a regression mask within the reference universe, preserve the explicit regression metadata table anyway so downstream users do not have to reconstruct it
- `chromosome_results` should together cover exactly the SNPs represented in the combined result, with no cross-chromosome duplication

Important structural rule:
- LD-score calculation is a chromosome-wise workflow
- combined `LDScoreResult` exists mainly as the post-aggregation object used by downstream regression and by aggregate output writing
- chromosome-specific `ld_scores`, `w_ld`, and SNP-count totals remain first-class intermediate results, not temporary implementation details

### Answer to your question about output attributes

**Yes, the final LD-score outputs should exist as attributes, but on a result object, not on the calculator itself.**

That gives you:
- easy reuse for regression workflows
- easy testing
- clear separation between “computation engine” and “result data”

Naming note:
- storing `M` and `M_5_50` as bare field names is too legacy-oriented
- prefer a named dictionary such as `snp_count_totals`
- example keys:
  - `all_reference_snp_counts`
  - `common_reference_snp_counts_maf_gt_0_05`

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
| `compute_chromosome(chrom, annotation_bundle, ref_panel, ldscore_config, common_config, regression_snps=None)` | backend-independent chromosome step producing one `ChromLDScoreResult` |
| `write_outputs(result, out_prefix)` | emit LDSC-compatible files |

### Suggested internal methods

| Method | Purpose |
| --- | --- |
| `_compute_from_plink(...)` | PLINK backend |
| `_compute_from_parquet(...)` | parquet backend |
| `_partition_annotation_rows_by_chrom(...)` | split aligned annotation rows into chromosome-specific bundles |
| `_partition_regression_snps_by_chrom(...)` | partition requested regression SNPs to the chromosome-specific retained SNP universes |
| `_build_reference_universe(...)` | intersect annotation SNPs with ref-panel presence |
| `_build_regression_mask(...)` | build a boolean mask on the retained reference SNP universe for regression SNPs |
| `_augment_annotation_matrix_with_regression_mask(...)` | append an internal mask column so `w_ld` is computed in the same LD-score traversal as the main annotation columns |
| `_build_regression_subset(...)` | enforce regression SNP subset rule |
| `_compute_counts(...)` | produce named SNP count arrays for all retained reference SNPs and optional common-SNP subsets |
| `_aggregate_chromosome_results(...)` | combine chromosome-specific outputs into the final `LDScoreResult` used by regression |

### Reference-SNP universe rule

Your statement is correct, but I recommend expressing it more strictly:

1. Start from the SNP universe defined by aligned annotation rows.
2. Intersect with SNPs actually present in the reference panel.
3. Apply optional global restriction set.
4. Define regression SNPs as a subset of the retained reference SNP universe.

This rule should be enforced by validation, not left as an informal convention.

Recommended algorithm detail:
- partition the aligned annotation rows by chromosome before LD-score computation
- if a global regression SNP list is supplied, partition it automatically to per-chromosome subsets after matching against the retained per-chromosome SNP universe
- add one internal annotation-like mask column corresponding to the retained regression SNP set
- compute reference LD-score columns and regression-weight LD scores in the same per-chromosome traversal
- do not require a second separate LD-score pass just to obtain `w_ld`
- emit one `ChromLDScoreResult` per chromosome
- after chromosome-wise computation, aggregate the ordinary annotation-derived LD-score columns, regression-weight LD scores, and SNP-count totals into the final combined `LDScoreResult`
- after computation, split the ordinary annotation-derived LD-score columns from the regression-mask-derived `w_ld` column in `LDScoreResult`

Current implementation note:
- the parquet backend already follows this design by appending a regression-mask column before one joint LD-score traversal
- the PLINK backend still computes `w_ld` in a second kernel call
- the refactor should unify both backends behind the same higher-level behavior
- the current pipeline already computes chromosome by chromosome and then aggregates results; the refactor should make that structure explicit in the public design rather than hiding it

Practical rule for regression SNP partitioning:
- if `snp_identifier == "chr_pos"`, chromosome partitioning is direct from the identifier itself
- if `snp_identifier == "rsid"`, partitioning should be done after matching the requested regression SNPs against chromosome-specific retained metadata/reference-panel rows
- users should pass one global regression-SNP input; the calculator should partition it automatically

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
| `has_alleles` | `bool` |
| `source_path` | `str \| None` |
| `trait_name` | `str \| None` |
| `provenance` | `dict[str, object]` |

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
- `subset_to(snps)`
- `align_to_metadata(metadata)`
- `summary()`

Validation rules should follow the old codebase behavior as closely as possible:
- preserve the existing LDSC-required columns and semantics
- preserve allele requirements for workflows that need them
- preserve the existing assumptions around `Z`, `N`, and optional `A1/A2`
- avoid introducing new sumstats features at this stage

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

This should follow the legacy LDSC workflow closely and deliberately avoid new feature additions in the first restructuring pass.

### Suggested public methods

| Method | Purpose |
| --- | --- |
| `run(raw_source, munge_config, common_config)` | full munging workflow using the legacy pipeline behavior |
| `write_output(sumstats, out_prefix)` | write LDSC-ready `.sumstats.gz` |
| `build_run_summary(sumstats)` | create a compact summary of what was removed and retained |

### Suggested internal methods

| Method | Purpose |
| --- | --- |
| `_infer_columns(...)` | header inference |
| `_filter_chunks(...)` | chunked QC logic |
| `_process_n(...)` | sample-size handling |
| `_restore_signed_z(...)` | sign logic |
| `_merge_alleles(...)` | allele harmonization |
| `_validate_against_legacy_schema(...)` | ensure outputs still match the current LDSC expectations |

Suggested behavior:
- keep the same logical pipeline as `munge_sumstats.py`
- preserve legacy-compatible output schema
- move complexity into internal helpers, not new user-visible options
- keep the CLI surface simple and mostly familiar

### Current-to-target mapping

| Target responsibility | Current implementation |
| --- | --- |
| full munging workflow | `munge_sumstats.munge_sumstats()` |
| header inference | `read_header()`, `get_cname_map()`, `parse_flag_cnames()` |
| row QC | `filter_pvals()`, `filter_info()`, `filter_frq()`, `filter_alleles()` |
| chunk processing | `parse_dat()` |
| sample-size handling | `process_n()` |
| allele merge | `allele_merge()` |

### Recommended supporting classes

To keep the munging workflow explicit without changing behavior, add these lightweight supporting objects.

#### `RawSumstatsSpec`

| Field | Type | Purpose |
| --- | --- | --- |
| `path` | `str` | raw summary-statistics file |
| `compression` | `Literal["auto", "gzip", "bz2", "none"]` | input compression handling |
| `trait_name` | `str \| None` | optional trait label for downstream summaries |
| `column_hints` | `dict[str, str]` | optional explicit column-name overrides corresponding to legacy CLI flags |

#### `MungeConfig`

| Field | Type | Purpose |
| --- | --- | --- |
| `out_prefix` | `str` | output prefix |
| `N` | `float \| None` | fixed sample size override |
| `N_cas` | `float \| None` | fixed case count override |
| `N_con` | `float \| None` | fixed control count override |
| `info_min` | `float` | INFO filter threshold |
| `maf_min` | `float` | MAF filter threshold |
| `n_min` | `float \| None` | minimum sample size threshold |
| `nstudy_min` | `float \| None` | minimum number-of-studies threshold |
| `chunk_size` | `int` | chunked reading size |
| `merge_alleles_path` | `str \| None` | optional allele-reference file |
| `signed_sumstats_spec` | `str \| None` | legacy signed-stat specification such as `OR,1` or `Z,0` |
| `ignore_columns` | `list[str]` | legacy ignore-list equivalent |
| `no_alleles` | `bool` | skip allele requirement |
| `a1_inc` | `bool` | treat `A1` as increasing allele |
| `keep_maf` | `bool` | preserve MAF column in output when available |
| `daner` | `bool` | legacy Daner format handling |
| `daner_n` | `bool` | legacy Daner-with-Nca/Nco handling |

#### `MungeRunSummary`

| Field | Type | Purpose |
| --- | --- | --- |
| `n_input_rows` | `int` | total rows read |
| `n_retained_rows` | `int` | rows remaining after munging |
| `drop_counts` | `dict[str, int]` | counts removed by each QC step |
| `inferred_columns` | `dict[str, str]` | final raw-to-canonical column mapping |
| `used_n_rule` | `str` | how sample size was determined |
| `output_paths` | `dict[str, str]` | written artifacts |

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
| `reference_snp_count_totals` | `dict[str, np.ndarray]` |
| `count_key_used_for_regression` | `str` |
| `retained_ld_columns` | `list[str]` |
| `dropped_zero_variance_ld_columns` | `list[str]` |
| `trait_names` | `list[str]` |
| `chromosomes_aggregated` | `list[str]` |

### Why it is useful

It gives a clean handoff between:
- IO and alignment logic
- regression kernels

Important scope rule:
- `RegressionDataset` is built only after all chromosome-level LD-score results are combined
- regression itself is not run chromosome by chromosome
- chromosome-specific computation belongs to `LDScoreCalculator`, while regression consumes the aggregated cross-chromosome result

Count-config rule:
- preserve the original LDSC default by using the common-SNP count vector corresponding to `.M_5_50` for regression unless the user explicitly requests otherwise
- `count_key_used_for_regression` should therefore default to something like `common_reference_snp_counts_maf_gt_0_05`
- using that count vector does **not** mean the LD scores themselves were computed only on common SNPs; it means the regression scaling/count input follows the legacy default

Naming correction:
- `novar_mask` is too implementation-oriented and too vague for this design document
- prefer explicit fields such as `retained_ld_columns` and `dropped_zero_variance_ld_columns`
- in the current code, `novar_cols` is a boolean mask/Series marking LD-score columns with zero variance that are removed before regression; the target design should expose the dropped column names directly instead

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
| `estimate_partitioned_h2(dataset, annotation_bundle, config)` | partitioned h2 for one query annotation or one already-selected query set |
| `estimate_partitioned_h2_batch(dataset, annotation_bundle, config)` | loop over multiple query annotations, run one partitioned LDSC model per query annotation, and combine/save the results |
| `estimate_rg(dataset1, dataset2, config)` | genetic correlation |
| `estimate_cell_type_specific(...)` | future extension |

### Important design recommendation

Keep the existing numerical estimator classes (`Hsq`, `Gencov`, `RG`) in the numerical layer. Do **not** replace them with high-level workflow objects.

Instead, `RegressionRunner` should:
- build `RegressionDataset`
- decide which estimator to call
- convert outputs into result objects or tables
- write or return summaries

Required new refactor feature:
- when users provide multiple query annotations, `RegressionRunner` should support a batch partitioned-h2 mode
- this mode should loop over query annotations one at a time against the fixed baseline model
- each loop iteration should run one partitioned LDSC regression
- the runner should then combine the per-query results into one clean summary table and save any requested detailed artifacts
- this is **not implemented in the old codebase** and should be treated as explicit new work during the restructuring

### Current-to-target mapping

| Target responsibility | Current implementation |
| --- | --- |
| h2 workflow | `ldscore.sumstats.estimate_h2()` |
| rg workflow | `ldscore.sumstats.estimate_rg()` |
| cts workflow | `ldscore.sumstats.cell_type_specific()` |
| estimators | `ldscore.regressions.Hsq`, `Gencov`, `RG` |

---

## 5.12 Output And Post-Processing Layer

Since you plan to build a richer output and post-processing suite, this layer should be part of the design now.

### Recommended classes

| Class | Type | Role |
| --- | --- | --- |
| `OutputSpec` | dataclass | user-facing output configuration for what artifacts to write and how to organize them |
| `ArtifactConfig` | optional dataclass | advanced per-artifact config for future summary/plot/report options without bloating `OutputSpec` |
| `RunSummary` | dataclass | compact derived summary for reporting and machine-readable metadata |
| `OutputManager` | service class | coordinates output writing and path construction |
| `ArtifactProducer` | protocol/service interface | pluggable artifact builder for one output type or one group of related outputs |
| `ResultWriter` | service/helper | writes primary artifacts such as LD-score tables, count files, manifests |
| `ResultFormatter` | service/helper | builds clean user-facing summary tables and manifests |
| `PostProcessor` | optional service layer | later-stage downstream organization, dashboard inputs, plots, combined reports |

Design rule:
- keep `OutputManager` stable and extensible
- future plots, summaries, and reports should be added through registered `ArtifactProducer` implementations rather than by turning `OutputManager` into a large hardcoded dispatcher

### `OutputSpec` recommended fields

| Field | Type | Purpose |
| --- | --- | --- |
| `out_prefix` | `str` | main output prefix |
| `output_dir` | `str \| None` | optional base directory |
| `artifact_layout` | `Literal["flat", "by_chrom", "run_dir"]` | output directory organization |
| `write_ldscore` | `bool` | write reference LD-score table |
| `write_w_ld` | `bool` | write regression-weight LD-score table |
| `write_counts` | `bool` | write SNP count artifacts |
| `write_annotation_manifest` | `bool` | write baseline/query annotation manifest |
| `write_per_chrom` | `bool` | emit chromosome-split artifacts |
| `aggregate_across_chromosomes` | `bool` | emit aggregated artifacts |
| `compression` | `Literal["gzip", "none"]` | output compression |
| `overwrite` | `bool` | whether existing files may be replaced |
| `log_path` | `str \| None` | optional explicit log path |
| `write_summary_json` | `bool` | machine-readable summary |
| `write_summary_tsv` | `bool` | compact human-readable summary |
| `write_run_metadata` | `bool` | provenance/config dump |
| `enabled_artifacts` | `list[str] \| None` | optional explicit list of registered artifacts to emit; `None` means use the clean built-in defaults |

Design defaults:
- include metadata columns such as `CM` and `MAF` whenever available
- fill unavailable metadata with `NA`
- keep the default artifact set clean and small for ordinary users
- do not put plot-specific or report-specific tuning options directly into `OutputSpec`
- if future artifacts need extra options, store them in a separate advanced `ArtifactConfig`

### `RunSummary` recommended contents

| Field | Purpose |
| --- | --- |
| `n_reference_snps` | retained reference SNP count |
| `n_regression_snps` | retained regression SNP count |
| `chromosomes_processed` | chromosomes included in the run |
| `count_artifacts_available` | whether common-SNP counts and total counts were emitted |
| `output_paths` | actual written artifacts |
| `config_snapshot` | compact provenance record |

### `ArtifactProducer` recommended interface

Each producer should define:
- `name`: stable artifact name used in config/registration
- `supports(result)`: whether the producer can handle the given result type
- `build(result, run_summary, output_spec, artifact_config=None)`: create the artifact content or artifact object
- `write(...)` or return a standardized artifact object for `ResultWriter` to serialize

Built-in producers should cover:
- LD-score tables
- regression-weight tables
- count files
- annotation manifest
- run summary TSV/JSON
- run metadata dump

Future producers can cover:
- partitioned-h2 query summary tables
- batch partitioned-h2 combined summary tables
- plots
- combined reports
- dashboard-ready derived tables

### Recommended output working streamline

1. A workflow service returns a raw result object such as `LDScoreResult` or a regression result.
2. `ResultFormatter` derives built-in user-facing tables from that raw result.
3. `RunSummary` is built from the raw result plus the formatter outputs.
4. `OutputManager` resolves which registered artifact producers are enabled by `OutputSpec`.
5. Each enabled `ArtifactProducer` builds one artifact or one related group of artifacts from the result, summary, and config.
6. `ResultWriter` writes primary machine-readable artifacts and any producer outputs that need file serialization.
7. `PostProcessor` can later register higher-level summary, plot, or report producers without changing the core computation classes.

For LD-score outputs specifically:
- `OutputManager` should be able to write either per-chromosome artifacts from `chromosome_results` or one combined artifact set from the aggregated `LDScoreResult`
- the default user-facing regression workflow should consume the combined cross-chromosome result
- future partitioned-h2 summaries and plots should be added as registered producers, not as special-case branches in the core writer

### Default output behavior

Keep the default output set small, readable, and directly useful:
- one main machine-readable result per feature
- one compact human-readable summary when applicable
- optional detailed artifacts for debugging, provenance, or downstream pipelines

Do not make users turn on many switches just to get the normal outputs they care about.

### Output design rule for downstream analyses

Default outputs should favor what most users care about first.

Example:
- if multiple query pathways are tested in partitioned heritability, the default user-facing summary should be one clean table with one row per query pathway and only the key query-pathway results
- baseline-category details can still be written as secondary artifacts, but should not dominate the default display

This same principle should guide the post-processing layer:
- default outputs should be concise
- detailed technical artifacts should still be available
- result formatting should be separated from raw result storage
- extensibility should come from registration and composition, not from hardcoded one-off branches in the output core

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
| global SNP restriction from `CommonConfig` | optional | applied after row-universe validation and treated as a global SNP-universe constraint |

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
| global SNP identifier mode from `CommonConfig` | yes | must match annotation workflow |
| genome build | conditional | needed when build-dependent files are used |
| optional MAF/frequency metadata | optional | required in parquet mode if you want the common-SNP count artifact analogous to `.M_5_50` |

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
- on disk through `OutputManager`: legacy-compatible LDSC artifacts by default, plus optional clean summaries and manifests

### Recommended working streamline

1. Load aligned annotation rows from `AnnotationBundle`.
2. Partition annotation rows by chromosome.
3. For each chromosome, intersect with reference-panel availability and apply the global SNP restriction from `CommonConfig`.
4. If a global regression SNP list is provided, partition it automatically to chromosome-specific retained SNP subsets.
5. For each chromosome, build one internal regression-mask column and compute chromosome-specific `ld_scores`, `w_ld`, and SNP-count totals.
6. Store those chromosome-level outputs as `ChromLDScoreResult` objects.
7. Aggregate all chromosome-level outputs into the final combined `LDScoreResult` used by regression and by aggregate output writing.

### Important design rules

1. The retained reference SNP universe is the intersection of:
   - annotation rows
   - SNPs present in the reference panel
   - optional global restriction

2. Regression SNPs must be a subset of retained reference SNPs.

3. If regression SNPs are provided as one global input, they should be partitioned automatically into chromosome-specific subsets during LD-score computation.
4. `ChromLDScoreResult` should carry chromosome-specific retained universes and file-ready outputs.
5. The final combined `LDScoreResult` should be built only after aggregating all chromosome-level results.
6. `w_ld` should be computed in the same chromosome traversal as the main LD-score columns by appending an internal regression-mask column.
7. `compute_m5_50` should only create an extra common-SNP count artifact when MAF is available; it should not redefine the main retained SNP universe or the main LD-score calculation.
8. Regression SNPs are not automatically redefined to `MAF > 0.05`; if a future workflow wants a different regression SNP config, it should be an explicit separate option.

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

This feature should preserve the **current legacy LDSC workflow and behavior**. No feature additions are planned in this stage.

### Inputs

| Input | Required | Notes |
| --- | --- | --- |
| raw summary-statistics file | yes | same role as current `--sumstats` input |
| optional column-name hints | optional | same role as current header-mapping flags |
| optional merge-alleles file | optional | preserve current allele harmonization behavior |
| QC thresholds | optional | preserve current INFO, MAF, p-value, N-related filtering behavior |
| shared `CommonConfig` | yes | provides the global SNP identifier config and logging defaults |
| `MungeConfig` | yes | workflow-specific munging options |

### Outputs

- `SumstatsTable`
- `.sumstats.gz`
- log/summary report
- optional `MungeRunSummary`

### Recommended public API

```python
sumstats = SumstatsMunger().run(raw_source, munge_config, common_config)
```

### Recommended working streamline

1. Read the raw header and infer canonical column names using legacy LDSC rules plus any explicit user hints.
2. Read the input in chunks.
3. Apply chunk-level QC and filtering for missing values, INFO, MAF, p-value bounds, SNP/allele validity, and optional allele-merge restrictions.
4. Concatenate retained rows and determine sample size using the same priority rules as the current codebase.
5. Restore or derive the final signed statistic / `Z` representation following the legacy munging behavior.
6. Build `SumstatsTable`, write `.sumstats.gz`, and emit `MungeRunSummary`.

### Mapping from current implementation

- `munge_sumstats.py`
- `ldscore.parse.sumstats()`

### Recommended design rule

Until later development stages, keep the munging feature intentionally conservative:
- preserve the old codebase pipeline
- preserve current output schema
- reorganize structure, not behavior

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
| partitioned heritability | query-focused summary table by default, plus optional detailed baseline-category artifacts; when multiple query annotations are provided, produce one combined results table with one row per query annotation |
| genetic correlation | `GeneticCorrelationResult` or equivalent summary object |

### Important practical recommendation

Keep one shared regression backbone.

Legacy-default count config:
- preserve the original LDSC default and use `M_5_50`-style common-SNP counts in regression by default
- only switch to all-SNP counts analogous to `.M` if the user explicitly requests that override
- this default should be represented explicitly in `RegressionConfig` and in the selected `count_key_used_for_regression`

Your note is correct:
- heritability and partitioned heritability should share the same base workflow
- partitioned h2 should loop over query annotations or query annotation sets against a fixed baseline model when needed
- cell-type-specific regression should be treated as a later specialization of the same pattern

Required new feature for the refactor:
- if users provide multiple query annotations, the code should automatically run one partitioned LDSC regression per query annotation
- each run should use the same baseline model and the same aggregated cross-chromosome LD-score inputs
- the outputs from those runs should then be combined into one saved summary table
- this feature is a planned addition for the refactored codebase and is **not present in the old implementation**

### Important rule on “all SNPs” baseline

Your note is also correct, but I recommend making it an explicit validation rule:

If the annotation design used for regression does not span the full SNP universe, require an explicit all-ones baseline category and raise a clear error if it is missing.

### Recommended public API

```python
runner = RegressionRunner(common_config, regression_config)
h2 = runner.estimate_h2(sumstats, ldscore_result)
part = runner.estimate_partitioned_h2(sumstats, ldscore_result, annotation_bundle)
part_batch = runner.estimate_partitioned_h2_batch(sumstats, ldscore_result, annotation_bundle)
rg = runner.estimate_rg(sumstats1, sumstats2, ldscore_result)
```

### Recommended working streamline

1. Build `RegressionDataset` by merging munged summary statistics with reference LD scores and regression-weight LD scores.
2. Remove zero-variance LD-score columns and record which columns were dropped.
3. Confirm that the LD-score inputs already represent the combined cross-chromosome result assembled from chromosome-specific LD-score runs.
4. Select the appropriate SNP-count vector from `reference_snp_count_totals` for the regression being run; by default this should be the common-SNP count vector corresponding to `.M_5_50`.
5. For ordinary partitioned h2 with one query annotation, call the existing numerical estimator class through `RegressionRunner`.
6. For batch partitioned h2 with multiple query annotations, loop over query annotations one at a time, run one partitioned LDSC model per query annotation, and collect the per-query outputs.
7. Combine the per-query outputs into one clean summary table with one row per query annotation, then optionally emit richer detailed artifacts through the output layer.
8. For h2 and rg, call the existing numerical estimator classes (`Hsq`, `Gencov`, `RG`) through `RegressionRunner`.

Important note on interpretation:
- in this refactor we intentionally preserve the legacy/default LDSC behavior of using `M_5_50` in regression
- this should be documented as a regression count input choice, not as a statement that LD scores were computed only from SNPs with `MAF > 0.05`

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
- `LDScore` should be split into **validated inputs + calculator + result**, not one mixed class.
- `config` should be replaced by **explicit immutable workflow config dataclasses**, not one global mutable object.
- `Sumstats` should not be fully skipped in the design; at least define the domain object and workflow boundary now.

### 9.2 Most practical architecture choice

If you want the easiest codebase to maintain and the easiest user experience:
- keep domain objects simple and validated
- keep computation in a small number of workflow services
- keep the statistical kernel mostly unchanged at first
- make user-facing arguments live in workflow config dataclasses and narrow CLI wrappers

That gives the cleanest bridge from the current implementation to the refactored package.
