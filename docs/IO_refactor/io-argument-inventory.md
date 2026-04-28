# IO Argument Inventory for Refactor Planning

Date: 2026-04-28

This document inventories user-facing input and output path arguments in the
refactored LDSC package. The goal is to support a later API cleanup where each
logical input or output object is represented by one consistent argument shape:
use `*_path` for one file-like input, `*_paths` for logical inputs that may
resolve to multiple files through globs or `@` chromosome tokens, and `*_dir`
for output directories. Avoid public output prefixes unless a legacy surface
still requires them.

## Scope

Included surfaces:

- Console entry point: `ldsc=ldsc.cli:main`.
- Public names re-exported by `ldsc.__all__`.
- Lazy public exports from `ldsc.__getattr__`.
- Public annotation objects re-exported by `ldsc.annotation_builder.__all__`.
- Public path-resolution helpers, because they define the package's current IO
  token semantics.
- Internal `_kernel.ref_panel` and `_kernel.annotation` classes only where they
  are re-exported through the public API.

Excluded surfaces:

- Private helpers and underscore-prefixed modules unless re-exported.
- Tests, tutorials, notebooks, and docs.
- Non-filesystem arguments such as thresholds, column names, trait names,
  model parameters, and genome-build switches, except where they affect output
  layout or path interpretation.

## Current Path Token Semantics

The current workflow layer already has a partial vocabulary in
`src/ldsc/path_resolution.py`.

| Function | Argument | Direction | Current semantic |
|---|---:|---:|---|
| `normalize_path_token(path)` | `path` | input/output token | Expands `~` and environment variables. Does not check the filesystem. |
| `normalize_path_tokens(values)` | `values` | input tokens | Accepts one token or a sequence and returns normalized strings. |
| `split_cli_path_tokens(value)` | `value` | input tokens | Splits CLI comma-delimited strings and `nargs` lists into normalized tokens. |
| `resolve_scalar_path(token, suffixes=(), label="input")` | `token` | input path | Exact path or glob that must resolve to exactly one file. |
| `resolve_file_group(tokens, ..., allow_chromosome_suite=False)` | `tokens` | input paths | Exact paths, globs, and optionally `@` chromosome-suite placeholders. Returns a deterministic deduplicated file list. |
| `resolve_chromosome_group(tokens, chrom, ...)` | `tokens` | input paths | Resolves group tokens for one chromosome, with filename-based chromosome filtering and `@` substitution. |
| `resolve_plink_prefix(token, chrom=None)` | `token` | input prefix | Exact PLINK prefix, `.bed/.bim/.fam` path, exact-one glob, or `@` suite token when a chromosome is supplied. |
| `resolve_plink_prefix_group(tokens, chrom=None, allow_chromosome_suite=False)` | `tokens` | input prefixes | PLINK prefix analogue of `resolve_file_group`; validates complete `.bed/.bim/.fam` trios. |
| `filter_paths_for_chromosome(paths, chrom)` | `paths` | input paths | Filters path names that encode `chrom`. |
| `ensure_output_directory(path)` | `path` | output dir | Creates the directory if missing; rejects existing non-directories. |
| `ensure_output_parent_directory(path)` | `path` | output path/prefix | Creates the parent directory for a path or prefix-like destination. |

Current terms are mixed:

- `path`: used for scalar files, path tokens, and sometimes output-like files.
- `paths`: used for file groups that may accept globs and `@`.
- `dir` / `output_dir`: used when a directory is the object.
- `prefix`: used for PLINK inputs and many legacy outputs.
- `out`: CLI shorthand that sometimes means output prefix and sometimes output
  directory.

## CLI Inventory

### `ldsc annotate`

Defined in `src/ldsc/cli.py` and dispatched to
`ldsc.annotation_builder.main_bed_to_annot`.

| Flag | Direction | Required | Current type | Current object | Resolution behavior | Notes |
|---|---:|---:|---|---|---|---|
| `--bed-files` | input | effectively yes | `nargs="+"`; comma-splittable tokens | BED file group | `split_cli_path_tokens` then `resolve_file_group(..., label="BED file")` | Accepts exact files and globs. In single-universe loading, BED basenames become query column names. |
| `--baseline-annot` | input | effectively yes | `nargs="+"`; comma-splittable tokens | baseline annotation file group | `split_cli_path_tokens` then `resolve_file_group(..., suffixes=ANNOTATION_SUFFIXES, allow_chromosome_suite=True)` | Accepts exact files, globs, and `@` suites. |
| `--output-dir` | output | effectively yes for standalone parser | scalar directory | generated query `.annot.gz` directory | `ensure_output_directory` | Unified CLI parser marks it optional, but `_run_annotate` requires `--bed-files` and `--baseline-annot`; standalone parser requires output dir. |

Output naming:

- `AnnotationBuilder.project_bed_annotations(..., output_dir=...)` writes
  `query.<chrom>.annot.gz` files below `output_dir`.
- Lower-level BED projection helpers may write `<baseline-derived>.annot.gz`
  into nested BED-name directories in compatibility paths.

Non-path flags that affect IO interpretation:

- `--snp-identifier`, `--genome-build`, `--no-batch`, `--log-level`.
- `--no-batch` is retained for compatibility and no longer changes the public
  output shape.

### `ldsc ldscore`

Defined in `src/ldsc/ldscore_calculator.py` and cloned into `ldsc.cli`.

| Flag | Direction | Required | Current type | Current object | Resolution behavior | Notes |
|---|---:|---:|---|---|---|---|
| `--output-dir` | output | yes | scalar directory | canonical LD-score result directory | `LDScoreOutputConfig(output_dir=...)`, then `ensure_output_directory` | Writes `manifest.json`, `baseline.parquet`, and optional `query.parquet`. |
| `--query-annot` | input | no | comma-delimited token string | query annotation file group | `split_cli_path_tokens` then annotation builder group resolution | Mutually exclusive with `--query-annot-bed`. |
| `--query-annot-bed` | input | no | comma-delimited token string | query BED file group | `split_cli_path_tokens` then `resolve_file_group(label="BED file")` | Projected in memory as query annotations. |
| `--baseline-annot` | input | no in parser, required by workflow | comma-delimited token string | baseline annotation file group | `split_cli_path_tokens` then annotation builder group resolution | Accepts exact paths, globs, and `@` suites. |
| `--bfile` | input | conditional | scalar PLINK prefix token | PLINK reference panel | `resolve_plink_prefix(..., chrom=...)` | Used when `--r2-table` is absent. |
| `--r2-table` | input | conditional | comma-delimited token string | parquet R2 file group | `resolve_chromosome_group(..., suffixes=PARQUET_SUFFIXES)` | If present, workflow selects parquet backend. |
| `--frqfile` | input | no | comma-delimited token string | frequency/metadata sidecar group | `resolve_chromosome_group(..., suffixes=FREQUENCY_SUFFIXES, required=False)` | Used for MAF and CM metadata. |
| `--ref-panel-snps-path` | input | no | scalar file path token | reference-panel SNP universe restriction | normalized, later `resolve_scalar_path` in backend/build flow | Despite `path` suffix, this is scalar exact-one. |
| `--regression-snps-path` | input | no | scalar file path token | regression SNP row-set restriction | `resolve_scalar_path(label="regression SNP list")` | Defines written LD-score row set. |
| `--keep` | input | no | scalar file path token | individual keep file | normalized and passed to backend namespace | PLINK mode only. |

Output naming:

- `--output-dir <dir>` writes a self-contained LD-score result directory.
- `manifest.json` stores relative parquet paths, baseline/query column lists,
  config metadata, chromosomes, row counts, and one count record per annotation
  column.
- `baseline.parquet` contains `CHR`, `SNP`, `BP`, `regr_weight`, then baseline
  LD-score columns.
- `query.parquet` is omitted unless query annotations exist. When present, it
  contains `CHR`, `SNP`, `BP`, then query LD-score columns.
- The new public path does not emit `.l2.ldscore.gz`, `.M`, `.M_5_50`,
  `.w.l2.ldscore.gz`, or a separate annotation manifest.

Non-path flags that affect IO interpretation:

- `--snp-identifier`, `--genome-build`: identifier and build interpretation.
- `--r2-bias-mode`, `--r2-sample-size`: parquet R2 interpretation.
- `--ld-wind-*`, `--maf`, `--chunk-size`, `--yes-really`: compute settings.

### `ldsc build-ref-panel`

Defined in `src/ldsc/ref_panel_builder.py`.

| Flag | Direction | Required | Current type | Current object | Resolution behavior | Notes |
|---|---:|---:|---|---|---|---|
| `--bfile` | input | yes | scalar PLINK prefix token | PLINK reference input | `resolve_plink_prefix_group((config.plink_prefix,), allow_chromosome_suite=True)` | Accepts prefix trios, globs, and `@` suites. |
| `--genetic-map-hg19` | input | yes | scalar/group-capable token | hg19 genetic map file group | `resolve_file_group(..., suffixes=_GENETIC_MAP_SUFFIXES, allow_chromosome_suite=True)` | May be one all-chromosome map or per-chromosome suite. |
| `--genetic-map-hg38` | input | yes | scalar/group-capable token | hg38 genetic map file group | Same as hg19 map | Same semantics. |
| `--liftover-chain-hg19-to-hg38` | input | conditional | scalar file path token | liftover chain | normalized in config; passed to `LiftOverTranslator` | Required when source build is hg19. |
| `--liftover-chain-hg38-to-hg19` | input | conditional | scalar file path token | liftover chain | normalized in config; passed to `LiftOverTranslator` | Required when source build is hg38. |
| `--ref-panel-snps-path` | input | no | scalar file path token | retained SNP universe restriction | `resolve_scalar_path(..., suffixes=_TABLE_SUFFIXES)` | Restricts retained panel SNPs before writing artifacts. |
| `--keep-indivs` | input | no | scalar file path token | individual keep file | `resolve_scalar_path(..., suffixes=_TABLE_SUFFIXES)` | Applied during PLINK loading. |
| `--out` | output | yes | scalar directory | reference-panel output root | `ensure_output_directory(config.output_dir)` | Here `--out` means directory, unlike LD-score and munging. |
| `--panel-label` | output naming | yes | string | filename stem prefix | used directly in emitted names | This is output identity, not a filesystem destination. |

Output naming:

- All outputs are written below `<out>/parquet`.
- Annotation parquet: `<out>/parquet/ann/<panel_label>_chr<chrom>_ann.parquet`.
- LD parquet: `<out>/parquet/ld/<panel_label>_chr<chrom>_LD.parquet`.
- Metadata sidecars:
  `<out>/parquet/meta/<panel_label>_chr<chrom>_meta_hg19.tsv.gz` and
  `<out>/parquet/meta/<panel_label>_chr<chrom>_meta_hg38.tsv.gz`.

### `ldsc munge-sumstats`

Parser is cloned from `src/ldsc/_kernel/sumstats_munger.py`; workflow wrapper is
`src/ldsc/sumstats_munger.py`.

| Flag | Direction | Required | Current type | Current object | Resolution behavior | Notes |
|---|---:|---:|---|---|---|---|
| `--sumstats` | input | kernel-dependent | scalar file path token | raw summary statistics | `resolve_scalar_path(label="raw sumstats")` | Exact-one path or glob. |
| `--merge-alleles` | input | no | scalar file path token | allele merge file | `resolve_scalar_path(label="merge-alleles file")` | Exact-one path or glob. |
| `--out` | output | kernel-dependent | scalar prefix | munged output prefix | `normalize_path_token`; `ensure_output_parent_directory` | Writes `<out>.sumstats.gz`; legacy kernel also uses `<out>.log`. |

Non-path flags that affect input interpretation:

- Column hints: `--snp`, `--N-col`, `--N-cas-col`, `--N-con-col`, `--a1`,
  `--a2`, `--p`, `--frq`, `--signed-sumstats`, `--info`, `--info-list`,
  `--nstudy`, `--ignore`.
- Format flags: `--daner`, `--daner-n`, `--no-alleles`, `--a1-inc`,
  `--keep-maf`.

### `ldsc h2`

Defined in `src/ldsc/regression_runner.py`.

| Flag | Direction | Required | Current type | Current object | Resolution behavior | Notes |
|---|---:|---:|---|---|---|---|
| `--ldscore-dir` | input | yes | scalar directory | canonical LD-score result directory | `load_ldscore_from_dir`; manifest paths resolved relative to the directory | Uses `baseline.parquet` only, including embedded `regr_weight`. |
| `--sumstats` | input | yes | scalar file path token | munged sumstats artifact | `load_sumstats`, then `resolve_scalar_path(label="munged sumstats")` | Exact-one path or glob. |
| `--out` | output | no | scalar prefix | summary output prefix | `normalize_path_token`; `ensure_output_parent_directory` | Writes `<out>.h2.tsv` if supplied. |

### `ldsc partitioned-h2`

Defined in `src/ldsc/regression_runner.py`.

| Flag | Direction | Required | Current type | Current object | Resolution behavior | Notes |
|---|---:|---:|---|---|---|---|
| `--ldscore-dir` | input | yes | scalar directory | canonical LD-score result directory | same as `h2` | Uses `baseline.parquet` plus `query.parquet`; loops over manifest `query_columns`. |
| `--sumstats` | input | yes | scalar file path token | munged sumstats artifact | same as `h2` | Exact-one path or glob. |
| `--out` | output | no | scalar prefix | summary output prefix | `normalize_path_token`; `ensure_output_parent_directory` | Writes `<out>.partitioned_h2.tsv` if supplied. |

### `ldsc rg`

Defined in `src/ldsc/regression_runner.py`.

| Flag | Direction | Required | Current type | Current object | Resolution behavior | Notes |
|---|---:|---:|---|---|---|---|
| `--ldscore-dir` | input | yes | scalar directory | canonical LD-score result directory | same as `h2` | Uses `baseline.parquet` only, including embedded `regr_weight`. |
| `--sumstats-1` | input | yes | scalar file path token | first munged sumstats artifact | `load_sumstats`, then `resolve_scalar_path` | Exact-one path or glob. |
| `--sumstats-2` | input | yes | scalar file path token | second munged sumstats artifact | `load_sumstats`, then `resolve_scalar_path` | Exact-one path or glob. |
| `--out` | output | no | scalar prefix | summary output prefix | `normalize_path_token`; `ensure_output_parent_directory` | Writes `<out>.rg.tsv` if supplied. |

## Public Python API Inventory

### Annotation API

Public exports:

- `AnnotationSourceSpec`
- `AnnotationBundle`
- `AnnotationBuilder`
- `run_bed_to_annot`
- `parse_bed_to_annot_args`
- `main_bed_to_annot`

| Object/function | Argument | Direction | Current type | Current object | Resolution behavior |
|---|---:|---:|---|---|---|
| `AnnotationSourceSpec` | `baseline_annot_paths` | input | path token or sequence | baseline annotation group | normalized in dataclass; resolved by `AnnotationBuilder.run` with `allow_chromosome_suite=True`. |
| `AnnotationSourceSpec` | `query_annot_paths` | input | path token or sequence | query annotation group | normalized in dataclass; resolved by `AnnotationBuilder.run` with `allow_chromosome_suite=True`. |
| `AnnotationSourceSpec` | `bed_paths` | input | path token or sequence | BED query annotation group | normalized in dataclass; resolved by `AnnotationBuilder.run` through `resolve_file_group`. |
| `AnnotationBuilder.run(source_spec, chrom=None)` | `source_spec` | input bundle spec | object | all annotation inputs | Delegates to file-group resolvers. |
| `AnnotationBuilder.run(source_spec, chrom=None)` | `chrom` | input selector | scalar | chromosome filter | Not a filesystem argument. |
| `AnnotationBuilder.parse_annotation_file(path, chrom=None)` | `path` | input | scalar file path | annotation table | Reads exact supplied path; no glob resolution in this method. |
| `AnnotationBuilder.project_bed_annotations(...)` | `bed_files` | input | path token or sequence | BED file group | `resolve_file_group(label="BED file")`. |
| `AnnotationBuilder.project_bed_annotations(...)` | `baseline_annot_paths` | input | path token or sequence | baseline annotation group | same semantics as `AnnotationSourceSpec.baseline_annot_paths`. |
| `AnnotationBuilder.project_bed_annotations(...)` | `output_dir` | output | directory path or `None` | generated query annotation directory | `ensure_output_directory`; writes `query.<chrom>.annot.gz` if provided. |
| `AnnotationBuilder.make_single_annotation_file(...)` | `bimfile` | input | scalar path token | PLINK `.bim` file | `resolve_scalar_path(label="PLINK BIM file")`. |
| `AnnotationBuilder.make_single_annotation_file(...)` | `annot_file` | output | scalar path | destination `.annot` or `.annot.gz` | `ensure_output_parent_directory`. |
| `run_bed_to_annot(...)` | `bed_files` | input | path token or sequence | BED file group | same as project method. |
| `run_bed_to_annot(...)` | `baseline_annot_paths` | input | path token or sequence | baseline annotation group | same as project method. |
| `run_bed_to_annot(...)` | `output_dir` | output | directory path or `None` | generated annotation directory | same as project method. |
| `parse_bed_to_annot_args(argv)` / `main_bed_to_annot(argv)` | `argv` | input | CLI arg vector | CLI path flags | See `ldsc annotate`. |

Related config:

| Dataclass | Field | Direction | Current object | Notes |
|---|---:|---:|---|---|
| `AnnotationBuildConfig` | `baseline_annotation_paths` | input | baseline annotation group | Normalized but not consistently used by workflow entry points. |
| `AnnotationBuildConfig` | `query_annotation_paths` | input | query annotation group | Same. |
| `AnnotationBuildConfig` | `query_bed_paths` | input | BED query group | Same. |
| `AnnotationBuildConfig` | `out_prefix` | output | output prefix | Present in config but current public BED projection uses `output_dir`. |

### Reference Panel API

Public exports:

- `ReferencePanelBuildConfig`
- `ReferencePanelBuildResult`
- `ReferencePanelBuilder`
- `run_build_ref_panel`
- `RefPanelConfig`
- `RefPanelSpec`
- `RefPanelLoader`
- `RefPanel`
- `PlinkRefPanel`
- `ParquetR2RefPanel`

| Object/function | Argument | Direction | Current type | Current object | Resolution behavior |
|---|---:|---:|---|---|---|
| `ReferencePanelBuildConfig` | `plink_prefix` | input | PLINK prefix token | source PLINK reference panel | normalized; `ReferencePanelBuilder.run` resolves with `resolve_plink_prefix_group(..., allow_chromosome_suite=True)`. |
| `ReferencePanelBuildConfig` | `genetic_map_hg19_path` | input | file/group token | hg19 genetic map | `resolve_file_group(..., allow_chromosome_suite=True)`. |
| `ReferencePanelBuildConfig` | `genetic_map_hg38_path` | input | file/group token | hg38 genetic map | same. |
| `ReferencePanelBuildConfig` | `output_dir` | output | directory | reference-panel artifact root | `ensure_output_directory`; artifacts written under `parquet/`. |
| `ReferencePanelBuildConfig` | `liftover_chain_hg19_to_hg38_path` | input | scalar path or `None` | liftover chain | normalized; required for hg19 source build. |
| `ReferencePanelBuildConfig` | `liftover_chain_hg38_to_hg19_path` | input | scalar path or `None` | liftover chain | normalized; required for hg38 source build. |
| `ReferencePanelBuildConfig` | `ref_panel_snps_path` | input | scalar path or `None` | SNP restriction file | `resolve_scalar_path(..., suffixes=_TABLE_SUFFIXES)`. |
| `ReferencePanelBuildConfig` | `keep_indivs_path` | input | scalar path or `None` | individual keep file | `resolve_scalar_path(..., suffixes=_TABLE_SUFFIXES)`. |
| `ReferencePanelBuilder.run(config)` | `config` | input/output spec | object | all build IO | Uses `ReferencePanelBuildConfig`. |
| `run_build_ref_panel(**kwargs)` | `bfile` | input | PLINK prefix token | source PLINK panel | CLI-equivalent alias for `plink_prefix`. |
| `run_build_ref_panel(**kwargs)` | `genetic_map_hg19`, `genetic_map_hg38` | input | file/group token | maps | CLI-equivalent aliases. |
| `run_build_ref_panel(**kwargs)` | `out` | output | directory | output root | CLI-equivalent alias for `output_dir`. |
| `RefPanelConfig` | `plink_prefix` | input | PLINK prefix token | PLINK backend source | normalized only; consumed by loader/specs. |
| `RefPanelConfig` | `parquet_r2_paths` | input | path token or sequence | parquet R2 group | normalized only; consumed by `RefPanelSpec`. |
| `RefPanelConfig` | `frequency_paths` | input | path token or sequence | metadata/frequency group | normalized only; consumed by `RefPanelSpec`. |
| `RefPanelSpec` | `bfile_prefix` | input | PLINK prefix token | PLINK backend source | resolved by `PlinkRefPanel` methods. |
| `RefPanelSpec` | `r2_table_paths` | input | path token or sequence | parquet R2 group | `resolve_chromosome_group(..., suffixes=PARQUET_SUFFIXES)`. |
| `RefPanelSpec` | `maf_metadata_paths` | input | path token or sequence | metadata/frequency group | `resolve_file_group` for discovery and `resolve_chromosome_group` for loads. |
| `RefPanelSpec` | `ref_panel_snps_path` | input | scalar path or `None` | SNP restriction file | normalized; `_apply_snp_restriction` reads it directly through identifier helper. |
| `PlinkRefPanel.build_reader(...)` | `keep_snps`, `keep_indivs` | input filters | in-memory collections | SNP/sample filters | Not filesystem arguments. |
| `ParquetR2RefPanel.build_reader(...)` | `metadata` | input data | DataFrame or `None` | metadata table | Not a path argument; path resolution uses spec fields. |

Output metadata:

- `ReferencePanelBuildResult.output_paths` is a dictionary keyed by artifact
  type (`ann`, `ld`, `meta_hg19`, `meta_hg38`) with per-chromosome path lists.

### LD-Score API

Public exports:

- `LDScoreConfig`
- `LDScoreOutputConfig`
- `LDScoreCalculator`
- `LDScoreResult`
- `ChromLDScoreResult`
- `run_ldscore`

| Object/function | Argument | Direction | Current type | Current object | Resolution behavior |
|---|---:|---:|---|---|---|
| `LDScoreConfig` | `keep_individuals_path` | input | scalar path or `None` | individual keep file | normalized in config; passed into PLINK backend namespace. |
| `LDScoreConfig` | `regression_snps_path` | input | scalar path or `None` | regression SNP list | normalized in config; CLI workflow resolves with `resolve_scalar_path`. |
| `LDScoreOutputConfig` | `output_dir` | output | directory | canonical LD-score result directory | normalized in config; writer creates if missing. |
| `LDScoreCalculator.run(...)` | `annotation_bundle` | input data | `AnnotationBundle` | annotation object | Not a path argument. |
| `LDScoreCalculator.run(...)` | `ref_panel` | input object | `RefPanel` | reference-panel object | Not a path argument; path resolution happened in the ref-panel object. |
| `LDScoreCalculator.run(...)` | `output_config` | output | `LDScoreOutputConfig` or `None` | canonical directory spec | If provided, `LDScoreDirectoryWriter.write` writes manifest/parquet files. |
| `LDScoreCalculator.write_outputs(result, output_config, ...)` | `output_config` | output | `LDScoreOutputConfig` | canonical directory spec | Delegates to `LDScoreDirectoryWriter`. |
| `run_ldscore(**kwargs)` | `output_dir` | output | directory | LD-score result directory | CLI-equivalent. |
| `run_ldscore(**kwargs)` | `query_annot`, `query_annot_bed`, `baseline_annot`, `bfile`, `r2_table`, `frqfile`, `ref_panel_snps_path`, `regression_snps_path`, `keep` | input | CLI-equivalent tokens | same objects as CLI | Same as `ldsc ldscore`. |

Output metadata:

- `LDScoreResult.output_paths` and `ChromLDScoreResult.output_paths` carry
  written artifact paths when outputs are produced.
- `LDScoreResult.baseline_table`, `query_table`, and `count_records` are the
  public in-memory equivalents of the directory payload.

### Output Layer API

Public exports:

- `OutputSpec`
- `LDScoreOutputConfig`
- `ArtifactConfig`
- `RunSummary`
- `ArtifactProducer`
- `LDScoreDirectoryWriter`
- `ResultFormatter`
- `ResultWriter`
- `OutputManager`
- `PostProcessor`

| Object/function | Argument | Direction | Current type | Current object | Resolution behavior |
|---|---:|---:|---|---|---|
| `OutputSpec` | `out_prefix` | output | scalar path-like prefix | artifact filename stem | normalized; only basename is used by `_artifact_filename`. |
| `OutputSpec` | `output_dir` | output | directory or `None` | artifact root | normalized; `None` means current directory. |
| `OutputSpec` | `artifact_layout` | output layout | enum string | layout rule | `flat`, `by_chrom`, or `run_dir`. |
| `OutputSpec` | `log_path` | output | scalar path or `None` | log destination | normalized, but no current writer uses it directly. |
| `OutputSpec` | `enabled_artifacts` | output selection | list or `None` | selected artifact producers | Not a path argument. |
| `LDScoreOutputConfig` | `output_dir` | output | directory | canonical LD-score result directory | normalized; writer creates if missing and preflights fixed filenames. |
| `Artifact` | `relative_path` | output | relative path string | artifact destination under root | Written by `ResultWriter`. |
| `LDScoreDirectoryWriter.write(result, output_config)` | `output_config` | output | `LDScoreOutputConfig` | LD-score directory writer config | Writes fixed `manifest.json`, `baseline.parquet`, and optional `query.parquet`. |
| `ResultWriter.resolve_path(artifact, root)` | `root` | output | `Path` | output root directory | Creates parent directory for `root / artifact.relative_path`. |
| `ResultWriter.prepare_writes(artifacts, root, overwrite=False)` | `root` | output | `Path` | output root directory | Preflights all output paths. |
| `ResultWriter.write(artifact, root, overwrite=False)` | `root` | output | `Path` | output root directory | Writes one artifact. |
| `ResultWriter.write_prepared(artifact, path)` | `path` | output | concrete `Path` | exact artifact destination | Assumes conflict checks already happened. |
| `OutputManager.write_outputs(result, output_spec, ...)` | `output_spec` | output | `OutputSpec` | output artifact spec | Creates root and writes enabled artifacts. |

Output naming rules:

- `artifact_layout="flat"`: `<prefix><suffix>[.gz]`.
- `artifact_layout="by_chrom"`: `chr<chrom>/<prefix>.<chrom><suffix>[.gz]`.
- `artifact_layout="run_dir"`: root becomes `<output_dir>/<basename(out_prefix)>`;
  artifact names still use basename prefix.

### Summary Statistics API

Public exports:

- `RawSumstatsSpec`
- `MungeConfig`
- `SumstatsMunger`
- `SumstatsTable`
- `MungeRunSummary`
- `load_sumstats`

| Object/function | Argument | Direction | Current type | Current object | Resolution behavior |
|---|---:|---:|---|---|---|
| `RawSumstatsSpec` | `path` | input | scalar path token | raw summary-statistics file | normalized in dataclass; `SumstatsMunger.run` resolves exact-one. |
| `MungeConfig` | `out_prefix` | output | scalar prefix | munged output prefix | normalized; `SumstatsMunger` writes `<out_prefix>.sumstats.gz` and tracks `<out_prefix>.log`. |
| `MungeConfig` | `merge_alleles_path` | input | scalar path or `None` | allele merge file | normalized; `SumstatsMunger` resolves exact-one. |
| `load_sumstats(path, trait_name=None)` | `path` | input | scalar path token | munged sumstats artifact | `resolve_scalar_path(label="munged sumstats")`. |
| `SumstatsMunger.run(raw_source, munge_config, ...)` | `raw_source.path` | input | scalar path token | raw summary-statistics file | `resolve_scalar_path(label="raw sumstats")`. |
| `SumstatsMunger.run(raw_source, munge_config, ...)` | `munge_config.out_prefix` | output | scalar prefix | munged artifact prefix | `ensure_output_parent_directory`; returned table is in memory but summary records output paths. |
| `SumstatsMunger.write_output(sumstats, out_prefix)` | `out_prefix` | output | scalar prefix | munged artifact prefix | Writes `<out_prefix>.sumstats.gz`. |
| `sumstats_munger.main(argv)` | `argv` | input | CLI arg vector | CLI path flags | See `ldsc munge-sumstats`. |

Output metadata:

- `MungeRunSummary.output_paths` records `sumstats_gz` and `log`.
- `SumstatsTable.source_path` records the resolved source path.

### Regression API

Public exports:

- `RegressionConfig`
- `RegressionDataset`
- `RegressionRunner`
- `load_ldscore_from_dir`

Most `RegressionRunner` methods consume in-memory `SumstatsTable`,
`LDScoreResult`, and `RegressionDataset` objects. Filesystem IO is concentrated
in CLI wrappers and `load_ldscore_from_dir`.

| Object/function | Argument | Direction | Current type | Current object | Resolution behavior |
|---|---:|---:|---|---|---|
| `load_ldscore_from_dir(...)` | `ldscore_dir` | input | scalar directory | LD-score result directory | Reads `manifest.json`; parquet file paths are resolved relative to the directory. |
| `run_h2_from_args(args)` | `args.ldscore_dir`, `args.sumstats` | input | directory plus scalar path token | h2 artifacts | CLI wrapper around loaders. |
| `run_h2_from_args(args)` | `args.out` | output | prefix or `None` | h2 summary | `_maybe_write_dataframe(..., ".h2.tsv")`. |
| `run_partitioned_h2_from_args(args)` | `args.ldscore_dir`, `args.sumstats` | input | directory plus scalar path token | partitioned h2 artifacts | CLI wrapper around loaders; requires query columns in the manifest. |
| `run_partitioned_h2_from_args(args)` | `args.out` | output | prefix or `None` | partitioned h2 summary | Writes `.partitioned_h2.tsv`. |
| `run_rg_from_args(args)` | `args.ldscore_dir`, `args.sumstats_1`, `args.sumstats_2` | input | directory plus scalar path tokens | rg artifacts | CLI wrapper around loaders. |
| `run_rg_from_args(args)` | `args.out` | output | prefix or `None` | rg summary | Writes `.rg.tsv`. |

### Global Config, Column Inference, and Build Inference

These public modules expose no user-facing filesystem output arguments.

| Object/function | IO-relevant argument | Direction | Notes |
|---|---:|---:|---|
| `GlobalConfig`, `set_global_config`, `get_global_config`, `reset_global_config` | none | none | Control identifier/build/log semantics only. |
| `ColumnSpec` and column inference helpers | none | none | Operate on in-memory column names. |
| `load_packaged_reference_table()` | none | input resource | Loads bundled package data, not a user-supplied path. |
| `resolve_chr_pos_table(df, ...)`, `infer_chr_pos_build(df, ...)` | none | none | Operate on in-memory data frames. |
| `normalize_chromosome`, `chrom_sort_key` | none | none | No filesystem IO. |

## Current Inconsistencies Relevant to Refactor

1. `--out` has different meanings across commands.
   - `munge-sumstats`, `h2`, `partitioned-h2`, and `rg`: output prefix.
   - `build-ref-panel`: output directory.
   - `ldscore` now uses `--output-dir`.

2. The package has both `path` and `paths` names for group-capable tokens.
   - `ref_panel_snps_path` and `regression_snps_path` are scalar exact-one
     paths.
   - `baseline_annot_paths`, `query_annot_paths`, `r2_table_paths`, and
     `maf_metadata_paths` are group-capable tokens that may include globs.

3. One logical object is sometimes split across multiple arguments.
   - Reference panel input can be `--bfile` plus `--frqfile`, or `--r2-table`
     plus `--frqfile`, plus optional `--r2-bias-*`.
   - Reference-panel build output is `--out` plus `--panel-label`; together
     they define the object's output identity.
   - Regression output summaries still use `--out` prefixes.

4. `prefix` is overloaded.
   - PLINK input prefix is a true input object identifier because it expands to
     `.bed/.bim/.fam`.
   - Output prefix is a legacy filename-stem convention.
   - `OutputSpec.out_prefix` can be path-like but only the basename is used
     when `output_dir` is supplied.

5. Directory-vs-prefix behavior is inconsistent at the CLI boundary.
   - `build-ref-panel --out` is a directory and `--panel-label` supplies the
     file stem.
   - `ldscore --output-dir` is a literal result directory.
   - Regression commands optionally write summary files only when a prefix is
     supplied.

6. Glob behavior is strong for inputs but intentionally absent for outputs.
   - Input scalar paths may be exact-one globs.
   - Input groups may be globs and `@` suites.
   - Output destinations are literal paths/directories; they should remain
     literal in the refactor.

## Candidate Target Vocabulary

This is not the full design, but it follows directly from the inventory.

| Target term | Use when | Examples |
|---|---|---|
| `*_path` | One file-like input. Input may allow exact-one glob. | `sumstats_path`, `keep_indivs_path`, `regression_snps_path`. |
| `*_paths` | A logical group of file inputs that may resolve to many files through globs or `@`. | `annotation_paths`, `r2_paths`, `metadata_paths`, `bed_paths`. |
| `*_dir` | A directory. All public workflow outputs should use directories because runs also need logs, metadata, and sidecar artifacts. Output directories are literal and created if missing. | `output_dir`, `ldscore_dir`. |

The main cleanup pressure points are `--out`, annotation input naming,
reference-panel source specification, and the relation between `output_dir`,
`out_prefix`, and `panel_label`.

## Accepted IO Refactor Rules

The target public API should use these rules for both Python and CLI surfaces:

1. Use `*_path` for one file-like input.
2. Use `*_paths` when one logical input may resolve to multiple files via
   globs or explicit `@` chromosome-suite tokens.
3. Use `*_dir` for directories.
4. Use `--output-dir` for workflow outputs that produce result directories.
   The LD-score workflow now follows this rule; regression summary output is a
   remaining prefix-style surface.
5. Avoid public output prefixes. If a workflow needs a stable run or panel
   identity, derive it from the basename of `output_dir` and record it in run
   metadata. If a file needs a name, use a canonical filename inside
   `output_dir`.
6. Use `--plink-path` for PLINK input instead of `--bfile`.
7. Rename annotation BED inputs to `query_annot_bed_paths` in Python and
   `--query-annot-bed-paths` in the CLI.
8. Rename LD-score PLINK sample filters to `keep_indivs_path` in Python and
   `--keep-indivs-path` in the CLI, unless the project intentionally accepts
   `--keep-indivs` as the single exception.

Under these rules, multi-artifact workflows should converge on `output_dir`.
`annotate` and `ldscore` now expose output directories. Regression commands
still use optional summary prefixes and are a remaining cleanup target.
