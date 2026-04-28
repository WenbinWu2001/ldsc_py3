# IO Argument Inventory and Implementation Plan

Date: 2026-04-28

This document records the final public input/output naming target after the
LD-score result-directory refactor. The LD-score workflow now uses a canonical
result directory as the baseline design:

```text
<ldscore_dir>/
  manifest.json
  baseline.parquet
  query.parquet        # omitted when no query annotations exist
```

Regression workflows consume this directory with `--ldscore-dir`; fragmented
inputs such as LD-score files, count vectors, regression-weight files, and
annotation manifests are no longer public inputs.

## Naming Rules

Public CLI flags and Python config fields follow these rules:

| Suffix | Meaning | Examples |
|---|---|---|
| `*_path` | one file-like input, exact-one glob allowed where the resolver supports it | `sumstats_path`, `merge_alleles_path`, `keep_indivs_path` |
| `*_paths` | one logical input that may resolve to many files via globs, comma lists, or `@` chromosome tokens | `baseline_annot_paths`, `query_annot_bed_paths`, `r2_paths` |
| `*_dir` | directory input or output location | `ldscore_dir`, `output_dir` |

Removed from public surfaces:

- `prefix`, `out`, `out_prefix`, `panel_label`
- `bfile`, `frqfile`
- legacy aliases for annotation BEDs, LD-score components, and regression
  outputs

Users customize run identity by choosing the `output_dir` name. Output filenames
inside that directory are fixed and workflow-specific.

## CLI Inventory

### `ldsc annotate`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--query-annot-bed-paths` | input | yes | BED interval files | Accepts exact files, globs, comma-separated tokens, and path-token lists. BED basenames become query annotation names. |
| `--baseline-annot-paths` | input | yes | baseline `.annot[.gz]` templates | Accepts exact files, globs, and `@` chromosome-suite tokens. |
| `--output-dir` | output | yes | generated query annotation directory | Writes `query.<chrom>.annot.gz` files or compatibility batch outputs. |

Removed flags: `--bed-files`, `--baseline-annot`.

### `ldsc ldscore`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--output-dir` | output | yes | canonical LD-score result directory | Writes `manifest.json`, `baseline.parquet`, and optional `query.parquet`. |
| `--baseline-annot-paths` | input | yes | baseline annotation files | Exact files, globs, comma lists, or `@` suites. |
| `--query-annot-paths` | input | no | prebuilt query annotation files | Mutually exclusive with `--query-annot-bed-paths`. |
| `--query-annot-bed-paths` | input | no | query BED interval files | Projected in memory onto the baseline SNP universe. |
| `--plink-path` | input | conditional | PLINK reference panel prefix | Required when not using `--r2-paths`; supports exact prefix, PLINK-prefix glob, or `@` suite. |
| `--r2-paths` | input | conditional | parquet R2 files | Required when not using `--plink-path`; supports exact files, globs, comma lists, or `@` suites. |
| `--metadata-paths` | input | no | reference metadata or frequency sidecars | Used for MAF and centiMorgan metadata where needed. |
| `--ref-panel-snps-path` | input | no | reference-panel SNP universe restriction | Scalar file-like input. |
| `--regression-snps-path` | input | no | persisted LD-score row-set restriction | Scalar file-like input. |
| `--keep-indivs-path` | input | no | PLINK individual keep file | PLINK mode only. |

Removed flags: `--bfile`, `--r2-table`, `--frqfile`, `--keep`,
`--baseline-annot`, `--query-annot`, `--query-annot-bed`, legacy `--out`.

LD-score output schema:

- `baseline.parquet`: `CHR`, `SNP`, `BP`, `regr_weight`, then baseline
  LD-score columns.
- `query.parquet`: `CHR`, `SNP`, `BP`, then query LD-score columns; omitted
  when there are no query annotations.
- `manifest.json`: format version, relative file paths, baseline/query column
  names, count records, config metadata, chromosomes, row counts, and writer
  metadata.

### `ldsc build-ref-panel`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--plink-path` | input | yes | PLINK reference panel prefix | Supports exact prefix, PLINK-prefix glob, or `@` suite. |
| `--source-genome-build` | input metadata | yes | source build selector | Not a filesystem argument. |
| `--genetic-map-hg19-path` | input | yes | hg19 genetic map file or suite | Exact file, exact-one glob, or `@` suite. |
| `--genetic-map-hg38-path` | input | yes | hg38 genetic map file or suite | Exact file, exact-one glob, or `@` suite. |
| `--liftover-chain-hg19-to-hg38-path` | input | conditional | liftover chain file | Required when source build is hg19. |
| `--liftover-chain-hg38-to-hg19-path` | input | conditional | liftover chain file | Required when source build is hg38. |
| `--ref-panel-snps-path` | input | no | retained SNP universe restriction | Scalar file-like input. |
| `--keep-indivs-path` | input | no | PLINK individual keep file | Applied during PLINK loading. |
| `--output-dir` | output | yes | reference-panel artifact directory | Run identity is `Path(output_dir).name`; no separate label is accepted. |

Removed flags: `--bfile`, `--out`, `--panel-label`, `--keep-indivs`,
`--genetic-map-hg19`, `--genetic-map-hg38`, old liftover-chain names without
`_path`.

Fixed output names:

```text
<output_dir>/parquet/ann/chr{chrom}_ann.parquet
<output_dir>/parquet/ld/chr{chrom}_LD.parquet
<output_dir>/parquet/meta/chr{chrom}_meta_hg19.tsv.gz
<output_dir>/parquet/meta/chr{chrom}_meta_hg38.tsv.gz
```

### `ldsc munge-sumstats`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--sumstats-path` | input | yes | raw summary-statistics file | Exact path or exact-one glob. |
| `--merge-alleles-path` | input | no | allele merge file | Exact path or exact-one glob. |
| `--output-dir` | output | yes | munged output directory | Internally uses `<output_dir>/sumstats` as the legacy kernel stem. |

Removed flags: `--sumstats`, `--merge-alleles`, `--out`.

Fixed output names:

```text
<output_dir>/sumstats.sumstats.gz
<output_dir>/sumstats.log
```

### `ldsc h2`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--ldscore-dir` | input | yes | canonical LD-score result directory | Reads baseline LD scores and embedded `regr_weight`. |
| `--sumstats-path` | input | yes | munged summary-statistics file | Exact path or exact-one glob. |
| `--output-dir` | output | no | result output directory | Writes `h2.tsv` when supplied. |

Removed flags: `--ldscore`, `--counts`, `--w-ld`, `--annotation-manifest`,
`--sumstats`, `--out`.

### `ldsc partitioned-h2`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--ldscore-dir` | input | yes | canonical LD-score result directory | Reads baseline plus query LD scores. |
| `--sumstats-path` | input | yes | munged summary-statistics file | Exact path or exact-one glob. |
| `--output-dir` | output | no | result output directory | Writes `partitioned_h2.tsv` when supplied. |

Removed flags: `--ldscore`, `--counts`, `--w-ld`, `--annotation-manifest`,
`--query-columns`, `--sumstats`, `--out`.

### `ldsc rg`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--ldscore-dir` | input | yes | canonical LD-score result directory | Reads baseline LD scores and embedded `regr_weight`. |
| `--sumstats-1-path` | input | yes | first munged summary-statistics file | Exact path or exact-one glob. |
| `--sumstats-2-path` | input | yes | second munged summary-statistics file | Exact path or exact-one glob. |
| `--output-dir` | output | no | result output directory | Writes `rg.tsv` when supplied. |

Removed flags: `--ldscore`, `--counts`, `--w-ld`, `--annotation-manifest`,
`--sumstats-1`, `--sumstats-2`, `--out`.

## Public Python API Inventory

### Annotation

| Object/function | Argument | Direction | Object |
|---|---:|---:|---|
| `AnnotationBuildConfig` | `baseline_annot_paths` | input | baseline annotation group |
| `AnnotationBuildConfig` | `query_annot_paths` | input | prebuilt query annotation group |
| `AnnotationBuildConfig` | `query_annot_bed_paths` | input | query BED group |
| `AnnotationBuildConfig` | `output_dir` | output | generated query annotation directory |
| `AnnotationBuilder.run(config=None, chrom=None)` | `config` | input/output | annotation workflow config; defaults to the builder config |
| `AnnotationBuilder.project_bed_annotations(...)` | `query_annot_bed_paths` | input | query BED group |
| `run_bed_to_annot(...)` | `query_annot_bed_paths` | input | query BED group |
| `run_bed_to_annot(...)` | `baseline_annot_paths` | input | baseline annotation templates |
| `run_bed_to_annot(...)` | `output_dir` | output | generated query annotation directory |

Removed Python names: `bed_paths`, `query_bed_paths`, `bed_files`,
`baseline_annot`, `out_prefix`.

### LD-score calculation

| Object/function | Argument | Direction | Object |
|---|---:|---:|---|
| `RefPanelConfig` | `plink_path` | input | PLINK reference-panel prefix token |
| `RefPanelConfig` | `r2_paths` | input | parquet R2 group |
| `RefPanelConfig` | `metadata_paths` | input | metadata/frequency sidecar group |
| `RefPanelConfig` | `ref_panel_snps_path` | input | reference SNP restriction |
| `LDScoreConfig` | `regression_snps_path` | input | regression SNP restriction |
| `LDScoreConfig` | `keep_indivs_path` | input | PLINK individual keep file |
| `LDScoreOutputConfig` | `output_dir` | output | canonical LD-score result directory |
| `run_ldscore(**kwargs)` | `baseline_annot_paths`, `query_annot_paths`, `query_annot_bed_paths` | input | annotation sources |
| `run_ldscore(**kwargs)` | `plink_path`, `r2_paths`, `metadata_paths` | input | reference-panel sources |
| `run_ldscore(**kwargs)` | `output_dir` | output | canonical result directory |

Removed Python names: `bfile`, `r2_table`, `frqfile`, `keep`,
`baseline_annot`, `query_annot`, `query_annot_bed`, `out`.

### Reference-panel building

| Object/function | Argument | Direction | Object |
|---|---:|---:|---|
| `ReferencePanelBuildConfig` | `plink_path` | input | PLINK reference-panel prefix token |
| `ReferencePanelBuildConfig` | `genetic_map_hg19_path` | input | hg19 genetic map |
| `ReferencePanelBuildConfig` | `genetic_map_hg38_path` | input | hg38 genetic map |
| `ReferencePanelBuildConfig` | `liftover_chain_hg19_to_hg38_path` | input | liftover chain |
| `ReferencePanelBuildConfig` | `liftover_chain_hg38_to_hg19_path` | input | liftover chain |
| `ReferencePanelBuildConfig` | `ref_panel_snps_path` | input | retained SNP restriction |
| `ReferencePanelBuildConfig` | `keep_indivs_path` | input | PLINK individual keep file |
| `ReferencePanelBuildConfig` | `output_dir` | output | artifact directory |
| `run_build_ref_panel(**kwargs)` | same config field names | input/output | CLI-equivalent wrapper |

Removed Python names: `plink_prefix`, `bfile`, `out`, `panel_label`,
`keep_indivs`, old genetic-map and liftover names without `_path`.

### Sumstats munging

| Object/function | Argument | Direction | Object |
|---|---:|---:|---|
| `MungeConfig` | `sumstats_path` | input | raw summary-statistics file |
| `MungeConfig` | `trait_name` | input metadata | optional trait label |
| `MungeConfig` | `column_hints` | input metadata | optional source-column hints |
| `MungeConfig` | `merge_alleles_path` | input | merge-alleles file |
| `MungeConfig` | `output_dir` | output | munged output directory |
| `SumstatsMunger.run(munge_config, ...)` | `munge_config` | input/output | normalized munging workflow |
| `SumstatsMunger.write_output(sumstats, output_dir)` | `output_dir` | output | writes fixed `sumstats.sumstats.gz` |

Removed Python names: legacy separate source-path object field,
`MungeConfig.out_prefix`, `write_output(..., out_prefix)`.

### Regression

| Object/function | Argument | Direction | Object |
|---|---:|---:|---|
| `run_h2_from_args(args)` | `ldscore_dir` | input | LD-score result directory |
| `run_h2_from_args(args)` | `sumstats_path` | input | munged summary-statistics file |
| `run_h2_from_args(args)` | `output_dir` | output | writes `h2.tsv` |
| `run_partitioned_h2_from_args(args)` | `ldscore_dir` | input | LD-score result directory |
| `run_partitioned_h2_from_args(args)` | `sumstats_path` | input | munged summary-statistics file |
| `run_partitioned_h2_from_args(args)` | `output_dir` | output | writes `partitioned_h2.tsv` |
| `run_rg_from_args(args)` | `ldscore_dir` | input | LD-score result directory |
| `run_rg_from_args(args)` | `sumstats_1_path`, `sumstats_2_path` | input | munged summary-statistics files |
| `run_rg_from_args(args)` | `output_dir` | output | writes `rg.tsv` |

Removed Python/public argparse names: `sumstats`, `sumstats_1`, `sumstats_2`,
`out`, `ldscore`, `counts`, `w_ld`, `annotation_manifest`, `query_columns`.

## Remaining Implementation Checklist

- [x] LD-score calculation writes the canonical result directory.
- [x] Regression commands read `--ldscore-dir` instead of fragmented LD-score
  inputs.
- [x] CLI parsers reject old public flags with abbreviation disabled on the
  unified surfaces.
- [x] Public dataclasses use `*_path`, `*_paths`, and `output_dir` names.
- [x] Annotation BED input is named `query_annot_bed_paths`.
- [x] PLINK input is named `plink_path` / `--plink-path`.
- [x] `keep` inputs are named `keep_indivs_path` / `--keep-indivs-path`.
- [x] Munging writes fixed files under `output_dir`.
- [x] Regression writes fixed TSV files under `output_dir`.
- [x] Build-ref-panel no longer accepts a separate panel label; output identity
  comes from the directory name.
- [x] Removed the low-level `ldsc.outputs.ArtifactOutputConfig.out_prefix`
  compatibility pipeline from the public output module; `ldsc.outputs` now exposes
  `LDScoreOutputConfig` and `LDScoreDirectoryWriter` for LD-score output writing.

## Clarifications and Recommendations

- Keep `--plink-path` rather than `--plink-prefix`. Although PLINK technically
  uses a file stem, the public rule is clearer if every filesystem input is
  visibly a path-like token.
- Do not add `--output-name` or `--panel-name`. Fixed output names make
  downstream automation simpler; users name runs by naming the output
  directory.
- Keep `--output-dir` required for artifact-writing workflows. Regression may
  continue to allow it as optional only when returning in-memory results is a
  supported Python/test path.
- Treat `_kernel` legacy namespace names (`bfile`, `frqfile`, `r2_table`,
  `keep`) as private adapter details until the numerical kernel itself is
  rewritten.
