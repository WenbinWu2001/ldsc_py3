# IO Argument Inventory

Date: 2026-04-28

This document records the current public input/output naming contract after the
LD-score result-directory refactor. The LD-score workflow uses a canonical
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
| `*_file` | one file-like input, exact-one glob allowed where the resolver supports it | `sumstats_file`, `sumstats_snps_file`, `keep_indivs_file` |
| `*_sources` | one logical input that may resolve to many files via globs, comma lists, or `@` chromosome tokens | `baseline_annot_sources`, `query_annot_bed_sources` |
| `*_dir` | directory input or output location | `ldscore_dir`, `output_dir` |

Removed from public surfaces:

- generic output/run identity names: `prefix`, `out`, `out_prefix`,
  `panel_label`
- `bfile`, `frqfile`
- legacy aliases for annotation BEDs, LD-score components, and regression
  outputs

Users customize run identity by choosing the `output_dir` name. Output filenames
inside that directory are fixed and workflow-specific.

## CLI Inventory

### `ldsc annotate`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--query-annot-bed-sources` | input | yes | BED interval files | Accepts exact files, globs, comma-separated tokens, and source-token lists. BED basenames become query annotation names. |
| `--baseline-annot-sources` | input | yes | baseline `.annot[.gz]` templates | Accepts exact files, globs, and `@` chromosome-suite tokens. |
| `--output-dir` | output | yes | generated query annotation directory | Writes combined `query.<chrom>.annot.gz` files, with all BED inputs represented as query columns. |
| `--overwrite` | output mode | no | collision policy | Permits replacement of generated `query.<chrom>.annot.gz` files. |
| `--snp-identifier`, `--genome-build` | config | no | coordinate interpretation | `--snp-identifier` defaults to `chr_pos`; `chr_pos` coordinate inputs need `--genome-build auto`, `hg19`, or `hg38`. |

Removed flags: `--bed-files`, `--baseline-annot`.

### `ldsc ldscore`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--output-dir` | output | yes | canonical LD-score result directory | Writes `manifest.json`, `baseline.parquet`, and optional `query.parquet`; parquet row groups are chromosome-aligned. |
| `--overwrite` | output mode | no | collision policy | Permits replacement of the fixed LD-score files in `output_dir`. |
| `--baseline-annot-sources` | input | no | baseline annotation files | Exact files, globs, comma lists, or `@` suites. If omitted with no query inputs, `ldscore` synthesizes an all-ones `base` column from retained reference-panel metadata. |
| `--query-annot-sources` | input | no | prebuilt query annotation files | Mutually exclusive with `--query-annot-bed-sources`; requires `--baseline-annot-sources`. |
| `--query-annot-bed-sources` | input | no | query BED interval files | Projected in memory onto the baseline SNP universe; requires `--baseline-annot-sources`. |
| `--plink-prefix` | input | conditional | PLINK reference panel prefix | Required when not using parquet reference-panel input; supports exact prefix, PLINK-prefix glob, or `@` suite. |
| `--r2-dir` | input | conditional | package-built parquet R2 directory | Parquet input; pass a build-specific directory such as `ref_panel/hg38`, containing `chr*_r2.parquet` plus optional `chr*_meta.tsv.gz`. |
| `--r2-bias-mode` | input metadata | conditional | parquet R2 bias declaration | Required with `--r2-dir`; choose `raw` or `unbiased`. |
| `--r2-sample-size` | input metadata | conditional | parquet R2 sample size | Required only when `--r2-bias-mode raw`. |
| `--ref-panel-snps-file` | input | no | reference-panel SNP universe restriction | Scalar file-like input. |
| `--regression-snps-file` | input | no | persisted LD-score row-set restriction | Scalar file-like input. |
| `--keep-indivs-file` | input | no | PLINK individual keep file | PLINK mode only. |
| `--maf-min` | input metadata | no | retained reference-panel MAF filter | Applied in the reference-panel layer when MAF is available. |
| `--common-maf-min` | input metadata | no | common-SNP count threshold | Default `0.05`; affects common count vectors only, using `MAF >= common_maf_min`. |
| `--chunk-size` | performance | no | PLINK block size | Default `128`; used by legacy PLINK block computations. |
| `--yes-really` | safety override | no | whole-chromosome LD windows | Allows whole-chromosome LD windows when no LD-window option is supplied. |

Removed flags: `--bfile`, `--r2-table`, `--frqfile`, `--r2-sources`,
`--metadata-sources`, `--keep`, `--maf`, `--baseline-annot`,
`--query-annot`, `--query-annot-bed`, legacy `--out`.

LD-score output schema:

- `baseline.parquet`: `CHR`, `POS`, `SNP`, `regr_weight`, then baseline
  LD-score columns. In no-annotation unpartitioned runs, the baseline column
  list is exactly `base`.
- `query.parquet`: `CHR`, `POS`, `SNP`, then query LD-score columns; omitted
  when there are no query annotations.
- `manifest.json`: format version, relative file paths, baseline/query column
  names, count records, `count_config`, config metadata, chromosomes, row
  counts, `row_group_layout`, `baseline_row_groups`, and
  `query_row_groups`.

### `ldsc build-ref-panel`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--plink-prefix` | input | yes | PLINK reference panel prefix | Supports exact prefix, PLINK-prefix glob, or `@` suite. |
| `--source-genome-build` | input metadata | no | source PLINK coordinate build | Not a filesystem argument; inferred from `.bim` before SNP restriction when omitted. |
| `--genetic-map-hg19-sources` | input | conditional | hg19 genetic map file or suite | Required when `--ld-wind-cm` is used and hg19 output is emitted; otherwise optional and missing hg19 CM values are written as `NA`. |
| `--genetic-map-hg38-sources` | input | conditional | hg38 genetic map file or suite | Required when `--ld-wind-cm` is used and hg38 output is emitted; otherwise optional and missing hg38 CM values are written as `NA`. |
| `--liftover-chain-hg19-to-hg38-file` | input | no | liftover chain file | Optional; enables hg38 outputs for hg19 source builds. |
| `--liftover-chain-hg38-to-hg19-file` | input | no | liftover chain file | Optional; enables hg19 outputs for hg38 source builds. |
| `--ref-panel-snps-file` | input | no | retained SNP universe restriction | Scalar file-like input; identifier mode comes from `GlobalConfig.snp_identifier`; `chr_pos` coordinates must align to the source PLINK build. |
| `--snp-identifier` | input metadata | conditional | global SNP identifier mode | CLI-only way to construct `GlobalConfig.snp_identifier` for this invocation. |
| `--keep-indivs-file` | input | no | PLINK individual keep file | Applied during PLINK loading. |
| `--maf-min` | input metadata | no | retained SNP MAF filter | Applied during PLINK loading. |
| `--output-dir` | output | yes | reference-panel artifact directory | Run identity is `Path(output_dir).name`; no separate label is accepted. |
| `--overwrite` | output mode | no | collision policy | Permits replacement of deterministic parquet and metadata outputs. |
| `--chunk-size` | performance | no | block size | Default `128`; used during block processing. |

Removed flags: `--bfile`, `--out`, `--panel-label`, `--keep-indivs`, `--maf`,
`--genetic-map-hg19`, `--genetic-map-hg38`, old liftover-chain names without
`_file` / `_sources`.

Fixed output names:

```text
<output_dir>/hg19/chr{chrom}_r2.parquet
<output_dir>/hg19/chr{chrom}_meta.tsv.gz
<output_dir>/hg38/chr{chrom}_r2.parquet
<output_dir>/hg38/chr{chrom}_meta.tsv.gz
```

When no usable source-to-target liftover chain is provided, the builder logs
an INFO message and emits source-build-only outputs. When a matching liftover
chain is provided, it emits both source and target build R2/metadata trees.

When SNP- or kb-window builds omit a genetic map for an emitted build, that
metadata sidecar is still written with `CM=NA`. cM-window builds require the
genetic map for every emitted build because each build's map defines that
build's LD window.

### `ldsc munge-sumstats`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--sumstats-file` | input | yes | raw summary-statistics file | Exact path or exact-one glob. |
| `--sumstats-snps-file` | input | no | summary-statistics SNP keep-list | Exact path or exact-one glob. Restricts munged rows only; does not allele-match or reorder rows. |
| `--output-dir` | output | yes | munged output directory | Internally uses `<output_dir>/sumstats` as the legacy kernel stem. |
| `--chr`, `--pos` | input metadata | no | raw column hints | Explicit chromosome and base-pair position columns; common aliases such as `#CHROM`, `CHROM`, `CHR`, `POS`, and `BP` are also inferred. |
| `--snp-identifier`, `--genome-build` | config | no | provenance | `--snp-identifier` defaults to `chr_pos`; `--genome-build` defaults to `hg38`; `--genome-build auto` can infer hg19/hg38 for complete `CHR`/`POS` rows. |
| `--overwrite` | output mode | no | collision policy | Permits replacement of the fixed sumstats outputs. |

Removed flags: `--sumstats`, `--merge-alleles`, `--merge-alleles-file`, `--out`.

Fixed output names:

```text
<output_dir>/sumstats.sumstats.gz
<output_dir>/sumstats.log
<output_dir>/sumstats.metadata.json
```

### `ldsc h2`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--ldscore-dir` | input | yes | canonical LD-score result directory | Reads baseline LD scores and embedded `regr_weight`. |
| `--sumstats-file` | input | yes | munged summary-statistics file | Exact path or exact-one glob. |
| `--output-dir` | output | no | result output directory | Writes `h2.tsv` when supplied. |
| `--count-kind` | model | no | count vector choice | `common` by default; `all` uses all-SNP counts. |
| `--overwrite` | output mode | no | collision policy | Permits replacement of `h2.tsv`. |

Removed flags: `--ldscore`, `--counts`, `--w-ld`, `--annotation-manifest`,
`--sumstats`, `--out`.

### `ldsc partitioned-h2`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--ldscore-dir` | input | yes | canonical LD-score result directory | Reads baseline plus query LD scores. |
| `--sumstats-file` | input | yes | munged summary-statistics file | Exact path or exact-one glob. |
| `--output-dir` | output | no | result output directory | Writes `partitioned_h2.tsv` when supplied. |
| `--count-kind` | model | no | count vector choice | `common` by default; `all` uses all-SNP counts. |
| `--write-per-query-results` | output mode | no | per-query result tree | Also writes `query_annotations/manifest.tsv` plus one sanitized query folder under `--output-dir`. |
| `--overwrite` | output mode | no | collision policy | Permits replacement of aggregate and requested per-query outputs. |

Removed flags: `--ldscore`, `--counts`, `--w-ld`, `--annotation-manifest`,
`--query-columns`, `--sumstats`, `--out`.

### `ldsc rg`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--ldscore-dir` | input | yes | canonical LD-score result directory | Reads baseline LD scores and embedded `regr_weight`. |
| `--sumstats-1-file` | input | yes | first munged summary-statistics file | Exact path or exact-one glob. |
| `--sumstats-2-file` | input | yes | second munged summary-statistics file | Exact path or exact-one glob. |
| `--output-dir` | output | no | result output directory | Writes `rg.tsv` when supplied. |
| `--count-kind` | model | no | count vector choice | `common` by default; `all` uses all-SNP counts. |
| `--overwrite` | output mode | no | collision policy | Permits replacement of `rg.tsv`. |

Removed flags: `--ldscore`, `--counts`, `--w-ld`, `--annotation-manifest`,
`--sumstats-1`, `--sumstats-2`, `--out`.

## Public Python API Inventory

### Annotation

| Object/function | Argument | Direction | Object |
|---|---:|---:|---|
| `AnnotationBuildConfig` | `baseline_annot_sources` | input | baseline annotation group |
| `AnnotationBuildConfig` | `query_annot_sources` | input | prebuilt query annotation group |
| `AnnotationBuildConfig` | `query_annot_bed_sources` | input | query BED group |
| `AnnotationBuildConfig` | `output_dir` | output | generated query annotation directory |
| `AnnotationBuilder.run(config=None, chrom=None)` | `config` | input/output | annotation workflow config; defaults to the builder config |
| `AnnotationBuilder.project_bed_annotations(...)` | `query_annot_bed_sources` | input | query BED group |
| `run_bed_to_annot(...)` | `query_annot_bed_sources` | input | query BED group |
| `run_bed_to_annot(...)` | `baseline_annot_sources` | input | baseline annotation templates |
| `run_bed_to_annot(...)` | `output_dir` | output | generated query annotation directory |

Removed Python names: `bed_paths`, `query_bed_paths`, `bed_files`,
`baseline_annot`, `out_prefix`.

### LD-score calculation

| Object/function | Argument | Direction | Object |
|---|---:|---:|---|
| `RefPanelConfig` | `plink_prefix` | input | PLINK reference-panel prefix token |
| `RefPanelConfig` | `r2_dir` | input | preferred package-built parquet reference-panel directory |
| `RefPanelConfig` | `ref_panel_snps_file` | input | reference SNP restriction |
| `RefPanelConfig` | `keep_indivs_file` | input | PLINK individual keep file |
| `RefPanelConfig` | `maf_min` | input metadata | retained reference-panel MAF filter |
| `LDScoreConfig` | `regression_snps_file` | input | regression SNP restriction |
| `LDScoreConfig` | `common_maf_min` | input metadata | common-SNP count threshold only |
| `LDScoreOutputConfig` | `output_dir` | output | canonical LD-score result directory |
| `run_ldscore(**kwargs)` | `baseline_annot_sources`, `query_annot_sources`, `query_annot_bed_sources` | input | optional annotation sources; query inputs require baseline sources, and no-annotation runs synthesize `base` |
| `run_ldscore(**kwargs)` | `plink_prefix`, `r2_dir` | input | reference-panel sources |
| `run_ldscore(**kwargs)` | `output_dir` | output | canonical result directory |

Removed Python names: `bfile`, `r2_table`, `frqfile`, `keep`, `maf`,
`baseline_annot`, `query_annot`, `query_annot_bed`, `out`,
`r2_ref_panel_dir`, `ref_panel_dir`, `r2_sources`, and `metadata_sources`.

### Reference-panel building

| Object/function | Argument | Direction | Object |
|---|---:|---:|---|
| `ReferencePanelBuildConfig` | `plink_prefix` | input | PLINK reference-panel prefix token |
| `ReferencePanelBuildConfig` | `source_genome_build` | input metadata | optional PLINK source build; inferred from `.bim` when omitted |
| `ReferencePanelBuildConfig` | `genetic_map_hg19_sources` | input | conditional hg19 genetic map |
| `ReferencePanelBuildConfig` | `genetic_map_hg38_sources` | input | conditional hg38 genetic map |
| `ReferencePanelBuildConfig` | `liftover_chain_hg19_to_hg38_file` | input | optional liftover chain |
| `ReferencePanelBuildConfig` | `liftover_chain_hg38_to_hg19_file` | input | optional liftover chain |
| `ReferencePanelBuildConfig` | `ref_panel_snps_file` | input | retained SNP restriction |
| `ReferencePanelBuildConfig` | `keep_indivs_file` | input | PLINK individual keep file |
| `ReferencePanelBuildConfig` | `maf_min` | input metadata | retained SNP MAF filter |
| `ReferencePanelBuildConfig` | `output_dir` | output | artifact directory |
| `run_build_ref_panel(**kwargs)` | same config field names except global settings | input/output | convenience wrapper; reads `snp_identifier` from the registered `GlobalConfig`; ignores `GlobalConfig.genome_build` |

Removed Python names: `plink_path`, `bfile`, `out`, `panel_label`,
`keep_indivs`, `maf`, old genetic-map and liftover names without `_file` /
`_sources`.

### Sumstats munging

| Object/function | Argument | Direction | Object |
|---|---:|---:|---|
| `MungeConfig` | `sumstats_file` | input | raw summary-statistics file |
| `MungeConfig` | `trait_name` | input metadata | optional trait label |
| `MungeConfig` | `column_hints` | input metadata | optional source-column hints |
| `MungeConfig` | `sumstats_snps_file` | input | summary-statistics SNP keep-list |
| `MungeConfig` | `output_dir` | output | munged output directory |
| `SumstatsMunger.run(munge_config, ...)` | `munge_config` | input/output | normalized munging workflow |
| `SumstatsMunger.write_output(sumstats, output_dir)` | `output_dir` | output | writes fixed `sumstats.sumstats.gz` |

Removed Python names: legacy separate source-path object field,
`MungeConfig.out_prefix`, `write_output(..., out_prefix)`.

### Regression

| Object/function | Argument | Direction | Object |
|---|---:|---:|---|
| `run_h2_from_args(args)` | `ldscore_dir` | input | LD-score result directory |
| `run_h2_from_args(args)` | `sumstats_file` | input | munged summary-statistics file |
| `run_h2_from_args(args)` | `output_dir` | output | writes `h2.tsv` |
| `run_partitioned_h2_from_args(args)` | `ldscore_dir` | input | LD-score result directory |
| `run_partitioned_h2_from_args(args)` | `sumstats_file` | input | munged summary-statistics file |
| `run_partitioned_h2_from_args(args)` | `output_dir` | output | writes `partitioned_h2.tsv` |
| `run_partitioned_h2_from_args(args)` | `write_per_query_results` | output mode | optionally writes `query_annotations/` |
| `run_rg_from_args(args)` | `ldscore_dir` | input | LD-score result directory |
| `run_rg_from_args(args)` | `sumstats_1_file`, `sumstats_2_file` | input | munged summary-statistics files |
| `run_rg_from_args(args)` | `output_dir` | output | writes `rg.tsv` |

Removed Python/public argparse names: `sumstats`, `sumstats_1`, `sumstats_2`,
`out`, `ldscore`, `counts`, `w_ld`, `annotation_manifest`, `query_columns`.

Current curated `.sumstats.gz` artifacts provide canonical `SNP`, `CHR`, `POS`,
`Z`, and `N` fields when written by `ldsc munge-sumstats`, plus a neighboring
metadata sidecar with the effective `snp_identifier` and `genome_build`.
Regression therefore merges on literal `SNP` in `rsid` mode and on normalized
`CHR:POS` coordinates in `chr_pos` mode.

## Remaining Implementation Checklist

- [x] LD-score calculation writes the canonical result directory.
- [x] Regression commands read `--ldscore-dir` instead of fragmented LD-score
  inputs.
- [x] CLI parsers reject old public flags with abbreviation disabled on the
  unified surfaces.
- [x] Public dataclasses use `*_file`, `*_sources`, and `output_dir` names.
- [x] Annotation BED input is named `query_annot_bed_sources`.
- [x] PLINK input is named `plink_prefix` / `--plink-prefix`.
- [x] `keep` inputs are named `keep_indivs_file` / `--keep-indivs-file`.
- [x] Munging writes fixed files under `output_dir`.
- [x] Regression writes fixed TSV files under `output_dir`; partitioned-h2 can
  additionally write per-query folders when requested.
- [x] Build-ref-panel no longer accepts a separate panel label; output identity
  comes from the directory name.
- [x] Removed the old prefix-based output compatibility pipeline from the
  public output module; `ldsc.outputs` now exposes
  `LDScoreOutputConfig`, `LDScoreDirectoryWriter`,
  `PartitionedH2OutputConfig`, and `PartitionedH2DirectoryWriter` for output
  writing.

## Clarifications and Recommendations

- Keep `--plink-prefix` rather than forcing it into the `*_file` or
  `*_sources` suffix rules. PLINK inputs are a file stem for a related
  `.bed`/`.bim`/`.fam` trio, so `prefix` is the most accurate public name.
- Prefer `--r2-dir` / `r2_dir` for package-built parquet
  panels. It mirrors `--plink-prefix` as the single reference-panel input for
  one backend. Direct R2-file and metadata-source groups are no longer part of
  the LD-score public API; standard panels should use the fixed directory
  layout with metadata sidecars next to their R2 parquet files.
- Do not add `--output-name` or `--panel-name`. Fixed output names make
  downstream automation simpler; users name runs by naming the output
  directory.
- Keep `--output-dir` required for artifact-writing workflows. Regression may
  continue to allow it as optional only when returning in-memory results is a
  supported Python/test path.
- Treat `_kernel` legacy namespace names (`bfile`, `frqfile`, `r2_table`,
  `keep`) as private adapter details until the numerical kernel itself is
  rewritten.
