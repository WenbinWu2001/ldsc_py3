# IO Argument Inventory

Date: 2026-04-28

This document records the current public input/output naming contract after the
LD-score result-directory refactor. The LD-score workflow uses a canonical
result directory as the baseline design:

```text
<ldscore_dir>/
  metadata.json
  ldscore.baseline.parquet
  ldscore.query.parquet        # omitted when no query annotations exist
  diagnostics/
    ldscore.log
```

Regression workflows consume this directory with `--ldscore-dir`; fragmented
inputs such as LD-score files, count vectors, regression-weight files, and
annotation manifests are no longer public inputs.

## Shared SNP Identity Contract

Public `--snp-identifier` values are exactly `rsid`, `rsid_allele_aware`,
`chr_pos`, and `chr_pos_allele_aware`; the default is
`chr_pos_allele_aware`. Mode names are exact. Column aliases apply only to
input headers.

Base modes are allele-blind: `rsid` uses only `SNP`, and `chr_pos` uses only
`CHR:POS`. Allele columns may be preserved as passive data, but they do not
affect identity, duplicate filtering, retention, or drop reasons. Allele-aware
modes require usable `A1/A2` for sumstats, reference-panel artifacts, R2 parquet
endpoints, and LD-score artifacts. Artifact duplicate filtering always computes
the effective merge key for the active mode, then drops all rows in duplicate
clusters.

Restriction files may omit alleles and then match by base key. Allele-bearing
restrictions, including the packaged HM3 map, match by the effective
allele-aware key in allele-aware modes. Annotation files may omit alleles in
allele-aware modes because they describe genomic membership; if
annotation alleles are present, they participate in allele-aware matching.
External raw R2 parquet inputs are supported only in `rsid` and `chr_pos`;
allele-aware modes require package-built canonical R2 parquet with
`A1_1/A2_1/A1_2/A2_2`.

## Format Policy

The package separates machine-consumed internal artifacts from
science-facing result tables. 
Internal artifacts are machine-consumed intermediate products and should use Parquet as their primary format. This applies to package-built R2 pair tables,
canonical LD-score result tables, and curated munged summary statistics. JSON
metadata is a downstream contract only for sumstats and LD-score artifacts;
non-consumed workflow metadata belongs under `diagnostics/`.

Science-facing result tables should use TSV. The public regression results are
`h2.tsv`, `partitioned_h2.tsv`, optional diagnostic per-query
`partitioned_h2.tsv` / `partitioned_h2_full.tsv`, and the rg family
(`rg.tsv`, `rg_full.tsv`, `h2_per_trait.tsv`, optional diagnostic `pairs/`).
Logs are audit files and should not be treated as output results.


### Current Adaptation Status

Adapted public paths:

- `ldsc ldscore` writes `ldscore.baseline.parquet` and optional
  `ldscore.query.parquet`.
- `ldsc munge-sumstats` writes `sumstats.parquet` by default.
- `ldsc h2`, `ldsc partitioned-h2`, and `ldsc rg` write result tables as TSV
  when an `output_dir` is supplied.
- `ldsc build-ref-panel` writes the primary R2 pair matrix as
  `chr{chrom}_r2.parquet`.

Not fully adapted or retained for compatibility:

- `ldsc annotate` writes generated query annotations as
  `query.<chrom>.annot.gz`; these are intermediate annotation artifacts rather
  than science-facing result tables.
- `ldsc build-ref-panel` still writes runtime metadata sidecars as
  `chr{chrom}_meta.tsv.gz` and liftover-stage drop audit files as
  `dropped_snps/chr{chrom}_dropped.tsv.gz`.
- `ldsc munge-sumstats` keeps `--output-format tsv.gz` and `both`, which write
  `sumstats.sumstats.gz` compatibility artifacts even though parquet is the
  default.
- Private `_kernel` emitters still write legacy `.sumstats.gz`,
  `.l2.ldscore.gz`, `.w.l2.ldscore.gz`, `.l2.M`, `.l2.M_5_50`, and
  `.annotation_groups.tsv` files for direct kernel or compatibility paths.

## Naming Rules

Public CLI flags and Python config fields follow these rules:

| Suffix | Meaning | Examples |
|---|---|---|
| `*_file` | one file-like input, exact-one glob allowed where the resolver supports it | `raw_sumstats_file`, `sumstats_file`, `sumstats_snps_file`, `keep_indivs_file` |
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

Workflow logs are fixed audit files under `output_dir`, preflighted with the
scientific outputs before they are opened. They are not included in workflow
`output_paths` mappings or thin metadata sidecars that downstream code
interprets as data artifacts.

For `munge-sumstats`, `ldscore`, `partitioned-h2`, and `annotate`, the fixed
outputs are coherent artifact families. Without overwrite, any existing owned
artifact in the family rejects the run, even when that artifact would not be
written by the current configuration. With overwrite, the workflow writes the
requested outputs and removes stale owned siblings not produced by the
successful run. Unrelated files in `output_dir` are preserved.

## CLI Inventory

### `ldsc annotate`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--query-annot-bed-sources` | input | yes | BED interval files | Accepts exact files, globs, comma-separated tokens, and source-token lists. BED basenames become query annotation names. |
| `--baseline-annot-sources` | input | yes | baseline `.annot[.gz]` templates | Accepts exact files, globs, and `@` chromosome-suite tokens. |
| `--output-dir` | output | yes | generated query annotation directory | Writes combined root `query.<chrom>.annot.gz` files, with all BED inputs represented as query columns, plus diagnostic `metadata.json`, `dropped_snps/dropped.tsv.gz`, and `annotate.log` under `diagnostics/`. |
| `--overwrite` | output mode | no | collision policy | Controls whether generated annotation files and diagnostics may be replaced; defaults to `False`, so any existing root-level `query.*.annot.gz` shard or owned diagnostic artifact is refused. With overwrite, stale query shards outside the current chromosome set are removed after a successful run. |
| `--log-level` | logging | no | workflow log verbosity | Controls ordinary LDSC logger records in console and `diagnostics/annotate.log`; lifecycle audit lines always appear in the file. |
| `--snp-identifier`, `--genome-build` | config | no | coordinate interpretation | Define how SNP coordinates are interpreted; `--snp-identifier` defaults to `chr_pos_allele_aware`, and `--genome-build` defaults to omitted/`None`, which is invalid for coordinate-family inputs unless the workflow can infer it. |

Removed flags: `--bed-files`, `--baseline-annot`.

### `ldsc ldscore`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--output-dir` | output | yes | canonical LD-score result directory | Writes root `metadata.json`, `ldscore.baseline.parquet`, optional `ldscore.query.parquet`, and `diagnostics/ldscore.log`; parquet row groups are chromosome-aligned. |
| `--overwrite` | output mode | no | collision policy | Controls whether fixed LD-score files and `diagnostics/ldscore.log` may be replaced; defaults to `False`, so any existing owned LD-score artifact in `output_dir` is refused. With overwrite, stale `ldscore.query.parquet` is removed after successful baseline-only runs. |
| `--log-level` | logging | no | workflow log verbosity | Controls ordinary LDSC logger records in console and `diagnostics/ldscore.log`; lifecycle audit lines always appear in the file. |
| `--baseline-annot-sources` | input | no | baseline annotation files | Supplies baseline annotation files; defaults to omitted/`None`, and if no query inputs are supplied `ldscore` synthesizes an all-ones `base` column. |
| `--query-annot-sources` | input | no | prebuilt query annotation files | Supplies prebuilt query annotation files; defaults to omitted/`None`, so no prebuilt query annotations are used. Mutually exclusive with `--query-annot-bed-sources` and requires `--baseline-annot-sources`. |
| `--query-annot-bed-sources` | input | no | query BED interval files | Supplies BED intervals to project as query annotations; defaults to omitted/`None`, so no BED query annotations are projected. Requires `--baseline-annot-sources`. |
| `--plink-prefix` | input | conditional | PLINK reference panel prefix | Selects PLINK reference-panel input; defaults to omitted/`None` and is required when `--r2-dir` is omitted. Supports exact prefix, PLINK-prefix glob, or `@` suite. |
| `--r2-dir` | input | conditional | package-built parquet R2 directory | Selects parquet reference-panel input; defaults to omitted/`None` and is required when `--plink-prefix` is omitted. Use a build-specific directory such as `ref_panel/hg38`. Allele-aware modes require package-built canonical R2 parquet with endpoint allele columns. |
| `--r2-bias-mode` | input metadata | optional | parquet R2 bias declaration | Declares whether parquet R2 values are raw or unbiased. Package-built panels auto-load this from `ldsc:r2_bias`; legacy files without metadata still default to `unbiased`. Choose `raw` only for external raw sample R2 values. |
| `--r2-sample-size` | input metadata | conditional | parquet R2 sample size | Provides the sample size for correcting raw R2. Package-built raw panels can auto-load this from `ldsc:n_samples`; legacy raw files still require an explicit value with `--r2-bias-mode raw`. |
| `--ref-panel-snps-file` | input | no | reference-panel SNP universe restriction | Restricts the retained reference-panel SNP universe using identity keys only; duplicate restriction keys collapse to one retained key, and non-identity columns such as `CM` or `MAF` are ignored. Defaults to omitted/`None`, so no additional restriction is applied. |
| `--use-hm3-ref-panel-snps` | input mode | no | packaged HM3 reference-panel SNP restriction | Restricts the retained reference-panel SNP universe to the packaged curated HM3 map. Mutually exclusive with `--ref-panel-snps-file`. |
| `--regression-snps-file` | input | no | persisted LD-score row-set restriction | Restricts the written LD-score row set using identity keys only; duplicate restriction keys collapse to one retained key, and non-identity columns such as `CM` or `MAF` are ignored. Defaults to omitted/`None`, so rows are not restricted by a persisted regression SNP set. |
| `--use-hm3-regression-snps` | input mode | no | packaged HM3 regression SNP restriction | Restricts the persisted LD-score row set to the packaged curated HM3 map. Mutually exclusive with `--regression-snps-file`. |
| `--keep-indivs-file` | input | no | PLINK individual keep file | Restricts PLINK individuals before LD calculation; defaults to omitted/`None`, so no individual keep filter is applied. PLINK mode only. |
| `--maf-min` | input metadata | no | retained reference-panel MAF filter | Filters retained reference-panel SNPs by MAF; defaults to omitted/`None`, so no retained-reference MAF filter is applied. |
| `--common-maf-min` | input metadata | no | common-SNP count threshold | Sets the MAF threshold for common-SNP count vectors; defaults to `0.05` and uses `MAF >= common_maf_min`. |
| `--snp-batch-size` | performance | no | LD-score SNP batch size | Number of SNPs processed per LD-score sliding batch; defaults to `128`. The canonical parquet reader sizes its decoded row-group cache automatically from this value and the chromosome LD window. |
| `--yes-really` | safety override | no | whole-chromosome LD windows | Allows whole-chromosome LD windows; defaults to `False`, so such windows are rejected unless this flag is supplied. |

Removed flags: `--bfile`, `--r2-table`, `--frqfile`, `--r2-sources`,
`--metadata-sources`, `--keep`, `--maf`, `--baseline-annot`,
`--query-annot`, `--query-annot-bed`, `--chunk-size`, legacy `--out`.

LD-score output schema:

- `ldscore.baseline.parquet`: `CHR`, `POS`, `SNP`, `regression_ld_scores`, then baseline
  LD-score columns. In no-annotation unpartitioned runs, the baseline column
  list is exactly `base`.
- `ldscore.query.parquet`: `CHR`, `POS`, `SNP`, then query LD-score columns; omitted
  when there are no query annotations.
- `metadata.json`: `schema_version`, `artifact_type`, relative file paths,
  baseline/query column names, count records, `count_config`, config metadata,
  chromosomes, row
  counts, `row_group_layout`, `baseline_row_groups`, and
  `query_row_groups`.

### `ldsc build-ref-panel`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--plink-prefix` | input | yes | PLINK reference panel prefix | Supports exact prefix, PLINK-prefix glob, or `@` suite. |
| `--source-genome-build` | input metadata | no | source PLINK coordinate build | Declares the PLINK coordinate build; defaults to omitted/`None`, so the build is inferred from `.bim` before SNP restriction. |
| `--genetic-map-hg19-sources` | input | conditional | hg19 genetic map file or suite | Supplies hg19 genetic-map CM values; defaults to omitted/`None` and is required when `--ld-wind-cm` is used and hg19 output is emitted. |
| `--genetic-map-hg38-sources` | input | conditional | hg38 genetic map file or suite | Supplies hg38 genetic-map CM values; defaults to omitted/`None` and is required when `--ld-wind-cm` is used and hg38 output is emitted. |
| `--liftover-chain-hg19-to-hg38-file` | input | no | liftover chain file | Enables hg38 outputs for hg19 source builds; defaults to omitted/`None`, so those target-build outputs are not emitted. |
| `--liftover-chain-hg38-to-hg19-file` | input | no | liftover chain file | Enables hg19 outputs for hg38 source builds; defaults to omitted/`None`, so those target-build outputs are not emitted. |
| `--ref-panel-snps-file` | input | no | retained SNP universe restriction | Restricts the emitted reference-panel SNP universe using identity keys only; duplicate restriction keys collapse to one retained key, and non-identity columns such as `CM` or `MAF` are ignored. Defaults to omitted/`None`, so no retained SNP restriction is applied. |
| `--use-hm3-snps` | input mode | no | packaged HM3 SNP restriction | Restricts the emitted reference-panel SNP universe to the packaged curated HM3 map. Mutually exclusive with `--ref-panel-snps-file`. |
| `--use-hm3-quick-liftover` | input mode | no | packaged HM3 coordinate map | Emits source and opposite-build R2/metadata for the HM3-restricted coordinate universe. Requires `--use-hm3-snps`, is valid only in `chr_pos`-family modes, and is mutually exclusive with chain-file liftover. |
| `--snp-identifier` | input metadata | conditional | global SNP identifier mode | Overrides the registered SNP identifier mode for this invocation; defaults to omitted/`None`, so `GlobalConfig.snp_identifier` is used. |
| `--keep-indivs-file` | input | no | PLINK individual keep file | Restricts PLINK individuals during panel building; defaults to omitted/`None`, so no individual keep filter is applied. |
| `--maf-min` | input metadata | no | retained SNP MAF filter | Filters retained SNPs by MAF during PLINK loading; defaults to omitted/`None`, so no retained-SNP MAF filter is applied. |
| `--output-dir` | output | yes | reference-panel artifact directory | Run identity is `Path(output_dir).name`; no separate label is accepted. |
| `--overwrite` | output mode | no | collision policy | Controls whether reference-panel artifacts, diagnostic metadata, always-written `diagnostics/dropped_snps/chr{chrom}_dropped.tsv.gz` audit files, and build-ref-panel workflow logs may be replaced; defaults to `False`, so existing owned outputs are refused. With overwrite enabled, stale owned target-build, out-of-scope chromosome, dropped-SNP, or log siblings are removed after a successful run. |
| `--log-level` | logging | no | workflow log verbosity | Controls ordinary LDSC logger records in console and the build-ref-panel workflow log; lifecycle audit lines always appear in the file. |
| `--snp-batch-size` | performance | no | SNP computation batch size | Number of SNPs loaded per pairwise-R2 computation batch; larger values may improve throughput but use more memory. Defaults to `128`. |

Removed flags: `--bfile`, `--out`, `--panel-label`, `--keep-indivs`, `--maf`,
`--genetic-map-hg19`, `--genetic-map-hg38`, old liftover-chain names without
`_file` / `_sources`, `--duplicate-position-policy`.

Fixed output names:

```text
<output_dir>/hg19/chr{chrom}_r2.parquet
<output_dir>/hg19/chr{chrom}_meta.tsv.gz
<output_dir>/hg38/chr{chrom}_r2.parquet
<output_dir>/hg38/chr{chrom}_meta.tsv.gz
<output_dir>/diagnostics/metadata.json
<output_dir>/diagnostics/dropped_snps/chr{chrom}_dropped.tsv.gz
<output_dir>/diagnostics/build-ref-panel.log
<output_dir>/diagnostics/build-ref-panel.chr<chrom>.log  # concrete single-chromosome PLINK prefix
```

Each `chr{chrom}_r2.parquet` stores Arrow schema metadata for
`ldsc:sorted_by_build`, `ldsc:row_group_size`, `ldsc:n_samples`, and
`ldsc:r2_bias`. Current package-built panels write unbiased R2 values and record
the PLINK sample count, so downstream LD-score runs can omit `--r2-bias-mode`
and `--r2-sample-size` for panels produced by this codebase.

When no usable source-to-target liftover method is provided, the builder logs an
INFO message and emits source-build-only outputs. When a matching liftover chain
is provided, it emits both source and target build R2/metadata trees. HM3 quick
liftover emits the same source and target tree shapes for the HM3-restricted
coordinate universe, backed by the packaged map rather than a chain file.
Reference-panel liftover is rejected when the active SNP identifier mode is
`rsid`-family modes; use a `chr_pos`-family mode for cross-build coordinate emission or omit liftover
for source-build-only rsID panels.

When SNP- or kb-window builds omit a genetic map for an emitted build, that
metadata sidecar is still written with `CM=NA`. cM-window builds require the
genetic map for every emitted build because each build's map defines that
build's LD window.

Use a fresh `build-ref-panel` output directory when changing emitted builds,
liftover/coordinate configuration, or chromosome scope.

### `ldsc munge-sumstats`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--raw-sumstats-file` | input | yes | raw summary-statistics file | Exact path or exact-one glob. |
| `--format` | input metadata | no | raw summary-statistics format profile | One of `auto`, `plain`, `daner-old`, `daner-new`, or `pgc-vcf`; defaults to `auto`, which detects common plain text, old DANER, new DANER, and PGC VCF-style headers. Conflicts with incompatible legacy DANER flags are rejected. |
| `--infer-only` | diagnostic | no | raw summary-statistics inference report | Reads the raw header and first data row, prints detected format, inferred hints, missing fields, notes, and a suggested minimal command. Missing `A1/A2` is reported only in allele-aware modes. Does not require `--output-dir` and writes no artifacts. |
| `--sumstats-snps-file` | input | no | summary-statistics SNP keep-list | Restricts munged summary-statistics rows using identity keys only; duplicate restriction keys collapse to one retained key, and non-identity columns such as `CM` or `MAF` are ignored. The keep-list is loaded before parsing and applied while chunks are streaming; defaults to omitted/`None`, so no keep-list restriction is applied. |
| `--use-hm3-snps` | input mode | no | packaged HM3 SNP restriction | Restricts munged summary-statistics rows to the packaged curated HM3 map while chunks are streaming. Mutually exclusive with `--sumstats-snps-file`. |
| `--trait-name` | input metadata | no | biological trait label | Optional label stored in root `metadata.json`; downstream regression uses it unless a regression CLI `--trait-name` override is supplied. |
| `--output-dir` | output | yes, except `--infer-only` | munged output directory | The workflow writes fixed `sumstats.*` artifacts under this directory and passes `<output_dir>/sumstats` as the kernel output stem. |
| `--output-format` | output mode | no | curated sumstats format | One of `parquet`, `tsv.gz`, or `both`; defaults to `parquet`. |
| `--chr`, `--pos` | input metadata | no | raw column hints | Identify raw chromosome and position columns; default to omitted/`None`, so common aliases such as `#CHROM`, `CHROM`, `CHR`, `POS`, and `BP` are inferred. |
| `--target-genome-build` | input metadata | no | optional liftover target build | Enables source-to-target coordinate liftover when paired with exactly one liftover method; valid only in `chr_pos`-family modes and never rewrites `SNP`. |
| `--liftover-chain-file` | input | no | optional munger liftover chain | Uses a source-to-target chain file for coordinate-only sumstats liftover. Mutually exclusive with `--use-hm3-quick-liftover`. |
| `--use-hm3-quick-liftover` | input mode | no | packaged HM3 coordinate map | Uses the curated dual-build HM3 map for HM3-only quick liftover. Requires `--use-hm3-snps` and is mutually exclusive with `--liftover-chain-file`. |
| `--daner-old`, `--daner-new` | input metadata | no | DANER schema interpretation | `--daner-old` parses case/control N from `FRQ_A_<Ncas>` and `FRQ_U_<Ncon>` headers; `--daner-new` parses exact `Nca` and `Nco` columns. |
| `--snp-identifier`, `--genome-build` | config | no | provenance | `--snp-identifier` defaults to `chr_pos_allele_aware`; `--genome-build` defaults to `hg38`; `--genome-build auto` can infer hg19/hg38 for complete `CHR`/`POS` rows. Allele-aware modes require usable `A1/A2`; rerun with `--snp-identifier chr_pos` or `--snp-identifier rsid` to run without allele-aware identity. |
| `--log-level` | logging | no | workflow log verbosity | Controls ordinary LDSC logger records in console and `diagnostics/sumstats.log`; lifecycle audit lines always appear in the file. |
| `--overwrite` | output mode | no | collision policy | Controls whether fixed sumstats outputs may be replaced; defaults to `False`, so any owned `sumstats.*` artifact is refused. With overwrite, stale sibling formats not produced by the current `--output-format` are removed after a successful run. |

Removed flags: `--sumstats`, `--sumstats-file` for raw munge input,
`--merge-alleles`, `--merge-alleles-file`, `--no-alleles`, `--out`.

Fixed output names:

```text
<output_dir>/sumstats.sumstats.gz
<output_dir>/sumstats.parquet
<output_dir>/metadata.json
<output_dir>/diagnostics/sumstats.log
<output_dir>/diagnostics/dropped_snps/dropped.tsv.gz
```

`sumstats.parquet` is the default curated artifact. `sumstats.sumstats.gz` is
written only for `--output-format tsv.gz` or `both`. `diagnostics/sumstats.log` is
preflighted and opened by the public workflow layer, but it is excluded from
`MungeRunSummary.output_paths`. Detailed provenance and output bookkeeping are
written to `diagnostics/sumstats.log`; row-level liftover drops are written to
`diagnostics/dropped_snps/dropped.tsv.gz`; root `metadata.json` stays limited
to the downstream contract. The kernel emits package logger records for QC progress
and preserves its direct legacy-compatible `.sumstats.gz` writer for
private/direct kernel calls.

Column-name flags should usually be omitted. The workflow infers safe aliases
for common text tables and format profiles, including `EA`/`EFFECT_ALLELE` as
`A1`, `NEA`/`NON_EFFECT_ALLELE` as `A2`, `PVAL` as `P`, and `IMPINFO` as
`INFO`. `A1` is the allele that the signed statistic is relative to, not the
genome reference allele; `A2` is its counterpart. Positive signed statistics
are interpreted relative to `A1`. `NEFF` is not treated as `N` automatically;
use `--N-col NEFF` only when that is analytically appropriate. Comma-separated
numeric/NA INFO values such as `IMPINFO=0.852,0.113,0.842,0.88,NA` are filtered
on the mean of numeric non-missing tokens; mixed values such as
`IMPINFO=0.95,LOW,0.88` raise with suggestions to use `--ignore` or explicit
INFO flags.

### `ldsc h2`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--ldscore-dir` | input | yes | canonical LD-score result directory | Reads baseline LD scores and embedded `regression_ld_scores`, the historical `w_ld` component used when final h2 weights are computed. |
| `--sumstats-file` | input | yes | munged summary-statistics file | Exact path or exact-one glob. |
| `--trait-name` | input metadata | no | output trait label | Optional label override. If omitted, regression uses sumstats root `metadata.json["trait_name"]` when present, then the filename fallback. |
| `--output-dir` | output | no | result output directory | Selects where to write h2 results; defaults to omitted/`None`, so the CLI prints the compact `h2.tsv` schema to stdout and writes no files. |
| `--count-kind` | model | no | count vector choice | Selects the count vector used by regression; defaults to `common`, while `all` uses all-SNP counts. |
| `--log-level` | logging | no | workflow log verbosity | Controls ordinary LDSC logger records in console and `diagnostics/h2.log` when `output_dir` is supplied; lifecycle audit lines always appear in the file. |
| `--overwrite` | output mode | no | collision policy | Controls whether `h2.tsv` and diagnostics may be replaced; defaults to `False`, so an existing owned artifact is refused. |

Removed flags: `--ldscore`, `--counts`, `--w-ld`, `--annotation-manifest`,
`--sumstats`, `--out`.

### `ldsc partitioned-h2`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--ldscore-dir` | input | yes | canonical LD-score result directory | Must contain baseline plus query LD scores; baseline-only directories are rejected. |
| `--sumstats-file` | input | yes | munged summary-statistics file | Exact path or exact-one glob. |
| `--trait-name` | input metadata | no | output trait label | Optional label override. If omitted, regression uses sumstats root `metadata.json["trait_name"]` when present, then the filename fallback. |
| `--output-dir` | output | no | result output directory | Selects where to write partitioned-h2 results; defaults to omitted/`None`, so the CLI prints the compact `partitioned_h2.tsv` schema to stdout and writes no files. |
| `--count-kind` | model | no | count vector choice | Selects the count vector used by regression; defaults to `common`, while `all` uses all-SNP counts. |
| `--write-per-query-results` | output mode | no | per-query result tree | Requests per-query output folders under `diagnostics/query_annotations`; defaults to `False`, so only the aggregate table is returned/written. |
| `--log-level` | logging | no | workflow log verbosity | Controls ordinary LDSC logger records in console and `diagnostics/partitioned-h2.log` when `output_dir` is supplied; lifecycle audit lines always appear in the file. |
| `--overwrite` | output mode | no | collision policy | Controls whether aggregate/per-query outputs and diagnostics may be replaced; defaults to `False`, so any owned partitioned-h2 artifact is refused. With overwrite, aggregate-only runs remove stale `diagnostics/query_annotations/` trees after successful writes. |

Removed flags: `--ldscore`, `--counts`, `--w-ld`, `--annotation-manifest`,
`--query-columns`, `--sumstats`, `--out`.

### `ldsc rg`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--ldscore-dir` | input | yes | canonical LD-score result directory | Reads baseline LD scores and embedded `regression_ld_scores`, the historical `w_ld` component used when final rg/gencov weights are computed. |
| `--sumstats-sources` | input | yes | two or more munged summary-statistics files | Accepts exact paths and glob patterns. With two files, computes one pair; with three or more files and no anchor, computes all unordered pairs in input order. |
| `--anchor-trait` | input selector | no | anchor trait label or path | When supplied, first matches a resolved trait name, then a resolved input path; computes anchor-vs-rest pairs only. |
| `--output-dir` | output | no | result output directory | Selects where to write the rg output family. Without it, Python returns `RgResultFamily` and the CLI prints only the concise `rg.tsv` schema to stdout. |
| `--write-per-pair-detail` | output mode | no | optional pair result tree | Requires `--output-dir`; writes `diagnostics/pairs/manifest.tsv` plus one `rg_full.tsv` and `metadata.json` per attempted pair. |
| `--count-kind` | model | no | count vector choice | Selects the count vector used by regression; defaults to `common`, while `all` uses all-SNP counts. |
| `--log-level` | logging | no | workflow log verbosity | Controls ordinary LDSC logger records in console and `diagnostics/rg.log` when `output_dir` is supplied; lifecycle audit lines always appear in the file. |
| `--overwrite` | output mode | no | collision policy | Controls whether `rg.tsv`, `rg_full.tsv`, `h2_per_trait.tsv`, optional `diagnostics/pairs/`, and diagnostics may be replaced; defaults to `False`, so an existing owned artifact is refused. |

Removed flags: `--ldscore`, `--counts`, `--w-ld`, `--annotation-manifest`,
`--sumstats-1`, `--sumstats-2`, `--sumstats-1-file`, `--sumstats-2-file`,
`--trait-name-1`, `--trait-name-2`, `--anchor-trait-file`, `--out`.

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
| `add_annotate_arguments(parser)` | `parser` | CLI surface | shared annotate argument registration for standalone and top-level parsers |
| `run_annotate_from_args(args)` / `main(argv)` | `query_annot_bed_sources` | input | query BED group |
| `run_annotate_from_args(args)` / `main(argv)` | `baseline_annot_sources` | input | baseline annotation templates |
| `run_annotate_from_args(args)` / `main(argv)` | `output_dir` | output | generated query annotation directory |
| `run_bed_to_annot(...)` | `query_annot_bed_sources` | input | query BED group |
| `run_bed_to_annot(...)` | `baseline_annot_sources` | input | baseline annotation templates |
| `run_bed_to_annot(...)` | `output_dir` | output | generated query annotation directory; convenience wrapper writes `diagnostics/annotate.log` |

Removed Python names: `bed_paths`, `query_bed_paths`, `bed_files`,
`baseline_annot`, `out_prefix`, `main_bed_to_annot`.

### LD-score calculation

| Object/function | Argument | Direction | Object |
|---|---:|---:|---|
| `RefPanelConfig` | `plink_prefix` | input | PLINK reference-panel prefix token |
| `RefPanelConfig` | `r2_dir` | input | preferred package-built parquet reference-panel directory |
| `RefPanelConfig` | `ref_panel_snps_file` | input | reference SNP restriction |
| `RefPanelConfig` | `use_hm3_ref_panel_snps` | input mode | packaged HM3 reference-panel SNP restriction |
| `RefPanelConfig` | `keep_indivs_file` | input | PLINK individual keep file |
| `RefPanelConfig` | `maf_min` | input metadata | retained reference-panel MAF filter |
| `LDScoreConfig` | `regression_snps_file` | input | regression SNP restriction |
| `LDScoreConfig` | `use_hm3_regression_snps` | input mode | packaged HM3 regression SNP restriction |
| `LDScoreConfig` | `snp_batch_size` | performance | LD-score SNP batch size |
| `LDScoreConfig` | `common_maf_min` | input metadata | common-SNP count threshold only |
| `LDScoreOutputConfig` | `output_dir` | output | canonical LD-score result directory |
| `run_ldscore(**kwargs)` | `baseline_annot_sources`, `query_annot_sources`, `query_annot_bed_sources` | input | optional annotation sources; query inputs require baseline sources, and no-annotation runs synthesize `base` |
| `run_ldscore(**kwargs)` | `plink_prefix`, `r2_dir` | input | reference-panel sources |
| `run_ldscore(**kwargs)` | `output_dir` | output | canonical result directory; convenience wrapper writes `diagnostics/ldscore.log` |

Removed Python names: `bfile`, `r2_table`, `frqfile`, `keep`, `maf`,
`baseline_annot`, `query_annot`, `query_annot_bed`, `out`,
`r2_ref_panel_dir`, `ref_panel_dir`, `r2_sources`, `metadata_sources`, and
LD-score `chunk_size`.

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
| `ReferencePanelBuildConfig` | `use_hm3_snps` | input mode | packaged HM3 retained SNP restriction |
| `ReferencePanelBuildConfig` | `use_hm3_quick_liftover` | input mode | packaged HM3-only coordinate liftover |
| `ReferencePanelBuildConfig` | `keep_indivs_file` | input | PLINK individual keep file |
| `ReferencePanelBuildConfig` | `maf_min` | input metadata | retained SNP MAF filter |
| `ReferencePanelBuildConfig` | `snp_batch_size` | performance | SNP computation batch size |
| `ReferencePanelBuildConfig` | `output_dir` | output | artifact directory |
| `run_build_ref_panel(**kwargs)` | same config field names except global settings | input/output | convenience wrapper; reads `snp_identifier` from the registered `GlobalConfig`; ignores `GlobalConfig.genome_build`; writes the build-ref-panel workflow log |

Removed Python names: `plink_path`, `bfile`, `out`, `panel_label`,
`keep_indivs`, `maf`, old genetic-map and liftover names without `_file` /
`_sources`.

### Sumstats munging

| Object/function | Argument | Direction | Object |
|---|---:|---:|---|
| `MungeConfig` | `raw_sumstats_file` | input | raw summary-statistics file |
| `MungeConfig` | `trait_name` | input metadata | optional trait label |
| `MungeConfig` | `column_hints` | input metadata | optional source-column hints |
| `MungeConfig` | `daner_old`, `daner_new` | input metadata | optional DANER schema interpretation switches |
| `MungeConfig` | `sumstats_snps_file` | input | summary-statistics SNP keep-list |
| `MungeConfig` | `use_hm3_snps` | input mode | packaged HM3 summary-statistics SNP restriction |
| `MungeConfig` | `target_genome_build` | input metadata | optional target build for `chr_pos`-family output coordinates |
| `MungeConfig` | `liftover_chain_file` | input | optional source-to-target chain file for munger liftover |
| `MungeConfig` | `use_hm3_quick_liftover` | input mode | use packaged curated HM3 dual-build map for coordinate-only liftover |
| `MungeConfig` | `output_dir` | output | munged output directory; `SumstatsMunger.run()` writes `diagnostics/sumstats.log` |
| `MungeConfig` | `output_format` | output mode | `parquet`, `tsv.gz`, or `both`; defaults to `parquet` |
| `SumstatsMunger.run(munge_config, ...)` | `munge_config` | input/output | normalized munging workflow; owns fixed output preflight, root `metadata.json`, diagnostics, always-written `diagnostics/dropped_snps/dropped.tsv.gz`, and result construction; summary `output_paths` excludes logs |
| `SumstatsMunger.write_output(sumstats, output_dir, output_format='parquet')` | `output_dir` | output | writes fixed `sumstats.parquet` and/or `sumstats.sumstats.gz` |

Removed Python names: legacy separate source-path object field,
`MungeConfig.sumstats_file`, `MungeConfig.out_prefix`,
`write_output(..., out_prefix)`.

### Regression

| Object/function | Argument | Direction | Object |
|---|---:|---:|---|
| `run_h2_from_args(args)` | `ldscore_dir` | input | LD-score result directory |
| `run_h2_from_args(args)` | `sumstats_file` | input | munged summary-statistics file |
| `run_h2_from_args(args)` | `output_dir` | output | writes `h2.tsv` plus diagnostics when supplied; otherwise prints compact TSV to stdout |
| `run_partitioned_h2_from_args(args)` | `ldscore_dir` | input | LD-score result directory |
| `run_partitioned_h2_from_args(args)` | `sumstats_file` | input | munged summary-statistics file |
| `run_partitioned_h2_from_args(args)` | `output_dir` | output | writes `partitioned_h2.tsv` plus diagnostics when supplied; otherwise prints compact TSV to stdout |
| `run_partitioned_h2_from_args(args)` | `write_per_query_results` | output mode | optionally writes `diagnostics/query_annotations/` |
| `run_rg_from_args(args)` | `ldscore_dir` | input | LD-score result directory |
| `run_rg_from_args(args)` | `sumstats_sources` | input | two or more munged summary-statistics files or glob patterns |
| `run_rg_from_args(args)` | `anchor_trait` | input selector | optional anchor trait label or path for anchor-vs-rest output |
| `run_rg_from_args(args)` | `output_dir` | output | writes `rg.tsv`, `rg_full.tsv`, `h2_per_trait.tsv`, and diagnostics when supplied; otherwise prints compact `rg.tsv` to stdout |
| `run_rg_from_args(args)` | `write_per_pair_detail` | output mode | optionally writes `diagnostics/pairs/manifest.tsv` and per-pair diagnostic folders when `output_dir` is supplied |

Removed Python/public argparse names: `sumstats`, `sumstats_1`, `sumstats_2`,
`out`, `ldscore`, `counts`, `w_ld`, `annotation_manifest`, `query_columns`.

Current curated `sumstats.parquet` and `.sumstats.gz` artifacts provide
canonical `SNP`, `CHR`, `POS`, `Z`, and `N` fields when written by
`ldsc munge-sumstats`, plus neighboring root `metadata.json` with
`schema_version`, `artifact_type`, `snp_identifier`, `genome_build`, and
optional `trait_name`. `load_sumstats()` reconstructs the config snapshot from
that identity provenance and resolves labels as explicit override, then sidecar
`trait_name`, then filename fallback.
Liftover reports, coordinate provenance, selected curated output files, and
Parquet row groups are log provenance, not metadata sidecar payload.
Regression therefore merges on the effective identity key: literal `SNP` in
`rsid`, `SNP:<allele_set>` in `rsid_allele_aware`, `CHR:POS` in `chr_pos`, and
`CHR:POS:<allele_set>` in `chr_pos_allele_aware`. Munger liftover is valid only in
`chr_pos`-family modes; it changes `CHR`/`POS`, never `SNP`, and runs after
`sumstats_snps_file` filtering. Liftover drop counts are written as readable log
records, examples appear only at `DEBUG`, and row-level dropped-SNP details are
written to `diagnostics/dropped_snps/dropped.tsv.gz`. The metadata sidecar stays limited
to current artifact provenance and must include the minimal identity fields
`schema_version`, `artifact_type`, `snp_identifier`, and `genome_build`.

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
  `PartitionedH2OutputConfig`, `PartitionedH2DirectoryWriter`,
  `RgOutputConfig`, and `RgDirectoryWriter` for output writing.

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
