# IO Argument Inventory

Last updated on: 2026-08-05

This document records the current public input/output naming contract after the
LD-score result-directory refactor. The LD-score workflow uses a canonical
result directory as the baseline design:

```text
<ldscore_dir>/
  metadata.json
  ldscore.baseline.parquet
  ldscore.query.parquet        # omitted when no query annotations exist
  ldscore.overlap.parquet      # overlap matrix; omitted for single-annotation (e.g. base-only) runs; required by partitioned-h2
  diagnostics/
    ldscore.log
    query_annotation_status.tsv  # BED/gene-list query runs
    gene_list_unresolved.tsv.gz  # gene-list runs
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
Package-built R2 parquets use the 4-column index format and serve all four
identifier modes. External R2 parquet formats are not supported.

## Format Policy

The package separates machine-consumed internal artifacts from
science-facing result tables. 
Internal artifacts are machine-consumed intermediate products and should use Parquet as their primary format. This applies to package-built R2 pair tables,
canonical LD-score result tables, and curated munged summary statistics. JSON
metadata is a downstream contract for LD-score artifacts; curated sumstats
store downstream identity metadata in the `sumstats.parquet` footer. Non-consumed
workflow metadata belongs under `diagnostics/`.

Science-facing result tables should use TSV. The public regression results are
`h2.tsv`, `partitioned_h2.tsv`, optional diagnostic per-query
`partitioned_h2.tsv` / `partitioned_h2_full.tsv`, and the rg family
(`rg.tsv`, `rg_full.tsv`, `h2_per_trait.tsv`, optional diagnostic
`diagnostics/pairs/`).
Logs are audit files and should not be treated as output results.


### Current Adaptation Status

Adapted public paths:

- `ldsc ldscore` writes `ldscore.baseline.parquet`, optional
  `ldscore.query.parquet`, and optional `ldscore.overlap.parquet`.
- `ldsc munge-sumstats` writes `sumstats.parquet` by default.
- `ldsc h2`, `ldsc partitioned-h2`, and `ldsc rg` write result tables as TSV
  when an `output_dir` is supplied.
- `ldsc build-ref-panel` writes the primary R2 pair matrix as
  `chr{chrom}_r2.parquet`.
- `ldsc build-gene-ldscore-index` writes one immutable gene LD-score index
  directory for later explicit indexed LD-score assembly.
- `ldsc convert-ldsc2-ldscores` is the only public reader of selected legacy
  LDSC2 LD-score fragments and writes a canonical LDSC3 LD-score directory.
- `ldsc query-r2` reads package-built index-format R2 panels and writes an
  annotated pair table to stdout or one explicit TSV path.

Not fully adapted or retained for compatibility:

- `ldsc annotate` writes generated query annotations as
  `query.<chrom>.annot.gz`; these are intermediate annotation artifacts rather
  than science-facing result tables.
- `ldsc build-ref-panel` still writes runtime metadata sidecars as
  `chr{chrom}_meta.tsv.gz` and liftover-stage drop audit files as
  `diagnostics/dropped_snps/chr{chrom}_dropped.tsv.gz`.
- `ldsc munge-sumstats` keeps `--output-format tsv.gz` and `both`, which write
  `sumstats.sumstats.gz` compatibility artifacts even though parquet is the
  default.

Private `_kernel` modules do not own LDSC2 artifact emission or regression
input compatibility. The workflow-layer sumstats writer retains the optional
metadata-free gzip output; legacy LD-score suites cross only the explicit
converter boundary.

## Naming Rules

Public CLI flags and Python config fields follow these rules:

| Suffix | Meaning | Examples |
|---|---|---|
| `*_file` | one file-like input, exact-one glob allowed where the resolver supports it | `raw_sumstats_file`, `sumstats_file`, `sumstats_snps_file`, `keep_indivs_file` |
| `*_sources` | one logical input that may resolve to many files via globs, comma lists, or supported `@` chromosome tokens | `baseline_annot_sources`, `query_annot_bed_sources`, `query_annot_gene_list_sources` |
| `*_dir` | directory input or output location | `ldscore_dir`, `output_dir` |

Removed from artifact-writing workflow surfaces:

- generic output/run identity names: `prefix`, `out`, `out_prefix`,
  `panel_label`
- `bfile`, `frqfile`
- legacy aliases for annotation BEDs, LD-score components, and regression
  outputs

Users customize run identity by choosing the `output_dir` name. Output filenames
inside that directory are fixed and workflow-specific. `ldsc query-r2` follows
the same `output_dir` model (`query_r2.tsv` plus the `diagnostics/` sidecar); it
additionally streams the result table to stdout when `--output-dir` is omitted,
for pipe-able interactive use.

Workflow logs are fixed audit files under `output_dir`, preflighted with the
scientific outputs before they are opened. They are not included in workflow
`output_paths` mappings or thin metadata sidecars that downstream code
interprets as data artifacts.

Fixed workflow outputs are coherent artifact families. Without overwrite, any
existing current-contract owned artifact rejects the run, even when that
artifact would not be written by the current configuration. With overwrite, the
workflow writes the requested outputs and removes stale current-contract owned
siblings not produced by the successful run. Legacy/root diagnostic names that
are no longer in the public layout are not blocked or cleaned as owned outputs.
Unrelated files in `output_dir` are preserved. Sharded workflows may narrow the
owned family to the shard selected by the current invocation. For
`build-ref-panel`, a concrete chromosome prefix owns only that chromosome's
artifact package, while a `@` chromosome-suite invocation owns the full panel
package.

Directory artifacts such as `diagnostics/query_annotations/` and
`diagnostics/pairs/` are owned as whole trees. They are staged under
`diagnostics/` and moved into final position as a unit when requested; if a
later overwrite run omits them, the stale tree is removed after the new outputs
are written.

## CLI Inventory

### `ldsc annotate`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--query-annot-bed-sources` | input | yes | BED interval files | Accepts exact files, globs, comma-separated tokens, and source-token lists. Each resolved BED file stem (`Path.stem`) becomes a query annotation name; duplicate stems are rejected. |
| `--baseline-annot-sources` | input | yes | baseline `.annot[.gz]` templates | Accepts exact files, globs, and `@` chromosome-suite tokens. |
| `--output-dir` | output | yes | generated query annotation directory | Writes combined root `query.<chrom>.annot.gz` files, with all BED inputs represented as query columns, plus diagnostic `metadata.json`, `dropped_snps/dropped.tsv.gz`, and `annotate.log` under `diagnostics/`. |
| `--padding-bp` | input transform | no | BED interval expansion | Adds this many base pairs to both sides of each BED interval before SNP projection; starts are clipped at zero. Defaults to `0`, so BED intervals are used as provided. |
| `--overwrite` | output mode | no | collision policy | Controls whether generated annotation files and diagnostics may be replaced; defaults to `False`, so any existing root-level `query.*.annot.gz` shard or owned diagnostic artifact is refused. With overwrite, stale query shards outside the current chromosome set are removed after a successful run. |
| `--log-level` | logging | no | workflow log verbosity | Controls ordinary LDSC logger record verbosity; defaults to `INFO`; these records go to `diagnostics/annotate.log` and the CLI console (stderr) shows only errors. Lifecycle audit lines always appear in the file. |
| `--snp-identifier` | config | no | SNP identity mode | Defines how SNPs are keyed. Defaults to `chr_pos_allele_aware`; valid values are `rsid`, `rsid_allele_aware`, `chr_pos`, and `chr_pos_allele_aware`. |
| `--genome-build` | config | no | coordinate interpretation | Coordinate build for `chr_pos`-family inputs. Defaults to omitted/`None`, which is invalid for coordinate-family inputs unless the workflow can infer it; accepted values include `auto`, `hg19`/`GRCh37`, and `hg38`/`GRCh38`. |

Removed flags: `--bed-files`, `--baseline-annot`, `--bed-padding-bp`.

### `ldsc ldscore`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--output-dir` | output | yes | canonical LD-score result directory | Writes root scientific artifacts and `diagnostics/ldscore.log`; BED/gene runs add `query_annotation_status.tsv`, and gene runs add `gene_list_unresolved.tsv.gz`. Parquet row groups are chromosome-aligned. |
| `--overwrite` | output mode | no | collision policy | Controls whether fixed LD-score files and `diagnostics/ldscore.log` may be replaced; defaults to `False`, so any existing owned LD-score artifact in `output_dir` is refused. With overwrite, stale `ldscore.query.parquet` is removed after successful baseline-only runs. |
| `--log-level` | logging | no | workflow log verbosity | Controls ordinary LDSC logger record verbosity; defaults to `INFO`; these records go to `diagnostics/ldscore.log` and the CLI console (stderr) shows only errors. Lifecycle audit lines always appear in the file. |
| `--baseline-annot-sources` | input | no | baseline annotation files | Supplies baseline annotation files; defaults to omitted/`None`, and if no query inputs are supplied `ldscore` synthesizes an all-ones `base` column. |
| `--query-annot-sources` | input | no | prebuilt query annotation files | Supplies prebuilt query annotations; defaults to omitted/`None`. Mutually exclusive with the BED and gene-list routes and requires `--baseline-annot-sources`. |
| `--query-annot-bed-sources` | input | no | query BED interval files | Supplies BED intervals projected in memory; defaults to omitted/`None`. Mutually exclusive with the prebuilt and gene-list routes and requires `--baseline-annot-sources`. Concrete-source failures are recorded and skipped while usable siblings continue. |
| `--query-annot-gene-list-sources` | input | no | one-column gene lists | Supplies exact/glob file groups resolved against the packaged protein-coding catalog; defaults to omitted/`None`. Explicit files retain order and globs expand lexically. Mutually exclusive with other query routes, does not use `@`, and requires `--baseline-annot-sources`. |
| `--padding-bp` | input transform | conditional | live BED or gene interval expansion | Adds this many base pairs to both sides of intervals derived from `--query-annot-bed-sources` or `--query-annot-gene-list-sources`; starts are clipped at zero. The effective default is `0`. Supplying the flag with prebuilt annotation queries, no query, or `--gene-ldscore-index-dir` is rejected, including an explicitly supplied value of `0`. |
| `--plink-prefix` | input | conditional | PLINK reference panel prefix | Selects PLINK reference-panel input; defaults to omitted/`None` and is required when `--r2-dir` is omitted. Supports exact prefix, PLINK-prefix glob, or `@` suite. |
| `--r2-dir` | input | conditional | package-built parquet R2 directory | Selects parquet reference-panel input; defaults to omitted/`None` and is required when `--plink-prefix` is omitted. Use a build-specific directory such as `ref_panel/hg38`. The directory must contain paired `chrN_r2.parquet` (4-column index format) and `chrN_meta.tsv.gz` sidecar files; the sidecar is mandatory. One parquet serves all identifier modes. |
| `--snp-identifier` | config | no | SNP identity mode | Defines how SNPs are keyed. Defaults to `chr_pos_allele_aware`; valid values are `rsid`, `rsid_allele_aware`, `chr_pos`, and `chr_pos_allele_aware`. |
| `--genome-build` | config | no | coordinate and gene-interval interpretation | Build for `chr_pos` identity and gene-list interval projection; defaults to omitted/`None`, while gene-list runs resolve omission to `auto`, including rsID modes. rsID artifact identity metadata remains build-independent. Accepted aliases normalize to hg19 or hg38. |
| `--ref-panel-snps-file` | input | no | reference-panel SNP universe restriction | Restricts the retained reference-panel SNP universe using identity keys only; duplicate restriction keys collapse to one retained key, and non-identity columns such as `CM` or `MAF` are ignored. Defaults to omitted/`None`, so no additional restriction is applied. |
| `--regression-snps-file` | input | no | regression row-set override | Replaces the bundled HM3 regression SNP set using identity keys only; defaults to omitted/`None`, so the bundled HM3 set is used. Duplicate restriction keys collapse to one retained key, and non-identity columns such as `CM` or `MAF` are ignored. |
| `--keep-indivs-file` | input | no | PLINK individual keep file | Restricts PLINK individuals before LD calculation; defaults to omitted/`None`, so no individual keep filter is applied. PLINK mode only. |
| `--maf-min` | input metadata | no | retained reference-panel MAF filter | Filters retained reference-panel SNPs by inclusive `MAF >= maf_min`; defaults to omitted/`None`. Applied identically in **both** backends (parquet via the sidecar MAF; PLINK via genotype-derived MAF). |
| `--exclude-regions` | input transform | no | named regression-region presets | After bundled HM3 or the custom regression list is selected, subtracts `none`, `mhc`, `centromeres`, or `mhc-and-centromeres` (`mhc-and-centromeres` default) from regression/output rows and `w_ld` contributors. It does not change LD-score contributors or reference count vectors. `--genome-build` selects named interval coordinates. |
| `--ld-wind-snps` | model | conditional | LD window in SNPs | Selects an LD window measured by SNP count; defaults to `None`. Exactly one of `--ld-wind-snps`, `--ld-wind-kb`, or `--ld-wind-cm` must be supplied. |
| `--ld-wind-kb` | model | conditional | LD window in kilobases | Selects an LD window measured by physical distance; defaults to `None`. Exactly one LD-window flag must be supplied. |
| `--ld-wind-cm` | model | conditional | LD window in centiMorgans | Selects an LD window measured by genetic distance; defaults to `None`. Exactly one LD-window flag must be supplied. Requires usable reference-panel `CM` (≥2 distinct finite values per chromosome); an all-zero/constant/missing `CM` raises a dedicated error that `--yes-really` does **not** bypass. PLINK panels with uninformative `.bim` `CM` can supply a genetic map (see below). |
| `--genetic-map-hg19-sources` | input | no | hg19 genetic map for PLINK cM windows | Genetic map (hg19) used to derive `CM` at `.bim` positions when the PLINK `.bim` `CM` column is uninformative; defaults to omitted/`None`. An explicit map always wins (used for all chromosomes). Ignored with a warning for the parquet backend (sidecar `CM` is authoritative). |
| `--genetic-map-hg38-sources` | input | no | hg38 genetic map for PLINK cM windows | Genetic map (hg38) used to derive `CM` at `.bim` positions when the PLINK `.bim` `CM` column is uninformative; defaults to omitted/`None`. The map build is selected from `--genome-build` or inferred (chr_pos modes); rsID modes require an explicit `--genome-build`. |
| `--export-ref-metadata` | output | no | opt-in reference-metadata sidecar | When set (PLINK backend), writes `ref_metadata/chrN_meta.tsv.gz` (`CHR POS SNP A1 A2 CM MAF`, matching the parquet panel sidecar) next to the LD-score output. Default `False`. Parquet panels already ship this sidecar. |
| `--common-maf-min` | input metadata | no | common-SNP count threshold | Sets the MAF threshold for common-SNP count vectors and the common-universe overlap matrix; defaults to `0.05` and uses inclusive `MAF >= common_maf_min` (deviates from legacy LDSC's strict `0.05 < FRQ < 0.95`). |
| `--snp-batch-size` | performance | no | LD-score SNP batch size | Genotype batch size for the PLINK reference-panel backend; defaults to `128`. The parquet-R2 backend streams stored pairs and ignores this value. |
| `--threads` | performance | no | cross-chromosome parallelism | Worker processes for cross-chromosome parallelism (joblib `n_jobs` convention): defaults to `1` (sequential), `N`=`N` workers, `-1`=all cores, `-2`=all but one. Capped at the chromosome count and CPU affinity. Workers are processes (the LD-score kernel is CPU-bound), each rebuilding its chromosome's panel, so peak RSS scales with the worker count. |
| `--yes-really` | safety override | no | whole-chromosome LD windows | Allows whole-chromosome LD windows; defaults to `False`, so such windows are rejected unless this flag is supplied. |

Removed flags: `--bfile`, `--r2-table`, `--frqfile`, `--r2-sources`,
`--metadata-sources`, `--keep`, `--maf`, `--baseline-annot`,
`--query-annot`, `--query-annot-bed`, `--chunk-size`, `--r2-bias-mode`,
`--r2-sample-size`, legacy `--out`. R2 bias mode and sample size are now read
solely from parquet schema metadata (`ldsc:r2_bias` / `ldsc:n_samples`).

LD-score output schema:

- `ldscore.baseline.parquet`: `CHR`, `POS`, `SNP`, `regression_ld_scores`, then baseline
  LD-score columns. In no-annotation unpartitioned runs, the baseline column
  list is exactly `base`.
- `ldscore.query.parquet`: `CHR`, `POS`, `SNP`, then query LD-score columns; omitted
  when there are no query annotations.
- `ldscore.overlap.parquet`: `row_annotation`, `col_annotation`,
  `overlap_all_snps`, `overlap_common_snps` (long-form annotation overlap matrix
  consumed by partitioned-h2); omitted for single-annotation (e.g. base-only) runs.
- `metadata.json`: `artifact_type`, relative file paths (`baseline`, optional
  `query`, optional `overlap`), baseline/query column names, count records, `count_config`,
  `overlap_config`, config metadata, chromosomes, row
  counts, `row_group_layout`, `baseline_row_groups`, and
  `query_row_groups`.

### `ldsc build-ref-panel`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--plink-prefix` | input | yes | PLINK reference panel prefix | Supports exact prefix, PLINK-prefix glob, or `@` suite. |
| `--source-genome-build` | input metadata | no | source PLINK coordinate build | Declares the PLINK coordinate build; defaults to `auto`, so the build is inferred from `.bim` before SNP restriction. |
| `--genetic-map-hg19-sources` | input | conditional | hg19 genetic map file or suite | Supplies hg19 genetic-map CM values; defaults to omitted/`None` and is required when `--ld-wind-cm` is used and hg19 output is emitted. |
| `--genetic-map-hg38-sources` | input | conditional | hg38 genetic map file or suite | Supplies hg38 genetic-map CM values; defaults to omitted/`None` and is required when `--ld-wind-cm` is used and hg38 output is emitted. |
| `--liftover-chain-hg19-to-hg38-file` | input | no | liftover chain file | Enables hg38 outputs for hg19 source builds; defaults to omitted/`None`, so those target-build outputs are not emitted. |
| `--liftover-chain-hg38-to-hg19-file` | input | no | liftover chain file | Enables hg19 outputs for hg38 source builds; defaults to omitted/`None`, so those target-build outputs are not emitted. |
| `--ref-panel-snps-file` | input | no | retained SNP universe restriction | Restricts the emitted reference-panel SNP universe using identity keys only; duplicate restriction keys collapse to one retained key, and non-identity columns such as `CM` or `MAF` are ignored. Defaults to omitted/`None`, so no retained SNP restriction is applied. |
| `--snp-identifier` | input metadata | conditional | global SNP identifier mode | Overrides the registered SNP identifier mode for this invocation; defaults to `chr_pos_allele_aware` (the parser and `GlobalConfig` default). |
| `--keep-indivs-file` | input | no | PLINK individual keep file | Restricts PLINK individuals during panel building; defaults to omitted/`None`, so no individual keep filter is applied. |
| `--maf-min` | input metadata | no | retained SNP MAF filter | Filters retained SNPs by MAF during PLINK loading; defaults to omitted/`None`, so no retained-SNP MAF filter is applied. |
| `--ld-wind-snps` | model | conditional | LD window in SNPs | Selects an LD window measured by SNP count; defaults to `None`. Exactly one of `--ld-wind-snps`, `--ld-wind-kb`, or `--ld-wind-cm` must be supplied. |
| `--ld-wind-kb` | model | conditional | LD window in kilobases | Selects an LD window measured by physical distance; defaults to `None`. Exactly one LD-window flag must be supplied. |
| `--ld-wind-cm` | model | conditional | LD window in centiMorgans | Selects an LD window measured by genetic distance; defaults to `None`. Exactly one LD-window flag must be supplied, and cM windows require genetic-map inputs for every emitted build. |
| `--output-dir` | output | yes | reference-panel artifact directory | Run identity is `Path(output_dir).name`; no separate label is accepted. |
| `--overwrite` | output mode | no | collision policy | Controls whether reference-panel artifacts, diagnostic metadata, always-written `diagnostics/dropped_snps/chr{chrom}_dropped.tsv.gz` audit files, and build-ref-panel workflow logs may be replaced; defaults to `False`, so existing owned outputs are refused. Concrete chromosome prefixes own and clean only that chromosome's package. `@` chromosome-suite runs own the full panel package and can remove stale target-build, out-of-scope chromosome, dropped-SNP, metadata, or log siblings after success. |
| `--log-level` | logging | no | workflow log verbosity | Controls ordinary LDSC logger record verbosity; defaults to `INFO`; these records go to the build-ref-panel workflow log and the CLI console (stderr) shows only errors. Lifecycle audit lines always appear in the file. |
| `--snp-batch-size` | performance | no | SNP computation batch size | Number of SNPs decoded per pairwise-R2 computation batch; larger values size the pairwise working set (decoded window and correlation block) and can improve throughput. It does **not** drive peak RSS: the genotype payload is read selectively (restricted builds) or streamed (unrestricted builds), so peak is governed by that bounded read plus the workflow/import floor, not by this batch size. Defaults to `128`. |
| `--min-r2` | output mode | no | pair-emission threshold | Optional unbiased-R2 floor for emitted pairs. Defaults to `0.0`, which writes every retained pair. Positive values reduce output size by omitting low-R2 pairs; downstream query and LD-score reads treat absent pairs as zero/absent according to their workflow contract. The threshold is recorded in parquet metadata as `ldsc:min_r2`. |

Removed flags: `--bfile`, `--out`, `--panel-label`, `--keep-indivs`, `--maf`,
`--genetic-map-hg19`, `--genetic-map-hg38`, old liftover-chain names without
`_file` / `_sources`, `--duplicate-position-policy`, `--chunk-size` (hidden alias;
use `--snp-batch-size`).

Fixed output names:

```text
<output_dir>/hg19/chr{chrom}_r2.parquet
<output_dir>/hg19/chr{chrom}_meta.tsv.gz
<output_dir>/hg38/chr{chrom}_r2.parquet
<output_dir>/hg38/chr{chrom}_meta.tsv.gz
<output_dir>/diagnostics/metadata.json
<output_dir>/diagnostics/metadata.chr<chrom>.json  # concrete single-chromosome PLINK prefix
<output_dir>/diagnostics/dropped_snps/chr{chrom}_dropped.tsv.gz
<output_dir>/diagnostics/build-ref-panel.log
<output_dir>/diagnostics/build-ref-panel.chr<chrom>.log  # concrete single-chromosome PLINK prefix
```

Each `chr{chrom}_r2.parquet` stores Arrow schema metadata for
`ldsc:sorted_by_build`, `ldsc:row_group_size`, `ldsc:n_samples`, and
`ldsc:r2_bias`. Current package-built panels write unbiased R2 values and record
the PLINK sample count. Downstream LD-score runs read R2 bias mode and sample
size from this metadata; there are no bias-related flags. External raw-R2 panels
must declare `ldsc:r2_bias=raw` and `ldsc:n_samples` in their parquet schema to
be corrected.

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

`build-ref-panel` owns only current-contract artifacts. The owned package depends
on the PLINK prefix scope: a concrete single-chromosome prefix owns only that
chromosome's R2, metadata sidecar, dropped-SNP audit, chromosome-scoped
diagnostic metadata, and chromosome-scoped log; a `@` chromosome-suite prefix
owns the full all-chromosome panel package. Existing owned artifacts block
without `--overwrite`; with `--overwrite`, stale current-contract siblings
inside the owned package are removed after the successful write.

### `ldsc build-gene-ldscore-index`

This advanced, deliberately closed command builds an immutable hg19 gene
LD-score index under an explicitly selected base `rsid` or `chr_pos` identity.
It is not a legacy compatibility command and performs no inference or liftover.

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--baseline-annot-sources` | input | yes | full baseline annotation suite | Canonical baseline used by later indexed LD-score assembly. |
| `--plink-prefix` | input | yes | PLINK reference panel | Exact prefix or supported chromosome suite. |
| `--output-dir` | output | yes | immutable index directory | Publishes one complete index artifact; no incremental update or component reuse. |
| `--genome-build` | config | yes | coordinate build assertion | Closed to explicit `hg19` in both modes; no default or inference. |
| `--snp-identifier` | config | yes | SNP identity | Exactly `rsid` or `chr_pos`; no default, `auto`, or allele-aware mode. |
| `--padding-bp` | model | no | gene interval padding | Padding used when constructing disjoint gene atoms; defaults to `0`. |
| `--gene-exclude-regions` | model | no | gene exclusion policy | `none` or `mhc`; defaults to `mhc`. |
| `--ld-wind-cm` | model | no | LD window | cM window used by the PLINK-backed operator; defaults to `1.0`. |
| `--maf-min`, `--common-maf-min` | QC/counts | no | reference MAF thresholds | Retained-reference and common-count thresholds stored in index identity; `--maf-min` defaults to `None`, and `--common-maf-min` defaults to `0.05`. |
| `--keep-indivs-file` | input | no | PLINK sample restriction | Restricts individuals used during index construction; defaults to omitted/`None`. |
| `--regression-snps-file` | input | no | regression SNP override | Identity-only list; defaults to omitted/`None`, so the bundled hg19 HM3 set is used. Repeated effective keys collapse. |
| `--exclude-regions` | input transform | no | regression-region subtraction | `none`, `mhc`, `centromeres`, or `mhc-and-centromeres`; defaults to `mhc-and-centromeres`. |
| `--genetic-map-hg19-sources` | input | no | genetic map | Optional explicit hg19 map; defaults to omitted/`None`; hg38 maps are not a builder CLI input. |
| `--chromosomes` | input selector | no | chromosome set | Explicit chromosome selection; defaults to `1-22`. |
| `--snp-batch-size`, `--atom-batch-size` | performance | no | computation batches | Bounds SNP and disjoint-atom matrix work; `--snp-batch-size` defaults to `128`, and `--atom-batch-size` defaults to `64`. |
| `--threads` | performance | no | chromosome workers | Cross-chromosome process count; defaults to `1`. |
| `--overwrite` | output mode | no | publication policy | Replaces only a complete valid owned index through staged publication; defaults to `False`. |
| `--log-level` | logging | no | workflow log verbosity | Controls the persistent index build-state log; defaults to `INFO`. |

### `ldsc convert-ldsc2-ldscores`

This rare, explicit migration utility accepts only reusable unpartitioned
reference/weight suites and complete baseline partitioned suites. It does not
accept query annotations, thin annotations, or arbitrary partitioned suites.

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--legacy-reference-dir` | input | yes | LDSC2 reference suite directory | Exactly one coherent chromosomes 1-22 `.l2.ldscore[.gz]` family. A complete full `.annot[.gz]` family selects baseline conversion; absence selects unpartitioned conversion only when the score table has one scientific column. |
| `--legacy-weight-dir` | input | yes | LDSC2 regression-weight suite directory | Exactly one coherent chromosomes 1-22 weight LD-score family. Rows are inner-joined to reference rows by rsID. |
| `--legacy-frequency-dir` | input | baseline only | LDSC2 `.frq[.gz]` suite directory | Defaults to omitted/`None`; required for baseline conversion to reconstruct strict legacy common counts and overlap, and omitted for unpartitioned conversion. |
| `--output-dir` | output | yes | canonical LDSC3 LD-score directory | Writes `metadata.json`, `ldscore.baseline.parquet`, optional baseline overlap, and conversion diagnostics without modifying source files. |
| `--snp-identifier` | config | no | output identity | `rsid` (default) or `chr_pos`; allele-aware modes are forbidden because the source suite has no authoritative alleles. All source-family joins still use rsID. |
| `--genome-build` | config | no | reference coordinate interpretation | `auto` (default), `hg19`, or `hg38`. Required inference applies only to `chr_pos`; rsID output remains buildless and records best-effort evidence as provenance. |
| `--overwrite` | output mode | no | collision policy | Replaces converter-owned outputs only; defaults to `False`; legacy inputs remain read-only. |
| `--log-level` | logging | no | workflow log verbosity | Defaults to `INFO`; writes `diagnostics/convert-ldsc2-ldscores.log`. |

There is intentionally no profile flag and no common-MAF threshold flag.
Profile selection is structural. `.l2.M_5_50` is required and always means the
strict legacy rule `MAF > 0.05` (or `0.05 < FRQ < 0.95`). `.l2.M` is optional:
it remains unavailable for an unpartitioned suite when any chromosome is
missing, while baseline conversion reconstructs missing values from full
annotations and frequencies. See
[`legacy-ldscore-conversion.md`](legacy-ldscore-conversion.md).

### `ldsc munge-sumstats`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--raw-sumstats-file` | input | yes | raw summary-statistics file | Exact path or exact-one glob. |
| `--format` | input metadata | no | raw summary-statistics format profile | One of `auto`, `plain`, `daner-old`, or `daner-new`; defaults to `auto`, which detects common plain text, including VCF-style headers and old DANER. This is the sole DANER selector (the legacy `--daner-old`/`--daner-new` booleans are removed). |
| `--infer-only` | diagnostic | no | raw summary-statistics inference report | Reads the raw header and first data row, prints detected format, inferred hints, missing fields, source/output genome-build status, liftover status, notes, and suggested commands. Defaults to `False`; missing `A1/A2` is reported only in allele-aware modes. Does not require `--output-dir` and writes no artifacts. |
| `--sumstats-snps-file` | input | no | summary-statistics SNP keep-list | Restricts munged summary-statistics rows using identity keys only; duplicate restriction keys collapse to one retained key, and non-identity columns such as `CM` or `MAF` are ignored. The keep-list is loaded before parsing and applied while chunks are streaming; defaults to omitted/`None`, so no keep-list restriction is applied. |
| `--use-hm3-snps` | input mode | no | packaged HM3 SNP restriction | Restricts munged summary-statistics rows to the packaged curated HM3 map while chunks are streaming; defaults to `False`. Mutually exclusive with `--sumstats-snps-file`. |
| `--trait-name` | input metadata | no | biological trait label | Optional label stored in the `sumstats.parquet` footer when parquet output is written; defaults to omitted/`None`. Downstream regression uses it unless a regression CLI `--trait-name` override is supplied. No root `metadata.json` sidecar is written. |
| `--output-dir` | output | yes, except `--infer-only` | munged output directory | Defaults to omitted/`None`; required for normal artifact-writing runs but not with `--infer-only`. The workflow writes fixed `sumstats.*` artifacts under this directory and passes `<output_dir>/sumstats` as the kernel output stem. |
| `--output-format` | output mode | no | curated sumstats format | One of `parquet`, `tsv.gz`, or `both`; defaults to `parquet`. |
| `--N`, `--N-cas`, `--N-con` | model/QC | no | sample-size overrides | Scalar total, case, and control sample-size overrides forwarded to the munging kernel; each defaults to `None`. |
| `--info-min` | QC | no | INFO threshold | Minimum INFO value retained by the munging kernel; defaults to `0.9`. |
| `--maf-min` | QC | no | MAF threshold | Minimum allele-frequency value retained by the munging kernel; defaults to `0.01`. |
| `--n-min` | QC | no | minimum sample size | Optional minimum sample-size row filter; defaults to `None`. |
| `--chunksize` | performance | no | raw chunk size | Number of raw rows streamed per munging chunk; defaults to `1000000`. |
| `--snp` | input metadata | no | raw SNP column hint | Identifies the raw SNP/variant identifier column; defaults to omitted/`None`, so common aliases are inferred. |
| `--chr` | input metadata | no | raw chromosome column hint | Identifies the raw chromosome column; defaults to omitted/`None`, so common aliases such as `#CHROM`, `CHROM`, and `CHR` are inferred. |
| `--pos` | input metadata | no | raw position column hint | Identifies the raw position column; defaults to omitted/`None`, so common aliases such as `POS` and `BP` are inferred. |
| `--N-col`, `--N-cas-col`, `--N-con-col` | input metadata | no | sample-size column hints | `--N-col` selects direct N and suppresses inferred case/control columns. `--N-cas-col` and `--N-con-col` must be supplied together, select case/control N, and suppress inferred direct N. The two strategies are mutually exclusive; fully automatic discovery of both is an actionable ambiguity error. |
| `--a1`, `--a2` | input metadata | no | allele column hints | Identify raw effect and other allele columns; each defaults to `None`. Allele-aware modes require usable alleles after parsing. |
| `--p` | input metadata | no | p-value column hint | Identifies the raw p-value column; defaults to `None`. |
| `--frq` | input metadata | no | allele-frequency column hint | Identifies the raw allele-frequency column used by munger QC and optional output preservation; defaults to `None`. |
| `--signed-sumstats` | input metadata | no | signed statistic specification | Identifies the signed statistic column and null value used to compute or validate `Z`; defaults to `None`. |
| `--info` | input metadata | no | INFO column hint | Identifies one raw INFO column for INFO filtering; defaults to `None`. |
| `--info-list` | input metadata | no | INFO column list hint | Identifies comma-separated INFO columns/tokens for INFO filtering; defaults to `None`. |
| `--nstudy` | input metadata | no | study-count column hint | Identifies the raw study-count column; defaults to `None`. |
| `--nstudy-min` | QC | no | minimum study count | Optional minimum study-count row filter; defaults to `None`. |
| `--ignore` | input metadata | no | ignored raw columns | Column names to ignore during munger detection and parsing; defaults to `None`. |
| `--a1-inc` | input mode | no | allele ordering compatibility | Legacy allele-ordering switch forwarded to the munging kernel; defaults to `False`. |
| `--keep-maf` | output mode | no | preserve allele frequency | Retains parsed allele frequency in compatibility outputs when available; defaults to `False`. |
| `--source-genome-build` | input metadata | no | raw coordinate source build | Defaults to `auto`, which infers hg19/hg38 from raw `CHR`/`POS` data in coordinate-family modes. May be set explicitly to `hg19` or `hg38`; rejected in rsid-family modes. |
| `--output-genome-build` | output metadata | yes in coordinate-family modes | final munged coordinate build | Defaults to omitted/`None`; required for `chr_pos`-family normal runs and `--infer-only`, and rejected in rsid-family modes. If it differs from the resolved source build, exactly one liftover method is required. |
| `--liftover-chain-file` | input | no | optional munger liftover chain | Uses a source-to-target chain file for coordinate-only sumstats liftover; defaults to omitted/`None`. Mutually exclusive with `--use-hm3-quick-liftover`. |
| `--use-hm3-quick-liftover` | input mode | no | packaged HM3 coordinate map | Uses the curated dual-build HM3 map for HM3-only quick liftover; defaults to `False`. Requires `--use-hm3-snps` and is mutually exclusive with `--liftover-chain-file`. |
| `--snp-identifier` | config | no | provenance | Defaults to `chr_pos_allele_aware`. Coordinate-family munger runs use `--source-genome-build` and `--output-genome-build`; rsid-family munger runs reject genome-build and liftover build flags and store `genome_build=None`. Allele-aware modes require usable `A1/A2`; rerun with `--snp-identifier chr_pos` or `--snp-identifier rsid` to run without allele-aware identity. |
| `--log-level` | logging | no | workflow log verbosity | Controls ordinary LDSC logger record verbosity; defaults to `INFO`; these records go to `diagnostics/sumstats.log` and the CLI console (stderr) shows only errors. Lifecycle audit lines always appear in the file. |
| `--overwrite` | output mode | no | collision policy | Controls whether fixed sumstats outputs may be replaced; defaults to `False`, so any owned `sumstats.*` artifact is refused. With overwrite, stale sibling formats not produced by the current `--output-format` are removed after a successful run. |

Removed flags: `--sumstats`, `--sumstats-file` for raw munge input,
`--merge-alleles`, `--merge-alleles-file`, `--no-alleles`, `--out`,
`--daner-old`, `--daner-new` (use `--format daner-old` / `--format daner-new`).

Fixed output names:

```text
<output_dir>/sumstats.sumstats.gz
<output_dir>/sumstats.parquet
<output_dir>/diagnostics/sumstats.log
<output_dir>/diagnostics/dropped_snps/dropped.tsv.gz
```

`sumstats.parquet` is the default curated artifact. `sumstats.sumstats.gz` is
written only for `--output-format tsv.gz` or `both`. `diagnostics/sumstats.log` is
preflighted and opened by the public workflow layer, but it is excluded from
`MungeRunSummary.output_paths`. Detailed provenance and output bookkeeping are
written to `diagnostics/sumstats.log`; row-level liftover drops are written to
`diagnostics/dropped_snps/dropped.tsv.gz`. There is no root metadata sidecar:
the downstream identity contract is embedded in the parquet footer when parquet
output is written. The kernel emits package logger records for QC progress and
returns an in-memory table; the workflow layer alone writes both Parquet and
the optional metadata-free gzip artifact.

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

### `ldsc query-r2`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--panel-dir` | input | yes | build-ref-panel output directory | Opens a package-built reference panel directory, optionally selecting a build subdirectory with `--genome-build`. This is the only panel input mode. |
| `--pairs` | input | yes | SNP-pair table | TSV by default, CSV when the filename ends with `.csv`, or stdin via `-`. Endpoint columns use `_1` and `_2` suffixes, with `SNP`, `CHR`, `POS`, `A1`, and `A2` supplied according to the active identifier mode. |
| `--output-dir` | output | no | result directory | Defaults to omitted/`None`; writes the canonical result directory (`query_r2.tsv` plus `diagnostics/metadata.json` and `diagnostics/query-r2.log`) when supplied. When omitted, the result table streams as a clean TSV to stdout. |
| `--overwrite` | output control | no | overwrite toggle | Replaces existing query-r2 output artifacts in `--output-dir`; defaults to `False`. No effect in stdout mode. |
| `--log-level` | output control | no | log verbosity | Verbosity of the directory-mode `diagnostics/query-r2.log`; defaults to `INFO`. No effect in stdout mode. |
| `--snp-identifier` | input metadata | no | query identity mode | Overrides panel metadata; defaults to omitted/`None`. If omitted, `R2Panel.open()` reads `ldsc:snp_identifier` from parquet metadata. |
| `--genome-build` | input selector | no | panel build selector | Concrete `hg19`/`hg38` selector for build-ref-panel directory layouts; defaults to omitted/`None`. |

Removed flags: `--out` (replaced by the `--output-dir` result directory with a
stdout fallback), `--meta`, `--parquet` (explicit single-chromosome input; use
`--panel-dir`), `--with-r` (signed `r` is now always emitted), `--strategy`,
`--strategy-threshold` (lookup strategy is fixed to the internal `auto` rule).

Output columns always append `r2`, nullable `sign`, signed Pearson `r`, and
`status` to the input pairs. `r` is computed by inverting unbiased R2 with panel
`ldsc:n_samples`; it is all-NaN in base/allele-blind modes (no sign) and in
panels lacking `ldsc:n_samples`. `status` is blank for found pairs and may report
`not_in_panel`, `cross_chromosome`, or `absent`.

### `ldsc h2`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--ldscore-dir` | input | yes | canonical LD-score result directory | Reads baseline LD scores and embedded `regression_ld_scores`, the historical `w_ld` component used when final h2 weights are computed. |
| `--sumstats-file` | input | yes | munged summary-statistics file | Exact path or exact-one glob. Current self-describing Parquet is native; legacy LDSC2 `.sumstats[.gz]` is automatically projected by rsID onto the canonical LD-score panel. Footerless Parquet is rejected. |
| `--trait-name` | input metadata | no | output trait label | Optional label override; defaults to omitted/`None`. If omitted, regression uses the sumstats parquet footer `ldsc:trait_name` when present, then the filename fallback. |
| `--output-dir` | output | no | result output directory | Selects where to write h2 results; defaults to omitted/`None`, so the CLI prints the compact `h2.tsv` schema to stdout and writes no files. |
| `--count-kind` | model | no | count vector choice | Selects the count vector used by regression; defaults to `common`, while `all` uses all-SNP counts. |
| `--n-blocks` | model | no | block jackknife partitions | Number of jackknife blocks used by the regression estimator; defaults to `200`. |
| `--no-intercept` | model | no | intercept policy | Fixes the LDSC intercept instead of estimating it; defaults to `False`. |
| `--allow-identity-downgrade` | model | no | identity compatibility override | Allows same-family allele-aware/base artifact mixes under the base identity mode; defaults to `False`. |
| `--intercept-h2` | model | no | fixed h2 intercept | Optional fixed h2 intercept supplied to the regression estimator; defaults to `None`. |
| `--two-step-cutoff` | model | no | two-step threshold | Optional cutoff for two-step regression fitting; defaults to `None`. Inclusive (`chi^2 <= cutoff` retained for step 1; deviates from legacy LDSC's strict `chi^2 < cutoff`). When unset, a single-annotation fit with a free intercept defaults to the legacy cutoff `30`. |
| `--chisq-max` | QC/model | no | chi-square filter | Optional maximum chi-square retained for regression fitting; defaults to `None`. Inclusive (`chi^2 <= chisq_max`; deviates from legacy LDSC's strict `chi^2 < chisq_max`). When unset, a single-annotation fit stays uncapped (outliers handled by the two-step estimator); a multi-annotation fit applies the legacy default cap `max(0.001 * N.max(), 80)`. |
| `--samp-prev` | model | no | sample (case) prevalence | Scalar sample prevalence `P` for liability-scale conversion of binary-trait h2; defaults to `None`. A probability in `(0, 1)`, or `nan` for a quantitative trait. Requires `--pop-prev`; omit both for observed scale. Output adds `total_h2_liab`/`total_h2_liab_se` and the applied prevalence columns. |
| `--pop-prev` | model | no | population prevalence | Scalar population prevalence `K` for liability-scale conversion; defaults to `None`. A probability in `(0, 1)`, or `nan`. Requires `--samp-prev`. Validated before inputs load. |
| `--log-level` | logging | no | workflow log verbosity | Controls ordinary LDSC logger record verbosity; defaults to `INFO`. With `output_dir`, records go to `diagnostics/h2.log` and the console (stderr) shows only errors; without it (console-only run) they print to the console. Lifecycle audit lines always appear in the file when one is created. |
| `--overwrite` | output mode | no | collision policy | Controls whether `h2.tsv` and diagnostics may be replaced; defaults to `False`, so an existing owned artifact is refused. |

Removed flags: `--ldscore`, `--counts`, `--w-ld`, `--annotation-manifest`,
`--sumstats`, `--out`.

### `ldsc partitioned-h2`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--ldscore-dir` | input | yes | canonical LD-score result directory | Requires an overlap artifact. With query columns, runs the cell-type regime (baseline plus one query per model); with no query columns, runs one functional-category model jointly over all baseline columns. This includes explicitly converted baseline-only LDSC2 suites. |
| `--sumstats-file` | input | yes | munged summary-statistics file | Exact path or exact-one glob. Accepts current Parquet or legacy LDSC2 text under the same projection rule as `h2`. |
| `--trait-name` | input metadata | no | output trait label | Optional label override; defaults to omitted/`None`. If omitted, regression uses the sumstats parquet footer `ldsc:trait_name` when present, then the filename fallback. |
| `--output-dir` | output | no | result output directory | Selects where to write partitioned-h2 results; defaults to omitted/`None`, so the CLI prints the `partitioned_h2.tsv` schema to stdout and writes no files. |
| `--count-kind` | model | no | count vector choice | Selects the count vector used by regression; defaults to `common`, while `all` uses all-SNP counts. |
| `--n-blocks` | model | no | block jackknife partitions | Number of jackknife blocks used by the regression estimator; defaults to `200`. |
| `--no-intercept` | model | no | intercept policy | Fixes the LDSC intercept instead of estimating it; defaults to `False`. |
| `--allow-identity-downgrade` | model | no | identity compatibility override | Allows same-family allele-aware/base artifact mixes under the base identity mode; defaults to `False`. |
| `--intercept-h2` | model | no | fixed h2 intercept | Optional fixed h2 intercept supplied to the regression estimator; defaults to `None`. |
| `--two-step-cutoff` | model | no | two-step threshold | Optional cutoff for two-step regression fitting; defaults to `None`. Inclusive (`chi^2 <= cutoff` retained for step 1). Not applied to partitioned (multi-annotation) models, which use the chi-square cap below instead. |
| `--chisq-max` | QC/model | no | chi-square filter | Optional maximum chi-square retained for regression fitting; defaults to `None`. Inclusive (`chi^2 <= chisq_max`; deviates from legacy LDSC's strict `chi^2 < chisq_max`). When unset, partitioned (multi-annotation) models apply the legacy default outlier cap `max(0.001 * N.max(), 80)` to keep extreme-chi-square SNPs from dominating the regression. |
| `--samp-prev` | model | no | sample (case) prevalence | Scalar sample prevalence `P` for liability-scale conversion; defaults to `None`. A probability in `(0, 1)`, or `nan` for a quantitative trait. Requires `--pop-prev`; omit both for observed scale. Adds the `*_liab` heritability columns (e.g. `category_h2_liab`/`category_h2_liab_se` and `total_h2_liab`/`total_h2_liab_se`) and the applied prevalence columns (proportions, enrichment, and coefficients are scale-invariant). |
| `--pop-prev` | model | no | population prevalence | Scalar population prevalence `K` for liability-scale conversion; defaults to `None`. A probability in `(0, 1)`, or `nan`. Requires `--samp-prev`. Validated before inputs load. |
| `--write-per-query-results` | output mode | no | per-query result tree | Requests per-query output folders under `diagnostics/query_annotations`; defaults to `False`, so only the aggregate table is returned/written. |
| `--summary-sort-by` | output mode | no | aggregate row sorting | Sort key for the partitioned-h2 table; defaults to `auto`, which resolves to `coefficient-p` in the cell-type regime (query annotations present) and `category` in the functional regime. Explicit choices: `category`, `prop-snps`, `prop-h2`, `enrichment`, `enrichment-p`, `coefficient`, and `coefficient-p`. |
| `--log-level` | logging | no | workflow log verbosity | Controls ordinary LDSC logger record verbosity; defaults to `INFO`. With `output_dir`, records go to `diagnostics/partitioned-h2.log` and the console (stderr) shows only errors; without it (console-only run) they print to the console. Lifecycle audit lines always appear in the file when one is created. |
| `--overwrite` | output mode | no | collision policy | Controls whether aggregate/per-query outputs and diagnostics may be replaced; defaults to `False`, so any owned partitioned-h2 artifact is refused. With overwrite, aggregate-only runs remove stale `diagnostics/query_annotations/` trees after successful writes. |

Removed flags: `--ldscore`, `--counts`, `--w-ld`, `--annotation-manifest`,
`--query-columns`, `--sumstats`, `--out`.

### `ldsc rg`

| Flag | Direction | Required | Object | Notes |
|---|---:|---:|---|---|
| `--ldscore-dir` | input | yes | canonical LD-score result directory | Reads baseline LD scores and embedded `regression_ld_scores`, the historical `w_ld` component used when final rg/gencov weights are computed. |
| `--sumstats-sources` | input | yes | two or more munged summary-statistics files | Accepts exact paths and glob patterns, including mixes of current Parquet and legacy LDSC2 text. Every legacy trait must have `A1/A2`; traits are projected independently and then pairwise harmonized. With two files, computes one pair; with three or more files and no anchor, computes all unordered pairs in input order. |
| `--anchor-trait` | input selector | no | anchor trait label or path | When supplied, first matches a resolved trait name, then a resolved input path; defaults to omitted/`None`, so all unordered pairs are computed. |
| `--output-dir` | output | no | result output directory | Selects where to write the rg output family; defaults to omitted/`None`. RG tables include nominal p-values only. Without it, Python returns `RgResultFamily` and the CLI prints only the concise `rg.tsv` schema to stdout. |
| `--write-per-pair-detail` | output mode | no | optional pair result tree | Requires `--output-dir`; defaults to `False`; writes `diagnostics/pairs/manifest.tsv` plus one `rg_full.tsv` and `metadata.json` per attempted pair when enabled. |
| `--count-kind` | model | no | count vector choice | Selects the count vector used by regression; defaults to `common`, while `all` uses all-SNP counts. |
| `--n-blocks` | model | no | block jackknife partitions | Number of jackknife blocks used by each regression estimator; defaults to `200`. |
| `--no-intercept` | model | no | intercept policy | Fixes LDSC intercepts instead of estimating them; defaults to `False`. |
| `--allow-identity-downgrade` | model | no | identity compatibility override | Allows same-family allele-aware/base artifact mixes under the base identity mode; defaults to `False`. |
| `--two-step-cutoff` | model | no | two-step threshold | Optional cutoff for two-step regression fitting; defaults to `None`. Inclusive (`chi^2 <= cutoff` retained for step 1; deviates from legacy LDSC's strict `chi^2 < cutoff`). When unset, a single-annotation rg fit with a free h2 intercept defaults to the legacy cutoff `30`. |
| `--chisq-max` | QC/model | no | chi-square filter | Optional opt-in rg filter on the product of statistics; defaults to `None`, so no default cap is applied. Inclusive (`Z1^2 * Z2^2 <= chisq_max^2`; deviates from legacy LDSC's strict `Z1^2 * Z2^2 < chisq_max^2`). |
| `--intercept-h2` | model | no | fixed h2 intercepts | Optional fixed h2 intercept value(s) supplied to per-trait h2 estimators; defaults to `None`. |
| `--intercept-gencov` | model | no | fixed genetic-covariance intercepts | Optional fixed intercept value(s) supplied to genetic-covariance estimators; defaults to `None`. |
| `--samp-prev` | model | no | per-trait sample prevalences | Comma-separated sample prevalences aligned to the resolved `--sumstats-sources` order, one per trait; defaults to `None`. Each a probability in `(0, 1)` or `nan` for a quantitative trait. Requires `--pop-prev`. Mutually exclusive with `--prevalence-manifest`. Adds `*_liab` and per-trait prevalence columns to `rg_full.tsv`/`h2_per_trait.tsv`; the `rg` ratio is unchanged. |
| `--pop-prev` | model | no | per-trait population prevalences | Comma-separated population prevalences aligned to the resolved order, one per trait; defaults to `None`. Each in `(0, 1)` or `nan`. Requires `--samp-prev`. Mutually exclusive with `--prevalence-manifest`. |
| `--prevalence-manifest` | model | no | prevalence lookup table | Whitespace/tab-delimited TSV with columns `trait_name`, `samp_prev`, `pop_prev` (`#` comment lines ignored); defaults to omitted/`None`. Looked up by exact munged trait name; may contain extra traits (every resolved trait must be present). Mutually exclusive with `--samp-prev`/`--pop-prev`. Duplicate resolved munged names abort the run. |
| `--log-level` | logging | no | workflow log verbosity | Controls ordinary LDSC logger record verbosity; defaults to `INFO`. With `output_dir`, records go to `diagnostics/rg.log` and the console (stderr) shows only errors; without it (console-only run) they print to the console. Lifecycle audit lines always appear in the file when one is created. |
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
| `AnnotationBuildConfig` | `query_annot_gene_list_sources` | input | one-column query gene-list group |
| `AnnotationBuildConfig` | `padding_bp` | input transform | BED or resolved-gene interval padding in base pairs; default `0` |
| `AnnotationBuildConfig` | `output_dir` | output | generated query annotation directory |
| `AnnotationBuilder.run(config=None, chrom=None)` | `config` | input/output | annotation workflow config; defaults to the builder config |
| `AnnotationBuilder.project_bed_annotations(...)` | `query_annot_bed_sources` | input | query BED group |
| `AnnotationBuilder.project_bed_annotations(...)` | `padding_bp` | input transform | BED interval padding in base pairs before projection; default `0` |
| `add_annotate_arguments(parser)` | `parser` | CLI surface | shared annotate argument registration for standalone and top-level parsers |
| `run_annotate_from_args(args)` / `main(argv)` | `query_annot_bed_sources` | input | query BED group |
| `run_annotate_from_args(args)` / `main(argv)` | `baseline_annot_sources` | input | baseline annotation templates |
| `run_annotate_from_args(args)` / `main(argv)` | `output_dir` | output | generated query annotation directory |
| `run_bed_to_annot(...)` | `query_annot_bed_sources` | input | query BED group |
| `run_bed_to_annot(...)` | `baseline_annot_sources` | input | baseline annotation templates |
| `run_bed_to_annot(...)` | `padding_bp` | input transform | BED interval padding in base pairs before projection; default `0` |
| `run_bed_to_annot(...)` | `output_dir` | output | generated query annotation directory; convenience wrapper writes `diagnostics/annotate.log` |

Removed Python names: `bed_paths`, `query_bed_paths`, `bed_files`,
`baseline_annot`, `bed_padding_bp`, `out_prefix`, `main_bed_to_annot`.

### LD-score calculation

| Object/function | Argument | Direction | Object |
|---|---:|---:|---|
| `RefPanelConfig` | `plink_prefix` | input | PLINK reference-panel prefix token |
| `RefPanelConfig` | `r2_dir` | input | preferred package-built parquet reference-panel directory |
| `RefPanelConfig` | `ref_panel_snps_file` | input | reference SNP restriction |
| `RefPanelConfig` | `keep_indivs_file` | input | PLINK individual keep file |
| `RefPanelConfig` | `maf_min` | input metadata | retained reference-panel MAF filter |
| `LDScoreConfig` | `regression_snps_file` | input | regression SNP override; omitted uses the bundled HM3 map |
| `LDScoreConfig` | `snp_batch_size` | performance | LD-score SNP batch size |
| `LDScoreConfig` | `threads` | performance | cross-chromosome worker processes |
| `LDScoreConfig` | `common_maf_min` | input metadata | common-SNP count threshold only |
| `LDScoreOutputConfig` | `output_dir` | output | canonical LD-score result directory |
| `run_ldscore(**kwargs)` | `baseline_annot_sources`, `query_annot_sources`, `query_annot_bed_sources`, `query_annot_gene_list_sources`, `padding_bp` | input | optional annotation sources and live BED/gene padding; query inputs require baseline sources, and no-annotation runs synthesize `base`. Explicit `padding_bp`, including `0`, is rejected unless a live BED or gene-list query is supplied and indexed mode is not used. |
| `run_ldscore(**kwargs)` | `plink_prefix`, `r2_dir` | input | reference-panel sources |
| `run_ldscore(**kwargs)` | `ref_panel_snps_file`, `regression_snps_file`, `exclude_regions` | input | independently selects the retained reference universe, regression rows, and regression-only named region subtraction |
| `run_ldscore(**kwargs)` | `output_dir` | output | canonical result directory; convenience wrapper writes `diagnostics/ldscore.log` |

Removed Python names: `bfile`, `r2_table`, `frqfile`, `keep`, `maf`,
`baseline_annot`, `query_annot`, `query_annot_bed`, `out`,
`r2_ref_panel_dir`, `ref_panel_dir`, `r2_sources`, `metadata_sources`, and
LD-score `chunk_size`; `use_hm3_ref_panel_snps`, `use_hm3_regression_snps`,
`exclude_regions_build`, and `exclude_regions_bed` are also removed.

### Gene LD-score index

| Object/function | Argument | Direction | Object |
|---|---:|---:|---|
| `GeneLDScoreIndexBuildConfig` | `baseline_annot_sources`, `plink_prefix` | input | baseline suite and PLINK panel |
| `GeneLDScoreIndexBuildConfig` | `output_dir` | output | complete immutable index directory |
| `GeneLDScoreIndexBuildConfig` | `genome_build`, `snp_identifier` | config | required explicit `hg19` and base `rsid|chr_pos` semantic identity |
| `GeneLDScoreIndexBuildConfig` | `padding_bp`, `gene_exclude_regions`, `ld_wind_cm` | model | closed index scientific identity; `padding_bp` defaults to `0` |
| `GeneLDScoreIndexBuildConfig` | `regression_snps_file`, `exclude_regions` | input | regression row universe |
| `GeneLDScoreIndexBuildConfig` | `snp_batch_size`, `atom_batch_size`, `threads` | performance | bounded index computation |
| `build_gene_ldscore_index(config)` | `config` | input/output | staged build and publication workflow |

### LDSC2 LD-score conversion

| Object/function | Argument | Direction | Object |
|---|---:|---:|---|
| `LegacyLDScoreConverter.convert(...)` | `legacy_reference_dir`, `legacy_weight_dir` | input | selected complete LDSC2 suite directories |
| `LegacyLDScoreConverter.convert(...)` | `legacy_frequency_dir` | input | required baseline frequency suite; omitted for unpartitioned conversion |
| `LegacyLDScoreConverter.convert(...)` | `snp_identifier`, `genome_build` | config | allele-unaware output identity and reference coordinate provenance |
| `LegacyLDScoreConverter.convert(...)` | `output_dir` | output | canonical LDSC3 LD-score directory plus diagnostics |
| `convert_ldsc2_ldscores(...)` | same names | input/output | convenience wrapper with no configurable common threshold |

### Reference-panel building

| Object/function | Argument | Direction | Object |
|---|---:|---:|---|
| `ReferencePanelBuildConfig` | `plink_prefix` | input | PLINK reference-panel prefix token |
| `ReferencePanelBuildConfig` | `source_genome_build` | input metadata | PLINK source build; defaults to `"auto"` and is inferred from `.bim` before SNP restriction |
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
| `ReferencePanelBuildConfig` | `min_r2` | output mode | optional unbiased-R2 emission floor; default `0.0` writes every retained pair |
| `ReferencePanelBuildConfig` | `exclude_regions`, `exclude_regions_bed` | input transform | packaged or custom region exclusions applied on source-build PLINK coordinates before panel emission |
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
| `MungeConfig` | `sumstats_format` | input metadata | `"auto"`/`"plain"`/`"daner-old"`/`"daner-new"`; sole DANER selector (the `daner_old`/`daner_new` boolean fields are removed) |
| `MungeConfig` | `sumstats_snps_file` | input | summary-statistics SNP keep-list |
| `MungeConfig` | `use_hm3_snps` | input mode | packaged HM3 summary-statistics SNP restriction |
| `MungeConfig` | `source_genome_build` | input metadata | raw source build for `chr_pos`-family coordinates; defaults to `"auto"` |
| `MungeConfig` | `output_genome_build` | output metadata | required final build for `chr_pos`-family output coordinates |
| `MungeConfig` | `liftover_chain_file` | input | optional source-to-target chain file for munger liftover |
| `MungeConfig` | `use_hm3_quick_liftover` | input mode | use packaged curated HM3 dual-build map for coordinate-only liftover |
| `MungeConfig` | `output_dir` | output | munged output directory; `SumstatsMunger.run()` writes `diagnostics/sumstats.log` |
| `MungeConfig` | `output_format` | output mode | `parquet`, `tsv.gz`, or `both`; defaults to `parquet` |
| `SumstatsMunger.run(munge_config, ...)` | `munge_config` | input/output | normalized munging workflow; owns fixed output preflight, the self-describing `sumstats.parquet` footer, diagnostics, always-written `diagnostics/dropped_snps/dropped.tsv.gz`, and result construction; summary `output_paths` excludes logs |
| `SumstatsMunger.write_output(sumstats, output_dir, output_format='parquet')` | `output_dir` | output | writes fixed `sumstats.parquet` and/or `sumstats.sumstats.gz` |

Removed Python names: legacy separate source-path object field,
`MungeConfig.sumstats_file`, `MungeConfig.out_prefix`,
`write_output(..., out_prefix)`.

### R2 pair query

| Object/function | Argument | Direction | Object |
|---|---:|---:|---|
| `R2Panel.open()` | `panel_dir` | input | package-built reference-panel directory (sole input mode) |
| `R2Panel.open()` | `snp_identifier`, `genome_build` | input metadata | active identity mode and build selector |
| `R2Panel.query_pairs(pairs)` | `pairs` | input | endpoint-suffixed SNP-pair table; always appends `r2`, `sign`, `r`, `status` |
| `query_r2(...)` | `pairs`, `panel_dir` | input | one-shot query wrapper around `R2Panel.open()` and `query_pairs()` |
| `unbiased_r2_to_pearson_r(r2_adj, n, sign)` | `r2_adj`, `n`, `sign` | transform | converts adjusted R2 plus sign to Pearson `r` |

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
| `run_rg_from_args(args)` | `output_dir` | output | writes `rg.tsv`, `rg_full.tsv`, `h2_per_trait.tsv`, and diagnostics when supplied; otherwise prints compact `rg.tsv` to stdout; rg tables include nominal p-values only |
| `run_rg_from_args(args)` | `write_per_pair_detail` | output mode | optionally writes `diagnostics/pairs/manifest.tsv` and per-pair diagnostic folders when `output_dir` is supplied |

Removed Python/public argparse names: `sumstats`, `sumstats_1`, `sumstats_2`,
`out`, `ldscore`, `counts`, `w_ld`, `annotation_manifest`, `query_columns`.

Current curated `sumstats.parquet` artifacts provide
canonical `SNP`, `CHR`, `POS`, `Z`, and `N` fields when written by
`ldsc munge-sumstats`, plus an embedded footer identity payload with
`ldsc:artifact_type`, `ldsc:snp_identifier`, `ldsc:genome_build`, and
optional `ldsc:trait_name`. `.sumstats.gz` carries no embedded metadata.
`load_sumstats()` reconstructs the config snapshot from that footer payload.
Metadata-free `.sumstats` and `.sumstats.gz` files are marked as LDSC2 legacy
inputs and projected at the regression boundary; footerless Parquet is rejected.
Labels resolve as explicit override, then footer `trait_name`, then filename
fallback.
Liftover reports, coordinate provenance, selected curated output files, and
Parquet row groups are log provenance, not sidecar payload.
Regression therefore merges on the effective identity key: literal `SNP` in
`rsid`, `SNP:<allele_set>` in `rsid_allele_aware`, `CHR:POS` in `chr_pos`, and
`CHR:POS:<allele_set>` in `chr_pos_allele_aware`. Munger liftover is valid only in
`chr_pos`-family modes; it changes `CHR`/`POS`, never `SNP`, and runs after
`sumstats_snps_file` filtering. Liftover drop counts are written as readable log
records, examples appear only at `DEBUG`, and row-level dropped-SNP details are
written to `diagnostics/dropped_snps/dropped.tsv.gz`. There is no root metadata
sidecar; the parquet footer carries the minimal identity fields
`ldsc:artifact_type`, `ldsc:snp_identifier`, `ldsc:genome_build`, and optional
`ldsc:trait_name`.

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
