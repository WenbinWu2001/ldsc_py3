# Data Flow

This document summarizes the user-visible file streams for each public workflow. The diagrams use Mermaid `flowchart LR` because it maps cleanly onto the package's left-to-right data movement and layered module boundaries.

## Layer Legend

| Layer | Public? | Modules | Role |
| --- | --- | --- | --- |
| CLI | yes | `ldsc.cli` | parse subcommands and dispatch |
| Preprocessing | yes | `ldsc.config`, `ldsc.path_resolution`, `ldsc.column_inference`, `ldsc.chromosome_inference` | normalize tokens, headers, identifiers, chromosome order |
| Workflow | yes | feature modules under `src/ldsc/` | build aligned in-memory tables |
| Kernel | no | `ldsc._kernel.*` | low-level readers and numerical work |
| Postprocessing | yes | `ldsc.outputs`, pandas writers in `ldsc.regression_runner`, `ldsc._logging` | preflight fixed output paths, emit files, summaries, and workflow audit logs |

## Package Overview

```mermaid
flowchart LR
  IN1[BED intervals]
  IN2[PLINK or parquet R2 reference]
  IN3[Raw GWAS sumstats]

  subgraph CLI[CLI Layer (public)<br/>ldsc.cli]
    C1[Dispatch subcommands]
  end

  subgraph WF[Workflow Layer (public)<br/>src/ldsc/*.py]
    W1[Build annotations]
    W2[Build parquet panel]
    W3[Compute LD scores]
    W4[Munge sumstats]
    W5[Run regression]
  end

  subgraph K[Compute Kernel (private)<br/>src/ldsc/_kernel/*.py]
    K1[Parse files]
    K2[Compute LD]
    K3[Run estimators]
  end

  subgraph OUT[Postprocessing (public)<br/>ldsc.outputs / regression_runner]
    O1[Write LDSC artifacts]
    O2[Write summary tables]
    O3[Write workflow logs]
  end

  IN1 --> C1 --> W1 --> K1
  IN2 --> C1 --> W2 --> K2
  IN1 --> C1 --> W3 --> K1
  IN2 --> C1 --> W3 --> K2
  IN3 --> C1 --> W4 --> K1
  W3 --> O1
  W4 --> O1
  W1 --> O3
  W2 --> O3
  W3 --> O3
  W4 --> O3
  O1 --> W5 --> K3 --> O2
  W5 --> O3
```

## Shared SNP Identity Contract

All workflows share one SNP identity contract. Public `snp_identifier` values
are exactly `rsid`, `rsid_allele_aware`, `chr_pos`, and
`chr_pos_allele_aware`; the default is `chr_pos_allele_aware`. Mode names are
exact. Column aliases apply only to input headers, not to mode values.

Base modes are allele-blind. `rsid` uses only `SNP`, and `chr_pos` uses only
`CHR:POS`; any allele columns present in base-mode inputs are passive data and
do not affect identity, duplicate filtering, retention, or drop reasons.
Allele-aware modes use unordered, strand-aware `A1/A2` allele sets only to make
merge keys safer. They require usable alleles on sumstats, reference-panel
artifacts, R2 parquet endpoints, and LD-score artifacts, and they drop missing,
invalid/non-SNP, identical, strand-ambiguous, multi-allelic base-key, and
duplicate effective-key clusters. Artifact duplicate filtering always computes
the effective key for the active mode, then drops all rows in duplicate-key
clusters.
Allele-free summary-statistics munging is selected by using the base
`--snp-identifier rsid` or `--snp-identifier chr_pos` mode, not by a separate
allele-skip flag.

Restriction files may omit alleles. Allele-free restrictions match by base key
and can retain multiple candidate rows before later artifact cleanup. Packaged
HM3 restrictions are allele-free base-key filters. Allele-bearing restrictions
in allele-aware modes match by the effective allele-aware key. Annotation files
may also omit alleles in allele-aware modes because they describe genomic
membership; when annotation alleles are present, they participate in
allele-aware matching.

## 1. `annotate`: BED Projection To SNP-Level `.annot.gz`

Output directories are created when missing and reused when present. Existing
root-level `query.*.annot.gz` files and the workflow `annotate.log` are refused
before any query shard is written unless the caller passes `--overwrite` or
`overwrite=True`. With overwrite enabled, stale query shards outside the
current chromosome set are removed after the current shards are written.
Annotation rows are cleaned before BED projection by computing the effective
SNP identity key for the active mode and dropping all rows in duplicate-key
clusters. Runs with an output directory also write
`dropped_snps/dropped.tsv.gz`, header-only when no rows are dropped.

`ldsc.annotation_builder` owns both command entry paths: `main(argv)` parses
standalone annotation arguments, and `run_annotate_from_args(args)` consumes
the namespace created by the unified `ldsc` CLI without reparsing.

### Required inputs

| File | Example | Notes |
| --- | --- | --- |
| baseline annotation shard | `CHR POS SNP CM base`<br/>`1 10583 rs58108140 0.0 1` | used as the SNP template; one shard per chromosome is typical |
| BED file | `chr1 10000 10100 enhancer_A` | one or more region files; BED basenames become query annotation column names |

### Flow

```mermaid
flowchart LR
  I1[Baseline .annot(.gz)]
  I2[BED intervals]

  subgraph P1[Preprocessing (public)<br/>config + path_resolution]
    A1[Resolve input shards]
    A2[Normalize BED and baseline tokens]
  end

  subgraph W1[Workflow (public)<br/>ldsc.annotation_builder]
    A3[Load baseline rows]
    A8[Drop duplicate identity clusters]
    A4[Project BEDs to SNP columns]
    A7[Preflight query shards + dropped sidecar + annotate.log<br/>Write outputs]
  end

  subgraph K1[Kernel (private)<br/>ldsc._kernel.annotation]
    A5[Normalize BED files]
    A6[Intersect regions with SNP grid]
  end

  I1 --> A1 --> A3 --> A8 --> A4 --> A5 --> A6 --> A7
  I2 --> A1
  A2 --> A4
  A7 --> O1[query.<chrom>.annot.gz + dropped_snps/dropped.tsv.gz + annotate.log]
```

### Outputs

| File | Example | Notes |
| --- | --- | --- |
| projected query annotation shard | `CHR POS SNP CM enhancer_A`<br/>`1 10583 rs58108140 0.0 1` | output name is `query.<chrom>.annot.gz` |
| dropped-SNP audit sidecar | `CHR SNP source_pos target_pos reason base_key identity_key allele_set stage` | always written as `dropped_snps/dropped.tsv.gz`; records annotation identity cleanup rows |
| workflow log | plain-text lifecycle and package records | `annotate.log` under `output_dir`; not included in returned data paths |

### Modules used

- Preprocessing: `ldsc.config`, `ldsc.path_resolution`, `ldsc.column_inference`
- Workflow: `ldsc.annotation_builder`
- Kernel: `ldsc._kernel.annotation`
- Postprocessing: gzip writer inside `ldsc.annotation_builder`

## 2. `build-ref-panel`: PLINK To Standard Parquet R2 Reference

Before chromosome processing starts, the builder precomputes flat candidate
paths under each emitted `{build}` directory, always-written per-chromosome
dropped-SNP audit files under `dropped_snps/`, plus `build-ref-panel.log`.
Existing candidates are refused unless `--overwrite` or
`ReferencePanelBuildConfig(overwrite=True)` is supplied; unrelated files in the
output directory are left untouched. Unlike the result-directory workflows,
`build-ref-panel` does not clean stale optional target-build or out-of-scope
chromosome siblings from earlier configurations; use a fresh output directory
when changing emitted builds, liftover/coordinate configuration, or chromosome
scope.
Reference-panel liftover is coordinate behavior: chain-file liftover and HM3
quick liftover are valid only when the active SNP identifier mode is in the
`chr_pos` family.
HM3 quick liftover requires the packaged HM3 SNP restriction flag and emits the
opposite build only for the retained HM3 coordinate universe. Duplicate-position
filtering also applies only in `chr_pos`-family modes and always drops all colliding
source or target coordinate groups. The sidecar also records unmapped and
cross-chromosome liftover drops; clean processed chromosomes get a header-only
sidecar.

### Required inputs

| File | Example | Notes |
| --- | --- | --- |
| PLINK prefix | `reference/genomes_30x_chr22` | resolves to `.bed`, `.bim`, `.fam` |
| `.bim` row | `22 rs123 0.0 16050075 A G` | variant metadata |
| `.fam` row | `fam1 iid1 0 0 0 -9` | sample metadata |
| genetic map, conditional | `chr position Genetic_Map(cM)`<br/>`22 16050000 0.42` | required for every emitted build when cM windows are used; optional for SNP/kb windows |
| liftover method, optional | `hg38ToHg19.over.chain.gz` or `--use-hm3-snps --use-hm3-quick-liftover` | matching source-to-target chain enables cross-build R2 and metadata in chr_pos-family modes (`chr_pos`, `chr_pos_allele_aware`); HM3 quick liftover uses the packaged map and requires HM3 restriction; omitted liftover produces source-build-only output; liftover is rejected in rsID-family modes (`rsid`, `rsid_allele_aware`) |
| keep or restrict file, optional | one IID per row, a headered SNP table, or `--use-hm3-snps` | filters individuals or variants; SNP restriction matching uses `GlobalConfig.snp_identifier`; `chr_pos`-family restrictions must match the source PLINK build; allele-free restrictions, including HM3, match by base key; allele-bearing restrictions in allele-aware modes match by effective allele-aware key |

### Flow

```mermaid
flowchart LR
  I1[PLINK prefix]
  I2[Genetic maps]
  I3[Liftover / keep / restrict]

  subgraph P2[Preprocessing (public)<br/>config + path_resolution]
    B1[Resolve PLINK prefixes]
    B0[Infer source build from .bim when omitted]
    B2[Resolve map and filter files]
    B8[Interpret SNP restrictions in source build]
  end

  subgraph W2[Workflow (public)<br/>ldsc.ref_panel_builder]
    B3[Discover chromosomes]
    B4[Prepare build state]
  end

  subgraph K2[Kernel (private)<br/>ldsc._kernel.ref_panel_builder]
    B5[Interpolate cM positions]
    B6[Emit pairwise R2 rows]
    B7[Format standard schemas]
  end

  I1 --> B1 --> B0 --> B3 --> B4 --> B5 --> B6 --> B7
  I2 --> B2 --> B4
  I3 --> B2 --> B8 --> B4
  B7 --> O2[{build}/chr*_r2.parquet + {build}/chr*_meta.tsv.gz + dropped_snps/ + build-ref-panel.log]
```

### Outputs

| File | Example | Notes |
| --- | --- | --- |
| build-specific R2 parquet | `hg38/chr22_r2.parquet` with columns `CHR`, `POS_1`, `POS_2`, `SNP_1`, `SNP_2`, `R2`, plus endpoint `A1_1/A2_1/A1_2/A2_2` when available | one row per unordered SNP pair inside the LD window; endpoint allele columns are required in allele-aware modes; row groups are sorted by that build's `POS_1`; schema metadata records minimal identity provenance plus `ldsc:n_samples` and `ldsc:r2_bias` |
| build-specific runtime metadata sidecar | `hg38/chr22_meta.tsv.gz` with leading `# ldsc:*` identity-provenance comments and `CHR POS SNP CM MAF A1 A2` when alleles are available | authoritative SNP universe for the matching R2 parquet when present; current package-written sidecars carry `schema_version`, `artifact_type`, `snp_identifier`, and `genome_build`; `A1/A2` are required in allele-aware modes |
| dropped-SNP audit sidecar | `dropped_snps/chr22_dropped.tsv.gz` with `CHR SNP source_pos target_pos reason base_key identity_key allele_set stage` | always written for each processed chromosome; header-only when no rows were dropped; identity reasons include `missing_allele`, `invalid_allele`, `strand_ambiguous_allele`, `multi_allelic_base_key`, and `duplicate_identity`; liftover reasons include `source_duplicate`, `unmapped_liftover`, `cross_chromosome_liftover`, and `target_collision` |
| workflow log | plain-text lifecycle and package records | `build-ref-panel.log` under `output_dir`; not included in `ReferencePanelBuildResult.output_paths` |

### Modules used

- Preprocessing: `ldsc.config`, `ldsc.path_resolution`
- Workflow: `ldsc.ref_panel_builder`
- Kernel: `ldsc._kernel.ref_panel_builder`, `ldsc._kernel.ldscore`, `ldsc._kernel.formats`
- Postprocessing: parquet and TSV writers in the kernel

## 3. `ldscore`: Reference Panel And Optional Annotations To LDSC Artifacts

The canonical LD-score workflow preflights `manifest.json`,
`ldscore.baseline.parquet`, `ldscore.query.parquet`, and `ldscore.log` as one
owned family before writing any of them. Use `--overwrite` or
`LDScoreOutputConfig(overwrite=True)` only for intentional reruns. With
overwrite enabled, a successful baseline-only run removes stale
`ldscore.query.parquet`.
The parquet payloads remain single flat files, but each row group contains rows
from exactly one chromosome. The manifest records the row-group layout and
per-chromosome offsets so readers can load one chromosome without scanning the
whole table.

For ordinary unpartitioned LD scores, callers may omit both baseline and query
inputs. The workflow then creates a synthetic baseline annotation named exactly
`base`, with value `1.0` for every row returned by the retained reference-panel
metadata. Query annotations are partitioned-LDSC inputs and require explicit
baseline annotations.

### Required inputs

| File | Example | Notes |
| --- | --- | --- |
| baseline annotation shard, optional | `CHR POS SNP CM base`<br/>`1 10583 rs58108140 0.0 1` | optional for unpartitioned runs; required when query annotations are supplied |
| query annotation shard, optional | `CHR POS SNP CM enhancer_A`<br/>`1 10583 rs58108140 0.0 1` | optional extra annotation columns; valid only with explicit baseline annotations |
| PLINK prefix or parquet R2 panel | `panel_chr@` or build directory `ref_panel/hg38` | choose one backend |
| frequency / metadata sidecar, optional | `CHR POS SNP CM MAF A1 A2` | used for MAF and runtime metadata; `A1/A2` are required for allele-aware modes |
| regression SNP list, optional | `rs123` or `CHR POS` table | restricts the weight-table SNP set; allele columns may be omitted and then match by base key; allele-bearing restrictions in allele-aware modes match by effective allele-aware key |

### Flow

```mermaid
flowchart LR
  I1[Optional .annot(.gz) shards]
  I2[PLINK or parquet R2 reference]
  I3[Metadata / restriction files]

  subgraph P3[Preprocessing (public)<br/>config + path_resolution + column_inference]
    C1[Resolve shards and tokens]
    C2[Normalize identifier mode]
  end

  subgraph W3[Workflow (public)<br/>ldsc.ldscore_calculator]
    C3[Bundle annotations or synthesize base]
    C4[Select backend]
    C5[Aggregate chromosome results]
  end

  subgraph K3[Kernel (private)<br/>ldsc._kernel.ldscore]
    C6[Compute chromosome LD scores]
    C7[Compute regression-universe LD score]
  end

  subgraph O3[Postprocessing (public)<br/>ldsc.outputs]
    C8[Write LDSC artifacts]
  end

  I1 --> C1 --> C3 --> C4 --> C6 --> C7 --> C5 --> C8
  I2 --> C1 --> C4
  I3 --> C1 --> C2 --> C6
```

### Outputs

| File | Example | Notes |
| --- | --- | --- |
| baseline LD-score table | `CHR POS SNP regression_ld_scores base`<br/>`1 10 rs1 1.7 1.2` | `ldscore.baseline.parquet` inside `output_dir`; `regression_ld_scores` is historical `w_ld`, not the final h2/rg regression weight; one row group per chromosome |
| query LD-score table | `CHR POS SNP enhancer_A`<br/>`1 10 rs1 0.4` | `ldscore.query.parquet` inside `output_dir`; one row group per chromosome; omitted when no query annotations exist |
| manifest | JSON metadata with files, columns, counts, chromosomes, config, row counts, and row-group metadata | `manifest.json` inside `output_dir` |
| workflow log | plain-text lifecycle and package records | `ldscore.log` inside `output_dir`; not included in `LDScoreResult.output_paths` |

### Modules used

- Preprocessing: `ldsc.config`, `ldsc.path_resolution`, `ldsc.column_inference`
- Workflow: `ldsc.ldscore_calculator`
- Kernel: `ldsc._kernel.ldscore`
- Postprocessing: `ldsc.outputs`

## 4. `munge-sumstats`: Raw GWAS Table To Curated Sumstats

The munging workflow preflights `sumstats.parquet`,
`sumstats.sumstats.gz`, `sumstats.log`, `sumstats.metadata.json`, and
`dropped_snps/dropped.tsv.gz` as one owned family before delegating to
`SumstatsMunger.run()` and then the
legacy-compatible munging kernel. The workflow owns the log file, metadata
sidecar, and curated output writing; the kernel keeps the low-level parsing and
QC. Before calling the kernel, the workflow runs format and column inference:
`--format auto` is the default and detects plain whitespace text, old DANER,
new DANER, and PGC VCF-style headers. `--infer-only` runs that inference pass
without requiring `--output-dir` and prints missing fields plus exact repair
suggestions. Default output is `sumstats.parquet`; `--output-format tsv.gz` or `both`
also supports the legacy `sumstats.sumstats.gz` artifact. With overwrite
enabled, stale sibling formats not produced by the current run are removed
after successful writes. Optional sumstats liftover is a `chr_pos`-family step:
the source build comes from `GlobalConfig.genome_build` or munger build
inference, SNP restrictions are interpreted in that source build before
liftover, and `--target-genome-build` must be paired with exactly one method
when the target differs from the source. Chain-file liftover uses
`--liftover-chain-file`; HM3 quick liftover requires `--use-hm3-snps`, uses the
packaged curated `hm3_curated_map.tsv.gz`, and is coordinate-only, so it never
rewrites `SNP`.
Missing coordinates, unmapped hits, cross-chromosome hits, and duplicate
source/target coordinate groups are dropped. Counts are readable audit records
in `sumstats.log`; row-level drops are written to
`dropped_snps/dropped.tsv.gz`; examples appear only at `DEBUG`.

### Required inputs

| File | Example | Notes |
| --- | --- | --- |
| raw sumstats | `#CHROM POS ID EA NEA PVAL BETA NEFF`<br/>`1 754182 rs3131969 A G 0.46 0.004 829249.58` | leading `##` metadata lines are skipped; header aliases are normalized in the workflow layer; `NEFF` is not inferred as `N` unless the user explicitly passes `--N-col NEFF` |
| DANER or PGC VCF-style raw sumstats, optional schema mode | old DANER: `FRQ_A_<Ncas>` and `FRQ_U_<Ncon>` headers<br/>new DANER: exact `Nca` and `Nco` columns<br/>PGC VCF-style: leading `##` metadata and `#CHROM` header | `--format auto` detects these profiles; explicit `--format daner-old`, `--format daner-new`, or `--format pgc-vcf` overrides auto-detection; legacy `--daner-old`/`--daner-new` remain supported |
| sumstats SNP keep-list, optional | headered `SNP` or `CHR`/`POS` restriction file, or `--use-hm3-snps` | optional row filter loaded once before parsing and applied inside each retained chunk; allele-free restrictions, including HM3, match by base key before later identity cleanup; allele-bearing restrictions in allele-aware modes match by effective allele-aware key |
| sumstats liftover method, optional | `--target-genome-build hg38 --liftover-chain-file hg19ToHg38.over.chain` or `--target-genome-build hg38 --use-hm3-snps --use-hm3-quick-liftover` | valid only in `chr_pos`-family modes; updates `CHR`/`POS` after SNP restriction and preserves `SNP` labels |
| column hints, optional | `--snp ID --chr '#CHROM' --pos POS --a1 EA --a2 NEA` | useful when headers are ambiguous; common aliases infer automatically, and `--infer-only` reports the hints it would apply |
| INFO lists, optional | `IMPINFO=0.852,0.113,0.842,0.88,NA` | numeric/NA comma-separated per-study values are filtered on their mean; mixed nonnumeric lists such as `0.95,LOW,0.88` are rejected with `--ignore` / `--info-list` suggestions |

### Flow

```mermaid
flowchart LR
  I1[Raw GWAS table]
  I2[Column hints / SNP keep-list]

  subgraph P4[Preprocessing (public)<br/>config + path_resolution + column_inference]
    D1[Resolve one raw file]
    D2[Skip leading ## lines<br/>Detect format + normalize aliases]
  end

  subgraph W4[Workflow (public)<br/>ldsc.sumstats_munger]
    D3[Normalize CLI/API config<br/>Build typed munging args]
    D4[Own log + metadata + dropped_snps<br/>Capture run summary]
  end

  subgraph K4[Kernel (private)<br/>ldsc._kernel.sumstats_munger]
    D5[QC and infer columns]
    D6[Compute Z/N<br/>Finalize CHR/POS]
    D7[Optional source-to-target liftover<br/>drop missing, unmapped, and colliding coordinates]
  end

  I1 --> D1 --> D2 --> D3 --> D5 --> D6 --> D7 --> D4
  I2 --> D2
  D4 --> O4[sumstats.parquet by default<br/>optional sumstats.sumstats.gz + log + metadata JSON + dropped_snps]
```

### Outputs

| File | Example | Notes |
| --- | --- | --- |
| curated sumstats | `SNP CHR POS A1 A2 Z N`<br/>`rs3131969 1 754182 A G 0.74 829249.58` | written as `sumstats.parquet` by default under `output_dir`; `--output-format tsv.gz` writes legacy `sumstats.sumstats.gz`, and `both` writes both; `CHR`/`POS` are present and may be missing when absent from raw input; optional `FRQ` may also be present |
| log file | plain-text lifecycle, QC log, coordinate provenance, readable liftover reports, HM3 provenance, output bookkeeping, and count-level drop summaries | workflow-owned `sumstats.log` under `output_dir`, populated from package logger messages emitted during workflow orchestration and kernel QC; excluded from `MungeRunSummary.output_paths` |
| metadata sidecar | thin JSON with `schema_version`, `artifact_type`, `snp_identifier`, `genome_build`, and optional `trait_name` | written as `sumstats.metadata.json` under `output_dir`; used by `load_sumstats()` to reconstruct config provenance and trait labels |
| dropped-SNP audit sidecar | `CHR SNP source_pos target_pos reason base_key identity_key allele_set stage` | always written as `dropped_snps/dropped.tsv.gz`; header-only when no rows were dropped; reasons may include identity drops (`missing_allele`, `invalid_allele`, `strand_ambiguous_allele`, `multi_allelic_base_key`, `duplicate_identity`) and liftover drops (`missing_coordinate`, `source_duplicate`, `unmapped_liftover`, `cross_chromosome_liftover`, `target_collision`) |

### Modules used

- Preprocessing: `ldsc.config`, `ldsc.path_resolution`, `ldsc.column_inference`
- Workflow: `ldsc.sumstats_munger`
- Kernel: `ldsc._kernel.sumstats_munger`, `ldsc._kernel.liftover`
- Postprocessing: workflow-owned Parquet/TSV writing, log, and metadata writing

Allele columns keep the LDSC-compatible names `A1` and `A2`. `A1` means the
allele that the signed statistic is relative to; `A2` is the counterpart
allele. This is a signed-statistic convention, not a genome reference-allele
claim. Positive `Z`, positive `BETA`, positive `LOG_ODDS`, and `OR > 1` are
interpreted relative to `A1`.

## 5. `h2`, `partitioned-h2`, and `rg`: Curated Artifacts To Regression Summaries

Regression summary commands write fixed result artifacts when `output_dir` is
supplied. `h2` writes `h2.tsv` plus `h2.log`; `partitioned-h2` writes
`partitioned_h2.tsv` plus `partitioned-h2.log`; `rg` writes `rg.tsv`,
`rg_full.tsv`, `h2_per_trait.tsv`, optional `pairs/`, plus workflow-owned
`rg.log`. Existing owned TSVs, optional trees, or logs raise before the new
table is written unless the command includes `--overwrite`. For
`partitioned-h2`, `partitioned_h2.tsv`, `query_annotations/`, and
`partitioned-h2.log` are treated as one owned family; aggregate-only overwrites
remove stale `query_annotations/` after the new summary is written. Without
`output_dir`, regression commands return in-memory results and do not create
log files; the rg CLI additionally prints the concise `rg.tsv` schema to
stdout for quick runs.
`partitioned-h2` can also write an opt-in per-query tree with
`--write-per-query-results`; the aggregate `partitioned_h2.tsv` remains the
stable summary entry point. It requires query LD scores in the LD-score
directory; baseline-only directories are valid for `h2` and `rg` but are
rejected by `partitioned-h2`.
For rg, `--sumstats-sources` accepts two or more files; three or more files
produce all unordered pairs unless `--anchor-trait` selects one trait label or
source path for anchor-vs-rest estimation. `--write-per-pair-detail` adds the
optional `pairs/` detail tree when an `output_dir` is supplied.

Regression merges on the effective key for the resolved mode: `SNP` in `rsid`,
`SNP:<allele_set>` in `rsid_allele_aware`, `CHR:POS` in `chr_pos`, and
`CHR:POS:<allele_set>` in `chr_pos_allele_aware`. `--allow-identity-downgrade`
is regression-only; it allows same-family allele-aware/base mixes to run under
the base mode and logs the original modes plus duplicate-key rows dropped before
merge. rsID-family and coordinate-family modes never mix.

### Required inputs

| File | Example | Notes |
| --- | --- | --- |
| munged sumstats | `SNP CHR POS A1 A2 Z N`<br/>`rs1 1 754182 A G 1.96 1000` | one file for `h2` and `partitioned-h2`, two or more files for `rg`; a neighboring `sumstats.metadata.json` recovers config provenance and `trait_name` when present |
| LD-score directory | `manifest.json`, `ldscore.baseline.parquet`, optional `ldscore.query.parquet` | produced by the LD-score workflow and supplied as `ldscore_dir`; `partitioned-h2` requires `ldscore.query.parquet` and non-empty `query_columns`; current parquet files have chromosome-aligned row groups; package-written directories without current manifest identity provenance are rejected and must be regenerated |

### Flow

```mermaid
flowchart LR
  I1[Curated sumstats parquet or .sumstats.gz]
  I2[LD-score artifacts]

  subgraph P5[Preprocessing (public)<br/>path_resolution + column_inference]
    E1[Resolve scalar files or rg file groups]
    E2[Reload canonical headers]
  end

  subgraph W5[Workflow (public)<br/>ldsc.regression_runner]
    E3[Rebuild LDScoreResult]
    E4[Merge on effective SNP identity key]
    E5[Drop zero-variance columns]
  end

  subgraph K5[Kernel (private)<br/>ldsc._kernel.regression]
    E6[Estimate h2]
    E7[Estimate partitioned h2]
    E8[Estimate rg]
  end

  I1 --> E1 --> E4 --> E5
  I2 --> E1 --> E2 --> E3 --> E4
  E5 --> E6 --> O5a[h2.tsv + h2.log]
  E5 --> E7 --> O5b[partitioned_h2.tsv + partitioned-h2.log]
  E7 --> O5d[query_annotations/]
  E5 --> E8 --> O5c[rg.tsv + rg_full.tsv + h2_per_trait.tsv + optional pairs/ + rg.log]
```

### Outputs

| Subcommand | Output columns | Example |
| --- | --- | --- |
| `h2` | `trait_name`, `n_snps`, `total_h2`, `total_h2_se`, `intercept`, `intercept_se`, `mean_chisq`, `lambda_gc`, `ratio`, `ratio_se` | `trait 105234 0.18 0.03 1.02 0.01 1.11 1.05 0.08 0.03` |
| `partitioned-h2` | `Category`, `Prop._SNPs`, `Prop._h2`, `Enrichment`, `Enrichment_p`, `Coefficient`, `Coefficient_p` | `enhancer_A 0.02 0.14 7.0 0.003 0.012 0.001` |
| `rg` | `trait_1`, `trait_2`, `rg`, `rg_se`, `z`, `p` | `trait_a trait_b 0.42 0.09 4.7 2.6e-06` |

When `output_dir` is supplied, the same directory also receives the matching
workflow log. The log is not part of any returned result `output_paths` mapping.

When `partitioned-h2 --write-per-query-results` is supplied, the command also
writes `query_annotations/manifest.tsv` plus one ordinal-prefixed sanitized
folder per query annotation. Each query folder contains `partitioned_h2.tsv`,
`partitioned_h2_full.tsv`, and `metadata.json`.
For column definitions and interpretation, see
[partitioned-h2-results.md](partitioned-h2-results.md).

### Modules used

- Preprocessing: `ldsc.path_resolution`, `ldsc.column_inference`
- Workflow: `ldsc.regression_runner`, `ldsc.sumstats_munger.load_sumstats()`
- Kernel: `ldsc._kernel.regression`, `ldsc._kernel._jackknife`, `ldsc._kernel._irwls`
- Postprocessing: pandas TSV writers in `ldsc.regression_runner`; partitioned-h2
  and rg directory writing in `ldsc.outputs`
