# Input Path Specification

This note explains how to specify filesystem inputs in the refactored package.
The goal is practical: help you choose the right path form for each workflow and
avoid ambiguous inputs.

## Supported Path Forms

Public input tokens support exactly three forms:

- Exact path
  Example: `annotations/baseline.1.annot.gz`
- Standard Python glob
  Example: `annotations/baseline.*.annot.gz`
- Explicit chromosome-suite token with `@`
  Example: `annotations/baseline.@.annot.gz`

Not supported:

- Bare prefixes without `@`
  Example: `annotations/baseline.`
- Automatic suffix guessing
  Example: `annotations/baseline.@` when the real files are `.annot.gz`
- Directory-only discovery for input files
  Example: `annotations/baseline_chr/`

Output paths are different:

- `output_dir` is treated as a literal directory
- if an output directory does not exist, the program warns once and creates it
- if the directory already exists, it is reused
- public workflows do not use output prefixes; output filenames inside
  `output_dir` are fixed by the workflow
- existing workflow-owned output artifacts raise `FileExistsError` before the
  workflow writes anything
- pass `--overwrite` on the CLI or `overwrite=True` in Python to intentionally
  replace those fixed files; for coherent result-directory workflows, a
  successful overwrite also removes stale owned siblings that the current run
  did not produce

### Coherent output artifact families

Several workflows write a fixed family of files that share one run identity
through `output_dir`. These families are treated as one coherent set, not as
independent optional files:

- `munge-sumstats`: `metadata.json`, `sumstats.parquet`,
  `sumstats.sumstats.gz`, `diagnostics/dropped_snps/dropped.tsv.gz`, and
  `diagnostics/sumstats.log`
  for CLI/workflow runs
- `ldscore`: `metadata.json`, `ldscore.baseline.parquet`,
  `ldscore.query.parquet`, and `diagnostics/ldscore.log` for CLI/workflow runs
- `build-ref-panel`: `{hg19,hg38}/chr*_r2.parquet`,
  `{hg19,hg38}/chr*_meta.tsv.gz`,
  `diagnostics/metadata.json`, `diagnostics/dropped_snps/chr*_dropped.tsv.gz`,
  and `diagnostics/build-ref-panel*.log` for CLI/workflow runs
- `partitioned-h2`: `partitioned_h2.tsv`,
  `diagnostics/metadata.json`, optional `diagnostics/query_annotations/`, and
  `diagnostics/partitioned-h2.log` for CLI/workflow runs
- `h2`: `h2.tsv`, `diagnostics/metadata.json`, and
  `diagnostics/h2.log` for CLI/workflow runs
- `annotate`: root-level `query.<chrom>.annot.gz` shards, plus diagnostic
  metadata, dropped-SNP audit, and `annotate.log` under `diagnostics/`
- `rg`: `rg.tsv`, `rg_full.tsv`, `h2_per_trait.tsv`,
  `diagnostics/metadata.json`, optional `diagnostics/pairs/`, and
  `diagnostics/rg.log` for CLI/workflow runs

Without overwrite, any existing current-contract owned sibling in the family
rejects the run, even if that sibling is not selected by the current output
mode. With overwrite enabled, the workflow writes the requested current outputs
and then removes stale current-contract owned siblings not produced by the
successful run. Removed legacy root diagnostic names are ignored by preflight
and cleanup, and unrelated files in the directory are preserved. Directory artifacts such as
`diagnostics/query_annotations/` and `diagnostics/pairs/` are owned as whole
trees: no-overwrite blocks if the root exists, and overwrite swaps or removes
the complete tree after the current run succeeds.

Direct Python data writers enforce the data artifact family they own. Workflow
wrappers add their workflow log to the preflight family.

## General Resolution Rules

### Exact path

Use this when one argument should point to one known file.

- For exact-one inputs, resolution must produce exactly one file
- For group inputs, an exact path just contributes one file to the group

### Glob

Use this when one argument should expand to multiple files.

- Matches are sorted and deduplicated while preserving first-seen order
- Group-style inputs may expand to many files
- Exact-one inputs fail if the glob expands to more than one file

### `@` chromosome suite

Use this when you have one file or prefix per chromosome.

- The workflow substitutes `@` with chromosome labels such as `1`, `2`, ..., `22`
- Only explicit `@` is supported
- This is the clearest way to express chromosome-sharded inputs

## How Group Inputs Are Handled

Group-style inputs include annotation files and BED inputs. The public
reference-panel parquet interface is directory-oriented (`r2_dir`); the loader
discovers fixed `chr*_r2.parquet` files and optional metadata sidecars inside
that directory.

During per-chromosome processing:

1. The workflow resolves the token or tokens to a concrete file list.
2. If a token expands to multiple files, the workflow tries to keep only files
   whose filenames encode the active chromosome, for example `.1.`, `_1_`, or `chr1`.
3. If filename-based filtering is not possible, the matched files are still used,
   and the file reader later filters rows by the `CHR` column.

This means both of the following are valid:

- one file per chromosome
  Example: `annotations/baseline.@.annot.gz`
- one or more shared multi-chromosome files
  Example: `annotations/baseline_joint*.annot.gz`

## Module-Specific Behavior

### Annotation loading

Relevant APIs:

- `AnnotationBuildConfig.baseline_annot_sources`
- `AnnotationBuildConfig.query_annot_sources`
- `AnnotationBuilder.run(...)`

Accepted path forms:

- exact path
- glob
- explicit `@`

How files are handled:

- all resolved files contributing to the same chromosome are combined column-wise
- after sorting, the SNP rows must match exactly for that chromosome
- annotation column names must be unique across files

Requirements:

- each file must contain metadata columns equivalent to `CHR`, `POS`, `SNP`, and `CM`
- query/baseline files for the same chromosome must refer to the same ordered SNP universe

Automatic inference:

- common header aliases for required columns are inferred
- chromosome labels are normalized
- multi-chromosome files are filtered by `CHR` during parsing

Example:

```python
AnnotationBuildConfig(
    baseline_annot_sources="annotations/baseline.*.annot.gz",
    query_annot_sources="annotations/query.@.annot.gz",
)
```

### BED to annotation projection

Relevant APIs:

- `annotation_builder.main(argv)`
- `run_annotate_from_args(args)`
- `run_bed_to_annot(...)`
- `AnnotationBuilder.project_bed_annotations(...)`

Accepted path forms:

- `query_annot_bed_sources`: exact path or glob
- `baseline_annot_sources`: exact path, glob, or explicit `@`

How files are handled:

- every resolved BED file becomes one annotation column
- every resolved baseline annotation file is used as a SNP template
- CLI dispatch through `ldsc annotate` calls the same workflow module directly;
  parsed namespaces are not converted back to argv and reparsed

Requirements:

- BED basenames must be unique because they become annotation names
- baseline templates must be `.annot` or `.annot.gz` files

Example:

```python
run_bed_to_annot(
    query_annot_bed_sources="beds/*.bed",
    baseline_annot_sources="annotations/baseline_chr/baseline.@.annot.gz",
    output_dir="annotations/query_from_beds",
    overwrite=True,
)
```

Output:

- `output_dir` is created when missing and reused when present.
- Projection writes `query.<chrom>.annot.gz` for every chromosome in the
  resulting bundle, plus diagnostics under `diagnostics/`.
- Existing root-level `query.*.annot.gz` files or owned diagnostics are refused
  before any annotation shard is written unless `overwrite=True` or CLI
  `--overwrite` is supplied. With overwrite enabled, stale query shards outside
  the current chromosome set are removed after the current shards are written.

### LD score calculation

Relevant APIs:

- `run_ldscore(...)`
- `ldsc ldscore`

Group-style inputs:

- `baseline_annot_sources`
- `query_annot_sources`
- `query_annot_bed_sources`

Directory-style inputs:

- `r2_dir`

Scalar-style inputs:

- `ref_panel_snps_file`
- `regression_snps_file`

Packaged HM3 convenience flags:

- `use_hm3_ref_panel_snps`
- `use_hm3_regression_snps`

PLINK prefix input:

- `plink_prefix`

How they are handled:

- group-style inputs may resolve to many files
- per chromosome, the workflow first tries filename-based filtering and then
  falls back to row-level `CHR` filtering
- `plink_prefix` must resolve to exactly one PLINK prefix per chromosome pass

Requirements:

- annotation files must align on SNP rows within each chromosome
- parquet R2 directories must use the fixed `chr{chrom}_r2.parquet` naming
  contract, with optional `chr{chrom}_meta.tsv.gz` sidecars
- package-built R2 parquet files record `ldsc:r2_bias` and `ldsc:n_samples` in
  schema metadata, so `--r2-bias-mode` and `--r2-sample-size` are only needed
  for external raw-R2 parquet files without LDSC metadata

Examples:

```bash
ldsc ldscore \
  --output-dir out/trait_ldscores \
  --r2-dir "r2_ref_panel_1kg30x_1cM_hm3/hg38" \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
```

```bash
ldsc ldscore \
  --output-dir out/trait_ldscores \
  --baseline-annot-sources "annotations/baseline_joint*.annot.gz" \
  --plink-prefix "resources/example_1kg_30x/genomes_30x_chr@" \
  --common-maf-min 0.05 \
  --ld-wind-kb 100
```

Output:

- `--output-dir` is a literal directory destination.
- LD-score calculation writes `metadata.json`, `ldscore.baseline.parquet`,
  optional `ldscore.query.parquet`, and `diagnostics/ldscore.log` inside that
  directory.
- `ldscore.baseline.parquet` and `ldscore.query.parquet` remain flat parquet files, but each
  row group contains exactly one chromosome. Root `metadata.json` records
  `row_group_layout`, `baseline_row_groups`, and `query_row_groups`.
- Existing canonical LD-score files or `diagnostics/ldscore.log` are refused before any of
  them are written unless `--overwrite` or
  `LDScoreOutputConfig(overwrite=True)` is supplied. With overwrite enabled, a
  successful baseline-only run removes any stale `ldscore.query.parquet` sibling.
- `ldscore.baseline.parquet` contains `CHR`, `POS`, `SNP`,
  `regression_ld_scores`, and baseline LD-score columns. `regression_ld_scores`
  is the historical `w_ld` LD score over the regression SNP universe, not the
  final h2/rg regression weight. When both baseline and query inputs are
  omitted, the LD-score workflow writes a synthetic all-ones baseline column
  named `base`.
- `ldscore.query.parquet` is present only when query annotations were supplied and
  contains `CHR`, `POS`, `SNP`, and query LD-score columns.
- Regression commands consume this directory via `--ldscore-dir`; users do not
  pass count vectors, weight files, or annotation manifests.

### Reference-panel building

Relevant APIs:

- `run_build_ref_panel(...)`
- `ldsc build-ref-panel`

Accepted path forms:

- `plink_prefix`: exact PLINK prefix, exact-one PLINK glob, or explicit `@` PLINK suite token
- map and chain inputs, when provided: exact path or exact-one glob
- `ref_panel_snps_file`, when provided: scalar file-like token interpreted
  using `GlobalConfig.snp_identifier`; restriction files may omit alleles and
  then match by base key, while allele-bearing restrictions in allele-aware
  modes match by the effective allele-aware key; packaged HM3 is allele-bearing
  and participates in allele-aware matching; duplicate restriction keys collapse
  to one retained key and non-identity columns such as `CM` or `MAF` are ignored;
  `chr_pos`-family coordinates must be aligned to the PLINK source build
- `use_hm3_snps`, when set: uses the packaged curated HM3 map instead of an
  explicit `ref_panel_snps_file`

How they are handled:

- `plink_prefix` is resolved at the PLINK prefix level, not at the individual `.bed/.bim/.fam` file level
- a chromosome suite such as `panel_chr@` is expanded one chromosome at a time
- liftover chains are optional; the matching source-to-target chain enables
  cross-build R2 and metadata outputs in `chr_pos`-family modes, while no matching
  chain produces source-build-only outputs; `use_hm3_quick_liftover` also emits
  the opposite build for the HM3-restricted coordinate universe and requires
  `use_hm3_snps`; all reference-panel liftover is rejected in `rsid`-family modes
- genetic maps are required for every emitted build when `--ld-wind-cm` is set;
  SNP- and kb-window builds may omit maps and write emitted metadata `CM` as
  `NA`
- `snp_identifier` for SNP restrictions comes from `GlobalConfig`; the CLI flag
  constructs a one-invocation identifier config, and the Python wrapper reads
  the registered config
- `build-ref-panel` ignores `GlobalConfig.genome_build`; in `chr_pos`-family modes, the
  restriction file must provide source-build coordinates, either through a
  source-specific column such as `hg19_POS` or through generic `POS` that
  infers to the source PLINK build

Example:

```bash
ldsc build-ref-panel \
  --plink-prefix data/reference/genomes_30x_chr@ \
  --source-genome-build hg38 \
  --genetic-map-hg19-sources maps/hg19.txt \
  --genetic-map-hg38-sources maps/hg38.txt \
  --liftover-chain-hg38-to-hg19-file chains/hg38ToHg19.over.chain \
  --ld-wind-cm 1.0 \
  --output-dir out/ref_panel
```

Output:

- `--output-dir` is created when missing and reused when present.
- Before chromosome processing starts, the builder checks both current-run
  deterministic paths and existing workflow-owned siblings under `hg19/`,
  `hg38/`, plus owned diagnostics under `diagnostics/`.
- Existing parquet, metadata, dropped-SNP audit, or workflow-log files are
  refused unless `--overwrite` or `ReferencePanelBuildConfig(overwrite=True)`
  is supplied.
- With overwrite enabled, a successful run removes stale owned source-build,
  target-build, out-of-scope chromosome, dropped-SNP, or log siblings that the
  current run did not produce.
- Dropped-SNP audit sidecars are always written for processed chromosomes
  (header-only when clean) and contain liftover-stage rows with reasons
  `source_duplicate`, `unmapped_liftover`, `cross_chromosome_liftover`, and
  `target_collision`.
- The ref-panel directory may contain many chromosome/build artifacts, but the
  workflow still treats them as one owned family for stale-output protection.

### Sumstats munging and regression

Relevant APIs:

- `MungeConfig.raw_sumstats_file`
- regression artifact paths such as `sumstats_file` and rg
  `sumstats_sources`
- `ldscore_dir` for the canonical LD-score result directory

Accepted path forms:

- exact path
- exact-one glob for scalar sumstats file inputs
- multi-file glob or mixed path/glob list for rg `--sumstats-sources`
- literal directory for `ldscore_dir`
- `ldsc munge-sumstats --infer-only` resolves only `--raw-sumstats-file`;
  `--output-dir` is not required because no artifacts are written

How they are handled:

- h2 and partitioned-h2 sumstats are scalar-style inputs; if a glob expands to
  more than one file, the workflow raises instead of combining them
- rg `--sumstats-sources` is a group-style input; all resolved files are
  deduplicated in first-seen order and then used for pair selection
- `ldscore_dir` is not glob-resolved; it is opened as a directory containing
  `metadata.json` plus parquet payload files
- `ldsc munge-sumstats` uses `--format auto` by default after the raw path is
  resolved. The inference layer can detect plain text, old DANER, new DANER,
  and PGC VCF-style headers before applying the usual column aliases and
  repair suggestions.

Output:

- `ldsc munge-sumstats` writes `sumstats.parquet` by default, plus root
  `metadata.json`, `diagnostics/sumstats.log`, and
  `diagnostics/dropped_snps/dropped.tsv.gz` under `output_dir`;
  `--output-format tsv.gz` writes legacy `sumstats.sumstats.gz`, and
  `--output-format both` writes both curated artifacts. Existing owned
  `sumstats.*` artifacts are refused unless `--overwrite` or
  `MungeConfig(overwrite=True)` is supplied. With overwrite enabled, a
  successful run removes stale sibling formats not produced by the current
  `--output-format`.
  `diagnostics/sumstats.log` is not recorded in
  `MungeRunSummary.output_paths`; detailed provenance and output bookkeeping
  are written to the log, row-level liftover drops are written to
  `diagnostics/dropped_snps/dropped.tsv.gz`, and root `metadata.json` stays
  thin.
- `--use-hm3-snps` uses the packaged curated HM3 map as the sumstats SNP
  restriction and conflicts with `--sumstats-snps-file`. HM3 quick liftover
  requires `--use-hm3-snps`.
- `ldsc h2`, `ldsc partitioned-h2`, and `ldsc rg` write fixed result families
  when `output_dir` is provided. Without `output_dir`, each command prints its
  compact TSV table to stdout and writes no diagnostics. For h2, the written
  family is `h2.tsv`, `diagnostics/metadata.json`, and workflow-owned
  `diagnostics/h2.log`. For rg, that family is `rg.tsv`, `rg_full.tsv`,
  `h2_per_trait.tsv`, optional `diagnostics/pairs/`, and workflow-owned
  `diagnostics/rg.log`; existing owned artifacts are refused unless
  `--overwrite` is supplied.
- `ldsc partitioned-h2` requires the LD-score directory to include
  `ldscore.query.parquet` and non-empty `query_columns`; baseline-only LD-score
  directories are valid for `h2` and `rg`.
- `ldsc partitioned-h2 --write-per-query-results` also writes a staged
  `diagnostics/query_annotations/` tree under `output_dir`. The tree contains
  `manifest.tsv` and one folder per query annotation, with per-query
  `partitioned_h2.tsv`, `partitioned_h2_full.tsv`, and `metadata.json`.
  Existing final per-query output is refused unless `--overwrite` is supplied;
  with overwrite enabled, an aggregate-only run removes a stale
  `diagnostics/query_annotations/` tree.
- Existing output directories are valid in every case. Only known files for the
  workflow-owned artifact family are checked; unrelated files are preserved.

## Automatic Inference

The package makes a limited amount of automatic inference after path resolution:

- column alias inference
  Example: `#CHROM` and `CHROM` may be accepted as `CHR`; `BP` may be accepted
  as `POS`
- raw sumstats metadata-line skipping
  Example: leading `##` lines are skipped before parsing the real header, while
  `#CHROM` remains a valid header column
- chromosome normalization
  Example: `chr1` may be normalized to `1`
- genome-build and coordinate-basis inference for `chr_pos` tables when a
  workflow is run with `genome_build="auto"` or `--genome-build auto`
  Example: an hg38 0-based `CHR`/`POS` restriction table can be inferred and
  converted to canonical 1-based `CHR:POS` identifiers
- row-level chromosome filtering inside parsed annotation and metadata files

The package does not infer:

- missing file suffixes
- directory contents from an input directory argument
- a hidden per-chromosome mode from a bare prefix
- target-build SNP restrictions for `ldsc build-ref-panel`; that workflow
  applies restrictions before liftover and requires restriction coordinates to
  align to the inferred or explicit source PLINK build

Programmatic build inference is available from the top-level Python API:
`from ldsc import infer_chr_pos_build, resolve_chr_pos_table`. The command-line
API exposes inference only through existing workflow flags; there is no
standalone `ldsc infer-build` command.

## Important Edge Cases

These are the main cases where users can get confused or introduce bugs.

- Ambiguous glob for an exact-one input
  Example: `data/*.parquet` or `data/*.sumstats.gz` matching two files will raise.
- Filename does not encode chromosome
  A multi-file glob may still work, but the workflow cannot pre-filter by name
  and must rely on `CHR` inside the file.
- Overlapping files for the same chromosome
  If multiple matched files contribute to the same chromosome, they are all
  combined. Their SNP rows must match exactly, or the run raises.
- Duplicate annotation column names
  If two files contribute the same annotation column name, the run raises.
- Mixed sharded and shared files
  This can be valid, but only if the resulting per-chromosome SNP rows still align.
- Existing output file in a reused directory
  The workflow fails before writing and tells you to pass `--overwrite` or
  `overwrite=True`. This prevents accidental reruns from silently replacing
  previous results.

## Recommended Usage

- Prefer explicit `@` suite tokens for chromosome-sharded inputs
- Prefer globs only when you really mean “all matching files”
- Keep chromosome labels visible in filenames if you want predictable per-chromosome file selection
- Use exact paths or exact-one globs for scalar inputs
- Use one dedicated output directory per reproducible run, and use
  `--overwrite` / `overwrite=True` only when replacing that run is intentional
