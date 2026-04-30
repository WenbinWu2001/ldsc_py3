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
- existing fixed output files raise `FileExistsError` before the workflow writes
  anything
- pass `--overwrite` on the CLI or `overwrite=True` in Python to intentionally
  replace those fixed files

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

Group-style inputs include annotation files, parquet R2 tables, BED inputs, and
frequency metadata. A single argument may resolve to multiple files.

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

- `run_bed_to_annot(...)`
- `AnnotationBuilder.project_bed_annotations(...)`

Accepted path forms:

- `query_annot_bed_sources`: exact path or glob
- `baseline_annot_sources`: exact path, glob, or explicit `@`

How files are handled:

- every resolved BED file becomes one annotation column
- every resolved baseline annotation file is used as a SNP template

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
  resulting bundle.
- Existing `query.<chrom>.annot.gz` files are refused before any annotation
  shard is written unless `overwrite=True` or CLI `--overwrite` is supplied.

### LD score calculation

Relevant APIs:

- `run_ldscore(...)`
- `ldsc ldscore`

Group-style inputs:

- `baseline_annot_sources`
- `query_annot_sources`
- `query_annot_bed_sources`
- `r2_sources`
- `metadata_sources`

Scalar-style inputs:

- `ref_panel_snps_file`
- `regression_snps_file`

PLINK prefix input:

- `plink_prefix`

How they are handled:

- group-style inputs may resolve to many files
- per chromosome, the workflow first tries filename-based filtering and then
  falls back to row-level `CHR` filtering
- `plink_prefix` must resolve to exactly one PLINK prefix per chromosome pass

Requirements:

- annotation files must align on SNP rows within each chromosome
- parquet and frequency inputs must have the schema expected by the downstream readers

Examples:

```bash
ldsc ldscore \
  --output-dir out/trait_ldscores \
  --r2-sources "r2/reference.*.parquet" \
  --metadata-sources "r2/reference_metadata.@.tsv.gz" \
  --r2-bias-mode unbiased \
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
- LD-score calculation writes `manifest.json`, `baseline.parquet`, and
  optional `query.parquet` inside that directory.
- `baseline.parquet` and `query.parquet` remain flat parquet files, but each
  row group contains exactly one chromosome. The manifest records
  `row_group_layout`, `baseline_row_groups`, and `query_row_groups`.
- Existing canonical LD-score files are refused before any of them are written
  unless `--overwrite` or `LDScoreOutputConfig(overwrite=True)` is supplied.
- `baseline.parquet` contains `CHR`, `POS`, `SNP`, `regr_weight`, and baseline
  LD-score columns. When both baseline and query inputs are omitted, the
  LD-score workflow writes a synthetic all-ones baseline column named `base`.
- `query.parquet` is present only when query annotations were supplied and
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
  using explicit CLI identifier/build flags or the registered `GlobalConfig`

How they are handled:

- `plink_prefix` is resolved at the PLINK prefix level, not at the individual `.bed/.bim/.fam` file level
- a chromosome suite such as `panel_chr@` is expanded one chromosome at a time
- liftover chains are optional; the matching source-to-target chain enables
  cross-build R2 and metadata outputs, while no matching chain produces
  source-build-only outputs
- genetic maps are required for every emitted build when `--ld-wind-cm` is set;
  SNP- and kb-window builds may omit maps and write emitted metadata `CM` as
  `NA`
- `snp_identifier` for SNP restrictions comes from the explicit argument or
  registered `GlobalConfig`; in `chr_pos` mode, `genome_build` selects the
  restriction-file coordinate build, falling back to the source build when
  omitted or `auto`

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
- Before chromosome processing starts, the builder checks the deterministic
  candidate paths under `{build}/r2` and `{build}/meta`.
- Existing candidate parquet or metadata files are refused unless
  `--overwrite` or `ReferencePanelBuildConfig(overwrite=True)` is supplied.
- The check covers source-build artifacts and covers target-build artifacts
  only when the matching liftover chain is configured.

### Sumstats munging and regression

Relevant APIs:

- `MungeConfig.sumstats_file`
- regression artifact paths such as `sumstats_file`, `sumstats_1_file`, and
  `sumstats_2_file`
- `ldscore_dir` for the canonical LD-score result directory

Accepted path forms:

- exact path
- exact-one glob for sumstats file inputs
- literal directory for `ldscore_dir`

How they are handled:

- these are scalar-style inputs
- if a glob expands to more than one file, the workflow raises instead of combining them
- `ldscore_dir` is not glob-resolved; it is opened as a directory containing
  `manifest.json` plus parquet payload files

Output:

- `ldsc munge-sumstats` writes `sumstats.sumstats.gz`, `sumstats.log`, and
  `sumstats.metadata.json` under `output_dir`; existing files are refused
  unless `--overwrite` or `MungeConfig(overwrite=True)` is supplied.
- `ldsc h2`, `ldsc partitioned-h2`, and `ldsc rg` write `h2.tsv`,
  `partitioned_h2.tsv`, and `rg.tsv`, respectively, when `output_dir` is
  provided; existing files are refused unless `--overwrite` is supplied.
- `ldsc partitioned-h2 --write-per-query-results` also writes a staged
  `query_annotations/` tree under `output_dir`; existing final per-query
  output is refused unless `--overwrite` is supplied.
- Existing output directories are valid in every case. Only known files for the
  active command are checked.

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
- the source build for `ldsc build-ref-panel`; that workflow still requires
  explicit `--source-genome-build`

Programmatic build inference is available from the top-level Python API:
`from ldsc import infer_chr_pos_build, resolve_chr_pos_table`. The command-line
API exposes inference only through existing workflow flags; there is no
standalone `ldsc infer-build` command.

## Important Edge Cases

These are the main cases where users can get confused or introduce bugs.

- Ambiguous glob for an exact-one input
  Example: `data/*.sumstats.gz` matching two files will raise.
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
