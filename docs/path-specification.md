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
- public workflows do not use output prefixes; output filenames inside
  `output_dir` are fixed by the workflow

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

- `AnnotationSourceSpec.baseline_annot_paths`
- `AnnotationSourceSpec.query_annot_paths`
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
AnnotationSourceSpec(
    baseline_annot_paths="annotations/baseline.*.annot.gz",
    query_annot_paths="annotations/query.@.annot.gz",
)
```

### BED to annotation projection

Relevant APIs:

- `run_bed_to_annot(...)`
- `AnnotationBuilder.project_bed_annotations(...)`

Accepted path forms:

- `query_annot_bed_paths`: exact path or glob
- `baseline_annot_paths`: exact path, glob, or explicit `@`

How files are handled:

- every resolved BED file becomes one annotation column
- every resolved baseline annotation file is used as a SNP template

Requirements:

- BED basenames must be unique because they become annotation names
- baseline templates must be `.annot` or `.annot.gz` files

Example:

```python
run_bed_to_annot(
    query_annot_bed_paths="beds/*.bed",
    baseline_annot_paths="annotations/baseline_chr/baseline.@.annot.gz",
    output_dir="annotations/query_from_beds",
)
```

### LD score calculation

Relevant APIs:

- `run_ldscore(...)`
- `ldsc ldscore`

Group-style inputs:

- `baseline_annot_paths`
- `query_annot_paths`
- `query_annot_bed_paths`
- `r2_paths`
- `metadata_paths`

Scalar-style inputs:

- `ref_panel_snps_path`
- `regression_snps_path`

PLINK prefix input:

- `plink_path`

How they are handled:

- group-style inputs may resolve to many files
- per chromosome, the workflow first tries filename-based filtering and then
  falls back to row-level `CHR` filtering
- `plink_path` must resolve to exactly one PLINK prefix per chromosome pass

Requirements:

- annotation files must align on SNP rows within each chromosome
- parquet and frequency inputs must have the schema expected by the downstream readers

Examples:

```bash
ldsc ldscore \
  --output-dir out/trait_ldscores \
  --baseline-annot-paths "annotations/baseline.@.annot.gz" \
  --r2-paths "r2/reference.*.parquet" \
  --metadata-paths "r2/reference_metadata.@.tsv.gz" \
  --ld-wind-cm 1.0
```

```bash
ldsc ldscore \
  --output-dir out/trait_ldscores \
  --baseline-annot-paths "annotations/baseline_joint*.annot.gz" \
  --plink-path "resources/example_1kg_30x/genomes_30x_chr@" \
  --ld-wind-kb 100
```

Output:

- `--output-dir` is a literal directory destination.
- LD-score calculation writes `manifest.json`, `baseline.parquet`, and
  optional `query.parquet` inside that directory.
- `baseline.parquet` contains `CHR`, `SNP`, `BP`, `regr_weight`, and baseline
  LD-score columns.
- `query.parquet` is present only when query annotations were supplied and
  contains `CHR`, `SNP`, `BP`, and query LD-score columns.
- Regression commands consume this directory via `--ldscore-dir`; users do not
  pass count vectors, weight files, or annotation manifests.

### Reference-panel building

Relevant APIs:

- `run_build_ref_panel(...)`
- `ldsc build-ref-panel`

Accepted path forms:

- `plink_path`: exact PLINK prefix, exact-one PLINK glob, or explicit `@` PLINK suite token
- map and chain inputs: exact path or exact-one glob

How they are handled:

- `plink_path` is resolved at the PLINK prefix level, not at the individual `.bed/.bim/.fam` file level
- a chromosome suite such as `panel_chr@` is expanded one chromosome at a time

Example:

```bash
ldsc build-ref-panel \
  --plink-path data/reference/genomes_30x_chr@ \
  --source-genome-build hg38 \
  --genetic-map-hg19-path maps/hg19.txt \
  --genetic-map-hg38-path maps/hg38.txt \
  --liftover-chain-hg38-to-hg19-path chains/hg38ToHg19.over.chain \
  --ld-wind-cm 1.0 \
  --output-dir out/ref_panel
```

### Sumstats munging and regression

Relevant APIs:

- `RawSumstatsSpec.sumstats_path`
- regression artifact paths such as `sumstats_path`, `sumstats_1_path`, and
  `sumstats_2_path`
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

## Automatic Inference

The package makes a limited amount of automatic inference after path resolution:

- column alias inference
  Example: `BP` may be accepted as `POS`
- chromosome normalization
  Example: `chr1` may be normalized to `1`
- row-level chromosome filtering inside parsed annotation and metadata files

The package does not infer:

- missing file suffixes
- directory contents from an input directory argument
- a hidden per-chromosome mode from a bare prefix

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

## Recommended Usage

- Prefer explicit `@` suite tokens for chromosome-sharded inputs
- Prefer globs only when you really mean “all matching files”
- Keep chromosome labels visible in filenames if you want predictable per-chromosome file selection
- Use exact paths or exact-one globs for scalar inputs
