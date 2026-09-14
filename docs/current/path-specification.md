# Input Path Specification

Last updated on: 2026-09-14

Munged data filenames use the filesystem-safe trait label when supplied: `<trait>.parquet` and optional `<trait>.sumstats.gz`. The `sumstats.parquet` and `sumstats.gz` names below describe runs without a trait label. See [munging output artifacts](munge-sumstats.md#output-artifacts) for naming and overwrite rules.

This note explains how to specify filesystem inputs in the refactored package.
The goal is practical: help you choose the right path form for each workflow and
avoid ambiguous inputs.

## Supported Path Forms

Public input file tokens support three forms (PLINK prefixes have the additional forms described below):

- Exact path
  Example: `annotations/baseline.1.annot.gz`
- Standard Python glob
  Example: `annotations/baseline.*.annot.gz`
- Explicit chromosome-suite token with `@`
  Example: `annotations/baseline.@.annot.gz`

Not supported:

- Bare prefixes without `@` for ordinary file inputs (PLINK is an exception)
  Example: `annotations/baseline.`
- Automatic suffix guessing
  Example: `annotations/baseline.@` when the real files are `.annot.gz`
- Directory-only discovery for input files
  Example: `annotations/baseline_chr/`

Output paths are different:

- `output_dir` is treated as a literal directory
- if an output directory does not exist, the program creates it and logs creation at INFO
- if the directory already exists, it is reused
- public workflows do not use output prefixes; output filenames inside
  `output_dir` are fixed by the workflow
- `ldsc plot` and `ldsc convert-h2-scale` are derived-result exceptions: they
  take a source result directory and write to its fixed nested `plots/` or
  `postprocessing/liability-scale/` root instead of accepting `--output-dir`
- existing workflow-owned output artifacts raise `FileExistsError` before the
  workflow writes anything
- pass `--overwrite` on the CLI or `overwrite=True` in Python to intentionally
  replace those fixed files; for coherent result-directory workflows, a
  successful overwrite also removes stale owned siblings that the current run
  did not produce

### PLINK prefix resolution

PLINK is an exception to the file-token rule against bare prefixes. `ldscore`, `build-gene-ldscore-index`, `build-r2-panel`, and direct Python PLINK panels use the shared `path_resolution.inspect_plink_inputs()` implementation.

- Accept an exact prefix (`panel.22`), a member path (`panel.22.bim`), a plain stem (`1000G.EUR.QC.`), a glob (`"panel.*"`), or an explicit `@` pattern (`"panel.@"`). Only `.bed`, `.bim`, and `.fam` are stripped as PLINK member extensions; `.22` remains part of a prefix.
- An exact selected prefix takes precedence over broader stem discovery. Otherwise a plain stem discovers matching PLINK members; a glob selects its actual matches. Every selected prefix must have a complete BED/BIM/FAM trio. Partially present trios are reported rather than silently omitted.
- Validated BIM `CHR` values establish chromosome identity. Filenames do not select chromosomes. A multi-chromosome trio can supply several chromosomes; different trios containing the same chromosome fail as ambiguous. An explicit `@` member must contain its declared chromosome only.
- Validation batches missing-member, malformed-trio, and chromosome-conflict issues before computation. BED header/size checks detect truncated or mismatched inputs. Workflows reuse the resulting chromosome-to-prefix mapping; workers do not reinterpret the original token.
- Coverage remains workflow-specific. Direct query LD-score `@` declarations require autosomes 1-22; plain/glob inputs select their actual validated scope and must match baseline coverage. Gene-index construction always requires autosomes 1-22. Generic group expansion and `build-r2-panel` retain existing `@` members, while rejecting incomplete selected trios.

[`inspect_plink_inputs()`](../../src/ldsc/path_resolution.py) returns prefixes, chromosome assignments, and issue records. `require_valid()` raises a combined error, optionally including missing required chromosomes. `resolve_plink_prefix_group()` retains the path-only group interface when no chromosome is requested; per-chromosome resolution uses the shared content validation. The direct workflow combines these issues with its annotation diagnostics. Failed gene-index attempts retain `plink_input_issues.tsv` under `.<output-name>.build-state/`; a later attempt archives that audit alongside the prior log. Failed `build-r2-panel` validation writes `diagnostics/plink_input_issues.tsv`.

### Pattern support in command help

Quote patterns so the package receives them intact. `*` matches filename text; `@`, where supported, substitutes chromosome numbers 1-22 by default. Use `*` to select available files for a chromosome subset where the workflow permits one. For example, `"baseline.*.annot.gz"` selects matching files without declaring that all 22 chromosome members must exist. The pattern must still match the intended files, and their validated contents determine the chromosome set.

The first suite-capable flag in each command explains this syntax. Later flags refer back to it only for compatible rules. Match cardinality and unsupported placeholders stay local to each description.

| Input family | Supported patterns and restrictions |
| --- | --- |
| Baseline/prebuilt-query annotations, PLINK prefixes, quantile target/reference metadata | Exact inputs, quoted `*` patterns, and `@` chromosome tokens. Direct query LD-score preflight and gene-index construction enforce their full-suite declarations. Generic group expansion tries chromosomes 1-22 and may retain only existing members; completeness depends on the consuming workflow. |
| Query BED files and focal gene-list sources | Exact files or quoted `*` patterns; `@` is not expanded. |
| Raw or curated scalar sumstats files; munging SNP-list and liftover-chain files | Exact file or quoted `*` pattern resolving to exactly one file; `@` is not expanded. |
| `rg --sumstats-sources` | Exact files or quoted `*` patterns selecting multiple files; at least two inputs are required. `@` is not expanded. |
| `ldscore --regr-snps-file` | Exact file or quoted `*` pattern resolving to one file; `@` is not expanded. |
| `build-r2-panel` SNP/individual restriction files | Exact file or quoted `*` pattern resolving to one file; `@` is not expanded. |
| `build-r2-panel` genetic-map sources | Exact files, quoted `*` patterns, or `@` suites. |
| `ldscore` and `build-gene-ldscore-index` genetic-map sources | Comma-separated exact paths. Their current readers do not expand `*` or `@`. |
| Gene catalogs, fixed control gene lists, `ldscore` reference-SNP/individual lists, gene-index regression-SNP/individual lists, and `build-r2-panel` liftover-chain files | Exact file paths. Pattern expansion is not implemented by these readers. |
| Input/output directories and `query-r2 --pairs` | Literal paths; no pattern expansion. `--pairs -` retains its standard-input meaning. |

Sources: [`path_resolution.resolve_scalar_path`, `resolve_file_group`, and PLINK resolvers](../../src/ldsc/path_resolution.py); [`_annotation_preflight.resolve_annotation_inputs`](../../src/ldsc/_annotation_preflight.py); [`_direct_annotation.prepare_direct_annotations`](../../src/ldsc/_direct_annotation.py); [`_quantile_inputs.prepare_quantile_statistics`](../../src/ldsc/_quantile_inputs.py); [`gene_list_resolver._expand_focal_gene_list_sources`](../../src/ldsc/gene_list_resolver.py); [`ref_panel_builder.ReferencePanelBuilder._prepare_build_state`](../../src/ldsc/ref_panel_builder.py); [`_kernel.ref_panel._resolve_genetic_map`](../../src/ldsc/_kernel/ref_panel.py); [`gene_ldscore_index._load_builder_genetic_map`](../../src/ldsc/gene_ldscore_index.py).

### Direct LD-score query scope

Direct gene-list, BED, and prebuilt-query `ldscore` runs use these rules for both PLINK and parquet-R² reference backends:

1. `@` declares all autosomes 1–22. Every required chromosome input must exist and validate; it does not mean “use whichever chromosomes are present.”
2. Ordinary globs select the files they actually match. Validated contents, not filenames, determine their chromosome set. For example, `"annotations/*.22.annot.gz"` does not guarantee chr22-only contents: any additional chromosomes in its matched files participate in scope.
3. Validated baseline and reference chromosome sets must match exactly. An `@` declaration paired with a chr22-only input group therefore fails alignment.
4. Every selected focal/control gene must lie within that shared scope. Selection means unique successfully resolved genes after explicit gene-region exclusions, before SNP-support filtering. A pathway need not contain genes on every covered chromosome. Any nonempty incompletely covered pathway or control fails the entire batch under both `strict` and `resolved-only`; it is never automatically truncated or skipped. Query BED regions and prebuilt-query SNPs must also be within scope.

Quote glob tokens on the command line so the package receives the pattern intact. Users own glob selection: a missing file may be undetectable if it disappears from the matches and the remaining baseline/reference artifacts consistently cover the same subset. Use `@` when complete-autosomal coverage is required.

For an intentional subset, use quoted `*` patterns or exact paths selecting matching baseline and PLINK chromosome sets. If a full suite was intended, restore missing files instead. With `--r2-dir`, use a literal directory containing the matching chromosome set. Every selected focal/control gene must still be covered; supply missing chromosome inputs or explicitly revise the gene lists. Replacing `@` with `*` does not filter genes or truncate pathways. Coverage errors and their workflow-log tracebacks include this repair guidance.

BED padding extends the coordinates already present in the input file. Use `--padding-bp 0` if the BED intervals already include the intended gene-region padding, so the regions are not padded twice.

Scope is logged and persisted in `diagnostics/chromosome_scope.json`; input failures use `diagnostics/input_issues.tsv`. See the [gene-list scope examples](gene-list-input-format.md#direct-query-chromosome-scope). Public indexes remain complete-autosomal. `_ldscore_preflight.inspect_direct_inputs()` enforces the direct-query contract before projection; other commands retain their own path contracts.

### Coherent output artifact families

Several workflows write a fixed family of files that share one run identity
through `output_dir`. These families are treated as one coherent set, not as
independent optional files:

- `munge-sumstats`: `sumstats.parquet` (self-describing footer),
  `sumstats.gz`, `diagnostics/dropped_snps/dropped.tsv.gz`, and
  `diagnostics/sumstats.log`
  for CLI/workflow runs
- `ldscore`: `metadata.json`, `ldscore.baseline.parquet`, optional single or numbered query Parquet files, optional `ldscore.overlap.parquet`, and `diagnostics/ldscore.log` for CLI/workflow runs
- `convert-ldsc2-ldscores`: canonical `metadata.json`,
  `ldscore.baseline.parquet`, optional `ldscore.overlap.parquet`,
  `diagnostics/conversion_issues.tsv.gz`, and
  `diagnostics/convert-ldsc2-ldscores.log`. This command deliberately accepts
  directory inputs and performs its own strict chromosomes 1-22 family
  discovery; that exception does not enable directory discovery elsewhere.
- `build-r2-panel`: `{hg19,hg38}/chr*_r2.parquet`,
  `{hg19,hg38}/chr*_meta.tsv.gz`,
  `diagnostics/metadata.json`, `diagnostics/metadata.chr*.json`,
  `diagnostics/dropped_snps/chr*_dropped.tsv.gz`,
  and `diagnostics/build-r2-panel*.log` for CLI/workflow runs
- `partitioned-h2`: `partitioned_h2.tsv`,
  `diagnostics/metadata.json`, optional `diagnostics/query_annotations/`, and
  `diagnostics/partitioned-h2.log` for CLI/workflow runs
- `h2`: `h2.tsv`, `diagnostics/ld_score_regression_bins.tsv`,
  `diagnostics/metadata.json`, and `diagnostics/h2.log` for CLI/workflow runs;
  the default `plots/` and `postprocessing/` roots are owned derived outputs
- `annotate`: root-level `query.<chrom>.annot.gz` shards, plus diagnostic
  metadata, dropped-SNP audit, and `annotate.log` under `diagnostics/`
- `rg`: `rg.tsv`, `rg_full.tsv`, `h2_per_trait.tsv`,
  `diagnostics/metadata.json`, optional `diagnostics/pairs/`, and
  `diagnostics/rg.log` for CLI/workflow runs; rg tables report nominal
  p-values only and do not include package-computed corrected p-value columns
- `plot`: one selected fixed PNG plus `diagnostics/metadata.json` and
  `diagnostics/plot.log` below `<result-dir>/plots/`
- `convert-h2-scale`: `h2_scale_conversion.tsv`, optional
  `h2_prevalence_sensitivity.png`, `diagnostics/metadata.json`, and
  `diagnostics/convert-h2-scale.log` below
  `<h2-result-dir>/postprocessing/liability-scale/`

Without overwrite, any existing current-contract owned sibling in the family
rejects the run, even if that sibling is not selected by the current output
mode. With overwrite enabled, the workflow writes the requested current outputs
and then removes stale current-contract owned siblings not produced by the
successful run. Removed legacy root diagnostic names are ignored by preflight
and cleanup, and unrelated files in the directory are preserved. Workflows that
can run independent shards into one directory may narrow the owned family to the
current shard. For `build-r2-panel`, a concrete chromosome prefix owns only
that chromosome's package, while an `@` chromosome-suite invocation owns the
full all-chromosome package. Directory artifacts such as
`diagnostics/query_annotations/` and `diagnostics/pairs/` are owned as whole
trees: no-overwrite blocks if the root exists, and overwrite swaps or removes
the complete tree after the current run succeeds.

Direct Python data writers enforce the data artifact family they own. Workflow
wrappers add their workflow log to the preflight family.

The default `plots/` root is an owned sibling of h2, partitioned-h2,
quantile-h2, and rg; the default `postprocessing/` root is an owned sibling of
h2. A core no-overwrite run rejects an orphaned owned root. A successful core
overwrite removes these derived roots after publishing the new core artifacts,
preventing a plot or conversion from silently describing an older result.
Python-only `output_dir=` overrides for plotting and conversion are unmanaged
and therefore outside this cleanup boundary.

An authorized overwrite that fails writes `RUN_FAILED.txt` in the applicable
output scope; concrete chromosome reference-panel attempts use
`RUN_FAILED.chr<chrom>.txt`. Marker creation adds no rollback, quarantine, or
action-order change. It does not enter scientific metadata or block a retry,
and a successful retry removes the applicable marker.

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

For chromosome-sharded query annotations, the canonical accepted layout is
**one query annotation file per chromosome, containing multiple annotation
columns; each column is one query annotation**:

```text
annotations/
  query.1.annot.gz   # CHR POS SNP CM query_A query_B ...
  query.2.annot.gz   # CHR POS SNP CM query_A query_B ...
  ...
```

Supply that suite with
`query_annot_sources="annotations/query.@.annot.gz"`. The query column names
and order must be identical in every chromosome shard. A layout with one
separately sharded suite per query is not accepted: after path expansion it
would contribute multiple query files for the same chromosome. Whole-genome,
unsharded query files may instead be supplied separately and are column-bound
after strict row-alignment validation.

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
- `run_annotate(...)`
- `parse_annotate_args(argv)`
- `AnnotationBuilder.run(source_spec, output_dir=...)`

Accepted path forms:

- `query_annot_bed_sources` or `query_annot_gene_list_sources`: exact paths or globs; standalone annotate accepts exactly one route
- `baseline_annot_sources`: exact path, glob, or explicit `@`

How files are handled:

- every resolved BED file becomes one annotation column
- the annotation column name is the resolved BED file basename with the final
  suffix removed (`pathlib.Path.stem`); directory names are ignored, so
  `/path1/annot.bed` and `/path2/annot.bed` both become `annot`
- every resolved baseline annotation file is used as a SNP template
- `padding_bp` / `--padding-bp` expands each BED interval on both
  sides before projection and clips starts at zero; the default `0` leaves
  intervals unchanged
- CLI dispatch through `ldsc annotate` calls the same workflow module directly;
  parsed namespaces are not converted back to argv and reparsed

Requirements:

- BED stems must be unique because they become annotation names; duplicate
  stems and clashes with existing annotation columns raise before projection or
  output writing
- baseline templates must be `.annot` or `.annot.gz` files

Example:

```python
run_annotate(
    query_annot_bed_sources="beds/*.bed",
    baseline_annot_sources="annotations/baseline_chr/baseline.@.annot.gz",
    output_dir="annotations/query_from_beds",
    padding_bp=0,
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
- `regr_snps_file`

These SNP restriction files are scalar identity filters, not chromosome-suite
inputs. They may resolve through exact-one globs where the workflow calls the
scalar resolver, but `@` is not a chromosome placeholder for these files. The
restriction universe is also materialized as an in-memory key set before
matching. `munge-sumstats` streams raw sumstats chunks and uses packed integer
keys for base `chr_pos` keep-lists, but it still keeps the complete restriction
key set in memory; `ldscore` and reference-panel workflows likewise do not
lazy-load restriction rows by chromosome.

Bundled HM3 is the default `ldscore` regression set. Supply
`regr_snps_file` only to replace that selected set; it never changes the
reference-panel universe.

PLINK prefix input:

- `plink_prefix`

How they are handled:

- group-style inputs may resolve to many files
- PLINK chromosome assignments come from validated BIM contents and are reused for every chromosome pass
- `plink_prefix` must assign each requested chromosome to exactly one complete trio

Requirements:

- annotation files must align on SNP rows within each chromosome
- parquet R2 directories must use the fixed `chr{chrom}_r2.parquet` naming
  contract, with a mandatory paired `chr{chrom}_meta.tsv.gz` sidecar for each
  chromosome (the sidecar defines the index space)
- package-built R2 parquet files record `ldsc:r2_bias` and `ldsc:n_samples` in
  schema metadata; R2 bias mode and sample size are read solely from this
  metadata (there are no bias-related flags), so external raw-R2 panels must
  declare `ldsc:r2_bias=raw` and `ldsc:n_samples` to be corrected

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
  --plink-prefix "resources/example_1kg_30x/genomes_30x_chr" \
  --common-maf-min 0.05 \
  --ld-wind-kb 100
```

Output:

- `--output-dir` is a literal directory destination.
- LD-score calculation writes `metadata.json`, `ldscore.baseline.parquet`, optional query Parquet files, optional `ldscore.overlap.parquet`, and `diagnostics/ldscore.log` inside that directory. One query batch uses `ldscore.query.parquet`; multiple batches use `ldscore.query.batch00001.parquet` and subsequent ordinals.
- Each Parquet file is genome-wide, with one row group per chromosome. Root metadata records `row_group_layout`, `baseline_row_groups`, and the required ordered `query_batches` manifest with each file's columns and row groups. `query_row_groups` is populated only for one query file; baseline-only output has an empty `query_batches` list. Current readers require the manifest; regenerate older directories.
- Existing canonical LD-score files or `diagnostics/ldscore.log` are refused before writing unless `--overwrite` or `LDScoreOutputConfig(overwrite=True)` is supplied. Successful overwrites remove stale owned single/numbered query files not produced by the new run. Query batches are staged under the output directory and published after successful computation, with metadata last. See `LDScoreDirectoryWriter.artifact_family` and `write_batches` in [outputs.py](../../src/ldsc/outputs.py).
- `ldscore.baseline.parquet` contains `CHR`, `POS`, `SNP`,
  `regression_ld_scores`, and baseline LD-score columns. `regression_ld_scores`
  is the historical `w_ld` LD score over the regression SNP universe, not the
  final h2/rg regression weight. When both baseline and query inputs are
  omitted, the LD-score workflow writes a synthetic all-ones baseline column
  named `base`.
- Query files contain `CHR`, `POS`, `SNP`, available alleles, and the query LD-score columns assigned to that batch. They are omitted for baseline-only runs.
- Regression commands consume this directory via `--ldscore-dir`; users do not
  pass count vectors, weight files, or annotation manifests.

### Reference-panel building

Relevant APIs:

- `run_build_ref_panel(...)`
- `ldsc build-r2-panel`

Accepted path forms:

- `plink_prefix`: exact PLINK prefix, plain chromosome-suite stem, PLINK-prefix glob, or explicit `@` suite token. A plain stem such as `panel_chr` discovers matching trios such as `panel_chr1.{bed,bim,fam}`; BIM contents establish chromosome identity
- genetic-map inputs, when provided: exact files, globs, or `@` chromosome suites; liftover-chain inputs require exact paths
- `ref_panel_snps_file`, when provided: scalar file-like token interpreted
  using `GlobalConfig.snp_identifier`; restriction files may omit alleles and
  then match by base key, while allele-bearing restrictions in allele-aware
  modes match by the effective allele-aware key; duplicate restriction keys collapse
  to one retained key and non-identity columns such as `CM` or `MAF` are ignored;
  `chr_pos`-family coordinates must be aligned to the PLINK source build; the
  resolved restriction is loaded into an in-memory key set before per-chromosome
  PLINK filtering, so very large custom keep-lists can become a memory input

How they are handled:

- `plink_prefix` is resolved at the PLINK prefix level, not at the individual `.bed/.bim/.fam` file level
- a chromosome suite such as `panel_chr@` is expanded one chromosome at a time
- liftover chains are optional; the matching source-to-target chain enables
  cross-build R2 and metadata outputs in `chr_pos`-family modes, while no matching
  chain produces source-build-only outputs; all reference-panel liftover is
  rejected in `rsid`-family modes
- genetic maps are required for every emitted build when `--ld-wind-cm` is set;
  SNP- and kb-window builds may omit maps and write emitted metadata `CM` as
  `NA`
- `snp_identifier` for SNP restrictions comes from `GlobalConfig`; the CLI flag
  constructs a one-invocation identifier config, and the Python wrapper reads
  the registered config
- `build-r2-panel` ignores `GlobalConfig.genome_build`; in `chr_pos`-family modes, the
  restriction file must provide source-build coordinates, either through a
  source-specific column such as `hg19_POS` or through generic `POS` that
  infers to the source PLINK build

Example:

```bash
ldsc build-r2-panel \
  --plink-prefix data/reference/genomes_30x_chr \
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
- For a concrete single-chromosome PLINK prefix such as `panel.1`, the owned
  package is restricted to chromosome 1: `chr1` R2 parquet, `chr1` metadata
  sidecars, `diagnostics/dropped_snps/chr1_dropped.tsv.gz`,
  `diagnostics/metadata.chr1.json`, and `diagnostics/build-r2-panel.chr1.log`.
  Sibling chromosomes are not collisions and are not stale cleanup targets.
- For a `@` chromosome-suite prefix such as `panel.@`, the owned package spans
  every discovered chromosome in the output directory, so full-panel overwrites
  can clean stale target-build, chromosome, dropped-SNP, metadata, and log
  siblings.
- Existing parquet, metadata, dropped-SNP audit, or workflow-log files are
  refused unless `--overwrite` or `ReferencePanelBuildConfig(overwrite=True)`
  is supplied.
- With overwrite enabled, a successful run removes stale owned artifacts inside
  the current package that the run did not produce.
- Dropped-SNP audit sidecars are always written for processed chromosomes
  (header-only when clean) and contain liftover-stage rows with reasons
  `source_duplicate`, `unmapped_liftover`, `cross_chromosome_liftover`, and
  `target_collision`.
- The ref-panel directory may contain many chromosome/build artifacts. Treat it
  as one owned family only for full `@` suite invocations; concrete
  per-chromosome invocations use chromosome-specific ownership so parallel
  array jobs can safely share the directory.

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
- `ldsc munge-sumstats --infer-only` resolves only `--raw-sumstats-file` and
  writes no artifacts; `--output-dir` remains syntactically required for the
  uniform command contract but is neither created nor used

How they are handled:

- h2 and partitioned-h2 sumstats are scalar-style inputs; if a glob expands to
  more than one file, the workflow raises instead of combining them
- rg `--sumstats-sources` is a group-style input; all resolved files are
  deduplicated in first-seen order and then used for pair selection
- `ldscore_dir` is not glob-resolved; it is opened as a directory containing
  `metadata.json` plus parquet payload files
- `ldsc munge-sumstats` uses `--input-format auto` by default after the raw path is
  resolved. The inference layer can detect plain text, including VCF-style
  headers, old DANER, and new DANER before applying the usual column aliases and
  repair suggestions.

Output:

- `ldsc munge-sumstats` writes a self-describing `sumstats.parquet` by default
  (identity in its footer; no `metadata.json`), plus
  `diagnostics/sumstats.log`, and
  `diagnostics/dropped_snps/dropped.tsv.gz` under `output_dir`;
  `--output-format tsv.gz` writes legacy `sumstats.gz`, and
  `--output-format both` writes both curated artifacts. Existing owned
  artifacts for the selected filename family and shared diagnostics are refused unless `--overwrite` or
  `MungeConfig(overwrite=True)` is supplied. With overwrite enabled, a
  successful run removes stale sibling formats for that same stem not produced by the current
  `--output-format`; other stems are preserved.
  `diagnostics/sumstats.log` is not recorded in
  `MungeRunSummary.output_paths`; detailed provenance and output bookkeeping
  are written to the log, row-level liftover drops are written to
  `diagnostics/dropped_snps/dropped.tsv.gz`. The Parquet footer carries the thin identity metadata; no root `metadata.json` is written.
- Packaged HM3 is the default SNP restriction. `--sumstats-snps-file FILE` replaces it; `--no-snp-restriction` disables it. These overrides are mutually exclusive. For different source/output builds, packaged HM3 uses automatic quick liftover, while custom-list or unrestricted runs require a chain file. An explicit chain overrides automatic mapping; matching builds need no liftover.
- `ldsc h2`, `ldsc partitioned-h2`, and `ldsc rg` require `output_dir` and write
  fixed result families. For h2, the written family is `h2.tsv`,
  `diagnostics/ld_score_regression_bins.tsv`, `diagnostics/metadata.json`, and
  workflow-owned `diagnostics/h2.log`. Its default `plots/` and
  `postprocessing/` roots are owned stale derivatives. For rg, that family is
  `rg.tsv`, `rg_full.tsv`,
  `h2_per_trait.tsv`, optional `diagnostics/pairs/`, and workflow-owned
  `diagnostics/rg.log`; rg outputs carry nominal p-values only; existing owned
  artifacts are refused unless `--overwrite` is supplied.
- `ldsc partitioned-h2` requires the LD-score directory to include
  `ldscore.overlap.parquet`. Baseline-only directories run the functional-category
  regime; directories with query columns run the cell-type-specific regime.
- Query-annotation `ldsc partitioned-h2` runs always write a staged
  `diagnostics/query_annotations/` tree under `output_dir`. The retired
  `--write-per-query-results` flag is rejected. The tree contains
  `manifest.tsv` and one folder per query annotation, with per-query
  `partitioned_h2.tsv`, `partitioned_h2_full.tsv`, and `metadata.json`.
  Existing final per-query output is refused unless `--overwrite` is supplied;
  baseline-only runs keep their complete fitted model at the result root and,
  with overwrite enabled, remove a stale `diagnostics/query_annotations/` tree.
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
- target-build SNP restrictions for `ldsc build-r2-panel`; that workflow
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
- For PLINK, verify BIM chromosome contents; filenames do not determine chromosome selection
- Use exact paths or exact-one globs for scalar inputs
- Use one dedicated output directory per reproducible run, and use
  `--overwrite` / `overwrite=True` only when replacing that run is intentional
