# LD Score Calculation

Last updated on: 2026-09-10

Goal: compute LDSC-compatible LD scores from a reference panel alone, from pre-built SNP-level annotations, or from raw BED/gene-list queries plus an explicit baseline.

## Chromosome scope for direct query runs

These rules apply to both PLINK and parquet-R² references, without an additional chromosome-testing flag:

- `@` requires all autosomes 1–22. Missing or invalid required chromosome inputs fail; `@` is not an instruction to use whatever files happen to exist.
- Ordinary globs select their actual matches. Validated contents establish chromosome scope, not filenames: `"annotations/*.22.annot.gz"` may match files containing chromosomes other than chr22, and those chromosomes participate.
- Validated baseline and reference chromosome sets must match exactly. Mixing `@` with a chr22-only group fails.
- Every selected focal/control gene must lie within the shared scope. Selected genes are unique successfully resolved genes after explicit gene-region exclusions, before SNP-support filtering. A chr22-only pathway is fully covered by matching chr21–22 baseline/reference inputs; a selected chr1 gene would fail the entire batch under both `strict` and `resolved-only`. Pathways are never automatically truncated or skipped for incomplete coverage. Query BED regions and prebuilt-query SNPs must also lie within scope.

Quote CLI glob patterns so the package receives them intact. Users own glob selection: a missing file can be undetectable if it disappears from the matches and the remaining required artifacts consistently cover the same subset. Check the chromosomes resolved and entering analysis in `diagnostics/ldscore.log` and `diagnostics/chromosome_scope.json`. See the [pass/fail examples](../docs/current/gene-list-input-format.md#direct-query-chromosome-scope) and [coverage diagnostics](../docs/current/gene-list-diagnostics-and-repair.md#chromosome-scope-and-pathway-coverage). Public indexed workflows still require complete autosomes 1–22.

## Reference inputs and conventions

The examples below assume chromosome-pattern annotation inputs such as `annotations/baseline.1.annot.gz` and a package-built R2 directory such as `r2_ref_panel_1kg30x_1cM_hm3/hg38`. Canonical Parquet files have exactly four columns: `IDX_1`, `IDX_2`, `R2`, and `SIGN`. The endpoint indices reference rows in the required `chr*_meta.tsv.gz` sidecar, which supplies SNP identities, alleles, `MAF`, and `CM`. Schema metadata binds the pair table to that sidecar. External raw R2 and the older identity-expanded format are unsupported; use `build-ref-panel` to construct a canonical panel. See the [format and read pipeline](../docs/current/parquet-r2-format-and-read-pipeline.md).

Input-token rules used below:

- exact path: one concrete file
- glob token: one token that expands to multiple files, for example `"beds/*.bed"`
- explicit chromosome suite: `annotations/baseline.@.annot.gz`

Output directories stay literal; only input fields are expanded.

SNP restriction files used for the reference-panel or regression universes are
identity-only filters. Duplicate restriction keys collapse to one retained key,
and non-identity columns such as `CM` or `MAF` are ignored rather than carried
into LD-score metadata.

The retained reference panel is the LD-score contributor and annotation-count
universe unless `--ref-panel-snps-file` explicitly restricts it. Regression
rows use the bundled HM3 map by default; `--regr-snps-file` replaces that
selection. Named `--regr-snps-exclude-regions` presets are then subtracted only from the
regression/output rows and `w_ld` contributors, not from LD-score contributors,
`M`, `M_5_50`, or overlap counts.

`CM` and `MAF` are population-specific and always come from the **reference
panel**, never the annotation: annotation `CM`/`MAF` are ignored. For the parquet
backend the `chr*_meta.tsv.gz` sidecar is authoritative; for the PLINK backend
`CM` comes from the `.bim` (or an interpolated genetic map via
`--genetic-map-hg19-sources` / `--genetic-map-hg38-sources`) and `MAF` from the
genotypes. `--ld-wind-cm` requires usable reference-panel `CM`; an all-zero /
constant / missing `CM` raises a dedicated error (not bypassable by
`--yes-really`). The two similarly named MAF flags have different effects:

- `--maf-min` filters reference-panel SNPs before LD-score calculation
  (inclusive `MAF >= maf_min`) in both backends. Filtered SNPs do not contribute
  LD scores or annotation counts.
- `--common-maf-min` does **not** filter reference-panel SNPs. It only sets the
  inclusive MAF threshold used for the common-SNP annotation-count vector
  (`M_5_50`); all retained reference SNPs still contribute to LD scores and the
  all-SNP count vector (`M`).

Resolution behavior:

- there is no separate `*_chr` argument anymore; annotation arguments accept exact paths, globs, or explicit `@` suite tokens
- `--plink-prefix` accepts one exact PLINK prefix or a plain chromosome-suite stem; for example, `panel_chr` discovers complete `panel_chr1.{bed,bim,fam}`, `panel_chr2.{bed,bim,fam}`, and so on. Globs and the older `panel_chr@` form remain supported
- group inputs such as `--baseline-annot-sources`, `--query-annot-sources`, `--query-annot-bed-sources`, and `--query-annot-gene-list-sources` may resolve to many files; package-built parquet panels are supplied as one build directory with `--r2-dir`
- direct query preflight validates all matched annotation/reference files before chromosome selection; chromosome-sharded annotation routing and PLINK trio selection use validated contents, not filename hints
- scalar inputs still must resolve to exactly one file

Genome-build behavior for `chr_pos` inputs:

- Python callers can import `infer_chr_pos_build()` and `resolve_chr_pos_table()` from `ldsc`
- workflow calls can set `GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="auto")`
- CLI calls can pass `--snp-identifier chr_pos_allele_aware --genome-build auto`
- auto mode infers hg19 or hg38 and converts 0-based `POS` values to canonical 1-based coordinates when enough reference SNPs are present
- the CLI has no standalone build-inference command; inference is part of existing workflows

Public SNP identifier modes are exactly `rsid`, `rsid_allele_aware`,
`chr_pos`, and `chr_pos_allele_aware`; the default is
`chr_pos_allele_aware`. Mode names are exact. Column aliases apply only to
input headers. Annotation files may omit alleles in allele-aware modes because
they describe genomic membership, but if annotation alleles are present they
participate in allele-aware matching.

Important output behavior:

- the in-memory result is one merged `LDScoreResult` with split `baseline_table` and optional `query_table`
- `--output-dir` writes a canonical LD-score result directory containing root `metadata.json`, `ldscore.baseline.parquet`, optional `ldscore.query.parquet`, and `diagnostics/ldscore.log`
- `ldscore.baseline.parquet` and `ldscore.query.parquet` are still single flat files, but each parquet row group contains exactly one chromosome
- root `metadata.json` records `row_group_layout`, `baseline_row_groups`, and `query_row_groups` for readers that want to load one chromosome by row-group index
- `LDScoreResult.output_paths` lists scientific data artifacts only; it does not include `diagnostics/ldscore.log`
- regression-universe LD scores live in the `regression_ld_scores` column of `ldscore.baseline.parquet`; there is no separate `.w.l2.ldscore.gz` output
- annotation counts are stored as metadata records, not as separate `.M` files
- if both baseline and query inputs are omitted, the workflow synthesizes an all-ones baseline column named exactly `base` over retained reference-panel metadata
- prebuilt, BED, and gene-list queries are mutually exclusive and require explicit baseline annotations; create an explicit all-ones `base` baseline yourself if you intentionally want that query universe
- evaluated BED/gene runs write `diagnostics/query_annotation_status.tsv`; gene runs also write the row-complete `gene_list_audit.tsv.gz` and per-source `gene_list_resolution_summary.tsv`
- missing output directories are created and existing directories are reused
- existing owned LD-score artifacts, including unselected siblings such as a
  stale `ldscore.query.parquet`, fail before writing starts; reruns that should
  replace them must pass `--overwrite` or `overwrite=True`
- with overwrite enabled, successful baseline-only runs remove stale
  `ldscore.query.parquet` so the result directory reflects the current
  metadata

Performance behavior for canonical parquet R2 input:

- `snp_batch_size` / `--snp-batch-size` controls the number of SNPs processed per LD-score sliding batch; the default is `128`
- the canonical parquet reader automatically sizes a decoded row-group cache once per chromosome from the actual LD window and `snp_batch_size`
- cache entries are decoded row groups with numeric endpoint indices and R2 values, not pandas DataFrames or SNP strings
- the cache only avoids rereading overlapping parquet row groups; every query still computes its required row-group set and filters to the current SNP window, so cache state does not affect correctness

## Case 1: Ordinary Unpartitioned LD Scores

For unpartitioned heritability, no baseline annotation input is needed. The
result directory still has the normal `ldscore.baseline.parquet` and manifest
metadata layout, with `baseline_columns == ["base"]`.

### Python API

```python
from ldsc import GlobalConfig, run_ldscore, set_global_config

set_global_config(
    GlobalConfig(
        snp_identifier="chr_pos_allele_aware",
        genome_build="hg38",
    )
)

result = run_ldscore(
    output_dir="tutorial_outputs/unpartitioned_ldscores",
    r2_dir="r2_ref_panel_1kg30x_1cM_hm3/hg38",
    common_maf_min=0.05,
    ld_wind_cm=1.0,
    # snp_batch_size=128,  # optional; also controls parquet cache sizing
)

print(result.baseline_columns)
print(result.baseline_table.loc[:, ["CHR", "SNP", "POS", "regression_ld_scores", "base"]].head())
```

### CLI

```bash
ldsc ldscore \
  --output-dir tutorial_outputs/unpartitioned_ldscores \
  --r2-dir "r2_ref_panel_1kg30x_1cM_hm3/hg38" \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg38 \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
# Add --snp-batch-size 128 only when intentionally tuning the sliding batch size.
```

## Case 2: Existing SNP-Level Annotation Files

### Python API

```python
from ldsc import GlobalConfig, run_ldscore, set_global_config

set_global_config(
    GlobalConfig(
        snp_identifier="chr_pos_allele_aware",
        genome_build="hg38",
    )
)

result = run_ldscore(
    output_dir="tutorial_outputs/r2_ldscores",
    baseline_annot_sources="annotations/baseline.@.annot.gz",
    r2_dir="r2_ref_panel_1kg30x_1cM_hm3/hg38",
    common_maf_min=0.05,
    ld_wind_cm=1.0,
    # overwrite=True,  # also removes stale LD-score siblings not produced by this run
)

print(result.baseline_table.columns.tolist())
print(result.baseline_table.head())
print(result.output_paths["baseline"])
print(result.output_paths["metadata"])
print(result.config_snapshot)
```

### CLI

```bash
ldsc ldscore \
  --output-dir tutorial_outputs/r2_ldscores \
  --baseline-annot-sources "annotations/baseline.@.annot.gz" \
  --r2-dir "r2_ref_panel_1kg30x_1cM_hm3/hg38" \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg38 \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
# Add --overwrite only when intentionally replacing the LD-score artifact family.
```

## Case 3: Use BED Files Directly During LD-Score Calculation

`ldsc ldscore` can now project BED files in memory through `--query-annot-bed-sources` without first materializing intermediate query `.annot.gz` files. Because BED inputs are query annotations, this mode requires explicit baseline annotations.

### Python API

```python
from ldsc import GlobalConfig, run_ldscore, set_global_config

set_global_config(
    GlobalConfig(
        snp_identifier="chr_pos_allele_aware",
        genome_build="hg38",
    )
)

result = run_ldscore(
    output_dir="tutorial_outputs/r2_ldscores_with_queries",
    baseline_annot_sources="annotations/baseline_chr/baseline.@.annot.gz",
    query_annot_bed_sources="beds/*.bed",
    r2_dir="r2_ref_panel_1kg30x_1cM_hm3/hg38",
    common_maf_min=0.05,
    ld_wind_cm=1.0,
)

print(result.query_columns)
print(result.baseline_table.head())
```

### CLI

```bash
ldsc ldscore \
  --output-dir tutorial_outputs/r2_ldscores_with_queries \
  --baseline-annot-sources "annotations/baseline_chr/baseline.@.annot.gz" \
  --query-annot-bed-sources "beds/*.bed" \
  --r2-dir "r2_ref_panel_1kg30x_1cM_hm3/hg38" \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg38 \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
```

## Case 4: Use Gene Lists Directly

Gene-list inputs are one-column, headerless plain or gzip files containing exact
authoritative gene IDs, exact case-sensitive gene names, or a mixture. The source
basename becomes the query name: `immune_genes.txt.gz` becomes `immune_genes`.

```bash
ldsc ldscore \
  --output-dir tutorial_outputs/gene_list_ldscores \
  --baseline-annot-sources "annotations/baseline_chr/baseline.@.annot.gz" \
  --query-annot-gene-list-sources "gene_lists/*.txt.gz" \
  --gene-coordinate-file "annotations/gene-coordinates.hg38.tsv.gz" \
  --padding-bp 100000 \
  --gene-list-resolution-policy strict \
  --r2-dir "r2_ref_panel_1kg30x_1cM_hm3/hg38" \
  --snp-identifier chr_pos_allele_aware \
  --genome-build auto \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
```

The required one-based coordinate catalog is the sole focal/control gene
universe. Its declared build must agree with baseline/reference-panel evidence;
`auto` reports the resolved build in `diagnostics/ldscore.log`. Gene-list
padding must be chosen explicitly, including `--padding-bp 0` for gene bodies.

Strict mode batches all focal/control identifier problems and stops before
LD-score work. Use `resolved-only` deliberately for exploratory subset analysis.
Start diagnosis with `gene_list_resolution_summary.tsv`, then filter
`gene_list_audit.tsv.gz`; after Gate A, `query_annotation_status.tsv` records
query-local viability. Valid siblings continue through SNP-support or variance
skips; when every focal query is skipped, only diagnostics are written.

See [Gene-List Query Input](../docs/current/gene-list-input-format.md) for exact
resolution, coordinate, naming, and partial-success rules.
For repair procedures and every diagnostic term, see
[Gene-list diagnostics and repair](../docs/current/gene-list-diagnostics-and-repair.md).

## Optional: Read One Chromosome From A Result Directory

Full-file readers such as `pd.read_parquet("ldscore.baseline.parquet")` still work.
For large outputs, use the root metadata row-group records to read exactly one
chromosome from `ldscore.baseline.parquet` or `ldscore.query.parquet`.

```python
import json
from pathlib import Path

import pyarrow.parquet as pq

ldscore_dir = Path("tutorial_outputs/r2_ldscores")
metadata = json.loads((ldscore_dir / "metadata.json").read_text())

if metadata["row_group_layout"] != "one_per_chromosome":
    raise ValueError(f"Unsupported row-group layout: {metadata['row_group_layout']!r}")

baseline_rg_by_chrom = {
    entry["chrom"]: entry["row_group_index"]
    for entry in metadata["baseline_row_groups"]
}

pf = pq.ParquetFile(ldscore_dir / "ldscore.baseline.parquet")
baseline_chr22 = pf.read_row_group(baseline_rg_by_chrom["22"]).to_pandas()

print(baseline_chr22["CHR"].unique())
print(baseline_chr22.head())
```

If query annotations were supplied, `metadata["query_row_groups"]` has the same
shape for `ldscore.query.parquet`. It is `None` for baseline-only LD-score
results.

## Case 5: Reuse an Exact Gene LD-Score Index

For many gene sets under one fixed hg19/PLINK configuration, build the
expensive reference calculation once. The builder requires an explicit base
identity mode and explicit hg19 assertion:

```bash
ldsc build-gene-ldscore-index \
  --baseline-annot-sources "annotations/baseline.@.annot.gz" \
  --plink-prefix "reference/1000G.EUR.QC." \
  --output-dir "indexes/baseline_100kb" \
  --gene-coordinate-file "annotations/gene-coordinates.hg19.tsv.gz" \
  --genome-build hg19 \
  --snp-identifier chr_pos \
  --padding-bp 100000 \
  --gene-exclude-regions mhc \
  --regr-snps-exclude-regions mhc-and-centromeres
```

Index construction defaults to `--padding-bp 0`. This example deliberately
requests a non-default 100 kb gene flank and names the index accordingly.

Use exactly `rsid` or `chr_pos`; there is no default, `auto`, inference,
liftover, hg38, or allele-aware index mode. `rsid` joins on `SNP`; `chr_pos`
joins on normalized positive 1-based `(CHR, POS)` and treats baseline SNP labels
as passive. PLINK publishes `CHR`, `POS`, `SNP`, `A1`, and `A2` after either
match. Advanced callers must ensure the baseline, PLINK, restriction, and map
coordinates all use hg19.

Mutable baseline or PLINK duplicate effective-key groups are removed in full,
with a warning and rows in `diagnostics/dropped_snps/`; no representative is
selected. Baseline-only and PLINK-only keys are dropped and counted. An empty
intersection fails. Repeated keys in an identity-only regression restriction
collapse because restrictions are sets. To replace bundled HapMap3 regression
candidates, add
`--regr-snps-file custom.snplist`. Region subtraction still follows
`--regr-snps-exclude-regions`.

Then assemble any number of gene-list query columns without the source PLINK or
baseline files:

```bash
ldsc ldscore \
  --gene-ldscore-index-dir "indexes/baseline_100kb" \
  --query-annot-gene-list-sources "gene_lists/*.txt" \
  --output-dir "tutorial_outputs/indexed_gene_ldscores"
```

Indexed assembly inherits the validated index identity and hg19 provenance.
It adds no gene control unless one existing one-column file is supplied with
`--control-gene-list-file`.
Do not pass live `--snp-identifier` or `--genome-build` options; either option
is rejected even when it equals the index. The resulting canonical directory is
self-contained and can be consumed by `h2`, `rg`, and `partitioned-h2` without
the source index.

One directory is one complete index. It has a single `index_id` and cannot be
extended with chromosomes or profiles. Rebuilding any input or setting requires
a complete replacement with `--overwrite`. The live build log is
`indexes/.baseline_100kb.build-state/build-gene-ldscore-index.log`. On success
the closed log moves to
`indexes/baseline_100kb/diagnostics/build-gene-ldscore-index.log`; failed logs
remain in hidden state and are archived on retry. This keeps the destination
absent or empty until a complete index is published.

Use an updated reader for newly built gene indexes: chromosome metadata now carries only the published-row fingerprint, which also covers the effective SNP key. Otherwise-valid existing indexes remain readable; the semantic `index_id` and six input fingerprints are unchanged. See the [artifact and identity contract](../docs/current/gene-ldscore-index.md#artifact-and-identity-contract).

## Optional: Materialize BED or Gene Projections for Reuse

If you want reusable query `.annot.gz` shards on disk, call `run_annotate(...)` or `ldsc annotate` explicitly.

### Python API

```python
from ldsc import GlobalConfig, run_annotate, set_global_config

set_global_config(GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"))

with run_annotate(
    query_annot_bed_sources="beds/*.bed",
    baseline_annot_sources="annotations/baseline_chr/baseline.@.annot.gz",
    output_dir="annotations/query_from_beds",
    # overwrite=True,  # also removes stale query shards outside the current chromosome set
) as bundle:
    print(bundle.query_columns)
    print(bundle.chromosomes)
```

### CLI

```bash
ldsc annotate \
  --query-annot-bed-sources "beds/*.bed" \
  --baseline-annot-sources "annotations/baseline_chr/baseline.@.annot.gz" \
  --output-dir annotations/query_from_beds
```

The generated query shards are named `query.<chrom>.annot.gz`, so downstream inputs should use a token such as `annotations/query_from_beds/query.@.annot.gz`.

To construct the same persistent annotation format from gene lists, select the gene route instead of BED sources and choose the padding explicitly:

```bash
ldsc annotate \
  --query-annot-gene-list-sources "gene_lists/*.txt" \
  --gene-coordinate-file references/genes_hg19.tsv.gz \
  --gene-list-resolution-policy strict \
  --gene-exclude-regions none \
  --padding-bp 0 \
  --genome-build hg19 \
  --baseline-annot-sources "annotations/baseline_chr/baseline.@.annot.gz" \
  --output-dir annotations/query_from_genes
```

Standalone gene annotation measures support on the cleaned baseline SNP grid and reports `annotation_snp_count`; it does not evaluate reference-panel or regression support. Empty and globally unsupported focal lists are skipped, usable siblings continue, and an all-skipped batch fails. All-one queries are retained. No control-gene route is accepted. See the [standalone guide](../docs/wiki/utility-functionalities/annotate.md) for coverage, diagnostics, and build rules.

The returned bundle references persistent query shards and the original baseline sources. Construction scratch is released before return; an explicit chromosome read can prepare private data again under the output directory. Close the bundle after use. Closure preserves the saved queries and original baselines.

If `--output-dir` does not exist yet, the workflow creates it automatically and logs creation at INFO. If any root-level `query.*.annot.gz` shard already exists, the
command fails before writing any shard, even if that shard is outside the
current chromosome set. Add `--overwrite` only for an intentional rerun; stale
query shards not produced by the successful run are removed.

For script-style annotation runs, `ldsc.annotation_builder.main(argv)` is the
parser entry point and returns the produced `AnnotationBundle`. If you already
have an argparse namespace from the unified `ldsc` parser, call
`ldsc.annotation_builder.run_annotate_from_args(args)`; it dispatches the
workflow directly without converting the namespace back to argv.

## Optional: Inspect `chr_pos` Genome Build Before A Workflow

Use the top-level Python API when you want to preflight a `CHR`/`POS` table in a
notebook or script. `infer_chr_pos_build()` reports the decision only;
`resolve_chr_pos_table()` also returns a normalized table with canonical
chromosome labels and 1-based positions.

```python
import pandas as pd
from ldsc import infer_chr_pos_build, resolve_chr_pos_table

# The table must contain enough HapMap3-overlapping SNPs to support inference.
restriction = pd.read_csv("filters/hapmap3_chr_pos.tsv.gz", sep="\t")

normalized, inference = resolve_chr_pos_table(
    restriction,
    context="tutorial restriction table",
)
decision_only = infer_chr_pos_build(
    normalized.loc[:, ["CHR", "POS"]],
    context="normalized tutorial restriction table",
)

print(inference.genome_build)
print(inference.coordinate_basis)
print(decision_only.summary_message)
print(normalized.head())
```

For command-line runs, use auto mode on the workflow itself:

```bash
ldsc ldscore \
  --output-dir tutorial_outputs/auto_build_ldscores \
  --baseline-annot-sources "annotations/baseline.@.annot.gz" \
  --r2-dir "r2_ref_panel_1kg30x_1cM_hm3/hg38" \
  --snp-identifier chr_pos_allele_aware \
  --genome-build auto \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
```

## Configuration Notes

For Python workflows, `GlobalConfig` now carries only shared runtime settings such as:

- `snp_identifier`
- `genome_build`
- `log_level`

Per-run SNP-universe controls are owned by the workflow-specific configs instead:

- `ref_panel_snps_file` optionally restricts the LD-score reference-panel input and is passed through `run_ldscore(...)` into `RefPanelConfig`; without it, the full retained reference panel remains the contributor universe
- the LD-score workflow intersects each chromosome bundle with `ref_panel.load_metadata(chrom)`, so reference-panel SNP restriction shrinks the sidecar-defined compute-time universe from `B` to `B ∩ A'`; in the no-annotation unpartitioned case, synthetic `B` is the retained reference-panel metadata itself
- `regr_snps_file` replaces the bundled HM3 regression-row default; named region exclusions are subsequently subtracted from those rows and from `w_ld` contributors without changing `B ∩ A'` LD-score contributors or annotation counts

Both explicit restriction files are interpreted only through their active SNP
identity keys. Repeated keys collapse to one retained key, while metadata-like
columns in the restriction table are ignored.

`run_annotate()` projects BED or resolved gene intervals onto the cleaned baseline rows and returns a bundle referencing saved query shards. Any reference-panel restriction is applied later, during LD-score calculation, when the workflow aligns each chromosome shard to the prepared reference-panel metadata.

In-process LD-score results carry a frozen `config_snapshot`. If you later change the registered `GlobalConfig`, existing results keep their original snapshot, and downstream merge points raise `ConfigMismatchError` if you try to combine artifacts produced under incompatible identifier or genome-build assumptions. Package-written LD-score directories whose root metadata is missing current identity provenance are rejected and must be regenerated with the current LDSC package.
