# Create reusable SNP annotations

Last updated on: 2026-10-26

`ldsc annotate` projects BED intervals or gene lists onto a BIM or baseline-annotation SNP grid and writes reusable binary query files.
Direct `ldsc ldscore` accepts BED files and gene lists too, so standalone annotation is optional for queries.
To prepare baseline categories from BED files, use the [baseline preparation recipe](#build-baseline-annotations-from-bed-files) below.


Auto build resolution checks every selected annotation or BIM grid shard together with the gene-catalog declaration, when supplied.
It accepts exact paths and globs as well as `@` suites and does not require chromosome 1 merely for inference.
Existing command-specific chromosome coverage requirements still apply.
A concrete build must agree with grid and catalog declarations; unlabeled coordinates use that assertion without inference or conversion.
Ordinary BED is unchecked and inherits the resolved projection build.
See the [shared policy](../main-functionalities/global-config.md#genome-build-inference-checking-and-overrides).

## Inputs

Supply `--output-dir` and exactly one SNP-grid option:

| Grid option | Output |
| --- | --- |
| `--ref-snp-annot-sources` | Supplied baseline values and generated queries, exported as aligned `baseline.<chrom>.annot.gz` and `query.<chrom>.annot.gz`. |
| `--ref-snp-bim-sources` | Generated `query.<chrom>.annot.gz` only; no baseline values are invented. |

BIM mode reads six-column headerless `.bim`/`.bim.gz` files directly, without BED/FAM companions.
Only SNP metadata is used: chromosome, SNP ID, one-based BP, and alleles under the selected identity policy; CM is ignored.
The two grid options are mutually exclusive; BIM plus supplied baseline annotations is unsupported.
The old standalone `--baseline-annot-sources` flag is removed without an alias.
Downstream `ldscore` retains that flag.

Select exactly one query route:

- `--query-annot-bed-sources`: user-created or externally obtained BED files. The first three fields are chromosome, 0-based start, and exclusive end. See [BED format](../../current/bed-input-format.md). Omitted padding is zero; explicit nonnegative `--padding-bp` expands both ends and clips starts at zero.
- `--query-annot-gene-list-sources`: one-column lists resolved exactly against `--gene-coordinate-file`. Explicit nonnegative `--padding-bp` is required; `0` selects gene bodies. Resolution is strict by default; add `--allow-unresolved-genes` to use the resolved subset. Use `--gene-exclude-regions none` or `mhc`. Standalone annotate has no control-gene-list option.

The selected source supplies the annotation grid.
Exact paths and globs select the actual input contents.
BIM `@` suites and annotation-grid gene-list runs require autosomes 1–22.
The existing BED-plus-annotation `@` route retains its available-member behavior; prefer explicit paths or globs for an intentionally partial scope.
Resolved genes outside this scope fail under either resolution policy after explicit exclusions.
Catalog and grid projection builds must agree; rsID identity records projection-build provenance separately.

```bash
ldsc annotate \
  --ref-snp-annot-sources 'baseline/baseline.@.annot.gz' \
  --query-annot-gene-list-sources 'pathways/*.txt' \
  --gene-coordinate-file genes.hg19.tsv.gz \
  --padding-bp 100000 \
  --gene-exclude-regions none \
  --snp-identifier rsid \
  --genome-build hg19 \
  --output-dir annotations/pathways
```

The example uses strict resolution. `--allow-unresolved-genes` takes no value and explicitly permits the audited resolved subset; structural input and coverage checks still apply. The retired `--gene-list-resolution-policy` option is rejected without an alias. Python `run_annotate()` still accepts `gene_list_resolution_policy="resolved-only"`; omission or `None` means strict. See [the resolution contract](../../current/gene-list-input-format.md#strict-and-exploratory-resolution); sources: `add_annotate_arguments()` in [annotation_builder.py](../../../src/ldsc/annotation_builder.py) and `run_annotate()` in [annotate_workflow.py](../../../src/ldsc/annotate_workflow.py).

For BED input, replace the gene-list/catalog/exclusion options with `--query-annot-bed-sources 'beds/*.bed'`; padding remains available.

## Outputs and validation

The output contains `query.<chrom>.annot.gz`, one file per surviving chromosome with all usable query columns.
Rows retain the selected grid identities and order after global identity cleanup.
The canonical header is `CHR BP SNP CM`, optional allele columns, then query columns; `CM` is explicitly `NA`, and generated values are integer `0`/`1`.
Ordinary prebuilt quantitative annotations remain supported by `ldscore` and are normalized to float32.

Annotation-grid mode exports the real baseline columns separately with identical metadata, row membership, row order, and chromosome shards.
Continuous baseline values retain float32 precision; original input files are unchanged.
Use the exported pair together in downstream `ldscore`:

```bash
ldsc ldscore \
  --baseline-annot-sources 'annotations/pathways/baseline.@.annot.gz' \
  --query-annot-sources 'annotations/pathways/query.@.annot.gz' \
  --plink-prefix 'panel/reference.@' \
  --snp-identifier rsid --genome-build hg19 \
  --ld-wind-kb 1000 --output-dir ldscores/pathways
```

Use a glob instead of `@` for an intentionally smaller chromosome scope.
Use the exported baseline files rather than the original baseline inputs: identity cleanup may remove rows, and pairing cleaned queries with the original rows can fail downstream alignment.
BIM mode exports only generated query columns.
Those files can be supplied to either downstream annotation flag according to their intended modeling role; supplying them as queries still requires actual aligned baseline annotations.
No automatic baseline synthesis or alignment with a separate baseline source is performed by BIM-mode `annotate`.
Do not split each query column into its own chromosome suite.

Gene Gate A reports safely discoverable source, naming, catalog, and identifier issues.
After header and catalog validation, naming collisions and unresolved identifiers are reported together before SNP-grid preparation.
Query names come from source filenames: rename a conflicting file and inspect both `gene_list_resolution_summary.tsv` and `gene_list_audit.tsv.gz` under `diagnostics/`.
Naming collisions remain fatal with `--allow-unresolved-genes`, including empty or fully unresolved colliding queries.
Gate B evaluates support only on grid rows surviving identity cleanup.
Empty and globally unsupported focal queries are skipped while usable siblings continue; an all-skipped batch fails.
All-one queries remain valid.
A query supported on another chromosome retains its zero-valued shards.
Standalone annotate performs no reference-panel or regression-design tests.

`diagnostics/` contains metadata, the command log, statuses, applicable input/scope/catalog diagnostics, a complete drop audit, and gene audits/summaries in gene mode.
Metadata `snp_rows` records input, output, and dropped logical SNP rows, in total and per chromosome.
Aligned column sources count the same SNP row once; all removed duplicate occurrences count toward drops.
The log reports totals even when zero rows were dropped and includes exported paths and the downstream baseline/query flags.
Gene support is labeled `annotation_snp_count`; reference-panel counts are unevaluated.
See the [complete condition–outcome table](../../current/annotate-gene-list-decisions.md#validation-gates-and-conditionoutcome-table).

At the default INFO level, `diagnostics/annotate.log` marks input reading, SNP identity checks and chromosome preparation, and preparation completion with retained counts and elapsed time. Gene projection and output writing follow preparation. Intentional exclusions use one line per gene set with line–gene pairs; CM/MAF notices appear once per annotation file read. Complete audit rows remain available. See [the messages and their meaning](../../current/workflow-logging.md#annotation-preparation).

## Build baseline annotations from BED files

This route is useful when rebuilding binary baseline categories for a new reference panel or assembling custom baseline regions. For example, the [East Asian partitioned-LD-score workflow](https://github.com/mkanai/eas_partitioned_ldscore#datasets) uses the authors' baseline BED files to annotate a different reference SNP set. BED intervals can be projected onto another panel in the same genome build; changing genome builds still requires matching interval coordinates. Existing annotation tables are sufficient when they already cover the intended reference SNPs.

The current command is named `ldsc annotate`, and its interval option is `--query-annot-bed-sources` even when preparing future baseline columns. The generated files are ordinary SNP annotation tables: their role in LD scoring is selected by passing them through `--baseline-annot-sources` or `--query-annot-sources`. The `query.<chrom>.annot.gz` filename does not permanently assign a regression role.

### 1. Define the target SNP rows

Supply the panel's BIM files directly with `--ref-snp-bim-sources 'panel/reference.@.bim'`.
Use the SNPs intended to contribute to LD scores, rather than restricting the grid to HapMap3 regression SNPs.
The example assumes a complete hg19 autosomal BIM suite; matching BED/FAM files are needed later for PLINK LD scoring, but not for annotation generation.
Use actual BIM filenames, not a PLINK prefix.
Identity cleanup applies globally before interval projection.

### 2. Project the baseline BED categories

Place one genome-wide BED file per category in `baseline_beds/`, with unique filenames such as `Coding_UCSC.bed` and `Coding_UCSC.extend.500.bed`.
Each filename stem becomes one column.
BED coordinates must match the panel build and use the [0-based, half-open convention](../../current/bed-input-format.md#standard-bed).

```bash
ldsc annotate \
  --ref-snp-bim-sources 'panel/reference.@.bim' \
  --query-annot-bed-sources 'baseline_beds/*.bed' \
  --padding-bp 0 \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg19 \
  --output-dir annotations/projected_baseline
```

**Use final category boundaries in the BED files.** Different categories can already encode different padding widths, so one projection command with `--padding-bp 0` preserves all of them. A filename containing `extend.500` does not cause automatic expansion. Likewise, flanking-only regions and intervals expanded to include their original region are different definitions; preserve the source model's definition. A nonzero `--padding-bp` applies the same additional expansion to every BED in this invocation.

Check `annotations/projected_baseline/diagnostics/query_annotation_status.tsv` and the output headers against the intended category list. Empty, malformed, or globally unsupported inputs can be skipped under the query-source rules; successful completion alone does not establish that every planned baseline category survived. A category supported on another chromosome retains its zero-valued chromosome rows.

### 3. Add an all-SNP column only if required by your model

BIM mode creates only the supplied region annotations and has no implicit `base` column.
For a baseline model that includes an all-SNP category, append `base = 1` to the retained rows and save a new suite:

```python
from pathlib import Path
import gzip

projected = Path("annotations/projected_baseline")
baseline = Path("annotations/custom_baseline")
baseline.mkdir(parents=True)
for chrom in range(1, 23):
    with gzip.open(projected / f"query.{chrom}.annot.gz", "rt") as source:
        with gzip.open(baseline / f"baseline.{chrom}.annot.gz", "wt") as output:
            output.write(source.readline().rstrip("\n") + "\tbase\n")
            for line in source:
                output.write(line.rstrip("\n") + "\t1\n")
```

This preserves the validated output row order and allele columns. Keep the projection diagnostics alongside the derived suite as its construction record. For a deliberately chromosome-limited run, replace the chromosome loops and `@` patterns with the same explicit chromosome scope throughout.

### 4. Use the tables as baseline inputs

```bash
ldsc ldscore \
  --plink-prefix 'panel/reference.@' \
  --baseline-annot-sources 'annotations/custom_baseline/baseline.@.annot.gz' \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg19 \
  --ld-wind-kb 1000 \
  --output-dir results/custom_baseline_ldscores
```

All categories in these tables now belong to the baseline group. A later gene-list or BED query can use its own `--padding-bp` in the `ldscore` invocation; that setting does not alter the prebuilt baseline values.

This recipe produces binary region-membership annotations. It does not reconstruct the complete baseline-LD v2.2 model, which also contains quantitative annotations and MAF-bin indicators. BED score fields and overlap multiplicities are not preserved by this projection route. Adding existing quantitative columns or combining outputs from separate padding groups requires alignment checks and column assembly into one baseline file per chromosome; passing two chromosome suites together is not a column-merge operation.

Implementation sources: `run_annotate()` in [annotate_workflow.py](../../../src/ldsc/annotate_workflow.py), `prepare_bed_queries()` and `build_query_shards()` in [_annotation_queries.py](../../../src/ldsc/_annotation_queries.py), and `_layout()` / `prepare_annotation_sources()` in [_annotation_sources.py](../../../src/ldsc/_annotation_sources.py).

## Memory and Python ownership

Preparation scans whole-genome and chromosome-sharded files in bounded chunks sized for the sources processed together.
It checks SNP identities globally using metadata, then copies retained numeric values directly into private chromosome artifacts.
Cross-chromosome duplicate identities are removed before projection.
Projection and writing proceed chromosome by chromosome; the returned bundle references exported files only, without original-input dependencies.
Binary values use the existing packed storage (eight SNP rows per byte), while continuous values use float32 stores.
Export reads bounded numeric row tiles; public reads return float32, and downstream arithmetic remains float64.
No whole-genome annotation matrix is returned.
New private files stay under the selected output directory.
Use `with run_annotate(...) as bundle:` or close a Python bundle after reading it.

This architecture supports large pathway batches, including 1,000 pathways subsequently tested separately against shared baseline categories. Standalone annotate is sequential; direct LD scoring additionally supports bounded chromosome workers and query batching. Released arrays become reusable memory, but process RSS may not immediately decrease. See the [developer memory design](../../current/annotation-memory-design.md) and [LD-score guide](../main-functionalities/ldscore.md). Existing artifacts require `--overwrite`; failure markers and stale-output cleanup retain their command contracts.

### Python grid selection

```python
from ldsc import GlobalConfig, run_annotate

with run_annotate(
    ref_snp_bim_sources="panel/reference.@.bim",
    query_annot_bed_sources="baseline_beds/*.bed",
    global_config=GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg19"),
    output_dir="annotations/from_bim",
) as bundle:
    print(bundle.query_columns)
    print(bundle.output_paths["query_annotations"])
```

For an existing annotation grid, replace `ref_snp_bim_sources` with `ref_snp_annot_sources`.
That result also exposes `bundle.output_paths["baseline_annotations"]`; pass both path lists to downstream `run_ldscore(baseline_annot_sources=..., query_annot_sources=..., ...)`.
The previous standalone Python keyword `baseline_annot_sources` is removed without an alias.
