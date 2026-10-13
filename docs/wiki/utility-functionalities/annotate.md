# Create reusable query annotations

Last updated on: 2026-09-15

`ldsc annotate` projects BED intervals or gene lists onto the SNP rows in baseline annotations and writes reusable binary query files. Direct `ldsc ldscore` accepts BED files and gene lists too, so standalone annotation is optional for queries. To prepare baseline categories from BED files, use the [baseline preparation recipe](#build-baseline-annotations-from-bed-files) below.

## Inputs

Supply `--baseline-annot-sources`, `--output-dir`, and exactly one query route:

- `--query-annot-bed-sources`: user-created or externally obtained BED files. The first three fields are chromosome, 0-based start, and exclusive end. See [BED format](../../current/bed-input-format.md). Omitted padding is zero; explicit nonnegative `--padding-bp` expands both ends and clips starts at zero.
- `--query-annot-gene-list-sources`: one-column lists resolved exactly against `--gene-coordinate-file`. Explicit nonnegative `--padding-bp` is required; `0` selects gene bodies. Resolution is strict by default; add `--allow-unresolved-genes` to use the resolved subset. Use `--gene-exclude-regions none` or `mhc`. Standalone annotate has no control-gene-list option.

Baseline rows supply the annotation grid. Exact paths and globs select the actual input contents; an `@` suite explicitly requires autosomes 1–22. Resolved genes outside this scope fail under either resolution policy after explicit exclusions. Catalog and baseline projection builds must agree; rsID identity records projection-build provenance separately.

```bash
ldsc annotate \
  --baseline-annot-sources 'baseline/baseline.@.annot.gz' \
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

The output contains `query.<chrom>.annot.gz`, one file per surviving chromosome with all usable query columns. Rows retain the baseline identities and order after global identity cleanup. The canonical header is `CHR BP SNP CM`, optional allele columns, then query columns; `CM` is explicitly `NA`, and generated values are integer `0`/`1`. Ordinary prebuilt quantitative annotations remain supported by `ldscore` and are normalized to float32.

The query files contain no baseline value columns. Reuse the original baseline sources for LD scoring and pass the generated query suite, for example `--query-annot-sources 'annotations/pathways/query.@.annot.gz'`. Use a glob instead of `@` for an intentionally smaller chromosome scope. Do not split each query column into its own chromosome suite.

Gene Gate A reports safely discoverable source, naming, catalog, and identifier issues. Gate B evaluates support only on baseline rows surviving identity cleanup. Empty and globally unsupported focal queries are skipped while usable siblings continue; an all-skipped batch fails. All-one queries remain valid. A query supported on another chromosome retains its zero-valued shards. Standalone annotate performs no reference-panel or regression-design tests.

`diagnostics/` contains metadata, the command log, statuses, applicable input/scope/catalog diagnostics, a complete drop audit, and gene audits/summaries in gene mode. Gene support is labeled `annotation_snp_count`; reference-panel counts are unevaluated. See the [complete condition–outcome table](../../current/annotate-gene-list-decisions.md#validation-gates-and-conditionoutcome-table).

At the default INFO level, `diagnostics/annotate.log` marks input reading, SNP identity checks and chromosome preparation, and preparation completion with retained counts and elapsed time. Gene projection and output writing follow preparation. Intentional exclusions use one line per gene set with line–gene pairs; CM/MAF notices appear once per annotation file read. Complete audit rows remain available. See [the messages and their meaning](../../current/workflow-logging.md#annotation-preparation).

## Build baseline annotations from BED files

This route is useful when rebuilding binary baseline categories for a new reference panel or assembling custom baseline regions. For example, the [East Asian partitioned-LD-score workflow](https://github.com/mkanai/eas_partitioned_ldscore#datasets) uses the authors' baseline BED files to annotate a different reference SNP set. BED intervals can be projected onto another panel in the same genome build; changing genome builds still requires matching interval coordinates. Existing annotation tables are sufficient when they already cover the intended reference SNPs.

The current command is named `ldsc annotate`, and its interval option is `--query-annot-bed-sources` even when preparing future baseline columns. The generated files are ordinary SNP annotation tables: their role in LD scoring is selected by passing them through `--baseline-annot-sources` or `--query-annot-sources`. The `query.<chrom>.annot.gz` filename does not permanently assign a regression role.

### 1. Define the target SNP rows

`annotate` requires `--baseline-annot-sources`; it does not accept a PLINK prefix as its SNP grid. If an existing annotation suite already has the desired rows and nonconflicting column names, it can supply that grid. To rebuild a baseline on a different panel, create a minimal all-ones grid from that panel's BIM files. Use the reference SNPs that will contribute to LD scores, rather than only the HapMap3 regression SNPs.

The following example assumes a complete hg19 autosomal PLINK panel named `panel/reference.1.bim` through `panel/reference.22.bim`, with matching `.bed` and `.fam` files for later LD scoring. Replace the panel path and genome build with your actual inputs. The Python snippets use only the standard library and create new output directories; use fresh destinations when repeating them.

```python
from pathlib import Path
import gzip

panel = Path("panel")
grid = Path("annotations/snp_grid")
grid.mkdir(parents=True)
for chrom in range(1, 23):
    with (panel / f"reference.{chrom}.bim").open() as source:
        with gzip.open(grid / f"grid.{chrom}.annot.gz", "wt") as output:
            output.write("CHR\tBP\tSNP\tCM\tA1\tA2\tbase\n")
            for line in source:
                chromosome, snp, cm, bp, a1, a2 = line.split()
                output.write(f"{chromosome}\t{bp}\t{snp}\t{cm}\t{a1}\t{a2}\t1\n")
```

The grid contains only `base`, avoiding name collisions with baseline categories being regenerated. Identity validation and cleanup still apply when `annotate` reads it.

### 2. Project the baseline BED categories

Place one genome-wide BED file per category in `baseline_beds/`, with unique filenames such as `Coding_UCSC.bed` and `Coding_UCSC.extend.500.bed`. Each filename stem becomes one column. Do not include a BED named `base.bed`, because the grid already declares that name. BED coordinates must match the panel build and use the [0-based, half-open convention](../../current/bed-input-format.md#standard-bed).

```bash
ldsc annotate \
  --baseline-annot-sources 'annotations/snp_grid/grid.@.annot.gz' \
  --query-annot-bed-sources 'baseline_beds/*.bed' \
  --padding-bp 0 \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg19 \
  --output-dir annotations/projected_baseline
```

**Use final category boundaries in the BED files.** Different categories can already encode different padding widths, so one projection command with `--padding-bp 0` preserves all of them. A filename containing `extend.500` does not cause automatic expansion. Likewise, flanking-only regions and intervals expanded to include their original region are different definitions; preserve the source model's definition. A nonzero `--padding-bp` applies the same additional expansion to every BED in this invocation.

Check `annotations/projected_baseline/diagnostics/query_annotation_status.tsv` and the output headers against the intended category list. Empty, malformed, or globally unsupported inputs can be skipped under the query-source rules; successful completion alone does not establish that every planned baseline category survived. A category supported on another chromosome retains its zero-valued chromosome rows.

### 3. Add the all-SNP baseline column

The output files contain the projected categories but omit the grid's `base` column. For a baseline model that includes an all-SNP category, append `base = 1` to the retained rows and save a new suite:

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

Preparation scans whole-genome and chromosome-sharded files in bounded chunks sized for the sources processed together. It checks SNP identities globally using metadata, then copies retained numeric values directly into private chromosome artifacts. Cross-chromosome duplicate identities are removed before projection. Projection and writing proceed chromosome by chromosome; the returned bundle references persistent outputs and may depend on original baselines. No whole-genome annotation matrix is returned. New private files stay under the selected output directory. Use `with run_annotate(...) as bundle:` or close a Python bundle after reading it.

This architecture supports large pathway batches, including 1,000 pathways subsequently tested separately against shared baseline categories. Standalone annotate is sequential; direct LD scoring additionally supports bounded chromosome workers and query batching. Released arrays become reusable memory, but process RSS may not immediately decrease. See the [developer memory design](../../current/annotation-memory-design.md) and [LD-score guide](../main-functionalities/ldscore.md). Existing artifacts require `--overwrite`; failure markers and stale-output cleanup retain their command contracts.
