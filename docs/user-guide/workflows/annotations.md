# Create reusable SNP annotations

Last updated on: 2026-10-27

`ldsc annotate` projects BED intervals or gene lists onto a SNP grid and writes reusable annotation files (`.annot.gz`), one query annotation per BED file or gene list.
`ldsc ldscore` reads BED files and gene lists itself, so this step is optional for query annotations; use it when you want annotation files on disk to inspect, share or reuse across LD-score runs.
It is also the way to build binary baseline annotations from BED files, as in the [baseline recipe](#build-baseline-annotations-from-bed-files) below.
The exact contract, including every validation outcome, is in the [`annotate` reference](../../reference/annotate.md).

## Inputs you need

A run takes `--output-dir`, exactly one SNP grid and exactly one query route.

| SNP grid option | Rows come from | Output |
| --- | --- | --- |
| `--ref-snp-annot-sources` | Existing annotation files, such as a baseline annotation suite | Aligned `baseline.<chrom>.annot.gz` (your baseline values) and `query.<chrom>.annot.gz` (the new query annotations). |
| `--ref-snp-bim-sources` | PLINK `.bim` or `.bim.gz` files | `query.<chrom>.annot.gz` only; no baseline values are created. |

The two grid options are mutually exclusive, and a BIM grid cannot be combined with separately supplied baseline annotations.
A BIM grid uses only the chromosome, SNP name, one-based position and, in allele-aware modes, the alleles of each six-column, headerless `.bim` row; the `CM` column is ignored, and no `.bed` or `.fam` files are needed.
Pass actual `.bim` file names, not a PLINK prefix.
`annotate` has no `--baseline-annot-sources` option; `ldscore` keeps that option.

| Query route | Option | Required with it |
| --- | --- | --- |
| BED files | `--query-annot-bed-sources` | Nothing else; `--padding-bp` is optional and defaults to 0. |
| Gene lists | `--query-annot-gene-list-sources` | `--gene-coordinate-file` and an explicit `--padding-bp` (0 for gene bodies). |

For the coordinate-family modes `chr_pos` and `chr_pos_allele_aware` (the default), pass `--genome-build` as `auto`, `hg19` or `hg38`; on the gene-list route an omitted genome build is treated as `auto`.
`auto` checks every selected grid file, and the gene coordinate catalog when present, together; a concrete genome build must agree with the genome builds those inputs declare.
BED files are not checked: they are assumed to use the resolved genome build, and nothing is lifted over.
See [genome build](../../reference/annotate.md#genome-build) and the [shared policy](../concepts/snp-identity-and-genome-build.md#genome-build-inference-checking-and-overrides).

## Prepare BED files

Use one text BED file per query annotation, with at least three tab-separated fields: chromosome, start and end.

```text
chr1	1000	1250
chr1	5000	5400
```

- Intervals are zero-based half-open intervals, so the first line covers the one-based SNP positions 1001 to 1250.
- Extra columns are allowed but ignored: each annotation is binary interval membership, not a weight from a score column.
- Use intervals in the same genome build as the SNP grid.
- Each file's name without its final extension becomes the query annotation name, so names must be unique and must not equal a baseline column; `enhancers.bed.gz` becomes `enhancers.bed`.
- Plain text and gzip are accepted; comment, `track`, `browser` and one leading header line are skipped.
- `--padding-bp` adds the same number of base pairs to both ends of every interval in the run; keep it at 0 when the intervals already contain the padding you want.

Text BED files are different from the binary PLINK `.bed` genotype files.
Pass BED files as exact paths or a quoted `*` pattern such as `"beds/*.bed"`; `@` is not expanded for BED sources.
An unreadable or malformed BED file stops the run, a BED file without intervals is skipped, and a BED file that overlaps no grid SNP is kept as an all-zero column.
The full format is in [BED input format](../../reference/formats/bed-intervals.md).

## Annotate from gene lists

```bash
ldsc annotate \
  --ref-snp-annot-sources "/path/to/baseline/baseline.@.annot.gz" \
  --query-annot-gene-list-sources "pathways/*.txt" \
  --gene-coordinate-file /path/to/genes/gene-coordinates.hg19.tsv.gz \
  --padding-bp 100000 \
  --gene-exclude-regions none \
  --snp-identifier rsid \
  --genome-build hg19 \
  --output-dir annotations/pathways
```

The gene list and gene coordinate catalog formats are described in [Calculate LD scores from gene lists](gene-list-ld-scores.md#gene-lists-and-coordinate-authority).
Resolution is strict: any unresolved or ambiguous identifier stops the run after the diagnostics are written.
Add `--allow-unresolved-genes`, which takes no value, to continue with the resolved genes and record the omissions; input structure, naming and chromosome coverage checks still apply ([strict and exploratory resolution](../../reference/formats/gene-lists.md#strict-and-exploratory-resolution)).
`--gene-exclude-regions` is `none` (default) or `mhc`, which removes genes whose unpadded interval overlaps the MHC.
`annotate` has no control gene list option; a control gene list exists only in `ldscore`, where it becomes the baseline annotation `gene_control`.

## Annotate from BED files

For BED input, replace the gene-list, gene coordinate catalog and exclusion options with `--query-annot-bed-sources`:

```bash
ldsc annotate \
  --ref-snp-annot-sources "/path/to/baseline/baseline.@.annot.gz" \
  --query-annot-bed-sources "beds/*.bed" \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg38 \
  --output-dir annotations/beds
```

## Chromosome scope

Exact paths and globs select the files they match, and the file contents decide which chromosomes are included.
A BIM chromosome suite (`@`) and an annotation-grid suite in a gene-list run require all autosomes 1 to 22.
An annotation-grid suite with BED files uses whichever chromosome files exist; use explicit paths or globs when you intend a partial chromosome set.
In gene-list runs, a selected gene outside the grid's chromosomes fails the run under either resolution policy, and gene lists are never truncated.

## Use the outputs

Each `query.<chrom>.annot.gz` holds all usable query annotations for one chromosome.
Its columns are `CHR`, `BP`, `SNP`, `CM` (always `NA`), the allele columns when present, and then the query annotations as integers 0 and 1.
With an annotation grid, `baseline.<chrom>.annot.gz` holds your baseline values on exactly the same rows, in the same order and chromosome files; continuous baseline values keep float32 precision, and your original files are unchanged.
Rows follow the grid after SNP identity cleanup, which may remove duplicate or invalid SNPs.
`annotate` applies no reference-panel restriction; `ldscore` aligns the files to the reference panel when it computes LD scores.

Pass the exported pair together to `ldscore`:

```bash
ldsc ldscore \
  --baseline-annot-sources "annotations/pathways/baseline.@.annot.gz" \
  --query-annot-sources "annotations/pathways/query.@.annot.gz" \
  --plink-prefix "/path/to/plink/1000G.EUR.QC.@" \
  --snp-identifier rsid \
  --ld-wind-kb 1000 \
  --output-dir results/ldscores/pathways
```

- Use the exported baseline files, not your original baseline inputs: SNP identity cleanup may have removed rows, and pairing the cleaned query files with the original rows can fail alignment in `ldscore`.
- Use a glob instead of `@` for a deliberately smaller chromosome set.
- Query files from a BIM grid can serve as either baseline or query annotations, according to their role in your model; used as query annotations, they still need actual baseline annotations on the same rows, which `annotate` does not create for a BIM grid.
- Pass each family as one chromosome suite; do not split query columns into separate suites.

The full output list is in [`annotate` outputs](../../reference/annotate.md#outputs) and [artifact metadata](../../reference/formats/artifact-metadata.md#annotate).

## Check the results

`diagnostics/query_annotation_status.tsv` records the outcome of every query annotation, and gene-list runs add `diagnostics/gene_list_resolution_summary.tsv` and `diagnostics/gene_list_audit.tsv.gz`.
Gene-list problems in the inputs, the gene coordinate catalog, query names and identifiers are all reported together before the SNP grid is read; a query name that collides with another query, a baseline column or a reserved name fails even with `--allow-unresolved-genes`, so rename the gene list file.
Gene SNP support is then measured on the grid rows that survive SNP identity cleanup and reported as `annotation_snp_count`; `annotate` makes no reference-panel or regression check.
Empty gene lists and gene-list query annotations that overlap no grid SNP are skipped while the others continue, whereas a BED file without overlap is kept as an all-zero column.
A run in which every query annotation is skipped fails, and an annotation that is 1 on every row is kept.
After a successful run with omitted genes, unsupported genes or skipped query annotations, the command prints a short notice naming them.
The complete condition-and-outcome table is in [validation and outcomes](../../reference/annotate.md#validation-and-outcomes).

`diagnostics/metadata.json` records, in `snp_rows`, the input, output and dropped SNP rows, in total and per chromosome; aligned columns count a SNP row once.
`diagnostics/annotate.log` reports the same totals, the exported paths and the `ldscore` options to use them, and at the default INFO level marks reading the inputs, checking SNP identities and preparing chromosomes, and completing preparation; projection and writing follow ([annotation preparation](../../reference/logs-and-outputs.md#annotation-preparation)).
Error remedies are in [troubleshooting](../troubleshooting.md#annotate).

## Build baseline annotations from BED files

This recipe rebuilds binary baseline annotations for a new reference panel or assembles custom baseline regions, for example by projecting another study's published baseline BED files onto your reference panel's SNPs.
BED intervals can be projected onto any reference panel in the same genome build; a different genome build needs intervals in that genome build.
If existing annotation files already cover your reference SNPs, you do not need this recipe.

The interval option is called `--query-annot-bed-sources` even when you prepare future baseline annotations.
The output files are ordinary annotation files; their role is chosen later, by passing them to `ldscore` through `--baseline-annot-sources` or `--query-annot-sources`, and the `query.<chrom>.annot.gz` file name does not fix that role.

### 1. Define the target SNPs

Use the reference panel's `.bim` files as the grid: `--ref-snp-bim-sources "panel/reference.@.bim"`.
Include all SNPs meant to contribute to LD scores, not only the HM3 regression SNPs.
The example assumes a complete hg19 autosomal `.bim` suite; the matching `.bed` and `.fam` files are needed later for PLINK LD scores but not for annotation.
SNP identity cleanup is applied across the whole grid before projection.

### 2. Project the baseline BED files

Put one genome-wide BED file per baseline annotation in `baseline_beds/`, with unique names such as `Coding_UCSC.bed` and `Coding_UCSC.extend.500.bed`; each name without its extension becomes one column.
The BED coordinates must use the reference panel's genome build.

```bash
ldsc annotate \
  --ref-snp-bim-sources "panel/reference.@.bim" \
  --query-annot-bed-sources "baseline_beds/*.bed" \
  --padding-bp 0 \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg19 \
  --output-dir annotations/projected_baseline
```

Write the final interval boundaries into the BED files.
Different annotations can already encode different padding, so one run with `--padding-bp 0` keeps all of them; a name containing `extend.500` does not expand anything by itself.
Flanking-only regions and regions expanded to include the original interval are different definitions, so keep the source model's definition.
A nonzero `--padding-bp` adds the same padding to every BED file in the run.

Compare `annotations/projected_baseline/diagnostics/query_annotation_status.tsv` and the output header with your intended list of annotations.
Empty BED files are skipped and malformed ones stop the run, but a BED file that overlaps no grid SNP is kept as an all-zero column, so a successful run alone does not prove that every planned annotation is usable ([BED queries](../../reference/annotate.md#bed-queries)).

### 3. Add an all-SNP column only if required by your model

A BIM grid creates only the annotations in your BED files and has no `base` column.
If your baseline model includes an all-SNP annotation, append a `base` column of ones and save a new suite:

```bash
mkdir -p annotations/custom_baseline
for chrom in $(seq 1 22); do
  gzip -dc "annotations/projected_baseline/query.${chrom}.annot.gz" \
    | awk 'BEGIN { OFS = "\t" } NR == 1 { print $0, "base"; next } { print $0, 1 }' \
    | gzip > "annotations/custom_baseline/baseline.${chrom}.annot.gz"
done
```

This keeps the validated row order and allele columns.
Keep the projection diagnostics with the new suite as the record of how it was built.
For a deliberately chromosome-limited run, use the same chromosome set in the loop and replace `@` with matching explicit paths throughout.

### 4. Use the files as baseline annotations

```bash
ldsc ldscore \
  --plink-prefix "panel/reference.@" \
  --baseline-annot-sources "annotations/custom_baseline/baseline.@.annot.gz" \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg19 \
  --ld-wind-kb 1000 \
  --output-dir results/ldscores/custom_baseline
```

Every annotation in these files is now a baseline annotation.
A later BED or gene-list query run can use its own `--padding-bp`, which does not change these baseline values.

This recipe produces binary region-membership annotations only.
It does not rebuild the complete baseline-LD v2.2 model, which also contains continuous annotations and MAF-bin indicators, and BED score fields and overlapping-interval counts are not kept.
Adding existing continuous columns, or combining outputs made with different padding, requires aligning the rows and assembling the columns into one baseline file per chromosome; passing two chromosome suites together does not merge their columns.

## Rerun and resources

If the output directory already holds any `query.*.annot.gz` or `baseline.*.annot.gz` file, even for a chromosome outside the current run, `annotate` stops before writing.
Add `--overwrite` only for an intended rerun; a successful overwrite removes stale annotation files that the new run did not produce, and a failed one leaves a failure marker ([logs and outputs](../../reference/logs-and-outputs.md#output-family-preflight)).
A missing output directory is created.

`annotate` runs serially; it reads inputs in bounded chunks and writes chromosome by chromosome, keeping its intermediate files under the output directory.
Released memory may not show as a drop in the process's resident memory right away.
Direct `ldscore` additionally offers chromosome workers and query batches for large pathway collections ([Calculate LD scores](ld-scores.md#many-query-annotations-and-workers)).

For the Python API, see [Python usage](../python-usage.md#annotations).
