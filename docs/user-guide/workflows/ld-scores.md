# Calculate LD scores

Last updated on: 2026-10-28

`ldsc ldscore` computes LD scores from a reference panel and writes them to an LD-score directory.
The regression commands `h2`, `rg` and `partitioned-h2` read that directory through `--ldscore-dir` and never recompute LD scores, so one LD-score directory serves every trait you analyze with the same reference panel and annotations.
Without annotations, the command computes the single all-SNP LD score used for SNP heritability and genetic correlation.
With baseline and query annotations, it computes the per-annotation LD scores used for partitioned heritability and cell-type enrichment.
For gene lists, follow [Calculate LD scores from gene lists](gene-list-ld-scores.md); to write reusable annotation files first, follow [Create reusable SNP annotations](annotations.md).
A worked example from summary statistics to regression is in the [guided tutorial](../guided-tutorial.md).

## Inputs you need

| Input | Options | Notes |
| --- | --- | --- |
| Reference panel | exactly one of `--plink-prefix` or `--r2-matrix-dir` | PLINK `.bed/.bim/.fam` genotypes, or an R² matrix written by `build-r2-matrix` (see [Build and query R² matrices](r2-matrices.md)). |
| LD window | exactly one of `--ld-wind-cm`, `--ld-wind-kb` or `--ld-wind-snps` | No default. |
| SNP identifier mode and genome build | `--snp-identifier`, `--genome-build` | The default mode is `chr_pos_allele_aware`, which needs `--genome-build` (`auto`, `hg19` or `hg38`). |
| Baseline annotations | `--baseline-annot-sources` | Optional; required when you add query annotations. |
| Query annotations | one of `--query-annot-sources`, `--query-annot-bed-sources` or `--query-annot-gene-list-sources` | Optional. |
| Output directory | `--output-dir` | Required. |

Where to find reference panels and baseline annotations is described in [reference data](../reference-data.md).

## Compute unpartitioned LD scores

This command computes LD scores for SNP heritability and genetic correlation from an R² matrix:

```bash
ldsc ldscore \
  --r2-matrix-dir /path/to/r2-matrix/hg38 \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg38 \
  --ld-wind-cm 1.0 \
  --output-dir results/ldscores/unpartitioned
```

With no baseline and no query annotations, LDSC3 creates one all-ones baseline annotation named `base` over the retained reference SNPs, so the output has the ordinary layout with `base` as its only LD-score column.
An R² matrix stores only the SNP pairs inside the LD window it was built with, so the requested LD window must use the same unit as the R² matrix's LD window and be no wider; this command needs a panel built with an LD window of at least 1 cM.

To compute from PLINK genotypes instead, replace `--r2-matrix-dir` with `--plink-prefix`:

```bash
ldsc ldscore \
  --plink-prefix /path/to/plink/1000G.EUR.QC. \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg19 \
  --ld-wind-cm 1.0 \
  --output-dir results/ldscores/unpartitioned
```

The plain stem `1000G.EUR.QC.` finds the complete per-chromosome trios `1000G.EUR.QC.1` to `1000G.EUR.QC.22`; see [PLINK input](#plink-input) for the other accepted forms.

## Read the outputs

| Path | Content |
| --- | --- |
| `metadata.json` | The LD-score directory manifest: SNP identifier mode, genome build, annotation columns, all-SNP and common-SNP counts per annotation, and the list of query files. Regression requires it. |
| `ldscore.baseline.parquet` | `CHR`, `POS`, `SNP`, `regression_ld_scores` and one LD-score column per baseline annotation. |
| `ldscore.query.parquet` | Query LD scores when the run fits in one query batch; larger runs write `ldscore.query.batch00001.parquet`, `ldscore.query.batch00002.parquet` and so on. Absent without query annotations. |
| `ldscore.overlap.parquet` | The annotation overlap matrix that `partitioned-h2` requires; written when the run has two or more LD-score columns. |
| `diagnostics/` | `ldscore.log`, per-chromosome logs under `chromosome-logs/`, and query status, scope and input-issue tables when they apply. |

The rows of the baseline and query files are the regression SNPs, by default the HM3 map after MHC and centromere exclusion.
Each of these files covers all computed chromosomes, with one Parquet row group per chromosome.
`regression_ld_scores` is the LD score computed over the regression SNPs only, which legacy LDSC wrote as the separate `.w.l2.ldscore` weight file.
Annotation counts are stored in `metadata.json` rather than as `.M` and `.M_5_50` files.
Pass the whole directory to the regression commands with `--ldscore-dir`; you never pass count vectors or weight files separately.
Every metadata field is listed in [artifact metadata](../../reference/formats/artifact-metadata.md#ldscore).

## Add annotations

### Baseline annotations

Pass prebuilt annotation tables (`.annot` or `.annot.gz`) with `--baseline-annot-sources`:

```bash
ldsc ldscore \
  --baseline-annot-sources "/path/to/baseline/baseline.@.annot.gz" \
  --r2-matrix-dir /path/to/r2-matrix/hg38 \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg38 \
  --ld-wind-cm 1.0 \
  --output-dir results/ldscores/baseline
```

The `@` in the path is a chromosome suite: it stands for chromosomes 1 to 22 and requires all 22 files.
Quote paths that contain `@` or `*` so the shell passes them unchanged; the path forms are described in [input paths](../../reference/input-paths.md).
Annotation values may be binary or continuous.
Annotation files may omit both allele columns even in allele-aware modes; when they contain `A1` and `A2`, the alleles take part in matching, and supplying only one of them is an error.
To build baseline annotations from BED files, follow [Build baseline annotations from BED files](annotations.md#build-baseline-annotations-from-bed-files).

### Query annotations

Query annotations are tested one at a time on top of the baseline in `partitioned-h2`.
Every query route requires explicit baseline annotations, and the three routes are mutually exclusive:

| Route | Option | Use when |
| --- | --- | --- |
| BED files | `--query-annot-bed-sources` | You have one BED file of intervals per query annotation; LDSC3 projects them onto SNPs during the run. |
| Gene lists | `--query-annot-gene-list-sources` | You have one gene list per query annotation; see [Calculate LD scores from gene lists](gene-list-ld-scores.md). |
| Prebuilt annotation files | `--query-annot-sources` | You already have query `.annot.gz` files aligned with the baseline, for example from `ldsc annotate`. |

BED queries are projected in memory, without writing intermediate annotation files:

```bash
ldsc ldscore \
  --baseline-annot-sources "/path/to/baseline/baseline.@.annot.gz" \
  --query-annot-bed-sources "beds/*.bed" \
  --r2-matrix-dir /path/to/r2-matrix/hg38 \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg38 \
  --ld-wind-cm 1.0 \
  --output-dir results/ldscores/bed_queries
```

Each BED file becomes one query annotation named after the file.
BED sources accept exact paths or quoted `*` patterns; `@` is not expanded for them.
`--padding-bp` extends every interval at both ends and defaults to 0; keep it at 0 when the intervals already include the padding you want.
BED file rules are in [Prepare BED files](annotations.md#prepare-bed-files) and the [BED input format](../../reference/formats/bed-intervals.md).

For prebuilt query files written by `ldsc annotate`, pass the exported baseline and query files together:

```bash
ldsc ldscore \
  --baseline-annot-sources "annotations/pathways/baseline.@.annot.gz" \
  --query-annot-sources "annotations/pathways/query.@.annot.gz" \
  --plink-prefix /path/to/plink/1000G.EUR.QC. \
  --snp-identifier rsid \
  --ld-wind-cm 1.0 \
  --output-dir results/ldscores/pathways
```

If you want query annotations tested against an all-SNP baseline only, create an explicit baseline file with a `base` column of ones; LDSC3 creates the synthetic `base` only when there are no query annotations.
Step 3 of the [baseline recipe](annotations.md#3-add-an-all-snp-column-only-if-required-by-your-model) shows how.
After a BED or gene-list run, `diagnostics/query_annotation_status.tsv` records which query annotations were kept or skipped.

## Choose which SNPs are used

LDSC3 keeps two SNP sets apart.
The reference SNPs, all retained reference-panel SNPs, contribute to every LD score and to the annotation counts and overlaps.
The regression SNPs are the rows written to the output and used in the regression.
The full definitions are in [reference and regression SNPs](../../reference/ldscore-snp-universes.md).

| Option | Effect |
| --- | --- |
| `--maf-min` | Removes reference SNPs with MAF below the value (keeps MAF at or above it); removed SNPs contribute nothing. Default: no filter. |
| `--common-maf-min` | Sets the MAF threshold, inclusive, of the common-SNP count; it does not remove any SNP. Default: 0.05. |
| `--ref-panel-snps-file` | Restricts the reference SNPs to a headered SNP identity table. |
| `--keep-indivs-file` | PLINK only: keeps the individuals listed, one IID per row. |
| `--regr-snps-file` | Replaces the HM3 map as the source of regression SNPs, using a headered SNP identity table. |
| `--regr-snps-exclude-regions` | Region-exclusion preset removed from the regression SNPs: `none`, `mhc`, `centromeres` or `mhc-and-centromeres` (default). Reference SNPs and counts are unchanged. |

SNP restriction files act only through the SNP identity columns of the current SNP identifier mode: repeated SNPs count once, and other columns such as `CM` or `MAF` are ignored.
`CM` and `MAF` always come from the reference panel, never from annotation files: an R² matrix takes them from its sidecar, and PLINK input takes `CM` from the `.bim` file (or a genetic map) and `MAF` from the genotypes.
The region presets are listed in [region-exclusion presets](../../reference/formats/region-exclusion-presets.md).

## Choose the LD window

Give exactly one of `--ld-wind-cm`, `--ld-wind-kb` or `--ld-wind-snps`; each is the maximum distance on either side of a SNP.

- `--ld-wind-cm` needs usable `CM` values.
  With PLINK input, supply `--genetic-map-hg19-sources` or `--genetic-map-hg38-sources` (comma-separated exact file paths) when the `.bim` `CM` column is uninformative; an all-zero or constant `CM` stops the run, and `--yes-really` does not override that check ([troubleshooting](../troubleshooting.md#ldscore-unusable-cm-for-ld-wind-cm)).
  With an R² matrix, `CM` comes from the sidecar and the genetic-map options do not apply.
- With an R² matrix, the LD window must use the unit the R² matrix was built with and be no wider than the R² matrix's LD window.
- `--ld-wind-snps` counts retained SNPs with PLINK input but sidecar rows with an R² matrix, so the two differ when a run keeps only part of the R² matrix; every R² matrix run with `--ld-wind-snps` writes the per-chromosome counts to `diagnostics/ld_window_snp_counts.tsv`.
- `--yes-really` allows an LD window that covers a whole chromosome, which is otherwise rejected.

LDSC3 sums exactly the SNP pairs inside the requested LD window for both PLINK and R² matrix input, so PLINK LD scores do not depend on `--snp-batch-size` and equal legacy LDSC run with `--chunk-size 1`.
PLINK and R² matrix LD scores still differ slightly, because the R² matrix stores R² at reduced precision and may use a different `CM` source.
The measured difference from published legacy LD scores is in the [legacy LDSC flag map](../../reference/legacy-flag-map.md#exact-plink-ld-windows), and the R² matrix window rules are in the [R² matrix format](../../reference/formats/r2-matrices.md#37-ld-windows-and-sidecar-values-in-ldscore).

## PLINK input

`--plink-prefix` accepts an exact prefix (`panel.22`), a plain stem (`1000G.EUR.QC.` or `panel_chr`), a quoted glob (`"panel.*"`) or a chromosome suite (`"panel.@"`).
Every selected prefix needs a complete `.bed/.bim/.fam` trio, and the chromosome of each trio comes from its `.bim` contents, not its file name ([PLINK prefix resolution](../../reference/input-paths.md#plink-prefix-resolution)).
`--snp-batch-size` (default 512) sets how many SNPs each genotype batch processes; it changes speed and memory but never the LD scores, and larger batches are not faster beyond about 512 to 1024 ([measurements](../../dev/design/memory-and-cpu.md#direct-ld-score-computation)).
It is ignored with `--r2-matrix-dir`.
`--export-ref-metadata` also writes `ref_metadata/chrN_meta.tsv.gz` with the effective `CM` and other reference metadata, without changing the `.bim` files.

## Genome build

The coordinate-family modes `chr_pos` and `chr_pos_allele_aware` need `--genome-build`; rsID-family runs without gene lists need none, and direct gene-list runs default to `auto`.
`--genome-build auto` checks all coordinate inputs together: annotation files, `.bim` files or R² matrix sidecars, the gene coordinate catalog and coordinate-based SNP restriction files must agree.
A concrete `hg19` or `hg38` must agree with any genome build those inputs declare, and it is used for unlabeled coordinates; `ldscore` never converts coordinates.
To check a file's genome build before a run, use `ldsc infer-genome-build`.
See [genome-build inference, checking and overrides](../concepts/snp-identity-and-genome-build.md#genome-build-inference-checking-and-overrides) and the [full policy](../../reference/snp-identity-and-genome-build.md#genome-build-policy).

## Chromosome scope of query runs

For a genome-wide query run, use `@` chromosome suites: `@` requires every autosome from 1 to 22, and a missing or invalid file fails the run.
For a deliberate chromosome subset, use quoted `*` patterns or exact paths that select the same chromosomes for the baseline annotations and the reference panel; file contents, not file names, decide which chromosomes are included, and the two chromosome sets must match exactly.
Every BED interval, prebuilt query SNP and selected gene must lie inside that shared scope.
A glob cannot reveal a file that is missing from all inputs consistently, so check the chromosomes reported in `diagnostics/ldscore.log` and `diagnostics/chromosome_scope.json`.
The rules and examples are in [direct LD-score query scope](../../reference/input-paths.md#direct-ld-score-query-scope) and [gene-list chromosome scope](../../reference/formats/gene-lists.md#direct-query-chromosome-scope); repairs are in [troubleshooting](../troubleshooting.md#ldscore-chromosome-coverage-preflight).

## Many query annotations and workers

`--query-batch-size` (default 100) sets the largest number of query annotations computed together in one pass; with 1,000 query annotations, the default runs ten query batches.
Each query batch is written to its own query file and released before the next starts, so a smaller value lowers memory, while direct computation repeats the reference-panel pass for every query batch.
The final output still grows with the total number of query annotations.

`--workers` sets how many chromosomes are computed at once; the default is one worker, and each worker adds memory for its own chromosome data.
Annotation inputs are prepared once, serially, before the workers start.
How workers are counted and limited is described in [parallel execution](../../reference/parallel-execution.md).
For indexed gene-list scoring, a measured starting point is given in [Calculate LD scores from gene lists](gene-list-ld-scores.md#step-2-assemble-gene-list-ld-scores), based on the [performance measurements](../../dev/validation/performance.md).

A run has no checkpoint or resume: if any chromosome fails, the run fails.
Intermediate files stay under the output directory and are removed after a handled success or failure, and `metadata.json` is written only after every query batch succeeds.
Memory released during a run may not show as a drop in the process's resident memory right away.

## Follow progress and fix input errors

Before preparing annotations, LDSC3 checks every declared input path, annotation header and reference-panel file, and collects all problems it finds; repair every row of `diagnostics/input_issues.tsv` together before rerunning.
At the default INFO level, `diagnostics/ldscore.log` records reading the annotation inputs, checking SNP identities and preparing chromosomes, completion of preparation with retained counts and elapsed time, and then each query batch and chromosome.
Completed preparation does not mean the LD scores are done; computation and writing follow.
Per-chromosome details and tracebacks are in `diagnostics/chromosome-logs/chrN.log`, and `--log-level DEBUG` adds per-source details.
The messages are described in [input validation](../../reference/logs-and-outputs.md#staged-input-validation), [annotation preparation](../../reference/logs-and-outputs.md#annotation-preparation) and [LD logs](../../reference/parallel-execution.md#ld-logs-and-failure-cleanup); error remedies are in [troubleshooting](../troubleshooting.md#ldscore).

## Rerun into an existing directory

If the output directory already holds LD-score files from an earlier run, including query files that the new run would not write, the command stops before writing.
Add `--overwrite` to replace them; a successful overwrite also removes stale query files so the directory matches the new `metadata.json`.
A missing output directory is created, unrelated files are kept, and a failed overwrite leaves a failure marker; see [logs and outputs](../../reference/logs-and-outputs.md#output-family-preflight).

## Caveats

- LD-score directories written by earlier versions of LDSC3 without the `query_batches` manifest or SNP identity provenance in `metadata.json` are rejected; regenerate them.
- To use published legacy LD scores, convert them with `ldsc convert-ldsc2-ldscores` rather than passing the legacy files to regression; see [Migrating from legacy LDSC](migrating-from-legacy-ldsc.md).
- The LD-score computation itself is specified in [LD-score computation from R² matrices](../../reference/methods/ld-score-computation.md).

For the Python API, see [Python usage](../python-usage.md#ld-scores).
