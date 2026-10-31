# Calculate LD scores from gene lists

Last updated on: 2026-10-31

Use gene-list mode to turn pathway, expression, proteomic, GO or other gene lists into query annotations for partitioned heritability, estimated with stratified LD score regression (S-LDSC; Finucane et al. 2015).
LDSC3 resolves each gene list against one gene coordinate catalog, marks the SNPs inside the listed genes' intervals, and writes an ordinary, self-contained LD-score directory that `partitioned-h2` reads.

There are two ways to compute these LD scores:

| Mode | Use when | What you supply |
| --- | --- | --- |
| Direct | You run a few gene-list batches, or need settings that change between runs. | Baseline annotations, a reference panel (PLINK or R² matrix), a gene coordinate catalog and the padding. |
| Indexed | You test many gene lists against one fixed baseline, PLINK reference panel and gene coordinate catalog. | A gene LD-score index, built once with `ldsc build-gene-ldscore-index`; it contains the gene coordinate catalog and all scientific settings. |

Both modes give the same kind of LD-score directory, and every general `ldscore` option and output is described in [Calculate LD scores](ld-scores.md).

## Gene lists and coordinate authority

### Gene coordinate catalog

The gene coordinate catalog is a headered TSV or `.tsv.gz` file with these six columns in any order:

```tsv
gene_id	gene_name	chrom	start	end	genome_build
GENE001	GENEA	1	100001	110000	hg19
GENE002	GENEB	2	200001	220000	hg19
```

The coordinates above are fictional.

- `start` and `end` are one-based and inclusive, so these are not BED coordinates.
- `gene_id` values must be unique; `gene_name` may be blank but the column is required.
- `chrom` is an autosome, written `1` or `chr1`.
- Every row declares the same genome build, `hg19`/`GRCh37` or `hg38`/`GRCh38`.
- Extra columns are ignored.

Build the gene coordinate catalog from one authoritative gene annotation release in one genome build: LDSC3 cannot detect rows that silently mix genome builds.
It defines the entire set of genes a gene list can name, and no gene coordinate catalog is packaged with LDSC3.
Direct mode reads it from `--gene-coordinate-file`; a gene LD-score index embeds it when the index is built.
The full contract is in [gene coordinate catalog](../../reference/formats/gene-lists.md#gene-coordinate-catalog).

### Gene lists

Each gene list is a headerless text file, plain or gzip, with one identifier per line; for example, `pathway_A.txt` selecting both genes above:

```text
GENE001
GENE002
```

- Each identifier must exactly match a `gene_id` or an unambiguous `gene_name`; prefer gene IDs, because names can be shared.
- Matching is case-sensitive, Ensembl version suffixes are not removed and synonyms are not inferred.
- Blank lines are ignored and surrounding whitespace is removed; there is no comment syntax, so `#TP53` is read as an identifier, and extra tab-separated fields make the line malformed.
- An identifier repeated in one file, or two names of the same gene, selects that gene once.
- The file name becomes the query annotation name after removing a final `.gz` and then one final `.txt`, `.tsv` or `.list`: `pathway_A.txt.gz` becomes `pathway_A`.
  Query names must be unique and must not equal a baseline column or `gene_control`.

Pass gene lists with `--query-annot-gene-list-sources` as exact paths or quoted `*` patterns such as `"/path/to/gene-lists/*.txt"`; `@` is not expanded.
The full parsing rules are in [gene-list source format](../../reference/formats/gene-lists.md#gene-list-source-format).

## Direct mode

```bash
ldsc ldscore \
  --snp-identifier rsid \
  --genome-build hg19 \
  --baseline-annot-sources "/path/to/baseline/baseline.@.annot.gz" \
  --query-annot-gene-list-sources "/path/to/gene-lists/*.txt" \
  --padding-bp 100000 \
  --gene-coordinate-file /path/to/genes/gene-coordinates.hg19.tsv.gz \
  --gene-exclude-regions mhc \
  --plink-prefix /path/to/plink/1000G.EUR.QC. \
  --ld-wind-cm 1.0 \
  --output-dir results/ldscores/gene-lists
```

- `--padding-bp` is required and must be chosen deliberately: `0` uses the gene bodies, and a positive value adds that many base pairs on both sides of each gene.
- `--gene-exclude-regions mhc` removes genes whose unpadded interval overlaps the MHC, before padding; the direct-mode default is `none`.
  This gene exclusion is recorded in the diagnostics and is separate from `--regr-snps-exclude-regions`, which removes regression SNPs.
- `--genome-build` defaults to `auto` for direct gene-list runs, and the gene coordinate catalog's declared genome build must agree with the baseline annotations and reference panel; there is no liftover.
- An R² matrix works as well: replace `--plink-prefix` with `--r2-matrix-dir`.
- Every selected gene must lie on a chromosome covered by both the baseline annotations and the reference panel; a gene list is never truncated to fit, and an uncovered gene fails the whole batch ([chromosome scope](../../reference/formats/gene-lists.md#direct-query-chromosome-scope)).

## Indexed mode

Indexed mode moves the expensive reference-panel calculation for one baseline, PLINK reference panel and gene coordinate catalog into a one-time build.
Later gene-list runs assemble exact LD scores from the stored index without the original PLINK or baseline files.

### Step 1: build a gene LD-score index

```bash
ldsc build-gene-ldscore-index \
  --snp-identifier rsid \
  --genome-build hg19 \
  --baseline-annot-sources "/path/to/baseline/baseline.@.annot.gz" \
  --gene-coordinate-file /path/to/genes/gene-coordinates.hg19.tsv.gz \
  --padding-bp 100000 \
  --gene-exclude-regions mhc \
  --plink-prefix /path/to/plink/1000G.EUR.QC. \
  --ld-wind-cm 1.0 \
  --output-dir /path/to/indexes/baseline_100kb \
  --regr-snps-exclude-regions mhc-and-centromeres
```

The build requires:

- an explicit `--snp-identifier` (any of the four modes) and an explicit `--genome-build hg19` or `hg38`, with no `auto`, inference or liftover;
- exactly one LD window, `--ld-wind-cm`, `--ld-wind-kb` or `--ld-wind-snps`, with no default; LD windows covering a whole chromosome are rejected and there is no `--yes-really`;
- a PLINK reference panel and baseline annotations covering all autosomes 1 to 22; partial builds are not supported;
- a gene coordinate catalog in the asserted genome build; every row is checked before the build starts, and the whole gene coordinate catalog, including MHC-excluded genes, is stored in the index.

The padding and gene exclusion are fixed in the index.
The index defaults are `--padding-bp 0` (gene bodies) and `--gene-exclude-regions mhc`; the example above requests 100 kb padding, which is why its directory is named `baseline_100kb`.
Regression SNPs default to the HM3 map minus `mhc-and-centromeres`.
To use your own regression SNPs, add `--regr-snps-file custom_regression_snps.tsv`, a headered table whose required columns depend on the SNP identifier mode (`SNP`, or `CHR` and `POS` in the asserted genome build); `--regr-snps-exclude-regions` still applies afterwards, and `none` keeps every listed SNP.
The file format is in [build one complete index](../../reference/gene-ldscore-index.md#build-one-complete-index).

How the inputs are combined:

- Baseline and PLINK SNPs are matched on the SNP identity of the chosen mode, as in PLINK-based `ldscore`; SNPs found on only one side are dropped and counted, every SNP of a duplicated SNP identity is dropped and listed in `diagnostics/dropped_snps/`, and an empty match fails.
- PLINK supplies the published `CHR`, `POS`, `SNP`, `A1` and `A2`; conflicting coordinates for one rsID, or incompatible alleles at one position, stop the build with the repair table `<index-dir>/.build/baseline_plink_issues.tsv`.
- The baseline LD scores, counts and regression weights are recomputed from the PLINK genotypes; precomputed `.l2.ldscore.gz`, `.M` and `.M_5_50` files distributed beside the baseline annotations are not used, so everything in the index shares one set of samples, filters and LD window.
- Gene intervals do not restrict the baseline LD scores: every retained reference SNP contributes according to its baseline annotation value, whether or not it lies in a gene.
  Do not prefilter the PLINK or baseline inputs to gene regions unless you intend that narrower reference panel.

The matching rules are in [SNP identity contract](../../reference/gene-ldscore-index.md#snp-identity-contract), the baseline recomputation in [baseline-suite component usage](../../reference/gene-ldscore-index.md#baseline-suite-component-usage), and the gene-region caveat in [baseline reference-SNP caveat](../../reference/gene-ldscore-index.md#baseline-reference-snp-caveat).

Resources:

- Keep `--workers 1` until you have profiled memory on your system, because each worker can multiply peak memory.
- Lower `--atom-batch-size` (default 64) or `--snp-batch-size` (default 512) to reduce memory; neither changes the results.
- `--keep-indivs-file` keeps the listed individuals, one IID per row.
- A cM LD window needs informative `.bim` `CM` values or `--genetic-map-hg19-sources` or `--genetic-map-hg38-sources` for the asserted genome build; kb and SNP-count LD windows need no genetic map.
  SNP-count LD windows count SNPs after SNP identity matching, genotype and MAF filtering.

While the build runs, follow its log in the destination's hidden build directory:

```bash
tail -f /path/to/indexes/baseline_100kb/.build/build-gene-ldscore-index.log
```

The log first reports annotation preparation (`Reading annotation inputs`, then SNP identity checks, then `Annotation preparation complete`) and then each chromosome.
Its final `Finished` or `Failed` line gives the outcome.
On success, the closed log moves to `/path/to/indexes/baseline_100kb/diagnostics/build-gene-ldscore-index.log` and `.build/` is removed; the build writes nothing outside the destination.
Input checks run before any computation, and failed checks write their repair tables, such as `plink_input_issues.tsv`, under `<index-dir>/.build/`.
If the gene coordinate catalog fails validation, repair it using `<index-dir>/.build/gene_coordinate_catalog_issues.tsv.gz` ([catalog repair](../../reference/results/gene-list-diagnostics.md#gene-coordinate-catalog-repair-and-index-builder-issues)).

One output directory holds one complete, immutable index with a single `index_id`; it cannot be extended with chromosomes or settings, so any change to an input or setting needs a complete rebuild.
The destination stays absent or empty until the whole index is published, and a failed build never leaves a partial index.
Replacing an existing valid index requires `--overwrite`, and the old index stays usable until its replacement has been validated.
A build never resumes from an earlier failed or interrupted attempt.
Indexes in earlier formats (before format version 4) must be rebuilt.
The directory layout, `index_id` and replacement rules are in [artifact and identity contract](../../reference/gene-ldscore-index.md#artifact-and-identity-contract) and [output preflight, logging and replacement](../../reference/gene-ldscore-index.md#output-preflight-logging-and-replacement); build errors are covered in [troubleshooting](../troubleshooting.md#build-gene-ldscore-index).

### Step 2: assemble gene-list LD scores

```bash
ldsc ldscore \
  --query-annot-gene-list-sources "/path/to/gene-lists/*.txt" \
  --output-dir results/ldscores/gene-lists \
  --gene-ldscore-index-dir /path/to/indexes/baseline_100kb
```

The index supplies the gene coordinate catalog, baseline, reference panel, genome build, SNP identifier mode, LD window, padding, gene exclusion and regression SNPs.
Do not pass any of them: options such as `--snp-identifier`, `--genome-build`, `--padding-bp` or `--baseline-annot-sources` are rejected in indexed mode, even when they equal the stored value.
Indexed mode accepts `--control-gene-list-file`, `--allow-unresolved-genes`, `--workers`, `--query-batch-size`, `--overwrite` and `--log-level`; the runtime options do not change the stored scientific settings.
A missing or damaged index stops the run; LDSC3 never falls back to direct computation.
The output records the inherited SNP identifier mode, genome build and `index_id`, and the regression commands use it without the index.

For indexed scoring of a workload like the one measured on the Longleaf cluster (1,000 gene lists across 22 chromosomes), start with `--workers 8` and a Slurm allocation of `--cpus-per-task=8 --mem=128G`, a provisional starting point rather than an established optimum, while the package default stays at one worker ([performance measurements](../../dev/validation/performance.md)).

## Optional control gene list

When your gene lists were selected from a specific set of genes, such as the genes measured by an assay, add that set as one control gene list:

```bash
ldsc ldscore \
  --query-annot-gene-list-sources "/path/to/gene-lists/*.txt" \
  --control-gene-list-file /path/to/gene-lists/control-genes.txt \
  --output-dir results/ldscores/gene-lists-with-control \
  --gene-ldscore-index-dir /path/to/indexes/baseline_100kb
```

The control gene list uses the gene-list file format but must be one exact file, not a pattern.
It becomes the baseline annotation `gene_control`, using the same gene coordinate catalog, gene exclusion and padding as the query annotations, so each query result is conditional on it.
A requested control gene list is never dropped silently: no resolved genes, no annotation SNPs or LD scores without variation stop the run.
Without this option, no control is added, in both modes.

## Strict versus exploratory resolution

Strict resolution is the default and suits final analyses.
LDSC3 checks every gene list and the control gene list together, and any unresolved or ambiguous identifier stops the run before the LD-score calculation, with the whole batch reported.

For a preliminary screen of hundreds of pathways, add `--allow-unresolved-genes`, which takes no value.
The run then uses only the resolved genes, records the resolved and total counts, and skips unusable query annotations without failing the usable ones.
Because this changes the submitted gene lists, the console, the result metadata, the audit and the summary all record it, and the command prints a short notice when genes have no SNP support or query annotations were skipped.
The option applies to direct and indexed runs; it does not relax malformed inputs, index integrity, chromosome coverage or the checks on a control gene list.
The retired option `--gene-list-resolution-policy` is rejected: remove `--gene-list-resolution-policy strict`, and replace `--gene-list-resolution-policy resolved-only` with `--allow-unresolved-genes`.
See [strict and exploratory resolution](../../reference/formats/gene-lists.md#strict-and-exploratory-resolution).

## Diagnose and curate

1. Open `diagnostics/gene_list_resolution_summary.tsv` to find the affected gene lists, including the control gene list.
2. Filter `diagnostics/gene_list_audit.tsv.gz` to the rejected rows of that gene list.
3. Fix your gene list first: replace ambiguous names with gene IDs, correct or remove unmatched values, and remove malformed or duplicate rows.
4. Rerun with `--overwrite` and compare the audit.
5. Change or rebuild the gene coordinate catalog only when its line-level evidence shows a defect in how it was produced.

After the identifiers are resolved, `diagnostics/query_annotation_status.tsv` records query annotations skipped for having no SNPs or no LD-score variation.
Genes without SNP support and genes removed by the MHC exclusion are recorded in the audit but are not unresolved identifiers.
`diagnostics/ldscore.log` summarizes intentional MHC exclusions on one INFO line per gene list and role, and direct runs also log the annotation-preparation steps ([annotation preparation](../../reference/logs-and-outputs.md#annotation-preparation)); indexed runs read the prepared index instead.
Every column, reason and repair step is in [gene-list diagnostics and repair](../../reference/results/gene-list-diagnostics.md), and common failures are in [troubleshooting](../troubleshooting.md#ldscore-gene-list-preflight-stopped-or-a-focal-query-is-missing).

## Continue to partitioned heritability

```bash
ldsc partitioned-h2 \
  --sumstats-file /path/to/trait/sumstats.parquet \
  --ldscore-dir results/ldscores/gene-lists \
  --output-dir results/partitioned-h2/trait
```

Each query result is conditional on the baseline annotations and on `gene_control` when present.
A query annotation missing from the output is a workflow status, not evidence of a biological null; check `diagnostics/query_annotation_status.tsv` before the regression.
Regression failures are handled separately from gene resolution: `--continue-on-query-error` in `partitioned-h2` skips query annotations whose regression fails, whereas `--allow-unresolved-genes` in `ldscore` concerns identifiers.
The regression options, result files and plots are covered in [Partitioned heritability](partitioned-heritability.md), [query failures and continuation](../../reference/results/partitioned-h2.md#query-failures-and-continuation) and [plotting](plotting.md).

## Many gene lists

For 1,000 pathways, `partitioned-h2` tests each pathway separately against the same baseline annotations.
`--query-batch-size` (default 100) sets how many query annotations `ldscore` computes together, in both modes.
A smaller value lowers memory, but direct mode then repeats the reference-panel pass and indexed mode rereads the stored chromosome data for each query batch.
More workers hold more chromosome data in memory at once ([parallel execution](../../reference/parallel-execution.md)).
The baseline and overlap files are shared, and each query batch writes its own query file; the batch widths of `ldscore` and `partitioned-h2` may differ.
Memory and query batches are covered further in [Calculate LD scores](ld-scores.md#many-query-annotations-and-workers).

For the Python API, see [Python usage](../python-usage.md#gene-list-ld-scores).
