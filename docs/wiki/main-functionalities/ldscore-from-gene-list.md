# Calculate LD scores for gene lists with an index

Last updated on: 2026-08-11

## Motivation

Gene expression, proteomic, GO, and SynGO analyses often produce many related
gene sets. The scientific question is whether common-variant heritability is
enriched near a focal set after accounting for ordinary functional annotations.

The original LDSC-SEG analysis ([Finucane et al., 2018](https://doi.org/10.1038/s41588-018-0081-4)) fitted a specifically expressed gene annotation with the baseline categories and a control annotation constructed from all genes in the analyzed gene-expression data set. This tutorial does not add a control-gene annotation, so each focal annotation is tested against the baseline categories only.

An exact index makes the LD-score stage fast to repeat across many focal gene lists; the downstream `partitioned-h2` model remains S-LDSC.

## Goal

Resolve one or more gene lists against the index's embedded catalog, assemble
their exact LD scores, and write the normal self-contained partitioned LD-score
directory used by `ldsc partitioned-h2`.

The full mathematical algorithm—from catalog resolution and Boolean atom union
through $Y Z$ assembly, counts, overlaps, chromosome aggregation, and output
artifacts—is documented in
[Exact gene LD-score index: mathematical algorithm](../../current/gene-ldscore-index-mathematics.md).

## Inputs

You need:

- one explicit complete index built with `ldsc build-gene-ldscore-index`;
- one or more plain-text or gzip-compressed one-column gene lists;
- munged GWAS summary statistics for the later `partitioned-h2` step.

Each nonempty gene-list row must contain either an exact Ensembl gene ID
(version suffixes are accepted) or a case-sensitive gene name:

```text
ENSG00000100393
EP300
SYN1
```

Duplicates and aliases resolving to the same gene do not double-count SNPs.
Overlapping and nested genes are combined by Boolean union.

## Calculate indexed LD scores

```bash
INDEX_DIR="/path/to/gene_ldscore_indexes/1000G_EUR_Phase3_baseline_100kb"
GENE_LIST_SOURCES="/path/to/gene_sets/*.txt"
LDSCORE_OUTPUT_DIR="/path/to/results/gene_set_ldscores"

ldsc ldscore \
  --gene-ldscore-index-dir "${INDEX_DIR}" \
  --query-annot-gene-list-sources "${GENE_LIST_SOURCES}" \
  --output-dir "${LDSCORE_OUTPUT_DIR}"
```

**Equivalent direct mode:** the same scientific analysis can be run without an
index by supplying the original baseline annotations and PLINK panel and
repeating the index's scientific settings. See the
[matched direct-mode command](#equivalent-direct-mode-and-when-to-use-it).
Direct mode recomputes the baseline and gene-list LD scores; indexed mode
reuses the stored baseline block and exact gene operator, so matching inputs
should produce equivalent canonical results.

The explicit index owns the baseline, reference panel, genome build, SNP
identity, LD window, MAF settings, padding, and region policies. Do not pass
live `--baseline-annot-sources`, `--plink-prefix`, `--r2-dir`, build/window/map,
padding, or SNP-filter options in indexed mode. In particular, remove
`--padding-bp` rather than passing `--padding-bp 0`; indexed assembly inherits
the padding stored in the index.

**SNP identity in fast mode:** an index stores either `rsid` or `chr_pos` and
explicit hg19 provenance. Fast assembly does not repeat source matching; it
loads the immutable mode/build, aligned rows, and operator. Therefore, do not
pass either `--snp-identifier` or `--genome-build` to an indexed command. Their
presence is rejected even when the values equal the index. A matching direct
command must pass both explicitly because it rebuilds the source alignment.

In `chr_pos` output, normalized `CHR/POS` are identity and PLINK `SNP` remains
a useful published label, not a merge key. Both supported base modes also retain
PLINK `A1/A2` as passive metadata. Indexed assembly performs no liftover or
output-build conversion and writes a self-contained canonical LD-score
directory. Downstream `h2`, `rg`, and `partitioned-h2` apply the same ordinary
summary-statistics identity, build, and explicit downgrade rules as direct
outputs; they do not need the source index after assembly.

Use `--overwrite` to replace an existing output family.

## Understand the optional control annotation

No gene control is added by default. To condition focal gene sets on a specific
background universe, supply one existing one-column gene-list file:

```bash
--control-gene-list-file /path/to/background_genes.txt
```

The selected genes create the fixed `gene_control` baseline column. A custom
background can be appropriate when focal sets were selected from a restricted
assay universe. The flag accepts one real file, not a list, glob, or sentinel;
an unusable control stops the run. Do not name a focal list `gene_control`;
that name is reserved.

## Outputs

```text
gene_set_ldscores/
    metadata.json
    ldscore.baseline.parquet
    ldscore.query.parquet
    ldscore.overlap.parquet
    diagnostics/
        ldscore.log
        query_annotation_status.tsv
        gene_list_unresolved.tsv.gz
```

- `ldscore.baseline.parquet` contains the indexed baseline columns followed by
  `gene_control` when enabled.
- `ldscore.query.parquet` contains one focal column per usable gene-list file,
  in requested source order.
- `ldscore.overlap.parquet` stores the all-SNP and common-SNP overlap blocks
  needed by partitioned regression.
- `metadata.json` records annotation counts, `index_id`, catalog and
  control provenance, and SNP-universe policies.

The canonical rows are filtered HM3 regression SNPs, but their baseline and
gene-set LD scores were accumulated from the broad retained PLINK reference
universe. Non-HM3 SNPs—including SNPs in MHC and centromeric regions—can
therefore contribute to an HM3 row's LD score even though they are not written
as regression rows. Indexed mode fixes this regression-row policy; it is not an
HM3-only LD-reference calculation.

> **Caveat: the protein-coding catalog does not gate supplied baseline
> columns.** Their stored LD scores were computed from the full retained
> baseline/PLINK reference intersection, weighted by each baseline annotation.
> Gene regions restrict only `gene_control` and focal gene-list annotations. A
> SNP outside every protein-coding interval may still contribute to a supplied
> baseline LD score.

The output is self-contained. The index does not need to remain installed for
downstream regression.

## Check query resolution before regression

Inspect `diagnostics/query_annotation_status.tsv`. A focal query can be omitted
from the scientific tables if its source is unreadable, no identifiers resolve,
it hits no retained SNP, or it has zero variance. For gene-level details,
inspect `diagnostics/gene_list_unresolved.tsv.gz`.

If every focal query is skipped, only diagnostics are written and the command
exits with an error. Do not infer a biological null result from a missing query
column.

## Run partitioned S-LDSC

```bash
SUMSTATS_FILE="/path/to/sumstats_processed/trait/sumstats.parquet"
PARTITIONED_H2_OUTPUT_DIR="/path/to/results/partitioned-h2/trait"

ldsc partitioned-h2 \
  --sumstats-file "${SUMSTATS_FILE}" \
  --ldscore-dir "${LDSCORE_OUTPUT_DIR}" \
  --output-dir "${PARTITIONED_H2_OUTPUT_DIR}" \
  --write-per-query-results
```

The workflow loads all focal columns, then fits one model per focal query. With
the command above, each model contains the supplied baseline annotations and
one focal gene-set annotation. If a control is enabled, its `gene_control`
column is also included. Interpret the focal coefficient and its p-value as
conditional on the fitted model. Correlated or overlapping gene sets can
therefore produce related results and should not be read as mutually
independent discoveries.

## Equivalent direct mode and when to use it

Use direct gene-list mode when no compatible index exists or when the desired
reference, baseline, build, padding, gene-region policy, or control design
differs from the available index:

```bash
ldsc ldscore \
  --baseline-annot-sources "/path/to/baseline.@.annot.gz" \
  --plink-prefix "/path/to/1000G.EUR.QC.@" \
  --query-annot-gene-list-sources "${GENE_LIST_SOURCES}" \
  --snp-identifier rsid \
  --genome-build hg19 \
  --ld-wind-cm 1.0 \
  --common-maf-min 0.05 \
  --padding-bp 100000 \
  --gene-exclude-regions mhc \
  --exclude-regions mhc-and-centromeres \
  --output-dir "${LDSCORE_OUTPUT_DIR}"
```

Direct mode is the scientific oracle and fallback. Indexed mode is an exact
acceleration for a matching immutable configuration, not an approximation.

## Scaling guidance

Indexed assembly does not batch user-requested query columns. This preserves
one canonical query table and stable ordering, but very wide gene-set
collections can make `ldscore.query.parquet` and `partitioned-h2` memory the
dominant cost. Split analyses only when separate output families and separate
multiple-testing accounting are scientifically acceptable; otherwise allocate
memory for the complete query table.

To build the index used here, see
[Build an exact gene LD-score index](../utility-functionalities/build-gene-ldscore-index.md).
