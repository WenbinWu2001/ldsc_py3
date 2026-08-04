# Calculate LD scores for gene lists with an index

Last updated on: 2026-08-04

## Motivation

Gene expression, proteomic, GO, and SynGO analyses often produce many related
gene sets. The scientific question is whether common-variant heritability is
enriched near a focal set after accounting for ordinary functional annotations
and the general tendency of SNPs near genes to carry heritability.

LDSC-SEG used S-LDSC for exactly this kind of conditional gene-set test: a
specifically expressed gene annotation was fitted with the baseline model and
an all-genes annotation ([Finucane et al., 2018](https://doi.org/10.1038/s41588-018-0081-4). An exact index makes the LD-score stage fast to repeat across many focal gene lists; the downstream `partitioned-h2` model and interpretation remain S-LDSC.

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
  --control-gene-list-source all-protein-coding \
  --output-dir "${LDSCORE_OUTPUT_DIR}"
```

The explicit index owns the baseline, reference panel, genome build, SNP
identity, LD window, MAF settings, padding, and region policies. Do not pass
live `--baseline-annot-sources`, `--plink-prefix`, `--r2-dir`, build/window/map,
padding, or SNP-filter options in indexed mode.

V1 indexes and their canonical LD-score rows are hg19. Indexed assembly does
not offer output-build conversion: it reuses the index's hg19 coordinates,
gene projection, regression masks, and map/window identity. Use direct mode for
a supported configuration outside the index compatibility domain; do not treat
rsID matching as permission to mix coordinate builds.

Use `--overwrite` to replace an existing output family.

## Understand the control annotation

The default `all-protein-coding` control creates the fixed `gene_control`
column and appends it to the baseline block. This plays the role of the
all-genes annotation in the LDSC-SEG design: a focal coefficient is interpreted
conditional on baseline annotations and general gene proximity, rather than as
a comparison of genic versus non-genic SNPs.

Alternative controls are explicit scientific choices:

```bash
# Disable the fixed gene control.
--control-gene-list-source none

# Use one custom background gene list.
--control-gene-list-source /path/to/assayed_protein_coding_genes.txt
```

A custom background can be appropriate when the focal sets were selected from
a restricted assay universe. An unusable custom control stops the run. Do not
name a focal list `gene_control`; that name is reserved.

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

The workflow loads all focal columns, then fits one model per focal query. Each
model contains the supplied baseline annotations, `gene_control` when enabled,
and one focal gene-set annotation. Interpret the focal coefficient and its
p-value as conditional on that model. Correlated or overlapping gene sets can
therefore produce related results and should not be read as mutually
independent discoveries.

## When to use direct mode instead

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
  --padding-bp 100000 \
  --gene-exclude-regions mhc \
  --control-gene-list-source all-protein-coding \
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
