# Gene-List Query Input

Last updated on: 2026-08-04

`ldsc ldscore` accepts one or more gene lists as query annotations through
`--query-annot-gene-list-sources`. Each list is resolved against the packaged
protein-coding catalog and projected directly onto the baseline SNP grid. No
intermediate BED or `.annot.gz` file is created.

## Input format

Each input is a plain-text or gzip-compressed, one-column, headerless file. Blank
lines are ignored. Every other line is an identifier; comments and headers have
no special syntax.

Accepted identifiers are:

- canonical Ensembl gene IDs such as `ENSG00000141510`;
- Ensembl IDs with a numeric version suffix, such as `ENSG00000141510.18`;
- exact, case-sensitive gene names such as `TP53`.

One file may mix IDs and names. Ensembl ID is always the canonical identity;
gene names are lookup aliases. Matching is not fuzzy or case-insensitive and
does not use synonym or historical-name tables. A name that maps to multiple
genes, or a token that conflicts across the ID and name namespaces, makes that
query ambiguous and causes it to be skipped. The gene-level audit lists the
conflicting Ensembl IDs.

Repeated tokens, blank rows, and different accepted aliases for the same gene
do not increase annotation values. A nonblank row with more than one
tab-delimited field is malformed.

## Query names

The query annotation name is derived from the source basename. An optional
`.gz` is removed first, followed by at most one of `.txt`, `.tsv`, or `.list`.

| Source | Query name |
| --- | --- |
| `immune_genes.txt.gz` | `immune_genes` |
| `immune.v2.tsv.gz` | `immune.v2` |
| `immune.custom.gz` | `immune.custom` |

Names must be unique and must not collide with baseline annotation columns.

## Catalog and coordinates

The package resource is `src/ldsc/data/protein_coding_genes.tsv.gz`. Its hg38
source is GENCODE Human release 49, based on
`gencode.v49.basic.annotation.gtf.gz`. Project curation keeps autosomal
protein-coding genes and excludes duplicate-name, novel, and read-through
genes. The exact derivation of the added hg19 coordinates is not documented by
the upstream preparation record, so LDSC reports the source limitation without
claiming a liftover method.

Catalog intervals are 0-based and half-open: `[start0, end)`. Baseline SNP
positions remain 1-based and a SNP at `POS=p` is treated as `[p-1, p)`. This is
the same coordinate contract as query BED files. `--padding-bp N` expands
both gene and BED intervals by `N` bases on each side, clips starts at zero, and
is applied exactly once. Overlapping intervals are unioned, so the annotation
remains binary. The effective default is `0`. This option belongs only to live
BED or gene-list projection; it is rejected for prebuilt annotation queries,
no-query LD-score runs, and indexed mode, where the index's stored padding is
authoritative.

## Genome build

Gene-list projection reuses `--genome-build`; aliases such as `GRCh37` and
`GRCh38` are accepted. If the flag is omitted for a gene-list run, its effective
default is `auto`. Automatic selection uses available baseline and R2-panel
evidence, requires the evidence to agree, and reports the inferred build in
`diagnostics/ldscore.log`.

In rsID modes the selected build controls coordinate-based gene projection and
named regression-region presets without changing SNP identity. Root LD-score
compatibility metadata keeps `genome_build: null`; the catalog provenance
separately records the projection build. A concrete automatically inferred
projection build is also used to select the default named exclusions. If
automatic inference has insufficient evidence, pass `--genome-build hg19` or
`--genome-build hg38` explicitly.

Gene-derived annotations follow the same SNP-universe contract as BED queries.
Named region exclusions remove only regression/output rows and `w_ld`
contributors; they do not remove SNPs from query LD-score calculation,
`n_annotation_snps`, `M`, `M_5_50`, or overlap counts. See
[LD-score SNP-universe contract](ldscore-snp-universe-contract.md).

## Example

```bash
ldsc ldscore \
  --output-dir results/gene_list_ldscores \
  --baseline-annot-sources "annotations/baseline.@.annot.gz" \
  --query-annot-gene-list-sources "gene_lists/*.txt.gz" \
  --r2-dir ref_panel/hg38 \
  --snp-identifier chr_pos_allele_aware \
  --genome-build auto \
  --ld-wind-cm 1.0
```

Gene-list, BED, and prebuilt query routes are mutually exclusive in one run.
Every query route requires explicit baseline annotations.

## Partial success and diagnostics

Each concrete BED or gene-list source is handled independently. A usable,
partially resolved gene list is retained with status `warning`; an empty,
malformed, ambiguous, unreadable, fully unresolved, zero-hit, or zero-variance
query is `skipped`. Valid siblings continue. Skipped queries do not receive
placeholder columns, counts, overlap entries, or downstream partitioned-h2
rows. If all queries are skipped, LDSC writes diagnostics and then exits with a
consolidated input error without writing canonical root LD-score artifacts.

The output diagnostics are:

- `diagnostics/query_annotation_status.tsv`: one row per requested BED or gene
  list, including `status`, `reason`, retained annotation-SNP count, and details;
- `diagnostics/gene_list_unresolved.tsv.gz`: unresolved/problematic gene rows
  only, with a header-only file for a clean gene-list run;
- `diagnostics/ldscore.log`: the inferred projection build, one warning per
  non-`ok` query, and a final status summary.

If an input query gene list or BED file is absent from the scientific results,
that query failed. Check `diagnostics/query_annotation_status.tsv` for the
reason, then `diagnostics/gene_list_unresolved.tsv.gz` for gene-level details.

By default, gene-list runs do not add a fixed gene control. To condition on a
specific background gene set, pass one existing one-column file with
`--control-gene-list-file`; it is projected as the baseline column
`gene_control`. `--gene-exclude-regions mhc` removes genes whose unpadded
transcribed intervals overlap the packaged MHC interval before padding; it is
independent of SNP `--exclude-regions` and applies to both focal and control
gene lists.

For an installed exact index, replace the live baseline/reference arguments
with `--gene-ldscore-index-dir <index-dir>`. Resolution then uses only the
index's embedded catalog and immutable padding/exclusion policy. See
[gene-ldscore-index.md](gene-ldscore-index.md) and the task-oriented
[indexed gene-list LD-score tutorial](../wiki/main-functionalities/ldscore.md).
