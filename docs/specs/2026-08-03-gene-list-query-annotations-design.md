# Gene-List Query Annotations Design

Last updated on: 2026-08-03

Status: implemented; the [exact-index specification](2026-08-03-exact-disjoint-atom-gene-ldscore-index-design.md) supersedes only its explicitly listed direct-workflow changes when implemented

## Problem and goal

`ldsc ldscore` accepts prebuilt SNP annotations and BED query annotations, but it does not accept the gene lists that users commonly maintain. Users must currently resolve identifiers, choose a genome build, convert genes to intervals, and construct BED files outside LDSC. That creates avoidable coordinate, naming, and provenance mistakes.

This design adds gene-list query annotations directly to `ldsc ldscore`. Each input list is resolved against one packaged protein-coding gene catalog, projected onto the same baseline SNP grid as a BED query, and processed as an ordinary exact query annotation. The feature must be scientifically equivalent to supplying the resolved gene intervals through `--query-annot-bed-sources`, while giving batch users a durable per-query status record instead of aborting because one query is unusable.

This specification supersedes the earlier Session 1 handoff wherever the approved batch partial-success behavior changes the previous BED or empty-list behavior.

## Scope

- one new `ldsc ldscore` input route, `--query-annot-gene-list-sources`;
- a packaged hg19/hg38 protein-coding gene catalog;
- exact Ensembl gene-ID and gene-name resolution;
- in-memory interval projection using the existing BED coordinate and padding rules;
- shared query-local status handling for gene-list and BED query sources;
- concise catalog/query provenance and row-level unresolved-gene diagnostics;
- ordinary exact LD scores, annotation counts, overlap artifacts, and downstream `partitioned-h2` compatibility.

## Out of scope

- the future sparse disjoint-atom gene LD-score index;
- a new command or a gene-list route for `ldsc annotate`;
- physical intermediate BED or `.annot.gz` files for gene-list inputs;
- fuzzy, case-insensitive, synonym, or historical gene-name matching;
- transcript, exon, strand-aware, TSS-only, or TES-only annotations;
- non-protein-coding or non-autosomal catalog entries;
- combining gene-list, BED, and prebuilt query source types in one invocation;
- an all-protein-coding-genes control annotation or a new partitioned-h2 mode.

## User-visible command contract

`ldsc ldscore` adds:

```text
--query-annot-gene-list-sources <source-group>
```

The source group accepts comma-separated exact file paths and glob tokens through the existing deterministic file-group resolver. Explicit files retain user order, each glob expands in lexical order, and duplicate resolved paths keep their first occurrence. A token that resolves to no files is a run-global usage error. Gene lists do not use chromosome-suite `@` expansion.

The new flag joins the existing mutually exclusive query-source group:

- `--query-annot-sources`
- `--query-annot-bed-sources`
- `--query-annot-gene-list-sources`

Any query route requires explicit `--baseline-annot-sources`. The synthetic all-ones baseline remains available only when no query source is requested. `--bed-padding-bp` applies symmetrically and exactly once to BED-derived and gene-derived intervals; starts are clipped at zero.

The Python workflow configuration adds the corresponding normalized source tuple to `AnnotationBuildConfig` and `run_ldscore(...)` accepts the matching CLI-style keyword. The gene resolver and its records remain internal and are not re-exported from `ldsc.__init__`.

## Packaged catalog contract

The durable package resource is:

```text
src/ldsc/data/protein_coding_genes.tsv.gz
```

The decompressed TSV schema is:

| Column | Contract |
| --- | --- |
| `ensgid` | Unique, nonempty canonical Ensembl gene ID. |
| `gene_name` | Nonempty exact lookup alias; duplicates are allowed and make that name ambiguous. |
| `hg38_chr`, `hg38_start0`, `hg38_end` | hg38 interval, either all missing or a valid autosomal 0-based half-open interval. |
| `hg19_chr`, `hg19_start0`, `hg19_end` | hg19 interval, either all missing or a valid autosomal 0-based half-open interval. |
| `hg38_strand`, `hg19_strand` | Optional retained source data; unused by symmetric padding and projection. |

Valid coordinates have a normalized autosomal chromosome, `start0 >= 0`, and `end > start0`. A coordinate triplet must be wholly present or wholly missing. A missing triplet is a valid catalog state that can produce a per-gene `build_missing` outcome. Duplicate Ensembl IDs, partial/invalid coordinates, incompatible columns, or invalid chromosomes make the catalog invalid and abort the whole run before query computation.

Intervals use the existing BED contract `[start0, end)`. Baseline SNP `POS` remains a canonical 1-based point coordinate and is projected as `[POS - 1, POS)`. This separation preserves exact BED equivalence and does not change SNP artifact schemas.

### Catalog provenance

The upstream gene release is GENCODE Human release 49, from `gencode.v49.basic.annotation.gtf.gz`. The preparation record is `/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/resources/gene_lists/protein-coding/12-pc-gene-list.R`, created in June 2026. It filters to autosomal protein-coding genes and excludes duplicate-name, novel, and read-through genes. It converts the GTF-derived hg38 start from 1-based to a 0-based BED start and preserves the end.

The preparation record says its `geneMatrix` source contains hg19 coordinates but does not document their exact derivation. Runtime and user documentation therefore identify GENCODE v49 as the release source without claiming an undocumented hg19 liftover method.

Runtime metadata records one concise catalog object per run: resource basename, release `GENCODE v49`, selected interval build, and a SHA-256 over the exact decompressed TSV bytes. The catalog checksum is not repeated for every query. It binds compact catalog indices to exact catalog content for future index compatibility.

## Gene-list parsing and identity resolution

A gene-list file is one-column and headerless.

- Leading and trailing whitespace is removed.
- Blank or whitespace-only lines are ignored and counted.
- Every nonblank line is treated as an identifier; headers and comments have no special syntax.
- A nonblank line containing more than one tab-delimited field is `malformed_input`.
- A readable file is scanned to completion so all problematic line records are collected before its final status is assigned.

One file may mix Ensembl gene IDs and exact gene names. Matching is case-sensitive. Resolution checks the exact Ensembl-ID namespace first and the exact gene-name namespace second. A terminal numeric Ensembl version suffix is stripped before lookup, so `ENSG00000123456.12` resolves as `ENSG00000123456`. An Ensembl-like token with a dotted nonnumeric suffix is `invalid_identifier`; an otherwise well-formed Ensembl ID absent from the catalog is `unmatched_identifier`.

Ensembl gene ID is the canonical identity. Gene names are lookup aliases only. No case folding, fuzzy matching, synonym table, or alias expansion is performed.

If one token maps through both namespaces to different genes, or an exact gene name maps to multiple Ensembl IDs, the token is `ambiguous_identifier` and the entire query is skipped. Diagnostics list all conflicting canonical Ensembl IDs. This protects the scientific annotation from an arbitrary alias choice.

### Duplicate and count semantics

Annotation construction uses unique canonical Ensembl genes. Repeated inputs and different accepted aliases for the same gene cannot increase annotation values.

Compact per-query accounting includes:

- nonblank input rows;
- unique normalized input tokens;
- repeated-token rows;
- matched input rows;
- unique resolved canonical genes;
- alias-collapsed rows, where different accepted identifiers resolve to the same canonical gene;
- blank rows;
- unresolved counts grouped by reason.

Duplicate identifiers and blank lines alone are benign and do not downgrade an otherwise usable query.

## Genome-build selection

Gene-list projection reuses `--genome-build`; no gene-list-specific build flag is added. Existing aliases normalize through shared rules: `hg37` and `GRCh37` select hg19, and `GRCh38` selects hg38.

When gene-list input is present and `--genome-build` is omitted, its effective default is `auto`, including in rsID-family modes. A concrete build selects that catalog interval set. In `auto`, the workflow uses existing LD-score evidence from baseline annotation coordinates and R2 parquet metadata; gene identifiers themselves provide no build evidence. All available evidence must agree. Insufficient or conflicting evidence is a run-global error with guidance to pass hg19 or hg38 explicitly.

For rsID-family workflows, the selected build controls interval projection and selection of named regression-region presets. SNP matching remains rsID-based and LD-score compatibility metadata keeps `genome_build=null`. The separately recorded catalog provenance carries the effective projection build. The declared or inferred projection build and the evidence used are reported in `diagnostics/ldscore.log`.

## Query naming

Each gene-list query name comes only from its resolved source basename. An optional `.gz` is removed, followed by at most one recognized gene-list suffix: `.txt`, `.tsv`, or `.list`. Other dots are retained.

Examples:

| Source basename | Query name |
| --- | --- |
| `immune_genes.txt.gz` | `immune_genes` |
| `immune.v2.tsv.gz` | `immune.v2` |
| `immune.custom.gz` | `immune.custom` |

BED naming remains backward compatible with its existing single-`Path.stem` behavior; for example, `immune.bed.gz` remains `immune.bed`.

Duplicate derived query names and collisions with baseline annotation columns are run-global preflight errors. The message names the conflicting query and source basenames. Directory components do not disambiguate names.

## Projection and scientific invariants

The resolver selects the requested build's intervals for unique canonical genes. `--bed-padding-bp` expands each interval on both sides once and clips starts to zero. Intervals are unioned before projection, so overlapping genes or padding never produce annotation values greater than one. No physical BED file is written.

Gene-derived and BED-derived queries use the same interval-overlap primitive and baseline SNP grid. For the same normalized interval union they must produce identical:

- binary SNP annotations;
- all-reference and common-reference annotation counts;
- query LD scores;
- overlap-matrix entries;
- retained regression SNP rows.

`n_annotation_snps`, `M`, `M_5_50`, and overlap entries are evaluated over `ld_reference_snps` after any explicit reference-SNP restriction and ordinary MAF retention. Named region exclusions do not remove LD-score contributors or change those counts; they are subtracted only from regression/output rows and `w_ld` contributors. Query LD-score variance is evaluated over the resulting written `ld_regression_snps` rows.

Only usable queries contribute columns to `AnnotationBundle`, `LDScoreResult`, `ldscore.query.parquet`, count records, or `ldscore.overlap.parquet`. Column order follows resolved source order after skipped queries are removed. Baseline columns and ordinary exact LD-score math are unchanged.

## Batch partial-success contract

BED and gene-list sources share query-local partial success. A problem in one concrete source does not interrupt usable siblings. Prebuilt query `.annot` inputs retain their existing contract and are not part of this batch-status change.

The workflow assigns one of these statuses to every requested BED or gene-list query:

| Status | Meaning |
| --- | --- |
| `ok` | Query is usable and has no non-benign resolution issue. |
| `warning` | Query is usable, but a gene list was only partially resolved. |
| `skipped` | Query is excluded from scientific outputs for a query-local failure. |

Allowed reasons are:

| Status | Reason | Trigger |
| --- | --- | --- |
| `ok` | empty | Usable clean query; duplicate/blank gene rows alone remain clean. |
| `warning` | `partial_resolution` | At least one selected-build gene is usable and at least one input is unmatched, invalid, or `build_missing`. |
| `skipped` | `empty_input` | No nonblank gene identifiers or no BED intervals. |
| `skipped` | `malformed_input` | Structurally invalid gene-list or BED content. |
| `skipped` | `fully_unresolved` | No input identifier resolves to a gene with selected-build coordinates. |
| `skipped` | `ambiguous_identifier` | A gene token has more than one possible canonical identity. |
| `skipped` | `unreadable_source` | A resolved concrete file cannot be opened, decoded, or decompressed. |
| `skipped` | `zero_annotation_snps` | Projection contains no retained LD-reference SNPs. |
| `skipped` | `zero_variance_ld_scores` | Query LD scores have zero variance on written regression rows. |

Malformed or ambiguous gene-list content skips the full query even if other rows resolve. Unmatched identifiers, invalid identifiers, and missing selected-build coordinates produce `warning` when at least one gene remains usable. A skipped query is never represented by an all-zero or placeholder scientific column.

The run succeeds if at least one requested query is `ok` or `warning`. If every requested query is skipped, the workflow writes the status manifest and log, then raises one consolidated input error and writes no root `metadata.json` or canonical LD-score parquet artifacts.

Run-global failures still abort the batch: an invalid packaged catalog, unresolved/conflicting build, missing baseline or reference panel, a source token matching no files, query-name collisions, baseline-column collisions, output collisions, and unexpected internal errors.

## Diagnostics and logging

### Query status manifest

Every BED/gene-list run writes `diagnostics/query_annotation_status.tsv` with one row per requested concrete source in deterministic source order:

| Column | Contract |
| --- | --- |
| `query` | Derived query annotation name. |
| `source` | Source basename only. |
| `input_type` | `bed` or `gene_list`. |
| `status` | `ok`, `warning`, or `skipped`. |
| `reason` | Exact reason above; empty for `ok`. |
| `n_annotation_snps` | Retained-reference-SNP annotation count before regression-region subtraction; null if not computable. |
| `details` | Concise explanation or relative diagnostic reference; empty for `ok`. |

`n_annotation_snps` is the query's `M_c` over retained `ld_reference_snps`. It is zero for a constructed query that hits no retained SNP and null when parsing/resolution/read failure prevents construction.

### Unresolved-gene audit

Every gene-list run writes `diagnostics/gene_list_unresolved.tsv.gz`, including a header-only file for a clean run. It contains unresolved/problematic entries only:

| Column | Contract |
| --- | --- |
| `query` | Derived query name. |
| `source` | Source basename only. |
| `line` | One-based input line number. |
| `input_gene` | Original trimmed token or malformed row text. |
| `reason` | `unmatched_identifier`, `invalid_identifier`, `build_missing`, `ambiguous_identifier`, or `malformed_input`. |
| `canonical_ensembl_id` | Canonical ID when one gene was identified but its selected-build interval is missing; otherwise empty. |
| `details` | Concise conflict or format details; otherwise empty. |

Absolute paths are excluded from metadata and sidecars. Resolved paths may appear in `diagnostics/ldscore.log` for troubleshooting.

Every non-`ok` query emits a WARNING in `diagnostics/ldscore.log` naming the query, status, reason, and relevant diagnostic file. The log includes a final status-count summary and the effective catalog build/evidence. Console behavior follows the existing file-authoritative logging contract.

The status manifest, unresolved audit, and LD-score log are owned output artifacts. They participate in family collision preflight. After an all-skipped failure, a corrected rerun in the same output directory requires `--overwrite`.

## Root metadata and provenance

For successful gene-list runs, root `metadata.json` retains the current LD-score schema and adds only fields that support reproducibility or diagnostics:

- one catalog object with resource, release, selected projection build, and decompressed-content SHA-256;
- an ordered compact query-provenance list containing source ordinal, query name, source basename, SHA-256 of the exact input file bytes, and the approved parsing/resolution counts;
- relative references to `diagnostics/query_annotation_status.tsv` and `diagnostics/gene_list_unresolved.tsv.gz`.

Successful canonical Ensembl IDs are not copied into metadata or a resolved-gene sidecar. The resolver keeps canonical IDs and compact catalog indices in memory. BED-only runs add the status-manifest reference but no catalog or gene-resolution metadata.

The shared `artifact_type` guard is unchanged by these additive LD-score metadata fields. Readers must continue to ignore additive fields they do not consume.

## Internal module boundary and performance

Gene-list parsing/resolution belongs in a workflow-layer module such as `src/ldsc/gene_list_resolver.py`. It owns catalog schema validation and pure resolution records, but performs no logging and writes no files. LD-score/annotation orchestration owns build selection, warnings, query status, preflight, sidecars, and metadata.

Per invocation, the implementation must:

- load and validate the packaged catalog once;
- build Ensembl-ID and gene-name maps once;
- parse each concrete list once;
- retain canonical Ensembl IDs, catalog indices, selected-build intervals, compact counts, and unresolved rows without expanding a dense gene-by-query matrix;
- preserve deterministic query and diagnostic row order.

Catalog indices are local to the exact catalog checksum. Any future persisted index must store and validate that checksum rather than treating row positions as stable across catalog revisions.

## Acceptance criteria and validation

Implementation is accepted when automated tests demonstrate all of the following:

1. The package resource is installed and importable, has the exact schema above, contains unique Ensembl IDs, and validates hg19/hg38 missing-coordinate rules.
2. Plain and gzip gene lists resolve identically; Ensembl versions, mixed IDs/names, blanks, duplicates, alias collapse, unmatched values, invalid values, build-missing genes, malformed rows, and both ambiguity forms follow this specification.
3. Source ordering, suffix removal, collision errors, and source-basename-only provenance are deterministic.
4. Explicit hg19/hg38 aliases and auto inference behave in chr_pos and rsID families; rsID artifact compatibility remains build-independent and inferred build evidence appears in the log.
5. A gene list and an equivalent BED interval union produce identical annotation vectors, counts, LD scores, overlap values, and regression rows with and without padding.
6. Mixed valid/invalid BED and gene-list batches continue with usable queries, emit exact status rows and warnings, and omit failed columns from all scientific artifacts.
7. Partial resolution writes every unresolved/build-missing row and compact counts; a clean run writes a header-only unresolved audit.
8. Zero-hit and zero-variance queries are skipped rather than persisted; all-skipped runs write owned diagnostics, raise, and leave no canonical result artifacts.
9. Existing prebuilt-query, BED naming, baseline-only, overlap-aware partitioned-h2, output collision, logging, and Python API tests remain green.
10. User documentation explains accepted formats, coordinate/build behavior, query naming, partial-success semantics, and that a missing requested query in scientific results means `diagnostics/query_annotation_status.tsv` must be checked.

## Risks and settled constraints

- Gene names are less stable and may be ambiguous; canonical Ensembl identity and fail-closed ambiguity handling are fixed safeguards.
- hg19 derivation is not documented by the source preparation script; the implementation must not invent more specific provenance.
- Query-local diagnostics change existing BED failure behavior intentionally. Tests must distinguish these failures from true global preflight/configuration errors.
- Filtering zero-variance LD scores occurs after computation, so status records must be preserved through chromosome aggregation without allowing result-column/count/overlap ordering to diverge.
- Gzip output must be deterministic enough for content checks where tested; catalog identity is based on decompressed TSV content, not gzip container metadata.

There are no unresolved product decisions in this specification. Any implementation discovery that would change a public flag, status/reason vocabulary, coordinate semantics, artifact layout, or scientific inclusion rule requires an explicit design update before coding continues.
