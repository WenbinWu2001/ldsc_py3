# `annotate`

Last updated on: 2026-10-27

`ldsc annotate` projects BED intervals or gene lists onto a SNP grid and writes reusable query annotation files.
This page is the exact contract; the task guide is [Create reusable SNP annotations](../user-guide/workflows/annotations.md), and every option is listed by `ldsc annotate --help`.
The Python entry point is `run_annotate()` (`src/ldsc/annotate_workflow.py`).

## Inputs

A run takes exactly one query route, exactly one SNP grid, and `--output-dir`.

| Option | Contract |
| --- | --- |
| `--query-annot-bed-sources` | BED query route; one query annotation per BED file. Mutually exclusive with gene lists. |
| `--query-annot-gene-list-sources` | Gene-list query route; one query annotation per gene list, named after its file. |
| `--ref-snp-annot-sources` | Annotation grid: existing annotation files define the SNP rows, and their baseline values are exported aligned with the query columns. |
| `--ref-snp-bim-sources` | BIM grid: PLINK `.bim` or `.bim.gz` files define the SNP rows from `CHR`, `SNP`, one-based `BP` and, under allele-aware modes, `A1` and `A2`; `CM` is ignored, no `.bed` or `.fam` files are read, and no baseline values are created. |
| `--gene-coordinate-file` | Required for gene lists; the gene coordinate catalog is the only authority for gene identifiers and intervals. |
| `--padding-bp` | Required for gene lists as an explicit non-negative integer (`0` means gene bodies only); defaults to 0 for BED files. |
| `--allow-unresolved-genes` | Gene lists only. Without it, resolution is strict; with it, allowlisted rejected rows are omitted and the resolved subset is used. Python: `gene_list_resolution_policy=None` (strict) or `"resolved-only"`. |
| `--gene-exclude-regions` | Gene lists only: `none` (default) or `mhc`; a gene is excluded by its unpadded interval. |

Gene-list options are rejected on the BED route.
There is no control gene list option: in `ldscore` a control gene list becomes the baseline column `gene_control`, but `annotate` writes query columns only.
The annotation grid and the BIM grid are mutually exclusive, and a BIM grid cannot be combined with separately supplied baseline files.

## Genome build

On the gene-list route, an omitted `--genome-build` is treated as `auto` in every SNP identifier mode.
The build is resolved from the catalog and the grid evidence without any reference panel; an explicit build must agree with the catalog, and no liftover is performed.
Under rsID-family modes, the SNP identity stays build-independent and the concrete projection build is recorded separately.
On the BED route, coordinate-family modes require `--genome-build`, and the BED files inherit that build without a check.

## Projection

The catalog has one declared build and one-based inclusive intervals `[start, end]`.
Each interval is converted once to `[start - 1, end)`, padded on both sides, and clipped at 0.
For example, an unpadded catalog interval `[101, 110]` marks the one-based SNP positions 101 and 110 but not 100 or 111, the same as the BED interval `[100, 110)` (`tests/test_annotate_grid.py`).
A query annotation is the binary union of its intervals over the grid rows that survive SNP identity cleanup.
Gene identifiers match exactly and case-sensitively, as described in the [gene-list format](formats/gene-lists.md).

A query name that equals a baseline column or a reserved name stops the run; queries are never renamed.
The reserved names are `gene_control`, `CHR`, `POS`, `BP`, `SNP`, `CM`, `A1`, `A2` and `MAF` (`src/ldsc/annotate_workflow.py`).

## Validation and outcomes

### BED queries

Every requested BED file is read before the grid.
Any unreadable or malformed BED file stops the run with one `BED query content validation failed` error that lists every such file, also written to `diagnostics/input_issues.tsv`.
A BED file without intervals is skipped as `empty_input`.
A non-empty BED file that overlaps no grid row is kept as an all-zero column with status `ok`; the zero-support skip below applies to gene lists only.
If every BED query is skipped, the run fails after writing `diagnostics/query_annotation_status.tsv`.

### Gene-list queries

Gate A validates the catalog and resolves every gene list together before projection, collecting every source, row, catalog and naming issue it can detect.
A structurally unusable catalog is reported without a gene audit, and an unreadable list gets a source-level error with unknown row counts.
Chromosome scope then comes from the grid contents: exact paths and globs select their actual files, and an `@` chromosome suite requires all autosomes 1 to 22.
Gate B measures support on the grid rows that survive SNP identity cleanup; unknown support is blank, never zero.
`annotate` has no reference panel, so it makes no baseline-versus-reference SNP comparison.

| Condition | Outcome |
| --- | --- |
| All identifiers resolve and the inputs validate | Check coverage, then grid support. |
| Empty readable gene list | Skipped as `empty_gene_list` under either policy; other queries continue. |
| Rejected identifiers, strict resolution | The batch fails at Gate A after diagnostics are written. |
| Rejected identifiers, `--allow-unresolved-genes` | Only allowlisted rejected rows are omitted; a retained partial query carries `partial_gene_resolution`. |
| No selected genes left after allowed omissions or exclusions | Skipped as `zero_resolved_genes`. |
| Unreadable list, bad compression or encoding, malformed multi-field row, or duplicate query name | Fails under either policy. |
| Query name equals a baseline column or a reserved name | Fails; no renaming. |
| Structurally unusable catalog or incompatible build | Fails with the available catalog and build diagnostics. |
| Catalog defect that no listed gene uses | Reported without failing the run. |
| Repeated rows or aliases select the same gene | The gene is selected once; duplicate rows stay in the audit. |
| Gene removed by `--gene-exclude-regions` | Recorded as `excluded`; coverage and support use the remaining genes. |
| Selected gene outside the grid's chromosome scope | The batch fails under either policy; the gene list is never truncated. |
| Grid input missing, unreadable, invalid or outside its declared scope | The input gate fails after collecting the issues it can detect. |
| One gene overlaps no grid row | Recorded as `unsupported` with reason `zero_annotation_snp_support`; not fatal on its own. |
| A whole query overlaps no grid row | Skipped as `zero_annotation_snps`; other queries continue. |
| A query is zero on one chromosome but supported elsewhere | Kept in every chromosome file, including the zeros. |
| Every query is skipped | Fails after diagnostics are written; no new query annotation files are published. |
| An annotation is 1 on every grid row | Kept; annotation generation has no variance or regression-design check. |

Precedence follows `gene_query_statuses()` in `src/ldsc/query_annotations.py`: incomplete coverage stops the run before support is assessed, a query with zero support is skipped, and `partial_gene_resolution` stays the primary warning when partial support also occurs.

### Command-line notices

After a successful run with omitted genes, unsupported genes or skipped queries, the command prints a short notice to stderr that names the affected queries and the diagnostics file to inspect.
Python callers receive the statuses and output paths instead; `run_annotate()` prints nothing.

## Outputs

| Path | Content |
| --- | --- |
| `query.<chrom>.annot.gz` | Query annotation columns, identical and in the same order in every chromosome file; skipped queries are omitted. |
| `baseline.<chrom>.annot.gz` | Annotation grid only: the supplied baseline values on the cleaned grid rows, as float32, including continuous values; the original files are not changed. |
| `diagnostics/metadata.json` | Provenance; fields are listed in [artifact metadata](formats/artifact-metadata.md#annotate). |
| `diagnostics/annotate.log` | Workflow log. |
| `diagnostics/dropped_snps/dropped.tsv.gz` | Grid rows dropped by SNP identity cleanup. |
| `diagnostics/query_annotation_status.tsv` | One status row per query. |
| `diagnostics/gene_list_audit.tsv.gz`, `diagnostics/gene_list_resolution_summary.tsv` | Gene lists only; columns are described in [gene-list diagnostics](results/gene-list-diagnostics.md). |
| `diagnostics/gene_catalog_issues.tsv`, `diagnostics/input_issues.tsv`, `diagnostics/chromosome_scope.json` | Written when the corresponding issue or scope record exists. |

Annotation files use the columns `CHR`, `BP`, `SNP` and `CM` with an explicit `CM=NA`, keep the allele columns, carry no `MAF`, and store query values as integers 0 and 1.
Pass the exported files to `ldscore` with `--baseline-annot-sources` and `--query-annot-sources`; columns from a BIM grid can also serve as baseline annotations when they define the intended baseline model.
Diagnostics are written before scientific outputs, and output preflight, `--overwrite` and failure markers follow the [logging and output rules](../current/workflow-logging.md#output-family-preflight).
