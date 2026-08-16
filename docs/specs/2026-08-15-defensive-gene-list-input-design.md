# Defensive Gene-List Input Specification

Last updated on: 2026-08-16

Status: approved for implementation

## Problem and goal

Legacy LDSC2 maps gene lists to coordinates with an inner join, so unmatched genes disappear silently. The current LDSC3 resolver records unresolved genes but normally continues with the resolved subset. Either behavior can produce a successful run whose scientific gene set differs from the submitted set.

This feature makes that behavior explicit and auditable:

- strict gene-list resolution is the default;
- an exploratory `resolved-only` policy deliberately permits partial resolution for large pathway batches;
- direct and indexed workflows use one authoritative, build-aware gene-coordinate catalog;
- focal and control lists are validated together before substantial LD-score work;
- all row-level outcomes, per-source counts, query viability, and catalog-builder defects have deterministic diagnostics; and
- users can repair their inputs from a dedicated diagnostics-and-repair guide without reconstructing resolver behavior from a log.

The governing operational principle is:

> Assume users do not inspect file logs after a successful run, especially under SLURM. An unrequested resolution, filtering, ambiguity, or compatibility condition that can change the scientific analysis must stop before substantial computation. An explicitly requested transformation may continue only when its effect is recorded in result diagnostics and provenance.

## Scope

This specification covers:

- direct and indexed `ldsc ldscore` gene-list query annotations;
- focal lists from `--query-annot-gene-list-sources`;
- the optional control list from `--control-gene-list-file`;
- the required live gene-coordinate catalog and the catalog embedded in an index;
- strict and exploratory resolution policies;
- gene exclusion, padding, retained-SNP support, and unusable-query behavior;
- index catalog validation and failure diagnostics;
- CLI mode separation, artifact provenance, documentation, and performance; and
- removal of the packaged protein-coding catalog and incompatible old index schema.

## Observable behavior

Implementation is acceptable only when these scenarios hold.

1. In strict mode, any focal or control identifier that cannot produce one valid catalog interval causes one consolidated Gate A failure. The error reports every safely discoverable problem across all sources, writes the complete row audit and per-source summary, and performs no LD-score calculation.
2. With `--gene-list-resolution-policy resolved-only`, approved row-level failures are omitted explicitly. Usable focal queries and a usable partial control continue, and results record the effective policy and exact resolved/total counts.
3. A malformed or unreadable list, duplicate query name, structurally invalid catalog, build conflict, corrupt index, or incompatible CLI mode remains fatal under both policies.
4. A focal list with zero resolved genes under `resolved-only`, zero retained annotation SNPs, or zero-variance LD scores is skipped for that run. Usable sibling queries continue. If every focal query is skipped, the run exits nonzero and does not publish misleading baseline-only scientific output.
5. A requested control is never silently dropped. An empty or zero-resolved control, zero-SNP control annotation, or zero-variance control LD-score column is fatal.
6. A gene with zero retained reference-SNP support is nonfatal. It is recorded in diagnostics and contributes to the per-source warning state, but it is not misclassified as a catalog-resolution failure.
7. Explicit MHC exclusion is an audited, non-error transformation. It is applied to the unpadded gene interval before padding and affects both focal and control lists.
8. Index construction validates the entire source catalog as canonical before atom construction. Any catalog defect aborts the build and writes every detectable defect to the catalog-issues artifact.
9. A public index always covers autosomes 1 through 22. Public indexed `ldscore` rejects partial coverage at the artifact boundary.
10. Live gene-derived and equivalent BED-derived interval unions remain scientifically equivalent after the required coordinate conversion and identical padding.

## CLI and mode contract

`ldsc ldscore` remains one command. It has four mutually exclusive query-input modes plus the existing no-query state.

| Mode | Required | Optional | Forbidden |
| --- | --- | --- | --- |
| Live gene list | `--query-annot-gene-list-sources`, `--gene-coordinate-file`, explicit `--padding-bp`, live baseline/reference inputs | `--control-gene-list-file`, `--gene-exclude-regions`, `--gene-list-resolution-policy` | Index, BED query, prebuilt query |
| Indexed gene list | `--query-annot-gene-list-sources`, `--gene-ldscore-index-dir` | `--control-gene-list-file`, `--gene-list-resolution-policy` | Live coordinate override, padding, live gene-exclusion override, live baseline/reference inputs, BED query, prebuilt query |
| Prebuilt query annotation | `--query-annot-sources`, live baseline/reference inputs | Ordinary LD-score controls | Gene-list, coordinate, control-gene, index, BED, explicit padding |
| Live query BED | `--query-annot-bed-sources`, live baseline/reference inputs | `--padding-bp` | Gene-list, coordinate, control-gene, gene exclusion, index, prebuilt query |

`--gene-ldscore-index-dir` is a backend selector, not a query source. It requires focal gene-list sources. Baseline-only/no-query LD-score calculation remains valid.

Padding has mode-specific meaning:

| Mode | Omitted | Explicit `0` | Explicit positive value |
| --- | --- | --- | --- |
| Live gene list | Error: the user must choose | Gene bodies | Padded genes |
| Live BED | Equivalent to `0` | Unpadded BED | Padded BED |
| Prebuilt annotation or no query | Accepted | Rejected | Rejected |
| Indexed gene list | Required omission; inherit index | Rejected | Rejected |

Gene-list modes add:

```text
--gene-list-resolution-policy {strict,resolved-only}
```

The default is `strict`. The option applies jointly to all focal lists and the optional control. It is invalid outside gene-list modes.

## Gene-coordinate catalog contract

### Required authority and format

Live gene-list mode and `ldsc build-gene-ldscore-index` require `--gene-coordinate-file`. There is no packaged fallback and no catalog overlay. The supplied catalog is the sole resolution authority for both focal and control lists.

The input is a headered TSV or TSV.GZ with these required columns in any order:

| Column | Contract |
| --- | --- |
| `gene_id` | Nonempty authoritative identity. |
| `gene_name` | Optional exact lookup alias; blank creates no alias. |
| `chrom` | Autosomal chromosome. |
| `start` | One-based inclusive start. |
| `end` | One-based inclusive end. |
| `genome_build` | Declared build; one catalog represents one build. |

The first nonblank, noncomment line is the header. Blank lines and lines beginning with `#` are ignored while physical line numbers are retained. Extra columns are accepted and scientifically ignored.

Catalog coordinates are one-based closed intervals `[start, end]`, with `start >= 1` and `end >= start`. Internal BED/projection coordinates are exactly:

```text
start0 = start - 1
end0_exclusive = end
```

Only the start is decremented.

### Normalization and build

The reader:

- trims surrounding whitespace from headers and cell values;
- preserves case, version suffixes, and internal characters in gene identities;
- never strips Ensembl or other identifier versions;
- normalizes `1` and `chr1` to chromosome `1`;
- normalizes `hg19`, `hg37`, and `GRCh37` to `hg19`;
- normalizes `hg38` and `GRCh38` to `hg38`; and
- accepts only integer text for `start` and `end`, rejecting values such as `1.0`.

The catalog build is mandatory evidence. Explicit `--genome-build` must agree with it. Under `--genome-build auto`, catalog evidence must agree with authoritative baseline/reference metadata. The workflow performs no implicit liftover.

### Identity and catalog validity

Gene matching is exact and case-sensitive. Users should prefer `gene_id`; `gene_name` is a secondary alias.

- Live analysis requires a submitted token to resolve to exactly one valid interval.
- A token matching one row by ID and another by name is an `identifier_namespace_conflict`.
- A name matching multiple rows is `ambiguous_gene_name`.
- A repeated input row or different aliases resolving to the same canonical ID is a benign list duplicate; the gene enters the source annotation once.
- A duplicate catalog ID is a catalog defect, not a benign list duplicate.
- The resolver never guesses among candidates, strips identifier versions, or silently chooses a row.

For live analysis, structural catalog defects abort immediately. A referenced row-level defect follows the selected resolution policy. An unreferenced row-level defect has no effect on the current scientific annotation, so it is logged and the live run may continue. Index construction is stricter: every structural or row-level defect anywhere in the catalog is fatal before atom construction.

## Gene-list source contract

A focal or control source is plain text or gzip-compressed text, is headerless, and contains one identifier per nonblank line.

- Leading and trailing whitespace is removed.
- Blank lines are counted but not audited.
- There is no comment syntax; a nonblank line beginning with `#` is an identifier.
- Multiple tab-separated fields are `malformed_input` and remain fatal under both resolution policies.
- All readable sources are scanned before Gate A reports problems.

A focal query name is its basename after removing an optional final `.gz` and then at most one final `.txt`, `.tsv`, or `.list`, case-insensitively. Names must be globally unique. A collision is reported with all other Gate A source and identifier issues rather than resolved by implicit renaming.

## Resolution policies and validation gates

### Strict and resolved-only policies

Under `strict`, every rejected gene-list row aborts Gate A.

Under `resolved-only`, these row-level failures may be omitted:

```text
unmatched_identifier
ambiguous_gene_name
identifier_namespace_conflict
catalog_duplicate_id
catalog_invalid_coordinates
outside_supported_chromosome
```

The policy does not relax structural catalog/build failures, missing columns, unparseable inputs, unreadable or incorrectly encoded/compressed sources, malformed multi-field list rows, duplicate query names, corrupt/partial production indexes, or incompatible CLI modes.

The audit outcome is policy-independent: an unusable row remains `disposition=rejected`. The policy controls whether rejected rows stop the run or are explicitly omitted. There is no arbitrary minimum resolution fraction.

### Gate A: catalog and identifier preflight

Before substantial LD-score work, Gate A validates all focal sources and the control together. The control cannot raise before focal issues are collected. Every safely discoverable source, catalog, and row problem is reported in one submission.

An unreadable source has no synthetic row-level audit record because its rows are unknowable. It does receive a source-summary row with blank counts and a source error. A structurally unusable catalog produces no gene-list audit because there is no reliable resolution authority.

In `resolved-only`, a focal source with no nonblank rows is skipped as `empty_gene_list`; a nonempty focal source with no usable genes is skipped as `zero_resolved_genes`. A partial control may continue, but an empty or zero-resolved control is fatal.

### Gate B: retained-SNP support

After the retained reference-panel SNP universe is loaded, but before LD-score calculation, Gate B evaluates all resolved, non-excluded genes and focal queries together.

- Per-gene zero retained-SNP support is nonfatal and audited.
- A focal query with zero annotation SNPs is skipped as `zero_annotation_snps`.
- A control with zero annotation SNPs is fatal.
- Every knowable Gate B condition is reported across the batch in one bounded console message and complete diagnostics.

Zero LD-score variance is evaluated only after scores exist. All zero-variance focal queries are pruned together as `zero_variance_ld_scores`; a zero-variance control is fatal. If every focal query is skipped across all stages, the run fails after writing available diagnostics.

### Explicit gene exclusion

`--gene-exclude-regions mhc` is an intentional transformation, not failed resolution. It is applied to the unpadded interval before padding and applies to focal and control genes. Excluded rows are audited as `excluded`; they do not trigger strict Gate A failure. Direct mode defaults to no gene exclusion. Indexed mode inherits the immutable index policy and rejects a live override.

## Diagnostics contract

The detailed field definitions and repair procedures will live in the dedicated current-behavior document:

```text
docs/current/gene-list-diagnostics-and-repair.md
```

That document is the single detailed user reference for interpreting the following artifacts and curating gene lists or the coordinate transformation.

### Gene-list audit

Every readable nonblank focal/control input row appears exactly once in:

```text
diagnostics/gene_list_audit.tsv.gz
```

Schema:

```text
argument
input_role
query
source
source_ordinal
line
input_gene
match_type
canonical_gene_id
catalog_lines
disposition
reason
chrom
start
end
details
```

`match_type` is `gene_id`, `gene_name`, `ambiguous`, or `unmatched`. Coordinates remain one-based. `catalog_lines` is one physical line for a unique match, a comma-separated ascending set for conflicting matches, or blank when unmatched. Ambiguous `details` lists conflicting gene IDs.

Disposition values are:

| Disposition | Meaning |
| --- | --- |
| `retained` | First usable occurrence; selected if the run proceeds. |
| `duplicate` | Later occurrence of the same canonical gene within that source. |
| `excluded` | Removed by explicit gene-region policy. |
| `unsupported` | Resolved gene with zero retained-SNP support. |
| `rejected` | No usable unique interval; strict stops, resolved-only may omit. |

Stable reasons are:

```text
duplicate_canonical_gene
excluded_gene_region
zero_reference_snp_support
unmatched_identifier
ambiguous_gene_name
identifier_namespace_conflict
catalog_duplicate_id
catalog_invalid_coordinates
outside_supported_chromosome
outside_index_chromosome_coverage
malformed_input
```

Disposition precedence within a source is rejected first; otherwise the first canonical occurrence is retained, excluded, or unsupported, and later occurrences are duplicates that point to the first line. Deduplication resets for every focal source and the control.

Audit order is deterministic: focal before control, focal `source_ordinal` ascending, input line ascending, and the control block last despite ordinal `0`. Ordering never depends on dictionary, filesystem, resolver-construction, or worker order. Gzip byte identity is not required because compression headers may contain timestamps.

### Per-source resolution summary

Every declared focal/control source appears once in:

```text
diagnostics/gene_list_resolution_summary.tsv
```

Schema:

```text
argument
input_role
query
source
source_ordinal
source_status
source_reasons
resolution_policy
nonblank_input_rows
uniquely_resolved_rows
rejected_rows
unique_resolved_genes
duplicate_rows
excluded_genes
zero_support_genes
genes_with_snp_support
resolution_fraction
```

`source_status` is `ok` or `error`. Stable source reasons are `unreadable_gene_list`, `invalid_gzip`, `invalid_utf8`, and `duplicate_query_name`; multiple reasons are comma-separated deterministically. Unreadable-source counts are blank.

Count definitions:

| Field | Definition |
| --- | --- |
| `nonblank_input_rows` | All nonblank rows, including rejected and duplicate rows. |
| `uniquely_resolved_rows` | Rows matching exactly one valid catalog gene before deduplication, exclusion, or SNP-support evaluation. |
| `rejected_rows` | Rows with `disposition=rejected`. |
| `unique_resolved_genes` | Distinct canonical IDs before explicit exclusion. |
| `duplicate_rows` | Later resolved occurrences of a canonical ID within the source. |
| `excluded_genes` | Distinct genes removed by explicit gene exclusion. |
| `zero_support_genes` | Distinct non-excluded genes with zero retained-SNP support. |
| `genes_with_snp_support` | Distinct non-excluded genes with at least one retained reference SNP. |
| `resolution_fraction` | `uniquely_resolved_rows / nonblank_input_rows`; blank for empty/unreadable sources. |

After Gate B:

```text
excluded_genes + zero_support_genes + genes_with_snp_support
    = unique_resolved_genes
```

If Gate B was not reached, SNP-support fields are blank rather than zero.

### Query annotation status

Gate A hard failures do not write `query_annotation_status.tsv`, because scientific viability was not evaluated. Once evaluation proceeds, the file contains one row per focal source with these final states:

| Status | Reason |
| --- | --- |
| `ok` | empty |
| `warning` | `partial_gene_resolution` |
| `warning` | `partial_snp_support` |
| `skipped` | `empty_gene_list` |
| `skipped` | `zero_resolved_genes` |
| `skipped` | `zero_annotation_snps` |
| `skipped` | `zero_variance_ld_scores` |

When partial gene resolution and partial SNP support both apply, `partial_gene_resolution` is the primary reason and `details` plus the resolution summary preserve both. The control has no query-status row.

### Console, log, and failure publication

Console output is bounded and actionable; logs and diagnostics are complete.

- `MAX_CONSOLE_GENE_ISSUES = 10` applies across the ordered batch, not per source.
- A strict hard stop reports rejected/submitted totals, per-source counts, the first 10 `(line, input_gene, reason)` records, omitted count, diagnostic paths, a source/role filter expression, and `--overwrite` guidance.
- A successful `resolved-only` run that omitted rows emits one bounded warning with aggregate resolved/total counts, the first 10 affected sources, and audit/summary paths.
- The complete log records every rejected gene with role and source.

A strict preflight failure publishes only available diagnostics under `<output-dir>/diagnostics/`; it publishes no root metadata or scientific Parquet artifact. A corrected rerun into the owned output directory requires `--overwrite`.

## Index catalog contract

The index embeds the complete canonical catalog used during construction, including intentionally region-excluded genes. The public embedded catalog preserves one-based `start` and `end`, original physical `catalog_line`, and source basename. Internal atoms use zero-based half-open coordinates.

The public builder:

- requires `--gene-coordinate-file`;
- always builds chromosomes 1 through 22;
- removes `--chromosomes` from the public parser, help, and documentation;
- validates every catalog row before atom work;
- never omits a defective row merely to finish construction; and
- publishes no index when catalog validation fails.

Partial coverage remains only through a private test seam. Internal partial builds still validate/embed the whole source catalog and distinguish out-of-coverage genes from zero-SNP-support genes. Public loading rejects any partial index before list resolution.

The indexed loader revalidates the complete embedded-catalog invariant through the same shared validator; it must not trust a weaker frame constructor. A corrupt, incomplete, ambiguous, partial, or incompatible index fails before query/control resolution. Indexed mode never consults an installed package catalog or accepts a live catalog override.

### Builder catalog-issues artifact

Catalog validation failure writes every detected defect to:

```text
.<index-name>.build-state/gene_coordinate_catalog_issues.tsv.gz
```

There is one row per defect occurrence; a physical catalog line may therefore appear more than once. Schema:

```text
source
catalog_line
gene_id
gene_name
chrom
start
end
genome_build
field
reason
related_catalog_lines
observed_value
details
repair
```

Duplicate groups produce an issue for each implicated line and list the complete ascending line group. Missing columns produce one issue per column at the header line. Raw values are preserved for repair. Rows sort by physical line, field, and reason.

Stable reasons are:

```text
catalog_unreadable
catalog_unparseable
malformed_catalog_row
missing_required_column
missing_genome_build
conflicting_genome_build
unsupported_genome_build
missing_gene_id
invalid_chromosome
outside_supported_chromosome
invalid_start
invalid_end
end_before_start
duplicate_gene_id
duplicate_gene_name
identifier_namespace_conflict
validation_incomplete
```

Validation collects every safely detectable schema, row, and global-canonicality defect. If a structural defect prevents later phases, `validation_incomplete` names the unevaluated phases and warns that more defects may appear after repair.

The fixed path represents the current failed attempt. Before another attempt, the prior issues file is archived with the corresponding historical builder log. A successful rebuild leaves no current issues file and publishes no empty issues file inside the index.

## Provenance and compatibility

Result metadata records:

- `gene_list_resolution_policy`;
- aggregate resolution counts;
- relative paths to the audit, resolution summary, and query-status artifact when applicable; and
- the index identity for indexed runs.

Gene-list and live catalog content SHA-256 values are removed from routine provenance. The catalog does not expose a separate checksum. The index retains its top-level `index_id` because it defines immutable artifact identity; normalized embedded catalog content participates directly in that identity calculation.

Implementation removes the packaged protein-coding catalog, all package constants/loaders for it, associated tests and documentation, gene-list input hashes, and catalog provenance hashes. There is no fallback, deprecation shim, migration layer, or compatibility support for old gene-index directories. Old indexes must be rebuilt.

## Documentation deliverables

Documentation changes ship with implementation, not before behavior exists.

1. Create `docs/current/gene-list-diagnostics-and-repair.md` as the detailed source of truth for:
   - every audit, summary, query-status, and builder-issues column;
   - every disposition, reason, status, count, and null-value rule;
   - strict versus resolved-only interpretation;
   - commands or filters for isolating one source and rejected rows;
   - deterministic diffing across reruns; and
   - detailed repair steps for user-owned lists and catalog transformations.
2. The repair workflow must prioritize editing focal/control lists: identify affected sources from the summary, filter audit rejections, replace ambiguous names with authoritative IDs, correct/remove unmatched tokens, and rerun. Catalog editing is secondary and is recommended only when `catalog_lines`, `details`, or builder issues demonstrate a faulty authoritative transformation.
3. The guide must explain how to repair/rebuild a coordinate catalog without changing its declared build, one-based coordinates, unique IDs, or authoritative source semantics. It must explain Gate B support and query-skip outcomes separately from identifier resolution.
4. `docs/troubleshooting.md` contains a concise entry and link to the dedicated guide rather than duplicating it.
5. Relevant `docs/current/` LD-score input/index documents define the four modes and padding matrix and link to the guide.
6. Relevant wiki pages remain task-oriented and high level. They explain the strict default, the exploratory policy, and the basic idea of inspecting the summary then fixing/filtering a list; they link to the detailed current guide instead of duplicating schemas, reason catalogs, or step-by-step repair details.

## Implementation and performance constraints

Gene-list resolution uses bulk tabular operations rather than Python-level iteration over genes.

1. Parse all readable focal/control sources into one table containing role, source, ordinal, and physical line.
2. Resolve identifiers through vectorized joins/indexed maps, masks, and group operations against the normalized catalog.
3. Derive ambiguity, duplication, exclusion, disposition, and summaries in bulk.
4. Sort only from declared deterministic fields.
5. Aggregate per-source summaries through vectorized grouping.

Bounded iteration over source files and per-chromosome/index work is acceptable. A per-gene Python resolver loop is not. Vectorization must not weaken exhaustive validation, actionable diagnostics, exact resolution, or deterministic ordering.

Direct resolution must execute inside the workflow output/logging boundary so Gate A failures can publish their owned diagnostics. Catalog and retained-panel inputs should each be read once per invocation; control and focal validation share one catalog pass and one Gate B panel pass.

## Validation strategy

Automated tests must cover:

1. Catalog parsing for TSV/TSV.GZ, flexible column order, comments/blanks, whitespace/build/chromosome normalization, exact identity preservation, invalid integers, duplicate IDs/names, namespace conflicts, invalid intervals, unsupported chromosomes, structural failures, and exhaustive deterministic builder issues.
2. Coordinate conversion at base 1, both gene boundaries, and padding boundaries, plus equality between gene-derived and equivalent BED-derived annotations.
3. Vectorized mixed focal/control resolution for IDs, names, unmatched and ambiguous tokens, referenced catalog defects, duplicates/aliases, MHC exclusion, strict failure, and resolved-only continuation.
4. Gate A batching across focal and control lists, including unreadable sources and duplicate query names, without the current early standalone control failure.
5. Gate B batching for per-gene zero support, zero-annotation focal skips, fatal zero-control annotations, and blank support counts when Gate B was not reached.
6. Post-computation batching for zero-variance focal skips, fatal zero-variance control, mixed usable/skipped batches, and all-focal-skipped failure.
7. Exact audit/summary/query-status schemas, ordering, counts, null semantics, bounded console samples, filter guidance, output ownership, and `--overwrite` reruns.
8. All four CLI modes, mutual exclusions, padding matrix, policy validity, required live catalog, and indexed inheritance.
9. Index build failure before atom construction, full issue-file lifecycle, embedded-catalog validation, full public chromosome coverage, internal partial-index seams, and rejection of old/corrupt/partial public indexes.
10. Direct/indexed equivalence for matching immutable inputs and successful downstream partitioned-LDSC consumption.
11. Removal of packaged-catalog resources, constants, hash fields, old diagnostic literals, stale documentation, and backward-compatibility paths.
12. Representative many-pathway workloads to confirm the resolution path uses vectorized table operations and avoids per-gene Python scaling.

Stable fixtures should include a small one-based catalog with unique IDs, duplicate IDs/names, namespace collisions, invalid coordinates, MHC overlap, unsupported chromosomes, and genes with/without retained SNP support. Expected audit and summary tables should be checked as complete tables, not only as individual exception messages.

## Scientific and operational constraints

- Gene intervals are one-based inclusive at input and converted once to zero-based half-open internal intervals.
- The catalog declares one build; the workflow performs no liftover.
- Gene identity is exact and case-sensitive after surrounding-whitespace cleanup.
- MHC exclusion precedes padding.
- Zero retained-SNP support is distinct from failed catalog resolution.
- A requested control defines the conditioning model and cannot disappear silently.
- Batch reporting is exhaustive at each gate to avoid repeated SLURM queue submissions that reveal one issue at a time.
- Audit row order is reproducible; compressed bytes need not be identical because gzip headers may contain timestamps.
- The user remains responsible for generating a catalog from one internally consistent authoritative genome build.

## Out of scope

- fuzzy, case-insensitive, synonym, or historical-name matching;
- automatic identifier-version stripping or identity inference;
- implicit liftover or row-level mixed-build detection;
- transcript-, exon-, strand-, TSS-, or TES-specific annotation semantics;
- non-autosomal gene LD-score analysis;
- partial production indexes;
- automatic catalog augmentation or fallback to package data;
- preserving compatibility with the packaged-catalog workflow or old indexes;
- changing ordinary BED or prebuilt annotation scientific semantics beyond the approved mode/padding validation; and
- imposing an arbitrary minimum resolution fraction under `resolved-only`.

## Risks and open questions

There are no unresolved product, scientific, or architectural decisions blocking implementation.

Accepted risks:

- A plausible catalog whose rows secretly mix coordinate builds cannot be detected cheaply when the declared build is internally consistent. This limitation must be documented.
- `resolved-only` deliberately changes a submitted gene set. Its explicit CLI policy, result metadata, per-source summary, row audit, and warning are the provenance boundary.
- Very wide pathway batches may still make query LD-score output or downstream partitioned regression, rather than identifier resolution, the dominant memory cost.
