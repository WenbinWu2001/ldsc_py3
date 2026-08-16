# Gene-list diagnostics and repair

Last updated on: 2026-08-16

This guide is the detailed reference for diagnosing direct and indexed
gene-list `ldscore` runs. The design assumes that users do not inspect a file
log after a successful SLURM job. Therefore, strict mode stops before LD-score
calculation when identifier resolution would silently change the submitted
gene set. The exploratory `resolved-only` policy continues only with an
explicit audit, summary, metadata record, and bounded console warning.

Start with `diagnostics/gene_list_resolution_summary.tsv`. It identifies the
affected focal or control files and shows their resolved/total counts. Then
filter `diagnostics/gene_list_audit.tsv.gz` to the source and reason you need to
repair. Edit the user-owned gene lists first. Change the coordinate catalog only
when catalog line evidence shows that its upstream transformation is defective.

## Diagnostic files by validation stage

| File | When written | Purpose |
| --- | --- | --- |
| `diagnostics/gene_list_audit.tsv.gz` | Every readable direct/indexed gene-list run that reaches catalog resolution, including Gate A failure | Exactly one row for every nonblank row in every readable focal/control list. |
| `diagnostics/gene_list_resolution_summary.tsv` | With the gene-list audit | Exactly one row for every declared focal/control source, including unreadable sources. |
| `diagnostics/query_annotation_status.tsv` | After Gate A succeeds and scientific viability is evaluated | Final outcome for each focal query. Controls are not query-status rows. |
| `diagnostics/ldscore.log` | Every ordinary `ldscore` workflow | Complete unbounded warnings, including every rejected and zero-support gene. |
| `.<index-name>.build-state/gene_coordinate_catalog_issues.tsv.gz` | A gene-index build whose catalog is noncanonical | Every detected source-catalog defect. It is outside the unpublished index. |

A structural catalog failure may prevent reliable gene-list matching, so it
does not create a synthetic gene-list audit. A Gate A identifier failure writes
the audit and summary but not query status or root scientific artifacts.
After Gate B, a successful CLI invocation prints one bounded notice naming the
affected query outcomes and the three diagnostic paths. The complete gene-level
record remains in the audit and log. Python callers receive the same structured
statuses and output paths in the result object without unsolicited console
output.

## Inspect and filter the files

Read a compressed audit with pandas:

```python
import pandas as pd

audit = pd.read_csv("results/diagnostics/gene_list_audit.tsv.gz", sep="\t")
summary = pd.read_csv("results/diagnostics/gene_list_resolution_summary.tsv", sep="\t")

print(summary.loc[summary["source_status"].ne("ok") | summary["rejected_rows"].fillna(0).gt(0)])
print(audit.loc[(audit["source"] == "set.txt") & (audit["disposition"] == "rejected")])
```

For a shell-only view:

```bash
gzip -dc results/diagnostics/gene_list_audit.tsv.gz | less -S
```

The console error names an audit predicate such as
`source == 'set.txt' & disposition == 'rejected'`. Use that source first. A
corrected rerun into the same owned output directory requires `--overwrite`.
Audit row order is deterministic, so decompressing two runs and diffing them is
a useful way to confirm that only intended rows changed. Compressed bytes need
not be identical because gzip headers can contain timestamps.

## Gene-list audit columns

Coordinates in this file remain one-based inclusive, matching the coordinate
catalog. Blank values mean that the field is not knowable or not applicable.

| Column | Meaning |
| --- | --- |
| `argument` | Originating CLI argument: focal `--query-annot-gene-list-sources` or control `--control-gene-list-file`. |
| `input_role` | `focal` or `control`. Focal rows always precede the control block. |
| `query` | Derived focal query name, or fixed `gene_control` for the control. |
| `source` | Source basename, suitable for filtering without exposing a long path. |
| `source_ordinal` | One-based focal declaration order. The control is `0` but is always sorted last. |
| `line` | Physical one-based line in the list, including intervening blank lines. |
| `input_gene` | Trimmed identifier exactly as submitted; case and version suffixes are preserved. |
| `match_type` | `gene_id`, `gene_name`, `ambiguous`, or `unmatched`. |
| `canonical_gene_id` | Authoritative catalog `gene_id` when one candidate gene is identifiable. |
| `catalog_lines` | One catalog physical line for a unique match; comma-separated ascending lines for conflicting matches; blank if unmatched. |
| `disposition` | Operational row outcome described below. |
| `reason` | Stable machine-readable explanation described below. Blank for an ordinary retained row. |
| `chrom`, `start`, `end` | Unique catalog interval when available, in one-based inclusive coordinates. |
| `details` | Human-readable conflicting IDs, first duplicate line, catalog defect, exclusion, or support context. |

Every readable nonblank input row appears exactly once. Rows sort by focal then
control, focal ordinal, and physical input line; the order never depends on
filesystem enumeration or dictionary construction.

### Dispositions

| Value | Meaning and action |
| --- | --- |
| `retained` | First usable occurrence of this canonical gene in this source. No repair is needed. |
| `duplicate` | A later row in the same source resolves to a gene already selected. The gene remains in the analysis once. Remove the redundant row if desired. |
| `excluded` | An explicit gene-region policy, currently MHC exclusion, intentionally removed the gene. This is not failed resolution. |
| `unsupported` | The gene resolved but overlaps zero retained reference-panel SNPs. This is nonfatal and does not change the annotation on this panel. |
| `rejected` | The row cannot provide one valid usable catalog interval. Strict mode stops; `resolved-only` may omit approved row-level failures. |

### Row reasons and repairs

| Reason | Meaning | Preferred repair |
| --- | --- | --- |
| `duplicate_canonical_gene` | This source already selected the same canonical gene on the line named in `details`. | Keep one authoritative-ID row in the list. |
| `excluded_gene_region` | The unpadded interval overlaps the requested gene-exclusion region. | Keep it if exclusion is intentional; otherwise choose a different explicit policy in direct mode or a different index. |
| `zero_reference_snp_support` | No retained panel SNP falls in the projected interval. | Usually no list repair is necessary. Check panel/build/padding if the result is surprising. |
| `unmatched_identifier` | Neither exact `gene_id` nor exact case-sensitive `gene_name` exists in the catalog. | Correct spelling/version/case, replace with the authoritative ID, or remove the row. Do not assume version suffixes are stripped. |
| `ambiguous_gene_name` | The exact name maps to multiple catalog genes. | Replace the name in the list with the intended authoritative `gene_id`. |
| `identifier_namespace_conflict` | The token is a `gene_id` on one catalog row and a `gene_name` on another. | Replace the list token with an unambiguous authoritative ID; report the catalog transformation conflict upstream. |
| `catalog_duplicate_id` | The authoritative ID occurs on multiple catalog lines. | Prefer removing that row from the list for an exploratory run; fix and rebuild the catalog/index for a durable analysis. Never choose a catalog row arbitrarily. |
| `catalog_invalid_coordinates` | The referenced row lacks a usable one-based autosomal interval. | Remove/replace the list token, or regenerate the catalog from the authoritative annotation. |
| `outside_supported_chromosome` | The referenced row is not on autosomes 1–22. | Remove it; this workflow does not analyze sex chromosomes or noncanonical contigs. |
| `outside_index_chromosome_coverage` | An internal partial test index does not cover the gene chromosome. | Rebuild/use a complete production index. Public indexes must cover 1–22. |
| `malformed_input` | A list row contains tab-separated extra fields rather than one identifier. | Keep exactly one identifier field on that line. This is fatal under both policies. |

## Resolution summary columns

| Column | Meaning |
| --- | --- |
| `argument`, `input_role`, `query`, `source`, `source_ordinal` | Same source identity fields as the audit. |
| `source_status` | `ok` or `error`. An unreadable/invalid source is an error even though no row audit can be produced. |
| `source_reasons` | Deterministic comma-separated source problems: `unreadable_gene_list`, `invalid_gzip`, `invalid_utf8`, or `duplicate_query_name`. An unmatched focal glob and a forbidden control glob are reported as unreadable declarations, with the console message naming the glob condition. |
| `resolution_policy` | `strict` or `resolved-only`. Resolver truth is identical under both policies. |
| `nonblank_input_rows` | Every nonblank submitted row, including duplicates and rejected rows. Blank when the source cannot be read. |
| `uniquely_resolved_rows` | Rows that match exactly one valid catalog gene before within-source deduplication, exclusion, and SNP support. |
| `rejected_rows` | Rows with `disposition=rejected`. |
| `unique_resolved_genes` | Distinct canonical IDs before explicit exclusion. |
| `duplicate_rows` | Later occurrences of a canonical ID within this source. Deduplication resets for each focal/control source. |
| `excluded_genes` | Distinct genes removed by the explicit gene exclusion. |
| `zero_support_genes` | Distinct non-excluded genes with zero retained-SNP support. Blank if Gate B was not reached. |
| `genes_with_snp_support` | Distinct non-excluded genes with at least one retained SNP. Blank if Gate B was not reached. |
| `resolution_fraction` | `uniquely_resolved_rows / nonblank_input_rows`; blank for empty or unreadable sources. |

After Gate B, each readable source satisfies:

```text
excluded_genes + zero_support_genes + genes_with_snp_support
    = unique_resolved_genes
```

Do not interpret blank Gate B counts as zero. They mean that an earlier hard
stop prevented SNP-support evaluation.

## Query annotation status

`query_annotation_status.tsv` has one row per focal source and no control row.
The common identity columns are `query`, `source`, and `input_type`; `status`,
`reason`, `n_annotation_snps`, and `details` describe the final run-specific
outcome.

| Status | Reason | Meaning |
| --- | --- | --- |
| `ok` | blank | The focal annotation is usable without a warning condition. |
| `warning` | `partial_gene_resolution` | `resolved-only` omitted one or more rejected list rows. |
| `warning` | `partial_snp_support` | Some resolved genes have no retained-SNP support, but the query remains usable. |
| `skipped` | `empty_gene_list` | The readable source contains no nonblank rows. |
| `skipped` | `zero_resolved_genes` | A nonempty `resolved-only` source has no usable catalog genes. |
| `skipped` | `zero_annotation_snps` | The resolved annotation has no retained panel SNPs in this run. |
| `skipped` | `zero_variance_ld_scores` | The computed query LD scores are constant on regression rows. |

When both partial gene resolution and partial SNP support occur,
`partial_gene_resolution` remains the primary reason and `details` plus the
summary record the SNP-support condition. Skipping one focal query does not
stop usable siblings. If every focal query is skipped, the command exits
nonzero after writing diagnostics and publishes no baseline-only scientific
result. A requested control instead fails if it has zero resolved genes, zero
annotation SNPs, or zero-variance LD scores because silently dropping it would
change the conditioning model.

## Strict versus resolved-only

Strict is the default and recommended policy for final analyses:

```bash
--gene-list-resolution-policy strict
```

Any rejected focal/control row stops Gate A after the complete batch has been
audited. Use the exploratory policy only when deliberate subset analysis is
appropriate, for example a preliminary screen of hundreds of pathways:

```bash
--gene-list-resolution-policy resolved-only
```

It may omit unmatched, ambiguous, namespace-conflicting, duplicate-ID,
invalid-coordinate, or unsupported-chromosome rows. It never relaxes malformed
or unreadable input, source-name collisions, structural catalog/build defects,
corrupt/partial indexes, or control viability. Successful subset runs emit a
bounded console warning because relying only on a file log would be unsafe on
batch systems.

## Recommended list-first repair workflow

1. Open the resolution summary and identify every source with an error,
   rejected rows, or a low resolution fraction.
2. Filter the audit by `source` and `disposition == 'rejected'`.
3. Replace ambiguous names with authoritative `gene_id` values. Prefer IDs for
   every curated list because names can be shared.
4. Correct case, spelling, or exact version suffixes for unmatched tokens. If
   the intended gene is outside the catalog universe, remove it or choose a
   catalog derived from an annotation that includes that gene.
5. Remove later duplicate rows and malformed extra fields. Duplicates are
   scientifically benign, but a clean list is easier to review and diff.
6. Treat `excluded` and `unsupported` separately from resolution failures.
   Decide whether the explicit exclusion, panel, build, and padding are the
   intended scientific configuration.
7. Rerun with `--overwrite`, then diff the decompressed audit. Confirm that only
   the curated rows changed and that the final query statuses are usable.

For exploratory pathway screening, preserve the original lists and run with
`resolved-only`; archive the summary and audit with the results. Before a final
analysis, curate the source lists and return to strict mode.

## Coordinate-catalog repair and index-builder issues

The coordinate TSV/TSV.GZ is the sole live catalog universe. It must contain
`gene_id`, `gene_name`, `chrom`, `start`, `end`, and `genome_build`. Coordinates
are one-based inclusive: internal projection uses `[start - 1, end)`, subtracting
one from the start only. The catalog represents one declared `hg19`/GRCh37 or
`hg38`/GRCh38 build and autosomes 1–22. Gene IDs must be unique; index building
also requires names and ID/name namespaces to be completely unambiguous.

When index construction fails, inspect
`.<index-name>.build-state/gene_coordinate_catalog_issues.tsv.gz`. Its columns
are:

| Column | Meaning |
| --- | --- |
| `source`, `catalog_line` | Catalog basename and physical source line. |
| `gene_id`, `gene_name`, `chrom`, `start`, `end`, `genome_build` | Raw source values retained for repair. |
| `field` | Field or validation phase implicated by this defect. |
| `reason` | Stable defect code. |
| `related_catalog_lines` | Complete ascending line group for duplicates/conflicts. |
| `observed_value` | Raw value that failed validation. |
| `details` | Additional conflict or incomplete-validation context. |
| `repair` | Concrete next action. |

The stable builder reasons are `catalog_unreadable`, `catalog_unparseable`,
`malformed_catalog_row`, `missing_required_column`, `missing_genome_build`,
`conflicting_genome_build`, `unsupported_genome_build`, `missing_gene_id`,
`invalid_chromosome`, `outside_supported_chromosome`, `invalid_start`,
`invalid_end`, `end_before_start`, `duplicate_gene_id`,
`duplicate_gene_name`, `identifier_namespace_conflict`, and
`validation_incomplete`. One source line can have multiple issue rows.
`validation_incomplete` means a structural defect prevented later checks; fix
the named prerequisite and rerun because more defects may then become visible.

Regenerate a defective catalog from its authoritative gene annotation rather
than hand-editing many coordinates. Keep the declared build unchanged, preserve
one-based starts/ends, select one intended gene-level interval per authoritative
ID, and remove ambiguity in IDs, names, and cross-namespace tokens. A successful
rebuild archives the prior issues file with the historical build log and leaves
no empty current issues file inside the published index.

No inexpensive validator can detect catalog rows that secretly mix coordinate
builds while all declarations and intervals remain plausible. Producing a
single-build catalog from one authoritative source remains the user's
responsibility; LDSC performs no implicit liftover.
