# Create reusable query annotations

Last updated on: 2026-09-10

`ldsc annotate` projects BED intervals or gene lists onto the SNP rows in baseline annotations and writes reusable binary query files. Direct `ldsc ldscore` accepts BED files and gene lists too, so standalone annotation is optional.

## Inputs

Supply `--baseline-annot-sources`, `--output-dir`, and exactly one query route:

- `--query-annot-bed-sources`: user-created or externally obtained BED files. The first three fields are chromosome, 0-based start, and exclusive end. See [BED format](../../current/bed-input-format.md). Omitted padding is zero; explicit nonnegative `--padding-bp` expands both ends and clips starts at zero.
- `--query-annot-gene-list-sources`: one-column lists resolved exactly against `--gene-coordinate-file`. Explicit nonnegative `--padding-bp` is required; `0` selects gene bodies. Use `--gene-list-resolution-policy strict` or `resolved-only`, and `--gene-exclude-regions none` or `mhc`. Standalone annotate has no control-gene-list option.

Baseline rows supply the annotation grid. Exact paths and globs select the actual input contents; an `@` suite explicitly requires autosomes 1–22. Resolved genes outside this scope fail under either resolution policy after explicit exclusions. Catalog and baseline projection builds must agree; rsID identity records projection-build provenance separately.

```bash
ldsc annotate \
  --baseline-annot-sources 'baseline/baseline.@.annot.gz' \
  --query-annot-gene-list-sources 'pathways/*.txt' \
  --gene-coordinate-file genes.hg19.tsv.gz \
  --padding-bp 100000 \
  --gene-exclude-regions none \
  --snp-identifier rsid \
  --genome-build hg19 \
  --output-dir annotations/pathways
```

For BED input, replace the gene-list/catalog/exclusion options with `--query-annot-bed-sources 'beds/*.bed'`; padding remains available.

## Outputs and validation

The output contains `query.<chrom>.annot.gz`, one file per surviving chromosome with all usable query columns. Rows retain the baseline identities and order after global identity cleanup. The canonical header is `CHR BP SNP CM`, optional allele columns, then query columns; `CM` is explicitly `NA`, and generated values are integer `0`/`1`. Ordinary prebuilt quantitative annotations remain supported by `ldscore` and are normalized to float32.

The query files contain no baseline value columns. Reuse the original baseline sources for LD scoring and pass the generated query suite, for example `--query-annot-sources 'annotations/pathways/query.@.annot.gz'`. Use a glob instead of `@` for an intentionally smaller chromosome scope. Do not split each query column into its own chromosome suite.

Gene Gate A reports safely discoverable source, naming, catalog, and identifier issues. Gate B evaluates support only on baseline rows surviving identity cleanup. Empty and globally unsupported focal queries are skipped while usable siblings continue; an all-skipped batch fails. All-one queries remain valid. A query supported on another chromosome retains its zero-valued shards. Standalone annotate performs no reference-panel or regression-design tests.

`diagnostics/` contains metadata, the command log, statuses, applicable input/scope/catalog diagnostics, a complete drop audit, and gene audits/summaries in gene mode. Gene support is labeled `annotation_snp_count`; reference-panel counts are unevaluated. See the [complete condition–outcome table](../../current/annotate-gene-list-decisions.md#validation-gates-and-conditionoutcome-table).

## Memory and Python ownership

Preparation scans whole-genome files in bounded chunks. Projection and writing proceed chromosome by chromosome; the returned bundle references persistent outputs and may depend on original baselines. No whole-genome annotation matrix is returned. New private files stay under the selected output directory. Use `with run_annotate(...) as bundle:` or close a Python bundle after reading it.

This architecture supports large pathway batches, including 1,000 pathways subsequently tested separately against shared baseline categories. Standalone annotate is sequential; direct LD scoring additionally supports bounded chromosome workers and query batching. Released arrays become reusable memory, but process RSS may not immediately decrease. See the [developer memory design](../../current/annotation-memory-design.md) and [LD-score guide](../main-functionalities/ldscore.md). Existing artifacts require `--overwrite`; failure markers and stale-output cleanup retain their command contracts.
