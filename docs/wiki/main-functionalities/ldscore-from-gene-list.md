# Calculate LD scores from gene lists

Last updated on: 2026-08-16

Use gene-list mode to turn pathway, expression, proteomic, GO, or other gene
sets into focal annotations for partitioned S-LDSC. LDSC resolves each list
against one explicit gene universe, projects its intervals, and writes the
ordinary self-contained LD-score directory consumed by `partitioned-h2`.

Two backends are available:

- direct mode recomputes LD scores from a live baseline/reference panel and
  requires your build-aware coordinate catalog;
- indexed mode assembles exact results from a curated distributed index and
  uses only that index's embedded catalog/configuration.

## Gene lists and coordinate authority

A gene-list file is headerless plain/gzip text with one exact `gene_id` or
case-sensitive `gene_name` per nonblank line. Prefer authoritative gene IDs:
names can be shared, and LDSC never guesses, strips identifier versions, or
uses fuzzy synonyms.

Focal arguments may use exact paths or deterministic glob patterns. The
control argument names one exact file.

Direct mode requires a headered TSV/TSV.GZ coordinate catalog containing
`gene_id`, `gene_name`, `chrom`, one-based inclusive `start`/`end`, and
`genome_build`. This file defines the entire focal and control gene universe;
there is no packaged fallback. Use a catalog generated from one authoritative
annotation/build.

## Direct mode

```bash
ldsc ldscore \
  --baseline-annot-sources "/path/to/baseline.@.annot.gz" \
  --plink-prefix "/path/to/1000G.EUR.QC." \
  --query-annot-gene-list-sources "/path/to/gene-sets/*.txt" \
  --gene-coordinate-file "/path/to/gene-coordinates.hg19.tsv.gz" \
  --padding-bp 100000 \
  --gene-exclude-regions mhc \
  --gene-list-resolution-policy strict \
  --snp-identifier rsid \
  --genome-build hg19 \
  --ld-wind-cm 1.0 \
  --output-dir "/path/to/results/gene-set-ldscores"
```

Live gene-list padding must be deliberate: pass `0` for gene bodies or a
positive flank. MHC gene exclusion is an intentional, audited transformation
applied before padding and is distinct from SNP `--exclude-regions`.

## Indexed mode

```bash
ldsc ldscore \
  --gene-ldscore-index-dir "/path/to/index" \
  --query-annot-gene-list-sources "/path/to/gene-sets/*.txt" \
  --gene-list-resolution-policy strict \
  --output-dir "/path/to/results/gene-set-ldscores"
```

The index owns its catalog, baseline, panel, build, identity, padding, and
exclusion policies. Do not pass live overrides. Production indexes always
cover autosomes 1–22 and old indexes must be rebuilt for the current catalog
schema.

## Optional control genes

Add one exact file when focal sets should be conditioned on a specific assay or
selection universe:

```bash
--control-gene-list-file /path/to/background-genes.txt
```

It becomes the fixed baseline column `gene_control`. It uses the same catalog,
exclusion, and padding as focal genes. A requested control is never silently
dropped: zero resolved genes, zero annotation SNPs, or zero-variance LD scores
stop the run.

## Strict versus exploratory batches

`strict` is the default for defensible final analyses. LDSC checks every focal
list and the control together; any unresolved or ambiguous row stops before
substantial calculation and reports the complete batch.

For preliminary screens of hundreds of pathways, explicitly use
`--gene-list-resolution-policy resolved-only`. The run uses only resolved genes
and records exact resolved/total counts. Unusable focal queries are skipped
without failing usable siblings. This policy deliberately changes submitted
sets, so the CLI console, result metadata, audit, and summary all identify it.
The CLI also prints a bounded successful-run notice for zero-SNP support or
skipped query outcomes.

## Diagnose and curate

The useful repair loop is:

1. inspect `diagnostics/gene_list_resolution_summary.tsv` to find affected
   focal/control sources;
2. filter `diagnostics/gene_list_audit.tsv.gz` to rejected rows in that source;
3. fix the user-owned list first—replace ambiguous names with authoritative
   IDs, correct/remove unmatched values, and remove malformed/duplicate rows;
4. rerun with `--overwrite` and compare the audit;
5. change/rebuild the catalog only when its physical line evidence shows an
   upstream catalog-transformation defect.

After identifier preflight,
`diagnostics/query_annotation_status.tsv` records focal zero-SNP or
zero-variance skips. Zero-support genes and explicit MHC exclusions are audited
but are not unresolved identifiers.

For every column/reason and step-by-step repair instructions, use the detailed
[Gene-list diagnostics and repair](../../current/gene-list-diagnostics-and-repair.md)
reference. The concise input contract is
[Gene-list query input](../../current/gene-list-input-format.md), and index
construction is covered by
[Build an exact gene LD-score index](../utility-functionalities/build-gene-ldscore-index.md).

## Continue to partitioned S-LDSC

```bash
ldsc partitioned-h2 \
  --sumstats-file "/path/to/trait/sumstats.parquet" \
  --ldscore-dir "/path/to/results/gene-set-ldscores" \
  --output-dir "/path/to/results/partitioned-h2/trait" \
  --write-per-query-results
```

Each focal result is conditional on the supplied baseline block and optional
`gene_control`. A missing focal output is a workflow status, not evidence for a
biological null; inspect diagnostics before regression.
