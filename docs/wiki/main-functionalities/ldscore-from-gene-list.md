# Calculate LD scores from gene lists

Last updated on: 2026-09-14

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

### Coordinate catalog

In direct mode, pass a headered TSV or `.tsv.gz` using `--gene-coordinate-file`. All six columns below must exist; column order does not matter.

| Column | Required content |
| --- | --- |
| `gene_id` | Unique, nonempty gene identifier. |
| `gene_name` | Exact gene-name alias; values may be blank, but the column is required. |
| `chrom` | Autosome `1`–`22`; `chr1` also works. |
| `start` | One-based, inclusive integer ≥ 1. |
| `end` | One-based, inclusive integer ≥ `start`. |
| `genome_build` | One consistent build: `hg19`/`GRCh37` or `hg38`/`GRCh38`. |

Illustrative catalog with **fictional coordinates**, separated by tabs:

```tsv
gene_id	gene_name	chrom	start	end	genome_build
GENE001	GENEA	1	100001	110000	hg19
GENE002	GENEB	2	200001	220000	hg19
```

Use actual gene coordinates from one authoritative annotation release/build. These are **not BED coordinates**: both `start` and `end` are one-based and inclusive. Extra columns are ignored. This catalog defines the entire focal and control gene universe; there is no packaged fallback. Indexed mode uses the catalog embedded in the index instead of a separate coordinate file.

### Gene-list files

Pass headerless text files using `--query-annot-gene-list-sources`, with one identifier per line. For example, `pathway_A.txt` selects both genes from the catalog above:

```text
GENE001
GENE002
```

- Plain text and gzip are supported; blank lines are ignored and surrounding whitespace is removed.
- Each identifier must exactly match a catalog `gene_id` or an unambiguous `gene_name`. Prefer authoritative gene IDs because names can be shared.
- Matching is **case-sensitive**. Ensembl version suffixes are **not removed**, and synonyms are not inferred.
- Do not include a header, extra columns, or comments: `#something` is treated as a gene identifier.
- Duplicate entries resolving to the same gene count once per file.
- Each file defines one query gene set. The filename determines its query name: remove a final `.gz`, then at most one final `.txt`, `.tsv`, or `.list` (extensions are case-insensitive). Query names must be unique and cannot collide with baseline columns or the reserved control name `gene_control`.

Focal arguments may use exact paths or quoted glob patterns, such as `"/path/to/gene-sets/*.txt"`. The optional `--control-gene-list-file` uses the same file format but accepts exactly one literal file, not a glob. By default, unresolved or ambiguous identifiers stop the run.

See [Gene-list query input](../../current/gene-list-input-format.md#gene-list-source-format) for the full parsing and resolution contract.

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
applied before padding and is distinct from SNP `--regr-snps-exclude-regions`.

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
schema. Runtime controls `--threads` and `--query-batch-size` are accepted in indexed mode; they do not override the stored scientific configuration.

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

Intentional MHC exclusions are summarized in `diagnostics/ldscore.log` on one INFO line per gene-list source and role, with the count and physical line–gene pairs. An alias includes its canonical ID when different; the full row audit remains unchanged. Direct runs also log the main annotation-preparation milestones before LD-score computation. See [preparation logging](../../current/workflow-logging.md#annotation-preparation) for examples; indexed runs reuse prepared index artifacts and do not scan live annotations.

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
  --output-dir "/path/to/results/partitioned-h2/trait"
```

Each focal result is conditional on the supplied baseline block and optional
`gene_control`. A missing focal output is a workflow status, not evidence for a
biological null; inspect diagnostics before regression.
Query runs write the per-query result tree automatically. To summarize nominal conditional-coefficient evidence across all query gene lists, pass the aggregate partitioned-h2 result root to `ldsc plot`; do not pass an individual query folder.

## Large pathway batches

For 1,000 pathways, regression tests each pathway separately against the same baseline categories. `--query-batch-size` defaults to `1000` in direct/indexed `ldscore` and `partitioned-h2`; lower it to reduce active query workspace. Direct LD scoring writes and releases one query batch before preparing the next, repeating reference/genotype work across batches. Indexed assembly keeps one chromosome operator per worker, writes and releases that chromosome's query batches, then releases the operator. `--threads` defaults to 1, accepts positive counts, `-1` for available cores and `-2` to leave one core free, and is capped at the chromosome count.

The saved baseline and overlap artifacts are shared. One query batch writes `ldscore.query.parquet`; multiple batches write numbered `ldscore.query.batchNNNNN.parquet` files, each with the same genome-wide SNP rows. Root `metadata.json.query_batches` records the ordered files, columns, and chromosome row groups. Python writing workflows return `LDScoreSource`; explicit `read_queries(names)` calls load columns without caching them. Generation and regression batch widths may differ. See the [LD-score guide](ldscore.md#memory-for-many-pathways) for the implementation sources and the [batch regression guide](partitioned-h2.md#testing-enrichment-for-a-large-batch-of-pathways) for result writing.
