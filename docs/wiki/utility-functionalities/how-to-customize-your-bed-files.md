# Prepare query BED files

Last updated on: 2026-10-20

Use one text BED file per query annotation, with at least three tab-separated fields: chromosome, 0-based inclusive start, and exclusive end. Extra columns do not supply annotation weights; LDSC uses binary interval membership. These files are distinct from PLINK binary `.bed` genotype inputs.

```text
chr1	1000	1250
chr1	5000	5400
```

The first interval covers 1-based SNP positions 1001–1250. Use actual intervals in the same build as the baseline projection grid. Filenames determine query names and must not collide with one another or baseline columns. Plain BED and gzip are supported; `--padding-bp` defaults to zero and can explicitly expand intervals.

Pass a quoted glob through `--query-annot-bed-sources` to [direct LD scoring](../main-functionalities/ldscore.md) or [standalone annotate](annotate.md).
`@` is not expanded for BED sources.
For reusable outputs, annotate writes one `query.<chrom>.annot.gz` file containing all usable query columns; in annotation-grid mode it also exports aligned `baseline.<chrom>.annot.gz` files.
Pass the exported baseline and query families together to `ldscore` with `--baseline-annot-sources` and `--query-annot-sources`.
For a BIM grid, use mutually exclusive `--ref-snp-bim-sources` instead of `--ref-snp-annot-sources`; that mode exports generated query columns only.

See the maintained [BED input format](../../current/bed-input-format.md), especially “Skipped Lines,” “Validation,” and “Workflow-Specific Handling,” for accepted headers, malformed-input outcomes, and naming rules.
