# Partitioned heritability

Last updated on: 2026-08-03

`ldsc partitioned-h2` reads a canonical LD-score directory and tests how its
annotations contribute to SNP heritability.

- With baseline annotations only, it reports a functional partitioning model.
- With focal query annotations, it fits one baseline-plus-one-query model per
  focal annotation.
- For indexed gene-list analyses, the baseline block normally also contains
  the fixed `gene_control` annotation.

For a complete gene-set workflow—from an explicit index profile through
LD-score calculation and regression—follow [Calculate LD scores for gene lists
with an index](ldscore.md). For the general end-to-end workflow and output
interpretation, see [the guided tutorial](../guided-tutorial.md) and the
[partitioned LDSC technical reference](../../current/partitioned-ldsc-workflow.md).
