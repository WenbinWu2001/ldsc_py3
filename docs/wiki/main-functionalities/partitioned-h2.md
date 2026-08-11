# Partitioned heritability

Last updated on: 2026-08-11

`ldsc partitioned-h2` reads a canonical LD-score directory and tests how its
annotations contribute to SNP heritability.

- With baseline annotations only, it reports a functional partitioning model.
- With focal query annotations, it fits one baseline-plus-one-query model per
  focal annotation.
- For gene-list analyses, the baseline block also contains the fixed
  `gene_control` annotation when a control is enabled. Pass
  `--control-gene-list-source none` during LD-score calculation to omit it.

For a complete gene-set workflow—from an explicit complete index through
LD-score calculation and regression—follow [Calculate LD scores for gene lists
with an index](ldscore-from-gene-list.md). For the general end-to-end workflow and output
interpretation, see [the guided tutorial](../guided-tutorial.md) and the
[partitioned LDSC technical reference](../../current/partitioned-ldsc-workflow.md).





three types of query inputs:

- gene list (may need coord table)
- bed files -- how are these usually obtained? created by users or downloaded online? Should this arg to be supported with padding?
- annotations (may deprecate)
