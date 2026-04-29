# Minimal External Resource Fixtures

These fixtures are a small, deterministic subset of the local 1KG 30x
chromosome 22 reference-panel resources. They exist so parity tests exercise the
same PLINK, genetic-map, liftover, parquet, and LD-score code paths without
loading the multi-gigabyte `resources/` tree.

Generation:

```bash
python tests/fixtures/generate_minimal_external_resources.py
```

Selection rule:

- take chromosome 22 positions from `src/ldsc/data/hm3_chr_pos_reference.tsv.gz`
- intersect those hg38 positions with `resources/example_1kg_30x/genomes_30x_chr22.bim`
- keep the first 32 matching SNPs in source order
- copy the matching SNP records from the SNP-major PLINK `.bed`
- write small hg19/hg38 genetic-map slices bracketing the selected positions

The fixture intentionally does not vendor the full liftover chain; tests still
use `resources/liftover/hg38ToHg19.over.chain`, which is small compared with
the original PLINK and genetic-map inputs.
