# Exact gene LD-score indexes

Last updated on: 2026-08-04

An exact gene LD-score index moves the repeated PLINK calculation for one
baseline, reference panel, regression-row policy, and gene projection offline.
Online gene-list runs resolve genes against the embedded catalog, assemble
Boolean unions from stored disjoint atoms, and write an ordinary canonical
LD-score directory.

## Build one complete index

```bash
ldsc build-gene-ldscore-index \
  --baseline-annot-sources annotations/baseline.@.annot.gz \
  --plink-prefix reference/1000G.EUR.QC.@ \
  --output-dir indexes/1000G_EUR_Phase3_baseline_100kb \
  --genome-build hg19 \
  --snp-identifier rsid \
  --ld-wind-cm 1.0 \
  --padding-bp 100000 \
  --gene-exclude-regions mhc \
  --exclude-regions mhc-and-centromeres
```

The builder is hg19/rsID/PLINK-only. The default persisted regression/output
rows are bundled HapMap3 SNPs with MHC and centromere regions removed. Supply
`--regression-snps-file custom.snplist` to replace the HapMap3 candidate set;
`--exclude-regions` is still applied afterward. Its choices are `none`, `mhc`,
`centromeres`, and `mhc-and-centromeres`.

The LD-reference universe is the inner intersection of baseline and PLINK SNPs
under the configured identifier mode (rsID in v1), matching ordinary PLINK
`ldscore` behavior. Duplicate effective identifiers in either source are
ambiguous and fail. Baseline-only and PLINK-only SNPs are dropped and counted.
PLINK supplies chromosome, position, alleles, cM, genotypes, output identity,
and gene projection coordinates. A matched rsID with different baseline and
PLINK coordinates is warned and counted, but PLINK coordinates win. This is
scientifically valid only when rsID is the intended identity contract and the
caller has independently verified that the PLINK panel is hg19; rsID matching
cannot prove genome build.

An individual chromosome may have zero regression rows after restriction and
region subtraction, which produces a warning. The complete build fails if the
selected chromosome set has zero regression rows in aggregate.

## Artifact and identity contract

One output directory contains one complete immutable index:

```text
<index-dir>/
    metadata.json
    gene_catalog.parquet
    diagnostics/
        build-gene-ldscore-index.log
        build-gene-ldscore-index.json
        history/
    chromosomes/
        chrN/
            metadata.json
            baseline_rows.parquet
            baseline_statistics.npz
            atoms.parquet
            gene_to_atom.npz
            ldscore_operator.npz
            atom_statistics.npz
```

`index_id` is a canonical SHA-256 identity over scientific content and settings,
not file paths or output names. It covers normalized baseline data, PLINK BED
content and BIM metadata, selected IID order, canonical regression keys and
region policy, genetic map, gene catalog, chromosome coverage, build/identifier,
window and MAF rules, padding, and gene-region policy. Batching, threads,
logging, output path, and overwrite are resource or publication controls and do
not change the identity.

Gene indexes do not support incremental updates, chromosome append, profile
addition, or common-layer reuse. Any input, configuration, or coverage change
requires a complete new build. A chromosome-22 prototype therefore belongs in
a different output directory from a chromosomes-1–22 production index.

## Output preflight, logging, and replacement

- A missing or empty directory starts a new build.
- A valid existing index requires `--overwrite`, even if its `index_id` would
  be unchanged.
- A nonempty invalid directory fails before chromosome computation, including
  with `--overwrite`.
- A diagnostics-only failed build is a controlled retry state and does not
  require `--overwrite`.
- Two builders cannot target the same absolute directory concurrently.

The stable live log is
`<index-dir>/diagnostics/build-gene-ldscore-index.log`, so `tail -f` works while
chromosomes are running. A retry archives the prior log under
`diagnostics/history/`. The lifecycle footer is the status authority: `Started`
without a terminal `Finished` or `Failed` line indicates abrupt termination;
there is no separate status file. The JSON diagnostic is a successful-build
summary, not a live status record.

Publication writes and reload-validates a complete sibling stage. During
overwrite, the old valid index remains loadable until the replacement passes
validation. A graceful failure keeps the old scientific index and its prior
success JSON while leaving the failed attempt in the stable log. Recognized
interrupted publication transactions are recovered on the next invocation;
ambiguous backup evidence fails rather than guessing.

## Assemble gene-list LD scores

```bash
ldsc ldscore \
  --gene-ldscore-index-dir indexes/1000G_EUR_Phase3_baseline_100kb \
  --query-annot-gene-list-sources 'gene_lists/*.txt' \
  --control-gene-list-source all-protein-coding \
  --output-dir results/gene_ldscores
```

The index directory must be explicit. Indexed mode accepts gene lists, the
control source, output/overwrite, and logging controls; it does not accept live
baseline, PLINK/R2, build, identity, window, map, region, or regression-SNP
overrides. A missing or corrupt index never triggers discovery or direct-mode
fallback.

`gene_control` is appended to the baseline block by default. Use `none` to
disable it or pass one custom control-list path. Overlapping, nested,
duplicated, and alias-selected genes use Boolean union. The stored operator is
float64, preserves negative adjusted-r-squared values, and is not clamped or
epsilon-pruned.

The broader intersected PLINK/baseline universe supplies LD-score contributors,
counts, and overlaps. Only the configured regression restriction and region
policy select persisted rows and `regression_ld_scores` contributors.

See also the task-oriented [build guide](../wiki/utility-functionalities/build-gene-ldscore-index.md)
and [indexed LD-score guide](../wiki/main-functionalities/ldscore.md).
