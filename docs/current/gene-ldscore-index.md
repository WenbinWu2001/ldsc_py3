# Exact gene LD-score indexes

Last updated on: 2026-08-03

An exact gene LD-score index moves the PLINK calculation for a fixed baseline,
reference panel, and gene profile offline. Online gene-list resolution still
uses Boolean interval union, but assembles LD scores from stored disjoint atoms
and writes the ordinary self-contained canonical LD-score directory.

## Why gene-set indexes matter

Many tissue, cell-type, proteomic, GO, and SynGO hypotheses are naturally gene
sets. LDSC-SEG showed how to test such hypotheses with genome-wide polygenic
signal: construct an annotation around specifically expressed genes and fit it
with S-LDSC conditional on the baseline model and an all-genes annotation
([Finucane et al., 2018](https://doi.org/10.1038/s41588-018-0081-4)). The index
preserves that conditional model while moving the repeated reference-panel LD
calculation offline. It accelerates repeated gene-set LD-score construction; it
does not change the downstream S-LDSC estimand.

## Build a profile

```bash
ldsc build-gene-ldscore-index \
  --baseline-annot-sources annotations/baseline.@.annot.gz \
  --plink-prefix reference/1000G.EUR.QC.@ \
  --output-dir indexes/1000G_EUR_Phase3_baseline \
  --genome-build hg19 \
  --snp-identifier rsid \
  --ld-wind-cm 1.0 \
  --padding-bp 100000 \
  --gene-exclude-regions mhc
```

V1 is PLINK-only and uses bundled HM3 regression rows minus
`mhc-and-centromeres`. Before genotype QC, every selected baseline shard must
have exactly the same `CHR/POS/SNP` identities as its PLINK BIM; reordering is
allowed, but missing, extra, duplicate, or conflicting rows abort the build.
The initial supported distributed configuration is hg19, rsID,
`1000G_EUR_Phase3`, a 1 cM window, no explicit retained-reference MAF filter,
inclusive `common_maf_min=0.05`, 100 kb padding, and MHC gene exclusion.

V1 performs no liftover and has no independent output-build setting. Baseline
annotation positions, PLINK BIM positions, gene projection, bundled
HM3/region-mask coordinates, and any explicit genetic map must all be hg19.
The published canonical rows inherit those coordinates. Strict baseline/BIM
`CHR/POS/SNP` equality detects disagreement between the two sources but cannot
prove that two mutually matching inputs were labeled with their true build;
source-build verification remains an input precondition. An explicit hg38 map
is rejected, while omitted maps fall back to informative BIM cM values tied to
the hg19 BIM positions.

The suite directory has immutable `common/` data and one or more
`profiles/<profile>/` directories. Semantic `suite_id` and `profile_id` values
bind scientific inputs. Publication stages and reload-validates all components;
`--overwrite` replaces only the targeted profile, preserves sibling profiles,
and reuses common data only under the same suite identity.

## Assemble a gene-list run

```bash
ldsc ldscore \
  --gene-ldscore-index-dir indexes/1000G_EUR_Phase3_baseline/profiles/padding-100000bp-mhc \
  --query-annot-gene-list-sources 'gene_lists/*.txt' \
  --control-gene-list-source all-protein-coding \
  --output-dir results/gene_ldscores
```

The profile path must be explicit. Indexed mode accepts focal gene lists, the
control source, output/overwrite, logging, and threading controls; live
baseline, PLINK/R², build, identity, window, region, padding, map, and SNP-filter
settings are rejected. A missing or corrupt profile never triggers discovery
or direct-mode fallback.

`gene_control` is appended to the baseline block by default. Use `none` to
disable it or pass one custom control-list path. Focal queries remain together
in `ldscore.query.parquet`; they are not batched online. `partitioned-h2`
continues to fit supplied baseline plus `gene_control` plus one focal query.

## Scientific and resource invariants

The index is not HM3-only. V1 fixes **persisted regression rows** to bundled
HM3 minus `mhc-and-centromeres`, while LD-score construction and sufficient
statistics retain the broader reference universe:

- Baseline/query contributors, counts, and overlaps use the broad retained
  PLINK universe, including MHC and centromeric SNPs.
- Persisted rows and `regression_ld_scores` use bundled HM3 minus
  `mhc-and-centromeres`.
- Overlapping, nested, duplicated, and alias-selected genes use Boolean union.
- The operator remains float64, retains negative adjusted-r² values, includes
  the diagonal once, and is neither clamped nor epsilon-pruned.
- SNP and internal atom columns are built in bounded batches. `--threads`
  defaults to one and can multiply chromosome-local memory.

A custom regression-SNP set is outside the v1 indexed domain and requires
direct mode. This restriction must not be described as restricting the PLINK
LD-reference universe to HM3.

The chromosome-22 validation and measured tolerance are recorded in
[the 2026-08-03 audit](../audits/2026-08-03-exact-gene-ldscore-index-chr22.md).

For task-oriented instructions, see:

- [Build an exact gene LD-score index](../wiki/utility-functionalities/build-gene-ldscore-index.md)
- [Calculate LD scores for gene lists with an index](../wiki/main-functionalities/ldscore.md)
