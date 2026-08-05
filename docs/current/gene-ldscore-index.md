# Exact gene LD-score indexes

Last updated on: 2026-08-05

An exact gene LD-score index moves the repeated PLINK calculation for one
baseline, reference panel, regression-row policy, and gene projection offline.
Online gene-list runs resolve genes against the embedded catalog, assemble
Boolean unions from stored disjoint atoms, and write an ordinary canonical
LD-score directory.

For the complete input-to-output derivation—including the builder factorization,
online matrix products, annotation counts, and overlap sufficient statistics—see
[Exact gene LD-score index: mathematical algorithm](gene-ldscore-index-mathematics.md).

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

The builder is hg19/PLINK-only and supports the base `rsid` and `chr_pos`
identity modes. Both `--genome-build hg19` and `--snp-identifier` are required;
there is no default, `auto`, inference, or liftover. The default persisted regression/output
rows are bundled HapMap3 SNPs with MHC and centromere regions removed. Supply
`--regression-snps-file custom.snplist` to replace the HapMap3 candidate set;
`--exclude-regions` is still applied afterward. Its choices are `none`, `mhc`,
`centromeres`, and `mhc-and-centromeres`.

The LD-reference universe is the inner intersection of baseline and PLINK SNPs
under the configured identifier mode, matching ordinary PLINK `ldscore`
behavior. Every row in a duplicate effective-identity group is dropped, warned,
and recorded under `diagnostics/dropped_snps/`; no representative is selected.
Baseline-only and PLINK-only SNPs are dropped and counted.
PLINK supplies chromosome, position, alleles, cM, genotypes, output identity,
and gene projection coordinates. A matched rsID with different baseline and
PLINK coordinates is warned and counted, but PLINK coordinates win. This is
scientifically valid only when rsID is the intended identity contract and the
caller has independently verified that the PLINK panel is hg19; rsID matching
cannot prove genome build.

### SNP identity contract

`rsid` uses `SNP` as the effective join key. `chr_pos` uses exactly normalized
`(CHR, POS)` with positive 1-based positions; a baseline `SNP` column is
optional and never participates in coordinate matching. After either join,
PLINK supplies the published `CHR`, `POS`, `SNP`, `A1`, and `A2`. Thus differing
labels at one matched coordinate are allowed and reported, and the PLINK label
is authoritative. `rsid_allele_aware` and `chr_pos_allele_aware` are not builder
choices.

The completed `index_id` binds the literal mode and explicit hg19 build. Online
indexed assembly does not rematch sources and rejects live `--snp-identifier`
and `--genome-build` options even when they equal the index. A matching direct
run must explicitly use the index's mode and hg19 build.

An individual chromosome may have zero regression rows after restriction and
region subtraction, which produces a warning. The complete build fails if the
selected chromosome set has zero regression rows in aggregate.

### Baseline LD-score contributor caveat

Protein-coding gene regions do not restrict recomputation of the supplied
baseline LD scores. For baseline column $c$ and persisted row $j$, the builder
computes

$$
l_c(j)=\sum_{k\in\text{retained reference SNPs within the LD window}}
\widetilde r^2_{jk}A_{kc}.
$$

Thus, all SNPs in the retained baseline/PLINK intersection are eligible
contributors; the value $A_{kc}$ determines a SNP's contribution to column
$c$. Protein-coding intervals define only the disjoint atoms used for
`gene_control` and focal gene-list annotations. They do not filter the supplied
baseline matrix or its LD-reference universe.

### Baseline-suite component usage

The builder consumes annotation values and reconstructs the numerical baseline
payload; it does not import the suite's precomputed LD-score or count files.

| Baseline-suite component | Used by index builder? |
| --- | --- |
| `baseline.N.annot.gz` | **Yes.** Supplies $A$, the baseline annotation matrix. |
| `baseline.N.l2.ldscore.gz` or another precomputed `.ldscore.gz` | **No.** The builder recomputes $L_A=PRA$. |
| `.M` and `.M_5_50` | **No.** Counts are recomputed over the retained reference and common-SNP universes. |
| PLINK BED/BIM/FAM | **Yes.** Supplies genotypes and PLINK-authoritative SNP metadata. |
| Existing regression-weight LD scores | **No.** The builder recomputes $w=PRp$. |

Recomputation guarantees that $L_A$, the atom operator $Y=PRH$, counts,
overlaps, and $w$ share the index's selected individuals, genotype and MAF
filters, reference intersection, map/window, adjusted-$r^2$ implementation,
and regression-row policy. Merely residing beside the annotation shards does
not make a legacy `.ldscore.gz`, `.M`, or `.M_5_50` file an index input.

## Artifact and identity contract

One output directory contains one complete immutable index:

```text
<index-dir>/
    metadata.json
    gene_catalog.parquet
    diagnostics/
        build-gene-ldscore-index.json
        build-gene-ldscore-index.log
        dropped_snps/
            chrN_dropped.tsv.gz
    chromosomes/
        chrN/
            metadata.json
            baseline_rows.parquet
            baseline_statistics.npz
            atoms.parquet
            gene_to_atom.npz
            ldscore_operator.npz
            atom_statistics.npz

.<index-name>.build-state/        # hidden operational state
    build-gene-ldscore-index.lock
    build-gene-ldscore-index.log # present while running or after failure
    history/                     # prior failed/interrupted attempts

.<index-name>.stage-<run-id>/     # hidden, private transaction while building
    .gene-index-publication.json
    <index-name>/chromosomes/chrN/
```

`index_id` is a canonical SHA-256 identity over scientific content and settings,
not file paths or output names. It covers normalized baseline data, PLINK BED
content and BIM metadata, selected IID order, canonical regression keys and
region policy, genetic map, gene catalog, chromosome coverage, build/identifier,
window and MAF rules, padding, and gene-region policy. Batching, threads,
logging, output path, and overwrite are resource or publication controls and do
not change the identity.

Root and chromosome metadata repeat the immutable identity mode and hg19 build.
Each chromosome also records ordered effective-identity and published-row
metadata digests. New PLINK-backed indexes require both allele columns and
persist canonical `CHR SNP POS A1 A2` metadata for regression rows only. The
strict loader rejects mode/build disagreement, duplicate effective identities,
row tampering, and older gene-index metadata contracts before indexed output is
published; rebuild older gene indexes with the current builder.

Gene indexes do not support incremental updates, chromosome append, profile
addition, or common-layer reuse. Any input, configuration, or coverage change
requires a complete new build. A chromosome-22 prototype therefore belongs in
a different output directory from a chromosomes-1–22 production index. The
incremental appearance of chromosome shards in a private run stage is only a
memory and durability strategy; it is not restart, resume, checkpoint reuse,
or incremental index-update support.

## Output preflight, logging, and replacement

- A missing destination stays absent, and an existing empty destination stays
  empty, until successful publication.
- A valid existing index requires `--overwrite`, even if its `index_id` would
  be unchanged.
- A nonempty invalid directory fails before chromosome computation, including
  with `--overwrite`.
- A failed build never creates a partial or diagnostics-only public index.
  Graceful failure removes its marked private transaction best-effort; an
  interrupted transaction is never reused scientifically on retry.
- Two builders cannot target the same absolute directory concurrently.

While the build is running, the live log is
`<parent>/.<index-name>.build-state/build-gene-ldscore-index.log`, so `tail -f`
works without placing an open file inside the replaceable artifact. After a
successful publication, the handler closes and the complete log moves atomically
to `<index-dir>/diagnostics/build-gene-ldscore-index.log`. Failed logs stay in
hidden state, and a retry archives them under `.<index-name>.build-state/history/`.
The lifecycle footer is the status authority: `Started`
without a terminal `Finished` or `Failed` line indicates abrupt termination;
there is no separate status file. The JSON diagnostic is a successful-build
summary, not a live status record.

Before chromosome computation, the builder creates one marked run-specific
sibling stage. Each worker writes its large payloads in a private temporary
shard directory, closes them, and atomically renames that directory to
`<stage>/<index-name>/chromosomes/chrN`. `Finished chromosome N` is logged only
after this rename, so it means the internal shard is durable; it does not mean
chromosome N is public. Workers return compact evidence, and the coordinator
writes shared metadata, the catalog, diagnostics, and final `index_id` metadata
after all chromosomes finish in canonical order.

The complete staged index is then reload-validated and moved into place without
rewriting or copying its chromosome payloads. During overwrite, the old valid
index remains loadable until the replacement passes validation. A graceful
failure keeps the old scientific index and its prior success diagnostics while
leaving the failed attempt in hidden build state. After the replacement is
reload-validated, transaction cleanup is best-effort garbage collection: a
cleanup error warns with the retained builder-owned path but does not turn the
completed publication into a failed command. Recognized computation-only
stages are discarded on the next invocation, never resumed. Ambiguous backup
evidence fails rather than guessing.

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

The assembled LD-score directory is self-contained and records the inherited
mode/build plus `index_id`; downstream `h2`, `rg`, and `partitioned-h2` use the
ordinary summary-statistics identity, genome-build, and explicit downgrade
rules. The source index is not needed after successful assembly.

**`gene_control` is a binary baseline annotation consists of the selected (padded) control genes. Its default is `all-protein-coding` (all eligible protein-coding genes); use `none` to omit it or provide one custom control-list path.** It captures broad protein-coding-gene background enrichment; the focal gene-list regression coefficient (`tau`) represents the incremental per-SNP heritability beyond that background and the other baseline annotations, not a standalone causal effect. Overlapping, nested, duplicated, and alias-selected genes use Boolean union. The stored operator is float64, preserves negative adjusted-r-squared values, and is not clamped or epsilon-pruned.

The broader intersected PLINK/baseline universe supplies LD-score contributors,
counts, and overlaps. Only the configured regression restriction and region
policy select persisted rows and `regression_ld_scores` contributors.

See also the task-oriented [build guide](../wiki/utility-functionalities/build-gene-ldscore-index.md)
and [indexed LD-score guide](../wiki/main-functionalities/ldscore.md).
