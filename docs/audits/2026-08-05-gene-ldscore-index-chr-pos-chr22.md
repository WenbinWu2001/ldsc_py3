# Explicit SNP-identity gene LD-score index chromosome-22 gate

Date: 2026-08-05
Last updated on: 2026-08-05

Status: passed locally for `rsid` and `chr_pos` with the Phase 3 baseline suite,
and for `chr_pos` with baselineLD v2.2. The earlier two-suite `rsid`
measurements remain recorded in
[the 2026-08-03 audit](2026-08-03-exact-gene-ldscore-index-chr22.md).

## Contract exercised

Both builders explicitly used `--genome-build hg19` and one of
`--snp-identifier rsid|chr_pos`. Indexed assembly supplied neither option and
inherited both immutable values from the validated index. The source baseline
and PLINK resources were caller-asserted to share hg19; no build inference,
overlap heuristic, liftover, or live indexed override was used.

The runs used chromosome 22, the local `1000G.EUR.QC.22` BED/BIM/FAM suite,
the Phase 3 baseline annotation, bundled HM3 regression candidates minus MHC
and centromeres, a 1 cM LD window, 100 kb gene padding, MHC gene exclusion,
all protein-coding genes as the control, and the overlapping-gene fixture as
the focal query.

## Build measurements

| Identity | Index ID | Wall time | Peak RSS | Retained reference | Regression rows | Atoms | `nnz(Y)` | Payload bytes |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `rsid` | `d6adf7f3a195...` | 109.97 s | 1,299,038,208 | 141,123 | 17,380 | 838 | 614,523 | 10,530,785 |
| `chr_pos` | `4b52ea06a025...` | 152.81 s | 1,523,449,856 | 141,123 | 17,380 | 838 | 614,523 | 10,544,779 |
| `chr_pos`, baselineLD v2.2 | `aba88288b84d...` | 116.52 s | 1,668,268,032 | 141,123 | 17,380 | 838 | 614,523 | 14,989,722 |

Both published indexes passed staged reload and strict validation. Their root
metadata records different full `index_id` values, the literal identity mode,
and `genome_build=hg19`. New PLINK-backed components retain `CHR`, `POS`, `SNP`,
`A1`, and `A2`; only the 17,380 regression rows are persisted. The broad
141,123-row PLINK intersection remains the LD contributor universe.

## Direct/index equivalence

For each mode, direct and indexed outputs had:

- baseline shape `(17380, 60)` and query shape `(17380, 6)`;
- exact pandas equality for baseline and query Parquet tables, including
  identity rows, ordering, column order, and persisted float32 LD scores;
- identical overlap shape `(2971, 4)` and labels;
- maximum absolute overlap difference `0.0`;
- identical gene control and focal query resolution.

The `chr_pos` output records `chr_pos+hg19`. Canonical `rsid` LD-score output
continues the package-wide convention that genome build is not part of public
row identity, while its source index still binds the explicit hg19 assertion in
metadata and `index_id`.

The separate baselineLD v2.2 coordinate run had baseline shape `(17380, 104)`,
query shape `(17380, 6)`, exact baseline/query equality, identical overlap shape
`(9703, 4)` and labels, and maximum absolute overlap difference
`7.654307410120964e-08`. This is inside the established absolute `1e-7`
overlap-only tolerance and exactly reproduces the earlier rsID accumulation
bound for that continuous baseline.

## Downstream gate

A deterministic coordinate-identity summary-statistics table containing all
17,380 rows was consumed successfully with indexed LD scores by `h2` and `rg`.
The `rg` pair completed with status `ok` and 17,380 SNPs. Matching direct and
indexed `partitioned-h2` runs produced byte-identical `partitioned_h2.tsv`
files, including coefficient, standard error, enrichment, p-value, and query
ordering fields.

## Duplicate and alignment evidence

Focused tests independently cover duplicate rsIDs at distinct coordinates,
duplicate coordinates with distinct labels, differing labels at matched
coordinates, empty intersections, repeated restriction keys, baseline-only and
PLINK-only identities, and immutable-artifact duplicate rejection. Mutable
duplicate effective-key groups are removed with vectorized drop-all cleanup,
warnings, and diagnostic rows. A deliberately unsorted BIM fixture verifies
that canonical sorting retains the original raw BIM/BED indices and gathers
the correct genotype columns; the implementation adds no Python per-SNP
alignment loop.

## Repository validation

The final repository-wide pytest run completed 1,206 tests successfully and
skipped one; its only five failures were multiprocessing checks stopped before
worker construction because the managed sandbox denied Python's semaphore-limit
`sysconf` query. Rerunning exactly those five tests outside the sandbox passed,
so all 1,211 executed pytest cases passed in an environment permitting the
required operating-system query. The independent standard-library compatibility
gate ran 962 tests successfully with one skip. Focused gene-index tests, CLI
help/omission checks, package compilation, and `git diff --check` also passed.

## Local output locations

- `rsid` index: `/private/tmp/codex_gene_rsid_20260805_index`
- `rsid` direct/indexed: `/private/tmp/codex_gene_rsid_20260805_direct` and
  `/private/tmp/codex_gene_rsid_20260805_indexed`
- `chr_pos` index: `/private/tmp/codex_gene_chrpos_20260805_index`
- `chr_pos` direct/indexed: `/private/tmp/codex_gene_chrpos_20260805_direct` and
  `/private/tmp/codex_gene_chrpos_20260805_indexed`
- baselineLD v2.2 `chr_pos`: `/private/tmp/codex_gene_chrpos_baselineLD_20260805_index`,
  `/private/tmp/codex_gene_chrpos_baselineLD_20260805_direct`, and
  `/private/tmp/codex_gene_chrpos_baselineLD_20260805_indexed`
- downstream direct/indexed: `/private/tmp/codex_gene_chrpos_h2_direct_20260805`
  and `/private/tmp/codex_gene_chrpos_h2_indexed_20260805`
- downstream smoke outputs: `/private/tmp/codex_gene_chrpos_h2_smoke_20260805`
  and `/private/tmp/codex_gene_chrpos_rg_smoke_20260805`

These are local acceptance artifacts, not shipped resources. Whole-genome
construction and publication remain separate operations.
