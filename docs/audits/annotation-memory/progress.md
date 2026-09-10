# Annotation memory implementation evidence

Last updated on: 2026-09-10

Implementation is in progress under the [confirmed plan](../../plans/2026-09-10-annotation-workflow-memory.md). This record reports completed checks, not claims about unfinished workflows.

## Baseline and storage pilot

Before replacing the original annotation workflow, full pytest passed **1,409 tests**, with **one skip**, **132 subtests**, and **156 warnings**, in **82.53 seconds**. Tests used the existing editable installation and `/Users/wenbinwu/miniforge3/envs/ldsc3-dev/bin/python` on the local Mac. Pytest and unittest are run sequentially.

The small resource fixture contains 6,000 SNPs on three chromosomes, two baseline columns, and 1,000 focal query columns. Chromosomes each contain 2,000 SNPs at positions `10 * (1 ... 2000)`; SNP identifiers are globally unique. The baseline columns are one and `row % 3 == 0`; query `q` is `(row + q) % 17 < 4`. Baseline and query values are written as separate whole-genome whitespace/gzip files. Measurement runs in a fresh Python process and excludes fixture generation. This is annotation preparation only, without reference or regression computation.

| Implementation | Peak process RSS | Preparation wall time | Private temporary disk |
| --- | --- | --- | --- |
| Original builder at `a505c45` | 477,347,840 bytes | 0.333 s | 0 bytes for this file-input route |
| Initial bounded preparation pilot | 367,542,272 bytes | 0.991 s | 41,334,521 bytes sampled peak |

Raw records: [baseline](baseline-annotation.json), [storage pilot](storage-pilot.json). The pilot additionally read 32 query columns at a time for every seventh row, taking 0.276 seconds across the three shards. RSS includes imports and runtime overhead; private disk was sampled every 20 ms using logical file sizes and can miss brief peaks. These are single observations, not a runtime guarantee or a production-scale benchmark. Later tuning changed the identity lookup implementation, so final measurements will rerun the completed workflows.

The pilot supports explicit column-major `.npy` reads as the initial private format. Coalesced reads are capped at 65,536 rows and stop before gaps greater than 64 rows. Preparation targets at most 16 MiB of parsed numeric cells across source columns per chunk, capped at 65,536 rows. Detached read results avoid permanent mappings. Exact identity bookkeeping uses a disk SQLite primary-key table with an 8 MiB page cache, adjacent delete journal, and no SQL sorting or temporary B-tree construction. SQLite temporary storage is memory-only; statements never direct scratch to system temporary directories.

## Completed checks

- New storage checks were observed failing before their implementations existed. Independent literal expectations verify row/column ordering, repeated and scattered indices, empty selections, float32 values, and close float64 target values.
- Whole-genome and chromosome-sharded source fixtures produce the same cleaned SNP rows and selected values. Aligned column files count as one logical observation. Cross-chunk/cross-chromosome rsID collisions remove every affected row; coordinate identities remain distinct when their coordinates differ.
- Allele checks preserve invalid-before-multiallelic-before-duplicate precedence. Diagnostic replay follows global policy-reason order and source row order.
- Workspace closure removes only owned private staging, preserves public files, is idempotent, and invalidates bundle reads. A released read array is not retained by the bundle.
- Focused annotation/storage/SNP-identity verification: **110 passed, 27 subtests passed**, in **2.10 seconds**. The original annotation parser now delegates normalization to the shared chunk normalizer; existing annotation behavior remains covered.

Public command migrations, direct/indexed LD batching, regression/quantile changes, complete resource comparisons, and final documentation checks remain pending.
