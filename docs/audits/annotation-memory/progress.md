# Annotation memory implementation evidence

Last updated on: 2026-09-10

The [confirmed plan](../../plans/2026-09-10-annotation-workflow-memory.md) is reopened after the [completion review](completion-review.md) identified R1–R3 at `8ada67c`. Core storage, projection, index, regression, and quantile changes are implemented; the remaining count, aggregation, and diagnostic requirements prevent full completion. The [developer design](../../current/annotation-memory-design.md) describes the final architecture; the [resource report](results.md) records dimensions, exact comparisons, peak process-tree RSS, elapsed time, private disk, persistent output size, and limitations separately.

## Verification before the completion review

- Full pytest on behavioral revision `91ad7bf`: **1,494 passed, one skipped, 132 subtests passed**, 189 warnings, **93.16 s**. Command: `python -m pytest -q --tb=short --show-capture=no`.
- Sequential compatibility run: **1,000 unittest tests run, one skipped, OK**, **51.699 s**. Command: `python -m unittest discover -s tests -p 'test*.py' -v`.
- Installed `ldsc --help`, `python -m ldsc --help`, and help for `annotate`, `ldscore`, `build-gene-ldscore-index`, `partitioned-h2`, and `quantile-h2` succeeded.
- All **21** measured current artifact sets matched their appropriate baseline or matching BED inputs. LD/count/overlap/membership/regression differences were zero; maximum quantile difference was **1.5210055437364645e-13**, within unchanged tolerances. Each regression comparison includes 3,002 tables plus model metadata.
- Source and benchmark modules compile. **44** concrete README/wiki/tutorial CLI examples parse; **17** Markdown Python blocks and **16** cells in the updated notebooks compile. The self-contained cell-specific notebook executed API and CLI workflows. The partitioned notebook requires user datasets and was checked without pretending its placeholder paths are executable inputs.
- Current/wiki/tutorial local file and anchor links, modified Markdown dates, cleared notebook outputs, and `git diff --check` were checked. No formatter, type checker, or documentation builder is configured. All work was local.

The final audit added independent mixed-allele/global-duplicate coverage, exhaustive direct gene naming failures alongside invalid siblings, persistent standalone return ownership after deleting original gene inputs, and source-summary/selection tuning checks. These are covered by the final suites. All final measured commands clean their owned scratch. Final HM3 LD tables remain aggregate; the new developer design and user guides explain 1,000 separately fitted pathways, default batching of 1000, sequential memory reuse, worker processes, and output-contained scratch.

## Milestone history

The following entries preserve earlier checkpoints and what remained at each stage. Their pending-work statements describe those earlier points, not the final status above.

### Baseline and storage pilot

Before replacing the original annotation workflow, full pytest passed **1,409 tests**, with **one skip**, **132 subtests**, and **156 warnings**, in **82.53 seconds**. Tests used the existing editable installation and `/Users/wenbinwu/miniforge3/envs/ldsc3-dev/bin/python` on the local Mac. Pytest and unittest are run sequentially.

The small resource fixture contains 6,000 SNPs on three chromosomes, two baseline columns, and 1,000 focal query columns. Chromosomes each contain 2,000 SNPs at positions `10 * (1 ... 2000)`; SNP identifiers are globally unique. The baseline columns are one and `row % 3 == 0`; query `q` is `(row + q) % 17 < 4`. Baseline and query values are written as separate whole-genome whitespace/gzip files. Measurement runs in a fresh Python process and excludes fixture generation. This is annotation preparation only, without reference or regression computation.

| Implementation | Peak process RSS | Preparation wall time | Private temporary disk |
| --- | --- | --- | --- |
| Original builder at `a505c45` | 477,347,840 bytes | 0.333 s | 0 bytes for this file-input route |
| Initial bounded preparation pilot | 367,542,272 bytes | 0.991 s | 41,334,521 bytes sampled peak |

Raw records: [baseline](baseline-annotation.json), [storage pilot](storage-pilot.json). The pilot additionally read 32 query columns at a time for every seventh row, taking 0.276 seconds across the three shards. RSS includes imports and runtime overhead; private disk was sampled every 20 ms using logical file sizes and can miss brief peaks. These are single observations, not a runtime guarantee or a production-scale benchmark. Later tuning changed the identity lookup implementation, so final measurements will rerun the completed workflows.

The pilot supports explicit column-major `.npy` reads as the initial private format. Coalesced reads are capped at 65,536 rows and stop before gaps greater than 64 rows. Preparation targets at most 16 MiB of parsed numeric cells across source columns per chunk, capped at 65,536 rows. Detached read results avoid permanent mappings. Exact identity bookkeeping uses a disk SQLite primary-key table with an 8 MiB page cache, adjacent delete journal, and no SQL sorting or temporary B-tree construction. SQLite temporary storage is memory-only; statements never direct scratch to system temporary directories.

### Completed checks

- New storage checks were observed failing before their implementations existed. Independent literal expectations verify row/column ordering, repeated and scattered indices, empty selections, float32 values, and close float64 target values.
- Whole-genome and chromosome-sharded source fixtures produce the same cleaned SNP rows and selected values. Aligned column files count as one logical observation. Cross-chunk/cross-chromosome rsID collisions remove every affected row; coordinate identities remain distinct when their coordinates differ.
- Allele checks preserve invalid-before-multiallelic-before-duplicate precedence. Diagnostic replay follows global policy-reason order and source row order.
- Workspace closure removes only owned private staging, preserves public files, is idempotent, and invalidates bundle reads. A released read array is not retained by the bundle.
- Focused annotation/storage/SNP-identity verification: **110 passed, 27 subtests passed**, in **2.10 seconds**. The original annotation parser now delegates normalization to the shared chunk normalizer; existing annotation behavior remains covered.

### Standalone annotation migration

`run_annotate()` and the real CLI now use staged baseline preparation, bounded BED parsing, staged gene-list resolution, chromosome-local projection, and incremental canonical writing. The old BED-only public wrapper and parser names are retired. One catalog lookup is shared; selections are loaded one source at a time and complete audits are replayed from private chunks. Query projection reuses a chromosome's sorted SNP positions, accumulates compact support counts, and releases that chromosome's working data before advancing. CLI-owned returned resources are closed on exit.

The focused annotation/gene/storage/region/configuration/layout/output/failure-marker check passed **288 tests and 26 subtests**, with 24 existing warnings, in **16.91 seconds**. Added checks cover independently expected interval boundaries, matching BED equivalence in both baseline layouts, support/coverage distinctions, empty and globally unsupported siblings, all-skipped failure, exact source audit ordering, complete missing-suite diagnostics, failure markers, persistent outputs, and a real module CLI subprocess. The legacy direct/indexed paths remain on their original builder until their coherent migration.

Direct/indexed LD batching, regression/quantile changes, complete resource comparisons, and final documentation checks remain pending.

### Direct workflow cutover

The direct command now prepares annotations once into owned chromosome artifacts, reuses that evidence for scope validation, and releases preparation state after each chromosome. Both numerical backends reuse LD blocks across query batches and retain only output SNP score rows. The process scheduler bounds submitted work by its worker count. Gene audit rows stream to canonical diagnostics; returned direct results reference those persistent diagnostics after scratch cleanup. A focused check passed 324 tests and 29 subtests, including real parallel runs, independent numerical expectations, output reloads, overwrite behavior, and a weak-reference lifetime check. End-to-end matched resource measurements remain pending; these correctness checks are not performance measurements.

### Shared annotation, direct LD, and index checkpoint

The index now retains paths and compact support summaries, validates and releases one chromosome at a time, and reloads each operator once for all assembly query batches. Construction uses shared prepared annotation shards inside its preserved transaction. The full pytest suite passed 1,468 tests and 132 subtests, with one skip and 195 warnings, in 85.89 seconds. This is scientific/artifact and lifetime evidence; matched index runtime, peak-memory, and temporary-disk measurements remain pending.

### Downstream and public API integration

Batch regression now prepares shared alignment once, reads explicit query batches, and stages each fit's complete details before releasing them. The focused regression/output integration check passed 242 tests and three subtests. Quantile reconstruction uses exact global float64 targets, disk-backed global identity joins, stable full-common moments, and a second annotation pass for quantile sums; 130 focused quantile/regression/output checks and four subtests passed. These are correctness checks; matched resource measurements remain pending.

The public annotation bundle now contains chromosome descriptors. The builder requires an explicit output workspace; eager fields, obsolete projection helpers, and calculator materialization fallbacks are removed. The migrated annotation/direct/reference/regression/quantile check passed 377 tests and 38 subtests. Eager gene batch resolution was then removed; 94 focused policy/index/output checks passed on staged audits.

Direct chromosome drop tables now stream to private files and become persistent references after publication. Source cleanup rows are partitioned once, including chromosomes with no retained rows, and precede reference-stage diagnostics. Prepared chromosome closure clears arrays and diagnostic references and is idempotent. The integration full pytest run passed 1,488 tests and 132 subtests with one skip, but exposed one retired private-wrapper test that omitted the new diagnostic stage. It was replaced by the real PLINK command test with the same header-only artifact expectation. The subsequent reference/direct/index focused check passed 104 tests. A final full pytest and sequential unittest check remain required after the complete audit and documentation pass.
