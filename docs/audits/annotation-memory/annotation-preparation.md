# Annotation preparation optimization

Last updated on: 2026-09-14

## Scope and implementation

This local change covers all four requested preparation opportunities. The starting tree was clean on `restructure` at `6373bdf85dddd73120b89d092dff91f6d4c4ab77`. A small baseline was measured before production edits, and the original Python sources were preserved for matched reruns. No HPC access, job changes, or additional preparation parallelism were involved.

- [`prepare_annotation_sources()`](../../../src/ldsc/_annotation_sources.py) sizes each independent scan from that file's header, then rechunks metadata using the combined width of the sources aligned together after layout discovery. Automatic tiles retain the 16 MiB cell estimate and 65,536-row ceiling. Different chromosomes no longer shrink each other's tiles. Explicit `chunk_rows` continues to apply to both stages.
- `_scan()` stores metadata separately from a raw row-major float32 stream. `_aligned_metadata()` never reads numeric staging. Alignment still combines baseline/query columns into one logical SNP observation before global indexing.
- [`DiskIdentityIndex`](../../../src/ldsc/_annotation_identity.py) uses parameterized primary-key lookups of at most 512 keys, bounded by SQLite's variable limit, and transaction batches of at most 262,144 valid rows. Lookup dictionaries cover only one annotation tile; the global table remains on disk with its 8 MiB page cache. Invalid alleles, unknown alleles, global duplicates, and allele conflicts retain their existing policy.
- `_select_rows()` writes final Parquet metadata and narrow retained-row locators. `_write_values()` uses the retained counts to create the existing column-major NumPy format, then copies each source directly to its final columns in bounded tiles. It avoids aligned numeric concatenations, retained numeric pickles, and their extra read/write pass. Sources are removed as their group finishes. [`FrameSpool`](../../../src/ldsc/_annotation_storage.py) now appends independently serialized frames to one file per spool, with no retained file handles.

The gene-index builder already calls this shared preparation before chromosome workers, so no change to its computation or scheduling code was necessary. See the updated [memory design](../../current/annotation-memory-design.md#bounded-source-preparation) and [builder documentation](../../current/gene-ldscore-index.md#baseline-suite-component-usage).

## Matched workloads and method

The [harness](../../../benchmarks/annotation_preparation.py) generates deterministic gzip annotations with signed, exactly representable fractional values and a duplicate SNP spanning every chromosome. It verifies every retained numeric value, float32 dtype, row and column selection order, duplicate diagnostics, and workspace cleanup against the input formula. Inputs contain no real research data or reference-panel calculations.

| Case | Chromosomes | Rows per chromosome | Baseline columns | Query columns | Layout |
| --- | ---: | ---: | ---: | ---: | --- |
| sharded | 22 | 8,000 | 96 | 0 | Separate chromosome files |
| wide-sharded | 3 | 6,000 | 8 | 1,000 | Separate baseline/query chromosome files |
| wide-whole | 3 | 6,000 | 8 | 1,000 | Two aligned whole-genome files |
| small | 2 | 100 | 4 | 2 | Separate baseline/query chromosome files |

Each measurement used a fresh Python process. Input generation and post-preparation numerical verification are outside elapsed time and the captured lifetime peak RSS. Three matched repetitions alternate before/after ordering; summaries below use their medians. The raw [JSON evidence](annotation-preparation.json) contains every repetition, phase timings, SQL statement counts, open counts, versions, and the SQLite pilots. Initial exploratory runs, including the noncontiguous-write pilot, are excluded from these matched summaries.

The environment was local macOS ARM64, Python 3.13.13, NumPy 2.4.5, pandas 2.3.3, and PyArrow 23.0.1. Measurements include the instrumented preparation and imports in RSS. Logical scratch bytes and simultaneous file counts are sampled every 25 ms under the owned output tree; short-lived peaks may be missed. Python file creation/open counts are also recorded separately; they do not count native SQLite/PyArrow file operations. File-system and SQLite caches were not flushed. These are local, cached synthetic preparation results, not cold-storage or HPC timings.

## Results

Each pair is **before → after**. MiB uses 1,048,576 bytes.

| Case | Median time (s) | Time reduction | Peak RSS (MiB) | Sampled peak scratch (MiB) | Sampled peak files |
| --- | ---: | ---: | ---: | ---: | ---: |
| sharded | 6.800 → 2.752 | 59.5% | 271.7 → 264.9 | 142.929 → 76.766 | 417 → 88 |
| wide-sharded | 3.915 → 2.160 | 44.8% | 274.4 → 357.9 | 139.019 → 93.736 | 84 → 20 |
| wide-whole | 2.521 → 2.130 | 15.5% | 493.4 → 375.2 | 134.867 → 139.871 | 32 → 12 |
| small | 0.049 → 0.043 | 12.9% | 114.6 → 114.2 | 0.027 → 0.030 | 6 → 11 |

Runtime ranges were 6.781–7.900 → 2.721–2.858 s for sharded, 3.744–4.172 → 2.104–2.720 s for wide-sharded, 2.425–2.642 → 2.069–2.194 s for wide-whole, and 0.046–0.076 → 0.040–0.043 s for small. The small result is only a control; millisecond differences do not establish a useful speedup.

The exact count of distinct scratch files created through Python fell from 440 to 89 for sharded, 87 to 19 for wide-sharded, and 35 to 9 for wide-whole; it increased from 10 to 13 for small. A spool's file count is now constant in its number of appended chunks, but separate metadata, numeric, and selection streams add fixed overhead for tiny inputs.

The sharded identity tile increased from at most 962 rows to one 8,000-row tile per chromosome; its computed ceiling is 21,183 rows. Wide-sharded identity tiles increased from at most 689 to 2,068 rows, matching the combined width of each baseline/query pair. Wide-whole retained the same 2,068-row identity ceiling, while its two files are scanned independently. Both declared and discovered chromosome layouts are covered by tests.

### Phase evidence

Times are medians in seconds. Alignment replay includes both identity passes; the old implementation also loaded and concatenated values there.

| Case | Scan/normalize/stage | Alignment replay | Identity add | Identity select | Final value copy |
| --- | ---: | ---: | ---: | ---: | ---: |
| sharded | 2.784 → 1.651 | 0.463 → 0.041 | 0.414 → 0.194 | 1.915 → 0.500 | 0.779 → 0.138 |
| wide-sharded | 2.655 → 1.663 | 0.337 → 0.092 | 0.048 → 0.022 | 0.202 → 0.060 | 0.422 → 0.192 |
| wide-whole | 1.651 → 1.582 | 0.200 → 0.091 | 0.028 → 0.020 | 0.162 → 0.059 | 0.293 → 0.309 |

The last column compares old `_finish_shard()` with new `_write_values()`, so it is not an identical unit of work: the old function also reread and wrote metadata; the new implementation writes metadata during selection. Old retained-frame serialization and other orchestration are outside the displayed component timings. In the JSON, `_select_rows` includes its alignment, identity selection, and Parquet writes; do not sum that enclosing phase with its nested timings.

An initial raw-write implementation slowed wide scans. A separate 2,090 × 1,000 float32 probe showed that `tofile()` on a column-major array took 68–78 ms, versus 4.5–4.9 ms after `np.ascontiguousarray()`. The final implementation explicitly makes a bounded contiguous row tile before writing. It does not move that copy into identity indexing.

### SQLite pilot

With 176,000 rows and explicit 962-row annotation tiles, isolated configuration pilots produced the following results. Each pilot is one run, so nearby timings should be treated as noise rather than a ranking.

| Keys per lookup | Rows per transaction | Select time (s) | Add time (s) | SELECT statements | Commits |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 262,144 | 2.264 | 0.259 | 176,000 | 1 |
| 128 | 262,144 | 0.774 | 0.253 | 1,474 | 1 |
| 512 | 262,144 | 0.798 | 0.259 | 374 | 1 |
| 2,048 | 262,144 | 0.812 | 0.262 | 198 | 1 |
| 512 | 962 | 0.828 | 0.390 | 374 | 183 |
| 512 | 65,536 | 0.830 | 0.278 | 374 | 3 |

The evidence supports bounded batched lookups and fewer commits. It does not justify an unbounded key dictionary, a much larger lookup batch, or an assertion that every SELECT was a physical disk read. The chosen 512-key limit is within the observed plateau; 262,144-row transactions amortize commit overhead while retaining a finite batch and output-owned journal. With automatic tiles, the final sharded case uses 352 SELECT statements and one commit, versus the original 176,000 SELECT statements and 198 commits.

## Regressions and remaining bottlenecks

- Wide-sharded peak RSS increased 30.4%, from 274.4 to 357.9 MiB. Corrected parser tiles are larger; the 16 MiB numeric-cell estimate is not a cap on pandas parser buffers, metadata, copies, imports, or allocator retention. Explicit smaller `chunk_rows` remains available at this private preparation seam. The automatic target was not increased.
- Wide-whole sampled scratch increased 3.7%, about 5.0 MiB, and its final value-copy phase was slightly slower. Whole-genome numeric streams coexist with final matrices until the aligned group finishes. Narrow retained-row locators add storage, and the 25 ms sampling can miss brief before/after peaks. There is no demonstrated whole-genome scratch reduction in this measurement.
- Small inputs incur extra fixed stream files and slightly more scratch bytes. No small-input fast path was added: it would duplicate preparation logic for a millisecond-scale case.
- Parsing, normalization, and initial staging now occupy roughly 60–77% of elapsed time in the three larger cases. Numeric column writes, metadata alignment, and individual SQLite upserts remain. Preparation parallelism and larger architecture changes remain separate follow-ups.

Reducing file creation, numeric serialization, SQL dispatch, and commits should reduce corresponding work on HPC, but the elapsed benefit depends on its filesystem, caches, input widths, chromosome sizes, and parser costs. These runs do not predict an HPC speedup or the total gene-index build time. All four requested opportunities were implemented; more aggressive batching or additional parallelism was not justified by this scope and evidence.

## Reproduction

From the repository root, use the configured `ldsc3-dev` environment. The harness imports whichever package is on `PYTHONPATH`, so use separate source trees to compare revisions. Preserve the original source before editing; the recorded baseline revision is above. Each output directory must be new.

```bash
python benchmarks/annotation_preparation.py generate /private/tmp/annotation-preparation-inputs
PYTHONPATH=/absolute/path/to/before/src python benchmarks/annotation_preparation.py measure /private/tmp/annotation-preparation-inputs /private/tmp/annotation-before-sharded --case sharded
PYTHONPATH=src python benchmarks/annotation_preparation.py measure /private/tmp/annotation-preparation-inputs /private/tmp/annotation-after-sharded --case sharded
```

Repeat with `wide-sharded`, `wide-whole`, and `small`, using a fresh output directory for every repetition and alternating variant order. No input regeneration is needed between variants. `--declared` exercises declared chromosome labels for sharded cases. To reproduce a SQL pilot on the final implementation:

```bash
PYTHONPATH=src python benchmarks/annotation_preparation.py measure /private/tmp/annotation-preparation-inputs /private/tmp/annotation-sql-pilot --case sharded --chunk-rows 962 --sql-batch-size 512 --transaction-rows 262144
```

The SQL options are benchmark-only overrides of private constants; they add no public API or CLI options. All generated inputs and measurement directories for this session were placed under `/private/tmp/ldsc-annotation-preparation-20260914/`.

## Validation

- Before production changes: 25 focused storage/streaming tests passed. New tests first demonstrated the excessive spool file count, per-key SQL dispatch, per-chunk commits, chromosome-dependent scan budget, and numeric deserialization during identity indexing.
- After implementation: 157 focused tests and 27 subtests passed across storage, preparation, streaming annotation, gene-index streaming, gene-resolution storage, annotation, and SNP identity. The final storage/preparation rerun passed 33 tests in 6.61 s, including two additional integral-float `chunk_rows` compatibility cases.
- Full pytest: 1,653 passed, one skipped, 132 subtests passed, 264 warnings, in 127.13 s. The skip is the existing PyArrow-unavailable branch when PyArrow is installed; warnings concern existing rsID build handling and pandas assignments outside this change.
- Both `ldsc --help` and `python -m ldsc --help` passed. Edited Python files compile and `git diff --check` passed.
- Standard-library unittest compatibility: 1,000 tests ran in 53.209 s, `OK (skipped=1)`, after pytest finished. These suites were run sequentially to avoid shared pybedtools temporary-file interference.

Focused guarantees include duplicates spanning chromosomes, aligned sources describing the same logical rows, different baseline/query widths with unequal scan boundaries, omitted/present alleles, strand-equivalent duplicates and allele conflicts, chromosome filtering after global cleanup, diagnostic reason/source ordering, exact signed/fractional values and float32 output, Fortran-order NumPy files, row/column ordering, invalid values on excluded chromosomes, row-count/order failures, empty retained results, numeric-read failure cleanup, bounded SQL parameters and transaction boundaries, spool replay, and preservation of unrelated files. The complete suite also exercises direct/indexed computation and canonical output writers through the shared preparation path.
