# Memory and CPU efficiency design

Last updated on: 2026-10-27

<a id="annotation-and-workflow-memory-design"></a>

This is the contributor design guide for implemented memory, CPU, and serialization optimizations across the package.
Its historical filename remains stable for existing links; its scope includes annotation storage, genotype and R² access, direct and indexed projection, worker transport, regression, and post-fit processing.
The intended workload includes thousands of pathways tested separately against shared baseline categories, as well as reference-panel construction and repeated use of saved artifacts.
It records mechanisms, owners, scientific invariants, costs, and remaining large allocations so future changes can preserve the design.
Munging-specific optimization remains outside the current work; shared identity helpers are included because other workflows use them.

Current code and the contracts linked below define implemented behavior.
Dated decisions and the measurements recorded below describe their original revisions; they are not evidence that a historical implementation is still current.
In particular, sequential query execution alone did not bound indexed fragments before `659738a`, and the June 1 panel-building spec predates implemented unrestricted BED streaming.
The older annotation benchmark scripts were removed from the active tree on September 15, 2026; their links below point to archived Git revisions.

## Design catalog

| Design | Main owner / section | Resource effect and accepted cost |
| --- | --- | --- |
| Two content passes and disk-backed global identity checks | [Source preparation](#bounded-source-preparation) | Avoid a dense whole-genome annotation matrix and numeric spool; repeat input parsing/decompression. |
| Exact binary classification and bit-packing | [Annotation format](#annotation-format-policy) | About 32 times smaller binary payload than float32; decode only selected tiles. |
| Shared baseline/control plus active query preparation | [BED and gene queries](#bed-and-gene-queries) | Retain one generated query batch; preserve complete source validation outside the batch loop. |
| Selective BED reads and unrestricted build streaming | [Genotypes and panel construction](#genotype-access-and-r2-panel-construction) | Avoid whole-source/sample-filter copies; streaming rereads genotype bytes. |
| Numeric pair emission and compact R² storage | [Genotypes and panel construction](#genotype-access-and-r2-panel-construction) | Avoid per-pair Python objects and repeated SNP identities; R² quantization is lossy. |
| Output-row accumulators, tiled float64 operands, fused weights | [Direct LD scoring](#direct-ld-score-computation) | Reuse each correlation block within a batch; another batch repeats reference traversal. |
| Bounded R² pair buffers and rectangular sparse products | [Sparse projection](#streamed-r2-projection) | Read stored pairs once per traversal; CSR construction and operand/product tiles still allocate. |
| Baseline-by-query overlaps and query diagonals | [Compact sufficient statistics](#compact-count-and-overlap-statistics) | Avoid a query-by-query matrix without changing separate-query models. |
| Persistent workers, bounded dispatch, controlled native threads | [Chromosome parallelism](#chromosome-parallelism) | Reuse startup and avoid thread oversubscription; each worker adds a live working set. |
| One active fragment batch and chromosome-order assembly | [Direct scoring](#direct-ld-score-computation), [index access](#immutable-gene-index-access) | Avoid genome-wide query transport/concatenation; write compressed fragments and staged final files. |
| Immutable atom index and Boolean union selectors | [Index construction](#gene-index-construction-and-reuse) | Amortize genotype work across later runs; pay construction, operator storage, and per-batch reload costs. |
| Shared regression inputs, read-only maps, owned model buffers | [Batch regression](#batch-regression-and-incremental-result-writing) | Avoid baseline IPC/copies; complete genome-wide models and private numerical work remain. |
| Block sufficient statistics and incremental fit artifacts | [Regression kernels](#regression-kernels-and-saved-result-reuse) | Avoid delete-block data reconstruction and retention of every estimator; keep exact model/jackknife semantics. |
| Exact quantile boundaries plus streamed aggregates | [Quantiles](#exact-global-quantiles) | Avoid whole fitted-annotation and SNP-by-quantile matrices; retain one eligible target vector for sorting. |
| Vectorized identities, reusable lookup state, streamed diagnostics | [Lookup reuse](#identity-and-pair-query-reuse), [diagnostics](#diagnostics-and-verification) | Reduce repeated Python work and audit retention; caches and some full diagnostics remain explicit costs. |

## Precision, universes, and resource accounting

Storage compression, arithmetic precision, and reduced live data are separate design choices.
Packed annotation bits represent exact 0/1 values losslessly; public annotation reads remain float32.
Continuous annotations retain the established float32 representation, while external quantile targets retain float64.
LD-score accumulation, annotation sums, overlap products, index operators, and saved coefficient delete values use float64.
Public LD-score Parquet values remain float32.
The existing int16 R² panel encoding is lossy and must not be described as equivalent to lossless annotation packing.

Preserve the [SNP-universe contract](../../reference/ldscore-snp-universes.md): reference contributors and all/common counts are broader than regression/output rows.
Moving the output mask before score allocation does not remove contributors.
Preserve physical BIM/BED mappings, genomic SNP order, query declaration order, Boolean interval unions, negative unbiased R², one diagonal contribution, and inclusive native common-MAF thresholds.
The direct route checks query variation at its existing post-kernel float32 precision; the indexed route checks float64 fragments before public narrowing.
Changing operation order may change floating-point rounding; algebraic equivalence is not a promise of bitwise equality across batch widths or backends.

Let $S_o$ be output rows, $B$ shared baseline columns, and $q$ active query columns for one chromosome.
For direct projection, the main score buffer contains $8 S_o(B+q+1)$ bytes when a weight column is present.
This is one array's payload, not a worker or job memory estimate.
Genotypes/operators, annotation tiles, index maps, CSR construction, conversion buffers, metadata, libraries, and concurrent workers add to the live set.
The fragment score payload for one active genome-wide batch scales with total output rows times $q$; direct fragments hold float32 and indexed fragments hold float64 before compression.
It does not scale with every requested query at once.
Privately staged final outputs still grow with total query count until publication.

Measure anonymous memory, file cache, process-tree memory, peak private scratch, final output size, CPU time, and wall time separately.
Neither logical I/O counters nor allocated file blocks measure resident anonymous memory.
Closing arrays or deleting fragments shortens their lifetime but does not force allocators or the OS page cache to release charged pages immediately.
A near-limit Slurm/cgroup peak alone establishes neither a Python leak nor a minimum RAM requirement.

## Identity and pair-query reuse

[`_kernel/snp_identity.py`](../../../src/ldsc/_kernel/snp_identity.py) normalizes allele sets with vectorized lookups while preserving identity and rejection reasons.
[`read_snp_restriction_keys()`](../../../src/ldsc/_kernel/identifiers.py) parses and collapses a restriction map once per process and file-identity/mode/build key, so shared reference and regression restrictions can reuse it.
The key includes resolved path, modification time, and size; this is a reuse key, not a content-integrity guarantee.
The cache retains full filters, has no automatic eviction, and can be reset with `clear_snp_restriction_keys_cache()` in long-lived callers.
Inputs are expected to remain unchanged during a workflow.

[`R2Panel`](../../../src/ldsc/r2_query.py) lazily retains per-chromosome sidecar/index/Parquet reader state for repeated pair queries.
[`lookup_pairs_in_parquet()`](../../../src/ldsc/_kernel/r2_query.py) canonicalizes query endpoint indices into numeric keys, deduplicates lookup work, and restores requested order and repeated queries in the result.
Small requests prune row groups using `IDX_1` footer ranges; absent statistics conservatively keep the group.
Larger requests scan the file sequentially, matching numeric keys without constructing an all-pairs DataFrame.
Returned arrays still scale with requested pairs, and a reused panel handle can retain metadata for every chromosome it touches.
This query access pattern is separate from the full sequential reader used for LD-score calculation.
See [the pair-query reference](../../reference/query-r2.md) and [`test_r2_query.py`](../../../tests/test_r2_query.py).

## Ownership and lifetime

[`AnnotationBundle`](../../../src/ldsc/_annotation_bundle.py) represents one complete dataset through shared metadata and chromosome descriptors.
It has ordered baseline/query names, source provenance, configuration, scope, and paths needed for access.
File-backed shard dictionaries contain descriptors.
The explicit small-input `AnnotationBundle.from_frames(metadata, baseline_annotations, query_annotations, config_snapshot=...)` route instead owns packed binary arrays, dense float32 continuous arrays, and complete diagnostics in memory, with no staging directory.
`bundle.read(chrom, rows=..., columns=...)` returns detached float32 values owned by the caller; `metadata_for_chromosome(chrom)` reads that chromosome's metadata without caching it.
No property implicitly concatenates every chromosome.

[`AnnotationWorkspace`](../../../src/ldsc/_annotation_storage.py) owns a unique `.ldsc-annotation-*` directory immediately below the supplied output directory.
Source-backed Python preparation requires `AnnotationBuilder.run(..., output_dir=...)`.
Use its returned bundle as a context manager or call `close()` after the last borrower finishes.
The calculator borrows a bundle and does not close it.
Closure is idempotent, deletes only owned private staging, and invalidates subsequent reads.
Original inputs must remain available and unchanged while a source-backed result depends on them.

Standalone `run_annotate()` uses [`_annotation_outputs.py`](../../../src/ldsc/_annotation_outputs.py) to return saved sources.
It writes persistent `query.<chrom>.annot.gz` files and, in annotation-grid mode, aligned `baseline.<chrom>.annot.gz` files incrementally.
Returned descriptors refer only to exported files; original baseline/BIM inputs are not dependencies.
BIM grids are scanned in bounded metadata chunks by [`_annotation_bim.py`](../../../src/ldsc/_annotation_bim.py), then use the shared disk identity index and retained-metadata writer.
They require no genotype data or placeholder baseline value stores.
Generated queries use packed binary stores; exported baselines retain the existing exact classification and float32 continuous-value contract.
The writer reads row tiles sized for the selected family width, preserving a bounded numeric working set.
It removes its construction scratch before returning.
A later explicit read can prepare output-owned private shards; the caller closes that owner when done.
Persistent outputs survive closure.
Direct and indexed `run_ldscore()` publish their diagnostics before releasing their private workspace, so returned diagnostic references remain valid.
A single-batch `LDScoreCalculator.run(..., output_config=None)` accepts already-prepared inputs, keeps diagnostics in memory, and creates no annotation, diagnostic, temporary, logging, or result files.
Multiple query batches require an output configuration; high-level `run_ldscore()` and the CLI always require an output directory.

Handled failures remove owned scratch while preserving the command's existing failure diagnostics and overwrite markers.
An interrupted process can leave its private directory behind.
A new run creates fresh staging and does not automatically reuse or remove abandoned directories.
The existing gene-index construction transaction is the explicit exception to output-contained staging: sibling build-state, locks, staging, publication, and recovery remain governed by the [index contract](../../reference/gene-ldscore-index.md).

## Bounded source preparation

[`prepare_annotation_sources()`](../../../src/ldsc/_annotation_sources.py) makes two content passes over each supplied text/gzip source, in addition to header preflight.
The first validates and classifies numeric columns before float32 conversion, and stages only row metadata.
Discovery and preflight use contents, so selecting a chromosome never requires first loading an entire annotation matrix.
Each file is scanned independently with an automatic tile targeting 16 MiB of parsed cells, estimated at eight bytes per header column, and at most 65,536 rows.
After discovering the layout, metadata streams are reblocked to a common tile using the combined header width of only the files aligned together.
Whole-genome column sources form one group; corresponding chromosome baseline/query shards form separate groups.
Adding unrelated chromosome shards therefore does not shrink the tile.
Explicit `chunk_rows` controls scanning, alignment, and the second content pass.
The target sizes numeric tiles; it is not a total process RSS limit and does not include variable-length metadata strings or parser/conversion copies.

Aligned baseline/query column files are checked positionally before chromosome partitioning.
Corresponding rows describe one logical SNP observation.
Both identity passes read only metadata: [`DiskIdentityIndex`](../../../src/ldsc/_annotation_identity.py) applies drop-all duplicate cleanup across the complete logical dataset, including cross-chromosome rsID collisions.
Its SQLite table stores narrow identity/count state with an 8 MiB page cache.
Parameterized lookups use at most 512 keys, reduced to the connection's variable limit when necessary; lookup dictionaries cover only the current bounded tile.
Upserts commit after at most 262,144 valid rows, with the final partial transaction committed before selection.
The journal stays in the workspace and any temporary IN key set is bounded by the lookup batch.
Fewer statements and commits reduce execution overhead; SQL statement counts do not measure physical disk reads.
Invalid allele rows are excluded from duplicate clusters under the existing rules.
A chromosome shard omitting both allele columns keeps the existing base-identity fallback; an unknown allele set cannot establish a multi-allelic conflict.

After global indexing, selection writes retained metadata directly to its final Parquet file and records bounded source-row/chromosome locators in one append stream per aligned group.
The retained counts size the final NumPy stores.
[`_write_values()`](../../../src/ldsc/_annotation_sources.py) rereads the original input at the aligned group's tile width and writes retained values directly into their final binary/continuous stores.
There is no numeric `values.bin` spool, replacement numeric staging file, or conversion copy on disk.
The repeated parsing/decompression is the accepted trade-off for this simpler preparation path.
Inputs must remain unchanged between passes.
Source metadata staging is removed after each source finishes; group selection locators are then removed.
[`FrameSpool`](../../../src/ldsc/_annotation_storage.py) uses one file of independently serialized frames per spool, including diagnostic spools; descriptors keep no open handles.
Selection preserves source order within each chromosome and applies requested chromosome filtering only after global identity cleanup.
Preparation remains serial before the gene-index builder starts chromosome workers.

A September 15 parallelism evaluation tested chromosome selection and writing after committing the shared global identity index.
Four process workers reduced local preparation time by 25.3% on 1.89 million real baseline rows, while increasing sampled process-tree RSS from 399 to 1,075 MiB.
Threads regressed, and a synthetic 1,000-query workload had no worthwhile total gain.
Numeric writing did not dominate either workload.
The user chose to keep preparation serial because the gains do not justify the memory and implementation costs.
The parallelism follow-up is closed, including further writer-only experiments.
Production preparation does not consume `--workers`; existing chromosome computation parallelism is unchanged.

## Annotation format policy

Type follows processing logic and complete value evidence, independently of baseline/query membership.
No dtype flag or user-supplied type declaration is required.
The [September 15 decisions](../../current/annotation-memory-decisions.md#confirmed-annotation-storage-scope-2026-09-15) govern this package-wide policy.

| Annotation source or case | Classification | Private representation |
| --- | --- | --- |
| Generated gene-list or BED membership, fixed gene control, synthetic all-ones `base` | Binary by construction, including all-zero and all-one columns | Packed bits written directly during projection/preparation |
| Supplied baseline or query annotation column containing only exact numeric 0/1 | Binary after exhaustive classification across every chunk and chromosome in the supplied dataset | Packed bits |
| Supplied column containing any other numeric value | Continuous for every chromosome of that dataset, even if individual shards contain only 0/1 or mostly zeros | Dense float32 |
| `AnnotationBundle.from_frames()` | Same exhaustive, pre-float32 classification across the input frames | Packed binary arrays and separate dense float32 arrays in RAM; no files |
| External quantile target | Separate target-value contract, not an ordinary annotation | Existing float64 storage and missing-token rules |

Classification uses parsed numeric values before float32 narrowing: for example, `1.00000001` keeps a supplied column continuous even though it rounds to one in float32.
There is no sampling, tolerance, rounding-based classification, or sparse-matrix choice.
Any nonbinary value in any scanned source row keeps that column continuous, including rows later excluded by identity or chromosome selection.
Ordinary annotation missing/nonnumeric input errors retain their existing behavior; zero is an observed value, not a missing-value code.
Existing advisory classification of fitted annotations still describes the values in its scientific analysis universe; it is separate from this source-wide storage choice.

### Encoding, identities, and order

[`ColumnStore`](../../../src/ldsc/_annotation_storage.py) uses private NPY v2 files in column-major order, with separate Parquet row metadata.
For a shard with $S$ logical SNP rows, $B$ binary columns and $C$ continuous columns:

- The binary file has physical shape `(ceil(S / 8), B)` and dtype `uint8`.
  Bits run along the SNP axis within each annotation column.
  `bitorder="little"` puts row zero in the least significant bit; row $r$ is bit `r % 8` of byte `r // 8`.
  Unused high bits of each column's final byte are zero and never become SNP rows.
- The continuous file has shape `(S, C)` and dtype `float32`.
  Empty type groups need no file.
  Supplied sources use `binary.npy` and `continuous.npy` within each chromosome directory; generated batches can add their own packed query store.
- The descriptor retains logical `n_rows`, ordered column names, logical dtype (`bool` for packed columns), NPY data offset, and explicit bit order.
  The NPY header describes physical shape/dtype/order.
  Descriptors carry this information to workers; a packed NPY file alone does not describe its logical SNP axis or names.
- All stores in one shard use exactly the same metadata row correspondence.
  Partial writes preserve neighboring bits at shared byte boundaries.
  Selected reads preserve requested rows and columns, including duplicates; padding is excluded by the logical row count.

Logical order remains baseline columns followed by queries, each in resolved input order.
Gene queries keep declaration order; filtering preserves surviving queries' relative order.
[`AnnotationBundle.read()`](../../../src/ldsc/_annotation_bundle.py) supplies that default list.
Private multi-store `AnnotationShard.read()` requires explicit ordered columns and reuses the name-to-store/position lookup to restore the requested order.
There is no additional shard-level logical-order list.
Physical binary/continuous grouping therefore cannot reorder default public reads, generated `.annot.gz` columns, or LD-score output columns.

### Bounded decoding and public compatibility

Column-major storage allows selected columns to be read without loading other queries.
File reads seek bounded, coalesced row runs, using at most `max_read_rows` logical rows plus up to fourteen covering-byte boundary bits per decoded run.
The default is 65,536 rows.
They return detached arrays and retain no memory map or value cache.
In-memory packed reads extract only the requested bits.
Internal binary-only reads may use Boolean arrays; `AnnotationBundle.read()` always returns float32, including binary-only and empty selections.
The numerical kernels retain float64 LD-score/count/overlap accumulation and the current float32 LD-score publication policy.

The returned selection itself is materialized.
Callers must select an LD block and query batch to bound its size; an explicit whole-chromosome read can still be large.
The package's computation tiles and execution batches supply those limits.
Packing alone does not shrink score accumulators or assembly copies; the separate projection and streaming designs below address those allocations.

These files are private, owner-scoped intermediate artifacts.
Canonical annotation text remains `.annot.gz` with generated 0/1 query values and supplied float32 baseline values, and canonical LD-score Parquet files, query-batch manifests, column order, and precision remain unchanged.
Reloading saved annotations uses the same automatic preparation policy.
No persistent packed-artifact format, migration, resume, or new user option is introduced.

### Size and memory evidence

Private value payload size per shard is `ceil(S / 8) * B + 4 * S * C` bytes, excluding NPY headers, Parquet metadata, and identity/diagnostic scratch.
The former all-float32 payload was `4 * S * (B + C)`.
Binary payload approaches a 32-fold reduction; row padding and file headers matter for small shards.
Binary density does not affect packed size.
Continuous density does not affect dense size.
File-free retained arrays use the same payload formula, excluding their Python/NumPy and metadata overhead.
Parsed input frames and returned float32 tiles remain larger transient allocations.

$S$ is the actual annotation shard row count after annotation identity cleanup and requested chromosome selection.
Supplied and generated query shards use the cleaned annotation grid, which can include SNPs absent from the retained reference panel.
During LD calculation, [`MappedAnnotations`](../../../src/ldsc/_kernel/ldscore_projection.py) maps the retained reference contributors onto those annotation rows; it does not create another complete dense annotation matrix.
LD output rows are a further regression/output subset.
Consequently, the reported 1,134,471 LD-score output rows cannot be used to size the intermediate annotation artifacts.

Local measurement on September 15 prepared `resources/baseline_ld_suites/baseline_v1.2/baseline.22.annot.gz` from the workspace and compared every decoded value against an independent pandas read aligned by SNP.
All values matched; all 53 columns were classified binary.
No source input was modified.

| Measured chromosome-22 quantity | Value |
| --- | ---: |
| Source and retained annotation rows | 141,123 |
| Packed value payload | 934,973 bytes |
| Packed NPY file, including header | 935,101 bytes |
| Separate Parquet metadata file | 2,125,835 bytes |
| Equivalent float32 value payload, calculated from those dimensions | 29,918,076 bytes |

This measurement covers annotation preparation for one local baseline shard.
It does not measure the 18,809-query Longleaf run, its intermediate row counts, or a new peak RSS.
Packing reduces annotation storage; it does not imply a 32-fold reduction in total job RSS.
The historical packing change left score accumulators and aggregation alone; later changes documented here independently reduce their width or lifetime.

## BED and gene queries

[`_annotation_queries.py`](../../../src/ldsc/_annotation_queries.py) and [`_gene_query_storage.py`](../../../src/ldsc/_gene_query_storage.py) prepare query sources separately from baseline artifacts.
BED intervals are normalized with bounded buffers.
Gene identifiers resolve once against shared catalog lookup state; selected catalog indices and complete row audits are staged per source.
Direct BED/gene projection prepares only the active execution batch, builds each query on each chromosome, and closes that private batch owner before advancing.
Global source resolution, duplicate cleanup, coverage/support gates, and the shared baseline/control remain outside that loop.
A fixed gene control, when requested by direct/indexed LD scoring, belongs to the shared baseline group.

Standalone [`annotate_workflow.py`](../../../src/ldsc/annotate_workflow.py) accepts exactly one BED or gene-list route.
Gene mode requires explicit nonnegative padding, uses catalog/build and coverage gates, and evaluates support on the cleaned baseline annotation grid.
It skips empty or globally unsupported focal queries while retaining usable siblings, all-one queries, and zero-valued chromosome shards for queries supported elsewhere.
It fails when all focal queries are skipped.
Its audits use `annotation_snp_count`; reference-panel support remains unevaluated.
Direct LD scoring evaluates support on retained reference rows and applies its additional regression-specific viability rules.
See the [gene-list condition–outcome table](../../reference/annotate.md#gene-list-queries).

## Genotype access and R2 panel construction

[`PlinkBEDFile`](../../../src/ldsc/_kernel/plink_bed.py) reads physical SNP records selectively and applies individual filtering during each read.
It avoids coexisting full-source and sample-filtered bitarrays and preserves the physical-row mapping when retained SNPs are reordered or dropped.
The ordinary reader retains the selected packed genotype payload in memory; standardized numerical genotype blocks are decoded as needed.
It does not imply that direct LD scoring or index construction streams all genotype payloads from disk.

[`ReferencePanelBuilder._build_chromosome()`](../../../src/ldsc/ref_panel_builder.py) enables disk streaming when there is no explicit reference-SNP restriction.
This path retains SNP/MAF bookkeeping and reads genotype batches again for computation instead of retaining the whole selected bitarray.
An explicitly restricted build retains its selected genotype payload.
The builder requests float32 standardized genotype blocks; the direct LD/index numerical path keeps its existing float64 computation.
When float32 builder genotypes were introduced, R² shifted by roughly 1e-6 to 1e-5 relative to float64, which was accepted for R² sums.
The build-side sliding window retains the current batch and window-spanning carry-over columns, so its cost depends on individuals, SNP density, window width, and SNP batch size.
It never needs a dense chromosome-wide correlation matrix, but a dense wide window can still be expensive.
The chromosome owner releases the reader after writing.

[`yield_pairwise_r2_rows()`](../../../src/ldsc/_kernel/ref_panel_builder.py) emits columnar numeric chunks rather than one dictionary per pair.
Vectorized triangle/window masks select cross-block and within-block pairs.
`_PendingPairs` retains numeric arrays for left-endpoint groups until all their right neighbors have been seen, then emits stable ordered pairs.
`write_r2_parquet()` converts these arrays directly to Arrow without an intervening pandas frame.
Its writer buffer flushes at its configured pair threshold and can overshoot by one incoming chunk; a row-group size controls file layout, not total process RAM.
See [`test_ref_panel_builder.py`](../../../tests/test_ref_panel_builder.py) and the [R² format reference](../../reference/formats/r2-panels.md).

The [R² format](../../reference/formats/r2-panels.md) stores two int32 sidecar indices, quantized int16 R², and a bit-packed correlation sign instead of repeating identity strings per pair.
The identity sidecar is stored once per SNP and bound to pair indices by the required identity hash.
Column-specific dictionary/delta/byte-stream encodings and compression reduce file size and I/O; compression also consumes CPU.
R² quantization uses scale 32767 and preserves signed unbiased estimates within the stored range, but introduces rounding error.
Aggregate LD-score error depends on the contributing pairs and annotation weights; historical measured differences are not universal bounds.
Positive `min_r2` filtering omits scientific contributions, so `ldscore` rejects such panels; they serve `query-r2` lookups only.
The default zero threshold preserves zero and negative unbiased pairs.

## Direct LD-score computation

[`PreparedChromosome`](../../../src/ldsc/_kernel/ldscore.py) owns the current chromosome's reference reader, contributor metadata, physical-row mappings, LD-window state, selected annotation accessor, and drop diagnostics.
`close()` clears these references and closes the reader, including on failure.
PLINK physical SNP indices remain distinct from annotation order.
The reference adapters in [`_kernel/ref_panel.py`](../../../src/ldsc/_kernel/ref_panel.py) do not retain a metadata cache containing all chromosomes.

[`_kernel/ldscore_projection.py`](../../../src/ldsc/_kernel/ldscore_projection.py) separates reference contributors from output rows.
Every retained reference SNP can contribute to LD scores, counts, and overlap calculations, even when that SNP is absent from the regression/output set.
The main accumulator is float64 and has `n_output_rows × (n_output_annotation_columns + 1)` entries; the extra column is the separately governed regression-weight projection.
Completed direct scores use the existing float32 output conversion.
Moving the output mask upstream does not restrict the reference contributor universe.

Within each active direct query batch, each PLINK correlation block is generated once.
The next execution batch repeats reference preparation and traversal.
Shared baseline columns and the weight mask are projected together, then query columns are read in batches for that same block.
For parquet R2, bounded pair batches are scattered through float64 sparse products.
Counts and overlap products are accumulated from bounded annotation tiles; the retained overlap summary consists of baseline-by-all-annotation products plus query diagonals, rather than a dense query-by-query matrix.
Preserve the existing accumulation and validation tolerances when changing these kernels.
Annotation sums reduce binary 0/1 and normalized float32 continuous values into float64, fixing batch/layout-dependent quantitative counts (completion-review R1).
Independent `math.fsum` checks cover signed and high-offset values at the existing tolerances; overlap products and LD-score accumulation retain their established precision.

`--query-batch-size` defaults to **100** for direct/indexed `ldscore` and batch `partitioned-h2`.
A smaller positive value reduces active query workspace.
It does not reduce the number of requested pathways or merge them into a joint regression model.
With exactly 1,000 queries, the default uses ten batches of up to 100.
The active output-row accumulator scales with batch width.
Only compact counts, baseline-by-query overlaps, and diagnostic/source descriptors scale with the total requested query count; saved query LD values do not accumulate in memory.

`--snp-batch-size` defaults to **512** for PLINK input in `ldscore`, `build-gene-ldscore-index`, and `build-r2-panel`; parquet R2 input ignores it.
LD windows are exact for every batch size, so the setting changes only speed and memory.
Each batch computes about one batch width of out-of-window pairs per SNP before dropping them, and holds a few `(window + batch) × batch` float64 blocks beside the `samples × window` genotype block.
Measured on 2026-10-27 at revision `38c98d8` (macOS, Apple Accelerate, one BLAS thread) by timing the PLINK kernel alone; LD scores were identical at every size, and memory is the peak above the loaded process:

| Batch | 1KG EUR chr22: 489 samples, 1 cM, about 2,470 SNPs per side | 1KG 30x chr22: 3,202 samples, 30,000-SNP stretch, 1 Mb, about 9,300 SNPs per side |
|---|---|---|
| 128 | 11.6 s, +0.59 GB | 25.0 s, +0.75 GB |
| 256 | 10.7 s, +0.57 GB | 21.0 s, +1.07 GB |
| 512 | 11.7 s, +0.80 GB | 17.7 s, +1.19 GB |
| 1024 | 13.7 s, +0.94 GB | 17.0 s, +1.35 GB |
| 2048 | 16.7 s, +1.25 GB | 17.6 s, +1.24 GB |
| 4096 | 23.6 s, +1.67 GB | not run |

512 was within about 4% of the fastest size on both panels; larger batches add time and memory without changing results.
Limits: one machine, one chromosome stretch per panel, and single runs; Linux BLAS libraries may favor somewhat larger batches.
Memory outside the kernel does not depend on the batch size, including each chromosome's in-memory genotypes (`SNPs × samples / 4` bytes) and annotation reads, and it repeats for each `--workers` process.
`build-r2-panel` uses the same chunk walk but was not timed; its stored pairs are identical for any batch size, with R2 differing by at most one int16 unit from float32 summation order (see `tests/test_plink_exact_windows.py`).

For direct file-writing workflows, [`stage_direct_chromosome()`](../../../src/ldsc/_direct_ldscore_batches.py) writes the existing post-kernel float32 query table to a private Parquet fragment before the worker returns.
The parent receives baseline data, counts, overlaps, variation summaries, and diagnostic descriptors, then combines them without a genome-wide query matrix.
[`_query_variation.py`](../../../src/ldsc/_query_variation.py) combines each chromosome's first values and within-chromosome variation flags, including columns that are constant within chromosomes but differ between chromosomes.
Direct filtering uses the same float32 values as the materialized API; indexed filtering continues to use float64.
The small prepared no-output Python route still returns a materialized `LDScoreResult`.

### Streamed R2 projection

[`SortedR2BlockReader.iter_all_pairs()`](../../../src/ldsc/_kernel/ldscore.py) visits row groups once in stored order and decodes bounded Arrow record batches of at most 65,536 rows, including when a stored row group is larger.
It reads only endpoint indices and R², with no sign column, decoded-row-group cache, or dense window matrix.
Mode-dependent sidecar identity validation and remapping happen once per chromosome reader.
Builder pair uniqueness plus an injective retained-SNP remap remove the need for a repeated per-pair identity join or deduplication set.

`ld_score_streaming_from_r2_reader()` splits decoded arrays into buffers of at most `_CSR_CHUNK_PAIRS` pairs, currently 16 million, independent of Parquet row-group size and PLINK SNP batch size.
[`ProjectionAccumulator.add_pairs()`](../../../src/ldsc/_kernel/ldscore_projection.py) filters the read-time LD window, retains each directed contribution whose target is an output row, and builds a rectangular float64 CSR over the involved targets and contributors.
It keeps contributions from non-output reference SNPs.
This is the selected-output implementation of the conceptual product $(RA)[O,:]$; production does not materialize a full-reference upper-triangular operator and full dense annotation matrix.

`ProjectionAccumulator._accumulate()` reads float64 contributor/column tiles and reuses each tile across output blocks.
Its default 16 MiB read target sizes a numerical operand tile, with at least one row; it is not a cap on all allocations.
Sparse slices, bidirectional pair arrays, concatenation/sorting, gathered accumulator entries, and multiplication outputs can coexist.
The pair-count limit likewise is not a 0.45 GiB worker-memory guarantee.
A narrower read-time window still scans the stored pairs, but can reduce the retained CSR and arithmetic; it may therefore change memory as well as runtime.
See the [mathematical specification](../../reference/methods/ld-score-computation.md) for selected-output equations and finite-precision caveats.

### Compact count and overlap statistics

[`annotation_statistics()`](../../../src/ldsc/_kernel/overlap.py) computes sums and products on bounded annotation tiles over the full retained contributor universe.
For $B$ baseline categories and $Q$ independently tested queries, it retains a $B\times(B+Q)$ baseline-by-all block and $Q$ query diagonals, for both all and common SNPs.
No query-query off-diagonal products are needed for the separate query-plus-baseline regressions.
This changes retained overlap storage from quadratic in all query columns to the sufficient statistics those models actually consume.
It must not be reused as a full joint-query overlap matrix.
Common counts use the inclusive configured MAF threshold, independent of output-row selection.

## Chromosome parallelism

`--workers` / `workers` defaults to one managed child process, including on HPC.
Every worker uses **1 BLAS thread**, in serial and parallel execution.
Managed children also limit OpenMP and Arrow CPU/I/O thread pools under the [execution policy](../../reference/parallel-execution.md).
The runtime checks actual native limits inside workers; changing only the parent's settings would not establish spawned-child behavior.
`joblib.cpu_count()` supplies the detected available CPU count; it does not measure idle CPUs or read the requested scheduler allocation.
Positive requests above that count fail before capping by available tasks; negative `-k` uses `max(1, detected_count + 1 - k)`.
LD scoring and index construction cap by chromosomes; regression caps by complete query-plus-baseline models and query batch size.
Zero, Boolean, and non-integer Python requests are invalid.
See [the execution policy](../../reference/parallel-execution.md) for native-limit verification, memory costs, process startup, and failure handling.

The sequential lifecycle is `load chromosome → prepare → compute/write diagnostics → close → next chromosome`.
Reference and annotation working arrays from the completed chromosome become available for reuse before the next preparation.
Python/NumPy allocators may keep freed pages, so RSS need not fall immediately.
A flat RSS trace alone is not evidence of retained chromosome objects; lifetime tests assert release directly.

Parallel submission keeps pending chromosome work bounded by the effective worker count.
Tasks receive descriptors, not eager annotation matrices for every chromosome.
Each worker owns its active chromosome state and releases it before advancing; the pool persists across batches and publication.
`LDScoreCalculator.run()` yields a `StagedLDScoreBatch` to `LDScoreDirectoryWriter.write_batches()` after all chromosome statistics are available.
The shared fragment writer reads, validates, and writes one chromosome at a time in genomic order, then releases its frame and deletes its fragment before reading another.
The entire active fragment directory is removed before the next batch, including when every query in a batch is skipped.
On failure, workers finish shutdown before borrowed annotation or fragment workspaces are removed.
Worker count and query batch size control different dimensions of memory; more workers can increase total process-tree memory and may be slower for small inputs.
Baseline tables, worker kernel arrays, and staged final outputs have separate lifetimes, so bounded query fragments do not establish a whole-job RAM limit.
See [the replay measurements](../../wiki/main-functionalities/ldscore.md#direct-assembly-measurements) and [`test_direct_batch_staging.py`](../../../tests/test_direct_batch_staging.py) for the evidence boundary.

**Public LD tables aggregate chromosomes within each query batch.** One batch uses `ldscore.query.parquet`; multiple batches use ordered `ldscore.query.batch00001.parquet` files.
Root `metadata.json.query_batches` records each file, ordered query columns, and chromosome row groups.
Baseline, counts, and overlap are shared once.
`_ldscore_batch_output.py` stages all query files inside the output directory, publishes after successful computation, and writes metadata last.
Its `_LDScorePublication` record retains only the shared baseline, counts, overlap, provenance, and diagnostic descriptors needed for publication.
Materialized query tables, producer SNP sets, chromosome results, and staged query paths remain with the active batch and can be released before the next batch is requested.
`ldscore_source.validate_ldscore_baseline_table()` checks the shared baseline before deferred reads; `validate_ldscore_query_table()` checks either a complete materialized query table or one fragment against the corresponding baseline rows.
Each fragment has one row-alignment check, with the existing allele-aware identity and row-order rules.
Handled failure closes the producer and removes owned staging while retaining available diagnostics.
Resume from the last completed batch is a future feature, not supported now.

Staging until computation succeeds does not make every filesystem publication step a rollback transaction.
Preserve the existing failure-marker and overwrite contracts; the gene-index builder's complete-directory transaction is a separate mechanism.
The managed scheduler allows at most one pending task per worker and collects a ready group before refilling, avoiding an unbounded queue of descriptors or completed results.
Process startup, IPC, and one native thread per worker have costs; extra workers are not assumed to improve small workloads.

## Immutable gene-index access

[`LoadedGeneLDScoreIndex`](../../../src/ldsc/gene_ldscore_index.py) retains catalog/settings metadata, chromosome component paths, and compact support summaries.
Integrity validation loads and releases one chromosome at a time.
Assembly uses `_indexed_ldscore_batches.py` with one active query batch and a persistent chromosome-worker pool.
Each task reloads its chromosome operator, writes lossless float64 Parquet for that batch, and returns only compact counts, overlap blocks, and exact variation summaries.
The coordinator reduces those summaries in index order and applies the shared query-status policy before narrowing values.
[`_ldscore_fragments.write_query_fragments`](../../../src/ldsc/_ldscore_fragments.py) validates, converts, and writes one chromosome at a time, deleting each consumed fragment.
There is no genome-wide query concatenation or retained-column copy in this path.
All active-batch scratch is removed before another batch is submitted, including when every query in a batch is skipped.
The score payload bound is proportional to output SNP rows times batch width, while each live query-score array is chromosome-sized.
Shared baseline values, compact metadata, and privately staged final outputs have separate lifetimes; the total job memory and scratch footprint is not just the raw score payload bound.
Canonical integrity checks remain intact, including validation when each operator is reopened.
`workers`/`--workers` has the direct-mode convention: default 1, positive worker counts, -1 detected available CPUs, -2 one fewer worker than that detected count, zero invalid, capped at chromosome count.
Negative offsets do not reserve or pin a CPU.

The index builder uses the same bounded annotation preparation and chromosome scheduling.
Shared baseline and regression-weight projection use one traversal; separate atom-batch traversals remain intentional.
Stored operators and common payloads remain float64.
Smaller query multiplications and separate validation/assembly passes can add overhead; performance claims require end-to-end measurements, not just a count of loaded operators.

### Gene-index construction and reuse

[`_kernel/gene_ldscore_index.py`](../../../src/ldsc/_kernel/gene_ldscore_index.py) splits the gene interval union into disjoint atoms and records a Boolean sparse gene-to-atom map.
Query selectors take the Boolean union of selected genes' atoms, so overlap, nesting, aliases, and duplicate gene selections do not count a SNP twice.
SNP-to-atom lookup uses sorted interval boundaries; atom counts and baseline-overlap summaries use bin counts and sparse membership rather than a complete dense SNP-by-atom matrix.
The [index mathematics](../../reference/methods/gene-ldscore-index-math.md) defines why multiplication by a union selector reproduces the requested query.

[`build_plink_index_chromosome()`](../../../src/ldsc/gene_ldscore_index.py) computes shared baseline and regression-weight projections in one traversal.
It then traverses the genotype blocks again for each atom batch, retaining the resulting float64 sparse operator blocks until their final horizontal assembly.
Atom batching bounds the active dense projection but does not bound the complete chromosome operator or the final assembly copies.
This builder also explicitly reads its chromosome baseline matrix into float64.
These are remaining large allocations, not eliminated by compressed NPZ storage.
Construction workers stage components and return `StagedIndexChromosome` descriptors rather than sending full operators through process IPC.

Later indexed query runs amortize that genotype work through sparse operator-by-selector multiplication.
The immutable index fixes the reference universe, baseline, gene geometry, output rows, and LD-window settings; live overrides would invalidate reuse and remain rejected.
Loading validates components one chromosome at a time without caching every operator.
Within scoring, the persistent process pool avoids repeated process startup, but each chromosome task reloads its operator for the active query batch.
This deliberate extra reading/validation trades CPU and I/O for a bounded fragment backlog and shorter operator/array lifetimes.
The complete [index contract](../../reference/gene-ldscore-index.md) and [`test_gene_index_streaming.py`](../../../tests/test_gene_index_streaming.py) cover ownership and numerical behavior.

## Batch regression and incremental result writing

Query concurrency defaults to one managed child.
`_partitioned_h2_parallel._QueryWorkers` reuses whole-model processes across loading batches, with 1 BLAS thread per worker.
Shared numerical columns and each active query batch use dtype-preserving read-only maps.
Nonnumeric columns and compact model metadata are loaded from private snapshots.
Each task carries its query name, column descriptor, and a unique fit directory.
Workers stage complete model results and return compact outcomes.
At most one task is in flight per process; workers exit before scratch removal.
See [query workers](../../reference/regression-options.md#43-query-workers-and-memory).

[`LDScoreSource`](../../../src/ldsc/ldscore_source.py), returned by `load_ldscore_from_dir()`, reads shared baseline/identity metadata and validates schemas without retaining all query values.
`read_queries(columns)` selects only requested Parquet columns and does not cache them.
Writing workflows return `LDScoreSource` directly; `read_queries(columns)` may span files and exceed generation width, preserves requested order, and imposes no RAM cap.
`LDScoreResult` remains the explicitly materialized result, including small prepared no-output calculations; writing workflows use private fragment descriptors instead.
Downstream `h2`, `rg`, and partitioned workflows consume only the columns they need.

[`PreparedRegressionInputs`](../../../src/ldsc/regression_runner.py) shares trait alignment, baseline LD scores, row mappings, and count lookup across a batch.
Each query batch is loaded/aligned once.
Each model still fits its complete retained genome-wide SNP set; variance filtering, outlier filtering, counts, regression weights, and jackknife work remain model-specific.
No chromosome-wise regression estimates are substituted for the genome-wide model.

`_dataset_from_query_values()` builds each model frame from borrowed aligned columns without concatenating or resetting copies of the baseline.
Zero-variance column removal preserves those shared values.
`_fit_h2_dataset()` slices only when the chi-square filter removes rows, then materializes the numerical design once.
The kernel reuses its own normalized/weighted design buffers only when the result dtype is preserved, including promotion from float32 inputs to float64 weights.
Shared inputs remain read-only throughout fitting.

`estimate_partitioned_h2_batch(..., output_dir=..., query_batch_size=100)` writes each completed fit's category table, float64 coefficient delete values, and metadata privately through [`stage_partitioned_h2_fit()`](../../../src/ldsc/outputs.py), then releases estimator details.
After final stable summary sorting, the existing writer assigns canonical folders and publishes the complete query tree.
The returned `PartitionedH2BatchResult` contains the aggregate summary, output paths, and `per_query_artifacts` referencing persistent files.
It does not contain dictionaries of all detailed DataFrames.
Baseline-only input remains one joint functional-category fit.

The writer moves explicitly owned staged scientific files into publication staging instead of copying them; persistent descriptors remain copy-only.
Run metadata is written once, while query metadata stores model fields and a relative reference to `diagnostics/metadata.json`.
The [result contract](../../reference/results/partitioned-h2.md) defines the ownership and complete-directory copying requirements.
[`test_regression_memory.py`](../../../tests/test_regression_memory.py) checks shared-input lifetime, filtering, and independent weighted least-squares agreement.
[`test_partitioned_publication.py`](../../../tests/test_partitioned_publication.py) checks exact file reuse, persistent-input preservation, metadata recovery after relocation, and failure cleanup.

Explicit `continue_on_query_error=True` permits failed query regressions while preserving private staging until the scan finishes.
Successful model artifacts publish with a compact attempted-query status ledger; strict mode remains the default and publishes diagnostics alone if any query fails.
Failed exception tracebacks are logged and detached from the retained first exception so model arrays do not accumulate through traceback references.
Panel coordinates retained with shared preparation allow post-filter block diagnostics even for rsID-only traits.
See [query failure handling](../../reference/results/partitioned-h2.md#query-failures-and-continuation).

### Regression kernels and saved-result reuse

[`_partitioned_h2_parallel._stage_frame()`](../../../src/ldsc/_partitioned_h2_parallel.py) groups contiguous numeric columns by dtype into private binary files and exposes column windows through read-only memory maps.
Labels and object-valued metadata use process-local snapshots; unsupported query dtypes retain their preserving fallback instead of silent numeric coercion.
This avoids copying the shared numeric baseline through every task message.
Maps can still populate charged file cache, and nonnumeric snapshots are loaded in each worker.
The parent owns the files until every borrowing worker has exited.

[`_kernel/regression.py`](../../../src/ldsc/_kernel/regression.py) and [`_kernel/_irwls.py`](../../../src/ldsc/_kernel/_irwls.py) reuse owned sample-size-scaled and weighted numerical buffers when the required result dtype permits it.
In-place reuse must not mutate caller-owned arrays or silently narrow promoted results.
Model filtering and weighted fitting remain private to each task.
The two IRWLS reweightings are the existing estimator, not a new early-stopping optimization.

[`LstsqJackknifeFast`](../../../src/ldsc/_kernel/_jackknife.py) uses existing block sufficient statistics: compute block $X^\top X$ and $X^\top y$ summaries once, then obtain delete-block fits by subtracting one block's summaries from the totals.
For $J$ blocks and $p$ predictors, the retained block summaries scale as $O(Jp^2)$ rather than storing $J$ copies of the SNP-by-predictor data.
Each delete fit still solves a small system; aligned block boundaries, weighting, singular-design errors, and estimator semantics remain unchanged.
[`test_jackknife.py`](../../../tests/test_jackknife.py) supplies numerical comparisons rather than asserting only that an allocation disappeared.

Legacy summary-statistic projection in [`regression_runner.py`](../../../src/ldsc/regression_runner.py) uses column-wise identity/allele masks rather than row-wise Python processing.
It preserves alignment and drop diagnostics; it still materializes the trait and projection-drop table.
[`test_legacy_sumstats_projection_performance.py`](../../../tests/test_legacy_sumstats_projection_performance.py) guards both agreement and a local performance budget, not a universal machine-independent runtime.
Saved fit coefficients, counts, overlaps, and delete values support quantile analysis, heritability-scale conversion, and plotting without refitting the model.
Figure-producing paths import Matplotlib lazily, avoiding its startup work for non-plotting calls.

Multi-trait `estimate_rg_pairs()` retains loaded trait tables and compact pair summaries, but builds each pair's aligned dataset separately.
It computes the standalone per-trait heritability summaries once per trait; those are distinct from the pair-specific fits on pair-specific SNP intersections.
There is no all-pairs merged-data cache or substitution of one trait's full-support fit for a pair-specific estimate.
Current control flow is in `RegressionRunner.estimate_rg_pairs()`.
A per-trait merge cache was rejected at design time: for fewer than ten traits it was estimated to save about 10% of wall time while holding 0.8–1.2 GB of extra state and roughly 3.5 times the code.
Reusing one heritability fit across pairs is statistically invalid, because `RatioJackknife` in [`_kernel/_jackknife.py`](../../../src/ldsc/_kernel/_jackknife.py) combines the `hsq1`, `hsq2`, and `gencov` delete-block vectors element by element, so all three must come from the same pair-specific SNP partition.
Parallelizing independent pairs across processes was identified as the larger potential gain and remains unimplemented.

## Exact global quantiles

[`_quantile_inputs.py`](../../../src/ldsc/_quantile_inputs.py) and [`_quantile_storage.py`](../../../src/ldsc/_quantile_storage.py) reconstruct the validated common SNP universe with bounded annotation/reference/target reads and narrow SQLite row locators.
Reference/target duplicates remain fatal globally, including base-key uniqueness when omitted alleles must be inferred.
External targets retain float64 values and raw missing-token semantics; they are separate from ordinary packed-binary/dense-float32 annotation storage.

The first annotation pass accumulates full-common sums, overlap products, stable centered mean/M2 statistics, and the permitted eligible target vector.
It validates original model aggregates, sorts that vector to obtain the existing **exact global** boundaries, and releases it.
The second pass reuses those boundaries for every chromosome and accumulates quantile sums/counts directly.
There is no whole-genome fitted-annotation matrix or dense SNP-by-quantile indicator matrix.
Boundary ties remain in the lower quantile.
Full-common standard deviations include target-excluded SNPs, as required for standardized coefficients.

Exact global boundaries deliberately retain a target vector proportional to eligible SNP count, including sorting workspace.
This workflow is bounded in annotation width, not constant in the number of eligible SNPs.
An approximate sketch would change the requested scientific result and is not implemented.

## Diagnostics and verification

Shared annotation preparation emits INFO milestones for input reading, global SNP identity checks and chromosome preparation, and successful completion with retained counts and elapsed time.
A storage summary reports binary/continuous column counts and actual value-file bytes, including NPY headers but excluding metadata and scratch.
There is no per-chunk progress output.
Intentional gene-region exclusions are summarized once per source and role, combining line–gene pairs across audit chunks while retaining only one source's message.
CM/MAF notices appear once per source during its classification pass.
Complete row audits and existing failure footers are preserved; see [workflow logging](../../current/workflow-logging.md#annotation-preparation).

Annotation identity audits, gene row audits, and quantile alignment diagnostics have complete streamed output routes.
Most result bookkeeping uses compact summaries or persistent paths.
`FrameSpool`, `IdentityDropSpool`, `TsvDiagnostics`, and quantile `AlignmentDiagnostics` preserve the stage-specific record order without an all-chromosome diagnostic DataFrame.
Gene audit rows are replayed in bounded chunks, but coverage failures still collect complete affected-gene lists in summary cells and error strings.
Batch regression also retains the full legacy-sumstats projection drop table through fitting.
These are remaining diagnostic-retention gaps (item R3 of the September 10 completion review); complete records should remain available on disk while subsequent work retains counts, bounded examples, and paths.
Direct results include source-cleanup drops even for a chromosome whose annotation rows were all removed.

The focused contract tests are [`test_annotation_format.py`](../../../tests/test_annotation_format.py), [`test_annotation_storage.py`](../../../tests/test_annotation_storage.py), [`test_annotation_preparation.py`](../../../tests/test_annotation_preparation.py), [`test_annotate_streaming.py`](../../../tests/test_annotate_streaming.py), [`test_ldscore_streaming.py`](../../../tests/test_ldscore_streaming.py), [`test_reference_preparation.py`](../../../tests/test_reference_preparation.py), [`test_gene_index_streaming.py`](../../../tests/test_gene_index_streaming.py), [`test_regression_streaming.py`](../../../tests/test_regression_streaming.py), and [`test_quantile_streaming.py`](../../../tests/test_quantile_streaming.py).
They combine independent numerical/membership expectations with selected-read, allocation-shape, lifetime, publication, and failure tests.
Packed-format checks include literal byte expectations, non-byte-aligned writes, padding, empty/duplicate/reordered selections, bounded decoding, worker descriptor transfer, global exact classification, and file-free preparation.
On September 15, the full suite passed 1,812 tests and 132 subtests, with one skip; unittest discovery ran 1,000 tests successfully with one skip.
Both `ldsc --help` and `python -m ldsc --help` passed.
The [local benchmark harness](https://github.com/WenbinWu2001/ldsc_py3/blob/dc023e80a2a9caf9f2c7a35aa327e24e31c65087/benchmarks/annotation_memory.py) measures simultaneous process-tree peak RSS, elapsed time, peak private logical disk bytes, and persistent output size separately.
These measurements complement correctness and ownership tests; they do not replace them.

Peak memory, runtime, private disk, and persistent output size are reported separately; small inputs do not necessarily benefit in RSS or runtime.

The September 14 sequential-query benchmarks ran at revision `cf1dfac` on chromosome 22 with 1,000 gene lists, 53 binary baselines, 489 individuals, 141,123 reference SNPs, and 17,380 output SNPs.
Reducing generation width from 1,000 to 100 used 82.2% less peak private scratch and took 1.903 times as long, with essentially unchanged peak RSS (about 1.17 GiB).
All 17,380,000 query LD values per run matched the preserved pre-change width-1,000 result exactly, as did baseline tables, row order, counts, and overlaps.
The numerical checks are in [`test_ldscore_query_batches.py`](../../../tests/test_ldscore_query_batches.py), [`test_ldscore_streaming.py`](../../../tests/test_ldscore_streaming.py), [`test_gene_index_streaming.py`](../../../tests/test_gene_index_streaming.py), and [`test_regression_streaming.py`](../../../tests/test_regression_streaming.py).
Whole-genome RAM and runtime figures remain conditional projections; the older benchmark report predates this implementation.

### Recent optimization sequence and evidence

| Implemented change | Revision / evidence | Interpretation boundary |
| --- | --- | --- |
| Sequential execution/output query batches | `cf1dfac`; [`test_ldscore_query_batches.py`](../../../tests/test_ldscore_query_batches.py) and the benchmark above | Smaller width reduced scratch in the measured local workload but repeated traversal increased wall time; this predates later assembly changes. |
| Lossless binary annotation storage | `59ee02b`; [storage decisions](../../current/annotation-memory-decisions.md#confirmed-annotation-storage-scope-2026-09-15) and the size measurement above | Exact classification and decoding verified; 32-fold binary payload reduction is not a job-RSS measurement. |
| Indexed active-batch staging and streamed assembly | `659738a`; [`test_indexed_batch_staging.py`](../../../tests/test_indexed_batch_staging.py), [LD-score guide](../../wiki/main-functionalities/ldscore.md) | Bounds outstanding raw fragments and parent arrays; per-batch operator reloads add work. |
| Regression metadata/file reuse and fewer model copies | `7a478e6`; [`test_partitioned_publication.py`](../../../tests/test_partitioned_publication.py), [`test_regression_memory.py`](../../../tests/test_regression_memory.py), [regression guide](../../wiki/main-functionalities/partitioned-h2.md) | Measures and preserves full-model fits; neither changes the estimator nor eliminates worker-private matrices. |
| Direct file-writing streamed assembly | `3a4d729`; [`test_direct_batch_staging.py`](../../../tests/test_direct_batch_staging.py), [replay measurements](../../wiki/main-functionalities/ldscore.md#direct-assembly-measurements) | Local end-to-end parity is separate from precomputed full-row assembly replays, which do not measure genotype computation. |

The indexed, regression, and direct rows summarize the optimizations integrated during the October session.
Their detailed evidence remains attached to the corresponding implementation or workflow guide rather than being generalized into a universal speedup or RAM recommendation.
Matched runs must control query width, worker/native-thread counts, selected chromosomes/SNPs, input artifacts, and output semantics.
Do not extrapolate a short indexed run by multiplying its elapsed time by the number of batches: startup, validation, operator reading, compression, and publication have different scaling.

### Deliberate limits and approaches not adopted

- Preparation stays serial after the measured process/thread tradeoff described above.
- Projection keeps explicit float64 operands.
  A September 16 evaluation (NumPy 2.4.5, SciPy 1.17.1) found that both libraries multiply a float64 operator by a Boolean operand through a hidden full float64 copy on every call: about 16 MiB extra per full-budget tile, for dense and CSR operators alike.
  The nominal 87.5% operand saving is therefore not realized, and wide CSR projections used more temporary memory.
  Boolean operands also carry a dtype hazard: on an adversarial CSR row of 100,000 copies of 0.1, a float32 operator times a Boolean operand lost 1.44 against `math.fsum`, versus 1.88e-8 for float64.
  All 243 native-operand measurements with float64 operators matched the explicit-float64 results exactly.
  Splitting mixed binary/continuous groups failed the `rtol=atol=1e-14` check (maximum difference 7.96e-12), so it was not adopted either.
- Generated queries are batch-owned, but explicit `AnnotationBundle.read()` and `LDScoreSource.read_queries()` selections remain caller-owned materializations with no automatic width cap.
- Genome-wide baseline data, chromosome genotype payloads or index operators, model-private regression arrays, and the exact quantile target vector remain substantial allocations.
- Restriction caches and reusable pair-panel metadata intentionally trade retained memory for avoided parsing; they are not constant-memory streams.
- All-query fragment staging and genome-wide direct query concatenation are superseded for writing workflows; the small no-output API intentionally retains a materialized result.
- There is no batch resume/checkpoint mechanism, public bit-packed annotation format, approximate quantile method, or automatic worker tuning in these designs.
- Complete affected-gene error lists and the legacy sumstats projection-drop table remain diagnostic-retention costs.

### Maintenance and validation requirements

Before changing a buffer, identify who owns it, whether another process borrows it, its dtype/shape, and the event that releases it.
Bound producers as well as consumers: deleting consumed files is insufficient if all batches can be produced first.
Preserve float64 calculation and each producer's established query-filtering precision; do not replace exact variation summaries with extrema or sums that mishandle missing values.
Avoid trading one removed allocation for hidden per-tile dtype conversion, sparse slicing, or IPC copies.

Validate both scientific results and resource lifetimes.
Use independent small numerical expectations for LD, counts, overlaps, and jackknife results, and matched-layout exact comparisons when isolating transport/publication changes.
Include negative R², one diagonal, contributor-only SNPs, mixed binary/continuous columns, reordered identities, all-skipped batches, cross-chromosome variation, and dtype promotion.
Exercise worker, reader, and writer failures to prove shutdown precedes scratch cleanup and unsuccessful computation leaves no published scientific result.
Record numerical tolerances, schemas, row/column order, diagnostics, and relocation behavior alongside measurements.

Useful entry tests beyond the annotation/regression list above are [`test_plink_bed_reader.py`](../../../tests/test_plink_bed_reader.py), [`test_ref_panel_builder.py`](../../../tests/test_ref_panel_builder.py), [`test_ldscore_projection.py`](../../../tests/test_ldscore_projection.py), [`test_direct_batch_staging.py`](../../../tests/test_direct_batch_staging.py), [`test_partitioned_query_workers.py`](../../../tests/test_partitioned_query_workers.py), and [`test_worker_policy.py`](../../../tests/test_worker_policy.py).
For performance claims, record package revision and dirty state, dependencies, input dimensions, concurrency/native limits, cache conditions, repeated matched wall/CPU measurements, anonymous/cache peaks, and peak scratch occupancy.
Distinguish local kernel tests, local integration, assembly-only replays, and full HPC runs.
Documentation-only updates do not establish a new performance result or justify reducing requested RAM.
