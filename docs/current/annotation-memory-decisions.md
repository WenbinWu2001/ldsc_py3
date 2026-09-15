# Annotation and Workflow Memory Decisions

Last updated on: 2026-09-15

Current status: the user confirmed the complete annotation-storage design and authorized implementation on September 15. Packed binary storage, exact automatic classification, two input content passes without numeric staging, and explicit private column selection are implemented. Public float32 reads, logical input order, and LD-score output behavior are preserved. The full suite passed 1,812 tests and 132 subtests, with one skip; an independent chromosome-22 baseline comparison verified every decoded value. See the [implemented format policy and measured file sizes](annotation-memory-design.md#annotation-format-policy). These decisions supersede the September 14 choice to retain dense annotation storage, while preserving the approved execution batching.

## Confirmed annotation storage scope (2026-09-15)

### Scope and scientific types

Optimize annotation handling across the package, including preparation and consumption in `annotate` and `ldscore`. The requested work is automatic binary detection, compact encoding during preparation, and bounded decoding for computation. The user accepts the current output-row LD-score accumulators and within-batch LD-score aggregation memory. Reducing their allocations, removing aggregation copies, and investigating completed-batch LD-table retention are outside this storage change.

Determine annotation type from processing logic and value evidence, independently of whether a column belongs to the baseline or query group. Gene-list and BED projection produce binary annotations by construction. For supplied annotation columns, use automatic classification with exact zero/one membership. Sampling and tolerance-based rounding do not establish that a supplied column is binary. Other numeric values must retain their numerical meaning.

Store binary annotation values as packed bits. Store continuous annotation values in dense matrices even when many entries are zero. Sparse annotation matrices, including SciPy CSC, are excluded from this change. Keep binary and continuous values in separate physical stores behind the shared selected-read interface; a continuous column must not force the complete binary annotation dataset into floating-point storage.

The user confirmed that the public Python `AnnotationBundle.read()` API continues returning float32 arrays, including binary-only selections. Packed storage and bounded decoding are internal details. The caller receives the requested numerical values in the established return dtype; a binary-only public read does not change to Boolean arithmetic. This decision accepts the larger returned tile while avoiding materialization beyond the explicitly requested selection.

### Logical order and physical correspondence

The logical annotation order is baseline columns followed by query columns, preserving resolved input order within each group. Gene-set queries follow their resolved declaration order. Existing filtering retains the relative order of surviving queries. Physical storage grouping must not change these logical sequences.

Reuse the existing lookup from annotation name to physical store and store-local column position. Every stored column must correspond to the correct annotation and the same ordered SNP metadata rows. Explicit reads restore the caller's requested column order across stores. Encoding and decoding must also preserve requested SNP row order and repeated row selections.

The user requested the simplest necessary ordering solution. Keep public default ordering at `AnnotationBundle.read()`, which already supplies `baseline_columns + query_columns`. Private multi-store shard reads take explicit ordered columns from their callers; migrate any remaining implicit calls instead of adding another logical column list to every shard or a new ordering abstraction. Physical store concatenation must not determine public read or output order. This replaces the earlier proposal to support logical-order defaults at the private shard level.

Downstream LD-score artifacts retain their current formats, floating-point policy, baseline/query grouping, and annotation-column order. This storage decision does not authorize reordering LD-score columns by annotation type or changing numerical accumulation. The existing metadata-column layout and query-batch manifest continue to govern output assembly.

### Packing and bounded access

Packed values must be accompanied by their logical row count and explicit bit order. Packing and unpacking must agree on the SNP axis, column correspondence, bit order, byte boundaries, and padding; padding bits are not SNP rows. Decode only the requested rows and columns through bounded reads. Apply compact encoding during preparation as well as in the final private annotation stores.

The user selected two content reads of supplied annotation files for simpler staging, accepting repeated parsing and decompression, and explicitly ruled out extra staging or temporary files for this optimization. The first bounded pass validates inputs and establishes column classifications and row metadata/selection information. The second bounded pass rereads the immutable original inputs and writes directly to the final private packed-binary and dense-continuous stores using those established identities and ordered columns. Remove the numeric `values.bin` staging for this route; do not replace it with a compact numeric spool, an on-disk conversion copy, or an additional temporary artifact stage. Reuse the existing metadata and identity bookkeeping where needed without introducing an extra staging layer. Generated gene-list/BED annotations are already binary by construction and can be packed directly when projected; the two-pass choice concerns supplied annotation tables.

### Pre-implementation code evidence

The ordering and storage seams were checked locally on `restructure` at `19381df` before recording these decisions:

- [`prepare_annotation_sources()`](../../src/ldsc/_annotation_sources.py) retains source/header order within each annotation group and forms `baseline_columns + query_columns`. Chromosome-sharded sources must have matching ordered headers.
- [`resolve_gene_lists_staged()`](../../src/ldsc/_gene_query_storage.py) creates declaration ordinals from resolved source order. [`gene_query_statuses()`](../../src/ldsc/query_annotations.py) and [`build_query_shards()`](../../src/ldsc/_annotation_queries.py) propagate the ordered usable queries. A gene control is appended to the baseline group under the existing policy.
- [`AnnotationBundle.read()`](../../src/ldsc/_annotation_bundle.py) already supplied the logical baseline-plus-query order for default reads. [`AnnotationShard.read()`](../../src/ldsc/_annotation_storage.py) restored explicitly requested column order across stores, but its former default used `AnnotationShard.columns`, which concatenates physical store columns. Public bundle reads, [`MappedAnnotations.read()`](../../src/ldsc/_kernel/ldscore_projection.py), and gene-index baseline preparation already passed explicit column lists to the shard. Implementation preserves that simple pattern and removes the private physical-order default.
- [`_split_ldscore_table()`](../../src/ldsc/ldscore_calculator.py) constructs baseline and query outputs from their supplied ordered name lists. [`AnnotationDirectoryWriter.write()`](../../src/ldsc/outputs.py) requests and writes generated queries in `bundle.query_columns` order.
- [`_scan()` and `_write_values()`](../../src/ldsc/_annotation_sources.py) formerly staged supplied values as float32 before writing dense column-major files. [`build_query_shards()`](../../src/ldsc/_annotation_queries.py) likewise wrote generated binary vectors through the default float32 store. The approved implementation replaces those paths with direct compact storage.
- Other annotation routes include the file-free [`AnnotationBundle.from_frames()` preparation](../../src/ldsc/_annotation_memory.py), [saved annotation reload](../../src/ldsc/_annotation_outputs.py), [gene-index baseline preparation](../../src/ldsc/gene_ldscore_index.py), and [quantile annotation reads](../../src/ldsc/_quantile_inputs.py). External quantile target vectors have their own float64/raw-missing-value contract, documented in [the current memory design](annotation-memory-design.md#exact-global-quantiles); sharing a low-level storage primitive does not make those targets ordinary annotation columns.

## Earlier batch implementation status (2026-09-14)

Status: on September 14, the user selected sequential query execution batches with separate query Parquet files and repeated genotype traversal. This replaces the proposed single-traversal streaming redesign and its planned second annotation-storage session. The subsequent integration interview confirmed the output layout, removal of old-directory compatibility requirements, saved-result handles, private publication, bounded indexed assembly, strictly write-free small Python calculations, and unrestricted explicit query-read widths below. The user then requested chromosome-parallel indexed generation using `--threads`, matching direct mode; the indexed ownership contract is now one operator per active worker. The user confirmed the final design and chromosome-count cap, then authorized implementation. The implementation and its verification follow the contracts below. The September 10 work is already implemented with gaps in the [completion review](../audits/annotation-memory/completion-review.md). Its [implemented design](annotation-memory-design.md), [specification](../specs/2026-09-10-annotation-workflow-memory-design.md), and [plan](../plans/2026-09-10-annotation-workflow-memory.md) describe existing behavior and historical scope, not a requirement to revive the abandoned optimization effort.

## Confirmed sequential query batches (2026-09-14)

The target workload is 18,000 queries across chromosomes 1–22, with 53 baseline annotations, query batch size 1,000, and 22 chromosome workers. Each query retains its separate baseline-plus-one-query regression model. After evaluating implementation complexity, memory, and runtime, the user explicitly chose the simpler batch strategy and cancelled the earlier broad optimization request.

The selected design is:

1. Execute direct query batches sequentially. For the target workload, this means 18 batches of 1,000 queries. The batch owns the lifecycle `prepare -> compute -> write -> release`; chromosome workers may run concurrently within that batch. Do not retain completed batch matrices or launch all query batches concurrently. Indexed generation has the confirmed chromosome/operator lifecycle below while preserving bounded query payloads.
2. Write query LD scores into separate Parquet files when there is more than one batch. Each file contains its query columns across the selected chromosomes in deterministic genomic order. This supersedes the single aggregate query-file requirement; it does not introduce public chromosome-sharded outputs. Preserve the shared baseline and scientific metadata. The integration interview confirmed the filenames and ordered batch manifest below.
3. Repeat genotype traversal and correlation computation for each query batch. Repeated reference work is an accepted runtime trade-off for simpler implementation and lower memory. The prior requirement for one numerical traversal across all queries is superseded. Preserve supported backends and migrate affected consumers without broad backend redesign.
4. Reuse the current numerical kernels and permit in-memory output-row-by-active-batch accumulators and completed tables for that batch. The memory target is bounded in total query count at fixed batch size and worker count, not independent of SNP count. Disk-backed numerical accumulators, active-window finalization, and fine-grained row streaming are not requirements of this design.
5. Keep completed batches file-backed, with explicit selected-query reads through an extended reader such as [`LDScoreSource`](../../src/ldsc/ldscore_source.py). Returning or opening the complete result must not reconstruct every query table. Materialized low-level per-batch result objects can remain; a wholesale public-class replacement is not required.
6. Prepare and release query annotations within the batch lifecycle to reduce peak scratch where the input route permits it. The September 14 decision retained the dense annotation representation and required no follow-up storage redesign; that representation choice is superseded by the [September 15 annotation storage scope](#confirmed-annotation-storage-scope-2026-09-15). The batch lifecycle remains applicable. Lower peak scratch does not imply fewer total annotation bytes written across all batches.
7. No hard RAM ceiling, temporary-disk limit, or maximum slowdown has been set. Measure process-tree memory, elapsed time, private scratch, and persistent output separately. Additional aggregation cleanup should follow a measured need; closing every historical audit item is not automatically part of this narrower task.

The current [`ProjectionAccumulator`](../../src/ldsc/_kernel/ldscore_projection.py) allocates for every column supplied to it, and [`LDScoreCalculator`](../../src/ldsc/ldscore_calculator.py) materializes and aggregates those supplied columns. The new outer batch boundary can reuse those mechanisms at a narrower width. Preserve reference contributors, regression-weight policy, numerical precision and tolerances, global identity/annotation validation, counts, overlaps, diagnostics, and failure/overwrite semantics. Do not undo already implemented annotation or regression improvements. The old nine-question streaming interview and two-session sequence are superseded.

## Confirmed batch integration (2026-09-14)

These decisions record the user's answers to both integration rounds and the subsequent request for indexed chromosome parallelism. They supersede conflicting compatibility or interface recommendations in the handoff and older sections; the user confirmed implementation on September 14. The current implementation follows these contracts.

### Execution batches and saved query files

Use the existing `query_batch_size` / `--query-batch-size`, with default `1000`, for the execution width and output query grouping. Preserve existing query order. A single batch uses `ldscore.query.parquet`; multiple batches use `ldscore.query.batch00001.parquet`, `ldscore.query.batch00002.parquet`, and subsequent deterministic ordinals. Smaller batches therefore increase both the number of query files and repeated direct computation.

Root `metadata.json` includes an ordered `query_batches` list recording each batch file, its ordered query columns, and its chromosome row groups. Keep one shared baseline file and the existing shared count/overlap semantics. The current file and metadata writer is [`LDScoreDirectoryWriter`](../../src/ldsc/outputs.py), `write()` and `build_metadata()`.

The user explicitly removed the requirement to read old LD-score directories: all affected writers and readers should use the latest canonical directory contract, following the repository's policy of replacing obsolete interfaces directly. Do not add an old-directory compatibility fallback. This decision concerns LD-score result directories; it does not retire supported PLINK/R2 inputs, immutable prebuilt gene indexes, or the explicit LDSC2 conversion workflow. Scripts that assume one query filename must migrate to the batch manifest when consuming multi-batch results.

### Saved and in-memory results

Calls that write return an object such as an extended [`LDScoreSource`](../../src/ldsc/ldscore_source.py) identifying saved artifacts and providing explicit selected-query access. This saved-result convention also applies to written baseline-only and single-batch results. Returning or opening it must not reconstruct all query values or retain completed chromosome/query tables.

A calculation with one execution batch may return an in-memory result without writing files, under the strict Python policy below. More than one execution batch requires an output directory. Materialized low-level batch results can remain.

Explicit `read_queries(names)` calls may request any number of query columns, including selections spanning several batch files or exceeding the generation batch width. Return exactly the requested columns in the requested order and retain no query-value cache. Do not enforce a read-width limit or automatically split the call into a different return type; callers accept the memory cost of their requested selection. Package regression workflows continue to use their own bounded query batches. The current selective-read seam is [`LDScoreSource.read_queries()`](../../src/ldsc/ldscore_source.py).

### Python calculations with no filesystem writes

The intended use is a small, single-batch Python calculation with prepared inputs. The user explicitly selected zero filesystem writes, including temporary files. Keep diagnostics in memory or return them directly, use already-prepared in-memory annotations, and add preparation that avoids disk staging where necessary for this route. The user accepts the additional RAM needed for data that otherwise lives in temporary files. This is a scoped exception to the large-run annotation and diagnostic storage rules below; it does not relax the output-directory requirement for multiple execution batches.

Omitting saved LD-score artifacts while still creating temporary diagnostics does not satisfy this policy. The calculation must not create scratch directories, staged annotation/diagnostic files, workflow file logs, output artifacts, or failure markers. Preserve scientific validation and diagnostic contents in memory, including failures; do not silently fall back to disk. Return the materialized scientific result and its diagnostics without depending on a temporary filesystem owner.

The current implementation requires an actual storage change: [`AnnotationBundle`](../../src/ldsc/_annotation_bundle.py) currently requires an `AnnotationWorkspace` and disk-backed shard descriptors, and [`_stage_chromosome_drops()`](../../src/ldsc/ldscore_calculator.py) writes a diagnostic file even without output configuration. Merely passing `output_config=None` is insufficient. The in-memory input/preparation route and diagnostic handling must be verified together through the Python calculation boundary.

### Publication and failure

Write batch files privately inside the supplied output directory and release their memory as they are completed. After every batch succeeds, move the files to their final names and write final root metadata last. Downstream workflows consume the complete published run.

On a handled computation failure, remove the current run's private batch files and retain diagnostics. Preserve the existing overwrite and failure-marker policy. Final publication carries no new rollback guarantee. The existing writer and artifact-family lifecycle are in [`outputs.py`](../../src/ldsc/outputs.py), `LDScoreDirectoryWriter` and `ArtifactFamily`; the current failure-marker contract is recorded in [path specification](path-specification.md).

### Indexed generation

Include bounded-memory result generation from a prebuilt gene LD-score index. Updating only filenames and the returned handle would leave the existing all-query `query_scores` allocation and retained `query_tables` in [`_run_indexed_ldscore_in_workspace()`](../../src/ldsc/gene_ldscore_index.py).

Support indexed chromosome parallelism through `ldscore --gene-ldscore-index-dir ... --threads N` and the corresponding `threads` argument on `run_indexed_ldscore()`. Use the same worker-count semantics as direct mode: default `1` runs sequentially, positive integers request that many chromosome workers, `-1` uses available cores, `-2` leaves one free, zero is invalid, and the effective count is capped at the chromosome count. This runtime control does not change immutable scientific inputs stored in the index. The current indexed dispatcher explicitly rejects `--threads`; remove that rejection as part of implementation. See [`_run_explicit_indexed_ldscore()`](../../src/ldsc/ldscore_calculator.py), the direct `_run_chromosomes()` process-pool boundary, and `_resolve_worker_count()`.

Each active worker owns one chromosome operator and processes that chromosome's query batches sequentially. Write its active batch to private disk storage and release the batch values before advancing; retain the operator until that chromosome's batches finish, then release it before taking another chromosome. With multiple workers, several chromosome operators and active query batches may coexist, bounded by the resolved worker count. Memory therefore grows with worker count as well as operator size, SNP rows, and query batch width; the guarantee is one operator per active worker, not one operator for the entire run.

Keep active and queued chromosome work bounded and pass paths or compact descriptors across the worker boundary rather than preloading every operator. Workers write distinct private chromosome/batch fragments. Assemble the same final query-batch files in deterministic chromosome and query order without reconstructing an all-query score table; completion order must not change scientific outputs. Preserve the approved complete-run publication and handled-failure cleanup policy. Existing exhaustive index validation remains a prerequisite; adding assembly workers does not weaken those gates or change index construction.

This adds indexed worker orchestration, writer integration, and corresponding lifetime tests to the approved work. Verification must compare one versus several chromosome workers, check the concurrency bound and operator/batch release, preserve numerical results and ordering across batch widths, and exercise worker failure before publication. No indexed parallel implementation or benchmark has run yet.

### Potential future feature: resume from a completed batch

The user requested documenting support for resuming from the last completed batch as a potential future feature. It is deferred and is not part of this implementation. Current retries start fresh under the existing scratch and overwrite rules; this note does not authorize checkpoint reuse, retention of handled-failure scratch, automatic recovery, or a new persistent cache. A later design would need to define what constitutes a reusable completed batch, particularly for indexed generation's chromosome-outer lifecycle.

## Intended use

This optimization is designed for testing many pathways in one run, for example enrichment testing for 1,000 pathways. LD-score generation shares reference-LD work within a query batch and may repeat it for later batches; subsequent batch partitioned heritability fits each pathway separately against the baseline categories. Query batching never turns the pathways into one joint regression model. README, current workflow documentation, examples, and the user wiki describe this use case.

## Evidence baseline

The active repository is the local `ldsc_py3_restructured` checkout on `restructure`. The pre-refactor evidence revision was `a505c45d61d35b39a3fc086c2b623dd7d6936d96` (`refactor: remove annotation fingerprints`). Unrelated working-tree changes are outside this task.

That commit removes annotation fingerprint generation, retained common-reference annotation values, result fields, and written fingerprint metadata. It preserves annotation-name validation and advisory value classification. Quantile reconstruction now checks alignment, recorded universe sizes, common-SNP annotation sums, and available overlap cross-products; these aggregate checks do not prove the original annotation value at every SNP. Independent artifact identity hashes remain. See [annotation semantics](../../src/ldsc/annotation_semantics.py), [LD-score calculation](../../src/ldsc/ldscore_calculator.py), and [quantile reconstruction and provenance checks](continuous-annotation-quantile-h2.md#reconstruction-and-provenance-checks).

At the pre-refactor evidence revision, the remaining eager paths included `AnnotationBuilder._detect_chromosome_shards()` and `parse_annotation_file()` reading complete source tables, `_run_sharded_inputs()` concatenating complete chromosome bundles, and `LDScoreCalculator._run_chromosomes()` eagerly slicing/submitting chromosome inputs. See [annotation builder](../../src/ldsc/annotation_builder.py) and [calculator](../../src/ldsc/ldscore_calculator.py). Those eager paths are now replaced as described in [the implemented design](annotation-memory-design.md); fingerprint aggregation remains removed.

## Confirmed scope and compatibility

The following records the earlier workflow coverage and compatibility contracts. Preserve existing supported behavior; these areas are not a fresh implementation checklist for the narrower sequential-query-batch task above.

- Optimize standalone `annotate` and direct `ldscore`, including input preflight, chromosome discovery, annotation loading, alignment, projection, computation, and writing.
- Cover chromosome-sharded and whole-genome annotation files, BED queries, gene-list queries, and synthetic `base`, with both PLINK and parquet-R2 reference backends where those routes apply. Standalone `annotate` will accept exactly one BED or focal gene-list query route, require baseline annotations and an output directory, and exclude control-gene lists. This functional expansion was approved separately; its validation boundary, condition–outcome table, diagnostics, and integration ownership are recorded in [standalone annotate gene-list decisions](annotate-gene-list-decisions.md).
- Backward compatibility with the existing Python API is not required. Replace obsolete interfaces directly and update affected internal callers, tests, and examples. This does not authorize unrelated interface changes.
- The September 14 integration interview also removes backward compatibility requirements for old LD-score result directories. All affected producers and consumers move to the latest batch-manifest contract.
- Preserve numerical results, SNP identities and ordering, filtering, annotation validation, counts, overlap semantics, and the reference-versus-regression SNP universes. Existing inconsistencies that would require a scientific behavior change must be surfaced separately.
- Annotation fingerprints remain removed under the policy in `a505c45`; the refactor does not restore their retained matrices or introduce replacement fingerprint staging.
- Implementation remains local. Do not launch a full 1,000-query HPC run for validation.
- The user subsequently expanded the optimization to batch `partitioned-h2` preparation and query loading, exact streaming `quantile-h2` statistics, gene-index loading and indexed LD-score assembly, incremental regression detail writing, and bounded shared metadata/diagnostics. These areas are now in scope; the earlier downstream-regression exclusion and proposed quantile materialization exception are superseded.
- Direct `ldscore` has query batching and score accumulation restricted to resolved output rows for both PLINK and parquet-R2. The September 14 decision permits repeating the numerical traversal for successive execution batches.

## Confirmed dataset contract

`AnnotationBundle` remains one object representing the complete annotation dataset. Its data are organized as chromosome shards loaded or constructed on demand. Shared metadata consists of chromosome scope, ordered baseline/query column names, configuration, provenance, and the information needed to access each chromosome. SNP row metadata belongs with the corresponding chromosome's annotation values.

For large, staged workflows, the bundle must not hold all chromosome DataFrames in a dictionary or list. Normal execution must not automatically cache every shard, concatenate shards into whole-genome annotation matrices, or rematerialize those matrices when returning results. The subsequently approved [small Python calculation](#python-calculations-with-no-filesystem-writes) may instead hold its prepared annotation inputs and diagnostics in memory; the user explicitly accepts that RAM cost to avoid all filesystem writes.

For sequential execution, the lifecycle is `load/build chromosome -> consume or write its shard -> release chromosome data -> process the next chromosome`. For parallel execution, active annotation data and queued work are bounded by the configured worker count; inputs must not be eagerly loaded or serialized for every chromosome.

Release means that chromosome-owned annotations, reader state, temporary compute buffers, mappings, and retaining views are no longer live before the next chromosome loads. Python/native allocators may keep freed pages available for reuse, so operating-system RSS need not fall immediately. Sequential query execution may retain completed chromosome LD tables for the active query batch; those tables must be released before advancing to the next query batch. Explicitly permitted shared state can remain live.

Standalone `annotate` writes chromosome shards incrementally and returns an object referencing those outputs. Final aggregation retains necessary shared metadata and diagnostics rather than complete annotation matrices.

`AnnotationBundle` is a resource-owning handle supporting context-manager use and explicit `close()`. It owns private staging. A calculator receiving an existing bundle borrows it without closing it. Closing removes owned private storage and prevents further shard loading, while leaving original inputs and persistent outputs untouched. Returning a bundle does not close resources that the returned handle still needs.

The returned complete annotation dataset may depend on original baseline sources. Generated query shards remain persistent; no portable copy of all baseline annotations is required. Private baseline staging remains owned by the live bundle.

New annotation/analysis disk staging stays inside the explicitly supplied output directory, rather than using the system temporary directory or `TMPDIR`. Workflows that require disk staging require an explicit output directory. The approved small, single-batch Python route instead holds preparation and diagnostics in memory and performs no filesystem writes; pure numerical functions operating on supplied arrays also remain file-free. Owned scratch in writing workflows is removed on closure or handled failure. Interrupted writing runs start fresh, without automatic reuse or deletion of abandoned directories.

The existing gene-index construction mechanism is an explicitly approved exception: retain its established sibling staging, build-state, locking, and atomic whole-directory publication behavior. This exception does not authorize new annotation or downstream-analysis scratch outside the supplied output directory.

## Confirmed duplicate identity policy

Apply the existing drop-all annotation identity rules across the complete logical dataset, independent of chunk boundaries and input layout. This includes cross-chromosome rsID-family collisions according to the selected identity mode. Corresponding rows in aligned baseline/query column sources represent one logical SNP row; they are not duplicate observations merely because several sources describe them.

This deliberately resolves an existing input-layout difference. A probe of the current `_apply_identity_cleanup()` retained two occurrences of the same rsID on different chromosomes when cleaned separately, but dropped both when supplied together. `_run_sharded_inputs()` currently uses separate cleanup. The user explicitly approved consistent dataset-wide cleanup, including the changed outcome for this malformed sharded-input edge case. See [identity cleanup](../../src/ldsc/_kernel/snp_identity.py) (`base_key_series`, `clean_identity_artifact_table`) and [annotation orchestration](../../src/ldsc/annotation_builder.py) (`_apply_identity_cleanup`, `_run_sharded_inputs`). This decision does not relax integrity validation of immutable index or canonical result artifacts.

## Confirmed whole-genome input strategy

The memory bound applies before numerical computation, including preflight and chromosome discovery. Reading a complete input and selecting one chromosome afterward does not meet the requirement.

Scan whole-genome annotation files in bounded chunks and, where needed, prepare private chromosome shards for subsequent access. Avoid rescanning and decompressing a whole-genome source once per chromosome. Use existing chromosome-sharded sources directly where practical. Chunking must preserve cross-chunk validation, duplicate handling, row alignment, and existing stage-specific validation requirements.

Shard access must also support selected query columns and needed SNP row ranges/indices. A query-batch request must not first load a complete dense chromosome-by-all-query matrix. Reuse normalized data and interval projections instead of reparsing whole-genome text or recomputing the same projections for every numerical LD block.

## Confirmed direct LD-score batching and output-row accumulation

Expose `--query-batch-size` in direct `ldscore`, indexed `ldscore`, and batch `partitioned-h2`, with a default of `1000`. It limits active query-column batches and is distinct from chromosome worker count, SNP block size, index-construction atom batches, and the number of models fitted together. A run with exactly 1,000 queries can use one batch at the default; smaller values permit a lower active query width.

Within one query execution batch, reuse each correlation block or decoded R2 chunk for baseline/control columns, the regression-weight mask, and that batch's query columns. The September 14 design deliberately permits wrapping successive query batches around the chromosome kernel, resetting genotype traversal, and rereading R2 inputs. Preserve the existing numerical projection and contributor rules within each traversal.

Batch annotation counts, advisory classification, baseline-query overlap, and query self-overlap too. Keep shared baseline statistics and the existing baseline-rows-plus-query-diagonal overlap representation. Do not introduce cross-query products or full-width float64 annotation/statistics temporaries.

Resolve the existing output/regression SNP selection and region exclusions before allocating scores. Accumulate only those output rows. Keep the supported all-output behavior when no output restriction is requested; do not hard-code HM3 into the numerical kernel. In-memory output-row-by-active-query-batch accumulators are permitted. Do not allocate scores for every query across all execution batches, or accumulate full-reference rows only to slice to output rows afterward.

For each output SNP, annotation contributions still come from the complete retained reference universe. Non-output endpoints must contribute to output endpoints in symmetric pair accumulation. Counts and overlap retain all/common reference universes. `w_ld` keeps its separate filtered-regression contributor universe, diagonal, and exclusion policy. See the [SNP-universe and traversal contract](ldscore-snp-universe-contract.md).

Preserve current accumulation precision and applicable numerical tolerances. No reduced-precision policy is approved. At fixed batch width, query-related accumulator and result RAM must not grow with the number of completed execution batches; total query count may increase persistent output and compact summary metadata.

## Confirmed result and resource boundary

The September 14 design requires completed query batches in writing workflows to remain file-backed with explicit selected-query access. It permits materialized tables for the active batch and does not require replacing every low-level `LDScoreResult` interface. A complete saved-result handle must not retain or reconstruct all query batches. Preserve the existing bounded annotation access and release each completed batch. The small Python exception returns one materialized batch with in-memory diagnostics and no filesystem writes. Explicit query reads may exceed generation batch width at the caller's memory cost.

Public LD-score files remain aggregated across chromosomes. Keep the shared `ldscore.baseline.parquet`, metadata, and overlap artifacts; query LD scores are split by query execution batch into multiple Parquet files when needed. The integration interview confirmed `ldscore.query.parquet` for one batch, numbered query files for multiple batches, and an ordered manifest with chromosome row groups. Old-directory readability is no longer required. See [confirmed batch integration](#confirmed-batch-integration-2026-09-14).

Require exact identities, ordering, inclusion decisions, annotation values, and quantile membership except for the explicitly approved duplicate-policy correction. Streaming reductions may differ in final floating-point digits within the applicable existing tolerances. Keep scientific validation thresholds unchanged, including quantile aggregate validation at `rtol=1e-6`, `atol=1e-8`; report observed differences instead of weakening checks to accommodate a regression.

The user has no hard RAM ceiling, runtime constraint, or maximum acceptable slowdown. Measure peak memory, runtime, and temporary disk usage separately across the expanded workflows; distinguish measurements from extrapolations and explain remaining bottlenecks. At the previously cited dimensions, one float32 HM3 matrix with 1,053 annotations contains approximately 4.45 GiB of values; this is an array-size estimate, not a measured process peak. The earlier approximately 23.4 GiB common-annotation matrix estimate no longer establishes a fingerprint storage requirement after `a505c45`.

## Confirmed downstream workflow scope

### Batch partitioned heritability

Prepare SNP alignment, trait data, and shared baseline LD scores once. Load one query or a bounded query batch on demand instead of retaining every query column. Each model fits its complete retained genome-wide SNP set. Preserve model-specific filtering, regression weights, and jackknife calculations; do not share quantities that depend on the fitted model. At the evidence revision, repeated assembly and detail accumulation occurred in [regression_runner.py](../../src/ldsc/regression_runner.py), `RegressionRunner.build_dataset()` and `estimate_partitioned_h2_batch()`; that directory loader read complete query tables in `load_ldscore_from_dir()`.

### Exact global quantile statistics

Determine the validated common SNP universe and the existing exact quantile boundaries from eligible target values. Process annotation shards using those same global boundaries, retaining running quantile sums, overlap-validation products, and numerically stable standard-deviation statistics instead of every SNP's fitted annotation values. Preserve ties, missing-target policy, aggregate validation, and full-common-universe statistics for standardized coefficients. Do not build whole-genome fitted-annotation matrices or dense SNP-by-quantile indicator matrices. The current numerical definitions and eager assembly are in [quantile_h2.py](../../src/ldsc/quantile_h2.py), `assign_legacy_quantiles()`, `_prepare_quantile_inputs()`, and `run_quantile_h2_from_args()`.

The exact boundary-selection stage may retain one genome-wide eligible numeric target vector and its sorting workspace. Release that vector when boundaries are established and it is no longer needed; fitted annotations remain chromosome-scoped. At 5,961,159 values, one float64 vector occupies approximately 45.5 MiB, excluding sorting workspace and other state.

### Gene-index loading and indexed assembly

Represent a loaded index using shared metadata and chromosome component paths. Validate with bounded memory. During assembly, each active chromosome worker loads one sparse operator, processes that chromosome's query batches sequentially, and releases the operator before taking another chromosome. Keep the operator loaded throughout its batches instead of rereading it for each batch. The September 14 integration decisions require writing and releasing each active query payload, removing all-query score allocations and retained query-table collections, and supporting `--threads` with direct-mode chromosome-worker semantics. Preserve index integrity checks and construction staging. Smaller multiplications, separate validation passes, and concurrent I/O can change throughput; measure rather than promise unchanged runtime. The relevant current seams are [gene_ldscore_index.py](../../src/ldsc/gene_ldscore_index.py), `LoadedGeneLDScoreIndex`, `_load_gene_ldscore_index()`, and `run_indexed_ldscore()`.

### Incremental regression detail writing

Write per-query category tables, coefficient delete values, and associated metadata after each fit completes. Retain aggregate summary rows and compact bookkeeping instead of all detailed results. Preserve output contents, ordering, diagnostics, and existing failure/overwrite contracts. Returned objects reference valid persistent outputs. Current [output writing](../../src/ldsc/outputs.py), `PartitionedH2DirectoryWriter._query_records()`, assigns folder ordinals from the final sorted summary; preserving this order requires final path assignment after summary ordering is known, even if detail contents are written earlier.

The September 14 decision kept incremental detail writes private until the complete batch succeeds. The September 15 query-failure decision preserves this strict default and adds explicit `--continue-on-query-error`: per-query model preparation, regression, and result calculation exceptions are logged and recorded, and successful fits publish after the scan finishes. Shared input/output errors and scans with no successful fits remain fatal. Canonical folders follow the final successful summary order; there is no early public access, resumable fitting, numerical fallback, or stronger rollback guarantee. See [query failure handling](partitioned-h2-results.md#query-failures-and-continuation).

### Shared metadata and diagnostics

Give chromosome metadata caches explicit ownership and lifetimes, releasing data when the consuming stage finishes. In writing workflows, stream potentially large row-level diagnostics to their output files while retaining summary counts and representative examples. Preserve complete diagnostic records and stage-specific validation behavior. This includes annotation identity drops, SNP alignment diagnostics, and gene-list row audits; returning saved results must not reload all diagnostic records. The approved small Python calculation keeps diagnostics in memory and returns them directly. Downstream commands load only needed LD-score columns; written generation results have the confirmed file-backed return requirement above.

## Decision map

The user confirmed the original map on September 10 and selected sequential query execution and output batches on September 14, then answered both integration rounds and requested indexed chromosome parallelism. The table records preserved existing contracts and the new batch boundary; it does not revive the abandoned broad optimization scope. The user confirmed the consolidated design and chromosome-count cap on September 14 and authorized implementation. Sequential direct batches, private publication, prepared-input zero-write calculation, and chromosome-parallel indexed assembly are implemented; validation evidence is recorded with the current memory design.

| Branch | Decision and boundary | State |
| --- | --- | --- |
| Objective and workflow scope | Many independent pathway fits; direct sequential query batches, bounded parallel indexed result generation, output/reader integration, and necessary consumer migrations | September 14 scope replaces the broad new optimization effort |
| Standalone gene-list behavior | Exact input route, catalog/build/padding/exclusion rules, Gate A/B outcomes, annotation-grid support, output reloadability | Confirmed; separate condition–outcome table |
| Dataset interface | One lazy complete-dataset handle; chromosome/column/row access; no whole-genome annotation materialization | Confirmed |
| Ownership and lifetime | Explicit owner/borrower/close; permitted original baseline dependencies; sequential release and bounded worker submissions | Confirmed |
| Scratch and interrupted runs | New scratch inside output directory; owned cleanup; fresh retry without abandoned-stage reuse/deletion; existing index-construction exception | Confirmed |
| Input preparation and identity | Bounded preflight/discovery and source preparation; dataset-wide annotation duplicate cleanup; stage-specific validation preserved | Confirmed |
| Direct computation | Sequential execution batches, target width 1000; repeated traversal and materialized active-batch scores permitted; contributor and weight policies preserved | Confirmed September 14 |
| Public LD output | Shared baseline; single or numbered query files; ordered batch manifest and chromosome row groups; latest contract without old-directory fallback | Confirmed September 14 integration round 1 |
| Result interface | Saved-artifact handle when writing; small single-batch Python calculation with prepared inputs, in-memory diagnostics, and no filesystem writes; multiple batches require output directory | Confirmed September 14 integration rounds |
| Explicit query reads | Requested columns may span files or exceed generation batch width; preserve requested order; no cache or enforced width limit; caller owns RAM cost | Confirmed September 14 integration round 2 |
| LD-score publication | Private batch files; release values after writing; final names and metadata after all batches succeed; existing failure/overwrite policy | Confirmed September 14 integration round 1 |
| Indexed assembly | `--threads` matches direct chromosome-worker semantics; one operator and one active query batch per worker; private fragments and deterministic final assembly; construction integrity/precision preserved | Confirmed September 14 integration round 1 and subsequent parallelism request |
| Batch regression | Shared preparation; each pathway fits its own full retained genome-wide model; per-model filtering/weights/jackknife preserved | Confirmed |
| Quantile statistics | Exact global boundaries; one numeric target vector/sort workspace permitted; sufficient statistics from shards | Confirmed |
| Regression publication | Private per-fit detail writing; final sorted canonical publication; existing failure/overwrite behavior | Confirmed |
| Diagnostics and shared state | Complete row records streamed; summaries/examples and explicit metadata lifetimes | Confirmed |
| Scientific equivalence | Exact identities/order/inclusion/values/membership; existing floating tolerances; no fingerprints restored | Confirmed |
| Implementation structure | Column-major private storage, disk-backed global identity bookkeeping, detached reads, explicit prepared/result split and accumulator owners | Confirmed |
| Evidence and delivery | Local numerical/artifact/lifetime checks and separate memory/runtime/disk measurements; post-implementation pathway-focused docs/wiki | Confirmed; validation details in plan |

## Closure audit and implementation checkpoints

The design-stage source audit identified preservation details, rather than new behavioral decisions. The following references describe the evidence revision; current replacements are documented in [annotation-memory-design.md](annotation-memory-design.md). At the evidence revision, gene-list resolution retained a combined audit and all per-source gene selections in `GeneListBatchResolution`; its replacement must stage row-level resolution and selection data while sharing catalog lookup state and compact source summaries. Controls remain the baseline column `gene_control` in direct/indexed workflows. LD-directory schema/identity checks remain separate from model-specific numeric checks and query-row alignment. Quantile external targets preserve float64/raw-token parsing and fatal duplicate-identity checks. See [gene resolution](../../src/ldsc/gene_list_resolver.py), `resolve_gene_lists()` and `GeneListBatchResolution`; [annotation construction](../../src/ldsc/annotation_builder.py), `_run_single_universe()`; [regression loading and assembly](../../src/ldsc/regression_runner.py), `load_ldscore_from_dir()` and `_assemble_regression_ldscore_table()`; and [quantile input preparation](../../src/ldsc/quantile_h2.py), `_read_target_annotation()` and `_prepare_quantile_inputs()`.

The September 10 interview closed the original design, and its subsequent implementation and completion gaps are recorded in the linked plan and audit. The user then rejected the proposed broad streaming redesign and approved sequential query batches. Existing stage-specific gates, immutable-input assumptions, diagnostic contents, canonical ordering within each output, and numerical tolerances remain fixed. The integration interview is closed. Implementation and validation now follow the confirmed contracts above. Do not resume the abandoned interview or make unrelated historical audit repair a prerequisite without a concrete interaction with this change.
