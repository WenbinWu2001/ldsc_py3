# Annotation and Workflow Memory Decisions

Last updated on: 2026-09-10

Status: the complete decision map and consolidated design, including the internal implementation choices, were confirmed on 2026-09-10. The refactor is implemented; final verification evidence is tracked by the plan. See the [implemented memory design](annotation-memory-design.md) for current module interfaces and ownership. The [specification](../specs/2026-09-10-annotation-workflow-memory-design.md) records the agreed contracts, and the [refactor plan](../plans/2026-09-10-annotation-workflow-memory.md) tracks delivery and validation. Private incremental regression-detail writing followed by final sorted publication is approved.

## Intended use

This optimization is designed for testing many pathways in one run, for example enrichment testing for 1,000 pathways. LD-score generation shares reference-LD work across pathways; subsequent batch partitioned heritability fits each pathway separately against the baseline categories. Query batching never turns the pathways into one joint regression model. README, current workflow documentation, examples, and the user wiki describe this use case.

## Evidence baseline

The active repository is the local `ldsc_py3_restructured` checkout on `restructure`. The pre-refactor evidence revision was `a505c45d61d35b39a3fc086c2b623dd7d6936d96` (`refactor: remove annotation fingerprints`). Unrelated working-tree changes are outside this task.

That commit removes annotation fingerprint generation, retained common-reference annotation values, result fields, and written fingerprint metadata. It preserves annotation-name validation and advisory value classification. Quantile reconstruction now checks alignment, recorded universe sizes, common-SNP annotation sums, and available overlap cross-products; these aggregate checks do not prove the original annotation value at every SNP. Independent artifact identity hashes remain. See [annotation semantics](../../src/ldsc/annotation_semantics.py), [LD-score calculation](../../src/ldsc/ldscore_calculator.py), and [quantile reconstruction and provenance checks](continuous-annotation-quantile-h2.md#reconstruction-and-provenance-checks).

At the pre-refactor evidence revision, the remaining eager paths included `AnnotationBuilder._detect_chromosome_shards()` and `parse_annotation_file()` reading complete source tables, `_run_sharded_inputs()` concatenating complete chromosome bundles, and `LDScoreCalculator._run_chromosomes()` eagerly slicing/submitting chromosome inputs. See [annotation builder](../../src/ldsc/annotation_builder.py) and [calculator](../../src/ldsc/ldscore_calculator.py). Those eager paths are now replaced as described in [the implemented design](annotation-memory-design.md); fingerprint aggregation remains removed.

## Confirmed scope and compatibility

- Optimize standalone `annotate` and direct `ldscore`, including input preflight, chromosome discovery, annotation loading, alignment, projection, computation, and writing.
- Cover chromosome-sharded and whole-genome annotation files, BED queries, gene-list queries, and synthetic `base`, with both PLINK and parquet-R2 reference backends where those routes apply. Standalone `annotate` will accept exactly one BED or focal gene-list query route, require baseline annotations and an output directory, and exclude control-gene lists. This functional expansion was approved separately; its validation boundary, condition–outcome table, diagnostics, and integration ownership are recorded in [standalone annotate gene-list decisions](annotate-gene-list-decisions.md).
- Backward compatibility with the existing Python API is not required. Replace obsolete interfaces directly and update affected internal callers, tests, and examples. This does not authorize unrelated interface changes.
- Preserve numerical results, SNP identities and ordering, filtering, annotation validation, counts, overlap semantics, and the reference-versus-regression SNP universes. Existing inconsistencies that would require a scientific behavior change must be surfaced separately.
- Annotation fingerprints remain removed under the policy in `a505c45`; the refactor does not restore their retained matrices or introduce replacement fingerprint staging.
- Implementation remains local. Do not launch a full 1,000-query HPC run for validation.
- The user subsequently expanded the optimization to batch `partitioned-h2` preparation and query loading, exact streaming `quantile-h2` statistics, gene-index loading and indexed LD-score assembly, incremental regression detail writing, and bounded shared metadata/diagnostics. These areas are now in scope; the earlier downstream-regression exclusion and proposed quantile materialization exception are superseded.
- Direct `ldscore` also gains query batching inside one shared numerical LD traversal and score accumulation restricted to resolved output rows. Both PLINK and parquet-R2 are covered. This addition preserves every previously agreed workflow area.

## Confirmed dataset contract

`AnnotationBundle` remains one object representing the complete annotation dataset. Its data are organized as chromosome shards loaded or constructed on demand. Shared metadata consists of chromosome scope, ordered baseline/query column names, configuration, provenance, and the information needed to access each chromosome. SNP row metadata belongs with the corresponding chromosome's annotation values.

The bundle must not hold all chromosome DataFrames in a dictionary or list. Normal execution must not automatically cache every shard, concatenate shards into whole-genome annotation matrices, or rematerialize those matrices when returning results.

For sequential execution, the lifecycle is `load/build chromosome -> consume or write its shard -> release chromosome data -> process the next chromosome`. For parallel execution, active annotation data and queued work are bounded by the configured worker count; inputs must not be eagerly loaded or serialized for every chromosome.

Release means that chromosome-owned annotations, reader state, temporary compute buffers, mappings, and retaining views are no longer live before the next chromosome loads. Python/native allocators may keep freed pages available for reuse, so operating-system RSS need not fall immediately. Final HM3 outputs and explicitly permitted shared state can remain live; the memory guarantee does not require those results to disappear between chromosomes.

Standalone `annotate` writes chromosome shards incrementally and returns an object referencing those outputs. Final aggregation retains necessary shared metadata and diagnostics rather than complete annotation matrices.

`AnnotationBundle` is a resource-owning handle supporting context-manager use and explicit `close()`. It owns private staging. A calculator receiving an existing bundle borrows it without closing it. Closing removes owned private storage and prevents further shard loading, while leaving original inputs and persistent outputs untouched. Returning a bundle does not close resources that the returned handle still needs.

The returned complete annotation dataset may depend on original baseline sources. Generated query shards remain persistent; no portable copy of all baseline annotations is required. Private baseline staging remains owned by the live bundle.

New annotation/analysis staging stays inside the explicitly supplied output directory, rather than using the system temporary directory or `TMPDIR`. Python workflows that may require staging must require an explicit output directory; pure numerical functions operating on supplied arrays may remain file-free. Owned scratch is removed on closure or handled failure. Interrupted runs start fresh, without automatic reuse or deletion of abandoned directories.

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

For each chromosome, compute/read one correlation block or decoded R2 chunk, project baseline/control columns and the regression-weight mask once, then project all query batches against that same block before releasing it. Do not wrap query batching around a complete chromosome kernel, reset genotype traversal per query batch, or reread the entire R2 file per query batch. Existing bounded validation reads remain allowed; one numerical traversal does not mean one total input read.

Batch annotation counts, advisory classification, baseline-query overlap, and query self-overlap too. Keep shared baseline statistics and the existing baseline-rows-plus-query-diagonal overlap representation. Do not introduce cross-query products or full-width float64 annotation/statistics temporaries.

Resolve the existing output/regression SNP selection and region exclusions before allocating scores. Accumulate only those output rows. Keep the supported all-output behavior when no output restriction is requested; do not hard-code HM3 into the numerical kernel. An output-row-by-all-query accumulator is permitted, but a full-reference-row-by-all-query score accumulator followed by slicing is not.

For each output SNP, annotation contributions still come from the complete retained reference universe. Non-output endpoints must contribute to output endpoints in symmetric pair accumulation. Counts and overlap retain all/common reference universes. `w_ld` keeps its separate filtered-regression contributor universe, diagonal, and exclusion policy. See the [SNP-universe and traversal contract](ldscore-snp-universe-contract.md).

Preserve current accumulation precision and applicable numerical tolerances. No reduced-precision policy is approved. Output-row accumulators and final materialized LD-score tables may grow with total query count; batching bounds active query workspace, not total process memory.

## Confirmed result and resource boundary

Final HM3 LD-score tables may remain materialized in memory. Changing `AnnotationBundle` does not require replacing `LDScoreResult` with a file-backed result interface. Eliminate whole-genome annotation retention and avoid unnecessary copies of completed LD scores.

Hard output rule: do not chromosome-shard the LD-score tables written by `ldscore`. Keep one aggregated `ldscore.baseline.parquet` and, when applicable, one aggregated `ldscore.query.parquet`, with existing metadata and overlap artifacts. Existing chromosome row groups inside those single files are not separate chromosome files. Private annotation shards and in-memory chromosome computation do not change this public output contract.

Require exact identities, ordering, inclusion decisions, annotation values, and quantile membership except for the explicitly approved duplicate-policy correction. Streaming reductions may differ in final floating-point digits within the applicable existing tolerances. Keep scientific validation thresholds unchanged, including quantile aggregate validation at `rtol=1e-6`, `atol=1e-8`; report observed differences instead of weakening checks to accommodate a regression.

The user has no hard RAM ceiling, runtime constraint, or maximum acceptable slowdown. Measure peak memory, runtime, and temporary disk usage separately across the expanded workflows; distinguish measurements from extrapolations and explain remaining bottlenecks. At the previously cited dimensions, one float32 HM3 matrix with 1,053 annotations contains approximately 4.45 GiB of values; this is an array-size estimate, not a measured process peak. The earlier approximately 23.4 GiB common-annotation matrix estimate no longer establishes a fingerprint storage requirement after `a505c45`.

## Confirmed downstream workflow scope

### Batch partitioned heritability

Prepare SNP alignment, trait data, and shared baseline LD scores once. Load one query or a bounded query batch on demand instead of retaining every query column. Each model fits its complete retained genome-wide SNP set. Preserve model-specific filtering, regression weights, and jackknife calculations; do not share quantities that depend on the fitted model. At the evidence revision, repeated assembly and detail accumulation occurred in [regression_runner.py](../../src/ldsc/regression_runner.py), `RegressionRunner.build_dataset()` and `estimate_partitioned_h2_batch()`; that directory loader read complete query tables in `load_ldscore_from_dir()`.

### Exact global quantile statistics

Determine the validated common SNP universe and the existing exact quantile boundaries from eligible target values. Process annotation shards using those same global boundaries, retaining running quantile sums, overlap-validation products, and numerically stable standard-deviation statistics instead of every SNP's fitted annotation values. Preserve ties, missing-target policy, aggregate validation, and full-common-universe statistics for standardized coefficients. Do not build whole-genome fitted-annotation matrices or dense SNP-by-quantile indicator matrices. The current numerical definitions and eager assembly are in [quantile_h2.py](../../src/ldsc/quantile_h2.py), `assign_legacy_quantiles()`, `_prepare_quantile_inputs()`, and `run_quantile_h2_from_args()`.

The exact boundary-selection stage may retain one genome-wide eligible numeric target vector and its sorting workspace. Release that vector when boundaries are established and it is no longer needed; fitted annotations remain chromosome-scoped. At 5,961,159 values, one float64 vector occupies approximately 45.5 MiB, excluding sorting workspace and other state.

### Gene-index loading and indexed assembly

Represent a loaded index using shared metadata and chromosome component paths. Validate with bounded memory; load one chromosome's sparse operator, process all query batches for that chromosome, and release it before advancing. Keep the operator loaded throughout those batches instead of rereading it for each batch. Preserve index integrity checks and construction staging. Smaller multiplications and separate validation passes can cost throughput or I/O; measure rather than promise unchanged runtime. The relevant current seams are [gene_ldscore_index.py](../../src/ldsc/gene_ldscore_index.py), `LoadedGeneLDScoreIndex`, `_load_gene_ldscore_index()`, and `run_indexed_ldscore()`.

### Incremental regression detail writing

Write per-query category tables, coefficient delete values, and associated metadata after each fit completes. Retain aggregate summary rows and compact bookkeeping instead of all detailed results. Preserve output contents, ordering, diagnostics, and existing failure/overwrite contracts. Returned objects reference valid persistent outputs. Current [output writing](../../src/ldsc/outputs.py), `PartitionedH2DirectoryWriter._query_records()`, assigns folder ordinals from the final sorted summary; preserving this order requires final path assignment after summary ordering is known, even if detail contents are written earlier.

The user explicitly approved keeping these incremental detail writes private until the complete batch succeeds. Assign canonical folders from the final stable summary order and publish through the existing writer lifecycle. This does not add early public access, resumable partial fits, a new skip-on-error policy, or stronger rollback guarantees.

### Shared metadata and diagnostics

Give chromosome metadata caches explicit ownership and lifetimes, releasing data when the consuming stage finishes. Stream potentially large row-level diagnostics to their output files while retaining summary counts and representative examples. Preserve complete diagnostic records and stage-specific validation behavior. This includes annotation identity drops, SNP alignment diagnostics, and gene-list row audits; returning results must not reload all diagnostic records. Downstream commands load only needed LD-score columns even though LD-score generation may return materialized final HM3 tables.

## Decision map

Every identified material branch is mapped below. The user confirmed the complete map and the consolidated implementation recommendations on 2026-09-10. Measured storage tuning and validation evidence remain implementation work within these contracts.

| Branch | Decision and boundary | State |
| --- | --- | --- |
| Objective and workflow scope | Many independent pathway fits; annotate, direct/indexed LD scores, batch regression, exact quantiles, diagnostics and writers all remain in scope | Confirmed |
| Standalone gene-list behavior | Exact input route, catalog/build/padding/exclusion rules, Gate A/B outcomes, annotation-grid support, output reloadability | Confirmed; separate condition–outcome table |
| Dataset interface | One lazy complete-dataset handle; chromosome/column/row access; no whole-genome annotation materialization | Confirmed |
| Ownership and lifetime | Explicit owner/borrower/close; permitted original baseline dependencies; sequential release and bounded worker submissions | Confirmed |
| Scratch and interrupted runs | New scratch inside output directory; owned cleanup; fresh retry without abandoned-stage reuse/deletion; existing index-construction exception | Confirmed |
| Input preparation and identity | Bounded preflight/discovery and source preparation; dataset-wide annotation duplicate cleanup; stage-specific validation preserved | Confirmed |
| Direct computation | Per-block query batching, default 1000; output-row accumulators retain full reference contributors and distinct regression-weight policy | Confirmed |
| Public LD output | One aggregate baseline/query file family; final HM3 tables may remain materialized | Confirmed |
| Indexed assembly | Validate boundedly; one chromosome operator retained through its query batches; construction integrity/precision preserved | Confirmed |
| Batch regression | Shared preparation; each pathway fits its own full retained genome-wide model; per-model filtering/weights/jackknife preserved | Confirmed |
| Quantile statistics | Exact global boundaries; one numeric target vector/sort workspace permitted; sufficient statistics from shards | Confirmed |
| Regression publication | Private per-fit detail writing; final sorted canonical publication; existing failure/overwrite behavior | Confirmed |
| Diagnostics and shared state | Complete row records streamed; summaries/examples and explicit metadata lifetimes | Confirmed |
| Scientific equivalence | Exact identities/order/inclusion/values/membership; existing floating tolerances; no fingerprints restored | Confirmed |
| Implementation structure | Column-major private storage, disk-backed global identity bookkeeping, detached reads, explicit prepared/result split and accumulator owners | Confirmed |
| Evidence and delivery | Local numerical/artifact/lifetime checks and separate memory/runtime/disk measurements; post-implementation pathway-focused docs/wiki | Confirmed; validation details in plan |

## Closure audit and implementation checkpoints

The design-stage source audit identified preservation details, rather than new behavioral decisions. The following references describe the evidence revision; current replacements are documented in [annotation-memory-design.md](annotation-memory-design.md). At the evidence revision, gene-list resolution retained a combined audit and all per-source gene selections in `GeneListBatchResolution`; its replacement must stage row-level resolution and selection data while sharing catalog lookup state and compact source summaries. Controls remain the baseline column `gene_control` in direct/indexed workflows. LD-directory schema/identity checks remain separate from model-specific numeric checks and query-row alignment. Quantile external targets preserve float64/raw-token parsing and fatal duplicate-identity checks. See [gene resolution](../../src/ldsc/gene_list_resolver.py), `resolve_gene_lists()` and `GeneListBatchResolution`; [annotation construction](../../src/ldsc/annotation_builder.py), `_run_single_universe()`; [regression loading and assembly](../../src/ldsc/regression_runner.py), `load_ldscore_from_dir()` and `_assemble_regression_ldscore_table()`; and [quantile input preparation](../../src/ldsc/quantile_h2.py), `_read_target_annotation()` and `_prepare_quantile_inputs()`.

No unresolved product, scientific, or architectural choice remains. The specification's storage layout, bounded reads, global bookkeeping, accumulator ownership, and prepared/result interfaces are confirmed. Tile sizes, read coalescing, and measured throughput are implementation experiments within those contracts, not undefined scientific policy. Existing stage-specific gates, immutable-input assumptions, complete diagnostics, canonical ordering, and numerical tolerances remain fixed. The user requested final specification documentation after confirming the design; implementation has not started.
