# Annotation and Workflow Memory Optimization Specification

Last updated on: 2026-09-10

Status: confirmed specification, approved on 2026-09-10. Implementation has not started. This document captures the complete agreed design, including the internal storage and resource interfaces. The [decision map](../current/annotation-memory-decisions.md#decision-map) records the decisions; the [implementation plan](../plans/2026-09-10-annotation-workflow-memory.md) sequences delivery and verification.

## Problem, goal, and scope

Current workflows can retain whole-genome annotations, all query columns, chromosome operators, or detailed results for every fitted model. These allocations become costly when one run tests many pathways. Annotation fingerprint retention was removed in commit `a505c45`; this optimization addresses the remaining preparation, computation, and result-lifetime allocations described in the [evidence baseline](../current/annotation-memory-decisions.md#evidence-baseline).

Enable large pathway batches, such as testing enrichment for 1,000 pathways in one run, while each pathway is fitted separately against the baseline categories. Share annotation preparation, reference-LD traversal, and model-independent regression preparation. Bound live annotation data, query workspace, metadata caches, and detailed regression results instead of retaining every chromosome and query throughout the workflow.

The scope includes standalone BED and gene-list `annotate`; direct PLINK and parquet-R2 `ldscore`; lazy gene-index loading and indexed LD-score assembly; batch `partitioned-h2`; exact global `quantile-h2`; and their discovery, preflight, diagnostics, and writers. Preserve each command's existing input capabilities except for the explicitly approved addition of focal gene lists to standalone `annotate`. The complete gene-list contract and condition–outcome table remain authoritative in [annotate gene-list decisions](../current/annotate-gene-list-decisions.md).

Keep the canonical LD-score directory: one aggregated `ldscore.baseline.parquet`, one optional aggregated `ldscore.query.parquet`, and the existing overlap and metadata artifacts. Do not chromosome-shard these LD tables. Existing chromosome row groups inside one Parquet file remain valid. Standalone annotation outputs continue to use `query.<chrom>.annot.gz`.

No Python compatibility layer is required. Replace obsolete interfaces and migrate internal callers, exports, tests, and examples. Annotation fingerprints remain removed under inspected HEAD `a505c45d61d35b39a3fc086c2b623dd7d6936d96`; independent index/reference integrity checks remain unchanged.

## Resource contract and ownership

| State | Owner and permitted lifetime |
| --- | --- |
| Dataset scope, ordered column names, configuration, provenance, source/shard descriptors | One `AnnotationBundle` representing the complete dataset |
| Shared gene catalog lookup and compact per-gene support | Resolution/projection or index context; no copy per pathway/chromosome; release after its consumers finish |
| Normalized private annotation storage, global identity bookkeeping, staged audits | Bundle or workflow owner, under its explicitly supplied output directory |
| One chromosome's row metadata, source-row mapping, selected annotation blocks, reference reader, LD blocks | Chromosome context; released when that chromosome finishes |
| Direct LD-score accumulator | Chromosome worker; float64, output rows by all requested output columns; released after existing direct-output conversion |
| Completed HM3 LD-score tables | Final result; may remain materialized in memory |
| Shared regression SNP alignment, traits, baseline LD scores | Prepared batch-regression context; retained for the batch, without duplicate chromosome caches |
| Query LD scores and fitted model state | Current query or bounded query batch; detailed fit state released after writing |
| Exact quantile target values and sort workspace | Boundary-selection stage only; one eligible genome-wide numeric vector is permitted |
| One index chromosome's sparse operator | Validation or assembly context; held across all assembly query batches, then released |
| Diagnostic counts, examples, summaries, artifact paths | Workflow/result; complete row-level records remain on disk |

`AnnotationBundle` provides chromosome, column, and row selection; it never automatically caches all chromosome frames, concatenates annotations over the genome, or materializes them on return. Loading a query batch must not first construct a dense chromosome-by-all-query annotation matrix. Selected blocks are detached bounded arrays, so a small view cannot retain a complete matrix or live mapping.

Use explicit context management and `close()`. A bundle owns its private storage; a calculator receiving a caller-owned bundle borrows it. Closing prevents further loads, removes owned private storage, and leaves original inputs and persistent outputs untouched. A returned bundle stays live and may depend on original baseline sources plus its owned staging. Source-backed workflows that may stage require an explicit output directory; pure numerical functions receiving supplied arrays may remain file-free. Original inputs must remain available and unchanged while used, consistent with the repository's immutable-input contract; do not add routine content hashes or promise a portable snapshot.

Each run gets a unique private directory inside the output directory. Route relevant library temporary work, including BED projection scratch, into that owned directory. Clean owned scratch on normal closure and handled failure. After interruption, a new run starts fresh and does not automatically reuse or delete abandoned private directories. Preserve the existing gene-index constructor's sibling staging, build-state, lock, and atomic publication as the explicitly approved exception. Existing output preflight, overwrite, failure-marker, and cleanup ordering remain authoritative; this refactor adds no rollback policy.

With one worker, finish consumption/writing and release the current chromosome's annotations, reader, mappings, retaining views, and compute buffers before loading the next chromosome. Freed allocations become available for reuse; Python/native allocators can retain pages, so process RSS may stay high. Final LD tables and other explicitly permitted shared state remain live. With multiple workers, bound active and pending chromosome tasks by worker count and submit lightweight descriptors instead of eager matrices. Consume completed results promptly without accumulating completed futures holding annotation workspaces.

## Bounded source preparation and access

The bound begins at chromosome discovery and input preflight. Whole-file parsing followed by chromosome selection fails this contract. Scan whole-genome text/gzip sources in bounded chunks, preserve source and row order, and create private chromosome storage when random access will be needed. Do not decompress and parse a whole-genome source once per chromosome or LD block. Existing chromosome-sharded inputs can be used directly where the consumer can read them with bounded memory; repeated numerical column/row access should use prepared seekable storage.

Preserve alignment validation at the same logical input boundary as today. Validate corresponding rows across aligned column sources before partitioning can hide a mismatch. Preserve active identity normalization, numeric validation, scope/build checks, and stage-specific error aggregation. Dataset-wide drop-all identity cleanup must span chunks and chromosomes. Rows from aligned baseline/query column files describe the same logical SNP row and must not count as duplicate observations. Cross-chromosome rsID collisions are removed consistently in sharded and whole-genome layouts, as explicitly approved. Do not transfer this cleanup behavior to immutable artifact integrity checks.

Use private chromosome-local float32 column-major value storage, with row metadata and ordering recorded separately. The current normalization boundary already uses float32 in `_kernel/annotation.py::_validate_annotation_values`; this is not a new precision reduction. Use explicit bounded reads into detached arrays for requested columns and row runs. Avoid permanent full-width mappings or row-major reads that fault unrelated query columns into memory. Coalesce nearby requested rows only within a bounded read budget; a sparse selection must not trigger reading its entire minimum-to-maximum span.

External quantile target sources use their existing distinct parsing contract: preserve raw-token exclusion semantics and float64 numeric values from `quantile_h2.py::_read_target_annotation()`. Do not narrow those targets through ordinary float32 annotation staging; doing so can collapse distinct values and change quantile ties. Their private staged storage must preserve that precision and missingness information.

Use bounded append spools while final row counts and global cleanup are unknown, then produce seekable chromosome storage once those dimensions are known. Keep genome-wide identity bookkeeping in a disk-backed table keyed by the existing exact identity representation, rather than an unbounded in-memory key set or new fingerprint. A standard-library SQLite table is the initial implementation choice; its scratch and bounded cache belong to the same private owner. A column-major `.npy` value file with explicit bounded I/O is the initial storage choice. Private format, tile sizes, and read coalescing are implementation details to check with the first local storage pilot; change them if evidence warrants without changing the access/ownership contract or adding a public file format.

Verify database journal and temporary-sort placement in that pilot; locating the database inside the output directory alone does not establish scratch containment. Use bounded operations and contained spill behavior. Likewise, enforce bounded Parquet decode buffers when an input row group exceeds the configured pair-chunk size, instead of reading the complete oversized row group before chunking it.

Do not normalize or project the same query repeatedly for every LD block. Build BED/gene values by chromosome and bounded query groups, store them once when repeated access is needed, and replay prepared values. Output writing uses bounded row buffers and preserves canonical full-row column order. Buffer sizing must account for column width; it must never select an entire chromosome merely because it is one output shard.

Bound query-source preparation as well as the resulting SNP values. Parse BED sources with bounded input buffers and keep prepared chromosome/source interval descriptors. Resolve gene sources once using shared catalog lookup state; stage complete row audits and deduplicated per-source selections instead of keeping the current all-source audit frame and tuples of every pathway's gene IDs/intervals live. Load selected-gene/interval data for the active chromosome/query group, retain compact source statuses, and preserve per-source duplicate order and exhaustive Gate A reporting. This is the application of the shared retention contract to `GeneListBatchResolution`, not a change to gene-resolution rules.

Stage complete row-level diagnostics with stable source/row/stage identifiers. Emit them in the existing deterministic order after global outcomes are known, retaining only counts and representative examples. Safely discoverable errors at each existing gate remain exhaustive. Global preparation may require bounded validation/replay passes; one numerical LD traversal is a separate requirement.

## Standalone annotation and gene lists

Replace BED-only Python wrappers with `run_annotate()` and `parse_annotate_args()`; retain `run_annotate_from_args()` for CLI dispatch. Require baseline sources, output directory, and exactly one BED or focal gene-list query route. No control-gene route is added.

Resolve the catalog and gene identifiers once. Preserve strict/resolved-only omission rules, explicit nonnegative gene padding, BED omitted padding of zero, `none`/`mhc` gene exclusions before padding, build agreement, and separate projection-build provenance for rsID identity. Gate A collects safely discoverable catalog, source, identifier, and naming failures. Coverage comes from actual baseline contents; `@` requires autosomes 1–22, and selected genes outside scope fail after explicit exclusions under either resolution policy.

Gate B measures support on baseline rows surviving global identity cleanup. Use bounded construction/diagnostic passes to establish global query status before deciding final columns. Skip empty and globally unsupported focal queries while continuing usable siblings; fail if all focal queries are skipped. Keep all-one annotations and zero-valued chromosome shards for queries supported elsewhere. Reference-panel and regression support remain unevaluated.

Write persistent annotation shards incrementally and return a bundle referencing them and the allowed baseline dependencies. Preserve SNP identity/order, values, duplicate policy, explicit `CM=NA`, absence of annotation MAF columns, reloadability, gene audits, summaries/statuses, provenance, diagnostics, and existing overwrite/failure behavior. Label support `annotation_snp_count`. The linked gene-list decision table governs detailed edge cases rather than a second divergent rule set here.

## Direct LD-score interfaces and traversal

Set `query_batch_size` to a positive integer, default **1000**, in direct `ldscore`, indexed `ldscore`, and batch `partitioned-h2`. Reject zero, negative, and nonintegral values at argument/config validation. A value at least the query count yields one query batch. This is separate from chromosome workers, SNP block/chunk size, gene-index construction atom batches, and the number of regression models fitted jointly. No new query worker pool is introduced.

Replace `PreparedChromosome.annotation_matrix` with a bounded accessor plus explicit chromosome state: full retained reference metadata, annotation-source row mapping, all/common masks, resolved output-row indices, and reference reader/window information. Output membership is the existing regression-SNP selection with declared exclusions, resolved before score allocation. Preserve supported unrestricted output selection; the kernel must not hard-code HM3 membership.

Allocate one chromosome-owned float64 score accumulator with shape `n_output_rows × (n_output_annotation_columns + 1)`; the final column denotes the separately governed regression-weight projection. Count each baseline, control, and focal query column exactly once. The optional control remains baseline column `gene_control` in direct/indexed output; it is projected with shared baseline work and is never moved into the focal query batch or output table. Preserve existing accumulation precision and direct float32 result conversion. Final output arrays may grow with total query count. Query batching bounds active projection/statistics workspace, not this explicitly permitted output memory. Avoid a full-reference-row-by-all-query score buffer followed by slicing.

For each PLINK correlation block or decoded parquet-R2 pair chunk, project baseline/control columns and the regression-weight mask once, then all query batches against the same block. Release the block only afterward. Preserve genotype normalization, filtering, windows, bias correction, pair directions, and diagonal handling. Do not reset the genotype cursor or rerun a chromosome kernel per query batch. An all-zero optimization must consider the correct block/columns; testing only the first query batch can incorrectly discard later contributions.

PLINK uses the existing within-block and cross-block correlations with output-row mappings applied to both contribution directions. Parquet-R2 constructs reusable output-destination/contributor sparse projections for both directions of each decoded chunk. Gather the corresponding reference contributor annotation blocks, cast only the active blocks to the required float64 multiplication precision, and accumulate only output destinations. This replaces the full-width float64 annotation copy without lowering multiplication precision.

For an output SNP \(j\) and annotation \(c\), compute \(\ell_{j,c}=\sum_{k\in\mathcal R}r_{jk}^{2}a_c(k)\) over the complete retained reference contributor universe \(\mathcal R\). An endpoint outside the output set must still contribute to an endpoint inside it. Output restriction does not shrink reference counts, annotation alignment, LD neighborhoods, or all/common overlap universes. The separate `w_ld` projection retains its filtered regression contributor set, diagonal, and exclusion policy. The [SNP-universe contract](../current/ldscore-snp-universe-contract.md) remains authoritative.

Return output-row metadata and scores separately from full/all/common reference counts, pre-exclusion selected-row counts, excluded-row counts, and compact sufficient statistics. Capture these existing provenance/count fields before releasing their full-reference metadata source. Audit every consumer that currently assumes score rows equal reference rows. Write required full-reference metadata sidecars while the prepared chromosome metadata is alive, then release it. Aggregate completed baseline/query outputs directly and retain compact chromosome summaries, avoiding join/concatenate/split copies plus a second collection of chromosome tables.

Batch annotation counts, classification, baseline-query overlap, and query self-overlap. Keep baseline statistics shared and preserve the baseline-rows-plus-query-diagonal overlap representation; cross-query products are unnecessary. Preserve existing operation-specific precision: current count reductions and float64 overlap reductions do not establish one universal dtype policy. Existing kernel callers for gene-index construction must retain their float64 stored operator/common payloads and intentional atom traversals; direct-result float32 conversion must not leak into those callers.

## Downstream workflows

### Gene-index loading and indexed assembly

Represent the index with shared metadata and chromosome component paths. Run existing integrity validation one chromosome at a time and release each operator after validation; accumulate compact gene-support summaries during that pass. Assembly may reload a validated chromosome once, keep its operator loaded across all query batches, and release it before advancing. Preserve query ordering, immutable-index checks, stored precision, exact atom semantics, construction staging, and aggregated LD-score outputs. Measure separate validation and assembly I/O and matrix throughput; batching can trade throughput for memory.

### Batch partitioned heritability

Use a selective LD-score source with shared schema/metadata, shared baseline values and identity alignment, and explicit query-column reads. Prepare trait data and common alignment once. Retain a query or bounded query batch, fit each model independently against the baseline over its complete retained genome-wide SNP set, and release query values when consumed.

Preserve loading-time checks of artifact provenance, required column schemas, allele metadata, and count/overlap configuration without loading every query value column. Preserve query-row alignment at the applicable fit boundary; loading currently uses `require_query_alignment=False`, while a query fit enforces aligned baseline/query identities and order. Validate numeric values and model-dependent counts/filters for the same selected columns as today. Do not introduce validation of unselected query values or skip existing validation merely to reduce I/O. No-query functional partitioned-h2 remains one joint baseline fit with its existing outputs and estimator behavior.

Share only model-independent work. Keep model-specific filtering, annotation counts/count-key behavior, design columns, variance/collinearity decisions, chi-square caps, weights, two-step behavior, jackknife calculations, and model outcomes inside the individual fit. Do not fit separate chromosome regressions or combine all pathways into one regression. Existing `RegressionRunner.build_dataset`, `_fit_h2_dataset`, `estimate_partitioned_h2_batch`, and selective loading in `regression_runner.py` are the migration seams.

Write each completed fit's category tables, coefficient delete values, and metadata privately, then release those detailed objects. Return persistent paths/manifests, summary rows, and compact bookkeeping instead of dictionaries retaining every detail matrix. Publish only after the complete batch succeeds, with canonical folders assigned after the final stable summary sort. `PartitionedH2DirectoryWriter._query_records()` currently derives ordinal folders from that sort. Preserve existing fit-abort and writer failure/overwrite behavior; private staging does not add skip-on-error fitting, public partial results, resumability, or whole-family rollback. Flush required failure diagnostics before cleaning owned scratch, preserving their existing persistent paths and schemas.

### Exact global quantile statistics

Use two bounded annotation passes. First validate the common SNP universe, annotation aggregates and available overlap products, and collect eligible numeric target values. Accumulate stable full-common-universe mean/variance statistics for standardized coefficients, including rows excluded from quantile assignment by target missingness. Determine the same exact global boundaries as `assign_legacy_quantiles()`; keep its boundary positions, lower-bin tie handling, and empty-bin failures. One eligible numeric target vector plus sorting workspace is permitted, then released when no longer needed.

In the second pass, load one chromosome's selected fitted annotations and target values, assign membership using the established global boundaries, and add quantile annotation sums to small sufficient-statistic arrays. Reuse `compute_quantile_h2()`'s statistical definitions. Do not construct whole-genome fitted annotations or a dense SNP-by-quantile indicator matrix. Preserve missing-target rules, same-name target checks, counts/overlap validation, coefficient ordering, and standardization over the complete validated common universe.

Preserve quantile-specific identity failures across the complete dataset. Duplicate effective identities in reference metadata or external targets remain fatal, and omitted-allele inference requires globally unique reference base identities. Do not replace these checks with ordinary annotation drop-all cleanup or chromosome-local uniqueness checks. See `quantile_h2.py::_prepare_quantile_inputs()` and `_effective_keys_with_reference_inference()`.

## Observable acceptance criteria

| Scenario | Required outcome |
| --- | --- |
| Equivalent whole-genome and chromosome-sharded annotation inputs | Same cleaned identities, ordering, values, scope, and diagnostics, including global duplicate handling |
| One worker advances to the next chromosome | Previous chromosome annotation/reference workspace has no retained owner or view; permitted output/shared state remains |
| Batch sizes 1, a non-divisor, and at least the query count | Equivalent outputs under existing tolerances and one shared numerical LD traversal per chromosome |
| Non-output reference SNP contributes to an output SNP | Contribution retained; output-only allocation does not change contributor or count universes |
| No output restriction | Existing all-output behavior retained |
| Mixed usable and unsupported standalone focal gene lists | Usable queries written with complete diagnostics; all skipped remains fatal |
| One index chromosome serves several query batches | Operator loaded once for that assembly chromosome, with existing integrity validation preserved |
| Many pathway regression models | Same per-model results as individual runs; shared preparation reused; details do not accumulate in RAM |
| Quantile ties span chromosomes or some target values are missing | Same exact global membership and full-common standardization as current numerical definitions |
| Successful workflow returns | Every returned artifact path is valid; no eager annotation or diagnostic rematerialization |
| Bundle closes after successful standalone annotation | Further bundle reads fail; original baselines and persistent query shards remain intact |
| A later regression fit fails after earlier fits completed | Earlier details remain private; required diagnostics and existing failure markers survive owned-scratch cleanup |
| Handled failure or interrupted prior run | Existing markers/overwrite behavior preserved; only current owner scratch cleaned automatically |

## Validation and measurement

Before structural refactors of files over 300 lines, complete the separately verified cleanup required by `AGENTS.md`. Implement through meaningful failing tests and stable scientific seams. Compare small independent dense LD calculations with `(R @ A)[output_rows]`, including off-output contributors, both directions, diagonal terms, binary and signed quantitative annotations, and the separate `w_ld` universe. Preserve index-construction precision with direct/indexed oracle checks. Add traversal counters, accessor read/allocation guards, and lifecycle tests so numerical success alone cannot hide full-width materialization or repeated LD work.

Require exact identities, row/column ordering, inclusion decisions, normalized annotation values, and quantile membership except for the approved global duplicate correction. Changed reduction order may produce ordinary final-digit differences within existing applicable tolerances; do not weaken those tolerances. Quantile aggregate validation stays at `rtol=1e-6`, `atol=1e-8`. Regression comparisons must include category tables, coefficient delete values, diagnostics, model-dependent exclusions, and canonical folder ordering, not only headline estimates.

Run matched local measurements against post-fingerprint-removal HEAD `a505c45`, separately reporting process-tree peak memory, runtime, peak temporary disk, and persistent output size. Vary chromosome count, query count, batch size, and worker count. Include a small 1,000-query workload where practical; no full production/HPC run is authorized. Use allocation/lifetime evidence to distinguish accepted final-output growth from annotation retention. No hard memory ceiling or slowdown limit is specified; do not claim unchanged runtime without measurements.

The implementation plan lists focused suites, real CLI/Python smoke workflows, artifact reload/failure checks, and final full pytest then unittest checks. Run these suites sequentially because shared BED temporary cleanup can interfere. Specification verification covers documentation and source contracts; numerical correctness and resource improvement require evidence from the implementation.

After implementation, update README, help, public docstrings, current architecture and workflow docs, examples/tutorials, and the user wiki. Explicitly describe the 1,000-pathway use case, independent baseline-adjusted models, default batch size, sequential release versus allocator RSS, output-contained scratch and its construction exception, and aggregated HM3 files.

## Risks and open questions

There are no unresolved product, scientific, or architectural decisions. The consolidated design and its internal implementation choices are confirmed. Tile sizes, read coalescing, and physical storage tuning remain measured implementation choices within the fixed contracts.

Staging I/O and smaller matrix calls may reduce throughput. Hidden full-width copies can defeat the memory bound. Source-row mistakes or restricting contributors together with output rows can change scientific results. The storage pilot, independent numerical oracles, allocation/lifetime checks, and matched resource measurements must address these risks. Revisit the design only if implementation evidence requires changing scientific behavior, public contracts, approved resource boundaries, or scope.

## Out of scope

Chromosome-sharded public LD tables; approximate quantiles; restored annotation fingerprints; reduced accumulation precision; new estimators, model universes, thresholds, or scientific filtering; persistent staging caches/resume; portable copies of every baseline source; redesigned index publication; extra standalone query/control modes; HPC benchmarks; and unrelated repository changes.
