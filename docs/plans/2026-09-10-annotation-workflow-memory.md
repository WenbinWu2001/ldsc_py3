# Annotation and Workflow Memory Implementation Plan

Last updated on: 2026-09-10

Status: implementation in progress. The complete design was confirmed on 2026-09-10 and implementation was authorized afterward. Design documentation was committed as `56d426a`.

Reference specification: [annotation and workflow memory optimization](../specs/2026-09-10-annotation-workflow-memory-design.md). Confirmed scope and decisions are recorded in [annotation memory decisions](../current/annotation-memory-decisions.md); standalone gene-list behavior is governed by [annotate gene-list decisions](../current/annotate-gene-list-decisions.md), especially its condition–outcome table. Keep this plan current as implementation evidence becomes available.

## Goal and success signal

Support large pathway batches, such as 1,000 pathways in one run with each pathway fitted separately against the shared baseline categories, without retaining whole-genome annotation matrices or every model's detailed regression results. Apply the improvement through annotation preparation, standalone `annotate`, direct and indexed `ldscore`, batch `partitioned-h2`, `quantile-h2`, writing, and diagnostics.

Success requires equivalent scientific selections, values, statistics, metadata, and canonical artifacts, together with measured reductions in avoidable active memory. Sequential execution releases each chromosome's annotation, reference, and temporary computation state before loading the next; parallel execution bounds active payloads and queued work by the configured worker count. Released allocations become available for reuse, but process RSS need not immediately fall because allocators may retain freed pages. Final materialized output tables and compact shared state remain live intentionally.

## Constraints and invariants

The inspected starting revision is `a505c45d61d35b39a3fc086c2b623dd7d6936d96`, after annotation fingerprint removal. Current eager seams are `AnnotationBuilder` source loading and bundle concatenation; `PreparedChromosome.annotation_matrix` and all-reference score buffers; `LoadedGeneLDScoreIndex.index_chromosomes`; repeated `RegressionRunner.build_dataset()` calls and retained per-query details; and `quantile_h2._prepare_quantile_inputs()` plus dense indicators. The slices below identify the source modules and tests for each seam.

- Work locally in `ldsc_py3_restructured` on `restructure`, preserving unrelated changes. No HPC execution is part of this plan.
- Keep public LD-score artifacts aggregate: `ldscore.baseline.parquet`, optional `ldscore.query.parquet`, and the existing overlap/metadata family. Do not publish chromosome-sharded HM3 LD tables. Standalone `query.<chrom>.annot.gz` annotation shards retain their existing format.
- Final HM3 LD-score tables and an output-row-by-all-query accumulator may remain materialized. Query batching bounds active query projection/statistics workspace; it does not make total memory independent of total query count.
- `--query-batch-size` defaults to **1000** for direct/indexed `ldscore` and batch `partitioned-h2`. Exactly 1,000 queries can therefore occupy one active batch at the default; a smaller setting reduces query workspace. The setting is independent of chromosome workers and does not combine separate query regression models. Validate it as a positive integer at argument/config boundaries.
- Keep float64 LD-score accumulation, the existing float32 output conversion, and applicable current numerical tolerances. SNP identities/order, annotation values, inclusion decisions, and quantile membership remain exact. Changed reduction order may change the final floating-point digits; do not loosen tests or scientific validation thresholds. Quantile aggregate validation remains `rtol=1e-6`, `atol=1e-8`.
- All new annotation/analysis scratch belongs to a uniquely named private directory inside the supplied output directory. Source-backed Python workflows that may stage require an explicit output directory. Preserve the existing gene-index construction sibling stage, build state, lock, and publication mechanism as the approved exception.
- A bundle owns private staging and supports context-manager use and `close()`. Borrowers do not close a caller-owned bundle. Persistent query outputs and original baseline files survive closure; returned bundles may depend on the original baselines. Handled failure/closure removes owned scratch; interrupted runs start fresh without automatically reusing or deleting abandoned directories.
- Replace obsolete Python interfaces and migrate internal callers, exports, tests, and examples without compatibility aliases. Preserve canonical formats, stage-specific validation, overwrite/failure-marker behavior, independent artifact integrity checks, and regression semantics. Annotation fingerprints remain removed.
- Before structurally refactoring a file over 300 lines, complete the repository-required dead-code/import/export/debug-log check and separate any necessary cleanup. Consult the completed [C1–C3 cleanup disposition](2026-09-10-c1-c3-cleanup.md#complete-disposition); do not repeat completed removals or reopen deliberately retained index locks, publication recovery, or scientific test oracles.

## Confirmed implementation choices

The specification defines the confirmed internal choices. Their performance and correctness need implementation evidence within the agreed scope:

1. Use normalized chromosome-local float32 column-major value storage, separate row metadata, and bounded detached reads by row and selected columns. Start with explicit I/O over `.npy` files and a bounded-cache SQLite table for global identity bookkeeping. Keep append spools and all database/library scratch inside the owned output directory. Measure read coalescing and tile sizes in slice 1; private format changes are permissible if the access/ownership contract remains unchanged.
2. Separate full-reference contributor metadata/mapping from output-row indexing. Allocate chromosome-owned in-memory float64 output-row-by-all-query accumulators; convert completed direct outputs using existing float32 policy and release the accumulation buffers. Keep index-construction float64 payloads. Write full-reference sidecars while their metadata owner remains live; final results retain output rows and compact reference summaries.
3. Keep original inputs available and immutable during source-backed access. Stream/replay diagnostics in existing deterministic source/row/stage order; preserve global gate consequences and omit new routine hashes.
4. Keep per-fit details private until the batch succeeds, then assign canonical folders from final stable summary ordering and publish through the existing writer lifecycle. The user explicitly approved this decision. Preserve existing fitting failures and writer publication behavior without adding public partial results or resume support.

The complete [decision map](../current/annotation-memory-decisions.md#decision-map) records the confirmed behavioral and architectural branches. Routine tile/coalescing measurements remain implementation checkpoints within the agreed contracts.

## Interfaces and migration strategy

| Boundary | Required final interface or invariant |
| --- | --- |
| `AnnotationBundle` | One complete-dataset handle with shared metadata and chromosome/column/row access; explicit ownership and `close()`; no whole-genome frame fields or automatic materialization |
| Prepared chromosome and kernel result | Full-reference contributor metadata/mappings remain chromosome-owned; returned score rows describe only resolved output SNPs; counts and exclusion provenance remain complete |
| Direct accumulator | Float64 `n_output_rows × (n_output_annotation_columns + 1)`; each annotation counted once, plus the separately governed regression-weight projection; existing direct float32 conversion |
| Loaded gene index | Shared catalog/settings/component paths, bounded validation, and one operator held through that chromosome's query batches; construction float64 payloads preserved |
| LD-directory reader and regression preparation | Shared metadata, identity alignment, traits, and baseline values; selected query-column reads; model-dependent filtering/weights/jackknife remain within each fit |
| Regression result and writer | Summary rows plus persistent detail paths/manifests; private per-fit writes, final sorted publication, and no collection retaining every detail matrix |
| Quantile inputs and statistics | Float64 external targets retain missing-token semantics; exact global boundaries and small sufficient statistics; no full fitted-annotation or SNP-by-quantile matrix |

Use expand–migrate–contract for the wide internal changes. Introduce the bounded storage and numerical primitives first, move each workflow and its tests through a complete path, then remove the old eager interfaces when their final callers are migrated. Keep unmigrated commands on their existing working path until their cutover. This is an implementation sequence, not a compatibility promise: do not add compatibility aliases, fallback materialization, or a second permanent storage framework. A migrated path must never call an adapter that rebuilds whole-genome annotations. Delay a shared public rename/removal until all affected callers can switch coherently; update the slice order if that dependency requires it.

Follow the numbered order below by default. Slice 1 establishes common storage/lifecycle primitives. Slices 2–4 migrate annotation and LD workflows; slice 3's numerical core can be developed from slice 1, but its gene-route verification also requires slice 2. Slice 5 can use existing canonical LD directories and supplies selective loading for slice 6. Slice 7 closes all migrations and verifies the complete workflow. No new task or handoff is needed.

## Implementation slices

### 1. Establish bounded annotation preparation and ownership

Status: shared storage and preparation primitives implemented and verified; public workflow migration remains in slices 2–4 and 6. The original whole-genome bundle remains only on unmigrated paths during the agreed expand–migrate–contract sequence.

Evidence: the independent source-to-selected-read tests cover both layouts, aligned column sources, cross-chunk/cross-chromosome duplicates, allele cleanup, detached reads, float64 target preservation, and owner closure. The focused annotation/storage/identity check passed 110 tests and 27 subtests. Before code changes, the full baseline passed 1,409 tests with one skip and 132 subtests. The [storage pilot](../audits/annotation-memory/progress.md) records memory, preparation runtime, and private disk separately; broader command benchmarks remain pending.

Likely areas: [annotation_builder.py](../../src/ldsc/annotation_builder.py), [config.py](../../src/ldsc/config.py), [_ldscore_preflight.py](../../src/ldsc/_ldscore_preflight.py), [snp_identity.py](../../src/ldsc/_kernel/snp_identity.py), annotation parsing and path/identity helpers, and a narrowly scoped shared storage module if the approved interface needs one.

Before replacing the eager paths, capture the small matched baseline workloads needed for the resource comparison in slice 7. Keep those measurements separate from correctness oracles and historical commit-reported test results.

Implement the shard-backed dataset primitives with shared scope/column/configuration/provenance metadata and on-demand access. Prove the complete source-to-private-shard-to-selected-read path before cutting over the public bundle and its consumers. Scan whole-genome sources in bounded chunks during discovery and preflight, normalize once where needed, and use existing chromosome shards directly where practical. Support row and query-column selection without whole-genome reassembly, eager serialization, or implicit caching of every shard. Apply drop-all identity cleanup over the complete logical dataset, including cross-chunk and cross-chromosome rsID collisions; matching logical rows across aligned column sources are not independent duplicate observations.

Introduce explicit owner/borrower closure and output-contained scratch. Stream row-level identity/alignment diagnostics while retaining counts and examples. Migrate callers as their dependent slices land and remove obsolete eager entry points before completion; do not add a second permanent storage framework.

Validation: independently expected duplicate/allele cases across chunks, chromosomes, and whole-genome versus sharded layouts; distinguishable column and row alignment fixtures; preflight scope failures; selected-column access; borrowed/owned cleanup on success and exception; returned-handle lifetime and missing-original-source behavior. Add allocation/load-count checks proving bounded discovery and no automatic whole-dataset materialization. Verify database journal/temporary-sort and library scratch containment rather than relying only on the database path. Use `tests/test_annotation.py`, `tests/test_snp_identity.py`, `tests/test_ldscore_chromosome_scope.py`, and focused tests for any new storage seam.

Exit checkpoint: a deterministic multi-chromosome fixture produces the independently expected selected rows/columns through both input layouts, including cross-chunk duplicates. A small storage pilot demonstrates bounded reads for scattered row indices, no permanent full-width mapping, valid handles, and contained scratch. Resolve storage layout or I/O problems here before migrating every consumer.

### 2. Deliver incremental standalone BED and gene-list annotation

Status: standalone BED/gene CLI and `run_annotate()` migration implemented and verified. Shared staged gene resolution is ready for direct/indexed callers; their migration and retirement of the remaining eager builder belong to slices 3–4. Focused annotation, gene, region, output, CLI/layout, configuration, and failure-marker checks passed 288 tests and 26 subtests. Real module CLI coverage verifies persistent output reload and scratch cleanup. The full repository milestone remains scheduled after the shared LD/index cutover.

Likely areas: [annotation_builder.py](../../src/ldsc/annotation_builder.py), [gene_list_resolver.py](../../src/ldsc/gene_list_resolver.py), [query_annotations.py](../../src/ldsc/query_annotations.py), CLI dispatch/exports, and output/logging helpers.

Expose general `run_annotate()` and `parse_annotate_args()` entry points while retaining `run_annotate_from_args()` as the CLI seam. Accept exactly one BED or focal gene-list route with baseline annotations and output directory; do not add controls. Reuse exact catalog resolution, explicit gene padding/exclusions, build agreement, and separate projection-build provenance for rsID identity.

Resolve genes once, collect Gate A issues, derive scope from baseline contents, and measure Gate B support on globally cleaned baseline rows. Replace all-source audit/selection retention with staged resolution rows and selected-gene descriptors, sharing catalog lookup state and loading only active chromosome/query selections. Apply the shared bounded preparation to BED interval sources too. Determine usable query columns over the complete scope using bounded passes/private staging, then emit canonical annotation shards incrementally and return a valid source-backed bundle. Preserve zero-valued chromosome shards for globally supported queries, all-one annotations, skipped-query policy, complete diagnostics, and `annotation_snp_count` terminology without reference-panel checks. Migrate direct/indexed gene-resolution callers to the same resource interface.

Validation: cover the complete approved condition–outcome table through real Python/CLI paths. Use independent SNP membership and interval-boundary expectations, matching BED/gene equivalence, both baseline layouts, mixed empty/usable queries, annotation-grid support, diagnostic ordering/provenance, overwrite/failure markers, explicit `CM=NA`, and downstream reloadability. Include duplicates crossing resolver chunks, streamed source failures, many-pathway selection lifetimes, and bounded BED preparation. Extend `tests/test_annotation.py`, `tests/test_gene_list_resolver.py`, and applicable output/workflow tests.

Exit checkpoint: real standalone BED and gene-list runs emit reloadable canonical shards and return live output-backed handles. Matching BED and gene intervals produce independently expected SNP memberships. Gate failures publish the required diagnostics with existing failure-marker behavior; the all-focal-skipped case publishes no new scientific query family.

### 3. Batch direct LD projection within one traversal and restrict score rows

Status: direct workflow cutover implemented and verified. File, BED, gene/control, and synthetic-base routes use output-contained source shards. Scope validation reuses prepared source evidence. Reference preparation maps selected reads; query batches reuse each correlation block; counts/overlap/classification are tiled; only output SNPs own score buffers. The scheduler bounds active plus pending work by worker count. Completed direct results detach gene diagnostics to persistent outputs before scratch closes. A one-worker lifetime test verifies release before the next chromosome. Focused workflow, scope, parallelism, gene integration, annotation, output, failure-marker, and configuration checks passed 324 tests and 29 subtests. Existing independent kernel/index checks passed in the preceding commits. Remaining old explicitly materialized Python bundle/test helpers await final shared API retirement after index and quantile migration. Resource trade-off measurements remain in slice 7.

Likely areas: [ldscore_calculator.py](../../src/ldsc/ldscore_calculator.py), [config.py](../../src/ldsc/config.py), [_kernel/ref_panel.py](../../src/ldsc/_kernel/ref_panel.py), [_kernel/ldscore.py](../../src/ldsc/_kernel/ldscore.py), [_kernel/plink_bed.py](../../src/ldsc/_kernel/plink_bed.py), [_kernel/overlap.py](../../src/ldsc/_kernel/overlap.py), and annotation semantics/count helpers.

Add and propagate positive-integer `query_batch_size=1000` through direct CLI/Python/configuration and worker arguments. Resolve the existing output/regression selection plus declared exclusions before score allocation, retaining supported unrestricted output behavior. Keep full contributor metadata and LD neighborhoods separately from output-row mappings. For every PLINK correlation block or decoded parquet-R2 chunk, project shared baseline/control columns and the regression-weight mask once, then project all query batches against that same block before releasing it. Do not wrap a batch loop around a whole-chromosome kernel or reread/reset the numerical traversal per query batch. Batch counts, advisory classification, baseline-query overlaps, and query diagonals without full-width float64 annotation/statistics temporaries.

Allocate scores only for output rows in float64; keep contributions from eligible non-output reference SNPs, both symmetric pair directions, diagonals, reference counts, and the distinct `w_ld` contributor universe. Preserve existing operation-specific count/overlap precision instead of imposing one new reduction dtype. Write full-reference sidecars while prepared metadata is alive, then return output metadata plus compact reference summaries. Bound active plus pending chromosome work by the resolved worker count, submit descriptors, and consume completed results promptly. Preserve ordered aggregation into the existing aggregate LD-score files without avoidable whole-result copies.

Validation: compare with an independent dense `(R @ A)[output_rows]` oracle containing a non-output SNP contributing to an output SNP. Cover both backends, quantitative/binary annotations, baseline-only/synthetic base, `gene_control` retaining its baseline role and output position, custom output restrictions, no restriction, region policies, `w_ld`, and all/common counts/overlaps plus pre/post-exclusion count metadata. Compare batch sizes 1, a non-divisor, and at least query count; compare sequential/parallel outputs. Instrument traversals and allocation shapes to detect repeated LD work, all-reference/all-query score buffers, and eager query conversions, including oversized Parquet row groups. Run focused reference, LD-score, parallelism, annotation-semantics, overlap, and affected index-kernel tests.

Exit checkpoint: each supported direct route works through the real entry points and reloads from the aggregate output directory. Both backends pass the independent oracle and traversal/allocation checks. Worker failure closes readers and owned scratch, and a one-worker run releases chromosome workspace before the next load. Migrate shared kernel callers and their tests with any changed signature; preserve the index constructor's precision before accepting that cutover.

### 4. Load gene-index chromosomes on demand and batch assembly

Status: not started. Depends on slices 1–3 for shared storage, query resolution, and kernel interfaces.

Likely areas: [gene_ldscore_index.py](../../src/ldsc/gene_ldscore_index.py), shared index kernels, query resolution/status helpers, and the canonical LD-score writer.

Represent the loaded index using shared metadata and component paths. Validate integrity through bounded reads and accumulate compact gene-support summaries. Assembly may reload each validated chromosome once; retain its sparse operator through all query batches for that chromosome and release it before proceeding. Apply the same positive-integer `query_batch_size=1000` contract to indexed CLI/Python assembly. Migrate construction to consume annotation shards while preserving intentional atom batching, float64 stored operator/common payloads, and the existing construction staging/publication exception. Avoid repeated operator loading per query batch and keep final LD-score files aggregate.

Validation: existing immutable-index corruption/identity cases still fail correctly; direct/indexed scientific equivalence, query batch invariance, operator-load counts, chromosome release, canonical schemas, and construction recovery/overwrite behavior remain covered by `tests/test_gene_ldscore_index.py` and `tests/test_gene_ldscore_index_kernel.py`. Measure validation I/O separately from assembly instead of assuming lazy loading is free.

Exit checkpoint: a multi-chromosome index is validated and assembled without a dictionary retaining every operator. Construction/recovery tests and direct/indexed scientific comparisons pass, including stored float64 payload checks. Run the broader repository suite at this shared annotation/kernel migration milestone.

### 5. Share batch-regression preparation and stream per-fit details

Status: not started. Depends on slice 1's shared lifecycle primitives. Existing canonical LD directories can exercise this slice independently of new LD-score generation.

Likely areas: [regression_runner.py](../../src/ldsc/regression_runner.py), [outputs.py](../../src/ldsc/outputs.py), selective LD-score loading, configuration, and partitioned-h2 result interfaces.

Prepare common SNP alignment, trait data, and shared baseline LD scores once; read only needed query columns in bounded batches from aggregate LD tables using positive-integer `query_batch_size=1000`. Preserve loading-time schema/provenance/allele/count-policy validation, fit-time query alignment, and selected-column numeric validation at their existing boundaries. Fit each query separately over its complete retained genome-wide SNP set; preserve the no-query functional regime as one joint baseline fit. Recompute model-dependent filtering, weights, regression design, and jackknife quantities. Write category tables, coefficient delete values, and per-query metadata privately as each model finishes, then release those detailed objects. Assign sorted canonical folders only after the batch succeeds. Retain summary rows, compact bookkeeping, and valid persistent output references.

Migrate every consumer of the changed LD-directory loader, including `h2` and `rg`, to request only its needed columns while preserving its current estimator and output behavior. This is integration of the selective reader, not a new statistical optimization scope for those commands.

Validation: compare every model with independent single-query execution, including model-specific zero/filtered columns and distinct SNP retention. Cover no-query functional fitting, unselected query values, missing schema columns, and query alignment at the existing validation boundary. Assert unchanged ordering, category tables, delete-value precision, diagnostics, canonical per-query paths, and overwrite/failure behavior. Inject a late fitting failure and publication failure to verify private details are cleaned, required diagnostics survive, and no stronger rollback promise is introduced. Check selective query reads and shared-preparation counts; verify no collection retains all detail matrices. Run `tests/test_regression_workflow.py`, `tests/test_regression_fit_outcomes.py`, `tests/test_kernel_regression.py`, and `tests/test_output.py` as applicable.

Exit checkpoint: a multi-query run agrees with the corresponding separate fits, preparation occurs once, and detailed objects do not survive their private writes. Published folder ordinals match every supported final-summary sort. `h2` and `rg` smoke checks confirm that changing the shared loader does not change their behavior.

### 6. Accumulate exact global quantile statistics from shards

Status: not started. Depends on slice 1's annotation access/lifetime work and the relevant selective-loading and persistent-fit interfaces from slice 5.

Likely areas: [quantile_h2.py](../../src/ldsc/quantile_h2.py), shared annotation/reference readers, numerical quantile/statistics helpers, and diagnostics/output integration.

Validate the common reference universe and retain the permitted eligible target-value vector plus sorting workspace to compute the existing exact global quantile boundaries. Preserve external target float64 values and raw-token exclusion semantics separately from fitted-annotation float32 storage. Keep external-target/reference duplicates fatal across the full dataset, including global uniqueness for omitted-allele inference. Reuse global boundaries over chromosome shards. Accumulate fitted annotation sums per quantile, overlap-validation products, and stable full-common-universe variance statistics; preserve target missingness, tie handling, empty-bin behavior, aggregate checks, and standardized coefficients. Do not create a whole-genome fitted-annotation matrix or dense SNP-by-quantile indicator matrix. Keep complete alignment diagnostics streamed and selected fitted columns in model order.

Validation: compare deterministic legacy quantile fixtures and an independent small dense calculation, including ties spanning chromosomes, float64 target differences that float32 would collapse, raw missing-token cases, cross-chromosome duplicate targets/reference identities, omitted-allele inference, empty bins, high-offset low-variance annotations, standardization over the full common universe, and deliberate count/overlap mismatches. Verify identical boundaries/membership and existing numeric tolerances, plus bounded shard lifetimes and no dense indicator allocation. Run `tests/test_quantile_h2.py` and affected annotation/output/regression tests.

Exit checkpoint: real baseline-only and single-query fitted-result inputs reproduce exact quantile membership and the expected quantile/standardized-coefficient tables. The target vector is released after boundary selection; annotation matrices remain chromosome-scoped, and full-common statistics still include target-excluded rows.

### 7. Verify workflow memory and finish user documentation

Status: not started. Depends on slices 1–6; collect targeted resource evidence earlier as each seam becomes available.

Complete the contract phase: remove migrated eager fields, obsolete BED-only wrappers, superseded readers/result collections, and temporary internal migration routes after a consumer search confirms they are unused. Update public exports and all remaining callers/tests/examples; do not retain compatibility aliases. Audit shared caches, query-resolution audits, alignment/drop diagnostics, writer conversions, completed futures, and returned objects across all workflows. Ensure each large payload has an explicit owner and release point and diagnostic completeness does not require reloading every row into result objects.

Finish matched local measurements against post-fingerprint-removal HEAD `a505c45` across chromosome counts, query counts, query batch sizes, and worker counts. Capture the baseline measurements before changing the corresponding paths, or reproduce them in an isolated local checkout without disturbing the working tree. Cover direct PLINK/parquet, standalone annotate, indexed assembly, batch regression, and quantile reconstruction. The baseline revision has no standalone gene-list command: compare that new route with a matching BED projection and label the differing resolution work explicitly; do not describe it as a same-command speed comparison. Existing commands use the same supported baseline arguments, with the new batching flag applied only to the new implementation.

Include a scaled 1,000-query local workload if practical without implying authorization for a full production/HPC run. Measure process-tree peak memory, elapsed runtime, and peak temporary disk separately from persistent output size; record dimensions, deterministic input/seed information, environment, stage boundaries, and the measurement method. Do not combine separate child-process peak measurements as though they occurred simultaneously. Compare fixed per-chromosome size as chromosome count grows, and fixed output dimensions as batch size changes, so accepted materialized-output memory is not mistaken for an annotation leak. There is no arbitrary RAM/runtime pass threshold: require bounded-retention evidence and report measured trade-offs honestly. Link the resulting local measurement report from this plan, clearly distinguishing measured results from projections.

After implementation, update affected public docstrings, exports/help, README, current architecture/data-flow/input/output contracts, troubleshooting, examples/tutorials, and user wiki. Explicitly explain the intended use: testing many pathways at once, for example 1,000 pathways in one run, with each pathway tested separately against baseline categories. Explain query batching, sequential memory release versus allocator RSS, output-contained scratch and its construction exception, and aggregate HM3 LD files. Preserve concurrent edits to existing wiki pages.

Final verification: run relevant focused suites and real CLI/Python smoke workflows, inspect generated/reloaded artifacts and failure diagnostics, then run the full `pytest` suite followed sequentially by `python -m unittest discover -s tests -p 'test*.py' -v`. Do not run the two suites concurrently because shared BED temporary cleanup can interfere. Verify `ldsc --help`, `python -m ldsc --help`, and affected command help, run `git diff --check`, and report actual results. Update this plan with evidence rather than carrying forward commit-reported prior outcomes.

Exit checkpoint: all migrated workflows satisfy the specification through production entry points; the final API has no eager compatibility path; scientific/artifact and memory-lifetime checks pass; resource measurements and the requested pathway-focused documentation are recorded.

## Validation milestones and progress tracking

For each slice, start with its smallest meaningful failing behavioral or numerical check and extend the affected focused suites as the path becomes runnable. Use existing independent oracles and immutable fixture resources; do not regenerate goldens or loosen tolerances to match a changed implementation. The focused numerical seams include `tests/test_reference_preparation.py`, `tests/test_plink_io.py`, `tests/test_ldscore_workflow.py`, `tests/test_ldscore_parallelism.py`, `tests/test_overlap_matrix.py`, and the slice-specific suites above.

Run broader repository validation after the shared annotation/kernel/index cutover and after all workflow migrations. Additional broad reruns need a new change, failure, or unresolved concern. Record the exact revision, commands, outcomes, and remaining limitations at each milestone. Use the existing editable development environment; this plan does not require a new dependency stack, external downloads, or benchmark jobs on HPC.

Keep the seven slice statuses current as implementation proceeds. Mark a slice complete only when its exit checkpoint passes, link any measurement evidence, and record newly discovered dependencies before continuing. Revise implementation details when evidence warrants, while preserving the confirmed goal, contracts, and scientific invariants.

## Risks and revision points

- Column batching can still retain full-width data through input parsing, projection, overlap helpers, pandas copies, or completed futures. Validate allocation shapes and ownership through the complete workflow.
- Moving the output mask upstream can silently remove reference contributors or alter `w_ld`; the independent dense oracle and explicit universe/result interfaces are required before kernel replacement.
- Global duplicate cleanup deliberately changes malformed sharded-input behavior; preserve every other stage-specific validation and source alignment rule.
- Smaller matrix calls and staging I/O may reduce throughput. Tune only after matched measurements, preserving the approved default of 1000 unless the user changes it.
- Source-backed access assumes original inputs remain available and unchanged. Verify global validation gate consequences and deterministic diagnostic ordering through bounded passes; do not add routine hashes or a persistent cache as an implicit remedy.
- Result-dependent folder numbering requires final path assignment after summary sorting. The approved private per-fit staging must not leak provisional paths into public manifests or returned results.

## Out of scope

Public chromosome-sharded LD-score tables; approximate quantiles; changed regression estimators, model universes, scientific thresholds, or accumulation precision; restored annotation fingerprints; portable copies of every baseline source; persistent scratch caches/resume; redesigned gene-index publication; unapproved standalone annotate input modes/controls; HPC benchmarks; unrelated repository changes.

## Verification record

Planning is grounded in local `restructure` at `a505c45`, the confirmed specification and decision records, current source/test consumers, the completed C1–C3 cleanup disposition, and repository guidance. The plan adds explicit internal migration/cutover rules, slice dependencies and exit checks, selective-reader integration for existing consumers, and a valid comparison for the new standalone gene-list route. No implementation, numerical tests, or resource benchmarks have run for this effort. On 2026-09-10, documentation checks passed for five related documents and 62 local links/anchors, update dates, trailing whitespace, and referenced test files. All seven slices have not-started status, validation instructions, and exit checkpoints. Source/test files remain unchanged; implementation evidence belongs here as slices complete.
