# Annotation memory completion review

Last updated on: 2026-09-10

The agreed work is **not fully complete** at `8ada67cdb89104c563b44937a144f64716f01d16` on `restructure`. The main architecture is implemented, but the earlier all-complete statement was premature. This review compares the conversation's final agreements, the [specification](../../specs/2026-09-10-annotation-workflow-memory-design.md), the [decision map](../../current/annotation-memory-decisions.md), the [standalone gene contract](../../current/annotate-gene-list-decisions.md), and the direct-query-batching handoff with the actual source and new probes. Later scope expansions take precedence over the original narrower request.

This audit changes documentation and adds reproducible evidence. It does not repair production code. The remaining implementation requirements below are already within the approved scope.

## Agreement coverage

| Approved requirement | Assessment at the audited revision | Evidence or remaining work |
| --- | --- | --- |
| One lazy dataset object; chromosome descriptors; no automatic whole-genome annotation matrix | Implemented | `_annotation_bundle.py`, `_annotation_storage.py`, and `test_annotation_storage.py`; detached reads and explicit closure are covered. |
| Bounded discovery/preflight for whole-genome and sharded sources; global drop-all duplicate cleanup | Implemented | `_annotation_sources.py` and `_annotation_identity.py`; chunk/alignment/allele/cross-chromosome cases are covered. No annotation fingerprints were restored. |
| Output-contained owned scratch; valid returned artifacts; original baseline dependencies allowed | Implemented | `AnnotationWorkspace` and `_annotation_outputs.py`; standalone construction scratch closes before returning persistent query descriptors. The approved index-construction staging exception remains. |
| Standalone BED and gene-list `annotate`, canonical shards, exact resolution, exclusions, Gate A/B, annotation-grid support, reloadability | Core behavior implemented; memory completion remains open | `annotate_workflow.py`, `_annotation_queries.py`, and `test_annotate_streaming.py`; shared coverage diagnostics still retain all affected identifiers (R3). |
| Sequential chromosome release; parallel active/submitted work bounded by workers | Implemented | `LDScoreCalculator._run_chromosomes`, `PreparedChromosome.close`, and sequential/parallel lifetime tests. `--threads` selects chromosome processes. Freed allocations can be reused even when RSS stays high. |
| Direct PLINK/R2 batching, default 1000, one LD traversal, shared baseline/control/weights, full reference contributors, output-row float64 accumulators | Projection implemented; statistics completion remains open | `_kernel/ldscore_projection.py` and independent dense-oracle tests; quantitative annotation counts depend on batching and can fail unchanged aggregate validation (R1). |
| Aggregate HM3 output files and permitted materialized final tables | Implemented | Canonical `ldscore.baseline.parquet` and `ldscore.query.parquet` remain aggregate. No chromosome-sharded public LD tables were added. |
| Direct aggregation avoids join/concatenate/split copies and a second retained chromosome-table collection | Not completed | `LDScoreCalculator._aggregate_chromosome_results` still uses both patterns (R2). |
| Lazy gene-index validation/assembly; one operator serves every batch for its chromosome; preserved integrity/construction policy | Operator design implemented | `LoadedGeneLDScoreIndex`, `_load_gene_ldscore_index`, and `assemble_indexed_ld_scores`; shared gene-coverage diagnostics still need R3. Runtime/I/O trade-offs were measured. |
| Batch `partitioned-h2` shares preparation, selectively loads queries, preserves separate whole-genome fits and model-dependent calculations | Implemented, except diagnostic retention | `RegressionRunner.prepare_inputs`, `_read_query_batch`, and `estimate_partitioned_h2_batch`; legacy projection drops remain resident across fitting (R3). |
| Per-fit details written/released incrementally; final sorted publication and persistent result paths | Implemented | `test_regression_streaming.py` covers detail lifetimes, output ordering, and late failures. |
| Exact global quantiles; bounded annotation passes; full-common stable statistics; no dense SNP-by-quantile matrix | Implemented; upstream compatibility checkpoint reopened | `_quantile_inputs.py` and `test_quantile_streaming.py`; the preserved count gate exposes R1. Do not weaken its tolerances. |
| Complete diagnostics streamed; only counts/examples/compact metadata retained | Partially implemented | Identity/gene row audits and quantile alignment diagnostics are streamed; two remaining paths are confirmed in R3. |
| Separate peak memory, runtime, temporary disk measurements; 1000-pathway use case; docs/wiki/tutorials/developer design | Substantial work completed; final sign-off reopened | [Resource report](results.md) includes matched local synthetic measurements; [developer design](../../current/annotation-memory-design.md) and tutorials explain separate baseline-adjusted models. Existing fixtures did not establish the missing guarantees above. |

## R1 — Quantitative annotation counts change beyond tolerance

Priority: high. In [`annotation_statistics()`](../../../src/ldsc/_kernel/overlap.py), each tile is summed in float32 and added into float32 count vectors. Tile height depends on query batch width. Detached selected reads also have a different memory layout from the earlier column-major annotation table. Those changes alter reduction order enough to exceed the agreed `rtol=1e-6`, `atol=1e-8` validation tolerance.

The public `run_ldscore()` probe used 6000 reference SNPs, 1200 output SNPs, two baseline columns, and 1000 queries on the same parquet-R2 reference. One query contained signed float32 values from `-np.random.default_rng(182).normal(size=6000)`; the other queries were unchanged binary inputs. All reference MAF values satisfy the common threshold.

| Quantity | Value |
| --- | ---: |
| Independent float64 sum of normalized query values | 25.222538138448726 |
| Original revision common count | 25.222549438476562 |
| Refactored default batch-size-1000 common count | 25.222614288330078 |
| Refactored relative error | 3.0191e-6 |
| Maximum absolute difference in written query LD scores | 0 |

SNP identities and order also match exactly. Passing the original artifact's counts and overlap to the existing `_validate_aggregates()` succeeds; passing the refactored artifact fails with an instruction to resupply the original annotations, even though the annotation values are unchanged. This probes the actual quantile validation seam; it is not a claim that a complete regression-plus-quantile run was executed during this audit.

A separate positive quantitative input, `1e6 + np.arange(6000) % 17`, gives all-reference sums of 6,000,048,128 before the refactor, 5,999,816,192 with batch size 32, and 6,000,163,328 with batch size 1000. The old common-count reduction was already inaccurate for this high-offset example; the signed example above isolates a newly failing compatibility case.

The hypotheses were distinguished as follows: staged-value corruption is ruled out by exact `ColumnStore` round-trip equality; output-row/reference selection changes are ruled out by identical rows and LD values on the same contributor universe; MAF changes are ruled out by the fixed all-common input; float32 reduction layout and batch-dependent partial sums reproduce the discrepancy with fixed values.

Reproduce the reduced case without the benchmark fixtures:

```bash
PYTHONPATH=src python docs/audits/annotation-memory/completion_probe.py OUTPUT_DIRECTORY
```

The [probe](completion_probe.py) creates and cleans scratch inside the specified output directory. [Recorded evidence](completion-evidence.json) includes the reduced case, public workflow results, and linked count/overlap validation. Completion requires stable count accumulation across query batches and row tiles, quantitative regression coverage, and unchanged scientific tolerances.

## R2 — Direct results retain duplicate LD tables

Priority: medium. [`_aggregate_chromosome_results()`](../../../src/ldsc/ldscore_calculator.py) still joins each chromosome's baseline/query tables, concatenates them, sorts the joined table, and splits it again. The returned `LDScoreResult.chromosome_results` also keeps every original chromosome table. `_replace_result_output_paths()` changes diagnostic references but retains these table payloads.

The public probe returned 1200 rows in the aggregate query table and another 1200 in retained chromosome query tables. `np.shares_memory()` was false for the same query column. This behavior existed before the refactor, but removing it was explicitly approved in the specification's “Direct LD-score interfaces and traversal” section. Permission to retain the final HM3 tables does not satisfy that separate requirement. Completion requires direct baseline/query aggregation and compact completed-chromosome bookkeeping while preserving aggregate public files.

## R3 — Two diagnostic paths retain full records

Priority: medium. [`StagedGeneListBatch.with_coverage()`](../../../src/ldsc/_gene_query_storage.py) retains every uncovered gene ID in each summary's `uncovered_gene_ids` string and again in a collected error string. A probe with 20 lists of 2000 out-of-scope genes retained 639,980 identifier characters in summary cells plus 684,050 characters in errors. Ten equivalent lists retained about half as much. This is retention proportional to all source-gene memberships during a failure gate, rather than counts and bounded examples. Complete per-gene records already have a streamed audit route.

[`RegressionRunner.prepare_inputs()`](../../../src/ldsc/regression_runner.py) stores the full legacy-sumstats projection drop DataFrame in the shared dataset. The CLI's `run_partitioned_h2_from_args()` separately keeps `legacy_drops` until all fits finish and then writes the audit. A preparation probe with only 50 aligned rows retained 5000 diagnostic rows occupying 1,749,132 bytes. This affects the explicitly included batch-regression workflow, even though optimizing unrelated single-trait estimators was outside scope.

Completion requires staging these complete diagnostics, preserving deterministic records and current output/failure contracts, and retaining only counts/examples/path references during subsequent work. Add failure-path and legacy-batch lifetime tests.

## Documentation and verification status

This audit corrects the plan/progress “all complete” claims, the specification's stale “implementation has not started” status, and the developer document's inaccurate claim that positive worker requests are capped by CPU availability. Positive requests are capped by chromosome count; negative requests derive their count from available CPUs. The developer design now identifies the unresolved implementation limitations.

The six focused streaming/storage suites passed **73 tests, 20 warnings, in 5.63 seconds** during this audit. Existing full-suite evidence remains **1494 pytest tests passed, one skipped, 132 subtests**, followed by **1000 unittest tests run, one skipped, OK** at behavioral revision `91ad7bf`. The full suites were not rerun for this documentation-only audit. The 21 prior matched artifact comparisons remain valid for their measured fixtures; they do not cover the counterexamples discovered here. No HPC run, external publication, or production source change occurred.

Finish R1 first, then R2 and R3, extend the relevant allocation/numerical/failure-path checks, repeat affected resource comparisons, and run the full suites sequentially before reinstating completion. These are implementation gaps under existing decisions, not new design questions.
