# Design Map

This file maps durable design documents to the implementation modules they
describe. Update it when workflow ownership, public entry points, or artifact
contracts move.

## Annotation Workflow

| Design document | Implementation |
| --- | --- |
| `docs/current/architecture.md` | `src/ldsc/annotation_builder.py`, `src/ldsc/_kernel/annotation.py`, `src/ldsc/cli.py` |
| `docs/current/code-structure.md` | public workflow boundary in `src/ldsc/annotation_builder.py`; private helper boundary in `src/ldsc/_kernel/annotation.py` |
| `docs/current/data-flow.md` | `AnnotationBuilder.run()`, `AnnotationBuilder.project_bed_annotations()`, `run_bed_to_annot()`, `run_annotate_from_args()`, `annotation_builder.main()` |
| `docs/current/io-argument-inventory.md` | annotate parser and namespace dispatch in `src/ldsc/annotation_builder.py`; top-level dispatch in `src/ldsc/cli.py` |
| `docs/current/layer-structure.md` | annotation row in the public workflow and kernel layer matrix |
| `docs/current/path-specification.md` | annotation path-token handling through `resolve_file_group()`, `split_cli_path_tokens()`, and output preflight helpers |

## Annotation Tutorials

| Tutorial | Covered implementation path |
| --- | --- |
| `tutorials/ld-score-calculation.md` | direct BED projection through `run_bed_to_annot()` and `ldsc annotate`; in-memory BED projection during `run_ldscore()` |
| `tutorials/partitioned-ldsc.md` | partitioned workflow using `AnnotationBuilder.run()` for in-memory query annotations and optional `ldsc annotate` materialization |
| `tutorials/cell-specific-ldsc.md` | cell-type BED annotations through `AnnotationBuilder.run()`, `ldsc ldscore`, and optional `run_bed_to_annot()` materialization |

## Current Boundary

- Public annotation workflow ownership is in `ldsc.annotation_builder`.
- `ldsc._kernel.annotation` contains only low-level table and BED helpers.
- `ldsc.cli` registers annotate flags from `annotation_builder` and dispatches
  parsed namespaces to `run_annotate_from_args()` without reparsing.
- `main_bed_to_annot()` is intentionally removed; `annotation_builder.main()`
  is the supported parser entry point.

## Workflow Logging

| Design document | Implementation |
| --- | --- |
| `docs/current/workflow-logging.md` | `src/ldsc/_logging.py`; workflow wrappers in `annotation_builder.py`, `ldscore_calculator.py`, `ref_panel_builder.py`, `sumstats_munger.py`, and `regression_runner.py` |
| `docs/superpowers/specs/2026-05-02-logging-harmonization-design.md` | final logging harmonization decisions and result-contract boundaries |
| `docs/superpowers/plans/2026-05-02-logging-harmonization-implementation-plan.md` | completed implementation checklist and verification commands |
| `docs/superpowers/plans/2026-06-04-logging-console-file-routing-plan.md` | console/file routing: file-authoritative sink, console error-only echo, no-output console fallback, run-aborting traceback capture; touches `src/ldsc/_logging.py`, `src/ldsc/cli.py`, `src/ldsc/_kernel/ldscore.py` |

## Liftover Harmonization

| Design document | Implementation |
| --- | --- |
| `docs/current/liftover-harmonization-decisions.md` | active liftover contracts for `munge-sumstats` and `build-ref-panel`; shared drop reasons and dropped-SNP sidecar schema |
| `docs/current/architecture.md` | public workflow responsibilities for `ldsc.sumstats_munger`, `ldsc.ref_panel_builder`, and shared `_kernel.liftover` helpers |
| `docs/current/data-flow.md` | sumstats liftover stage, reference-panel liftover/drop flow, and dropped-SNP audit outputs |
| `docs/current/io-argument-inventory.md` | removed `--duplicate-position-policy`, remaining liftover flags, output paths, and Python config fields |
| `docs/current/path-specification.md` | always-written sumstats dropped-SNP sidecar and reference-panel stale class-2 output-directory contract |
| `docs/superpowers/specs/2026-05-10-liftover-harmonization.md` | final harmonization design: drop-all only, unified sidecars, DEBUG-only examples, and class-2 hygiene contract |
| `docs/superpowers/plans/2026-05-10-liftover-harmonization.md` | implemented task plan and verification checklist for liftover harmonization |
| `docs/superpowers/specs/2026-05-11-ref-panel-stale-class2-warning.md` | additive stale class-2 reference-panel warning design |
| `docs/superpowers/plans/2026-05-11-ref-panel-stale-class2-warning.md` | implemented task plan for stale class-2 warning tests and helper |

## Liftover Tutorials

| Tutorial | Covered implementation path |
| --- | --- |
| `tutorials/build-parquet-reference-panel-from-plink.md` | reference-panel chain liftover, rsid-mode rejection, drop-all collision behavior, stale class-2 warning, and per-chromosome dropped-SNP sidecars |
| `tutorials/heritability-estimates.md` | `munge-sumstats` liftover outputs, count-only logs, and `dropped_snps/dropped.tsv.gz` audit sidecar |
| `tutorials/partitioned-ldsc.md` | `chr_pos` sumstats liftover identity and dropped-SNP audit sidecar in partitioned workflows |
| `tutorials/cell-specific-ldsc.md` | liftover audit behavior before cell-specific regression |
| `tutorials/cross-trait-genetic-correlation.md` | multi-trait `chr_pos` liftover audit behavior and coherent sumstats output family |

## Regression And RG Results

| Design document | Implementation |
| --- | --- |
| `docs/current/class-and-features.md` | public regression API surface, including `RegressionRunner.estimate_rg_pairs()`, `RgResultFamily`, `RgOutputConfig`, and `RgDirectoryWriter` |
| `docs/current/data-flow.md` | regression data flow from curated sumstats and canonical LD-score directories to h2, partitioned-h2, and rg output families |
| `docs/current/io-argument-inventory.md` | current sumstats trait-label metadata contract, `ldsc rg --sumstats-sources` CLI contract, `--anchor-trait` selection, and output-family flags |
| `docs/current/layer-structure.md` | ownership split between `regression_runner.py`, `outputs.py`, CLI dispatch, and workflow logging |
| `docs/current/path-specification.md` | scalar h2/partitioned-h2 sumstats path handling vs. group-style rg `sumstats_sources` resolution |
| `docs/superpowers/specs/2026-05-09-batch-rg-design.md` | final multi-trait rg design decisions and output schemas |
| `docs/superpowers/plans/2026-05-09-batch-rg-implementation-plan.md` | implemented rg refactor checklist and verification plan |

## Regression Tutorials

| Tutorial | Covered implementation path |
| --- | --- |
| `tutorials/cross-trait-genetic-correlation.md` | Python `estimate_rg_pairs()` flow, CLI `--sumstats-sources`, glob input, anchor mode, and rg output family |
| `tutorials/cross-trait-genetic-correlation.ipynb` | executable toy `estimate_rg_pairs()` and `ldsc rg --sumstats-sources` smoke path |

## Reference Panel (R²) Build & Read

| Design document | Implementation |
| --- | --- |
| `docs/current/parquet-r2-format-and-read-pipeline.md` | canonical R² parquet schema/metadata (incl. `ldsc:min_r2`, `ldsc:r2_encoding`/`ldsc:r2_scale`), int16 R² quantization + `BYTE_STREAM_SPLIT` (`_quantize_r2` in `_kernel/ref_panel_builder.py`), zstd-L9 column-chunk compression, and `SortedR2BlockReader` decode/window read path with int16→float32 dequant (`_resolve_r2_scale`) and vectorized endpoint-index lookup in `_kernel/ldscore.py` |
| `docs/superpowers/specs/2026-06-01-build-ref-panel-memory-optimization-design.md` | build-ref-panel memory/speed optimization: `--min-r2`, compact numpy pair buffers, vectorized pair extraction, builder-scoped float32 genotypes, direct PyArrow writes, BED bitarray releases (MAF-filter and per-build), zstd-L9, vectorized canonical decode lookup, and (#11) int16 R² quantization + `BYTE_STREAM_SPLIT` |
| `docs/superpowers/specs/2026-06-02-build-ref-panel-reader-streaming-design.md` | shared PLINK genotype reader extracted to `_kernel/plink_bed.py`; lazy header read + per-SNP selective decode with fused individual filter (O1/O3) and opt-in disk streaming for unrestricted builds (O2). Benefits `build-ref-panel`, the legacy `compute_chrom_from_plink` LD-score path, and `PlinkRefPanel`; never materializes the whole-chromosome bitarray |
| `docs/superpowers/specs/2026-06-03-pair-emission-columnar-and-doc-sync-design.md` | columnar `PairColumns` pair emission end-to-end (O6): `_pop_pair_rows_before`/`yield_pairwise_r2_rows` yield `(i, j, r2, sign)` chunks and `write_r2_parquet`/`_standard_r2_index_table` consume them, removing the per-pair dict round-trip; parquet format/encoding unchanged. Also the O7 memory doc sync across `docs/current/` and the build tutorial |

## LD-Score Calculation Parallelism

| Design document | Implementation |
| --- | --- |
| `docs/superpowers/specs/2026-06-05-ldscore-chromosome-parallelism-design.md` | single-machine cross-chromosome process-pool parallelism (parallel by default): `LDScoreConfig.parallel` (master switch, default `True`) + `num_workers` (default `0`=auto/all cores) with CLI `--no-parallel`/`--num-workers`, `_resolve_worker_count`, `LDScoreCalculator._run_chromosomes` (inline vs spawn `ProcessPoolExecutor`), module-level `_compute_one_chromosome` worker, `_ChromOutcome`, `_init_worker` BLAS-thread guard, deterministic chromosome-ordered aggregation. Output contract unchanged. All in `src/ldsc/ldscore_calculator.py` |
| `docs/superpowers/plans/2026-06-05-ldscore-chromosome-parallelism-plan.md` | TDD implementation checklist and the two-chromosome canonical index-panel fixture in `tests/test_ldscore_parallelism.py` |
