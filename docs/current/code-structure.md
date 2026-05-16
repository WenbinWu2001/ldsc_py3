# Code Structure

This is the contributor-facing module map for `ldsc_py3_Jerry`.

## Repository Map

```text
ldsc_py3_Jerry/
├── docs/
├── src/ldsc/
│   ├── __init__.py
│   ├── __main__.py
│   ├── cli.py
│   ├── config.py
│   ├── path_resolution.py
│   ├── _logging.py
│   ├── column_inference.py
│   ├── chromosome_inference.py
│   ├── genome_build_inference.py
│   ├── annotation_builder.py
│   ├── ref_panel_builder.py
│   ├── ldscore_calculator.py
│   ├── sumstats_munger.py
│   ├── regression_runner.py
│   ├── outputs.py
│   └── _kernel/
├── tests/
└── tutorials/
```

## Dependency Direction

- `ldsc.cli` -> public workflow modules
- public workflow modules -> shared helpers in `config.py`, `path_resolution.py`, `column_inference.py`, `chromosome_inference.py`, `genome_build_inference.py`
- public workflow modules -> private `_kernel` modules
- `ldsc.outputs` is called from the workflow layer, not from `_kernel`
- regression reloads written LD-score artifacts; it does not depend on annotation or reference-panel kernels directly

## Module Map

| Module | Responsibility |
| --- | --- |
| `ldsc.cli` | unified `ldsc` command and subcommand dispatch |
| `ldsc.config` | frozen public config dataclasses and basic validation |
| `ldsc.path_resolution` | normalize path tokens, resolve concrete input files, create output directories, preflight fixed output paths, and enforce coherent output artifact families |
| `ldsc._logging` | shared workflow logging context, LDSC logger level handling, lifecycle audit lines, and log-only formatting helpers |
| `ldsc.column_inference` | resolve header aliases and normalize identifier/build tokens |
| `ldsc.chromosome_inference` | canonical chromosome normalization and ordering |
| `ldsc.genome_build_inference` | public `chr_pos` build and coordinate-basis inference helpers |
| `ldsc.hm3` | public packaged curated HM3 map loader and installed map path helper for workflow internals |
| `ldsc.hm3_reference` | reproducible builder for the compact HM3 coordinate reference used by genome-build inference |
| `ldsc._kernel.liftover` | shared hg19/hg38 liftover helpers, chain-file translation, curated HM3 dual-build coordinate conversion, drop-all coordinate collision helpers, and readable drop reports |
| `ldsc.annotation_builder` | public annotation workflow: CLI args, parser entry point, path resolution, bundle loading, BED projection, and query `.annot.gz` writing |
| `ldsc.ref_panel_builder` | parquet reference-panel build workflow |
| `ldsc.ldscore_calculator` | LD-score orchestration, optional synthetic `base` annotation construction, aggregation, and output routing |
| `ldsc.sumstats_munger` | raw-sumstats CLI/API orchestration, `--format auto` / `--infer-only` header inference, Parquet/TSV curated output writing, root downstream `metadata.json`, diagnostics under `diagnostics/`, canonical `CHR`/`POS` sumstats output, and curated sumstats loader |
| `ldsc.regression_runner` | file-driven regression dataset assembly, active effective identity-key merging (`SNP`, `SNP:<allele_set>`, `CHR:POS`, or `CHR:POS:<allele_set>`), h2/partitioned-h2/rg estimator dispatch, and rg result-family writing |
| `ldsc.outputs` | artifact naming, LD-score parquet layout, partitioned-h2 per-query layout, rg result-family layout, metadata JSON payloads, and serialization |
| `ldsc._kernel.annotation` | low-level annotation table reading and BED intersection helpers |
| `ldsc._kernel.ref_panel_builder` | optional genetic-map parsing, optional liftover, parquet schemas, pairwise LD emission |
| `ldsc._kernel.ref_panel` | runtime PLINK/parquet reference-panel adapters |
| `ldsc._kernel.ldscore` | LD-score math and legacy-compatible computation helpers |
| `ldsc._kernel.sumstats_munger` | legacy-compatible raw summary-statistics QC, normalization, optional coordinate liftover, and `.sumstats.gz` writing without public CLI or log-file ownership |
| `ldsc._kernel.regression` | LDSC estimators for `Hsq` and `RG` |
| `ldsc._kernel._jackknife`, `ldsc._kernel._irwls` | supporting numerical routines used by regression |
| `ldsc._kernel.formats`, `ldsc._kernel.identifiers` | file-format readers and SNP identifier helpers |

## Where To Change Code

| Goal | Start here |
| --- | --- |
| change CLI flags or subcommand wiring | `src/ldsc/cli.py` |
| change path-token behavior | `src/ldsc/path_resolution.py` |
| change output collision policy | `src/ldsc/path_resolution.py`, then the workflow writer that owns the artifact |
| change workflow log files, lifecycle lines, or log-level handling | `src/ldsc/_logging.py`, then the workflow orchestration function that owns the output directory |
| change header aliases or identifier/build normalization | `src/ldsc/column_inference.py` |
| change automatic `chr_pos` genome-build inference | `src/ldsc/genome_build_inference.py` |
| change annotation loading, `ldsc annotate` behavior, or BED projection | `src/ldsc/annotation_builder.py`, then `src/ldsc/_kernel/annotation.py` |
| change parquet reference-panel build logic | `src/ldsc/ref_panel_builder.py`, then `src/ldsc/_kernel/ref_panel_builder.py` |
| change runtime PLINK/parquet reference access | `src/ldsc/_kernel/ref_panel.py` |
| change LD-score orchestration, optional-baseline behavior, or output packaging | `src/ldsc/ldscore_calculator.py`, `src/ldsc/outputs.py` |
| change LD-score math | `src/ldsc/_kernel/ldscore.py` |
| change raw sumstats ingestion, format inference, `CHR`/`POS` handling, sumstats SNP keep-list filtering, liftover drop audit sidecars, sidecar provenance, or curated loading | `src/ldsc/sumstats_munger.py`, then `src/ldsc/_kernel/sumstats_munger.py` |
| change regression dataset assembly or CLI summaries | `src/ldsc/regression_runner.py`, then `src/ldsc/outputs.py`, `docs/current/partitioned-h2-results.md` for partitioned-h2 output layout, and `docs/current/partitioned-ldsc-workflow.md` for rg output contracts |
| change LDSC estimators | `src/ldsc/_kernel/regression.py` |
| change LD-score result-directory files, parquet row-group layout, partitioned-h2 per-query layout, rg result-family layout, or metadata JSON payloads | `src/ldsc/outputs.py` |

## Architectural Rules That Matter In Practice

- Treat `src/ldsc/` as the only supported Python import surface.
- Do not add user-facing path discovery to `_kernel`; pass concrete files in.
- Keep public file contracts for `.annot(.gz)`, Parquet-first munged sumstats plus optional `.sumstats.gz` compatibility output, root sumstats `metadata.json`, canonical LD-score result directories, and regression summary directories stable unless the change is intentional and coordinated. LD-score parquet files remain flat files, with chromosome-aligned row groups documented through root `metadata.json`; diagnostic logs, dropped-SNP reports, and diagnostic metadata live under `diagnostics/`. Legacy `.l2.ldscore(.gz)`, `.w.l2.ldscore(.gz)`, `.l2.M`, and `.l2.M_5_50` files are compatibility concerns rather than the public LD-score output surface.
- Keep optional-baseline behavior in the public LD-score workflow layer: no baseline and no query means a synthetic all-ones `base`; query annotations require explicit baseline annotations. Regression treats that synthetic `base` path as an `h2`/`rg` input only; `partitioned-h2` requires explicit query LD-score columns.
- Keep `ldsc annotate` orchestration in `ldsc.annotation_builder`. The CLI
  registers annotation flags from that module and dispatches parsed namespaces
  to `run_annotate_from_args()`; `_kernel.annotation` must not own parser
  functions, result objects, or public workflow aliases.
- Use the shared output preflight helpers for fixed output files, including
  workflow log files. Public workflows should create missing output
  directories, reuse existing directories, fail on existing owned artifacts by
  default, and require `--overwrite` or `overwrite=True` for replacement.
- Keep coherent output families consistent. For `munge-sumstats`,
  `build-ref-panel`, `ldscore`, `partitioned-h2`, `rg`, and `annotate`,
  no-overwrite mode rejects any owned sibling from the workflow family.
  Overwrite mode writes the requested current outputs and then removes stale
  owned siblings that were not produced by the successful run. Unrelated files
  in the output directory must be preserved.
- Preflight deterministic output paths before expensive or multi-file writes.
  This is especially important for `build-ref-panel`, `munge-sumstats`,
  `annotate`, and summary-table regression commands.
- Treat workflow logs as audit artifacts. Do not include log paths in
  `output_paths` mappings or thin metadata sidecars that downstream code
  interprets as scientific data artifacts.
- Keep regression file-driven: it should be able to rebuild state from written artifacts without recomputing LD scores.
- Prefer extending shared helpers or the workflow-owned writer over duplicating local parsing or writing logic.

## Test Map

| Area | Main tests |
| --- | --- |
| output layer | `tests/test_output.py` |
| annotation workflow | `tests/test_annotation.py` |
| reference-panel builder | `tests/test_ref_panel_builder.py` |
| LD-score workflow | `tests/test_ldscore_workflow.py` |
| sumstats munging | `tests/test_sumstats_munger.py` |
| regression workflow | `tests/test_regression_workflow.py` |
| path and config contracts | `tests/test_path_resolution.py`, `tests/test_config_identifiers.py`, `tests/test_column_inference.py`, `tests/test_genome_build_inference.py` |
