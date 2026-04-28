# Code Structure

This is the contributor-facing module map for `ldsc_py3_restructured`.

## Repository Map

```text
ldsc_py3_restructured/
├── docs/
├── src/ldsc/
│   ├── __init__.py
│   ├── __main__.py
│   ├── cli.py
│   ├── config.py
│   ├── path_resolution.py
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
| `ldsc.path_resolution` | normalize path tokens, resolve concrete input files, create output directories, and preflight fixed output paths |
| `ldsc.column_inference` | resolve header aliases and normalize identifier/build tokens |
| `ldsc.chromosome_inference` | canonical chromosome normalization and ordering |
| `ldsc.genome_build_inference` | public `chr_pos` build and coordinate-basis inference helpers |
| `ldsc.annotation_builder` | public annotation bundle loading and BED projection surface |
| `ldsc.ref_panel_builder` | parquet reference-panel build workflow |
| `ldsc.ldscore_calculator` | LD-score orchestration, optional synthetic `base` annotation construction, aggregation, and output routing |
| `ldsc.sumstats_munger` | raw-sumstats munging wrapper, canonical `CHR`/`POS` sumstats output, metadata sidecar handling, and curated sumstats loader |
| `ldsc.regression_runner` | file-driven regression dataset assembly, `SNP` or `CHR:POS` merging, and estimator dispatch |
| `ldsc.outputs` | artifact naming, layout, and serialization |
| `ldsc._kernel.annotation` | annotation parsing, validation, BED intersection |
| `ldsc._kernel.ref_panel_builder` | optional genetic-map parsing, optional liftover, parquet schemas, pairwise LD emission |
| `ldsc._kernel.ref_panel` | runtime PLINK/parquet reference-panel adapters |
| `ldsc._kernel.ldscore` | LD-score math and legacy-compatible computation helpers |
| `ldsc._kernel.sumstats_munger` | legacy-compatible raw summary-statistics QC and normalization |
| `ldsc._kernel.regression` | LDSC estimators for `Hsq` and `RG` |
| `ldsc._kernel._jackknife`, `ldsc._kernel._irwls` | supporting numerical routines used by regression |
| `ldsc._kernel.formats`, `ldsc._kernel.identifiers` | file-format readers and SNP identifier helpers |

## Where To Change Code

| Goal | Start here |
| --- | --- |
| change CLI flags or subcommand wiring | `src/ldsc/cli.py` |
| change path-token behavior | `src/ldsc/path_resolution.py` |
| change output collision policy | `src/ldsc/path_resolution.py`, then the workflow writer that owns the artifact |
| change header aliases or identifier/build normalization | `src/ldsc/column_inference.py` |
| change automatic `chr_pos` genome-build inference | `src/ldsc/genome_build_inference.py` |
| change annotation loading or BED projection | `src/ldsc/annotation_builder.py`, then `src/ldsc/_kernel/annotation.py` |
| change parquet reference-panel build logic | `src/ldsc/ref_panel_builder.py`, then `src/ldsc/_kernel/ref_panel_builder.py` |
| change runtime PLINK/parquet reference access | `src/ldsc/_kernel/ref_panel.py` |
| change LD-score orchestration, optional-baseline behavior, or output packaging | `src/ldsc/ldscore_calculator.py`, `src/ldsc/outputs.py` |
| change LD-score math | `src/ldsc/_kernel/ldscore.py` |
| change raw sumstats ingestion, `CHR`/`POS` handling, sidecar provenance, or curated loading | `src/ldsc/sumstats_munger.py`, then `src/ldsc/_kernel/sumstats_munger.py` |
| change regression dataset assembly or CLI summaries | `src/ldsc/regression_runner.py` |
| change LDSC estimators | `src/ldsc/_kernel/regression.py` |
| change LD-score result-directory files | `src/ldsc/outputs.py` |

## Architectural Rules That Matter In Practice

- Treat `src/ldsc/` as the only supported Python import surface.
- Do not add user-facing path discovery to `_kernel`; pass concrete files in.
- Keep public file contracts for `.annot(.gz)`, `.sumstats.gz` plus `sumstats.metadata.json`, and canonical LD-score result directories stable unless the change is intentional and coordinated. Legacy `.l2.ldscore(.gz)`, `.w.l2.ldscore(.gz)`, `.l2.M`, and `.l2.M_5_50` files are compatibility concerns rather than the public LD-score output surface.
- Keep optional-baseline behavior in the public LD-score workflow layer: no baseline and no query means a synthetic all-ones `base`; query annotations require explicit baseline annotations.
- Reuse `ensure_output_paths_available()` for fixed output files. Public
  workflows should create missing output directories, reuse existing
  directories, fail on existing known files by default, and require
  `--overwrite` or `overwrite=True` for replacement.
- Preflight deterministic output paths before expensive or multi-file writes.
  This is especially important for `build-ref-panel`, `munge-sumstats`,
  `annotate`, and summary-table regression commands.
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
