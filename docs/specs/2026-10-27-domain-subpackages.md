# Domain subpackage reorganization of `src/ldsc`

Last updated on: 2026-10-27

Status: design confirmed by the maintainer on 2026-10-27, including the dead-code policy D1.
Implementation is not yet authorized.
This is a temporary in-flight specification; the commit sequence is in the [implementation plan](../plans/2026-10-27-domain-subpackages.md).
After implementation, move the durable contract into [code-structure.md](../dev/code-structure.md) and [architecture.md](../dev/architecture.md), then delete this file and the plan.

## Goal and scope

Reorganize the 87 modules of `src/ldsc/` from one flat directory plus a global `_kernel/` layer into private domain subpackages whose imports point in one direction only.
The public interface does not change: users import supported objects from `ldsc`, and the `ldsc` command keeps its subcommands, flags, outputs, log content and messages.

The reorganization does not:

- change numerical results, artifact formats, column names, CLI behavior or public names;
- add backward-compatibility aliases for old private module paths;
- split large workflow modules other than `gene_ldscore_index.py`;
- move CLI parsers out of the workflow modules that own them.

## Evidence

The evidence was collected at `118b398` and rechecked at `ee2886e`.

- `src/ldsc/` holds 66 root modules, 38 of them private helpers, plus 19 modules in `_kernel/` and 2 in `plotting/`.
- Filename prefixes act as namespaces: 12 `_annotation_*`, 5 `_ldscore_*`, 2 `_quantile_*` and 2 `_genome_build_*` modules.
- Four names exist at two levels with different meanings: `ref_panel_builder`, `sumstats_munger`, `r2_query` and `gene_ldscore_index`.
- No import cycle exists among module-level imports, but function-level imports form one 21-module cycle (for example `outputs` and `ldscore_calculator`) and a second one between `regression_runner` and `_partitioned_h2_parallel`.
- Nine of the 19 `_kernel` modules import root helpers, and [`_kernel/ref_panel.py`](../../src/ldsc/_kernel/ref_panel.py) imports `config`, `path_resolution` and `genome_build_inference`.
- [`_kernel/ldscore.py`](../../src/ldsc/_kernel/ldscore.py) mixes reference-panel preparation code (about 770 lines: the `PreparedChromosome` contract, LD-window helpers, identity and sample helpers, R² panel readers) with about 175 lines of LD-score math.
- The R² panel format is spread over four files: the writer in `_kernel/ref_panel_builder.py`, the readers in `_kernel/ldscore.py`, and the pair lookup in `r2_query.py` and `_kernel/r2_query.py`.
- [architecture.md](../dev/architecture.md) already defines the preparation boundary: `RefPanel.prepare_chromosome()` returns a `PreparedChromosome`, which `compute_chromosome()` consumes without reopening sources.

## Tiers and import rules

| Tier | Package or module | Contents |
| --- | --- | --- |
| 0 | `errors.py` | Exception classes |
| 1 | `_core/` | Scientific primitives: chromosomes, columns, coordinates, SNP identity, regions, PLINK companion files |
| 2 | `_runtime/` | Infrastructure: paths, input preflight gate, logging, progress, workers, artifact publication base |
| 3 | `config.py` | Public configuration dataclasses |
| 4 | `_genome/` | Genome-build inference and policy, liftover, HM3 map |
| 5 | `_annotation/`, `_refpanel/`, `_sumstats/` | Peer domains |
| 6 | `_ldscore/` | LD-score calculation, LD-score directory, gene LD-score index, legacy LD-score conversion |
| 7 | `_regression/` | Estimators, regression workflows, post-fit conversion and quantiles |
| 8 | `_plotting/` | Result plotting |
| 9 | `cli.py`, `__init__.py`, `__main__.py` | Command dispatch and the public namespace |

1. A module may import from its own package and from lower tiers only.
   Peers in tier 5 may not import each other.
2. The rule covers every import inside the package: module-level, function-level and `TYPE_CHECKING` imports.
   There is no allowlist.
3. Function-level imports remain allowed for optional dependencies such as Matplotlib and for SciPy-backed workflows that the public namespace loads lazily, but never to avoid an import cycle.
4. A domain `__init__.py` holds only a docstring that states the domain's purpose and tier; it imports nothing.
   Python runs a package's `__init__.py` before any of its modules, so a re-exporting `__init__.py` would import the whole domain and break the lazy loading added in `f321e0e`.
5. Imports between domains name the module, for example `from .._ldscore.source import LDScoreSource`.
6. Modules inside private packages have no leading underscore; the package underscore marks them private.
   Role names (`outputs.py`, `preflight.py`, `builder.py`) repeat across domains with the same meaning; different meanings get different names.
   Names that read as standard modules are avoided: `workflow_logging.py`, not `logging.py`, and `munge_input.py`, not `input.py`.
7. The package root keeps only `__init__.py`, `__main__.py`, `cli.py`, `config.py`, `errors.py` and `data/`.

A layering test parses every module under `src/ldsc/`, derives each module's tier from its package path, and fails on any import that breaks rules 1 and 2.
While modules are still being moved, a module inside a tiered package may not import a module that has not yet moved, because that module's tier is not yet defined.
The final version also fails if any module other than the root shell sits at the package root.

## Module map

```text
src/ldsc/
├── __init__.py  __main__.py  cli.py  config.py  errors.py  data/   (stay; imports updated)
│
├── _core/
│   ├── chromosomes.py        ← chromosome_inference.py
│   ├── columns.py            ← column_inference.py
│   ├── coordinates.py        ← _coordinates.py
│   ├── identifiers.py        ← _kernel/identifiers.py
│   ├── plink_files.py        ← _kernel/formats.py (PLINK BIM/FAM and keep-list readers)
│   ├── regions.py            ← _kernel/regions.py
│   ├── row_alignment.py      ← _row_alignment.py (+ sort_frame_by_genomic_position)
│   └── snp_identity.py       ← _kernel/snp_identity.py (+ VALID_SNPS)
│
├── _runtime/
│   ├── artifacts.py          ← outputs.py shared base: ArtifactFamily, _declare_artifacts,
│   │                           atomic JSON/table writes, _to_serializable, _result_metadata
│   ├── cli_help.py           ← _cli_help.py
│   ├── parallelism.py        ← _parallelism.py
│   ├── paths.py              ← path_resolution.py
│   ├── preflight.py          ← _input_preflight.py generic gate (+ input_issue, INPUT_ISSUE_COLUMNS)
│   ├── progress.py           ← _progress.py
│   ├── result_files.py       ← _result_files.py
│   ├── serial_calculation.py ← _serial_calculation.py
│   ├── worker_runtime.py     ← _worker_runtime.py
│   ├── workflow_kwargs.py    ← _workflow_kwargs.py
│   └── workflow_logging.py   ← _logging.py
│
├── _genome/
│   ├── build_inference.py    ← genome_build_inference.py
│   ├── build_inspection.py   ← _genome_build_inspection.py
│   ├── build_policy.py       ← _genome_build_policy.py
│   ├── hm3.py                ← hm3.py
│   └── liftover.py           ← _kernel/liftover.py
│
├── _annotation/
│   ├── annotate.py           ← annotate_workflow.py
│   ├── bim_grid.py           ← _annotation_bim.py
│   ├── builder.py            ← annotation_builder.py
│   ├── bundle.py             ← _annotation_bundle.py
│   ├── gene_lists.py         ← gene_list_resolver.py
│   ├── gene_query_storage.py ← _gene_query_storage.py
│   ├── identity_cleanup.py   ← _annotation_identity.py
│   ├── in_memory.py          ← _annotation_memory.py
│   ├── interval_projection.py← _annotation_projection.py
│   ├── interval_queries.py   ← _annotation_queries.py
│   ├── loading.py            ← _annotation_loading.py
│   ├── outputs.py            ← outputs.py: AnnotationDirectoryWriter
│   ├── overlap.py            ← _kernel/overlap.py
│   ├── parsing.py            ← _annotation_parsing.py
│   ├── preflight.py          ← _annotation_preflight.py + annotation header checks from _input_preflight.py
│   ├── query_outcomes.py     ← query_annotations.py
│   ├── saved_results.py      ← _annotation_outputs.py
│   ├── semantics.py          ← annotation_semantics.py
│   ├── sources.py            ← _annotation_sources.py
│   ├── storage.py            ← _annotation_storage.py
│   └── tables.py             ← _kernel/annotation.py
│
├── _refpanel/
│   ├── adapters.py           ← _kernel/ref_panel.py (+ identity and sample helpers from _kernel/ldscore.py)
│   ├── build_kernel.py       ← _kernel/ref_panel_builder.py: pair emission, genetic maps, restriction masks
│   ├── builder.py            ← ref_panel_builder.py
│   ├── outputs.py            ← outputs.py: QueryR2DirectoryWriter
│   ├── plink.py              ← _kernel/plink_bed.py, genotype reading only
│   ├── preflight.py          ← inspect_r2_headers + genetic-map header check from _input_preflight.py
│   ├── prepared.py           ← PreparedChromosome and ChromosomeAnnotations from _kernel/ldscore.py;
│   │                           ArrayAnnotations and MappedAnnotations from _kernel/ldscore_projection.py
│   ├── query.py              ← r2_query.py
│   ├── r2_format.py          ← R² readers, SortedR2BlockReader and build_index_remap from _kernel/ldscore.py;
│   │                           writer and constants from _kernel/ref_panel_builder.py; pair lookup from _kernel/r2_query.py
│   └── windows.py            ← both LD-window implementations, under distinct names
│
├── _sumstats/
│   ├── munge_input.py        ← _sumstats_input.py
│   ├── munge_kernel.py       ← _kernel/sumstats_munger.py
│   └── munger.py             ← sumstats_munger.py
│
├── _ldscore/
│   ├── batch_publication.py  ← _ldscore_batch_output.py
│   ├── calculator.py         ← ldscore_calculator.py
│   ├── compute.py            ← LD-score math from _kernel/ldscore.py, PLINK block sums from
│   │                           _kernel/plink_bed.py, ProjectionAccumulator from _kernel/ldscore_projection.py
│   ├── diagnostics.py        ← _ldscore_diagnostics.py
│   ├── direct_annotation.py  ← _direct_annotation.py
│   ├── direct_batches.py     ← _direct_ldscore_batches.py
│   ├── fragments.py          ← _ldscore_fragments.py
│   ├── legacy_converter.py   ← legacy_ldscore_converter.py
│   ├── outputs.py            ← outputs.py: LDScoreDirectoryWriter and its Parquet helpers
│   ├── overlap_table.py      ← overlap_matrix.py: LDScoreOverlap and long-frame read/write
│   ├── preflight.py          ← _ldscore_preflight.py
│   ├── query_variation.py    ← _query_variation.py
│   ├── source.py             ← ldscore_source.py
│   └── gene_index/
│       ├── artifact.py       ← gene_ldscore_index.py: on-disk format and identity hashes (+ inspect_index_paths)
│       ├── atoms.py          ← _kernel/gene_ldscore_index.py
│       ├── batches.py        ← _indexed_ldscore_batches.py
│       ├── build.py          ← gene_ldscore_index.py: index building and publication
│       └── scoring.py        ← gene_ldscore_index.py: index loading and indexed scoring
│
├── _regression/
│   ├── estimators.py         ← _kernel/regression.py
│   ├── h2_scale.py           ← h2_scale.py
│   ├── irwls.py              ← _kernel/_irwls.py
│   ├── jackknife.py          ← _kernel/_jackknife.py
│   ├── outputs.py            ← outputs.py: H2, PartitionedH2, Rg and QuantileH2 writers and column constants
│   ├── overlap_model.py      ← overlap_matrix.py: overlap_aware_category_table, assemble_model_overlap,
│   │                           model_collinearity_error
│   ├── partitioned_parallel.py ← _partitioned_h2_parallel.py
│   ├── prevalence.py         ← prevalence.py
│   ├── quantile_h2.py        ← quantile_h2.py
│   ├── quantile_inputs.py    ← _quantile_inputs.py
│   ├── quantile_storage.py   ← _quantile_storage.py
│   └── runner.py             ← regression_runner.py
│
└── _plotting/
    ├── builders.py           ← plotting/_builders.py
    └── dispatch.py           ← plotting/__init__.py
```

## Changes beyond file moves

Each change keeps behavior identical; the decision that requires it is in parentheses.

1. **Reference-panel and LD-score boundary (Q1, option A2).**
   `_refpanel` owns reading, the preparation contract and the R² format; `_ldscore` owns all LD-score math.
   The preparation code of `_kernel/ldscore.py` moves to `_refpanel/` and its LD-score math to `_ldscore/compute.py`.
2. **PLINK LD sums leave the reader class (Q1, A2).**
   `ldScoreVarBlocks`, `__corSumVarBlocks__`, `__l2_unbiased__` and `_zero_pairs_outside_window` (about 110 lines in [`_kernel/plink_bed.py`](../../src/ldsc/_kernel/plink_bed.py)) become a function in `_ldscore/compute.py` that takes the reader's `nextSNPs`, `m` and `n`.
   The loop body is unchanged.
   Legacy LDSC defines the same algorithm as the method `ldScoreVarBlocks` of `__GenotypeArrayInMemory__` in `ldscore/ldscore.py:121` of the legacy repository; legacy comparisons cite the new function next to that method.
   Callers: `compute_chromosome()`, two calls in `gene_ldscore_index.py`, and tests including two `monkeypatch.setattr(geno_type, "ldScoreVarBlocks", ...)` calls.
3. **Projection module split (1f).**
   `ArrayAnnotations` and `MappedAnnotations`, the accessor types of `PreparedChromosome.annotations`, move to `_refpanel/prepared.py`.
   `ProjectionAccumulator` moves to `_ldscore/compute.py`; it reads accessors only through `shape`, `columns` and `read()`.
4. **One R² format module (1c).**
   `_refpanel/r2_format.py` holds the schema constants (`R2_ENCODING`, `R2_QUANT_SCALE`, `SIDECAR_FLOAT_FORMAT`), the contract readers and checks, sidecar binding, `SortedR2BlockReader`, `build_index_remap`, the index pair lookup and the writer.
5. **Renames.**
   The kernel `AnnotationBundle` (`_kernel/ldscore.py`) becomes `ChromosomeAnnotations`, because the public `ldsc.AnnotationBundle` keeps its name (1g).
   The two `build_window_coordinates` functions get distinct names in `_refpanel/windows.py`; the proposed names are `ldscore_window_coordinates()` for the `ldscore` version and `r2_panel_window_coordinates()` for the `build-r2-panel` version (1e).
   Merging them is a separate, behavior-changing follow-up.
6. **Input preflight split (Q2).**
   `_runtime/preflight.py` keeps `InputGate`, the issue records and the generic declaration, PLINK-companion and artifact-path checks.
   The annotation header checks move to `_annotation/preflight.py` and reach `inspect_declared_inputs()` through its `checks=` hook or a header callback.
   The genetic-map header check and `inspect_r2_headers()` move to `_refpanel/preflight.py`, and `inspect_index_paths()` to `_ldscore/gene_index/artifact.py`.
7. **`overlap_matrix.py` split (Q3).**
   `overlap_aware_category_table()`, `assemble_model_overlap()` and `model_collinearity_error()` are used only by regression and quantile code, and the first is the file's only user of the regression estimators and of `scipy.stats`.
   They move to `_regression/overlap_model.py`; `LDScoreOverlap` and its long-frame read/write go to `_ldscore/overlap_table.py`.
   `_kernel/overlap.py` goes to `_annotation/`, because `query_annotations.py` builds `OverlapContribution` objects and the overlap matrix is a product of annotations over reference SNPs without LD.
8. **Writer split (Q5).**
   Each writer moves to the domain that owns its artifact; the helpers used by more than one domain go to `_runtime/artifacts.py`.
   `_cast_parquet_floats`, `_write_chromosome_aligned_parquet`, `DEFAULT_COUNT_CONFIG`, `QUERY_STATUS_COLUMNS` and `REGRESSION_LD_SCORE_COLUMN` go to `_ldscore/outputs.py`; the result-table column constants go to `_regression/outputs.py`.
9. **`config` above `_runtime` (Q4).**
   `run_from_kwargs()` receives the default log level from its caller instead of calling `get_global_config()` ([`_workflow_kwargs.py:102`](../../src/ldsc/_workflow_kwargs.py)).
10. **Small relocations (Q5).**
    `VALID_SNPS` moves from `_kernel/regression.py` to `_core/snp_identity.py`; `input_issue` and `INPUT_ISSUE_COLUMNS` from `_annotation_preflight.py` to `_runtime/preflight.py`; `sort_frame_by_genomic_position` to `_core/row_alignment.py`, because placing it in `_core/chromosomes.py` would create the import cycle `chromosomes -> columns -> snp_identity -> coordinates -> chromosomes`; `_LDSCORE_INTERSECTION_DOC` to `_refpanel`; `_direct_annotation.py` to `_ldscore`.
11. **Gene LD-score index split (Q6).**
    `gene_ldscore_index.py` (3,298 lines) splits into index building, the on-disk index format and indexed scoring; no other workflow module is split.
12. **Dead code (Q8, 1d, D1).**
    A first commit removes unreferenced definitions and unused imports, repeating the scan until it finds nothing, before any structural change.

## Dead-code policy D1

Five functions have no production caller and are referenced only by tests: `parse_annotation_file()` and `combine_annotation_groups()` with their private helpers and aliases, `block_left_to_right()`, and `_accumulate_pair_contributions()`, all in `_kernel/ldscore.py`.
They are deleted.
Their 15 tests are handled one by one:

- a test that checks only the deleted function, or a rule already tested on the production path, is deleted;
- a test whose rule has no production-path test is ported to the production code (`normalize_annotation_chunk()`, the annotation source loader or `ProjectionAccumulator`);
- `test_annotate_output_round_trips_through_ldscore_parser` is ported so that it reads `ldsc annotate` output through the production loader.

Rules at risk of losing their only test include short rows that must not shift columns and allele-free annotations that must match on the base identity.
A table of each test's disposition is shown to the maintainer before the commit.

## Verified constraints

- `setup.py` uses `find_packages(where="src")`, so new subpackages are installed without configuration changes; the `package_data` path `ldsc/data` is unchanged.
- The lazy loading of SciPy-backed workflows (`f321e0e`) is guarded by `test_top_level_help_does_not_import_heavy_workflows` ([`tests/test_package_layout.py:179`](../../tests/test_package_layout.py)), which names old module paths.
  Each move commit updates those names, and the updated test is checked once to fail when the rule is broken; otherwise it would pass on module names that no longer exist.
- Logger names are explicit strings such as `"LDSC"` and `"LDSC.preflight"`, so log content does not depend on module paths.
- No artifact records a module path; pickles carry data to worker processes within one run and to scratch files only.
- 89 of the 102 test files and the maintenance scripts under `tools/` import or patch modules that move; tests, string patch targets and tools are updated in the commit that moves each module.
- [status.md](../../status.md) holds a separate annotation-grid design handoff whose implementation seams cite current module paths; the two efforts must not run concurrently, and the handoff's paths are updated when this reorganization lands.

## Verification

- Every commit: `pytest`, then `python -m unittest discover -s tests -p 'test*.py' -v`, run one after the other; `ldsc --help` and `python -m ldsc --help`; the layering test.
- The commit that lifts the PLINK LD sums: the existing known-value tests pass unchanged, and LD-score outputs computed on the repository's test fixtures before and after the commit are bitwise identical.
- The commit that introduces the layering test: a deliberate upward import makes the test fail.
- The final commit: no `_kernel/` directory and no module at the package root other than the shell.

## Documentation

- [AGENTS.md](../../AGENTS.md) (with its `CLAUDE.md` link): module paths only, no rule changes.
- [code-structure.md](../dev/code-structure.md) and [architecture.md](../dev/architecture.md): rewritten for the new layout, receiving the contract from this specification.
- Other developer, reference and user-guide pages: module paths updated in the commit that moves each module.
  User-guide pages include the "Developer reference" lines of [troubleshooting.md](../user-guide/troubleshooting.md), [workflows/munge-sumstats.md](../user-guide/workflows/munge-sumstats.md), [concepts/snp-identity-and-genome-build.md](../user-guide/concepts/snp-identity-and-genome-build.md) and the private-module example in [python-usage.md](../user-guide/python-usage.md).
- The stale citations of `_kernel.ldscore.compute_chrom_from_parquet()` and `compute_chrom_from_plink()` in [troubleshooting.md](../user-guide/troubleshooting.md) are corrected when that line is edited.
- New contributor terms, such as the tier names, are recorded under the proposed terms of [DOCS_VOCABULARY.md](../dev/DOCS_VOCABULARY.md) when the developer pages are rewritten.

## Decision record

| ID | Decision |
| --- | --- |
| Layout | Private domain subpackages with per-domain writers and one-directional layering enforced by a test. |
| Q1 | Option A2: the preparation side owns reading, the contract and the R² format; all LD-score math lives in `_ldscore`. |
| 1b | Resolved by A2: `ProjectionAccumulator` belongs to `_ldscore`. |
| 1c | One R² format module for reading and writing. |
| 1d | Delete the production-unused annotation-parsing code; refined by D1. |
| 1e | Distinct names for the two LD-window functions; merge later. |
| 1f | Split the projection module between `_refpanel` (accessors) and `_ldscore` (accumulator). |
| 1g | Rename the kernel `AnnotationBundle` to `ChromosomeAnnotations`. |
| Q2 | Split `_input_preflight.py` between `_runtime` and the domains. |
| Q3 | Split `overlap_matrix.py` minimally; `_kernel/overlap.py` goes to `_annotation`. |
| Q4 | `config` sits above `_runtime`; `run_from_kwargs()` receives the default log level. |
| Q5 | Placement corrections and the writer split as listed above. |
| Q6 | Split only `gene_ldscore_index.py`. |
| Q7 | The root keeps only the shell; `plotting/` becomes `_plotting/`. |
| Q8 | The test checks every import with no allowlist; commits run per domain, bottom-up, after a dead-code commit. |
| Q9 | `tests/` stays flat. |
| Q10 | Update `AGENTS.md` paths and all affected documentation; fix the stale troubleshooting citations. |
| R1 | Domain `__init__.py` files hold only a docstring. |
| R2 | No leading underscore on modules inside private packages. |
| R3 | The module map above. |
| R4 | Record the design in this specification and its plan, with a pointer in `status.md`. |
| D1 | Delete the five production-unused functions; port each of their 15 tests whose rule has no production-path test, delete the rest, and show the per-test table before the commit. |

## Out of scope

- Merging the two LD-window implementations, which differ in `--ld-wind-snps` behavior and error messages.
- Unifying `_result_files.atomic_write_json()` with the writers' `_atomic_write_json()`; [code-structure.md](../dev/code-structure.md) records that the writers keep their separate JSON serialization.
- The migration code for the older gene-index build-state location (`_migrate_legacy_gene_index_build_state`).
- Per-domain test folders.
- Splitting the remaining large workflow modules or moving CLI parsers out of them.
