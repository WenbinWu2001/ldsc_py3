# Implementation plan: domain subpackage reorganization

Last updated on: 2026-10-27

Status: not started; implementation is not yet authorized.
This plan implements the [domain subpackage specification](../specs/2026-10-27-domain-subpackages.md), which holds the decisions, the module map and the import rules.
Delete this plan together with the specification once the final commit lands.

## Ground rules

- Work in `ldsc_py3_restructured` on branch `restructure`, starting from a clean worktree.
- Run no other source change on this branch while the plan is in progress; the annotation-grid design in [status.md](../../status.md) waits until the final commit.
- Each commit leaves the package importable and every test passing.
- Move files with `git mv` so that history follows them; a file that splits keeps its history in the part that keeps most of its lines.
- A module's imports, tests, string patch targets, `tools/` scripts and documentation paths are updated in the commit that moves the module.
- Commit messages follow `<type>(<scope>): <description>` with a subject of at most 50 characters and body lines of at most 72; they carry no AI tool names or attribution lines.
- Before each commit, run the checks below and report their results; after ten or more messages or a gap, reread the specification, this plan, `status.md` and `lessons.md`.

Checks for every commit:

```bash
pytest
python -m unittest discover -s tests -p 'test*.py' -v
ldsc --help
python -m ldsc --help
```

Run `pytest` and `unittest` one after the other, never in parallel.
From commit 3 on, `pytest` includes the layering test.

## Commits

### 0. Baseline

Record the test counts and results of both runners at the starting commit, so later commits can be compared against them.
No commit.

### 1. `refactor: remove unused definitions`

- Remove unreferenced definitions and unused imports across `src/ldsc/`, repeating the scan until it finds nothing new.
- Apply the dead-code policy D1 after the maintainer approves the per-test disposition table.

### 2. `refactor(core): move shared helpers down`

- Move `VALID_SNPS` from `_kernel/regression.py` to `_kernel/snp_identity.py`.
- Move `sort_frame_by_genomic_position()` from `_kernel/ldscore.py` to `chromosome_inference.py`.

### 3. `refactor(core): move primitives into _core`

- Move the eight `_core` modules.
- Add `tests/test_layering.py`, which derives tiers from package paths and checks every import inside the package.
- Confirm once that a deliberate upward import makes the new test fail.

### 4. `refactor(runtime): split input preflight`

- Keep the generic gate in `_input_preflight.py` and move `input_issue` and `INPUT_ISSUE_COLUMNS` into it.
- Move the annotation header checks to `_annotation_preflight.py`.
- Create the root module `_ref_panel_preflight.py` for `inspect_r2_headers()` and the genetic-map header check; it moves into `_refpanel/` in commit 11.
- Move `inspect_index_paths()` into `gene_ldscore_index.py`.

### 5. `refactor(runtime): pass default log level in`

- `run_from_kwargs()` receives the default log level from each of its six callers.

### 6. `refactor(runtime): move infrastructure into _runtime`

- Move the eleven `_runtime` modules, including the new `artifacts.py` that receives the shared writer base from `outputs.py`.
- Update the names in `test_top_level_help_does_not_import_heavy_workflows` and confirm once that the updated test fails when a heavy module is imported eagerly.

### 7. `refactor(genome): move build and liftover code`

- Move the five `_genome` modules.

### 8. `refactor(annotation): move annotation domain`

- Move the annotation modules, `_kernel/overlap.py` and `_kernel/annotation.py`.
- Move `AnnotationDirectoryWriter` from `outputs.py` to `_annotation/outputs.py`.
- `_direct_annotation.py` stays at the root until commit 14, because it belongs to `_ldscore`.

### 9. `refactor(ldscore): lift PLINK LD sums from reader`

- Turn `ldScoreVarBlocks()` and its helpers in `_kernel/plink_bed.py` into a function in `_kernel/ldscore.py` that takes the reader's `nextSNPs`, `m` and `n`.
- Update `compute_chromosome()`, the two calls in `gene_ldscore_index.py` and the tests.
- Compute LD scores on the repository's test fixtures before and after the change and require bitwise-identical outputs, in addition to the existing known-value tests.

### 10. `refactor(refpanel): extract panel contract and format`

- Create `_refpanel/prepared.py`, `_refpanel/windows.py` and `_refpanel/r2_format.py` from the preparation code of `_kernel/ldscore.py`, the accessors of `_kernel/ldscore_projection.py`, the format code of `_kernel/ref_panel_builder.py` and all of `_kernel/r2_query.py`.
- Rename the kernel `AnnotationBundle` to `ChromosomeAnnotations` and the two window functions as the specification proposes.

### 11. `refactor(refpanel): move reference-panel domain`

- Move `ref_panel_builder.py`, `r2_query.py`, `_kernel/ref_panel.py`, `_kernel/ref_panel_builder.py`, `_kernel/plink_bed.py` and `_ref_panel_preflight.py`.
- Move `QueryR2DirectoryWriter` to `_refpanel/outputs.py`.

### 12. `refactor(sumstats): move munging domain`

- Move the three `_sumstats` modules.

### 13. `refactor(ldscore): split index and overlap modules`

- Split `gene_ldscore_index.py` at the root into index building, the on-disk format and indexed scoring.
- Split `overlap_matrix.py` into the LD-score container and the regression-only functions.

### 14. `refactor(ldscore): move LD-score domain`

- Move the LD-score modules, the gene LD-score index modules, `_direct_annotation.py` and the container part of `overlap_matrix.py`.
- Merge the LD-score math left in `_kernel/ldscore.py` and `ProjectionAccumulator` into `_ldscore/compute.py`.
- Move `LDScoreDirectoryWriter` and its Parquet helpers to `_ldscore/outputs.py`.

### 15. `refactor(regression): move regression domain`

- Move the regression modules and the regression-only part of `overlap_matrix.py`.
- Move what remains of `outputs.py` to `_regression/outputs.py` and remove the empty `_kernel/` directory.

### 16. `refactor(plotting): make plotting private`

- Move `plotting/__init__.py` to `_plotting/dispatch.py` and `plotting/_builders.py` to `_plotting/builders.py`.
- Extend the layering test to fail if any module other than the shell sits at the package root.

### 17. `docs: describe the domain layout`

- Update the module paths in `AGENTS.md`.
- Rewrite [code-structure.md](../dev/code-structure.md) and [architecture.md](../dev/architecture.md) with the `architecture-doc` skill, moving the durable contract out of the specification.
- Record new contributor terms under the proposed terms of [DOCS_VOCABULARY.md](../dev/DOCS_VOCABULARY.md).
- Correct the stale `compute_chrom_from_parquet()` and `compute_chrom_from_plink()` citations in [troubleshooting.md](../user-guide/troubleshooting.md).
- Update the module paths in the annotation-grid handoff in `status.md` and remove this plan's pointer.
- Delete the specification and this plan.

## Risks

- Commits 9 to 11 change the most code; review them one at a time and stop on any numerical difference.
- A string patch target or module name that a test no longer resolves can make a test pass without testing anything; after each move, search tests for the old dotted paths and require zero matches.
- Test files that patch `ldsc.ref_panel_builder.kernel_ldscore.PlinkBEDFile` and similar paths depend on module aliases that change; update them to the new owner of each name.
