# Plotting and Post-processing Integration Plan

Last updated on: 2026-09-07

## Status

Status: complete.

Reference specification: [`docs/specs/2026-09-07-plotting-module-integration-design.md`](../specs/2026-09-07-plotting-module-integration-design.md).

This is a living execution guide. Update slice status, validation evidence, and discovered constraints as implementation proceeds. Evidence that changes a public, scientific, or lifecycle contract returns to design review rather than being absorbed silently.

## Goal and success signal

Integrate the approved detached plotting prototype into LDSC3 as an optional Matplotlib-only post-processing layer, add exact saved data for the binned unpartitioned h2 diagnostic, expose post-fit observed-to-liability conversion, and add durable failed-overwrite markers without changing existing workflow action order.

Completion is observable when a current canonical result can be passed to `ldsc plot --result-dir ...` and produces the single approved figure for its analysis regime; `ldsc convert-h2-scale` supports exact and prevalence-range conversion without refitting; every new h2 result carries the exact bin table required by its diagnostic plot; the core package remains installable and runnable without Matplotlib; and failed authorized overwrites leave the specified marker while successful retries remove it.

## Context and constraints

- Work in `ldsc_py3_restructured` on branch `restructure`.
- Preserve the unrelated existing edit in `docs/wiki/continuous-annotation-partitioned-ldsc.md`.
- The detached `ldsc_plotting_sandbox` supplies approved visual behavior and focused prototype tests, but its standalone package, separate commands, output-directory flags, duplicated liability formula, and current metadata names must not be copied blindly.
- `src/ldsc/cli.py` owns parsing and lazy dispatch only. Scientific calculations and result validation remain in workflow or kernel modules.
- `src/ldsc/outputs.py` and `src/ldsc/path_resolution.py` own canonical artifact writing, collision preflight, and stale-owned output cleanup. Reuse those abstractions and preserve produced-path versus owned-path distinctions.
- `src/ldsc/_logging.py` already retains a failed workflow log and traceback. Failure markers supplement that log; they do not replace it or add rollback.
- The existing `--samp-prev`/`--pop-prev` behavior and vectorized `liability_conversion_factor` are implemented and numerically authoritative. Migration verifies and reuses them rather than rewriting the estimator.
- Matplotlib is optional and imported only on an explicit plot-producing path. Seaborn is not added.
- Plotting is never automatic. No plotting failure may alter a source result or prevent a core numerical workflow from completing.
- Apply focused TDD in each slice: write the smallest failing behavioral or numerical test, implement, then refactor while green.

## Interfaces and invariants

### Public seams

- `ldsc plot --result-dir RESULT_DIR [--overwrite] [--log-level LEVEL]` dispatches one plot from canonical metadata and has no `--output-dir`.
- `ldsc convert-h2-scale --h2-result-dir H2_RESULT --samp-prev P (--pop-prev K | --pop-prev-range MIN MAX) [--num-points N] [--overwrite] [--log-level LEVEL]` writes to the fixed nested conversion root and has no `--output-dir`.
- `plot_result(result_dir, output_dir=None, overwrite=False)` returns an immutable `PlotArtifact` with kind, path, live figure, and axes.
- `convert_h2_scale(..., output_dir=None, overwrite=False)` returns an immutable conversion artifact with table, metadata, log, and optional plot paths.
- Public Python exports use the package's lazy-export mechanism so importing `ldsc` does not import Matplotlib.

### Scientific invariants

- Every plotted estimate, SE, p-value, fitted expectation, and regression weight comes from a canonical saved artifact or the exact h2 fit that creates that artifact. Plotting never refits.
- The h2 diagnostic bins exactly the final-slope SNP population after the active chi-square filter. Under two-step estimation, this is the second-step slope population.
- The h2 bin table has the fixed columns and at-most-50 deterministic rank-bin behavior specified in the reference specification.
- Liability conversion always starts from observed-scale h2 and scales the block-jackknife SE by the same factor while treating prevalence as fixed.
- Functional and quantile enrichment use one SE; anchor rg uses one SE; the query suite alone displays nominal one-sided p-values; no multiple-testing correction is introduced.

### Artifact and lifecycle invariants

- Plot input dispatch requires `diagnostics/metadata.json` and the plotting-relevant `files` entry. It does not inspect `schema_version` or infer from filenames.
- Older h2 results without `diagnostics/ld_score_regression_bins.tsv` fail with rerun guidance; there is no compatibility layer.
- Default plot and conversion roots are package-owned. An explicit Python destination is unmanaged.
- Existing workflow publication, stale-cleanup, and failure order is unchanged. There is no new quarantine, rollback, or whole-family transaction.
- A failed authorized overwrite writes the applicable `RUN_FAILED` marker while leaving the artifacts exactly as the command's existing order left them. A successful materializing retry removes its marker after normal success and cleanup.

## Implementation slices

### Slice 1: Persist exact h2 regression-bin diagnostics

Status: complete.

Goal: make every new unpartitioned h2 result independently plot-ready without adding a plotting dependency.

Likely areas: `src/ldsc/regression_runner.py`, a small pure diagnostic helper in the regression workflow or private kernel, `src/ldsc/outputs.py`, `tests/test_kernel_regression.py`, `tests/test_regression_workflow.py`, and `tests/test_output.py`.

Work:

- Factor the current h2 input preparation so one internal fit path returns the `Hsq` estimate plus the exact post-filter arrays needed for diagnostics while preserving `RegressionRunner.estimate_h2()` and its existing `Hsq` return contract.
- Build stable ascending rank bins from the fitted unpartitioned LD Score, preserving fitted row order for ties and using `min(50, n_snps)` nonempty, nearly equal-count bins.
- Evaluate per-SNP fitted chi-square from the final slope/intercept and exact per-SNP sample size. Evaluate the existing `Hsq.weights` function with the final model and exact LD Score, regression-weight LD Score, sample size, and reference-SNP count, then aggregate the required table fields.
- Extend `H2DirectoryWriter` to validate and atomically write `diagnostics/ld_score_regression_bins.tsv`, record `files.ld_score_regression_bins`, and include the table in both produced and owned h2 artifact lists.
- Carry the bin table through the h2 workflow without changing `h2.tsv`, regression estimates, exact-K reporting, or Matplotlib availability.

Validation:

- Hand-computed arrays verify bin membership, stable tie behavior, counts, LD Score bounds and means, chi-square means and sample SDs, sample-size means, fitted expectations, and weights.
- Cover free and constrained intercepts, default two-step estimation, explicit chi-square filtering, variable sample sizes, fewer than 50 SNPs, and ordinary 50-bin output.
- End-to-end h2 tests assert `sum(n_snps)` equals metadata `n_snps`, metadata points to the table, overwrite replaces the table, and public `estimate_h2()` still returns the same estimates and type.
- Run the focused regression-kernel, regression-workflow, and h2-writer tests.

Checkpoint: independently compare the saved two-step fitted expectations and final-model weights with the estimator formula before any plot consumes them. Revise the plan if exact fit-input reuse would require retaining a large SNP table beyond output construction.

### Slice 2: Add post-fit liability-scale conversion

Status: complete.

Goal: provide exact-K and prevalence-sensitivity conversion as a canonical derived result family while reusing the existing numerical primitive.

Likely areas: a new focused workflow module such as `src/ldsc/h2_scale.py`, `src/ldsc/cli.py`, lazy exports in `src/ldsc/__init__.py`, package dependency metadata, and focused conversion/CLI tests adapted from the sandbox.

Work:

- Add strict loading for a one-row canonical h2 result using `artifact_type`, `files.summary`, `total_h2_obs`, and `total_h2_obs_se`; never select populated liability fields.
- Implement exact and inclusive linear-range modes through `_kernel.regression.liability_conversion_factor`, with scalar validation aligned to existing `--samp-prev`/`--pop-prev` semantics.
- Write the fixed conversion table, derived metadata with `artifact_type="h2_scale_conversion_result"`, workflow log, and conditional sensitivity figure below `postprocessing/liability-scale/` by default; support only the approved Python `output_dir=` override.
- Add the lazy CLI command and public Python export. Preserve the existing h2, partitioned-h2, and rg exact-K shortcuts and verify the spelling `--samp-prev` throughout.
- Add the `plot` optional dependency extra with Matplotlib as its only member. Keep exact mode free of Matplotlib imports; range mode checks the optional dependency before writing and uses the approved sensitivity style.
- Apply the fixed-family overwrite rules, including exact-mode removal of a stale range plot after successful publication.

Validation:

- Compare exact and vectorized range results to direct calls of the kernel primitive, including the existing numerical anchor and scalar/array parity.
- Test mutual exclusion, bounds, one-row source enforcement, missing/invalid observed fields, inclusive endpoints, default/custom grid sizes, provenance, default and Python-override destinations, collisions, and stale sensitivity removal.
- In a no-Matplotlib environment or import-isolation test, verify exact mode succeeds while range mode fails before partial output with `ldsc[plot]` guidance.
- Exercise both `ldsc convert-h2-scale --help` and `python -m ldsc convert-h2-scale --help` without loading unrelated heavy workflows.

Checkpoint: confirm dependency checking happens before range artifacts are written but after enough source validation to produce an accurate user error; do not duplicate the conversion formula from the sandbox.

### Slice 3: Deliver plotting infrastructure, h2 diagnostic, and rg plots

Status: complete.

Goal: establish the optional dispatcher and artifact lifecycle end to end with the structurally distinct h2 and rg plots.

Likely areas: a small `src/ldsc/plotting` implementation area or equivalently focused plotting modules, `src/ldsc/cli.py`, lazy exports in `src/ldsc/__init__.py`, package metadata, and new plotting/CLI tests adapted from the sandbox.

Work:

- Add strict canonical metadata/table loading and dispatch for `h2_result` and `rg_result`; require the specified metadata fields and `files` paths and reject unsupported, contradictory, loose, per-pair, or incomplete inputs.
- Add the immutable `PlotArtifact`, shared Matplotlib style helpers, noninteractive backend selection, fixed 300-dpi saving, plot-family metadata/log writing, fixed default destination, and advanced Python output override.
- Implement the binned h2 plot directly from saved means, fitted expectations, and weights. Normalize only the displayed relative-weight color scale and do not fit a trend line.
- Migrate the all-pairs lower-triangular rg heatmap and anchor-trait forest plot from the prototype using current `rg.tsv` and metadata contracts. Preserve failed cells, estimate/SE labels, heatmap colorbar omission, trait order, and anchor validation.
- Add one lazy `plot` CLI dispatcher. The CLI closes its figure after saving; Python callers receive the live figure and axes.

Validation:

- Dispatch tests cover both rg regimes and h2, exact fixed filenames, metadata/log contents, returned objects, default and override roots, collisions, unsupported artifact types, missing files/columns, duplicate rg pairs, and old h2 results without bins.
- Numerical artist tests assert h2 point coordinates, fitted path, normalized weight ordering, rg cell/forest coordinates, one-SE intervals, lower-triangle masking, failed labels, and absence of a heatmap color bar.
- Verify `import ldsc`, top-level help, and non-plot commands do not import Matplotlib; verify `ldsc plot` reports the optional-extra installation command when unavailable.
- Run focused plotting, CLI, output, and import-isolation tests under a headless backend.

Checkpoint: inspect one h2, all-pairs rg, and anchor rg example image against the approved sandbox outputs before expanding dispatch to additional result types.

### Slice 4: Migrate partitioned-h2 and quantile-h2 plots

Status: complete.

Goal: complete metadata-driven plotting for functional, cell-type/query, and continuous-annotation result suites without exposing scientifically discouraged plots.

Likely areas: the plotting implementation from slice 3, `tests/test_plotting.py` or focused plot modules, existing partitioned/quantile fixtures, and selected sandbox tests.

Work:

- Add functional-category dispatch from `analysis_type=functional_category` and `headline_metric=enrichment`; read the declared root summary and draw dark-grey horizontal enrichment bars with one SE and a medium-grey dashed null behind the bars.
- Add cell-type dispatch from the aggregate root only. Plot one query per row as nominal one-sided `-log10(P)` from separate baseline-conditional fits; preserve unavailable values and recover underflowed zero p-values from finite coefficient z values rather than clipping scientific evidence arbitrarily.
- Add quantile dispatch using `artifact_type=quantile_h2_result`, `target_annotation`, and the declared quantile table. Preserve ascending quantile order, realized bounds, darker colors for larger quantiles, one SE, and the same enrichment null behind bars.
- Keep raw functional coefficients/p-values, individual baseline-plus-query results, query enrichments, corrected significance, and thresholds inaccessible through the public dispatcher.
- Apply shared figure sizing, concise labels, estimate/SE annotations, and explicit unavailable-result handling without reintroducing Seaborn.

Validation:

- Assert metadata regime selection, aggregate-root enforcement, required source columns, query ordering, `-log10(P)` values including p-value underflow recovery, quantile order/bounds, bar orientation/colors, error magnitudes, null-line color/style/z-order, and explicit unavailable labels.
- Run the complete optional plotting test group and regenerate all seven example figures for visual review.
- Run the relevant partitioned-h2 and quantile-h2 source-workflow tests to confirm plotting additions do not alter numerical outputs.

Checkpoint: compare the functional and quantile null-line z-order directly in rendered figures and confirm the cell-type title states that queries come from separate baseline-conditional fits.

### Slice 5: Integrate derived-output cleanup and failed-overwrite markers

Status: complete.

Goal: apply the approved nested-output ownership and durable failure visibility across all materializing workflows without changing their existing action order.

Likely areas: a small shared helper near `src/ldsc/path_resolution.py` and `src/ldsc/_logging.py`, each public workflow boundary, existing writer/preflight code in `src/ldsc/outputs.py`, and lifecycle-focused tests.

Work:

- Expand first with a narrowly scoped marker helper that writes the approved plain-text content on an escaping failure from an authorized overwrite and removes the applicable marker after successful completion. The helper must not move, restore, delete, or classify scientific artifacts.
- Apply the helper first to the new plot and conversion workflows and h2 as a pilot. Verify failures before log creation and after log opening, then migrate annotate, ldscore, build-ref-panel, build-gene-ldscore-index, convert-ldsc2-ldscores, munge-sumstats, partitioned-h2, quantile-h2, rg, and query-r2.
- Preserve command-specific output scopes. Implement chromosome-scoped `build-ref-panel` marker names without allowing one shard to remove another shard's marker.
- Add default `plots/` to the owned stale set of h2, partitioned-h2, rg, and quantile-h2, and default `postprocessing/` to h2. Ensure every early workflow preflight and final writer agrees on produced versus owned paths; remove these derived roots only through existing post-success cleanup.
- Treat orphaned reserved derived roots as no-overwrite collisions. Keep `RUN_FAILED` outside scientific metadata and collision families so it cannot prevent a retry.
- Contract after migration with one inventory-level test or explicit workflow matrix proving every public materializing CLI/high-level Python path uses the policy; exclude in-memory APIs, `munge-sumstats --infer-only`, and no-overwrite preflight failures.

Validation:

- Inject failures before a workflow log opens, during output publication, and before stale cleanup. Assert the exact files left by existing order are untouched by marker handling and the marker truthfully warns about incomplete or mixed results.
- Verify the canonical failed log retains its `Failed` footer and traceback when opened, the marker points to it, and a successful retry removes the marker only after ordinary success and stale cleanup.
- Test plot/conversion marker roots, Python overrides, ordinary roots, chromosome and full-suite build-ref-panel scopes, and no marker for no-overwrite collisions.
- Re-run existing staging/backup tests unchanged to prove the marker layer did not replace or generalize their command-specific behavior.
- Run logging, path-resolution, writer, and representative workflow suites after the pilot and again after the full migration.

Checkpoint: use the repository lesson on produced/owned path parity as a release gate. A newly produced path must never also appear in the stale set, and marker integration must introduce no `try` block that reorders publication or cleanup.

### Slice 6: Documentation, examples, and release validation

Status: complete.

Goal: make the implemented behavior discoverable to developers and scientifically interpretable to users, then validate the optional dependency and full repository boundary.

Likely areas: `docs/current/plotting-module.md`, `tutorials/plotting-results.md`, `README.md`, `AGENTS.md`, `docs/current/architecture.md`, `docs/current/code-structure.md`, `docs/current/layer-structure.md`, `docs/current/data-flow.md`, `docs/current/path-specification.md`, `docs/current/io-argument-inventory.md`, `docs/current/artifact-metadata-field-inventory.md`, `docs/current/workflow-logging.md`, relevant LDSC tutorials, and `docs/troubleshooting.md` where multi-cause abort guidance is needed.

Work:

- Write the dedicated developer module/seam document and scientist-facing plotting manual required by the specification; explain compatible result roots, each scientific question, uncertainty semantics, and interpretation limits.
- Update package navigation, command/output inventories, metadata fields, dependency installation, nested ownership, failure markers, and the two CLI exceptions to the former universal `--output-dir` rule.
- Link the user manual from the README and relevant h2, partitioned-h2, quantile-h2, and rg tutorials. Keep the existing liability tutorial consistent with exact shortcuts and the new post-fit command.
- Add or adapt small deterministic examples that generate all seven approved plots without bundling large study data. Keep the detached sandbox as development history rather than a runtime dependency.
- Use scientific Python documentation conventions for the new public functions, artifacts, and modules.

Validation:

- Run every documented CLI example against actual `--help` output and inspect one successful artifact tree for each supported source regime.
- Verify generated PNG dimensions/DPI and visually review all seven example plots.
- Run the core test suite in an environment without Matplotlib, then the full optional plotting tests with `ldsc[plot]` installed.
- Run `pytest`, `python -m unittest discover -s tests -p 'test*.py' -v`, `ldsc --help`, `python -m ldsc --help`, both new command help paths, focused end-to-end h2/plot/conversion workflows, `git diff --check`, and `git status --short`.

Checkpoint: do not mark completion until a core-only install and test pass proves that plotting dependencies and imports are isolated, and until every metadata-listed file in each new artifact family exists after overwrite.

## Completion evidence

- The core-only pytest run completed with 1,319 passed and 4 expected skips; Matplotlib was absent, while package import, exact conversion, and all core workflows remained available.
- The plotting-enabled focused matrix completed with 37 passed across plotting, h2 conversion, h2 diagnostics, derived-output lifecycle, and failure markers.
- Standard-library discovery ran 1,015 tests successfully with 2 expected skips after the optional plotting module was made compatible with both test runners.
- Installed and module CLI help succeeded for the full command surface, `plot`, and `convert-h2-scale`; the core-only installed `ldsc plot` command failed cleanly with `ldsc[plot]` installation guidance.
- `tools/generate_plot_examples.py` produced all seven approved figures. Each PNG recorded approximately 300 dpi, and all 13 generated metadata documents pointed only to files that existed.
- Python compilation, documentation/contract review, `git diff --check`, and the final repository status review completed successfully; the pre-existing wiki edit remained untouched.

## Risks and checkpoints

- The highest scientific risk is using a subtly different SNP population or intermediate IRWLS state for the h2 diagnostic. Resolve this in slice 1 with explicit fit-input reuse and independent formula checks.
- The broadest operational risk is the package-wide failure-marker migration. Keep the shared helper limited to marker I/O, migrate workflow boundaries incrementally, and retain every command's existing publication and cleanup behavior.
- Stale-output cleanup can delete a freshly written artifact when produced and owned lists diverge. Test complete write-then-overwrite paths and metadata-listed file existence for every modified writer.
- Lazy loading can be defeated by an import in `ldsc.__init__`, `cli.py`, a type annotation, or conversion exact mode. Treat import-isolation tests as a milestone gate rather than relying on code inspection.
- Very large trait/category suites can create unreadable figures. Preserve data completeness and deterministic dynamic sizing; do not introduce hidden truncation or a new public customization surface without design review.
- Matplotlib version differences can shift layout. Test semantic artists and saved metadata rather than pixel identity, while retaining human review of representative outputs.
- Revise the plan if implementation evidence changes a plotted estimand, uncertainty definition, artifact ownership boundary, or public command. Ordinary module placement and private helper names do not require design review.

## Out of scope

- Automatic plots from h2, partitioned-h2, rg, or quantile-h2.
- Seaborn, interactive plotting, publication themes, exhaustive customization, or a separate plotting package.
- A CLI output-directory override or public one-command-per-plot surface.
- Manhattan plots, legacy LDSC2 plotting, or compatibility for older h2 results without bins.
- The legacy textual LD-score construction summary as a plotting source.
- Functional coefficient/p-value plots, individual baseline-plus-query plots, query enrichment plots, or multiple-testing correction.
- Changes to LDSC estimators, regression defaults, exact-K shortcut semantics, or rg scale invariance.
- Quarantine, rollback, whole-family transactions, or reordering existing workflow actions.
