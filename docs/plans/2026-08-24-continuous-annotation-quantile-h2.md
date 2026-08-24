# Continuous Annotation and Quantile Heritability Implementation Plan

Last updated on: 2026-08-24

## Status

Status: complete.

Reference specification: `docs/specs/2026-08-24-continuous-annotation-quantile-h2-design.md`.

This is a living execution guide. Update slice status, validation evidence, and discovered constraints as implementation proceeds. A discovery that changes a scientific or public contract returns to design review rather than being absorbed into the plan.

## Goal and success signal

Add advisory quantitative-annotation semantics to native LD-score and partitioned-h2 artifacts, persist per-model coefficient delete values, and deliver `ldsc quantile-h2` as a verified post-fit projection workflow.

Completion is observable when a native baseline-only or per-query fitted model can be passed to `ldsc quantile-h2` with resupplied fitted annotations, target annotations, and reference metadata; the command verifies the common reference-SNP universe, reproduces legacy quantile statistics on a deterministic fixture, reports paper-compatible `tau_star`, writes the specified diagnostics and metadata, and leaves existing binary partitioned-h2 numerics unchanged.

## Context and constraints

- Work in `ldsc_py3_restructured` on branch `restructure`.
- The current tree already contained unrelated changes in `tests/test_snp_identity.py` and `docs/current/clarification-heritability-from-ldsc-vs-pldsc.md`; they were preserved unchanged by this feature work. The later implementation request explicitly asked to commit all uncommitted repository changes, bringing them into the final commit scope.
- `AnnotationBuilder` owns heterogeneous annotation loading and projection; `ldscore_calculator.py` owns LD-score workflow assembly; `regression_runner.py` owns fitted-model orchestration; `outputs.py` owns directory contracts; `cli.py` dispatches but contains no scientific logic.
- Pure quantile, projection, standardized-coefficient, and jackknife transformations belong in a private numerical kernel. Input resolution, identity alignment, provenance checks, logging, and failure diagnostics belong in a public workflow module.
- The regression estimator and existing partitioned-h2 numerical tables do not change. Annotation classification is advisory.
- Full SNP-by-annotation matrices are never persisted. Coefficient delete matrices are small and remain float64.
- The common reference-SNP universe, coefficient order, identity mode, genome build, common-MAF operator/threshold, and prevalence are immutable fitted-model provenance.
- SHA-256 is used only for the approved semantic identity check over common-reference-SNP identities and fitted annotation values. It is not general file provenance.
- Existing dependencies are sufficient. Do not add packages.
- The user-approved `aggregate_only` path for older LDSC3/converted LD-score artifacts is an explicit compatibility exception.
- Apply TDD per vertical slice: create the smallest failing behavioral or numerical test, implement, then refactor while green.

## Interfaces and invariants

### Existing surfaces to expand

- `LDScoreResult` and LD-score `metadata.json` gain annotation classifications and optional `annotation_fingerprints` without changing LD-score parquet schemas.
- Annotation names are globally unique across baseline and query groups in every in-memory and loaded LD-score artifact.
- `PartitionedH2BatchResult` carries coefficient delete values for the baseline-only model or each written per-query model.
- `PartitionedH2DirectoryWriter` writes `coefficient_delete_values.parquet` beside the model that owns it and exposes the path, block count, and retained order through metadata.
- Partitioned-h2 logs the persisted annotation classifications once per run but does not branch numerical calculation on them.

### New public surface

- `ldsc quantile-h2` uses the exact flags and output schemas in the specification, requires `--output-dir`, and follows standard output preflight, overwrite, log-level, and diagnostic metadata conventions.
- The fitted result directory is authoritative for model identity and coefficient order; resupplied sources provide values only.
- The command accepts a baseline-only result root or one per-query result directory and rejects aggregate multi-query roots.
- Failure diagnostics are durable even when headline results are not committed.

### Numerical shapes and invariants

- `S`: fitted annotation sums by quantile, shape `C x Q`.
- `tau`: point coefficient vector, shape `C`.
- `T_delete`: coefficient delete values, shape `B x C`.
- `h = S.T @ tau`, shape `Q`; `H_delete = T_delete @ S`, shape `B x Q`.
- No stage may infer coefficient order from source-file column order.
- No jackknife stage may discard an invalid block selectively.
- Quantile assignment, tie behavior, normal P values, and fixed-scale `tau_star` follow the specification and legacy anchors exactly.

## Implementation slices

### Slice 1: Expand annotation semantics and LD-score metadata

Status: completed.

Goal: native LD-score construction emits validated, reusable annotation semantics and exact common-value fingerprints while preserving all current numerical outputs.

Likely areas: `src/ldsc/_kernel/annotation.py`, `src/ldsc/annotation_builder.py`, `src/ldsc/ldscore_calculator.py`, `src/ldsc/outputs.py`, LD-score loader validation in `src/ldsc/regression_runner.py`, and focused annotation/LD-score/output tests.

Work:

- Reject nonfinite fitted annotation values in the shared annotation validator.
- Add exact `binary`/`quantitative` classification over retained reference-SNP values and carry it through chromosome aggregation without rescanning baseline annotations per query regression.
- Enforce global annotation-name uniqueness across baseline and query groups at construction, result validation, write, and load boundaries.
- Implement one shared, versioned canonical-frame fingerprint helper for sorted effective identities and normalized float32 values over the common reference-SNP universe.
- Persist annotation classifications and fingerprints in LD-score metadata; reconstruct `unknown` and missing-fingerprint states for supported older artifacts.
- Add the compact LD-score classification log block.

Validation:

- Unit tests cover exact binary values, discrete quantitative values, signed/centered values, NaN/infinity rejection, duplicate names across groups, row-order-independent fingerprints, value-sensitive fingerprints, and excluded rare-SNP changes that do not alter common-value fingerprints.
- LD-score writer/loader round-trip tests assert metadata schema and unchanged LD-score/count/overlap numerics.
- Run `pytest tests/test_annotation.py tests/test_ldscore_workflow.py tests/test_output.py tests/test_overlap_matrix.py -q`.

Checkpoint: verify fingerprint canonicalization is stable across whole-genome versus `@`-sharded inputs before later slices depend on it. Revise the plan if the current chromosome aggregation cannot produce identical canonical content without retaining a large matrix.

### Slice 2: Persist per-model coefficient delete values and interpretation logs

Status: completed.

Goal: every model eligible for later quantile processing owns a self-describing coefficient delete artifact.

Likely areas: `src/ldsc/regression_runner.py`, `src/ldsc/outputs.py`, `tests/test_regression_workflow.py`, and `tests/test_output.py`.

Work:

- Extend the partitioned-h2 result carrier with one baseline-only delete matrix and a per-query mapping, keeping model column order explicit.
- Extract `Hsq.part_delete_values` immediately after each fit and pair it with `dataset.retained_ld_columns`.
- Write float64 wide parquet with zero-based `delete_block`; do not route this artifact through helpers that narrow float64 to float32.
- Add baseline-only and staged per-query file paths and metadata fields without changing current headline TSV schemas.
- Add the once-per-run quantitative interpretation warning using persisted classifications, including quantile-h2 and binary-bin-refit guidance.
- Preserve aggregate-only cell-type behavior when per-query results are not requested; document that such roots are not quantile-h2 model inputs.

Validation:

- Hand-built regression-result tests assert exact point/delete coefficient scaling, block rows, column order, and float64 persistence.
- Writer collision, overwrite, stale-tree, and staged per-query failure tests include the new files.
- Regression-equivalence tests compare all pre-existing binary and quantitative summary fields before and after the result-carrier expansion.
- Run `pytest tests/test_regression_workflow.py tests/test_output.py tests/test_kernel_regression.py -q`.

Checkpoint: confirm two-step and ordinary Hsq fits expose compatible coefficient delete semantics before declaring this slice complete.

### Slice 3: Resolve one fitted model and align resupplied inputs

Status: completed.

Goal: establish a verified, deterministic common-reference-SNP table and fitted annotation matrix for one model, with complete row-level diagnostics.

Likely areas: a new public workflow module `src/ldsc/quantile_h2.py`, `src/ldsc/config.py`, existing annotation/path/identity helpers, and a new `tests/test_quantile_h2.py`.

Work:

- Add a configuration/result vocabulary for fitted-model directory, fitted annotation sources, target source/name, reference metadata sources, missing token, quantile count, and output policy.
- Resolve baseline-only versus per-query model metadata and reject aggregate multi-query roots or incomplete result families.
- Reuse existing annotation-source resolution and BED/gene projection behavior, selecting only retained fitted columns and reordering them to model order.
- Load whole-genome or `@`-sharded reference metadata, inherit the exact common-MAF rule, and build effective identities with existing identity helpers.
- Align by effective identity rather than row position; implement baseline-grid intersection semantics, unambiguous allele inference, duplicate collection, target coverage validation, and numeric/raw missing-token behavior.
- Collect all safely discoverable alignment issues at each validation gate into the stable issue schema.
- Perform exact common-value fingerprint checks when present and fixed-tolerance count/overlap checks in every case; record `exact_common_values` or `aggregate_only`.
- Design failure output so the issue sidecar and log survive fatal validation while headline results and successful metadata remain absent.

Validation:

- Fixture matrices cover reordered rows, `@` shards, extra source columns, extra reference/annotation rows, missing target rows, duplicate identities, absent alleles, missing MAF, numeric and string missing tokens, global-name collisions, exact fingerprint matches/mismatches, and aggregate-only fallback.
- Tests assert exhaustive issue rows and `excluded` versus `fatal` actions, including header-only clean diagnostics.
- Run `pytest tests/test_quantile_h2.py tests/test_annotation.py tests/test_snp_identity.py -q` without modifying unrelated user changes in `tests/test_snp_identity.py`.

Checkpoint: inspect peak memory and number of annotation scans on a representative multi-column source. If the workflow retains more than the selected fitted columns plus identity metadata, revise before adding numerical projection.

### Slice 4: Implement the numerical quantile and standardized-coefficient kernel

Status: completed.

Goal: reproduce legacy quantile calculations and paper-compatible `tau_star` through pure array/data-frame inputs.

Likely areas: a small private kernel such as `src/ldsc/_kernel/quantile_h2.py`, shared jackknife/liability helpers, and numerical tests in `tests/test_quantile_h2.py`.

Work:

- Implement legacy rounded-order-statistic boundaries, interval closure, lower-bin tie handling, and empty-bin rejection.
- Accumulate `n_snps`, target boundaries, and the `C x Q` fitted annotation-sum matrix in one grouped pass.
- Compute point and delete quantile heritability by matrix multiplication; derive total, proportion, enrichment, SEs, and inside-versus-complement normal P values using the existing pseudovalue convention.
- Implement negative/zero total and invalid-delete behavior without dropping blocks.
- Compute population annotation SD, point total over the complete common reference-SNP universe, fixed-scale `tau_star`, duplicated two-sided zero-tests, and nonpositive-total behavior.
- Reuse the existing liability conversion factor for absolute quantile estimates and SEs.

Validation:

- Hand-computed small matrices verify `S.T @ tau`, `T_delete @ S`, quantile totals, ratios, contrasts, and standardized coefficients.
- Port deterministic fixtures from documented `quantile_M.pl` and `quantile_h2g.r`, including boundary ties and non-200 block counts, and compare all seven legacy statistics at strict justified tolerances.
- Add tests for external targets, fitted targets, target-name/vector mismatch, negative totals, zero totals, invalid delete denominators, nonpositive standardization totals, and liability scaling.
- Run `pytest tests/test_quantile_h2.py tests/test_kernel_regression.py -q`.

Checkpoint: require exact agreement with the legacy fixture and independently hand-computed matrix examples before exposing the CLI.

### Slice 5: Deliver the quantile-h2 writer and CLI end to end

Status: completed.

Goal: expose the verified workflow through the canonical package CLI and fixed result family.

Likely areas: `src/ldsc/quantile_h2.py`, `src/ldsc/outputs.py`, `src/ldsc/cli.py`, public exports if appropriate, `tests/test_output.py`, and CLI/integration tests.

Work:

- Add the exact parser surface, mutual-exclusion rules, shared projection options, required output directory, and CLI dispatch without moving logic into `cli.py`.
- Add a `QuantileH2DirectoryWriter` that validates stable column order, preflights the whole family, writes headline TSVs and successful diagnostic metadata atomically, and cooperates with the fatal-diagnostic path from slice 3.
- Populate all required metadata and log fields, including fitted-model provenance, boundary/missingness rules, P-value conventions, classifications, fingerprint canonicalization, and verification level.
- Ensure overwrite removes stale owned siblings only after a successful replacement and never erases useful failure diagnostics prematurely.
- Add baseline-only, per-query, fitted-target, and external-target real entry-point tests.

Validation:

- Writer schema/collision/stale-output tests and CLI help/error tests pass.
- Run real minimal commands for baseline-only and per-query fixtures and inspect every TSV, parquet, JSON, diagnostic, and log artifact.
- Run `pytest tests/test_quantile_h2.py tests/test_output.py tests/test_regression_workflow.py -q`, `ldsc quantile-h2 --help`, and `python -m ldsc quantile-h2 --help`.

Checkpoint: manually inspect a successful and fatal output directory to confirm the intended atomicity boundary is understandable and machine-readable.

### Slice 6: Documentation, compatibility audit, and repository validation

Status: completed.

Goal: make the implemented contracts navigable, reproducible, and safe for release.

Likely areas: `README.md`, `docs/current/partitioned-ldsc-workflow.md`, `docs/current/partitioned-h2-results.md`, `docs/current/data-flow.md`, `docs/current/artifact-metadata-field-inventory.md`, `docs/current/io-argument-inventory.md`, `docs/current/code-structure.md`, `docs/current/layer-structure.md`, `docs/troubleshooting.md`, and a focused tutorial/example.

Work:

- Document annotation classifications, why quantitative-row weighted summaries remain visible but are not category enrichments, and the distinction between post-fit quantile projection and refitted binary bins.
- Document coefficient delete artifacts, quantile-h2 inputs/outputs, missing tokens, fingerprints, aggregate-only verification, diagnostics, common reference-SNP universe terminology, and all numerical edge cases.
- Explicitly document Parquet-R2 `chr*_meta.tsv.gz` inputs and PLINK `ldscore --export-ref-metadata` preparation.
- Update architecture/data-flow/artifact inventories and troubleshooting anchors for run-aborting multi-cause errors.
- Add a reproducible small example that shows a fitted target and an external target without bundling large reference matrices.
- Audit legacy references against the original Perl/R scripts and record any intentional naming or artifact-layout differences.

Validation:

- Run focused documentation-linked examples and verify every shown command against `--help`.
- Run `pytest`.
- Run `python -m unittest discover -s tests -p 'test*.py' -v`.
- Run `ldsc --help`, `python -m ldsc --help`, and the end-to-end baseline/per-query smoke workflows.
- Run `git diff --check` and inspect `git status --short` to confirm unrelated user changes remain untouched.

## Risks and checkpoints

- Fingerprint stability is the earliest high-risk checkpoint; do not build downstream compatibility around an unverified canonicalization.
- Annotation/reference intersection must reproduce the original LD-score universe exactly. Aggregate counts and overlaps are required even after hashes pass because they make mismatches diagnosable.
- Global annotation-name uniqueness may expose previously latent invalid artifacts. Fail with actionable source/group information rather than auto-renaming.
- Failure diagnostics intentionally survive an aborted workflow. Keep the successful-result metadata boundary unambiguous so partial output is never mistaken for a valid result.
- Large text annotation suites can dominate runtime and memory. Read only identities and selected model columns where parsers permit, avoid repeated SNP scans inside jackknife loops, and measure before optimizing further.
- Preserve coefficient and delete matrices in float64; narrowing them would undermine jackknife and aggregate comparisons.
- The plan must be revised if legacy fixtures reveal a difference in boundary rounding, jackknife pseudovalue convention, or P-value calculation.

## Completion evidence

- `conda run -n ldsc3-dev pytest -q`: 1,280 passed, 1 skipped, 113 subtests passed.
- `conda run -n ldsc3-dev python -m unittest discover -s tests -p 'test*.py' -v`: 1,005 passed, 1 skipped.
- `python -m ldsc quantile-h2 --help` and unified `python -m ldsc --help` expose the approved command and flag names.
- `git diff --check` passes.
- Focused numerical and workflow coverage includes legacy boundary ties, empty quantiles, matrix/jackknife projection, fixed-scale `tau_star`, negative/nonpositive totals, explicit `NaN` missing-token exclusion, baseline-only and per-query model loading, aggregate-root rejection, delete-value persistence, global annotation-name uniqueness, and deterministic SHA256 fingerprints.

## Out of scope

- Changes to the regression estimator, current partitioned-h2 numerical rows, overlap formulas, or coefficient tests.
- Automatic normalization, MAF adjustment, binary-bin construction/refitting, plotting, or meta-analysis.
- Raw LDSC2 regression-prefix support in `ldsc quantile-h2`.
- Full annotation-matrix persistence, PLINK MAF recomputation, or a post-fit common-MAF override.
- Unrelated cleanup of the current dirty worktree; pre-existing changes were preserved and included only because the user explicitly requested one commit of all uncommitted changes.
