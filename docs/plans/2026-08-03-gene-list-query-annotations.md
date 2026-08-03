# Gene-List Query Annotations Implementation Plan

Last updated on: 2026-08-03

Status: ready for implementation; this plan does not itself authorize code changes

## Goal and governing specification

Implement direct gene-list query annotations for `ldsc ldscore`, with BED-equivalent interval projection and batch-safe per-query diagnostics, as specified in `docs/specs/2026-08-03-gene-list-query-annotations-design.md`.

The implementation stays in the `ldsc_py3_restructured` worktree on the `restructure` branch. It preserves the existing split baseline/query LD-score artifacts, overlap-aware partitioned-h2 contracts, BED naming, public exports, and SNP identity metadata.

## Architecture summary

The workflow gains one internal resolver seam. It loads and validates the packaged catalog once, maps each list to canonical Ensembl genes and selected-build intervals, and returns pure in-memory resolution records. `AnnotationBuilder` uses the existing interval-overlap primitive to project usable intervals. LD-score orchestration owns query-local status transitions, filters zero-hit/zero-variance queries, and carries ordered provenance/status data into the output layer. `LDScoreDirectoryWriter` preflights and writes the two new diagnostics plus additive root metadata; scientific tables and overlap data contain usable queries only.

The main implementation risk is not gene parsing itself. It is preserving one ordered query identity across source resolution, projection, chromosome computation, zero-variance filtering, count records, overlap labels, metadata, and diagnostics. Each slice below adds boundary assertions/tests before extending the next layer.

## Repository conventions and constraints

- Use the workflow layer for path resolution, build inference, catalog loading, status classification, logging, and output ownership. Keep numeric overlap and LD-score kernels unaware of gene identifiers and files.
- Reuse `resolve_file_group`, shared genome-build normalization/inference, BED `[start0, end)` projection, output-family preflight, and workflow logging rather than creating parallel conventions.
- Preserve `GlobalConfig.genome_build=None` in rsID artifacts; catalog projection build is separate provenance.
- Extend frozen dataclasses with defaults where compatibility requires it and use `dataclasses.replace(...)` when carrying new fields through transformations.
- Preserve explicit/glob query order. Assert alignment among query columns, count records, overlap labels, query provenance, and status records at assembly/writer boundaries.
- Do not bump the shared artifact `schema_version` for additive LD-score metadata.
- Do not modify the user’s unrelated existing edits, including the current `docs/current/io-argument-inventory.md` work, without reconciling them during the documentation slice.
- Use focused red-green-refactor cycles within each slice and run broader regressions at the slice checkpoint.

## Delivery slices

### Slice 1 — Package and validate the catalog; implement pure resolution

**Outcome:** A tested internal resolver converts one readable list into canonical catalog indices/intervals, counts, and unresolved records without logging or I/O side effects beyond reading inputs.

**Primary files:**

- Add `src/ldsc/data/protein_coding_genes.tsv.gz`.
- Add `src/ldsc/gene_list_resolver.py`.
- Update `src/ldsc/data/readme.txt` and package-data checks if needed; `setup.py` already includes `ldsc/data/*.tsv.gz`.
- Add focused tests, preferably `tests/test_gene_list_resolver.py`, and extend `tests/test_package_layout.py`.

**Work:**

1. Convert the prepared catalog into the approved compressed TSV resource without changing row identity or coordinates. Validate its 18,401-row source expectation, unique `ensgid`, build triplet rules, autosomal chromosomes, and coordinate bounds. Record the decompressed-content SHA-256 in a test fixture/expectation rather than duplicating it in runtime source constants.
2. Define small frozen internal records for catalog identity, per-query counts, unresolved rows, and resolved query data. Keep record field order explicit and serializable, but do not expose these types from `ldsc.__init__`.
3. Implement one catalog loader that normalizes missing cells, validates the full schema before returning any query result, builds the exact Ensembl and gene-name maps once, and computes the stable content checksum.
4. Implement complete-file parsing and resolution: whitespace/blanks, one-field enforcement, numeric Ensembl version stripping, exact/case-sensitive namespace order, mixed IDs/names, cross-namespace/name ambiguity, selected-build missing intervals, canonical deduplication, and all approved counts.
5. Keep error classification data-oriented. Catalog corruption raises a global input error; query problems are returned as unresolved/structural records for orchestration to classify later.

**Validation checkpoint:**

- Catalog load/installability and checksum tests.
- Unit tests for plain/gzip parity, duplicate and alias collapse, both genome builds, all unresolved reasons, full-file diagnostic collection, empty input, malformed rows, ambiguous names/namespaces, and deterministic record order.
- Memory-shape assertion that the resolver returns compact indices/intervals and does not build a gene-by-query dense matrix.

### Slice 2 — Add CLI/config source resolution and catalog-build selection

**Outcome:** CLI and Python entrypoints accept the new mutually exclusive source group, resolve deterministic query names, enforce baseline/name preflight, and select a concrete projection build in every gene-list run.

**Primary files:**

- Update `src/ldsc/config.py`.
- Update `src/ldsc/ldscore_calculator.py`.
- Update `src/ldsc/path_resolution.py` only if a small reusable basename helper belongs there; otherwise keep gene suffix handling with the resolver/workflow.
- Extend `tests/test_config_identifiers.py`, `tests/test_path_resolution.py`, and `tests/test_ldscore_workflow.py`.

**Work:**

1. Add `query_annot_gene_list_sources` to `AnnotationBuildConfig`, its normalization/docs, `run_ldscore(...)` accepted keywords/defaults, CLI argument normalization, and the three-way mutually exclusive parser group.
2. Update the explicit-baseline guard so all three query routes require `--baseline-annot-sources`. Preserve the no-query synthetic-base path.
3. Resolve comma-separated exact/glob gene-list tokens deterministically, without `@` suite semantics. Distinguish a token matching no files (global usage error) from a resolved concrete file that later proves unreadable (query-local status).
4. Derive gene-list names by removing optional `.gz` plus one `.txt`/`.tsv`/`.list`; retain current BED `Path.stem` behavior. Preflight duplicate query names and baseline-column collisions before calculation.
5. Refactor build resolution narrowly so gene-list runs treat omitted `--genome-build` as `auto` in all SNP identity families. Reuse baseline/R2 evidence and agreement rules. Store the concrete projection build separately from `GlobalConfig` identity metadata in rsID modes.
6. Add INFO log data for declared/inferred catalog build and evidence, leaving warning emission to the later status slice.

**Validation checkpoint:**

- Parser mutual-exclusion and explicit-baseline tests for CLI and Python wrapper.
- Exact/glob ordering, path deduplication, suffix examples, and collision-message tests.
- Explicit alias and auto-build tests in all four SNP identifier modes, including insufficient/conflicting evidence and rsID metadata remaining build-independent.
- Help smoke tests for both `ldsc ldscore --help` and `python -m ldsc ldscore --help`.

### Slice 3 — Project gene intervals and introduce shared BED/gene query statuses

**Outcome:** BED and gene-list files are handled independently, usable sources produce ordered binary annotations, and query-local parse/read/zero-hit failures do not abort usable siblings.

**Primary files:**

- Update `src/ldsc/annotation_builder.py`.
- Update `src/ldsc/ldscore_calculator.py`.
- Reuse `src/ldsc/_kernel/annotation.py` without adding gene concepts unless a general interval-union primitive is genuinely missing.
- Extend `tests/test_annotation.py` and `tests/test_ldscore_workflow.py`.

**Work:**

1. Introduce an internal ordered query-status/provenance record carried with the annotation workflow. Add compatible default fields to `AnnotationBundle` only if that is the narrowest way to preserve records through chromosome sharding/concatenation.
2. Resolve gene lists once, select unique build intervals, apply symmetric padding once, union intervals, and call the same binary SNP-overlap path used by BED projection. Do not materialize BED or generated annotation files.
3. Refactor multi-BED loading so a malformed, empty, or unreadable concrete BED can become one skipped status while siblings continue. Preserve existing BED format rules and basename derivation.
4. Classify empty, malformed, fully unresolved, ambiguous, unreadable, and partially resolved sources. Retain all parsed problem rows. Keep partially resolved gene queries usable with `warning`.
5. After reference-SNP restrictions and region exclusions, compute `n_annotation_snps` from the exact count universe and skip zero-hit queries before scientific result assembly.
6. Assert that the remaining annotation columns and ordered usable status records agree before chromosome dispatch and after chromosome bundle concatenation.

**Validation checkpoint:**

- A hand-built list and equivalent BED interval union produce identical annotation vectors with padding 0 and nonzero padding, including start clipping and overlapping genes.
- Mixed clean/warning/skipped BED and gene sources preserve source order and allow usable siblings to proceed.
- Empty, malformed, ambiguous, fully unresolved, unreadable, and zero-hit cases have exact statuses/reasons and null-versus-zero count semantics.
- Existing annotation and BED naming/projection tests remain green.

### Slice 4 — Filter zero-variance queries and preserve LD-score/count/overlap alignment

**Outcome:** Scientific results contain only queries that remain estimable on regression rows; all dependent result structures are pruned consistently.

**Primary files:**

- Update `src/ldsc/ldscore_calculator.py` result assembly.
- Update `src/ldsc/annotation_builder.py` only for status handoff if needed.
- Update `src/ldsc/overlap_matrix.py` only through existing label-based selection helpers if pruning support is needed.
- Extend `tests/test_ldscore_workflow.py`, `tests/test_overlap_matrix.py`, and `tests/test_regression_workflow.py`.

**Work:**

1. Evaluate query LD-score variance on the final written regression SNP rows after chromosome aggregation. Mark constant queries `skipped/zero_variance_ld_scores`.
2. Prune skipped query columns from query tables, `query_columns`, count records, overlap labels/blocks, and any per-chromosome result views using label-aware selection. Use `dataclasses.replace(...)` to retain new status/provenance fields when replacing result objects.
3. Add explicit boundary validation that every final query column has one count record, one usable status/provenance record, and the expected overlap representation in the same order.
4. Confirm `partitioned-h2` needs no fabricated failure rows: it sees only the final usable query columns and continues to use current baseline-plus-one-query models.
5. Implement the all-skipped decision at the workflow boundary: prepare diagnostics, suppress canonical output writing, and raise one consolidated user-facing error only after diagnostic ownership is established in Slice 5.

**Validation checkpoint:**

- Zero-variance and zero-hit tests prove no query parquet/count/overlap/metadata leakage.
- Multi-query tests prove removing a middle query preserves the order and values of surviving columns across all structures.
- Gene-list/BED equivalence covers exact LD scores, count records, overlap values, and regression row alignment.
- Partitioned-h2 regression tests prove only usable queries generate result rows.

### Slice 5 — Own diagnostics, warnings, provenance, and partial-failure writes

**Outcome:** Successful and all-skipped runs emit the exact status/audit/log contract, with atomic/collision-safe ownership and concise root metadata.

**Primary files:**

- Update `src/ldsc/outputs.py`.
- Update `src/ldsc/ldscore_calculator.py` orchestration.
- Update `src/ldsc/_logging.py` only if the existing workflow-log context cannot preserve an intentional post-diagnostic error.
- Extend `tests/test_output.py`, `tests/test_ldscore_workflow.py`, and `tests/test_logging_refactor.py`.

**Work:**

1. Add serializers for `diagnostics/query_annotation_status.tsv` and gzip `diagnostics/gene_list_unresolved.tsv.gz` with fixed columns, deterministic row order, empty-field semantics, and a header-only clean audit.
2. Extend LD-score output-family preflight/overwrite ownership to include both diagnostics when applicable. Preserve the current rule that the workflow log is owned by orchestration rather than returned in `output_paths`.
3. Extend `LDScoreResult` or the writer input contract with defaulted ordered status, catalog, and query-provenance fields. Add only the approved concise root metadata and relative diagnostic references.
4. Emit one WARNING per non-`ok` query and a status-count summary. Include the effective catalog build/evidence in the file log without dumping unbounded unresolved identifiers to the console.
5. Implement a diagnostics-only all-skipped commit path: preflight the owned family, write status/audit/log, raise the consolidated input error, and ensure no root metadata/parquet is present. With `--overwrite`, remove stale canonical outputs only after the diagnostics-only result has been written successfully; without it, reject collisions before mutation.
6. Test failure injection around diagnostic and canonical writes so a failed run does not silently leave a mixed current/stale artifact family.

**Validation checkpoint:**

- Exact schema/content tests for clean, partial, malformed, ambiguous, unreadable, zero-hit, and zero-variance diagnostics.
- Metadata tests for one catalog record, ordered compact query records, basename-only sources, stable checksum definitions, relative diagnostic paths, and no successful gene list.
- Log capture tests for every non-ok warning, build evidence, summary counts, and consolidated all-skipped error.
- Collision/overwrite/stale-artifact tests for successful and diagnostics-only failure paths.

### Slice 6 — End-to-end compatibility, documentation, and handoff

**Outcome:** The feature is documented, installed-package safe, regression compatible, and ready for the separate sparse-index design session.

**Primary files:**

- Update `docs/current/io-argument-inventory.md` after reconciling existing user edits.
- Update `docs/current/bed-input-format.md`, `docs/current/config-design.md`, `docs/current/code-structure.md`, `docs/current/layer-structure.md`, `docs/current/data-flow.md`, `docs/current/artifact-metadata-field-inventory.md`, and `docs/current/snp-identifier-genome-build-defaults.md` where their authoritative contracts change.
- Update `docs/wiki/guided-tutorial.md` and `docs/troubleshooting.md`.
- Update packaging/release documentation only where the catalog resource or user-facing flag requires it.
- Extend end-to-end tests in `tests/test_ldscore_workflow.py`, `tests/test_regression_workflow.py`, and `tests/test_package_layout.py`.

**Work:**

1. Add end-to-end CLI fixtures for multiple gene lists, mixed BED validity, both builds, rsID auto inference, partial resolution, all-skipped behavior, output collision handling, and a downstream partitioned-h2 run.
2. Run installed-resource tests so catalog loading does not depend on a repository-relative path.
3. Document the one-column format, exact identifier rules, canonical Ensembl identity, coordinates/padding, build inference, basename naming, diagnostics, compact provenance, and BED-equivalence guarantee.
4. Add the approved tutorial warning: if an input query gene list or BED file is absent from scientific results, that query failed; inspect `diagnostics/query_annotation_status.tsv` for the reason and `diagnostics/gene_list_unresolved.tsv.gz` for gene-level details.
5. Document that no all-protein-coding control is added in this design and that a future persisted gene index must match the catalog checksum.
6. Perform a final leakage audit: no physical BEDs, no dense gene-by-query matrix, no skipped query columns/counts/overlaps/result rows, no absolute paths in metadata/sidecars, and no rsID projection build in shared identity metadata.
7. Write the requested next-session handoff in an OS temporary directory, summarizing the stable resolver result seam and catalog checksum binding for the future sparse disjoint-atom index. The handoff must not scaffold or implement that future command.

**Validation checkpoint:**

- Focused feature suite plus all touched workflow/output/regression suites.
- Full `pytest` suite and repository-standard `unittest` compatibility run.
- CLI help and representative end-to-end commands in both console entrypoints.
- `git diff --check`, package-resource import smoke test, and documentation link/path check.

## Suggested validation commands

Use the repository development environment for every command:

```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh
conda activate ldsc3-dev

pytest -q tests/test_gene_list_resolver.py tests/test_annotation.py tests/test_ldscore_workflow.py tests/test_overlap_matrix.py tests/test_output.py tests/test_regression_workflow.py
pytest -q
python -m unittest discover -s tests
ldsc ldscore --help
python -m ldsc ldscore --help
git diff --check
```

If `tests/test_gene_list_resolver.py` is named differently during implementation, substitute the chosen focused test module without changing the coverage checkpoint.

## Completion definition

The implementation is complete only when every acceptance criterion in the governing specification has an automated test or a documented manual check, all repository tests pass, current docs describe the shipped behavior, and representative gene-list/BED-equivalent runs produce matching scientific artifacts. A passing resolver unit suite alone is not sufficient; the status/diagnostic output family and downstream partitioned-h2 behavior are part of the feature.

## Risks and checkpoints requiring design escalation

- If existing logging/output preflight cannot safely support diagnostics-only all-skipped failures, pause before inventing a new partial-commit policy and update the specification.
- If catalog data do not satisfy the approved schema or the 12-pc preparation record cannot justify a field, preserve the source limitation and request a catalog decision rather than silently transforming identity/coordinates.
- If zero-variance pruning cannot keep counts and overlap labels aligned without changing public result objects, prefer additive defaulted fields and label-based selection; any public schema break requires review.
- If rsID auto inference lacks sufficient baseline/R2 evidence in a supported backend, retain the specified actionable error. Do not infer build from gene identifiers or region presets.
- If implementation suggests accepting mixed query source types, changing BED names, or emitting placeholder failed rows, treat that as out of scope and seek an explicit design amendment.

No product ambiguities remain. Implementation may choose private helper names and exact test-fixture organization as long as the public contracts and layer ownership above remain unchanged.
