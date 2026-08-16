# Defensive Gene-List Input Implementation Plan

Last updated on: 2026-08-16

Status: implemented and repository-validated on 2026-08-16

Governing specification: [Defensive Gene-List Input Specification](../specs/2026-08-15-defensive-gene-list-input-design.md)

## Goal and success signal

Replace the packaged-catalog, per-source partial resolver with one defensive, build-aware gene-list workflow shared by direct and indexed `ldscore` modes.

The implementation succeeds when strict mode stops before LD-score calculation on any rejected focal/control gene, `resolved-only` deliberately continues with the audited usable subset, all validation gates report the complete batch, the index builder rejects any noncanonical catalog before atom construction, and direct/indexed outputs remain scientifically equivalent for matching inputs. The feature is not complete until its diagnostics can be used to repair lists/catalogs, the packaged protein-coding catalog and obsolete schemas are gone, and representative many-pathway resolution is vectorized rather than implemented as a Python loop over genes.

## Progress

| Slice | State |
| --- | --- |
| 1. Shared catalog and batch-resolution foundation | Complete |
| 2. Direct CLI, build authority, and Gate A | Complete |
| 3. Gate B, result diagnostics, and direct scientific completion | Complete |
| 4. Canonical index construction and builder diagnostics | Complete |
| 5. Indexed resolution, assembly, and direct/index parity | Complete |
| 6. Contract cleanup, documentation, and repository validation | Complete |

Implementation evidence is recorded below. Do not preserve stale implementation detail when code inspection or tests show a better internal boundary.

## Context and constraints

- Work only in `ldsc_py3_restructured` on branch `restructure`.
- Preserve the unrelated existing modification in `tests/test_snp_identity.py` unless the user separately brings it into scope.
- Before this implementation, `gene_list_resolver.py` owned the packaged `GeneCatalog`, SHA-256 provenance, per-file parsing, and a Python per-row resolver. `AnnotationBuilder` invoked it separately for focal lists and the control and could fail on the control before focal issues were examined.
- Direct gene resolution previously occurred before `run_ldscore_from_args()` created/preflighted the output directory and entered `workflow_logging`; Gate A therefore had to move inside that owned diagnostic boundary.
- `LDScoreDirectoryWriter` previously wrote query status and `gene_list_unresolved.tsv.gz`, while orchestration duplicated the LD-score output-family prediction. The replacement diagnostic family required one authoritative path/ownership rule so overwrite cleanup could not delete newly written artifacts.
- The prior indexed builder loaded the packaged catalog, exposed public `--chromosomes`, embedded a 0-based catalog, and recorded a catalog checksum. The online loader reconstructed a weaker catalog through `from_index_frame`; both paths had to migrate to the shared validator.
- The CLI console handler normally emits only ERROR when a workflow log exists. Bounded successful gene-list notices therefore required a narrow one-shot console route; globally lowering the console threshold would leak complete per-gene warnings and violate the bounded-console contract.
- Preserve LD-score numerical kernels, adjusted-\(r^2\), reference/regression SNP universes, overlap semantics, row order, and downstream partitioned-h2 behavior. This work changes gene resolution and orchestration, not the estimand.
- Preserve current index constraints not changed by the new specification, including explicit hg19 and base `rsid`/`chr_pos` construction. Public index coverage becomes fixed to autosomes 1–22; partial coverage survives only in a private test seam.
- Use red-green-refactor within each slice. The plan uses expand–migrate–contract so the branch remains testable while direct and indexed consumers move to the new internal model; no legacy path remains at completion.

## Interfaces and invariants

### Public inputs

- Add live/indexed `--gene-list-resolution-policy {strict,resolved-only}`, default `strict`.
- Add required live/builder `--gene-coordinate-file` and remove every packaged fallback.
- Remove public builder `--chromosomes`; public construction and loading require chromosomes 1–22.
- Preserve one `ldsc ldscore` command and enforce the specification's four mutually exclusive query modes.
- Preserve omission versus explicit zero for `--padding-bp`: live gene lists require an explicit value, live BED defaults to zero, and prebuilt/no-query/indexed modes reject an explicit value.

### Shared internal result seam

Replace per-source resolution as the orchestration boundary with one batch result containing:

- a normalized catalog authority and its build/source metadata;
- the complete ordered gene-list audit frame;
- one ordered source-summary frame, including unreadable sources;
- selected unique catalog rows/intervals for each focal source and the control;
- source-level fatal issues and strict/resolved-only policy outcome; and
- enough stable source/query identity to derive `QueryAnnotationStatus` without rereading files.

The exact private class/helper names are not contractual. `AnnotationBundle`/`LDScoreResult` should carry batch frames and selected data without forcing the writer to understand resolver internals. Use `dataclasses.replace(...)` when enriching or pruning frozen results.

### Scientific and artifact invariants

- Catalog input coordinates are 1-based inclusive and convert exactly once to `[start - 1, end)` internally.
- ID/name matching is exact and case-sensitive after surrounding-whitespace cleanup; identifier versions are never stripped.
- MHC exclusion uses unpadded intervals and precedes padding.
- Strict and resolved-only modes produce the same row audit; policy changes the run consequence, not resolver truth.
- Gate A completes across focal and control sources before any LD-score computation. Gate B completes across all surviving genes/queries after retained-panel load and before LD-score computation.
- Zero reference-SNP support is not failed resolution. Query zero-hit and zero-variance states are query-local; corresponding control states are fatal.
- Diagnostic schemas, disposition/reason vocabularies, source order, null semantics, and console cap are exactly those in the governing specification.
- Routine gene-list/catalog hashes disappear. Only the top-level semantic `index_id` remains; normalized embedded catalog content participates directly in that identity.
- No old index, old embedded-catalog schema, packaged-catalog route, or old diagnostic filename is accepted after contraction.

## Implementation slices

### Slice 1 — Expand a shared catalog and vectorized batch resolver

**Outcome:** A pure workflow-layer service can validate a user/embedded catalog and resolve all focal/control rows in bulk, producing the approved audit/summary shapes without logging or artifact writes.

**Likely areas:** `src/ldsc/gene_list_resolver.py`; a focused internal catalog module if separating catalog validation from list resolution keeps concerns smaller; shared build/chromosome normalization helpers; new deterministic fixtures under `tests/fixtures/`; `tests/test_gene_list_resolver.py`.

**Work:**

1. Add the canonical fixture catalog and focal/control sources first. Include unique IDs, duplicate IDs/names, ID/name namespace collision, invalid coordinates, unsupported chromosome, MHC overlap, duplicates/aliases, unmatched tokens, malformed rows, and genes with/without retained SNP support. Expected audit/summary tables are complete golden tables rather than exception-only assertions.
2. Implement one path-required catalog reader for TSV/TSV.GZ that preserves physical lines/raw values, performs approved lexical normalization, and returns normalized rows plus deterministic structural/row/global issue frames. Support consequence-scoped live validation and all-or-nothing index validation through one validator, not duplicated parsers.
3. Build all readable focal/control sources into one row frame with role, argument, query, basename, ordinal, physical line, and trimmed token. Iterate only over source files; use pandas/NumPy string operations for rows.
4. Resolve the row frame through vectorized ID/name joins, group operations, and masks. Derive ambiguity, catalog-line groups, canonical identity, interval validity, duplicate precedence, MHC exclusion, dispositions, details, and initial per-source counts without a Python per-gene loop or a gene-by-query dense matrix.
5. Keep source errors separate from row audit records, but emit one source-summary row per declared source. Ensure strict/resolved-only evaluation is a policy pass over the same audit.
6. Establish a compact per-source selected-gene representation (catalog indices/interval arrays or sparse membership) for direct projection and indexed selector assembly.

**Validation checkpoint:**

- Focused catalog tests cover every builder issue reason, `validation_incomplete`, raw observed values, related line groups, stable ordering, and live versus index consequence differences.
- Focused batch tests cover strict/resolved-only parity of audit rows, all approved list dispositions/reasons, duplicate precedence/reset, unreadable source summaries, and exact count identities.
- A representative hundreds-of-lists test or profile demonstrates resolution uses a bounded set of joins/group operations and has no Python call path once per submitted gene.
- Existing direct/index consumers may still use the old seam at this checkpoint; the new seam is independently green.

### Slice 2 — Migrate direct CLI, catalog build authority, and Gate A

**Outcome:** Direct gene-list mode uses the required catalog and new batch resolver, enforces the mode/padding matrix, writes complete Gate A diagnostics, and never begins LD-score calculation after a fatal preflight.

**Likely areas:** `src/ldsc/config.py`, `src/ldsc/ldscore_calculator.py`, `src/ldsc/annotation_builder.py`, `src/ldsc/outputs.py`, `src/ldsc/genome_build_inference.py` or existing shared inference seams, and focused config/annotation/workflow/output tests.

**Work:**

1. Add/normalize `gene_coordinate_file` and `gene_list_resolution_policy` in CLI/Python configuration. Retain an omission sentinel for padding until query mode is known; normalize an allowed omission only after mode validation.
2. Centralize the four-mode argument matrix so CLI and `run_ldscore()` reject the same invalid combinations. Preserve ordinary no-query behavior and BED/prebuilt semantics.
3. Reorder direct initialization: perform cheap usage/path normalization, establish the predictable gene diagnostic family and `ldscore.log`, enter workflow logging, then load the catalog and combine its build evidence with baseline/reference evidence. Preserve build-independent shared identity metadata for rsID artifacts.
4. Resolve focal and control sources together. Remove the current standalone control raise. Detect duplicate query names while continuing every safely possible source/catalog check.
5. Give the output writer one authoritative diagnostic path/family description used by early preflight, diagnostics-only writes, canonical writes, and stale cleanup. Include audit and summary even on strict Gate A failure when their inputs are knowable; do not write query status when scientific viability was not evaluated.
6. Under strict, emit the consolidated bounded error after diagnostic writes. Under resolved-only, pass only selected genes to annotation projection, create skipped empty/zero-resolved focal statuses, and enforce nonempty usable control.
7. Project resolved intervals through the existing BED-equivalent interval primitive. Gene files/catalogs are read once even when baseline annotations are chromosome-sharded.

**Validation checkpoint:**

- CLI/Python tests cover the complete mode and padding matrices, required catalog, policy scope/default, query-name collision, and build agreement/auto inference.
- Gate A integration tests combine problems in multiple focal files and the control and assert one console error, complete log, complete audit/summary, no query-status file, and no scientific output.
- Diagnostics-only write/overwrite tests assert every metadata-listed/current artifact survives replacement and stale outputs are removed only after a successful diagnostic commit.
- Annotation tests prove exact 1-based-to-0-based boundary behavior and gene/BED projection equality at padding 0 and positive padding.

### Slice 3 — Add Gate B, control viability, final statuses, and direct output provenance

**Outcome:** Direct mode evaluates retained-SNP support before LD-score work, prunes unusable focal queries in batches, protects the control model, and publishes the final diagnostics/metadata contract.

**Likely areas:** `src/ldsc/ldscore_calculator.py`, `src/ldsc/annotation_builder.py`, reference-panel preparation seams, `src/ldsc/query_annotations.py`, `src/ldsc/outputs.py`, `src/ldsc/_logging.py`, overlap/result tests, and downstream regression tests.

**Work:**

1. Prepare each chromosome's retained panel metadata once and reuse it for Gate B and the LD-score kernel. Compute per-gene support vectorially from sorted SNP positions and interval start/end arrays; avoid a gene-by-SNP matrix and avoid rereading the panel after preflight.
2. Update the audit/summary with `unsupported`, support counts, and query totals. Aggregate every zero-support gene and zero-annotation query before deciding which focal columns proceed. Enforce fatal zero-control annotation.
3. Run LD-score calculation only for surviving focal queries/control. After aggregation, prune all zero-variance focal columns together and enforce fatal zero-variance control while preserving count/overlap/column alignment.
4. Replace old per-resolution result provenance with the audit/summary/policy fields needed by `LDScoreDirectoryWriter`. Write `gene_list_audit.tsv.gz`, `gene_list_resolution_summary.tsv`, final query status, and relative metadata references; remove input/catalog SHA fields from this path.
5. Implement bounded messages centrally. Add a narrow console-notice helper for the approved resolved-only success warning so complete per-gene WARNING records remain file-only; do not lower the global console threshold.
6. Complete all-focal-skipped handling: write available audit, summary, query status, and log; raise one consolidated error; publish no root scientific artifacts.

**Validation checkpoint:**

- Gate B tests cover multiple zero-support genes/queries at once, blank support fields when Gate B is not reached, mixed usable/skipped focal batches, partial-gene plus partial-SNP status precedence, and fatal control outcomes.
- Post-computation tests prune a middle zero-variance query consistently from tables, counts, overlaps, metadata, and downstream partitioned-h2 inputs.
- Console/log tests prove the console is bounded and actionable while the file log contains every rejected/unsupported record exactly once.
- A direct end-to-end result reloads through the public LD-score loader and produces unchanged downstream partitioned-h2 behavior for surviving queries.

### Slice 4 — Migrate index construction to the canonical catalog contract

**Outcome:** The builder accepts only an explicit fully canonical catalog, validates it before atom work, writes the repairable issue artifact on failure, embeds the new catalog schema, and publishes only complete 1–22 indexes.

**Likely areas:** `src/ldsc/gene_ldscore_index.py`, `src/ldsc/config.py`, index tests/fixtures, build-state/log lifecycle helpers, and index metadata/identity tests.

**Work:**

1. Require `--gene-coordinate-file`; remove public `--chromosomes` from parser/help and public build configuration. Keep a private explicit-coverage seam solely for fast tests/prototypes, and make the public loader require exact 1–22 coverage.
2. Extend build-state preparation to archive the previous `gene_coordinate_catalog_issues.tsv.gz` with its historical log. Preflight the destination, initialize logging/build state, validate the entire catalog, and only then create a transaction or invoke chromosome/atom workers.
3. Serialize every detectable issue with the approved schema, raw values, repair guidance, ordering, and `validation_incomplete`. On catalog failure, retain the current issue file/log and prove atom construction was never called.
4. Embed the complete canonical source catalog with one-based coordinates, physical `catalog_line`, source basename, inclusion/exclusion state, and chromosome-local row mapping. Convert to zero-based intervals only for atom construction; retain intentionally excluded genes with no selected atoms.
5. Make the normalized catalog content part of the semantic `index_id` payload without storing a separate catalog checksum. Preserve unrelated identity/integrity hashes already justified by the existing index contract.
6. Preserve transactional replacement/rollback, existing valid index on failed overwrite, chromosome worker behavior, and numerical atom/operator construction.

**Validation checkpoint:**

- Invalid-catalog tests assert exhaustive issue rows, bounded console output/path, no atom invocation, no public index, archive-on-retry, and no current issues file after success.
- Valid builder tests assert exact embedded schema/basis/source lines, MHC inclusion flags, full public coverage, and private partial-fixture behavior.
- Index identity tests show normalized catalog content changes `index_id`, row-order normalization is deterministic where specified, and no catalog checksum field is exposed.
- Existing atom/operator numerical tests remain unchanged and green.

### Slice 5 — Migrate indexed online resolution and prove direct/index parity

**Outcome:** Indexed `ldscore` uses only the validated embedded catalog, shares both resolution policies/diagnostics with direct mode, performs Gate B from index statistics, and rejects every incompatible old/corrupt/partial index.

**Likely areas:** `src/ldsc/gene_ldscore_index.py`, explicit-index dispatch in `src/ldsc/ldscore_calculator.py`, shared resolver/output modules, index/output/workflow/regression tests.

**Work:**

1. Route embedded `gene_catalog.parquet` through the full shared canonical validator. Remove `from_index_frame` weak reconstruction and reject old columns, missing lines/build/source metadata, ambiguity, invalid intervals, noncanonical order/mapping, and non-1–22 public coverage before list resolution.
2. Forward `gene_list_resolution_policy` into indexed execution and resolve focal/control sources in one batch. Gate A writes the same audit/summary and uses the same strict/resolved-only consequence rules as direct mode.
3. Compute per-gene and per-query SNP support before operator multiplication from the gene-to-atom mapping and stored atom counts using sparse/vectorized operations. Enforce the same focal skip and control failure rules as direct Gate B.
4. Assemble all surviving focal selectors and optional control without looping over genes. Reuse the existing vectorized operator multiplication, float64 accumulation, count/overlap statistics, zero-variance pruning, and canonical writer.
5. Remove indexed catalog/input hash provenance and record only the selected policy, diagnostic paths/counts, and top-level `index_id`.
6. Ensure explicit index errors never trigger catalog discovery, direct fallback, partial scientific publication, or a baseline-only success.

**Validation checkpoint:**

- Direct/index tests compare surviving query/control columns, rows, LD scores, counts, overlaps, statuses, audit/summary tables, metadata, and downstream partitioned-h2 results for the same catalog/policy.
- Strict and resolved-only indexed batches cover unmatched genes, zero-resolved focal lists, partial controls, MHC-excluded genes, zero-support genes, zero-hit/zero-variance queries, and all-focal-skipped failure.
- Loader tests reject old schema, catalog corruption, partial coverage, and identity mismatch before any list source is read.
- CLI/Python explicit-index mode accepts only the approved gene-list controls and rejects all live scientific overrides.

### Slice 6 — Contract obsolete paths, write user documentation, and validate the repository

**Outcome:** Only the new contract remains in code, package data, tests, help, and active documentation; users have one detailed repair guide and concise high-level tutorials.

**Likely areas:** `src/ldsc/gene_list_resolver.py`, `src/ldsc/data/`, `setup.py`, result/config/docstrings, repository-wide tests, `README.md`, `tutorials/`, `docs/current/`, `docs/troubleshooting.md`, and relevant `docs/wiki/` pages.

**Work:**

1. Remove `protein_coding_genes.tsv.gz`, packaged-catalog constants/loaders, Ensembl-version stripping, old per-source resolution records no longer consumed, input/catalog hash fields, `gene_list_unresolved.tsv.gz`, `partial_resolution`/old reason names, and public chromosome selection. Remove fallback/compatibility branches rather than retaining aliases.
2. Update or retire obsolete tests/fixtures. Add a repository-wide negative search so only explicit historical/supersession references may mention removed names; active code, help, current docs, tutorials, and wiki must have none.
3. Create `docs/current/gene-list-diagnostics-and-repair.md` with every diagnostic term/null rule and detailed list-first/catalog-second repair workflow. Add a concise link in `docs/troubleshooting.md`.
4. Update relevant current contracts, including gene-list input, gene index, argument inventory, configuration, class/features, data flow/layer structure, artifact metadata, workflow logging, and SNP-universe documentation.
5. Keep wiki pages high-level: explain strict versus exploratory use, show the summary-then-audit repair idea, and link to the detailed current guide. Do not duplicate schemas/reason catalogs. Update representative README/tutorial commands to require the catalog and explicit live padding.
6. Update old dated design/plan documents that would otherwise present the packaged catalog or partial production index as governing behavior with an unambiguous superseded status/link, or retire them according to repository convention. Preserve mathematical index documentation that remains valid while replacing its obsolete catalog/CLI sections.
7. Run focused, full pytest, standard-library unittest, CLI help, package-resource, artifact reload, documentation-link, old-literal, and diff checks. Record any unavoidable historical references explicitly rather than allowing accidental leftovers.

**Validation checkpoint:**

- No active code/package resource/current documentation exposes a packaged gene catalog, old diagnostic filename, old hash field, version stripping, or public partial-index flag.
- Both console entrypoints show the final flags/modes only.
- Detailed diagnostics documentation and high-level wiki links agree with generated fixtures and output schemas.
- Full `pytest` and standard-library unittest compatibility suites pass in `ldsc3-dev`; representative direct and indexed outputs reload and run partitioned-h2.

## Validation commands

Use the repository development environment and adjust only focused filenames if tests are reorganized:

```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh
conda activate ldsc3-dev

pytest -q tests/test_gene_list_resolver.py tests/test_annotation.py tests/test_ldscore_workflow.py tests/test_output.py tests/test_gene_ldscore_index.py tests/test_logging_refactor.py tests/test_regression_workflow.py tests/test_package_layout.py
pytest -q
python -m unittest discover -s tests -p 'test*.py' -v

ldsc ldscore --help
ldsc build-gene-ldscore-index --help
python -m ldsc ldscore --help
python -m ldsc build-gene-ldscore-index --help

git diff --check
```

At contraction, also use repository-wide `rg` checks for the exact obsolete resource, diagnostic, reason, hash-field, and public-flag names enumerated in Slice 6. Review each remaining match rather than assuming every historical occurrence is acceptable.

## Completion evidence

- Focused gene-list/index/workflow/output/regression/package suite: 523 tests passed with 45 subtests.
- Full pytest suite: 1,259 tests passed, one expected dependency-conditioned test skipped, and 113 subtests passed.
- Standard-library compatibility suite: 985 tests passed with one expected skip.
- Installed and module CLI help both expose the final four-mode gene-list contract; the index builder has no public `--chromosomes` option.
- Package compilation completed successfully for `src/ldsc` and `tests`.
- Changed-document link validation checked 54 relative links with no missing target.
- Repository searches found no active packaged-catalog resource/path, old gene-list diagnostic/reason name, or routine gene-list/catalog hash. Remaining `--chromosomes` mentions are documentation of its removal and tests enforcing that removal; remaining `start0` uses are internal interval/atom representations after the required one-based input conversion.
- `git diff --check` completed cleanly. The pre-existing unrelated modification to `tests/test_snp_identity.py` was preserved and excluded from this implementation.

## Risks and revision checkpoints

- **Vectorized ambiguity expansion:** ID/name candidate joins can multiply rows for ambiguous aliases. Check memory/cardinality with the canonical fixture and a many-pathway sample before locking the internal frame shape. Do not replace joins with per-gene loops to avoid the issue; aggregate candidate groups deliberately.
- **Physical line preservation:** catalog/list blank/comment handling must retain physical source lines. Verify the chosen bulk reader against CRLF, gzip, blank lines, and comments before building diagnostics on top of it.
- **Output-family ownership:** early preflight, diagnostics-only writes, canonical writes, and overwrite cleanup must share one produced/owned-path rule. Run write-then-overwrite/reload tests before proceeding past Slice 2.
- **Build inference timing:** catalog evidence is required for auto build, but catalog resolution must be inside the log/diagnostic boundary. If current global-config initialization cannot be reordered cleanly, revise the workflow boundary rather than reading the catalog twice.
- **Gate B reuse:** per-gene support must be knowable before LD-score computation without doubling panel I/O. Verify prepared chromosome metadata can be reused by both PLINK and Parquet backends; pause if a backend would require a second expensive read.
- **Bounded console routing:** do not lower the whole workflow console to WARNING. Successful resolved-only and Gate B notices use narrow bounded console summaries; complete row warnings remain in the log/audit.
- **Control integrity:** test the control at every stage. A convenient shared query-pruning helper must not silently prune `gene_control` or turn a control failure into query-local skip behavior.
- **Index identity:** removing the catalog checksum must not make catalog content disappear from semantic identity. Compare `index_id` under one-field catalog changes before removing old metadata.
- **Private partial indexes:** public builder and loader must reject partial coverage, while tests still need cheap private fixtures. Keep the escape hatch private and assert it is absent from CLI/help/public wrapper signatures.
- **Contract blast radius:** old names appear across outputs, tests, specs, current docs, wiki, tutorials, and metadata. Use the contraction search as a release gate; do not leave compatibility aliases.
- **Working tree safety:** preserve unrelated dirty files, especially `tests/test_snp_identity.py`, and reassess overlap before editing any file already changed by the user.

Any evidence requiring a change to coordinate basis, build authority, strict/resolved-only consequence classes, Gate A/B timing, control viability, diagnostic schemas, public index coverage, or the scientific LD-score operator requires a specification amendment before continuing.

## Out of scope

- fuzzy, synonym, case-insensitive, or historical-name matching;
- identifier-version stripping;
- liftover or row-level mixed-build detection;
- transcript/exon/strand/TSS/TES annotation semantics;
- non-autosomal gene LD-score analysis;
- partial production indexes or public chromosome selection;
- catalog augmentation/fallback or old-index migration;
- a minimum resolution-fraction threshold;
- changes to BED/prebuilt scientific semantics beyond the approved mode/padding validation;
- new numerical LD-score kernels or altered SNP universes; and
- publishing indexes or external artifacts.

## Completion definition

The work is complete only when every observable behavior in the specification has an automated test or recorded check; direct and indexed modes share the same catalog/resolution/diagnostic semantics; invalid catalogs cannot reach atom construction; all output families survive write/overwrite/reload validation; the packaged catalog and obsolete contracts have no active leftovers; the dedicated repair guide and high-level wiki links are shipped; vectorized many-pathway resolution is evidenced; and the focused, full pytest, unittest, CLI, package, documentation, and diff checks pass.
