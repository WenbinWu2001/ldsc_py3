# Exact disjoint-atom gene LD-score index implementation plan

Last updated on: 2026-08-03

Status: implementation complete; final validation recorded below

## Implementation progress

- [x] Slice 1 — migrated the public padding vocabulary and completed direct gene-region/control semantics. Focused annotation/resolver tests: 82 passed; affected LD-score/overlap/output/regression/region checkpoint: 343 passed.
- [x] Slice 2 — implemented the pure disjoint-atom model, bounded SNP-by-atom block iteration, exact sufficient statistics, structural validation, and float64 indexed assembly. Focused kernel checkpoint passes, including duplicate-CSR rejection.
- [x] Slice 3 — exact PLINK-backed common/profile chromosome construction. Focused pair/batching/strict-identity tests pass; both real chromosome-22 inputs retained 141,123 reference rows and produced 17,380 fixed regression rows.
- [x] Slice 4 — serialization, validation, semantic identities, and transactional publication. Round-trip, corruption, common reuse, sibling preservation, targeted overwrite, and injected rollback tests pass.
- [x] Slice 5 — explicit indexed online assembly and canonical output integration. Indexed routing, incompatibility, all-skipped diagnostic-only, canonical output, control/query/count/overlap, and provenance tests pass.
- [x] Slice 6 — chromosome-22 empirical gate, measured tolerance/resources, and documentation completion. Both final-code baseline suites produced 838 atoms and 614,523 operator nonzeros, selected the same 489 individuals, and reloaded under strict validation; canonical LD scores were exactly equal, maximum continuous overlap difference was `7.654307410120964e-08`, and downstream `partitioned-h2` was exactly equal. Final validation: 34 focused index tests passed; the affected checkpoint passed 428 tests plus 22 subtests; full pytest passed 1,176 tests plus 85 subtests with one environment-dependent skip; unittest passed 977 tests with one skip; all four required help commands and `git diff --check` passed.

## Goal and success signal

Implement the closed design in `docs/specs/2026-08-03-exact-disjoint-atom-gene-ldscore-index-design.md`: an offline PLINK-backed `ldsc build-gene-ldscore-index` workflow and an explicit indexed adapter on `ldsc ldscore` that assemble exact gene-list LD scores from disjoint genomic atoms.

The implementation succeeds when a matching direct and indexed run has structurally identical binary annotations, SNP rows and ordering, baseline/query grouping, annotation counts, overlap blocks, control semantics, statuses, and canonical output schemas; LD scores and downstream `partitioned-h2` results must agree within a tolerance measured by the required chromosome-22 validation rather than chosen in advance. A default single-worker chromosome build must remain within the approximately 4–8 GB target through SNP batching and internal atom-column batching.

## Current baseline and governing constraints

- The `restructure` branch already implements the Design 1 gene-list resolver, direct interval projection, query diagnostics, canonical split LD-score outputs, overlap-aware `partitioned-h2`, corrected regression-only named-region exclusions, and the PLINK LD-score path. Preserve these as the direct scientific oracle.
- The Design 2 specification is authoritative where older repository guidance still names obsolete LD-score files or describes the former SNP-region behavior.
- Keep workflow responsibilities in public `src/ldsc/` modules and numerical primitives in `src/ldsc/_kernel/`. `src/ldsc/cli.py` only registers and dispatches commands.
- Preserve the corrected SNP-universe split: the broad retained PLINK universe, including MHC and centromeric SNPs, contributes to baseline/query LD scores, counts, and overlaps; bundled HM3 minus `mhc-and-centromeres` defines persisted rows and the identical `w_ld` contributor set.
- Reuse the existing PLINK genotype decoding, adjusted-\(r^2\), cM-window, individual filtering, inclusive MAF, genetic-map, HM3, region, query-status, overlap, output-preflight, and workflow-logging seams. Do not create a numerically separate LD implementation for the index.
- Accumulate sparse operator values in float64, preserve negative adjusted-\(r^2\), include the diagonal exactly once, and do not clamp or epsilon-prune. Cast only at the existing canonical Parquet output boundary.
- User-requested focal columns are always assembled together. Only offline atom columns and existing SNP work may be batched.
- Indexes are separate distribution artifacts, not package data. There is no generic index schema/software version and no redundant checksum for every payload; semantic IDs, structural validation, staged publication, and optional archive transport checksums provide integrity.
- The initial supported profile domain is hg19, rsID identity, the `1000G_EUR_Phase3` PLINK suite, 1 cM, bundled HM3 minus `mhc-and-centromeres`, no explicit retained-reference MAF threshold, `common_maf_min=0.05`, no individual keep file, informative BIM cM, 100 kb padding, and MHC gene exclusion. Both approved baseline suites must pass the chromosome-22 gate.
- Use focused red-green-refactor cycles within each slice. The plan is a living guide: record completed slices and revise implementation details if evidence changes, without changing scientific or artifact contracts unless the user approves a specification amendment.

## Interfaces and invariants

### Public command and Python seams

- Add `ldsc build-gene-ldscore-index` with a workflow entrypoint in a new `src/ldsc/gene_ldscore_index.py` module and a dedicated validated config record, preferably in `src/ldsc/config.py` with the other workflow configs.
- Add explicit indexed mode to `ldsc ldscore` through `--gene-ldscore-index-dir <profile-dir>`. It must never discover an index, silently fall back, or accept live scientific settings owned by the profile.
- Rename `--bed-padding-bp` / `bed_padding_bp` to `--padding-bp` / `padding_bp` across `annotate`, direct/indexed `ldscore`, the builder, Python entrypoints, configuration snapshots, tests, and current docs. The old name must be rejected; no public compatibility alias remains.
- Add `--gene-exclude-regions {none,mhc}` and singular `--control-gene-list-source`. Direct gene-list mode defaults to gene exclusion `none`; initial index profiles use `mhc`. The control defaults to reserved token `all-protein-coding`, also accepts `none` or a path, and is valid only for gene-list workflows.
- `gene_control` is a reserved fixed baseline column, not a query. It is appended after supplied baseline columns, participates in the baseline rows of the overlap artifact, and leaves the current baseline-plus-one-focal-query `partitioned-h2` loop unchanged.

### Exact atom model

For one chromosome, the implementation must preserve the specification's identity

$$
z(a)=\mathbf{1}[B^\mathsf{T}a>0],\qquad q(a)=Hz(a),\qquad \ell(a)=Yz(a),\qquad Y=PRH.
$$

`B` is Boolean gene-to-atom CSR, each retained reference SNP belongs to at most one atom in `H`, and selecting multiple overlapping genes is a Boolean union. Duplicate identifiers, aliases, nested genes, padding-created overlaps, and repeated genes must never double-count a SNP. Gene exclusion uses the unpadded transcribed interval and occurs before padding; excluded genes remain in the embedded catalog with `included=false` and `exclusion_reason=excluded_gene_region`.

### Artifact and loading boundary

- A suite owns immutable common inputs and one or more profile directories. Root/common metadata repeat `suite_id`; profile/chromosome metadata repeat `suite_id` and `profile_id` and declare chromosome coverage, dimensions, ordering, sparse format, and dtypes.
- Common chromosome payloads are `baseline_rows.parquet` and `baseline_statistics.npz`. Profile payloads are `gene_catalog.parquet`, `atoms.parquet`, `gene_to_atom.npz`, `ldscore_operator.npz`, and `atom_statistics.npz`, with the exact members, shapes, and dtypes in the specification.
- `baseline_rows.parquet` follows the canonical public baseline table's identity columns, `regression_ld_scores`, ordered baseline score columns, and dtypes, but remains an internal chromosome shard with its distinct filename.
- Loader validation must complete before canonical output publication. It verifies artifact type, semantic IDs, coverage, component presence, row keys/order, dimensions, CSR invariants, dtypes, and finite values without relying on per-payload checksums.
- Existing `common/` content is reusable only when the calculated `suite_id` matches. Publication stages and validates a new common layer/profile before replacement. `--overwrite` replaces only the targeted profile and cannot mutate a shared common layer under other profiles.
- Indexed online output remains an ordinary self-contained directory written through `LDScoreDirectoryWriter`; no output may depend on the index remaining installed.

## Delivery slices

### Slice 1 — Complete direct-mode semantics and migrate the padding vocabulary

**Outcome:** Direct gene-list LD score calculation implements every behavior that indexed mode must reproduce, and all public padding surfaces use the final name.

**Likely areas:** `src/ldsc/config.py`, `src/ldsc/annotation_builder.py`, `src/ldsc/gene_list_resolver.py`, `src/ldsc/query_annotations.py`, `src/ldsc/ldscore_calculator.py`, `src/ldsc/outputs.py`, `src/ldsc/overlap_matrix.py`, `src/ldsc/cli.py`, and their existing tests.

**Work:**

1. Use an expand–migrate–contract sequence for the wide padding rename: introduce the internal `padding_bp` vocabulary, migrate every caller, persisted field, test, and current document, then remove the old parser/config/Python keyword in the same slice. The checkpoint must prove both the new name works and the old name is an ordinary unrecognized/invalid argument.
2. Extend catalog projection with a pure pre-padding MHC gene filter driven by the packaged hg19/hg38 MHC interval. Preserve excluded catalog identities and add `excluded_gene_region` counts/audit rows without conflating this with SNP `--exclude-regions`.
3. Add control-source normalization and resolution. The `all-protein-coding` token selects all profile-eligible catalog genes, `none` disables the control, and a path uses the same resolver, build, gene exclusion, and padding as focal lists. A partially resolved custom control warns and retains usable genes; an unusable control aborts the run.
4. Project the control through the same direct binary interval-union primitive as focal queries, append `gene_control` to the supplied baseline matrix, and carry its counts and overlap rows through the existing result and writer contracts. Reject collisions with a supplied baseline or focal query of that name.
5. Preserve focal query ordering, status/pruning behavior, gene-level diagnostics, broad-reference counts, regression-row filtering, and canonical output semantics. Add concise control and gene-region provenance without storing absolute runtime dependencies.

**Validation checkpoint:**

- Unit tests cover MHC overlap at interval boundaries, filtering before padding, HLA genes as ordinary MHC-overlapping genes, non-MHC chromosome-6 genes, and no centromeric gene filter.
- Equivalent gene-list and hand-built BED controls/queries produce identical direct annotations, LD scores, counts, and overlap values.
- `all-protein-coding`, `none`, valid/partial/unusable custom controls, and `gene_control` collisions have exact behavior.
- Existing chromosome-6 tests continue to prove that MHC/centromere SNPs contribute to LD scores/counts while remaining absent from persisted rows and `w_ld`.
- CLI, Python, metadata, and documentation searches contain no remaining public `bed_padding_bp` contract.

### Slice 2 — Implement and prove the pure disjoint-atom model

**Outcome:** A private numerical module can deterministically construct atoms and assemble annotations/statistics from catalog selectors without PLINK or filesystem concerns.

**Likely areas:** add `src/ldsc/_kernel/gene_ldscore_index.py`; add a focused `tests/test_gene_ldscore_index_kernel.py`; reuse catalog and region records from the workflow layer only through primitive arrays/frames passed into the kernel.

**Work:**

1. Define compact internal records for chromosome atom geometry, catalog row ordering, Boolean `gene_to_atom` CSR, and atom-statistic arrays. Keep public file and gene identifiers out of the numerical functions.
2. Build maximal nonempty half-open atoms from included padded gene intervals using deterministic genomic ordering. Permit omission of SNP-empty atoms only after proving that SNP-level annotations and diagnostics do not change.
3. Map retained SNP positions to at most one atom without materializing a dense SNP-by-atom matrix. Provide an iterator that emits dense Boolean SNP-by-atom blocks only for the builder's internal atom batch.
4. Assemble a focal/control gene selector into `z` by Boolean OR, then assemble `q`, atom counts, baseline-atom overlap vectors, control-query intersections, and query self-overlaps. Keep all count universes explicit rather than inferring them from output rows.
5. Add structural validators for sorted/nonoverlapping atoms, catalog and CSR ordering, excluded-gene empty rows, shapes, dtypes, and Boolean membership.

**Validation checkpoint:**

- Synthetic tests cover individual, overlapping, adjacent, nested, duplicate/alias-selected, padded, chromosome-start-clipped, excluded, and SNP-empty genes.
- For every synthetic selector, `H @ z` exactly equals the direct interval-union annotation and never exceeds one.
- All/common atom counts and baseline-atom overlaps reproduce direct dense calculations, including continuous baseline columns.
- Empty chromosomes and chromosomes with genes but no retained atom SNPs produce valid zero-dimensional records.

### Slice 3 — Build exact common and profile chromosome data from PLINK

**Outcome:** The new builder workflow produces validated in-memory common/profile chromosome records using the same reference, window, and adjusted-\(r^2\) behavior as direct PLINK `ldscore`, within the memory target.

**Likely areas:** add `src/ldsc/gene_ldscore_index.py`; extend `src/ldsc/config.py`, `src/ldsc/cli.py`, `src/ldsc/_kernel/ldscore.py`, and possibly `src/ldsc/_kernel/plink_bed.py` only to expose a reusable prepared-PLINK seam; add `tests/test_gene_ldscore_index.py` and focused PLINK fixture tests.

**Work:**

1. Register `build-gene-ldscore-index` in both lightweight and full CLI paths, add its public workflow/config entrypoints, and enforce the v1 compatibility matrix before heavy work: PLINK only, one cM window, fixed bundled HM3 and `mhc-and-centromeres`, no R2/custom SNP universes/alternate windows/whole-chromosome override/reference metadata export.
2. Resolve baseline annotation suites, PLINK chromosome prefixes, optional keep individuals, optional inclusive `maf_min`, `common_maf_min`, explicit matching-build genetic maps, chromosome coverage, and catalog/profile inputs through existing helpers. An explicit map overrides BIM cM; otherwise BIM cM must be informative.
3. Before genotype-derived filtering, compare canonical genomic rows across every ordered baseline source and the PLINK BIM. Require exact `CHR/POS/SNP` identity for the initial allele-free rsID suites after canonical sorting; reject missing, extra, duplicate, or conflicting rows with actionable diagnostics. Do not reuse direct mode's permissive intersection at this boundary.
4. Prepare the retained PLINK chromosome once after individual selection, monomorphic/unusable genotype removal, and optional MAF filtering. Derive the broad contributor universe, common mask, fixed regression/output selector, and `w_ld` mask as separate named arrays and record their counts.
5. Compute common baseline LD scores, `w_ld`, counts, and overlap statistics through the existing PLINK adjusted-\(r^2\) path. Compute `Y` in internal atom-column batches by reusing that same path, selecting only fixed regression rows from each float64 result block, converting exact nonzeros to CSR, and releasing dense blocks promptly. Never allocate the full SNP-by-atom matrix or dense full `Y`.
6. Keep the default worker count at one. Make chromosome parallelism opt-in through `--threads`, log its multiplicative memory consequence, and record effective SNP/atom batch sizes and chromosome dimensions.

**Validation checkpoint:**

- A deterministic small PLINK fixture proves strict pre-QC identity equality, reordering allowance, genotype-derived filtering, keep-file behavior, inclusive optional MAF, map precedence, and uninformative-CM errors.
- Pair-level numerical tests prove window boundaries, one diagonal, symmetric contributions, negative adjusted-\(r^2\), and contribution from non-HM3 and excluded-region reference SNPs to retained HM3 rows.
- Batched and single-block atom construction agree in float64; varying only the internal atom batch does not change structural output and differs numerically no more than the underlying accumulation order permits.
- Memory-oriented tests or instrumentation prove no dense whole-chromosome SNP-by-atom allocation; a focused benchmark records peak RSS for the largest locally available chromosome fixture/prototype.

### Slice 4 — Serialize, validate, and transactionally publish suite/profile artifacts

**Outcome:** The builder writes the approved directory tree, can safely reuse a matching common layer, and the loader rejects incomplete or incompatible profiles before scientific outputs are exposed.

**Likely areas:** module-local index reader/writer and identity helpers in `src/ldsc/gene_ldscore_index.py` unless implementation evidence justifies a separate workflow I/O module; reuse atomic tree-replacement patterns from `src/ldsc/outputs.py`; extend output, logging, and package-layout tests.

**Work:**

1. Define canonical JSON normalization and content-identity calculation for `suite_id` and `profile_id`. Include every specified ordered source identity and scientific setting, selected-individual identity/count, chromosome coverage, removal counts, and decompressed catalog SHA; exclude diagnostics and human-readable directory names.
2. Serialize common and profile metadata plus Parquet/NPZ payloads with the exact names, members, shapes, ordering, and dtypes in the specification. Use standard SciPy CSR NPZ representation with float64/int32 `Y` and Boolean gene-to-atom data; do not add generic version fields or payload checksums.
3. Implement a reader that validates the whole selected profile and common layer as one compatibility unit. Treat missing metadata, wrong IDs, undeclared chromosomes, shape/dtype/order disagreement, invalid CSR, nonfinite values, or row-identity mismatches as corruption/incompatibility.
4. Stage new payloads under the destination filesystem, reload and validate the staged tree, then publish it with backup/rollback replacement semantics consistent with existing whole-tree writers. On failure, retain an actionable profile log but no loadable profile.
5. Reuse existing common data only after recalculating and matching `suite_id`. Refuse a conflicting common layer even with `--overwrite`; overwrite may replace only the requested profile. Do not maintain a mutable root profile list.
6. Write `profiles/<profile>/diagnostics/build-gene-ldscore-index.log` with input resolution, strict alignment, filtering, IDs, common creation/reuse, per-chromosome dimensions, batches/resources, validation, publication, warnings, and elapsed time.

**Validation checkpoint:**

- Round-trip tests inspect every JSON/Parquet/NPZ field and prove deterministic IDs across equivalent path spellings and different IDs after any scientific input changes.
- Corruption tests independently alter IDs, coverage, rows, shapes, dtypes, CSR indices, and required members and confirm failure before publication/use.
- Failure-injection tests around common/profile writes prove that an old valid profile survives a failed overwrite and that a failed new build cannot be discovered as valid.
- Common reuse, conflicting-suite refusal, targeted-profile overwrite, subset-chromosome identity, and sibling-profile preservation are covered.

### Slice 5 — Add explicit indexed online assembly and canonical output integration

**Outcome:** `ldsc ldscore --gene-ldscore-index-dir <profile>` resolves focal/control lists against the embedded catalog, assembles all requested columns exactly, and writes an ordinary canonical LD-score directory without source PLINK/baseline files.

**Likely areas:** `src/ldsc/gene_ldscore_index.py`, `src/ldsc/ldscore_calculator.py`, `src/ldsc/gene_list_resolver.py`, `src/ldsc/outputs.py`, `src/ldsc/overlap_matrix.py`, `src/ldsc/regression_runner.py`, and their workflow/output/regression tests.

**Work:**

1. Add a strict indexed-mode argument branch. Accept only the profile path, focal sources, control source, output/overwrite/logging/threading/diagnostic controls; reject every live scientific option listed in the specification with guidance to use direct mode. Validate the full index before preflighting canonical scientific output.
2. Let `GeneCatalog` load the embedded `gene_catalog.parquet` as the sole online authority while reusing the exact Design 1 identifier resolver and diagnostics. Apply the profile's stored inclusion policy; no installed-catalog discovery, build inference, padding, or region projection occurs online.
3. For each chromosome, assemble every surviving focal selector and the optional control together: Boolean-union gene rows into atom selectors, multiply the float64 CSR operator by the selector matrix, compute all/common counts and overlap blocks from stored sufficient statistics, and cast through the canonical writer only. Do not loop/batch user query columns.
4. Combine stored common baseline rows with the assembled control and focal columns in the existing `LDScoreResult`/`LDScoreOverlap` structures. Preserve source order, zero-hit/zero-variance pruning, count-record order, overlap labels, statuses, all-skipped diagnostics, one row group per chromosome, and downstream baseline-plus-one-query semantics.
5. Extend result/writer metadata additively with `suite_id`, `profile_id`, profile catalog provenance, control provenance, and index scientific policy. Use `dataclasses.replace(...)` when enriching results so overlap/status/provenance fields cannot be dropped.
6. Keep explicit-index failures fail-closed: no discovery, direct fallback, partial thin result, or canonical publication after a corrupt/incompatible component.

**Validation checkpoint:**

- Synthetic end-to-end direct/index tests compare annotations, retained rows/order, baseline/query columns, `gene_control`, counts, all/common overlap values, statuses, diagnostics, metadata, and output schemas exactly.
- Indexed tests cover overlapping/nested/duplicate/alias genes, partial resolution, MHC-filtered genes, zero-hit/zero-variance queries, custom/default/disabled controls, subset chromosome profiles, and all-skipped batches.
- The argument compatibility matrix has positive and negative CLI/Python tests, including explicit proof that missing/corrupt profiles never fall back.
- Matching canonical direct/index directories feed unchanged `load_ldscore_from_dir` and `partitioned-h2`; all focal columns are loaded together and each fitted model still contains supplied baseline plus `gene_control` plus one focal query.

### Slice 6 — Run the empirical gate, document the shipped workflow, and complete regression validation

**Outcome:** Numerical tolerance and resource behavior are supported by recorded evidence, current documentation matches the implementation, and production whole-genome index construction is unblocked but remains a separate distribution operation.

**Likely areas:** add a dated validation report under `docs/audits/`; update the governing specification with the measured tolerance; update affected `docs/current/`, `README.md`, `docs/wiki/guided-tutorial.md`, and `docs/troubleshooting.md`; extend end-to-end and package-layout tests.

**Work:**

1. Build a chromosome-22 prototype for each approved baseline suite using the local `1000G_EUR_Phase3` PLINK chromosome-22 files. Confirm the strict 141,123-row baseline/BIM identity prerequisite and record resolved inputs and IDs.
2. Compare representative direct and indexed runs containing overlapping/nested genes, MHC-adjacent/excluded genes, aliases/duplicates, and non-HM3 contributors. Measure differences after canonical float32 output casting for LD scores and downstream `partitioned-h2` coefficients, standard errors, enrichment, p-values, and ordering. Record the observed maxima and establish the acceptance tolerance; do not loosen it to hide an unexplained discrepancy.
3. Record build wall time, peak RSS, atom count, `nnz(Y)`, payload sizes, and batch settings for chromosome 22. Run/record chromosome-6 correctness and resource evidence when its full PLINK shard is available on HPC; keep deterministic local chromosome-6 tests mandatory meanwhile. Stop before production whole-genome builds if exactness or the 4–8 GB single-worker target fails.
4. Document the strict `baseline identities == PLINK identities` requirement, supported v1 configuration, suite/profile selection, direct fallback domain, immutable scientific identity, common/profile tree, overwrite behavior, control model, gene versus SNP region exclusions, no-query-batching decision, internal atom batching, and separately distributed index policy.
5. Update the partitioned-h2 manual with the approved complete-query-Parquet memory table, its HM3-row/float32 assumptions, and the recommendation for approximately 1.5–2 times table-size headroom. Make clear that more memory is the supported path for very wide focal sets.
6. Audit current architecture/data-flow/configuration/argument/artifact metadata docs and CLI help for the new workflow and public rename. Do not retrofit obsolete generic version or old output-layout language into the new index contract.

**Validation checkpoint:**

- The chromosome-22 report is reproducible and linked from the specification; production full-genome construction remains blocked until the report passes.
- Focused kernel, builder, loader, direct/index equivalence, output, regression, CLI, logging, region, config, and package-layout tests pass.
- The full `pytest` suite and standard-library `unittest` compatibility run pass in `ldsc3-dev`.
- Both `ldsc --help` and `python -m ldsc --help`, plus command-specific help, expose the final command/flags only.
- `git diff --check`, package-resource import checks, and documentation path/link checks pass.

## Suggested validation commands

Use the repository development environment:

```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh
conda activate ldsc3-dev

pytest -q tests/test_gene_ldscore_index_kernel.py tests/test_gene_ldscore_index.py
pytest -q tests/test_gene_list_resolver.py tests/test_annotation.py tests/test_ldscore_workflow.py tests/test_overlap_matrix.py tests/test_output.py tests/test_regression_workflow.py tests/test_regions.py
pytest -q
python -m unittest discover -s tests -p 'test*.py' -v
ldsc --help
ldsc build-gene-ldscore-index --help
ldsc ldscore --help
python -m ldsc --help
git diff --check
```

Test filenames may be adjusted if implementation evidence favors a different focused split, but the kernel, builder/loader, online adapter, canonical writer, and regression coverage must remain independently runnable.

## Risks and revision checkpoints

- **Public rename blast radius:** the old padding name appears in configuration, internal kernel helpers, metadata, tests, and docs. Use the explicit expand–migrate–contract checkpoint and a final repository search; do not leave a public alias accidentally.
- **Direct-oracle drift:** implementing control/gene filtering separately in indexed mode would make equivalence ambiguous. Complete and test the shared direct semantics before index integration; reuse resolver/status definitions rather than duplicating them.
- **Strict identity versus current permissive intersection:** the builder must inspect pre-QC baseline and BIM identities before invoking ordinary PLINK filtering. If current annotation cleanup hides duplicates or missing rows, add a builder-specific strict boundary rather than weakening the requirement.
- **Sparse numerical cancellation:** adjusted-\(r^2\) includes negative values, so atom batches and CSR conversion can expose accumulation-order differences. Retain float64 through construction/assembly and use the chromosome-22 gate to distinguish ordinary final-cast noise from algorithmic error.
- **Sparse payload size:** `Y` may be less sparse than genomic intuition suggests because negative adjusted-\(r^2\) values are retained. Record `nnz(Y)` and bytes early; if the prototype is impractical, pause for a design review rather than prune values or change the scientific operator.
- **Transactional shared state:** adding or replacing a profile beside shared common data has a larger failure surface than a flat output family. Prove rollback and sibling preservation with failure injection before allowing overwrite on real suites.
- **Memory multiplication:** `--threads` multiplies PLINK/genotype, baseline, dense atom-batch, and sparse-builder working sets. Benchmark one worker first, log the effective worker count, and avoid an automatic all-core default.
- **Catalog authority:** direct mode uses the installed catalog while indexed mode uses the embedded catalog. Their content SHA, release, build, inclusion flags, and ordering must match for equivalence; a mismatch is incompatibility, not a warning.
- **Evidence availability:** only chromosome 22 of the approved PLINK suite is local. Record local chr22 evidence first and treat full chromosome-6/full-genome measurements as an HPC validation checkpoint, not a reason to substitute R2 Parquet or relax acceptance.

Any result that suggests changing the exact operator, SNP universes, regression rows, baseline equality rule, public output family, control model, or artifact identity requires a specification amendment before implementation proceeds.

## Out of scope

- R2-Parquet index construction;
- custom indexed reference/regression SNP lists or SNP exclusion policies other than fixed `mhc-and-centromeres`;
- SNP-count, kb, or whole-chromosome indexed LD windows;
- automatic index discovery, compatibility search, or silent direct fallback;
- query-column batching, a public query-batch flag, or a query-count cap;
- concurrent profile append, incremental chromosome finalization, or mutation of a profile's coverage;
- packaging or installing the built indexes with the Python distribution;
- thin canonical outputs that refer to index payloads;
- per-payload checksums, generic schema/software version fields, or a mutable root profile registry;
- approximate gene windows, summing overlapping per-gene LD scores, clamping, or epsilon sparsification;
- production construction and publication of both whole-genome distributed suites before the empirical gate passes.

## Completion definition

The feature is complete when every acceptance criterion in the governing specification has an automated test or a recorded empirical check, the chromosome-22 gate establishes and documents the numerical tolerance for both approved baseline suites, the single-worker resource target is evidenced, current docs describe only the final interfaces, and matching direct/index outputs remain interchangeable for downstream `partitioned-h2`. Passing unit tests for atomization alone, or producing a loadable index without direct/downstream equivalence, is not sufficient.
