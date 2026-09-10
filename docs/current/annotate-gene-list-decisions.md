# Standalone Annotate Gene-list Decisions

Last updated on: 2026-09-10

Status: standalone BED/gene-list CLI and `run_annotate()` are implemented using staged chromosome access. Direct/indexed workflow migration continues under the [memory implementation plan](../plans/2026-09-10-annotation-workflow-memory.md). The decision table remains the behavior contract; implementation evidence is recorded in the [memory audit](../audits/annotation-memory/progress.md).

## Problem and scope

Standalone `ldsc annotate` accepts BED or focal gene-list queries over the explicitly supplied baseline annotation SNP grid. The implementation is local to `ldsc_py3_restructured`; no HPC execution is part of this refactor.

The inspected implementation is commit `a505c45` (`refactor: remove annotation fingerprints`). Relevant seams are [annotation_builder.py](../../src/ldsc/annotation_builder.py), `add_annotate_arguments()`, `run_annotate_from_args()`, `AnnotationBuilder.run()`, and `_project_intervals_to_metadata()`; [ldscore_calculator.py](../../src/ldsc/ldscore_calculator.py), `build_parser()` and `_normalize_run_args()`; and [gene_list_resolver.py](../../src/ldsc/gene_list_resolver.py), `GeneCatalog.load()` and `resolve_gene_lists()`.

## CLI and Python contract

Accept exactly one of `--query-annot-bed-sources` and `--query-annot-gene-list-sources`. Continue requiring `--baseline-annot-sources` and `--output-dir`. Preserve existing BED behavior, SNP identities and row ordering, duplicate handling, and annotation values.

| Argument | Approved standalone behavior |
| --- | --- |
| `--query-annot-bed-sources` | Existing BED route; mutually exclusive with gene-list sources. |
| `--query-annot-gene-list-sources` | Focal gene-list exact paths or globs, using existing path-token normalization, deterministic expansion, source naming, and resolution rules. |
| `--gene-coordinate-file` | Required for gene lists; the user-supplied TSV/TSV.GZ is the sole coordinate and gene-resolution authority. |
| `--padding-bp` | Gene lists require an explicit nonnegative integer, including explicit `0` for gene bodies. BED omission remains effective zero. Preserve omission until the route is known. |
| `--gene-list-resolution-policy` | Gene-only; `strict` by default or explicit `resolved-only`, using the existing omission allowlist. |
| `--gene-exclude-regions` | Gene-only; `none` by default or explicit `mhc`. Exclude genes using their unpadded intervals before applying padding. |
| `--snp-identifier`, `--genome-build` | Preserve identity modes and reuse applicable direct-ldscore build-resolution rules described below. |
| `--overwrite`, `--log-level` | Existing output-family and logging controls. |

Do not expose `--control-gene-list-file`. A control becomes baseline column `gene_control` in direct ldscore, whereas annotate writes query columns. No control or copied-baseline output suite is added. Existing shared control support remains available to its current workflows.

Do not add prebuilt-query or indexed query routes, synthetic baseline generation, PLINK/R² inputs, reference-panel filters, regression-SNP selectors, LD windows, or regression-specific exclusions. The baseline annotation grid remains required; a BIM-based replacement is outside this decision.

Use a general public `run_annotate()` entry point and rename the parser helper to `parse_annotate_args()`. Keep `run_annotate_from_args()` as the CLI dispatch seam. Replace obsolete BED-only wrappers and update exports, internal callers, tests, and examples without compatibility aliases. The return value follows the shared shard-based `AnnotationBundle` interface and references persistent outputs; its precise resource interface belongs to the memory refactor.

## Catalog, build, and projection semantics

Reuse the existing [gene-list input contract](gene-list-input-format.md#gene-list-source-format): exact case-sensitive identifiers or names, no implicit alias repair or Ensembl-version stripping, one identifier per nonblank row, and within-source deduplication by canonical gene ID. Preserve source and row audit ordering and established annotation-name collision rules.

The catalog contains one declared build and one-based inclusive intervals. Convert `[start, end]` once to internal `[start - 1, end)`, then pad both sides and clip the start at zero. The annotation is the binary union of selected intervals on surviving baseline SNP rows. Explicit MHC exclusion is an audited transformation, not rejected resolution; there is no implicit gene exclusion.

In gene-list mode, omitted `--genome-build` behaves as `auto`. Combine catalog and available baseline evidence through existing inference helpers, with no reference-panel dependency. Conflicting or insufficient evidence fails according to the direct workflow's applicable rules; an explicit build must agree with the catalog. No implicit liftover occurs. For rsID-family identity, retain build-independent identity configuration and separately record the concrete projection build. Preserve BED build behavior when no catalog is involved. See [ldscore_calculator.py](../../src/ldsc/ldscore_calculator.py), `_normalize_run_args()` and `_resolve_ldscore_chr_pos_genome_build()`, and [genome_build_inference.py](../../src/ldsc/genome_build_inference.py).

## Validation gates and condition–outcome table

The table applies to the new standalone gene-list route. Existing BED acceptance and scientific behavior remain unchanged.

Gate A validates the catalog and resolves all focal sources together before substantial annotation projection. Collect every safely discoverable source, row, catalog, and naming issue at the gate where its prerequisites are available. A structurally unusable catalog may prevent reliable matching; report its issues without inventing a gene-resolution audit. An unreadable list receives a source-summary error with unknown row counts. Unreferenced latent catalog defects are diagnosed without failing an otherwise valid run. Reuse the current resolver's `RESOLVED_ONLY_REASONS` rather than creating another omission policy.

After Gate A, establish chromosome coverage from validated baseline contents before interpreting absent support. Exact paths and globs select their actual files; file contents establish scope. An `@` declaration requires all autosomes 1–22. Collect safely discoverable missing, unreadable, invalid, and inconsistent baseline inputs. Every selected gene must be inside the validated scope after resolution, deduplication, and explicit gene exclusion. A pathway need not contain genes on every covered chromosome. Do not import direct ldscore's baseline/reference-set equality check: standalone annotate has no reference panel. Preserve the documented glob-selection caveat from [path specification](path-specification.md#direct-ld-score-query-scope).

Gate B measures support on baseline SNP rows remaining after existing annotation identity cleanup. Count surviving SNP rows intersecting each padded gene interval; query counts measure the interval union without double-counting overlap. Determine query viability across the complete validated chromosome scope. Unknown or unevaluated support is blank, never zero.

| Condition | Outcome |
| --- | --- |
| All identifiers resolve and required inputs validate | Assess coverage and then annotation-grid support. |
| Empty readable focal list | Skip as `empty_gene_list` under either policy; usable siblings may continue. Strict policy rejects invalid identifiers, not blank files. |
| Rejected identifiers under `strict` | Fail the complete batch at Gate A after collecting available issues and writing diagnostics. |
| Permitted identifier rejections under `resolved-only` | Omit only allowlisted rejected rows; use the audited resolved subset. Retained partial queries carry `partial_gene_resolution`. |
| Nonempty focal list has no selected genes after allowed omissions or explicit exclusions | Skip as `zero_resolved_genes`; under `strict`, rejected rows would already have failed Gate A. |
| Unreadable or incorrectly compressed/encoded list, malformed multi-field row, or duplicate query name | Fail under either policy, collecting independent issues at the applicable gate. |
| Query name collides with a baseline column or established reserved name | Fail at the first gate with enough information to detect the collision; no implicit renaming. |
| Structurally unusable catalog or incompatible build | Fail with available catalog/build diagnostics; do not fabricate resolution or support results. |
| Unreferenced latent catalog defect | Diagnose without failing an otherwise valid run; referenced defects follow existing resolution policy. |
| Repeated rows or aliases select the same canonical gene | Select that gene once within the source and retain duplicate audit rows. |
| Explicit gene-region exclusion removes a gene | Record `excluded`; assess coverage and support only for the remaining selected genes. |
| Selected gene lies outside validated baseline chromosome scope | Fail the complete batch under either policy as incomplete chromosome coverage; do not truncate or skip the pathway. |
| Required baseline input is missing, unreadable, invalid, or violates its declared scope | Fail the input gate after collecting safely discoverable issues; support remains unevaluated. |
| Individual resolved gene overlaps no surviving annotation-grid SNPs | Diagnose measured zero support; this alone is nonfatal. A supported sibling gene can keep the focal query usable. |
| Entire focal query overlaps no surviving annotation-grid SNPs | Skip globally as `zero_annotation_snps`; retain usable focal siblings. |
| A usable query is zero on one chromosome but supported elsewhere | Retain the column in every emitted chromosome shard, including its zero values. |
| Every focal query is skipped | Fail after writing available diagnostics; do not publish a new scientific query family. |
| Annotation is all one | Retain it. Annotation generation has no LD-score variance or regression-design test. |

Preserve the shared precedence in [query_annotations.py](../../src/ldsc/query_annotations.py), `gene_query_statuses()`: incomplete coverage aborts before support; zero-query support causes a skip; partial resolution remains the primary warning when partial support also occurs. All-focal-skipped handling is shared policy. Reference preparation, MAF/sample/SNP restrictions, LD-score variance checks, and control viability are outside standalone validation.

## Outputs, diagnostics, and failure behavior

Keep root `query.<chrom>.annot.gz` files with identical ordered query columns across chromosomes. Omit globally skipped queries from every shard. Preserve the canonical `CHR/BP/SNP/CM` layout, explicit `CM=NA`, existing allele-column behavior, absence of annotation MAF, and generated integer `0/1` values. Preserve reloadability through the shared annotation readers and downstream prebuilt-query ldscore route. See [annotation_builder.py](../../src/ldsc/annotation_builder.py), `_write_bundle_query_as_annot_files()`, and [lessons.md](../../lessons.md#missing-metadata-must-be-explicit-in-reusable-whitespace-parsed-tables).

Persist the row-complete `diagnostics/gene_list_audit.tsv.gz`, per-source `diagnostics/gene_list_resolution_summary.tsv`, and post-Gate-A `diagnostics/query_annotation_status.tsv`, plus applicable chromosome-scope and input/catalog-issue diagnostics. Retain the existing `diagnostics/metadata.json`, `diagnostics/annotate.log`, and dropped-SNP audit. Structural catalog failure may produce catalog issues without a list audit; Gate A rejection writes the available audit/summary without fabricating later-stage statuses. Reuse existing diagnostic schemas and writers where applicable.

Annotation-grid support must be explicit: use `annotation_snp_count` for measured per-gene grid counts and an annotate-specific zero-support reason, leaving reference-panel counts unevaluated. Summary support counts describe the annotation grid; `n_annotation_snps` counts marked surviving rows over all emitted chromosomes. Preserve reference-support semantics for direct/indexed ldscore. The support-update helper and shared status text must distinguish these universes rather than calling baseline-grid SNPs reference-panel SNPs.

Provenance records query and baseline sources, catalog source/build, concrete projection build, SNP identity configuration, padding, explicit exclusions, resolution policy, counts, chromosome scope, and applicable diagnostic paths. Do not add routine content hashes or restore removed annotation fingerprints.

Fatal validation publishes available diagnostics before new scientific outputs. Successful partial runs emit bounded CLI notices for deliberate resolution omissions and zero-support/skipped-query outcomes, with detailed diagnostics retained. Python callers receive structured statuses and paths without direct console notices. Follow [workflow logging](workflow-logging.md#output-family-preflight) for family preflight, logs, successful stale-owned-file cleanup, and authorized-overwrite failure markers. Do not add rollback, restoration, quarantine, or a new failure-time deletion policy. A failed overwrite is not promised to restore a prior successful family.

## Memory-refactor handoff and ownership

The [annotation memory decisions](annotation-memory-decisions.md#confirmed-dataset-contract) govern bundle storage and resource lifetime. The memory refactor owns the shared shard interface, bounded input scanning and preparation, private-storage ownership, incremental output machinery, and resolution of existing input-layout-dependent identity cleanup. The standalone feature owns the annotate arguments/API, applicable validation gates, diagnostics/provenance, and behavioral tests, using those shared facilities.

Resolve gene lists once and reuse the batch across chromosomes. Accumulate coverage/support/status summaries without whole-genome annotation matrices. Globally determine which query columns are retained before publishing a complete scientific family; the refactor determines the necessary bounded passes or private staging. Standalone annotate writes shards incrementally and returns an object referencing persistent outputs. Do not build a second storage framework, eagerly retain all chromosome DataFrames, reread a whole-genome source once per chromosome, or rematerialize whole-genome annotations on return.

The memory-refactor discussion subsequently approved dataset-wide annotation identity cleanup, resource-owning bundle handles with explicit closure and borrowing, and dependence on original baseline sources. It also requires new annotation/analysis staging inside the supplied output directory; remaining lifecycle details are tracked in [annotation memory decisions](annotation-memory-decisions.md). Those decisions govern this feature's integration, including the deliberate correction of layout-dependent duplicate handling. No portable copies of all baseline annotations are required. Standalone query-mode scope is settled and should not be reopened as an unresolved memory-refactor question.

## Implementation acceptance checks

The standalone implementation uses `annotate_workflow.run_annotate()`, shared chunk normalization, output-contained shard storage, staged gene resolution, and `AnnotationDirectoryWriter`. Affected parser/export tests now use `run_annotate()` and `parse_annotate_args()`. The complete documentation/tutorial migration and full-suite milestone remain part of the continuing memory-refactor plan.

Validate independently expected SNP membership, interval boundaries, overlap unions, padding, clipping, exact aliases and duplicates, and explicit exclusions. For example, an unpadded catalog interval `[101, 110]` must mark one-based SNP positions 101 and 110 but not 100 or 111, matching BED `[100, 110)`. Repeat with matching BED intervals and padding using independently derived expected values.

Cover chromosome-sharded and whole-genome baselines, identity and row-order preservation, globally supported queries with zero-only shards, every condition in the table, batched failures across independent sources, output reloadability, diagnostic/provenance contents, collision/overwrite handling, and stale-owned-artifact cleanup. Verify actual CLI and Python paths and that support does not consult a reference panel. Coordinate bounded-memory checks with the refactor, including discovery/preflight and return behavior. Run appropriate focused tests, the full pytest suite, CLI help, and the standard-library unittest compatibility check; run pytest and unittest sequentially because they share BED temporary-file cleanup. Report actual results after implementation; documentation-only verification does not establish runtime correctness.
