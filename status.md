# Handoff: annotation grids and LD-score annotation coverage

Last updated on: 2026-10-28

Status: high-level design confirmed by the user; CLI details and any optional intersection policy remain unresolved.
This handoff records the next design task, not an implemented feature or a completed specification.
The current request authorizes this handoff and a continuation prompt, not implementation.

## Completed: domain subpackage reorganization of `src/ldsc`

- The reorganization landed on `restructure` on 2026-10-28, from `9236b80` (dead-code cleanup) to the documentation commit `docs: describe the domain layout`.
  [Code structure](docs/dev/code-structure.md) and [architecture](docs/dev/architecture.md) describe the layout, and `tests/test_layering.py` enforces its import tiers.
- The module paths in this handoff follow the new layout.
- Follow-up cleanup landed on 2026-10-28 in the `refactor` commits from `cabe537` to `5b3cede`:
  - `window_coordinates` and `get_block_lefts` in `src/ldsc/_refpanel/windows.py` serve both `ldscore` and `build-r2-panel`; each command keeps its own missing-CM error.
  - `_runtime.artifacts.atomic_write_json` is the one atomic JSON writer, and every file it writes ends with a newline; `h2_scale` uses the shared `atomic_write_dataframe`.
  - The gene-index build no longer migrates `<index>.build/` logs or takes the standalone lock of earlier releases.
  - `_core.identifiers.open_text` is the one text-file opener; it decodes UTF-8 and treats `.gz` in any letter case as gzip.
  - Four dead or test-only helpers were removed, `_genome.liftover.load_hm3_curated_map(path)` became `load_hm3_map_file`, and three modules stopped re-exporting names they do not define.
  - The 24 underscore-prefixed names that other domains imported lost the underscore (`_runtime.preflight._readable` became `check_readable`), so an underscore again means module-local.
- Look-alike helpers deliberately kept, because each carries a domain-specific error message or return type: the four `_require_columns`, the two `_read_json`, and the two `get_compression`.
- Remaining follow-up, not started: consider per-domain test folders.

## Workspace and checkpoint

- Work locally in `ldsc_py3_restructured` on branch `restructure`.
- Read [AGENTS.md](AGENTS.md) before starting.
- At handoff creation, HEAD is `4e2fa3c` (`refactor: clarify LD-score batch publication`).
- The worktree was clean before this handoff was added.
- Recheck HEAD and the dirty worktree in the next session; preserve any subsequent user changes.
- An August session had an unrelated modification in `tests/test_snp_identity.py`; it is no longer dirty at this checkpoint.
- No source changes, commits, HPC runs, or external publication were performed in this session.

## Confirmed design from the user

These five decisions are settled and should not be reopened as unanswered preferences.
They supersede older baseline-only grid assumptions where those conflict with the new design.
Other existing scientific and artifact contracts remain applicable unless separately revised.

1. **Separate annotation definitions from the SNP grid.**
   Gene lists, together with the coordinate catalog and padding, define annotation columns.
   The SNP-grid input defines which SNPs receive those annotations.
2. **`ldsc annotate` accepts exactly one of two mutually exclusive SNP-grid inputs: PLINK `.bim` files or baseline annotation files.**
   Baseline annotation files are temporarily retained as an alternative grid source.
   Either input can explicitly define the candidate reference SNP universe.
   Supplying both as grid sources is invalid.
3. **Baseline annotations have two distinct roles.**
   They may temporarily supply the grid to `annotate`.
   Downstream, they supply the conditioning annotations for LD scoring and regression.
   Choosing BIM for the annotation grid does not eliminate downstream baseline covariates.
4. **Preserve the distinction between SNP universes.**
   The selected grid defines candidate SNPs.
   Reference-panel QC and explicit restrictions determine the effective computational universe.
   Regression SNPs form a separate subset.
   For parquet-R² workflows, reference SNP identities and coordinates come from the panel's metadata.
5. **Validate annotation coverage during `ldscore` reference preparation, before numerical LD-score computation.**
   Required annotation values must exist for the selected, retained computational SNPs.
   Missing annotations must not be filled with zero or silently shrink that universe.
   In baseline-grid mode, deliberately selecting baseline rows is an explicit universe choice.
   In BIM-grid mode, baseline coverage must not redefine the BIM-selected universe.

Exact CLI flag names and any optional intersection policy remain to be specified.
The parquet-R² metadata authority in decision 4 does not add a third standalone `annotate` grid-input mode.

### Consequences to preserve

- Distinguish a known annotation value of zero from a missing annotation row or value.
  A SNP outside a supplied gene/BED interval can legitimately receive zero during projection.
  An absent required annotation for a retained computational SNP cannot be imputed to zero.
- Validate coverage against the retained computational SNPs after applicable reference preparation and restrictions, before the numerical projection kernel runs.
  Do not equate every raw candidate SNP with a SNP retained after QC.
- Preserve the broader LD-contributor/count/overlap universe and the distinct regression/output/weight universe.
  Do not move regression-only filtering into contributor selection.
- The prior discussion of BIM as a complete replacement is superseded by the user's explicit temporary two-grid-mode design.
  No retirement date or removal mechanism for the baseline-grid alternative has been approved.

Concrete coverage example for future tests: a BIM-selected candidate grid contains SNPs A, B, and C, and all three remain after reference QC and restrictions.
If a required baseline annotation covers only A and B, LD-score preparation must fail rather than silently remove C or invent its annotation values.
If reference QC first removes C, coverage is assessed on the retained set instead.
Deliberately selecting baseline-grid rows A and B is a different, explicit candidate-universe choice.

## Current implementation versus the confirmed direction

The current standalone command already accepts BED or gene-list query definitions, but requires `--baseline-annot-sources` for the SNP grid.
It does not yet accept BIM files directly.
The [standalone decision record](docs/reference/annotate.md#inputs) explicitly excluded BIM replacement when it was written; that exclusion must not override the new user decision above.
The [annotate wiki](docs/user-guide/workflows/annotations.md#build-baseline-annotations-from-bed-files) currently documents an external BIM-to-minimal-annotation-grid recipe.

Existing standalone behavior to preserve unless a specific change is confirmed:

- BED and gene-list query definitions are mutually exclusive.
- The shared catalog/resolver provides exact case-sensitive matching, audited duplicates/exclusions, and strict versus explicitly permitted subset resolution.
- The current CLI uses the value-free `--allow-unresolved-genes` flag.
  The retired CLI `--gene-list-resolution-policy` has no alias; Python retains `gene_list_resolution_policy="strict"` or `"resolved-only"`.
- Gene-list padding must be explicit, including zero; BED omission means zero.
- Catalog intervals are one-based inclusive and convert once to `[start - 1, end)` before padding.
- Standalone annotate has no control-gene-list option.
  Direct/indexed LD-score workflows retain their existing optional control behavior.
- Standalone support is measured on the cleaned annotation grid and recorded as `annotation_snp_count`, separately from reference-panel support.
- Globally unsupported focal gene queries are skipped while usable siblings continue; all-skipped gene batches fail, and all-one annotations are retained.
- Root outputs are `query.<chrom>.annot.gz`; `CM=NA` and existing allele behavior preserve whitespace-reader compatibility.
- `run_annotate()` returns a resource-owning chromosome-backed `AnnotationBundle` with persistent outputs and structured paths/statuses.
  Callers use a context manager or `close()`; original baseline sources are currently dependencies for later reads.
- Keep bounded source scans, shared interval projection, output-contained scratch, and the existing writer ownership model.
- Failed authorized overwrites use `RUN_FAILED.txt`.
  The current contract does not promise rollback, restoration, or a new failure-time deletion policy.

The current [LD-score SNP-universe contract](docs/reference/ldscore-snp-universes.md#chromosome-coverage-versus-snp-support) describes a prepared baseline/reference intersection.
The next session must trace actual reference preparation and alignment before specifying how the new grid authority changes that path.
Do not assume the newly confirmed strict annotation-coverage behavior is already implemented.

## Remaining decisions and investigation

1. Choose exact CLI/Python names and mutually exclusive grid arguments.
   Define BIM path/suite resolution, supported input forms, and the representation of the temporary baseline-grid alternative.
   Do not silently require PLINK genotype companion files just to read a BIM grid; resolve any such proposed requirement explicitly.
2. Specify how selected grid identity and provenance reach `ldscore` when generated annotations are supplied as prebuilt inputs.
   The consumer must distinguish an explicit candidate-universe choice from accidental absence of baseline/query rows.
   No metadata schema or new downstream grid selector has yet been approved.
3. Specify exact identity, duplicate, ordering, allele, build-evidence, and chromosome-scope behavior for BIM grids using existing shared policies where applicable.
   Verify how those policies combine with PLINK reference rows and parquet-R² panel metadata.
4. Define the coverage gate's required annotation columns, diagnostic schema, batching, failure message, and position relative to QC/restrictions.
   Check baseline, prebuilt query, and live gene/BED paths as applicable.
   Determine whether coverage can be checked before expensive genotype/R² work without changing the selected computational universe.
5. Resolve whether any explicit optional intersection policy should exist.
   None is approved by this handoff; do not add automatic intersection, silent truncation, or a permissive fallback.
6. Identify affected downstream readers and documentation before finalizing the interface.
   Keep changes to gene-index behavior, numerical kernels, or regression models outside the default scope unless investigation establishes a necessary change and the user approves it.

Inspect repository facts before asking questions.
Resume a focused design discussion about these remaining details; obtain explicit implementation authorization before modifying source/tests.
The user has not asked this session to produce a full specification or implementation plan.

## Carryover findings from the previous read-only audit

The audit checked `restructure` at the same `4e2fa3c` checkpoint.
These are observed gaps, not additional approved feature scope; recheck them before deciding how to address them alongside the grid work.

| Finding | Reproduction and location |
| --- | --- |
| Missing successful partial-run console notices | A CLI run omitted one of two gene rows under `--allow-unresolved-genes` and exited zero without the promised LDSC terminal notice; another successful run skipped a zero-support query without a notice. Inspect [_annotation/annotate.py](src/ldsc/_annotation/annotate.py), `run_annotate()` final logging, against the standalone decision record's output policy. |
| Incomplete batched validation | `base.txt` containing a valid gene and an unknown identifier reported the identifier failure but missed its collision with baseline column `base`, because the early Gate A exit precedes naming checks. A missing sibling file also stops the input gate before readable sibling lists receive resolution audits. Inspect `run_annotate()` and [the shared input gate](src/ldsc/_runtime/preflight.py). |
| Incomplete strict error guidance | Twelve rejected rows were fully audited, with ten terminal examples, but `_gate_a_error()` omitted per-source rejected/submitted counts, an explicit omitted-example count, an audit filter, and overwrite guidance. |
| Stale workflow documentation | [data-flow.md](docs/dev/data-flow.md) still presents standalone annotate as BED-only. Main standalone decisions, wiki, and tutorials already describe the gene route. |

No fixes were made during that audit.
The August list of unresolved questions is historical: most were subsequently settled by the standalone and memory-refactor decisions.
Do not reopen control output, standalone LD-score variance checks, rollback, or obsolete Python wrappers solely because the August conversation considered them.

## Evidence and verification limits

The previous audit ran this focused command in the installed `ldsc3-dev` environment against this checkout:

```bash
/Users/wenbinwu/miniforge3/envs/ldsc3-dev/bin/python -m pytest -q tests/test_annotate_streaming.py tests/test_annotation.py tests/test_annotation_preparation.py tests/test_gene_list_resolver.py tests/test_cli_help.py
```

Result: **148 tests passed and 3 subtests passed in 15.30 seconds**.
Six additional real `python -m ldsc annotate` probes exercised partial resolution, partial/zero support, a missing file, strict rejection, a query-name collision, and all-one output.
The all-one probe retained `[1, 1, 1, 1]`; interval-boundary probes produced independently expected `[0, 1, 1, 0]` memberships.
The loaded package was verified as this repository's `src/ldsc`, using Python 3.13.13.
Full LD-score/index/regression suites and HPC workflows were not rerun.
These checks do not validate the newly confirmed BIM/coverage design.
No scientific tests were rerun merely to write this handoff.

Optional local probe evidence remains at `/private/tmp/ldsc-annotate-audit.whc642/evidence.json` on the user's Mac, with fixtures, commands, and outputs alongside it.
That temporary location may disappear; this handoff contains the consequential results needed for continuation.

## Reading map for the next session

Start with the current instructions, this file, and the relevant sections of:

- [Standalone annotate decisions](docs/reference/annotate.md).
- [LD-score SNP-universe contract](docs/reference/ldscore-snp-universes.md).
- [Architecture](docs/dev/architecture.md), especially annotation and reference/LD-score workflow ownership.
- [Annotation memory design](docs/dev/design/memory-and-cpu.md).
- [Gene-list input](docs/reference/formats/gene-lists.md), [diagnostics and repair](docs/reference/results/gene-list-diagnostics.md), [path specification](docs/reference/input-paths.md), and [workflow logging](docs/reference/logs-and-outputs.md).
- [Lessons](lessons.md), particularly reference preparation, physical BIM/BED row mapping, output ownership, and resource lifetime.

Implementation seams:

- [_annotation/builder.py](src/ldsc/_annotation/builder.py): public parser and CLI dispatch.
- [_annotation/annotate.py](src/ldsc/_annotation/annotate.py): standalone orchestration, builds, validation gates, and publication.
- [_annotation/sources.py](src/ldsc/_annotation/sources.py), [_annotation/parsing.py](src/ldsc/_annotation/parsing.py), and [_annotation/bundle.py](src/ldsc/_annotation/bundle.py): grid loading/cleanup, normalized values, and dataset access.
- [_annotation/interval_queries.py](src/ldsc/_annotation/interval_queries.py) and [_annotation/interval_projection.py](src/ldsc/_annotation/interval_projection.py): shared BED/gene interval projection.
- [_annotation/gene_lists.py](src/ldsc/_annotation/gene_lists.py), [_annotation/gene_query_storage.py](src/ldsc/_annotation/gene_query_storage.py), and [_annotation/query_outcomes.py](src/ldsc/_annotation/query_outcomes.py): shared gene identity, coverage, support, and query viability.
- [_ldscore/direct_annotation.py](src/ldsc/_ldscore/direct_annotation.py), [_ldscore/preflight.py](src/ldsc/_ldscore/preflight.py), [_ldscore/calculator.py](src/ldsc/_ldscore/calculator.py), and [_refpanel/adapters.py](src/ldsc/_refpanel/adapters.py): annotation/reference alignment, scope, retained reference preparation, and computation entry.
- [_annotation/outputs.py](src/ldsc/_annotation/outputs.py), [_runtime/paths.py](src/ldsc/_runtime/paths.py), and [_runtime/preflight.py](src/ldsc/_runtime/preflight.py): artifact ownership, path tokens, and staged input diagnostics.

Relevant existing tests include `tests/test_annotate_streaming.py`, `tests/test_annotation.py`, `tests/test_annotation_preparation.py`, `tests/test_gene_list_resolver.py`, `tests/test_ref_panel.py`, `tests/test_ldscore_workflow.py`, `tests/test_ldscore_chromosome_scope.py`, `tests/test_ldscore_projection.py`, `tests/test_staged_input_preflight.py`, and `tests/test_cli_help.py`.

## Ready-to-paste new-session prompt

```text
Continue the annotation-grid and LD-score coverage design in:
/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured

Use branch restructure. Read AGENTS.md and status.md first, then inspect the relevant current implementation and contracts.
The five decisions in status.md are finalized: annotation definitions are separate from the grid; annotate accepts exactly one of BIM files or temporarily retained baseline annotation files; downstream baseline covariates remain required where applicable; candidate, retained computational, and regression SNP universes remain distinct; required annotation coverage is validated during ldscore reference preparation before numerical computation, without zero-filling or silent universe shrinkage.

Treat those decisions as settled. Resolve only the remaining CLI/Python details, grid provenance, BIM identity/build/scope rules, coverage diagnostics/timing, and whether an explicit optional intersection policy should exist.
Do not assume an intersection option is approved. Recheck the carryover audit findings without automatically broadening implementation scope.
Preserve the shared gene-list policy, bounded-memory architecture, existing scientific invariants, and unrelated user changes.
Start with repository-backed findings and a focused round of remaining design questions. Do not implement until I explicitly authorize implementation.
```
