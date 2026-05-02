# Design Map

This file maps durable design documents to the implementation modules they
describe. Update it when workflow ownership, public entry points, or artifact
contracts move.

## Annotation Workflow

| Design document | Implementation |
| --- | --- |
| `docs/current/architecture.md` | `src/ldsc/annotation_builder.py`, `src/ldsc/_kernel/annotation.py`, `src/ldsc/cli.py` |
| `docs/current/code-structure.md` | public workflow boundary in `src/ldsc/annotation_builder.py`; private helper boundary in `src/ldsc/_kernel/annotation.py` |
| `docs/current/data-flow.md` | `AnnotationBuilder.run()`, `AnnotationBuilder.project_bed_annotations()`, `run_bed_to_annot()`, `run_annotate_from_args()`, `annotation_builder.main()` |
| `docs/current/io-argument-inventory.md` | annotate parser and namespace dispatch in `src/ldsc/annotation_builder.py`; top-level dispatch in `src/ldsc/cli.py` |
| `docs/current/layer-structure.md` | annotation row in the public workflow and kernel layer matrix |
| `docs/current/path-specification.md` | annotation path-token handling through `resolve_file_group()`, `split_cli_path_tokens()`, and output preflight helpers |

## Annotation Tutorials

| Tutorial | Covered implementation path |
| --- | --- |
| `tutorials/ld-score-calculation.md` | direct BED projection through `run_bed_to_annot()` and `ldsc annotate`; in-memory BED projection during `run_ldscore()` |
| `tutorials/partitioned-ldsc.md` | partitioned workflow using `AnnotationBuilder.run()` for in-memory query annotations and optional `ldsc annotate` materialization |
| `tutorials/cell-specific-ldsc.md` | cell-type BED annotations through `AnnotationBuilder.run()`, `ldsc ldscore`, and optional `run_bed_to_annot()` materialization |

## Current Boundary

- Public annotation workflow ownership is in `ldsc.annotation_builder`.
- `ldsc._kernel.annotation` contains only low-level table and BED helpers.
- `ldsc.cli` registers annotate flags from `annotation_builder` and dispatches
  parsed namespaces to `run_annotate_from_args()` without reparsing.
- `main_bed_to_annot()` is intentionally removed; `annotation_builder.main()`
  is the supported parser entry point.
