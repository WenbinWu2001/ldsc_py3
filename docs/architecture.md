# Architecture

`ldsc_py3_restructured` is the refactored Python 3 LDSC package. It reads SNP-level annotations, PLINK or parquet LD references, and GWAS summary statistics; resolves user-facing path and header conventions in the public workflow layer; delegates numerical work to `ldsc._kernel`; and writes LDSC-compatible artifacts that can be chained into later runs.

Related docs:

- [data_flow.md](data_flow.md): user-visible file streams and flowcharts
- [class-and-features.md](class-and-features.md): public API surface and major types
- [code_structure.md](code_structure.md): module map and change guide

![Package overview](assets/ldsc-package-overview.png)

## Bird's-Eye View

- **Build query annotations**: project BED or gene-set regions onto a baseline SNP grid. Entry points: `ldsc annotate`, `ldsc.AnnotationBuilder`
- **Build parquet reference panels**: convert PLINK genotype panels into standard parquet LD artifacts. Entry points: `ldsc build-ref-panel`, `ldsc.ReferencePanelBuilder`
- **Compute LD scores**: align annotations to a reference panel and emit LDSC-compatible LD-score artifacts. Entry points: `ldsc ldscore`, `ldsc.LDScoreCalculator`
- **Munge raw summary statistics**: normalize raw GWAS tables into curated `.sumstats.gz` artifacts. Entry points: `ldsc munge-sumstats`, `ldsc.SumstatsMunger`
- **Run LDSC regression**: consume munged sumstats and LD-score artifacts to estimate `h2`, partitioned `h2`, or `rg`. Entry points: `ldsc h2`, `ldsc partitioned-h2`, `ldsc rg`, `ldsc.RegressionRunner`

## Layer Structure

- **CLI Layer**: public command dispatch in `ldsc.cli`
- **Workflow And Preprocessing Layer**: public services in `ldsc.annotation_builder`, `ldsc.ref_panel_builder`, `ldsc.ldscore_calculator`, `ldsc.sumstats_munger`, and `ldsc.regression_runner`, plus shared normalization in `ldsc.config`, `ldsc.path_resolution`, `ldsc.column_inference`, and `ldsc.chromosome_inference`
- **Compute Kernel**: private file-format and numerical code in `ldsc._kernel.*`
- **Output Layer**: artifact writing in `ldsc.outputs` plus regression summary writers in `ldsc.regression_runner`

## File Tree

```text
ldsc_py3_restructured/
├── docs/
│   ├── architecture.md      # package-level architecture
│   ├── class-and-features.md
│   ├── code_structure.md
│   ├── data_flow.md         # workflow-level file streams
│   └── assets/
├── src/ldsc/
│   ├── __init__.py          # public Python exports
│   ├── __main__.py          # `python -m ldsc` entry point
│   ├── cli.py               # unified `ldsc` command
│   ├── config.py            # validated public config dataclasses
│   ├── path_resolution.py   # input-token normalization and discovery
│   ├── column_inference.py  # header and identifier normalization
│   ├── chromosome_inference.py
│   ├── annotation_builder.py
│   ├── ref_panel_builder.py
│   ├── ldscore_calculator.py
│   ├── sumstats_munger.py
│   ├── regression_runner.py
│   ├── outputs.py           # LDSC artifact writer
│   └── _kernel/
│       ├── annotation.py    # annotation parsing and BED projection
│       ├── ref_panel_builder.py
│       ├── ref_panel.py
│       ├── ldscore.py
│       ├── sumstats_munger.py
│       ├── regression.py
│       ├── _jackknife.py
│       ├── _irwls.py
│       ├── formats.py
│       └── identifiers.py
├── tests/                   # workflow and parity checks
└── tutorials/               # end-to-end usage examples
```

## Code Map

### `ldsc.cli`

This is the only command-line entry point. It parses subcommands, preserves the user-facing CLI contract, and hands control to the public workflow modules. Architecture invariant: CLI code should dispatch only; it should not contain numerical logic.

### `ldsc.config`, `ldsc.path_resolution`, `ldsc.column_inference`, `ldsc.chromosome_inference`

These modules define the package-wide contracts that every workflow shares. They normalize path tokens, chromosome ordering, genome-build aliases, SNP identifier modes, and input header aliases before the kernel is called. Architecture invariant: path discovery and header inference happen here or in the workflow layer, never inside `_kernel`.

### `ldsc.annotation_builder`

This is the public interface for annotation loading and BED projection. It re-exports `AnnotationBuilder`, `AnnotationBundle`, and the legacy-compatible projection helpers while the implementation lives in `ldsc._kernel.annotation`. Public interface: users should start here rather than importing the kernel directly.

### `ldsc.ref_panel_builder`

This module builds standard parquet LD reference artifacts from PLINK inputs. It handles genetic-map loading, liftover selection, restriction filtering, and output-path construction, then delegates pairwise-LD emission to the kernel. Architecture invariant: emitted parquet schemas are part of the public file contract for parquet-backed LDSC workflows.

### `ldsc.ldscore_calculator`

This module orchestrates chromosome-wise LD-score computation. It resolves annotation and reference-panel inputs, builds per-chromosome runs, aggregates them into `LDScoreResult`, and routes artifact writing through `ldsc.outputs`. Architecture invariant: computation stays chromosome-wise; the aggregate result is assembled only after all chromosome runs finish.

### `ldsc.sumstats_munger`

This module wraps the historical munging behavior in typed public objects such as `RawSumstatsSpec`, `SumstatsTable`, and `SumstatsMunger`. It keeps raw-input parsing flexible but reloads curated `.sumstats(.gz)` artifacts through stricter internal column contracts. Public interface: downstream workflows should consume `SumstatsTable` or curated `.sumstats.gz`, not raw heterogeneous GWAS tables.

### `ldsc.regression_runner`

This module rebuilds an `LDScoreResult` from on-disk artifacts, merges it with munged sumstats, drops zero-variance LD-score columns, and dispatches to the regression kernel for `h2`, partitioned `h2`, and `rg`. Architecture invariant: regression only consumes aggregated LD-score artifacts; it does not recompute LD scores.

### `ldsc.outputs`

This is the artifact writer for LD-score runs. It owns normalized per-chromosome
`.l2.ldscore.gz` files with embedded `regr_weight`, plus `.l2.M`, `.l2.M_5_50`,
manifests, and summary files. Architecture invariant: new postprocessing should
be added through producers here instead of feature-specific ad hoc file writes.

### `ldsc._kernel.*`

The kernel layer contains the actual numerical methods and low-level readers. It includes annotation parsing, PLINK/parquet reference-panel access, LD-score math, legacy-compatible summary-statistics munging, and regression estimators. Private boundary: `_kernel` is not the supported import surface and should receive only resolved primitive inputs.

## Cross-Cutting Concerns

- **Input token language**: public inputs accept exact paths, globs, `@` chromosome suites, and some legacy bare prefixes. Output paths are normalized but remain literal destinations.
- **Column and identifier normalization**: `column_inference.py` is the single source of truth for raw-input aliases, internal artifact headers, SNP identifier modes, and genome-build aliases.
- **Artifact compatibility**: LDSC-compatible outputs remain the public contract for downstream chaining: `.annot.gz`, `.sumstats.gz`, `.l2.ldscore(.gz)`, `.l2.M`, and `.l2.M_5_50`. Legacy separate `.w.l2.ldscore(.gz)` files remain loadable for regression, but the refactored writer embeds weights as `regr_weight` instead of emitting that artifact by default.
- **Chromosome ordering**: chromosome-sharded inputs are validated and reassembled in stable genomic order by the workflow layer.
- **Testing approach**: tests under `tests/` cover file contracts, workflow behavior, and legacy compatibility expectations.

## Architectural Invariants

- Import from `ldsc`, not from `ldsc._kernel`, for supported public use.
- Kernel code must not resolve globs, `@` suites, or other user-facing path tokens.
- `column_inference.py` owns alias families and internal artifact header strictness.
- LD-score computation remains chromosome-wise; regression consumes only the aggregated artifacts.
- Output extensibility goes through `ldsc.outputs`, not through feature-specific one-off writers.
