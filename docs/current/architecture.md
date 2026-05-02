# Architecture

`ldsc_py3_Jerry` is the refactored Python 3 LDSC package. It reads optional SNP-level annotations, PLINK or parquet R2 references, and GWAS summary statistics; resolves user-facing path and header conventions in the public workflow layer; delegates numerical work to `ldsc._kernel`; and writes LDSC-compatible artifacts that can be chained into later runs.

Related docs:

- [data-flow.md](data-flow.md): user-visible file streams and flowcharts
- [class-and-features.md](class-and-features.md): public API surface and major types
- [code-structure.md](code-structure.md): module map and change guide
- [layer-structure.md](layer-structure.md): layer-by-function object matrix

![Package overview](../assets/ldsc-package-overview.png)

## Bird's-Eye View

- **Build query annotations**: project BED intervals onto a baseline SNP grid. Entry points: `ldsc annotate`, `ldsc.AnnotationBuilder`
- **Build parquet reference panels**: convert PLINK genotype panels into standard parquet R2 artifacts. Entry points: `ldsc build-ref-panel`, `ldsc.ReferencePanelBuilder`
- **Compute LD scores**: align annotations to a reference panel and emit LDSC-compatible LD-score artifacts; ordinary unpartitioned runs may omit annotations and receive a synthetic all-ones `base` annotation. Entry points: `ldsc ldscore`, `ldsc.run_ldscore()`, `ldsc.LDScoreCalculator`
- **Munge raw summary statistics**: normalize raw GWAS tables into curated `.sumstats.gz` artifacts. Entry points: `ldsc munge-sumstats`, `ldsc.SumstatsMunger`
- **Run LDSC regression**: consume munged sumstats and LD-score artifacts to estimate `h2`, partitioned `h2`, or `rg`. Entry points: `ldsc h2`, `ldsc partitioned-h2`, `ldsc rg`, `ldsc.RegressionRunner`

## Layer Structure

- **CLI Layer**: public command dispatch in `ldsc.cli`
- **Workflow And Preprocessing Layer**: public services in `ldsc.annotation_builder`, `ldsc.ref_panel_builder`, `ldsc.ldscore_calculator`, `ldsc.sumstats_munger`, and `ldsc.regression_runner`, plus shared normalization in `ldsc.config`, `ldsc.path_resolution`, `ldsc.column_inference`, `ldsc.chromosome_inference`, and `ldsc.genome_build_inference`
- **Compute Kernel**: private file-format and numerical code in `ldsc._kernel.*`
- **Output Layer**: canonical LD-score and partitioned-h2 artifact writing in `ldsc.outputs`, plus fixed h2/rg summary writers in `ldsc.regression_runner`

## File Tree

```text
ldsc_py3_Jerry/
├── docs/
│   ├── current/
│   │   ├── architecture.md      # package-level architecture
│   │   ├── class-and-features.md
│   │   ├── code-structure.md
│   │   └── data-flow.md         # workflow-level file streams
│   └── assets/
├── src/ldsc/
│   ├── __init__.py          # public Python exports
│   ├── __main__.py          # `python -m ldsc` entry point
│   ├── cli.py               # unified `ldsc` command
│   ├── config.py            # validated public config dataclasses
│   ├── path_resolution.py   # input-token normalization and discovery
│   ├── column_inference.py  # header and identifier normalization
│   ├── chromosome_inference.py
│   ├── genome_build_inference.py
│   ├── annotation_builder.py
│   ├── ref_panel_builder.py
│   ├── ldscore_calculator.py
│   ├── sumstats_munger.py
│   ├── regression_runner.py
│   ├── outputs.py           # LD-score and partitioned-h2 artifact writers
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

### `ldsc.config`, `ldsc.path_resolution`, `ldsc.column_inference`, `ldsc.chromosome_inference`, `ldsc.genome_build_inference`

These modules define the package-wide contracts that every workflow shares. They normalize path tokens, chromosome ordering, genome-build aliases, SNP identifier modes, input header aliases, and `chr_pos` genome-build inference before the kernel is called. `ldsc.path_resolution` also owns output-directory creation and fixed-artifact collision preflight through `ensure_output_paths_available()`. `ldsc.genome_build_inference` is public for Python callers through the top-level `ldsc` exports, while the CLI exposes it only through `--genome-build auto` on existing workflows. Architecture invariant: path discovery, header inference, output preflight, and user-facing build inference happen here or in the workflow layer, never inside `_kernel`.

### `ldsc.annotation_builder`

This is the public interface for annotation loading and BED projection. It re-exports `AnnotationBuilder`, `AnnotationBundle`, and `run_bed_to_annot()` while the implementation lives in `ldsc._kernel.annotation`. Public interface: users should start here rather than importing the kernel directly.

### `ldsc.ref_panel_builder`

This module builds standard parquet R2 reference artifacts from PLINK inputs. It handles optional genetic-map loading, optional liftover selection, restriction filtering, and output-path construction, then delegates pairwise-R2 emission to the kernel. Architecture invariant: emitted parquet schemas are part of the public file contract for parquet-backed LDSC workflows.

### `ldsc.ldscore_calculator`

This module orchestrates chromosome-wise LD-score computation. It resolves annotation and reference-panel inputs, synthesizes the all-ones `base` annotation when an unpartitioned run omits baseline/query inputs, builds per-chromosome runs, aggregates them into `LDScoreResult`, and routes artifact writing through `ldsc.outputs`. Architecture invariant: computation stays chromosome-wise; the aggregate result is assembled only after all chromosome runs finish.

### `ldsc.sumstats_munger`

This module wraps the historical munging behavior in typed public objects such as `MungeConfig`, `SumstatsTable`, and `SumstatsMunger`. Raw heterogeneous GWAS input is named explicitly as `raw_sumstats_file` in Python and `--raw-sumstats-file` in the CLI. `main()` parses only, `run_munge_sumstats_from_args()` maps CLI arguments into `MungeConfig`, and `SumstatsMunger.run()` owns path resolution, fixed-output preflight, `sumstats.log`, metadata sidecars, and result construction before delegating primitive munging work to the kernel. It keeps raw-input parsing flexible, skips leading `##` raw metadata lines, writes canonical `CHR` and `POS` columns alongside `SNP`, `Z`, and `N`, optionally applies a headered `sumstats_snps_file` keep-list through the shared identifier parser, and persists `snp_identifier`/`genome_build` provenance in `sumstats.metadata.json`. Public interface: downstream workflows should consume `SumstatsTable` or curated `.sumstats.gz` plus its sidecar, not raw heterogeneous GWAS tables.

### `ldsc.regression_runner`

This module rebuilds an `LDScoreResult` from on-disk artifacts, merges it with munged sumstats, drops zero-variance LD-score columns, and dispatches to the regression kernel for `h2`, partitioned `h2`, and `rg`. In `rsid` mode the merge uses `SNP`; in `chr_pos` mode it builds a private normalized `CHR:POS` key from sumstats and LD-score coordinates. Architecture invariant: regression only consumes aggregated LD-score artifacts; it does not recompute LD scores.

### `ldsc.outputs`

This is the canonical LD-score result-directory writer. It owns fixed files
inside `output_dir`: `manifest.json`, `baseline.parquet`, and optional
`query.parquet`. The parquet files stay flat for compatibility, but are written
with one row group per chromosome; `manifest.json` records row-group layout and
per-chromosome offsets for chromosome-scoped reads. It reuses existing
directories but refuses existing canonical files unless
`LDScoreOutputConfig(overwrite=True)` is supplied. Architecture invariant:
public output customization chooses the directory name and explicit overwrite
policy, not per-run filename prefixes.

### `ldsc._kernel.*`

The kernel layer contains the actual numerical methods and low-level readers. It includes annotation parsing, PLINK/parquet reference-panel access, LD-score math, legacy-compatible summary-statistics munging, and regression estimators. Private boundary: `_kernel` is not the supported import surface and should receive only resolved primitive inputs.

## Cross-Cutting Concerns

- **Input token language**: public inputs accept exact paths, globs, `@` chromosome suites, and some legacy bare prefixes. Output paths are normalized but remain literal destinations.
- **Column and identifier normalization**: `column_inference.py` is the single source of truth for raw-input aliases, including `#CHROM`/`CHROM` as `CHR`, internal artifact headers, SNP identifier modes, and genome-build aliases. `genome_build_inference.py` owns automatic hg19/hg38 and 0-based/1-based inference for `chr_pos` tables.
- **Artifact compatibility**: Public downstream chaining uses `.annot.gz`, `.sumstats.gz` with optional `sumstats.metadata.json`, and canonical LD-score result directories. LD-score parquet payloads are flat files with chromosome-aligned row groups, so full-file readers still work while chromosome-specific readers can skip unrelated row groups. Legacy `.l2.ldscore(.gz)`, `.l2.M`, `.l2.M_5_50`, and separate `.w.l2.ldscore(.gz)` files remain internal/legacy file-format concerns, not the public LD-score writer contract.
- **Output collision handling**: output directories are literal destinations.
  Missing directories are created, existing directories are reused, and fixed
  output files are checked before writing. By default an existing artifact
  raises `FileExistsError`; `--overwrite` or `overwrite=True` makes replacement
  explicit without deleting unrelated files or cleaning the directory.
- **Dependency split**: base package dependencies cover core pandas/numpy/SciPy
  workflows and parquet I/O. PLINK-backed LD computation requires the
  `plink` extra (`bitarray`), BED projection requires the `bed` extra
  (`pybedtools` plus the external `bedtools` executable), and cross-build
  reference-panel output requires the `liftover` extra (`pyliftover`).
- **Chromosome ordering**: chromosome-sharded inputs are validated and reassembled in stable genomic order by the workflow layer.
- **Testing approach**: tests under `tests/` cover file contracts, workflow behavior, and legacy compatibility expectations.

## Architectural Invariants

- Import from `ldsc`, not from `ldsc._kernel`, for supported public use.
- Kernel code must not resolve globs, `@` suites, or other user-facing path tokens.
- `column_inference.py` owns alias families and internal artifact header strictness.
- LD-score computation remains chromosome-wise; regression consumes only the aggregated artifacts.
- Public LD-score output layout is fixed by `LDScoreDirectoryWriter`; optional
  per-query partitioned-h2 output layout is fixed by
  `PartitionedH2DirectoryWriter`.
- Query annotations are valid only with explicit baseline annotations; the synthetic `base` path is for ordinary unpartitioned LD-score generation.
- Every workflow that writes fixed artifacts must precompute expected output
  paths and call the shared output preflight before the first write.
