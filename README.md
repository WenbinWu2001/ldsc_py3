# ldsc3_Jerry

Last updated on: 2026-08-04

This repository is the active refactored LDSC package.

## Structure

- `src/ldsc/`: public package surface
- `src/ldsc/_kernel/`: internal compute and file-format modules
- `tests/`: local parity and workflow tests
- `tutorials/`: package-level usage examples
- `docs/current/architecture.md`, `docs/current/code-structure.md`, `docs/current/class-and-features.md`, `docs/current/workflow-logging.md`, `docs/current/liftover-harmonization-decisions.md`: active design and navigation docs

## Install

Choose the stable branch for routine use, or the `restructure` branch for the
current development version. Both installs use an editable package install so
local source changes are picked up immediately.

### Stable version

Clone the default `main` branch to a directory `ldsc3`:

```bash
git clone https://github.com/WenbinWu2001/ldsc_py3.git ldsc3
cd ldsc3
```

Create the conda environment and install the package with development extras:

```bash
conda env create -f environment.yml
conda activate ldsc3
python -m pip install -e ".[dev]"
```

**For later runs,** activate the environment before launching LDSC commands or scripts:

```bash
conda activate ldsc3
```

To update the package, enter the repository and pull the latest changes from the `main` branch:

```bash
cd ldsc3
git pull origin main
conda activate ldsc3
```

Because the package is installed in editable mode, source-code updates take effect immediately. 
If `environment.yml` or the package dependencies have changed, also update the environment and reinstall the package:

```bash
conda env update -f environment.yml --prune
python -m pip install -e ".[dev]"
```


### Development version

Clone the `restructure` branch to a directory `ldsc3-dev/`:

```bash
git clone --branch restructure https://github.com/WenbinWu2001/ldsc_py3.git ldsc3-dev
cd ldsc3-dev
```

Create a separate development environment and install the package with
development extras:

```bash
conda env create -f environment.yml -n ldsc3-dev
conda activate ldsc3-dev
python -m pip install -e ".[dev]"
```

**For later runs,** activate the environment before launching LDSC commands or scripts:

```bash
conda activate ldsc3-dev
```

To update the development version, enter the repository and pull the latest changes from the `restructure` branch:

```bash
cd ldsc3-dev
git pull origin restructure
conda activate ldsc3-dev
```

Because the package is installed in editable mode, source-code updates take effect immediately.
If `environment.yml` or the package dependencies have changed, also update the development environment and reinstall the package:

```bash
conda env update -f environment.yml -n ldsc3-dev --prune
python -m pip install -e ".[dev]"
```

The package supports Python 3.11 through 3.13. The base install includes the
core NumPy/pandas/SciPy/PyArrow stack. Optional extras are split by workflow:
`.[plink]` installs `bitarray` for PLINK-backed LD computation, `.[bed]`
installs `pybedtools` for BED projection, and `.[liftover]` installs
`pyliftover` for chain-file liftover in sumstats munging and cross-build
reference-panel output. `.[dev]` installs all of
those extras plus pytest. BED-based annotation projection also requires the
external `bedtools` executable, which `environment.yml` installs from bioconda.
For non-conda installs, make sure `bedtools` is available on `PATH` before
running BED annotation workflows.

## CLI

```bash
ldsc --help
```

Equivalent development entrypoint:

```bash
python -m ldsc --help
```

Subcommands:

- `ldsc annotate`
- `ldsc build-gene-ldscore-index`
- `ldsc build-ref-panel`
- `ldsc convert-ldsc2-ldscores`
- `ldsc ldscore`
- `ldsc munge-sumstats`
- `ldsc h2`
- `ldsc partitioned-h2`
- `ldsc rg`

`ldsc munge-sumstats` writes `sumstats.parquet` by default, with optional
legacy `sumstats.sumstats.gz` output through `--output-format tsv.gz` or
`--output-format both`. Package-written artifacts use canonical `SNP`, `CHR`,
`POS`, `Z`, and `N` columns when possible and always include `CHR`/`POS`. Raw
files may provide chromosome and position columns through common aliases such as
`#CHROM`, `CHROM`, `CHR`, `POS`, or `BP`, or explicitly through `--chr` and
`--pos`; pass the raw input as `--raw-sumstats-file`. Leading `##`
metadata/comment lines are skipped before the real header is parsed. The
`sumstats.parquet` is self-describing: its `snp_identifier`, `genome_build`, and
optional `--trait-name` provenance ride in the Parquet footer, so later
regression commands need only that one file -- no `metadata.json` sidecar is
written. The legacy `sumstats.sumstats.gz` carries no embedded metadata and is
treated as an rsID lookup artifact at regression time. Footerless Parquet is
rejected rather than guessed.
Detailed coordinate and liftover bookkeeping is written to `sumstats.log`. The default
`snp_identifier` is `chr_pos_allele_aware`, which requires usable `A1/A2`; rerun
with `--snp-identifier chr_pos` to use coordinate identity without
allele-aware matching. The legacy `--no-alleles` escape hatch is no longer
accepted; choose the base `rsid` or `chr_pos` identity mode instead. Use `--sumstats-snps-file`
when the munged artifact should be restricted to a headered SNP keep-list, or
`--use-hm3-snps` to use the packaged curated HM3 map. Restriction files may omit
alleles; allele-free restrictions, including packaged HM3 restrictions, match by
base key. Allele-bearing restrictions in allele-aware modes match by the
effective allele-aware key. Restriction files are identity-only filters:
duplicate restriction keys collapse to one retained key, and non-identity
columns such as `CM` or `MAF` are ignored. These filters are loaded once before
raw chunk parsing and applied while chunks are still streaming, after canonical
columns and coordinate normalization are available. They keep matching rows only
and do not reorder the output. In `chr_pos`-family modes, `SNP` is a label; matching uses
source-build `CHR` and `POS`, plus the allele set in
`chr_pos_allele_aware`, represented internally as effective coordinate-family
keys. Rows with missing or invalid coordinates are dropped and counted at
coordinate match/map stages. The base `chr_pos` mode uses coordinate identity
without allele-aware matching. To convert coordinates after QC and after SNP
restriction, pass `--target-genome-build` with `--liftover-chain-file`, or pair
`--target-genome-build --use-hm3-snps --use-hm3-quick-liftover` for the HM3
coordinate shortcut. Liftover is invalid in `rsid`-family modes because
positions are not the row identity there.

`ldsc build-ref-panel` keeps a separate source-build contract for PLINK input:
provide or infer `--source-genome-build`, and a matching chain file emits the
opposite build. Deliberate reference-universe restriction uses an explicit
`--ref-panel-snps-file`; the builder has no HM3-only restriction or quick-liftover
mode. Chain-file liftover is invalid when the active
`snp_identifier` is in the `rsid` family; omit liftover for source-build-only rsID panels. In
`chr_pos`-family modes,
duplicate source or target coordinate groups are dropped by default
(`--duplicate-position-policy drop-all`), with details in `build-ref-panel.log`
and duplicate-only sidecars under `dropped_snps/`.

Artifact-writing workflows write completed per-run logs under their output
directories. During `build-gene-ldscore-index`, the open log temporarily lives
under hidden `.<index-name>.build-state/` so it is never part of an atomic index
replacement; after its handler closes, a successful log moves into the index's
`diagnostics/`. `munge-sumstats` keeps the historical `sumstats.log` name; other
commands use `annotate.log`, `ldscore.log`, `build-ref-panel.log`,
`build-gene-ldscore-index.log`,
`h2.log`, `partitioned-h2.log`, or `rg.log`. Concrete single-chromosome
`build-ref-panel` runs use `build-ref-panel.chr<chrom>.log` so parallel
per-chromosome jobs can share an output directory without sharing one log file.
Logs are audit artifacts, so Python result objects and `output_paths` mappings
only list scientific data outputs.

`ldsc ldscore` supports ordinary unpartitioned LD-score generation without
baseline annotations:

```bash
ldsc ldscore --output-dir ldscores --plink-prefix panel --ld-wind-snps 10
ldsc h2 --sumstats-file trait/sumstats.parquet --ldscore-dir ldscores --output-dir h2_out
```

When no baseline and no query annotations are supplied, the workflow writes a
synthetic all-ones baseline column named exactly `base` in
`ldscore.baseline.parquet`.
Query annotation inputs still require explicit `--baseline-annot-sources`.
Use this synthetic `base` directory for `ldsc h2` or `ldsc rg`. A baseline-only
directory is also accepted by `ldsc partitioned-h2` in its functional-category
regime, although a single all-ones `base` column is a degenerate one-category
fit rather than a meaningful partitioned analysis.
The LD-score parquet files remain flat `ldscore.baseline.parquet` and
`ldscore.query.parquet` files, but they are written with one row group per chromosome. The metadata
records `row_group_layout`, `baseline_row_groups`, and `query_row_groups` so
callers can load a single chromosome by row-group index when needed.

Gene-list LD scores may also use an explicitly installed exact profile:

```bash
ldsc ldscore \
  --gene-ldscore-index-dir indexes/1000G_EUR_baseline_100kb \
  --query-annot-gene-list-sources gene_lists/immune.txt \
  --output-dir gene_ldscores
```

Build complete indexes offline with `ldsc build-gene-ldscore-index`. Indexed mode is explicit and fail-closed: it validates the complete index, accepts no live baseline/reference/window settings, and writes the same self-contained canonical LD-score directory. During the offline build, each completed chromosome is atomically persisted in a hidden run-specific stage, allowing its large in-memory payload to be released while later chromosomes continue. A `Finished chromosome N` log line reports that private durability boundary; the public destination remains absent, empty, or at its prior complete version until every chromosome and shared metadata pass reload validation. These private shards are not resumable checkpoints and never support chromosome append or incremental index updates. See
[the exact gene-index guide](docs/current/gene-ldscore-index.md) and its
[mathematical algorithm](docs/current/gene-ldscore-index-mathematics.md).
Task-oriented walkthroughs cover [building the index](docs/wiki/utility-functionalities/build-gene-ldscore-index.md)
and [using it for gene-list LD scores](docs/wiki/main-functionalities/ldscore.md).

## LDSC2 compatibility boundary

Backward compatibility is deliberately asymmetric and limited to artifacts
that are costly or impractical for users to recreate:

- Regression accepts genuine LDSC2 `.sumstats` and `.sumstats.gz` text files
  automatically. Their `SNP` values are treated as rsID lookup keys and are
  projected onto the canonical LDSC3 LD-score panel. `A1` and `A2` are required;
  allele orientation is validated and `Z` is negated when a swap is needed.
  `FRQ` is optional, is never imputed from the panel, and does not control SNP
  retention. Footerless Parquet is not treated as a legacy artifact.
- LDSC2 LD-score fragments are never accepted directly by regression. Run
  `ldsc convert-ldsc2-ldscores` explicitly with a complete reference directory
  and a complete regression-weight directory. The converter supports a
  one-column unpartitioned suite or a complete baseline partitioned suite; the
  latter also requires a frequency directory. Query/cell-type suites and thin
  annotations are intentionally unsupported.
- Converted suites remain allele-unaware (`rsid` by default, or `chr_pos`).
  `.l2.M_5_50` is required and retains the fixed strict LDSC2 common-frequency
  rule. There is no converter threshold flag. Missing `.l2.M` is tolerated only
  under the documented count policy; requesting unavailable all-SNP counts
  later is an error.
- LDSC3 does not promise that its outputs can be fed back into LDSC2. The
  `BP` header written in text `.annot.gz` files is a narrow interoperability
  convenience, and LDSC3 accepts either `BP` or `POS` when reading them.

See the complete policies for
[legacy sumstats](docs/current/legacy-sumstats-compatibility.md) and
[legacy LD-score conversion](docs/current/legacy-ldscore-conversion.md).

## Python API

```python
from ldsc import (
    AnnotationBuilder,
    ChrPosBuildInference,
    GeneLDScoreIndexBuildConfig,
    LDScoreCalculator,
    ReferencePanelBuilder,
    RegressionRunner,
    SumstatsMunger,
    build_gene_ldscore_index,
    infer_chr_pos_build,
    load_gene_ldscore_index,
    load_sumstats,
    resolve_chr_pos_table,
)
```

Genome-build inference for `chr_pos` tables is public through the Python API:
use `infer_chr_pos_build()` when you only need the `hg19`/`hg38` decision, and
use `resolve_chr_pos_table()` when you also want 0-based inputs converted to
canonical 1-based coordinates. The CLI exposes this behavior inside workflows
with `--genome-build auto`, including during `munge-sumstats`; there is no
standalone inference subcommand.

## Input Path Tokens

Public workflow APIs accept normalized string tokens for inputs:

- exact paths
- Python glob patterns such as `annotations/*.annot.gz`
- explicit chromosome suites using `@`, for example `baseline.@`
- PLINK prefix tokens for reference-panel inputs, such as `panel_chr@`

Inputs are resolved before the internal kernel runs. Public outputs use fixed
filenames inside the selected `output_dir`.

## Output Collision Policy

Every workflow treats `--output-dir` or `output_dir` as a directory:

- missing output directories are created with a warning
- existing directories are reused
- known output files fail the run before writing starts
- reruns that intentionally replace known files must pass `--overwrite` on the
  CLI or `overwrite=True` in Python
- per-run log files are preflighted with the scientific outputs, so a collision
  fails before a new log is opened

The overwrite flag applies only to the fixed files owned by that workflow. It
does not remove unrelated files and never cleans a whole directory.

## Verification

```bash
pytest
```

The suite also remains compatible with the standard-library runner while the
project transitions fully to pytest:

```bash
python -m unittest discover -s tests -p 'test*.py' -v
```
