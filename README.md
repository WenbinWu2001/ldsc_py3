# ldsc_py3_restructured

This repository is the active refactored LDSC package.

## Structure

- `src/ldsc/`: public package surface
- `src/ldsc/_kernel/`: internal compute and file-format modules
- `tests/`: local parity and workflow tests
- `tutorials/`: package-level usage examples
- `docs/architecture.md`, `docs/code_structure.md`, `docs/class-and-features.md`: active design and navigation docs

## Install

```bash
pip install -e .
```

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
- `ldsc build-ref-panel`
- `ldsc ldscore`
- `ldsc munge-sumstats`
- `ldsc h2`
- `ldsc partitioned-h2`
- `ldsc rg`

`ldsc munge-sumstats` writes canonical `SNP`, `CHR`, `POS`, `Z`, and `N`
columns when possible and always includes `CHR`/`POS` in package-written
artifacts. Raw files may provide chromosome and position columns through common
aliases such as `#CHROM`, `CHROM`, `CHR`, `POS`, or `BP`, or explicitly through
`--chr` and `--pos`. Leading `##` metadata/comment lines are skipped before the
real header is parsed. Each munged run also writes `sumstats.metadata.json`
beside `sumstats.sumstats.gz` so later regression commands can recover
`snp_identifier` and `genome_build` provenance.

`ldsc ldscore` supports ordinary unpartitioned LD-score generation without
baseline annotations:

```bash
ldsc ldscore --output-dir ldscores --plink-prefix panel --ld-wind-snps 10
ldsc h2 --sumstats-file trait.sumstats.gz --ldscore-dir ldscores --output-dir h2_out
```

When no baseline and no query annotations are supplied, the workflow writes a
synthetic all-ones baseline column named exactly `base` in `baseline.parquet`.
Query annotation inputs still require explicit `--baseline-annot-sources`.

## Python API

```python
from ldsc import (
    AnnotationBuilder,
    ChrPosBuildInference,
    LDScoreCalculator,
    ReferencePanelBuilder,
    RegressionRunner,
    SumstatsMunger,
    infer_chr_pos_build,
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
- legacy bare prefixes such as `baseline.`

Inputs are resolved before the internal kernel runs. Public outputs use fixed
filenames inside the selected `output_dir`.

## Output Collision Policy

Every workflow treats `--output-dir` or `output_dir` as a directory:

- missing output directories are created with a warning
- existing directories are reused
- known output files fail the run before writing starts
- reruns that intentionally replace known files must pass `--overwrite` on the
  CLI or `overwrite=True` in Python

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
