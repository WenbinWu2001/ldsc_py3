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

## Python API

```python
from ldsc import (
    AnnotationBuilder,
    LDScoreCalculator,
    ReferencePanelBuilder,
    RegressionRunner,
    SumstatsMunger,
    load_sumstats,
)
```

## Input Path Tokens

Public workflow APIs accept normalized string tokens for inputs:

- exact paths
- Python glob patterns such as `annotations/*.annot.gz`
- explicit chromosome suites using `@`, for example `baseline.@`
- legacy bare prefixes such as `baseline.`

Inputs are resolved before the internal kernel runs. Public outputs use fixed
filenames inside the selected `output_dir`.

## Verification

```bash
pytest
```

The suite also remains compatible with the standard-library runner while the
project transitions fully to pytest:

```bash
python -m unittest discover -s tests -p 'test*.py' -v
```
