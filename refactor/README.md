# ldsc_py3_restructured

This repository is the active refactored LDSC package.

## Structure

- `src/ldsc/`: public package surface
- `src/ldsc/_kernel/`: internal compute and file-format modules
- `tests/`: local parity and workflow tests
- `tutorials/`: package-level usage examples
- `architecture.md`, `code_structure.md`, `class-and-features.md`: active design and navigation docs

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
python -m ldsc.cli --help
```

Subcommands:

- `ldsc annotate`
- `ldsc ldscore`
- `ldsc munge-sumstats`
- `ldsc h2`
- `ldsc partitioned-h2`
- `ldsc rg`

## Python API

```python
from ldsc import AnnotationBuilder, LDScoreCalculator, RegressionRunner, SumstatsMunger
```

## Input Path Tokens

Public workflow APIs accept normalized string tokens for inputs:

- exact paths
- Python glob patterns such as `annotations/*.annot.gz`
- explicit chromosome suites using `@`, for example `baseline.@`
- legacy bare prefixes such as `baseline.`

Inputs are resolved before the internal kernel runs. Output paths such as `out_prefix` stay literal destinations.

## Verification

```bash
python -m unittest discover -s tests -p 'test*.py' -v
```
