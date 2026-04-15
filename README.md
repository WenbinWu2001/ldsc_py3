# ldsc_py3_Jerry

This repository now contains two tracked codebases under the same Git root:

- `refactor/`: the active refactored codebase
- `legacy/`: the preserved old codebase

Use each subtree as its own working directory.

## Refactor Codebase

Path:

- `refactor/`

Contains:

- the active workflow/domain refactor
- current design docs
- current `PLANS.md`
- current test suite

Typical commands:

```bash
cd refactor
pip install -e .
python ldsc_new.py -h
python -m unittest discover -s tests -p 'test*.py' -v
```

## Legacy Codebase

Path:

- `legacy/`

Contains:

- the preserved old Python 3 LDSC port
- its original docs, resources, and CLI layout

Typical commands:

```bash
cd legacy
pip install -e .
python ldsc.py -h
python munge_sumstats.py -h
```

## Repo-Level Rule

Do not add active code back to the repository root. Keep code, docs, tests, and plans inside either `refactor/` or `legacy/`.
