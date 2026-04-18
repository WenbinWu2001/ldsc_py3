# ldsc_py3_Jerry

This worktree is organized around the active refactored codebase under `refactor/`.

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

## Repo-Level Rule

Do not add active code back to the repository root. Keep code, docs, tests, and plans inside `refactor/`.
