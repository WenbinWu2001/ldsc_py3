# Repository Layout

This Git repository is organized around one active codebase:

- `refactor/`: active refactored codebase

Do not place active package code at the repository root.

# Build Commands

Refactor:

```bash
cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/refactor
pip install -e .
python ldsc_new.py -h
python -m unittest discover -s tests -p 'test*.py' -v
```

# Key Rules

- Keep refactor work inside `refactor/`.
- Update design docs together with structural changes in `refactor/`.
- Preserve legacy behavior and file contracts unless the refactor docs explicitly change them.
