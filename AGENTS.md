# Repository Layout

This Git repository now contains two tracked codebases:

- `refactor/`: active refactored codebase
- `legacy/`: preserved old codebase

Do not place active package code at the repository root.

# Build Commands

Refactor:

```bash
cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_Jerry_workspace/ldsc_py3_Jerry/.worktrees/restructure/refactor
pip install -e .
python ldsc_new.py -h
python -m unittest discover -s tests -p 'test*.py' -v
```

Legacy:

```bash
cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_Jerry_workspace/ldsc_py3_Jerry/.worktrees/restructure/legacy
pip install -e .
python ldsc.py -h
python munge_sumstats.py -h
```

# Key Rules

- Keep refactor work inside `refactor/`.
- Keep old-code preservation work inside `legacy/`.
- Update design docs together with structural changes in `refactor/`.
- Preserve legacy behavior and file contracts unless the refactor docs explicitly change them.
