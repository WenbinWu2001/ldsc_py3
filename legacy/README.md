# Legacy Codebase

This directory contains the preserved old Python 3 LDSC codebase.

It is kept as a separate tracked subtree so the refactor can evolve without mixing old and new layouts in the same directory.

Typical setup:

```bash
cd legacy
pip install -e .
python ldsc.py -h
python munge_sumstats.py -h
```

Use this subtree only for legacy reference, parity checks, and old-code maintenance.
