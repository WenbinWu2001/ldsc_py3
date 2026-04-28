"""Module entry point for ``python -m ldsc``."""

from .cli import main, run_cli


if __name__ == "__main__":
    raise SystemExit(run_cli())
