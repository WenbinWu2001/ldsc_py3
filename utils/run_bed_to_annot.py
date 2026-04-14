"""Thin wrapper around the shared annotation projection workflow."""

from __future__ import annotations

import sys

try:
    from ldscore.annotation import main_bed_to_annot, run_bed_to_annot
except ImportError:  # pragma: no cover - package import fallback
    from ..ldscore.annotation import main_bed_to_annot, run_bed_to_annot


def main(argv=None) -> int:
    return main_bed_to_annot(argv)


if __name__ == "__main__":  # pragma: no cover
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        sys.exit(1)
