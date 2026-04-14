#!/usr/bin/env python
"""Thin wrapper around the shared single-annotation build workflow."""

from __future__ import annotations

try:
    from ldscore.annotation import gene_set_to_bed, main_make_annot, make_annot_files
except ImportError:  # pragma: no cover - package import fallback
    from .ldscore.annotation import gene_set_to_bed, main_make_annot, make_annot_files


def main(argv=None) -> int:
    return main_make_annot(argv)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
