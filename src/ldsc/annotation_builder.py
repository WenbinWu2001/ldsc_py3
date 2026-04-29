"""Public workflow module for SNP-level annotation building.

This package-facing module re-exports the annotation workflow objects while the
actual file discovery and bundle assembly remain in the workflow layer. Public
annotation inputs may be exact paths, standard Python glob patterns, or
explicit chromosome-suite tokens using ``@``. Those tokens are resolved to
concrete files before the internal kernel code runs.
"""

from __future__ import annotations

from .config import AnnotationBuildConfig
from ._kernel.annotation import (
    AnnotationBuilder,
    AnnotationBundle,
    main_bed_to_annot,
    parse_bed_to_annot_args,
    run_bed_to_annot,
)

__all__ = [
    "AnnotationBuildConfig",
    "AnnotationBuilder",
    "AnnotationBundle",
    "main_bed_to_annot",
    "parse_bed_to_annot_args",
    "run_bed_to_annot",
]
