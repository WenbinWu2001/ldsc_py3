"""Public workflow module for SNP-level annotation building.

This package-facing module re-exports the annotation workflow objects while the
actual file discovery and bundle assembly remain in the workflow layer. Public
annotation inputs may be exact paths, standard Python glob patterns, explicit
chromosome-suite tokens using ``@``, or legacy bare prefixes such as
``baseline.``. Those tokens are resolved to concrete files before the internal
kernel code runs.
"""

from __future__ import annotations

from ._kernel.annotation import (
    AnnotationBuilder,
    AnnotationBundle,
    AnnotationSourceSpec,
    gene_set_to_bed,
    main_bed_to_annot,
    main_make_annot,
    make_annot_files,
    parse_bed_to_annot_args,
    parse_make_annot_args,
    run_bed_to_annot,
)

__all__ = [
    "AnnotationBuilder",
    "AnnotationBundle",
    "AnnotationSourceSpec",
    "gene_set_to_bed",
    "main_bed_to_annot",
    "main_make_annot",
    "make_annot_files",
    "parse_bed_to_annot_args",
    "parse_make_annot_args",
    "run_bed_to_annot",
]
