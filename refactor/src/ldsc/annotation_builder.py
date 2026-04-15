"""Public workflow module for SNP-level annotation building."""

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
