"""Per-SNP annotations: reading, building and checking annotation values.

Annotation file parsing and value checks, the lazily read annotation
dataset, gene-list and BED query projection onto an annotation SNP grid,
query viability, the annotation overlap computation, and the ``ldsc annotate``
workflow and its writer. This is tier 5 of the package layout, a peer of
``ldsc._refpanel`` and ``ldsc._sumstats``; peers do not import each other.
"""
