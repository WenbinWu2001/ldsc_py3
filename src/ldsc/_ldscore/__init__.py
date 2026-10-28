"""LD scores: calculation, the LD-score directory, gene LD-score indexes and legacy conversion.

``compute`` holds the LD-score math over a prepared reference chromosome;
``calculator`` orchestrates the ``ldscore`` command; ``source`` reads an
LD-score directory and ``outputs`` writes one; ``gene_index`` builds and
scores exact gene LD-score indexes; ``legacy_converter`` imports legacy
LD-score files. This is tier 6 of the package layout: it may import the
reference-panel, annotation and summary-statistics domains below it.
"""
