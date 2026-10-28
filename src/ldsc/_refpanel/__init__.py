"""Reference panels: PLINK genotypes and package-built R2 matrices.

Reading PLINK data, the R2 matrix file format, LD-window geometry, the
runtime adapters that prepare one chromosome for LD-score computation, R2
matrix building and pair queries. This is tier 5 of the package layout, a
peer of ``ldsc._annotation`` and ``ldsc._sumstats``; peers do not import
each other.
"""
