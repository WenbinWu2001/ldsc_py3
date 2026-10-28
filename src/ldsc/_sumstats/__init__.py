"""GWAS summary statistics: reading, munging and loading curated tables.

Raw-file format inference, the munging workflow and its quality-control
kernel, liftover of summary-statistic coordinates, and the loader for
curated summary statistics that regression reads. This is tier 5 of the
package layout, a peer of ``ldsc._annotation`` and ``ldsc._refpanel``;
peers do not import each other.
"""
