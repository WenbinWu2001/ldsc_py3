"""Exact gene LD-score indexes: building, the on-disk format and scoring.

``build`` builds and publishes an index from a PLINK reference panel and a
gene catalog; ``artifact`` owns the on-disk format and identity hashes;
``scoring`` loads an index and scores gene lists with it; ``atoms`` holds
the sparse disjoint-atom primitives; ``batches`` runs indexed query batches
in bounded workers. See ``docs/reference/gene-ldscore-index.md``.
"""
