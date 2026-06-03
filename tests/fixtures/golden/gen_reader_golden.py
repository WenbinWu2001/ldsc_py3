"""Regenerate reader golden arrays. Run only on a known-good commit:

    python tests/fixtures/golden/gen_reader_golden.py

Captures the in-memory PLINK reader's decoded standardized genotypes plus the
retained-SNP bookkeeping for the small ``tests/fixtures/plink`` fixture. The
selective-read and streaming refactors must reproduce these arrays exactly.
"""
from pathlib import Path

import numpy as np

from ldsc._kernel import formats as legacy_parse
from ldsc._kernel.plink_bed import PlinkBEDFile

FIX = Path("tests/fixtures/plink/plink")
OUT = Path("tests/fixtures/golden/reader_golden.npz")


def decode_all(geno):
    geno._currentSNP = 0
    cols = [geno.nextSNPs(1, dtype=np.float64).ravel() for _ in range(geno.m)]
    return np.column_stack(cols) if cols else np.zeros((geno.n, 0))


KEEP_INDIVS = [0, 1, 2, 3]


def build(keep_indivs=None):
    bim = legacy_parse.PlinkBIMFile(str(FIX) + ".bim")
    fam = legacy_parse.PlinkFAMFile(str(FIX) + ".fam")
    return PlinkBEDFile(str(FIX) + ".bed", len(fam.IDList), bim, keep_indivs=keep_indivs)


def main():
    base = build()
    ki = build(keep_indivs=KEEP_INDIVS)
    np.savez(
        OUT,
        decoded=decode_all(base),
        kept_snps=np.asarray(base.kept_snps),
        freq=np.asarray(base.freq, dtype=np.float64),
        maf=np.asarray(base.maf, dtype=np.float64),
        m=np.asarray(base.m),
        n=np.asarray(base.n),
        ki_decoded=decode_all(ki),
        ki_kept_snps=np.asarray(ki.kept_snps),
        ki_freq=np.asarray(ki.freq, dtype=np.float64),
        ki_m=np.asarray(ki.m),
        ki_n=np.asarray(ki.n),
    )
    print(f"wrote {OUT} base m={base.m} n={base.n}; keep_indivs m={ki.m} n={ki.n}")


if __name__ == "__main__":
    main()
