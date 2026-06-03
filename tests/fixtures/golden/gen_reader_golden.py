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


def build():
    bim = legacy_parse.PlinkBIMFile(str(FIX) + ".bim")
    fam = legacy_parse.PlinkFAMFile(str(FIX) + ".fam")
    return PlinkBEDFile(str(FIX) + ".bed", len(fam.IDList), bim)


def main():
    geno = build()
    np.savez(
        OUT,
        decoded=decode_all(geno),
        kept_snps=np.asarray(geno.kept_snps),
        freq=np.asarray(geno.freq, dtype=np.float64),
        maf=np.asarray(geno.maf, dtype=np.float64),
        m=np.asarray(geno.m),
        n=np.asarray(geno.n),
    )
    print(f"wrote {OUT} m={geno.m} n={geno.n}")


if __name__ == "__main__":
    main()
