"""PLINK genotype reader and in-memory LD-score block sums.

Moved verbatim from ``ldscore.py`` so the reader has a focused home; behavior
is unchanged. Selective-read and streaming optimizations build on this module.
"""
from __future__ import annotations

import numpy as np

try:  # pragma: no cover - optional dependency
    import bitarray as ba
except ImportError:  # pragma: no cover - optional dependency
    ba = None

from ..errors import LDSCDependencyError


class __GenotypeArrayInMemory__(object):
    """Parent class for in-memory genotype matrices."""

    def __init__(self, fname, n, snp_list, keep_snps=None, keep_indivs=None, mafMin=None):
        """Load, filter, and normalize one in-memory genotype matrix."""
        self.m = len(snp_list.IDList)
        self.n = n
        self.keep_snps = keep_snps
        self.keep_indivs = keep_indivs
        self.df = np.array(snp_list.df[["CHR", "SNP", "BP", "CM"]])
        self.colnames = ["CHR", "SNP", "POS", "CM"]
        self.mafMin = mafMin if mafMin is not None else 0
        self._currentSNP = 0
        (self.nru, self.geno) = self.__read__(fname, self.m, n)
        if keep_indivs is not None:
            keep_indivs = np.array(keep_indivs, dtype="int")
            if np.any(keep_indivs > self.n):
                raise ValueError("keep_indivs indices out of bounds")
            (self.geno, self.m, self.n) = self.__filter_indivs__(self.geno, keep_indivs, self.m, self.n)
            if self.n <= 0:
                raise ValueError("After filtering, no individuals remain")
        if keep_snps is not None:
            keep_snps = np.array(keep_snps, dtype="int")
            if np.any(keep_snps > self.m):
                raise ValueError("keep_snps indices out of bounds")
        (self.geno, self.m, self.n, self.kept_snps, self.freq) = self.__filter_snps_maf__(
            self.geno, self.m, self.n, self.mafMin, keep_snps
        )
        if self.m <= 0:
            raise ValueError("After filtering, no SNPs remain")
        self.df = self.df[self.kept_snps, :]
        self.maf = np.minimum(self.freq, np.ones(self.m) - self.freq)
        self.sqrtpq = np.sqrt(self.freq * (np.ones(self.m) - self.freq))
        self.df = np.c_[self.df, self.maf]
        self.colnames.append("MAF")

    def __read__(self, fname, m, n):
        """Read the backend-specific genotype representation into memory."""
        raise NotImplementedError

    def __filter_indivs__(geno, keep_indivs, m, n):
        """Apply backend-specific sample filtering to the genotype matrix."""
        raise NotImplementedError

    def __filter_maf_(geno, m, n, maf):
        """Apply backend-specific SNP and MAF filtering to the genotype matrix."""
        raise NotImplementedError

    def ldScoreVarBlocks(self, block_left, c, annot=None):
        """Compute LD-score block sums using the unbiased :math:`r^2` transform."""
        func = lambda x: self.__l2_unbiased__(x, self.n)
        snp_getter = self.nextSNPs
        return self.__corSumVarBlocks__(block_left, c, func, snp_getter, annot)

    def ldScoreBlockJackknife(self, block_left, c, annot=None, jN=10):
        """Compute block-jackknife LD-score summaries using squared correlations."""
        func = lambda x: np.square(x)
        snp_getter = self.nextSNPs
        return self.__corSumBlockJackknife__(block_left, c, func, snp_getter, annot, jN)

    def __l2_unbiased__(self, x, n):
        """Convert correlation values to the unbiased LD-score contribution."""
        denom = n - 2 if n > 2 else n
        sq = np.square(x)
        return sq - (1 - sq) / denom

    def __corSumVarBlocks__(self, block_left, c, func, snp_getter, annot=None):
        """Accumulate transformed correlation sums over LDSC-style LD blocks."""
        m, n = self.m, self.n
        block_sizes = np.array(np.arange(m) - block_left)
        block_sizes = np.ceil(block_sizes / c) * c
        if annot is None:
            annot = np.ones((m, 1))
        else:
            annot_m = annot.shape[0]
            if annot_m != self.m:
                raise ValueError("Incorrect number of SNPs in annot")

        n_a = annot.shape[1]
        cor_sum = np.zeros((m, n_a))
        b = np.nonzero(block_left > 0)
        if np.any(b):
            b = b[0][0]
        else:
            b = m
        b = int(np.ceil(b / c) * c)
        if b > m:
            c = 1
            b = m
        l_A = 0
        A = snp_getter(b)
        rfuncAB = np.zeros((b, c))
        rfuncBB = np.zeros((c, c))
        for l_B in range(0, b, c):
            B = A[:, l_B:l_B + c]
            np.dot(A.T, B / n, out=rfuncAB)
            rfuncAB = func(rfuncAB)
            cor_sum[l_A:l_A + b, :] += np.dot(rfuncAB, annot[l_B:l_B + c, :])
        b0 = b
        md = int(c * np.floor(m / c))
        end = md + 1 if md != m else md
        for l_B in range(b0, end, c):
            old_b = b
            b = int(block_sizes[l_B])
            if l_B > b0 and b > 0:
                A = np.hstack((A[:, old_b - b + c:old_b], B))
                l_A += old_b - b + c
            elif l_B == b0 and b > 0:
                A = A[:, b0 - b:b0]
                l_A = b0 - b
            elif b == 0:
                A = np.array(()).reshape((n, 0))
                l_A = l_B
            if l_B == md:
                c = m - md
                rfuncAB = np.zeros((b, c))
                rfuncBB = np.zeros((c, c))
            if b != old_b:
                rfuncAB = np.zeros((b, c))
            B = snp_getter(c)
            p1 = np.all(annot[l_A:l_A + b, :] == 0)
            p2 = np.all(annot[l_B:l_B + c, :] == 0)
            if p1 and p2:
                continue
            np.dot(A.T, B / n, out=rfuncAB)
            rfuncAB = func(rfuncAB)
            cor_sum[l_A:l_A + b, :] += np.dot(rfuncAB, annot[l_B:l_B + c, :])
            cor_sum[l_B:l_B + c, :] += np.dot(annot[l_A:l_A + b, :].T, rfuncAB).T
            np.dot(B.T, B / n, out=rfuncBB)
            rfuncBB = func(rfuncBB)
            cor_sum[l_B:l_B + c, :] += np.dot(rfuncBB, annot[l_B:l_B + c, :])
        return cor_sum


if ba is not None:
    class PlinkBEDFile(__GenotypeArrayInMemory__):
        """Interface for PLINK .bed format."""

        def __init__(self, fname, n, snp_list, keep_snps=None, keep_indivs=None, mafMin=None,
                     streaming=False):
            """Initialize the PLINK reader and configure BED bit-pattern decoding.

            When ``streaming`` is true the retained genotype bitarray is never
            materialized; ``nextSNPs`` decodes each batch directly from disk. Only
            ``build-ref-panel`` opts in, and only for unrestricted builds where the
            kept set spans the whole chromosome.
            """
            self._streaming = bool(streaming)
            self._bedcode = {
                2: ba.bitarray("11"),
                9: ba.bitarray("10"),
                1: ba.bitarray("01"),
                0: ba.bitarray("00"),
            }
            __GenotypeArrayInMemory__.__init__(
                self,
                fname,
                n,
                snp_list,
                keep_snps=keep_snps,
                keep_indivs=keep_indivs,
                mafMin=mafMin,
            )

        def __read__(self, fname, m, n):
            """Open the BED file, validate the header, and record geometry.

            The payload is decoded on demand from disk (see ``_source_snp_bits``),
            so the whole-chromosome bitarray is never materialized. The retained,
            individual-filtered bitarray that ``nextSNPs`` reads is built lazily by
            ``__filter_snps_maf__``.
            """
            if not fname.endswith(".bed"):
                raise ValueError(".bed filename must end in .bed")
            fh = open(fname, "rb")
            magicNumber = ba.bitarray(endian="little")
            magicNumber.fromfile(fh, 2)
            bedMode = ba.bitarray(endian="little")
            bedMode.fromfile(fh, 1)
            if magicNumber != ba.bitarray("0011011011011000"):
                raise IOError("Magic number from Plink .bed file not recognized")
            if bedMode != ba.bitarray("10000000"):
                raise IOError("Plink .bed file must be in default SNP-major mode")
            e = (4 - n % 4) if n % 4 != 0 else 0
            nru = n + e
            self._fh = fh
            self._data_start = 3
            self._nru_source = nru
            self._source_bytes_per_snp = 2 * nru // 8
            self._pending_keep_indivs = None
            fh.seek(0, 2)
            file_size = fh.tell()
            expected = self._data_start + m * self._source_bytes_per_snp
            if file_size != expected:
                raise IOError(
                    "Plink .bed file has {n1} bytes, expected {n2}".format(n1=file_size, n2=expected)
                )
            return (nru, None)

        def _source_snp_bits(self, source_index):
            """Return one SNP's bits, individual-subset when keep_indivs is active.

            Reads only this SNP's bytes from disk. With no individual filter the
            full source bits (length ``2 * nru_source``) are returned; otherwise the
            result has length ``2 * self.nru`` for the kept individuals, packed as
            the former ``__filter_indivs__`` did.
            """
            self._fh.seek(self._data_start + source_index * self._source_bytes_per_snp)
            raw = ba.bitarray(endian="little")
            raw.fromfile(self._fh, self._source_bytes_per_snp)
            keep = self._pending_keep_indivs
            if keep is None:
                return raw
            nru_src = self._nru_source
            z = ba.bitarray(2 * self.nru, endian="little")
            z.setall(0)
            for e, i in enumerate(keep):
                z[2 * e::2 * self.nru] = raw[2 * i::2 * nru_src]
                z[2 * e + 1::2 * self.nru] = raw[2 * i + 1::2 * nru_src]
            return z

        def __filter_indivs__(self, geno, keep_indivs, m, n):
            """Record the individual filter and compute kept-sample geometry.

            No genotype copy is made; individuals are subset per SNP during the
            on-demand read in ``__filter_snps_maf__``, so the raw and filtered
            bitarrays never coexist.
            """
            self._pending_keep_indivs = [int(i) for i in keep_indivs]
            n_new = len(self._pending_keep_indivs)
            e = (4 - n_new % 4) if n_new % 4 != 0 else 0
            self.nru = n_new + e
            return (geno, m, n_new)

        def __filter_snps_maf__(self, geno, m, n, mafMin, keep_snps):
            """Select SNPs by keep list and MAF, reading each SNP from disk.

            Each candidate SNP is read from source (and individual-subset) on
            demand, then counted. In the default (in-RAM) mode survivors are packed
            into the retained bitarray; in streaming mode only the kept-SNP/MAF
            bookkeeping is built and the genotypes stay on disk. Either way the
            whole-chromosome payload is never held in memory.
            """
            n_eff = n
            if keep_snps is None:
                keep_snps = range(m)
            m_poly = 0
            y = None if self._streaming else ba.bitarray()
            kept_snps = []
            freq = []
            for j in keep_snps:
                z = self._source_snp_bits(int(j))
                A = z[0::2]
                a = A.count()
                B = z[1::2]
                b = B.count()
                c = (A & B).count()
                major_ct = b + c
                n_nomiss = n_eff - a + c
                f = major_ct / (2 * n_nomiss) if n_nomiss > 0 else 0
                het_miss_ct = a + b - 2 * c
                if np.minimum(f, 1 - f) > mafMin and het_miss_ct < n_eff:
                    freq.append(f)
                    if not self._streaming:
                        y += z
                    m_poly += 1
                    kept_snps.append(int(j))
            return (y, m_poly, n_eff, kept_snps, freq)

        def nextSNPs(self, b, dtype=np.float64):
            """Return the next ``b`` standardized SNP columns from the BED stream.

            ``dtype`` selects the working precision of the decoded and standardized
            genotype matrix. It defaults to ``float64`` so the in-PLINK LD-score
            path is bit-for-bit unchanged; the reference-panel builder passes
            ``float32`` to halve the carry-over block's memory footprint.
            """
            try:
                b = int(b)
                if b <= 0:
                    raise ValueError("b must be > 0")
            except TypeError:
                raise TypeError("b must be an integer")
            if self._currentSNP + b > self.m:
                raise ValueError("{b} SNPs requested, {k} SNPs remain".format(b=b, k=(self.m - self._currentSNP)))
            c = self._currentSNP
            n = self.n
            nru = self.nru
            if self._streaming:
                snp_bits = ba.bitarray(endian="little")
                for k in range(c, c + b):
                    snp_bits += self._source_snp_bits(self.kept_snps[k])
            else:
                snp_bits = self.geno[2 * c * nru:2 * (c + b) * nru]
            X = np.array(list(snp_bits.decode(self._bedcode)), dtype=dtype).reshape((b, nru)).T
            X = X[0:n, :]
            Y = np.zeros(X.shape, dtype=dtype)
            for j in range(0, b):
                newsnp = X[:, j]
                ii = newsnp != 9
                avg = np.mean(newsnp[ii])
                newsnp[np.logical_not(ii)] = avg
                denom = np.std(newsnp)
                if denom == 0:
                    denom = 1
                Y[:, j] = (newsnp - avg) / denom
            self._currentSNP += b
            return Y
else:
    class PlinkBEDFile:  # pragma: no cover - dependency-gated fallback
        """Fallback PLINK reader that raises when `bitarray` is unavailable."""
        def __init__(self, *args, **kwargs):
            """Raise an informative import error for dependency-gated PLINK support."""
            raise LDSCDependencyError("PLINK LD-score support requires the optional dependency 'bitarray'.")
