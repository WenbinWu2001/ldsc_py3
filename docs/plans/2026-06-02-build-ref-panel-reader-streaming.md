# build-ref-panel genotype-reader streaming — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stop `build-ref-panel` from loading the whole-chromosome PLINK `.bed` into RAM — read only kept SNP blocks for restricted builds, and stream a sliding window for unrestricted builds — without changing any LD-score output or adding user flags.

**Architecture:** Extract the shared genotype reader into its own module, then make `__read__` lazy and `__filter_snps_maf__` seek-and-decode per kept SNP (folding individual filtering in). Add an opt-in `streaming` mode that determines kept SNPs by a bit-count pass and serves `nextSNPs` directly from the file. `build-ref-panel` opts into streaming only when no SNP restriction was supplied; every other consumer keeps the in-RAM default.

**Tech Stack:** Python 3.13, numpy, pandas, pyarrow, `bitarray`; pytest.

**Reference spec:** `docs/superpowers/specs/2026-06-02-build-ref-panel-reader-streaming-design.md`

**Environment for every command:**
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev
```
All paths are relative to the repo root
`/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured`.

**Bit-for-bit rule:** Tasks 3–5 must not change the values asserted by the golden tests written in Task 2. If a golden test fails, the change is wrong — fix the change, never the golden.

---

## File structure

| File | Responsibility | Task |
|---|---|---|
| `src/ldsc/_kernel/plink_bed.py` (create) | The PLINK genotype reader (`__GenotypeArrayInMemory__`, `PlinkBEDFile`, `bitarray`-absent fallback) | 1 |
| `src/ldsc/_kernel/ldscore.py` (modify) | Re-export the reader; keep coordinate helpers + LD-score workflow | 1, 5 |
| `tests/fixtures/golden/gen_reader_golden.py` (create) | One-off script to snapshot current reader/workflow outputs | 2 |
| `tests/fixtures/golden/reader_golden.npz` (create) | Committed golden arrays | 2 |
| `tests/test_plink_bed_reader.py` (create) | Golden bit-for-bit + mode-equivalence + filter-matrix tests | 2, 4 |
| `src/ldsc/ref_panel_builder.py` (modify) | Pass `streaming=` at the reader construction site | 4 |

---

## Task 1: Extract the genotype reader (pure move, zero logic change)

**Files:**
- Create: `src/ldsc/_kernel/plink_bed.py`
- Modify: `src/ldsc/_kernel/ldscore.py` (remove the two classes at lines 363–639; add a re-export import)

- [ ] **Step 1: Establish the green baseline**

Run:
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest -q
```
Expected: all tests pass. Record the pass count.

- [ ] **Step 2: Create `plink_bed.py` and move the reader classes into it**

Cut these spans from `src/ldsc/_kernel/ldscore.py` and paste them verbatim into a new `src/ldsc/_kernel/plink_bed.py`:
- `class __GenotypeArrayInMemory__` (currently `ldscore.py:363-494`)
- the `if ba is not None: class PlinkBEDFile(...) ... else: class PlinkBEDFile:` block (currently `ldscore.py:497-639`)

Add the imports the moved code needs to the top of `plink_bed.py`:
```python
"""PLINK genotype reader and in-memory LD-score block sums.

Moved verbatim from ``ldscore.py``; behavior is unchanged. Streaming and
selective-read optimizations land in later tasks.
"""
from __future__ import annotations

import numpy as np

try:
    import bitarray as ba
except ImportError:  # pragma: no cover - dependency-gated
    ba = None

from .errors import LDSCDependencyError  # adjust to the actual LDSCDependencyError import path
```

Verify the actual import location of `LDSCDependencyError` and any other names the moved code references (`np`, `ba`) by grepping `ldscore.py`'s header:
```bash
grep -nE "LDSCDependencyError|^import |^from |bitarray as ba" src/ldsc/_kernel/ldscore.py | head
```
Match `plink_bed.py`'s imports to whatever `ldscore.py` actually uses.

- [ ] **Step 3: Re-export from `ldscore.py`**

At the point where the classes were removed (or with the other intra-package imports near the top of `ldscore.py`), add:
```python
from .plink_bed import __GenotypeArrayInMemory__, PlinkBEDFile
```
This keeps `sys.modules["ldsc._kernel.ldscore"].PlinkBEDFile` resolvable, so `get_legacy_ld_module()` (`ldscore.py:301`) and `compute_chrom_from_plink`'s `legacy_ld.PlinkBEDFile` are unaffected. Leave `getBlockLefts`, `get_block_lefts`, and `block_left_to_right` in `ldscore.py`.

- [ ] **Step 4: Verify nothing changed**

Run:
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest -q
```
Expected: identical pass count to Step 1, zero failures.

Also confirm the public import surface is intact:
```bash
python -c "from ldsc._kernel.ldscore import PlinkBEDFile; from ldsc._kernel.plink_bed import PlinkBEDFile as P2; print(PlinkBEDFile is P2)"
```
Expected: `True`.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/plink_bed.py src/ldsc/_kernel/ldscore.py
git commit -m "refactor(kernel): extract PLINK genotype reader to plink_bed.py"
```

---

## Task 2: Golden bit-for-bit snapshot (safety net before any behavior change)

**Files:**
- Create: `tests/fixtures/golden/gen_reader_golden.py`
- Create: `tests/fixtures/golden/reader_golden.npz`
- Create: `tests/test_plink_bed_reader.py`

This task captures the *current* reader and LD-score outputs as committed golden
arrays, then writes tests that assert against them. The tests pass now and must
keep passing through Tasks 3–5.

- [ ] **Step 1: Write the snapshot generator**

Create `tests/fixtures/golden/gen_reader_golden.py`:
```python
"""Regenerate reader/workflow golden arrays. Run only on a known-good commit:

    python tests/fixtures/golden/gen_reader_golden.py
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


def build(keep_snps=None, keep_indivs=None, maf_min=None):
    bim = legacy_parse.PlinkBIMFile(str(FIX) + ".bim")
    fam = legacy_parse.PlinkFAMFile(str(FIX) + ".fam")
    return PlinkBEDFile(
        str(FIX) + ".bed", len(fam.IDList), bim,
        keep_snps=keep_snps, keep_indivs=keep_indivs, mafMin=maf_min,
    )


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
```

- [ ] **Step 2: Generate and commit the golden artifact**

Run (on the current, known-good code):
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && python tests/fixtures/golden/gen_reader_golden.py
```
Expected: prints `wrote tests/fixtures/golden/reader_golden.npz m=<M> n=<N>` with the fixture's SNP/sample counts. Confirm the file exists:
```bash
ls -l tests/fixtures/golden/reader_golden.npz
```

- [ ] **Step 3: Write the golden test**

Create `tests/test_plink_bed_reader.py`:
```python
import importlib.util
from pathlib import Path

import numpy as np
import pytest

_HAS_BA = importlib.util.find_spec("bitarray") is not None
pytestmark = pytest.mark.skipif(not _HAS_BA, reason="bitarray required")

from ldsc._kernel import formats as legacy_parse
from ldsc._kernel.plink_bed import PlinkBEDFile

FIX = Path("tests/fixtures/plink/plink")
GOLDEN = Path("tests/fixtures/golden/reader_golden.npz")


def _build(keep_snps=None, keep_indivs=None, maf_min=None):
    bim = legacy_parse.PlinkBIMFile(str(FIX) + ".bim")
    fam = legacy_parse.PlinkFAMFile(str(FIX) + ".fam")
    return PlinkBEDFile(
        str(FIX) + ".bed", len(fam.IDList), bim,
        keep_snps=keep_snps, keep_indivs=keep_indivs, mafMin=maf_min,
    )


def _decode_all(geno):
    geno._currentSNP = 0
    cols = [geno.nextSNPs(1, dtype=np.float64).ravel() for _ in range(geno.m)]
    return np.column_stack(cols) if cols else np.zeros((geno.n, 0))


def test_in_ram_reader_matches_golden():
    g = np.load(GOLDEN)
    geno = _build()
    assert int(geno.m) == int(g["m"])
    assert int(geno.n) == int(g["n"])
    np.testing.assert_array_equal(np.asarray(geno.kept_snps), g["kept_snps"])
    np.testing.assert_allclose(np.asarray(geno.freq), g["freq"], rtol=0, atol=0)
    np.testing.assert_allclose(np.asarray(geno.maf), g["maf"], rtol=0, atol=0)
    np.testing.assert_allclose(_decode_all(geno), g["decoded"], rtol=0, atol=0)
```

- [ ] **Step 4: Run the golden test**

Run:
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest tests/test_plink_bed_reader.py -v
```
Expected: `test_in_ram_reader_matches_golden` PASSES (asserts current behavior against the just-captured golden).

- [ ] **Step 5: Add a build-ref-panel parquet golden test**

Append to `tests/test_plink_bed_reader.py`. This runs an unrestricted `build-ref-panel` on the fixture and snapshots the parquet on first run, then asserts equality on subsequent runs (so it guards Tasks 3–5):
```python
def _run_build_ref_panel(tmp_path, prefix):
    from ldsc.ref_panel_builder import ReferencePanelBuilder, ReferencePanelBuildConfig
    from ldsc.config import GlobalConfig
    cfg = ReferencePanelBuildConfig(
        plink_prefix=str(prefix),
        source_genome_build="hg38",
        output_dir=str(tmp_path / "panel"),
        overwrite=True,
        ld_wind_snps=8,
    )
    ReferencePanelBuilder(cfg, GlobalConfig(snp_identifier="rsid")).run()
    import glob
    paths = sorted(glob.glob(str(tmp_path / "panel" / "**" / "*_r2.parquet"), recursive=True))
    assert paths, "no parquet produced"
    return paths[0]


def _read_r2(parquet_path):
    import pyarrow.parquet as pq
    t = pq.read_table(parquet_path)
    return {name: t.column(name).to_numpy() for name in t.column_names}


def test_build_ref_panel_parquet_golden(tmp_path):
    gp = Path("tests/fixtures/golden/build_ref_panel_r2.npz")
    out = _read_r2(_run_build_ref_panel(tmp_path, FIX))
    if not gp.exists():
        np.savez(gp, **out)
        pytest.skip("captured build-ref-panel golden; rerun to assert")
    g = np.load(gp)
    for name in out:
        np.testing.assert_array_equal(out[name], g[name], err_msg=f"column {name} drifted")
```
Adjust `ReferencePanelBuildConfig`/`GlobalConfig` keyword names to the real signatures (grep `src/ldsc/config.py`). If the fixture lacks the columns `build-ref-panel` needs (e.g. allele info), use `tests/fixtures/minimal_external_resources/plink/hm3_chr22_subset` instead and set `snp_identifier`/window accordingly.

- [ ] **Step 6: Capture then assert the parquet golden**

Run twice:
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev \
  && pytest tests/test_plink_bed_reader.py::test_build_ref_panel_parquet_golden -v \
  && pytest tests/test_plink_bed_reader.py::test_build_ref_panel_parquet_golden -v
```
Expected: first run SKIPS (captures `build_ref_panel_r2.npz`), second run PASSES.

- [ ] **Step 7: Commit**

```bash
git add tests/fixtures/golden/ tests/test_plink_bed_reader.py
git commit -m "test(kernel): golden snapshot of PLINK reader and build-ref-panel R2"
```

---

## Task 3: O1 + O3 — lazy, selective read with fused individual filter (in-RAM mode)

**Files:**
- Modify: `src/ldsc/_kernel/plink_bed.py` (`PlinkBEDFile.__read__`, `__filter_snps_maf__`; the base `__init__` filtering order)
- Test: `tests/test_plink_bed_reader.py`

Goal: never materialize the whole-chromosome bitarray for restricted builds.
`nextSNPs` is unchanged (still reads the small retained `self.geno`). MAF must be
computed on kept individuals, preserving today's order.

- [ ] **Step 1: Write the failing filter-matrix test**

Append to `tests/test_plink_bed_reader.py`:
```python
@pytest.mark.parametrize("maf_min", [0.0, 0.05])
@pytest.mark.parametrize("with_keep_indivs", [False, True])
def test_selective_read_matches_whole_read(maf_min, with_keep_indivs):
    full = _build()  # decode reference via current path semantics
    keep_snps = list(range(0, int(full.m), 2)) or [0]
    keep_indivs = list(range(0, max(1, full.n - 1))) if with_keep_indivs else None
    geno = _build(keep_snps=keep_snps, keep_indivs=keep_indivs, maf_min=maf_min)
    # Reader must still expose a usable, finite standardized matrix and contract.
    assert geno.m >= 0 and geno.n >= 1
    if geno.m:
        decoded = _decode_all(geno)
        assert decoded.shape == (geno.n, geno.m)
        assert np.isfinite(decoded).all()
    assert hasattr(geno, "kept_snps") and hasattr(geno, "freq")
    assert hasattr(geno, "maf") and hasattr(geno, "sqrtpq")
    assert geno.df.shape[0] == geno.m and "MAF" in geno.colnames


def test_no_whole_chromosome_bitarray_after_init():
    # After a restricted init, the retained bitarray must hold only kept SNPs,
    # not the whole chromosome: len(geno.geno) == 2 * nru * m.
    keep_snps = [0]
    geno = _build(keep_snps=keep_snps)
    assert len(geno.geno) == 2 * geno.nru * geno.m
```

- [ ] **Step 2: Run to confirm the contract tests pass on current code, then drive the change with the golden**

Run:
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest tests/test_plink_bed_reader.py -v
```
Expected: all PASS on current code (these assert the contract O1 must preserve). The real regression guard for O1 is `test_in_ram_reader_matches_golden`.

- [ ] **Step 3: Make `__read__` lazy**

In `plink_bed.py`, replace `PlinkBEDFile.__read__` (currently reads the whole payload) with a version that opens the handle, validates the header, stores geometry, and returns no payload:
```python
def __read__(self, fname, m, n):
    """Open the BED file, validate the header, and record geometry.

    The payload is no longer loaded here; SNPs are decoded on demand from the
    retained, individual-filtered bitarray built in __filter_snps_maf__.
    """
    if not fname.endswith(".bed"):
        raise ValueError(".bed filename must end in .bed")
    fh = open(fname, "rb")
    magicNumber = ba.bitarray(endian="little")
    magicNumber.fromfile(fh, 2)
    bedMode = ba.bitarray(endian="little")
    bedMode.fromfile(fh, 1)
    e = (4 - n % 4) if n % 4 != 0 else 0
    nru = n + e
    self.nru = nru
    if magicNumber != ba.bitarray("0011011011011000"):
        raise IOError("Magic number from Plink .bed file not recognized")
    if bedMode != ba.bitarray("10000000"):
        raise IOError("Plink .bed file must be in default SNP-major mode")
    self._fh = fh
    self._data_start = 3                 # bytes: 2 magic + 1 mode
    self._nru_source = nru               # source padding (full individuals)
    self._source_bytes_per_snp = 2 * nru // 8
    return (nru, None)                   # payload built later; self.geno stays None for now
```

- [ ] **Step 4: Add a per-SNP source reader helper**

Add to `PlinkBEDFile`:
```python
def _source_snp_bits(self, source_index, keep_indivs):
    """Return one SNP's bits, individual-subset when keep_indivs is given.

    Reads only this SNP's bytes from the file. When keep_indivs is None the
    full source bits (length 2*nru_source) are returned; otherwise the result
    has length 2*self.nru for the kept individuals, packed exactly as the old
    __filter_indivs__ did.
    """
    self._fh.seek(self._data_start + source_index * self._source_bytes_per_snp)
    raw = ba.bitarray(endian="little")
    raw.fromfile(self._fh, self._source_bytes_per_snp)
    if keep_indivs is None:
        return raw
    nru_src = self._nru_source
    z = ba.bitarray(2 * self.nru, endian="little")
    z.setall(0)
    for e, i in enumerate(keep_indivs):
        z[2 * e::2 * self.nru] = raw[2 * i::2 * nru_src]
        z[2 * e + 1::2 * self.nru] = raw[2 * i + 1::2 * nru_src]
    return z
```

- [ ] **Step 5: Fold individual filtering into `__filter_snps_maf__` and read per SNP**

Replace `PlinkBEDFile.__filter_snps_maf__` so it (a) iterates source SNP indices, (b) reads each SNP via `_source_snp_bits` with the kept individuals applied, (c) computes MAF on the kept-individual bits using the exact current arithmetic, and (d) packs survivors into `self.geno`:
```python
def __filter_snps_maf__(self, geno, m, n, mafMin, keep_snps):
    """Select SNPs by keep list and MAF, reading each SNP from disk on demand.

    Individual filtering is fused here: each SNP is read from source and
    subset to the kept individuals before counting, so the whole-chromosome
    bitarray and a second filtered copy are never materialized.
    """
    keep_indivs = self._pending_keep_indivs   # set by __init__ (see Step 6)
    if keep_snps is None:
        keep_snps = range(m)
    nru = self.nru
    n_eff = n
    m_poly = 0
    y = ba.bitarray()
    kept_snps = []
    freq = []
    for j in keep_snps:
        z = self._source_snp_bits(int(j), keep_indivs)  # length 2*nru
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
            y += z
            m_poly += 1
            kept_snps.append(int(j))
    return (y, m_poly, n_eff, kept_snps, freq)
```
Notes:
- `nru` is already `nru_new` (kept padded) because `__init__` sets it before this call (Step 6).
- The arithmetic (`major_ct`, `n_nomiss`, `f`, `het_miss_ct`, the `> mafMin` and `< n` tests) is copied unchanged from the current implementation to preserve drop decisions and `freq` exactly.

- [ ] **Step 6: Adjust the base `__init__` to set geometry before filtering and drop the standalone individual-filter pass**

In `__GenotypeArrayInMemory__.__init__` (now in `plink_bed.py`), compute the kept-individual geometry up front and stash the kept list for `__filter_snps_maf__`, replacing the `__filter_indivs__` call:
```python
    self._currentSNP = 0
    (self.nru, self.geno) = self.__read__(fname, self.m, n)   # nru = source padding; geno is None
    self._pending_keep_indivs = None
    if keep_indivs is not None:
        keep_indivs = np.array(keep_indivs, dtype="int")
        if np.any(keep_indivs > self.n):
            raise ValueError("keep_indivs indices out of bounds")
        self._pending_keep_indivs = [int(i) for i in keep_indivs]
        n_new = len(self._pending_keep_indivs)
        e = (4 - n_new % 4) if n_new % 4 != 0 else 0
        self.nru = n_new + e          # nru_new used by _source_snp_bits packing
        self.n = n_new
        if self.n <= 0:
            raise ValueError("After filtering, no individuals remain")
    if keep_snps is not None:
        keep_snps = np.array(keep_snps, dtype="int")
        if np.any(keep_snps > self.m):
            raise ValueError("keep_snps indices out of bounds")
    (self.geno, self.m, self.n, self.kept_snps, self.freq) = self.__filter_snps_maf__(
        self.geno, self.m, self.n, self.mafMin, keep_snps
    )
```
Leave the post-filter block (`self.df = self.df[self.kept_snps, :]`, `maf`, `sqrtpq`,
`df` MAF column, `colnames.append("MAF")`) unchanged. `__filter_indivs__` is now
unused for `PlinkBEDFile`; leave the base stub as-is (it is still part of the
abstract contract) but it is no longer called.

- [ ] **Step 7: Confirm `nextSNPs` is untouched and reads from the retained `self.geno`**

No change to `nextSNPs`: after Step 5 `self.geno` is the retained, kept-individual,
kept-SNP bitarray with `self.nru` = kept padding — exactly the shape `nextSNPs`
already expects. Verify the slice math (`self.geno[2*c*nru : 2*(c+b)*nru]`) still
holds with `nru = self.nru`.

- [ ] **Step 8: Run golden + matrix tests**

Run:
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest tests/test_plink_bed_reader.py -v
```
Expected: `test_in_ram_reader_matches_golden`, `test_build_ref_panel_parquet_golden`, `test_selective_read_matches_whole_read[*]`, and `test_no_whole_chromosome_bitarray_after_init` all PASS.

- [ ] **Step 9: Run the full suite (LD-score path regression)**

Run:
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest -q
```
Expected: all pass, including `tests/test_ldscore_workflow.py` and `tests/test_ref_panel*.py` (these exercise `compute_chrom_from_plink` and `PlinkRefPanel`, proving the LD-score and API paths are unaffected).

- [ ] **Step 10: Commit**

```bash
git add src/ldsc/_kernel/plink_bed.py tests/test_plink_bed_reader.py
git commit -m "perf(kernel): selective per-SNP BED read with fused individual filter"
```

---

## Task 4: O2 — streaming mode (opt-in; build-ref-panel unrestricted only)

**Files:**
- Modify: `src/ldsc/_kernel/plink_bed.py` (`PlinkBEDFile.__init__` signature, pass-1 bit-count, streaming `nextSNPs`)
- Modify: `src/ldsc/ref_panel_builder.py:802` (pass `streaming=`)
- Test: `tests/test_plink_bed_reader.py`

Goal: for unrestricted builds, determine kept SNPs by a no-retention bit-count
pass and serve `nextSNPs` from the file, so `self.geno` is never packed.

- [ ] **Step 1: Write the mode-equivalence failing test**

Append to `tests/test_plink_bed_reader.py`:
```python
def test_streaming_matches_in_ram_unrestricted():
    in_ram = _build()                       # streaming=False default
    streamed = _build_streaming()           # helper defined below
    assert int(streamed.m) == int(in_ram.m)
    assert int(streamed.n) == int(in_ram.n)
    np.testing.assert_array_equal(np.asarray(streamed.kept_snps),
                                  np.asarray(in_ram.kept_snps))
    np.testing.assert_allclose(np.asarray(streamed.freq),
                               np.asarray(in_ram.freq), rtol=0, atol=0)
    np.testing.assert_allclose(_decode_all(streamed), _decode_all(in_ram),
                               rtol=0, atol=0)


def _build_streaming(keep_indivs=None, maf_min=None):
    bim = legacy_parse.PlinkBIMFile(str(FIX) + ".bim")
    fam = legacy_parse.PlinkFAMFile(str(FIX) + ".fam")
    return PlinkBEDFile(
        str(FIX) + ".bed", len(fam.IDList), bim,
        keep_snps=None, keep_indivs=keep_indivs, mafMin=maf_min, streaming=True,
    )
```

- [ ] **Step 2: Run to verify it fails**

Run:
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest tests/test_plink_bed_reader.py::test_streaming_matches_in_ram_unrestricted -v
```
Expected: FAIL with `TypeError: __init__() got an unexpected keyword argument 'streaming'`.

- [ ] **Step 3: Add the `streaming` flag and branch the filter pass**

In `plink_bed.py`, thread `streaming` through both `__init__`s:
```python
# __GenotypeArrayInMemory__.__init__ signature:
def __init__(self, fname, n, snp_list, keep_snps=None, keep_indivs=None,
             mafMin=None, streaming=False):
    ...
    self._streaming = bool(streaming)
    ...
# PlinkBEDFile.__init__ signature mirrors it and forwards streaming to super().__init__.
```
In `__filter_snps_maf__`, when `self._streaming` is true, do the bit-count pass
without packing genotypes and record the kept→source map:
```python
def __filter_snps_maf__(self, geno, m, n, mafMin, keep_snps):
    keep_indivs = self._pending_keep_indivs
    if keep_snps is None:
        keep_snps = range(m)
    n_eff = self.n
    m_poly = 0
    kept_snps = []
    freq = []
    y = None if self._streaming else ba.bitarray()
    for j in keep_snps:
        z = self._source_snp_bits(int(j), keep_indivs)
        A = z[0::2]; a = A.count()
        B = z[1::2]; b = B.count()
        c = (A & B).count()
        major_ct = b + c
        n_nomiss = n_eff - a + c
        f = major_ct / (2 * n_nomiss) if n_nomiss > 0 else 0
        het_miss_ct = a + b - 2 * c
        if np.minimum(f, 1 - f) > mafMin and het_miss_ct < n_eff:
            freq.append(f)
            kept_snps.append(int(j))
            m_poly += 1
            if not self._streaming:
                y += z
    return (y, m_poly, n_eff, kept_snps, freq)
```
`kept_snps` already serves as the kept→source index map (source indices, in kept order).

- [ ] **Step 4: Implement the streaming `nextSNPs` branch**

Refactor `nextSNPs` so the decode/standardize logic is shared and only the byte
source differs:
```python
def nextSNPs(self, b, minorRef=None, dtype=np.float64):
    """Return the next b standardized SNP columns (in-RAM or streaming)."""
    b = int(b)
    if b <= 0:
        raise ValueError("b must be > 0")
    if self._currentSNP + b > self.m:
        raise ValueError(f"{b} SNPs requested, {self.m - self._currentSNP} remain")
    nru = self.nru
    n = self.n
    if self._streaming:
        slice_bits = ba.bitarray(endian="little")
        for k in range(self._currentSNP, self._currentSNP + b):
            slice_bits += self._source_snp_bits(self.kept_snps[k], self._pending_keep_indivs)
    else:
        c = self._currentSNP
        slice_bits = self.geno[2 * c * nru:2 * (c + b) * nru]
    X = np.array(list(slice_bits.decode(self._bedcode)), dtype=dtype).reshape((b, nru)).T
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
        if minorRef is not None and self.freq[self._currentSNP + j] > 0.5:
            denom = denom * -1
        Y[:, j] = (newsnp - avg) / denom
    self._currentSNP += b
    return Y
```
This is the current `nextSNPs` body with only the `slice_bits` source branched; the
decode and standardization are byte-identical, so in-RAM and streaming produce the
same values. Rewind works because `nextSNPs` re-reads `self.kept_snps[k]` from
`_currentSNP`.

- [ ] **Step 5: Run the mode-equivalence test**

Run:
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest tests/test_plink_bed_reader.py::test_streaming_matches_in_ram_unrestricted -v
```
Expected: PASS.

- [ ] **Step 6: Wire the trigger in build-ref-panel**

In `src/ldsc/ref_panel_builder.py` at the `PlinkBEDFile(...)` construction
(currently line 802), pass `streaming` computed from whether a SNP restriction was
resolved. Use the resolved build-state restriction (preferred) or the config:
```python
            unrestricted = config.ref_panel_snps_file is None and not config.use_hm3_snps
            geno = kernel_ldscore.PlinkBEDFile(
                prefix + ".bed",
                len(fam.IDList),
                bim,
                keep_snps=list(build_keep_snps),
                keep_indivs=None if keep_indivs is None else list(keep_indivs),
                mafMin=config.maf_min,
                streaming=unrestricted,
            )
```
If a resolved `restriction_values`/`restriction_keys` is available on the build
state in this scope, prefer `unrestricted = build_state.restriction_values is None`
(grep the surrounding method to confirm the variable name; `ref_panel_builder.py:130,646`).

- [ ] **Step 7: Run the parquet golden under streaming**

The build-ref-panel golden from Task 2 was captured unrestricted, so it now exercises
the streaming path. Run:
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest tests/test_plink_bed_reader.py -v
```
Expected: `test_build_ref_panel_parquet_golden` PASSES (streaming output equals the
pre-refactor parquet), all other reader tests PASS.

- [ ] **Step 8: Run the full suite**

Run:
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest -q
```
Expected: all pass.

- [ ] **Step 9: Commit**

```bash
git add src/ldsc/_kernel/plink_bed.py src/ldsc/ref_panel_builder.py tests/test_plink_bed_reader.py
git commit -m "perf(build-ref-panel): stream genotypes for unrestricted builds"
```

---

## Task 5: Targeted cleanup — remove the dead `minorRef`/`freq`-in-`nextSNPs` branch

**Files:**
- Modify: `src/ldsc/_kernel/plink_bed.py`
- Test: `tests/test_plink_bed_reader.py`

No consumer passes `minorRef`; the branch reading `self.freq` inside `nextSNPs` is
dead. Remove it to keep the decode path clean.

- [ ] **Step 1: Confirm `minorRef` is unused**

Run:
```bash
grep -rn --include='*.py' "minorRef" src/ldsc
```
Expected: only the `nextSNPs` definition and its internal branch — no call site
passes it.

- [ ] **Step 2: Remove the `minorRef` parameter and branch**

In `nextSNPs`, drop the `minorRef=None` parameter and delete:
```python
        if minorRef is not None and self.freq[self._currentSNP + j] > 0.5:
            denom = denom * -1
```
Keep everything else identical.

- [ ] **Step 3: Run reader + full suite**

Run:
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest -q
```
Expected: all pass; golden values unchanged.

- [ ] **Step 4: Commit**

```bash
git add src/ldsc/_kernel/plink_bed.py
git commit -m "refactor(kernel): drop dead minorRef branch from nextSNPs"
```

---

## Self-review notes (author)

- **Spec coverage:** §3.1 → Task 1; §6 golden/equivalence → Task 2 + Task 4; §3.2 (O1) + O3 → Task 3; §3.3 (O2) + trigger → Task 4; §3.4 single-class shape → Tasks 3–4 (two `if self._streaming` sites); cleanup → Task 5. LD-score regression coverage → Task 3 Step 9 / Task 4 Step 8 (full suite).
- **Type/name consistency:** `_source_snp_bits`, `_pending_keep_indivs`, `_streaming`, `_nru_source`, `_source_bytes_per_snp`, `_data_start`, `_fh` are introduced in Task 3 and reused in Task 4. `streaming` kwarg name matches the trigger in Task 4 Step 6.
- **Known follow-ups (out of scope):** separating the LD-score math from the reader class; compact `.bim` parse (O5); chr-only genetic map (O4); columnar pair emission (O6).
- **Verification caveats flagged for the implementer:** exact `LDSCDependencyError` import path (Task 1 Step 2); exact `ReferencePanelBuildConfig`/`GlobalConfig` kwargs and fixture column adequacy (Task 2 Step 5); exact build-state restriction variable name (Task 4 Step 6).
