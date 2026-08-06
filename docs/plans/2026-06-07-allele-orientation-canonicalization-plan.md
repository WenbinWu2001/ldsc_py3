# Allele Orientation Canonicalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `A1` the minor allele in every non-sumstats artifact (so `MAF = freq(A1) <= 0.5`) and give the reference-panel `SIGN`/`query-r2 --with-r` a stable minor-allele orientation, while leaving sumstats (`A1` = effect allele) untouched.

**Architecture:** Approach A from the design doc — the flip happens at `build-ref-panel` standardization: per-SNP, where the `.bim` `A2` frequency `< 0.5`, swap `A1`/`A2` in the meta sidecar and negate the standardized genotype column before the correlation/`SIGN` computation. `nextSNPs` stays orientation-free. R²/LD-scores/h2/rg are sign-invariant and unchanged. A shared `canonicalize_alleles()` helper normalizes externally-supplied frequency sidecars on read; the column registry gains an oriented `FRQ`-family vs folded `MAF` distinction.

**Tech Stack:** Python 3, NumPy, pandas, pyarrow, pytest. Spec: `docs/superpowers/specs/2026-06-07-allele-orientation-canonicalization-design.md`.

**Environment for every command below:**
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev
cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured
```
Confirm branch is `restructure` before starting: `git branch --show-current`.

---

## File structure

| File | Responsibility | Action |
|---|---|---|
| `src/ldsc/_kernel/ref_panel_builder.py` | `orientation_flip_sign`, oriented-getter wrapper, `A1/A2` swap in metadata frame | Modify |
| `src/ldsc/ref_panel_builder.py` | Compute flip from `geno.freq`; wire swap + oriented getter into the build loop | Modify |
| `src/ldsc/_kernel/snp_identity.py` | `canonicalize_alleles()` defensive-read helper | Modify |
| `src/ldsc/column_inference.py` | Split oriented `FRQ`-family from folded `MAF` for non-sumstats sidecars | Modify |
| `src/ldsc/ldscore_calculator.py` | Assert `MAF <= 0.5` invariant on emitted LD-score metadata | Modify |
| `src/ldsc/r2_query.py` | Document signed r as minor-allele (`A1`) orientation | Modify |
| `src/ldsc/_kernel/sumstats_munger.py` | One-line carve-out comment (no behavior change) | Modify |
| `docs/current/column-schema.md` | Per-artifact `A1`/`A2`/`FRQ`/`MAF` definitions, sign rule, oriented-vs-folded | Modify |
| `docs/current/data-flow.md` | Note build-stage flip point; `nextSNPs` orientation-free | Modify |
| `tests/test_ref_panel_builder.py` | Unit tests for flip helpers + metadata swap + sign invariance | Modify |
| `tests/test_snp_identity.py` | Unit tests for `canonicalize_alleles` | Modify |
| `tests/test_column_inference.py` | Tests for the `FRQ`/`MAF` split | Modify |
| `tests/test_ref_panel.py` | Integration: build chr22 example, assert invariants | Modify |

---

## Task 1: Kernel orientation primitives (flip sign + oriented getter)

**Files:**
- Modify: `src/ldsc/_kernel/ref_panel_builder.py`
- Test: `tests/test_ref_panel_builder.py`

These are the pure, testable core of the refactor: the canonical rule (`freq < 0.5`)
and the genotype-column sign wrapper.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_ref_panel_builder.py` (add `import numpy as np` and the
kernel import at the top if not already present:
`from ldsc._kernel import ref_panel_builder as kb`):

```python
def test_orientation_flip_sign_rule():
    # freq is the .bim A2 allele frequency; A1 must end up minor (freq(A1) <= 0.5).
    freq = np.array([0.9, 0.5, 0.49999, 0.1, 0.5000001])
    sign = kb.orientation_flip_sign(freq)
    # Flip (sign -1) exactly where A2 freq < 0.5 (strict; 0.5 keeps PLINK order).
    np.testing.assert_array_equal(sign, np.array([1.0, 1.0, -1.0, -1.0, 1.0], dtype=np.float32))
    assert sign.dtype == np.float32


def test_oriented_snp_getter_negates_flagged_columns_sequentially():
    # Base getter returns deterministic standardized-like columns in two batches.
    rng = np.random.default_rng(0)
    full = rng.standard_normal((6, 5)).astype(np.float32)  # n=6 individuals, m=5 SNPs

    class _Base:
        def __init__(self):
            self.pos = 0

        def __call__(self, b):
            out = full[:, self.pos:self.pos + b].copy()
            self.pos += b
            return out

    sign = np.array([1.0, -1.0, 1.0, -1.0, 1.0], dtype=np.float32)
    getter = kb.make_oriented_snp_getter(_Base(), sign)
    first = getter(3)   # columns 0,1,2
    second = getter(2)  # columns 3,4
    got = np.hstack((first, second))
    expected = full * sign  # broadcast over rows
    np.testing.assert_allclose(got, expected, rtol=0, atol=0)


def test_oriented_getter_is_identity_when_no_flips():
    rng = np.random.default_rng(1)
    full = rng.standard_normal((4, 3)).astype(np.float32)

    class _Base:
        def __init__(self):
            self.pos = 0

        def __call__(self, b):
            out = full[:, self.pos:self.pos + b].copy()
            self.pos += b
            return out

    getter = kb.make_oriented_snp_getter(_Base(), np.ones(3, dtype=np.float32))
    np.testing.assert_allclose(getter(3), full, rtol=0, atol=0)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/test_ref_panel_builder.py -k "orientation or oriented_snp_getter or oriented_getter" -v`
Expected: FAIL with `AttributeError: module ... has no attribute 'orientation_flip_sign'`.

- [ ] **Step 3: Implement the primitives**

In `src/ldsc/_kernel/ref_panel_builder.py`, add these functions just above
`build_plink_metadata_frame` (around line 210):

```python
def orientation_flip_sign(freq_values: Sequence[float]) -> np.ndarray:
    """Per-SNP sign that re-orients genotypes so ``A1`` becomes the minor allele.

    ``freq_values`` is the frequency of the ``.bim`` second allele (the package's
    ``A2``). A SNP needs flipping when ``A2`` is the minor allele (``freq < 0.5``);
    at ``freq == 0.5`` PLINK order is kept (no flip), so the enforced invariant is
    ``freq(A1) <= 0.5``. Returns ``+1.0`` (keep) / ``-1.0`` (flip) as ``float32``
    so it broadcasts against ``float32`` standardized genotype columns without
    upcasting.
    """
    freq = np.asarray(freq_values, dtype=np.float64)
    return np.where(freq < 0.5, np.float32(-1.0), np.float32(1.0)).astype(np.float32)


def make_oriented_snp_getter(base_getter, flip_sign: np.ndarray):
    """Wrap a sequential standardized-SNP getter to apply ``flip_sign`` per column.

    The reference-panel windowed builder reads each retained SNP column exactly
    once, in increasing order, via successive ``base_getter(width)`` calls. This
    wrapper multiplies each returned batch by the matching slice of ``flip_sign``
    (``-1`` negates a standardized column, i.e. recodes it to the other allele),
    tracking the sequential position. Negation flips the correlation sign for
    pairs that include exactly one flipped SNP and leaves R2 (``sign**2``)
    unchanged.
    """
    flip_sign = np.asarray(flip_sign, dtype=np.float32)
    state = {"pos": 0}

    def getter(width):
        columns = base_getter(width)
        ncols = columns.shape[1]
        signs = flip_sign[state["pos"]:state["pos"] + ncols]
        state["pos"] += ncols
        if signs.size:
            columns = columns * signs  # (n, ncols) * (ncols,) broadcast over rows
        return columns

    return getter
```

Confirm `Sequence` is already imported at the top of the file (it is used by
`build_plink_metadata_frame`); if not, add `from typing import Sequence`.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/test_ref_panel_builder.py -k "orientation or oriented_snp_getter or oriented_getter" -v`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/ref_panel_builder.py tests/test_ref_panel_builder.py
git commit -m "feat(ref-panel): add A1=minor orientation primitives"
```

---

## Task 2: Swap A1/A2 in the metadata frame

**Files:**
- Modify: `src/ldsc/_kernel/ref_panel_builder.py:210-238` (`build_plink_metadata_frame`)
- Test: `tests/test_ref_panel_builder.py`

`build_plink_metadata_frame` currently copies `A1`/`A2` verbatim. Make it swap the
pair (using the same `orientation_flip_sign` rule) so the sidecar is canonical.
`MAF` already equals `min(freq, 1-freq)`, which after the swap is exactly
`freq(A1)`.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_ref_panel_builder.py`:

```python
def test_build_plink_metadata_frame_swaps_to_minor_a1():
    import pandas as pd

    class _Bim:
        df = pd.DataFrame(
            {
                "CHR": ["22", "22", "22"],
                "SNP": ["rs1", "rs2", "rs3"],
                "CM": [0.0, 0.0, 0.0],
                "BP": [100, 200, 300],
                "A1": ["G", "C", "A"],   # .bim allele 1
                "A2": ["A", "T", "G"],   # .bim allele 2 (package A2; freq is its frequency)
            }
        )

    bim = _Bim()
    kept_snps = [0, 1, 2]
    a2_freq = [0.8, 0.3, 0.5]            # SNP rs2 must flip; rs3 (==0.5) keeps order
    maf = np.minimum(a2_freq, 1 - np.asarray(a2_freq))

    out = kb.build_plink_metadata_frame(
        bim=bim, kept_snps=kept_snps, maf_values=maf, freq_values=a2_freq
    )

    # rs1: A2 major -> no swap; rs2: A2 minor -> swap; rs3: tie -> keep order.
    assert list(out["A1"]) == ["G", "T", "A"]
    assert list(out["A2"]) == ["A", "C", "G"]
    # MAF is unchanged numerically and is now freq(A1) (<= 0.5 for all rows).
    np.testing.assert_allclose(out["MAF"].to_numpy(), maf, rtol=0, atol=1e-12)
    assert (out["MAF"] <= 0.5).all()
    # Unordered allele set is preserved per SNP.
    assert {out["A1"][1], out["A2"][1]} == {"C", "T"}
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/test_ref_panel_builder.py::test_build_plink_metadata_frame_swaps_to_minor_a1 -v`
Expected: FAIL with `TypeError: build_plink_metadata_frame() got an unexpected keyword argument 'freq_values'`.

- [ ] **Step 3: Implement the swap**

Replace `build_plink_metadata_frame` body in
`src/ldsc/_kernel/ref_panel_builder.py` (signature gains `freq_values`):

```python
def build_plink_metadata_frame(
    *,
    bim,
    kept_snps: Sequence[int],
    maf_values: Sequence[float],
    freq_values: Sequence[float],
) -> pd.DataFrame:
    """Build the retained SNP metadata table after PLINK filtering.

    ``A1``/``A2`` are emitted in canonical orientation: ``A1`` is the minor allele
    (``freq(A1) <= 0.5``). ``freq_values`` is the ``.bim`` ``A2`` allele frequency;
    rows with ``freq < 0.5`` have ``A1``/``A2`` swapped. ``maf_values`` is the
    folded ``min(freq, 1-freq)`` and, after the swap, equals ``freq(A1)``.
    """

    kept_snps = np.asarray(kept_snps, dtype=int)
    maf_values = np.asarray(maf_values, dtype=float)
    freq_values = np.asarray(freq_values, dtype=float)
    if not (len(kept_snps) == len(maf_values) == len(freq_values)):
        raise LDSCInternalError(
            "build-ref-panel metadata assembly failed: kept_snps, maf_values, and "
            f"freq_values have different lengths ({len(kept_snps)} vs "
            f"{len(maf_values)} vs {len(freq_values)}). Most likely PLINK filtering "
            "bookkeeping desynchronized. Re-run with DEBUG logging and report the traceback."
        )
    kept = bim.df.iloc[kept_snps].reset_index(drop=True).copy()
    a1 = kept["A1"].astype(str).to_numpy()
    a2 = kept["A2"].astype(str).to_numpy()
    flip = orientation_flip_sign(freq_values) < 0
    canonical_a1 = np.where(flip, a2, a1)
    canonical_a2 = np.where(flip, a1, a2)
    out = pd.DataFrame(
        {
            "CHR": kept["CHR"].map(_normalize_map_chromosome),
            "SNP": kept["SNP"].astype(str),
            "CM": pd.to_numeric(kept["CM"], errors="coerce"),
            "POS": pd.to_numeric(kept["BP"], errors="raise").astype(np.int64),
            "A1": canonical_a1,
            "A2": canonical_a2,
            "MAF": maf_values.astype(float),
        }
    )
    return out.reset_index(drop=True)
```

- [ ] **Step 4: Run to verify it passes**

Run: `pytest tests/test_ref_panel_builder.py::test_build_plink_metadata_frame_swaps_to_minor_a1 -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/ref_panel_builder.py tests/test_ref_panel_builder.py
git commit -m "feat(ref-panel): canonicalize A1=minor in metadata frame"
```

---

## Task 3: Wire the flip into the build loop

**Files:**
- Modify: `src/ldsc/ref_panel_builder.py:928-932` (metadata call) and `:999-1007` (getter)
- Test: `tests/test_ref_panel.py` (integration, Task 8 covers the heavy assertions; here add a focused unit-ish test)

Pass `geno.freq` into the metadata frame and replace the bare `nextSNPs` lambda
with the oriented getter so `SIGN` is computed in canonical orientation.

- [ ] **Step 1: Update the metadata call**

In `src/ldsc/ref_panel_builder.py`, change the `build_plink_metadata_frame` call
(currently lines 928-932):

```python
            metadata = kernel_builder.build_plink_metadata_frame(
                bim=bim,
                kept_snps=geno.kept_snps,
                maf_values=geno.maf,
                freq_values=geno.freq,
            )
```

- [ ] **Step 2: Compute the flip sign once and wrap the getter**

Still in `src/ldsc/ref_panel_builder.py`, immediately before the
`kernel_builder.write_r2_parquet(` call (line 999), insert:

```python
            flip_sign = kernel_builder.orientation_flip_sign(geno.freq)
            oriented_snp_getter = kernel_builder.make_oriented_snp_getter(
                lambda b: geno.nextSNPs(b, dtype=np.float32),
                flip_sign,
            )
```

Then replace the `standardized_snp_getter=` argument inside
`yield_pairwise_r2_rows(...)` (line 1003) with:

```python
                    standardized_snp_getter=oriented_snp_getter,
```

- [ ] **Step 3: Run the existing ref-panel suite to verify nothing breaks structurally**

Run: `pytest tests/test_ref_panel_builder.py tests/test_ref_panel.py -v`
Expected: PASS (existing tests still pass; new orientation tests pass). If an
existing test asserted verbatim PLINK `A1`/`A2` order on a flipped SNP, update its
expectation to the canonical order and note it in the commit body.

- [ ] **Step 4: Commit**

```bash
git add src/ldsc/ref_panel_builder.py
git commit -m "feat(ref-panel): apply A1=minor flip in build loop (meta + SIGN)"
```

---

## Task 4: Sign-invariance + canonical-SIGN integration check

**Files:**
- Test: `tests/test_ref_panel.py`

Prove on the real example panel that (a) R² is unchanged by orientation and
(b) the meta sidecar is canonical. R²-invariance is the safety property for
h2/rg.

- [ ] **Step 1: Write the pure orientation-invariance test**

Append to `tests/test_ref_panel_builder.py` (it already imports `numpy as np`):

```python
def test_r2_is_orientation_invariant():
    # Unit-level proof that negating standardized columns leaves R2 unchanged
    # and flips SIGN iff exactly one endpoint of a pair was flipped.
    rng = np.random.default_rng(7)
    n, m = 50, 8
    base = rng.standard_normal((n, m)).astype(np.float32)
    base = (base - base.mean(0)) / base.std(0)
    sign = np.array([1, -1, 1, 1, -1, -1, 1, -1], dtype=np.float32)
    corr_plain = base.T @ (base / n)
    flipped = base * sign
    corr_flipped = flipped.T @ (flipped / n)
    np.testing.assert_allclose(corr_plain ** 2, corr_flipped ** 2, rtol=0, atol=1e-6)
    outer = np.outer(sign, sign)
    np.testing.assert_allclose(np.sign(corr_flipped), np.sign(corr_plain) * outer, atol=1e-6)
```

- [ ] **Step 2: Write the end-to-end canonical-sidecar test**

Append to `tests/test_ref_panel.py`. Mirror the existing build tests in that file
— they construct `RefPanelConfig(backend="plink", plink_prefix=...)` and drive the
real builder; the chr22 PLINK fixture and skip guard already exist as
`_CHR22 = Path(__file__).resolve().parent / "fixtures" / "minimal_external_resources" / "plink" / "hm3_chr22_subset"`.
Reuse the exact build-invocation idiom from the nearest passing build test in the
file (e.g. the one near line 715 that writes a panel to `tmp_path`); do not invent
a new builder API. After the panel is written to `tmp_path`, assert:

```python
import glob
import pandas as pd

def _assert_meta_canonical(out_dir):
    meta_paths = glob.glob(str(out_dir / "**" / "chr*_meta.tsv.gz"), recursive=True)
    assert meta_paths, "no meta sidecar written"
    meta = pd.read_csv(meta_paths[0], sep="\t")
    assert (meta["MAF"] <= 0.5 + 1e-12).all()                      # MAF = freq(A1) <= 0.5
    assert (meta["A1"].astype(str) != meta["A2"].astype(str)).all()  # valid biallelic pair
```

Call `_assert_meta_canonical(tmp_path / <panel-subdir>)` from a test that builds
the `_CHR22` fixture panel, guarded by the same `pytest.skip(...)` the
neighboring fixture-dependent tests use when `_CHR22` is unavailable.

- [ ] **Step 3: Run**

Run: `pytest tests/test_ref_panel_builder.py -k orientation_invariant -v` and
`pytest tests/test_ref_panel.py -k canonical -v`
Expected: PASS. (`test_r2_is_orientation_invariant` is pure and must pass
immediately; the `_CHR22` build test validates end-to-end canonicalization, or is
skipped if the fixture is unavailable.)

- [ ] **Step 4: Commit**

```bash
git add tests/test_ref_panel_builder.py tests/test_ref_panel.py
git commit -m "test(ref-panel): assert A1=minor sidecar and R2 orientation-invariance"
```

---

## Task 5: Defensive-read helper `canonicalize_alleles`

**Files:**
- Modify: `src/ldsc/_kernel/snp_identity.py`
- Test: `tests/test_snp_identity.py`

Normalizes any externally-supplied table that carries `A1`/`A2` plus an *oriented*
`A1` frequency to canonical `A1`=minor. Idempotent; ties (`== 0.5`) kept.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_snp_identity.py` (ensure `import numpy as np`,
`import pandas as pd`, and `from ldsc._kernel import snp_identity as si`):

```python
def test_canonicalize_alleles_flips_oriented_rows():
    df = pd.DataFrame(
        {
            "A1": ["G", "C", "A"],
            "A2": ["A", "T", "G"],
            "FRQ": [0.20, 0.80, 0.50],  # oriented freq of A1
        }
    )
    out = si.canonicalize_alleles(df, freq_col="FRQ")
    # Row 0: freq(A1)=0.2 <=0.5 keep. Row 1: 0.8 >0.5 flip. Row 2: 0.5 tie keep.
    assert list(out["A1"]) == ["G", "T", "A"]
    assert list(out["A2"]) == ["A", "C", "G"]
    np.testing.assert_allclose(out["FRQ"].to_numpy(), [0.20, 0.20, 0.50])
    assert (out["FRQ"] <= 0.5 + 1e-12).all()


def test_canonicalize_alleles_is_idempotent():
    df = pd.DataFrame({"A1": ["G", "C"], "A2": ["A", "T"], "FRQ": [0.7, 0.3]})
    once = si.canonicalize_alleles(df, freq_col="FRQ")
    twice = si.canonicalize_alleles(once, freq_col="FRQ")
    pd.testing.assert_frame_equal(once.reset_index(drop=True), twice.reset_index(drop=True))


def test_canonicalize_alleles_does_not_mutate_input():
    df = pd.DataFrame({"A1": ["C"], "A2": ["T"], "FRQ": [0.9]})
    _ = si.canonicalize_alleles(df, freq_col="FRQ")
    assert list(df["A1"]) == ["C"] and list(df["FRQ"]) == [0.9]
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_snp_identity.py -k canonicalize_alleles -v`
Expected: FAIL with `AttributeError: module ... has no attribute 'canonicalize_alleles'`.

- [ ] **Step 3: Implement the helper**

Add to `src/ldsc/_kernel/snp_identity.py` (it already imports `pandas as pd`; add
`import numpy as np` if absent):

```python
def canonicalize_alleles(
    df: pd.DataFrame,
    *,
    a1_col: str = "A1",
    a2_col: str = "A2",
    freq_col: str = "MAF",
) -> pd.DataFrame:
    """Re-orient ``A1``/``A2`` so ``A1`` is the minor allele (``freq(A1) <= 0.5``).

    ``freq_col`` must hold the *oriented* frequency of the allele in ``a1_col``
    (may exceed 0.5). Rows with ``freq(A1) > 0.5`` have their alleles swapped and
    ``freq`` replaced by ``1 - freq``; ``freq == 0.5`` rows are left as-is
    (tie-break: keep order). The function is idempotent and returns a copy without
    mutating the input. It cannot recover orientation from a folded value, so the
    caller must pass an oriented column (see design §2.1).
    """
    out = df.copy()
    freq = pd.to_numeric(out[freq_col], errors="coerce").to_numpy(dtype=float)
    flip = freq > 0.5
    a1 = out[a1_col].to_numpy()
    a2 = out[a2_col].to_numpy()
    out[a1_col] = np.where(flip, a2, a1)
    out[a2_col] = np.where(flip, a1, a2)
    out[freq_col] = np.where(flip, 1.0 - freq, freq)
    return out
```

- [ ] **Step 4: Run to verify pass**

Run: `pytest tests/test_snp_identity.py -k canonicalize_alleles -v`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/snp_identity.py tests/test_snp_identity.py
git commit -m "feat(identity): add canonicalize_alleles defensive-read helper"
```

---

## Task 6: Split oriented `FRQ`-family from folded `MAF` in the registry

**Files:**
- Modify: `src/ldsc/column_inference.py:112-118,177-181`
- Test: `tests/test_column_inference.py`

Non-sumstats sidecars must be able to declare orientation by column name. Make
`MAF` mean folded-only, add an oriented `FRQ`-family spec for sidecar reads.
Sumstats specs (`RAW_SUMSTATS_*`, `INTERNAL_SUMSTATS_*`) are unchanged.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_column_inference.py` (mirror the existing resolution-helper
usage in that file; the assertions below use the public alias families directly):

```python
from ldsc import column_inference as ci


def test_maf_alias_family_is_folded_only():
    # MAF (folded) no longer absorbs FRQ/FREQ; those are the oriented family.
    assert ci.MAF_COLUMN_ALIASES == ("MAF",)


def test_oriented_frq_alias_family_exists():
    fam = set(ci.ORIENTED_FRQ_COLUMN_ALIASES)
    assert {"FRQ", "FREQ", "FREQUENCY", "EAF", "A1_FRQ", "FRQ_A1", "FREQ_A1"} <= fam
    assert "MAF" not in fam  # MAF stays folded, never oriented
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_column_inference.py -k "folded_only or oriented_frq" -v`
Expected: FAIL (`MAF_COLUMN_ALIASES` still contains FRQ/FREQ; `ORIENTED_FRQ_COLUMN_ALIASES` undefined).

- [ ] **Step 3: Implement the split**

In `src/ldsc/column_inference.py`, change line 112 and add the oriented family +
spec right after `MAF_COLUMN_SPEC` (line 118):

```python
MAF_COLUMN_ALIASES = ("MAF",)
```

```python
MAF_COLUMN_SPEC = ColumnSpec("MAF", MAF_COLUMN_ALIASES, "folded minor-allele frequency")
# Oriented frequency of the A1 allele (may exceed 0.5). Distinct from folded MAF
# so external sidecars can signal orientation by column name (design §2.1).
ORIENTED_FRQ_COLUMN_ALIASES = ("FRQ", "FREQ", "FREQUENCY", "EAF", "A1_FRQ", "FRQ_A1", "FREQ_A1")
ORIENTED_FRQ_COLUMN_SPEC = ColumnSpec("FRQ", ORIENTED_FRQ_COLUMN_ALIASES, "oriented A1-allele frequency")
```

Leave `ANNOTATION_METADATA_SPECS` / `REFERENCE_METADATA_SPECS` keyed on
`MAF_COLUMN_SPEC` (those internal artifacts use folded `MAF`). The oriented spec
is for external sidecar readers that opt into normalization; export it for that
use. Do **not** touch `RAW_SUMSTATS_REQUIRED_OR_OPTIONAL_SPECS` (line 158) or the
`INTERNAL_SUMSTATS_*` specs — sumstats keep the existing conflation.

- [ ] **Step 4: Run to verify pass + no regressions**

Run: `pytest tests/test_column_inference.py -v`
Expected: PASS. If a pre-existing test asserted that `FRQ`/`FREQ` resolve to
`MAF` in an annotation/reference context, update it to reflect the split (folded
`MAF` only) and record it in the commit body.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/column_inference.py tests/test_column_inference.py
git commit -m "feat(columns): split oriented FRQ family from folded MAF"
```

---

## Task 7: Assert the `MAF <= 0.5` invariant on LD-score output

**Files:**
- Modify: `src/ldsc/ldscore_calculator.py` (helper + call from `_aggregate_chromosome_results`, `:590`)
- Test: `tests/test_ldscore_workflow.py`

LD-score artifacts inherit `MAF`/`A1`/`A2` from the canonical ref panel. Add a
cheap invariant assertion so a non-canonical panel fails loudly instead of
silently mis-tying `MAF` to `A1`.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_ldscore_workflow.py` (add imports if absent):

```python
import pandas as pd
import pytest
from ldsc import ldscore_calculator as lc
from ldsc.errors import LDSCInternalError


def test_ldscore_rejects_non_minor_maf():
    bad = pd.DataFrame({"CHR": ["22"], "SNP": ["rs1"], "POS": [1],
                        "A1": ["A"], "A2": ["G"], "MAF": [0.73]})
    with pytest.raises(LDSCInternalError, match="MAF"):
        lc._assert_canonical_maf(bad)


def test_assert_canonical_maf_passes_for_minor():
    ok = pd.DataFrame({"MAF": [0.0, 0.2, 0.5]})
    lc._assert_canonical_maf(ok)  # no raise
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_ldscore_workflow.py -k canonical_maf -v`
Expected: FAIL with `AttributeError: module ... has no attribute '_assert_canonical_maf'`.

- [ ] **Step 3: Implement the assertion**

Add the module-level helper in `src/ldsc/ldscore_calculator.py` (module already
imports `pandas as pd`; ensure `from .errors import LDSCInternalError` is present)
and call `_assert_canonical_maf(<aggregated metadata frame>)` from
`_aggregate_chromosome_results` (`:590`) immediately before it returns the
aggregated result, passing the DataFrame that carries the `MAF` column:

```python
def _assert_canonical_maf(metadata: pd.DataFrame) -> None:
    """Verify the canonical invariant that MAF = freq(A1) is the minor allele."""
    if "MAF" not in metadata.columns:
        return
    maf = pd.to_numeric(metadata["MAF"], errors="coerce")
    over = maf[maf > 0.5 + 1e-9]
    if len(over):
        raise LDSCInternalError(
            "LD-score metadata carries MAF > 0.5, violating the canonical "
            f"A1=minor invariant ({len(over)} rows; first={float(over.iloc[0]):.4f}). "
            "Most likely the reference panel was built before allele-orientation "
            "canonicalization. Rebuild the reference panel with the current package."
        )
```

Ensure `LDSCInternalError` is imported in the module (it is used elsewhere in the
package; add the import from `.errors` if missing).

- [ ] **Step 4: Run to verify pass**

Run: `pytest tests/test_ldscore_workflow.py -k canonical_maf -v`
Expected: PASS. Then confirm canonical panels still pass end-to-end:
`pytest tests/test_ldscore_workflow.py -v`.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/ldscore_calculator.py tests/test_ldscore_workflow.py
git commit -m "feat(ldscore): assert canonical A1=minor MAF invariant"
```

---

## Task 8: Document signed r and the sumstats carve-out (no behavior change)

**Files:**
- Modify: `src/ldsc/r2_query.py` (docstrings + argparse help around `:323,376,489`)
- Modify: `src/ldsc/_kernel/sumstats_munger.py` (one carve-out comment)

- [ ] **Step 1: Update `query-r2` signed-r documentation**

In `src/ldsc/r2_query.py`, update the `_signed_r` docstring (line 379) and the
`--with-r` argparse help (line 489) to state the orientation:

`_signed_r` docstring →
```python
        """Compute the signed Pearson r from adjusted r2 and the sign column.

        Sign is defined in the canonical reference-panel orientation: r is the
        correlation between the two SNPs' A1 (minor) allele dosages. A positive
        value means the minor alleles co-occur on haplotypes more than chance
        (positive LD between minor alleles). A1 is guaranteed minor by the
        canonical panel (design 2026-06-07 allele-orientation-canonicalization).
        """
```

`--with-r` help →
```python
    parser.add_argument(
        "--with-r",
        action="store_true",
        help="Add a signed Pearson r column (correlation of A1/minor allele dosages; "
        "positive = positive LD between minor alleles). Allele-aware panels only.",
    )
```

- [ ] **Step 2: Add the sumstats carve-out comment**

In `src/ldsc/_kernel/sumstats_munger.py`, just above the `frq = np.minimum(frq, 1 - frq)`
line (line 252), add:

```python
    # Carve-out: sumstats A1 is the EFFECT allele, not the minor allele. The fold
    # below is a local mask for the --maf-min threshold only; never reorient
    # A1/A2 or overwrite FRQ here (FRQ stays freq(A1), may exceed 0.5).
```

- [ ] **Step 3: Verify nothing broke**

Run: `pytest tests/test_sumstats_munger.py tests/test_r2_query.py -q`
(documentation-only change; all tests should still pass). Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add src/ldsc/r2_query.py src/ldsc/_kernel/sumstats_munger.py
git commit -m "docs(code): document A1=minor signed-r and sumstats carve-out"
```

---

## Task 9: User-facing documentation

**Files:**
- Modify: `docs/current/column-schema.md`
- Modify: `docs/current/data-flow.md`

- [ ] **Step 1: Update `column-schema.md`**

Add a section "Allele orientation and frequency" containing: the unifying rule
(design §1), the per-artifact `A1`/`A2`/`FRQ`/`MAF` table (design §2), the
`MAF <= 0.5` guarantee and its sumstats exception (design §3), the sign rule
(design §4: `A1`=minor; signed r = correlation of `A1`/minor allele dosages;
positive = positive LD between minor alleles), and the oriented-vs-folded
external-sidecar naming convention (design §2.1: `FRQ`-family oriented vs `MAF`
folded). Copy the tables from the spec verbatim and link to the spec path.

- [ ] **Step 2: Update `data-flow.md`**

In the `build-ref-panel` flow description, add a note: the A1=minor flip is
applied at standardization in the builder (negate standardized columns where
`.bim` A2 freq < 0.5; swap A1/A2 in the meta sidecar) so `SIGN` and `MAF` are
canonical; `nextSNPs` stays orientation-free and the r²-based in-PLINK LD-score
path is unaffected (sign-invariant).

- [ ] **Step 3: Verify docs cross-references resolve**

Run: `grep -rn "allele-orientation-canonicalization" docs/current/`
Expected: the new references appear and point to the spec filename.

- [ ] **Step 4: Commit**

```bash
git add docs/current/column-schema.md docs/current/data-flow.md
git commit -m "docs: document package-wide A1/A2/FRQ/MAF orientation rules"
```

---

## Task 10: Full-suite verification

**Files:** none (verification only)

- [ ] **Step 1: Run the whole test suite**

Run: `pytest -q`
Expected: all tests pass. Investigate and fix any failure before proceeding;
record any pre-existing test whose expectation legitimately changed (e.g. a
flipped-SNP `A1`/`A2` order) in that fix's commit body.

- [ ] **Step 2: Build the example panel end-to-end and spot-check the invariant**

Run:
```bash
ldsc build-ref-panel \
  --plink-prefix /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/resources/example_1kg_30x/genomes_30x_chr22 \
  --output-dir /tmp/ldsc_canon_check --ld-wind-snps 100
python - <<'PY'
import pandas as pd, glob
f = glob.glob('/tmp/ldsc_canon_check/**/chr22_meta.tsv.gz', recursive=True)[0]
m = pd.read_csv(f, sep='\t')
print('rows', len(m), 'MAF<=0.5 all?', bool((m['MAF']<=0.5+1e-12).all()),
      'A1!=A2 all?', bool((m['A1'].astype(str)!=m['A2'].astype(str)).all()))
PY
```
Expected: `MAF<=0.5 all? True` and `A1!=A2 all? True`. (Use the build's actual
emitted-build subdir; adjust `--ld-wind-snps`/flags to match the example's
required options if the CLI reports missing inputs.)

- [ ] **Step 3: Final commit (if any verification fix was needed)**

```bash
git add -A && git commit -m "test: verify package-wide allele canonicalization"
```

---

## Self-review notes (for the implementer)

- **R²/h2/rg safety:** the flip only negates standardized columns; R² = `corr**2`
  is provably unchanged (Task 4 `test_r2_is_orientation_invariant`). Do not add
  any orientation logic to `nextSNPs` or to the r²-based LD-score path.
- **Ordering alignment:** `geno.freq`, `geno.maf`, `geno.kept_snps`, the metadata
  rows, and the sequential `nextSNPs` output are all in the same kept-SNP order;
  `flip_sign` indexes that order. This must hold per emitted build (hg19/hg38),
  which is why `flip_sign` is recomputed from `geno.freq` inside the per-build
  loop (Task 3), not once globally.
- **Streaming and in-RAM modes** share `nextSNPs`/`standardized_snp_getter`, so
  the wrapper covers both with no extra work.
- **Legacy artifacts** are out of scope — they are rebuilt; the LD-score
  assertion (Task 7) is the guardrail that catches a stale panel.
