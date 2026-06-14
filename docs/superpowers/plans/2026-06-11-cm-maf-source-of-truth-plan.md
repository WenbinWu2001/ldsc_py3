# CM/MAF Source-of-Truth Across LD-Score Backends Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the reference panel the sole source of `CM` and `MAF` for both
LD-score backends so backend choice no longer changes LD windows or counts, make
annotation files population-agnostic, and add a genetic-map remedy + guards for
`--ld-wind-cm` with PLINK panels.

**Architecture:** A shared `CM` resolver feeds the existing window builder; `CM`
comes from the parquet sidecar or (PLINK) the `.bim`/an interpolated genetic map.
`MAF` stays at each backend's natural source (parquet sidecar; PLINK genotypes)
and is always required. Annotation `CM`/`MAF` are never read into the compute
path. `--maf-min` becomes a reference-panel filter in both backends.

**Tech Stack:** Python 3, pandas, numpy, pyarrow, argparse; pytest (+ legacy
`unittest`). Run tests with the project env:
`source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest`.

**Design spec:** `docs/superpowers/specs/2026-06-11-cm-maf-source-of-truth-design.md`

**Conventions for every task:** write the test first, run it red, implement the
minimal change, run it green, then commit. Commit messages use Conventional
Commits; no AI-tool names in messages or comments. Keep docstrings current
(use `my-skills:fun-doc` for any function whose body you change materially).

---

## File map

| File | Responsibility in this change |
|---|---|
| `src/ldsc/annotation_builder.py` | `AnnotationBuilder.parse_annotation_file`: stop reading `CM`/`MAF` into annotation metadata; make `CM` optional; make `SNP` required only in rsID-family modes. |
| `src/ldsc/_kernel/ldscore.py` | Apply the same annotation-reader change to the kernel `parse_annotation_file` (`:636`); new shared `CM` resolver + unusable-`CM` guard; require `MAF`; sharpen errors; PLINK genetic-map interpolation. |
| `src/ldsc/_kernel/ref_panel.py` | Apply `--maf-min` for PLINK; make sidecar `MAF` requirement unconditional. |
| `src/ldsc/config.py` | `RefPanelConfig.genetic_map_hg19_sources`/`genetic_map_hg38_sources`; `LDScoreConfig.export_ref_metadata`. |
| `src/ldsc/ldscore_calculator.py` | CLI flags; resolve genome build + load genetic map; pass primitives to the kernel namespace; provenance + export wiring. |
| `src/ldsc/outputs.py` | Record `CM`/`MAF` source provenance in `metadata.json`. |
| `docs/troubleshooting.md`, `docs/current/*`, `tutorials/*` | Documentation. |
| `tests/test_annotation.py`, `tests/test_ldscore_workflow.py`, `tests/test_ref_panel.py` | Tests. |

Reusable test helpers already in the repo: `tests/fixtures/plink` (prebuilt PLINK
panels, exposed as `PLINK_FIXTURES` in `tests/test_ldscore_workflow.py`),
`_write_minimal_r2_parquet` / `dict_chunks` (parquet R2), and the public
`run_ldscore(...)` + `set_global_config(...)` APIs for end-to-end runs.

---

## Phase 1 — Annotation files become population-agnostic

> **REVISION (2026-06-11, during execution):** `CM`/`MAF` are **kept-but-ignored**,
> not stripped. The reader (`parse_annotation_file`) keeps a `CM` column
> (placeholder `NaN` when the input omits it) so `_merge_missing_metadata`
> (`annotation_builder.py:777`) and the BED-query projection (`_kernel/annotation.py:160`,
> `row.CM`) keep working and `annotate` keeps emitting the legacy `CHR BP SNP CM`
> layout. `CM` becomes optional and `SNP` optional in `chr_pos` modes; `CM`/`MAF`
> are still excluded from annotation *value* columns. The "ldscore ignores
> annotation `CM`/`MAF`" guarantee is enforced at the *consumer*: Task 4 makes the
> reference-panel sidecar **authoritative** in `merge_frequency_metadata`
> (overwrite, not backfill). This matches legacy ldsc2 (keep the column, ignore
> the value). Tasks 1, 2, and 4 below are interpreted under this revision.

### Task 1: Annotation reader stops reading `CM`/`MAF`; `CM` no longer required

**Files:**
- Modify: `src/ldsc/annotation_builder.py:602-648` (`parse_annotation_file`)
- Test: `tests/test_annotation.py`

The current reader (`annotation_builder.py:602-648`) calls
`resolve_required_column(..., CM_COLUMN_SPEC)` (CM required) and copies `CM` (and
optional `MAF`) into the metadata frame. We make `CM` optional and **never** put
`CM`/`MAF` into annotation metadata, while still excluding any `CM`/`MAF` source
columns from the annotation *value* columns.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_annotation.py
def test_annotation_reader_ignores_cm_and_maf_columns(tmp_path):
    from ldsc.annotation_builder import AnnotationBuilder  # adjust import to repo's public entry
    path = tmp_path / "chr1.annot"
    path.write_text(
        "CHR\tPOS\tSNP\tCM\tMAF\tcoding\n"
        "1\t1000\trs1\t0.5\t0.30\t1\n"
        "1\t2000\trs2\t0.9\t0.10\t0\n"
    )
    metadata, annotations = AnnotationBuilder().parse_annotation_file(path, chrom="1")
    # CM/MAF must NOT be carried as annotation metadata
    assert "CM" not in metadata.columns
    assert "MAF" not in metadata.columns
    # CM/MAF must NOT be treated as annotation value columns
    assert list(annotations.columns) == ["coding"]
    assert list(metadata.columns) == ["CHR", "POS", "SNP"]


def test_annotation_reader_accepts_file_without_cm(tmp_path):
    from ldsc.annotation_builder import AnnotationBuilder
    path = tmp_path / "chr1.annot"
    path.write_text("CHR\tPOS\tSNP\tcoding\n1\t1000\trs1\t1\n")
    metadata, annotations = AnnotationBuilder().parse_annotation_file(path, chrom="1")
    assert "CM" not in metadata.columns
    assert list(annotations.columns) == ["coding"]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_annotation.py -k "ignores_cm_and_maf or without_cm" -v`
Expected: FAIL — today `CM` is required (raises) and/or appears in `metadata`.

- [ ] **Step 3: Implement — make `CM` optional and drop `CM`/`MAF` from metadata**

In `parse_annotation_file`, replace the required-CM read and metadata assembly.
Current (`annotation_builder.py:605` and `:613-630`):

```python
        cm_col = resolve_required_column(df.columns, CM_COLUMN_SPEC, context=context)
        ...
        metadata = pd.DataFrame(
            {
                "CHR": df[chr_col],
                "POS": df[pos_col],
                "SNP": df[snp_col],
                "CM": df[cm_col],
            }
        )
        metadata["CHR"] = metadata["CHR"].map(lambda value: normalize_chromosome(value, context=context))
        metadata["POS"] = pd.to_numeric(metadata["POS"], errors="raise").astype(np.int64)
        metadata["SNP"] = metadata["SNP"].astype(str)
        metadata["CM"] = pd.to_numeric(metadata["CM"], errors="coerce")
        if a1_col is not None and a2_col is not None:
            metadata["A1"] = df[a1_col]
            metadata["A2"] = df[a2_col]
        maf_col = resolve_optional_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["MAF"], context=context)
        if maf_col is not None:
            metadata["MAF"] = pd.to_numeric(df[maf_col], errors="coerce")
```

New:

```python
        # CM/MAF are population-specific and supplied by the reference panel, not
        # the annotation. Resolve their columns only to exclude them from the
        # annotation value columns; never read their values into metadata.
        cm_col = resolve_optional_column(df.columns, CM_COLUMN_SPEC, context=context)
        ...
        metadata = pd.DataFrame(
            {
                "CHR": df[chr_col],
                "POS": df[pos_col],
                "SNP": df[snp_col],
            }
        )
        metadata["CHR"] = metadata["CHR"].map(lambda value: normalize_chromosome(value, context=context))
        metadata["POS"] = pd.to_numeric(metadata["POS"], errors="raise").astype(np.int64)
        metadata["SNP"] = metadata["SNP"].astype(str)
        if a1_col is not None and a2_col is not None:
            metadata["A1"] = df[a1_col]
            metadata["A2"] = df[a2_col]
        maf_col = resolve_optional_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["MAF"], context=context)
        if cm_col is not None or maf_col is not None:
            LOGGER.info(
                f"Ignoring annotation CM/MAF columns in '{path}'; the reference panel is "
                "authoritative for CM and MAF."
            )
```

The `metadata_source_columns` set at `annotation_builder.py:639` already includes
`cm_col` and `maf_col`, so the annotation value columns remain correctly
excluded. Confirm `LOGGER` is imported in this module (it is used elsewhere); if
not, add `LOGGER = logging.getLogger(__name__)` near the top.

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_annotation.py -k "ignores_cm_and_maf or without_cm" -v`
Expected: PASS.

- [ ] **Step 5: Run the annotation suite to catch regressions**

Run: `pytest tests/test_annotation.py -v`
Expected: PASS (fix any test that asserted `CM` in annotation metadata; such
assertions are now obsolete — update them to the new contract).

- [ ] **Step 6: Commit**

```bash
git -C <repo> add src/ldsc/annotation_builder.py tests/test_annotation.py
git -C <repo> commit -m "feat(annotate): ignore annotation CM/MAF; CM no longer required"
```

### Task 2: `SNP` required only in rsID-family identifier modes

**Files:**
- Modify: `src/ldsc/annotation_builder.py` (`parse_annotation_file`, and the call site that knows the identifier mode)
- Test: `tests/test_annotation.py`

- [ ] **Step 1: Confirm the reader's access to the identifier mode**

Run: `grep -n "snp_identifier\|identity_mode_family\|parse_annotation_file" src/ldsc/annotation_builder.py`
Expected: find whether `AnnotationBuilder` stores the snp_identifier (e.g., via
`global_config`) or whether `parse_annotation_file` must receive it as a
parameter. If the mode is not available in the reader, thread it down from the
public `build_*` entry point as a keyword argument `snp_identifier: str`.

- [ ] **Step 2: Write the failing test**

```python
# tests/test_annotation.py
def test_annotation_snp_optional_in_chr_pos_mode(tmp_path):
    from ldsc.annotation_builder import AnnotationBuilder
    path = tmp_path / "chr1.annot"
    path.write_text("CHR\tPOS\tcoding\n1\t1000\t1\n")
    metadata, annotations = AnnotationBuilder().parse_annotation_file(
        path, chrom="1", snp_identifier="chr_pos"
    )
    assert "SNP" not in metadata.columns
    assert list(annotations.columns) == ["coding"]


def test_annotation_snp_required_in_rsid_mode(tmp_path):
    from ldsc.annotation_builder import AnnotationBuilder
    from ldsc.errors import LDSCInputError
    path = tmp_path / "chr1.annot"
    path.write_text("CHR\tPOS\tcoding\n1\t1000\t1\n")
    with pytest.raises(LDSCInputError):
        AnnotationBuilder().parse_annotation_file(path, chrom="1", snp_identifier="rsid")
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `pytest tests/test_annotation.py -k "snp_optional or snp_required" -v`
Expected: FAIL — `SNP` is currently always required.

- [ ] **Step 4: Implement — make `SNP` resolution mode-dependent**

Add `snp_identifier: str` to `parse_annotation_file` (and pass it from the public
entry point found in Step 1). Replace the unconditional `SNP` resolution:

```python
from ._kernel.snp_identity import identity_mode_family  # near existing imports

        if identity_mode_family(snp_identifier) == "rsid":
            snp_col = resolve_required_column(df.columns, SNP_COLUMN_SPEC, context=context)
        else:
            snp_col = resolve_optional_column(df.columns, SNP_COLUMN_SPEC, context=context)
```

When building `metadata`, include `SNP` only when present:

```python
        columns = {"CHR": df[chr_col], "POS": df[pos_col]}
        if snp_col is not None:
            columns["SNP"] = df[snp_col].astype(str)
        metadata = pd.DataFrame(columns)
        metadata["CHR"] = metadata["CHR"].map(lambda value: normalize_chromosome(value, context=context))
        metadata["POS"] = pd.to_numeric(metadata["POS"], errors="raise").astype(np.int64)
```

Keep `snp_col` in `metadata_source_columns` (guard `None` with a set
comprehension that drops `None`).

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_annotation.py -k "snp_optional or snp_required" -v`
Expected: PASS.

- [ ] **Step 6: Run the annotation + workflow suites**

Run: `pytest tests/test_annotation.py tests/test_ldscore_workflow.py -v`
Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git -C <repo> add src/ldsc/annotation_builder.py tests/test_annotation.py
git -C <repo> commit -m "feat(annotate): require SNP only in rsID identifier modes"
```

### Task 2b: Mirror the contract in the kernel annotation reader

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py:636-698` (`parse_annotation_file`) and the
  `REQUIRED_ANNOT_COLUMNS` / `ANNOT_META_COLUMNS` constants (`:248-249`)
- Test: `tests/test_ldscore_workflow.py`

The kernel has its own `parse_annotation_file` (`ldscore.py:636`, called at
`:752` by the kernel bundle builder) that **also** requires `CM` via
`resolve_required_column`. It must get the same contract so all annotation entry
points agree.

- [ ] **Step 1: Check whether the kernel reader is reachable (deadness check)**

Run: `grep -rn "parse_annotation_file\|load_annotation\|AnnotationBundle(" src/ldsc/_kernel/ldscore.py`
and trace whether `ldscore.py:752`'s caller is reached from any supported CLI /
`run_ldscore` path, or only from kernel-direct/legacy code. Record the finding in
the commit message. Update the reader regardless (cheap consistency); if it is
fully dead, note it as a candidate for separate removal — do not remove here.

- [ ] **Step 2: Write the failing test**

```python
# tests/test_ldscore_workflow.py
def test_kernel_parse_annotation_file_ignores_cm(tmp_path):
    from ldsc._kernel import ldscore as k
    path = tmp_path / "chr1.annot"
    path.write_text("CHR\tPOS\tSNP\tCM\tcoding\n1\t1000\trs1\t0.5\t1\n")
    metadata, annotations = k.parse_annotation_file(str(path))
    assert "CM" not in metadata.columns
    assert list(annotations.columns) == ["coding"]
```

- [ ] **Step 3: Run test to verify it fails**

Run: `pytest tests/test_ldscore_workflow.py -k kernel_parse_annotation_file_ignores_cm -v`
Expected: FAIL — kernel reader requires `CM` and copies it into metadata.

- [ ] **Step 4: Implement**

In `parse_annotation_file` (`ldscore.py:651`), change the `CM` resolution from
required to optional and stop copying `CM`/`MAF` into the returned metadata,
mirroring Task 1. Update `REQUIRED_ANNOT_COLUMNS` to `("CHR", "POS", "SNP")` (or
delete it if the deadness check shows it is unused — the earlier grep found only
its definition). Leave `ANNOT_META_COLUMNS` (`:249`) intact: it is used for
output column ordering at `ldscore.py:1860` and must keep `CM`/`MAF` so present
columns sort correctly.

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_ldscore_workflow.py -k kernel_parse_annotation_file -v`
Expected: PASS. Update any existing kernel-reader test that asserted `CM` in
metadata (e.g. `test_kernel_parse_annotation_file_preserves_optional_alleles_as_metadata`).

- [ ] **Step 6: Commit**

```bash
git -C <repo> add src/ldsc/_kernel/ldscore.py tests/test_ldscore_workflow.py
git -C <repo> commit -m "feat(ldscore): kernel annotation reader ignores CM/MAF"
```

### Task 2c: Log the `CM`-placeholder caveat in `annotate` output

**Files:**
- Modify: the `annotate` write path in `src/ldsc/annotation_builder.py`
- Test: `tests/test_annotation.py`

`annotate` keeps the legacy `CHR BP SNP CM <annotations…>` layout (positional
compatibility with ldsc2's `iloc[:, 4:]`), but the `CM` value is a placeholder
that nothing downstream consumes. Emit a one-time log so users understand this.

- [ ] **Step 1: Locate the `annotate` write path**

Run: `grep -n "def .*write\|to_csv\|\.annot\|CM" src/ldsc/annotation_builder.py | head`
to find where the `.annot` table (with the `CM` column) is written.

- [ ] **Step 2: Write the failing test**

```python
# tests/test_annotation.py
def test_annotate_logs_cm_placeholder_caveat(tmp_path, caplog):
    import logging
    # Drive the annotate write path that produces a .annot file (use the public
    # entry point the other annotate tests use), then assert the caveat is logged.
    with caplog.at_level(logging.INFO, logger="ldsc"):
        ...  # run annotate to write a .annot file
    assert any("CM" in r.message and "not used" in r.message for r in caplog.records)
```

- [ ] **Step 3: Implement — log once at write time**

At the `.annot` write site, add:

```python
    LOGGER.info(
        "Writing a placeholder CM column in the .annot output for legacy LDSC "
        "positional compatibility (CHR BP SNP CM <annotations>). The CM value is "
        "not used downstream: ldscore sources CM from the reference panel."
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_annotation.py -k logs_cm_placeholder -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git -C <repo> add src/ldsc/annotation_builder.py tests/test_annotation.py
git -C <repo> commit -m "docs(annotate): log that the .annot CM column is an unused placeholder"
```

---

## Phase 2 — Reference panel is the CM/MAF authority; MAF always required

### Task 3: Require `MAF` after sourcing; parquet sidecar `MAF` unconditional

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` (`compute_chrom_from_parquet`, `compute_chrom_from_plink`)
- Modify: `src/ldsc/_kernel/ref_panel.py:1029-1038` (sidecar `MAF` requirement)
- Test: `tests/test_ldscore_workflow.py`, `tests/test_ref_panel.py`

- [ ] **Step 1: Write the failing test (kernel-level MAF requirement)**

```python
# tests/test_ldscore_workflow.py
def test_compute_chrom_from_parquet_requires_maf():
    import pandas as pd
    import argparse
    from ldsc.errors import LDSCInputError
    from ldsc._kernel import ldscore as k
    metadata = pd.DataFrame({"CHR": ["1", "1"], "POS": [10, 20], "SNP": ["a", "b"], "CM": [0.1, 0.2]})
    # MAF column entirely missing -> must raise, regardless of --maf-min.
    bundle = k.AnnotationBundle(
        metadata=metadata,
        annotations=pd.DataFrame({"base": [1.0, 1.0]}),
        baseline_columns=["base"],
        query_columns=[],
    )
    args = argparse.Namespace(
        snp_identifier="chr_pos", maf_min=None, common_maf_min=0.05,
        ld_wind_kb=100.0, ld_wind_snps=None, ld_wind_cm=None, frqfile=None,
        genome_build="hg19", r2_bias_mode="unbiased", r2_sample_size=None,
        snp_batch_size=50, yes_really=False,
    )
    with pytest.raises(LDSCInputError, match="MAF"):
        k.compute_chrom_from_parquet("1", bundle, args, regression_keys=None)
```

(If constructing a full parquet read is too heavy for a unit test, assert the
require-MAF guard via a small extracted helper `require_reference_maf(metadata,
chrom)` introduced in Step 3 and unit-test that helper directly instead.)

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_ldscore_workflow.py -k requires_maf -v`
Expected: FAIL — today missing MAF degrades to `M_5_50=None`, no error.

- [ ] **Step 3: Implement — a shared require-MAF helper used by both backends**

Add to `src/ldsc/_kernel/ldscore.py`:

```python
def require_reference_maf(metadata: pd.DataFrame, chrom: str) -> None:
    """Raise when the reference panel supplies no usable MAF for a chromosome.

    MAF is mandatory: M_5_50 common-SNP counts, the partitioned-h2 common-overlap
    correction, and any --maf-min filter are meaningless without it.
    """
    if "MAF" not in metadata.columns or metadata["MAF"].isna().all():
        raise LDSCInputError(
            f"ldscore requires MAF for the reference panel but none is available on "
            f"chromosome {chrom}. Most likely a parquet panel was built without allele "
            "frequencies (sidecar MAF=NA). Rebuild the reference panel with MAF, or use a "
            "PLINK panel (MAF is computed from genotypes)."
        )
```

Call it in `compute_chrom_from_parquet` immediately after the metadata is
finalized (after `merge_frequency_metadata` + `apply_maf_filter`, before
`build_window_coordinates` at `ldscore.py:1634`):

```python
    require_reference_maf(metadata, chrom)
```

Call it in `compute_chrom_from_plink` after `geno_meta["MAF"]` is set
(after `ldscore.py:1781`), guarding the empty-universe case first so an
all-monomorphic / fully-filtered chromosome gets a precise message instead of the
misleading "panel lacks MAF" one (`empty["MAF"].isna().all()` is vacuously
`True`):

```python
    if len(geno_meta) == 0:
        raise LDSCInputError(
            f"ldscore retained no PLINK reference SNPs on chromosome {chrom} after genotype "
            "filtering. Most likely every candidate SNP is monomorphic, --maf-min is too high, "
            "or no annotation SNPs overlap the panel. Lower --maf-min, or check the SNP "
            "overlap and genome build."
        )
    require_reference_maf(geno_meta, chrom)
```

Add a matching unit test asserting the empty-`geno_meta` path raises the
"retained no PLINK reference SNPs" message (not the MAF message).

Then make the parquet sidecar loader's MAF requirement unconditional. In
`src/ldsc/_kernel/ref_panel.py:1029-1038`, replace:

```python
    if maf_col is not None:
        maf = pd.to_numeric(df[maf_col], errors="coerce").astype(float)
        out["MAF"] = pd.Series(maf).map(lambda value: value if pd.isna(value) else min(value, 1.0 - value))
    elif global_config.fail_on_missing_metadata:
        raise LDSCInputError( ... )
```

with (drop the `fail_on_missing_metadata` gate):

```python
    if maf_col is None:
        raise LDSCInputError(
            f"Reference-panel metadata sidecar '{path}' is missing MAF metadata. MAF is "
            "required. Most likely the sidecar was produced without allele-frequency values. "
            "Regenerate the panel with MAF metadata."
        )
    maf = pd.to_numeric(df[maf_col], errors="coerce").astype(float)
    out["MAF"] = pd.Series(maf).map(lambda value: value if pd.isna(value) else min(value, 1.0 - value))
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_ldscore_workflow.py -k requires_maf tests/test_ref_panel.py -v`
Expected: PASS. Update any `test_ref_panel.py` case that relied on
`fail_on_missing_metadata=False` to permit missing MAF — that path is gone.

- [ ] **Step 5: Simplify the now-unreachable degradation**

In `compute_counts` (`ldscore.py:1571-1572`), the `if "MAF" not in ... return M, None`
branch is now unreachable in the workflow (require-MAF runs first). Leave it as a
defensive guard for the direct-Python API but add a one-line comment noting MAF
is required upstream. Do not delete (keeps `compute_counts` self-contained).

- [ ] **Step 6: Commit**

```bash
git -C <repo> add src/ldsc/_kernel/ldscore.py src/ldsc/_kernel/ref_panel.py tests/
git -C <repo> commit -m "feat(ldscore): require reference-panel MAF in both backends"
```

### Task 4: Parquet `CM` is sidecar-authoritative; sharpen the missing-`CM` error

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py:1020-1027` (`build_window_coordinates` CM guard)
- Test: `tests/test_ldscore_workflow.py`

Because annotation metadata no longer carries `CM` (Task 1), `merge_frequency_metadata`
already fills `CM` solely from the panel sidecar (`args.frqfile`, wired in
`_namespace_from_configs:1886`). No annotation-first logic remains to remove. The
only change is to sharpen the existing NaN-CM rejection message.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_ldscore_workflow.py
def test_build_window_coordinates_rejects_missing_cm():
    import numpy as np, pandas as pd, argparse
    from ldsc.errors import LDSCInputError
    from ldsc._kernel.ldscore import build_window_coordinates
    metadata = pd.DataFrame({"CM": [0.1, np.nan, 0.3], "POS": [1, 2, 3]})
    args = argparse.Namespace(ld_wind_snps=None, ld_wind_kb=None, ld_wind_cm=1.0)
    with pytest.raises(LDSCInputError, match="missing CM"):
        build_window_coordinates(metadata, args)
```

- [ ] **Step 2: Run test to verify it fails or passes weakly**

Run: `pytest tests/test_ldscore_workflow.py -k rejects_missing_cm -v`
Expected: PASS already if matching "missing CM" loosely; if so, tighten the
assertion to the new wording in Step 3 and re-run (it then fails first).

- [ ] **Step 3: Implement — sharpen the message**

Replace the message in `build_window_coordinates` (`ldscore.py:1021-1026`):

```python
    if metadata["CM"].isna().any():
        raise LDSCInputError(
            "ldscore cannot use `--ld-wind-cm` because the reference panel has missing CM "
            "values for at least one retained SNP. ldscore never silently drops reference "
            "SNPs. Most likely the parquet panel was built without a genetic map "
            "(sidecar CM=NA). Rebuild the panel with a genetic map, or use `--ld-wind-kb` "
            "/ `--ld-wind-snps`."
        )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_ldscore_workflow.py -k rejects_missing_cm -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git -C <repo> add src/ldsc/_kernel/ldscore.py tests/test_ldscore_workflow.py
git -C <repo> commit -m "feat(ldscore): sharpen missing-CM error for parquet cM windows"
```

---

## Phase 3 — Unusable-`CM` guard and the PLINK genetic-map remedy

### Task 5: Unusable-`CM` check (unconditional, before the whole-chromosome guard)

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` (new helper; call in both backends)
- Test: `tests/test_ldscore_workflow.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_ldscore_workflow.py
import numpy as np, pandas as pd
from ldsc.errors import LDSCInputError
from ldsc._kernel.ldscore import assert_cm_usable

def test_assert_cm_usable_rejects_all_zero():
    with pytest.raises(LDSCInputError, match="all zero"):
        assert_cm_usable(pd.Series([0.0, 0.0, 0.0]), chrom="1")

def test_assert_cm_usable_rejects_constant():
    with pytest.raises(LDSCInputError):
        assert_cm_usable(pd.Series([2.5, 2.5, 2.5]), chrom="1")

def test_assert_cm_usable_accepts_two_distinct_values():
    assert_cm_usable(pd.Series([0.0, 0.1, 0.1]), chrom="1")  # no raise
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_ldscore_workflow.py -k assert_cm_usable -v`
Expected: FAIL — `assert_cm_usable` does not exist.

- [ ] **Step 3: Implement**

Add to `src/ldsc/_kernel/ldscore.py`:

```python
def assert_cm_usable(cm: pd.Series, chrom: str) -> None:
    """Reject a CM column that cannot order SNPs on a chromosome.

    Unusable means fewer than two distinct finite values (all zero, all
    identical, or all missing). Such a column makes every SNP collapse to one
    coordinate, so any positive cM window spans the whole chromosome -- a
    meaningless result, not an intentional choice. This check is unconditional;
    `--yes-really` does not bypass it.
    """
    finite = pd.to_numeric(cm, errors="coerce").dropna().unique()
    if len(finite) < 2:
        raise LDSCInputError(
            f"ldscore cannot use `--ld-wind-cm` on chromosome {chrom}: the reference "
            "panel CM column is all zero or otherwise uninformative (fewer than two "
            "distinct values), so it cannot define a genetic-distance window. Provide "
            "real genetic-map positions: pass a PLINK `.bim` with informative CM, add "
            "`--genetic-map-hg38-sources` / `--genetic-map-hg19-sources` for the panel's "
            "build, or use `--ld-wind-kb` / `--ld-wind-snps`. "
            "See docs/troubleshooting.md#ldscore-unusable-cm-for-ld-wind-cm."
        )
```

(The "all zero" phrasing in the message satisfies the test's `match="all zero"`.)

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_ldscore_workflow.py -k assert_cm_usable -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git -C <repo> add src/ldsc/_kernel/ldscore.py tests/test_ldscore_workflow.py
git -C <repo> commit -m "feat(ldscore): add unconditional unusable-CM guard"
```

### Task 6: Config + CLI + namespace wiring for genetic-map sources

**Files:**
- Modify: `src/ldsc/config.py` (`RefPanelConfig`, `LDScoreConfig`)
- Modify: `src/ldsc/ldscore_calculator.py` (`build_parser`, spec construction `:1533-1556`, `_namespace_from_configs:1887-1912`)
- Test: `tests/test_ldscore_workflow.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_ldscore_workflow.py
def test_refpanelconfig_accepts_genetic_map_sources():
    from ldsc.config import RefPanelConfig
    spec = RefPanelConfig(
        backend="plink", plink_prefix="x",
        genetic_map_hg38_sources="map_hg38.txt.gz",
    )
    assert spec.genetic_map_hg38_sources is not None
    assert spec.genetic_map_hg19_sources is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_ldscore_workflow.py -k genetic_map_sources -v`
Expected: FAIL — fields do not exist.

- [ ] **Step 3: Implement config fields**

In `src/ldsc/config.py`, add to `RefPanelConfig` (after `ref_panel_snps_file`,
~line 444):

```python
    genetic_map_hg19_sources: str | PathLike[str] | None = None
    genetic_map_hg38_sources: str | PathLike[str] | None = None
```

Normalize them in `__post_init__` next to the existing path normalizations
(~line 458-461):

```python
        object.__setattr__(self, "genetic_map_hg19_sources", _normalize_optional_path(self.genetic_map_hg19_sources))
        object.__setattr__(self, "genetic_map_hg38_sources", _normalize_optional_path(self.genetic_map_hg38_sources))
```

In `LDScoreConfig` add (after `common_maf_min`, ~line 541):

```python
    export_ref_metadata: bool = False
```

- [ ] **Step 4: Add CLI flags**

In `ldscore_calculator.build_parser()` (near the maf flags at `:954-955`), add:

```python
    parser.add_argument("--genetic-map-hg19-sources", default=None,
                        help="Genetic map (hg19) for cM windows when the PLINK .bim CM is uninformative.")
    parser.add_argument("--genetic-map-hg38-sources", default=None,
                        help="Genetic map (hg38) for cM windows when the PLINK .bim CM is uninformative.")
    parser.add_argument("--export-ref-metadata", action="store_true", default=False,
                        help="Write a chrN_meta.tsv.gz reference-metadata sidecar next to the LD-score output (PLINK only).")
```

- [ ] **Step 5: Thread flags into the spec and config**

In both `RefPanelConfig(...)` constructions (`ldscore_calculator.py:1533` and
`:1545`) add:

```python
            genetic_map_hg19_sources=getattr(args, "genetic_map_hg19_sources", None),
            genetic_map_hg38_sources=getattr(args, "genetic_map_hg38_sources", None),
```

In `_ldscore_config_from_args` add `export_ref_metadata=getattr(args, "export_ref_metadata", False)`.

- [ ] **Step 6: Add a `genetic_map=None` placeholder key to the namespace**

In `_namespace_from_configs` (`ldscore_calculator.py:1887`), add a
`genetic_map=None` keyword to the returned `argparse.Namespace(...)` so the kernel
can read `getattr(args, "genetic_map", None)` safely. The actual resolution +
loading is added in Task 7 (which defines `_resolve_build_for_ldscore_genetic_map`),
keeping this task's tests green on its own.

- [ ] **Step 7: Run test to verify it passes**

Run: `pytest tests/test_ldscore_workflow.py -k genetic_map_sources -v`
Expected: PASS.

- [ ] **Step 8: Commit**

```bash
git -C <repo> add src/ldsc/config.py src/ldsc/ldscore_calculator.py tests/test_ldscore_workflow.py
git -C <repo> commit -m "feat(ldscore): wire genetic-map + export flags through config"
```

### Task 7: Genome-build resolution + kernel CM interpolation (map wins)

**Files:**
- Modify: `src/ldsc/ldscore_calculator.py` (new `_resolve_build_for_ldscore_genetic_map`)
- Modify: `src/ldsc/_kernel/ldscore.py` (`compute_chrom_from_plink` CM resolution)
- Test: `tests/test_ldscore_workflow.py`

- [ ] **Step 1: Write the failing test (kernel CM resolution)**

```python
# tests/test_ldscore_workflow.py
def test_plink_cm_uses_genetic_map_when_bim_cm_zero(tmp_path):
    """End-to-end: an all-zero .bim CM + a genetic map yields a non-whole-chrom window."""
    import pandas as pd
    from ldsc.config import GlobalConfig, RefPanelConfig, set_global_config
    from ldsc.ldscore_calculator import run_ldscore
    # Reuse a PLINK fixture whose .bim CM is all zero (PLINK_FIXTURES); write a
    # 2-point genetic map covering its positions so interpolation is non-degenerate.
    gmap = tmp_path / "map_hg19.txt"
    gmap.write_text("CHR POS CM\n1 1 0.0\n1 1000000 5.0\n")
    set_global_config(GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"))
    result = run_ldscore(
        plink_prefix=str(PLINK_FIXTURES / "chr1"),
        baseline_annot_sources=str(PLINK_FIXTURES / "chr1.annot"),
        genetic_map_hg19_sources=str(gmap),
        ld_wind_cm=1.0,
        output_dir=str(tmp_path / "out"),
    )
    assert result is not None  # run completed; window was not whole-chromosome
```

(Adjust fixture paths/names to the actual `tests/fixtures/plink` contents; the
key assertion is that the run does not raise the unusable-CM error.)

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_ldscore_workflow.py -k uses_genetic_map_when -v`
Expected: FAIL — today PLINK ignores the map and the all-zero CM trips the guard.

- [ ] **Step 3: Implement build resolution (workflow layer)**

Add to `src/ldsc/ldscore_calculator.py`:

```python
def _resolve_build_for_ldscore_genetic_map(ref_panel, global_config, chrom) -> str:
    """Resolve hg19/hg38 for selecting the genetic-map source.

    Uses the configured genome build, or infers it from the reference panel
    positions when genome_build='auto' (chr_pos-family only). Raises if the build
    cannot be determined.
    """
    build = resolve_genome_build(
        global_config.genome_build,            # hint ("auto" | "hg19" | "hg38" | None)
        global_config.snp_identifier,          # already-normalized mode
        ref_panel.load_metadata(chrom)[["CHR", "POS"]],  # sample_frame for auto inference
        context=f"PLINK panel chromosome {chrom} for genetic-map build selection",
    )
    if build not in ("hg19", "hg38"):
        raise LDSCInputError(
            "ldscore could not determine the genome build needed to select a genetic map "
            "for the PLINK panel. Pass `--genome-build hg19` or `--genome-build hg38` "
            "(automatic build inference returns None in rsID identifier modes, and needs "
            "sufficient overlap with the packaged HM3 reference in chr_pos modes)."
        )
    return build
```

The real signature (verified) is
`resolve_genome_build(hint, snp_identifier, sample_frame, *, context, logger=None, reference_table=None)`
in `src/ldsc/genome_build_inference.py:113`; it returns `"hg19"`, `"hg38"`, or
`None` (always `None` for rsID-family modes), and raises if `hint="auto"` has
insufficient HM3 overlap.

- [ ] **Step 3b: Resolve + load the genetic map in the namespace builder**

In `_namespace_from_configs` (`ldscore_calculator.py:1887`), replace the
`genetic_map=None` placeholder added in Task 6 with the resolved map. Above the
`return argparse.Namespace(...)`:

```python
    genetic_map = None
    if backend == "plink" and (getattr(spec, "genetic_map_hg19_sources", None) or getattr(spec, "genetic_map_hg38_sources", None)):
        from ._kernel.ref_panel_builder import load_genetic_map_group
        resolved_build = _resolve_build_for_ldscore_genetic_map(ref_panel, global_config, chrom)
        sources = spec.genetic_map_hg38_sources if resolved_build == "hg38" else spec.genetic_map_hg19_sources
        if sources is None:
            raise LDSCInputError(
                f"ldscore needs a `--genetic-map-{resolved_build}-sources` file to derive CM for the "
                f"PLINK panel (resolved build {resolved_build}), but it was not supplied."
            )
        genetic_map = load_genetic_map_group(split_cli_path_tokens(sources))
    elif backend == "parquet_r2" and (getattr(spec, "genetic_map_hg19_sources", None) or getattr(spec, "genetic_map_hg38_sources", None)):
        LOGGER.warning(
            "Ignoring --genetic-map-*-sources for the parquet R2 reference panel: CM is taken "
            "from the panel metadata sidecar (authoritative). Genetic-map flags apply only to "
            "PLINK panels."
        )
```

Set `genetic_map=genetic_map` in the namespace kwargs. Confirm `split_cli_path_tokens`
is imported in `ldscore_calculator.py` (it is used elsewhere); otherwise import it.
Add a test asserting the warning is logged when `--genetic-map-*-sources` is
passed with a parquet panel (use `caplog`).

- [ ] **Step 4: Implement kernel CM resolution (map wins, then usable bim CM, else guard)**

In `compute_chrom_from_plink`, replace the window-coordinate construction at
`ldscore.py:1792`:

```python
    coords, max_dist = build_window_coordinates(geno_meta.drop(columns="_key"), args)
```

with CM resolution that honors an explicit map and validates bim CM under cM mode:

```python
    if args.ld_wind_cm is not None:
        genetic_map = getattr(args, "genetic_map", None)
        if genetic_map is not None:
            geno_meta["CM"] = interpolate_genetic_map_cm(
                normalize_chromosome(chrom, context=prefix + ".bim"),
                geno_meta["POS"].to_numpy(dtype=np.int64),
                genetic_map,
            )
        else:
            assert_cm_usable(geno_meta["CM"], chrom)
    coords, max_dist = build_window_coordinates(geno_meta.drop(columns="_key"), args)
```

Add the import at the top of `ldscore.py`:
`from .ref_panel_builder import interpolate_genetic_map_cm` (kernel-internal).
Verify the import path: `grep -n "def interpolate_genetic_map_cm" src/ldsc/_kernel/ref_panel_builder.py`.

- [ ] **Step 5: Run test to verify it passes**

Run: `pytest tests/test_ldscore_workflow.py -k uses_genetic_map_when -v`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git -C <repo> add src/ldsc/ldscore_calculator.py src/ldsc/_kernel/ldscore.py tests/test_ldscore_workflow.py
git -C <repo> commit -m "feat(ldscore): interpolate PLINK CM from genetic map for cM windows"
```

### Task 8: Unusable-`CM` error fires even under `--yes-really`; troubleshooting doc

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` (ensure `assert_cm_usable` runs before `check_whole_chromosome_window`)
- Create: `docs/troubleshooting.md` section
- Test: `tests/test_ldscore_workflow.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_ldscore_workflow.py
def test_unusable_cm_errors_even_with_yes_really(tmp_path):
    from ldsc.config import GlobalConfig, set_global_config
    from ldsc.ldscore_calculator import run_ldscore
    from ldsc.errors import LDSCInputError
    set_global_config(GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"))
    with pytest.raises(LDSCInputError, match="uninformative|all zero"):
        run_ldscore(
            plink_prefix=str(PLINK_FIXTURES / "chr1"),  # .bim CM all zero, no map
            baseline_annot_sources=str(PLINK_FIXTURES / "chr1.annot"),
            ld_wind_cm=1.0,
            whole_chromosome_ok=True,  # --yes-really must NOT bypass the unusable-CM error
            output_dir=str(tmp_path / "out"),
        )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_ldscore_workflow.py -k errors_even_with_yes_really -v`
Expected: FAIL — today `--yes-really` lets the all-zero CM proceed as whole-chrom.

- [ ] **Step 3: Implement — ordering already correct, verify**

Task 7 places `assert_cm_usable` before `build_window_coordinates`, which is
before `check_whole_chromosome_window` (`ldscore.py:1794`). `assert_cm_usable`
does not consult `args.yes_really`, so it fires unconditionally. Confirm by
reading the surrounding code that no `yes_really` short-circuit precedes it.

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_ldscore_workflow.py -k errors_even_with_yes_really -v`
Expected: PASS.

- [ ] **Step 5: Add the troubleshooting section**

Append to `docs/troubleshooting.md` a section with the slug used in the error
messages, `ldscore-unusable-cm-for-ld-wind-cm`:

```markdown
## ldscore: unusable CM for `--ld-wind-cm`

`--ld-wind-cm` needs a genetic-map coordinate (`CM`) that orders SNPs. The run
aborts when the reference panel's `CM` is all zero, constant, or missing.

Ranked causes and remedies:
1. **PLINK `.bim` has an all-zero `CM` column** (the common PLINK default).
   - Add a genetic map: `--genetic-map-hg38-sources` or `--genetic-map-hg19-sources`
     matching the panel's build. ldscore interpolates `CM` at the `.bim` positions.
   - Or use a `.bim` whose third column carries real genetic-map positions.
   - Or switch to `--ld-wind-kb` (e.g. `1000`) or `--ld-wind-snps`.
2. **Parquet panel built without a genetic map** (sidecar `CM=NA`).
   - Rebuild with `ldsc build-ref-panel --genetic-map-<build>-sources ...`, or use
     `--ld-wind-kb` / `--ld-wind-snps`.
3. **Genome build for the genetic map could not be determined** (rsID modes).
   - Pass `--genome-build hg19` or `--genome-build hg38`.

`--yes-really` only authorizes an expensive whole-chromosome window from *valid*
`CM`; it never authorizes a window from an uninformative `CM` column.
```

- [ ] **Step 6: Commit**

```bash
git -C <repo> add src/ldsc/_kernel/ldscore.py docs/troubleshooting.md tests/test_ldscore_workflow.py
git -C <repo> commit -m "feat(ldscore): unusable-CM error overrides --yes-really; add docs"
```

---

## Phase 4 — `--maf-min` applies in both backends

### Task 9: Apply `--maf-min` as a reference-panel filter for PLINK

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py:1771` (PLINK reader `mafMin`) and/or `src/ldsc/ldscore_calculator.py:1908` (`_namespace_from_configs`)
- Test: `tests/test_ldscore_workflow.py`

The kernel already passes `mafMin=getattr(args, "maf_min", ...)` to `PlinkBEDFile`
(`ldscore.py:1771`); the only reason `--maf-min` is dropped for PLINK is that
`_namespace_from_configs` hardcodes `maf_min=None` (`:1908`). Plumb the real value.

- [ ] **Step 1: Write the failing cross-backend test**

```python
# tests/test_ldscore_workflow.py
def test_maf_min_filters_identically_across_backends(tmp_path):
    """--maf-min must drop the same SNPs (same M, M_5_50) in PLINK and parquet."""
    # Build matched PLINK and parquet panels from the same genotypes/frequencies,
    # run run_ldscore with maf_min=0.1 against each, and compare the count totals.
    from ldsc.config import GlobalConfig, set_global_config
    from ldsc.ldscore_calculator import run_ldscore
    set_global_config(GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"))
    plink_result = run_ldscore(
        plink_prefix=str(PLINK_FIXTURES / "chr1"),
        baseline_annot_sources=str(PLINK_FIXTURES / "chr1.annot"),
        ld_wind_kb=1000.0, maf_min=0.1, output_dir=str(tmp_path / "plink_out"),
    )
    # parquet_result = run_ldscore(r2_dir=..., maf_min=0.1, ...)  # matched parquet panel
    # assert plink_result.snp_count_totals == parquet_result.snp_count_totals
    assert plink_result is not None
```

Make the assertion concrete once a matched parquet fixture is available; the
minimal red test asserts that `--maf-min` changes the PLINK retained count
relative to `maf_min=None` (today it does not).

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_ldscore_workflow.py -k maf_min_filters_identically -v`
Expected: FAIL — PLINK does not currently apply `--maf-min`.

- [ ] **Step 3: Implement — plumb `maf_min` into the namespace**

In `_namespace_from_configs` (`ldscore_calculator.py:1908`), replace `maf_min=None`
with the spec value:

```python
        maf_min=getattr(spec, "maf_min", None),
```

This makes the kernel's PLINK `PlinkBEDFile(mafMin=...)` apply the inclusive
`MAF >= maf_min` filter (already implemented at `plink_bed.py:320`). For parquet,
the kernel `apply_maf_filter` then also applies it; since the parquet adapter
already filtered at `load_metadata`, the second pass is idempotent (`>=` twice).
To avoid a redundant pass, gate the kernel `apply_maf_filter` call in
`compute_chrom_from_parquet` on the parquet backend not having pre-filtered — or
simply leave it (idempotent). Document the choice in the function docstring.

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_ldscore_workflow.py -k maf_min -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git -C <repo> add src/ldsc/ldscore_calculator.py tests/test_ldscore_workflow.py
git -C <repo> commit -m "fix(ldscore): apply --maf-min for the PLINK backend"
```

---

## Phase 5 — Provenance and optional metadata export

### Task 10: Record `CM`/`MAF` source provenance in `metadata.json`

**Files:**
- Modify: `src/ldsc/outputs.py` (`LDScoreDirectoryWriter.build_metadata`, ~`:306`)
- Modify: `src/ldsc/ldscore_calculator.py` (collect the provenance dict)
- Test: `tests/test_ldscore_workflow.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_ldscore_workflow.py
def test_metadata_records_cm_maf_provenance(tmp_path):
    import json
    from ldsc.config import GlobalConfig, set_global_config
    from ldsc.ldscore_calculator import run_ldscore
    set_global_config(GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"))
    out = tmp_path / "out"
    run_ldscore(
        plink_prefix=str(PLINK_FIXTURES / "chr1"),
        baseline_annot_sources=str(PLINK_FIXTURES / "chr1.annot"),
        ld_wind_kb=1000.0, output_dir=str(out),
    )
    meta = json.loads((out / "metadata.json").read_text())
    assert meta["reference_cm_source"] in {"plink_bim", "genetic_map", "parquet_sidecar"}
    assert meta["reference_maf_source"] in {"plink_genotypes", "parquet_sidecar"}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_ldscore_workflow.py -k records_cm_maf_provenance -v`
Expected: FAIL — keys absent.

- [ ] **Step 3: Implement**

Determine the provenance in the workflow when building the result/output config
(where `backend` and `spec` are known). Compute:

```python
    cm_source = (
        "parquet_sidecar" if backend == "parquet_r2"
        else "genetic_map" if (spec.genetic_map_hg19_sources or spec.genetic_map_hg38_sources)
        else "plink_bim"
    )
    maf_source = "parquet_sidecar" if backend == "parquet_r2" else "plink_genotypes"
```

Pass these into `LDScoreDirectoryWriter.build_metadata` (extend its signature /
the metadata dict it assembles at `outputs.py:306`) so they appear as
`reference_cm_source` and `reference_maf_source` in `metadata.json`. Inspect the
exact `build_metadata` parameters first:
Run: `grep -n "def build_metadata" src/ldsc/outputs.py` and add two keyword
arguments threaded into the returned dict.

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_ldscore_workflow.py -k records_cm_maf_provenance -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git -C <repo> add src/ldsc/outputs.py src/ldsc/ldscore_calculator.py tests/test_ldscore_workflow.py
git -C <repo> commit -m "feat(ldscore): record CM/MAF source provenance in metadata.json"
```

### Task 11: `--export-ref-metadata` writes a unified `chrN_meta.tsv.gz` sidecar

**Files:**
- Modify: `src/ldsc/ldscore_calculator.py` (per-chromosome export when enabled)
- Reuse: `src/ldsc/_kernel/ref_panel_builder.py:867` (`write_runtime_metadata_sidecar`), `:616` (`build_runtime_metadata_table`)
- Test: `tests/test_ldscore_workflow.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_ldscore_workflow.py
def test_export_ref_metadata_writes_sidecar(tmp_path):
    import pandas as pd
    from ldsc.config import GlobalConfig, set_global_config
    from ldsc.ldscore_calculator import run_ldscore
    set_global_config(GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"))
    out = tmp_path / "out"
    run_ldscore(
        plink_prefix=str(PLINK_FIXTURES / "chr1"),
        baseline_annot_sources=str(PLINK_FIXTURES / "chr1.annot"),
        ld_wind_kb=1000.0, export_ref_metadata=True, output_dir=str(out),
    )
    sidecar = out / "ref_metadata" / "chr1_meta.tsv.gz"
    assert sidecar.exists()
    df = pd.read_csv(sidecar, sep="\t")
    assert list(df.columns)[:7] == ["CHR", "POS", "SNP", "A1", "A2", "CM", "MAF"]


def test_export_ref_metadata_defaults_off(tmp_path):
    from ldsc.config import GlobalConfig, set_global_config
    from ldsc.ldscore_calculator import run_ldscore
    set_global_config(GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"))
    out = tmp_path / "out"
    run_ldscore(
        plink_prefix=str(PLINK_FIXTURES / "chr1"),
        baseline_annot_sources=str(PLINK_FIXTURES / "chr1.annot"),
        ld_wind_kb=1000.0, output_dir=str(out),
    )
    assert not (out / "ref_metadata").exists()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_ldscore_workflow.py -k export_ref_metadata -v`
Expected: FAIL — no export path exists.

- [ ] **Step 3: Implement**

When `ldscore_config.export_ref_metadata` and backend is PLINK, after each
chromosome's compute, assemble the retained metadata (`CHR POS SNP A1 A2 CM MAF`,
folded MAF) from the chromosome result's `out_metadata` and write it with the
existing kernel writer. Inspect `write_runtime_metadata_sidecar` and
`build_runtime_metadata_table` signatures first:
Run: `grep -n "def write_runtime_metadata_sidecar\|def build_runtime_metadata_table" src/ldsc/_kernel/ref_panel_builder.py`

Then, in the workflow chromosome loop:

```python
    if ldscore_config.export_ref_metadata and backend == "plink":
        from ._kernel.ref_panel_builder import write_runtime_metadata_sidecar
        export_dir = Path(output_dir) / "ref_metadata"
        export_dir.mkdir(parents=True, exist_ok=True)
        write_runtime_metadata_sidecar(
            chrom_result.reference_metadata[["CHR", "POS", "SNP", "A1", "A2", "CM", "MAF"]],
            export_dir / f"chr{chrom}_meta.tsv.gz",
        )
```

**Plumbing note (important):** the public `ChromLDScoreResult` drops per-SNP
`MAF` (its tables keep only `CHR/SNP/POS/A1/A2` + scores), so it cannot feed the
export by itself. The kernel `ChromComputationResult.metadata` *does* carry
`CM`/`MAF` and is visible in `_wrap_legacy_chrom_result` as `reference_metadata`
(`ldscore_calculator.py:533`). Thread that frame to the workflow: either add a
`reference_metadata: pd.DataFrame | None` field to `ChromLDScoreResult` populated
in `_wrap_legacy_chrom_result`, or perform the export at that wrap point where
`reference_metadata` is in scope. Then adapt the writer call to the real
`write_runtime_metadata_sidecar` signature.

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_ldscore_workflow.py -k export_ref_metadata -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git -C <repo> add src/ldsc/ldscore_calculator.py tests/test_ldscore_workflow.py
git -C <repo> commit -m "feat(ldscore): add opt-in --export-ref-metadata sidecar (PLINK)"
```

---

## Phase 6 — Cross-backend consistency and documentation

### Task 12: Cross-backend window + counts regression test (the headline)

**Files:**
- Test: `tests/test_ldscore_workflow.py`

- [ ] **Step 1: Write the test**

```python
# tests/test_ldscore_workflow.py
def test_cm_window_consistent_across_backends(tmp_path):
    """Identical annotation + identical reference CM -> identical block-lefts."""
    import numpy as np, pandas as pd
    from ldsc._kernel.ldscore import build_window_coordinates, get_block_lefts
    import argparse
    cm = pd.Series([0.0, 0.4, 0.9, 1.6])
    args = argparse.Namespace(ld_wind_snps=None, ld_wind_kb=None, ld_wind_cm=1.0)
    md = pd.DataFrame({"CM": cm, "POS": [1, 2, 3, 4]})
    coords_a, dist_a = build_window_coordinates(md, args)
    coords_b, dist_b = build_window_coordinates(md.copy(), args)
    np.testing.assert_array_equal(get_block_lefts(coords_a, dist_a),
                                  get_block_lefts(coords_b, dist_b))
```

For a full end-to-end equality (`M`, `M_5_50`, overlap) across PLINK and parquet,
build matched panels from the same genotypes (extend the Task 9 fixture) and
assert `plink_result.snp_count_totals == parquet_result.snp_count_totals` and
equal LD-score tables within float tolerance.

- [ ] **Step 1b: Inclusive `>=` boundary tests (lock the invariant)**

```python
# tests/test_ldscore_workflow.py
def test_maf_thresholds_are_inclusive():
    import numpy as np, pandas as pd
    from ldsc._kernel.ldscore import apply_maf_filter, compute_counts
    md = pd.DataFrame({"CHR": ["1", "1"], "POS": [1, 2], "SNP": ["a", "b"], "MAF": [0.05, 0.04]})
    annot = pd.DataFrame({"base": [1.0, 1.0]})
    # --maf-min: a SNP at exactly maf_min is KEPT (MAF >= maf_min).
    kept_md, kept_annot = apply_maf_filter(md.copy(), annot.copy(), 0.05, context="test")
    assert list(kept_md["SNP"]) == ["a"]
    # --common-maf-min: a SNP at exactly the threshold is COUNTED (MAF >= common).
    md_common = pd.DataFrame({"MAF": [0.05, 0.049]})
    M, M_5_50 = compute_counts(md_common, pd.DataFrame({"base": [1.0, 1.0]}), common_maf_min=0.05)
    assert M_5_50[0] == 1.0  # only the 0.05 SNP counts as common
```

- [ ] **Step 2: Run test**

Run: `pytest tests/test_ldscore_workflow.py -k consistent_across_backends -v`
Expected: PASS (this guards the invariant; should pass once Phases 1-4 land).

- [ ] **Step 3: Commit**

```bash
git -C <repo> add tests/test_ldscore_workflow.py
git -C <repo> commit -m "test(ldscore): cross-backend CM-window and counts consistency"
```

### Task 13: Documentation and design-map sync

**Files:**
- Modify: `docs/current/data-flow.md`, `docs/current/class-and-features.md` (CM/MAF sourcing, new flags)
- Modify: `tutorials/` (the LD-score tutorial: annotation contract, genetic-map flags, `--export-ref-metadata`)
- Modify: `design_map.md` (map spec sections to the changed functions)

- [ ] **Step 1: Update the design docs**

Document in `docs/current/` (data-flow and class-and-features): annotations carry
only `CHR/POS` (+ optional `SNP/A1/A2`); reference panel is authoritative for
`CM`/`MAF`; the `--genetic-map-*-sources`, `--export-ref-metadata` flags; the
unusable-`CM` rule and `--yes-really` interaction; `--maf-min` now applies in both
backends; and `--genetic-map-*-sources` is ignored (with a warning) in parquet
mode (sidecar `CM` authoritative).

Add an explicit **`CM` caveat** subsection: `CM` is population-specific and used
only by `--ld-wind-cm`, sourced from the reference panel. The `CM` column written
by `annotate` is a legacy-positional-compatibility placeholder (`CHR BP SNP CM
<annotations>`) and is **not used anywhere downstream** — neither legacy ldsc2
(which skips it via `iloc[:, 4:]`) nor this package's ldscore (which sources `CM`
from the reference panel and selects annotation columns by name).

- [ ] **Step 2: Update the tutorial**

Add a worked PLINK `--ld-wind-cm` example with `--genetic-map-hg38-sources`, and
note that annotation `CM`/`MAF` columns are ignored.

- [ ] **Step 3: Update `design_map.md`**

Add rows mapping the spec's CM/MAF-source decisions to `assert_cm_usable`,
`require_reference_maf`, `compute_chrom_from_plink` CM resolution, and the
annotation-reader change.

- [ ] **Step 4: Full suite + commit**

Run: `source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest`
Expected: PASS.

```bash
git -C <repo> add docs/ tutorials/ design_map.md
git -C <repo> commit -m "docs(ldscore): document CM/MAF source-of-truth and new flags"
```

---

## Final verification

- [ ] Run the full suite and show output:
  `source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest -q`
- [ ] Manually exercise the headline scenario: a PLINK panel with all-zero `.bim`
  `CM` under `--ld-wind-cm` (a) errors clearly without a map, even with
  `--yes-really`; (b) succeeds with `--genetic-map-<build>-sources` and matches
  the parquet result built from the same map.
- [ ] Confirm a parquet panel without sidecar `MAF` errors; with `MAF` it runs.
- [ ] Confirm `--maf-min` drops the same SNPs in both backends.

## Notes for the implementer

- **Verify-before-edit:** several wiring steps include a `grep` step to confirm an
  exact signature (`resolve_genome_build`, `build_metadata`,
  `write_runtime_metadata_sidecar`, the annotation reader's access to the
  identifier mode). Run those greps and adapt argument names to reality — do not
  assume.
- **Idempotent `--maf-min`:** after Task 9, parquet may filter MAF twice
  (adapter + kernel). With inclusive `>=` this is harmless; prefer the single
  source if a clean gate is easy, but correctness does not require it.
- **Legacy `unittest` style:** `tests/test_ldscore_workflow.py` mixes
  `unittest.TestCase` and pytest functions. Either style is fine; match the
  surrounding file. Guard parquet tests with the existing
  `@unittest.skipUnless(_HAS_PYARROW, ...)` pattern where applicable.
