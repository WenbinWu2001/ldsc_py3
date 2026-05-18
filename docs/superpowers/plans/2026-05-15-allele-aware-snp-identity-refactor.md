# Allele-Aware SNP Identity Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement package-wide SNP identity harmonization with four exact SNP identifier modes, allele-aware merge keys, drop-all duplicate policy, conservative multi-allelic filtering, artifact provenance validation, and regression-only identity downgrade.

**Architecture:** Add one shared internal identity service in `src/ldsc/_kernel/snp_identity.py`, then route every workflow merge, restriction, artifact validation, and row-alignment check through that service. Keep low-level file parsing in `src/ldsc/_kernel/identifiers.py`; the new identity module owns mode semantics, allele normalization, effective keys, drop reports, and compatibility decisions.

**Tech Stack:** Python 3.11, pandas, numpy, pyarrow parquet metadata, pytest/unittest, existing LDSC workflow wrappers and kernel modules.

---

## Source Design

Primary spec:

- `docs/superpowers/specs/2026-05-15-allele-aware-snp-identity-design.md`

Clarified decisions applied by this plan:

- Public SNP identifier modes are exactly `rsid`, `rsid_allele_aware`, `chr_pos`, `chr_pos_allele_aware`.
- Default mode is `chr_pos_allele_aware`.
- Mode aliases such as `rsID`, `SNPID`, `rsid_alleles`, and `chr_pos_alleles` are not accepted as `snp_identifier` values.
- Column aliases remain accepted for input headers.
- In base modes, alleles are never inspected. If rows duplicate the active base merge key, drop reason is `duplicate_identity`.
- Base modes are not "allele-aware with relaxed validation"; they are fully allele-blind. `rsid` uses only `SNP`; `chr_pos` uses only `CHR:POS`. Missing, invalid, non-SNP, strand-ambiguous, or multi-base allele columns must not change retention or drop reasons in base modes.
- External raw R2 parquet support remains available only for `rsid` and `chr_pos`; allele-aware modes require package-built canonical R2 parquet with endpoint allele columns.
- Latest Q4 decision: external raw R2 support intentionally stays base-mode only.
- Latest Q6 decision: base-mode duplicate cleanup uses `duplicate_identity` and never checks whether allele columns exist.
- External annotation inputs do not require alleles. If allele columns are present, they can participate in allele-aware matching; otherwise annotation matching uses the base key after the reference universe is cleaned.
- `--allow-identity-downgrade` is regression-only and is represented in the Python API as `RegressionConfig.allow_identity_downgrade`.
- Artifact schema version for this refactor is integer `1`.
- Minimal identity provenance on reloadable identity artifacts is `schema_version`, `artifact_type`, `snp_identifier`, and `genome_build`; runtime events such as downgrade decisions are logged.

## Global Invariants

These invariants must be true after every task:

1. Drop-duplication policy always applies on the effective merge key for the active mode.
2. Artifact-like tables use drop-all for duplicate effective-key clusters.
3. Restriction files are set-like: duplicate effective restriction keys collapse to one key.
4. `rsid` and `chr_pos` modes preserve allele-free support and do not inspect allele columns.
5. In base modes, `missing_allele`, `invalid_allele`, `strand_ambiguous_allele`, and `multi_allelic_base_key` are unreachable drop reasons.
6. `rsid_allele_aware` and `chr_pos_allele_aware` require usable `A1/A2` for sumstats, reference-panel artifacts, R2 parquet endpoints, and LD-score artifacts.
7. Allele-aware identity modes drop missing, invalid, non-SNP, identical, and strand-ambiguous allele pairs.
8. Allele-aware identity modes drop all rows in multi-allelic base-key clusters, where the base key is `SNP` for rsID-family modes and `CHR:POS` for coordinate-family modes.
9. Liftover collision filtering ignores allele set and checks only `CHR:POS`.
10. Coordinates lift; allele sets stay unchanged.
11. Ordered `A1/A2` affect signed-statistic orientation only, never the unordered identity token.

## Base-Mode Wall

The implementation must keep a hard branch between base and allele-aware modes:

- `base_key_series()` is the only identity-key input used by `rsid` and `chr_pos`.
- `allele_set_series()`, `_allele_failure_reason()`, `normalize_allele_set()`, strand-ambiguity checks, and multi-allelic base-key cleanup are called only after `is_allele_aware_mode(mode)` is true.
- In restriction cleanup, base modes return base keys before any allele cleanup branch, even if the restriction file contains `A1/A2`, `REF/ALT`, or other allele aliases.
- In regression downgrade, the effective base mode stays allele-blind; it must not reuse the recorded allele-aware mode for duplicate cleanup.
- Base modes may preserve allele columns for downstream orientation or output convenience, but those columns are passive data and never identity inputs.

## File Map

| File | Responsibility |
| --- | --- |
| `src/ldsc/_kernel/snp_identity.py` | New central identity service: exact mode validation, mode families, base/effective keys, allele QC, artifact cleanup, restriction cleanup, compatibility/downgrade decisions, sidecar schemas. |
| `src/ldsc/_kernel/identifiers.py` | Keep restriction file parsing and delimiter/header handling; delegate key semantics and allele-aware restriction cleanup to `snp_identity.py`. |
| `src/ldsc/column_inference.py` | Keep column alias registry; remove `snp_identifier` mode aliases; re-export exact mode normalizer for transitional imports. Add allele column specs for external restrictions/annotations and parquet endpoints. |
| `src/ldsc/config.py` | Expand `SNPIdentifierMode`, update defaults, add `RegressionConfig.allow_identity_downgrade`, update compatibility validation. |
| `src/ldsc/_row_alignment.py` | Use effective identity keys and optional `A1/A2` columns for baseline/query row alignment. |
| `src/ldsc/sumstats_munger.py` | Enforce allele-aware inputs, write/read minimal identity metadata, run global identity cleanup after liftover, extend dropped-SNP sidecar schema. |
| `src/ldsc/_kernel/sumstats_munger.py` | Use four modes in parser and chunk filtering; remove rsID-only duplicate keep-first behavior and call shared cleanup at final identity point. |
| `src/ldsc/ref_panel_builder.py` | Apply shared identity cleanup after restrictions and after liftover synchronization; reject rsID-family liftover; add allele-aware metadata and sidecar rows. |
| `src/ldsc/_kernel/ref_panel_builder.py` | Emit endpoint allele columns in R2 parquet; write identity schema metadata to parquet and runtime metadata sidecars; write `A1/A2` in runtime metadata sidecars. |
| `src/ldsc/_kernel/ref_panel.py` | Validate package-built R2 and metadata artifacts; reject old package artifacts; reject external raw R2 in allele-aware modes. |
| `src/ldsc/_kernel/ldscore.py` | Use identity keys in annotation/reference/R2 matching; validate canonical R2 endpoint alleles in allele-aware modes; write LD-score rows with `A1/A2`. |
| `src/ldsc/ldscore_calculator.py` | Preserve `A1/A2` through `ChromLDScoreResult`/`LDScoreResult`, manifest, and split baseline/query tables. |
| `src/ldsc/outputs.py` | Add minimal LD-score manifest fields and allele-aware table validation. |
| `src/ldsc/regression_runner.py` | Exact-mode compatibility by default; implement regression-only downgrade; orient h2 sumstats to LD-score alleles; keep rg ordered harmonization after merge. |
| `src/ldsc/annotation_builder.py` | Accept annotation allele columns without requiring them; validate annotation rows by available precision. |
| `tests/test_snp_identity.py` | New focused identity-service tests. |
| Existing workflow tests | Update expectations for modes, defaults, duplicate drops, sidecars, artifact metadata, and regression compatibility. |

## Pre-Flight

- [ ] **Step 1: Confirm branch and local edits**

Run:

```bash
git status --short --branch
```

Expected: the first line starts with `## restructure`; extra modified or untracked files are allowed only when they are unrelated to the current task and are not overwritten.

If `docs/superpowers/specs/2026-05-15-allele-aware-snp-identity-design.md` is modified, treat it as user-owned unless the current task explicitly edits the spec.

- [ ] **Step 2: Re-read the design spec**

Run:

```bash
sed -n '1,430p' docs/superpowers/specs/2026-05-15-allele-aware-snp-identity-design.md
```

Expected: the file states the four public modes, duplicate policy, allele-aware drops, no old artifact compatibility, regression-only downgrade, and raw R2 base-mode support.

- [ ] **Step 3: Baseline targeted test inventory**

Run:

```bash
pytest -q tests/test_config_identifiers.py tests/test_column_inference.py tests/test_sumstats_munger.py tests/test_sumstats_liftover.py tests/test_ref_panel_builder.py tests/test_ref_panel.py tests/test_ldscore_workflow.py tests/test_regression_workflow.py tests/test_plink_io.py
```

Expected before implementation: current branch behavior may pass or fail independently of this plan. Record failures before editing so new failures are distinguishable.

---

## Task 1: Add The Shared SNP Identity Service

**Files:**

- Create: `src/ldsc/_kernel/snp_identity.py`
- Create: `tests/test_snp_identity.py`
- Modify: `src/ldsc/column_inference.py`
- Modify: `src/ldsc/_kernel/identifiers.py`

- [ ] **Step 1: Write identity-service tests**

Create `tests/test_snp_identity.py` with these test cases:

```python
import pandas as pd
import pytest

from ldsc._kernel.snp_identity import (
    clean_identity_artifact_table,
    collapse_restriction_identity_keys,
    effective_merge_key_series,
    identity_mode_family,
    normalize_allele_set,
    normalize_snp_identifier_mode,
)


def test_accepts_only_four_public_modes():
    assert normalize_snp_identifier_mode("rsid") == "rsid"
    assert normalize_snp_identifier_mode("rsid_allele_aware") == "rsid_allele_aware"
    assert normalize_snp_identifier_mode("chr_pos") == "chr_pos"
    assert normalize_snp_identifier_mode("chr_pos_allele_aware") == "chr_pos_allele_aware"
    for bad in ["rsID", "SNPID", "snp_id", "chrpos", "rsid_alleles", "chr_pos_alleles"]:
        with pytest.raises(ValueError, match="snp_identifier"):
            normalize_snp_identifier_mode(bad)


def test_mode_family_and_base_mode():
    assert identity_mode_family("rsid") == "rsid"
    assert identity_mode_family("rsid_allele_aware") == "rsid"
    assert identity_mode_family("chr_pos") == "chr_pos"
    assert identity_mode_family("chr_pos_allele_aware") == "chr_pos"


def test_unordered_strand_aware_allele_set_normalization():
    assert normalize_allele_set("A", "C") == "A:C"
    assert normalize_allele_set("C", "A") == "A:C"
    assert normalize_allele_set("T", "G") == "A:C"
    assert normalize_allele_set("G", "T") == "A:C"


def test_allele_set_rejects_invalid_pairs():
    for a1, a2 in [(None, "A"), ("A", ""), ("A", "A"), ("A", "AT"), ("A", "T"), ("C", "G")]:
        with pytest.raises(ValueError):
            normalize_allele_set(a1, a2)


def test_effective_keys_for_all_modes():
    frame = pd.DataFrame(
        {
            "CHR": ["1", "1"],
            "POS": [100, 200],
            "SNP": ["rs1", "rs2"],
            "A1": ["A", "C"],
            "A2": ["C", "T"],
        }
    )
    assert effective_merge_key_series(frame, "rsid").tolist() == ["rs1", "rs2"]
    assert effective_merge_key_series(frame, "chr_pos").tolist() == ["1:100", "1:200"]
    assert effective_merge_key_series(frame, "rsid_allele_aware").tolist() == ["rs1:A:C", "rs2:A:G"]
    assert effective_merge_key_series(frame, "chr_pos_allele_aware").tolist() == ["1:100:A:C", "1:200:A:G"]


def test_base_modes_drop_duplicate_identity_without_inspecting_alleles():
    frame = pd.DataFrame(
        {
            "CHR": ["1", "1", "1"],
            "POS": [100, 100, 200],
            "SNP": ["rs1", "rs1", "rs2"],
            "A1": ["A", "A", "A"],
            "A2": ["C", "G", "T"],
        }
    )
    result = clean_identity_artifact_table(frame, "rsid", context="test", stage="unit")
    assert result.cleaned["SNP"].tolist() == ["rs2"]
    assert result.dropped["reason"].tolist() == ["duplicate_identity", "duplicate_identity"]
    assert result.dropped["base_key"].tolist() == ["rs1", "rs1"]


def test_base_modes_keep_singletons_even_when_alleles_are_invalid_or_ambiguous():
    frame = pd.DataFrame(
        {
            "CHR": ["1", "1", "1"],
            "POS": [100, 200, 300],
            "SNP": ["rs1", "rs2", "rs3"],
            "A1": ["A", None, "C"],
            "A2": ["T", "AT", "G"],
        }
    )
    result = clean_identity_artifact_table(frame, "chr_pos", context="test", stage="unit")
    assert result.cleaned["SNP"].tolist() == ["rs1", "rs2", "rs3"]
    assert result.dropped.empty


def test_allele_aware_drops_multi_allelic_base_key_cluster():
    frame = pd.DataFrame(
        {
            "CHR": ["1", "1", "1"],
            "POS": [100, 100, 200],
            "SNP": ["rs1", "rs1", "rs2"],
            "A1": ["A", "A", "A"],
            "A2": ["C", "G", "C"],
        }
    )
    result = clean_identity_artifact_table(frame, "chr_pos_allele_aware", context="test", stage="unit")
    assert result.cleaned["SNP"].tolist() == ["rs2"]
    assert result.dropped["reason"].tolist() == ["multi_allelic_base_key", "multi_allelic_base_key"]


def test_restriction_duplicate_keys_collapse_instead_of_drop_all():
    frame = pd.DataFrame(
        {
            "CHR": ["1", "1", "1"],
            "POS": [100, 100, 200],
            "SNP": ["rs1", "rs1", "rs2"],
        }
    )
    keys = collapse_restriction_identity_keys(frame, "chr_pos", context="restriction")
    assert keys.match_kind == "base"
    assert keys.keys == {"1:100", "1:200"}
```

- [ ] **Step 2: Run the new tests and verify they fail**

Run:

```bash
pytest -q tests/test_snp_identity.py
```

Expected: import failure for `ldsc._kernel.snp_identity`.

- [ ] **Step 3: Implement `snp_identity.py` public constants and mode helpers**

Add these exact public names:

```python
from __future__ import annotations

from dataclasses import dataclass
import logging
from typing import Literal

import numpy as np
import pandas as pd

from .._coordinates import CHR_POS_KEY_COLUMN, build_chr_pos_key_frame

SNPIdentifierMode = Literal["rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"]
SNP_IDENTIFIER_MODES: tuple[str, str, str, str] = ("rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware")
RSID_FAMILY_MODES = frozenset({"rsid", "rsid_allele_aware"})
CHR_POS_FAMILY_MODES = frozenset({"chr_pos", "chr_pos_allele_aware"})
ALLELE_AWARE_MODES = frozenset({"rsid_allele_aware", "chr_pos_allele_aware"})
BASE_MODES = frozenset({"rsid", "chr_pos"})
SCHEMA_VERSION = 1

IDENTITY_DROP_COLUMNS = [
    "CHR",
    "SNP",
    "source_pos",
    "target_pos",
    "reason",
    "base_key",
    "identity_key",
    "allele_set",
    "stage",
]

IDENTITY_DROP_REASONS = frozenset(
    {
        "missing_allele",
        "invalid_allele",
        "strand_ambiguous_allele",
        "multi_allelic_base_key",
        "duplicate_identity",
    }
)

COMPLEMENT = str.maketrans({"A": "T", "T": "A", "C": "G", "G": "C"})


@dataclass(frozen=True)
class IdentityCleanupResult:
    cleaned: pd.DataFrame
    dropped: pd.DataFrame


@dataclass(frozen=True)
class RestrictionIdentityKeys:
    keys: set[str]
    match_kind: Literal["base", "identity"]
    dropped: pd.DataFrame
    n_input_rows: int
    n_retained_keys: int


def normalize_snp_identifier_mode(value: str) -> str:
    if value in SNP_IDENTIFIER_MODES:
        return value
    allowed = ", ".join(SNP_IDENTIFIER_MODES)
    raise ValueError(f"Unsupported snp_identifier mode {value!r}; expected one of: {allowed}.")


def identity_mode_family(mode: str) -> Literal["rsid", "chr_pos"]:
    mode = normalize_snp_identifier_mode(mode)
    return "rsid" if mode in RSID_FAMILY_MODES else "chr_pos"


def identity_base_mode(mode: str) -> Literal["rsid", "chr_pos"]:
    return identity_mode_family(mode)


def is_allele_aware_mode(mode: str) -> bool:
    return normalize_snp_identifier_mode(mode) in ALLELE_AWARE_MODES
```

- [ ] **Step 4: Implement allele normalization and error classification**

Rules:

- Convert alleles to uppercase stripped strings.
- Accept only single-base `A/C/G/T`.
- Reject identical pairs as `invalid_allele`.
- Reject strand-ambiguous unordered pairs `A:T` and `C:G` as `strand_ambiguous_allele`.
- Convert both input alleles and their complements to the lexicographically smaller unordered representation. This makes `A:C`, `C:A`, `T:G`, and `G:T` all become `A:C`.

Implementation shape:

```python
def _clean_allele(value: object) -> str | None:
    if value is None or pd.isna(value):
        return None
    token = str(value).strip().upper()
    return token or None


def _allele_failure_reason(a1: object, a2: object) -> str | None:
    left = _clean_allele(a1)
    right = _clean_allele(a2)
    if left is None or right is None:
        return "missing_allele"
    if len(left) != 1 or len(right) != 1 or left not in "ACGT" or right not in "ACGT" or left == right:
        return "invalid_allele"
    unordered = frozenset({left, right})
    if unordered in {frozenset({"A", "T"}), frozenset({"C", "G"})}:
        return "strand_ambiguous_allele"
    return None


def normalize_allele_set(a1: object, a2: object) -> str:
    reason = _allele_failure_reason(a1, a2)
    if reason is not None:
        raise ValueError(reason)
    left = _clean_allele(a1)
    right = _clean_allele(a2)
    observed = ":".join(sorted((left, right)))
    complemented = ":".join(sorted((left.translate(COMPLEMENT), right.translate(COMPLEMENT))))
    return min(observed, complemented)
```

- [ ] **Step 5: Implement base/effective key builders**

Expected behavior:

- `rsid`: key from canonical `SNP`.
- `chr_pos`: key from normalized `CHR:POS`.
- allele-aware modes append normalized allele set to base key.
- The function requires canonical columns; external alias inference stays outside this module.

Implementation shape:

```python
def base_key_series(frame: pd.DataFrame, mode: str, *, context: str) -> pd.Series:
    mode = normalize_snp_identifier_mode(mode)
    if identity_mode_family(mode) == "rsid":
        if "SNP" not in frame.columns:
            raise ValueError(f"{context} is missing SNP required for snp_identifier={mode!r}.")
        return frame["SNP"].astype("string").astype(object)
    missing = [column for column in ("CHR", "POS") if column not in frame.columns]
    if missing:
        raise ValueError(f"{context} is missing {missing} required for snp_identifier={mode!r}.")
    keyed, _report = build_chr_pos_key_frame(frame.loc[:, ["CHR", "POS"]].copy(), context=context, drop_missing=True)
    keys = pd.Series(pd.NA, index=frame.index, dtype="object")
    keys.loc[keyed.index] = keyed[CHR_POS_KEY_COLUMN].astype(str)
    return keys


def allele_set_series(frame: pd.DataFrame, *, context: str) -> tuple[pd.Series, pd.Series]:
    missing = [column for column in ("A1", "A2") if column not in frame.columns]
    if missing:
        raise ValueError(
            f"{context} requires A1/A2 allele columns for allele-aware SNP identity; missing {missing}."
        )
    reasons = pd.Series(pd.NA, index=frame.index, dtype="object")
    values = pd.Series(pd.NA, index=frame.index, dtype="object")
    for idx, a1, a2 in zip(frame.index, frame["A1"], frame["A2"]):
        reason = _allele_failure_reason(a1, a2)
        if reason is None:
            values.loc[idx] = normalize_allele_set(a1, a2)
        else:
            reasons.loc[idx] = reason
    return values, reasons


def effective_merge_key_series(frame: pd.DataFrame, mode: str, *, context: str = "table") -> pd.Series:
    mode = normalize_snp_identifier_mode(mode)
    base = base_key_series(frame, mode, context=context)
    if not is_allele_aware_mode(mode):
        return base
    allele_set, reasons = allele_set_series(frame, context=context)
    if bool(reasons.notna().any()):
        bad_reason = str(reasons.dropna().iloc[0])
        raise ValueError(f"{context} contains allele rows invalid for {mode}: {bad_reason}.")
    return base.astype(str) + ":" + allele_set.astype(str)
```

- [ ] **Step 6: Implement drop-frame and cleanup helpers**

The cleanup helper must return a cleaned table plus a dropped table with `IDENTITY_DROP_COLUMNS`. Use `duplicate_identity` for all duplicate effective-key clusters, including base-mode duplicate clusters. Use `multi_allelic_base_key` only in allele-aware modes after allele sets are inspected.

Base-mode cleanup must not call `allele_set_series()` and must not inspect `A1/A2`. The implementation branch is:

```python
if is_allele_aware_mode(mode):
    # allele validation, allele-set construction, and multi-allelic cleanup live here
else:
    # duplicate_identity cleanup on base keys only
```

Implementation shape:

```python
def empty_identity_drop_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "CHR": pd.Series(dtype="string"),
            "SNP": pd.Series(dtype="string"),
            "source_pos": pd.Series(dtype="Int64"),
            "target_pos": pd.Series(dtype="Int64"),
            "reason": pd.Series(dtype="string"),
            "base_key": pd.Series(dtype="string"),
            "identity_key": pd.Series(dtype="string"),
            "allele_set": pd.Series(dtype="string"),
            "stage": pd.Series(dtype="string"),
        },
        columns=IDENTITY_DROP_COLUMNS,
    )


def coerce_identity_drop_frame(frame: pd.DataFrame | None) -> pd.DataFrame:
    if frame is None:
        return empty_identity_drop_frame()
    out = frame.copy()
    for column in IDENTITY_DROP_COLUMNS:
        if column not in out.columns:
            out[column] = pd.NA
    out = out.loc[:, IDENTITY_DROP_COLUMNS]
    for column in ("CHR", "SNP", "reason", "base_key", "identity_key", "allele_set", "stage"):
        out[column] = out[column].astype("string")
    out["source_pos"] = pd.to_numeric(out["source_pos"], errors="coerce").astype("Int64")
    out["target_pos"] = pd.to_numeric(out["target_pos"], errors="coerce").astype("Int64")
    return out.reset_index(drop=True)


def _identity_drop_rows(frame: pd.DataFrame, *, reason: str, stage: str) -> pd.DataFrame:
    out = empty_identity_drop_frame().reindex(range(len(frame))).copy()
    if "CHR" in frame.columns:
        out["CHR"] = frame["CHR"].astype("string").reset_index(drop=True)
    if "SNP" in frame.columns:
        out["SNP"] = frame["SNP"].astype("string").reset_index(drop=True)
    if "POS" in frame.columns:
        out["source_pos"] = pd.to_numeric(frame["POS"], errors="coerce").astype("Int64").reset_index(drop=True)
    if "_ldsc_base_key" in frame.columns:
        out["base_key"] = frame["_ldsc_base_key"].astype("string").reset_index(drop=True)
    if "_ldsc_identity_key" in frame.columns:
        out["identity_key"] = frame["_ldsc_identity_key"].astype("string").reset_index(drop=True)
    if "_ldsc_allele_set" in frame.columns:
        out["allele_set"] = frame["_ldsc_allele_set"].astype("string").reset_index(drop=True)
    out["reason"] = reason
    out["stage"] = stage
    return coerce_identity_drop_frame(out)
```

The artifact cleanup algorithm:

```python
def clean_identity_artifact_table(
    frame: pd.DataFrame,
    mode: str,
    *,
    context: str,
    stage: str,
    logger: logging.Logger | None = None,
) -> IdentityCleanupResult:
    mode = normalize_snp_identifier_mode(mode)
    work = frame.copy()
    dropped_frames: list[pd.DataFrame] = []
    work["_ldsc_base_key"] = base_key_series(work, mode, context=context)

    if is_allele_aware_mode(mode):
        work["_ldsc_allele_set"], reasons = allele_set_series(work, context=context)
        for reason in ("missing_allele", "invalid_allele", "strand_ambiguous_allele"):
            mask = reasons == reason
            if bool(mask.any()):
                dropped_frames.append(_identity_drop_rows(work.loc[mask], reason=reason, stage=stage))
        work = work.loc[reasons.isna()].copy()
        multi = work.groupby("_ldsc_base_key", dropna=False)["_ldsc_allele_set"].transform("nunique") > 1
        if bool(multi.any()):
            dropped_frames.append(_identity_drop_rows(work.loc[multi], reason="multi_allelic_base_key", stage=stage))
            work = work.loc[~multi].copy()
        work["_ldsc_identity_key"] = work["_ldsc_base_key"].astype(str) + ":" + work["_ldsc_allele_set"].astype(str)
    else:
        work["_ldsc_allele_set"] = pd.NA
        work["_ldsc_identity_key"] = work["_ldsc_base_key"]

    dup = work["_ldsc_identity_key"].duplicated(keep=False)
    if bool(dup.any()):
        dropped_frames.append(_identity_drop_rows(work.loc[dup], reason="duplicate_identity", stage=stage))
        work = work.loc[~dup].copy()

    drop_cols = ["_ldsc_base_key", "_ldsc_allele_set", "_ldsc_identity_key"]
    cleaned = work.drop(columns=[column for column in drop_cols if column in work.columns]).reset_index(drop=True)
    dropped = coerce_identity_drop_frame(pd.concat(dropped_frames, ignore_index=True) if dropped_frames else None)
    if logger is not None and len(dropped):
        counts = dropped["reason"].value_counts(sort=False).to_dict()
        logger.warning("%s identity cleanup dropped %d rows: %s", context, len(dropped), counts)
    return IdentityCleanupResult(cleaned=cleaned, dropped=dropped)
```

- [ ] **Step 7: Implement restriction key collapsing**

Behavior:

- If the restriction file lacks allele columns, return `match_kind="base"` and base keys.
- If the restriction file has allele columns and mode is allele-aware, drop invalid allele rows, drop multi-allelic base-key clusters, build effective keys, and collapse duplicate keys.
- If the restriction file has allele columns but active mode is base, ignore alleles and collapse base keys.

Function signature:

```python
def collapse_restriction_identity_keys(
    frame: pd.DataFrame,
    mode: str,
    *,
    context: str,
    has_allele_columns: bool | None = None,
    logger: logging.Logger | None = None,
) -> RestrictionIdentityKeys:
    mode = normalize_snp_identifier_mode(mode)
    has_alleles = {"A1", "A2"}.issubset(frame.columns) if has_allele_columns is None else has_allele_columns
    work = frame.copy()
    work["_ldsc_base_key"] = base_key_series(work, mode, context=context)

    if not has_alleles or not is_allele_aware_mode(mode):
        keys = set(work["_ldsc_base_key"].dropna().astype(str))
        return RestrictionIdentityKeys(
            keys=keys,
            match_kind="base",
            dropped=empty_identity_drop_frame(),
            n_input_rows=len(frame),
            n_retained_keys=len(keys),
        )

    dropped_frames: list[pd.DataFrame] = []
    work["_ldsc_allele_set"], reasons = allele_set_series(work, context=context)
    for reason in ("missing_allele", "invalid_allele", "strand_ambiguous_allele"):
        mask = reasons == reason
        if bool(mask.any()):
            dropped_frames.append(_identity_drop_rows(work.loc[mask], reason=reason, stage="restriction_cleanup"))
    work = work.loc[reasons.isna()].copy()

    multi = work.groupby("_ldsc_base_key", dropna=False)["_ldsc_allele_set"].transform("nunique") > 1
    if bool(multi.any()):
        dropped_frames.append(_identity_drop_rows(work.loc[multi], reason="multi_allelic_base_key", stage="restriction_cleanup"))
        work = work.loc[~multi].copy()

    work["_ldsc_identity_key"] = work["_ldsc_base_key"].astype(str) + ":" + work["_ldsc_allele_set"].astype(str)
    keys = set(work["_ldsc_identity_key"].dropna().astype(str))
    dropped = coerce_identity_drop_frame(pd.concat(dropped_frames, ignore_index=True) if dropped_frames else None)
    if logger is not None and len(dropped):
        logger.warning("%s restriction cleanup dropped %d rows before matching.", context, len(dropped))
    return RestrictionIdentityKeys(
        keys=keys,
        match_kind="identity",
        dropped=dropped,
        n_input_rows=len(frame),
        n_retained_keys=len(keys),
    )
```

- [ ] **Step 8: Re-export the exact mode normalizer from `column_inference.py`**

Replace the old alias-based function body with an import:

```python
from ._kernel.snp_identity import SNP_IDENTIFIER_MODES, normalize_snp_identifier_mode
```

Remove `_RSID_IDENTIFIER_ALIASES` and `_CHR_POS_IDENTIFIER_ALIASES` from mode normalization. Keep all column alias constants.

- [ ] **Step 9: Keep `identifiers.py` as parser, not policy owner**

Update `src/ldsc/_kernel/identifiers.py` so:

- `build_snp_id_series()` delegates to `effective_merge_key_series()` for canonical artifact-like tables.
- Restriction readers parse aliases into canonical frames and return `RestrictionIdentityKeys` when allele-aware matching is needed.
- Existing `read_global_snp_restriction()` remains available for base-mode callers during transition, but new workflow code uses the richer restriction object where allele precision matters.

- [ ] **Step 10: Run focused tests**

Run:

```bash
pytest -q tests/test_snp_identity.py tests/test_column_inference.py tests/test_config_identifiers.py
```

Expected: new identity tests pass; column/config tests fail only where old alias/default expectations still need Task 2 changes.

- [ ] **Step 11: Commit Task 1**

```bash
git add src/ldsc/_kernel/snp_identity.py src/ldsc/column_inference.py src/ldsc/_kernel/identifiers.py tests/test_snp_identity.py tests/test_column_inference.py tests/test_config_identifiers.py
git commit -m "feat: add shared SNP identity service"
```

---

## Task 2: Update Public Modes, Defaults, Config, And CLI Surfaces

**Files:**

- Modify: `src/ldsc/config.py`
- Modify: `src/ldsc/annotation_builder.py`
- Modify: `src/ldsc/ldscore_calculator.py`
- Modify: `src/ldsc/ref_panel_builder.py`
- Modify: `src/ldsc/sumstats_munger.py`
- Modify: `src/ldsc/_kernel/sumstats_munger.py`
- Modify: `src/ldsc/_kernel/ldscore.py`
- Modify: `src/ldsc/__init__.py`
- Test: `tests/test_config_identifiers.py`
- Test: `tests/test_column_inference.py`

- [ ] **Step 1: Update tests for exact mode values**

Change tests that currently expect alias normalization. Required expectations:

```python
def test_normalize_snp_identifier_mode_accepts_only_public_modes(self):
    for value in ("rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"):
        self.assertEqual(normalize_snp_identifier_mode(value), value)
    for value in ("rsID", "SNPID", "snp_id", "snp", "chrpos", "rsid_alleles", "chr_pos_alleles"):
        with self.assertRaises(ValueError):
            normalize_snp_identifier_mode(value)
```

Global config default expectation:

```python
config = GlobalConfig()
self.assertEqual(config.snp_identifier, "chr_pos_allele_aware")
self.assertEqual(config.genome_build, "auto")
```

- [ ] **Step 2: Expand `SNPIdentifierMode` and defaults**

In `src/ldsc/config.py`:

```python
SNPIdentifierMode = Literal["rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"]
```

Update `GlobalConfig` defaults:

```python
snp_identifier: SNPIdentifierMode = "chr_pos_allele_aware"
```

Update `__init__` default:

```python
snp_identifier: SNPIdentifierMode = "chr_pos_allele_aware"
```

Use identity helpers:

```python
mode = normalize_snp_identifier_mode(self.snp_identifier)
family = identity_mode_family(mode)
```

Validation rules:

- `chr_pos` and `chr_pos_allele_aware` require `genome_build` not `None`.
- `rsid` and `rsid_allele_aware` set `genome_build` to `None`.
- `genome_build="auto"` is invalid for both rsID-family modes.
- The warning for ignored genome build applies to both rsID-family modes.

- [ ] **Step 3: Add `RegressionConfig.allow_identity_downgrade`**

Add:

```python
allow_identity_downgrade: bool = False
```

Document:

```python
allow_identity_downgrade : bool, optional
    If True, same-family allele-aware/base regression inputs may run under the
    base identity mode. Cross-family mixes remain rejected. Default is False.
```

- [ ] **Step 4: Update CLI `--snp-identifier` choices and defaults**

Every public parser that exposes `--snp-identifier` must use:

```python
choices=("rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware")
default="chr_pos_allele_aware"
```

Files:

- `src/ldsc/_kernel/sumstats_munger.py`
- `src/ldsc/sumstats_munger.py`
- `src/ldsc/ref_panel_builder.py`
- `src/ldsc/annotation_builder.py`
- `src/ldsc/ldscore_calculator.py`
- `src/ldsc/_kernel/ldscore.py`

- [ ] **Step 5: Update mode-family conditions**

Replace direct checks:

```python
if mode == "chr_pos":
```

with:

```python
if identity_mode_family(mode) == "chr_pos":
```

Replace direct rsID checks:

```python
if mode == "rsid":
```

with:

```python
if identity_mode_family(mode) == "rsid":
```

Use exact mode checks only when behavior genuinely differs between base and allele-aware.

- [ ] **Step 6: Add regression CLI downgrade flag**

In `src/ldsc/regression_runner.py`, add `--allow-identity-downgrade` to h2, partitioned-h2, and rg parsers:

```python
parser.add_argument(
    "--allow-identity-downgrade",
    action="store_true",
    default=False,
    help=(
        "Allow same-family allele-aware/base identity mixes in regression by "
        "running with the base identity mode. rsID-family and chr_pos-family "
        "inputs still cannot be mixed."
    ),
)
```

Wire it into `_runner_from_args()`:

```python
allow_identity_downgrade=getattr(args, "allow_identity_downgrade", False),
```

- [ ] **Step 7: Run focused mode/config tests**

Run:

```bash
pytest -q tests/test_config_identifiers.py tests/test_column_inference.py
```

Expected: tests pass with exact modes and new default.

- [ ] **Step 8: Commit Task 2**

```bash
git add src/ldsc tests/test_config_identifiers.py tests/test_column_inference.py
git commit -m "feat: expose four SNP identifier modes"
```

---

## Task 3: Implement Minimal Identity Provenance And Artifact Validation

**Files:**

- Modify: `src/ldsc/_kernel/snp_identity.py`
- Modify: `src/ldsc/sumstats_munger.py`
- Modify: `src/ldsc/outputs.py`
- Modify: `src/ldsc/regression_runner.py`
- Modify: `src/ldsc/_kernel/ref_panel.py`
- Modify: `src/ldsc/_kernel/ref_panel_builder.py`
- Test: `tests/test_sumstats_munger.py`
- Test: `tests/test_output.py`
- Test: `tests/test_ref_panel.py`
- Test: `tests/test_regression_workflow.py`

- [ ] **Step 1: Add artifact metadata helpers**

In `snp_identity.py`, add:

```python
ARTIFACT_TYPES = frozenset({"sumstats", "ref_panel_r2", "ref_panel_metadata", "ldscore"})


def identity_artifact_metadata(*, artifact_type: str, snp_identifier: str, genome_build: str | None) -> dict[str, object]:
    if artifact_type not in ARTIFACT_TYPES:
        raise ValueError(f"Unknown LDSC artifact_type {artifact_type!r}.")
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_type": artifact_type,
        "snp_identifier": normalize_snp_identifier_mode(snp_identifier),
        "genome_build": genome_build,
    }


def validate_identity_artifact_metadata(metadata: dict[str, object], *, expected_artifact_type: str) -> str:
    if metadata.get("schema_version") != SCHEMA_VERSION or metadata.get("artifact_type") != expected_artifact_type:
        raise ValueError("This artifact was not written with the current LDSC schema/provenance contract. Regenerate it with the current LDSC package.")
    return normalize_snp_identifier_mode(str(metadata.get("snp_identifier")))
```

- [ ] **Step 2: Write sumstats metadata tests**

Add tests that assert:

- New `sumstats.metadata.json` has top-level `schema_version`, `artifact_type`, `snp_identifier`, and `genome_build`.
- `load_sumstats()` rejects a metadata sidecar missing `schema_version`.
- Allele-aware metadata without `A1/A2` in the artifact raises a malformed-artifact error.

- [ ] **Step 3: Update sumstats metadata writer and loader**

Change `_write_sumstats_metadata()` payload to:

```python
payload = {
    "schema_version": 1,
    "artifact_type": "sumstats",
    "snp_identifier": config_snapshot.snp_identifier,
    "genome_build": config_snapshot.genome_build,
    "trait_name": trait_name,
}
```

Do not persist `log_level` or `fail_on_missing_metadata`.

Change `_global_config_from_sumstats_metadata()` to reconstruct:

```python
mode = validate_identity_artifact_metadata(metadata, expected_artifact_type="sumstats")
return GlobalConfig(
    snp_identifier=mode,
    genome_build=metadata.get("genome_build"),
    log_level="INFO",
)
```

In `load_sumstats()`, after loading data and metadata:

```python
if is_allele_aware_mode(config_snapshot.snp_identifier) and not {"A1", "A2"}.issubset(df.columns):
    raise ValueError(
        "Munged sumstats artifact is malformed: snp_identifier="
        f"{config_snapshot.snp_identifier!r} requires A1/A2 columns. "
        "Regenerate it with the current LDSC package."
    )
```

- [ ] **Step 4: Update LD-score manifest writer**

In `outputs.py::LDScoreDirectoryWriter.build_manifest()`, add top-level:

```python
"schema_version": 1,
"artifact_type": "ldscore",
```

Keep existing workflow-specific fields: `files`, `chromosomes`, `baseline_columns`, `query_columns`, `counts`, `count_config`, row groups, and row counts.

Remove `config_snapshot` from newly written manifest. Keep reading old manifests rejected by `schema_version` check in Task 9.

- [ ] **Step 5: Add R2 parquet schema metadata**

In `_kernel/ref_panel_builder.py::write_r2_parquet()`, add Arrow schema metadata:

```python
b"ldsc:schema_version": b"1",
b"ldsc:artifact_type": b"ref_panel_r2",
b"ldsc:snp_identifier": snp_identifier.encode("utf-8"),
b"ldsc:genome_build": genome_build.encode("utf-8"),
```

Update the function signature:

```python
def write_r2_parquet(
    *,
    pair_rows: Iterable[dict[str, float | int | str]],
    reference_snp_table: pd.DataFrame,
    path: str | PathLike[str],
    genome_build: str,
    snp_identifier: str,
    n_samples: int,
    batch_size: int = 100_000,
    row_group_size: int = 50_000,
) -> str:
```

Pass `self.global_config.snp_identifier` from `ref_panel_builder.py`.

- [ ] **Step 6: Validate R2 parquet schema metadata on load**

In `_kernel/ref_panel.py`, add reader helper:

```python
def _read_identity_schema_meta(path: str, *, expected_artifact_type: str) -> dict[str, object]:
    import pyarrow.parquet as pq
    raw = pq.read_schema(path).metadata or {}
    if b"ldsc:schema_version" not in raw or b"ldsc:artifact_type" not in raw:
        raise ValueError("This artifact was not written with the current LDSC schema/provenance contract. Regenerate it with the current LDSC package.")
    return {
        "schema_version": int(raw[b"ldsc:schema_version"].decode("utf-8")),
        "artifact_type": raw[b"ldsc:artifact_type"].decode("utf-8"),
        "snp_identifier": raw[b"ldsc:snp_identifier"].decode("utf-8"),
        "genome_build": raw[b"ldsc:genome_build"].decode("utf-8"),
    }
```

For external raw R2 parquet, do not call this package-artifact validation in base modes; call it for package-built canonical R2. In allele-aware modes, raw layout is rejected in Task 7.

- [ ] **Step 7: Run artifact metadata tests**

Run:

```bash
pytest -q tests/test_sumstats_munger.py tests/test_output.py tests/test_ref_panel.py tests/test_regression_workflow.py
```

- [ ] **Step 8: Commit Task 3**

```bash
git add src/ldsc tests/test_sumstats_munger.py tests/test_output.py tests/test_ref_panel.py tests/test_regression_workflow.py
git commit -m "feat: validate identity artifact provenance"
```

---

## Task 4: Refactor Restriction File Semantics

**Files:**

- Modify: `src/ldsc/column_inference.py`
- Modify: `src/ldsc/_kernel/identifiers.py`
- Modify: `src/ldsc/_kernel/ref_panel_builder.py`
- Modify: `src/ldsc/_kernel/ldscore.py`
- Modify: `src/ldsc/ldscore_calculator.py`
- Modify: `src/ldsc/sumstats_munger.py`
- Test: `tests/test_config_identifiers.py`
- Test: `tests/test_sumstats_munger.py`
- Test: `tests/test_ref_panel_builder.py`
- Test: `tests/test_ldscore_workflow.py`

- [ ] **Step 1: Add allele column specs for external restriction files**

In `column_inference.py`, add reusable allele specs:

```python
A1_COLUMN_SPEC = ColumnSpec(
    "A1",
    ("A1", "ALLELE1", "ALLELE_1", "EFFECT_ALLELE", "REFERENCE_ALLELE", "REF", "EA"),
    "allele 1",
)
A2_COLUMN_SPEC = ColumnSpec(
    "A2",
    ("A2", "ALLELE2", "ALLELE_2", "OTHER_ALLELE", "NON_EFFECT_ALLELE", "ALT", "NEA"),
    "allele 2",
)
RESTRICTION_ALLELE_SPECS = (A1_COLUMN_SPEC, A2_COLUMN_SPEC)
```

Do not treat `REF/ALT` as package output columns. They are external input aliases only.

- [ ] **Step 2: Add tests for restriction precision**

Required cases:

- No allele columns in allele-aware mode returns `match_kind="base"`.
- Allele columns in allele-aware mode returns `match_kind="identity"`.
- Allele-bearing restriction rows with missing alleles are dropped and logged at `WARNING`.
- Duplicate restriction keys collapse to one key.
- Multi-allelic restriction clusters are dropped in allele-aware modes.
- Base modes ignore allele columns and use `duplicate_identity` only during later artifact cleanup.
- Base-mode restriction files with invalid, missing, strand-ambiguous, or multi-base allele values still produce base-key restriction sets; allele cleanup functions are not called.

- [ ] **Step 3: Refactor `read_global_snp_restriction()` or add a richer reader**

Add:

```python
def read_snp_restriction_keys(
    path: str | Path,
    snp_identifier: str,
    genome_build: str | None = None,
    logger=None,
) -> RestrictionIdentityKeys:
    raw_frame = _read_snp_restriction_table(path)
    canonical_frame, has_allele_columns = _canonicalize_snp_restriction_table(
        raw_frame,
        snp_identifier=snp_identifier,
        genome_build=genome_build,
    )
    return collapse_restriction_identity_keys(
        canonical_frame,
        snp_identifier,
        context=f"SNP restriction file {path}",
        has_allele_columns=has_allele_columns,
        logger=logger,
    )
```

Implementation:

- Parse the file into canonical columns `SNP`, `CHR`, `POS`, optional `A1`, `A2`.
- Use the active mode family to require `SNP` or `CHR/POS`.
- Detect allele columns only when both `A1` and `A2` resolve.
- Pass canonical frame to `collapse_restriction_identity_keys()`.
- Keep the packed CHR/POS streaming helper for munger base-mode chunk filtering; add an allele-aware fallback that uses string keys because allele-aware filtering requires `A1/A2`.

- [ ] **Step 4: Update restriction matching masks**

Add a shared matcher:

```python
def restriction_membership_mask(
    frame: pd.DataFrame,
    restriction: RestrictionIdentityKeys,
    mode: str,
    *,
    context: str,
) -> pd.Series:
    if restriction.match_kind == "base":
        keys = base_key_series(frame, mode, context=context)
    else:
        keys = effective_merge_key_series(frame, mode, context=context)
    return keys.isin(restriction.keys)
```

Use this matcher in:

- `ref_panel_builder.py` for `--ref-panel-snps-file` and HM3 restrictions.
- `ldscore_calculator.py` for regression SNP restrictions.
- `_kernel/ldscore.py` for kernel path regression key masks.
- `sumstats_munger.py` for non-packed allele-aware keep-list filtering.

- [ ] **Step 5: Preserve fast packed filtering for `chr_pos` base mode**

Keep `read_global_chr_pos_restriction_key_set()` and `build_packed_chr_pos_series()` for `chr_pos` base mode in munger chunk filtering. Use the richer restriction reader for:

- `chr_pos_allele_aware`
- `rsid_allele_aware`
- external restriction files with allele columns when active mode is allele-aware

- [ ] **Step 6: Run restriction-focused tests**

Run:

```bash
pytest -q tests/test_config_identifiers.py tests/test_sumstats_munger.py tests/test_ref_panel_builder.py tests/test_ldscore_workflow.py
```

- [ ] **Step 7: Commit Task 4**

```bash
git add src/ldsc tests/test_config_identifiers.py tests/test_sumstats_munger.py tests/test_ref_panel_builder.py tests/test_ldscore_workflow.py
git commit -m "feat: support allele-aware SNP restrictions"
```

---

## Task 5: Refactor Sumstats Munging And Loading

**Files:**

- Modify: `src/ldsc/_kernel/sumstats_munger.py`
- Modify: `src/ldsc/sumstats_munger.py`
- Modify: `src/ldsc/_kernel/liftover.py`
- Test: `tests/test_sumstats_munger.py`
- Test: `tests/test_sumstats_liftover.py`

- [ ] **Step 1: Write failing sumstats tests**

Add tests for:

- Default `chr_pos_allele_aware` rejects raw input without `A1/A2` and suggests `--snp-identifier chr_pos`.
- `chr_pos` runs without alleles when `--no-alleles` is passed.
- `rsid` runs without alleles when `--no-alleles` is passed.
- `chr_pos` and `rsid` keep singleton SNP rows even if allele columns are present but invalid, missing, strand-ambiguous, or multi-base.
- `rsid_allele_aware` requires alleles.
- Global identity cleanup after liftover drops duplicate effective-key clusters.
- Liftover source/target duplicate detection ignores allele set.
- Dropped sidecar has columns `CHR`, `SNP`, `source_pos`, `target_pos`, `reason`, `base_key`, `identity_key`, `allele_set`, `stage`.
- Base-mode duplicate drops use `duplicate_identity`.

- [ ] **Step 2: Enforce allele-aware input requirements**

In `_kernel/sumstats_munger.py`, replace:

```python
if not args.no_alleles and not all(x in cname_translation.values() for x in ['A1', 'A2']):
```

with mode-aware logic:

```python
mode = normalize_snp_identifier_mode(getattr(args, "snp_identifier", "chr_pos_allele_aware"))
requires_alleles = is_allele_aware_mode(mode)
if requires_alleles and not all(x in cname_translation.values() for x in ["A1", "A2"]):
    raise ValueError(
        f"This run is using snp_identifier={mode!r}, which requires A1/A2 allele columns. "
        f"No usable allele columns were found. To run without allele-aware SNP identity, rerun with --snp-identifier {identity_base_mode(mode)}."
    )
if args.no_alleles and requires_alleles:
    raise ValueError(f"--no-alleles cannot be used with snp_identifier={mode!r}. Rerun with --snp-identifier {identity_base_mode(mode)}.")
```

Keep allele-free support for `rsid` and `chr_pos`.

- [ ] **Step 3: Remove old rsID keep-first duplicate policy**

Replace:

```python
dat = dat.drop_duplicates(subset='SNP')
```

with the shared cleanup after liftover and after `P` to `Z` conversion. Duplicate identity must be drop-all, not keep-first.

- [ ] **Step 4: Run shared identity cleanup after liftover**

In `_apply_liftover_if_requested()` or immediately after it, call:

```python
cleanup = clean_identity_artifact_table(
    dat,
    normalize_snp_identifier_mode(getattr(args, "snp_identifier", "chr_pos_allele_aware")),
    context="munged sumstats",
    stage="post_liftover_identity_cleanup",
    logger=LOGGER,
)
dat = cleanup.cleaned
```

Append `cleanup.dropped` to the existing liftover drop frame before writing the sidecar.

- [ ] **Step 5: Extend sumstats dropped-SNP sidecar schema**

Replace the current `LIFTOVER_DROP_COLUMNS` sidecar schema in `sumstats_munger.py` with `IDENTITY_DROP_COLUMNS`. Existing liftover frames must be coerced with empty `base_key`, `identity_key`, `allele_set`, and `stage="liftover"`.

Implementation:

```python
drop_frame = coerce_identity_drop_frame(coordinate_metadata.pop("liftover_drop_frame", None))
identity_drop_frame = coerce_identity_drop_frame(getattr(args, "_identity_drop_frame", None))
drop_frame = coerce_identity_drop_frame(pd.concat([drop_frame, identity_drop_frame], ignore_index=True))
```

- [ ] **Step 6: Update `SumstatsTable` identity helpers**

Change `SumstatsTable.snp_identifiers()` and `align_to_metadata()` to use:

```python
effective_merge_key_series(self.data, mode, context="sumstats table")
```

For regression downgrade, Task 9 will pass an effective base mode.

- [ ] **Step 7: Validate loaded sumstats artifacts**

In `load_sumstats()`:

- Reject missing metadata sidecar as old package-written artifact.
- Validate `schema_version`, `artifact_type`, `snp_identifier`, and `genome_build`.
- Reject allele-aware modes missing `A1/A2`.
- Run duplicate effective-key drop-all on loaded artifacts or raise malformed. Prefer raising for package-written artifacts because they should already be clean.

Error text:

```text
This artifact was not written with the current LDSC schema/provenance contract. Regenerate it with the current LDSC package.
```

- [ ] **Step 8: Run sumstats tests**

Run:

```bash
pytest -q tests/test_sumstats_munger.py tests/test_sumstats_liftover.py
```

- [ ] **Step 9: Commit Task 5**

```bash
git add src/ldsc/_kernel/sumstats_munger.py src/ldsc/sumstats_munger.py src/ldsc/_kernel/liftover.py tests/test_sumstats_munger.py tests/test_sumstats_liftover.py
git commit -m "feat: apply identity cleanup to munged sumstats"
```

---

## Task 6: Refactor Reference-Panel Builder And R2 Output Contract

**Files:**

- Modify: `src/ldsc/ref_panel_builder.py`
- Modify: `src/ldsc/_kernel/ref_panel_builder.py`
- Modify: `src/ldsc/config.py`
- Test: `tests/test_ref_panel_builder.py`
- Test: `tests/test_plink_io.py`

- [ ] **Step 1: Write failing ref-panel builder tests**

Required cases:

- `rsid` and `rsid_allele_aware` reject matching chain-file liftover.
- `rsid` and `rsid_allele_aware` reject HM3 quick liftover.
- `chr_pos_allele_aware` requires BIM allele columns and drops invalid allele rows.
- Multi-allelic base-key clusters are dropped package-wide in allele-aware modes.
- Base modes drop duplicate effective keys with reason `duplicate_identity`.
- Cross-build synchronized filtering drops a variant from both emitted builds if it collides in either source or target build.
- R2 parquet includes `A1_1`, `A2_1`, `A1_2`, `A2_2`.
- Metadata sidecars include `A1`, `A2`.
- Per-chromosome dropped sidecar has the extended identity schema.

- [ ] **Step 2: Include alleles in runtime metadata**

In `_kernel/ref_panel_builder.py::build_runtime_metadata_table()`, add:

```python
"A1": metadata["A1"].astype(str),
"A2": metadata["A2"].astype(str),
```

Place `A1/A2` after `SNP` and before `CM`.

- [ ] **Step 3: Include endpoint alleles in R2 parquet**

Update `_STANDARD_R2_COLUMNS`:

```python
_STANDARD_R2_COLUMNS = [
    "CHR",
    "POS_1",
    "POS_2",
    "SNP_1",
    "SNP_2",
    "A1_1",
    "A2_1",
    "A1_2",
    "A2_2",
    "R2",
]
```

Update `_empty_standard_r2_table()` and `build_standard_r2_table()` to populate endpoint allele columns from `left["REF"]`, `left["ALT"]`, `right["REF"]`, and `right["ALT"]`, or rename `REF/ALT` to `A1/A2` inside `build_reference_snp_table()` to avoid genome-reference terminology in package outputs.

Preferred implementation:

```python
reference_snp_table = pd.DataFrame(
    {
        "chr": chromosomes.astype(str),
        "hg19_pos": pd.to_numeric(metadata["hg19_pos"], errors="raise").astype("int64"),
        "hg38_pos": pd.to_numeric(metadata["hg38_pos"], errors="raise").astype("int64"),
        "rsID": metadata["SNP"].astype(str),
        "MAF": pd.to_numeric(metadata["MAF"], errors="raise").astype(float),
        "A1": metadata["A1"].astype(str),
        "A2": metadata["A2"].astype(str),
    }
)
```

Then:

```python
"A1_1": left["A1"].astype(str),
"A2_1": left["A2"].astype(str),
"A1_2": right["A1"].astype(str),
"A2_2": right["A2"].astype(str),
```

- [ ] **Step 4: Apply restriction then shared identity cleanup**

In `ReferencePanelBuilder._build_chromosome()`:

1. Load BIM rows and normalize `CHR/BP`.
2. Apply `ref_panel_snps_file` or HM3 restriction.
3. Build metadata with `CHR`, `SNP`, `POS`, `A1`, `A2`.
4. Call `clean_identity_artifact_table(metadata, self.global_config.snp_identifier, context=f"reference panel chromosome {chrom}", stage="ref_panel_source_identity_cleanup", logger=LOGGER)`.
5. Convert surviving rows back to PLINK row indices for genotype loading.

Do not run allele-aware multi-allelic retention. V1 drops every multi-allelic base-key cluster.

- [ ] **Step 5: Keep liftover duplicate policy coordinate-only**

Existing `_resolve_unique_snp_set()` checks source and target `CHR:POS`. Keep that behavior for liftover only, but apply it for both coordinate-family modes:

```python
if identity_mode_family(self.global_config.snp_identifier) == "chr_pos":
    unified_drop_frame = _resolve_unique_snp_set(source_frame, target_frame, chrom=chrom)
```

When liftover drops rows, use reasons `source_duplicate` and `target_collision`, not identity reasons.

- [ ] **Step 6: Synchronize cross-build drops**

After source/target liftover collision detection, compute the union of dropped PLINK row indices from either emitted build and remove those indices from all emitted build outputs. Log:

```text
Reference-panel liftover synchronized cross-build drops: dropped N SNPs from all emitted builds because they collided in at least one emitted build.
```

- [ ] **Step 7: Reject rsID-family liftover**

Update checks:

```python
if matching_chain is not None and identity_mode_family(self.global_config.snp_identifier) == "rsid":
    raise ValueError("Reference-panel chain liftover is only valid for chr_pos-family modes; omit the matching liftover chain in rsID-family modes.")
if config.use_hm3_quick_liftover and identity_mode_family(self.global_config.snp_identifier) == "rsid":
    raise ValueError("Reference-panel HM3 quick liftover is only valid for chr_pos-family modes.")
```

- [ ] **Step 8: Extend ref-panel dropped sidecar schema**

Replace local five-column `_empty_unified_drop_frame()` and `_coerce_unified_drop_frame()` with `empty_identity_drop_frame()` and `coerce_identity_drop_frame()`.

For liftover rows, fill:

- `reason`: liftover reason
- `stage`: `ref_panel_liftover`
- `base_key`, `identity_key`, `allele_set`: empty unless the row was dropped by identity cleanup

For identity cleanup rows, fill all identity fields and `stage`.

- [ ] **Step 9: Run ref-panel tests**

Run:

```bash
pytest -q tests/test_ref_panel_builder.py tests/test_plink_io.py
```

- [ ] **Step 10: Commit Task 6**

```bash
git add src/ldsc/ref_panel_builder.py src/ldsc/_kernel/ref_panel_builder.py src/ldsc/config.py tests/test_ref_panel_builder.py tests/test_plink_io.py
git commit -m "feat: emit allele-aware reference-panel artifacts"
```

---

## Task 7: Refactor Reference-Panel Loading And R2 Readers

**Files:**

- Modify: `src/ldsc/_kernel/ref_panel.py`
- Modify: `src/ldsc/_kernel/ldscore.py`
- Modify: `src/ldsc/column_inference.py`
- Test: `tests/test_ref_panel.py`
- Test: `tests/test_plink_io.py`
- Test: `tests/test_ldscore_workflow.py`

- [ ] **Step 1: Add parquet endpoint allele specs**

In `column_inference.py`, extend `PARQUET_R2_CANONICAL_SPECS` with:

```python
ColumnSpec("A1_1", ("A1_1", "A11", "ALLELE1_1", "ALLELE_1_1"), "canonical parquet R2 left A1", allow_suffix_match=False),
ColumnSpec("A2_1", ("A2_1", "A21", "ALLELE2_1", "ALLELE_2_1"), "canonical parquet R2 left A2", allow_suffix_match=False),
ColumnSpec("A1_2", ("A1_2", "A12", "ALLELE1_2", "ALLELE_1_2"), "canonical parquet R2 right A1", allow_suffix_match=False),
ColumnSpec("A2_2", ("A2_2", "A22", "ALLELE2_2", "ALLELE_2_2"), "canonical parquet R2 right A2", allow_suffix_match=False),
```

Use these columns as required only when active mode is allele-aware.

- [ ] **Step 2: Write R2 loading tests**

Required cases:

- Canonical package-built R2 with endpoint allele columns loads in all four modes.
- Canonical R2 missing endpoint allele columns loads only in `rsid` and `chr_pos`.
- External raw R2 parquet loads only in `rsid` and `chr_pos`.
- External raw R2 parquet in `rsid_allele_aware` or `chr_pos_allele_aware` raises:

```text
allele-aware SNP identity requires package-built canonical R2 parquet with A1_1/A2_1/A1_2/A2_2 endpoint allele columns
```

- [ ] **Step 3: Validate metadata sidecars**

In `_kernel/ref_panel.py::_read_metadata_table()` or the nearest loader helper:

- Require `A1/A2` for allele-aware modes.
- Preserve `A1/A2` in returned metadata.
- Reject metadata sidecar missing `A1/A2` under allele-aware mode as malformed.

- [ ] **Step 4: Update `SortedR2BlockReader` index maps**

For allele-aware modes:

- Metadata index map key is `effective_merge_key_series(metadata, mode)`.
- Canonical parquet endpoint key for each endpoint is:
  - `SNP_1:<allele_set_1>` for `rsid_allele_aware`
  - `CHR:POS_1:<allele_set_1>` for `chr_pos_allele_aware`

Implementation helper:

```python
def _endpoint_identity_keys(table: pd.DataFrame, *, side: int, mode: str) -> pd.Series:
    if side == 1:
        frame = table.rename(columns={"SNP_1": "SNP", "POS_1": "POS", "A1_1": "A1", "A2_1": "A2"})
    else:
        frame = table.rename(columns={"SNP_2": "SNP", "POS_2": "POS", "A1_2": "A1", "A2_2": "A2"})
    return effective_merge_key_series(frame, mode, context=f"R2 endpoint {side}")
```

For base modes, keep current `SNP` or `POS` index maps and do not inspect alleles.

- [ ] **Step 5: Reject raw layout in allele-aware modes**

In `SortedR2BlockReader.__init__()` after layout detection:

```python
if layout == "raw" and is_allele_aware_mode(self.identifier_mode):
    raise ValueError(
        "allele-aware SNP identity requires package-built canonical R2 parquet with "
        "A1_1/A2_1/A1_2/A2_2 endpoint allele columns; external raw R2 parquet is supported only for rsid and chr_pos."
    )
```

- [ ] **Step 6: Run R2 reader tests**

Run:

```bash
pytest -q tests/test_ref_panel.py tests/test_plink_io.py tests/test_ldscore_workflow.py
```

- [ ] **Step 7: Commit Task 7**

```bash
git add src/ldsc/_kernel/ref_panel.py src/ldsc/_kernel/ldscore.py src/ldsc/column_inference.py tests/test_ref_panel.py tests/test_plink_io.py tests/test_ldscore_workflow.py
git commit -m "feat: load allele-aware R2 endpoints"
```

---

## Task 8: Refactor Annotation And LD-Score Workflows

**Files:**

- Modify: `src/ldsc/annotation_builder.py`
- Modify: `src/ldsc/_kernel/annotation.py`
- Modify: `src/ldsc/ldscore_calculator.py`
- Modify: `src/ldsc/_kernel/ldscore.py`
- Modify: `src/ldsc/outputs.py`
- Test: `tests/test_annotation.py`
- Test: `tests/test_ldscore_workflow.py`
- Test: `tests/test_output.py`

- [ ] **Step 1: Write annotation/LD-score tests**

Required cases:

- Annotation files without `A1/A2` are accepted in `chr_pos_allele_aware`.
- Annotation files with `A1/A2` participate in allele-aware matching.
- Reference-panel cleanup drops multi-allelic clusters before annotation matching.
- LD-score outputs in allele-aware modes include `A1/A2`.
- LD-score baseline/query row alignment includes `A1/A2` when required.
- LD-score manifest includes minimal identity provenance.

- [ ] **Step 2: Parse optional annotation alleles**

In `annotation_builder.py::parse_annotation_file()` and `_kernel/ldscore.py::parse_annotation_file()`:

- Resolve optional `A1/A2` with the external allele specs.
- Add them to metadata if both resolve.
- If only one resolves, raise:

```text
Annotation file has only one allele column; provide both A1 and A2 or neither.
```

Do not require alleles for annotation inputs.

- [ ] **Step 3: Validate annotation rows by available precision**

In `AnnotationBundle.validate()`:

- If active mode is allele-aware and metadata contains both `A1/A2`, validate effective allele-aware keys.
- If active mode is allele-aware and metadata lacks `A1/A2`, validate base keys.
- If active mode is base, validate base keys and ignore allele columns.

Do not call `clean_identity_artifact_table()` on annotation bundles; annotations are external inputs and do not write reloadable identity artifacts.

- [ ] **Step 4: Align annotation bundles to reference panel by available precision**

In `_align_annotation_bundle_to_ref_panel()`:

```python
annotation_has_alleles = {"A1", "A2"}.issubset(annotation_bundle.metadata.columns)
mode_for_annotation_match = global_config.snp_identifier if annotation_has_alleles else identity_base_mode(global_config.snp_identifier)
```

Reference metadata is already cleaned and contains alleles in allele-aware modes. If annotation lacks alleles, match on base key.

- [ ] **Step 5: Preserve `A1/A2` in LD-score result tables**

In `LDScoreCalculator._wrap_legacy_chrom_result()`:

```python
metadata_cols = ["CHR", "SNP", pos_column]
for allele_col in ("A1", "A2"):
    if allele_col in reference_metadata.columns:
        metadata_cols.append(allele_col)
```

Then rename `pos_column` to `POS`.

In `_split_ldscore_table()`, metadata columns become:

```python
metadata_columns = ["CHR", "SNP", "POS", *[col for col in ("A1", "A2") if col in ldscore_table.columns]]
```

For allele-aware modes, require `A1/A2`.

- [ ] **Step 6: Update `assert_same_snp_rows()`**

In `_row_alignment.py`, use `effective_merge_key_series()` for row identity. For allele-aware modes, require `A1/A2` in both frames. For base modes, do not compare allele columns.

- [ ] **Step 7: Update output validation**

In `outputs.py::_validate_tables()`:

- Required LD-score columns include `A1/A2` when `is_allele_aware_mode(result.config_snapshot.snp_identifier)` is true.
- Manifest writes `schema_version`, `artifact_type`, `snp_identifier`, `genome_build`.

- [ ] **Step 8: Run LD-score tests**

Run:

```bash
pytest -q tests/test_annotation.py tests/test_ldscore_workflow.py tests/test_output.py
```

- [ ] **Step 9: Commit Task 8**

```bash
git add src/ldsc/annotation_builder.py src/ldsc/_kernel/annotation.py src/ldsc/ldscore_calculator.py src/ldsc/_kernel/ldscore.py src/ldsc/outputs.py src/ldsc/_row_alignment.py tests/test_annotation.py tests/test_ldscore_workflow.py tests/test_output.py
git commit -m "feat: propagate allele-aware LD-score identity"
```

---

## Task 9: Refactor Regression Compatibility, Downgrade, And Orientation

**Files:**

- Modify: `src/ldsc/regression_runner.py`
- Modify: `src/ldsc/config.py`
- Modify: `src/ldsc/_kernel/regression.py`
- Test: `tests/test_regression_workflow.py`

- [ ] **Step 1: Write regression compatibility tests**

Required cases:

- Exact same modes merge by effective identity.
- Cross-family modes reject even with `allow_identity_downgrade=True`.
- `chr_pos_allele_aware` sumstats with `chr_pos` LD-score rejects by default.
- The same pair runs with `allow_identity_downgrade=True` under effective `chr_pos`.
- `rsid_allele_aware` plus `rsid` runs with flag under effective `rsid`.
- Downgrade runs duplicate effective-key drop-all on each input under the base mode before merge.
- Downgrade log names original modes, effective mode, and rows dropped.
- h2 allele-aware runs orient sumstats `Z` to LD-score `A1/A2`.
- rg keeps current ordered allele harmonization after identity merge.

- [ ] **Step 2: Add identity compatibility helper**

In `snp_identity.py`:

```python
@dataclass(frozen=True)
class IdentityCompatibility:
    effective_mode: str
    downgrade_applied: bool


def resolve_regression_identity_mode(
    left_mode: str,
    right_mode: str,
    *,
    allow_identity_downgrade: bool,
) -> IdentityCompatibility:
    left = normalize_snp_identifier_mode(left_mode)
    right = normalize_snp_identifier_mode(right_mode)
    if left == right:
        return IdentityCompatibility(effective_mode=left, downgrade_applied=False)
    if identity_mode_family(left) != identity_mode_family(right):
        raise ValueError(f"Cannot mix SNP identifier families in regression: {left!r} vs {right!r}.")
    if not allow_identity_downgrade:
        raise ValueError(
            f"snp_identifier mismatch: {left!r} vs {right!r}. "
            "Use --allow-identity-downgrade to run same-family allele-aware/base inputs under the base mode."
        )
    return IdentityCompatibility(effective_mode=identity_base_mode(left), downgrade_applied=True)
```

- [ ] **Step 3: Add pre-merge downgrade filtering**

In `regression_runner.py`, add:

```python
def _prepare_regression_identity_table(frame: pd.DataFrame, mode: str, *, context: str, logger: logging.Logger) -> tuple[pd.DataFrame, pd.DataFrame]:
    cleanup = clean_identity_artifact_table(frame, mode, context=context, stage="regression_identity_downgrade", logger=logger)
    return cleanup.cleaned, cleanup.dropped
```

When downgrade applies, call it on each input with the effective base mode before merging.

Do not run allele-aware validation during base-mode downgrade.

- [ ] **Step 4: Replace `_effective_snp_identifier_mode()`**

Use the new compatibility helper and `RegressionConfig.allow_identity_downgrade`.

For h2:

- Compare sumstats mode and LD-score mode.
- Store effective mode on `RegressionDataset.config_snapshot` or add `effective_snp_identifier` field to `RegressionDataset`.

For rg:

- Compare trait 1, trait 2, and LD-score modes pairwise.
- Use one effective mode for all three inputs.

- [ ] **Step 5: Log downgrade**

When downgrade applies:

```python
LOGGER.warning(
    "Identity downgrade enabled: LD-score mode %s, sumstats mode %s; running regression with effective snp_identifier=%r. Dropped %d duplicate effective-key rows before merge.",
    ldscore_mode,
    sumstats_mode,
    effective_mode,
    dropped_count,
)
```

For rg, include both trait modes:

```text
Identity downgrade enabled: LD-score mode chr_pos, trait 1 mode chr_pos_allele_aware, trait 2 mode chr_pos_allele_aware; running rg with effective snp_identifier='chr_pos'. Dropped N duplicate effective-key rows before merge.
```

- [ ] **Step 6: Orient h2 sumstats to LD-score alleles**

After identity merge in `build_dataset()` and before regression:

```python
if is_allele_aware_mode(identifier_mode) and {"A1", "A2", "A1_ld", "A2_ld"}.issubset(merged.columns):
    merged["Z"] = _orient_sumstats_z_to_reference_alleles(merged)
```

Implementation:

```python
def _orient_sumstats_z_to_reference_alleles(merged: pd.DataFrame) -> pd.Series:
    same = (merged["A1"] == merged["A1_ld"]) & (merged["A2"] == merged["A2_ld"])
    flipped = (merged["A1"] == merged["A2_ld"]) & (merged["A2"] == merged["A1_ld"])
    incompatible = ~(same | flipped)
    if bool(incompatible.any()):
        LOGGER.warning("Dropping %d SNPs with incompatible sumstats and LD-score allele order.", int(incompatible.sum()))
    z = merged.loc[~incompatible, "Z"].copy()
    z.loc[flipped.loc[~incompatible]] *= -1
    return z
```

When implementing, preserve row filtering and index reset so every regression column remains aligned.

- [ ] **Step 7: Keep rg ordered harmonization after merge**

Do not replace the existing `reg._filter_alleles()` and `reg._align_alleles()` flow. Update it only enough to handle allele-aware effective merge keys and suffix column names cleanly.

The rg rule remains:

- Merge by effective identity.
- Then compare ordered trait alleles.
- Flip trait 2 `Z` when allele order is reversed.
- Drop incompatible rows.

- [ ] **Step 8: Update output metadata for regression detail writers**

When writing partitioned-h2 or rg detail metadata, include only:

```python
"effective_snp_identifier": effective_mode,
"identity_downgrade_applied": downgrade_applied,
```

Put original modes and dropped counts in the log, not JSON metadata.

- [ ] **Step 9: Run regression tests**

Run:

```bash
pytest -q tests/test_regression_workflow.py
```

- [ ] **Step 10: Commit Task 9**

```bash
git add src/ldsc/regression_runner.py src/ldsc/config.py src/ldsc/_kernel/regression.py tests/test_regression_workflow.py
git commit -m "feat: add regression identity downgrade"
```

---

## Task 10: Update Documentation, Docstrings, And User-Facing Errors

**Files:**

- Modify: `docs/current/architecture.md`
- Modify: `docs/current/data_flow.md`
- Modify: `docs/current/file_formats.md` or nearest current file-format doc
- Modify: `docs/current/inference-genome-build.md`
- Modify: `tutorials/munge_sumstats.md`
- Modify: public module docstrings touched by Tasks 1-9
- Test: `tests/test_logging_refactor.py`

- [ ] **Step 1: Document four modes and defaults**

Add a short section to current docs:

```text
LDSC supports exactly four SNP identifier modes: rsid, rsid_allele_aware,
chr_pos, and chr_pos_allele_aware. The default is chr_pos_allele_aware.
Mode names are exact; column aliases apply to input headers only.
```

- [ ] **Step 2: Document base-mode allele-free support**

Add:

```text
The base modes rsid and chr_pos do not inspect alleles. If allele columns are
present, they are preserved when useful but do not affect identity, duplicate
filtering, or merge keys.
```

- [ ] **Step 3: Document allele-aware drop policy**

Add:

```text
Allele-aware modes use alleles only to make merge keys safer. They do not
retain multi-allelic same-base-key clusters in v1. Rows with the same rsID or
CHR:POS but more than one allele set are dropped before LD-score or regression
artifacts are written.
```

- [ ] **Step 4: Document restrictions and annotations**

Add:

```text
Restriction files may omit alleles. Allele-free restriction files match by the
base key and can retain multiple candidate rows before later artifact cleanup.
Annotation files may omit alleles even in allele-aware modes because annotations
describe genomic membership, not variant alleles.
```

- [ ] **Step 5: Document external raw R2 limitation**

Add:

```text
External raw R2 parquet inputs are supported only in rsid and chr_pos modes.
Allele-aware modes require package-built canonical R2 parquet with endpoint
allele columns A1_1/A2_1/A1_2/A2_2.
```

- [ ] **Step 6: Document regression downgrade flag**

Add:

```text
--allow-identity-downgrade is regression-only. It allows same-family
allele-aware/base mixes to run under the base mode. rsID-family and
coordinate-family modes never mix.
```

- [ ] **Step 7: Update error message tests**

Ensure tests assert these repair messages:

```text
This run is using snp_identifier='chr_pos_allele_aware', which requires A1/A2 allele columns.
To run without allele-aware SNP identity, rerun with --snp-identifier chr_pos.
```

```text
This artifact was not written with the current LDSC schema/provenance contract. Regenerate it with the current LDSC package.
```

- [ ] **Step 8: Run docs-adjacent tests**

Run:

```bash
pytest -q tests/test_logging_refactor.py tests/test_config_identifiers.py
```

- [ ] **Step 9: Commit Task 10**

```bash
git add docs/current tutorials src/ldsc tests/test_logging_refactor.py tests/test_config_identifiers.py
git commit -m "docs: document allele-aware SNP identity"
```

---

## Task 11: End-To-End Verification And Cleanup

**Files:**

- Modify tests only if failures expose stale expectations.
- Do not change implementation unless a verification failure identifies a concrete defect.

- [ ] **Step 1: Scan for old mode aliases**

Run:

```bash
rg -n "rsid_alleles|chr_pos_alleles|SNPID|snp_id|chrpos|rsID\"|\"rsID|choices=\\(\" src tests docs tutorials
```

Expected:

- No old mode aliases accepted as `snp_identifier`.
- `SNPID`, `rsID`, and related strings may remain only as input column aliases or fixture column names.

- [ ] **Step 2: Scan for direct base-mode conditionals**

Run:

```bash
rg -n "== ['\\\"]chr_pos['\\\"]|== ['\\\"]rsid['\\\"]|!= ['\\\"]chr_pos['\\\"]|!= ['\\\"]rsid['\\\"]" src/ldsc
```

Expected:

- Remaining exact comparisons are only places where base and allele-aware modes intentionally differ.
- Family decisions use `identity_mode_family()`.

- [ ] **Step 3: Run focused workflow suite**

Run:

```bash
pytest -q tests/test_snp_identity.py tests/test_config_identifiers.py tests/test_column_inference.py tests/test_sumstats_munger.py tests/test_sumstats_liftover.py tests/test_ref_panel_builder.py tests/test_ref_panel.py tests/test_ldscore_workflow.py tests/test_regression_workflow.py tests/test_plink_io.py tests/test_output.py
```

Expected: all tests pass.

- [ ] **Step 4: Run full suite**

Run:

```bash
pytest -q
```

Expected: all tests pass.

- [ ] **Step 5: Run CLI help smoke tests**

Run:

```bash
python -m ldsc munge-sumstats --help
python -m ldsc build-ref-panel --help
python -m ldsc ldscore --help
python -m ldsc h2 --help
python -m ldsc partitioned-h2 --help
python -m ldsc rg --help
```

Expected:

- Help text shows exact four modes where choices are printed.
- Regression commands show `--allow-identity-downgrade`.

- [ ] **Step 6: Verify package artifact rejection**

Create small old-style artifacts in a temporary directory:

- sumstats metadata without `schema_version`
- LD-score manifest without `schema_version`
- canonical R2 parquet without `ldsc:schema_version`

Expected:

```text
Regenerate it with the current LDSC package.
```

- [ ] **Step 7: Commit verification-only expectation fixes**

If only stale tests/docs were changed:

```bash
git add tests docs tutorials
git commit -m "test: align identity refactor expectations"
```

If implementation defects were fixed, use a message naming the defect:

```bash
git add src tests
git commit -m "fix: handle allele-aware identity edge cases"
```

---

## Final Acceptance Checklist

- [ ] `GlobalConfig().snp_identifier == "chr_pos_allele_aware"`.
- [ ] Only the four exact public modes are accepted.
- [ ] Base modes never inspect allele columns.
- [ ] Base modes never emit allele-derived drop reasons.
- [ ] Base-mode duplicate effective keys drop all rows with `duplicate_identity`.
- [ ] Allele-aware invalid/missing/ambiguous allele rows are dropped and reported.
- [ ] Allele-aware multi-allelic base-key clusters are dropped package-wide.
- [ ] Restrictions without alleles match by base key.
- [ ] Restrictions with alleles match by effective key in allele-aware modes.
- [ ] Raw external R2 parquet works only in `rsid` and `chr_pos`.
- [ ] Canonical R2 parquet writes `A1_1`, `A2_1`, `A1_2`, `A2_2`.
- [ ] Reference metadata sidecars write identity provenance plus `A1`, `A2`.
- [ ] LD-score outputs include `A1`, `A2` in allele-aware modes.
- [ ] Old package-written artifacts are rejected with a regeneration message.
- [ ] Regression exact-mode compatibility is enforced by default.
- [ ] `--allow-identity-downgrade` works only within a mode family.
- [ ] h2 and partitioned-h2 orient `Z` to LD-score/reference alleles in allele-aware modes.
- [ ] rg still harmonizes ordered trait alleles after identity merge.
- [ ] Logs explain downgrade and dropped-row counts not stored in minimal manifests.
- [ ] `pytest -q` passes.

## Execution Strategy

This refactor is broad enough that each task should be implemented and committed independently. The safest order is exactly the task order above:

1. Build the shared identity service.
2. Update public modes and config defaults.
3. Add artifact provenance validation.
4. Refactor restrictions.
5. Refactor sumstats.
6. Refactor reference-panel builder.
7. Refactor R2 loading.
8. Refactor annotation and LD-score.
9. Refactor regression compatibility and orientation.
10. Update docs and error messages.
11. Run end-to-end verification.

Do not start workflow rewrites before Task 1 is complete. The package currently has identity logic spread across multiple modules; replacing it piecemeal without a shared service will leave inconsistent duplicate behavior.
