# PLANS

## Active: GlobalConfig default alignment — chr_pos + auto

**Design doc:** `docs/current/snp-identifier-genome-build-defaults.md`

### Background

`GlobalConfig` currently defaults to `snp_identifier="chr_pos"`, `genome_build=None`.
That pair is invalid (chr_pos requires a genome build). The package works around this
by registering `GlobalConfig(snp_identifier="rsid")` as the process default, creating
a mismatch between the class declaration and actual runtime behavior.

Goal: make `GlobalConfig()` produce a valid state that encodes the preferred SNP
identification convention (`chr_pos`) with the honest uncertainty sentinel (`"auto"`).

---

### Contract After This Change

- `GlobalConfig()` → `chr_pos + auto`. Always valid.
- `GlobalConfig(snp_identifier="chr_pos", genome_build=None)` → raises with clear message.
- `GlobalConfig(snp_identifier="rsid")` → `rsid + None`. Still valid.
- `GlobalConfig(snp_identifier="rsid", genome_build="auto")` → still raises (no change).
- Process default (`_GLOBAL_CONFIG`) → `chr_pos + auto`.

---

### Step 1 — `src/ldsc/config.py`

**1a. Change field default**

```python
# Line 108-109 — before
snp_identifier: SNPIdentifierMode = "chr_pos"
genome_build: GenomeBuildInput | None = None

# after
snp_identifier: SNPIdentifierMode = "chr_pos"
genome_build: GenomeBuildInput | None = "auto"
```

**1b. Add explicit chr_pos + None validation in `__post_init__`**

Insert after the existing `normalize_genome_build` call (around line 117), before
the `rsid + auto` check:

```python
if self.snp_identifier == "chr_pos" and self.genome_build is None:
    raise ValueError(
        "genome_build is required when snp_identifier='chr_pos'. "
        "Pass genome_build='auto' to infer from data, or 'hg19'/'hg38' explicitly."
    )
```

**1c. Simplify process-default and reset**

```python
# Line 130 — before
_GLOBAL_CONFIG: GlobalConfig = GlobalConfig(snp_identifier="rsid")

# after
_GLOBAL_CONFIG: GlobalConfig = GlobalConfig()

# reset_global_config (line ~153) — before
_GLOBAL_CONFIG = GlobalConfig(snp_identifier="rsid")

# after
_GLOBAL_CONFIG = GlobalConfig()
```

**1d. Update docstring** for `genome_build` parameter to say default is `"auto"`.

---

### Step 2 — `src/ldsc/ref_panel_builder.py`

`config_from_args` (line ~902) constructs `GlobalConfig` without `genome_build`.
With the new validation, if `snp_identifier="chr_pos"` (the default), this raises.

```python
# before
global_config = GlobalConfig(
    snp_identifier=snp_identifier,
    log_level=args.log_level,
)

# after
global_config = GlobalConfig(
    snp_identifier=snp_identifier,
    genome_build="auto",   # build-ref-panel ignores GlobalConfig.genome_build; "auto" satisfies the invariant
    log_level=args.log_level,
)
```

Note: `build-ref-panel` intentionally ignores `GlobalConfig.genome_build` (uses
`ReferencePanelBuildConfig.source_genome_build` instead). The `"auto"` value here
simply satisfies the construction invariant; it is never consumed by this workflow.

---

### Step 3 — `src/ldsc/sumstats_munger.py`

`_effective_sumstats_config` (line ~593) reads `genome_build` from
`coordinate_metadata`. If the metadata dict contains `"genome_build": None`
(from an older artifact), the current code passes `None` through and would now
raise on `chr_pos + None`.

```python
# before
genome_build = coordinate_metadata.get("genome_build", config.genome_build)

# after
genome_build = coordinate_metadata.get("genome_build") or config.genome_build
```

This makes `None` (or missing key) fall back to `config.genome_build`, which
is now `"auto"` for chr_pos mode. Old artifacts with an explicit `None` do not
corrupt the recovered config.

---

### Step 4 — `src/ldsc/regression_runner.py`

**4a. `_global_config_from_manifest` (line ~818)**

The snapshot recovery uses `snapshot.get("genome_build") or manifest.get("genome_build")`,
which evaluates to `None` for old manifests. If the recovered `snp_identifier` is
`"chr_pos"`, `GlobalConfig(snp_identifier="chr_pos", genome_build=None)` now raises
and falls into the existing `except Exception` block — returning `None` (unknown
provenance). This is acceptable, but more explicit:

```python
# before (inside GlobalConfig(...) construction)
genome_build=snapshot.get("genome_build") or manifest.get("genome_build"),

# after
genome_build=snapshot.get("genome_build") or manifest.get("genome_build") or (
    "auto" if (snapshot.get("snp_identifier") or manifest.get("snp_identifier") or "rsid") == "chr_pos"
    else None
),
```

This recovers old chr_pos manifests with a `"auto"` fallback instead of silently
dropping provenance. Rsid manifests continue to get `None`.

**4b. Empty-merge error message (line ~161-166)**

```python
# before
raise ValueError(
    f"No overlapping {identifier_mode} SNPs remain after merging sumstats '{source}' "
    f"with {len(ldscore_frame)} LD-score rows. Check that snp_identifier and genome_build match."
)

# after
raise ValueError(
    f"No overlapping {identifier_mode} SNPs remain after merging sumstats '{source}' "
    f"with {len(ldscore_frame)} LD-score rows. Check that snp_identifier and genome_build match. "
    f"Active config: {self.global_config!r}."
)
```

---

### Step 5 — Tests

Run the full suite before and after:

```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3
pytest
```

Specific tests to check or add:

1. **`test_global_config_registry.py`**: Add tests:
   - `GlobalConfig()` succeeds and produces `snp_identifier="chr_pos"`, `genome_build="auto"`.
   - `GlobalConfig(snp_identifier="chr_pos", genome_build=None)` raises `ValueError` with
     message containing "Pass genome_build='auto'".
   - `GlobalConfig(snp_identifier="rsid")` still produces `genome_build=None`.
   - `get_global_config()` returns `chr_pos + auto` after module import.
   - `reset_global_config()` returns `chr_pos + auto`.

2. **`test_config_identifiers.py`**: Verify existing chr_pos + None construction sites
   that were previously valid now raise.

3. **`test_regression_workflow.py`**: Check that empty-merge error message now includes
   `Active config:`.

---

### Step 6 — Verify no regressions in CLI

```bash
# build-ref-panel should still work
ldsc build-ref-panel --help

# ldscore should still require --genome-build in chr_pos mode
ldsc ldscore --snp-identifier chr_pos  # should raise without --genome-build

# regression workflows should pick up chr_pos+auto from registry
ldsc h2 --help
```

---

### Completion criteria

- [ ] `GlobalConfig()` valid, produces `chr_pos + auto`
- [ ] `GlobalConfig(snp_identifier="chr_pos", genome_build=None)` raises with clear message
- [ ] `_GLOBAL_CONFIG` and `reset_global_config()` use `GlobalConfig()`
- [ ] `ref_panel_builder.config_from_args` passes `genome_build="auto"` in construction
- [ ] `sumstats_munger._effective_sumstats_config` guards against `None` override
- [ ] `regression_runner._global_config_from_manifest` handles chr_pos + None recovery
- [ ] Empty-merge error includes `Active config: {self.global_config!r}`
- [ ] Full test suite passes
