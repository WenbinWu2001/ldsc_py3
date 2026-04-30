# PLANS

## Active: GlobalConfig default alignment — chr_pos + auto

**Design doc:** `docs/current/snp-identifier-genome-build-defaults.md`

### Goal

Make `GlobalConfig()` always produce a valid object that encodes the preferred
SNP identification convention. Current state: `GlobalConfig()` produces
`chr_pos + None`, which is an invalid pair that several workflows silently reject
later. The process registry works around this by registering `rsid + None`
instead — creating a mismatch between the class declaration and actual runtime.

Target state:

- `GlobalConfig()` → `chr_pos + auto`. Always valid. No override needed in registry.
- `GlobalConfig(snp_identifier="chr_pos", genome_build=None)` → raises immediately
  with: `"genome_build is required when snp_identifier='chr_pos'. Pass genome_build='auto' to infer from data, or 'hg19'/'hg38' explicitly."`
- `GlobalConfig(snp_identifier="rsid")` → `rsid + None`. Unchanged.

---

### Step 1 — `src/ldsc/config.py`

**1a. Change `genome_build` field default**

```python
# before
genome_build: GenomeBuildInput | None = None

# after
genome_build: GenomeBuildInput | None = "auto"
```

**1b. Add explicit `chr_pos + None` validation in `__post_init__`**

Insert after `object.__setattr__(self, "genome_build", normalize_genome_build(self.genome_build))`
and before the existing `rsid + auto` check:

```python
if self.snp_identifier == "chr_pos" and self.genome_build is None:
    raise ValueError(
        "genome_build is required when snp_identifier='chr_pos'. "
        "Pass genome_build='auto' to infer from data, or 'hg19'/'hg38' explicitly."
    )
```

#### 1c. Simplify process default and reset

```python
# _GLOBAL_CONFIG line — before
_GLOBAL_CONFIG: GlobalConfig = GlobalConfig(snp_identifier="rsid")

# after
_GLOBAL_CONFIG: GlobalConfig = GlobalConfig()

# reset_global_config body — before
_GLOBAL_CONFIG = GlobalConfig(snp_identifier="rsid")

# after
_GLOBAL_CONFIG = GlobalConfig()
```

**1d. Update the `genome_build` docstring parameter line** to say default is `"auto"`.

---

### Step 2 — `src/ldsc/ref_panel_builder.py`

`config_from_args` (around line 902) constructs `GlobalConfig` without
`genome_build`. With the new validation, `chr_pos + None` now raises.

```python
# before
global_config = GlobalConfig(
    snp_identifier=snp_identifier,
    log_level=args.log_level,
)

# after
global_config = GlobalConfig(
    snp_identifier=snp_identifier,
    genome_build="auto",
    log_level=args.log_level,
)
```

Note: `build-ref-panel` intentionally ignores `GlobalConfig.genome_build`; it
uses `ReferencePanelBuildConfig.source_genome_build` instead. The `"auto"` here
only satisfies the construction invariant.

---

### Step 3 — `src/ldsc/sumstats_munger.py`

`_effective_sumstats_config` (around line 593) reads `genome_build` from
`coordinate_metadata`. If that dict contains `"genome_build": None` (from an
older artifact on disk), the current code passes `None` through and the new
validation would raise on `chr_pos + None`.

```python
# before
genome_build = coordinate_metadata.get("genome_build", config.genome_build)

# after
genome_build = coordinate_metadata.get("genome_build") or config.genome_build
```

`None` (explicit) and missing key both fall back to `config.genome_build`, which
is `"auto"` under the new default.

---

### Step 4 — `src/ldsc/regression_runner.py`

**4a. `_global_config_from_manifest` (around line 818)**

The recovered `genome_build` can be `None` for old chr_pos manifests that did
not record a build. Rather than letting that fall into the `except` block and
silently returning `None` provenance, fall back to `"auto"`:

```python
# before
genome_build=snapshot.get("genome_build") or manifest.get("genome_build"),

# after — determine the recovered snp_identifier first, then pick the right fallback
_recovered_snp = (
    snapshot.get("snp_identifier") or manifest.get("snp_identifier") or "rsid"
)
genome_build=(
    snapshot.get("genome_build")
    or manifest.get("genome_build")
    or ("auto" if _recovered_snp == "chr_pos" else None)
),
```

#### 4b. Empty-merge error message (around line 161–166)

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

Activate environment first:

```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3
```

Run full suite:

```bash
pytest
```

Specific cases to add in `tests/test_global_config_registry.py`:

1. `GlobalConfig()` succeeds → `snp_identifier="chr_pos"`, `genome_build="auto"`.
2. `GlobalConfig(snp_identifier="chr_pos", genome_build=None)` → raises `ValueError`
   with message containing `"Pass genome_build='auto'"`.
3. `GlobalConfig(snp_identifier="rsid")` → `genome_build=None`. Still valid.
4. `get_global_config()` after module import → `snp_identifier="chr_pos"`, `genome_build="auto"`.
5. `reset_global_config()` → `snp_identifier="chr_pos"`, `genome_build="auto"`.

In `tests/test_regression_workflow.py`: verify the empty-merge error string
contains `"Active config:"`.

---

### Completion Criteria

- [ ] `GlobalConfig()` valid; produces `chr_pos + auto`
- [ ] `GlobalConfig(snp_identifier="chr_pos", genome_build=None)` raises with fix-it message
- [ ] `_GLOBAL_CONFIG` and `reset_global_config()` use `GlobalConfig()` with no args
- [ ] `ref_panel_builder.config_from_args` passes `genome_build="auto"`
- [ ] `sumstats_munger._effective_sumstats_config` guards against `None` override from old metadata
- [ ] `regression_runner._global_config_from_manifest` falls back to `"auto"` for chr_pos manifests with no recorded build
- [ ] Empty-merge error includes `Active config: {self.global_config!r}`
- [ ] Full `pytest` suite passes with no new failures
