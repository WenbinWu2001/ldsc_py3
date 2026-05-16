# SNP Identifier / Genome Build: Default Design

## Design Decisions

### The invariant pair

`GlobalConfig` has two fields that form a constraint pair:

| Field | Purpose |
|---|---|
| `snp_identifier` | one of exactly `"rsid"`, `"rsid_allele_aware"`, `"chr_pos"`, `"chr_pos_allele_aware"` — how SNPs are identified |
| `genome_build` | `"hg19"`, `"hg38"`, `"auto"`, or `None` — coordinate system for CHR:POS |

`genome_build` is only meaningful in the `chr_pos` family. In the `rsid` family
it is irrelevant and is silently coerced to `None`.

### Why `chr_pos_allele_aware` is the package default

rsIDs have well-known anomalies: retired IDs, cross-build merges, and
platform-specific artefacts. Coordinate identity makes the coordinate system
explicit, while allele-aware identity reduces accidental merges at the same
rsID or `CHR:POS`. The package therefore defaults to
`chr_pos_allele_aware`.

### Why `"auto"` is the default companion for coordinate-family modes

`"auto"` is the only honest value when the user has not yet specified a genome
build. It encodes "I'll infer the build from data when a workflow needs it"
without silently committing to `hg19` or `hg38`.

`None` is not appropriate as a coordinate-family default because it is
semantically ambiguous — in `rsid`-family modes `None` means "irrelevant"; in
coordinate-family modes `None` would mean "missing required information."
Using the same literal for two opposite semantics is confusing and error-prone.

### Identity semantics

Mode names are exact. Column aliases apply only to input headers, not to
`snp_identifier` values.

Base modes are fully allele-blind:

- `rsid` uses only `SNP`.
- `chr_pos` uses only `CHR:POS`.

Allele columns may be preserved in base modes, but they never affect identity,
duplicate filtering, retention, or drop reasons.

Allele-aware modes require usable `A1/A2` on package-written identity
artifacts: munged sumstats, reference-panel metadata, R2 parquet endpoints, and
LD-score artifacts. They drop missing, invalid/non-SNP, identical,
strand-ambiguous allele pairs, multi-allelic base-key clusters, and duplicate
effective-key clusters. Artifact duplicates are handled by computing the
effective merge key for the active mode and dropping all rows in duplicate-key
clusters.

### Validity table

| `snp_identifier` | `genome_build` | Valid? | Notes |
|---|---|---|---|
| `rsid` / `rsid_allele_aware` | `None` | yes | only valid rsID-family state |
| `rsid` / `rsid_allele_aware` | `"auto"` | **raise** | auto-inference is meaningless without coordinates |
| `rsid` / `rsid_allele_aware` | `"hg19"` / `"hg38"` | warn + coerce to `None` | build is irrelevant in rsID-family modes |
| `chr_pos` / `chr_pos_allele_aware` | `None` | **raise** | fix-it message: "Pass genome_build='auto' to infer from data." |
| `chr_pos` / `chr_pos_allele_aware` | `"auto"` | yes | `GlobalConfig()` default companion |
| `chr_pos` / `chr_pos_allele_aware` | `"hg19"` / `"hg38"` | yes | concrete coordinate mode |

`__post_init__` enforces this table. The coordinate-family `None` raise
converts a latent error into an early, actionable failure.

---

## Default State

```python
GlobalConfig()
# → snp_identifier="chr_pos_allele_aware", genome_build="auto"
```

The process-wide registry and `reset_global_config()` both return `GlobalConfig()`.
Earlier designs used `GlobalConfig(snp_identifier="rsid")` or
`GlobalConfig(snp_identifier="chr_pos", genome_build="auto")` as defaults.
The current class default is valid and allele-aware, so the registry no longer
needs an override.

---

## Per-Workflow Behavior

### CLI: stricter than the Python API

CLI `--genome-build` defaults to `None` in most workflow parsers. A CLI user in
a coordinate-family mode must type `--genome-build auto` (or a concrete build)
explicitly unless that workflow can infer the build from another current
artifact.
This is intentional: CLI users benefit from the choice being visible in their
shell history and scripts. The API and CLI have different ergonomic defaults
because they serve different user models.

| Workflow | Identifier default | Genome-build default | Notes |
|---|---|---|---|
| `annotate` | `chr_pos_allele_aware` | `None` (CLI requires explicit for coordinate-family modes) | Raises unless `--genome-build` is supplied or inferable by that workflow. |
| `munge-sumstats` | `chr_pos_allele_aware` | workflow default/inference | Requires usable `A1/A2`; rerun with `--snp-identifier chr_pos` or `--snp-identifier rsid` to run without allele-aware identity. The removed `--no-alleles` flag is not accepted. |
| `ldscore` | `chr_pos_allele_aware` | `None` (CLI requires explicit for coordinate-family modes) | Allele-aware parquet mode requires package-built canonical R2 endpoint alleles. |
| `build-ref-panel` | registry or `--snp-identifier` | **ignored** | Uses `--source-genome-build` (separate field); `GlobalConfig.genome_build` is never consulted. |
| `h2`, `partitioned-h2`, `rg` | registry at construction | registry at construction | On-disk provenance from LD-score and sumstats root `metadata.json` contracts takes precedence over the runner's live config. |
| Python `run_ldscore()` wrapper | registry | registry (`"auto"`) | Inherits `chr_pos_allele_aware + auto` from the registry; `auto` is resolved to `hg19`/`hg38` during inference. |

### `build-ref-panel` isolation

`build-ref-panel` is intentionally isolated from `GlobalConfig.genome_build`.
The PLINK source build is owned by `ReferencePanelBuildConfig.source_genome_build`
and is inferred from `.bim` coordinates when omitted. This keeps the panel-building
workflow self-contained: it neither reads nor writes the global genome-build
assumption. The `GlobalConfig` that `build-ref-panel` constructs internally holds
`genome_build="auto"` only to satisfy the construction invariant; that value is
never consumed by the workflow.

---

## Regression Provenance Hierarchy And Downgrade

Regression workflows (`h2`, `partitioned-h2`, `rg`) can see three independent
`GlobalConfig` sources for any given run:

1. The LD-score root metadata snapshot — `LDScoreResult.config_snapshot`
2. The reconstructed sumstats artifact snapshot — `SumstatsTable.config_snapshot`
3. The runner's live `self.global_config` — from the process registry at construction time

**Pair validation** (`build_dataset`, lines 132–137) fires when both on-disk
snapshots are present. It compares the two artifacts against each other.
The runner's live config is **not** a participant in this check.

Effective identifier resolution selects the mode in priority order:

```
LD-score snapshot → sumstats snapshot → runner config → "chr_pos_allele_aware" fallback
```

The runner's live config is only the fallback for manually constructed in-memory
objects that do not carry provenance. Package-written disk artifacts must carry
current identity provenance, and persisted provenance always dominates.

By default, regression requires exact mode compatibility when both artifacts
carry current provenance. `--allow-identity-downgrade` is regression-only. It
allows same-family allele-aware/base mixes to run under the base mode, with
drop-all duplicate filtering recomputed under that base mode. rsID-family and
coordinate-family modes never mix. Original modes and dropped duplicate counts
are logged, not persisted in detail metadata.

### Failure modes by provenance state

| LD metadata | Sumstats snapshot | Outcome |
|---|---|---|
| `rsid` | `rsid` | Pair check passes; merge uses `SNP` column. Runner's `chr_pos_allele_aware+auto` default is ignored. |
| `rsid_allele_aware` | `rsid` | Raises unless `--allow-identity-downgrade`; with downgrade, merge uses `SNP` and base-mode duplicate cleanup. |
| `rsid` | `chr_pos` | Error: rsID-family and coordinate-family modes cannot mix. |
| `chr_pos_allele_aware` | `chr_pos_allele_aware` | Pair check passes; merge uses `CHR:POS:<allele_set>`. |
| `chr_pos_allele_aware` | missing current provenance | Package-written sumstats are rejected; regenerate them with the current LDSC package. |
| missing current provenance | missing current provenance | Package-written artifacts are rejected; manually constructed in-memory objects without provenance fall back to runner config and still undergo column/allele validation. |

The key implication: on-disk provenance is transparent to the
`chr_pos_allele_aware + auto` registry default. An rsID-keyed pipeline whose
artifacts carry current production snapshots continues to work correctly
without any user intervention.

### Empty-merge diagnostic

When scenario D (no provenance on either side) produces an empty merge, the
error message includes the active `GlobalConfig`:

```
No overlapping chr_pos_allele_aware SNPs remain after merging sumstats 'sumstats.parquet'
with 45321 LD-score rows. Check that snp_identifier and genome_build match.
Active config: GlobalConfig(snp_identifier='chr_pos_allele_aware', genome_build='auto', ...).
```

---

## Call Sites Affected By The Coordinate-Family `None` Invariant

Any `GlobalConfig` construction that omits `genome_build` while
`snp_identifier` is in the `chr_pos` family (explicitly or by class default)
must supply `genome_build="auto"` or a concrete build.

| File | Site | Change |
|---|---|---|
| `config.py` | field defaults | `genome_build = "auto"` |
| `config.py` | `__post_init__` | Raise for coordinate-family `genome_build=None` with fix-it message |
| `config.py` | `_GLOBAL_CONFIG` and `reset_global_config()` | `GlobalConfig()` — no explicit args needed |
| `ref_panel_builder.py` | `config_from_args` | Add `genome_build="auto"` (satisfies invariant; value is ignored by the workflow) |
| `sumstats_munger.py` | `_effective_sumstats_config` | Guard: `genome_build = coordinate_metadata.get("genome_build") or config.genome_build` so the effective in-memory snapshot uses inferred or lifted coordinates when available; new sidecars serialize only minimal identity provenance used to reconstruct that snapshot on reload |
| `regression_runner.py` | LD-score metadata loading | Current root metadata must carry minimal identity provenance; old package-written artifacts should be regenerated |
| `regression_runner.py` | `build_dataset` empty-merge error | Append `f"Active config: {self.global_config!r}."` |

Call sites in `ldscore_calculator.py` already resolve `genome_build` to a
concrete value or `"auto"` before constructing `GlobalConfig`. They are not
affected.
