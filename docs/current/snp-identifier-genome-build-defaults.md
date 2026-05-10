# SNP Identifier / Genome Build: Default Design

## Design Decisions

### The invariant pair

`GlobalConfig` has two fields that form a constraint pair:

| Field | Purpose |
|---|---|
| `snp_identifier` | `"rsid"` or `"chr_pos"` — how SNPs are identified |
| `genome_build` | `"hg19"`, `"hg38"`, `"auto"`, or `None` — coordinate system for CHR:POS |

`genome_build` is only meaningful in `chr_pos` mode. In `rsid` mode it is
irrelevant and is silently coerced to `None`.

### Why `chr_pos` is the package default

rsIDs have well-known anomalies: retired IDs, cross-build merges, and
platform-specific artefacts. `chr_pos` makes the coordinate system explicit and
unambiguous. The package actively encourages users toward coordinate-based SNP
identification by making it the zero-argument default.

### Why `"auto"` is the default companion for `chr_pos`

`"auto"` is the only honest value when the user has not yet specified a genome
build. It encodes "I'll infer the build from data when a workflow needs it"
without silently committing to `hg19` or `hg38`.

`None` is not appropriate as a `chr_pos` default because it is semantically
ambiguous — in `rsid` mode `None` means "irrelevant"; in `chr_pos` mode `None`
would mean "missing required information." Using the same literal for two
opposite semantics is confusing and error-prone.

### Validity table

| `snp_identifier` | `genome_build` | Valid? | Notes |
|---|---|---|---|
| `rsid` | `None` | yes | only valid rsid state |
| `rsid` | `"auto"` | **raise** | auto-inference is meaningless without coordinates |
| `rsid` | `"hg19"` / `"hg38"` | warn + coerce to `None` | build is irrelevant in rsid mode |
| `chr_pos` | `None` | **raise** | fix-it message: "Pass genome_build='auto' to infer from data." |
| `chr_pos` | `"auto"` | yes | **new package default** |
| `chr_pos` | `"hg19"` / `"hg38"` | yes | concrete coordinate mode |

`__post_init__` enforces this table. The `chr_pos + None` raise is new; it
converts a previously latent error into an early, actionable failure.

---

## Default State

```python
GlobalConfig()
# → snp_identifier="chr_pos", genome_build="auto"
```

The process-wide registry and `reset_global_config()` both return `GlobalConfig()`.
Before this design, the registry was `GlobalConfig(snp_identifier="rsid")` to
work around the invalid `chr_pos + None` class default. After this design, the
class default is valid so the registry no longer needs an override.

---

## Per-Workflow Behavior

### CLI: stricter than the Python API

CLI `--genome-build` defaults to `None` in every workflow parser. A CLI user in
`chr_pos` mode must type `--genome-build auto` (or a concrete build) explicitly.
This is intentional: CLI users benefit from the choice being visible in their
shell history and scripts. The API and CLI have different ergonomic defaults
because they serve different user models.

| Workflow | Identifier default | Genome-build default | Notes |
|---|---|---|---|
| `annotate` | `chr_pos` | `None` (CLI requires explicit) | Raises unless `--genome-build` is supplied. |
| `munge-sumstats` | `chr_pos` | `None` (CLI requires explicit) | Raises unless `--genome-build` is supplied. |
| `ldscore` | `chr_pos` | `None` (CLI requires explicit) | Raises unless `--genome-build` is supplied. |
| `build-ref-panel` | registry or `--snp-identifier` | **ignored** | Uses `--source-genome-build` (separate field); `GlobalConfig.genome_build` is never consulted. |
| `h2`, `partitioned-h2`, `rg` | registry at construction | registry at construction | On-disk provenance from LD-score manifests and sumstats sidecars takes precedence over the runner's live config. |
| Python `run_ldscore()` wrapper | registry | registry (`"auto"`) | Inherits `chr_pos + auto` from the registry; `auto` is resolved to `hg19`/`hg38` during inference. |

### `build-ref-panel` isolation

`build-ref-panel` is intentionally isolated from `GlobalConfig.genome_build`.
The PLINK source build is owned by `ReferencePanelBuildConfig.source_genome_build`
and is inferred from `.bim` coordinates when omitted. This keeps the panel-building
workflow self-contained: it neither reads nor writes the global genome-build
assumption. The `GlobalConfig` that `build-ref-panel` constructs internally holds
`genome_build="auto"` only to satisfy the construction invariant; that value is
never consumed by the workflow.

---

## Regression Provenance Hierarchy

Regression workflows (`h2`, `partitioned-h2`, `rg`) can see three independent
`GlobalConfig` sources for any given run:

1. The LD-score manifest snapshot — `LDScoreResult.config_snapshot`
2. The sumstats sidecar snapshot — `SumstatsTable.config_snapshot`
3. The runner's live `self.global_config` — from the process registry at construction time

**Pair validation** (`build_dataset`, lines 132–137) fires when both on-disk
snapshots are present. It compares the two artifacts against each other.
The runner's live config is **not** a participant in this check.

**Effective identifier resolution** (`_effective_snp_identifier_mode`) selects
the mode in priority order:

```
LD-score snapshot → sumstats snapshot → runner config → "rsid" fallback
```

The runner's live config is only the fallback for fully un-instrumented inputs.
Persisted provenance always dominates.

### Failure modes by provenance state

| LD manifest | Sumstats snapshot | Outcome |
|---|---|---|
| `rsid` | `rsid` | Pair check passes; merge uses `SNP` column. Runner's `chr_pos+auto` default is ignored. |
| `rsid` | `chr_pos` | `ConfigMismatchError("snp_identifier mismatch ... 'rsid' vs 'chr_pos'")`. Clear, actionable. |
| `rsid` | none | Pair check skipped; mode resolution returns `rsid` from manifest. Runner default is irrelevant. |
| none | none | Pair check skipped; mode resolution falls through to runner → `chr_pos`. If sumstats lack CHR/POS columns: column inference raises. If merge is empty: runtime error includes active config. |

The key implication: on-disk provenance is transparent to the new `chr_pos+auto`
registry default. An rsid-keyed pipeline whose artifacts carry their production
snapshots continues to work correctly without any user intervention.

### Empty-merge diagnostic

When scenario D (no provenance on either side) produces an empty merge, the
error message includes the active `GlobalConfig`:

```
No overlapping chr_pos SNPs remain after merging sumstats 'sumstats.parquet'
with 45321 LD-score rows. Check that snp_identifier and genome_build match.
Active config: GlobalConfig(snp_identifier='chr_pos', genome_build='auto', ...).
```

---

## Call Sites Affected by the New `chr_pos + None` Invariant

Any `GlobalConfig` construction that omits `genome_build` while `snp_identifier`
is `chr_pos` (explicitly or by class default) must supply `genome_build="auto"` or
a concrete build.

| File | Site | Change |
|---|---|---|
| `config.py` | field defaults | `genome_build = "auto"` |
| `config.py` | `__post_init__` | Add raise for `chr_pos + None` with fix-it message |
| `config.py` | `_GLOBAL_CONFIG` and `reset_global_config()` | `GlobalConfig()` — no explicit args needed |
| `ref_panel_builder.py` | `config_from_args` | Add `genome_build="auto"` (satisfies invariant; value is ignored by the workflow) |
| `sumstats_munger.py` | `_effective_sumstats_config` | Guard: `genome_build = coordinate_provenance.get("genome_build") or config.genome_build` so a `None` in old metadata does not override a valid config build; new sidecars serialize `coordinate_provenance`, while old `coordinate_metadata` sidecars are normalized on load |
| `regression_runner.py` | `_global_config_from_manifest` | For old chr_pos manifests with no recorded build, fall back to `"auto"` rather than `None` |
| `regression_runner.py` | `build_dataset` empty-merge error | Append `f"Active config: {self.global_config!r}."` |

Call sites in `ldscore_calculator.py` already resolve `genome_build` to a
concrete value or `"auto"` before constructing `GlobalConfig`. They are not
affected.
