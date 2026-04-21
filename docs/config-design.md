# Config Design: Immutable Config + Provenance-Carrying Results

## Implementation Status

This design is now implemented in the package. The shipped behavior matches the
core design: `GlobalConfig` remains immutable, workflow result objects carry
`config_snapshot`, and critical compatibility checks are enforced at combination
points.

Two implementation details are important to know:

- Compatibility checks run when both inputs carry real snapshots. Legacy objects
  with `config_snapshot=None` are still accepted for backward compatibility.
- `load_sumstats()` cannot recover the original munge-time config from disk, so
  it emits a warning and attaches the current global config as a proxy snapshot.

## The Problem This Design Solves

Genome build (`hg19`/`hg38`) and SNP identifier mode (`rsid`/`chr_pos`) are global
analysis assumptions that every step of the pipeline depends on. If these settings
differ between the annotation, reference-panel, and regression steps — even
temporarily, because of a notebook cell re-run — the mismatch will silently corrupt
results rather than raise an error.

The original design allowed mutating the process-wide default via `set_global_config()`
at any time. Any result computed before or after a call to `set_global_config()` would
be indistinguishable to downstream steps. This is the core danger.

---

## Design Principles

### 1. `GlobalConfig` is structurally immutable

`GlobalConfig` is a `frozen=True` dataclass. "Changing" the config always produces a
new object; nothing that already holds a reference to the old config is affected.

```python
cfg = GlobalConfig(genome_build="hg38", snp_identifier="chr_pos")
# cfg.genome_build = "hg19"  → raises FrozenInstanceError immediately
cfg2 = GlobalConfig(genome_build="hg19", snp_identifier="chr_pos")  # correct way
```

### 2. Every result object carries a config snapshot

Each computed artifact — `LDScoreResult`, `ChromLDScoreResult`, `AnnotationBundle`,
`SumstatsTable`, `RegressionDataset` — embeds a `config_snapshot: GlobalConfig` frozen
at the moment of computation. This snapshot is the authoritative record of what
assumptions were active when that object was produced.

This means reproducibility is structural: you can always interrogate a result object
to find out what config it was computed under, regardless of what the global registry
currently holds.

### 3. Compatibility is validated at combination points

The highest-leverage moment to catch a mismatch is when two result objects are
**combined** — not when they are constructed. Combination points include:

- `RegressionRunner.build_dataset(sumstats_table, ldscore_result)` — merges a
  `SumstatsTable` with an `LDScoreResult`
- `LDScoreCalculator.run(annotation_bundle, ref_panel, ...)` — combines an
  `AnnotationBundle` with a `RefPanel`

At each of these points, the package calls `validate_config_compatibility()` on the
config snapshots of both inputs. Mismatches on critical fields (`genome_build`,
`snp_identifier`) raise `ConfigMismatchError` with a precise diagnostic message.

The implementation also validates two adjacent consistency boundaries:

- `LDScoreCalculator.run()` checks `AnnotationBundle.config_snapshot` against the
  active runtime `GlobalConfig`
- `LDScoreCalculator.run()` checks `RefPanelSpec.genome_build` against the active
  runtime `GlobalConfig.genome_build`
- LD-score aggregation checks that per-chromosome `ChromLDScoreResult` snapshots
  agree before constructing the final `LDScoreResult`

```
ConfigMismatchError: genome_build mismatch when combining SumstatsTable and
LDScoreResult: 'hg38' vs 'hg19'. These objects were computed under different
assumptions and cannot be safely merged.
```

### 4. The global registry is a convenience default only

`set_global_config()` and `get_global_config()` are retained for notebook UX:
setting the config once at the top of a notebook is ergonomic. But the registry is
demoted: workflow functions read from it only as a fallback when no explicit config
is passed, and the snapshot captured at computation time — not the registry — is what
travels with results.

Calling `set_global_config()` mid-analysis will not corrupt previously computed
results, because those results carry their own frozen snapshots.

---

## Critical vs. Advisory Fields

| Field | Class | Mismatch behavior |
|---|---|---|
| `genome_build` | Critical | Hard error (`ConfigMismatchError`) |
| `snp_identifier` | Critical | Hard error (`ConfigMismatchError`) |
| `log_level` | Advisory | Ignored during compatibility checks |
| `fail_on_missing_metadata` | Advisory | Ignored during compatibility checks |
| `restrict_snps_path` | Advisory | Warning only |

---

## Usage Patterns

### Recommended: set once at the top of a notebook

```python
import ldsc
ldsc.set_global_config(ldsc.GlobalConfig(genome_build="hg38", snp_identifier="chr_pos"))

annot = ldsc.AnnotationBuilder(...).run(source_spec)
ldscore = ldsc.LDScoreCalculator().run(annot, ref_panel, ldscore_cfg)
regression = ldsc.RegressionRunner().build_dataset(sumstats, ldscore)
```

Each result carries `config_snapshot = GlobalConfig(genome_build="hg38", ...)`.
If you accidentally re-run a cell with a different `set_global_config()` call, the
results you already have are unaffected — and combining them with new results will
raise `ConfigMismatchError` at merge time.

### Explicit passing (safest for library code)

```python
cfg = ldsc.GlobalConfig(genome_build="hg38", snp_identifier="chr_pos")
annot = ldsc.AnnotationBuilder(global_config=cfg).run(source_spec)
ldscore = ldsc.LDScoreCalculator().run(annot, ref_panel, ldscore_cfg, global_config=cfg)
```

### What you cannot do

```python
# This will raise ConfigMismatchError at build_dataset time:
ldscore_hg38 = compute_with(genome_build="hg38")
sumstats_hg19 = munge_with(genome_build="hg19")
regression = runner.build_dataset(sumstats_hg19, ldscore_hg38)  # ← hard error here
```

---

## New Public API Surface

| Name | Location | Purpose |
|---|---|---|
| `ConfigMismatchError` | `ldsc.config` | Raised on incompatible config combination |
| `validate_config_compatibility(a, b)` | `ldsc.config` | Public compatibility helper used by workflow merge points |
| `config_snapshot` field | All result classes | Frozen config from computation time |

---

## Caveats and Known Limitations

**`RefPanelSpec.genome_build` is not a full `GlobalConfig`.**
The `RefPanelSpec` dataclass carries its own `genome_build` string field (set during
panel construction). Compatibility checks between a `RefPanel` and other results use
this field directly, not a full `GlobalConfig` snapshot. A future refactor could
promote this to a full snapshot, but for now the field is checked explicitly.

**`AnnotationBundle` does not carry chromosome-level genome-build metadata.**
Genome build for annotations is inferred from the reference metadata embedded in
`AnnotationBundle.metadata`. Compatibility validation trusts the `GlobalConfig`
snapshot rather than re-inferring the build from coordinates.

**`restrict_snps_path` is advisory.**
Two results computed with different `restrict_snps_path` values will not raise a hard
error. A warning is emitted. This is intentional: restricting SNPs is a filter, not
a coordinate-space change, so the results remain combinable with a caveat.

**`load_sumstats()` attaches a proxy snapshot, not original provenance.**
Curated `.sumstats(.gz)` artifacts on disk do not embed the `GlobalConfig` that was
active when they were munged. The loader therefore emits a warning and attaches the
current global config as a best-effort proxy snapshot. This supports notebook
ergonomics, but users should treat the snapshot on a loaded table as caller-time
context rather than guaranteed historical provenance.

**Legacy objects with `config_snapshot=None` skip strict compatibility checks.**
The package preserves backward compatibility for objects constructed before this
design shipped, or for file-reloaded artifacts that do not carry full provenance.
When either side of a compatibility boundary has `config_snapshot=None`, the merge
is allowed to proceed rather than raising immediately. This is intentional: the
package prefers explicit checks when provenance exists, without fabricating a false
history for older objects.

**Notebook re-execution order still matters for the registry.**
The global registry is process-wide. If a user calls `set_global_config()` in a cell,
re-runs an earlier cell that did not capture a config explicitly, and then combines the
results, they will get a `ConfigMismatchError`. This is the intended behavior — the
error message explains the mismatch clearly. It is not a silent corruption.

**The `run_ldscore()` and `run_build_ref_panel()` convenience wrappers read the global
registry at call time.** Explicit `global_config=` arguments are always preferred.
