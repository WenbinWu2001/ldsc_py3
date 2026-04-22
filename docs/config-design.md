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

```text
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
| --- | --- | --- |
| `genome_build` | Critical | Hard error (`ConfigMismatchError`) |
| `snp_identifier` | Critical | Hard error (`ConfigMismatchError`) |
| `log_level` | Advisory | Ignored during compatibility checks |
| `fail_on_missing_metadata` | Advisory | Ignored during compatibility checks |

The compatibility helper compares only `GlobalConfig` snapshots, so the table
above lists only `GlobalConfig` fields. Workflow-specific LD-score controls such
as `RefPanelSpec.ref_panel_snps_path` and `LDScoreConfig.regression_snps_path`
still matter to the result, but they are owned by those workflow objects rather
than by `GlobalConfig`.

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
| --- | --- | --- |
| `ConfigMismatchError` | `ldsc.config` | Raised on incompatible config combination |
| `validate_config_compatibility(a, b)` | `ldsc.config` | Public compatibility helper used by workflow merge points |
| `config_snapshot` field | All result classes | Frozen config from computation time |

---

## Workflow-Specific SNP Controls

Standard LDSC still needs two distinct SNP-set controls, but they no longer live
on `GlobalConfig`. The refactored package assigns them to the workflow objects
that actually own those decisions.

### The hierarchy

```text
Annotation bundle rows B                      (AnnotationBuilder)
  intersects prepared reference panel A'     (RefPanelSpec.ref_panel_snps_path)
    then optional regression row subset C    (LDScoreConfig.regression_snps_path)
```

| Concept | Owner | CLI flag | When applied | Artifact effect |
| --- | --- | --- | --- | --- |
| Reference-panel restriction | `RefPanelSpec.ref_panel_snps_path` | `--ref-panel-snps-path` | `RefPanel.load_metadata()`; then `LDScoreCalculator.compute_chromosome()` aligns `B_chrom` to the restricted panel before the kernel call | Shrinks the compute-time universe to `ld_reference_snps = B ∩ A'`; affects LD scores, `.M`, and `.M_5_50` |
| Regression row restriction | `LDScoreConfig.regression_snps_path` | `--regression-snps-path` | After LD computation, when normalized/public rows are selected | Shrinks written rows to `ld_regression_snps = B ∩ A' ∩ C`; `regr_weight` is embedded in the same row table |

### What each control does

**`RefPanelSpec.ref_panel_snps_path`** — the *reference-panel SNP universe*

- `AnnotationBuilder.run()` still builds the full annotation universe `B`.
- `run_bed_to_annot()` and `ldsc annotate` do **not** apply this restriction.
- `RefPanel.load_metadata()` applies the restriction to the raw panel `A`, producing `A'`.
- `LDScoreCalculator.compute_chromosome()` then intersects the chromosome-local annotation bundle with that prepared metadata so the kernel sees `B_chrom ∩ A'_chrom`.

When `None`, the workflow uses the full reference panel `A`.

**`LDScoreConfig.regression_snps_path`** — the *regression row set*

- The restriction is loaded once into the `regression_snps` set `C`.
- LD scores and count totals are still computed over `ld_reference_snps = B ∩ A'`.
- Only the normalized/public rows are reduced to `ld_regression_snps = B ∩ A' ∩ C`.
- There is no separate `.w.l2.ldscore.gz` artifact in the new format; the selected rows keep their regression weights in the embedded `regr_weight` column.

When `None`, the normalized/public row table uses all of `ld_reference_snps`.

### Typical configuration

```python
cfg = ldsc.GlobalConfig(genome_build="hg38", snp_identifier="chr_pos")
ldsc.set_global_config(cfg)

ref_panel = ldsc.RefPanelLoader(cfg).load(
    ldsc.RefPanelSpec(
        backend="parquet_r2",
        r2_table_paths="r2/reference.@.parquet",
        maf_metadata_paths="r2/reference_metadata.@.tsv.gz",
        ref_panel_snps_path="filters/reference_universe.txt",
    )
)

ldscore = ldsc.LDScoreCalculator().run(
    annot,
    ref_panel,
    ldsc.LDScoreConfig(
        ld_wind_cm=1.0,
        regression_snps_path="filters/hapmap3.txt",
    ),
    global_config=cfg,
)
```

### Artifact contract

- Materialized query `.annot.gz` files and in-memory `AnnotationBundle` objects stay on the annotation universe `B`.
- `ref_panel_snps_path` becomes visible only during LD-score calculation, when the workflow aligns `B_chrom` to `ref_panel.load_metadata(chrom)`.
- `.M` and `.M_5_50` are accumulated over `ld_reference_snps = B ∩ A'`.
- Normalized/public `.l2.ldscore.gz` rows and `LDScoreResult.ldscore_table` rows are `ld_regression_snps = B ∩ A' ∩ C`.

### Compatibility with older code

Stop passing these controls to `GlobalConfig()`:

- move reference-panel restriction to `RefPanelSpec(ref_panel_snps_path=...)`
- move regression row restriction to `LDScoreConfig(regression_snps_path=...)` or `run_ldscore(...)`

The old `--regression-snps` and `--print-snps` behavior is unified under
`--regression-snps-path`. There is no longer a separate default `.w.l2.ldscore.gz`
output artifact.

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

**Workflow-specific SNP controls are not part of `config_snapshot`.**
`config_snapshot` records shared assumptions such as `genome_build` and
`snp_identifier`. Per-run LD-score controls such as
`RefPanelSpec.ref_panel_snps_path` and `LDScoreConfig.regression_snps_path`
still materially affect the outputs, but callers who need to preserve that
provenance should persist the `RefPanelSpec` / `LDScoreConfig` they used
alongside the written artifacts.

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
