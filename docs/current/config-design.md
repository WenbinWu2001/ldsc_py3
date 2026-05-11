# Config Design: Immutable Config + Provenance-Carrying Results

## Implementation Status

This design is now implemented in the package. The shipped behavior matches the
core design: `GlobalConfig` remains immutable, workflow result objects carry
`config_snapshot`, and critical compatibility checks are enforced at combination
points.

Two implementation details are important to know:

- Compatibility checks run when both inputs carry real snapshots. Legacy objects
  with `config_snapshot=None` are still accepted for backward compatibility.
- Current `ldsc munge-sumstats` writes `sumstats.metadata.json` beside
  `sumstats.parquet` by default, or beside legacy `sumstats.sumstats.gz` when
  `--output-format tsv.gz` is selected; `load_sumstats()` recovers the original
  munge-time `GlobalConfig` from that thin sidecar. Row-level liftover drops are
  audited separately in the always-written `dropped_snps/dropped.tsv.gz` file.
  Older `.sumstats.gz` files without the metadata sidecar still warn and load
  with `config_snapshot=None`. Metadata sidecars that exist but lack
  `config_snapshot` are invalid and are not migrated.
- `load_ldscore_from_dir()` keeps strict format checks but treats missing or
  invalid manifest config provenance as unknown, warning and returning
  `config_snapshot=None`.

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

### 2. Computed result objects carry a config snapshot

Each computed artifact — `LDScoreResult`, `ChromLDScoreResult`, `AnnotationBundle`,
`SumstatsTable`, `RegressionDataset` — has a `config_snapshot` field. For
in-process results, that field stores the `GlobalConfig` frozen at the moment
of computation. This snapshot is the authoritative record of what assumptions
were active when that object was produced.

File-reloaded artifacts may instead use `config_snapshot=None` when the on-disk
format cannot prove its original settings. For example, `load_sumstats()` can
recover provenance from `sumstats.metadata.json` written by the current munger,
but warns and returns unknown provenance for older `.sumstats.gz` files that do
not have that sidecar.

This means reproducibility is structural when provenance exists: interrogating a
result object tells you the frozen config it was computed under, or tells you
explicitly that the provenance is unknown, regardless of what the global
registry currently holds.

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

The implementation also validates adjacent consistency boundaries:

- `LDScoreCalculator.run()` checks `AnnotationBundle.config_snapshot` against the
  active runtime `GlobalConfig`
- reference-panel loaders receive the active `GlobalConfig`, so parquet metadata
  build checks and `chr_pos` interpretation use the same runtime build that is
  later captured in LD-score result snapshots
- LD-score aggregation checks that per-chromosome `ChromLDScoreResult` snapshots
  agree before constructing the final `LDScoreResult`

```text
ConfigMismatchError: genome_build mismatch when combining SumstatsTable and
LDScoreResult: 'hg38' vs 'hg19'. These objects were computed under different
assumptions and cannot be safely merged.
```

For regression merges, the active SNP identifier mode is resolved in this order:
`LDScoreResult.config_snapshot.snp_identifier`, then
`SumstatsTable.config_snapshot.snp_identifier`, then the runner's active
`GlobalConfig.snp_identifier`, then the package default `chr_pos`. Disk-loaded
LD-score artifacts therefore use their manifest provenance when present; if old
manifests omit SNP-identifier provenance entirely, regression defaults to
coordinate identity rather than `SNP` labels.

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
as `RefPanelConfig.ref_panel_snps_file` and `LDScoreConfig.regression_snps_file`
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

Each in-process result carries `config_snapshot = GlobalConfig(genome_build="hg38", ...)`.
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
| `config_snapshot` field | Result classes | Frozen config from computation time, or `None` when provenance is genuinely unknown |

---

## Workflow-Specific SNP Controls

Standard LDSC still needs two distinct SNP-set controls, but they no longer live
on `GlobalConfig`. The refactored package assigns them to the workflow objects
that actually own those decisions.

### The hierarchy

```text
Annotation bundle rows B                      (AnnotationBuilder)
  intersects prepared reference panel A'     (RefPanelConfig ref-panel filters)
    then optional regression row subset C    (LDScoreConfig.regression_snps_file)
```

| Concept | Owner | CLI flag | When applied | Artifact effect |
| --- | --- | --- | --- | --- |
| Reference-panel SNP restriction | `RefPanelConfig.ref_panel_snps_file` | `--ref-panel-snps-file` | `RefPanel.load_metadata()`; then `LDScoreCalculator.compute_chromosome()` aligns `B_chrom` to the restricted panel before the kernel call | Shrinks the compute-time universe to `ld_reference_snps = B ∩ A'`; affects LD scores and count records |
| Reference-panel MAF/sample filters | `RefPanelConfig.maf_min`, `RefPanelConfig.keep_indivs_file` | `--maf-min`, `--keep-indivs-file` | `RefPanel.load_metadata()` and PLINK reader construction | Affects the prepared panel A' and LD computation; separate from `LDScoreConfig.common_maf_min` |
| Regression row restriction | `LDScoreConfig.regression_snps_file` | `--regression-snps-file` | After LD computation, when normalized/public rows are selected | Shrinks written rows to `ld_regression_snps = B ∩ A' ∩ C`; `regression_ld_scores` is embedded in the same row table |
| Common-count threshold | `LDScoreConfig.common_maf_min` | `--common-maf-min` | During count-vector computation after LD scores are computed | Affects only `common_reference_snp_count` / `common_reference_snp_counts`; does not change LD rows, LD scores, or the stored regression-universe LD score |

### What each control does

**`RefPanelConfig.ref_panel_snps_file`** — the *reference-panel SNP universe*

- `AnnotationBuilder.run()` still builds the full annotation universe `B`.
- `run_bed_to_annot()` and `ldsc annotate` do **not** apply this restriction.
- `RefPanel.load_metadata()` applies the restriction to the raw panel `A`, producing `A'`.
- `LDScoreCalculator.compute_chromosome()` then intersects the chromosome-local annotation bundle with that prepared metadata so the kernel sees `B_chrom ∩ A'_chrom`.

When `None`, the workflow uses the full reference panel `A`.

**`LDScoreConfig.regression_snps_file`** — the *regression row set*

- The restriction is loaded once into the `regression_snps` set `C`.
- LD scores and count totals are still computed over `ld_reference_snps = B ∩ A'`.
- Only the persisted LD-score rows are reduced to `ld_regression_snps = B ∩ A' ∩ C`.
- There is no separate public weight artifact. The selected baseline rows keep
  `regression_ld_scores`, the historical `w_ld` LD score computed over the
  regression SNP universe, in `ldscore.baseline.parquet`. Final h2 and rg
  regression weights are computed later in the regression kernel from this
  value plus model-specific quantities.

When `None`, the normalized/public row table uses all of `ld_reference_snps`.

### Typical configuration

```python
cfg = ldsc.GlobalConfig(genome_build="hg38", snp_identifier="chr_pos")
ldsc.set_global_config(cfg)

ref_panel = ldsc.RefPanelLoader(cfg).load(
    ldsc.RefPanelConfig(
        backend="parquet_r2",
        r2_dir="r2_ref_panel_1kg30x_1cM_hm3/hg38",
        ref_panel_snps_file="filters/reference_universe.txt",
    )
)

ldscore = ldsc.LDScoreCalculator().run(
    annot,
    ref_panel,
    ldsc.LDScoreConfig(
        ld_wind_cm=1.0,
        regression_snps_file="filters/hapmap3.txt",
    ),
    global_config=cfg,
)
```

For package-built parquet R2 panels, `RefPanelConfig.r2_bias_mode` and
`sample_size` can usually stay omitted because the loader reads
`ldsc:r2_bias` and `ldsc:n_samples` from the parquet schema. Set them explicitly
only when loading legacy or external raw-R2 parquet files that lack LDSC schema
metadata.

### Artifact contract

- Materialized query `.annot.gz` files and in-memory `AnnotationBundle` objects stay on the annotation universe `B`.
- Ordinary unpartitioned `run_ldscore()` calls may omit baseline/query annotation inputs; the workflow creates a synthetic all-ones `base` annotation after the reference panel has applied `ref_panel_snps_file`.
- Query annotations are valid only with explicit baseline annotations, so the synthetic `base` path is not used for partitioned/query LDSC.
- `ref_panel_snps_file` becomes visible only during LD-score calculation, when the workflow aligns `B_chrom` to `ref_panel.load_metadata(chrom)`.
- Count records are accumulated over `ld_reference_snps = B ∩ A'` and stored in `manifest.json`.
- Public `ldscore.baseline.parquet` and optional `ldscore.query.parquet` rows are `ld_regression_snps = B ∩ A' ∩ C`.

### Migration Notes

Stop passing these controls to `GlobalConfig()`:

- move reference-panel restriction to `RefPanelConfig(ref_panel_snps_file=...)`
- move regression row restriction to `LDScoreConfig(regression_snps_file=...)` or `run_ldscore(...)`

The old `--regression-snps` and `--print-snps` behavior is unified under
`--regression-snps-file`. LD-score outputs are fixed files under `output_dir`;
legacy `.l2.*` and `.w.l2.*` filenames are not emitted by the public writer.
Existing owned LD-score family artifacts are refused by default and require
`--overwrite` or `overwrite=True`, so reruns cannot silently replace artifacts
produced under a different configuration snapshot. With overwrite enabled, a
successful run removes stale owned siblings that the current run did not
produce, such as `ldscore.query.parquet` after switching from query LD scores
to baseline-only LD scores.

The same coherent-family rule applies to `munge-sumstats`, `partitioned-h2`,
and `annotate`: no-overwrite mode rejects any owned sibling, while successful
overwrites delete stale owned siblings and preserve unrelated files. The
`build-ref-panel` workflow is the documented expert exception; its overwrite
mode replaces current candidate panel artifacts but does not clean stale
optional target-build or `dropped_snps` siblings.

---

## Caveats and Known Limitations

**Reference-panel runtime state does not carry its own full `GlobalConfig`.**
`RefPanelConfig` owns backend and filtering options such as `r2_dir`,
`maf_min`, and `ref_panel_snps_file`. The shared identifier and genome-build
assumptions come from the `GlobalConfig` passed to `RefPanelLoader` and then
captured by the LD-score workflow result snapshots.

**`AnnotationBundle` does not carry chromosome-level genome-build metadata.**
Genome build for annotations is inferred from the reference metadata embedded in
`AnnotationBundle.metadata`. Compatibility validation trusts the `GlobalConfig`
snapshot rather than re-inferring the build from coordinates.

**Workflow-specific SNP controls are not part of `config_snapshot`.**
`config_snapshot` records shared assumptions such as `genome_build` and
`snp_identifier`. Per-run LD-score controls such as
`RefPanelConfig.ref_panel_snps_file` and `LDScoreConfig.regression_snps_file`
still materially affect the outputs, but callers who need to preserve that
provenance should persist the `RefPanelConfig` / `LDScoreConfig` they used
alongside the written artifacts.

**`load_sumstats()` recovers provenance for current disk artifacts.**
Curated `sumstats.parquet` and `.sumstats(.gz)` artifacts written by the current
munger have a neighboring `sumstats.metadata.json` sidecar that records a thin
metadata payload: schema marker, optional trait label, and the active or
inferred `GlobalConfig` snapshot. Detailed coordinate provenance, liftover
reports, HM3 provenance, output bookkeeping, and row counts are written to
`sumstats.log`; row-level liftover drops are written to
`dropped_snps/dropped.tsv.gz`. Neither belongs in the metadata sidecar. The
loader uses the sidecar to populate `config_snapshot`.
Older artifacts without the sidecar still emit a warning and load with
`config_snapshot=None`; pre-`config_snapshot` sidecars are treated as invalid
metadata rather than a backward-compatible format.

**Legacy LD-score directories may also have unknown provenance.**
Canonical LD-score directories written by the current workflow include a manifest
config snapshot. Older or malformed manifests may not. `load_ldscore_from_dir()`
warns and returns `LDScoreResult.config_snapshot=None` in that case while still
loading the baseline/query tables if the rest of the manifest is usable.

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
