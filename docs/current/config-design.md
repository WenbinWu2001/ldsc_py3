# Config Design: Immutable Config + Provenance-Carrying Results

## Implementation Status

This design is now implemented in the package. The shipped behavior matches the
core design: `GlobalConfig` remains immutable, workflow result objects carry
`config_snapshot`, and critical compatibility checks are enforced at combination
points.

Three implementation details are important to know:

- Compatibility checks run against the current provenance recorded on
  package-written artifacts. Old package-written sumstats and LD-score artifacts
  without the current schema/provenance contract are rejected and must be
  regenerated with the current LDSC package.
- Current `ldsc munge-sumstats` writes root `metadata.json` beside
  `sumstats.parquet` by default, or beside legacy `sumstats.sumstats.gz` when
  `--output-format tsv.gz` is selected; `load_sumstats()` recovers the original
  downstream compatibility `GlobalConfig` from that thin sidecar. In
  coordinate-family modes, that snapshot stores the final output genome build;
  in rsid-family modes, it stores `genome_build=None`. Row-level liftover drops are
  audited separately in the always-written
  `diagnostics/dropped_snps/dropped.tsv.gz` file.
  Package-written sumstats artifacts without current identity provenance are not
  migrated.
- `load_ldscore_from_dir()` keeps strict metadata checks and rejects missing or
  invalid package-written root metadata identity provenance with a regeneration
  message.

## The Problem This Design Solves

Genome build (`hg19`/`hg38`) and SNP identifier mode are global analysis
assumptions that every step of the pipeline depends on. LDSC supports exactly
four SNP identifier modes: `rsid`, `rsid_allele_aware`, `chr_pos`, and
`chr_pos_allele_aware`; the default is `chr_pos_allele_aware`. Mode names are
exact. Column aliases apply to input headers only. If these settings differ
between the annotation, reference-panel, and regression steps — even
temporarily, because of a notebook cell re-run — the mismatch will silently
corrupt results rather than raise an error.

The original design allowed mutating the process-wide default via `set_global_config()`
at any time. Any result computed before or after a call to `set_global_config()` would
be indistinguishable to downstream steps. This is the core danger.

---

## Design Principles

### 1. `GlobalConfig` is structurally immutable

`GlobalConfig` is a `frozen=True` dataclass. "Changing" the config always produces a
new object; nothing that already holds a reference to the old config is affected.

```python
cfg = GlobalConfig(genome_build="hg38", snp_identifier="chr_pos_allele_aware")
# cfg.genome_build = "hg19"  → raises FrozenInstanceError immediately
cfg2 = GlobalConfig(genome_build="hg19", snp_identifier="chr_pos_allele_aware")  # correct way
```

### 2. Computed result objects carry a config snapshot

Each computed artifact — `LDScoreResult`, `ChromLDScoreResult`, `AnnotationBundle`,
`SumstatsTable`, `RegressionDataset` — has a `config_snapshot` field. For
in-process results, that field stores the `GlobalConfig` frozen at the moment
of computation. This snapshot is the authoritative record of what assumptions
were active when that object was produced.

File-reloaded package-written artifacts must prove their original settings with
current identity provenance. For example, `load_sumstats()` recovers provenance
from root `metadata.json` written by the current munger, while older
package-written `.sumstats.gz` files without the current sidecar are rejected
and must be regenerated.

This means reproducibility is structural for package-written artifacts:
interrogating a loaded result object tells you the frozen config it was computed
under, regardless of what the global registry currently holds.

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
`GlobalConfig.snp_identifier`, then the package default
`chr_pos_allele_aware`. Disk-loaded LD-score artifacts therefore use their
metadata provenance when present; old package-written artifacts without the
current schema/provenance contract must be regenerated.

`--allow-identity-downgrade` is regression-only. It allows same-family
allele-aware/base mixes to run under the base mode and logs the original modes
plus dropped duplicate effective-key rows. rsID-family modes and
coordinate-family modes never mix.

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
ldsc.set_global_config(ldsc.GlobalConfig(genome_build="hg38", snp_identifier="chr_pos_allele_aware"))

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
cfg = ldsc.GlobalConfig(genome_build="hg38", snp_identifier="chr_pos_allele_aware")
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
    then optional regression row subset C    (LDScoreConfig regression SNP controls)
```

| Concept | Owner | CLI flag | When applied | Artifact effect |
| --- | --- | --- | --- | --- |
| Reference-panel SNP restriction | `RefPanelConfig.ref_panel_snps_file` or `use_hm3_ref_panel_snps` | `--ref-panel-snps-file` or `--use-hm3-ref-panel-snps` | `RefPanel.load_metadata()`; then `LDScoreCalculator.compute_chromosome()` aligns `B_chrom` to the restricted panel before the kernel call | Shrinks the compute-time universe to `ld_reference_snps = B ∩ A'`; affects LD scores and count records |
| Reference-panel MAF/sample filters | `RefPanelConfig.maf_min`, `RefPanelConfig.keep_indivs_file` | `--maf-min`, `--keep-indivs-file` | `RefPanel.load_metadata()` and PLINK reader construction | Affects the prepared panel A' and LD computation; separate from `LDScoreConfig.common_maf_min` |
| Regression row restriction | `LDScoreConfig.regression_snps_file` or `use_hm3_regression_snps` | `--regression-snps-file` or `--use-hm3-regression-snps` | After LD computation, when normalized/public rows are selected | Shrinks written rows to `ld_regression_snps = B ∩ A' ∩ C`; `regression_ld_scores` is embedded in the same row table |
| Common-count threshold | `LDScoreConfig.common_maf_min` | `--common-maf-min` | During count-vector computation after LD scores are computed | Affects only `common_reference_snp_count` / `common_reference_snp_counts`; does not change LD rows, LD scores, or the stored regression-universe LD score |

### What each control does

**`RefPanelConfig.ref_panel_snps_file` / `use_hm3_ref_panel_snps`** — the *reference-panel SNP universe*

- `AnnotationBuilder.run()` still builds the full annotation universe `B`.
- `run_bed_to_annot()` and `ldsc annotate` do **not** apply this restriction.
- `RefPanel.load_metadata()` applies the explicit restriction or packaged HM3
  restriction to the raw panel `A`, producing `A'`.
- SNP restriction files are identity-only. Duplicate restriction keys collapse
  to one retained key, and non-identity columns such as `CM`, `MAF`, or other
  metadata are ignored rather than carried into LD-score or regression outputs.
- `LDScoreCalculator.compute_chromosome()` then intersects the chromosome-local annotation bundle with that prepared metadata so the kernel sees `B_chrom ∩ A'_chrom`.

When `None`, the workflow uses the full reference panel `A`.

**`LDScoreConfig.regression_snps_file` / `use_hm3_regression_snps`** — the *regression row set*

- The explicit restriction file or packaged HM3 map is loaded once into the
  `regression_snps` set `C`.
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
cfg = ldsc.GlobalConfig(genome_build="hg38", snp_identifier="chr_pos_allele_aware")
ldsc.set_global_config(cfg)

ref_panel = ldsc.RefPanelLoader(cfg).load(
    ldsc.RefPanelConfig(
        backend="parquet_r2",
        r2_dir="r2_ref_panel_1kg30x_1cM_hm3/hg38",
        use_hm3_ref_panel_snps=True,
    )
)

ldscore = ldsc.LDScoreCalculator().run(
    annot,
    ref_panel,
    ldsc.LDScoreConfig(
        ld_wind_cm=1.0,
        use_hm3_regression_snps=True,
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
- Ordinary unpartitioned `run_ldscore()` calls may omit baseline/query annotation inputs; the workflow creates a synthetic all-ones `base` annotation after the reference panel has applied any explicit or HM3 reference-panel restriction.
- Query annotations are valid only with explicit baseline annotations, so the synthetic `base` path is not used for partitioned/query LDSC.
- Reference-panel SNP restrictions become visible only during LD-score calculation, when the workflow aligns `B_chrom` to `ref_panel.load_metadata(chrom)`.
- Count records are accumulated over `ld_reference_snps = B ∩ A'` and stored in LD-score root `metadata.json`.
- Public `ldscore.baseline.parquet` and optional `ldscore.query.parquet` rows are `ld_regression_snps = B ∩ A' ∩ C`.

### Migration Notes

Stop passing these controls to `GlobalConfig()`:

- move reference-panel restriction to `RefPanelConfig(ref_panel_snps_file=...)`
- move regression row restriction to `LDScoreConfig(regression_snps_file=...)` or `run_ldscore(...)`
- for the packaged HM3 map, prefer `RefPanelConfig(use_hm3_ref_panel_snps=True)`
  and `LDScoreConfig(use_hm3_regression_snps=True)` instead of wiring a custom
  HM3 path.

The old `--regression-snps` and `--print-snps` behavior is unified under
`--regression-snps-file`. LD-score outputs are fixed files under `output_dir`;
legacy `.l2.*` and `.w.l2.*` filenames are not emitted by the public writer.
Existing owned LD-score family artifacts are refused by default and require
`--overwrite` or `overwrite=True`, so reruns cannot silently replace artifacts
produced under a different configuration snapshot. With overwrite enabled, a
successful run removes stale owned siblings that the current run did not
produce, such as `ldscore.query.parquet` after switching from query LD scores
to baseline-only LD scores.

The same coherent-family rule applies to `munge-sumstats`, `build-ref-panel`,
`partitioned-h2`, `rg`, and `annotate`: no-overwrite mode rejects any owned
sibling in the current public layout, while successful overwrites delete stale
current-contract owned siblings and preserve unrelated files. Removed legacy
root diagnostic names are ignored by preflight and cleanup. Sharded workflows
may narrow the owned package to the current shard; for `build-ref-panel`,
concrete chromosome PLINK prefixes own only that chromosome's package, while
`@` chromosome-suite prefixes own the full all-chromosome panel package.

---

## Caveats and Known Limitations

**Reference-panel runtime state does not carry its own full `GlobalConfig`.**
`RefPanelConfig` owns backend and filtering options such as `r2_dir`,
`maf_min`, `ref_panel_snps_file`, and `use_hm3_ref_panel_snps`. The shared identifier and genome-build
assumptions come from the `GlobalConfig` passed to `RefPanelLoader` and then
captured by the LD-score workflow result snapshots.

**`AnnotationBundle` does not carry chromosome-level genome-build metadata.**
Genome build for annotations is inferred from the reference metadata embedded in
`AnnotationBundle.metadata`. Compatibility validation trusts the `GlobalConfig`
snapshot rather than re-inferring the build from coordinates.

**Workflow-specific SNP controls are not part of `config_snapshot`.**
`config_snapshot` records shared assumptions such as `genome_build` and
`snp_identifier`. Per-run LD-score controls such as
`RefPanelConfig.ref_panel_snps_file`, `RefPanelConfig.use_hm3_ref_panel_snps`,
`LDScoreConfig.regression_snps_file`, and
`LDScoreConfig.use_hm3_regression_snps` still materially affect the outputs, but
callers who need to preserve that provenance should persist the
`RefPanelConfig` / `LDScoreConfig` they used alongside the written artifacts.

**`load_sumstats()` recovers provenance for current disk artifacts.**
Curated `sumstats.parquet` and `.sumstats(.gz)` artifacts written by the current
munger have a neighboring root `metadata.json` sidecar that records a thin
metadata payload: `schema_version`, `artifact_type`, `snp_identifier`,
`genome_build`, and optional `trait_name`. Detailed coordinate provenance, liftover
reports, HM3 provenance, output bookkeeping, and row counts are written to
`diagnostics/sumstats.log`; row-level liftover drops are written to
`diagnostics/dropped_snps/dropped.tsv.gz`. Neither belongs in the metadata sidecar. The
loader reconstructs `config_snapshot` from the sidecar identity fields.
Older package-written artifacts without the current sidecar are rejected with a
regeneration message; sidecars missing the current identity provenance are
treated as invalid metadata rather than a migrated older format.

**LD-score directories must carry current identity provenance.**
Canonical LD-score directories written by the current workflow include root
metadata with a config snapshot and minimal identity provenance. Older or
malformed package-written metadata is rejected with a regeneration message instead of
loading with inferred provenance.

**Package-written artifacts do not use missing snapshots as a compatibility
bypass.**
Objects constructed manually with missing provenance are not a promise of disk
artifact compatibility. Package-written sumstats and LD-score artifacts that lack
the current schema/provenance contract must be regenerated before they can be
combined in downstream workflows.

**Notebook re-execution order still matters for the registry.**
The global registry is process-wide. If a user calls `set_global_config()` in a cell,
re-runs an earlier cell that did not capture a config explicitly, and then combines the
results, they will get a `ConfigMismatchError`. This is the intended behavior — the
error message explains the mismatch clearly. It is not a silent corruption.

**The `run_ldscore()` and `run_build_ref_panel()` convenience wrappers read the global
registry at call time.** Explicit `global_config=` arguments are always preferred.
