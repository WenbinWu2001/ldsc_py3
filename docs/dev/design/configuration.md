# Config Design: Immutable Config + Provenance-Carrying Results

Last updated on: 2026-11-06

`--workers` sets the number of managed workers (default one, each with one BLAS thread); CPU budgets, limits and failure handling are in [parallel execution](../../reference/parallel-execution.md).

For `partitioned-h2`, `--workers` belongs to the batch workflow, not to `RegressionConfig`.

Examples below use the munged file names written without a trait label (`sumstats.parquet`, `sumstats.gz`); naming rules are in [artifact metadata](../../reference/formats/artifact-metadata.md).
Memory controls such as `query_batch_size` are separate from these statistical settings; see the [memory design](memory-and-cpu.md).

## The Problem This Design Solves

Genome build (`hg19`/`hg38`) and SNP identifier mode are global analysis
assumptions that every step of the pipeline depends on (the modes are defined in
[column schema](../../reference/formats/column-schema.md#snp-identifier-modes)). If these settings differ
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
from the `sumstats.parquet` footer written by the current munger. A legacy
`.sumstats` or `.sumstats.gz` input has no embedded provenance and is explicitly
projected by rsID onto the LD-score panel; footerless Parquet is rejected.

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

Regression resolves the effective SNP identifier mode from saved provenance first, and `--allow-identity-downgrade` permits only same-family allele-aware/base mixes; see [regression provenance hierarchy and downgrade](../../reference/snp-identity-and-genome-build.md#regression-provenance-hierarchy-and-downgrade).

### 4. The global registry is a convenience default only

`set_global_config()` and `get_global_config()` are retained for notebook UX: setting the config once at the top of a notebook is ergonomic.
But the registry is demoted: a workflow function that accepts an explicit config reads the registry only as a fallback when none is passed, and the snapshot captured at computation time — not the registry — is what travels with results.
The keyword wrappers that mirror one command are the exception: `run_ldscore()`, `run_build_r2_matrix()`, `run_h2()`, `run_partitioned_h2()`, `run_rg()` and `run_quantile_h2()` take no explicit config and read the fields they need from the registry at every call (see [Caveats and Known Limitations](#caveats-and-known-limitations)).

Calling `set_global_config()` mid-analysis will not corrupt previously computed
results, because those results carry their own frozen snapshots.

---

## Critical vs. Advisory Fields

| Field | Class | Mismatch behavior |
| --- | --- | --- |
| `genome_build` | Critical | Hard error (`ConfigMismatchError`) |
| `snp_identifier` | Critical | Hard error (`ConfigMismatchError`) |
| `log_level` | Advisory | Ignored during compatibility checks |

The compatibility helper compares only `GlobalConfig` snapshots, so the table
above lists only `GlobalConfig` fields. Workflow-specific LD-score controls such
as `RefPanelConfig.ref_panel_snps_file` and `LDScoreConfig.regr_snps_file`
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

## Public API Surface

| Name | Location | Purpose |
| --- | --- | --- |
| `ConfigMismatchError` | `ldsc.config` | Raised on incompatible config combination |
| `validate_config_compatibility(a, b)` | `ldsc.config` | Compatibility helper used by workflow merge points; not exported from `ldsc` |
| `config_snapshot` field | Result classes | Frozen config from computation time, or `None` when provenance is genuinely unknown |

---

## LD window

`LDScoreConfig`, `R2MatrixBuildConfig` and `GeneLDScoreIndexBuildConfig` apply one LD-window rule when they are constructed (`_ld_window()` in `src/ldsc/config.py`).
Exactly one of `ld_wind_snps`, `ld_wind_kb` and `ld_wind_cm` is set, to a finite number greater than 0 that is not a boolean, and `ld_wind_snps` is an integer.
The value is stored as `int` for SNP counts and `float` for kilobases and centiMorgans, so one window has one spelling in R² matrix footers, gene-index identities and LD-score `metadata.json`.
Each of the three configs returns it as an `LDWindow` (`unit`, `value`) from its `ld_window` property; `LDWindow` is internal and not exported from `ldsc`.
The window helpers in `src/ldsc/_refpanel/windows.py`, the R² matrix footer writer and reader, and the gene-index identity take that value.
A violation raises `LDSCConfigError` naming both the command-line option and the Python field.
On the command line, `ldscore` and `build-r2-matrix` reject two windows while parsing (exit 2), and `build-r2-matrix` also requires one; `ldscore` cannot, because indexed mode takes its window from the index.

---

## Workflow-Specific SNP Controls

Standard LDSC still needs two distinct SNP-set controls, but they no longer live
on `GlobalConfig`. LDSC3 assigns them to the workflow objects
that actually own those decisions.

### The hierarchy

```text
Annotation bundle rows B                      (AnnotationBuilder)
  intersects prepared reference panel A'     (RefPanelConfig ref-panel filters)
    then optional regression SNP subset C    (LDScoreConfig regression SNP controls)
```

| Concept | Owner | CLI flag | When applied | Artifact effect |
| --- | --- | --- | --- | --- |
| Reference-panel SNP restriction | `RefPanelConfig.ref_panel_snps_file` | `--ref-panel-snps-file` | `RefPanel.prepare_chromosome()` restricts the panel and aligns `B_chrom` before the kernel call | Deliberately shrinks the reference SNPs: after all panel filtering and intersection they are the reference SNPs $B \cap A'$, the annotation rows $B$ that are also in the filtered panel $A'$; affects LD scores and count records |
| Reference-panel MAF/sample filters | `RefPanelConfig.maf_min`, `RefPanelConfig.keep_indivs_file` | `--maf-min`, `--keep-indivs-file` | `RefPanel.prepare_chromosome()` applies configured sample selection before PLINK genotype QC/MAF filtering, or filters parquet sidecar rows by MAF | Affects the prepared panel A' and LD computation with inclusive `MAF >= maf_min`; separate from `LDScoreConfig.common_maf_min` |
| PLINK genetic map (cM windows) | `RefPanelConfig.genetic_map_hg19_sources`, `genetic_map_hg38_sources` | `--genetic-map-hg19-sources`, `--genetic-map-hg38-sources` | `PlinkRefPanel.prepare_chromosome()` resolves the build/map and interpolates `CM` at retained `.bim` positions (map always wins). `RefPanelConfig(backend="r2_matrix")` rejects a genetic map when it is constructed, because the R² matrix sidecars hold `CM`. | Defines `CM` for `--ld-wind-cm` when the `.bim` `CM` is uninformative; no effect on SNP/kb windows |
| Reference-metadata export | `LDScoreConfig.export_ref_metadata` | `--export-ref-metadata` | After each chromosome's compute (PLINK worker), writes `ref_metadata/chrN_meta.tsv.gz`, which `LDScoreDirectoryWriter.write_batches()` publishes with the directory | Opt-in provenance/QC artifact (`CHR POS SNP A1 A2 CM MAF`, with the minor allele as `A1`); rejected for an R² matrix, whose sidecars hold the same metadata; accepted later by `quantile-h2 --ref-metadata-sources` |
| Regression SNP restriction | `LDScoreConfig.regr_snps_file` or bundled HM3 default | `--regr-snps-file` | After LD computation, when normalized/public rows are selected; `--regr-snps-exclude-regions` subtracts named intervals here | Written rows are the regression SNPs $B \cap A' \cap C$ outside the excluded regions; `regression_ld_scores` (`w_ld`) is computed over that same SNP set |
| Common-SNP-count threshold | `LDScoreConfig.common_maf_min` | `--common-maf-min` | During count-vector and overlap-matrix computation after LD scores are computed | Defines the common reference SNPs with inclusive `MAF >= common_maf_min` (deviates from legacy LDSC's strict `0.05 < FRQ < 0.95`); affects `common_reference_snp_count(s)` and the common-SNP overlap matrix, but not LD rows, LD scores, or the stored regression-SNP LD score |

### What each control does

**`RefPanelConfig.ref_panel_snps_file`** — the *reference SNPs*

- `AnnotationBuilder.run()` still builds annotations on the full annotation SNP grid `B`.
- `run_annotate()` and `ldsc annotate` do **not** apply this restriction.
- `RefPanel` applies the explicit restriction to the raw panel `A`, producing `A'`, alongside configured sample/MAF filters. Both metadata inspection and chromosome preparation honor it. Bundled HM3 is never a convenience way to shrink the reference SNPs.
- SNP restriction files are identity-only. Duplicate restriction keys collapse
  to one retained key, and non-identity columns such as `CM`, `MAF`, or other
  metadata are ignored rather than carried into LD-score or regression outputs.
- `RefPanel.prepare_chromosome()` aligns the chromosome-local annotation bundle to the filtered reference metadata so the kernel sees `B_chrom ∩ A'_chrom`.

When `None`, the workflow uses the full reference panel `A`.

**`LDScoreConfig.regr_snps_file`** — the *regression SNPs*

- The explicit restriction file or packaged HM3 map is loaded once into the
  `regression_snps` set `C`.
- LD scores and count totals are still computed over the reference SNPs $B \cap A'$.
- Only the persisted LD-score rows are reduced to the regression SNPs $B \cap A' \cap C$.
- There is no separate public weight artifact. The selected baseline rows keep
  `regression_ld_scores`, the historical `w_ld` LD score computed over the
  regression SNPs, in `ldscore.baseline.parquet`. Final `h2` and `rg`
  regression weights are computed later in the regression kernel from this
  value plus model-specific quantities.

When `None` on the public CLI, the normalized/public row table uses the bundled
HM3 map. Direct calculator calls may still provide no set explicitly.

The file format of `--regr-snps-file` and `--ref-panel-snps-file`, with an example per SNP identifier mode, is in [column schema](../../reference/formats/column-schema.md#snp-restriction-files).

### Typical configuration

```python
cfg = ldsc.GlobalConfig(genome_build="hg38", snp_identifier="chr_pos_allele_aware")
ldsc.set_global_config(cfg)

ref_panel = ldsc.RefPanelLoader(cfg).load(
    ldsc.RefPanelConfig(
        backend="r2_matrix",
        r2_matrix_dir="r2_ref_panel_1kg30x_1cM_hm3/hg38",
    )
)

ldscore = ldsc.LDScoreCalculator().run(
    annot,
    ref_panel,
    ldsc.LDScoreConfig(
        ld_wind_cm=1.0,
    ),
    global_config=cfg,
)
```

For R² matrices, R² bias and sample size are recorded only in the
Parquet schema (`ldsc:r2_bias` and `ldsc:n_samples`); there are no
`RefPanelConfig` bias or sample-size fields.
Readers accept only already-unbiased R² matrices written by the current `build-r2-matrix`; other R² matrices must be rebuilt.

### Regression-region exclusion

`--regr-snps-exclude-regions` belongs to the regression SNP controls: it subtracts named regions from the regression SNPs after `regr_snps_file` or HM3 selection and leaves the reference SNPs, annotation LD scores and counts untouched.
Its presets, command scope and build rule are in [region-exclusion presets](../../reference/formats/region-exclusion-presets.md).

### Artifact contract

- Materialized query `.annot.gz` files and in-memory `AnnotationBundle` objects stay on the annotation SNP grid `B`.
- Ordinary unpartitioned `run_ldscore()` calls may omit baseline/query annotation inputs; the workflow creates a synthetic all-ones `base` annotation after the reference panel has applied any explicit reference-panel restriction.
- Prebuilt, BED, and gene-list query annotations are valid only with explicit baseline annotations, so the synthetic `base` path is not used for partitioned heritability.
- BED and gene-list sources are query-local failure units. Only `ok`/`warning` queries enter scientific result objects; every source remains represented in diagnostic status records.
- Reference-panel SNP restrictions become visible only during LD-score calculation, when `ref_panel.prepare_chromosome()` aligns `B_chrom` to the filtered reference rows.
- Count records are accumulated over the reference SNPs $B \cap A'$ and stored in LD-score root `metadata.json`.
- Public `ldscore.baseline.parquet` and optional query batch file rows are the regression SNPs $B \cap A' \cap C$ outside the excluded regions; `w_ld` uses that same final set.
- LD-score `metadata.json` records `snp_universe_policy`, including independent reference and regression policies plus the regression-region removal count.
- LD-score `metadata.json` records `ld_window`: the run's window unit and value, and whether windows spanning a whole chromosome were allowed.

---

## Caveats and Known Limitations

**Reference-panel runtime state does not carry its own full `GlobalConfig`.**
`RefPanelConfig` owns backend and filtering options such as `r2_matrix_dir`,
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
`RefPanelConfig.ref_panel_snps_file` and `LDScoreConfig.regr_snps_file`
still materially affect the outputs, but
callers who need to preserve that provenance should persist the
`RefPanelConfig` / `LDScoreConfig` they used alongside the written artifacts.

**`load_sumstats()` recovers provenance for current disk artifacts.**
It rebuilds `config_snapshot` from the identity footer of munged `sumstats.parquet` and rejects missing or invalid footer provenance; `.sumstats.gz` carries none.
The footer keys and what the munger writes elsewhere are in [embedded identity metadata](../../reference/munge-sumstats.md#embedded-identity-metadata).

**LD-score directories must carry current identity provenance.**
LD-score directories written by the current workflow include root
metadata with a config snapshot and minimal identity provenance. Older or
malformed package-written metadata is rejected with a regeneration message instead of
loading with inferred provenance.

**Package-written artifacts do not use missing snapshots as a compatibility
bypass.**
Objects constructed manually with missing provenance are not a promise of disk
artifact compatibility. Package-written munged summary statistics and LD-score directories that lack
the current schema/provenance contract must be regenerated before they can be
combined in downstream workflows.

**Notebook re-execution order still matters for the registry.**
The global registry is process-wide. If a user calls `set_global_config()` in a cell,
re-runs an earlier cell that did not capture a config explicitly, and then combines the
results, they will get a `ConfigMismatchError`. This is the intended behavior — the
error message explains the mismatch clearly. It is not a silent corruption.

**Stored snapshots are read by named field.**
Loaders recover the required identity fields from a stored snapshot explicitly and ignore any other keys, rather than constructing a config from arbitrary snapshot keys.

**The `run_ldscore()` and `run_build_r2_matrix()` convenience wrappers read the global registry at call time.**
They take no `global_config=` argument.
`_runtime.workflow_kwargs.namespace_from_kwargs()` is the one place they read the registered SNP identifier mode, genome build and log level, so a parser default never stands in for them; `build-r2-matrix` records `genome_build="auto"` (`None` in rsID-family modes) and reads the PLINK build from `source_genome_build`.
