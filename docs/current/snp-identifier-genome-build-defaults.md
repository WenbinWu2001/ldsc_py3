# SNP Identifier / Genome Build: Default Design

Last updated on: 2026-10-26

## Standalone build inspection

```bash
ldsc infer-genome-build PATH
```

Supply one exact file or recognized LDSC artifact directory; `--help` is the only option.
The helper supports standalone BIM files without BED/FAM companions, headered text tables (including gzip), Parquet, and recognized LDSC metadata JSON.
It recognizes annotation, munged-sumstats, LD-score, R²-panel (root or build child), gene-index, and saved-result layouts.
It reads their own metadata and relevant members without following upstream input paths or recursively searching unrelated files.
It does not accept batches, globs, chromosome-suite expansion, build overrides, or column overrides.
No new Python API is exported.

Concrete metadata or build-specific position columns establish a declared build, using the existing build aliases.
Otherwise, recognizable SNP coordinates are sampled from every relevant file/chromosome, using at most 5,000 usable rows per file/chromosome and the existing HM3 criteria of at least 200 informative matches and at least 99% support.
The inferred output includes the input coordinate basis; no coordinates are changed.
Declared metadata is not independently verified, and inference does not examine every coordinate row for agreement.
The helper reports the artifact's own coordinates rather than substituting an upstream source build.
For example, a panel's emitted builds describe its outputs; `source_genome_build` does not replace them.

```text
Genome build: hg19 (inferred from SNP coordinates; 1-based).
This artifact’s workflow uses rsIDs for SNP matching.
```

The second sentence appears only when applicable metadata records `rsid` or `rsid_allele_aware`.
An rsID column alone does not establish the workflow's matching mode.
The sentence may also accompany an unknown result.

Insufficient or inconclusive evidence produces `Genome build: unknown (...)` with a short reason.
This includes rsID-only inputs and ordinary BED intervals without build metadata; interval endpoints are not SNP positions.
An unresolved chromosome or member prevents a single inferred build for the complete artifact.
Conflicting declarations or confident chromosome/member inferences for the same coordinate set cause an error.
Explicit alternatives, such as separate hg19/hg38 position columns or R² build directories, receive separate labeled results ordered hg19 then hg38.
Independent alternatives are not treated as conflicts.

| Exit code | Meaning |
| --- | --- |
| 0 | Inspection completed, including `unknown` or deliberately available multiple builds. |
| 1 | Missing/unreadable input, malformed or unsupported format, unrecognized directory, or conflicting build evidence. |
| 2 | Invalid command syntax or an unexpected internal failure. |

Reports go to stdout and errors go to stderr.
No artifact files, directories, logs, or failure markers are created, including when invalid options are supplied.
This informational helper does not weaken the workflow consistency rules below.
See [inspection troubleshooting](../troubleshooting.md#infer-genome-build-input-cannot-be-inspected) for input failures.

## Genome-build policy

This contract governs inference, checking, and resolution before coordinate-dependent computation.
The shared collector is `GenomeBuildCheck` in [_genome_build_policy.py](../../src/ldsc/_genome_build_policy.py); statistical evidence remains in `genome_build_inference.py`.
No new flags or per-input build declarations are required by this policy.

### One build for inputs used together

Inputs whose coordinates are used together must share one workflow genome build.
A concrete build supplied through the command's existing option must agree with every applicable input declaration after equivalent names are normalized, including GRCh37/hg19 and GRCh38/hg38.
An existing metadata label is a declaration regardless of whether an earlier workflow inferred it or a user supplied it.
Unlabeled inputs are admissible under an explicit build assertion; fresh statistical inference cannot override that assertion.
For example, `--genome-build hg19` declares hg19 for the applicable workflow, while `--source-genome-build hg19` declares the source build for panel building or munging.
An explicit assertion establishes how unlabeled coordinates are interpreted; it never overwrites an input label, converts positions, or replaces an existing liftover operation.
Source and target builds in a workflow that performs liftover remain separate coordinate stages, so an intentional source-to-target conversion is not itself a build conflict.
Errors must suggest only options supported by the affected command; workflows without a build override direct the user to matching upstream artifacts.

Without a concrete assertion, inspect build evidence across all relevant inputs and their selected shards before deciding whether the workflow can proceed.
No input wins a disagreement because it was listed or read first.
Conflicting hg19 and hg38 evidence causes an error before substantial computation or publication of scientific outputs.
If neither the inputs nor a concrete assertion establish a workflow build, resolution fails instead of choosing a default build.
Path syntax does not determine the build: exact paths, globs, and `@` declarations use the same build policy for the files they select.
Existing chromosome-coverage requirements remain separate from build consistency.

### Build declarations and bounded inference

Valid artifact build metadata or a build-specific position column establishes that input's build without coordinate re-inference.
Compare the declared builds across inputs using the same consistency rule as inferred builds; a declaration on one input does not override conflicting evidence from another input.
A concrete user-supplied build bypasses fresh coordinate inference while retaining declaration checks.
A disagreement stops the run before substantial computation or scientific output publication.
Repair the inputs by supplying matching coordinates, correcting erroneous metadata only after verification, or performing an explicit coordinate conversion.

For inputs that require coordinate inference, sample at most 5,000 rows from every relevant file/chromosome.
Coordinate columns are streamed to discover chromosomes; inference uses only those bounded samples.
The current criteria require at least 200 uniquely informative packaged HM3 matches and at least 99% support for the winning build/coordinate-basis hypothesis.
See `MIN_INSPECTED_REFERENCE_SNPS` and `DOMINANCE_THRESHOLD` in [genome_build_inference.py](../../src/ldsc/genome_build_inference.py).
Successful inference for one input must not stop evidence collection for the remaining inputs or shards.
An exhaustive build check of every coordinate row is not performed.
Mixtures confined to unsampled rows can escape detection.

### SNP keep-list rules

Build validation depends on the SNP identifiers actually used for filtering, not merely on which columns happen to be present in the keep-list.
Apply these rules consistently to user-supplied reference-panel and regression SNP restrictions and other keep-list readers:

| Keep-list type and filtering identity | Genome-build validation |
| --- | --- |
| rsIDs only, with filtering by rsID | Not required; no coordinate build is used by this restriction. |
| Explicit build-specific positions, such as `hg19_POS` | Select the column matching the workflow's build; its declaration supplies the build without coordinate inference. |
| Generic `CHR/POS`, with coordinate filtering and no build declaration | Under `auto`, attempt inference and raise an error if evidence is insufficient. |

Existing valid artifact build metadata remains sufficient evidence under the shared declaration rule above.
When both generic `POS` and a position column matching the resolved build are present, prefer the matching build-specific column.
For example, an hg19 workflow uses `hg19_POS` even when `POS` is also present.
Use generic `POS` only when the matching build-specific column is absent.
An unused generic column does not create a genome-build conflict.
The presence of both hg19 and hg38 coordinate columns represents available alternatives; select the column for the resolved workflow build rather than treating unused alternatives as conflicting inputs.

A generic coordinate keep-list with insufficient inference evidence remains `unresolved` and causes an automatic run to fail, even when the other inputs agree on one build.
It does not inherit the build from those other inputs.
Finish examining the remaining inputs and include the keep-list's flag, path, and inference failure reason in the grouped error and workflow log.
The actionable message asks the user to check the keep-list's coordinate build and explicitly specify `--genome-build hg19` or `--genome-build hg38` for commands supporting that option.
Use `--source-genome-build` for source-build workflows instead; do not suggest an unsupported flag.
A concrete assertion permits an unlabeled generic coordinate keep-list without successful coordinate-based inference, while ordinary input-format and identifier validation still apply.
See the [insufficient-evidence example](../troubleshooting.md#common-keep-list-build-could-not-be-resolved).

### Ordinary BED files inherit the workflow build

Ordinary BED interval files are not independently checked for genome build.
Report them as `unchecked`, then interpret their intervals in the workflow's resolved or explicitly declared build.
This is an accepted assumption about compatibility, not an inferred or verified build for the BED file.
An ordinary BED file does not supply build evidence, and its unchecked status alone does not prevent a run whose build is otherwise resolved.
BED's zero-based, half-open interval convention remains unchanged.
Unlike an unresolved generic coordinate keep-list, an ordinary BED file is deliberately unchecked and does not need inference evidence.

### Complete, grouped validation errors

Finish examining every relevant input before raising a genome-build consistency error.
Record a missing, unreadable, or otherwise unusable input and continue examining the other inputs that can be inspected; do not stop at the first mismatch or failed inference.
The report must account for every relevant input, including inputs that agree and ordinary BED files that are deliberately unchecked.
Group the findings by `hg19`, `hg38`, `unresolved`, and `unchecked`.
These are user-facing report labels, not a persisted artifact status enum.

| Report group | Meaning |
| --- | --- |
| `hg19` or `hg38` | Available evidence identifies that build; distinguish coordinate inference from a metadata declaration. |
| `unresolved` | The input could not be assigned a build; include the reason, such as insufficient informative coordinates or an unreadable file. |
| `unchecked` | Build checking is intentionally omitted; ordinary BED files inherit the resolved workflow build, while rsID-only restrictions used for rsID filtering need no coordinate build. |

List the CLI input flag and the concrete file or shard under each group, with chromosome and evidence source where useful.
When members supplied by one flag disagree, show that flag in each affected group and identify the individual members.
The same flags identify input roles in reports from Python workflows; Python callers set the corresponding configuration fields.
Use a group heading such as "Build evidence across inputs" because not every known build was inferred from coordinates.
State the observed conflict, identify unresolved evidence without pretending it was checked successfully, and explain what the user can do next.
The actionable guidance asks the user to check the listed inputs for a consistent build and gives a command-appropriate explicit assertion example for verified unlabeled inputs when available.
The message must identify the affected input, declared and requested or conflicting builds, and actionable remedies.
The grouped findings belong in the run-aborting message and workflow log, so users do not need `DEBUG` or repeated failed runs to discover additional input conflicts.
See the [grouped report example](../troubleshooting.md#common-genome-build-consistency-check-failed).

### Command scope

| Command | Evidence and interpretation |
| --- | --- |
| `annotate` | Selected annotation/BIM grids and the supplied gene-catalog declaration; ordinary BED inherits the resolved projection build. |
| Direct `ldscore` | Baseline and prebuilt query annotations, PLINK BIM shards or R2 Parquet and sidecar declarations, the gene catalog when used, and coordinate keep lists. |
| `build-r2-panel` | All source BIM shards and the source-coordinate keep list, before filtering; genetic maps are selected by the resolved source/emitted build. cM map validation waits until source inference. |
| `munge-sumstats` | Raw source coordinates and the source-coordinate keep list, before QC; explicit `--output-genome-build` remains a separate target. `--infer-only` checks the same source-build declarations and reports conflicts without creating output files. |
| `convert-ldsc2-ldscores` | Reference declarations or per-chromosome inference establish coordinate identity; weight/annotation declarations must agree. Other unlabeled legacy coordinates remain descriptive because rsID joins transfer reference coordinates to the output. rsID output requires no build inference. |
| `h2`, `partitioned-h2`, `rg` | Coordinate-family saved artifacts must agree; all RG traits participate before fitting. No CLI build override. Footerless Parquet remains invalid; legacy LDSC2 text uses panel-owned coordinate projection. |
| `quantile-h2` | Saved LD-score build and resupplied coordinate artifacts must agree; BED inherits that build. No CLI build override. |
| `build-gene-ldscore-index` | Required explicit hg19 or hg38 must agree with the catalog and applicable input declarations; no auto or liftover. Indexed scoring later inherits the immutable index and rejects live configuration overrides. |
| `query-r2` | Existing concrete build selects/interprets a panel; otherwise coordinate-panel shard metadata must agree. Pair endpoints are lookup requests interpreted in the selected panel's build, with no separate pair-build inference. |
| `convert-h2-scale`, `plot` | Saved-result consumers; no coordinate-build inference. |

Build metadata describes emitted coordinates.
After munging liftover it is the output build; each R2 output records its own emitted build.
Neither should be read as an assertion about the original pre-liftover coordinates.
SNP identity and file-format validation remain separate: an explicit build does not permit invalid alleles, missing required columns, broken sidecar bindings, or incompatible identity modes.

## Design Decisions

### The invariant pair

`GlobalConfig` has two fields that form a constraint pair:

| Field | Purpose |
|---|---|
| `snp_identifier` | one of exactly `"rsid"`, `"rsid_allele_aware"`, `"chr_pos"`, `"chr_pos_allele_aware"` — how SNPs are identified |
| `genome_build` | `"hg19"`, `"hg38"`, `"auto"`, or `None` — coordinate system for CHR:POS |

`GlobalConfig.genome_build` is only meaningful in the `chr_pos` family. In the
`rsid` family it is irrelevant and is silently coerced to `None`. A gene-list
LD-score run is the narrow exception at the workflow level: it separately
selects hg19 or hg38 for interval projection without changing rsID identity
metadata.

`build-gene-ldscore-index` is a second narrow workflow-level exception. New
construction requires an explicit choice among all four identity modes and an explicit hg19 or hg38 assertion in its own `GeneLDScoreIndexBuildConfig`.
Exactly one cM, kb, or SNP-count window is required, without a default.
In rsID-family modes, ordinary `GlobalConfig` still carries `genome_build=None`; the index independently binds its projection build because gene intervals, region masks, maps, and PLINK positions use that coordinate system.

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

Multi-allelic base-key clusters are dropped because allele-aware modes exist for safer merging, not for retaining multi-allelic sites.
Computing LD scores over retained same-position clusters would need its own design for same-position pair semantics, LD-window traversal, and R2 endpoint indexing.

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
artifact. `munge-sumstats` is the exception: its raw source build defaults to
`--source-genome-build auto`, and coordinate-family runs require the user to
choose the final build with `--output-genome-build`.
This is intentional: CLI users benefit from the choice being visible in their
shell history and scripts. The API and CLI have different ergonomic defaults
because they serve different user models.

| Workflow | Identifier default | Genome-build default | Notes |
|---|---|---|---|
| `annotate` | `chr_pos_allele_aware` | `None` (CLI requires explicit for coordinate-family modes) | Raises unless `--genome-build` is supplied or inferable by that workflow. |
| `munge-sumstats` | `chr_pos_allele_aware` | source: `auto`; output: required for coordinate-family modes | Pass `--output-genome-build hg19` or `--output-genome-build hg38`; the raw source build is inferred unless `--source-genome-build hg19/hg38` is supplied. rsid-family modes reject source/output/liftover build flags and store `genome_build=None`. Requires usable `A1/A2`; rerun with `--snp-identifier chr_pos` or `--snp-identifier rsid` to run without allele-aware identity. The removed `--no-alleles` flag is not accepted. |
| `ldscore` | `chr_pos_allele_aware` | `None`; gene-list route treats omission as `auto` | Coordinate modes require a concrete/inferred build. Gene-list projection also requires a build in rsID modes, but shared rsID artifact metadata remains `genome_build=None`. |
| `build-r2-panel` | registry or `--snp-identifier` | **ignored** | Uses `--source-genome-build` (separate field); `GlobalConfig.genome_build` is never consulted. |
| `h2`, `partitioned-h2`, `rg` | registry at construction | registry at construction | Current Parquet uses footer provenance. Legacy LDSC2 `.sumstats[.gz]` is treated as rsID lookup input and projected onto the LD-score panel's identity; footerless Parquet is rejected. |
| `convert-ldsc2-ldscores` | `rsid` | `auto` evidence | Output identity is restricted to allele-unaware `rsid` or `chr_pos`. rsID output remains buildless; chr_pos output requires inferred or explicit reference build. |
| Python `run_ldscore()` wrapper | registry | registry (`"auto"`) | Inherits `chr_pos_allele_aware + auto` from the registry; `auto` is resolved to `hg19`/`hg38` during inference. |

### `build-r2-panel` isolation

`build-r2-panel` is intentionally isolated from `GlobalConfig.genome_build`.
The PLINK source build is owned by `ReferencePanelBuildConfig.source_genome_build`
and defaults to `"auto"`, which checks every `.bim` shard together with coordinate-based SNP restrictions. This keeps the panel-building
workflow self-contained: it neither reads nor writes the global genome-build
assumption. The `GlobalConfig` that `build-r2-panel` constructs internally holds
`genome_build="auto"` only to satisfy the construction invariant; that value is
never consumed by the workflow.

---

## Regression Provenance Hierarchy And Downgrade

Regression workflows (`h2`, `partitioned-h2`, `rg`) can see three independent
`GlobalConfig` sources for any given run:

1. The LD-score root metadata snapshot — `LDScoreResult.config_snapshot`
2. The reconstructed sumstats artifact snapshot — `SumstatsTable.config_snapshot`
3. The runner's live `self.global_config` — from the process registry at construction time

Coordinate-family regression checks the LD-score artifact and all participating traits before fitting.
Known artifact snapshots establish builds; suitable in-memory objects without a concrete build use bounded coordinate inference.
Conflicting or unresolved builds produce one grouped error.
The runner's live configuration does not override saved artifact provenance.

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
| `ldscore_source.py` | LD-score metadata loading | Current root metadata must carry minimal identity provenance; old package-written artifacts should be regenerated |
| `regression_runner.py` | `build_dataset` empty-merge error | Append `f"Active config: {self.global_config!r}."` |

Call sites in `ldscore_calculator.py` already resolve `genome_build` to a
concrete value or `"auto"` before constructing `GlobalConfig`. They are not
affected.
