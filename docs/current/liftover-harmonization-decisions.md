# Liftover Harmonization Decisions

Date: 2026-05-10
Branch: `codex/liftover-sumstats-munger`

> **STATUS — 2026-05-11: implemented harmonization contract.** This document
> records the liftover contracts for the harmonized workflows. The
> implementation plan is
> [`docs/superpowers/plans/2026-05-10-liftover-harmonization.md`](../superpowers/plans/2026-05-10-liftover-harmonization.md).

This note records the agreed liftover contracts after harmonizing summary
statistics munging and PLINK reference-panel building. It is intended as the
handoff document for future liftover harmonization work.

## Scope

The implementation harmonizes mechanics, not workflows. The shared layer lives
inside `src/ldsc/_kernel/liftover.py` and is still internal. It is not a public
Python API. Public workflow modules keep their own configuration, validation,
output, and metadata contracts.

## Shared Mechanics

- Shared chain-file mapping uses 1-based input and output coordinates.
- When `pyliftover` returns multiple hits, the first same-chromosome hit wins.
- If a position has both cross-chromosome and same-chromosome hits, the
  same-chromosome hit is retained.
- Positions with no same-chromosome hit are dropped and counted separately as
  unmapped or cross-chromosome drops when the reason is knowable.
- Shared duplicate-coordinate helpers implement the `drop-all` primitive:
  every row in a duplicate coordinate group is dropped. In other words, if multiple SNPs are mapped to the same coordinate in target build, all of them are dropped.
- Shared drop reports carry readable counts plus up to five examples with core
  fields: `SNP`, `CHR`, `source_POS`, optional `target_POS`, and `reason`.
  Normal log levels are count-only; examples are emitted only at `DEBUG`.
- Missing chain-file and missing optional dependency errors include
  workflow-aware labels and flag hints.

## Summary-Statistics Munger Contract

- Liftover is meaningful only in `chr_pos` mode.
- `SNP` is always a label field in `chr_pos` mode. Sumstats liftover updates
  only `CHR` and `POS`; it never rewrites `SNP`.
- Sumstats liftover requires explicit source/target/method validation:
  `--target-genome-build` plus exactly one of `--liftover-chain-file` or
  `--use-hm3-quick-liftover` when source and target differ. HM3 quick liftover
  additionally requires `--use-hm3-snps`, so HM3 filtering is explicit.
- Missing `CHR`/`POS` rows are dropped before mapping and logged.
- Source duplicate `CHR/POS` groups are dropped before mapping.
- Target duplicate `CHR/POS` groups are dropped after mapping.
- Duplicate dropping happens only when a liftover request reaches the map
  stage; ordinary no-liftover `chr_pos` munging does not drop duplicate
  coordinates just because they are duplicated.
- If missing-coordinate, unmapped, cross-chromosome, source-duplicate, or
  target-duplicate filtering removes every row, sumstats munging raises a hard
  error.
- `dropped_snps/dropped.tsv.gz` is always written under the sumstats output
  directory. Clean runs produce a header-only gzip TSV. It uses columns `CHR`,
  `SNP`, `source_pos`, `target_pos`, `reason` with nullable `string`/`Int64`
  dtypes and may contain all five reasons: `missing_coordinate`,
  `source_duplicate`, `unmapped_liftover`, `cross_chromosome_liftover`, and
  `target_collision`.
- Newly written `sumstats.metadata.json` sidecars are intentionally thin:
  `format`, optional `trait_name`, and `config_snapshot`.
- Detailed coordinate provenance, liftover counts, HM3 map provenance, output
  paths, and row-group details are readable `sumstats.log` entries rather than
  sidecar payloads.
- Old pre-`config_snapshot` sidecars are not supported. Completely missing
  sidecars remain legacy unknown-provenance inputs for existing loader
  compatibility.
- HM3 quick liftover is a packaged-map coordinate shortcut, not a general
  chain-liftover replacement.

## Reference-Panel Builder Contract

- PLINK reference-panel building keeps its existing UX: the source build is
  explicit or inferred, and a matching chain file emits the opposite build.
- `--use-hm3-snps` restricts the retained reference-panel universe to the
  packaged curated HM3 map. `--use-hm3-quick-liftover` requires that restriction
  and emits the opposite build using the packaged map.
- No new public `target_genome_build` field was added for reference-panel
  building.
- Matching chain-file liftover and HM3 quick liftover are invalid when the active
  `GlobalConfig.snp_identifier` is `rsid`. In rsID mode, omit the matching
  chain or quick-liftover flag and build source-genome coordinates only.
- Coordinate duplicate groups are dropped with `drop-all`. There is no public
  duplicate-position policy knob.
- Duplicate-position handling applies only in `chr_pos` mode.
- In `rsid` source-only builds, coordinate duplicate filtering is skipped and
  logged once per run.
- Source duplicate groups are dropped before liftover mapping.
- Target duplicate groups are dropped after target-build positions are known.
- If all SNPs are dropped for one chromosome, that chromosome is skipped.
  The run fails only if no chromosome emits artifacts.
- Existing runtime metadata TSV and parquet schemas remain stable.
- Liftover/drop provenance goes to workflow logs and
  `dropped_snps/chr*_dropped.tsv.gz` sidecars.
- Per-chromosome `dropped_snps` sidecars are always written for every
  chromosome processed by the run. Clean chromosomes produce header-only gzip
  TSVs. Ref-panel sidecars use the shared columns `CHR`, `SNP`, `source_pos`,
  `target_pos`, `reason` and cover the four ref-panel-applicable reasons:
  `source_duplicate`, `unmapped_liftover`, `cross_chromosome_liftover`, and
  `target_collision`. `missing_coordinate` is not produced by ref-panel because
  PLINK BIM `BP` is structurally non-null.
- Output preflight includes the always-written dropped-SNP sidecar paths so
  overwrite behavior stays deterministic.
- Before chromosome processing starts, the builder warns once when existing
  class-2 parquet or metadata artifacts under `hg19/` or `hg38/` are outside
  the current run's expected output set. The warning is informational; it does
  not delete files or abort the run.

## GlobalConfig And Identity

- `GlobalConfig.snp_identifier='chr_pos'` means coordinate identity is `CHR/POS`
  under a concrete genome build.
- `GlobalConfig.snp_identifier='rsid'` means row identity is `SNP`; genome build
  is ignored and normalized to `None`.
- Liftover is coordinate behavior and therefore valid only when the relevant
  workflow is operating in `chr_pos` mode.
- Build aliases continue to normalize through the existing rules:
  `hg37`/`GRCh37` to `hg19`, `GRCh38` to `hg38`, and `auto` only where a
  workflow supports inference.

## Backward Compatibility Rules

- Sumstats sidecars written by the current workflow must include
  `config_snapshot`.
- Old package-written sumstats metadata sidecars without current identity
  provenance do not need migration support.
- Missing sumstats sidecars on package-written artifacts are rejected; those
  artifacts must be regenerated with the current LDSC package.
- Reference-panel runtime metadata TSV and R2 parquet schemas must not change
  for this harmonization.
- Dropped-SNP sidecar column names are stable: `CHR`, `SNP`, `source_pos`,
  `target_pos`, `reason`. Consumers should read with explicit nullable
  `string`/`Int64` dtypes because empty/header-only files and nullable
  `target_pos` cells cannot be represented safely with plain NumPy integer
  dtypes.
- Existing multi-build reference-panel loading ambiguity rules remain for old
  or external panels.

## Open Follow-Up Areas

- Decide whether the internal liftover helper should become a fuller
  package-level internal service object with workflow adapters, or remain a
  function/dataclass helper module.
- Consider a common typed drop-report object that can be surfaced through
  Python results without expanding disk sidecars.
- HM3 quick liftover remains an HM3-restricted coordinate shortcut. Revisit only
  with a separate design if broader non-HM3 liftover behavior is needed.
- Audit whether downstream loaders should expose liftover provenance from logs
  or leave logs as human audit artifacts only.
- Consider adding a small compatibility validator for old external reference
  panels whose metadata predates the current build/identifier contracts.

## Current Maintenance Checklist

When touching liftover code, preserve these contracts unless a new design
explicitly changes them: liftover is `chr_pos`-only, `SNP` is a label in
`chr_pos`, sumstats liftover updates only `CHR/POS`, reference-panel liftover is
rejected in `rsid`, HM3 quick liftover requires HM3 SNP restriction, dropped-SNP
sidecars stay schema-stable, and detailed provenance belongs in workflow logs.
