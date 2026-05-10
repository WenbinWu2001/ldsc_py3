# Liftover Harmonization Decisions

Date: 2026-05-10
Branch: `codex/liftover-sumstats-munger`

This note records the current liftover contracts after harmonizing summary
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
  every row in a duplicate coordinate group is dropped.
- Shared drop reports carry readable counts plus up to five examples with core
  fields: `SNP`, `CHR`, `source_POS`, optional `target_POS`, and `reason`.
- Missing chain-file and missing optional dependency errors include
  workflow-aware labels and flag hints.

## Summary-Statistics Munger Contract

- Liftover is meaningful only in `chr_pos` mode.
- `SNP` is always a label field in `chr_pos` mode. Sumstats liftover updates
  only `CHR` and `POS`; it never rewrites `SNP`.
- Sumstats liftover requires explicit source/target/method validation:
  `--target-genome-build` plus exactly one of `--liftover-chain-file` or
  `--use-hm3-quick-liftover` when source and target differ.
- Missing `CHR`/`POS` rows are dropped before mapping and logged.
- Source duplicate `CHR/POS` groups are dropped before mapping.
- Target duplicate `CHR/POS` groups are dropped after mapping.
- Duplicate dropping happens only when a liftover request reaches the map
  stage; ordinary no-liftover `chr_pos` munging does not drop duplicate
  coordinates just because they are duplicated.
- If missing-coordinate, unmapped, cross-chromosome, source-duplicate, or
  target-duplicate filtering removes every row, sumstats munging raises a hard
  error.
- Newly written `sumstats.metadata.json` sidecars are intentionally thin:
  `format`, optional `trait_name`, and `config_snapshot`.
- Detailed coordinate provenance, liftover counts, HM3 map provenance, output
  paths, and row-group details are readable `sumstats.log` entries rather than
  sidecar payloads.
- Old pre-`config_snapshot` sidecars are not supported. Completely missing
  sidecars remain legacy unknown-provenance inputs for existing loader
  compatibility.
- HM3 quick liftover remains sumstats-only unless a later design explicitly
  expands it.

## Reference-Panel Builder Contract

- PLINK reference-panel building keeps its existing UX: the source build is
  explicit or inferred, and a matching chain file emits the opposite build.
- No new public `target_genome_build` or liftover-method config field was added
  for reference-panel building.
- Matching chain-file liftover is invalid when the active
  `GlobalConfig.snp_identifier` is `rsid`. In rsID mode, omit the matching
  chain and build source-genome coordinates only.
- `duplicate_position_policy` now defaults to `drop-all`; `error` remains
  available.
- Duplicate-position handling applies only in `chr_pos` mode.
- In `rsid` source-only builds, duplicate-position policy is ignored and logged
  once per run.
- Source duplicate groups are dropped before liftover mapping.
- Target duplicate groups are dropped after target-build positions are known.
- If all SNPs are dropped for one chromosome, that chromosome is skipped.
  The run fails only if no chromosome emits artifacts.
- Existing runtime metadata TSV and parquet schemas remain stable.
- Liftover/drop provenance goes to workflow logs and existing
  `dropped_snps/chr*_dropped.tsv.gz` sidecars.
- `dropped_snps` sidecars contain duplicate drops only. Liftover unmapped and
  cross-chromosome drops are logged but not written to those sidecars.
- Output preflight still includes conditional duplicate sidecar paths so
  overwrite behavior stays deterministic.

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
- Old sumstats metadata sidecars without `config_snapshot` do not need
  backward compatibility.
- Missing sumstats sidecars remain tolerated as unknown-provenance legacy
  artifacts, because older disk outputs may have no sidecar at all.
- Reference-panel runtime metadata TSV and R2 parquet schemas must not change
  for this harmonization.
- Duplicate-drop sidecar schema remains `CHR`, `SNP`, `source_pos`,
  `target_pos`, `reason`.
- Existing multi-build reference-panel loading ambiguity rules remain for old
  or external panels.

## Open Follow-Up Areas

- Decide whether the internal liftover helper should become a fuller
  package-level internal service object with workflow adapters, or remain a
  function/dataclass helper module.
- Consider a common typed drop-report object that can be surfaced through
  Python results without expanding disk sidecars.
- Decide whether HM3 quick liftover should ever be allowed for reference-panel
  building. The current answer is no.
- Audit whether downstream loaders should expose liftover provenance from logs
  or leave logs as human audit artifacts only.
- Consider adding a small compatibility validator for old external reference
  panels whose metadata predates the current build/identifier contracts.

## Prompt For The Next Session

Work in:

`/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_Jerry/.worktrees/liftover-sumstats-munger`

Start by reading this decision document:

`docs/current/liftover-harmonization-decisions.md`

Then inspect:

- `src/ldsc/_kernel/liftover.py`
- `src/ldsc/sumstats_munger.py`
- `src/ldsc/ref_panel_builder.py`
- `src/ldsc/config.py`
- `docs/superpowers/specs/2026-05-09-munge-sumstats-liftover.md`
- `docs/superpowers/plans/2026-05-09-munge-sumstats-liftover.md`

Goal for a follow-up session: continue liftover harmonization without changing
the public workflow contracts unless explicitly approved. Preserve these
decisions: liftover is `chr_pos`-only, `SNP` is a label in `chr_pos`, sumstats
liftover updates only `CHR/POS`, reference-panel chain liftover is rejected in
`rsid`, HM3 quick liftover is sumstats-only, and detailed provenance belongs in
workflow logs while compatibility sidecars stay thin or schema-stable.
