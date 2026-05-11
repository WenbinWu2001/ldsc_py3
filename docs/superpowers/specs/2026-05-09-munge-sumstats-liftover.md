# Design - `munge-sumstats` Liftover Feature

## Context

Raw GWAS summary statistics may come in hg19 or hg38. Downstream LDSC steps
(`h2`, `partitioned-h2`, `rg`) need coordinate-keyed artifacts to use the same
genome build. The package already supports chain-file liftover for reference
panel construction and already rejects incompatible concrete genome builds at
regression time via `GlobalConfig` snapshots. The missing package-level feature
is coordinate liftover inside `munge-sumstats`.

This design adds opt-in liftover for `snp_identifier = chr_pos` only. In
`chr_pos` mode, SNP identity is always `CHR` and `POS` jointly. The `SNP` column
is treated as a label, usually an rsID, even when its values look like
coordinates. Liftover updates `CHR` and `POS`; it never rewrites `SNP`.

## Core Decisions

- Liftover is valid only when `snp_identifier = chr_pos`.
- Any liftover request in `rsid` mode errors before input IO.
- `SNP` is never the identity key in `chr_pos` mode and must not be rewritten by
  liftover.
- All `chr_pos` matching in the package uses `CHR/POS` keys, not `SNP`.
- Source build remains `GlobalConfig.genome_build`; do not add
  `MungeConfig.genome_build`.
- `--sumstats-snps-file` is interpreted in the source build because SNP
  filtering runs before liftover.
- Newly written sidecars stay thin: only schema, trait label, and
  `config_snapshot` are stored. Liftover and coordinate provenance are logged.
- Shared liftover mechanics live in an internal helper module, not a public API.
  Sumstats and reference-panel workflows share chain mapping, drop reports,
  duplicate-coordinate detection, and example formatting while retaining
  separate workflow contracts.
- Reference-panel PLINK builds keep the source-build plus optional opposite-build
  UX. Matching chain-file liftover is rejected in `rsid` mode; duplicate
  coordinate policy applies only in `chr_pos` mode and defaults to `drop-all`.

## CLI And Config

Add three `munge-sumstats` flags:

| Flag | Type | Meaning |
|---|---|---|
| `--target-genome-build` | `{hg19,hg37,GRCh37,hg38,GRCh38}` | Desired output build. Normalized to `hg19` or `hg38`; `auto` is invalid. |
| `--liftover-chain-file` | path | Explicit chain file for general liftover. |
| `--use-hm3-quick-liftover` | bool flag | Use the packaged HM3 curated coordinate map. |

Add only these fields to `MungeConfig`:

```python
target_genome_build: GenomeBuildInput | None = None
liftover_chain_file: Path | None = None
use_hm3_quick_liftover: bool = False
```

Validation rules:

- `liftover_chain_file` and `use_hm3_quick_liftover` are mutually exclusive.
- A method flag without `target_genome_build` is an error.
- If target build differs from resolved source build, one method is required.
- If target build equals resolved source build with no method, this is a no-op.
- If target build equals resolved source build with a method, this is an error.
- If source build is `auto` and cannot be inferred, liftover errors and asks the
  user to pass `--genome-build hg19` or `--genome-build hg38`.

Recommended rsID-mode error:

```text
munge-sumstats liftover requires snp_identifier='chr_pos'. In rsid mode,
regression merges by SNP ID and genome build is not a compatibility key.
Remove the liftover flags or run with --snp-identifier chr_pos.
```

## Pipeline

The munger stage order is:

```text
1. Existing parsing, QC, allele handling, and coordinate finalization
2. --sumstats-snps-file filtering in source-build coordinates
3. Optional liftover in chr_pos mode
4. Output writing in target-build coordinates
```

If liftover is requested, rows missing `CHR` or `POS` after existing QC and
filtering are dropped before mapping and counted in the run log. Missing
coordinates mean the SNP has no usable `chr_pos` identity for the map/match
stage; they are not a fatal input error. Malformed non-missing coordinates still
raise clear validation errors. Liftover may also drop rows for unmapped
coordinates, cross-chromosome chain hits, duplicated source coordinates before
mapping, or duplicated target coordinates after mapping, but it must not write
an empty artifact.

After successful liftover:

- `CHR` and `POS` are target-build coordinates.
- `SNP` is unchanged.
- logged coordinate provenance records the target build.
- `config_snapshot.genome_build` is the target build.
- logged coordinate provenance keeps coordinate basis as `"1-based"`.

## Liftover Methods

### Chain File

Move the existing `LiftOverTranslator` from `_kernel/ref_panel_builder.py` to
shared `_kernel/liftover.py`; reference-panel code imports the shared class with
no behavior change. Chain-file liftover uses `pyliftover`; if the optional
dependency is missing, raise the existing LDSC dependency error with an install
hint for the liftover extra.

Chain behavior:

- No-hit mappings are dropped and counted as unmapped.
- Cross-chromosome hits are dropped, counted separately, and included in total
  dropped rows.
- If a chain hit list contains both cross-chromosome and same-chromosome hits,
  the first same-chromosome hit is kept.
- Logs show count summaries by default; up to 5 examples per drop category are
  available only at `DEBUG`.

### HM3 Quick Liftover

Use the provided curated map packaged as:

```text
src/ldsc/data/hm3_curated_map.tsv.gz
```

Source file:

```text
/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/resources/SNP_lists/hm3_curated_map.tsv.gz
```

Runtime loader requirements:

- Required logical columns: chromosome, hg19 position, hg38 position, and SNP
  label. Resolve them with the package's existing column-alias inference
  helpers instead of exact-only names.
- Canonical preferred names are `CHR`, `hg19_POS`, `hg38_POS`, and `SNP`.
- Extra columns are ignored.
- Duplicate `(CHR, hg19_POS)` or duplicate `(CHR, hg38_POS)` in the map are a
  hard package-data validation error.

Lookup is coordinate-only:

- hg19 to hg38: key `(CHR, hg19_POS)`, write `hg38_POS`.
- hg38 to hg19: key `(CHR, hg38_POS)`, write `hg19_POS`.
- `SNP` from the map is used only for diagnostics/provenance, never for lookup.

Rows missing from the HM3 map are dropped and reported. Detailed HM3 provenance
is recorded in the workflow log; the CLI flag remains
`--use-hm3-quick-liftover`.

## Duplicate Coordinates

For sumstats liftover, duplicated source `CHR/POS` groups are dropped before
mapping. After either liftover method, if multiple original rows map to the same
target `CHR/POS`, every row in that duplicated target-coordinate group is
dropped. The munger does not keep the first row. It records separate
`n_duplicate_source_dropped` and `n_duplicate_target_dropped` counts and exposes
examples only at `DEBUG`.

If source or target duplicate removal leaves no rows, liftover errors instead
of writing an empty artifact.

For reference-panel PLINK builds, the same duplicate-coordinate mechanics apply
only when the active `GlobalConfig.snp_identifier` is `chr_pos`. There is no
public duplicate-position policy knob: source duplicate groups are dropped
before chain mapping and target duplicate groups after mapping with `drop-all`.
Per-chromosome `dropped_snps/chr{chrom}_dropped.tsv.gz` sidecars are always
written for processed chromosomes and cover the ref-panel-applicable reasons:
`source_duplicate`, `unmapped_liftover`, `cross_chromosome_liftover`, and
`target_collision`. If all SNPs drop for one chromosome, that chromosome is
skipped; the full run errors only if no chromosome emits artifacts. In
source-only `rsid` builds, coordinate duplicate filtering is logged once as not
applicable.

## Metadata And Logging

The sidecar is intentionally thin. It stores only fields needed for downstream
workflow behavior: the sidecar marker, optional trait label, and the full
`GlobalConfig` snapshot used for compatibility checks and reloads.

```jsonc
{
  "format": "ldsc.sumstats.v1",
  "trait_name": "Trait label",
  "config_snapshot": {
    "snp_identifier": "chr_pos",
    "genome_build": "hg38",
    "log_level": "INFO",
    "fail_on_missing_metadata": false
  }
}
```

Top-level entries:

- `format`: sidecar schema identifier.
- `trait_name`: optional trait label used by downstream summaries.

`config_snapshot` is the compatibility block used by downstream modules:

- `snp_identifier`: `rsid` or `chr_pos`; controls identity-key semantics.
- `genome_build`: concrete build for `chr_pos`, `null` for `rsid`; used by
  `validate_config_compatibility()`.
- `log_level`: log-level setting captured from the run.
- `fail_on_missing_metadata`: missing-metadata compatibility policy captured
  from the run.

The log file records detailed provenance that does not belong in the sidecar:

- coordinate provenance: coordinate source columns, resolved/final build,
  inference status, coordinate basis, missing-coordinate counts, and optional
  build-inference details;
- liftover report: `applied`, source/target builds, method, method resource path
  when present, input/retained/drop counts, missing-coordinate drop count,
  unmapped count, cross-chromosome count, duplicate-source drop count, and
  duplicate-target drop count;
- HM3 provenance: method token and packaged HM3 map path when HM3 quick liftover
  is used;
- output bookkeeping: selected output format, written output files, Parquet
  compression, and Parquet row-group layout.

`n_lifted` means the final retained row count after all liftover-specific drops,
so it should match the post-liftover output row count. Log messages are readable
key/value text with count summaries at normal verbosity and examples only at
`DEBUG`; they are not JSON payloads. In no-op reports, `source_build` and
`target_build` are both `null`; `null` means the field is not applicable because
liftover did not run.

## Package-Wide chr_pos Identity Rule

This feature depends on a package-wide invariant:

- In `rsid` mode, identity and merges use `SNP`.
- In `chr_pos` mode, identity and merges use normalized `CHR/POS`.
- In `chr_pos` mode, rows missing `CHR` or `POS` are dropped and logged at
  merge, match, map, and restriction boundaries. They may still be preserved by
  non-matching stages such as basic sumstats QC and artifact writing.
- Invalid non-missing coordinates, such as non-numeric POS, zero/negative POS,
  non-integer POS, or unsupported chromosome labels, remain errors.
- `SNP` is always a label and should be assumed to contain rsIDs, even when it
  does not.

Implementation must audit identity-sensitive operations and make them
mode-aware. Required changes include:

- `SumstatsTable.snp_identifiers()` uses `CHR/POS` in `chr_pos` mode.
- `SumstatsTable.subset_to()` uses `CHR/POS` in `chr_pos` mode.
- `SumstatsTable.align_to_metadata()` aligns by `CHR/POS` in `chr_pos` mode.
- `_row_alignment.assert_same_snp_rows()` accepts identifier mode; in `chr_pos`
  mode it compares `CHR/POS` and ignores `SNP` label mismatches.
- Regression `h2` and `rg` already use private `CHR/POS` keys in `chr_pos`
  mode; tests must lock this in.
- Annotation, reference-panel, LD-score, and regression merge/filter sites must
  be audited for unconditional `SNP` identity use and replaced with mode-aware
  keys where applicable.

Operations that are explicitly label/allele based, such as `--merge-alleles`,
may remain `SNP` based if their contract is rsID/label matching rather than
coordinate identity.

## Files To Modify

| File | Change |
|---|---|
| `src/ldsc/_kernel/liftover.py` | Shared internal liftover module with chain translator, HM3 curated-map loader/lifter, duplicate-coordinate helpers, and readable drop reports. |
| `src/ldsc/_kernel/ref_panel_builder.py` | Import relocated translator. |
| `src/ldsc/_kernel/sumstats_munger.py` | Insert liftover after SNP filtering and before output writing. |
| `src/ldsc/config.py` | Add munger target/method fields only; reference-panel duplicate handling has no public policy field. |
| `src/ldsc/sumstats_munger.py` | Workflow validation, sidecar metadata, dropped-SNP audit sidecar, and mode-aware `SumstatsTable` helpers. |
| `src/ldsc/_row_alignment.py` | Make row-alignment checks identifier-mode aware. |
| `src/ldsc/data/hm3_curated_map.tsv.gz` | Packaged provided HM3 curated map. |
| `src/ldsc/cli.py` | Register new CLI flags. |
| `tests/` | Add focused config, liftover, metadata, and package-wide identity tests. |

## Verification

Required test coverage:

1. Config and CLI flag validation.
2. HM3 curated map load, column-alias inference, and duplicate-key rejection.
3. HM3 coordinate-only lookup with mismatching `SNP` labels.
4. Chain and HM3 row dropping.
5. Cross-chromosome drop counts.
6. Duplicate target-coordinate group removal.
7. All-dropped error.
8. Missing coordinate error.
9. `SNP` unchanged after liftover.
10. No-op and applied metadata schema.
11. Thin sidecar metadata and detailed log provenance.
12. Package-wide `chr_pos` identity helpers and merge behavior.
13. Lifted hg38 `chr_pos` sumstats reject hg19 LD scores through existing
    `ConfigMismatchError`.
14. Lifted hg38 `chr_pos` sumstats proceed with hg38 LD scores.

## 2026-05-10/11 Harmonization Update

The follow-up liftover harmonization design supersedes the reference-panel
duplicate-policy details in this original sumstats-liftover spec. See
[`2026-05-10-liftover-harmonization.md`](2026-05-10-liftover-harmonization.md).

Updated contract:

- Reference-panel `duplicate_position_policy` / `--duplicate-position-policy`
  is removed. `drop-all` is the only coordinate duplicate behavior in
  `chr_pos` mode.
- Sumstats writes an always-present `dropped_snps/dropped.tsv.gz` audit sidecar
  with unified columns `CHR`, `SNP`, `source_pos`, `target_pos`, and `reason`.
  It may contain `missing_coordinate`, `source_duplicate`,
  `unmapped_liftover`, `cross_chromosome_liftover`, and `target_collision`.
- Reference-panel writes always-present per-chromosome
  `dropped_snps/chr{chrom}_dropped.tsv.gz` sidecars for processed chromosomes.
  It uses the same column names and nullable dtype contract, but never emits
  `missing_coordinate` because PLINK BIM `BP` is structurally non-null.
- Default logs are count-only with sidecar pointers. Row examples are available
  only at `DEBUG`.
