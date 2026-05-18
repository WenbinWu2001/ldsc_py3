# Allele-Aware SNP Identity Design

Date: 2026-05-15
Branch: `restructure`

## General Policies

The package supports exactly four public SNP identifier modes, with no mode
aliases:

- `rsid`
- `rsid_allele_aware`
- `chr_pos`
- `chr_pos_allele_aware`

`chr_pos_allele_aware` is the package default. Column-name aliases such as
`RSID`, `BP`, or `POSITION` remain input-column aliases; they are not aliases for
the `snp_identifier` mode value.

Effective merge keys are:

| Mode | Effective merge key |
| --- | --- |
| `rsid` | `rsID` |
| `rsid_allele_aware` | `rsID:<allele_set>` |
| `chr_pos` | `CHR:POS` |
| `chr_pos_allele_aware` | `CHR:POS:<allele_set>` |

Duplicate policy: compute the effective merge key for the active mode, then drop all rows in duplicate-key clusters. For any artifact-like table, compute the merge key for the active SNP identifier mode and drop all rows in duplicate-key clusters.

`<allele_set>` is an unordered, strand-aware biallelic SNP allele token. For
identity, `A:C`, `C:A`, `T:G`, and `G:T` normalize to the same allele set.
Ordered `A1/A2` still matter for signed-statistic orientation and harmonization,
but not for the identity token itself.

Allele-aware modes drop rows with missing alleles, non-ACGT alleles, non-SNP
alleles, identical allele pairs, and strand-ambiguous pairs (`A/T`, `T/A`,
`C/G`, `G/C`). Allele-aware modes also drop multi-allelic base-key clusters:
rows with the same base key (`rsID` in rsID-family modes, `CHR:POS` in
coordinate-family modes) and more than one allele set. This conservative v1
policy means allele-aware modes are used for safer merging, not for retaining
multi-allelic sites in LDSC analyses.

Artifact row cleanup and restriction cleanup differ intentionally. Artifact rows
with duplicate effective keys use drop-all. Restriction files are set-like:
duplicate effective restriction keys collapse to one key, while invalid allele
rows and multi-allelic restriction clusters are dropped and reported.

## Scope

This design covers package-wide SNP identity and merge-key harmonization for:

- summary-statistics munging;
- PLINK reference-panel building and R2 parquet artifacts;
- reference-panel loading;
- annotation and LD-score construction;
- h2, partitioned-h2, and rg regression preprocessing;
- SNP restriction files used by the above workflows.

The design intentionally does not provide backward compatibility for old
package-written artifacts. Existing artifacts from previous schema generations
must be regenerated.

## Semantics Of `A1/A2`

`A1/A2` remain the package-written allele column names everywhere. They are not
genome `REF/ALT` columns.

Context-specific meanings:

- In raw or munged summary statistics before harmonization, `A1` is the allele
  that the signed statistic is relative to, and `A2` is the counterpart allele.
- In reference-panel and LD-score artifacts, `A1/A2` are the reference-panel
  orientation columns.
- In h2 and partitioned-h2 preprocessing under allele-aware modes, sumstats are
  oriented to LD-score/reference-panel `A1/A2`.
- In rg preprocessing, trait 2 is oriented to trait 1 after the two sumstats are
  merged by identity.

The unordered allele set is used for identity. Ordered `A1/A2` are used only
when signed-statistic orientation matters.

## Shared Internal Service

Create a shared internal module:

```text
src/ldsc/_kernel/snp_identity.py
```

The module owns:

- exact mode validation for the four public modes;
- mode-family helpers (`rsid` family versus `chr_pos` family);
- base-key construction;
- unordered strand-aware allele-set normalization;
- effective merge-key construction;
- allele QC for identity;
- multi-allelic base-key cluster detection in allele-aware modes;
- duplicate effective-key detection for artifact tables;
- restriction matching and restriction cleanup semantics;
- drop reports with controlled reason vocabulary.

Low-level restriction file parsing, delimiter handling, and raw row extraction
remain in `src/ldsc/_kernel/identifiers.py`. Mode-aware matching delegates to
`snp_identity.py`.

Controlled identity drop reasons:

- `missing_allele`
- `invalid_allele`
- `strand_ambiguous_allele`
- `multi_allelic_base_key`
- `duplicate_identity`

Existing liftover reasons remain separate:

- `missing_coordinate`
- `source_duplicate`
- `unmapped_liftover`
- `cross_chromosome_liftover`
- `target_collision`

## Minimal Artifact Provenance

Reloadable package-written identity artifacts record minimal provenance:

- `schema_version`
- `artifact_type`
- `snp_identifier`
- `genome_build`

Workflow-specific technical metadata remains where already required, for
example R2 sorted build, row-group metadata, R2 bias/sample-size metadata, count
records, and LD-score output layout metadata.

Runtime events that are not stable artifact identity belong in logs, not in the
minimal metadata. Examples include identity downgrade, original input modes in a
downgraded regression run, effective regression mode, duplicate rows dropped
during downgrade filtering, and identity cleanup counts.

Reloading rules:

- Hard validation of recorded identity mode against required allele fields applies
  to reloadable identity artifacts: munged sumstats, reference-panel artifacts,
  and LD-score artifacts. It does not apply to external annotations, external
  restriction files, or regression summary outputs.
- Old package-written artifacts without current schema/provenance are rejected
  with a direct "regenerate with current LDSC package" message.
- Current artifacts whose recorded mode requires fields that are missing are
  rejected as malformed.
- Raw sumstats, external annotations, and user restriction files are external
  inputs, not package-written artifacts; they use their own parsing and repair
  rules.

## Artifact Output Contract

### Munged Sumstats

In `rsid_allele_aware` and `chr_pos_allele_aware`, munged sumstats require
`A1/A2`. In `rsid` and `chr_pos`, `A1/A2` remain optional so allele-free h2 and
partitioned-h2 workflows continue to work.

Package-written sumstats artifacts record the active `snp_identifier` and
`genome_build` in current metadata. Allele-aware sumstats artifacts without
`A1/A2` are malformed.

Dropped-SNP audit sidecars remain always written:

```text
dropped_snps/dropped.tsv.gz
```

Clean runs write a header-only sidecar.

### Reference-Panel Builder Outputs

Package-built R2 parquet files write endpoint allele columns whenever available:

```text
A1_1, A2_1, A1_2, A2_2
```

These columns are required in allele-aware modes and ignored for identity in
base modes.

Reference-panel metadata sidecars write current identity provenance plus
`A1/A2` whenever available. Alleles are required in allele-aware modes, and old
package-written metadata sidecars without current provenance must be
regenerated.

Reference-panel dropped-SNP sidecars remain always written per processed
chromosome:

```text
dropped_snps/chr{chrom}_dropped.tsv.gz
```

Clean chromosomes write header-only sidecars.

### LD-Score Outputs

LD-score baseline and query artifacts include `A1/A2` whenever available. They
are required in allele-aware modes and ignored for identity in base modes.

LD-score artifacts written in allele-aware modes must be reloadable with their
recorded mode and allele columns. Missing allele columns are malformed.

### Dropped-SNP Sidecar Columns

The dropped-SNP sidecars extend the current audit shape with identity fields:

```text
CHR
SNP
source_pos
target_pos
reason
base_key
identity_key
allele_set
stage
```

`source_pos` and `target_pos` remain nullable because not every drop reason has
both positions. `base_key`, `identity_key`, and `allele_set` are nullable for
coordinate-only liftover reasons or base-mode drops where the field is not
applicable.

## Sequential Workflow Steps

### Shared Identity Cleanup For Artifact Tables

For a candidate artifact-like table:

1. Build the base key for the active mode.
2. In allele-aware modes, validate and normalize `A1/A2` into an unordered,
   strand-aware allele set.
3. In allele-aware modes, drop rows with missing, invalid, non-SNP, identical,
   or strand-ambiguous allele pairs.
4. In allele-aware modes, drop all rows in base-key clusters that contain more
   than one allele set.
5. Build the effective merge key for the active mode.
6. Drop all rows in duplicate effective-key clusters.
7. Return the cleaned table plus a drop report.

### Restriction Files

Restriction files define candidate membership before artifact identity cleanup.

1. Parse the restriction file in `identifiers.py`.
2. Detect whether the file has usable allele columns.
3. If no allele columns are present, match by base key. A base-key restriction
   can retain multiple allele-aware candidates; later identity cleanup decides
   whether those candidates survive.
4. If allele columns are present, rows with missing, invalid, non-SNP,
   identical, or strand-ambiguous alleles are dropped from the restriction set
   and reported with a warning.
5. In allele-aware modes, allele-bearing restriction files drop multi-allelic
   base-key clusters before matching.
6. Duplicate effective restriction keys collapse to one set key rather than
   using drop-all.
7. Match retained candidates by either base keys or effective allele-aware keys,
   depending on restriction precision.

### `munge-sumstats`

1. Parse and QC raw chunks.
2. Apply `--sumstats-snps-file` or HM3 restriction chunk-wise as an early
   candidate-universe filter for memory efficiency.
3. Keep existing early allele QC so invalid raw allele rows do not flow through
   expensive downstream steps.
4. Concatenate retained candidate chunks.
5. Apply sample-size processing and signed-statistic conversion.
6. If requested, apply coordinate-family liftover. Liftover is valid only for
   `chr_pos` and `chr_pos_allele_aware`.
7. Around liftover, apply the special coordinate-collision policy:
   source duplicate `CHR:POS` clusters are dropped before mapping, and target
   duplicate `CHR:POS` clusters are dropped after mapping. Allele set is ignored
   for liftover collision detection.
8. Once output identity is final, run shared identity cleanup globally.
9. Write curated sumstats, metadata, logs, and the always-written dropped-SNP
   sidecar.

### `ref-panel-builder`

1. Resolve PLINK source build and optional SNP restrictions.
2. Apply `--ref-panel-snps-file` or HM3 restriction to define the candidate
   universe.
3. Run shared identity cleanup on the retained candidate universe.
4. Reject matching chain-file liftover and HM3 quick liftover for `rsid` and
   `rsid_allele_aware`; reference-panel liftover is valid only for `chr_pos` and
   `chr_pos_allele_aware`.
5. If liftover is used, apply source and target coordinate-collision drops on
   `CHR:POS`, ignoring allele set.
6. If liftover emits both hg19 and hg38 outputs, synchronize drops across builds:
   a variant dropped because of source or target coordinate collision in either
   emitted build is dropped from all emitted builds.
7. Emit R2 parquet and metadata sidecars with allele columns when available.
8. Write the always-written per-chromosome dropped-SNP sidecar.

### `ldscore`

1. Load annotation bundle and reference-panel metadata/R2 inputs.
2. Apply reference-panel and regression SNP restrictions as candidate filters.
3. Run shared identity cleanup on the reference-panel universe before LD-score
   calculation.
4. Match annotations to the cleaned reference universe.
5. Annotation inputs do not require alleles. If annotations include alleles, they
   can participate in allele-aware matching. If they do not include alleles, they
   match by base key after the reference universe has already been cleaned.
6. Compute LD scores and regression-weight LD scores on the cleaned universe.
7. Write LD-score artifacts with `A1/A2` when available and required in
   allele-aware modes.

### `h2`, `partitioned-h2`, And `rg`

1. Load sumstats and LD-score artifacts.
2. By default, require exact `snp_identifier` mode compatibility.
3. If `--allow-identity-downgrade` is set, same-family allele-aware/base pairs
   may run under the base mode:
   - `chr_pos_allele_aware` plus `chr_pos` uses effective `chr_pos`;
   - `rsid_allele_aware` plus `rsid` uses effective `rsid`.
4. Cross-family mixes are always rejected: rsID-family modes never mix with
   coordinate-family modes.
5. In downgrade mode, run duplicate effective-key drop-all on each input under
   the effective base mode before merging.
6. Record minimal regression metadata:
   - `effective_snp_identifier`
   - `identity_downgrade_applied`
7. Log a readable downgrade message naming the original input modes, effective
   mode, and rows dropped during downgrade filtering.
8. In h2 and partitioned-h2 under allele-aware modes, orient sumstats `Z` to the
   LD-score/reference-panel `A1/A2` before regression.
9. In rg, merge by effective identity, then keep the current ordered-allele
   harmonization step to filter incompatible rows and flip trait 2's `Z` when
   needed.

## Error Handling And User Repair

When the default `chr_pos_allele_aware` mode requires alleles but none are
available, errors should include the repair path:

```text
This run is using snp_identifier='chr_pos_allele_aware', which requires A1/A2
allele columns. No usable allele columns were found. To run without allele-aware
SNP identity, rerun with --snp-identifier chr_pos.
```

When old package-written artifacts are loaded:

```text
This artifact was not written with the current LDSC schema/provenance contract.
Regenerate it with the current LDSC package.
```

When downgrade is used, logs must be explicit:

```text
Identity downgrade enabled: LD-score mode chr_pos, sumstats mode
chr_pos_allele_aware; running regression with effective
snp_identifier='chr_pos'. Dropped 42 duplicate effective-key rows before merge.
```

## Testing Requirements

The implementation should include focused tests for:

- exact mode validation: only the four public modes are accepted;
- removal of old `snp_identifier` mode aliases;
- allele-set normalization including unordered and strand-complement equivalent
  pairs;
- missing, invalid, identical, non-SNP, and strand-ambiguous allele drops;
- multi-allelic base-key cluster drops in allele-aware modes;
- duplicate effective-key drop-all for artifact rows;
- restriction files without alleles matching by base key;
- allele-bearing restriction rows with invalid alleles being dropped and warned;
- allele-bearing restriction multi-allelic clusters being dropped;
- duplicate restriction keys collapsing to one set key;
- munger final identity cleanup after liftover;
- ref-panel liftover rejection for `rsid` and `rsid_allele_aware`;
- ref-panel cross-build synchronized coordinate-collision drops;
- LD-score output requiring `A1/A2` in allele-aware modes;
- regression exact-mode compatibility by default;
- regression `--allow-identity-downgrade` within a family;
- regression cross-family rejection;
- h2/partitioned-h2 sumstats orientation to LD-score/reference alleles;
- rg ordered-allele harmonization after identity merge.

## Non-Goals For V1

V1 does not compute LDSC over retained multi-allelic same-base-key clusters.
Those clusters are dropped package-wide in allele-aware modes. A future design
can revisit full multi-allelic LD-score support, including same-position pair
semantics, window traversal assumptions, and R2 endpoint indexing.

V1 does not provide compatibility loaders for old package-written artifacts.
Existing outputs from earlier schema generations should be regenerated.
