# LDSC2 Artifact Compatibility Design

Last updated on: 2026-08-04

Status: approved for implementation

## Problem and goal

Large consortiums have already munged, quality-controlled, distributed, and
curated substantial collections of LDSC2 summary statistics. Requiring those
users to recover raw GWAS inputs and re-run LDSC3 munging would make otherwise
usable datasets inaccessible. Selected precomputed LDSC2 LD-score suites are
also expensive or operationally awkward to reproduce, especially published
European weight scores and baseline models.

LDSC3 nevertheless has a new canonical identity and result-directory contract.
Allowing every LDSC2 fragment directly into every workflow would make that
contract optional, spread legacy parsing through regression and numerical
kernels, and make provenance-dependent behavior difficult to audit.

This design establishes two deliberately different compatibility boundaries:

1. LDSC2 munged `.sumstats` and `.sumstats.gz` files are accepted automatically
   at the regression boundary and projected onto the canonical LDSC3 LD-score
   panel.
2. Selected reusable LDSC2 LD-score suites are accepted only by an explicit,
   one-time converter that writes a canonical LDSC3 LD-score directory.

All other regression inputs continue to follow current LDSC3 contracts. The
successful outcome is that approved legacy artifacts can be used without
recomputation while regression, output writers, and `_kernel` operate on one
canonical in-memory and on-disk representation.

## Scope

- automatic legacy munged-sumstats compatibility for `h2`, `partitioned-h2`,
  and `rg`;
- an explicit `convert-ldsc2-ldscores` command and public Python workflow;
- unpartitioned reference/weight suites such as `eur_w_ld_chr`;
- one self-contained, baseline-only partitioned reference suite plus a separate
  one-column regression-weight suite and frequency suite;
- allele-aware and allele-unaware panel projection for legacy sumstats;
- current canonical LD-score Parquet, metadata, count, and overlap artifacts;
- deterministic discovery, validation, provenance, logging, and issue audits;
- removal of obsolete private LDSC2 LD-score emitters and unused regression
  readers from `_kernel`;
- rejection of footerless curated Parquet sumstats;
- `BP` serialization for `.annot.gz` while retaining canonical in-memory
  `POS`; and
- reconciliation of current architecture, I/O inventory, workflow, user, and
  troubleshooting documentation.

## Out of scope

- direct use of LDSC2 `.l2.ldscore`, `.M`, `.M_5_50`, `.annot`, or frequency
  fragments by regression;
- automatic LD-score conversion during `h2`, `partitioned-h2`, or `rg`;
- conversion of query, cell-type-specific, or composed multi-prefix reference
  models;
- thin annotations;
- accepting a legacy annotation-overlap output as authoritative;
- attaching or imputing alleles onto legacy LD-score suites;
- imputing frequency values between sumstats and LD-score/reference artifacts;
- liftover of legacy suites or legacy sumstats during compatibility matching;
- configurable common-MAF thresholds in the converter;
- footerless Parquet compatibility;
- a general LDSC3-to-LDSC2 interoperability promise; and
- changing PLINK `.bed/.bim/.fam` support, ordinary annotation inputs, or other
  domain-standard source formats merely because LDSC2 also used them.

## Compatibility matrix

| Artifact | Entry point | Compatibility | Canonical boundary |
| --- | --- | --- | --- |
| LDSC2 `.sumstats` / `.sumstats.gz` | `h2`, `partitioned-h2`, `rg`, public `RegressionRunner` | Automatic | Projected to the supplied canonical LD-score panel before ordinary dataset assembly |
| Current `sumstats.parquet` with complete footer metadata | Regression | Native | Footer identity remains authoritative |
| Footerless `sumstats.parquet` | Regression | Unsupported | Reject; never reinterpret as legacy |
| Unpartitioned LDSC2 reference and weight suites | `convert-ldsc2-ldscores` | Explicit conversion | New canonical LD-score directory with `base` and `regression_ld_scores` |
| Self-contained baseline partitioned suite | `convert-ldsc2-ldscores` | Explicit conversion | New canonical baseline table, counts, and overlap artifact |
| Query/cell-type-specific LD-score suite | None | Unsupported | Recompute through current LDSC3 query contracts |
| Thin annotation | None | Unsupported by converter | Recreate a full annotation or recompute |
| `.M`, `.M_5_50`, `.annot`, `.frq` | Converter components only | Not standalone regression inputs | Validated and represented in canonical metadata/overlap |
| LDSC3 `.annot.gz` output | `ldscore`/`annotate` input | Native LDSC3 text artifact | Writes `BP`, reads either `BP` or `POS` |
| LDSC3 gzip sumstats output | Regression | Metadata-free text compatibility path | Reloaded through the legacy projection boundary; no reverse guarantee |

No public regression or kernel function accepts a legacy LD-score path or
prefix after this work.

## Public command contract

The unified CLI adds:

```text
ldsc convert-ldsc2-ldscores \
    --legacy-reference-dir REF_DIR \
    --legacy-weight-dir WEIGHT_DIR \
    --output-dir OUTPUT_DIR
```

Baseline partitioned conversion additionally requires:

```text
--legacy-frequency-dir FRQ_DIR
```

Optional controls are:

```text
--snp-identifier rsid|chr_pos       # default: rsid
--genome-build auto|hg19|hg38       # default: auto
--overwrite
--log-level DEBUG|INFO|WARNING|ERROR
```

The converter has no `--common-maf-min`, profile, prefix, glob, column-map, or
allele option. Users reusing one unpartitioned suite for both reference and
weight roles pass the same path to both directory arguments.

The converter always processes autosomes 1 through 22 as one complete suite.

## Legacy sumstats classification

`load_sumstats()` classifies only files ending in `.sumstats` or
`.sumstats.gz` as legacy munged sumstats. It reads them as whitespace-delimited
text and marks the returned `SumstatsTable` with explicit legacy source
provenance. A missing `config_snapshot` alone is not the classifier because
current Parquet without identity metadata is invalid rather than legacy.

Legacy tables require `SNP`, `A1`, `A2`, `Z`, and `N`. `FRQ` is optional.
Allele-less legacy tables are rejected for `h2`, `partitioned-h2`, and `rg`.
This stricter uniform rule keeps compatibility behavior independent of whether
a particular estimator happens to square `Z`.

Current self-describing Parquet sumstats retain their footer provenance and are
never reinterpreted by the compatibility layer. An `rg` run may mix current
Parquet and legacy text traits when every current trait is normally compatible
with the supplied panel.

## Legacy sumstats projection algorithm

Projection occurs in `RegressionRunner` before ordinary identity resolution,
configuration compatibility checks, or regression dataset merging. It returns
a new working `SumstatsTable`; it never mutates the source file or caller-owned
DataFrame.

### Source cleanup

For each legacy trait:

1. normalize `SNP` as a nonempty string and uppercase `A1`/`A2`;
2. reject missing required file columns;
3. apply drop-all duplicate handling to source rsIDs, so every row in a
   duplicated source-rsID cluster is excluded;
4. validate alleles as single-base `A`, `C`, `G`, or `T` values;
5. treat missing, invalid, and strand-ambiguous `A/T` or `C/G` source pairs as
   row-level drops; and
6. parse optional `FRQ` numerically without making it a retention field.

Missing, nonnumeric, infinite, or out-of-range `FRQ` becomes missing and emits
an aggregate warning. It does not drop a SNP because regression does not
consume frequency.

### Allele-aware panel

When the canonical LD-score panel has an allele-aware identifier mode and
usable `A1`/`A2`:

1. find every panel candidate with the source rsID;
2. retain candidates compatible through direct, complement, swap, or
   swapped-complement orientation;
3. accept the source row only when exactly one compatible panel candidate
   remains;
4. keep `Z` for direct or complement orientation;
5. negate `Z` for swap or swapped-complement orientation;
6. transform a valid source `FRQ` to `1 - FRQ` on a swap;
7. replace `CHR`, `POS`, `SNP`, `A1`, and `A2` with panel values; and
8. run subsequent regression under the panel's recorded identifier mode.

A panel rsID need not initially be unique because alleles may select one
candidate. Zero compatible candidates are `incompatible_alleles`; multiple
compatible candidates are `ambiguous_panel_mapping`.

### Allele-unaware panel

When the canonical panel uses `rsid` or `chr_pos` and has no authoritative
alleles:

1. find candidates by source rsID;
2. accept only a unique panel row for that rsID;
3. replace source `CHR`, `POS`, and `SNP` with panel values;
4. retain the source `A1`, `A2`, `Z`, and valid `FRQ` orientation; and
5. run subsequent regression under the panel's recorded identifier mode.

Multiple panel candidates cannot be disambiguated without alleles and are
dropped as `ambiguous_panel_mapping`. `h2` and `partitioned-h2` require no
further sign operation because their response uses `Z` squared.

For `rg`, project both traits independently, then harmonize trait 2 to trait 1
using the same direct/complement/swap/swapped-complement rules. Negate trait 2
`Z` and transform its valid `FRQ` to `1 - FRQ` when swapped. Incompatible or
strand-ambiguous trait pairs are dropped. Trait 1 is the pairwise orientation
anchor; neither trait is oriented to an allele-unaware panel.

Missing `FRQ` must not enter a blanket `dropna()` condition anywhere in h2 or
rg assembly. Missingness filtering is limited to fields actually consumed by
the estimator and allele harmonizer.

### Legacy sumstats drop audit

When a regression output directory is supplied and at least one legacy trait
is used, the workflow owns:

```text
diagnostics/dropped_snps/legacy_sumstats.tsv.gz
```

The file is written header-only when clean and contains:

| Column | Meaning |
| --- | --- |
| `trait_name` | Resolved trait label |
| `source_path` | Legacy source path |
| `SNP` | Source rsID when available |
| `A1`, `A2` | Source alleles when available |
| `reason` | Stable drop reason |
| `panel_candidate_count` | Number of same-rsID panel candidates before final acceptance |

Stable reasons include `duplicate_source_rsid`, `missing_panel_rsid`,
`missing_allele`, `invalid_allele`, `strand_ambiguous`,
`incompatible_alleles`, and `ambiguous_panel_mapping`. File-level missing
columns and zero retained rows abort the affected regression. No new
retention-percentage threshold is introduced.

Logs report aggregate counts by reason at `INFO` or `WARNING`; example rsIDs
are `DEBUG`-only. Stdout-only regression reports aggregates without writing a
sidecar.

## Converter suite discovery

Discovery is directory-based and deterministic. LD-score families match either
of these layouts:

```text
<chrom>.l2.ldscore[.gz]
<prefix><chrom>.l2.ldscore[.gz]
```

The chromosome token is the numeric token immediately before `.l2.ldscore`.
Examples include `1.l2.ldscore.gz`, `weights.1.l2.ldscore.gz`, and
`baseline.1.l2.ldscore.gz`. Files such as `6_old.l2.ldscore.gz` do not belong
to a numeric chromosome family.

Files are grouped by the exact inferred prefix. Each source role must resolve
to exactly one coherent 1--22 family. Zero or multiple coherent families are
errors. Annotation and count files use the selected reference prefix;
frequency files resolve one analogous `<prefix><chrom>.frq[.gz]` family from
the frequency directory.

When both plain and gzip representations exist for one logical shard, compare
their decompressed bytes. Accept and select the gzip representation only when
they agree; otherwise fail. Ignored non-family files and stale alternatives are
recorded in diagnostics and provenance but never selected heuristically.

## Converter profile detection

There is no user profile flag.

- Exactly one reference LD-score value column and no annotation family selects
  the unpartitioned profile.
- A complete full annotation family with a bijective annotation/LD-score column
  mapping selects the baseline partitioned profile, including a deliberately
  single-column full annotation suite.
- Multiple reference LD-score columns without a complete annotation family are
  invalid.
- A partially present annotation family is invalid and cannot silently fall
  back to unpartitioned conversion.

The converter assumes an advanced caller supplies a baseline suite. It does not
try to infer biological intent from directory names or annotation values and
does not require an all-ones root annotation.

## Legacy converter schemas

### Reference and weight LD-score tables

Every shard requires `CHR`, `SNP`, and `BP`. Optional legacy `CM` and `MAF`
columns are parsed for diagnostics or source-build evidence but are not
scientific LD-score columns.

After excluding metadata:

- the weight suite must contain exactly one numeric LD-score column;
- an unpartitioned reference suite must contain exactly one numeric LD-score
  column; and
- a baseline reference suite may contain one or more numeric LD-score columns.

Scientific values must be finite and nonmissing. `CHR`, `SNP`, and `BP` must be
valid for every reference row. Duplicate rsIDs within a shard or across the
complete family are errors rather than row-level cleanup.

### Full annotations

Every baseline annotation shard requires `CHR`, `BP`, `SNP`, and `CM`, followed
by one or more numeric annotation columns. Thin annotations and non-finite
scientific values are errors. Annotation rsIDs must be globally unique.

For each LD-score column, column matching tries:

1. exact equality with an annotation name; then
2. if the LD-score name has one terminal `L2`, equality after removing that
   terminal suffix.

The mapping must be complete, one-to-one, and unambiguous in both directions.
The exact LD-score table name is the canonical output name; `L2` is never added
or removed in output. Count records and overlap axes use those exact output
names.

### Counts

Count files are headerless whitespace-delimited numeric vectors in LD-score
column order.

- `.l2.M_5_50` is required for every chromosome and both profiles.
- `.l2.M` is optional per chromosome.
- Vector length must equal the reference LD-score column count.
- Nonnumeric, missing, negative, or non-finite count values are errors.

For unpartitioned conversion, source counts are authoritative because the
reference SNP annotation matrix is absent. If any chromosome lacks `.l2.M`,
the aggregate all-SNP count is unavailable rather than partially summed.

For baseline conversion, counts are reconstructed from annotations and
frequencies as described below. Present legacy values are validation values;
missing `.l2.M` values are reconstructed.

### Frequencies

Every baseline frequency shard requires `SNP` and exactly one usable `MAF` or
`FRQ` column. Optional allele and sample-size columns are ignored. Frequency
rsIDs must be globally unique.

Annotations are left-joined to frequency by rsID. Every annotation SNP must
have exactly one finite frequency in its allowed range: `[0, 0.5]` for `MAF`
or `[0, 1]` for `FRQ`. Missing annotation frequencies are errors. Extra
frequency SNPs are ignored and audited.

## Converter identity and coordinates

Every relationship among reference scores, weights, annotations, frequencies,
and counts is resolved by rsID. Cross-family `CHR`/`BP` equality is never a join
or retention requirement. Coordinate disagreements warn and produce issue rows
but do not discard an otherwise unique rsID match.

The reference LD-score suite supplies canonical output `CHR`, `POS`, and `SNP`.
The converter never attaches `A1`/`A2`; only allele-unaware output modes are
accepted.

### `rsid` output

- canonical `snp_identifier` is `rsid`;
- canonical `genome_build` is null;
- `auto` performs best-effort reference, annotation, and weight build inference
  for provenance;
- unresolved reference inference warns and succeeds;
- an explicit reference build is validated when inference is decisive and is
  otherwise retained as a source declaration; and
- inferred annotation/weight disagreement is diagnostic only.

### `chr_pos` output

- all source-family joins still use rsID;
- canonical identity uses the reference LD-score `CHR/POS` after the join;
- `auto` must resolve the reference build or conversion fails;
- an explicit reference build that contradicts decisive reference evidence is
  an error;
- detected zero-based reference coordinates are shifted to canonical one-based
  `POS` in output only; and
- duplicate canonical reference `CHR/POS` clusters are conversion errors.

Source files are never modified and no liftover is performed.

## Reference/weight row assembly

For each chromosome:

1. validate reference and weight shards independently;
2. inner-join their unique rows by rsID;
3. preserve reference `CHR`, `BP`, `SNP`, and scientific LD-score values;
4. rename only the weight suite's single scientific value to
   `regression_ld_scores`;
5. record reference-only, weight-only, and retained rsIDs; and
6. sort output by canonical chromosome and position.

After chromosome assembly, validate global rsID uniqueness and, for `chr_pos`,
global coordinate uniqueness. No retention-percentage cutoff is applied. An
empty aggregate intersection is an error. Low final regression SNP counts are
reported later by regression after sumstats have also been intersected.

Unpartitioned output renames the reference scientific LD-score column to
`base`. Baseline output preserves all exact reference LD-score names.

## Common-count and overlap mathematics

The converter's common threshold is fixed at `0.05`; neither CLI nor Python API
accepts another value.

For a frequency column already expressed as MAF, the legacy common mask is:

```text
MAF > 0.05
```

For an allele-frequency column, the equivalent mask is:

```text
0.05 < FRQ < 0.95
```

Let `A_c` be one chromosome's full annotation matrix in canonical output-column
order and `C_c` its common-row restriction. Baseline conversion reconstructs:

```text
M_c       = column_sum(A_c)
M_5_50_c  = column_sum(C_c)
O_all_c   = A_c.T @ A_c
O_common_c = C_c.T @ C_c
```

Genome-wide values are chromosome sums. Marginal `M` values are column sums,
not overlap diagonals; those quantities differ for continuous annotations.

Present `.M` and `.M_5_50` values are checked per chromosome and annotation
before aggregation. Binary/integer annotation counts must agree exactly.
Continuous counts use `numpy.isclose` with `rtol=1e-8` and `atol=1e-6`.
Agreement preserves the exact legacy value. A missing baseline `.M` uses the
reconstructed value. A disagreement beyond tolerance fails conversion with no
override because marginal counts and overlaps would otherwise describe
different SNP universes.

The canonical count configuration records:

```json
{
  "common_reference_snp_maf_min": 0.05,
  "common_reference_snp_maf_operator": ">",
  "common_reference_snp_semantics": "legacy_ldsc2"
}
```

For baseline output, `ldscore.overlap.parquet` stores the full
baseline-by-baseline all/common block and no query diagonal. `query_columns` is
empty and there is no `ldscore.query.parquet`.

For unpartitioned output there is no overlap artifact. When any unpartitioned
chromosome lacks `.M`, its `base` count record contains
`all_reference_snp_count: null`; JSON `NaN` is forbidden. Common regression
continues to work. An all-count request raises a targeted error before dataset
assembly.

The LD-score loader validates that `count_config` and `overlap_config` agree on
the fixed threshold and operator. Inconsistent internal metadata is rejected,
so the threshold cannot drift through downstream use.

## Canonical converter output

Successful unpartitioned conversion writes:

```text
<output_dir>/metadata.json
<output_dir>/ldscore.baseline.parquet
<output_dir>/diagnostics/convert-ldsc2-ldscores.log
<output_dir>/diagnostics/conversion_issues.tsv.gz
```

The baseline table schema is:

```text
CHR POS SNP regression_ld_scores base
```

Successful baseline conversion additionally writes:

```text
<output_dir>/ldscore.overlap.parquet
```

Its baseline table schema is:

```text
CHR POS SNP regression_ld_scores <exact legacy LD-score columns...>
```

Annotations are validation and overlap inputs; raw annotation values are not
copied into the canonical result directory.

The root artifact remains ordinary `artifact_type="ldscore"`. Native regression
therefore needs no imported-artifact file reader. Root metadata adds a
`legacy_ldsc2_import` object containing:

- conversion profile and package/converter version;
- source directory tokens and selected prefixes;
- selected, ignored, duplicate-representation, and missing files;
- streaming SHA-256 hashes for selected source files;
- per-chromosome reference-only, weight-only, and retained row counts;
- source/effective build and coordinate-basis inference;
- count origin per chromosome/column: `legacy_validated`,
  `legacy_unvalidated`, `reconstructed`, or `missing`;
- fixed legacy common-frequency semantics;
- coordinate disagreement summaries; and
- relative diagnostic paths.

For an unpartitioned suite, `legacy_unvalidated` means an available source count
could not be independently reconstructed; it does not mean validation failed.

## Converter diagnostics and failure behavior

`diagnostics/conversion_issues.tsv.gz` has a stable schema:

```text
severity source_role file chromosome SNP annotation reason observed expected details
```

It contains warnings and errors, and is header-only after a clean conversion.
Stable row-level reasons cover duplicate rsIDs, missing frequencies, invalid
scientific values, coordinate disagreements, discarded non-family files, and
reference/weight membership differences. Count mismatches use one row per
chromosome/annotation because scalar source files cannot reveal exact causal
SNP membership. Frequency-boundary SNPs are reported when they could explain a
common-count discrepancy; LDSC3 never fabricates a precise causal-SNP claim.

Raised conversion errors name the violated invariant, role and file,
chromosome or annotation, affected count, and a bounded rsID sample. Complete
available rows remain in the issue table. Diagnostics may be written for the
current failing invariant without promising collection of unrelated errors
that occur later in the pipeline.

The converter preflights its full owned family before reading large inputs.
During validation it may create only `diagnostics/`. Root metadata and Parquet
files are written only after all validation succeeds. A failed diagnostics-only
directory is not loadable by regression. A corrected rerun requires
`--overwrite`, consistent with other workflow-owned failure diagnostics.
Original legacy files remain read-only.

## Public Python and module boundary

One workflow-layer module owns conversion. It exposes:

- `LegacyLDScoreConverter`;
- `LegacyLDScoreConversionResult`; and
- `convert_ldsc2_ldscores(...)`.

The convenience function mirrors the CLI directory, identity, build, overwrite,
and logging arguments. There is no separate public configuration hierarchy and
no common-threshold parameter. The implementation reuses `LDScoreResult`,
`LDScoreOverlap`, the pure overlap mathematics, and
`LDScoreDirectoryWriter`.

`ldsc.__init__` re-exports the converter, result, and convenience function.
`ldsc.cli` owns command registration and dispatch.

Legacy filename discovery, text parsing, validation, and diagnostics remain in
the workflow module. They are not added to `_kernel`. Regression continues to
load only the newly written canonical directory through
`load_ldscore_from_dir()`.

## Canonical metadata and reader changes

`LDScoreResult` gains one optional additive import-provenance field consumed by
`LDScoreDirectoryWriter.build_metadata()`. Native LD-score results leave it
unset. Readers ignore the additive field for computation.

The writer must use the actual `count_config` operator in `overlap_config`
rather than hard-coding LDSC3's native `>=`. The reader validates internal
agreement among:

- `count_config.common_reference_snp_maf_min`;
- `count_config.common_reference_snp_maf_operator`;
- `overlap_config.common_maf_min`; and
- `overlap_config.common_maf_operator`.

Count assembly treats null all-SNP counts as unavailable. It builds an
all-count vector only when every requested model column has a finite all count.
An explicit all-count request against an imported unpartitioned suite with null
counts raises a message naming the directory and missing columns. Common-count
selection never falls back silently to all counts when a malformed common
record is present.

## Kernel and obsolete compatibility cleanup

The current public `ldscore` workflow is the only LD-score computation and
output entry point. Remove from `_kernel.ldscore` the obsolete LDSC2 output
assembly and standalone output path, including writers for:

- `.l2.ldscore.gz`;
- `.w.l2.ldscore.gz`;
- `.l2.M`;
- `.l2.M_5_50`; and
- `.annotation_groups.tsv`.

Retain pure computational functions still called by `ldscore_calculator`.
Delete or privatize the obsolete standalone parser/run functions once their
call graph is empty.

From `_kernel.formats`, remove legacy sumstats, LD-score, count, and annotation
readers proven unused after obsolete kernel-path removal. Retain PLINK and other
format primitives still used by reference-panel workflows. The converter does
not reuse `_kernel.formats`. Update tests and module documentation so private
modules no longer claim regression directly consumes LDSC2 artifacts.

The metadata-free gzip sumstats writer remains in `sumstats_munger`, the
workflow layer. Correct stale `_kernel.sumstats_munger` documentation that
claims the kernel owns it.

## Footerless Parquet behavior

For a `.parquet` curated sumstats path, `load_sumstats()` requires the complete
current identity footer. Absence of `ldsc:artifact_type`,
`ldsc:snp_identifier`, or the required genome-build field is an artifact error
with regeneration guidance. The loader never assigns legacy provenance to
Parquet. Footerless-Parquet support and inference are intentionally absent.

## Annotation `BP` serialization

Annotation workflows retain canonical `POS` in memory. Immediately before
writing `query.<chrom>.annot.gz`, the serialization copy renames `POS` to `BP`
and writes the legacy positional metadata order:

```text
CHR BP SNP CM <annotation columns...>
```

The operation does not mutate `AnnotationBundle`. LDSC3 annotation readers
continue accepting both `BP` and `POS` and normalize to in-memory `POS`. This
small interoperability improvement is not a general promise that LDSC3 outputs
can be supplied to LDSC2.

## Documentation changes

After implementation and tests:

- add `convert-ldsc2-ldscores` to `io-argument-inventory.md`, command help,
  data flow, code structure, and user utility documentation;
- add the actual `build-gene-ldscore-index` command section to
  `io-argument-inventory.md`;
- replace the private legacy-emitter statement with the compatibility matrix;
- update regression docs to describe automatic legacy projection, allele
  behavior, mixed rg inputs, FRQ missingness, and audits;
- update LD-score docs with converter profiles, fixed common threshold,
  required frequency inputs, and missing unpartitioned `.M` behavior;
- update metadata inventories for import provenance, nullable all counts, and
  legacy common-count semantics;
- update troubleshooting for every stable converter and sumstats drop/error
  reason; and
- avoid describing gzip sumstats or `BP` annotations as a reverse LDSC2
  compatibility guarantee.

## Validation strategy and acceptance criteria

### Converter unit and workflow tests

1. Discover plain, gzip, empty-prefix, and named-prefix 22-chromosome families
   deterministically.
2. Prefer identical gzip/plain duplicates and reject conflicting duplicates.
3. Ignore and audit `6_old`-style files; reject missing chromosomes and
   ambiguous coherent families.
4. Auto-detect unpartitioned and baseline profiles and reject partial or
   unsupported profiles.
5. Validate exact/terminal-`L2` column mappings, ambiguity, count lengths,
   duplicate rsIDs, finite values, and frequency coverage.
6. Reconstruct binary and continuous all/common counts and overlaps, validate
   the approved tolerance, preserve agreeing legacy values, reconstruct missing
   baseline `.M`, and fail real conflicts.
7. Represent missing unpartitioned `.M` as JSON null and reject downstream
   all-count regression without affecting common-count regression.
8. Inner-join reference and weights by rsID, preserve reference coordinates,
   audit asymmetric membership, and accept coordinate disagreements with
   warnings.
9. Produce both `rsid` and `chr_pos` directories, including best-effort and
   required build inference, zero-based normalization, and coordinate-identity
   duplicate failures.
10. Validate current Parquet schemas, row-group metadata, count/overlap config
    consistency, hashes, provenance, and clean/nonclean issue sidecars.
11. Leave only diagnostics after failure and require overwrite on rerun.
12. Convert the provided local `eur_w_ld_chr` and `weights_hm3_no_hla` as a
    non-CI smoke test without modifying either source directory.

### Legacy sumstats tests

1. Mark `.sumstats` and `.sumstats.gz` as legacy and reject footerless Parquet.
2. Drop all duplicated source-rsID clusters.
3. Against an allele-aware panel, cover direct, complement, swap,
   swapped-complement, invalid, ambiguous, missing, incompatible, zero-candidate,
   and multi-candidate cases.
4. Verify `Z` and valid `FRQ` transformations and prove missing/invalid `FRQ`
   never drops an otherwise usable SNP.
5. Against an allele-unaware `rsid` or `chr_pos` panel, preserve source allele
   orientation and use panel coordinates.
6. For rg, harmonize two projected traits pairwise and cover mixed current and
   legacy inputs.
7. Verify h2, partitioned-h2, and rg emit stable aggregate logs and a complete
   header-only/nonempty legacy drop audit as appropriate.
8. Verify the same projection occurs through public `RegressionRunner`, not
   only CLI wrappers.

### End-to-end regression tests

1. Run `h2` and `rg` with a converted unpartitioned suite and legacy sumstats.
2. Run baseline-only `partitioned-h2` with converted counts and overlap.
3. Verify common-count metadata consistency is enforced downstream.
4. Verify missing all counts fail only an all-count request.
5. Verify regression never reads a legacy LD-score fragment directly.

### Cleanup and documentation tests

1. Confirm top-level help and command parsers include the converter.
2. Confirm private kernel emitters and obsolete reader tests are absent.
3. Confirm annotation output writes `BP` and can be read back into canonical
   `POS` by LDSC3.
4. Confirm `io-argument-inventory.md` lists both the converter and
   `build-gene-ldscore-index` and no longer claims private legacy emitters are a
   supported path.

## Scientific and operational constraints

- Legacy sumstats rsIDs are lookup keys, not canonical identity assertions.
- The canonical LD-score panel owns regression identity and output coordinates.
- No allele or frequency value crosses artifact ownership boundaries by
  imputation.
- LDSC2 common counts use strict `0.05` semantics and are never described as
  LDSC3-native inclusive counts.
- Annotation overlaps and marginal counts must describe the same reference SNP
  universe.
- Source suites are read-only and conversion is deterministic.
- Large annotation shards are processed chromosome-wise; overlap and count
  accumulation must not require concatenating the full genome-wide annotation
  matrix in memory.
- Hashing is streaming and deterministic.
- Diagnostics use stable reasons and bounded console/error examples; full
  available issue sets go to compressed sidecars.
- Canonical Parquet row order remains genomic so block-jackknife behavior is
  unchanged.

## Risks and resolved edge cases

- **Published suites may omit `.M`.** Unpartitioned all counts remain
  unavailable; baseline counts are reconstructed from full annotations.
- **Reference and weight SNP sets differ.** Their rsID inner join is expected
  and audited.
- **Coordinates can reflect different builds.** They are diagnostic outside
  the reference source because joins are rsID-based.
- **Annotation names are not semantic labels.** The converter trusts the
  advanced caller's baseline assertion and validates only structural
  completeness.
- **Continuous annotations break diagonal=count intuition.** Marginal sums and
  overlap products are computed and validated separately.
- **Count conflicts can indicate mixed releases or filters.** They fail rather
  than produce an internally inconsistent canonical artifact.
- **Footerless Parquet is theoretically recoverable.** It remains rejected
  because the approved compatibility need is legacy text, not provenance-free
  current artifacts.

There are no unresolved product, scientific, or architectural questions. Exact
error wording and private helper decomposition may be refined during
implementation without changing these contracts.
