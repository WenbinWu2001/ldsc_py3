# Legacy LDSC2 LD-Score Suite Conversion

Last updated on: 2026-08-04

This document defines the interoperability boundary for explicitly converting
selected reusable LDSC2 LD-score suites into the canonical LDSC3 LD-score
result-directory contract.

## Command

Conversion is exposed as:

```text
ldsc convert-ldsc2-ldscores \
    --legacy-reference-dir REF_DIR \
    --legacy-weight-dir WEIGHT_DIR \
    --output-dir OUTPUT_DIR
```

Baseline partitioned conversion additionally requires
`--legacy-frequency-dir FRQ_DIR`. Users who reuse one unpartitioned suite for
both reference scores and regression weights pass the same directory to both
arguments. Optional `--snp-identifier rsid|chr_pos` and
`--genome-build auto|hg19|hg38` controls select the allele-unaware converted
identity; `rsid` and `auto` are the defaults. The command normally needs no
glob, prefix, column-map, or suite-profile argument.

## Compatibility Boundary

Legacy LD-score files are not accepted directly by `h2`, `partitioned-h2`, or
`rg`. A user explicitly runs a one-time conversion command, then supplies the
new output directory to LDSC3 regression commands. Conversion never changes,
renames, or deletes the source files.

Conversion accepts a reference LD-score suite and a one-column regression-weight
LD-score suite. The two inputs need not contain identical SNP rows. The
converted regression rows are their inner join by rsID. Users may supply the
same suite for both roles, as is customary for `eur_w_ld_chr` in
non-partitioned regression.

Two reference-suite profiles are supported:

1. An ordinary one-column LDSC2 suite such as `eur_w_ld_chr`. It consists of
   chromosome-sharded `.l2.ldscore` or `.l2.ldscore.gz` tables and count files.
   The converted directory has one baseline LD-score column named `base` and
   is consumed by `h2` and `rg`.
2. A baseline-only partitioned suite such as `baseline_v1.2`. It consists of
   matching chromosome-sharded `.annot` or `.annot.gz`, `.l2.ldscore` or
   `.l2.ldscore.gz`, `.l2.M`, and `.l2.M_5_50` families. The converted directory
   retains the input baseline annotation/LD-score column family and is consumed
   by `partitioned-h2` in its baseline-only functional-category regime. All
   imported baseline columns are fitted jointly; the converted result contains
   no query columns and does not invoke the baseline-plus-one-query regime.

The converter detects these profiles without a profile flag. One reference
LD-score value column with no annotation family is unpartitioned. A complete
full annotation family with a bijective annotation/LD-score column mapping is
baseline partitioned. Multiple reference LD-score columns without annotations,
or an incomplete annotation family, are errors. LDSC3 validates structural
completeness and internal consistency but does not attempt to infer whether the
annotation columns are biologically a baseline rather than a query model; that
advanced-use assertion belongs to the caller. No all-ones root annotation is
required.

The inputs are suite directories and the output is a new canonical LDSC3
LD-score directory. Routine conversion must not require users to enumerate
chromosome files or provide low-level column mappings.

The following artifacts and workflows are outside this compatibility
boundary:

- cell-type, query, or other non-baseline annotation-specific LD scores;
- thin annotations;
- annotation overlap artifacts; and
- automatic, implicit conversion during a regression command.

Conversion never accepts query annotations. It does not turn the LDSC2 input
directory into a general LDSC3 annotation source or permit later query
extension of the imported suite.

Consequently, “baseline partitioned conversion” does not mean conversion of a
complete historical cell-type S-LDSC analysis. It imports only the reusable
baseline reference model. Users who need a query/cell-type model must compute
that query through the current LDSC3 annotation and LD-score contracts rather
than combining it with the converted directory.

## Counts

Every converted suite requires `.l2.M_5_50` files. For an unpartitioned suite,
these files are the authoritative common-reference-SNP counts because the suite
does not contain the LD-reference SNP annotation matrix from which they could
be reconstructed. `.l2.M` files are optional. When they are absent from an
unpartitioned suite, the converted metadata represents the all-reference-SNP
count as missing, conversion emits a warning, and downstream regression rejects
an explicit request for all-SNP counts. LDSC3 never substitutes the common-SNP
count for the missing all-SNP count.

A baseline partitioned suite contains the full SNP-level annotation matrix, and
its required frequency suite supplies the common-SNP mask. Baseline conversion
therefore reconstructs overlap information and can independently calculate
annotation count vectors. Source `.l2.M` and `.l2.M_5_50` files participate in
validation of those reconstructed values under the baseline-suite policy.

Baseline counts are reconstructed and validated per chromosome and annotation
before aggregation. When a source count is present and agrees exactly for
integer/binary annotations or within the documented floating-point tolerance
for continuous annotations, the converted metadata preserves the exact legacy
value. When `.l2.M` is absent, conversion uses the reconstructed all-SNP count.
A disagreement beyond tolerance is a conversion error: LDSC3 does not combine
legacy marginal counts with an overlap matrix calculated from a conflicting
annotation/frequency universe, and there is no mismatch override. This
reconstruction and validation are unavailable for unpartitioned suites, whose
available legacy count files remain authoritative.

Reconstruction uses LDSC2 common-variant semantics: `MAF > 0.05` for a minor
allele frequency column, or `0.05 < FRQ < 0.95` for an allele-frequency column.
The imported metadata records this legacy operator rather than describing the
counts as LDSC3-native `MAF >= 0.05` counts. The compatibility threshold is
fixed at `0.05` for both conversion profiles. The converter exposes no
`--common-maf-min` option and its Python API has no configurable threshold.

## Identity and Coordinates

Legacy suites are allele-unaware. Conversion does not infer, impute, or attach
reference and alternate alleles. The converted result uses an allele-unaware
SNP identifier mode. The reference LD-score suite supplies the converted
`CHR`, `POS`, and `SNP` fields; the weight suite supplies regression LD-score
values after the rsID inner join.

All relationships among legacy reference LD-score, annotation, frequency, and
weight families are matched by rsID. Reference LD-score `CHR` and `BP` become
canonical output `CHR` and `POS`. Coordinates in the other families are not
join keys and disagreement does not discard an otherwise unique rsID match.
Instead, coordinate disagreement produces a warning and a diagnostic record.

The converter may infer a genome build separately from each legacy source
family for validation and provenance. It does not alter the legacy directory
or perform a liftover. In `chr_pos` output mode, reference LD-score coordinates
must form unique canonical identities and the reference build must resolve;
there is no strict cross-family `CHR`/`POS` equality requirement.

For `rsid` output, canonical `genome_build` remains null. Automatic reference
build inference is best-effort provenance and an unresolved result warns but
does not fail; an explicit build validates or declares source provenance rather
than changing rsID identity. For `chr_pos` output, the reference build must be
resolved automatically or explicitly, detected zero-based coordinates are
normalized to canonical one-based `POS`, and the resolved build becomes
canonical identity metadata.

## Column Names

For a baseline partitioned suite, annotation and LD-score columns must have an
unambiguous one-to-one correspondence. An LD-score column may either equal its
annotation column or append one terminal `L2` suffix. The converted scientific
column name is always the name present in the LD-score table; conversion does
not remove or add `L2`. Counts and overlap axes use those same converted
LD-score column names.

## Legacy Table Schemas

Reference and weight LD-score tables require `CHR`, `SNP`, and `BP`. Legacy
`CM` and `MAF` metadata columns are accepted and excluded from the scientific
LD-score columns. Weight tables require exactly one numeric LD-score column;
unpartitioned reference tables likewise require exactly one and emit it as
canonical `base`. Baseline reference tables may contain multiple numeric
LD-score columns.

Full baseline annotations require `CHR`, `BP`, `SNP`, and `CM` followed by
numeric annotation columns. Thin annotations are rejected. Missing, nonnumeric,
or non-finite scientific values are conversion errors.

The baseline name resolver first tries exact LD-score-to-annotation equality,
then, for an LD-score name with one terminal `L2`, tries the stripped name. The
complete mapping must be unambiguous and bijective. Converted output always
preserves the exact LD-score table name.

The frequency family requires `SNP` and exactly one usable `MAF` or `FRQ`
column; optional allele and sample metadata are ignored. Frequencies are
matched to annotations by rsID. Every annotation SNP must have exactly one
usable frequency. Missing or duplicate annotation-frequency matches are
conversion errors; extra frequency SNPs are ignored and audited.

## Suite Discovery

The converter discovers `CHR.l2.*` and `PREFIX.CHR.l2.*` chromosome families
inside each supplied directory and requires one coherent autosomal 1--22
family for each role. Matching annotations and counts use the selected
reference prefix. Duplicate compressed and uncompressed representations are
accepted only when their decompressed contents agree. Non-family alternatives
such as `6_old.*` are ignored and recorded; ambiguous families are rejected.

Reference and weight LD-score rows are inner-joined by rsID within chromosome,
then validated for global uniqueness. Reference values and coordinates are
preserved, and the weight suite contributes only canonical
`regression_ld_scores`. Reference-only, weight-only, and retained rsIDs are
audited per chromosome. There is no conversion retention-percentage cutoff;
an empty final intersection is an error and ordinary low-SNP scientific
warnings remain regression concerns.

## Conversion Diagnostics

Conversion errors state the violated invariant, source role and file,
chromosome or annotation when applicable, affected-SNP count, and a bounded
sample of causal rsIDs. The complete row-level set is written to
`diagnostics/conversion_issues.tsv.gz`. Count-vector conflicts, whose scalar
source files do not identify causal SNP membership, instead report source and
reconstructed values per chromosome and annotation; frequency-boundary SNPs
are additionally reported when they could explain a common-count discrepancy.

The command creates the requested output directory's `diagnostics/` area while
it validates the sources. It writes `convert-ldsc2-ldscores.log` and
`conversion_issues.tsv.gz`, including a header-only issue table after a clean
conversion. Root `metadata.json` and parquet artifacts are written only after
all validation succeeds. A diagnostics-only failed directory is not a
consumable LD-score result and a corrected rerun into it requires
`--overwrite`.

## Module Boundary

LDSC2 discovery, parsing, validation, and conversion belong to a small public
workflow-layer `LegacyLDScoreConverter`, with a convenience
`convert_ldsc2_ldscores()` function and structured conversion result. It reuses
the pure overlap calculation and canonical LD-score directory writer. No
legacy LD-score reader is added to regression and no legacy file parsing or
emission is added to `_kernel`; after conversion, the ordinary canonical
`load_ldscore_from_dir()` path consumes the result.

## Provenance

The converted directory must identify itself as an imported LDSC2 suite and
record enough source and discovery provenance to make the conversion
auditable. It remains a normal current `artifact_type="ldscore"` directory and
adds a `legacy_ldsc2_import` metadata block containing the profile, converter
version, source directories and selected prefixes, selected/ignored/missing
files and hashes, intersection counts, source/effective build inference, count
origins, fixed legacy common-frequency rule, coordinate disagreement summaries,
and diagnostic paths. Original legacy files remain untouched.

For unpartitioned conversion, the `base` count record stores
`all_reference_snp_count: null` when any chromosome lacks `.l2.M`, plus the
aggregated required common count. Provenance names the chromosomes missing
`.l2.M`. Common-count regression remains available; all-count regression raises
a targeted error before dataset assembly. JSON `NaN` is never emitted.

For baseline conversion, all- and common-universe overlap matrices are computed
as `A.T @ A`, with the common matrix restricted by the fixed legacy frequency
mask. Both axes are relabeled with exact LD-score names and stored as the full
baseline-by-baseline block. Marginal counts remain annotation column sums, not
overlap diagonals. Baseline conversion writes no query columns or query parquet;
unpartitioned conversion writes no overlap artifact.

## Related Contracts

- [`legacy-sumstats-compatibility.md`](legacy-sumstats-compatibility.md)
  defines the separate automatic compatibility boundary for legacy munged
  sumstats.
- [`io-argument-inventory.md`](io-argument-inventory.md) defines public command
  inputs and outputs.
- [`data-flow.md`](data-flow.md) defines the canonical LD-score result-directory
  contract consumed by regression.
