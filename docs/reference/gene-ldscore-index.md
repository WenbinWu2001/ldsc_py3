# Exact gene LD-score indexes

Last updated on: 2026-10-27

`--workers` / `workers` defaults to one managed child process, including on HPC.
Every worker uses **1 BLAS thread**, in serial and parallel execution.
`joblib.cpu_count()` supplies the detected available CPU count; it does not measure idle CPUs or read the requested scheduler allocation.
Positive requests above that count fail before capping by available tasks; negative `-k` uses `max(1, detected_count + 1 - k)`.
LD scoring and index construction cap by chromosomes; regression caps by complete query-plus-baseline models and query batch size.
Zero, Boolean, and non-integer Python requests are invalid.
See [the execution policy](parallel-execution.md) for native-limit verification, memory costs, process startup, and failure handling.

An exact gene LD-score index moves the repeated PLINK calculation for one
baseline, reference panel, regression-row policy, and gene projection offline.
Online gene-list runs resolve genes against the embedded catalog, assemble
Boolean unions from stored disjoint atoms, and write an ordinary canonical
LD-score directory.

For the complete input-to-output derivation—including the builder factorization,
online matrix products, annotation counts, and overlap sufficient statistics—see
[Exact gene LD-score index: mathematical algorithm](methods/gene-ldscore-index-math.md).

## Build one complete index

`--workers` controls a reusable chromosome process pool during index building.
The default of one still uses a managed child; every worker uses 1 BLAS thread.
A [focused transport probe](../decisions/2026-10-11-parallel-execution.md#final-decisions-and-remaining-boundaries) measured the additional process memory and verified compact task transport before migration.

```bash
ldsc build-gene-ldscore-index \
  --baseline-annot-sources annotations/baseline.@.annot.gz \
  --plink-prefix reference/1000G.EUR.QC. \
  --output-dir indexes/1000G_EUR_Phase3_baseline_100kb \
  --gene-coordinate-file annotations/gene-coordinates.hg19.tsv.gz \
  --genome-build hg19 \
  --snp-identifier rsid \
  --ld-wind-cm 1.0 \
  --padding-bp 100000 \
  --gene-exclude-regions mhc \
  --regr-snps-exclude-regions mhc-and-centromeres
```

The supplied coordinate catalog is required and is the complete gene universe
embedded in the index. It uses the one-based format in
[Gene-list query input](formats/gene-lists.md). The builder validates every
catalog row as canonical before creating a transaction or constructing atoms;
any defect stops the build and writes the full repair table described in
[Gene-list diagnostics and repair](results/gene-list-diagnostics.md).

Index construction defaults to unpadded gene intervals (`--padding-bp 0`).
The example above explicitly builds a non-default 100 kb index, which is why
both its destination name and command record `100kb`/`100000`. Omit the flag
for the default unpadded index.

PLINK prefixes use the [shared content-based resolver](input-paths.md#plink-prefix-resolution). A plain stem such as `1000G.EUR.QC.` is sufficient; `@` is optional. The builder validates all selected trios and requires chromosomes 1-22 before annotation preparation and genotype computation. Failed PLINK validation writes a complete `plink_input_issues.tsv` repair audit under `.<output-name>.build-state/`, alongside the live build log.

The builder supports PLINK inputs in hg19 or hg38 and all four shared SNP identity modes.
`--genome-build` and `--snp-identifier` are required without defaults; there is no `auto`, inference, or liftover.
The gene catalog and applicable input metadata must agree with the explicit build after alias normalization.
Conflicting labels are never overwritten, and coordinates are never converted.
The default persisted regression/output
rows are bundled HapMap3 SNPs with MHC and centromere regions removed. Supply
`--regr-snps-file custom.snplist` to replace the HapMap3 candidate set;
`--regr-snps-exclude-regions` is still applied afterward. Its choices are `none`, `mhc`,
`centromeres`, and `mhc-and-centromeres`.
The former `--exclude-regions` alias is removed; use `--regr-snps-exclude-regions`.

The custom regression SNP file must be a **headered text table containing SNP
identities**. Its required columns depend on the builder's explicit
`--snp-identifier`:

| `--snp-identifier` | Required columns | Optional columns |
|---|---|---|
| `rsid` | `SNP` | anything else |
| `chr_pos` | `CHR`, `POS` | anything else |
| `rsid_allele_aware` | `SNP` | both `A1` and `A2`, or neither |
| `chr_pos_allele_aware` | `CHR`, `POS` | both `A1` and `A2`, or neither |

Restriction files without allele columns select by base identity in allele-aware modes.
With both allele columns, they use the shared allele-aware restriction rules; supplying only one is invalid.
Duplicate restriction keys collapse, and columns outside the active identity schema are ignored.
For coordinate restrictions, a matching build-specific position column takes precedence over generic `POS`.

For example, `custom_regression_snps.tsv` for
`--snp-identifier chr_pos --genome-build hg19` can contain:

```tsv
CHR	POS
1	10583
1	13302
2	21537
```

Then pass it as `--regr-snps-file custom_regression_snps.tsv`. For an
`rsid` build, use the same headered layout with one `SNP` column instead.

The LD-reference universe is the inner intersection of baseline and PLINK SNPs
under the configured identifier mode, matching ordinary PLINK `ldscore`
behavior. Every row in a duplicate effective-identity group is dropped, warned,
and recorded under `diagnostics/dropped_snps/`; no representative is selected.
Baseline-only and PLINK-only SNPs are dropped and counted.
PLINK supplies chromosome, position, alleles, cM, genotypes, output identity, and gene projection coordinates for consistent matches.
Conflicting coordinates for a shared rsID or incompatible usable allele pairs at a shared base key stop the build before chromosome workers run.
The repair table is `.<index-name>.build-state/baseline_plink_issues.tsv`, with affected keys, input coordinates/alleles, reasons, and sources.
Cross-chromosome rsID conflicts are checked using a temporary disk-backed join, preserving chromosome-bounded metadata memory.
Ordinary unmatched SNPs remain an inner intersection; rsID matching alone does not establish a genome build.

### SNP identity contract

`rsid` uses `SNP` as the effective join key. `chr_pos` uses exactly normalized
`(CHR, POS)` with positive 1-based positions; a baseline `SNP` column is
optional and never participates in coordinate matching. After either join,
PLINK supplies the published `CHR`, `POS`, `SNP`, `A1`, and `A2`. Thus differing
labels at one matched coordinate are allowed and reported, and the PLINK label
is authoritative.
`rsid_allele_aware` and `chr_pos_allele_aware` additionally use shared unordered, strand-complement-normalized allele pairs.
Missing, invalid/non-SNP, identical, or palindromic A/T and C/G pairs are removed with diagnostics; multiallelic base-key groups and duplicate effective identities are dropped in full.
Baseline annotations may omit both allele columns, in which case base-key matching recovers unambiguous reference alleles.
These rules can change contributors, counts, and SNP-count windows.

The completed `index_id` binds the literal mode and explicit hg19/hg38 build.
Online
indexed assembly does not rematch sources and rejects live `--snp-identifier`
and `--genome-build` options even when they equal the index. A matching direct
run must explicitly use the index's mode and coordinate build.

An individual chromosome may have zero regression rows after restriction and
region subtraction, which produces a warning. Public construction always
builds autosomes 1 through 22; `--chromosomes` is not a public option. Public
loading rejects partial coverage. Smaller coverage exists only as a private
test seam. Online coverage is assessed after resolution and explicit gene exclusions; an uncovered selected gene fails the batch under both policies. Uncovered genes are not rejected identifiers or zero-support genes. Independent missing/invalid chromosome shards are collected before aborting index preflight.

### LD windows

Specify exactly one of `--ld-wind-cm`, `--ld-wind-kb`, or `--ld-wind-snps`; there is no default.
Each describes the maximum distance on either side of a SNP using the same PLINK kernel as direct `ldscore`.
SNP-count windows use the retained SNP order after identity cleanup, the baseline/PLINK intersection, and genotype/MAF filtering.
cM windows require informative BIM CM values or matching `--genetic-map-hg19-sources` / `--genetic-map-hg38-sources`; kb and SNP-count windows need no map.
The builder rejects whole-chromosome windows and exposes no override.
Window unit and value are stored in the immutable index identity.

### Baseline LD-score contributor caveat

Gene-catalog regions do not restrict recomputation of the supplied
baseline LD scores. For baseline column $c$ and persisted row $j$, the builder
computes

$$
l_c(j)=\sum_{k\in\text{retained reference SNPs within the LD window}}
\widetilde r^2_{jk}A_{kc}.
$$

Thus, all SNPs in the retained baseline/PLINK intersection are eligible
contributors; the value $A_{kc}$ determines a SNP's contribution to column
$c$. Catalog intervals define only the disjoint atoms used for focal
gene-list annotations and an optional custom `gene_control`. They do not filter
the supplied baseline matrix or its LD-reference universe.

### Baseline-suite component usage

The builder consumes annotation values and reconstructs the numerical baseline
payload; it does not import the suite's precomputed LD-score or count files.

Before launching chromosome workers, [`_run_gene_ldscore_index_build()`](../../src/ldsc/gene_ldscore_index.py) calls shared [`prepare_annotation_sources()`](../../src/ldsc/_annotation_sources.py). This serial phase scans and validates the complete baseline suite, applies global duplicate cleanup using metadata alone, and copies retained values into private NumPy/Parquet shards inside the existing build transaction. Automatic tiles use each aligned group's width, so unrelated chromosome shards do not shrink them. `--workers` controls the subsequent chromosome calculations. The parallelism evaluation is complete: preparation remains serial because the measured speed gains did not justify the memory and implementation costs. See the [storage design](../dev/design/memory-and-cpu.md#bounded-source-preparation) for the measurements.

At INFO, the live build log marks annotation input reading, SNP identity checking and chromosome preparation, and completion with retained counts and elapsed preparation time. These milestones precede `Starting chromosome N`; preparation completion means annotation shards are ready, while computation and index publication still follow. The CM/MAF notice appears once per source file read. See [annotation preparation logging](../current/workflow-logging.md#annotation-preparation) for the exact messages and failure interpretation.

| Baseline-suite component | Used by index builder? |
| --- | --- |
| `baseline.N.annot.gz` | **Yes.** Supplies $A$, the baseline annotation matrix. |
| `baseline.N.l2.ldscore.gz` or another precomputed `.ldscore.gz` | **No.** The builder recomputes $L_A=PRA$. |
| `.M` and `.M_5_50` | **No.** Counts are recomputed over the retained reference and common-SNP universes. |
| PLINK BED/BIM/FAM | **Yes.** Supplies genotypes and PLINK-authoritative SNP metadata. |
| Existing regression-weight LD scores | **No.** The builder recomputes $w=PRp$. |

Recomputation guarantees that $L_A$, the atom operator $Y=PRH$, counts,
overlaps, and $w$ share the index's selected individuals, genotype and MAF
filters, reference intersection, map/window, adjusted-$r^2$ implementation,
and regression-row policy. Merely residing beside the annotation shards does
not make a legacy `.ldscore.gz`, `.M`, or `.M_5_50` file an index input.

The fixed baseline and weight payload is projected in one PLINK traversal: the
builder appends the binary regression selector $p$ after the supplied baseline
matrix $A$, evaluates $R[A\;p]$ once, then splits $L_A=PRA$ from $w=PRp$ on the
persisted rows. Construction of $Y=PRH$ remains separate and intentionally
evaluates bounded atom-column batches, resetting the genotype cursor for each
batch to keep peak memory independent of the total atom count.

## Artifact and identity contract

One output directory contains one complete immutable index:

```text
<index-dir>/
    metadata.json
    gene_catalog.parquet
    diagnostics/
        build-gene-ldscore-index.json
        build-gene-ldscore-index.log
        dropped_snps/
            chrN_dropped.tsv.gz
    chromosomes/
        chrN/
            metadata.json
            baseline_rows.parquet
            baseline_statistics.npz
            atoms.parquet
            gene_to_atom.npz
            ldscore_operator.npz
            atom_statistics.npz

.<index-name>.build-state/        # hidden operational state
    build-gene-ldscore-index.lock
    build-gene-ldscore-index.log # present while running or after failure
    gene_coordinate_catalog_issues.tsv.gz # present after catalog failure
    baseline_plink_issues.tsv    # present after SNP consistency failure
    history/                     # prior failed/interrupted attempts

.<index-name>.stage-<run-id>/     # hidden, private transaction while building
    .gene-index-publication.json
    <index-name>/chromosomes/chrN/
```

`index_id` is a canonical SHA-256 identity over scientific content and settings,
not file paths or output names. It covers normalized baseline data, PLINK BED
content and BIM metadata, selected IID order, canonical regression keys and
region policy, genetic map, gene catalog, chromosome coverage, build/identifier,
window and MAF rules, padding, and gene-region policy. Batching, threads,
logging, output path, and overwrite are resource or publication controls and do
not change the identity.

Root and chromosome metadata record `format_version: 3` and repeat the immutable identity mode and hg19/hg38 build.
Version 3 marks exact PLINK LD windows; version 2 operators used LDSC2's batch-rounded windows, so the loader rejects them and asks for a rebuild.
Each chromosome records `published_row_metadata_sha256` over ordered `CHR SNP POS A1 A2` metadata for regression rows only.
Both allele columns are required.
This one digest covers the effective SNP key in all four identity modes; the redundant `effective_identity_sha256` field is no longer written or required.
The strict loader retains its mode/build, duplicate-effective-identity, canonical-order, schema, and row-metadata checks.
These fingerprints identify declared inputs/settings and verify row metadata; they do not checksum numerical LD-score or sparse-matrix payloads.

Earlier index formats are unsupported and must be rebuilt completely with the current package.
There is no compatibility adapter, migration, or partial component reuse.
The semantic `index_id` retains content-based identity, now including the selected window unit/value and allele-aware identity where applicable.
See the [fingerprint rationale](../current/sha256-usage-decisions.md) for the retained hash roles.

Gene indexes do not support incremental updates, chromosome append, profile
addition, or common-layer reuse. Any input or configuration change requires a
complete new 1–22 build. Partial production indexes are unsupported. The
incremental appearance of chromosome shards in a private run stage is only a
memory and durability strategy; it is not restart, resume, checkpoint reuse,
or incremental index-update support.

Online assembly and direct gene-list workflows share `query_annotations.assess_gene_coverage()`, `gene_query_statuses()`, `gene_viability_errors()`, and `finalize_query_statuses()`. Indexed support comes from gene-to-atom counts; direct support comes from the prepared reference/annotation intersection. Shared finalization prunes query tables, counts, overlap blocks, chromosome results, and classification metadata together. Online outputs include `diagnostics/chromosome_scope.json`; invalid indexes produce `diagnostics/input_issues.tsv` in the requested LD-score output directory without altering the index.

## Output preflight, logging, and replacement

- A missing destination stays absent, and an existing empty destination stays
  empty, until successful publication.
- A valid existing index requires `--overwrite`, even if its `index_id` would
  be unchanged.
- A directory containing only the regular root `RUN_FAILED.txt` marker, optionally alongside recognized legacy build diagnostics, is reusable. Preflight leaves the marker in place; successful publication replaces the directory and clears it. Marker presence does not authorize replacing unrelated files, directories, or incomplete scientific artifacts.
- Diagnostic recognition checks all entries, including empty subdirectories and symlinks. Interrupted-publication recovery also accepts an owned marker-only destination when exactly one validated matching backup exists, and carries the marker into the recovered index until the current build succeeds. Recovery still refuses unrelated destination contents.
- Other nonempty invalid directories fail before chromosome computation, including with `--overwrite`.
- A failed build never publishes a partial scientific index. An authorized overwrite failure may leave a marker-only public directory. Graceful failure removes its marked private transaction best-effort; an interrupted transaction is never reused scientifically on retry.
- Two builders cannot target the same absolute directory concurrently.

While the build is running, the live log is
`<parent>/.<index-name>.build-state/build-gene-ldscore-index.log`, so `tail -f`
works without placing an open file inside the replaceable artifact. After a
successful publication, the handler closes and the complete log moves atomically
to `<index-dir>/diagnostics/build-gene-ldscore-index.log`. Failed logs stay in
hidden state, and a retry archives them under `.<index-name>.build-state/history/`.
The lifecycle footer is the status authority: `Started`
without a terminal `Finished` or `Failed` line indicates abrupt termination;
there is no separate status file. The JSON diagnostic is a successful-build
summary, not a live status record.

Before chromosome computation, the builder creates one marked run-specific
sibling stage. Each worker writes its large payloads in a private temporary
shard directory, closes them, and atomically renames that directory to
`<stage>/<index-name>/chromosomes/chrN`. `Finished chromosome N` is logged only
after this rename, so it means the internal shard is durable; it does not mean
chromosome N is public. Workers return compact evidence, and the coordinator
writes shared metadata, the catalog, diagnostics, and final `index_id` metadata
after all chromosomes finish in canonical order.

The complete staged index is then reload-validated and moved into place without
rewriting or copying its chromosome payloads. During overwrite, the old valid
index remains loadable until the replacement passes validation. A graceful
failure keeps the old scientific index and its prior success diagnostics while
leaving the failed attempt in hidden build state. After the replacement is
reload-validated, transaction cleanup is best-effort garbage collection: a
cleanup error warns with the retained builder-owned path but does not turn the
completed publication into a failed command. Recognized computation-only
stages are discarded on the next invocation, never resumed. Ambiguous backup
evidence fails rather than guessing.

## Assemble gene-list LD scores

```bash
ldsc ldscore \
  --gene-ldscore-index-dir indexes/1000G_EUR_Phase3_baseline_100kb \
  --query-annot-gene-list-sources 'gene_lists/*.txt' \
  --output-dir results/gene_ldscores
```

The index directory must be explicit. Indexed mode accepts gene lists, the
optional `--control-gene-list-file`, output/overwrite, and logging controls; it
does not accept live baseline, PLINK/R2, build, identity, window, map, region,
or regression-SNP overrides. A missing or corrupt index never triggers
discovery or direct-mode fallback.

The assembled LD-score directory is self-contained and records the inherited
mode/build plus `index_id`; downstream `h2`, `rg`, and `partitioned-h2` use the
ordinary summary-statistics identity, genome-build, and explicit downgrade
rules. The source index is not needed after successful assembly.

By default, indexed assembly adds no gene control. Supplying one existing
one-column gene list with `--control-gene-list-file` creates the binary baseline
annotation `gene_control` from those selected, padded genes. A focal gene-list
coefficient (`tau`) is then conditional on that chosen background and the other
baseline annotations. Overlapping, nested, duplicated, and alias-selected genes
use Boolean union. The stored operator is float64, preserves negative adjusted-r-squared
values, and is not clamped or epsilon-pruned.

The broader intersected PLINK/baseline universe supplies LD-score contributors,
counts, and overlaps. Only the configured regression restriction and region
policy select persisted rows and `regression_ld_scores` contributors.

See also the task-oriented [build guide](../wiki/utility-functionalities/build-gene-ldscore-index.md)
and [indexed LD-score guide](../user-guide/workflows/gene-list-ld-scores.md).

The build lifecycle retains preceding/current lock coordination in `_gene_index_build_lock` and recognized old-log archival in `_migrate_legacy_gene_index_build_state`. This protects destinations that may still have an earlier writer and preserves operational history. Active build-state locking, validated publication, and interrupted-publication recovery remain required; these are not unused runtime helpers.

## Bounded validation and query assembly

`LoadedGeneLDScoreIndex` holds shared catalog/settings, chromosome component paths, and compact support counts.
Validation releases each chromosome before loading the next.
Each assembly worker loads one operator for the active query batch, writes a lossless float64 Parquet fragment, and releases the operator before its next task.
`run_indexed_ldscore(..., query_batch_size=100, workers=1)` bounds outstanding raw query fragments and concurrent chromosomes; the worker pool persists across batches.
The parent combines compact counts, overlaps, and exact float64 variation summaries in chromosome order, then validates and writes one chromosome row group at a time without concatenating genome-wide query values.
Constant-query filtering precedes float32 output conversion, including queries that vary across chromosomes or only below float32 precision.
The active batch's scratch is removed before the next batch starts; privately staged final Parquet files accumulate until all required calculations succeed and publication begins.
Repeated index reads are the trade-off for the bounded batch lifecycle and must be included in end-to-end performance measurements.
One batch uses `ldscore.query.parquet`; multiple batches use the ordered `metadata.json.query_batches` manifest.
Saved calls return `LDScoreSource` with explicit uncached `read_queries(names)` reads.
Every worker uses 1 BLAS thread.
See [parallel execution](parallel-execution.md); no checkpoint or resume workflow is provided.
