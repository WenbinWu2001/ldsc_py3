# Exact gene LD-score indexes

Last updated on: 2026-10-31

`--workers` sets the number of workers (default one, each with one BLAS thread); CPU budgets, limits and failure handling are in [parallel execution](parallel-execution.md).

An exact gene LD-score index moves the repeated PLINK calculation for one baseline, reference panel, regression-SNP policy, and gene projection offline.
Online gene-list runs resolve genes against the embedded gene coordinate catalog, assemble Boolean unions from stored disjoint atoms, and write an ordinary LD-score directory.

For the complete input-to-output derivation—including the builder factorization,
online matrix products, annotation counts, and overlap sufficient statistics—see
[Exact gene LD-score index: mathematical algorithm](methods/gene-ldscore-index-math.md).

## Build one complete index

`--workers` controls a reusable pool of chromosome workers during index building.
The default of one still uses a worker; every worker uses 1 BLAS thread.
A [focused transport probe](../dev/validation/performance.md#gene-ld-score-index-transport-probe) measured the additional process memory and verified compact task transport before migration.

```bash
ldsc build-gene-ldscore-index \
  --snp-identifier rsid \
  --genome-build hg19 \
  --baseline-annot-sources annotations/baseline.@.annot.gz \
  --gene-coordinate-file annotations/gene-coordinates.hg19.tsv.gz \
  --padding-bp 100000 \
  --gene-exclude-regions mhc \
  --plink-prefix reference/1000G.EUR.QC. \
  --ld-wind-cm 1.0 \
  --output-dir indexes/1000G_EUR_Phase3_baseline_100kb \
  --regr-snps-exclude-regions mhc-and-centromeres
```

The supplied gene coordinate catalog is required and is the complete set of genes embedded in the index.
It uses the one-based format in [Gene-list query input](formats/gene-lists.md).
The builder validates every catalog row as canonical before creating the staged index or constructing atoms; any defect stops the build and writes the full repair table described in [Gene-list diagnostics and repair](results/gene-list-diagnostics.md).

Index construction defaults to unpadded gene intervals (`--padding-bp 0`).
The example above explicitly builds a non-default 100 kb index, which is why
both its destination name and command record `100kb`/`100000`. Omit the flag
for the default unpadded index.

PLINK prefixes use the [shared content-based resolver](input-paths.md#plink-prefix-resolution). A plain stem such as `1000G.EUR.QC.` is sufficient; `@` is optional. The builder validates all selected trios and requires chromosomes 1-22 before annotation preparation and genotype computation. Failed PLINK validation writes a complete `plink_input_issues.tsv` repair audit under `<index-dir>/.build/`, alongside the live build log.

The builder supports PLINK inputs in hg19 or hg38 and all four SNP identifier modes.
`--genome-build` and `--snp-identifier` are required without defaults; there is no `auto`, inference, or liftover.
The gene coordinate catalog and applicable input metadata must agree with the explicit build after alias normalization.
Conflicting labels are never overwritten, and coordinates are never converted.
The default regression SNPs, which are also the persisted output rows, are the bundled HapMap3 SNPs with MHC and centromere regions removed.
Supply `--regr-snps-file custom.snplist` to replace the HapMap3 candidate set; `--regr-snps-exclude-regions` is still applied afterward.
Its choices are `none`, `mhc`, `centromeres`, and `mhc-and-centromeres`.
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

The reference SNPs are the inner intersection of baseline and PLINK SNPs under the configured SNP identifier mode, matching ordinary PLINK `ldscore` behavior.
PLINK SNPs off the baseline's annotation SNP grid are excluded and counted in each chromosome's reference-SNP funnel, which the build stores in the index, records in `diagnostics/build-gene-ldscore-index.json` and reports with a WARNING ([reference SNPs on the annotation SNP grid](ldscore-snp-universes.md#reference-snps-on-the-annotation-snp-grid)).
Every row in a duplicate effective-identity group is dropped, warned, and recorded under `diagnostics/dropped_snps/`; no representative is selected.
Baseline-only and PLINK-only SNPs are dropped and counted.
PLINK supplies chromosome, position, alleles, cM, genotypes, output identity, and gene projection coordinates for consistent matches.
Conflicting coordinates for a shared rsID or incompatible usable allele pairs at a shared base key stop the build before chromosome workers run.
The repair table is `<index-dir>/.build/baseline_plink_issues.tsv`, with affected keys, input coordinates/alleles, reasons, and sources.
Cross-chromosome rsID conflicts are checked using a temporary disk-backed join, preserving chromosome-bounded metadata memory.
Its database is in `<index-dir>/.build/`, and each table's `SNP` index exists before rows are inserted, so the join writes no temporary sort file to the system temporary directory.
Ordinary unmatched SNPs remain an inner intersection; rsID matching alone does not establish a genome build.

### SNP identity contract

`rsid` uses `SNP` as the effective join key. `chr_pos` uses exactly normalized
`(CHR, POS)` with positive one-based positions; a baseline `SNP` column is
optional and never participates in coordinate matching. After either join,
PLINK supplies the published `CHR`, `POS`, `SNP`, `A1`, and `A2`. Thus differing
labels at one matched coordinate are allowed and reported, and the PLINK label
is authoritative.
`rsid_allele_aware` and `chr_pos_allele_aware` additionally use shared unordered, strand-complement-normalized allele pairs.
Missing, invalid/non-SNP, identical, or palindromic A/T and C/G pairs are removed with diagnostics; multiallelic base-key groups and duplicate effective identities are dropped in full.
Baseline annotations may omit both allele columns, in which case base-key matching recovers unambiguous alleles from the reference panel.
These rules can change the reference SNPs, counts, and SNP-count windows.

The completed `index_id` binds the literal mode and explicit hg19/hg38 build.
Online
indexed assembly does not rematch sources and rejects live `--snp-identifier`
and `--genome-build` options even when they equal the index. A matching direct
run must explicitly use the index's mode and coordinate build.
In Python, `run_ldscore(gene_ldscore_index_dir=...)` ignores the registered `GlobalConfig.snp_identifier` and `GlobalConfig.genome_build` and records the index values in the result's `config_snapshot`.

An individual chromosome may have zero regression SNPs after restriction and region subtraction, which produces a warning.
Public construction always builds autosomes 1 through 22; `--chromosomes` is not a public option.
Public loading rejects partial coverage.
Smaller coverage exists only as a private test seam.
Online coverage is assessed after resolution and explicit gene exclusions; an uncovered selected gene fails the batch under both policies.
Uncovered genes are not rejected identifiers or zero-support genes.
Independent missing/invalid chromosome shards are collected before aborting index preflight.

### LD windows

Specify exactly one of `--ld-wind-cm`, `--ld-wind-kb`, or `--ld-wind-snps`; there is no default.
Each describes the maximum distance on either side of a SNP using the same PLINK kernel as direct `ldscore`.
SNP-count windows use the retained SNP order after identity cleanup, the baseline/PLINK intersection, and genotype/MAF filtering.
cM windows require informative BIM CM values or matching `--genetic-map-hg19-sources` / `--genetic-map-hg38-sources`; kb and SNP-count windows need no map.
The builder rejects whole-chromosome windows and exposes no override.
Window unit and value are stored in the immutable index identity.

### Baseline reference-SNP caveat

Gene coordinate catalog regions do not restrict recomputation of the supplied baseline LD scores.
For baseline column $c$ and persisted row $j$, the builder computes

$$
l_c(j)=\sum_{k\in\text{retained reference SNPs within the LD window}}
\widetilde r^2_{jk}A_{kc}.
$$

Thus, all SNPs in the retained baseline/PLINK intersection are eligible reference SNPs; the value $A_{kc}$ determines a SNP's contribution to column $c$.
Catalog intervals define only the disjoint atoms used for gene-list query annotations and an optional custom `gene_control`.
They do not filter the supplied baseline matrix or its reference SNPs.

### Baseline-suite component usage

The builder consumes annotation values and reconstructs the numerical baseline
payload; it does not import the suite's precomputed LD-score or count files.

Before launching chromosome workers, [`_run_gene_ldscore_index_build()`](../../src/ldsc/_ldscore/gene_index/build.py) calls shared [`prepare_annotation_sources()`](../../src/ldsc/_annotation/sources.py). This serial phase scans and validates the complete baseline suite, applies global duplicate cleanup using metadata alone, and copies retained values into private NumPy/Parquet shards inside the existing build transaction. Automatic tiles use each aligned group's width, so unrelated chromosome shards do not shrink them. `--workers` controls the subsequent chromosome calculations. The parallelism evaluation is complete: preparation remains serial because the measured speed gains did not justify the memory and implementation costs. See the [storage design](../dev/design/memory-and-cpu.md#bounded-source-preparation) for the measurements.

At INFO, the live build log marks annotation input reading, SNP identity checking and chromosome preparation, and completion with retained counts and elapsed preparation time. These milestones precede `Starting chromosome N`; preparation completion means annotation shards are ready, while computation and index publication still follow. The CM/MAF notice appears once per source file read. See [annotation preparation logging](logs-and-outputs.md#annotation-preparation) for the exact messages and failure interpretation.

| Baseline-suite component | Used by index builder? |
| --- | --- |
| `baseline.N.annot.gz` | **Yes.** Supplies $A$, the baseline annotation matrix. |
| `baseline.N.l2.ldscore.gz` or another precomputed `.ldscore.gz` | **No.** The builder recomputes $L_A=PRA$. |
| `.M` and `.M_5_50` | **No.** All-SNP and common-SNP counts are recomputed over the retained reference SNPs. |
| PLINK BED/BIM/FAM | **Yes.** Supplies genotypes and PLINK-authoritative SNP metadata. |
| Existing regression-weight LD scores | **No.** The builder recomputes $w=PRp$. |

Recomputation guarantees that $L_A$, the atom operator $Y=PRH$, counts, overlaps, and $w$ share the index's selected individuals, genotype and MAF filters, reference intersection, map/window, unbiased-R² implementation, and regression-SNP policy.
Merely residing beside the annotation shards does not make a legacy `.ldscore.gz`, `.M`, or `.M_5_50` file an index input.

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
            reference_snp_funnel.parquet
    .build/                       # hidden; present while building or after a failure
        build-gene-ldscore-index.lock
        build-gene-ldscore-index.log # live log
        gene_coordinate_catalog_issues.tsv.gz # after a catalog failure
        baseline_plink_issues.tsv    # after an SNP consistency failure
        history/                     # earlier failed or interrupted attempts
        stage/                       # the index being built
```

`index_id` is a canonical SHA-256 identity over scientific content and settings, not file paths or output names.
It covers normalized baseline data, PLINK BED content and BIM metadata, selected IID order, canonical regression keys and region policy, genetic map, gene coordinate catalog, chromosome coverage, build/identifier, window and MAF rules, padding, and gene-region policy.
Batching, workers, logging, output path, and overwrite are resource or publication controls and do not change the identity.
The input fingerprints behind `index_id` (normalized baseline content, whole PLINK BED file bytes, parsed BIM metadata, selected individual IDs, sorted regression keys, and the genetic map) only define the identity; indexed scoring never re-reads those inputs to check them again.
Identical inputs and settings must produce identical fingerprints and the same `index_id`; changing what a fingerprint covers changes the identity of every index.

Root and chromosome metadata record `format_version: 4` and repeat the immutable SNP identifier mode and hg19/hg38 build.
Version 4 adds each chromosome's reference-SNP funnel, `reference_snp_funnel.parquet`, which indexed `ldscore` reports; version 3 introduced exact PLINK LD windows, and version 2 operators used legacy LDSC's batch-rounded windows.
The loader rejects every earlier version and asks for a rebuild.
Each chromosome records `published_row_metadata_sha256` over ordered `CHR SNP POS A1 A2` metadata for regression SNPs only.
Both allele columns are required.
This one digest covers the effective SNP identity in all four SNP identifier modes; the redundant `effective_identity_sha256` field is no longer written or required.
The strict loader retains its mode/build, duplicate-effective-identity, canonical-order, schema, and row-metadata checks.
These fingerprints identify declared inputs/settings and verify row metadata; they do not checksum numerical LD-score or sparse-matrix payloads.

Earlier index formats are unsupported and must be rebuilt completely with the current package.
There is no compatibility adapter, migration, or partial component reuse.
The semantic `index_id` retains content-based identity, now including the selected window unit/value and allele-aware identity where applicable.
Together with the R² matrix sidecar binding ([R² matrix format](formats/r2-matrices.md#4-index--sidecar-binding)), these are the package's only hashes; routine provenance hashes are not computed.

Gene LD-score indexes do not support incremental updates, chromosome append, profile addition, or common-layer reuse.
Any input or configuration change requires a complete new 1–22 build.
Partial production indexes are unsupported.
The incremental appearance of chromosome shards in a private run stage is only a memory and durability strategy; it is not restart, resume, checkpoint reuse, or incremental index-update support.

Online assembly and direct gene-list workflows share `_annotation.query_outcomes.assess_gene_coverage()`, `gene_query_statuses()`, `gene_viability_errors()`, and `finalize_query_statuses()`. Indexed support comes from gene-to-atom counts; direct support comes from the prepared reference/annotation intersection. Shared finalization prunes query tables, counts, overlap blocks, chromosome results, and classification metadata together. Online outputs include `diagnostics/chromosome_scope.json`; invalid indexes produce `diagnostics/input_issues.tsv` in the requested LD-score output directory without altering the index.

## Output preflight, logging, and replacement

Every file the build creates, writes, renames or deletes is inside the index directory.
The lock, the live log, issue tables, earlier failed attempts and the index being built are in its hidden `.build/` directory, which a successful build removes.

- A missing destination is created when the build starts; until a build succeeds it holds only `.build/` and, after a failed authorized overwrite, `RUN_FAILED.txt`.
- A valid existing index requires `--overwrite`, even if its `index_id` would be unchanged.
- A directory holding only `.build/` and a regular `RUN_FAILED.txt` is reused without `--overwrite`.
  Preflight leaves the marker in place; a successful build removes it.
- Index components beside `.build/` without a loadable index are the remains of an interrupted build or publication; rebuilding over them needs `--overwrite`.
- Any other content fails before chromosome computation, including with `--overwrite`: unrelated files or directories, a symlink, a marker that is not a regular file, or index component names without `.build/`.
- Two builders cannot target the same directory concurrently.

While the build is running, the live log is `<index-dir>/.build/build-gene-ldscore-index.log`, so `tail -f` works.
After a successful publication, the handler closes, the complete log moves to `<index-dir>/diagnostics/build-gene-ldscore-index.log`, and `.build/` is removed.
Failed logs stay in `.build/`, and a retry archives them under `.build/history/`.
The lifecycle footer is the status authority: `Started` without a terminal `Finished` or `Failed` line indicates abrupt termination; there is no separate status file.
The JSON diagnostic is a successful-build summary, not a live status record.

Before chromosome computation, the builder creates the staged index in `.build/stage/`.
Each worker writes its large payloads in a private temporary shard directory, closes them, and atomically renames that directory to `.build/stage/chromosomes/chrN`.
`Finished chromosome N` is logged only after this rename, so it means the internal shard is durable; it does not mean chromosome N is public.
Workers return compact evidence, and the coordinator writes shared metadata, the catalog, diagnostics, and final `index_id` metadata after all chromosomes finish in canonical order.

The complete staged index is then reload-validated and published without copying its chromosome payloads.
Publication removes the old root `metadata.json`, which records `index_id`, then replaces each other component with its staged counterpart, merging `diagnostics/` entry by entry and keeping unrelated files, and writes the new `metadata.json` last.
The published index is reload-validated again.
The loader rejects a directory without `metadata.json`, so an interrupted publication is never read as an index.
A failure before publication begins keeps the earlier index complete and loadable; a failure during publication leaves no loadable index, which the failure marker and the loader report.
Leftovers of a failed attempt are discarded on the next invocation, never resumed.

## Assemble gene-list LD scores

```bash
ldsc ldscore \
  --query-annot-gene-list-sources 'gene_lists/*.txt' \
  --output-dir results/gene_ldscores \
  --gene-ldscore-index-dir indexes/1000G_EUR_Phase3_baseline_100kb
```

The index directory must be explicit.
Indexed mode accepts gene lists, the optional `--control-gene-list-file`, output/overwrite, and logging controls; it does not accept live baseline, PLINK/R², build, identity, window, map, region, or regression-SNP overrides.
A missing or corrupt index never triggers discovery or direct-mode fallback.

The assembled LD-score directory is self-contained and records the inherited
mode/build plus `index_id`; downstream `h2`, `rg`, and `partitioned-h2` use the
ordinary summary-statistics identity, genome-build, and explicit downgrade
rules. The source index is not needed after successful assembly.

By default, indexed assembly adds no gene control.
Supplying one existing one-column gene list with `--control-gene-list-file` creates the binary baseline annotation `gene_control` from those selected, padded genes.
The coefficient (column `coefficient` in `partitioned_h2.tsv`) of a gene-list query annotation is then conditional on that chosen background and the other baseline annotations.
Overlapping, nested, duplicated, and alias-selected genes use Boolean union.
The stored operator is float64, preserves negative unbiased R² values, and is not clamped or epsilon-pruned.

The broader intersection of PLINK and baseline SNPs supplies the reference SNPs for LD scores, counts, and overlaps.
Only the configured regression restriction and region policy select the regression SNPs, which are the persisted rows and the SNPs that contribute to `regression_ld_scores`.

See also the task-oriented [build guide](../user-guide/workflows/gene-list-ld-scores.md)
and [indexed LD-score guide](../user-guide/workflows/gene-list-ld-scores.md).

`_gene_index_build_lock` holds one nonblocking lock, `<index-dir>/.build/build-gene-ldscore-index.lock`, for each destination.
A successful build removes `.build/` while still holding the lock, so taking the lock checks that the locked file is still the one at the lock path and retries otherwise.
The builder neither reads nor migrates the hidden sibling `.<index-name>.build-state/` and `.<index-name>.stage-*/` directories, the visible `<index-name>.build/` log directory, or the standalone `.<index-name>.build-gene-ldscore-index.lock` of earlier releases.
Build-directory locking and validated publication remain required; these are not unused runtime helpers.

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
