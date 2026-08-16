# Build an exact gene LD-score index

Last updated on: 2026-08-16

For the mathematical construction of the disjoint atoms, stored operator, and
sufficient statistics—and the full downstream indexed-assembly derivation—see
[Exact gene LD-score index: mathematical algorithm](../../current/gene-ldscore-index-mathematics.md).

## Goal

Build one complete reusable index that contains the fixed baseline LD scores,
an embedded caller-supplied gene catalog, exact disjoint-gene atoms, and the sparse
operator needed to assemble later gene-list annotations. The index is an
offline artifact; it does not run regression.

## Supported scientific contract

- explicit hg19 and either base `rsid` or base `chr_pos` identity, with a
  PLINK BED/BIM/FAM reference;
- a 1 cM window by default;
- unpadded gene intervals (`--padding-bp 0`) and MHC gene exclusion by default;
- bundled HapMap3 regression SNP candidates by default, or one custom
  `--regression-snps-file`;
- `--exclude-regions mhc-and-centromeres` applied after candidate selection by
  default.

Baseline and PLINK rows need not be identical. The builder inner-joins them by
the configured effective identifier, as ordinary PLINK-backed `ldscore` does.
It drops every row in a duplicate effective-key group with a warning and audit
record, drops and reports rows found on only one side, fails on an empty
intersection, and uses PLINK metadata and
genotypes for matched rows. Coordinate disagreement under a shared rsID is a
warning; PLINK coordinates win. Verify independently that PLINK and gene/map
inputs are hg19 because rsID matching cannot establish build.

### Choose the explicit identity mode

The builder has no identity or genome-build default. Always pass
`--snp-identifier rsid --genome-build hg19` or
`--snp-identifier chr_pos --genome-build hg19`. It does not accept `auto`, hg38,
or allele-aware modes and performs no liftover. Advanced users are responsible
for ensuring every baseline, PLINK, restriction, map, and catalog coordinate is
hg19.

`rsid` joins on `SNP`. `chr_pos` joins on normalized positive 1-based
`(CHR, POS)`; baseline and restriction SNP labels are ignored and a baseline
SNP column may be absent. PLINK supplies the published `CHR`, `POS`, `SNP`,
`A1`, and `A2`, so differing labels at a matched coordinate are reported but
the PLINK label is retained. Repeated restriction keys collapse because a
restriction is a set; duplicate mutable variant-source groups instead use
drop-all. Inspect `diagnostics/dropped_snps/chrN_dropped.tsv.gz`.

Fast indexed assembly inherits the stored identity mode and build and therefore
accepts neither live option, even when a supplied value would match.

### Regression SNP file format

`--regression-snps-file` must be a **headered text table containing SNP
identities**. Its required columns depend on the builder's explicit
`--snp-identifier`:

| `--snp-identifier` | Required columns | Optional columns |
|---|---|---|
| `rsid` | `SNP` | anything else |
| `chr_pos` | `CHR`, `POS` | anything else |

The gene-index builder supports only these two base identity modes. The reader
accepts whitespace-, tab-, or comma-delimited plain text and `.gz` files.
Duplicate restriction keys collapse, and columns outside the active identity
schema are ignored. In `chr_pos` mode, `CHR` and `POS` must be hg19 coordinates.

For example, the `custom_regression_snps.tsv` used by a
`--snp-identifier chr_pos --genome-build hg19` build can contain:

```tsv
CHR	POS
1	10583
1	13302
2	21537
```

For an `rsid` build, use a headered one-column file instead:

```tsv
SNP
rs123
rs456
```

The intersected baseline/PLINK SNPs are the LD-reference contributor, count,
and overlap universe. Regression candidates and `--exclude-regions` select
only persisted output rows and `regression_ld_scores` contributors.

> **Caveat: gene-catalog regions do not restrict baseline LD scores.** For
> each supplied baseline column, the builder recomputes LD scores over the full
> retained baseline/PLINK reference intersection. A retained SNP contributes
> according to its value in that baseline column, whether or not it lies in a
> catalog gene region. For an all-ones baseline column, every retained
> reference SNP within the LD window contributes. Catalog intervals are
> used only to construct the disjoint atoms for focal gene-list annotations
> and an optional custom `gene_control`; they never redefine the baseline LD-reference
> universe. Do not prefilter the PLINK or baseline inputs to gene
> regions unless that narrower reference universe is intentionally the desired
> scientific input.

### Which baseline-suite files are used?

The builder uses baseline annotation shards as the matrix to which the new
PLINK LD operator is applied. 

At build time **it recomputes baseline LD scores** from the these baseline annotation shards using the PLINK data. It does not reuse or copy the precomputed LD scores (`.l2.ldscore.gz`) or counts (`.M`, `.M_5_50`) previously distributed beside those annotations.

| Baseline-suite component | Used by index builder? |
| --- | --- |
| `baseline.N.annot.gz` | **Yes.** Supplies the baseline annotation matrix. |
| `baseline.N.l2.ldscore.gz` or another precomputed `.ldscore.gz` | **No.** Baseline LD scores are recomputed from the selected PLINK genotypes and current index settings. |
| `.M` and `.M_5_50` | **No.** Annotation counts are recomputed over the retained reference universe. |
| PLINK BED/BIM/FAM | **Yes.** Supplies genotypes, variant identities, authoritative coordinates, alleles, and BIM cM values when used. |
| Existing regression-weight LD scores | **No.** `regression_ld_scores` is recomputed for the selected regression rows. |

This is intentional. Recalculation keeps the supplied baseline block, gene
operator, counts, overlaps, and regression weights consistent with the same
PLINK samples, genotype/MAF filtering, baseline/PLINK intersection, genetic
map, LD window, adjusted-$r^2$ implementation, and regression-row policy. A
legacy `.ldscore.gz` file may have been produced under different choices even
when it is distributed in the same baseline suite.

## Build the production index

```bash
RESOURCE_ROOT="/path/to/ldsc_resources"
INDEX_ROOT="/path/to/gene_ldscore_indexes"

ldsc build-gene-ldscore-index \
  --baseline-annot-sources "${RESOURCE_ROOT}/baseline/baseline.@.annot.gz" \
  --plink-prefix "${RESOURCE_ROOT}/plink/1000G.EUR.QC.@" \
  --output-dir "${INDEX_ROOT}/baseline_100kb" \
  --gene-coordinate-file "${RESOURCE_ROOT}/genes/gene-coordinates.hg19.tsv.gz" \
  --genome-build hg19 \
  --snp-identifier rsid \
  --ld-wind-cm 1.0 \
  --padding-bp 100000 \
  --gene-exclude-regions mhc \
  --exclude-regions mhc-and-centromeres \
  --threads 1
```

This command explicitly requests non-default 100 kb padding. Omit
`--padding-bp` to build the default unpadded index.

The catalog is required, one-based, hg19, and must be fully canonical. The
builder checks every row before atom work and embeds the complete catalog,
including intentionally MHC-excluded genes. Production construction always
covers autosomes 1–22; partial builds and public `--chromosomes` selection are
unsupported.

If validation fails, use the full issue table under the hidden build-state
directory to repair the upstream catalog transformation. See
[Gene-list diagnostics and repair](../../current/gene-list-diagnostics-and-repair.md)
for every column/reason and the recommended workflow.

To use custom regression SNPs while retaining the standard region subtraction:

```bash
ldsc build-gene-ldscore-index \
  --baseline-annot-sources "${RESOURCE_ROOT}/baseline/baseline.@.annot.gz" \
  --plink-prefix "${RESOURCE_ROOT}/plink/1000G.EUR.QC.@" \
  --output-dir "${INDEX_ROOT}/custom_regression_index" \
  --gene-coordinate-file "${RESOURCE_ROOT}/genes/gene-coordinates.hg19.tsv.gz" \
  --genome-build hg19 \
  --snp-identifier chr_pos \
  --regression-snps-file custom_regression_snps.tsv \
  --exclude-regions mhc-and-centromeres
```

Duplicate keys collapse; ordering and nonidentity columns do not affect the
scientific selection. Use `--exclude-regions none` only when intentionally
retaining all candidate rows.

Keep `--threads 1` until the target system has been profiled because each
chromosome worker can multiply peak memory. Use `--keep-indivs-file` for a
one-IID-per-row sample restriction, and `--genetic-map-hg19-sources` when BIM
cM values are uninformative.

## Output and identity

```text
<index-dir>/
    metadata.json
    gene_catalog.parquet
    diagnostics/
        build-gene-ldscore-index.json
        build-gene-ldscore-index.log
        dropped_snps/chr1_dropped.tsv.gz ... chr22_dropped.tsv.gz
    chromosomes/
        chr1/ ... chr22/

.<index-name>.build-state/        # hidden operational state
    build-gene-ldscore-index.lock
    build-gene-ldscore-index.log # running or failed only
    gene_coordinate_catalog_issues.tsv.gz # current catalog failure only
    history/
```

Each chromosome directory contains `baseline_rows.parquet`,
`baseline_statistics.npz`, `atoms.parquet`, `gene_to_atom.npz`,
`ldscore_operator.npz`, `atom_statistics.npz`, and component `metadata.json`.

One canonical `index_id` binds all scientific content and settings. Paths,
output names, compression, harmless restriction-file ordering, threads, and
batch sizes do not define identity. There are no suite or profile IDs.

## Rerun and failure behavior

- missing output: leave it absent until a complete index is published;
- empty output: leave it empty until a complete index is published;
- failed build: preserve the hidden live log without creating a partial index;
- valid existing index: require `--overwrite` and rebuild everything;
- nonempty invalid output: fail even with `--overwrite`;
- same absolute target already building: fail immediately and point to the
  stable live log.

`--overwrite` never updates components in place. The old valid index remains
loadable until a complete staged replacement passes reload validation. Failed
overwrites preserve the old scientific index. The next invocation recovers a
single recognized valid backup after an interrupted publication, but refuses
ambiguous recovery evidence. Once a replacement has been reload-validated,
failure to remove its transaction directory is a warning rather than a build
failure; the warning reports the retained path for later cleanup.

The hidden run transaction is created before chromosome computation. A worker
writes its payloads to a private temporary shard and atomically renames it to
`<stage>/<index-name>/chromosomes/chrN`; only then does the log say `Finished
chromosome N`. That message means the chromosome is durable inside the private
stage, while the public destination is still missing, empty, or serving its old
complete version. After all chromosomes finish, the coordinator writes shared
metadata in canonical chromosome order, reload-validates the complete stage,
and moves the already-written tree into place without recopying its payloads.

These durable internal shards are not resumable checkpoints. A failed or killed
run is never continued, and a retry never reuses its scientific payloads. A
recognized stage with no publication backup is removed best-effort; build a new
complete transaction for every retry or configuration change.

Monitor a running build with:

```bash
tail -f "${INDEX_ROOT}/.production.build-state/build-gene-ldscore-index.log"
```

On success, the closed log moves to
`${INDEX_ROOT}/production/diagnostics/build-gene-ldscore-index.log`. The log is
the lifecycle status authority. The JSON file is written for a successful
publication summary; no separate status file is needed.

## Validate from Python

```python
from ldsc import load_gene_ldscore_index

index = load_gene_ldscore_index("/path/to/gene_ldscore_indexes/production")
print(index.index_id)
print(index.chromosomes)
```

## Next step

Pass the same directory to `ldsc ldscore --gene-ldscore-index-dir`; see
[Calculate LD scores](../main-functionalities/ldscore-from-gene-list.md).
