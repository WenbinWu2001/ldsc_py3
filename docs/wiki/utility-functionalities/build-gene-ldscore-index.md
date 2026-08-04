# Build an exact gene LD-score index

Last updated on: 2026-08-04

For the mathematical construction of the disjoint atoms, stored operator, and
sufficient statistics—and the full downstream indexed-assembly derivation—see
[Exact gene LD-score index: mathematical algorithm](../../current/gene-ldscore-index-mathematics.md).

## Goal

Build one complete reusable index that contains the fixed baseline LD scores,
an embedded protein-coding catalog, exact disjoint-gene atoms, and the sparse
operator needed to assemble later gene-list annotations. The index is an
offline artifact; it does not run regression.

## Supported scientific contract

- hg19, rsID identity, and a PLINK BED/BIM/FAM reference;
- a 1 cM window by default;
- 100 kb gene padding and MHC gene exclusion by default;
- bundled HapMap3 regression SNP candidates by default, or one custom
  `--regression-snps-file`;
- `--exclude-regions mhc-and-centromeres` applied after candidate selection by
  default.

Baseline and PLINK rows need not be identical. The builder inner-joins them by
the configured effective identifier (rsID in v1), as ordinary PLINK-backed
`ldscore` does. It fails on duplicate effective IDs or an empty intersection,
drops and reports rows found on only one side, and uses PLINK metadata and
genotypes for matched rows. Coordinate disagreement under a shared rsID is a
warning; PLINK coordinates win. Verify independently that PLINK and gene/map
inputs are hg19 because rsID matching cannot establish build.

The intersected baseline/PLINK SNPs are the LD-reference contributor, count,
and overlap universe. Regression candidates and `--exclude-regions` select
only persisted output rows and `regression_ld_scores` contributors.

> **Caveat: protein-coding regions do not restrict baseline LD scores.** For
> each supplied baseline column, the builder recomputes LD scores over the full
> retained baseline/PLINK reference intersection. A retained SNP contributes
> according to its value in that baseline column, whether or not it lies in a
> protein-coding gene region. For an all-ones baseline column, every retained
> reference SNP within the LD window contributes. Protein-coding intervals are
> used only to construct the disjoint atoms for `gene_control` and focal
> gene-list annotations; they never redefine the baseline LD-reference
> universe. Do not prefilter the PLINK or baseline inputs to protein-coding
> regions unless that narrower reference universe is intentionally the desired
> scientific input.

## Prototype chromosome 22

```bash
RESOURCE_ROOT="/path/to/ldsc_resources"
INDEX_ROOT="/path/to/gene_ldscore_indexes"

ldsc build-gene-ldscore-index \
  --baseline-annot-sources "${RESOURCE_ROOT}/baseline/baseline.@.annot.gz" \
  --plink-prefix "${RESOURCE_ROOT}/plink/1000G.EUR.QC.@" \
  --output-dir "${INDEX_ROOT}/prototype_chr22" \
  --chromosomes 22 \
  --genome-build hg19 \
  --snp-identifier rsid \
  --ld-wind-cm 1.0 \
  --padding-bp 100000 \
  --gene-exclude-regions mhc \
  --exclude-regions mhc-and-centromeres \
  --threads 1
```

The prototype is a complete chromosome-22 index. It cannot be extended in
place; use a different output directory for the production chromosomes-1–22
index.

To use custom regression SNPs while retaining the standard region subtraction:

```bash
ldsc build-gene-ldscore-index \
  --baseline-annot-sources "${RESOURCE_ROOT}/baseline/baseline.@.annot.gz" \
  --plink-prefix "${RESOURCE_ROOT}/plink/1000G.EUR.QC.@" \
  --output-dir "${INDEX_ROOT}/custom_regression_index" \
  --regression-snps-file custom_regression_snps.tsv \
  --exclude-regions mhc-and-centromeres
```

Restriction files are identity-only. Duplicate keys collapse; ordering and
nonidentity columns do not affect the scientific selection. Use
`--exclude-regions none` only when intentionally retaining all candidate rows.

## Production build

Omit `--chromosomes` to build autosomes 1–22. Keep `--threads 1` until the
target system has been profiled because each chromosome worker can multiply
peak memory. Use `--keep-indivs-file` for a one-IID-per-row sample restriction,
and `--genetic-map-hg19-sources` when BIM cM values are uninformative.

## Output and identity

```text
<index-dir>/
    metadata.json
    gene_catalog.parquet
    diagnostics/
        build-gene-ldscore-index.json
        build-gene-ldscore-index.log
    chromosomes/
        chr1/ ... chr22/

.<index-name>.build-state/        # hidden operational state
    build-gene-ldscore-index.lock
    build-gene-ldscore-index.log # running or failed only
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
[Calculate LD scores](../main-functionalities/ldscore.md).
