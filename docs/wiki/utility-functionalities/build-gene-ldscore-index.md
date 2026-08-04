# Build an exact gene LD-score index

Last updated on: 2026-08-03

## Motivation

Many biological hypotheses arrive as gene sets: genes differentially expressed
in a tissue or cell type, genes prioritized by proteomics, or genes grouped by
GO or SynGO terms. Stratified LD score regression (S-LDSC) tests whether SNP
heritability is enriched near those genes while controlling for broader
functional annotations.

LDSC-SEG established this pattern for gene-expression data: it converted each
set of specifically expressed genes into a genomic annotation with 100 kb
windows, then tested that annotation conditional on the baseline model and an
all-genes annotation. This allowed genome-wide polygenic signal—not only
genome-wide significant loci—to identify relevant tissues and cell types
([Finucane et al., 2018](https://doi.org/10.1038/s41588-018-0081-4); [local
paper](../../../../docs/ldsc_papers/paper_ldsc-seg.pdf)).

The expensive part of repeating this analysis for many gene sets is computing
LD scores from the same reference panel. An exact gene LD-score index performs
that PLINK calculation once for a fixed baseline/reference/profile combination.
Later gene lists reuse the stored operator while preserving the direct S-LDSC
calculation exactly.

## Goal

Build a reusable index suite containing:

- immutable common data for one baseline suite and PLINK reference panel;
- an embedded protein-coding gene catalog and one projection profile;
- exact disjoint-atom operators used to assemble arbitrary gene-set unions.

This is an offline utility step. It does not run a GWAS regression and it does
not package the resulting index with the Python distribution.

## Initial supported configuration

The v1 builder is intentionally strict:

- hg19 coordinates and rsID SNP identity;
- the `1000G_EUR_Phase3` PLINK suite;
- a 1 cM LD window;
- bundled HapMap3 regression rows minus MHC and centromeres;
- no retained-reference MAF filter by default;
- common-SNP counts at inclusive MAF 0.05;
- 100 kb gene padding and MHC gene exclusion by default.

Before genotype QC, every baseline chromosome must contain exactly the same
`CHR/POS/SNP` identities as the corresponding PLINK BIM. Reordering is allowed;
missing, extra, duplicate, or conflicting identities stop the build.

## Keep every coordinate-bearing input on hg19

The v1 index has one coordinate build: **hg19**. It does not lift coordinates
between builds, and it cannot write an hg38 index from hg19 inputs. The
following components must all use the same hg19 coordinate system:

- baseline annotation `CHR/POS` rows;
- PLINK BIM positions;
- the embedded gene-catalog projection;
- bundled HM3 regression rows and MHC/centromere masks;
- an explicit genetic map, when supplied.

An explicit map must therefore use `--genetic-map-hg19-sources`; an hg38 map is
rejected. Without an explicit map, informative BIM cM values are used and must
correspond to the hg19 BIM positions.

Exact baseline/BIM `CHR/POS/SNP` equality detects mismatched coordinate sets,
but it cannot prove that two mutually matching files are truly hg19. Confirm
the documented source build of both suites before construction. Two hg38 files
misdeclared as hg19 could agree with each other while projecting hg19 genes and
region masks onto the wrong coordinates.

## Distinguish reference SNPs from regression SNPs

The v1 builder fixes the **written regression rows** to bundled HapMap3 (HM3)
SNPs minus MHC and centromeric regions. It does **not** build LD scores from
HM3 SNPs alone.

- LD-score contributors are all retained SNPs in the PLINK reference panel,
  including non-HM3, MHC, and pericentromeric SNPs.
- `M`, `M_5_50`, annotation counts, and overlaps use that same broad retained
  PLINK universe; only the common-MAF rule further restricts `M_5_50`.
- Canonical baseline/query output rows are the bundled HM3 set minus MHC and
  centromeres, intersected with retained PLINK SNPs.
- `regression_ld_scores` (`w_ld`) uses that filtered regression set as both its
  rows and contributors.

Therefore, “the index supports only HapMap3 SNPs” is incorrect. The precise v1
restriction is that the regression-row policy is fixed to bundled filtered HM3;
a custom regression-SNP set is not supported in indexed mode. Use direct mode
when a different regression-row policy is required.

## Set input paths

```bash
RESOURCE_ROOT="/path/to/ldsc_resources"
BASELINE_ANNOT_SOURCES="${RESOURCE_ROOT}/1000G_EUR_Phase3_baseline/baseline.@.annot.gz"
PLINK_PREFIX="${RESOURCE_ROOT}/1000G_EUR_Phase3_plink/1000G.EUR.QC.@"

INDEX_ROOT="/path/to/gene_ldscore_indexes"
INDEX_SUITE_DIR="${INDEX_ROOT}/1000G_EUR_Phase3_baseline"
```

`@` is replaced by the chromosome number. The PLINK prefix must resolve matching
`.bed`, `.bim`, and `.fam` files.

## Run a chromosome-22 prototype first

**Recommended memory allocation:** 8 GB with one worker

**Observed running time:** about 2 minutes for chromosome 22 in the local
validation; larger chromosomes may take longer and use more memory.

```bash
INDEX_PROTOTYPE_DIR="${INDEX_ROOT}/prototype_chr22_baseline"

ldsc build-gene-ldscore-index \
  --baseline-annot-sources "${BASELINE_ANNOT_SOURCES}" \
  --plink-prefix "${PLINK_PREFIX}" \
  --output-dir "${INDEX_PROTOTYPE_DIR}" \
  --chromosomes 22 \
  --genome-build hg19 \
  --snp-identifier rsid \
  --ld-wind-cm 1.0 \
  --padding-bp 100000 \
  --gene-exclude-regions mhc \
  --threads 1
```

A chromosome subset is a complete suite with its own identity. It cannot later
be extended in place, so use a separate output directory for the prototype.
Check:

```text
<prototype>/profiles/padding-100000bp-mhc/diagnostics/
    build-gene-ldscore-index.log
    build-gene-ldscore-index.json
```

The JSON records wall time, peak RSS, batch settings, row counts, atom count,
operator nonzeros, and payload bytes.

## Build the chromosomes 1–22 suite

Omit `--chromosomes` to select autosomes 1–22:

```bash
ldsc build-gene-ldscore-index \
  --baseline-annot-sources "${BASELINE_ANNOT_SOURCES}" \
  --plink-prefix "${PLINK_PREFIX}" \
  --output-dir "${INDEX_SUITE_DIR}" \
  --genome-build hg19 \
  --snp-identifier rsid \
  --ld-wind-cm 1.0 \
  --padding-bp 100000 \
  --gene-exclude-regions mhc \
  --threads 1
```

Keep `--threads 1` for the first production run. Each additional chromosome
worker can multiply chromosome-local memory. Increase workers only after
measuring peak RSS on the target system.

If BIM cM coordinates are absent or uninformative, provide an explicit matching
hg19 map through `--genetic-map-hg19-sources`. To restrict samples, provide a
one-IID-per-row file through `--keep-indivs-file`. Both choices become part of
the immutable suite identity.

## Output directory

```text
1000G_EUR_Phase3_baseline/
    metadata.json
    common/
        metadata.json
        chr1/ ... chr22/
    profiles/
        padding-100000bp-mhc/
            metadata.json
            gene_catalog.parquet
            chr1/ ... chr22/
            diagnostics/
                build-gene-ldscore-index.log
                build-gene-ldscore-index.json
```

`suite_id` binds the baseline, PLINK files, selected individuals, map/window,
SNP policies, and chromosome coverage. `profile_id` additionally binds the gene
catalog, projection build, padding, and gene-region policy.

The builder stages and reload-validates scientific payloads before publication.
`--overwrite` replaces only the named profile, preserves sibling profiles, and
reuses `common/` only when the suite identity matches. It never converts a
conflicting suite into a new one.

## Validate a distributed profile

The indexed LD-score command validates the whole profile automatically. For a
standalone installation check:

```bash
python - <<'PY'
from ldsc import load_gene_ldscore_index

profile = load_gene_ldscore_index(
    "/path/to/gene_ldscore_indexes/1000G_EUR_Phase3_baseline/"
    "profiles/padding-100000bp-mhc"
)
print(profile.suite_id)
print(profile.profile_id)
print(profile.chromosomes)
PY
```

A missing, corrupt, or identity-mismatched component is an error. There is no
automatic profile discovery or direct-mode fallback.

## Next step

Use the explicit profile to calculate LD scores for focal gene sets:
[Calculate LD scores for gene lists with an index](../main-functionalities/ldscore.md).
