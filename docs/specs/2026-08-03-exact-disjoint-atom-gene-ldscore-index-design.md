# Exact disjoint-atom gene LD-score index design

Last updated on: 2026-08-03

Status: Design interrogation in progress. Round 1--3 decisions are recorded below;
the open decisions remain intentionally unresolved.

## Objective

Design an offline-built index that answers gene-list LD-score queries exactly under
the selected gene catalog, padding, baseline annotation suite, and PLINK reference
panel. A query should not require access to the source PLINK suite, precomputed
pairwise-R2 Parquet files, or baseline reference files.

The index uses exact disjoint genomic atoms to recover the same query annotation
and LD scores as the corresponding direct computation. Approximation is not part
of the contract.

## Confirmed decisions

### Artifact boundary and online contract

- One index directory is self-contained and is bound to one immutable scientific
  configuration. At minimum, that configuration includes the padding value, gene
  catalog and genome build, baseline annotation suite, PLINK reference panel and
  retained LD-reference universe, genetic map and LD-window rule, SNP identity
  policy, filters and exclusions, and the regression-row subset.
- The online command receives an explicit index directory and a gene list. It
  does not discover a compatible index, fall back to direct computation, or
  accept flags that silently alter the index's scientific configuration.
- A conflicting or incompatible request fails with a specific error.
- Indexes are distributed separately from the Python package.

### Reference data and SNP universes

- PLINK is the primary source for offline index construction. Pairwise-R2 Parquet
  files are not required.
- The initial reference configuration uses a 1 cM LD window and the
  `1000G_EUR_Phase3` PLINK suite. The local chromosome-22 fixture is:
  `resources/example_1kg_phase3_PLINK`.
- HapMap3 SNPs select the output/regression rows only. The LD-reference universe
  is not restricted to HapMap3: every retained non-HapMap3 reference SNP can
  still contribute to the LD score of a HapMap3 output SNP.
- In partitioned heritability regression, the effective regression rows are the
  intersection of the index's output rows and the SNPs present in the munged
  summary statistics, under the effective SNP identity key.

### Initial baseline configurations

The first supported configurations use the 1 cM PLINK setup above with each of:

1. `resources/baseline_ld_suites/1000G_EUR_Phase3_baseline`
2. `resources/baseline_ld_suites/1000G_Phase3_baselineLD_v2.2_ldscores`

The baseline annotation shards define the annotation inputs. Existing derived
legacy LD-score files do not replace construction of the index's own canonical
baseline artifacts.

Both initial distributed configurations use 100-kb padding. This is the named
LDSC-SEG-compatible padding profile. Padding remains a configurable index-build
parameter so that a separate 0-kb gene-body index can be built; 20 kb is not a
second default.

### All-protein-coding-gene control

- Every padding-specific index contains a fixed all-protein-coding-gene control
  annotation constructed from the same catalog and padding as the query atoms.
- Indexed gene-list partitioned-heritability regressions fit the baseline
  annotations, the all-protein-coding-gene control, and one focal gene-list query
  at a time. This intentionally expands the earlier direct gene-list workflow's
  baseline-plus-query model to match the LDSC-SEG conditioning structure.

### Artifact identity metadata

- The index uses `artifact_type: gene_ldscore_index` as its format guard. It does
  not store a generic schema version or software version.
- One semantic `index_id` is derived from a canonical manifest containing the
  scientific configuration and content identities for the gene catalog, baseline
  annotation inputs, and PLINK `.bed`, `.bim`, and `.fam` inputs.
- Component headers repeat the `index_id`, dimensions, and dtypes so components
  from different builds cannot be mixed silently.
- Large derived matrices do not each receive an additional payload checksum.
  Staged publication protects against partial workflow output, and distribution
  archives may use a separate transport checksum.

### Physical suite and profile layout

- One physical index suite is bound to the shared baseline annotation suite,
  PLINK panel and individuals, LD window, SNP identity and filtering policy, and
  region-exclusion policy. Its payload is sharded by autosomal chromosome.
- The suite contains a shared `common/` layer and one or more logical padding
  profiles. Each padding profile is bound to exactly one catalog and padding
  value and is itself sharded by chromosome.
- The `common/` layer contains canonical baseline-derived artifacts produced by
  the index builder. It is not a byte-for-byte copy of the input baseline suite.
  Padding-dependent atom and query artifacts are kept in the profile layer.
- The complete suite is relocatable and self-contained. Baseline annotation
  sources are offline builder inputs; an online query does not require the user
  to supply `--baseline-annot-sources`.
- Sharing `common/` refines the earlier physical-directory statement without
  changing scientific identity: one logical profile still represents exactly
  one padding value and one immutable common configuration.

### Control-gene universe

- The catalog-wide all-protein-coding annotation is the stable default control.
  It is an LDSC-SEG-style standardized control, not necessarily a literal copy
  of a paper analysis whose all-genes universe depended on the expression
  dataset.
- Exact replication may replace the default with one explicitly supplied,
  dataset-specific control-gene list. The control is evaluated through the same
  exact index as focal queries and does not require rebuilding the index.

### Lifecycle and resource constraints

- Index creation follows the repository's existing output-family preflight,
  overwrite, removal, and staged-publication conventions.
- Construction and query execution must be explicitly bounded for a practical
  single-run memory budget of 4--8 GB. The design must stream or batch over SNPs,
  atoms, and queries rather than materializing dense SNP-by-atom or
  SNP-by-query matrices of unbounded size.
- There is no arbitrary semantic limit on the number of query gene lists. Memory
  bounds are achieved through batching.

## Historical padding evidence

LDSC-SEG extended each transcribed gene region by 100 kb on both sides. Its
supplement compared 20 kb and 100 kb windows, crossed with three gene-set size
choices, in two trait/tissue analyses and selected 100 kb because it produced the
strongest significance in that comparison. This supports 100 kb as a historical
LDSC-SEG replication setting; it does not establish 100 kb as a universal margin
derived from gene length.

The repository's current direct gene-list implementation uses the requested
padding value and therefore can serve as the exact oracle for any chosen value.

## Open decisions

The following items remain on the grilling frontier:

1. LDSC-SEG region-exclusion semantics: HLA-gene filtering, regression-row
   exclusion, LD-reference-universe exclusion, and whether the repository's
   centromere preset belongs in this profile.
2. Exact filenames and sufficient-statistic layout for common baseline overlap,
   profile-specific baseline-to-atom overlap, and all/common SNP universes.
3. The optional dataset-specific control-list CLI and how its fixed-covariate
   role is encoded for `partitioned-h2`.
4. Exact atom-index representation, numeric dtypes, batching, and the acceptance
   tolerance against the current direct-computation oracle.
5. Detailed treatment of panel filters, allele identity, duplicate variants,
   genetic maps, and baseline-annotation row alignment during construction.
6. Command names, output filenames, and the exact boundary between reusable
   canonical baseline artifacts and query-specific outputs.
