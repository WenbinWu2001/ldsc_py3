# Exact disjoint-atom gene LD-score index specification

Last updated on: 2026-08-03

Status: approved specification; design closed

## Problem and goal

The direct gene-list route in `ldsc ldscore` resolves every requested list, projects its padded gene intervals onto the full retained reference-SNP grid, and recomputes annotation LD scores from PLINK or pairwise-R2 data. This is exact but unnecessarily expensive when many users query the same gene catalog, baseline suite, reference panel, LD window, and regression-row policy.

Design 2 adds an expensive offline command, `ldsc build-gene-ldscore-index`, and an explicit indexed mode on the existing `ldsc ldscore` command. The index represents overlapping padded genes as disjoint genomic atoms and stores the exact linear operator and sufficient statistics needed to assemble focal gene-list annotations. Online queries require only the selected index profile and gene-list files; they do not require the source PLINK, R2, baseline annotation, or genetic-map files.

Success means that indexed and matching direct runs produce the same binary annotations, retained rows and ordering, annotation counts, overlap entries, query statuses, canonical LD-score artifact structure, and downstream `partitioned-h2` interpretation. Floating-point LD scores and regression results must agree within a tolerance measured by the chromosome-22 validation gate. Approximation, whole-gene LD windows, and summing overlapping per-gene LD scores are outside the contract.

## Relationship to the current direct workflow

The direct gene-list implementation remains the scientific oracle and fallback. This specification preserves the resolution, naming, build-selection, partial-success, pruning, provenance, and diagnostics contracts in [Gene-List Query Annotations Design](2026-08-03-gene-list-query-annotations-design.md), subject only to these approved changes:

- rename public `--bed-padding-bp` / `bed_padding_bp` to `--padding-bp` / `padding_bp` everywhere, with no compatibility alias;
- add pre-padding gene-region filtering through `--gene-exclude-regions {none,mhc}`;
- add the fixed control option `--control-gene-list-source`, defaulting to the reserved token `all-protein-coding`;
- allow the explicit indexed adapter `--gene-ldscore-index-dir`.

The rename covers `annotate`, direct and indexed `ldscore`, the index builder, `AnnotationBuildConfig`, public Python entry points, persisted configuration snapshots, tests, and current user documentation. The operation is unchanged: it expands each BED or gene interval symmetrically, once, and clips the lower coordinate at zero.

The current corrected SNP-universe behavior in [LD-score SNP-universe contract](../current/ldscore-snp-universe-contract.md) is authoritative. In particular, named MHC and centromere exclusions affect regression rows and regression weights, not the broad LD-reference contributor universe.

## Scientific model and exactness

For one chromosome, define:

- \(m\): retained PLINK reference SNPs after genotype usability checks, optional individual selection, and optional MAF filtering;
- \(r\): persisted regression/output rows, fixed in v1 to bundled HapMap3 SNPs minus the `mhc-and-centromeres` regions and intersected with the retained reference rows;
- \(g\): catalog genes on the chromosome;
- \(t\): nonempty disjoint atoms induced by retained padded gene intervals;
- \(R\in\mathbb{R}^{m\times m}\): the exact within-window adjusted-\(r^2\) operator with its unit diagonal;
- \(P\in\{0,1\}^{r\times m}\): the persisted-row selector;
- \(H\in\{0,1\}^{m\times t}\): reference-SNP-to-atom membership, with at most one nonzero per row;
- \(B\in\{0,1\}^{g\times t}\): gene-to-atom membership;
- \(a\in\{0,1\}^{g}\): one resolved gene-list selector.

The selected-atom vector is the Boolean union

$$
z(a)=\mathbf{1}[B^\mathsf{T}a>0].
$$

Because atoms are maximal nonoverlapping intervals with a constant covering-gene set, the direct binary-union annotation is exactly

$$
q(a)=Hz(a).
$$

The offline index stores

$$
Y=PRH,
$$

so the online query score is

$$
\ell(a)=Yz(a)=PRq(a).
$$

This identity is the core correctness invariant. Duplicate identifiers, overlapping genes, nested genes, and overlaps created by padding select atoms by Boolean OR and never count a SNP twice. Atoms containing no retained reference SNP may be omitted only if this SNP-level equality and all query diagnostics remain unchanged.

Construction must reuse the current PLINK LD-window, adjusted-\(r^2\), and genotype filtering behavior. It must preserve the unit diagonal exactly once and preserve negative adjusted-\(r^2\) values. It must not clamp assembled scores or epsilon-prune sparse values.

### Counts and overlap sufficient statistics

Let \(A_B\in\mathbb{R}^{m\times b}\) be the ordered supplied baseline annotations. Store, for all retained SNPs and separately for common SNPs selected by `MAF >= common_maf_min`:

$$
d=H^\mathsf{T}\mathbf{1},\qquad D=A_B^\mathsf{T}H.
$$

Then a selected atom set gives its annotation count \(d^\mathsf{T}z\), baseline-query overlaps \(Dz\), and query self-overlap \(d^\mathsf{T}z\). The same operations assemble the all-protein-coding or custom control. Control-query overlap is computed from the atom intersection of their selectors, so no extra source annotation is needed online.

## Observable behavior

### Offline builder

`ldsc build-gene-ldscore-index` builds one profile into a suite directory. Its scientific inputs are baseline annotation sources, a PLINK prefix, SNP identity, genome build, `--padding-bp`, `--gene-exclude-regions`, and a cM LD window (default `--ld-wind-cm 1.0`). It also accepts the existing optional `--maf-min`, `--common-maf-min`, `--keep-indivs-file`, and build-specific genetic-map inputs, plus chromosome selection, SNP/atom construction tuning, threading, output, overwrite, and logging controls.

V1 supports only PLINK-backed cM construction. It rejects an R2 backend, explicit reference-panel or regression SNP files, configurable SNP-region policy, SNP-count or kb windows, whole-chromosome override, and reference-metadata export. The regression policy is always bundled HapMap3 minus `mhc-and-centromeres`. Internal atom-column batch sizing is a builder implementation/resource control and is never reused as a public online query-column batch control.

V1 has one coordinate build, hg19, and no liftover or separate output-build concept. Baseline annotation positions, PLINK BIM positions, gene-catalog projection, bundled HM3 and named-region coordinates, and any explicit genetic map must all refer to hg19. Canonical index and assembled LD-score rows inherit that coordinate system. The public `--genome-build` choice is therefore restricted to hg19. The legacy-named `--genetic-map-hg38-sources` argument is rejected when supplied; only an explicit hg19 map or informative BIM cM values are valid.

### Direct gene-list mode

Without `--gene-ldscore-index-dir`, `ldsc ldscore` follows the ordinary Design 1 path. It accepts the normal PLINK or R2 scientific inputs and can represent configurations outside the first index compatibility domain. `--gene-exclude-regions` defaults to `none`.

The direct path must implement the same gene filtering and control-annotation semantics as an otherwise matching index profile. It is the acceptance oracle for indexed computation.

### Indexed gene-list mode

Indexed mode is selected only by:

```text
ldsc ldscore \
  --gene-ldscore-index-dir <profile-dir> \
  --query-annot-gene-list-sources <sources> \
  --output-dir <ldscore-dir>
```

The profile directory, not the containing suite directory, is supplied. One selected profile represents exactly one padding value, catalog/build, gene-region policy, and immutable common suite.

Indexed mode accepts focal gene-list sources, `--control-gene-list-source`, output/overwrite/logging/threading controls, and ordinary diagnostic controls. It rejects live baseline sources, PLINK or R2 inputs, reference or regression restriction files, SNP or gene region settings, LD-window settings, `--padding-bp`, genetic maps, genome build, and SNP-identity settings. Those properties have one authority: the profile and suite metadata.

An explicitly supplied index never triggers discovery or fallback. A conflicting or unsupported option is an actionable usage error directing the user to remove `--gene-ldscore-index-dir` and use the direct route. Missing, corrupt, structurally invalid, or identity-mismatched components fail before canonical scientific outputs are published.

### Control gene list

`--control-gene-list-source` is singular and applies to direct and indexed gene-list runs:

| Value | Meaning |
| --- | --- |
| `all-protein-coding` | Default. Select every profile-included protein-coding gene from the same catalog. |
| `none` | Fit the former baseline-plus-focal-query model without a fixed gene control. |
| path | Resolve a dataset-specific control list through the same catalog, build, gene-region, and padding rules as focal lists. |

When present, the fixed column is named `gene_control`. It is appended to the supplied baseline columns in `ldscore.baseline.parquet`, is recorded in `baseline_columns`, and is not listed in `query_columns`. Consequently, the unchanged `partitioned-h2` cell-type regime fits `supplied baseline + gene_control + one focal query`. The overlap artifact stores `gene_control` as part of its baseline-rows block. A focal query or supplied baseline column named `gene_control` is a collision error.

An unusable custom control is a run-level error because the requested conditioning model cannot be formed. Partial control resolution follows the ordinary warning/audit behavior and uses the retained genes.

### Gene-region filtering

`--gene-exclude-regions {none,mhc}` is separate from SNP `--exclude-regions`. The `mhc` choice removes a gene when its unpadded transcribed interval overlaps the packaged 25-35 Mb MHC interval on chromosome 6; filtering occurs before padding. The HLA genes are located within the broader MHC region, but HLA and MHC are not treated as synonyms.

The generic direct default is `none`. Initial distributed LDSC-SEG profiles use `mhc`. Centromeric genes are not removed: the paper-supported centromere rule is a regression-row exclusion, not a gene-set definition.

Excluded genes are recorded at gene level with reason `excluded_gene_region` in the existing gene-list diagnostic audit. Partial exclusion warns and continues. A list with no retained genes is skipped, and an all-skipped batch follows the existing consolidated-error contract. When resolution and region exclusion both affect one list, both conditions remain visible in diagnostic counts/details even though the query-status record has one primary reason.

## Reference, regression, and count universes

The first indexed mode fixes the following separation:

| Quantity | SNP universe |
| --- | --- |
| Baseline and focal LD-score contributors | All retained PLINK reference SNPs, including MHC and pericentromeric SNPs |
| `M`, `M_5_50`, and annotation overlaps | The same retained reference universe; `M_5_50` additionally uses the common-MAF rule |
| Persisted baseline/query rows | Bundled HapMap3 rows minus `mhc-and-centromeres`, intersected with retained PLINK rows |
| `regression_ld_scores` / `w_ld` contributors and rows | The identical filtered regression set |
| Final regression observations | Intersection of persisted rows, weights, and munged summary statistics |

A description of v1 as “supporting only HapMap3 SNPs” is scientifically
incorrect. Only the persisted regression-row policy is fixed to bundled,
region-filtered HapMap3. The reference/contributor universe remains the broad
retained PLINK set, and non-HapMap3 or excluded-region reference SNPs can
contribute LD to a persisted HapMap3 row. Implementations and user
documentation must preserve this distinction.

A custom `--regression-snps-file` or an SNP `--exclude-regions` choice other than `mhc-and-centromeres` is outside the v1 indexed domain and remains available through the direct path. There is no v1 indexed override for the reference universe.

## Initial supported configurations

The first separately distributed suites use:

- hg19;
- `snp_identifier=rsid`, compatible with the allele-free legacy annotation shards;
- the `1000G_EUR_Phase3` PLINK suite;
- a 1 cM LD window;
- bundled HapMap3 regression rows minus `mhc-and-centromeres`;
- no explicit retained-reference MAF threshold;
- `common_maf_min=0.05` with the package's inclusive `>=` rule;
- no individual keep file;
- informative PLINK `.bim` CM values rather than an external genetic map;
- 100 kb padding and `gene_exclude_regions=mhc`.

One suite is built for each baseline annotation source:

1. `1000G_EUR_Phase3_baseline`;
2. `1000G_Phase3_baselineLD_v2.2_ldscores`.

The 100 kb profile is the historical LDSC-SEG-compatible default. A separate 0 bp gene-body profile is supported but is not a second initially distributed default; 20 kb is not a default. The LDSC-SEG paper extends gene regions by 100 kb, and its supplement reports that 100 kb was chosen from an empirical comparison of 20 kb and 100 kb settings, not from gene-length distribution theory.

### PLINK filtering and genetic maps

No `--maf-min` means no explicit frequency threshold. The existing PLINK reader still removes monomorphic, zero-variance, or unusable genotype rows; it does not apply a special singleton exclusion. Optional `--maf-min` is genotype-derived after `--keep-indivs-file` and is inclusive. The suite identity records the selected-individual content identity, selected sample count, MAF policy and source, and per-chromosome removal counts.

An explicit hg19 map is interpolated at PLINK positions and overrides `.bim` CM. V1 rejects a supplied hg38 map. Without an explicit map, informative `.bim` CM is used; an uninformative CM column is an error with guidance to supply an hg19 map. The effective hg19 map/BIM-CM source is part of suite identity.

### Strict baseline-to-PLINK alignment

All baseline inputs contributing to one chromosome must agree on row identities and ordering after canonical genomic sorting. Before genotype-derived filtering, the baseline identities must equal the PLINK identities for every selected chromosome; reordering is allowed, but missing, extra, duplicated, or conflicting identities are build errors. The initial allele-free suites compare the full `CHR/POS/SNP` triple while using `rsid` as the public matching mode. Only subsequent PLINK genotype usability checks and optional MAF filtering may reduce this common universe.

This is deliberately stricter than direct `ldscore`, which permits a baseline/PLINK intersection. Builder documentation and troubleshooting must state `baseline identities == PLINK identities` prominently. The supplied chromosome-22 files already satisfy this precondition: both approved baseline suites and the PLINK `.bim` have the same 141,123 `CHR/POS/SNP` rows.

Identity equality is a consistency check, not independent genome-build
inference. Mutually matching baseline and PLINK rows that are both hg38 but
misdeclared as hg19 could pass equality while receiving hg19 gene intervals,
HM3 rows, and region masks. The caller must therefore verify the documented
source build; a matching `CHR/POS/SNP` table alone is insufficient evidence.

## Index artifact contract

Indexes are distributed separately from the package. A physical suite holds one shared common layer and one or more selectable profiles:

```text
<suite-dir>/
  metadata.json
  common/
    metadata.json
    chr<chrom>/
      metadata.json
      baseline_rows.parquet
      baseline_statistics.npz
  profiles/
    padding-100000-mhc/
      metadata.json
      gene_catalog.parquet
      diagnostics/
        build-gene-ldscore-index.log
      chr<chrom>/
        metadata.json
        atoms.parquet
        gene_to_atom.npz
        ldscore_operator.npz
        atom_statistics.npz
```

Root metadata does not maintain a mutable profile list; loaders discover profiles beneath `profiles/`. Human-readable directory names are descriptive only. Semantic identities establish compatibility.

### Identity and publication

- `artifact_type: gene_ldscore_index` is the format guard. The index does not add a generic schema-version or software-version field.
- `suite_id` is derived from canonical content identities and settings for the ordered baseline suite, PLINK files and selected individuals, chromosome coverage, SNP/build/allele policy, genotype and MAF filters, genetic map, LD window, regression rows, and SNP-region policy.
- `profile_id` is derived from `suite_id` plus the exact gene-catalog resource/release and decompressed-content SHA-256, projection build, padding, and gene-region policy.
- Common component metadata repeats `suite_id`; profile components repeat both IDs and declare dimensions, ordering, sparse format, and dtypes.
- Large payloads do not receive redundant per-file checksums. Structural validation, semantic IDs, and staged publication protect local use; separately distributed archives may use transport checksums.

The builder stages and validates scientific payloads before atomic publication. An existing `common/` layer is reused only when its calculated `suite_id` matches. `--overwrite` replaces only the targeted profile and never changes a shared common layer beneath other profiles. Diagnostics are excluded from semantic identity and may survive a failed build without making the profile loadable.

Profiles default to chromosomes 1-22. `--chromosomes` may build a prototype subset such as chromosome 22, but that coverage is part of `suite_id` and a published profile's chromosome set is immutable. Indexed output contains only the declared chromosomes; genes that hit no retained SNP in that coverage follow the ordinary zero-annotation query behavior. V1 does not support concurrent shard append or separate finalization; extending a partial profile requires rebuilding the suite/profile under the resulting full-coverage identity.

### Common payload

`baseline_rows.parquet` is a chromosome shard with the same identity, `regression_ld_scores`, ordered baseline LD-score columns, and dtypes as the canonical public `ldscore.baseline.parquet`. It has a distinct filename because it is an internal shard, not a public unsharded LD-score result.

`baseline_statistics.npz` contains:

| Member | Shape and dtype |
| --- | --- |
| `baseline_count_all` | `(b,)`, float64 |
| `baseline_count_common` | `(b,)`, float64 |
| `baseline_overlap_all` | `(b,b)`, float64 |
| `baseline_overlap_common` | `(b,b)`, float64 |
| `total_reference_snps_all` | scalar, int64 |
| `total_reference_snps_common` | scalar, int64 |

The common layer is produced by the builder; it is not copied byte-for-byte from the legacy baseline LD-score files. Source baseline annotation files are provenance-bearing offline inputs and are not required online.

### Profile payload

`gene_catalog.parquet` is self-contained and includes all genes, including region-excluded genes. Its fields are `gene_index`, `canonical_ensembl_id`, `gene_name`, `CHR`, `start0`, `end`, `included`, `exclusion_reason`, and `chromosome_gene_row`. Coordinates are unpadded 0-based half-open intervals in the profile build. No duplicate TSV is written.

`atoms.parquet` contains deterministic chromosome-local `atom_id`, `CHR`, `start0`, and `end`, ordered by genomic start. `gene_to_atom.npz` is Boolean CSR with rows ordered by `chromosome_gene_row` and columns by `atom_id`; excluded genes have no selected atoms. A chromosome with no retained genes still has valid empty components with declared zero dimensions.

`ldscore_operator.npz` stores \(Y\) as SciPy CSR with float64 data and int32 indices. `gene_to_atom.npz` is Boolean CSR. Construction accumulates in float64 and performs no epsilon pruning.

`atom_statistics.npz` contains:

| Member | Shape and dtype |
| --- | --- |
| `atom_count_all` | `(t,)`, int64 |
| `atom_count_common` | `(t,)`, int64 |
| `baseline_atom_overlap_all` | `(b,t)`, float64 |
| `baseline_atom_overlap_common` | `(b,t)`, float64 |

## Canonical online output

An indexed run materializes the ordinary self-contained output family through the existing `LDScoreDirectoryWriter` contract:

```text
<ldscore-dir>/
  metadata.json
  ldscore.baseline.parquet
  ldscore.query.parquet
  ldscore.overlap.parquet
  diagnostics/
    ldscore.log
    query_annotation_status.tsv
    gene_list_unresolved.tsv.gz
```

The common baseline shards are reused without recomputing baseline LD scores. `gene_control`, when enabled, is assembled from atoms and appended to the baseline table/block. Focal scores are written to `ldscore.query.parquet` in resolved source order after skipped queries are removed. The overlap artifact stores the baseline-rows block, now including the control, plus focal query self-overlaps exactly as expected by `assemble_model_overlap`. Root metadata records `suite_id`, `profile_id`, control provenance, catalog/profile provenance, and the ordinary query diagnostics. A thin result that refers back to the index is outside v1.

Public Parquet files retain one row group per chromosome. All user-requested focal columns are evaluated together for a chromosome; neither direct nor indexed `ldscore` exposes a query-column batch option or an arbitrary query-count cap. Implementations may stream chromosome payloads into those row groups, while a public in-memory API may still materialize its final `LDScoreResult` to preserve the existing return contract.

## Resource constraints

Offline construction targets approximately 4-8 GB peak RSS for the default single-chromosome worker by batching SNP work and internal atom columns. The atom columns are an implementation basis, not user-requested focal query columns, so this internal batching does not violate the no-query-batching decision. The 100 kb prototype yielded about 34,546 atoms genome-wide; chromosome 1 had 3,823 atoms and 779,354 reference variants, making dense float32 annotations about 11.9 GB and a dense float64 accumulator about 23.8 GB before other working arrays. Dense SNP-by-atom construction is therefore forbidden.

`--threads` may parallelize chromosomes, defaulting to one worker. Higher concurrency is opt-in and may multiply chromosome-local memory; the build log records effective workers, batching, dimensions, and resource settings.

`partitioned-h2` deliberately loads the complete `ldscore.query.parquet` and all focal columns into memory. It does not lazily read or batch them. For the packaged HM3 map, the approved hg19 masks leave approximately 1,157,301 candidate rows before PLINK intersection. At float32, each focal column is about 4.42 MiB. Empirical resident table floors are:

| Focal columns | Baseline suite (~53 columns) | baselineLD v2.2 (~97 columns) |
| ---: | ---: | ---: |
| 100 | 1.3 GB | 1.5 GB |
| 500 | 3.1 GB | 3.3 GB |
| 1,000 | 5.5 GB | 5.7 GB |
| 2,000 | 10.1 GB | 10.3 GB |

These are DataFrame floors, not peak RSS. User documentation must identify the row-count/dtype assumptions and recommend roughly 1.5-2 times the table estimate for Parquet decoding, pandas copies, overlap data, and regression work arrays. Machines with more memory are the supported solution for very wide query sets.

## Diagnostics and operational behavior

The builder uses the package's output preflight, overwrite, staged-publication, and logging conventions. Each successful profile contains `diagnostics/build-gene-ldscore-index.log`, recording resolved inputs, common creation or reuse, suite/profile identities, baseline/PLINK alignment, sample and SNP filtering, per-chromosome dimensions, batching/resources, warnings, publication, and elapsed time. Failures preserve an actionable log without publishing a loadable scientific profile.

Online query warnings and skips reuse the direct path's `QueryAnnotationStatus` and gene-level audit. The profile catalog, rather than the installed package catalog, is the online resolution authority. Absolute source paths are not persisted as runtime dependencies.

## Validation strategy and acceptance criteria

Implementation is accepted only after the following evidence exists:

1. Atomization reproduces direct binary SNP annotations for individual, overlapping, nested, duplicate, alias-selected, padded, MHC-filtered, and chromosome-start-clipped genes.
2. Synthetic pair tests prove exact LD-window boundaries, one diagonal contribution, symmetric upper-triangle handling, preservation of negative adjusted-\(r^2\), and contributions from non-HapMap3 and excluded-region reference SNPs to retained HapMap3 rows.
3. \(Yz(a)\), counts, common counts, baseline-query overlaps, control overlaps, query self-overlaps, row keys/order, column grouping/order, statuses, and pruning match the direct oracle.
4. Compatibility tests reject mismatched catalog, build, padding, gene-region policy, baseline order/content, PLINK content or selected individuals, SNP identity, MAF/common-MAF settings, map/window policy, regression rows, chromosome coverage, shapes, dtypes, and component IDs before output publication.
5. Chromosome 6 tests prove that MHC and centromere SNPs remain LD-score contributors and count/overlap members while being absent from persisted rows and `w_ld`, and that MHC genes are filtered independently before padding.
6. A chromosome-22 prototype runs both initial baseline suites against the supplied PLINK panel with overlapping/nested genes, excluded-region-adjacent genes, and non-HapMap3 contributors. It measures direct-versus-index accumulation differences after canonical float32 output casting and establishes the acceptance tolerance before whole-genome construction.
7. Matching direct and indexed canonical directories produce downstream `partitioned-h2` coefficients, standard errors, enrichments, p-values, and query ordering within the measured tolerance.
8. Builder benchmarks record wall time, peak RSS, atom count, `nnz(Y)`, payload bytes, and batch settings for chromosome 22 and chromosome 6 before whole-genome extrapolation. Online benchmarks record cold/warm load, resolution, assembly, output time, and peak memory.
9. Existing no-index gene-list, BED, prebuilt-query, baseline-only, PLINK/R2, canonical output, overlap-aware regression, CLI, and full repository tests remain green after the intentional public rename.

Structural equality is required for annotations, identities, order, counts, overlaps, statuses, and schemas. Numerical tolerance is not guessed in this specification. A discrepancy materially larger than canonical output floating-point granularity requires an accumulation fix rather than a looser threshold.

The local chromosome-22 gate passed on 2026-08-03 for both approved baseline
suites. Canonical LD-score columns were exactly equal after output casting;
continuous baselineLD overlap accumulation differed by at most
`7.654307410120964e-08`, establishing an absolute overlap tolerance of `1e-7`.
The unchanged downstream `partitioned-h2` comparison was exactly equal. Full
measurements and inputs are recorded in
[`docs/audits/2026-08-03-exact-gene-ldscore-index-chr22.md`](../audits/2026-08-03-exact-gene-ldscore-index-chr22.md).

### Scientific and repository evidence

The region split follows the LDSC paper's Online Methods: HapMap3 and long-range/pericentromeric rules constrain regression variants, while LD from the long-range regions remains part of LD-score estimation. The stratified-LDSC supplement separately distinguishes HapMap3 regression SNPs from the broad 1000 Genomes reference/causal universe. The LDSC-SEG paper and Supplementary Note support the historical 100 kb profile and show that it was selected empirically from tested window choices.

Repository seams that the implementation must preserve include `GeneCatalog` and `GeneListResolution` in `gene_list_resolver.py`, the direct projection oracle in `AnnotationBuilder`, PLINK accumulation and regression masks in `_kernel/ldscore.py`, the canonical artifact contract in `LDScoreDirectoryWriter`, full-file loading in `load_ldscore_from_dir`, and baseline-plus-one-query overlap assembly in `assemble_model_overlap`.

## Out of scope

- R2-Parquet index construction in v1;
- custom indexed reference or regression SNP sets;
- indexed SNP exclusion choices other than `mhc-and-centromeres`;
- SNP-count or kb LD windows, or a whole-chromosome override;
- automatic index discovery or silent fallback;
- query-column batching or a public query-batch flag;
- concurrent profile append or incremental chromosome finalization;
- shipping indices inside the Python package;
- thin canonical outputs that depend on the index remaining installed;
- approximate gene-level LD windows, per-gene LD-score addition, clamping, or epsilon sparsification.

## Risks and remaining evidence gate

- Strict baseline/PLINK identity equality intentionally rejects inputs that direct `ldscore` would intersect permissively.
- A profile is scientifically rigid; changing padding, catalog, gene exclusion, sample selection, MAF, map, window, baseline, or row policy requires another compatible profile or suite.
- Wide focal sets can make canonical output and downstream `partitioned-h2` memory dominate even when indexed score assembly is fast.
- The initial hg19 catalog coordinates retain the documented upstream provenance limitation from Design 1; the implementation must not claim an unrecorded liftover method.
- Omitting per-payload checksums trades fine-grained corruption detection for simpler structural validation and archive-level transport checksums.

There are no unresolved product, scientific, or architectural decisions. The
local chromosome-22 gate has passed; production whole-genome construction and
full chromosome-6 resource measurement remain separate distribution/HPC work.
