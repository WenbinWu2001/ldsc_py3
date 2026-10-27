# Genome-Build and Coordinate-Basis Inference

Last updated on: 2026-10-27

## Standalone artifact inspection

Use `ldsc infer-genome-build PATH` to inspect one exact file or recognized LDSC artifact directory without running a workflow.
It reports the artifact's own build independently of its SNP matching mode, including recorded rsID workflows.
Concrete metadata and build-specific columns are reported as declared; otherwise the helper samples at most 5,000 usable rows per file/chromosome using the HM3 decision rule below.
Insufficient evidence, rsID-only inputs, and ordinary BED intervals without metadata yield `unknown` with exit 0.
Interval endpoints are not SNP positions.
All relevant members participate: unresolved evidence prevents a single inferred result, incompatible builds are errors, and deliberately available hg19/hg38 alternatives receive separate labels.
The helper reports coordinate basis without converting positions and creates no logs or artifacts.
See the [complete inspection contract](../reference/snp-identity-and-genome-build.md#standalone-build-inspection) for supported inputs, the exact rsID notice, and exit codes.

## Workflow inference

Coordinate-dependent workflows resolve builds through their existing `--genome-build auto` or source-build options.
Their current metadata precedence, bounded sampling, and all-input consistency rules are defined in the [shared build policy](../reference/snp-identity-and-genome-build.md#genome-build-policy).
The workflow-specific descriptions below concern workflow processing, not the read-only inspection helper.

The package infers the build by comparing a subset of input (CHR, POS) pairs against the packaged HapMap3 coordinate map, then selecting the hypothesis that best explains those positions:

- `hg19`, 1-based positions
- `hg19`, 0-based positions
- `hg38`, 1-based positions
- `hg38`, 0-based positions

If a 0-based input is detected, positions are converted to the package's
canonical 1-based coordinates before downstream matching.

---

## Reference Maps

Genome-build inference uses the compact packaged reference
`src/ldsc/data/hm3_chr_pos_reference.tsv.gz` (~200 KB compressed). It contains
22,000 autosomal HapMap3 SNP positions with columns:

| Column | Meaning |
|---|---|
| `CHR` | chromosome |
| `hg19_POS` | 1-based hg19 position |
| `hg38_POS` | 1-based hg38 position |

The reference has up to 1,000 SNPs per autosome.
SNPs were selected after filtering to non-strand-ambiguous `A1/A2` sites that are unique in both `(CHR, hg19_POS)` and `(CHR, hg38_POS)` and have different hg19 and hg38 positions, and are evenly spaced by hg38 position after filtering.
No allele-frequency filter is applied; the curated map carries no frequencies.
Rebuild it from the curated map with `python tools/hm3/build_hm3_chr_pos_reference.py --curated-map src/ldsc/data/hm3_curated_map.tsv.gz --output src/ldsc/data/hm3_chr_pos_reference.tsv.gz`.

The reference is loaded by `load_packaged_reference_table()` and cached once
per Python process with `lru_cache`.

The separate full curated map
`src/ldsc/data/hm3_curated_map.tsv.gz` is the package-backed HM3 restriction and
quick-liftover map. Public code loads it through `ldsc.load_hm3_curated_map()`,
which returns its six columns `CHR`, `hg19_POS`, `hg38_POS`, `SNP`, `A1`, and
`A2` in normalized form. The same file is parsed by SNP-restriction
readers, so allele-aware modes use the packaged `A1/A2` values for HM3
keep-list identity matching.

---

## Decision Rule

The inference engine is `infer_chr_pos_build()` in
`src/ldsc/genome_build_inference.py`.

1. Deduplicate valid input `(CHR, POS)` pairs.
2. Compare each pair against the four hypothesis sets above.
3. Count only uniquely informative matches: a pair is informative only if it
   matches exactly one hypothesis.
4. Require at least 200 informative HapMap3 matches.
5. Require the best hypothesis to explain at least 99% of informative matches.

If either threshold is not met, standalone inspection reports `unknown` successfully.
Coordinate-dependent workflows instead require a concrete build through their supported build option or matching upstream artifacts; see the shared policy for command scope.

---

## Input-Specific Sampling and Workflow Resolution

`GenomeBuildCheck` in [_genome_build_policy.py](../../src/ldsc/_genome_build_policy.py) collects applicable declarations and bounded coordinate evidence across relevant files and chromosomes before resolving a workflow build.
Concrete workflow build options must agree with existing declarations and supply the build for unlabeled inputs without fresh inference.
Without a concrete assertion, incompatible evidence or unresolved required coordinates stop the workflow.
The read-only inspector reuses table readers, aliases, and the inference engine but keeps its successful `unknown` reporting separate from workflow resolution.

| Input type | Evidence |
| --- | --- |
| Raw sumstats text/gzip or PLINK BIM | Stream coordinate columns and retain at most 5,000 usable rows per file/chromosome for inference. |
| Annotation chromosome suites | Check all selected shards, using applicable metadata before coordinate samples. |
| Parquet artifacts | Prefer concrete build declarations in schema metadata; otherwise use recognized SNP coordinate columns. |
| Build-specific coordinate columns | Treat named hg19/hg38 positions as declared alternatives, selecting the workflow's required build for computation. |
| Ordinary BED or rsID-only filtering lists | No SNP coordinate inference; BED inherits the workflow build and rsID filtering does not require coordinate validation. |

Path expansion belongs to workflows; standalone inspection requires one exact path.
A panel directory name alone is not evidence of its coordinate build.
Use a build-specific R² directory when a workflow requires one build; the inspection helper can report both explicitly available builds at the panel root.

Munging resolves source coordinates before requested conversion and records the final output build in its Parquet identity metadata.
Source-build evidence must not be substituted for that output build during artifact inspection.
For current source/output options and conversion behavior, see [munging](../reference/munge-sumstats.md) and the [shared policy](../reference/snp-identity-and-genome-build.md#genome-build-policy).

## Parquet R2 Reference Panels

Canonical parquet R2 files written by `ldsc build-r2-panel` store the build in
Arrow schema metadata under `ldsc:sorted_by_build`. Reading that metadata costs
only a parquet footer read and no row data.

Index-format R2 files **require** `ldsc:sorted_by_build`: a file lacking it is a
hard error in `SortedR2BlockReader._init_index_path()` (there is no row-data
inference fallback). Regenerate such files with `ldsc build-r2-panel`.

If the inferred or declared R2 build disagrees with the runtime
`--genome-build`, the reader raises an error rather than mixing coordinate
systems.

---

## Workflow Output Behavior

When inference succeeds, logs include the selected build, whether the input was
1-based or 0-based, and the support for all four hypotheses. For sumstats, the
root metadata records only the final compatibility build through the minimal
`genome_build` identity field; inferred-build evidence, coordinate basis, and
count-level liftover provenance stay in the log. Row-level liftover drops are
audited in `diagnostics/dropped_snps/dropped.tsv.gz`, while example SNPs remain
`DEBUG`-only.

All normalized `chr_pos` coordinates used downstream are 1-based.
