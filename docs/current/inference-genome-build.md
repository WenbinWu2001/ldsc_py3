# Genome-Build and Coordinate-Basis Inference

Automatic genome-build inference is used only for `chr_pos` workflows when the
user passes `--genome-build auto`. It is not used in `rsid` mode.

The package infers the build by comparing a subset of input (CHR, POS) pairs against the packaged HapMap3 coordinate map, then selecting the hypothesis that best explains those positions:

- `hg19`, 1-based positions
- `hg19`, 0-based positions
- `hg38`, 1-based positions
- `hg38`, 0-based positions

If a 0-based input is detected, positions are converted to the package's
canonical 1-based coordinates before downstream matching.

---

## Reference Map

The packaged reference is
`src/ldsc/data/hm3_chr_pos_reference.tsv.gz` (~104 KB compressed). It contains
11,000 autosomal HapMap3 SNP positions with columns:

| Column | Meaning |
|---|---|
| `CHR` | chromosome |
| `hg19_POS` | 1-based hg19 position |
| `hg38_POS` | 1-based hg38 position |

The reference has 500 SNPs per autosome. SNPs were selected after filtering to
common, non-strand-ambiguous, conflict-free sites that are unique in both
`(CHR, hg19_POS)` and `(CHR, hg38_POS)`, have different hg19 and hg38
positions, and are evenly spaced by position after filtering.

The reference is loaded by `load_packaged_reference_table()` and cached once
per Python process with `lru_cache`.

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

If either threshold is not met, inference fails and the user should pass
`--genome-build hg19` or `--genome-build hg38` explicitly.

---

## Input-Specific Sampling

Different inputs are sampled differently so inference reads only the data it
needs.

| Input type | Strategy |
|---|---|
| Packaged HM3 reference | Load full 11,000-SNP map once, then reuse from cache |
| Raw sumstats text / `.gz` | Kernel-side normalization of munged `CHR` + `POS`, with build inference when requested |
| Annotation chromosome-suite inputs | Read a small head sample from the first resolvable `@` chromosome file |
| Canonical parquet R2 reference panel | Prefer schema metadata; otherwise inspect the first row group |

## Raw Sumstats

Raw sumstats can be large, compressed, and sequential-access. `ldsc
munge-sumstats` keeps `--genome-build auto` unresolved until the munging kernel
has parsed the raw table, applied standard QC, and translated coordinate aliases
into canonical `CHR` and `POS` columns. The kernel then calls
`resolve_chr_pos_table()` so build inference and coordinate-basis normalization
operate on the same coordinate table that will be written to the curated
artifact.

The build resolver uses `GenomeBuildEvidenceAccumulator` from
`genome_build_inference.py`:

```text
normalize CHR spellings
coerce POS to integer base-pair positions
compare unique valid (CHR, POS) pairs to the packaged HM3 hypotheses
require:
  - at least 200 informative HM3 matches, AND
  - best hypothesis explains >= 99% of informative matches
```

If those thresholds are not met, inference raises an actionable error that asks
the user to pass `--genome-build hg19` or `--genome-build hg38` explicitly. The
same normalized coordinate metadata is also used by munge-time
`--sumstats-snps-file` filtering in `chr_pos` mode, so a keep-list with both
`hg19_POS` and `hg38_POS` selects the position column matching the effective raw
sumstats build.

## Annotation Inputs

For annotation workflows that use chromosome-suite path tokens, such as
`baseline.@.annot.gz`, inference samples the first existing chromosome-specific
file. Chromosome `1` is tried first, followed by any configured chromosome
order. The sampler reads a small head sample and infers the `CHR` and `POS`
columns using the shared column-alias registry.

If the sample does not contain enough informative HapMap3 matches, pass the
build explicitly.

## Parquet R2 Reference Panels

Canonical parquet R2 files written by `ldsc build-ref-panel` store the build in
Arrow schema metadata under `ldsc:sorted_by_build`. Reading that metadata costs
only a parquet footer read and no row data.

For externally generated canonical R2 files that lack this metadata,
`SortedR2BlockReader._init_canonical_path()` reads row group 0 with only
`CHR` and `POS_1`, then runs the same coordinate inference. Dense reference
panels usually have enough HapMap3 overlap in the first row group.

If the inferred or declared R2 build disagrees with the runtime
`--genome-build`, the reader raises an error rather than mixing coordinate
systems.

---

## Output Behavior

When inference succeeds, logs include the selected build, whether the input was
1-based or 0-based, and the support for all four hypotheses. For sumstats,
the metadata sidecar also records the inferred build and coordinate basis.

All normalized `chr_pos` coordinates used downstream are 1-based.
