# Genome-Build and Coordinate-Basis Inference

Automatic genome-build inference is used for coordinate-family workflows when
the user passes `--genome-build auto` or workflow-specific source build
`auto`, including `munge-sumstats --source-genome-build auto` and
`build-ref-panel --source-genome-build auto`. It is not used in `rsid`-family modes.

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
`src/ldsc/data/hm3_chr_pos_reference.tsv.gz` (~104 KB compressed). It contains
11,000 autosomal HapMap3 SNP positions with columns:

| Column | Meaning |
|---|---|
| `CHR` | chromosome |
| `hg19_POS` | 1-based hg19 position |
| `hg38_POS` | 1-based hg38 position |

The reference has 500 SNPs per autosome. SNPs were selected after filtering to
common, non-strand-ambiguous `A1/A2` sites that are unique in both `(CHR,
hg19_POS)` and `(CHR, hg38_POS)`, have different hg19 and hg38 positions, and
are evenly spaced by position after filtering. Rebuild it from the curated map
with `python -m ldsc.hm3_reference --curated-map src/ldsc/data/hm3_curated_map.tsv.gz --output src/ldsc/data/hm3_chr_pos_reference.tsv.gz`.

The reference is loaded by `load_packaged_reference_table()` and cached once
per Python process with `lru_cache`.

The separate full curated map
`src/ldsc/data/hm3_curated_map.tsv.gz` is the package-backed HM3 restriction and
quick-liftover map. Public code loads it through `ldsc.load_hm3_curated_map()`,
which preserves packaged columns while normalizing `CHR`, `hg19_POS`,
`hg38_POS`, `SNP`, `A1`, and `A2`. The same file is parsed by SNP-restriction
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

If either threshold is not met, inference fails and the user should pass
`--genome-build hg19` or `--genome-build hg38` explicitly.

---

## Input-Specific Sampling

Different inputs are sampled differently so inference reads only the data it
needs.

| Input type | Strategy |
|---|---|
| Packaged inference reference | Load the compact 11,000-SNP HM3 coordinate reference once, then reuse from cache |
| Raw sumstats text / `.gz` | Stream lightweight `CHR`/`POS` chunks into the HM3 evidence accumulator and stop once evidence is sufficient |
| Annotation chromosome-suite inputs | Read a small head sample from the first resolvable `@` chromosome file |
| `ldscore --r2-dir` directory | Locate candidate R2 parquet files, then infer from `ldsc:sorted_by_build` schema metadata |
| PLINK `.bim` source panel | Stream `.bim` `CHR/BP` chunks into the HM3 evidence accumulator before SNP restriction when `source_genome_build` is `auto` |
| build-ref-panel SNP restriction generic `POS` | Infer the restriction file's local build and require it to match the source PLINK build |
| Canonical parquet R2 reference panel | Prefer schema metadata; otherwise inspect the first row group |

Adaptive evidence collection is intentionally limited to large sequential
inputs: raw summary statistics and PLINK `.bim` source panels. Small
chromosome-suite artifacts keep their fixed head-sample behavior, and
metadata-backed parquet artifacts remain metadata-first.

## Raw Sumstats

Raw sumstats can be large, compressed, and sequential-access. In
`chr_pos`-family modes (`chr_pos` and `chr_pos_allele_aware`),
`ldsc munge-sumstats` defaults `--source-genome-build` to `auto` and requires
`--output-genome-build hg19` or `--output-genome-build hg38` in `chr_pos`-family
modes. When source `auto` is requested, the workflow resolves it before chunk parsing.
After the raw header has been mapped to canonical columns, the workflow streams
lightweight `CHR`/`POS` chunks into the shared build-evidence accumulator and
stops once the inference thresholds are met. The resolved source build and
coordinate basis are then reused while each chunk is parsed.

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
the user to pass `--source-genome-build hg19` or `--source-genome-build hg38` explicitly. The
same resolved source-build metadata is used by munge-time SNP filtering in
`chr_pos`-family modes. A keep-list, including the packaged map selected by
`--use-hm3-snps`, with both `hg19_POS` and `hg38_POS` selects the position
column matching the effective raw sumstats build. Filtering happens inside each
chunk after canonical coordinate normalization and before the retained chunks are
concatenated.

Optional munger liftover runs after this source-build resolution and after
SNP filtering. It is valid for chr_pos-family modes because rsID-family modes
do not use positions for row identity. Chain-file liftover and HM3 quick liftover both
update `CHR`/`POS` only; `SNP` remains a label. HM3 quick liftover requires
`--use-hm3-snps`, so HM3 filtering and HM3 coordinate conversion are explicit
separate steps. The metadata sidecar records the final output build in its
`genome_build` identity field; the source/target/method/drop counts,
duplicate-coordinate drops, and coordinate inference details are written to the
run log.

## Build-Ref-Panel Source Inputs

PLINK `.bim` files can also be large across chromosome suites. When
`build-ref-panel --source-genome-build auto` is used, the workflow reads only
the `.bim` chromosome and base-pair columns in chunks, feeds them to the shared
HM3 evidence accumulator, and stops once the same inference thresholds are met.
The final build decision still goes through `resolve_genome_build()`, so
confidence thresholds and error messages stay aligned with other workflows.

SNP restriction files are not converted to adaptive streaming. They are small
control artifacts in typical use and are read for filtering after source-build
resolution, so their generic `POS` build check continues to use the materialized
restriction frame.

## Annotation Inputs

For annotation workflows that use chromosome-suite path tokens, such as
`baseline.@.annot.gz`, inference samples the first existing chromosome-specific
file. Chromosome `1` is tried first, followed by any configured chromosome
order. The sampler reads a small head sample and infers the `CHR` and `POS`
columns using the shared column-alias registry.

If the sample does not contain enough informative HapMap3 matches, pass the
build explicitly.

## LD-Score Runtime Auto Build

The `ldsc ldscore` CLI requires `--genome-build` in coordinate-family modes.
Passing `--genome-build auto` defers the choice until runtime, where the workflow
collects build evidence from the supplied LD-score inputs and requires all
available evidence to agree.

There is no silent precedence between baseline annotations and the R2 reference
panel. The workflow checks both sources when present:

1. If `--baseline-annot-sources` is supplied, sample a chromosome-suite
   annotation file and infer its coordinate build from `CHR` and `POS`.
2. If `--r2-dir` is supplied, infer the reference-panel build from the
   parquet schema metadata key `ldsc:sorted_by_build`.
   - Directory names such as `hg19` or `hg38` are not build evidence.
   - The directory layout is used only to find candidate files: first
     `chr*_r2.parquet` directly under `--r2-dir`, otherwise files under
     conventional `hg19/` and `hg38/` children.
   - If no candidate parquet has `ldsc:sorted_by_build`, the R2 directory
     contributes no build evidence.
   - If candidate parquets report conflicting builds, the workflow raises.
3. If both annotation and R2-directory evidence are available, they must infer
   the same build. A disagreement raises an error instead of selecting one over
   the other.

Examples:

```text
--genome-build auto --r2-dir /panels/my_panel/hg38
  -> the build recorded in /panels/my_panel/hg38/chr*_r2.parquet metadata

--genome-build auto --r2-dir /panels/my_panel
  where /panels/my_panel contains only hg19/
  -> the build recorded in /panels/my_panel/hg19/chr*_r2.parquet metadata

--genome-build auto --baseline-annot-sources baseline.@.annot.gz --r2-dir /panels/my_panel/hg38
  -> succeeds only if the baseline annotation sample agrees with the parquet metadata

--genome-build auto --r2-dir /panels/my_panel
  where /panels/my_panel contains both hg19/ and hg38/
  -> succeeds only if all candidate parquets report the same metadata build;
     otherwise errors with conflicting R2 parquet genome-build metadata
```

For durable scripts, prefer a concrete build plus a build-specific R2
directory, for example `--genome-build hg38 --r2-dir /panels/my_panel/hg38`.

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
1-based or 0-based, and the support for all four hypotheses. For sumstats, the
root metadata records only the final compatibility build through the minimal
`genome_build` identity field; inferred-build evidence, coordinate basis, and
count-level liftover provenance stay in the log. Row-level liftover drops are
audited in `diagnostics/dropped_snps/dropped.tsv.gz`, while example SNPs remain
`DEBUG`-only.

All normalized `chr_pos` coordinates used downstream are 1-based.
