# Genome-Build and Coordinate-Basis Inference

Automatic genome-build inference is used for `chr_pos` workflows when the user
passes `--genome-build auto`, and by `build-ref-panel` when
`--source-genome-build` is omitted. It is not used in `rsid` mode.

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
| `ldscore --r2-dir` directory | Locate candidate R2 parquet files, then infer from `ldsc:sorted_by_build` schema metadata |
| PLINK `.bim` source panel | `build-ref-panel` reads `.bim` `CHR/BP` rows before SNP restriction when `source_genome_build` is omitted |
| build-ref-panel SNP restriction generic `POS` | Infer the restriction file's local build and require it to match the source PLINK build |
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

Optional munger liftover runs after this source-build resolution and after
`--sumstats-snps-file` filtering. It is valid only in `chr_pos` mode because
`rsid` mode does not use positions for row identity. Chain-file liftover and
HM3 quick liftover both update `CHR`/`POS` only; `SNP` remains a label. The
metadata sidecar records the final output build through `config_snapshot`; the
source/target/method/drop counts, duplicate-coordinate drops, and coordinate
inference details are written to the run log.

## Annotation Inputs

For annotation workflows that use chromosome-suite path tokens, such as
`baseline.@.annot.gz`, inference samples the first existing chromosome-specific
file. Chromosome `1` is tried first, followed by any configured chromosome
order. The sampler reads a small head sample and infers the `CHR` and `POS`
columns using the shared column-alias registry.

If the sample does not contain enough informative HapMap3 matches, pass the
build explicitly.

## LD-Score Runtime Auto Build

The `ldsc ldscore` CLI requires `--genome-build` in `chr_pos` mode. Passing
`--genome-build auto` defers the choice until runtime, where the workflow
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
1-based or 0-based, and the support for all four hypotheses. For sumstats,
the metadata sidecar records only the final compatibility build through
`config_snapshot`; inferred-build evidence and coordinate basis stay in the log.

All normalized `chr_pos` coordinates used downstream are 1-based.
