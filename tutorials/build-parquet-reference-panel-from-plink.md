# Build a Standard R2 Parquet Reference Panel from PLINK

Goal: start from a PLINK reference panel (`.bed/.bim/.fam`) and build a standard per-chromosome parquet representation of pairwise R2, plus sidecar metadata that can be used in downstream LDSC workflows.

This tutorial is written for first-time users. The examples use the chr22 1000 Genomes 30x example files bundled in this repository, but the same workflow applies to your own PLINK reference panel.

All example paths below are relative to the workspace root that contains both `resources/` and `ldsc_py3_restructured/`.

## What This Builder Produces

The `build-ref-panel` workflow converts PLINK genotypes into build-specific R2 reference-panel files for each chromosome:

- one R2 parquet per emitted build (`hg19/chr*_r2.parquet` and/or `hg38/chr*_r2.parquet`): one row per unordered SNP pair within the chosen LD window
- one runtime metadata sidecar per emitted build (`hg19/chr*_meta.tsv.gz` and/or `hg38/chr*_meta.tsv.gz`): one row per retained SNP, used by LDSC-style downstream tools
- one dropped-SNP audit file under `diagnostics/dropped_snps/` for each processed chromosome, header-only when no liftover-stage rows are dropped
- `diagnostics/build-ref-panel.log` when run through the CLI or convenience wrapper; a
  concrete single-chromosome PLINK prefix writes
  `diagnostics/build-ref-panel.chr<chrom>.log` so parallel per-chromosome jobs do not
  share one log file

The R2 output is a long pairwise table, not a dense square matrix on disk. That is usually the practical format for large reference panels.

Each package-built R2 parquet records the LD reference sample size and R2 bias
state plus minimal identity provenance in Arrow schema metadata. Downstream
`ldsc ldscore` runs can therefore read `ldsc:n_samples`, `ldsc:r2_bias`,
`ldsc:snp_identifier`, and `ldsc:genome_build` directly and do not need
`--r2-bias-mode` or `--r2-sample-size` for panels built by this workflow.
Allele-aware modes require package-built canonical R2 parquet with endpoint
allele columns `A1_1/A2_1/A1_2/A2_2`.

By default, the builder keeps all SNPs in the PLINK panel after:

- optional user-requested filters
- automatic liftover sanity filtering, when a usable source-to-target chain is provided

When a matching chain is provided in a chr_pos-family mode, SNPs are dropped if they fail hg19/hg38 liftover or liftover onto a different chromosome in the target build, then both build-specific R2 and metadata outputs are written. Matching chain liftover is not used in rsID-family modes because row identity is the SNP label there. When no usable matching chain is provided, the builder logs that it is skipping liftover and writes source-build-only outputs.

## Input Files

### Standard PLINK files

The builder expects a standard PLINK trio:

- `.bed`: binary genotype matrix
- `.bim`: variant metadata
- `.fam`: sample metadata

You pass the file prefix, not the individual filenames. For example, if you have:

- `genomes_30x_chr22.bed`
- `genomes_30x_chr22.bim`
- `genomes_30x_chr22.fam`

then the PLINK prefix is:

```text
resources/example_1kg_30x/genomes_30x_chr22
```

The `.bim` file is especially important here because it supplies:

- chromosome
- SNP identifier
- genetic-position placeholder from PLINK
- base-pair position
- allele labels

The builder copies the PLINK `.bim` `SNP` field directly into the output SNP
labels. If the `.bim` uses dbSNP IDs, downstream metadata and R2 endpoint
columns will show dbSNP IDs. If the `.bim` uses labels like
`22:10684250:C:G`, those exact strings are preserved; `build-ref-panel` does
not convert them to rsIDs. The `--snp-identifier` option only controls how a
`--ref-panel-snps-file` restriction is interpreted (`rsid` versus `chr_pos`);
it does not rewrite the emitted `SNP`, `SNP_1`, or `SNP_2` values.

### Chromosome selection

There is currently no separate `--chrom` flag.

Chromosome selection comes from the PLINK prefix token you provide through `--plink-prefix`:

- use `--plink-prefix <single-prefix>` when the prefix already points to one chromosome or one specific PLINK dataset
- use `--plink-prefix <suite-token>` when you have one PLINK prefix per chromosome, typically with an explicit `@` token
- use `--plink-prefix <glob-like-token>` only when that token resolves cleanly to one or more complete PLINK prefixes; unlike scalar file inputs elsewhere, a PLINK prefix token is resolved at the prefix level rather than at the individual `.bed/.bim/.fam` file level

Examples:

```text
--plink-prefix resources/example_1kg_30x/genomes_30x_chr22
--plink-prefix data/reference/genomes_30x_chr@
```

If the input resolves to multiple chromosomes, the builder writes one output set per chromosome.

In other words, the old split between `--plink-prefix` and a dedicated per-chromosome flag is gone. The unified `--plink-prefix` argument now handles either one concrete prefix or an explicit chromosome suite such as `panel_chr@`.

### Genetic map files

Genetic maps are optional unless the LD window is defined in centiMorgans.
They serve two separate purposes:

- every emitted build requires its matching map for `--ld-wind-cm`, because each build-specific R2 file is windowed and sorted in that build's coordinates
- any provided map populates `CM` in the matching metadata sidecar
- omitted maps are allowed for `--ld-wind-snps` and `--ld-wind-kb`; emitted metadata sidecars store `CM=NA` for builds without a map

The recommended file format is a plain text table with columns:

- `chr`
- `position`
- `Genetic_Map(cM)`

The bundled Alkes-group maps in `resources/genetic_maps/genetic_map_alkesgroup/` already use an accepted format.

## Parameters and Configuration

### Core arguments

- `--plink-prefix`
  Plain-English meaning: where the PLINK panel lives.
  Recommended usage: use `--plink-prefix` for a single chromosome, a single prefix, or a chromosome-split suite such as `panel_chr@`.

- `--source-genome-build`
  Plain-English meaning: which genome build the input PLINK coordinates already use.
  Accepted values: `hg19`, `hg37`, `GRCh37`, `hg38`, `GRCh38`.
  Recommended usage: be explicit for production runs when you know the source
  panel build. If omitted, `build-ref-panel` infers the source build from the
  PLINK `.bim` `CHR/BP` rows before applying any SNP restriction.

- `--genetic-map-hg19-sources`
  Plain-English meaning: genetic map aligned to hg19 coordinates.
  Recommended usage: provide it when hg19 output is emitted and you use
  `--ld-wind-cm`, or when you want hg19 metadata `CM` values. For SNP- or
  kb-window builds, omitting it writes hg19 `CM` as `NA` when hg19 metadata is
  emitted.

- `--genetic-map-hg38-sources`
  Plain-English meaning: genetic map aligned to hg38 coordinates.
  Recommended usage: provide it when hg38 output is emitted and you use
  `--ld-wind-cm`, or when you want hg38 metadata `CM` values. For SNP- or
  kb-window builds, omitting it writes hg38 `CM` as `NA` when hg38 metadata is
  emitted.

- `--liftover-chain-hg19-to-hg38-file` or `--liftover-chain-hg38-to-hg19-file`
  Plain-English meaning: explicit chain file used to translate positions into the other genome build.
  Recommended usage: pass the chain that matches `--source-genome-build` when you need both hg19 and hg38 outputs in `chr_pos`-family modes. If you omit it, the build completes with source-build-only outputs and logs an informational message. Do not pass a matching chain in `rsid`-family modes.

- `--output-dir`
  Plain-English meaning: output root directory.
  Recommended usage: point this to a dedicated directory for the new reference
  panel build. The run identity is the directory name; output filenames are
  fixed directly under each `{build}/` directory. Missing directories are created
  and existing directories are reused, but existing candidate parquet or
  metadata files and build-ref-panel workflow logs are refused before
  chromosome processing starts.

- `--overwrite`
  Plain-English meaning: allow replacement of the fixed parquet and metadata
  files that this build may write.
  Recommended usage: omit it for reproducible first runs. Add it only when you
  intentionally want to replace an existing panel build. `build-ref-panel`
  does not delete unrelated files or clean the output directory. It owns only
  current-contract artifacts: emitted `{build}/chr*` files, diagnostic
  metadata, build-ref-panel logs under `diagnostics/`, and
  `diagnostics/dropped_snps/` audit files. With overwrite enabled, stale owned
  siblings from earlier target-build, coordinate, or chromosome configurations
  are removed after the successful write. Removed legacy root diagnostic names
  are ignored.

### Choose exactly one LD-window option

Exactly one of the following must be set:

- `--ld-wind-cm`
  Window size in centiMorgans.
  Recommended usage: this is the most natural choice for LDSC-style reference panels. A common starting point is `1.0`.

- `--ld-wind-kb`
  Window size in kilobases.
  Recommended usage: useful when you want a purely physical window rather than a recombination-based window.

- `--ld-wind-snps`
  Window size in number of SNPs.
  Recommended usage: mostly useful for debugging or for reproducing workflows that are defined in terms of SNP counts.

### Optional filters and controls

- `--maf-min`
  Plain-English meaning: minimum minor allele frequency among retained SNPs.
  Optional: yes.
  Recommended usage: leave unset to keep all panel SNPs by default. Use values like `0.01` or `0.05` only if your downstream analysis truly wants a frequency-restricted panel.

- `--ref-panel-snps-file`
  Plain-English meaning: restrict the retained reference-panel SNP rows to a
  user-supplied list.
  Optional: yes.
  Recommended usage: use this when you want to build a panel only on HM3 SNPs,
  a curated common-SNP list, or another pre-defined reference universe.

  When this path is supplied, the restriction key comes from the invocation
  `GlobalConfig.snp_identifier`. Restriction files may omit alleles and then
  match by base key. Allele-bearing restrictions in allele-aware modes match by
  the effective allele-aware key. The packaged HM3 restriction is allele-free and
  matches by base key. Restriction files are identity-only filters: duplicate
  restriction keys collapse to one retained key, and non-identity columns such
  as `CM` or `MAF` are ignored. In `chr_pos`-family modes, the restriction file
  must be aligned to the PLINK source build. Target-build restriction files are
  not lifted over by the builder.

- `--use-hm3-snps`
  Plain-English meaning: restrict the retained reference-panel SNP rows to the
  packaged curated HM3 map.
  Optional: yes.
  Recommended usage: prefer this over a custom HM3 file path when the desired
  universe is HapMap3. It is mutually exclusive with `--ref-panel-snps-file`.

- `--use-hm3-quick-liftover`
  Plain-English meaning: for an HM3-restricted `chr_pos` build, emit the
  opposite genome build using the packaged curated HM3 coordinate map.
  Optional: yes.
  Recommended usage: use this only with `--use-hm3-snps` when a chain-file
  liftover is unnecessary for an HM3-only reference panel.

- `--snp-identifier`
  Plain-English meaning: how to interpret the SNP restriction file supplied by
  `--ref-panel-snps-file`.
  Optional: conditional.
  Accepted values: `rsid`, `chr_pos`.
  Recommended usage: use `rsid` for one-column rsID/dbSNP-style lists, and use
  `chr_pos` for one-column `CHR:POS` lists or tables with `CHR` and `POS`
  columns.

- `--keep-indivs-file`
  Plain-English meaning: restrict the individuals used to compute LD.
  Optional: yes.
  Recommended usage: use a PLINK-style keep file when you want population-specific LD or want to exclude certain samples.

- `--snp-batch-size`
  Plain-English meaning: number of SNPs loaded per pairwise-R2 computation batch.
  Optional: yes; default is `128`.
  Recommended usage: keep the default unless you are tuning memory and throughput on a large machine. Larger values may improve throughput but use more memory.

Coordinate duplicate handling has no public flag. In `chr_pos`-family modes,
source and target coordinate collision groups are always handled with `drop-all`,
and the dropped rows are recorded under `diagnostics/dropped_snps/`. In source-only
`rsid`-family builds, coordinate duplicate filtering is not applicable because
SNP labels define row identity.

- `--log-level`
  Plain-English meaning: how much progress logging to print.
  Optional: yes; default is `INFO`.
  Recommended usage: use `INFO` for routine runs and `DEBUG` only when you are diagnosing behavior.

## Canonical Run: CLI

If `ldsc` is installed as a command-line tool, run:

```bash
ldsc build-ref-panel \
  --plink-prefix resources/example_1kg_30x/genomes_30x_chr22 \
  --source-genome-build hg38 \
  --genetic-map-hg19-sources resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg19_withX.txt \
  --genetic-map-hg38-sources resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg38_withX.txt \
  --liftover-chain-hg38-to-hg19-file resources/liftover/hg38ToHg19.over.chain \
  --ld-wind-cm 1.0 \
  --output-dir tutorial_outputs/ref_panel_chr22
```

If you are running directly from a source checkout instead of an installed CLI, use:

```bash
PYTHONPATH=src python -m ldsc build-ref-panel \
  --plink-prefix resources/example_1kg_30x/genomes_30x_chr22 \
  --source-genome-build hg38 \
  --genetic-map-hg19-sources resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg19_withX.txt \
  --genetic-map-hg38-sources resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg38_withX.txt \
  --liftover-chain-hg38-to-hg19-file resources/liftover/hg38ToHg19.over.chain \
  --ld-wind-cm 1.0 \
  --output-dir tutorial_outputs/ref_panel_chr22
```

What this command is doing:

- reads the chr22 PLINK panel
- treats the input coordinates as hg38
- uses a 1 cM LD window
- uses the explicit hg38->hg19 liftover chain to populate the hg19 coordinates
- interpolates cM values from the provided hg19 and hg38 genetic maps
- writes a build-separated R2 parquet panel rooted at `tutorial_outputs/ref_panel_chr22`
- writes `tutorial_outputs/ref_panel_chr22/build-ref-panel.chr22.log`

This example uses `--ld-wind-cm` and emits both hg38 and hg19 because a matching
liftover chain is provided, so both genetic maps are required.

## Canonical Run: Python API

The most explicit Python API uses the public config object and builder class.
This direct class API writes only data artifacts; use `run_build_ref_panel(...)`
or the CLI when you want the workflow log as well:

```python
from ldsc import GlobalConfig, ReferencePanelBuildConfig, ReferencePanelBuilder, set_global_config

GLOBAL_CONFIG = GlobalConfig(snp_identifier="chr_pos_allele_aware", log_level="INFO")
set_global_config(GLOBAL_CONFIG)

config = ReferencePanelBuildConfig(
    plink_prefix="resources/example_1kg_30x/genomes_30x_chr22",
    source_genome_build="hg38",
    genetic_map_hg19_sources="resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg19_withX.txt",
    genetic_map_hg38_sources="resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg38_withX.txt",
    liftover_chain_hg38_to_hg19_file="resources/liftover/hg38ToHg19.over.chain",
    output_dir="tutorial_outputs/ref_panel_chr22",
    ld_wind_cm=1.0,
)

result = ReferencePanelBuilder(global_config=GLOBAL_CONFIG).run(config)

print(result.chromosomes)
print(result.output_paths["r2_hg38"][0])
print(result.output_paths["r2_hg19"][0])
print(result.output_paths["meta_hg38"][0])
```

When you use the lower-level `ReferencePanelBuilder` API with
`ReferencePanelBuildConfig(ref_panel_snps_file=...)`, put the restriction-file
identifier mode on the injected `GlobalConfig`. For the packaged HapMap3 map,
use `ReferencePanelBuildConfig(use_hm3_snps=True)` instead of a custom file
path. In `chr_pos`-family modes, `GlobalConfig.genome_build` is ignored by this
builder; the restriction file must be in the same build as
`source_genome_build`:

```python
GLOBAL_CONFIG = GlobalConfig(
    snp_identifier="chr_pos_allele_aware",
    log_level="INFO",
)
```

The convenience wrapper shown later reads `snp_identifier` from the registered
`GlobalConfig`; set it before calling the wrapper.

## Output Format

### File tree

For the chr22 example above, the output tree looks like:

```text
tutorial_outputs/ref_panel_chr22/
├── diagnostics/
│   ├── metadata.json          # diagnostic provenance only
│   ├── build-ref-panel.chr22.log
│   └── dropped_snps/          # audit files, header-only when no SNPs dropped
│       └── chr22_dropped.tsv.gz
├── hg19/
│   ├── chr22_r2.parquet
│   └── chr22_meta.tsv.gz
└── hg38/
    ├── chr22_r2.parquet
    └── chr22_meta.tsv.gz
```

If you omit the matching liftover chain, only the source-build directory is
written. If you use `--ld-wind-snps` or `--ld-wind-kb`, you may omit genetic
maps entirely; emitted metadata sidecars will keep their rows and write missing
`CM` values.

For a genome-wide build, the same pattern repeats once per chromosome.

### `r2`: pairwise R2 parquet

This file has one row per unordered SNP pair that falls inside the chosen LD window.
It is written in the canonical row-group-prunable format used by the parquet R2
backend. The writer requires PyArrow so it can set schema metadata and row-group
size explicitly.

Columns:

- `CHR`
- `POS_1`
- `POS_2`
- `SNP_1`
- `SNP_2`
- `R2`

The physical schema is:

- `CHR`, `SNP_1`, `SNP_2`: string
- `POS_1`, `POS_2`: int64
- `R2`: float32

The parquet schema metadata includes:

- `ldsc:sorted_by_build`: the emitted genome build used for `POS_1` and `POS_2`
- `ldsc:row_group_size`: the intended row-group size, defaulting to `50000`
- `ldsc:n_samples`: the number of PLINK samples used to compute LD
- `ldsc:r2_bias`: currently `unbiased` for package-built panels

Example rows from the same chr22 build:

```text
CHR	POS_1	POS_2	SNP_1	SNP_2	R2
22	10684250	10684299	22:10684250:C:G	22:10684302:C:A	-0.0002415361
22	10684250	10685981	22:10684250:C:G	22:10685981:G:A	-0.0002331376
```

Interpretation notes:

- `R2` is the unbiased estimator used by LDSC-style workflows
- because it is unbiased, very weak LD can produce slightly negative values near zero
- downstream readers auto-load the bias state and sample size from schema metadata
- `POS_1` and `POS_2` are in the emitted build recorded in `ldsc:sorted_by_build`
- rows are sorted by non-decreasing `POS_1`; `POS_2` ordering within equal `POS_1` is not required
- legacy columns such as `hg19_pos_1`, `hg38_pos_1`, `Dprime`, and `+/-corr` are intentionally not written

### `meta`: LDSC runtime sidecars

These are gzip-compressed tab-separated files with one row per retained SNP.
Only emitted builds have sidecars; source-build-only runs write just the source
build sidecar.

Schema:

- `CHR`
- `POS`
- `SNP`
- `CM`
- `MAF`

Example `hg19/chr22_meta.tsv.gz` rows:

```text
CHR	POS	SNP	CM	MAF
22	17383676	22:10684250:C:G	2.98554	0.0065584
22	17383725	22:10684302:C:A	2.98571	0.00265459
22	17385403	22:10685981:G:A	2.99174	0.0029669
```

Example `hg38/chr22_meta.tsv.gz` rows:

```text
CHR	POS	SNP	CM	MAF
22	10684250	22:10684250:C:G	0	0.0065584
22	10684299	22:10684302:C:A	0	0.00265459
22	10685981	22:10685981:G:A	0	0.0029669
```

Use the build directory from the coordinate system you want downstream. For
parquet-backed LD-score calculation, the sidecar is the authoritative
reference-panel SNP universe when present. If it is absent, the loader falls
back to the SNPs that appear as R2 endpoints, but cannot recover SNPs with no
off-diagonal R2 pairs and cannot provide CM or MAF metadata.

## Downstream Use

### Load the outputs in Python

For small files, `pandas.read_parquet(...)` is enough. For large R2 files, `pyarrow.parquet.ParquetFile(...)` is often a better first step because it lets you read row groups incrementally.

```python
import pandas as pd
import pyarrow.parquet as pq

meta_hg38 = pd.read_csv(
    "tutorial_outputs/ref_panel_chr22/hg38/chr22_meta.tsv.gz",
    sep="\t",
)

r2_file = pq.ParquetFile("tutorial_outputs/ref_panel_chr22/hg38/chr22_r2.parquet")
first_row_group = r2_file.read_row_group(0).to_pandas()

print(meta_hg38.head())
print(first_row_group.head())
```

### Use the panel in LDSC-style downstream analysis

After you have built a genome-wide parquet panel, feed the build-specific
directory into `ldsc ldscore`.

Example:

```bash
ldsc ldscore \
  --output-dir tutorial_outputs/ldscores_from_parquet_panel \
  --r2-dir "tutorial_outputs/ref_panel/hg38" \
  --snp-identifier rsid \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
```

This ordinary unpartitioned example omits baseline and query annotations, so
`ldsc ldscore` writes a synthetic all-ones `base` baseline column. Add explicit
baseline annotations when computing partitioned or query LD scores.

Use the `hg19/` files instead of the `hg38/` files when your downstream
coordinate system is hg19/GRCh37. In source-build-only runs, choose the single
emitted build directory and keep downstream coordinate settings aligned with
that build.

### Other downstream tasks

These files are also useful outside of LDSC proper. Common examples include:

- inspecting local LD around lead GWAS hits
- building custom SNP neighborhoods inside a physical or genetic window
- comparing hg19 and hg38 coordinates for the same reference panel
- computing custom summaries on retained SNPs, allele labels, or MAF
- feeding pairwise R2 edges into custom fine-mapping or graph-style workflows

## Appendix: Other Supported Options

### Build a chromosome suite instead of one chromosome

If your PLINK files are split by chromosome, pass the suite token through `--plink-prefix`:

```bash
ldsc build-ref-panel \
  --plink-prefix data/reference/genomes_30x_chr@ \
  --source-genome-build hg38 \
  --genetic-map-hg19-sources resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg19_withX.txt \
  --genetic-map-hg38-sources resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg38_withX.txt \
  --liftover-chain-hg38-to-hg19-file resources/liftover/hg38ToHg19.over.chain \
  --ld-wind-cm 1.0 \
  --output-dir tutorial_outputs/ref_panel
```

If `--output-dir` does not exist yet, the workflow warns once and creates it automatically.
If candidate files such as `hg38/chr22_r2.parquet` or
`hg38/chr22_meta.tsv.gz` already exist, the build fails before processing
chromosomes unless `--overwrite` is supplied.

### Restrict to a predefined SNP universe

Use `--use-hm3-snps` when you want the packaged HapMap3 universe. Use
`--ref-panel-snps-file` when you want to keep another specific SNP set.

Accepted forms include:

- tables with an `SNP`/`rsID`-style column in `rsid`-family modes
- tables with `CHR` and `POS` columns in `chr_pos`-family modes
- tables with build-specific columns such as `hg19_POS` and `hg38_POS`; in
  `chr_pos`-family modes, the builder reads the column matching the source PLINK
  build

Allele-free restrictions match by base key, including in allele-aware modes.
Allele-bearing restrictions in allele-aware modes match by the effective
allele-aware key after restriction cleanup. Duplicate restriction keys collapse
to one retained key, and non-identity columns such as `CM` or `MAF` are ignored
instead of being carried into emitted metadata.

Example:

```bash
ldsc build-ref-panel \
  --plink-prefix data/reference/genomes_30x_chr@ \
  --source-genome-build hg38 \
  --use-hm3-snps \
  --snp-identifier rsid \
  --ld-wind-kb 1000 \
  --output-dir tutorial_outputs/ref_panel_hm3
```

For `chr_pos` restriction files, provide source-build coordinates. If both
`hg19_POS` and `hg38_POS` are present, the builder selects the column matching
the explicit or inferred source PLINK build. If only generic `POS` is present,
the builder infers that restriction file's build and errors if it differs from
the source PLINK build.

If you use `--snp-identifier rsid`, omit matching liftover chains. The builder
will emit source-build-only outputs and will not apply coordinate duplicate
filtering because rsID labels, not coordinates, define row identity.

### Restrict the sample set

Use `--keep-indivs-file` when you want LD computed from only a subset of individuals, for example one ancestry group or a QC-passed subset.

### Apply an MAF filter

Use `--maf-min 0.01` or another threshold when you explicitly want to exclude very rare variants from the built panel.

Remember that the canonical examples in this tutorial intentionally do not set `--maf-min`, so they use the full retained SNP set by default.

### Choose a different LD window

All three window modes are supported:

- `--ld-wind-cm`
- `--ld-wind-kb`
- `--ld-wind-snps`

For LDSC-like reference panels, `--ld-wind-cm 1.0` is the most natural starting point.

### Use the convenience Python wrapper

If you prefer a thinner Python wrapper around the CLI-style arguments, the package also exports `run_build_ref_panel(...)`:

```python
from ldsc import GlobalConfig, run_build_ref_panel, set_global_config

GLOBAL_CONFIG = GlobalConfig(snp_identifier="chr_pos", log_level="INFO")
set_global_config(GLOBAL_CONFIG)

result = run_build_ref_panel(
    plink_prefix="resources/example_1kg_30x/genomes_30x_chr22",
    source_genome_build="hg38",
    genetic_map_hg19_sources="resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg19_withX.txt",
    genetic_map_hg38_sources="resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg38_withX.txt",
    liftover_chain_hg38_to_hg19_file="resources/liftover/hg38ToHg19.over.chain",
    output_dir="tutorial_outputs/ref_panel_chr22",
    ld_wind_cm=1.0,
    # overwrite=True,  # replaces current panel candidates but does not clean stale optional siblings
)
```

If you add `ref_panel_snps_file=...` or `use_hm3_snps=True` to this wrapper
call, set one of the exact public modes (`rsid`, `rsid_allele_aware`,
`chr_pos`, or `chr_pos_allele_aware`) on the registered `GlobalConfig`. In
`chr_pos`-family modes, explicit restriction files must be in the source PLINK
build; `GlobalConfig.genome_build` is ignored by
`build-ref-panel`.

### A note on coordinate systems

The builder keeps an in-memory reference SNP table while constructing the R2
rows, but it does not persist an annotation parquet. When a matching liftover
chain is provided in a `chr_pos`-family mode, it emits one R2 parquet and one metadata sidecar per build;
otherwise it emits only the source build. Each R2 parquet stores positions from
its own build in `POS_1`/`POS_2`, with that build recorded in
`ldsc:sorted_by_build`. That makes the panel easier to reuse across projects,
but it also means you should stay intentional about which build you use
downstream:

- use the R2 parquet and metadata sidecar from the same build directory
- if you match SNPs by chromosome-position instead of by `rsID`, make sure the build is consistent all the way through
- if you later load this panel through `RefPanelLoader`, keep the downstream
  `GlobalConfig.genome_build` aligned with the chosen sidecar build so LD-score
  calculation does not fail on an explicit build mismatch
