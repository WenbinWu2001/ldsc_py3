# Build a Standard LD Parquet Reference Panel from PLINK

Goal: start from a PLINK reference panel (`.bed/.bim/.fam`) and build a standard per-chromosome parquet representation of pairwise LD, plus sidecar metadata that can be used in downstream LDSC workflows.

This tutorial is written for first-time users. The examples use the chr22 1000 Genomes 30x example files bundled in this repository, but the same workflow applies to your own PLINK reference panel.

All example paths below are relative to the workspace root that contains both `resources/` and `ldsc_py3_restructured/`.

## What This Builder Produces

The `build-ref-panel` workflow converts PLINK genotypes into three kinds of files for each chromosome:

- one SNP annotation parquet (`ann`): one row per retained SNP
- one LD parquet (`ld`): one row per unordered SNP pair within the chosen LD window
- two runtime metadata sidecars (`meta_hg19` and `meta_hg38`): one row per retained SNP, used by LDSC-style downstream tools

The LD output is a long pairwise table, not a dense square matrix on disk. That is usually the practical format for large reference panels.

By default, the builder keeps all SNPs in the PLINK panel after:

- optional user-requested filters
- automatic liftover sanity filtering

In particular, SNPs are dropped if they fail hg19/hg38 liftover or liftover onto a different chromosome in the target build.

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

The builder copies the PLINK `SNP` field directly into the output `rsID` column. If your BIM uses dbSNP IDs, you will see dbSNP IDs. If your BIM uses IDs like `22:10684250:C:G`, those exact strings are preserved.

### Chromosome selection

There is currently no separate `--chrom` flag.

Chromosome selection comes from the PLINK prefix token you provide through `--bfile`:

- use `--bfile <single-prefix>` when the prefix already points to one chromosome or one specific PLINK dataset
- use `--bfile <suite-token>` when you have one PLINK prefix per chromosome, typically with an explicit `@` token
- use `--bfile <glob-like-token>` only when that token resolves cleanly to one or more complete PLINK prefixes; unlike scalar file inputs elsewhere, a PLINK prefix token is resolved at the prefix level rather than at the individual `.bed/.bim/.fam` file level

Examples:

```text
--bfile resources/example_1kg_30x/genomes_30x_chr22
--bfile data/reference/genomes_30x_chr@
```

If the input resolves to multiple chromosomes, the builder writes one output set per chromosome.

In other words, the old split between `--bfile` and a dedicated per-chromosome flag is gone. The unified `--bfile` argument now handles either one concrete prefix or an explicit chromosome suite such as `panel_chr@`.

### Genetic map files

You must supply:

- one hg19-aligned genetic map
- one hg38-aligned genetic map

The recommended file format is a plain text table with columns:

- `chr`
- `position`
- `Genetic_Map(cM)`

The bundled Alkes-group maps in `resources/genetic_maps/genetic_map_alkesgroup/` already use an accepted format.

## Parameters and Configuration

### Required arguments

- `--bfile`
  Plain-English meaning: where the PLINK panel lives.
  Recommended usage: use `--bfile` for a single chromosome, a single prefix, or a chromosome-split suite such as `panel_chr@`.

- `--panel-label`
  Plain-English meaning: short label that becomes part of the output filenames.
  Recommended usage: keep it short, stable, and filesystem-safe, for example `1KG30X`, `EUR`, or `UKB_EUR`.

- `--source-genome-build`
  Plain-English meaning: which genome build the input PLINK coordinates already use.
  Accepted values: `hg19`, `hg37`, `GRCh37`, `hg38`, `GRCh38`.
  Recommended usage: be explicit. Build auto-inference is not part of this workflow.

- `--genetic-map-hg19`
  Plain-English meaning: genetic map aligned to hg19 coordinates.
  Recommended usage: use the bundled Alkes-group map unless you have a strong reason to substitute your own.

- `--genetic-map-hg38`
  Plain-English meaning: genetic map aligned to hg38 coordinates.
  Recommended usage: same as above; provide the hg38 mate of your hg19 map.

- `--liftover-chain-hg19-to-hg38` or `--liftover-chain-hg38-to-hg19`
  Plain-English meaning: explicit chain file used to translate positions into the other genome build.
  Recommended usage: pass the chain that matches `--source-genome-build`.
  Current behavior: one liftover chain is required for every build-ref-panel run.

- `--out`
  Plain-English meaning: output root directory.
  Recommended usage: point this to a dedicated directory for the new reference panel build.

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

- `--maf`
  Plain-English meaning: minimum minor allele frequency among retained SNPs.
  Optional: yes.
  Recommended usage: leave unset to keep all panel SNPs by default. Use values like `0.01` or `0.05` only if your downstream analysis truly wants a frequency-restricted panel.

- `--ref-panel-snps-path`
  Plain-English meaning: restrict the retained reference-panel SNP rows to a
  user-supplied list.
  Optional: yes.
  Recommended usage: use this when you want to build a panel only on HM3 SNPs,
  a curated common-SNP list, or another pre-defined reference universe.

  The builder auto-detects two common identifier styles:

  - one-column rsID-style lists
  - chromosome-position identifiers such as `22:10684250`, or tables with `CHR` and `POS` columns

- `--keep-indivs`
  Plain-English meaning: restrict the individuals used to compute LD.
  Optional: yes.
  Recommended usage: use a PLINK-style keep file when you want population-specific LD or want to exclude certain samples.

- `--chunk-size`
  Plain-English meaning: block size used by the sliding LD writer.
  Optional: yes; default is `50`.
  Recommended usage: keep the default unless you are tuning memory and throughput on a large machine.

- `--log-level`
  Plain-English meaning: how much progress logging to print.
  Optional: yes; default is `INFO`.
  Recommended usage: use `INFO` for routine runs and `DEBUG` only when you are diagnosing behavior.

## Canonical Run: CLI

If `ldsc` is installed as a command-line tool, run:

```bash
ldsc build-ref-panel \
  --bfile resources/example_1kg_30x/genomes_30x_chr22 \
  --panel-label 1KG30X \
  --source-genome-build hg38 \
  --genetic-map-hg19 resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg19_withX.txt \
  --genetic-map-hg38 resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg38_withX.txt \
  --liftover-chain-hg38-to-hg19 resources/liftover/hg38ToHg19.over.chain \
  --ld-wind-cm 1.0 \
  --out tutorial_outputs/ref_panel_chr22
```

If you are running directly from a source checkout instead of an installed CLI, use:

```bash
PYTHONPATH=src python -m ldsc build-ref-panel \
  --bfile resources/example_1kg_30x/genomes_30x_chr22 \
  --panel-label 1KG30X \
  --source-genome-build hg38 \
  --genetic-map-hg19 resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg19_withX.txt \
  --genetic-map-hg38 resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg38_withX.txt \
  --liftover-chain-hg38-to-hg19 resources/liftover/hg38ToHg19.over.chain \
  --ld-wind-cm 1.0 \
  --out tutorial_outputs/ref_panel_chr22
```

What this command is doing:

- reads the chr22 PLINK panel
- treats the input coordinates as hg38
- uses a 1 cM LD window
- uses the explicit hg38->hg19 liftover chain to populate the hg19 coordinates
- interpolates cM values from both hg19 and hg38 genetic maps
- writes a standard parquet panel rooted at `tutorial_outputs/ref_panel_chr22`

## Canonical Run: Python API

The most explicit Python API uses the public config object and builder class:

```python
from ldsc import GlobalConfig, ReferencePanelBuildConfig, ReferencePanelBuilder, set_global_config

GLOBAL_CONFIG = GlobalConfig(log_level="INFO")
set_global_config(GLOBAL_CONFIG)

config = ReferencePanelBuildConfig(
    panel_label="1KG30X",
    plink_prefix="resources/example_1kg_30x/genomes_30x_chr22",
    source_genome_build="hg38",
    genetic_map_hg19_path="resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg19_withX.txt",
    genetic_map_hg38_path="resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg38_withX.txt",
    liftover_chain_hg38_to_hg19_path="resources/liftover/hg38ToHg19.over.chain",
    output_dir="tutorial_outputs/ref_panel_chr22",
    ld_wind_cm=1.0,
)

result = ReferencePanelBuilder(global_config=GLOBAL_CONFIG).run(config)

print(result.chromosomes)
print(result.output_paths["ann"][0])
print(result.output_paths["ld"][0])
print(result.output_paths["meta_hg38"][0])
```

## Output Format

### File tree

For the chr22 example above, the output tree looks like:

```text
tutorial_outputs/ref_panel_chr22/
└── parquet/
    ├── ann/
    │   └── 1KG30X_chr22_ann.parquet
    ├── ld/
    │   └── 1KG30X_chr22_LD.parquet
    └── meta/
        ├── 1KG30X_chr22_meta_hg19.tsv.gz
        └── 1KG30X_chr22_meta_hg38.tsv.gz
```

For a genome-wide build, the same pattern repeats once per chromosome.

### `ann`: SNP annotation parquet

This file has one row per retained SNP.

Columns:

- `chr`
- `hg19_pos`
- `hg38_pos`
- `hg19_Uniq_ID`
- `hg38_Uniq_ID`
- `rsID`
- `MAF`
- `REF`
- `ALT`

Example rows from a real chr22 build:

```text
chr	hg19_pos	hg38_pos	hg19_Uniq_ID	hg38_Uniq_ID	rsID	MAF	REF	ALT
22	17383676	10684250	22:17383676:G:C	22:10684250:G:C	22:10684250:C:G	0.0065584009993754355	G	C
22	17383725	10684299	22:17383725:A:C	22:10684299:A:C	22:10684302:C:A	0.0026545908806995255	A	C
22	17385403	10685981	22:17385403:A:G	22:10685981:A:G	22:10685981:G:A	0.002966895690193594	A	G
```

Interpretation notes:

- `rsID` is copied directly from the PLINK BIM `SNP` field
- `REF` and `ALT` are currently copied from the PLINK BIM allele columns used by this pipeline
- the `hg19_*` and `hg38_*` columns let you work in either coordinate system downstream

### `ld`: pairwise LD parquet

This file has one row per unordered SNP pair that falls inside the chosen LD window.

Columns:

- `chr`
- `rsID_1`
- `rsID_2`
- `hg38_pos_1`
- `hg38_pos_2`
- `hg19_pos_1`
- `hg19_pos_2`
- `hg38_Uniq_ID_1`
- `hg38_Uniq_ID_2`
- `hg19_Uniq_ID_1`
- `hg19_Uniq_ID_2`
- `R2`
- `Dprime`
- `+/-corr`

Example rows from the same chr22 build:

```text
chr	rsID_1	rsID_2	hg38_pos_1	hg38_pos_2	hg19_pos_1	hg19_pos_2	hg38_Uniq_ID_1	hg38_Uniq_ID_2	hg19_Uniq_ID_1	hg19_Uniq_ID_2	R2	Dprime	+/-corr
22	22:10684250:C:G	22:10684302:C:A	10684250	10684299	17383676	17383725	22:10684250:G:C	22:10684299:A:C	22:17383676:G:C	22:17383725:A:C	-0.00024153611420225257		-
22	22:10684250:C:G	22:10685981:G:A	10684250	10685981	17383676	17385403	22:10684250:G:C	22:10685981:A:G	22:17383676:G:C	22:17385403:A:G	-0.00023313758664503438		-
```

Interpretation notes:

- `R2` is the unbiased estimator used by LDSC-style workflows
- because it is unbiased, very weak LD can produce slightly negative values near zero
- `+/-corr` stores the sign of the underlying correlation
- `Dprime` is written as `NA` in the current implementation

### `meta_hg19` and `meta_hg38`: LDSC runtime sidecars

These are gzip-compressed tab-separated files with one row per retained SNP.

Schema:

- `CHR`
- `POS`
- `SNP`
- `CM`
- `MAF`

Example `meta_hg19` rows:

```text
CHR	POS	SNP	CM	MAF
22	17383676	22:10684250:C:G	2.98554	0.0065584
22	17383725	22:10684302:C:A	2.98571	0.00265459
22	17385403	22:10685981:G:A	2.99174	0.0029669
```

Example `meta_hg38` rows:

```text
CHR	POS	SNP	CM	MAF
22	10684250	22:10684250:C:G	0	0.0065584
22	10684299	22:10684302:C:A	0	0.00265459
22	10685981	22:10685981:G:A	0	0.0029669
```

Use the sidecar that matches the coordinate system you want downstream.

## Downstream Use

### Load the outputs in Python

For small files, `pandas.read_parquet(...)` is enough. For large LD files, `pyarrow.parquet.ParquetFile(...)` is often a better first step because it lets you read row groups incrementally.

```python
import pandas as pd
import pyarrow.parquet as pq

ann = pd.read_parquet("tutorial_outputs/ref_panel_chr22/parquet/ann/1KG30X_chr22_ann.parquet")
meta_hg38 = pd.read_csv(
    "tutorial_outputs/ref_panel_chr22/parquet/meta/1KG30X_chr22_meta_hg38.tsv.gz",
    sep="\t",
)

ld_file = pq.ParquetFile("tutorial_outputs/ref_panel_chr22/parquet/ld/1KG30X_chr22_LD.parquet")
first_row_group = ld_file.read_row_group(0).to_pandas()

print(ann.head())
print(meta_hg38.head())
print(first_row_group.head())
```

### Use the panel in LDSC-style downstream analysis

After you have built a genome-wide parquet panel, you can feed the LD parquet and matching metadata sidecar into `ldsc ldscore`.

Example:

```bash
ldsc ldscore \
  --out tutorial_outputs/ldscores_from_parquet_panel \
  --baseline-annot "annotations/baseline.@.annot.gz" \
  --r2-table "tutorial_outputs/ref_panel/parquet/ld/EUR_chr@_LD.parquet" \
  --frqfile "tutorial_outputs/ref_panel/parquet/meta/EUR_chr@_meta_hg38.tsv.gz" \
  --snp-identifier rsid \
  --ld-wind-cm 1.0
```

Use `meta_hg19` instead of `meta_hg38` when your downstream coordinate system is hg19/GRCh37.

### Other downstream tasks

These files are also useful outside of LDSC proper. Common examples include:

- inspecting local LD around lead GWAS hits
- building custom SNP neighborhoods inside a physical or genetic window
- comparing hg19 and hg38 coordinates for the same reference panel
- computing custom summaries on retained SNPs, allele labels, or MAF
- feeding pairwise LD edges into custom fine-mapping or graph-style workflows

## Appendix: Other Supported Options

### Build a chromosome suite instead of one chromosome

If your PLINK files are split by chromosome, pass the suite token through `--bfile`:

```bash
ldsc build-ref-panel \
  --bfile data/reference/genomes_30x_chr@ \
  --panel-label 1KG30X \
  --source-genome-build hg38 \
  --genetic-map-hg19 resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg19_withX.txt \
  --genetic-map-hg38 resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg38_withX.txt \
  --liftover-chain-hg38-to-hg19 resources/liftover/hg38ToHg19.over.chain \
  --ld-wind-cm 1.0 \
  --out tutorial_outputs/ref_panel
```

If `--out` does not exist yet, the workflow warns once and creates it automatically.

### Restrict to a predefined SNP universe

Use `--ref-panel-snps-path` when you want to keep only a specific SNP set, for
example HapMap3 or another curated reference-panel universe.

Accepted forms include:

- plain one-column SNP lists
- one-column `chr:pos` lists
- tables with `CHR` and `POS` columns

### Restrict the sample set

Use `--keep-indivs` when you want LD computed from only a subset of individuals, for example one ancestry group or a QC-passed subset.

### Apply an MAF filter

Use `--maf 0.01` or another threshold when you explicitly want to exclude very rare variants from the built panel.

Remember that the canonical examples in this tutorial intentionally do not set `--maf`, so they use the full retained SNP set by default.

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

GLOBAL_CONFIG = GlobalConfig(log_level="INFO")
set_global_config(GLOBAL_CONFIG)

result = run_build_ref_panel(
    bfile="resources/example_1kg_30x/genomes_30x_chr22",
    panel_label="1KG30X",
    source_genome_build="hg38",
    genetic_map_hg19="resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg19_withX.txt",
    genetic_map_hg38="resources/genetic_maps/genetic_map_alkesgroup/genetic_map_hg38_withX.txt",
    liftover_chain_hg38_to_hg19="resources/liftover/hg38ToHg19.over.chain",
    out="tutorial_outputs/ref_panel_chr22",
    ld_wind_cm=1.0,
)
```

### A note on coordinate systems

The builder stores both hg19 and hg38 positions in the output parquet files. That makes the panel easier to reuse across projects, but it also means you should stay intentional about which build you use downstream:

- use the metadata sidecar that matches your downstream coordinate system
- if you match SNPs by chromosome-position instead of by `rsID`, make sure the build is consistent all the way through
- if you later load this panel through `RefPanelLoader`, keep the downstream
  `GlobalConfig.genome_build` aligned with the chosen sidecar build so LD-score
  calculation does not fail on an explicit build mismatch
