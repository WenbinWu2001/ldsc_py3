# Build and query R² matrices

Last updated on: 2026-10-28

`ldsc build-r2-matrix` computes the unbiased R² between pairs of SNPs within an LD window from the genotypes of a PLINK reference panel and saves them as an R² matrix.
An R² matrix is a reusable reference panel: build it once, then use it with `ldsc ldscore --r2-matrix-dir` to compute LD scores for any annotations and traits, or with `ldsc query-r2-matrix` to look up the LD between specific SNP pairs.
The R² matrix describes LD in the reference individuals you supply, so choose individuals from the population whose LD you need.

## What you need

- PLINK genotypes: a complete `.bed`, `.bim`, `.fam` trio per prefix, one prefix or a chromosome suite.
- The genome build of the PLINK coordinates, or `--source-genome-build auto` to infer it.
- One LD window, in SNPs, kb or cM.
- Optional: genetic maps (required for a cM window), a liftover chain file (to also write the other genome build), a SNP restriction file, and a list of individuals to keep.

## Try it on the bundled fixture

A source checkout of the repository includes a 32-SNP hg38 chromosome-22 PLINK fixture.
It demonstrates the commands; it is not a reference panel for scientific analysis.
Run this from the repository root:

```bash
ldsc build-r2-matrix \
  --source-genome-build hg38 \
  --plink-prefix tests/fixtures/minimal_external_resources/plink/hm3_chr22_subset \
  --ld-wind-snps 10 \
  --output-dir tutorial_outputs/r2_matrix_chr22
```

The same arguments work with `python -m ldsc build-r2-matrix`.
This example supplies no genetic map and no liftover chain, so it writes hg38 only and stores missing `CM` values in the sidecar.
Use a new output directory, or add `--overwrite` to replace this command's earlier outputs.

## Build an R² matrix for your analysis

A genome-wide R² matrix with a 1 cM LD window, written in both genome builds, looks like this:

```bash
ldsc build-r2-matrix \
  --source-genome-build hg38 \
  --plink-prefix "/path/to/resources/1000G_EUR/chr@" \
  --ld-wind-cm 1.0 \
  --genetic-map-hg19-sources "/path/to/resources/genetic_maps/hg19/chr@.txt" \
  --genetic-map-hg38-sources "/path/to/resources/genetic_maps/hg38/chr@.txt" \
  --liftover-chain-hg38-to-hg19-file /path/to/resources/liftover/hg38ToHg19.over.chain \
  --output-dir results/r2_matrix
```

`--plink-prefix` accepts an exact prefix, a quoted glob, a quoted chromosome suite such as `chr@`, or a plain stem such as `1000G.EUR.QC.`.
The chromosome of each trio comes from its `.bim` contents, not its file name, and numeric suffixes such as `.22` stay part of the prefix.
Missing trio members, malformed files and chromosomes supplied by more than one trio are reported together before any computation, with the full list in `diagnostics/plink_input_issues.tsv`.
See [PLINK prefix resolution](../../reference/input-paths.md#plink-prefix-resolution) for every accepted form.

Chromosomes are built one after another in one run.
To spread the work over a cluster, run one command per chromosome prefix with the same `--output-dir`: per-chromosome runs keep their own diagnostics, so they do not collide.
Worker and BLAS-thread settings are described in [parallel execution](../../reference/parallel-execution.md).

### Choose the LD window

Give exactly one of `--ld-wind-snps N`, `--ld-wind-kb KB` or `--ld-wind-cm CM`.
`--ld-wind-snps` counts positions in the retained SNP sequence.
`--ld-wind-cm` needs a genetic map for every genome build the run writes; SNP and kb windows need no map.
A supplied genetic map fills the sidecar `CM` column; the cM column of the `.bim` file is never used as the genetic map.

The R² matrix stores no pairs beyond its LD window.
`ldsc ldscore` later accepts only a window in the same unit and no wider than the build window, so build with the widest window you will need.
See [LD windows in `ldscore`](../../reference/formats/r2-matrices.md#37-ld-windows-and-sidecar-values-in-ldscore).

`--r2-min` defaults to 0, which keeps every pair in the window, including negative unbiased R² values.
A positive `--r2-min` drops smaller values; `ldsc ldscore` then rejects the R² matrix, because it would count the dropped pairs as zero.
Such an R² matrix serves only `query-r2-matrix`, which reports a dropped pair as `r2=NaN` with `status=absent`; an absent pair is not measured zero LD.

### Choose the reference SNPs and individuals

By default every usable SNP in the PLINK files is kept.
Three options restrict the R² matrix:

- `--keep-indivs-file` lists one individual ID (IID) per row.
- `--maf-min` keeps SNPs with MAF at or above the value, in [0, 0.5].
- `--ref-panel-snps-file` is a headered SNP table in the PLINK source genome build.
  Use a `SNP` column in rsID-family modes, or `CHR` and `POS` columns in coordinate-family modes, with `A1` and `A2` when available.
  A restriction without alleles matches by base identity even in an allele-aware mode; duplicate keys collapse and extra columns are ignored.

`--snp-identifier` (default `chr_pos_allele_aware`) sets the SNP identifier mode used for duplicate removal, for matching the restriction file, and for the SNP identity written into the R² matrix.
In coordinate-family modes, SNPs that share a position are all dropped; there is no option to keep one of them.
`ldsc ldscore` later requires the same SNP identifier mode as the R² matrix; `query-r2-matrix` accepts the R² matrix's mode or its base form.
See [SNP identity and genome build](../concepts/snp-identity-and-genome-build.md).

`--source-genome-build auto` (the default) checks every `.bim` file and a coordinate-based restriction file before any SNP is removed.
A small restriction file with generic `CHR`/`POS` columns may give too little evidence; the run then stops and asks for `--source-genome-build hg19` or `hg38`.
A concrete build must agree with any build declared by the inputs; it labels unlabeled coordinates and never converts them.
In a restriction file, a column named for the build, such as `hg38_POS`, takes precedence over generic `POS`, and rsID restrictions need no build check.
See [genome-build inference, checking and overrides](../concepts/snp-identity-and-genome-build.md#genome-build-inference-checking-and-overrides).

A liftover chain that starts from the source build (`--liftover-chain-hg19-to-hg38-file` or `--liftover-chain-hg38-to-hg19-file`) adds the other genome build as a second R² matrix.
Liftover is available only in coordinate-family modes, and each written build needs its own genetic map for a cM window.
SNPs whose source or target positions collide are dropped; see [liftover rules](../../reference/snp-identity-and-genome-build.md#liftover).

### When a chromosome loses all its SNPs

Before reading any genotypes, `build-r2-matrix` applies SNP selection to every chromosome: the SNP restriction, identity cleanup, liftover and duplicate-position filtering.
If this leaves any autosome (chromosomes 1 to 22) without SNPs, the run stops before computing any R².
The error reports every such autosome together, naming the stage that removed its SNPs, and `diagnostics/input_issues.tsv` lists them with a suggested repair.
Genotype filtering comes later: if MAF filtering, monomorphic SNPs or missing genotypes among the kept individuals remove every SNP of an autosome, the run stops when it reaches that chromosome.
Chromosomes X, Y and MT that lose all their SNPs are skipped with a warning and are left out of the R² matrix.

The run stops because an R² matrix missing an autosome would silently change every LD score computed from it.
Check that the restriction file, `--snp-identifier`, `--source-genome-build` and liftover chain match the PLINK input, and that `--maf-min` and `--keep-indivs-file` are what you intended.
`diagnostics/dropped_snps/chr<N>_dropped.tsv.gz` lists the identity-cleanup, liftover and duplicate-position drops for each chromosome.
Do not loosen an intended filter only to make the command finish; a different filter builds a different R² matrix.
See [troubleshooting](../troubleshooting.md#build-r2-matrix-no-reference-panel-artifacts-were-produced).

## Read the outputs

The fixture example writes:

```text
tutorial_outputs/r2_matrix_chr22/
├── hg38/
│   ├── chr22_r2.parquet
│   └── chr22_meta.tsv.gz
└── diagnostics/
    ├── metadata.chr22.json
    ├── build-r2-matrix.chr22.log
    └── dropped_snps/
        └── chr22_dropped.tsv.gz
```

Each genome build gets its own folder (`hg19/`, `hg38/`) holding one pair table and one sidecar per chromosome.
Pass that build folder to `ldscore --r2-matrix-dir` or `query-r2-matrix --r2-matrix-dir`.

**Pair table (`chr<N>_r2.parquet`).**
Each row is one unordered pair of different SNPs within the LD window.
`IDX_1` and `IDX_2` are zero-based row numbers in the matching sidecar, with `IDX_1 < IDX_2`.
`R2` is the unbiased R², $r^2 - (1 - r^2)/(N - 2)$ for $N$ reference individuals, capped at 1 and stored as a 16-bit integer: divide it by 32767 (the footer key `ldsc:r2_scale`) to recover the value.
Negative values are valid estimates and are kept.
`SIGN_R` is true when the Pearson correlation is zero or positive, with alleles oriented as in the sidecar.
The R² of a SNP with itself is 1 and is not stored.
The 16-bit storage rounds each value slightly; LD-score calculation adds the decoded values in double precision but cannot recover the rounding.
How much the rounding moves an LD score depends on the number of pairs and the annotation values, so there is no general per-SNP error bound.

**Sidecar (`chr<N>_meta.tsv.gz`).**
A gzip TSV with `# ldsc:*` comment lines and the columns `CHR`, `POS`, `SNP`, `A1`, `A2`, `CM`, `MAF`, one row per SNP in position order.
Positions are one-based; `A1` is the minor allele and `MAF` is computed in the retained reference individuals.
`SNP` keeps the `.bim` label, which need not be an rsID.
The sidecar defines the pair indices and is bound to its pair table by an identity hash, so keep the two files together and never sort, filter or edit either one.
A pair table without its sidecar cannot be used.
The complete format is in the [R² matrix specification](../../reference/formats/r2-matrices.md).

**Diagnostics.**
A run on one chromosome prefix writes `diagnostics/metadata.chr<N>.json` and `diagnostics/build-r2-matrix.chr<N>.log`.
A run over several chromosomes writes `diagnostics/metadata.json` and `diagnostics/build-r2-matrix.log`.
The metadata JSON records provenance only.
The dropped-SNP audit is written for each processed chromosome and is header-only when nothing was dropped; SNPs removed by the restriction file or by genotype filtering are not listed there.
`--log-level` sets the detail of the log file; it does not add terminal progress output.
With `--overwrite`, a single-chromosome run replaces only that chromosome's files, and a chromosome-suite run replaces the whole R² matrix.
Logs, overwrite rules and the `RUN_FAILED.txt` failure marker are described in [logs and outputs](../../reference/logs-and-outputs.md); the metadata fields are in [artifact metadata](../../reference/formats/artifact-metadata.md#build-r2-matrix).

## Query SNP pairs

`ldsc query-r2-matrix` looks up specific SNP pairs in an R² matrix and returns their unbiased R² and, in allele-aware modes, the signed Pearson correlation in your allele orientation.
It lets you inspect local LD without loading a whole pair table.

Write the pairs as a TSV (or `.csv`) with one row per pair and endpoint columns suffixed `_1` and `_2`.
For the fixture R² matrix, in the default `chr_pos_allele_aware` mode:

```text
CHR_1	POS_1	A1_1	A2_1	CHR_2	POS_2	A1_2	A2_2
22	16406147	G	A	22	16805059	C	T
```

Use positions in the R² matrix's genome build, and the alleles whose dosage correlation you want.
rsID-family modes use `SNP_1` and `SNP_2`; `rsid_allele_aware` also needs both allele pairs.
Common column aliases are accepted, and `--pairs-file -` reads a TSV from standard input.

```bash
ldsc query-r2-matrix \
  --r2-matrix-dir tutorial_outputs/r2_matrix_chr22/hg38 \
  --genome-build hg38 \
  --pairs-file pairs.tsv \
  --output-dir tutorial_outputs/r2_matrix_query_chr22
```

`--r2-matrix-dir` takes a build folder, or an R² matrix root together with `--genome-build`; single Parquet files are not accepted.
`--snp-identifier` defaults to the mode the R² matrix was built with; the only other accepted value is its base form, `rsid` or `chr_pos`, which ignores alleles.
`--genome-build` must agree with the R² matrix's declared build and does not convert positions; without it, all chromosome files of a coordinate-based R² matrix must declare the same build.
Pair positions are read in the R² matrix's build; their build is not inferred separately.

`query_r2_matrix.tsv` keeps your input columns and adds:

| Column | Meaning |
| --- | --- |
| `r2` | Unbiased R² as stored (decoded); can be negative; `NaN` when unavailable |
| `sign_r` | +1 or -1 in your query allele orientation; missing in base modes or when `r2` is `NaN` |
| `r` | Signed Pearson correlation, recovered using the R² matrix's reference sample size |
| `status` | Empty when `r2` is available; otherwise `not_in_matrix`, `cross_chromosome` or `absent` |

`not_in_matrix` means an endpoint did not match an R² matrix SNP, including invalid or ambiguous query alleles in allele-aware modes.
`cross_chromosome` means the endpoints are on different chromosomes.
`absent` means the pair was not stored, because it lies outside the build window or below a positive `--r2-min`; it does not mean zero correlation.
A pair of a SNP with itself returns `r2 = 1`.
Swapping the alleles of exactly one endpoint flips `sign_r` and `r` and leaves `r2` unchanged.
`sign_r` is aligned to your alleles, unlike the stored `SIGN_R`, which follows the R² matrix's alleles.
Base modes (`rsid`, `chr_pos`) ignore alleles and return missing `sign_r` and `r`.

The command also writes `diagnostics/metadata.json`, with provenance and a count of each status, and `diagnostics/query-r2-matrix.log`.
`--overwrite` replaces this command's earlier outputs, and `--log-level` sets the detail of the log file.
See the [query reference](../../reference/query-r2-matrix.md) for allele matching and lookup details, and the [command-line reference](../../reference/cli.md#ldsc-query-r2-matrix) for every option.

## Compute LD scores from the R² matrix

Pass the build folder to `ldsc ldscore --r2-matrix-dir` with the R² matrix's SNP identifier mode and a window of the same unit, no wider than the build window:

```bash
ldsc ldscore \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg38 \
  --r2-matrix-dir tutorial_outputs/r2_matrix_chr22/hg38 \
  --ld-wind-snps 10 \
  --output-dir tutorial_outputs/ldscores_chr22
```

Without annotations, `ldscore` computes ordinary unpartitioned LD scores from an all-ones `base` annotation; query annotations require explicit baseline annotations.
A cM window in `ldscore` needs sidecar `CM` values, so build the R² matrix with a genetic map and `--ld-wind-cm`.
For scientific analyses, build the R² matrix from a genome-wide reference population with the window you intend to use.
See [compute LD scores](ld-scores.md).

## Options

| Option | Purpose and default |
| --- | --- |
| `--plink-prefix` | Required PLINK prefix, plain stem, quoted glob or chromosome suite |
| `--output-dir` | Required output directory |
| `--source-genome-build` | Build of the PLINK coordinates: `auto` (default), `hg19` (also `hg37`, `GRCh37`) or `hg38` (also `GRCh38`) |
| `--snp-identifier` | SNP identifier mode: `chr_pos_allele_aware` (default), `chr_pos`, `rsid` or `rsid_allele_aware` |
| `--ld-wind-snps`, `--ld-wind-kb`, `--ld-wind-cm` | Exactly one LD window is required |
| `--genetic-map-hg19-sources`, `--genetic-map-hg38-sources` | Genetic maps for sidecar `CM`; required for a cM window in each written build |
| `--liftover-chain-hg19-to-hg38-file`, `--liftover-chain-hg38-to-hg19-file` | Chain file (exact path) that adds the other build; coordinate-family modes only |
| `--ref-panel-snps-file` | Optional headered SNP restriction in the source build |
| `--keep-indivs-file` | Optional list of individual IDs to keep, one per row |
| `--maf-min` | Optional MAF threshold in [0, 0.5], inclusive |
| `--r2-min` | Unbiased R² below which pairs are dropped; default 0 keeps every pair; a positive value makes the R² matrix usable by `query-r2-matrix` only |
| `--snp-batch-size` | SNPs per genotype batch; default 512; the stored pairs do not depend on it |
| `--overwrite` | Replace this command's existing outputs and remove its obsolete outputs after success |
| `--log-level` | Log-file detail: `DEBUG`, `INFO` (default), `WARNING` or `ERROR` |

The [command-line reference](../../reference/cli.md#ldsc-build-r2-matrix) and `ldsc build-r2-matrix --help` list every option with its accepted path forms.

## Memory and run time

Without `--ref-panel-snps-file`, genotypes are streamed from disk; with it, the selected genotypes are held in memory.
Memory also grows with the number of retained individuals and SNPs, the density of SNPs in the window and `--snp-batch-size`.
Wider windows store more pairs and take longer.
There is no general run-time estimate: time one representative chromosome before scheduling a genome-wide build.

For the Python API, see [Python usage](../python-usage.md#r2-matrices).
