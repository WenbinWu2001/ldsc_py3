# Munge-Sumstats

This document explains the public shape of `ldsc munge-sumstats`: what it does,
what it writes, how genome builds are handled, and how to use `--infer-only`
before spending time on a full run.

## Workflow

`munge-sumstats` converts one raw GWAS summary-statistics file into an
LDSC-ready artifact family. The workflow accepts one exact file path, or a glob
that resolves to exactly one file.

The normal workflow is:

1. Resolve the raw input path and preflight fixed output artifacts.
2. Infer or apply the raw format profile (`plain`, `daner-old`, `daner-new`,
   or `auto`). VCF-style headers are handled by `plain`.
3. Infer safe column hints for SNP, chromosome, position, alleles, p-values,
   sample size, INFO, and signed statistics.
4. In coordinate-family modes, resolve `source_genome_build`. The public
   default is `auto`, which infers hg19 versus hg38 from raw `CHR`/`POS`
   coordinates.
5. Load the optional SNP keep-list once. `--use-hm3-snps` uses the packaged
   curated HM3 map; `--sumstats-snps-file` uses a custom headered keep-list.
6. Stream the raw file in chunks, apply core LDSC QC, normalize coordinates,
   apply the SNP restriction, and compute or validate `Z` and `N`.
7. Clean identity keys for the active `--snp-identifier` mode. Duplicate
   effective identities are dropped as whole duplicate groups.
8. If the resolved source build differs from the requested output build, apply
   exactly one liftover method.
9. Write fixed outputs, root metadata, and diagnostic audit files.

The default identity mode is `chr_pos_allele_aware`, so `A1/A2` are required by
default. Use `--snp-identifier chr_pos` for coordinate identity without
allele-aware matching, or `--snp-identifier rsid` for rsID-only identity.

## Genome-Build Contract

Coordinate-family modes separate the raw input build from the final output
build:

| Concept | CLI flag | Python config | Meaning |
| --- | --- | --- | --- |
| Source build | `--source-genome-build` | `source_genome_build` | Build of raw `CHR`/`POS`; defaults to `auto`. |
| Output build | `--output-genome-build` | `output_genome_build` | Final build written to munged sumstats; required in coordinate-family modes. |

Accepted concrete builds are `hg19` and `hg38`; common aliases such as `hg37`,
`GRCh37`, and `GRCh38` normalize to canonical values. `auto` is accepted only
for the source build.

In rsID-family modes, genome build is not part of the merge identity. The
munger rejects source-build, output-build, and liftover flags in those modes.
`--use-hm3-snps` remains valid because it is a row filter. The written
`metadata.json` and in-memory `SumstatsTable.config_snapshot` store
`genome_build=null` for rsID-family artifacts.

## Output Artifacts

The output directory is a coherent artifact family:

```text
<output_dir>/
  metadata.json
  sumstats.parquet
  sumstats.sumstats.gz
  diagnostics/
    sumstats.log
    dropped_snps/
      dropped.tsv.gz
```

`sumstats.parquet` is the default data artifact. `sumstats.sumstats.gz` is
written only with `--output-format tsv.gz` or `--output-format both`.

Without `--overwrite`, any existing owned artifact blocks the run before output
files are opened. With `--overwrite`, the workflow replaces current-run outputs
and removes stale owned siblings not produced by the successful run. Unrelated
files in the output directory are preserved.

## Root Metadata Schema

Root `metadata.json` is a downstream compatibility contract. It is intentionally
small:

| Field | Meaning |
| --- | --- |
| `schema_version` | Artifact metadata schema version. |
| `artifact_type` | Always `sumstats`. |
| `files` | Relative paths for written data artifacts, usually `sumstats.parquet` and optionally `sumstats.sumstats.gz`. |
| `snp_identifier` | Identity mode used when writing the artifact. |
| `genome_build` | Final output build in coordinate-family modes; `null` in rsID-family modes. |
| `trait_name` | Optional biological trait label supplied by the user. |

The same compatibility state is attached to in-memory results as
`SumstatsTable.config_snapshot`. Downstream regression checks `snp_identifier`
and `genome_build` compatibility against LD-score artifacts before merging.

Detailed source-build inference, liftover decisions, method names, and drop
counts live in `diagnostics/sumstats.log`, not root metadata.

## Dropped-SNP Audit

`diagnostics/dropped_snps/dropped.tsv.gz` is an always-owned audit sidecar. It
records identity-cleanup and liftover drops when they occur, and is still
written for clean runs so a missing file is never confused with zero drops.

The audit table uses fields such as `CHR`, `SNP`, `source_pos`, `target_pos`,
`reason`, `base_key`, `identity_key`, `allele_set`, and `stage` when those
values are available. Downstream LDSC commands do not consume this file.

## `--infer-only`

`--infer-only` is the quick pre-run diagnosis mode. It reads the raw header and
a small coordinate sample, prints the inferred configuration, and writes no
artifacts. It does not require `--output-dir`.

In coordinate-family modes, `--infer-only` also requires
`--output-genome-build` because liftover status cannot be diagnosed without the
requested final build.

Example:

```bash
ldsc munge-sumstats \
  --raw-sumstats-file data/trait.tsv.gz \
  --snp-identifier chr_pos_allele_aware \
  --output-genome-build hg38 \
  --infer-only
```

The report includes:

- detected raw format
- inferred column hints
- INFO-list handling
- missing required fields
- resolved source genome build, if inference succeeds
- requested output genome build
- whether liftover is required
- whether the command is runnable
- suggested command fragments for repair

If source-build inference fails, the report is non-runnable and suggests
rerunning with `--source-genome-build hg19` or `--source-genome-build hg38`.

If liftover is required but no method is supplied, the report suggests both
valid options:

```bash
--use-hm3-snps --use-hm3-quick-liftover
```

or:

```bash
--liftover-chain-file <hg19ToHg38.over.chain>
```

The chain-file label changes with direction. If a chain file is supplied to
`--infer-only`, the workflow checks that the path exists and reports the
expected source -> output direction; it does not parse or apply the chain.

## HM3 Filter and Quick Liftover

`--use-hm3-snps` is a filter only. It restricts the summary-statistics rows to
the packaged curated HM3 map while chunks are streaming. It does not imply
liftover by itself.

`--use-hm3-quick-liftover` uses the packaged dual-build HM3 coordinate map for a
coordinate-only liftover. It always requires `--use-hm3-snps`, because the quick
map is defined for that filtered HM3 universe.

Example:

```bash
ldsc munge-sumstats \
  --raw-sumstats-file data/trait.tsv.gz \
  --snp-identifier chr_pos \
  --source-genome-build auto \
  --output-genome-build hg38 \
  --use-hm3-snps \
  --use-hm3-quick-liftover \
  --output-dir outputs/trait_hg38
```

Quick liftover changes only `CHR` and `POS`. The `SNP` column remains a label,
so rsIDs are not rewritten.

## Liftover Rules

Liftover is valid only in coordinate-family modes.

Allowed specifications:

| Situation | Required behavior |
| --- | --- |
| Resolved source build equals output build | No liftover is applied. If a liftover method was supplied, the workflow warns and ignores it. |
| Resolved source build differs from output build | Supply exactly one of `--use-hm3-quick-liftover` or `--liftover-chain-file`. |
| Both liftover methods are supplied | Hard error. |
| `--use-hm3-quick-liftover` is supplied without `--use-hm3-snps` | Hard error. |
| rsID-family mode receives source/output/liftover flags | Hard error. |

Chain-file mode expects a chain in the resolved source -> requested output
direction. For example, source `hg19` and output `hg38` expects an
`hg19ToHg38` chain.

During liftover, the workflow drops rows that cannot be mapped cleanly,
including missing coordinates, duplicate source coordinate groups, unmapped
coordinates, cross-chromosome mappings, and duplicate target coordinate groups.
Counts are written to `diagnostics/sumstats.log`; row-level examples live in
the dropped-SNP audit sidecar.

## Common Commands

Diagnose first:

```bash
ldsc munge-sumstats \
  --raw-sumstats-file data/trait.tsv.gz \
  --output-genome-build hg38 \
  --infer-only
```

Run with inferred source build and no liftover if the raw file is already hg38:

```bash
ldsc munge-sumstats \
  --raw-sumstats-file data/trait.tsv.gz \
  --source-genome-build auto \
  --output-genome-build hg38 \
  --output-dir outputs/trait_hg38
```

Run HM3-filtered quick liftover:

```bash
ldsc munge-sumstats \
  --raw-sumstats-file data/trait.tsv.gz \
  --source-genome-build auto \
  --output-genome-build hg38 \
  --use-hm3-snps \
  --use-hm3-quick-liftover \
  --output-dir outputs/trait_hg38
```

Run chain-file liftover:

```bash
ldsc munge-sumstats \
  --raw-sumstats-file data/trait.tsv.gz \
  --source-genome-build hg19 \
  --output-genome-build hg38 \
  --liftover-chain-file resources/liftover/hg19ToHg38.over.chain \
  --output-dir outputs/trait_hg38
```

Run rsID identity:

```bash
ldsc munge-sumstats \
  --raw-sumstats-file data/trait.tsv.gz \
  --snp-identifier rsid \
  --use-hm3-snps \
  --output-dir outputs/trait_rsid
```

Do not pass source/output/liftover flags in rsID modes.

## Python API

```python
from ldsc import MungeConfig, SumstatsMunger

sumstats = SumstatsMunger().run(
    MungeConfig(
        raw_sumstats_file="data/trait.tsv.gz",
        source_genome_build="auto",
        output_genome_build="hg38",
        use_hm3_snps=True,
        output_dir="outputs/trait_hg38",
    )
)
```

If the inferred source build differs from `output_genome_build`, add exactly one
liftover method:

```python
MungeConfig(
    raw_sumstats_file="data/trait.tsv.gz",
    source_genome_build="auto",
    output_genome_build="hg38",
    use_hm3_snps=True,
    use_hm3_quick_liftover=True,
    output_dir="outputs/trait_hg38",
)
```

or:

```python
MungeConfig(
    raw_sumstats_file="data/trait.tsv.gz",
    source_genome_build="hg19",
    output_genome_build="hg38",
    liftover_chain_file="resources/liftover/hg19ToHg38.over.chain",
    output_dir="outputs/trait_hg38",
)
```
