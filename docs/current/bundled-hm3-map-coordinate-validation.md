# Bundled HapMap3 Coordinate Validation

Last validated: 2026-08-05

## Material Passport

- Origin skill: `academic-research-suite` / `experiment-agent`
- Origin mode: `validate`
- Verification status: `ANALYZED`
- Version label: `hm3_coordinate_validation_v1`
- Repository commit at validation: `3e3622376d551b9dee8e2ec7827e13a1302df3f3`

## Purpose

This document records the external validation of the chromosome and base-pair
positions in the bundled HapMap3 map. It also defines a repeatable method for
validating this and other bundled coordinate resources before a release.

This validation covers `CHR`, `hg19_POS`, and `hg38_POS`. It does **not**
validate genetic-map positions (`hg19_CM` or `hg38_CM`), alleles, MAF, HapMap3
membership, or the suitability of a SNP for regression.

## Validated Artifact

| Property | Value |
|---|---|
| Repository path | `src/ldsc/data/hm3_curated_map.tsv.gz` |
| Compressed size | 34,262,175 bytes |
| SHA-256 | `c803754b24b3fbd6bd250242ecd88f00f9090743071e13b3722e57ac67297b43` |
| Rows | 1,215,390 |
| SNP identifiers | 1,215,390 unique rsIDs |
| Chromosomes | Autosomes 1-22 |

The gzip stream passed `gzip -t`. The artifact had no missing or nonpositive
coordinates, duplicate rsIDs, duplicate `(CHR, hg19_POS)` keys, or duplicate
`(CHR, hg38_POS)` keys.

## Official External Sources

### UCSC dbSNP build 151 common tracks

The first pass used the fixed UCSC `snp151Common` database tables:

- [hg19 table](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151Common.txt.gz)
  and [schema](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151Common.sql)
- [hg38 table](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151Common.txt.gz)
  and [schema](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151Common.sql)

| Source | Compressed size | SHA-256 | Rows read |
|---|---:|---|---:|
| hg19 `snp151Common.txt.gz` | 784,441,382 bytes | `1dcaf21d5f85a477b81e29025c8c86537202a5de81bbc8ccda956c04bff39c95` | 14,831,956 |
| hg38 `snp151Common.txt.gz` | 804,963,492 bytes | `439fa417416766485c0962fd83825e987aeee22e0b8798f27b2de307ef94373a` | 15,175,044 |

Both downloaded gzip streams passed `gzip -t`. These checksums identify the
exact fixed-release evidence used in this run.

### Ensembl assembly-specific variation APIs

Exceptions from the UCSC pass were checked against the official Ensembl
variation endpoints:

- GRCh37: `https://grch37.rest.ensembl.org/variation/homo_sapiens`
- GRCh38: `https://rest.ensembl.org/variation/homo_sapiens`
- [Variation endpoint documentation](https://rest.ensembl.org/documentation/info/variation_id)

Requests were sent in batches of 200 rsIDs. Only mappings on chromosomes 1-22
for the requested assembly were used.

### NCBI RefSNP API

Records unresolved or conflicting after Ensembl were adjudicated with the
official NCBI RefSNP endpoint:

`https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{numeric_rsid}`

The comparison selected top-level, non-alt, non-patch chromosome placements
for these exact assembly accessions:

| Build label | Assembly | NCBI assembly accession |
|---|---|---|
| hg19 | GRCh37.p13 | `GCF_000001405.25` |
| hg38 | GRCh38.p14 | `GCF_000001405.40` |

See the [NCBI RefSNP JSON API announcement](https://ncbiinsights.ncbi.nlm.nih.gov/2018/06/15/dbsnp-updates-json-refsnp-report-api/)
for the API and report format.

## Tools and Runtime

- Python 3.9.10
- Python standard-library `csv`, `gzip`, and `json` modules
- `requests` 2.27.1 for Ensembl and NCBI HTTPS requests
- Apple `gzip` 448.80.1 for compressed-stream integrity checks
- OpenSSL SHA-256 for artifact and source identification
- Git for recording the repository commit

No liftover chain was used as validation evidence. Liftover would not be a
fully independent check because liftover was one of the strategies used when
the bundled map was constructed.

## Comparison Method

### 1. Artifact integrity and structural checks

1. Test the bundled gzip stream.
2. Record its size, SHA-256, row count, and column schema.
3. Require unique rsIDs and autosomal chromosomes.
4. Reject missing or nonpositive coordinates.
5. Check each build independently for duplicate `(CHR, POS)` keys.

### 2. Exhaustive UCSC comparison

Every bundled rsID was joined to every record with the same rsID in the
corresponding UCSC table. UCSC BED-style `chromStart` is zero-based, so the
comparison coordinate was:

```text
one_based_position = chromStart + 1
```

A SNP was classified as exact when **any** same-rsID UCSC record had the same
autosomal chromosome and one-based position. This handles rsIDs with multiple
UCSC records without depending on row order. Missing rsIDs and present rsIDs
without an exact placement were retained as exceptions.

The UCSC pass produced:

| Build | Exact | Missing from common track | Present but position differed |
|---|---:|---:|---:|
| hg19 | 1,209,864 | 5,526 | 0 |
| hg38 | 1,209,832 | 5,547 | 11 |

The common track is a fixed dbSNP151 subset, not an exhaustive current dbSNP
catalog. Absence from this track was therefore treated as an exception for
further investigation, not as an error.

### 3. Ensembl exception comparison

The union of build-specific UCSC exceptions contained 5,560 rsIDs. Ensembl
records were compared by rsID, requested assembly, chromosome, and one-based
position. Alternative loci and mappings outside chromosomes 1-22 were not
accepted as primary coordinate evidence.

Ten of the eleven hg38 disagreements in UCSC dbSNP151 were confirmed by both
current Ensembl and current NCBI records to match the bundled map:

`rs12420650`, `rs12574060`, `rs9325911`, `rs9709158`, `rs3751916`,
`rs9475973`, `rs9476035`, `rs7748620`, `rs1861014`, and `rs12668140`.

These were classified as stale positions in the fixed dbSNP151 track rather
than bundled-map errors.

### 4. RefSNP merge and identity resolution

Forty-one rsIDs still required NCBI adjudication. Resolution followed these
rules:

1. Query the original numeric rsID.
2. If the response contains `merged_snapshot_data.merged_into`, follow the
   official target rsID recursively until reaching the current RefSNP record.
3. Detect and stop merge loops defensively.
4. Never infer identity from proximity, matching coordinates, matching
   alleles, or similar names. Only an explicit NCBI merge establishes identity.
5. Select placements whose `seq_id_traits_by_assembly` entry matches the exact
   assembly accession and has `is_top_level=true`, `is_chromosome=true`,
   `is_alt=false`, and `is_patch=false`.
6. Convert the RefSNP SPDI position from zero-based to one-based with
   `position + 1`.
7. Classify a record without a qualifying placement as externally
   unverifiable for that build, not automatically incorrect.

This process fetched 77 RefSNP records after merge targets and cache reuse were
accounted for. Among the 41 queried original rsIDs, 38 hg19 and 36 hg38
placements matched exactly. Some exact resolutions involved retired rsIDs; for
example, `rs10493050` resolved through the official merge to `rs2377924`.

### 5. Final classification precedence

For each build, evidence was applied in this order:

1. Exact UCSC same-rsID placement.
2. Exact Ensembl requested-assembly placement for a UCSC exception.
3. Exact NCBI primary placement, after official merge resolution, for a record
   unresolved or conflicting in Ensembl.
4. `confirmed_error` only when current authoritative placement evidence
   disagreed with the bundle.
5. `externally_unverifiable` when no current qualifying placement existed.

No nearest-position matching or chain liftover was used to make an exception
appear exact.

## Validation Results

| Build | Exact | Externally unverifiable | Confirmed errors |
|---|---:|---:|---:|
| hg19 / GRCh37 | 1,215,388 | 2 | 0 |
| hg38 / GRCh38 | 1,215,386 | 3 | 1 |

All externally verifiable chromosome assignments were correct. The sole
confirmed error was a base-pair position, not a chromosome mismatch.

### Confirmed error

| SNP | Build | Bundled coordinate | Authoritative coordinate |
|---|---|---|---|
| `rs28542093` | hg38 | chr9:67,107,264 | chr9:42,080,196 |

NCBI maps this record to `NC_000009.12` position 42,080,196 on GRCh38.p14.
Ensembl and the fixed UCSC dbSNP151 hg38 track agree. The bundled hg19
coordinate, chr9:43,734,193, matches the GRCh37.p13 RefSNP placement. No other
bundled row currently occupies hg38 chr9:42,080,196.

Correcting this base-pair position should be accompanied by separate review of
the row's `hg38_CM`; the genetic-map coordinate was outside this validation's
scope.

### Externally unverifiable build-specific records

| SNP | Build | Bundled coordinate | Current-source limitation |
|---|---|---|---|
| `rs12361431` | hg19 | chr11:66,845,209 | Current RefSNP record has unsupported snapshot data and no primary assembly placement. |
| `rs9331726` | hg19 | chr9:136,368,685 | Current RefSNP placement is on obsolete contig `NT_079538.1`, without a qualifying GRCh37 chromosome placement. |
| `rs12361431` | hg38 | chr11:67,077,738 | Current RefSNP record has unsupported snapshot data and no primary assembly placement. |
| `rs16980132` | hg38 | chr22:25,327,333 | Current RefSNP record has unsupported snapshot data and no primary assembly placement. |
| `rs28665796` | hg38 | chr9:41,103,190 | Current RefSNP has a GRCh37 chromosome placement but no primary GRCh38 placement. |

These entries are not proven wrong. They require an explicit release decision:
retain them as legacy coordinates with provenance, replace or retire them, or
find another authoritative archived source.

## Release Revalidation Protocol

Before validating bundled sources after development is complete:

1. Freeze the candidate artifact and record its Git commit, compressed size,
   SHA-256, schema, and row count.
2. Download fixed external source files and retain their URLs, sizes, hashes,
   release identifiers, and access date.
3. Save raw or normalized responses from mutable APIs so the adjudication can
   be audited after those APIs change.
4. Run build-specific comparisons using explicit assembly accessions and
   explicit coordinate-system conversions.
5. Write machine-readable exception tables for every stage, including the
   original identifier, resolved identifier, merge chain, bundled coordinate,
   external placements, evidence source, and status.
6. Maintain a reviewed allowlist for build-specific records that remain
   externally unverifiable. A new unverifiable record should require review;
   it should not silently pass.
7. Treat any confirmed chromosome or base-pair mismatch as a release-blocking
   error unless a documented source-priority decision overrides it.
8. Re-run duplicate-coordinate checks after applying corrections.
9. Validate alleles, genetic-map coordinates, MAF, membership, and biological
   filtering in separate reports; coordinate agreement does not validate those
   fields.
10. Recompute the final counts from build-specific exception sets. Do not use
    counts from the union of exceptions as if every exception affected both
    builds.

## Reproducibility and Limitations

- The fixed UCSC inputs are reproducible from the recorded URLs and hashes.
- Ensembl and NCBI are mutable services. Their responses describe current
  records as accessed on 2026-08-05 and may change in later releases.
- The Ensembl and UCSC variation data ultimately share dbSNP lineage with NCBI;
  agreement among them is useful release/version corroboration but is not
  equivalent to three biologically independent observations.
- The UCSC `snp151Common` table excludes variants outside its common subset, so
  missing records required secondary adjudication.
- The temporary comparison tables and raw API responses from this run were not
  committed to the repository. A final release-validation implementation
  should preserve those artifacts or provide a deterministic script that
  regenerates them.
- Because the current API evidence is mutable and this run was not produced by
  a committed validation pipeline, the Material Passport remains `ANALYZED`
  rather than `VERIFIED`.
