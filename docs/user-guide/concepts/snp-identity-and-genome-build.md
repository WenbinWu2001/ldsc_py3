# SNP identity and genome build

Last updated on: 2026-11-06

Every command that matches SNPs across files uses one SNP identifier mode, and every command that reads positions uses one genome build.
Choose both once, before you munge summary statistics or compute LD scores, and use the same choice for every input of an analysis.
The exact contract is in [SNP identity and genome build](../../reference/snp-identity-and-genome-build.md).

## SNP identifier modes

The SNP identifier mode (`--snp-identifier`) decides which fields make up a SNP identity; allele-aware modes match the two alleles as an unordered, strand-aware pair:

| Mode | SNP identity |
| --- | --- |
| `rsid` | rsID |
| `rsid_allele_aware` | rsID and the two alleles |
| `chr_pos` | chromosome and one-based position |
| `chr_pos_allele_aware` (default) | chromosome, one-based position and the two alleles |

The coordinate-family modes (`chr_pos`, `chr_pos_allele_aware`) need a genome build; the rsID-family modes do not.
The guided tutorial selects the base mode `chr_pos`.
Column requirements and the duplicate policy for each mode are in the [column schema](../../reference/formats/column-schema.md#snp-identifier-modes).

## Inspect an artifact's genome build

Before choosing inputs, check the build that a file or LDSC3 output directory records or implies:

```bash
ldsc infer-genome-build reference.chr22.bim
ldsc infer-genome-build results/ldscore
```

The command takes one exact file or recognized artifact directory, writes nothing, and reports `unknown` (exit 0) when the evidence is insufficient.
It does not replace the checks that each workflow runs; see the [inspection contract](../../reference/snp-identity-and-genome-build.md#standalone-build-inspection).

## Genome-build inference, checking, and overrides

With `--genome-build auto` (`--source-genome-build auto` for `build-r2-matrix` and `munge-sumstats`), each command reads build metadata where inputs have it and infers the build from a bounded sample of positions where they do not.
A concrete `hg19` or `hg38` must agree with every input that declares a build; it labels unlabeled inputs but never converts coordinates.
Coordinate-family munging also needs `--output-genome-build`, and indexed LD scoring inherits the index's stored build.
Inputs that disagree stop the run before computation, with every checked input grouped by build in one error.
Repair the inputs (matching coordinates, verified metadata, or an explicit liftover); the build option cannot override a conflicting declaration.

See the [full policy and command scopes](../../reference/snp-identity-and-genome-build.md#genome-build-policy), the [keep-list rules](../../reference/snp-identity-and-genome-build.md#snp-keep-list-rules), and the [grouped error examples](../troubleshooting.md#common-genome-build-consistency-check-failed).
Path patterns and output directories are covered in [inputs and outputs](inputs-and-outputs.md).
