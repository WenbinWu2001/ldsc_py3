# Global configuration and paths

Last updated on: 2026-10-26

Choose SNP identity and genome build consistently before materializing inputs. The general default is `chr_pos_allele_aware`; the guided tutorial explicitly selects base `chr_pos`. Coordinate-based munging requires an explicit output build. Direct annotation projection also needs coordinates in the intended build. Indexed LD scoring inherits its stored configuration and rejects live overrides. See [configuration design](../../current/config-design.md) and [identity/build defaults](../../current/snp-identifier-genome-build-defaults.md).

## Inspect an artifact's genome build

Before choosing workflow inputs, inspect their recorded or inferable build:

```bash
ldsc infer-genome-build reference.chr22.bim
ldsc infer-genome-build raw-gwas.tsv.gz
ldsc infer-genome-build results/ldscore
ldsc infer-genome-build panels/reference
```

Each invocation takes one exact file or recognized LDSC artifact directory; `--help` is the only option.
Standalone BIM needs no BED/FAM companions.
Headered text (including gzip), Parquet, and recognized LDSC metadata JSON are supported, as are annotation, sumstats, LD-score, R²-panel, gene-index, and saved-result directory layouts.
The helper inspects only the artifact's own metadata and relevant members, without following upstream paths or searching unrelated neighboring files.
It does not expand globs or chromosome suites, accept build/column overrides, or create files, directories, logs, or failure markers.

```text
Genome build: hg19 (inferred from SNP coordinates; 1-based).
```

Metadata and build-specific coordinate columns produce a declared result without independent coordinate verification.
Otherwise bounded HM3 inference uses at most 5,000 usable rows per file/chromosome, at least 200 informative matches, and at least 99% support.
The reported build belongs to the inspected artifact's output coordinates, even when an upstream source used another build.
No coordinates are converted.

Insufficient or inconclusive evidence produces `Genome build: unknown (...)` with exit 0, including rsID-only tables and ordinary BED intervals without build metadata.
An unresolved member prevents a single inferred build; confidently conflicting evidence is an error.
Deliberately available hg19/hg38 alternatives, such as an R² panel with both build directories, are reported separately in that order and succeed.

Inspection uses available build evidence regardless of SNP matching mode.
When applicable metadata records `rsid` or `rsid_allele_aware`, either a known or unknown result also includes:

```text
This artifact’s workflow uses rsIDs for SNP matching.
```

An rsID column alone does not trigger this sentence.
Reports go to stdout; actionable errors go to stderr.
Missing/unreadable inputs, malformed/unsupported formats, unrecognized directories, and conflicts exit 1; invalid syntax and unexpected internal failures exit 2.
A successful inspection does not replace the workflow consistency checks below.
See the [full inspection contract](../../current/snp-identifier-genome-build-defaults.md#standalone-build-inspection).

Source: `inspect_genome_build()` in [_genome_build_inspection.py](../../../src/ldsc/_genome_build_inspection.py), CLI dispatch in [cli.py](../../../src/ldsc/cli.py), and [inspection acceptance tests](../../../tests/test_genome_build_inspection.py).

## Genome-build inference, checking, and overrides

Use the existing build flag for the coordinate stage you are running.
`--genome-build auto` checks all relevant inputs and selected shards; `build-r2-panel` and `munge-sumstats` use `--source-genome-build auto`.
Valid artifact metadata and build-specific position columns establish builds without coordinate re-inference.
Otherwise LDSC samples up to 5,000 rows per input/chromosome and requires at least 200 informative HM3 matches with at least 99% support.
Every relevant input participates before a consistency error is raised, regardless of input order.
This is bounded inference, not an exhaustive check for within-file mixtures.

A concrete `--genome-build hg19` or `hg38` is authoritative, including over conflicting metadata.
Use `--source-genome-build` for source-build workflows.
An override changes interpretation and does not lift coordinates over.
Munging still requires a separate `--output-genome-build` for coordinate-family output.
Saved output build metadata describes the coordinates after any requested conversion.

| SNP keep-list used for filtering | Rule |
| --- | --- |
| rsID keys | No coordinate-build validation. |
| Build-specific positions | Use the column matching the workflow build, ahead of generic `POS`; no coordinate inference is needed. |
| Generic `CHR/POS` without a declaration | Infer in auto mode; insufficient evidence is an error even if the panel is resolved. |

A list containing both `hg19_POS` and `hg38_POS` supplies alternatives; unused columns do not conflict.
If the matching named column is absent, generic `POS` is the fallback and must supply enough evidence under `auto`.
For a small generic list, check its build and explicitly pass the command's concrete build option.
Ordinary BED is unchecked and inherits the resolved workflow build; its compatibility is assumed.

On failure, the terminal error and workflow log group every checked input by `hg19`, `hg38`, or `unresolved`, and list BED/rsID-only inputs as `unchecked`.
Each entry identifies its flag, file or shard, and evidence or failure reason.
Check those inputs and use the suggested existing override only when you know the intended build.
Commands without a build override, such as regression and quantile analysis, require matching upstream artifacts.
See the [full policy and command scopes](../../current/snp-identifier-genome-build-defaults.md#genome-build-policy) and [grouped error examples](../../troubleshooting.md#common-genome-build-consistency-check-failed).

## Input paths

Quote patterns so the shell passes them intact. Pattern support depends on the command and option, not just the `-file` or `-sources` suffix.

| Input | Examples and rules |
| --- | --- |
| `infer-genome-build PATH` | One exact file or recognized artifact directory; no glob or `@` expansion. |
| Baseline/prebuilt query annotation sources | `"baseline/baseline.@.annot.gz"` declares a complete autosomal suite in direct query LD scoring; `"baseline/*.annot.gz"` selects actual matches. Validated contents determine chromosome scope. |
| BED and gene-list sources | `"beds/*.bed"` or `"pathways/*.txt"`; `@` is not expanded. |
| Scalar sumstats | A literal file or a quoted glob resolving to exactly one file. |
| `rg --sumstats-sources` | Exact files or globs selecting at least two inputs; `@` is not expanded. |
| Genetic maps for `ldscore`/`build-gene-ldscore-index` | Comma-separated exact file paths; neither `*` nor `@` is expanded. `build-r2-panel` has its own map resolver that supports both patterns. |
| Gene-coordinate catalog and control-gene file | One literal file per option; no glob expansion. |
| Input and output directories | Literal directory paths; no glob expansion. |

Consult the [complete command-specific table](../../current/path-specification.md#pattern-support-in-command-help) for restriction files, PLINK prefixes, reference metadata, and other exceptions. An exact-one glob that matches multiple files fails; an ordinary multi-file glob does not establish that all required chromosomes are present.

## Outputs and diagnostics

Core materializing commands require `--output-dir`. Existing command-owned artifacts require `--overwrite`; unrelated files are preserved. `munge-sumstats --infer-only` still requires the argument but writes nothing. Derived `plot` writes under `<result-dir>/plots/`, and `convert-h2-scale` writes under `<h2-result-dir>/postprocessing/liability-scale/`, without CLI output overrides.

Use the named scientific artifacts as downstream inputs. LD-score and gene-index roots contain required scientific `metadata.json`; most workflows also write logs and provenance under `diagnostics/`. Munging stores its configuration in Parquet metadata and has no JSON metadata sidecar. Some downstream workflows consume specific diagnostics: plotting reads declared sources, and `quantile-h2` uses saved complete-model coefficients and metadata. Keep complete result directories together. See [workflow logging](../../current/workflow-logging.md) and the [artifact metadata inventory](../../current/artifact-metadata-field-inventory.md).

Sources: `resolve_scalar_path` and `resolve_file_group` in [path_resolution.py](../../../src/ldsc/path_resolution.py), command parsers in [ldscore_calculator.py](../../../src/ldsc/ldscore_calculator.py) and [sumstats_munger.py](../../../src/ldsc/sumstats_munger.py), and the shared writers in [outputs.py](../../../src/ldsc/outputs.py).
