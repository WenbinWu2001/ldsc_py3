# Global configuration and paths

Last updated on: 2026-10-20

Choose SNP identity and genome build consistently before materializing inputs. The general default is `chr_pos_allele_aware`; the guided tutorial explicitly selects base `chr_pos`. Coordinate-based munging requires an explicit output build. Direct annotation projection also needs coordinates in the intended build. Indexed LD scoring inherits its stored configuration and rejects live overrides. See [configuration design](../../current/config-design.md) and [identity/build defaults](../../current/snp-identifier-genome-build-defaults.md).

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
