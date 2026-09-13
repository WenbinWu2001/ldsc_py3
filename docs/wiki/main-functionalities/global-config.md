# Global configuration and paths

Last updated on: 2026-09-13

Choose SNP identity and genome build consistently before materializing inputs. The general default is `chr_pos_allele_aware`; the guided tutorial explicitly selects base `chr_pos`. Coordinate-based munging requires an explicit output build. Direct annotation projection also needs coordinates in the intended build. Indexed LD scoring inherits its stored configuration and rejects live overrides. See [configuration design](../../current/config-design.md) and [identity/build defaults](../../current/snp-identifier-genome-build-defaults.md).

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
