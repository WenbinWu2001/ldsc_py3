# Output-directory retry audit

Last updated on: 2026-09-14

## Scope and conclusion

Audited the local `restructure` branch starting at `214e891`, following the marker-only gene-index retry repair. The scope is all 13 materializing CLI commands, their workflow entry points, seven shared directory writers, and the output helpers and private writers they call. No HPC files or existing analysis outputs were changed.

The shared artifact-family checks already ignored root failure markers. The whole-directory gene-index validator was the exception. This audit additionally reproduced an interrupted-publication retry blocker, incomplete recognition of diagnostic-directory contents, and inconsistent destination handling for markers, logs, audits, and scratch. The fixes reuse existing path normalization and ownership checks; they do not change numerical computations, scientific schemas, or chromosome coverage requirements.

## Findings and corrections

| Finding | Reproduction and correction |
| --- | --- |
| Marker scope used an unexpanded path | With an output token such as `$LDSC_RETRY_OUTPUT/result`, scientific files went to the expanded destination but `_logging.overwrite_failure_marker` used a literal directory. The marker helper now calls `normalize_path_token` for both failure publication and success cleanup. |
| Repeated CLI options selected different destinations | argparse used the last destination, while `cli._raw_option_value` selected the first for the outer failure guard. Raw marker-scope extraction now uses the last value for both `--option value` and `--option=value`. Tests exercise core and derived commands through actual CLI failures. |
| Some logs, audits, and scratch bypassed normalization | Indexed LD-score logging, direct indexed-query scratch, direct partitioned-regression scratch, and h2/rg legacy-input audits used raw output arguments. They now use normalized destinations. Plotting and liability conversion normalize their source directory and Python output override consistently with preflight and markers. Regression and indexed-query tests check that no literal variable-named directory remains. |
| Diagnostic recognition ignored non-file entries | The gene-index check filtered `rglob` results to files, overlooking empty unrelated directories and dangling/directory symlinks. It now examines every entry, permits only owned regular files and their containing directories, and rejects symlinks and special files. Direct publication delegates to the same preflight used by the build workflow. |
| Existing interrupted-publication recovery rejected a failure-only destination | After a valid index moved to its owned transaction backup, a failure marker could recreate the destination. Recovery then demanded a publication marker in that diagnostic-only directory. Recovery now recognizes the same owned diagnostic contents as preflight, still requires exactly one validated matching backup, and carries the failure marker into the recovered index until the current build succeeds. Unrecognized destination contents still stop recovery and preserve the backup. |

Implementation: [`overwrite_failure_marker`](../../src/ldsc/_logging.py), [`_raw_option_value`](../../src/ldsc/cli.py), [`_preflight_gene_index_output`, `_is_diagnostics_only_gene_index`, `_recover_gene_index_publication`, and `run_indexed_ldscore`](../../src/ldsc/gene_ldscore_index.py), [`_run_explicit_indexed_ldscore`](../../src/ldsc/ldscore_calculator.py), [`RegressionRunner.estimate_partitioned_h2_batch`, `run_h2_from_args`, and `run_rg_from_args`](../../src/ldsc/regression_runner.py), [`plot_result`](../../src/ldsc/plotting/__init__.py), and [`convert_h2_scale`](../../src/ldsc/h2_scale.py).

## Workflow coverage

| Command | Output check and ownership | Retry verification |
| --- | --- | --- |
| `annotate` | `annotate_workflow.run_annotate` and `AnnotationDirectoryWriter`; annotation shards and named diagnostics | Marker-only destination, saved query reads, and scratch release through the real workflow |
| `ldscore` | `LDScoreDirectoryWriter`, direct preflight, indexed preflight, and `_ldscore_batch_output.write_ldscore_batches` | Direct PLINK forms and indexed CLI/API; marker removal, normalized logs/scratch, and numerical scores |
| `build-r2-panel` | `ReferencePanelBuilder._run`, `_ref_panel_output_family`; full-suite or chromosome-specific artifacts | Actual missing-BED failure followed by a successful retry for both scopes; applicable marker removed, other marker retained |
| `build-gene-ldscore-index` | `_preflight_gene_index_output`, `publish_gene_ldscore_index`, and existing transaction recovery | Repeated failure then success, marker-only publication, unrelated-content rejection, owned-backup recovery, and backup preservation on rejection |
| `convert-ldsc2-ldscores` | `LegacyLDScoreConverter.convert`; LD-score family plus conversion audit/log | Marker-only destination converted to a reloadable 22-chromosome suite |
| `munge-sumstats` | `SumstatsMunger.run`/`write`; named Parquet/text family plus workflow diagnostics | Marker-only destination through real compressed/plain input munging; row accounting and reload assertions retained |
| `h2` | `_preflight_regression_outputs` and `H2DirectoryWriter`; summary, bins, metadata, audit, log, and reserved derived roots | Marker-only destination and expanded-path artifact/audit publication; estimator isolated with the existing kernel stub |
| `partitioned-h2` | Shared regression preflight and `PartitionedH2DirectoryWriter`; full models, optional per-query trees, derived plots | Real CLI fits, sorted output trees, marker removal, overwrite cleanup, and direct API scratch normalization |
| `quantile-h2` | `run_quantile_h2_from_args` and `QuantileH2DirectoryWriter`; summaries, coefficients, alignment audit, metadata/log, derived plots | Marker-only destination through the real post-fit workflow; hand-worked heritability assertions retained |
| `rg` | Shared regression preflight and `RgDirectoryWriter`; summaries, optional pair tree, legacy audit, metadata/log, derived plots | Marker-only destination and expanded-path publication using the real writer; pair estimation isolated by the existing stub |
| `query-r2` | `_write_query_r2_directory` and `QueryR2DirectoryWriter` | Real CLI pair query from a marker-only destination; R2 value, schema, and metadata assertions retained |
| `convert-h2-scale` | Explicit table/optional plot/metadata/log family under `postprocessing/liability-scale` or Python override | Missing-metadata failure then successful retry; normalized source/default/custom destination; source metadata and unrelated file preserved |
| `plot` | Explicit plot/metadata/log family under `plots` or Python override | Missing-metadata failure then successful retry; normalized source/default/custom destination; source metadata and unrelated file preserved |

The common writer-family test independently checks all seven declarations with overwrite enabled and disabled: root markers and unrelated notes do not collide, owned files still require overwrite, stale cleanup preserves unowned files, and markers never enter declared metadata. Workflow-owned diagnostics from a failed attempt may still require `--overwrite`; this is an intentional collision, unlike rejection caused solely by a marker.

## Supporting output paths

Inspected directory enumeration (`iterdir`, `rglob`, `glob`, `listdir`, and `walk`), directory creation/removal, file writes, shared output-preflight callers, and every failure-marker scope across `src/ldsc`. Input-only scans in legacy conversion, PLINK resolution, and reference-panel loading are not output-directory validators.

[`path_resolution`](../../src/ldsc/path_resolution.py) creates/reuses directories and checks explicitly enumerated output families. [`outputs`](../../src/ldsc/outputs.py) adds only declared scientific/diagnostic paths and narrowly matched annotation/query/drop-report siblings; its globs do not match `RUN_FAILED` files. Reserved derived roots and optional per-query/per-pair trees retain their existing whole-tree ownership and cleanup contracts. [`_result_files`](../../src/ldsc/_result_files.py) writes an explicit destination via a temporary sibling and separately validates declared derived-workflow inputs.

Private annotation, gene-query, quantile, LD-score batch, and regression workspaces write under the resolved workflow root or an explicitly supplied owned temporary directory. Their leaf writers receive resolved paths and do not require a globally empty result directory. Direct low-level writers do not declare a whole workflow successful and therefore leave workflow failure-marker lifecycle management to their entry-point guard. Read-only CLI help and munging inference retain their no-write behavior, including preservation of existing markers.

## Verification

- Before repair, focused reproductions failed for each corrected behavior; shared writer-family marker checks already passed.
- Final focused workflow set: **252 passed**, including 3 subtests. This set covers the repaired branches; the full suite additionally covers the other seeded workflow tests and existing overwrite/collision protections.
- Full pytest suite: **1,729 passed, 1 skipped** in 119.51 seconds. The existing dependency-absence test is skipped because PyArrow is installed.
- Sequential unittest compatibility run: **1,000 tests, OK, 1 skipped**.
- Both `ldsc --help` and `python -m ldsc --help` passed; `git diff --check` passed and all 15 audit links resolve locally.
- The audit uses local deterministic fixtures. It does not simulate every filesystem failure or rerun the user's HPC analysis.

Primary regression evidence: [`test_output_directory_retries.py`](../../tests/test_output_directory_retries.py), [`test_gene_ldscore_index.py`](../../tests/test_gene_ldscore_index.py), [`test_plink_workflow_resolution.py`](../../tests/test_plink_workflow_resolution.py), [`test_regression_workflow.py`](../../tests/test_regression_workflow.py), and [`test_regression_streaming.py`](../../tests/test_regression_streaming.py). Existing successful workflow tests for annotation, munging, conversion, quantiles, and R2 queries now also start with an old failure marker and verify its removal.
