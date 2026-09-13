# Wiki consistency audit

Last updated on: 2026-09-13

## Remediation status

The user authorized corrections after the initial audit. Findings 1–9 below have been addressed in the wiki: the tutorial keeps the focal trait explicit, stale contracts and destinations are corrected, empty pages now provide workflow or resource guidance, and performance guidance refers to measured evidence. The optional Parquet preview improvement, removed-flag wording, and link to the existing quantitative-count limitation were also added. Package implementation was unchanged; the separately recorded quantitative-count issue remains open.

Post-edit verification: all **36 complete Bash-fenced LDSC commands** parse with the current CLI, all local file links and heading anchors resolve, and all 21 wiki pages have content and update dates. Executing the actual guided-tutorial Bash blocks with scientific commands replaced by an argument recorder confirms that MDD, SCZ, and ADHD are prepared and the final partitioned analysis selects MDD. Analysis 3 also selects MDD when started with stale ADHD variables. Both index build variants use matching monitoring and exported Python input paths. The Parquet preview was executed against a synthetic file and read only 1,000 rows of one selected column. `git diff --check` passed. No production analysis was run; the full suite was not repeated for these documentation changes.

The findings, line references, page dispositions, and test results below record the **pre-fix audit snapshot**, not outstanding work. See the [wiki update log](../wiki/log.md) for the correction record.

## Initial audit scope and result

Reviewed all 21 Markdown files under `docs/wiki/` against the local `restructure` checkout at `55e127b`, current workflow documentation, CLI parsers, and relevant tests. The working copies of `main-functionalities/ldscore-from-gene-list.md` and `main-functionalities/partitioned-h2.md` already contained user edits; those contents were included in the review and preserved. No wiki or package source files were changed.

The wiki needs updates. The most consequential issue is a guided-tutorial variable leak that runs the final partitioned analysis on ADHD instead of the advertised MDD trait. Other confirmed discrepancies concern identity defaults, required baselines, path expansion, private disk use, index monitoring paths, genetic-map handling, and legacy compatibility coverage. Several pages remain empty or skeletal. Recent HM3, R² naming, gene-list catalog, and rg plotting updates are substantially reflected in the reviewed pages.

Priorities below distinguish an example that silently selects a different analysis (P1), misleading operational or contract guidance (P2), and completeness/editorial work (P3).

## Findings and recommended corrections

### 1. P1 — The guided tutorial changes the focal trait before Analysis 3

**Location:** [guided tutorial](../wiki/guided-tutorial.md), “Munge two additional sumstats,” lines 218–231, and Analysis 3 Step 2, lines 379–382.

The Analysis 2 loop assigns `TRAIT_NAME`, `RAW_SUMSTATS_FILE`, and `SUMSTATS_OUT_DIR` for each additional trait. After it finishes, the trait is `adhd2019`. Analysis 3 reuses those variables without restoring `mdd2025`, although the motivating example and its output tree describe MDD. The command succeeds on the wrong intended trait when the instructions are executed sequentially.

**Evidence:** A Bash reproduction using only the documented variable assignments resolved the final input to `/example/tutorial_output/sumstats_processed/adhd2019/adhd2019.parquet` and the output to `/example/tutorial_output/partitioned-h2/adhd2019`. [`run_partitioned_h2_from_args`](../../src/ldsc/regression_runner.py), lines 2769–2791, loads exactly the supplied `sumstats_file`.

**Correction:** Give the additional-trait loop its own variable names and explicitly set the focal trait and input path at the start of Analysis 3. Verify the tutorial sections together as well as individually.

### 2. P2 — The LDSC2 changes page states the wrong default and suggests removing required baselines

**Location:** [changes-from-LDSC2](../wiki/main-functionalities/changes-from-LDSC2.md), lines 13–15.

“Defaults to chr::pos” omits allele awareness and is not the public mode spelling. “deprecate the baseline annotations” contradicts the current direct query workflow, which requires explicit baseline annotations. If this note meant precomputed legacy baseline LD scores rather than annotations, it needs to say so; explicit legacy conversion also remains supported.

**Evidence:** [`GlobalConfig`](../../src/ldsc/config.py), lines 176–182, and [`ldscore_calculator.build_parser`](../../src/ldsc/ldscore_calculator.py), line 993, default to `chr_pos_allele_aware`. The same source, lines 89–95 and 945, explicitly requires baseline annotations for direct queries. The index supplies its own stored baseline in indexed mode. See [legacy conversion](../current/legacy-ldscore-conversion.md) for the supported alternative.

**Correction:** Replace the draft notes with a current contract table, including command-specific identity exceptions and the distinction between annotations and precomputed LD scores.

### 3. P2 — The general path rule promises expansion that some options do not support

**Location:** [guided tutorial](../wiki/guided-tutorial.md), Remarks item 2, line 49; [global-config](../wiki/main-functionalities/global-config.md), lines 13–15.

The tutorial says all `-sources` options support both `*` and `@`. In fact, BED/gene-list sources and rg sumstats sources do not expand chromosome placeholders. Genetic-map sources for `ldscore` and `build-gene-ldscore-index` accept comma-separated exact paths and expand neither wildcard nor chromosome tokens. The global-config draft's suffix-only shorthand is similarly insufficient: scalar inputs can sometimes use exact-one globs, while the control-gene file is literal.

**Evidence:** [Path specification, “Command-specific path support” table](../current/path-specification.md), lines 57–67; [`ldscore_calculator.build_parser`](../../src/ldsc/ldscore_calculator.py), lines 1153–1168; [`_expand_focal_gene_list_sources`](../../src/ldsc/gene_list_resolver.py).

**Correction:** Link the command-specific table and explain that `@` requests complete autosomal suites only for supported inputs. Avoid using flag suffixes as a universal parser contract.

### 4. P2 — The gene-list tutorial incorrectly rules out disk-backed annotation preparation

**Location:** [LDSC-SEG-PC-genes](../wiki/LDSC-SEG-PC-genes.md), direct-mode description, line 35.

The page says query annotations are constructed in memory and “not written to disk.” Direct mode now writes private chromosome query matrices as part of bounded annotation preparation. It still does not publish reusable `.annot.gz` query files by default, but temporary disk space and I/O are real requirements.

**Evidence:** [`build_query_shards`](../../src/ldsc/_annotation_queries.py), lines 94–115, creates `query-<chrom>.npy` stores and writes projected values. [`ColumnStore.create`](../../src/ldsc/_annotation_storage.py), lines 99–114, creates the backing file. See [annotation memory design, “Workspace ownership and lifetime”](../current/annotation-memory-design.md).

**Correction:** Explain private on-disk staging, cleanup, and the distinction from published annotation outputs. Add the current query-batch and chromosome-worker controls or link their existing wiki descriptions.

### 5. P2 — Index monitoring and Python validation use a different destination from the build example

**Location:** [build-gene-ldscore-index](../wiki/utility-functionalities/build-gene-ldscore-index.md), build command line 140, monitoring lines 253–257, Python example line 266.

The primary build writes `${INDEX_ROOT}/baseline_100kb`, but the monitoring command tails `.production.build-state`, the success-log path points into `production`, and Python loads `production`. Neither build example creates that index. Copying these follow-up commands fails to monitor or load the index just built.

**Evidence:** The same page's “Output and identity” tree correctly defines the hidden state directory as `.<index-name>.build-state`. [`gene_ldscore_index.py`](../../src/ldsc/gene_ldscore_index.py) derives build-state and publication paths from the actual output target.

**Correction:** Define one `INDEX_DIR` and reuse it in building, monitoring, validation, and indexed LD scoring. For the existing primary example, use `.baseline_100kb.build-state` and `baseline_100kb`.

### 6. P2 — Legacy interoperability is still labeled untested

**Location:** [guided tutorial](../wiki/guided-tutorial.md), “Backward compatibility with the legacy ldsc python2 codebase,” lines 469–473.

The “not tested” note and TODO contradict the implemented compatibility boundary and dedicated tests. Legacy munged summary statistics are accepted directly by regression; legacy LD-score suites require explicit conversion into the canonical directory. The beginning of the same tutorial already describes that conversion route correctly.

**Evidence:** [Legacy sumstats compatibility](../current/legacy-sumstats-compatibility.md), [legacy LD-score conversion](../current/legacy-ldscore-conversion.md), [`tests/test_legacy_sumstats_compatibility.py`](../../tests/test_legacy_sumstats_compatibility.py), and [`tests/test_legacy_ldscore_converter.py`](../../tests/test_legacy_ldscore_converter.py). Both test modules passed during this audit.

**Correction:** Replace the note with those two supported routes and their restrictions. Describe actual tested coverage without implying universal scientific equivalence across all inputs.

### 7. P2 — Genetic-coordinate guidance needs the PLINK/R² distinction

**Location:** [guided tutorial](../wiki/guided-tutorial.md), cM caveat at line 359 and annotation remarks at line 452; [ldscore](../wiki/main-functionalities/ldscore.md), line 9.

The tutorial's unconditional error claim for all-zero BIM CM overlooks the supported explicit genetic-map route. The ldscore page describes map-or-BIM behavior without limiting it to PLINK; R² input instead uses its required metadata sidecar. Finally, “CM column is preserved” can be read as preserving input values, but generated query annotations contain an `NA` placeholder, as the newer annotate page correctly explains.

**Evidence:** [`_kernel/ref_panel.py`](../../src/ldsc/_kernel/ref_panel.py), lines 465–470 and `_resolve_genetic_map` at lines 1095–1113, applies PLINK map interpolation and uses authoritative sidecar CM for R². [`_annotation_parsing.py`](../../src/ldsc/_annotation_parsing.py), lines 52–55, discards input CM values and sets `NaN`.

**Correction:** Explain PLINK's informative-BIM-or-matching-map requirement, R² sidecar authority, and the generated annotation placeholder separately.

### 8. P3 — The main workflow navigation is incomplete

**Locations:** [h2](../wiki/main-functionalities/h2.md) and [how-to-customize-your-bed-files](../wiki/utility-functionalities/how-to-customize-your-bed-files.md) are empty; [resources](../wiki/utility-functionalities/resources.md) is whitespace-only. [global-config](../wiki/main-functionalities/global-config.md) and [changes-from-LDSC2](../wiki/main-functionalities/changes-from-LDSC2.md) are unfinished notes. [ldscore](../wiki/main-functionalities/ldscore.md) contains only CM and memory subsections despite being a main entry point.

There is also no wiki introduction or navigation entry for `convert-h2-scale`, which is exposed by the current CLI. Plotting and quantile-h2 do have coverage through other wiki pages, so their lack of a same-named page is not itself an omission. The [partitioned-h2 page](../wiki/main-functionalities/partitioned-h2.md) could use a concise baseline-only example and output interpretation; its detailed query-batch coverage does not fill that introductory role.

**Evidence:** Full file inventory below; [`cli.build_parser`](../../src/ldsc/cli.py), lines 76–136, exposes all 13 commands. [`h2_scale.py`](../../src/ldsc/h2_scale.py) implements saved-result conversion and prevalence sensitivity.

**Correction:** Complete or explicitly redirect the empty/stub pages to maintained guides. Add a short h2-scale route to the main tutorial. Use the existing [writing guide](../wiki/main-functionalities/writing-guide.md) for inputs, minimal examples, outputs, and interpretation. Six existing pages lack an update date: the three empty/blank pages and `changes-from-LDSC2`, `global-config`, and `writing-guide`. Add dates when editing them.

### 9. P3 — Release links and performance wording need maintenance

**Location:** [guided tutorial](../wiki/guided-tutorial.md), lines 37, 55, 313, 353–354 and final TODOs; [LDSC-SEG-PC-genes](../wiki/LDSC-SEG-PC-genes.md), lines 37–38 and 88–89.

The guide asks readers to install `main` but links two local documentation topics through hard-coded `ldsc3-beta` URLs. Those links are not version-consistent with this checkout. Their remote availability was not checked. Prefer relative links for the two wiki topics and specify the intended release/install target deliberately.

The statement that globbing BED files makes their LD scores compute “in parallel” needs to distinguish query batching from chromosome workers, whose default is one. The claim that the codebase has not been benchmarked is also too broad: local synthetic benchmark evidence now exists, although it does not validate the production GB/hour estimates shown in these tutorials.

**Evidence:** [`LDScoreConfig`](../../src/ldsc/config.py), lines 532–535, and [`_resolve_worker_count`](../../src/ldsc/ldscore_calculator.py); [local resource results, “Matched command measurements” and “Inputs, environment, and measurement method”](annotation-memory/results.md).

**Correction:** Link the benchmark evidence with its workload and limitations; qualify or remove unattributed production estimates. Explain that grouping sources shares computation, while `--threads` controls chromosome concurrency.

## Additional coverage worth adding

- **Known quantitative-annotation limitation:** The [memory completion review, R1](annotation-memory/completion-review.md) and [active memory plan, Status](../plans/2026-09-10-annotation-workflow-memory.md) still record batch-dependent float32 annotation-count discrepancies. The [continuous-annotation wiki](../wiki/continuous-annotation-partitioned-ldsc.md) would benefit from linking that limitation until it is resolved. This is an existing recorded issue, not a new failure reproduced in this audit; the passing focused tests below do not close it.
- **Large Parquet inspection:** [examine-parquet-and-npz](../wiki/utility-functionalities/examine-parquet-and-npz.md), line 11, reads the entire file before showing examples. That API usage is valid, but a row-group/column-first example would better suit large R² panels. This is a usability improvement, not a stale schema finding.
- **Removed flag wording:** The guided and LDSC-SEG tutorials call `--write-per-query-results` “deprecated.” It is absent from the current parser. Say “removed; per-query results are automatic” so readers do not infer that the old flag remains accepted.

## Page-by-page disposition

“No discrepancy found” means no concrete mismatch was identified in the reviewed claims; it is not a guarantee of exhaustive scientific validation.

| Wiki file | Disposition |
| --- | --- |
| `guided-tutorial.md` | Update: findings 1, 3, 6, 7, 9; complete remaining draft sections. |
| `LDSC-SEG-PC-genes.md` | Update private-disk description and resource estimates; clarify removed flag. |
| `continuous-annotation-partitioned-ldsc.md` | No new contract discrepancy found; link the recorded quantitative-count limitation. |
| `log.md` | Historical entries remain consistent with their described changes; not a complete maintenance history. |
| `main-functionalities/build-r2-panel.md` | No discrepancy found in current naming, schema, threshold, and output guidance. |
| `main-functionalities/changes-from-LDSC2.md` | Rewrite stale draft; finding 2. |
| `main-functionalities/global-config.md` | Complete draft and qualify path rules; findings 3 and 8. |
| `main-functionalities/h2.md` | Empty; add guide or redirect. |
| `main-functionalities/ldscore-from-gene-list.md` | No discrepancy found in current working copy's input, control, catalog, and batching guidance. |
| `main-functionalities/ldscore.md` | Expand stub and scope CM explanation to backend; findings 7 and 8. |
| `main-functionalities/munge-sumstats.md` | No discrepancy found in HM3, liftover, naming, or sample-size guidance. |
| `main-functionalities/partitioned-h2.md` | Current batch guidance is consistent; add introductory baseline-only example and interpretation. |
| `main-functionalities/rg.md` | No discrepancy found in current h2 annotation and plotting guidance. |
| `main-functionalities/writing-guide.md` | Editorial guidance, not obsolete API documentation; add update date when maintained. |
| `utility-functionalities/annotate.md` | No discrepancy found in the reviewed gene-route and output contracts. |
| `utility-functionalities/build-gene-ldscore-index.md` | Fix follow-up destination mismatch; finding 5. |
| `utility-functionalities/convert-ldsc2-ldscores.md` | No discrepancy found in supported filenames, input families, and flags. |
| `utility-functionalities/examine-parquet-and-npz.md` | Current schema examples; improve large-file read example. |
| `utility-functionalities/how-to-customize-your-bed-files.md` | Empty; add guide or redirect to maintained BED format documentation. |
| `utility-functionalities/query-r2.md` | No discrepancy found in current `SIGN_R`/`sign_r`, status, and panel guidance. |
| `utility-functionalities/resources.md` | Whitespace-only; add resource navigation or remove from intended navigation. |

## Validation and limits

- Read all 21 Markdown files, including both pre-existing modified working copies.
- Verified the imported package path is this checkout's `src/ldsc/__init__.py` using the existing `ldsc3-dev` Python environment.
- Parsed all **31 complete `ldsc` commands in Bash fences** with `ldsc.cli.build_parser()`: all accepted. Skipped two explicitly incomplete examples containing ellipses. Parsing verifies syntax and required arguments, not resource existence, shell state, or numerical validity.
- Checked relative Markdown file links and local heading anchors: no missing targets found. Remote URLs and external resource availability were not checked.
- Reproduced the guided tutorial's trait-variable leak in Bash without reading datasets or invoking scientific workflows.
- Ran `PYTHONPATH=src /Users/wenbinwu/miniforge3/envs/ldsc3-dev/bin/python -m pytest -q tests/test_cli_help.py tests/test_munge_hm3_defaults.py tests/test_legacy_sumstats_compatibility.py tests/test_legacy_ldscore_converter.py tests/test_annotation_storage.py tests/test_plotting.py`: **112 passed, 7 warnings in 33.29 seconds**. Warnings were pandas dtype-assignment FutureWarnings in `sumstats_munger.py:1800`.
- Did not run the full suite, the production-scale tutorial, external-data validation, or HPC commands. No package behavior was changed.
