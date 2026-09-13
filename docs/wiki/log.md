# Wiki Update Log

Last updated on: 2026-09-13

## 2026-09-13 update | wiki consistency corrections

Reason: fix the guided tutorial's MDD-to-ADHD variable leak; correct identity defaults, baseline requirements, path expansion, private annotation disk use, index destinations, genetic-coordinate handling, legacy compatibility, and removed-flag wording. Replace branch-specific wiki URLs with relative links, replace unsupported resource estimates with benchmark context, complete empty navigation pages, and add baseline-only and saved-h2 conversion examples. Existing gene-list edits were preserved.

Validation: all 36 complete Bash-fenced LDSC examples parse with the current CLI; all local wiki file links and heading anchors resolve; all pages contain content and update dates. The actual guided-tutorial Bash blocks, with scientific commands replaced by an argument recorder, select MDD after the three-trait preparation sequence and when Analysis 3 starts with stale ADHD variables. Both index variants resolve matching build/monitor/Python paths. The Parquet preview reads 1,000 rows from one selected column in a synthetic file. `git diff --check` passed. No production datasets or HPC jobs were used.

Source: [audit and remediation record](../audits/2026-09-13-wiki-consistency.md), current CLI parsers, `build_query_shards`, `_resolve_genetic_map`, `normalize_annotation_chunk`, shared output writers, and the current path and result-schema documentation. The separate quantitative-count limitation remains open and is now linked from the continuous-annotation page.

## 2026-09-11 update | build-r2-panel, query-r2, examine-parquet-and-npz

Reason: replace the build/query placeholders with current usage and artifact contracts, rename the build page, explain SIGN_R versus sign_r and required sidecars, and correct inspection examples for quantized R².

Source: `src/ldsc/ref_panel_builder.py`, `src/ldsc/_kernel/ref_panel_builder.py`, `src/ldsc/r2_query.py`, `src/ldsc/outputs.py`, and the current Parquet and query format documentation.

## 2026-09-11 update | munge-sumstats, guided-tutorial

Reason: remove obsolete HM3 switches, explain default restriction and automatic mapping with explicit chain precedence, document stdout/log count summaries, require an output directory for inference previews, and use trait-named data files in downstream commands.

Source: `src/ldsc/sumstats_munger.py` (`build_parser`, `run_munge_sumstats_from_args`, `_render_munge_summary`), `src/ldsc/config.py` (`MungeConfig`), and [the current munging guide](../current/munge-sumstats.md#liftover-rules).
