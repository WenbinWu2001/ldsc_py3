# Approved Flag Cleanup Follow-up

Last updated on: 2026-09-11

Current names and legacy behavior are recorded in [`docs/current/legacy-cli-flag-map.md`](../../current/legacy-cli-flag-map.md). The original audit, inventory, and probe results in this directory remain historical evidence of the pre-cleanup interface.

## Implemented decisions

- Removed the hidden `--exclude-regions` alias from `ldscore` and `build-gene-ldscore-index`; retained `--regr-snps-exclude-regions` and its defaults.
- Retained `--N`, `--N-cas`, `--N-con`, and `--chunksize`. Corrected sample-size help to match LDSC2: input columns precede constants, omitted/zero `--n-min` uses N90/1.5, and constant N bypasses that filter.
- Removed `--keep-maf` and its configuration fields. Recognized frequency is always preserved as FRQ without MAF conversion and is omitted when absent. The MAF filter is unchanged.
- Preserved frequency through Parquet and gzip TSV reloads. The new test exposed an existing TSV blank-field problem; explicit `NA` serialization prevents missing coordinates from shifting Z/N/FRQ. FRQ export also avoids three-decimal rounding.
- Renamed `--format` to `--input-format`; retained `daner-new` and removed its separate preparation branch. Automatic and explicit new-DANER runs now share aliases, optional frequency handling, validation, and reading.
- Retained `--threads` for both `ldscore` and `build-gene-ldscore-index`, following the user's subsequent decision to revert the proposed `--workers` rename. `LDScoreConfig.threads`, the `run_ldscore(..., threads=...)` keyword, and `GeneLDScoreIndexBuildConfig.threads` remain. Relevant user docs explain that direct LD-score computation uses chromosome worker processes and gene-index construction uses a chromosome thread pool.
- Updated README, current docs, wiki guides, tutorials, relevant docstrings/comments, tests, and the annotation-memory benchmark. Historical audit/plan/spec/archive references and deliberate removed-name documentation remain unchanged in meaning.

## Verification

Behavior tests were run failing before their corresponding implementation changes, then rerun successfully. The new tests cover automatic frequency retention/absence, unchanged MAF filtering, Parquet and TSV frequency precision/alignment, DANER case aliases and optional fields, shared missing-N validation, removed flag rejection, retained sample-size fallback/threshold behavior, and concurrency option handling.

- Full pytest run before restoring `--threads`: 1520 passed, 1 skipped, 189 warnings, 132 subtests passed.
- Standard-library compatibility run, sequentially after pytest: 1000 tests, OK, 1 skipped.
- All 13 real subcommand help invocations exited successfully.
- Changed Python sources compiled; Markdown dates, new mapping links, notebook JSON, and `git diff --check` passed.
- Current documentation and examples use `--threads` for both commands. The benchmark retains its independent orchestration `--workers` option and forwards its value to `ldscore --threads`.

After restoring `--threads` on 2026-09-11, the combined LD-score parallelism, streaming, and gene-index suites passed: **113 tests**, with 68 existing rsid/genome-build warnings. Both real command help pages expose `--threads` and omit `--workers`; changed Python sources compile and `git diff --check` passes. This focused run exposed an existing test-order issue: initializer tests left the package logger at WARNING. Those tests now scope their logger changes with `caplog.at_level`, restoring the prior level before subsequent tests.

Other audit recommendations were not implemented without a corresponding decision. In particular, increasing-allele naming, plural INFO-list support, and general infer-only validation remain outside this approved cleanup. Unrelated edits appeared in the shared workspace during verification and were not modified or committed by this work. Test results describe their respective run snapshots.
