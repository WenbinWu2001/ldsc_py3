# Plan: multi-trait genetic-correlation refactor for `ldsc rg`

## Context

The current `ldsc rg` command accepts exactly two summary-statistics files via
`--sumstats-1-file` / `--sumstats-2-file` and emits a single-row `rg.tsv`.
Geneticists running rg on more than two traits today have to script a loop
themselves and concatenate per-pair outputs, which is error-prone and hides the
heritability and intercept context that helps interpret rg.

This refactor turns `ldsc rg` into one unified multi-trait command that takes
N≥2 sumstats files (paths and/or globs), computes either all pairs or
anchor-vs-rest, and writes or returns a long-form, interpretation-ready result
family. Pairwise rg is just the N=2 case of the same workflow. The numerical
kernel stays untouched; the change is at the workflow/UX layer.

Worktree:
`/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/.worktrees/refactor/rg/`

## Confirmed design decisions

1. **CLI surface (unified replacement, hard removal of old flags).**
   - New required flag: `--sumstats-sources <PATH|GLOB> [...]`. Multi-token,
     glob-aware, dedup-preserving-first-seen-order, must resolve to ≥2 distinct
     concrete paths.
   - New optional flag: `--anchor-trait-file <PATH>`. When set, the path must
     match exactly one resolved entry from `--sumstats-sources`; only
     anchor-vs-rest pairs (N−1 rows) are computed. When unset, all C(N,2) pairs.
   - New optional flag: `--write-per-pair-detail` (default off). Writes a
     diagnostic `pairs/{idx:04d}_{safe_trait_a}_vs_{safe_trait_b}/` subtree with a
     per-pair `rg_full.tsv` and `metadata.json`. This tree is not part of the
     main scientific result; the root `rg.tsv` and `rg_full.tsv` already
     contain every pair. This flag requires `--output-dir`; reject otherwise
     before loading inputs.
   - Removed outright: `--sumstats-1-file`, `--sumstats-2-file`,
     `--trait-name-1`, `--trait-name-2`. No deprecated stubs, no
     migration messages — the refactored version is unpublished, so
     unknown-flag rejection from argparse is sufficient.
   - Intercept flags: `--intercept-h2 <FLOAT>` is now a single scalar
     broadcast across all traits; `--intercept-gencov <FLOAT>` is a single
     scalar broadcast across all pairs; `--no-intercept` constrains h2 to 1
     and gencov to 0. `--no-intercept` is rejected when combined with
     `--intercept-h2` or `--intercept-gencov`; the current silent precedence
     behavior is intentionally not preserved. The legacy `nargs=2` form is
     removed outright.
   - Other CLI options (`--ldscore-dir`, `--output-dir`, `--overwrite`,
     `--count-kind`, `--n-blocks`, `--two-step-cutoff`, `--chisq-max`,
     `--log-level`) are unchanged. `--output-dir` remains optional, but its
     help text and `run_rg_from_args` docstring should strongly recommend it.
     Without an output directory, the Python workflow still returns the full
     result family, but CLI invocations only print the concise table and write
     no diagnostic files or log. `run_rg_from_args` itself must not print.

2. **Per-pair SNP support: per-pair intersection** (rigorous; matches existing
   pairwise rg). Each output row carries its own `n_snps_used`.

3. **Computational strategy: pure pairwise loop.**
   - `estimate_rg_pairs` builds one `RGRegressionDataset` per pair, then fits
     that dataset through the same private helper used by public
     `estimate_rg`. No factor-out of `build_rg_dataset`, no merge cache, and
     no alternate preprocessing path.
   - Per-pair preprocessing goes through `build_rg_dataset` exactly as
     committed in `f36982c` (`regression_runner.py:265-379`), which
     produces an `RGRegressionDataset` (`regression_runner.py:123-142`)
     on the per-pair trait_1 ∩ trait_2 ∩ LD-ref SNP set. Allele
     harmonization, zero-variance drop, and `Hsq` fits inside kernel
     `reg.RG` all happen on that same per-pair SNP set, so jackknife
     delete-blocks for `hsq1`, `hsq2`, and `gencov` line up by
     construction. **Zero kernel touch, statistically identical to
     today's single-pair `ldsc rg`.**
   - Per-trait h2 (for `h2_per_trait.tsv`) is computed independently via
     the existing `build_dataset(trait_i, ldscore_result)` +
     `estimate_h2(dataset)` path; that produces an `Hsq` fit on each
     trait's full LD-ref intersection (matches `ldsc h2`), separate from
     the per-pair `Hsq` inside `RG`. The two h2 numbers will differ
     because they are fit on different SNP sets — intentional and
     documented in the schema section below.

   **Compute strategy rationale.** A merge-cache design that hoisted the
   trait + LD-score merge out of the per-pair loop and stored it per
   trait was considered. Estimated savings for N<10 traits with one
   baseline LD column: ~10% wall-time reduction (~1–2 min on a ~15-min
   run for N=9 in chr_pos mode; less in rsid mode), at a cost of ~0.8–1.2
   GB additional cached state and ~3.5× the implementation surface
   (~230 vs ~85 LOC plus a permanent byte-identity test against
   `build_rg_dataset`). The savings are below the user-noticeable
   threshold at this scale; the maintenance liability is permanent.
   Rejected pending real-workload profiling at N≥20. If preprocessing
   ever becomes the bottleneck, parallelizing the per-pair loop across
   processes (independent pairs, shared read-only `LDScoreResult`) is a
   bigger and more direct win than caching, and is a more natural future
   addition.

4. **Output layout (canonical rg result directory).**
   When `--output-dir` is supplied:
   ```
   <output-dir>/
     rg.tsv                                # concise, publication-ready
     rg_full.tsv                           # comprehensive diagnostic
     h2_per_trait.tsv                      # one row per input trait
     rg.log                                # from workflow_logging
     pairs/                                # only when --write-per-pair-detail
       0001_<safe_trait_a>_vs_<safe_trait_b>/
         rg_full.tsv                       # one-row, comprehensive view
         metadata.json                     # per-pair reproducibility record
       0002_...
   ```
   This mirrors the established `partitioned_h2.tsv` (compact) /
   `partitioned_h2_full.tsv` (detailed) convention. The default `rg.tsv` is
   the file 90% of users will open and report from; `rg_full.tsv` is for
   diagnostic deep-dives.
   When `--output-dir` is omitted, no files are written. The Python workflow
   still constructs and returns the full `RgResultFamily`; CLI dispatch prints
   only the concise `rg.tsv`-schema table to stdout. Users who need persisted
   diagnostics should rerun with `--output-dir`.

5. **`rg.tsv` schema (concise, publication-ready, ALWAYS long-form):**
   8 columns: `trait_1`, `trait_2`, `n_snps_used`, `rg`, `rg_se`, `p`,
   `p_fdr_bh`, `note`. `p_fdr_bh` is Benjamini-Hochberg over finite
   p-values only. `note` is empty for successful pairs and gives failed
   pairs a short pointer such as
   `Failed; see rg_full.tsv error column; use --output-dir for rg.log`.
   This is the headline table — it matches what genetic-correlation papers
   typically report verbatim.

5b. **`rg_full.tsv` schema (comprehensive diagnostic, ALWAYS long-form):**
    32 columns in this order, grouped by purpose:
    - **Identifiers and sample size** (3): `trait_1`, `trait_2`, `n_snps_used`
    - **Headline rg result** (6): `rg`, `rg_se`, `z`, `p`, `p_fdr_bh`,
      `p_bonferroni`
    - **Per-trait heritability fit on the pair's SNP intersection** (4):
      `h2_1`, `h2_1_se`, `h2_2`, `h2_2_se`
    - **Genetic covariance** (2): `gencov`, `gencov_se`
    - **Intercepts (sample-overlap diagnostics)** (6): `intercept_h2_1`,
      `intercept_h2_1_se`, `intercept_h2_2`, `intercept_h2_2_se`,
      `intercept_gencov`, `intercept_gencov_se`
    - **Confounding-vs-polygenic ratio diagnostics** (4): `ratio_1`,
      `ratio_1_se`, `ratio_2`, `ratio_2_se` —
      `ratio_t = (intercept_h2_t − 1) / (mean_chisq_t − 1)`
    - **Inflation diagnostics** (4): `lambda_gc_1`, `lambda_gc_2`,
      `mean_chisq_1`, `mean_chisq_2`
    - **Provenance** (1): `pair_kind` ∈ {`all_pairs`, `anchor`}
    - **Failure status** (2): `status` ∈ {`ok`, `failed`}, `error`
    Notes:
    - `p_bonferroni`: `min(1.0, p * n_valid_pairs)` where
      `n_valid_pairs` = the number of finite p-values being adjusted.
    - `p_fdr_bh` is identical to the value in `rg.tsv` for the same row.
    - Benjamini-Hochberg is implemented locally with NumPy: sort finite
      p-values ascending, compute `p_sorted * m / rank`, apply reverse
      cumulative minima for monotonicity, clip to `[0, 1]`, and scatter back
      to original row order. NaN/non-finite p-values are ignored during
      adjustment and remain NaN in adjusted columns.
    - When a pair has zero SNP overlap or any ordinary exception occurs during
      that pair's build+fit+summary phase: the row is still emitted in both
      `rg.tsv` and `rg_full.tsv` with `n_snps_used=0` and NaN-filled stats;
      `status="failed"` and `error` contains `<ExceptionClass>: <message>`;
      pair iteration continues.
    - Do not catch `BaseException` subclasses such as `KeyboardInterrupt` and
      `SystemExit`.
    - Numeric rg estimates are preserved even outside `[-1, 1]` or
      `[-1.2, 1.2]`; non-numeric kernel `NA` values are converted to failed
      NaN rows with an explanatory `error`.
    - All TSV and stdout table writes use `na_rep="NaN"`.
    - `intercept_gencov_se`, `ratio_1`, `ratio_1_se`, `ratio_2`, `ratio_2_se`
      are sourced from kernel attributes already produced by `RG.gencov` and
      `RG.hsq1` / `RG.hsq2` — no kernel changes needed.

6. **`h2_per_trait.tsv` schema:** reuse `summarize_total_h2` columns —
   `trait_name`, `n_snps`, `total_h2`, `total_h2_se`, `intercept`,
   `intercept_se`, `mean_chisq`, `lambda_gc`, `ratio`, `ratio_se`. One row
   per resolved sumstats source.

7. **Self-pairs:** skip (no rg-vs-self diagonal rows).

8. **Trait-name collision policy:** `load_sumstats(path, trait_name=None)`
   currently defaults `trait_name` to `Path(path).name`. If the resolved
   sources produce duplicate basenames, disambiguate by appending
   `@<parent_dir_basename>`; if collisions persist, append a deterministic
   `@<8-char SHA1 prefix of full path>`. Emit one INFO log per
   disambiguation.

   **Trait-name UX TODO:** keep filename-derived names for v1. Future work
   should design propagation of scientific trait labels from sumstats munging
   metadata or a manifest-style input. Keep the flag name
   `--anchor-trait-file` for now; revisit it with the broader trait-label
   design.

9. **Edge-case handling:**
   - Glob resolves to <2 distinct files → workflow raises `ValueError` with a
     message pointing at `--sumstats-sources`.
   - `--anchor-trait-file` doesn't match any resolved source (compared by
     normalized absolute path AND by resolved trait_name) → hard error
     before any compute.
   - Anchor file matched both directly and by another glob → dedupe so the
     anchor never self-pairs.
   - Trait with zero SNP overlap to the LD-score reference: the per-trait
     h2 path will hard-fail at `build_dataset` (preserves current
     `estimate_h2` behavior). The run should hard-fail at trait load /
     per-trait h2 fitting rather than silently dropping the trait —
     traits that cannot be analyzed should not be quietly excluded from
     the run.
   - `--write-per-pair-detail` without `--output-dir` → reject before input
     loading. The flag requests filesystem artifacts, so ignoring it would
     mislead users.
   - Per-pair failures (including the three `ValueError`s raised by
     `build_rg_dataset` — "No overlapping SNPs remain after merging
     both sumstats tables" at `:333`, "No allele-compatible SNPs remain
     after harmonizing" at `:344`, "All LD-score columns have zero
     variance" at `:356` — plus ordinary numerical/kernel exceptions) →
     caught in `estimate_rg_pairs` around that pair's build+fit+summary
     block. Emit a NaN row with `n_snps_used=0`, `status="failed"`, an
     `error` message in `rg_full.tsv`, a concise `note` in `rg.tsv`, log a
     WARNING with traceback, and continue to later pairs. Do not catch
     `BaseException` subclasses such as `KeyboardInterrupt` and `SystemExit`.
   - Kernel non-numeric `NA` values for rg outputs → convert to a failed NaN
     row. Numeric out-of-range rg estimates are preserved.
   - `--no-intercept` combined with either fixed-intercept flag →
     reject before compute with a clear `ValueError`. The current helper
     silently lets `--no-intercept` force defaults over fixed values; this
     refactor replaces that ambiguity with an error.

## Files modified

- `src/ldsc/regression_runner.py`
  - Replace `add_rg_arguments`: register the new flags above. Old flags
    (`--sumstats-{1,2}-file`, `--trait-name-{1,2}`, `--intercept-h2`
    nargs=2) are simply dropped — argparse's default unknown-argument
    error is sufficient.
  - Add argument validation for intercept conflicts:
    `--no-intercept` may not be combined with `--intercept-h2` or
    `--intercept-gencov`.
  - Drop the `include_h2_intercept: bool` parameter on
    `_add_common_regression_arguments` (`:1046-1060`) — always register
    `--intercept-h2` as a scalar `type=float, default=None`. Update the
    three callers (`add_h2_arguments:864`,
    `add_partitioned_h2_arguments:871`, `add_rg_arguments:884`)
    accordingly.
  - Drop the `index: int` parameter and the
    `isinstance(value, list)` branch in `_select_intercept`
    (`:851-859`). The helper becomes
    `_select_intercept(value, use_intercept, default_when_disabled)`.
    Update the three call sites in `estimate_rg:560-562` (intercept_hsq1,
    intercept_hsq2 collapse to a single resolution; intercept_gencov
    likewise).
  - Drop the `not isinstance(config.intercept_h2, list)` guard in
    `estimate_h2:403`. The condition becomes plain
    `elif config.intercept_h2 is not None:`.
  - Keep the existing `RGRegressionDataset` (`:123-142`) and
    `build_rg_dataset` (`:265-379`) semantically unchanged.
  - Refactor the current `estimate_rg` body into two layers:
    public `estimate_rg(sumstats_table_1, sumstats_table_2, ldscore_result,
    config)` still builds one `RGRegressionDataset` and returns the same
    kernel `reg.RG` object as today; private `_fit_rg_dataset(dataset,
    config)` contains the current array extraction, intercept resolution,
    `n_blocks=min(len(dataset.merged), config.n_blocks)`, and `reg.RG(...)`
    call. This helper is the only fit path used by both public `estimate_rg`
    and the new multi-trait workflow.
  - Add `RegressionRunner.estimate_rg_pairs(sumstats_tables: Sequence[
    SumstatsTable], ldscore_result: LDScoreResult, *, anchor_index: int |
    None = None, config: RegressionConfig | None = None) -> RgResultFamily`.
    Validates all traits' config snapshots against the LD-score result up
    front (one loop, mirroring `build_rg_dataset:282-288`). For
    `h2_per_trait.tsv`, fits `Hsq` once per trait via the existing
    `build_dataset(trait_i, ldscore_result)` + `estimate_h2(dataset)`
    path (independent of the per-pair Hsq inside `RG`). Iterates pairs by
    all-pairs or anchor policy; for each pair `(i, j)` builds the
    `RGRegressionDataset`, fits it with `_fit_rg_dataset(dataset, config)`,
    and emits rows from the resulting `reg.RG` plus the exact dataset
    metadata. `n_snps_used` is `len(dataset.merged)` from the dataset that was
    actually fit. Wrap each pair's build+fit+summary block in
    `try/except Exception`; failed pairs emit a NaN row with `n_snps_used=0`,
    log a WARNING with traceback (`LOGGER.warning(..., exc_info=True)`),
    record `status="failed"` and `error="<ExceptionClass>: <message>"`, and
    pair iteration continues. Do not catch `BaseException`.
  - Replace `run_rg_from_args`: resolve `--sumstats-sources` via
    `resolve_file_group(...)` (label `"munged sumstats source"`); validate
    arity ≥2; handle `--anchor-trait-file` resolution (compare on absolute
    paths AND resolved trait_name); reject `--write-per-pair-detail` without
    `--output-dir`; call `estimate_rg_pairs`. Always return the full
    `RgResultFamily` to Python and never print directly. If `--output-dir`
    is supplied, preflight the full workflow family
    before opening `workflow_logging`: data artifacts
    `[rg.tsv, rg_full.tsv, h2_per_trait.tsv, pairs/]` plus `rg.log`. Then
    pass the result to `RgDirectoryWriter`.
  - Add `RgResultFamily` dataclass mirroring `PartitionedH2BatchResult`
    (`:145-157`) but with explicit output-family fields: `rg` (8-column
    concise frame), `rg_full` (32-column diagnostic frame), `h2_per_trait`
    (h2_per_trait.tsv frame), and `per_pair_metadata` (ordered metadata
    records for the optional `pairs/` tree).
  - Add small private helpers: `_validate_rg_inputs(traits,
    ldscore_result)` (the up-front config-compat loop),
    `_iter_rg_pairs(n, anchor_index)` (yields input-order all-pairs, e.g.
    `(0,1), (0,2), (1,2)`, or `(anchor, other)` in input order),
    `_fit_rg_dataset(dataset, config)`, `_summarize_rg_pair(rg_result,
    dataset, trait_a, trait_b, pair_kind)` (builds one diagnostic row from a
    fitted `reg.RG` and its actual dataset), `_nan_rg_row(trait_a, trait_b,
    pair_kind, error_message)` (builds a NaN row for failed pairs),
    `_format_pair_error(exc)` (returns `<ExceptionClass>: <message>`),
    `_coerce_rg_numeric(value, field_name)` (preserves numeric values,
    including numeric out-of-range rg, and raises a pair-local exception for
    kernel string `NA`), and `_apply_multiple_testing(rg, rg_full)` (local
    NumPy BH and Bonferroni scatter-back).
- `src/ldsc/config.py`
  - Simplify the `RegressionConfig.intercept_h2` annotation
    (`config.py:629`) from `float | list[float] | None` to
    `float | None`. The list form was only reachable via the
    now-removed `--intercept-h2 nargs=2` CLI; no Python-API caller
    relies on it. Update the field docstring at `:614` correspondingly.
- `src/ldsc/outputs.py`
  - Add `RgDirectoryWriter` and `RgOutputConfig`, modeled exactly on
    `PartitionedH2DirectoryWriter` and `PartitionedH2OutputConfig`. The
    writer owns only scientific data artifacts:
    `[rg.tsv, rg_full.tsv, h2_per_trait.tsv, pairs/]`. `rg.log` remains
    owned by `run_rg_from_args()` and `workflow_logging()`. The
    `pairs/` subtree must be staged in a temp dir and `os.replace`d into
    place at the end (mirrors the partitioned-h2 staging at
    `outputs.py:373-405`). Inside each `pairs/<idx>_<safe_a>_vs_<safe_b>/`, write
    only `rg_full.tsv` (the one-row comprehensive view) plus
    `metadata.json`; do not duplicate the concise table at the per-pair
    level. All TSV writes use `na_rep="NaN"`.
  - Pair directory names use filesystem-safe slugs derived from resolved trait
    names, with the numeric prefix as the authoritative order/collision
    identifier. Reuse the partitioned-h2 slug policy: Unicode normalize to
    ASCII, lowercase, replace runs of non `[a-z0-9._-]` characters with `_`,
    trim leading/trailing punctuation, and fall back to `trait`.
  - Define a constant `RG_RESULT_FORMAT = "ldsc.rg_result_family.v1"` for
    the per-pair `metadata.json` files.
- `src/ldsc/cli.py`
  - Update the `rg` command description/help text to "Multi-trait genetic
    correlation".
  - After dispatching `run_rg_from_args`, if the command is `rg` and no
    `--output-dir` was supplied, call `result.rg.to_csv(sys.stdout,
    sep="\t", index=False, na_rep="NaN")`. With `--output-dir`, do not print
    the table by default.

## New files

- None. The new types (`RgResultFamily`, `RgOutputConfig`, `RgDirectoryWriter`)
  live in the same modules as their h2 / partitioned-h2 cousins to match the
  existing layout.

## Files to read (no edits expected)

- `src/ldsc/_kernel/regression.py` (RG/Gencov/Hsq, allele helpers — read-only
  callers).
- `src/ldsc/path_resolution.py` (`resolve_file_group`, preflight helpers).
- `src/ldsc/sumstats_munger.py` (`load_sumstats`, `SumstatsTable`).
- `src/ldsc/_logging.py` (`workflow_logging`, `log_inputs`, `log_outputs`).
- `src/ldsc/cli.py` (the `rg` subparser registration plus CLI-only stdout
  printing when no output directory is supplied).

## Reused existing helpers

- `resolve_file_group(tokens, suffixes=(".sumstats", ".sumstats.gz",
  ".parquet"), label="munged sumstats source")` —
  `src/ldsc/path_resolution.py:130-190`. Multi-path + glob + dedup. The
  workflow itself enforces arity ≥2.
- `load_sumstats(path, trait_name=None)` —
  `src/ldsc/sumstats_munger.py:176-252`. Per-source loader.
- `RGRegressionDataset` — `src/ldsc/regression_runner.py:123-142`. New
  dataclass produced by `build_rg_dataset`. **Unchanged.**
- `RegressionRunner.build_rg_dataset` —
  `src/ldsc/regression_runner.py:265-379`. **Unchanged.** Called once per
  pair from inside the multi-trait loop.
- `RegressionRunner.estimate_rg` —
  `src/ldsc/regression_runner.py:541-576`. Public behavior unchanged. Its
  internals delegate to `_fit_rg_dataset`, which the multi-trait path also
  uses.
- `RegressionRunner.build_dataset` —
  `src/ldsc/regression_runner.py:172-263`. Used once per trait for the
  per-trait h2 path that feeds `h2_per_trait.tsv` (separate from the rg
  preprocessing).
- `RegressionRunner.estimate_h2` and `summarize_total_h2` —
  `src/ldsc/regression_runner.py:381-419` and `:683-700` (approx; offsets
  shifted post-refactor). For `h2_per_trait.tsv`.
- Kernel `reg.RG`, `reg._filter_alleles`, `reg._align_alleles` —
  `src/ldsc/_kernel/regression.py:712`, `:780`, `:785`. Untouched.
  Invoked once per pair through `build_rg_dataset` → `_fit_rg_dataset` →
  kernel, exactly as today's single-pair flow.
- `_with_chr_pos_key`, `_assemble_regression_ldscore_table`,
  `_count_totals_for_columns`, `_select_count_key` —
  `src/ldsc/regression_runner.py`. Untouched. Used internally by
  `build_rg_dataset` and `build_dataset`; the multi-trait path does not call
  them directly.
- `_select_intercept` — `src/ldsc/regression_runner.py` (helper around the
  estimate_rg body). Adapted for scalar-broadcast.
- `_preflight_regression_outputs`, `_runner_from_args`,
  `_load_sumstats_table`, `_maybe_write_dataframe` —
  `src/ldsc/regression_runner.py` (line offsets shifted post-refactor).
  Reused where appropriate. Add an rg-specific preflight helper in the
  workflow layer so `rg.log` is preflighted before `workflow_logging()` opens
  it, while data-artifact writing remains delegated to `RgDirectoryWriter`.
- `PartitionedH2DirectoryWriter` and `PartitionedH2BatchResult` —
  `src/ldsc/outputs.py:295-410` and `regression_runner.py:145-157`. Pattern
  models only for `RgDirectoryWriter` and `RgResultFamily`.
- `workflow_logging`, `log_inputs`, `log_outputs` —
  `src/ldsc/_logging.py`.

## Implementation order (TDD)

Each step writes the failing test first, then the implementation.

1. **Glob/arity at the workflow boundary.** Test that `--sumstats-sources`
   resolving to a single file errors clearly. Implement the arity check.
2. **`add_rg_arguments` rewrite.** Test that `--sumstats-sources` accepts
   a glob and multiple tokens. (No tests for the removed flags — argparse
   rejects unknown arguments by default; that's already covered by the
   stdlib.)
3. **Intercept conflict validation.** Test that `--no-intercept` with
   `--intercept-h2` or `--intercept-gencov` raises a clean workflow error
   before compute.
4. **Trait-name collision disambiguation.** Test the `@<parent>` and
   `@<sha>` fallbacks with synthetic file layouts.
5. **`estimate_rg_pairs` skeleton (all-pairs N=3).** Test that it produces 3
   rows with the full diagnostic schema and correct trait pairings (1-2,
   1-3, 2-3) in input order.
6. **N=2 numerical equivalence.** Test that the unified path with two
   sources produces a row whose `rg`, `rg_se`, `z`, `p`, `n_snps_used`
   match exactly what `RegressionRunner.estimate_rg(t1, t2, ldscore)`
   produces on the same inputs. Because the unified path uses the same
   `_fit_rg_dataset` helper as `estimate_rg` — no alternate preprocessing
   surface — equivalence is structural rather than numerical, but the test
   still pins the contract: the row that a multi-trait run emits for a single
   pair must be the row a single-pair run would have emitted.
7. **Anchor mode.** Test N=4 with `--anchor-trait-file` produces 3 rows
   (anchor vs each of the other three) with `pair_kind="anchor"`. Test
   anchor not in sources → clean error. Test anchor matched twice via
   overlapping globs → dedupe, no self-pair.
8. **Per-pair zero overlap.** Test that one bad pair produces a NaN row but
   pair iteration continues; the bad pair has `status="failed"`, an `error`
   message in `rg_full.tsv`, a non-empty `note` in `rg.tsv`, and is logged
   as WARNING.
8b. **Arbitrary pair-local exception.** Monkeypatch the pair fit or summary
    helper to raise a generic `RuntimeError` for one pair. Assert the result
    has a failed NaN row with `error="RuntimeError: ..."`, a WARNING with
    traceback, and later pairs still run.
8c. **Kernel invalid values.** Simulate a fitted `reg.RG` with string `NA`
    in headline rg fields and assert it becomes a failed NaN row. Simulate a
    numeric rg outside `[-1, 1]` and assert the numeric value is preserved
    with `status="ok"`.
9. **`h2_per_trait.tsv`.** Test that values match running `ldsc h2`
   separately on each input.
10. **Multiple-testing columns.** Test local NumPy Benjamini-Hochberg
   against hand-computed expected values; finite p-values are adjusted with
   reverse cumulative minima and clipping, NaN/non-finite p-values stay NaN,
   `p_bonferroni` uses the finite p-value count and is capped at 1.0, and
   `p_fdr_bh` is identical between `rg.tsv` and `rg_full.tsv` for the same
   `(trait_1, trait_2)` row.
11. **Concise vs full schema parity.** Test that for every row, the eight
    concise columns in `rg.tsv` are byte-identical to the same columns
    sliced from `rg_full.tsv`. (Catches column-formatting drift.)
12. **`RgDirectoryWriter` + per-pair detail tree.** Test that the
    no-detail writer run writes `rg.tsv`, `rg_full.tsv`, and
    `h2_per_trait.tsv`; that `--write-per-pair-detail` additionally
    produces the `pairs/` tree (each subdir contains exactly `rg_full.tsv`
    and `metadata.json`); that pair directory names use safe slugs with
    numeric prefixes; that staging temp dirs are cleaned up on success; and
    that a forced failure leaves no partial subtree behind. Separately test
    the workflow writes `rg.log` only when `--output-dir` is supplied.
12b. **Per-pair detail requires output dir.** Test
    `--write-per-pair-detail` without `--output-dir` raises before loading
    sumstats or LD scores.
13. **Preflight + overwrite.** Test that `--overwrite=False` refuses on
    pre-existing artifacts (any of `rg.tsv`, `rg_full.tsv`,
    `h2_per_trait.tsv`, `pairs/`, `rg.log`); that `--overwrite=True`
    cleans stale family siblings (e.g. a stale `pairs/` from a prior
    anchor-mode run when rerunning in all-pairs mode without
    `--write-per-pair-detail`).
14. **No-output-dir stdout mode.** Test that `run_rg_from_args` returns the
    full `RgResultFamily` to Python and does not print. Test that CLI dispatch
    prints only the concise `rg.tsv`-schema table to stdout when
    `--output-dir` is absent, and does not create `rg_full.tsv`,
    `h2_per_trait.tsv`, `pairs/`, or `rg.log`.
14b. **Literal NaN rendering.** Test writer output and CLI stdout render
    missing values as literal `NaN`.
15. **Cross-validation against legacy LDSC rg.** Slow-marker test: take
    one fixture that was previously analyzed by the original `ldsc.py`
    rg (or by the existing pairwise path on N=3 traits) and assert that
    the multi-trait results agree to within 1e-6 on rg/rg_se for each pair.
    (If a legacy fixture is not on disk, generate it once via the
    current pairwise `RegressionRunner.estimate_rg` over the 3 pairs and
    pin the output as a JSON fixture; mark with `@pytest.mark.slow`.)

## Verification

Run the full test suite from the worktree root:
```
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && \
cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured && \
pytest tests/test_regression_workflow.py -v
```
End-to-end smoke test using existing test fixtures (the `write_ldscore_dir`
and `make_sumstats_table` helpers in `tests/test_regression_workflow.py`):
```
ldsc rg --sumstats-sources <fixture-dir>/trait_*.sumstats.gz \
        --ldscore-dir <fixture-ldscore-dir> \
        --output-dir /tmp/ldsc-rg-multitrait \
        --write-per-pair-detail
```
Expected: `rg.tsv` (8-column concise), `rg_full.tsv` (32-column
diagnostic), `h2_per_trait.tsv`, `rg.log`, and (with the flag)
`pairs/0001_*..0003_*` each containing `rg_full.tsv` and `metadata.json`.
Row count = C(N,2); all-pairs row order follows input order (`A-B`, `A-C`,
`B-C`); columns match the schemas in §5 and §5b;
`rg_full.tsv` has `pair_kind = all_pairs` for every row. Repeat with
`--anchor-trait-file` and confirm row count = N−1 and `pair_kind` =
`anchor`.

Repeat without `--output-dir`; expected Python return value is the full
`RgResultFamily`, `run_rg_from_args` emits no stdout, CLI stdout is the
concise 8-column table only, literal missing values appear as `NaN`, and no
diagnostic files or log are created.

## Documentation updates (implemented)

- Updated `docs/current/class-and-features.md`,
  `docs/current/data-flow.md`, `docs/current/io-argument-inventory.md`,
  `docs/current/layer-structure.md`, and
  `docs/current/path-specification.md` to describe the multi-trait rg flow,
  `RgResultFamily`, and `RgDirectoryWriter` family.
- Updated `tutorials/cross-trait-genetic-correlation.md` and
  `tutorials/cross-trait-genetic-correlation.ipynb` to use
  `--sumstats-sources`. The Markdown tutorial includes glob and
  `--anchor-trait-file` examples.
- Updated `design_map.md` to map `RgDirectoryWriter`, `RgOutputConfig`,
  `RgResultFamily`, and `estimate_rg_pairs` to their implementation modules.
- Updated rg docstrings and CLI argument help to document the no-stdout
  workflow boundary and the strongly recommended `--output-dir` diagnostics.

## Out of scope for this PR

- Hsq caching across pairs (would require `--shared-snp-set` mode and a
  small kernel-adjacent helper). Keep as a follow-up if profiling on real
  N=20+ runs justifies it.
- A shared `IndexedSubtreeWriter` extracted from `PartitionedH2DirectoryWriter`
  + `RgDirectoryWriter`. Two call sites is the right threshold; revisit
  after this PR lands.
- Liability-scale conversion of rg or h2 (kernel already has
  `gencov_obs_to_liab` / `h2_obs_to_liab`; CLI exposure is a separate
  feature).
