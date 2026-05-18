# Multi-Trait Genetic-Correlation Refactor Design

**Date:** 2026-05-09
**Implementation status:** implemented in the `refactor/rg` worktree.
**Documentation status:** runtime docs and tutorials updated to describe
`--sumstats-sources`, anchor mode, `RgResultFamily`, `RgDirectoryWriter`, and
the rg output family.
**Scope:** `src/ldsc/regression_runner.py`, `src/ldsc/outputs.py`,
`src/ldsc/cli.py`, `src/ldsc/config.py`, `tests/test_regression_workflow.py`.

---

## Problem

The current `ldsc rg` command takes exactly two summary-statistics files
through `--sumstats-1-file` / `--sumstats-2-file` and writes a one-row
`rg.tsv`. Geneticists who want pairwise genetic correlations across a panel
of N traits have to drive `ldsc rg` from a shell loop and concatenate
per-pair outputs. That workflow is error-prone (allele alignment varies
per pair, file naming conflicts, no consistent multiple-testing column),
and the resulting tables omit the heritability and intercept context that
non-statistical users need to interpret rg.

Centralizing the preprocessing step in `RGRegressionDataset` and
`RegressionRunner.build_rg_dataset()` (commit `f36982c`) made the public
`estimate_rg` short and uniform — every rg fit goes through one merge +
allele harmonization + zero-variance drop pipeline. That uniformity is
what makes the unified multi-trait command straightforward: the workflow
iterates trait pairs through the existing `estimate_rg`, with the
orchestration layer handling source resolution, anchor mode,
multiple-testing adjustment, and output aggregation.

The deliverable is a refactored `ldsc rg` that accepts N≥2 summary
statistics, computes either all C(N,2) pairs or anchor-vs-rest, and
writes a coherent output directory with one publication-ready table, one
diagnostic table, and one per-trait heritability table. The numerical
kernel is unchanged. The single-pair preprocessing pipeline is unchanged.
Only the workflow and output layers grow.

---

## Design Decisions

| Question | Decision |
|---|---|
| CLI shape for N traits | Single flag `--sumstats-sources` accepting multiple paths and/or globs. The old `--sumstats-{1,2}-file`, `--trait-name-{1,2}` flags are removed outright (no deprecated stubs). The refactored package is unpublished, so argparse's default unknown-flag rejection is sufficient. |
| Pair selection | One unified `ldsc rg` path. With two resolved sources and no anchor, compute one pair. With three or more sources and no anchor, compute all unordered C(N,2) pairs in input order (`A-B`, `A-C`, `B-C`). With `--anchor-trait-file`, compute anchor-vs-rest pairs in input order, skipping the anchor itself. There is no separate mode or flag. |
| Anchor mode | First-class flag `--anchor-trait-file <PATH>`. The token must match exactly one resolved source path or one post-disambiguation trait label; the run yields N−1 anchor-vs-rest pairs. |
| Per-pair SNP support | Per-pair intersection (each pair's `RGRegressionDataset` is built on its own trait_1 ∩ trait_2 ∩ LD-ref SNP set). Matches single-pair `ldsc rg` semantics. |
| Compute strategy | Pure pairwise loop. `estimate_rg_pairs` builds one `RGRegressionDataset` per pair, then fits it through the same private helper used by public `estimate_rg`. No factor-out of `build_rg_dataset`, no merge caching layer. The 3.5× implementation surface and ~1 GB extra cache memory of a merge-cache design were not justified by the ~10% wall-time savings at N<10. Rejected pending real-workload profiling at N≥20. |
| Intercept controls | `--intercept-h2 <FLOAT>` is a single scalar broadcast across all traits; `--intercept-gencov <FLOAT>` is a single scalar broadcast across all pairs. `--no-intercept` constrains h2 to 1 and gencov to 0, and is rejected if combined with either fixed-intercept flag. The pre-refactor `nargs=2` form is removed; the `RegressionConfig.intercept_h2` field becomes `float | None` (the `list[float]` branch is dropped). |
| Output dual-table | `rg.tsv` is the concise, publication-ready 8-column table; `rg_full.tsv` is the 32-column diagnostic. Mirrors `partitioned_h2.tsv` / `partitioned_h2_full.tsv`. |
| Per-pair detail | Opt-in `--write-per-pair-detail`. Requires `--output-dir`; reject otherwise because the flag explicitly requests filesystem artifacts. Writes optional diagnostic directories at `diagnostics/pairs/{idx:04d}_{safe_trait_a}_vs_{safe_trait_b}/{rg_full.tsv, metadata.json}`. The root tables are complete without this tree. |
| Per-trait h2 | Always emit `h2_per_trait.tsv` from `summarize_total_h2`. Fit on each trait's full LD-ref intersection (matches `ldsc h2`); independent of the per-pair Hsq inside `RG`. |
| Self-pairs | Skipped. |
| Trait-name collision policy | If two resolved sources share `Path(path).name`, append `@<parent_dir>`; if collision persists, append `@<sha1[:8] of full path>`. INFO-log each disambiguation. |
| Trait-name UX TODO | Keep filename-derived names for v1. Future design should propagate scientific trait labels from sumstats munging metadata or a manifest-style input; `--anchor-trait-file` naming is kept for now and can be revisited with that broader trait-label design. |
| Per-pair failure | Caught in `estimate_rg_pairs` with `except Exception`; `BaseException` subclasses such as `KeyboardInterrupt` and `SystemExit` still stop the process. Emits a NaN row with `n_snps_used=0`, `status="failed"`, an `error` message formatted as `<ExceptionClass>: <message>` in `rg_full.tsv`, a concise `note` in `rg.tsv`, logs WARNING with traceback, and pair iteration continues. |
| Whole-trait failure | Hard-fail at trait load / per-trait h2 fitting time. Traits that cannot be analyzed are not silently dropped. |
| Multiple-testing columns | `p_fdr_bh` (Benjamini-Hochberg over finite p-values only) in both tables; `p_bonferroni = min(1.0, p * n_valid_pairs)` only in `rg_full.tsv`. Failed/NaN rows stay NaN in adjusted columns. |
| Result object | `RgResultFamily` has explicit fields: `rg`, `rg_full`, `h2_per_trait`, and `per_pair_metadata`. Python callers should not have to slice one ambiguous full table to recover the public output family. |
| Stdout boundary | `run_rg_from_args` returns `RgResultFamily` and does not print. `cli.py` prints `result.rg` to stdout only for CLI invocations without `--output-dir`. |
| Missing-value rendering | All TSV and stdout table writes use `na_rep="NaN"` so failed or unavailable numeric values are visible and round-trip through pandas/R readers. |
| Kernel invalid values | Numeric rg estimates are preserved even outside `[-1, 1]` or `[-1.2, 1.2]`. Non-numeric kernel `NA` values are converted to failed NaN rows with an explanatory `error`. |

---

## CLI Surface

```
ldsc rg \
    --sumstats-sources <PATH|GLOB> [<PATH|GLOB> ...]
    --ldscore-dir <DIR>
    [--anchor-trait-file <PATH>]
    [--output-dir <DIR>]
    [--write-per-pair-detail]
    [--overwrite]
    [--count-kind {common,all}]
    [--intercept-h2 <FLOAT>] [--intercept-gencov <FLOAT>] [--no-intercept]
    [--n-blocks <INT>] [--two-step-cutoff <FLOAT>] [--chisq-max <FLOAT>]
    [--log-level {DEBUG,INFO,WARNING,ERROR}]
```

Removed flags (`--sumstats-1-file`, `--sumstats-2-file`,
`--trait-name-1`, `--trait-name-2`) are dropped outright. The refactored
package is unpublished, so the default argparse "unrecognized arguments"
error is sufficient feedback. Likewise the pre-refactor `--intercept-h2
nargs=2` form is removed; the new `--intercept-h2` is registered as a
single scalar through the shared `_add_common_regression_arguments`
helper (whose `include_h2_intercept: bool` parameter is dropped — every
regression subcommand registers `--intercept-h2` the same way now).

`--sumstats-sources` resolution uses
`path_resolution.resolve_file_group(...)` — multi-token, glob-aware,
deduplicated in first-seen order. The workflow itself enforces arity ≥ 2.

`--anchor-trait-file` is resolved by both normalized absolute path **and**
the trait_name produced by `load_sumstats(...)` (after collision
disambiguation), so the user can name the anchor file by either spelling
and it still matches.

`--output-dir` is optional but strongly recommended in the argument help and
workflow docstring. The workflow always constructs the same result family in
memory: concise rg table, full diagnostic rg table, per-trait h2 table, and
per-pair metadata. When `--output-dir` is omitted, no files are written; the
Python-facing return value still carries the full result family, while the CLI
prints only the concise `rg.tsv`-schema table to stdout so quick runs remain
pipeable.

`--write-per-pair-detail` is valid only with `--output-dir`. Without an
output directory there is nowhere to place the requested `diagnostics/pairs/` diagnostic
tree, so the workflow rejects that argument combination before loading inputs.

---

## Workflow Layer

`src/ldsc/regression_runner.py`:

- **`RGRegressionDataset` (`:123-142`) and `build_rg_dataset`
  (`:265-379`)** — unchanged. Every pair still gets its own per-pair
  trait_1 ∩ trait_2 ∩ LD-ref SNP set.
- **`estimate_rg` (`:541-576`)** — public behavior unchanged. Its internals
  are refactored so it builds one `RGRegressionDataset` and delegates fitting
  to a private helper that accepts an already-built dataset. This lets
  `estimate_rg_pairs` reuse the same fit path while recording `n_snps_used`
  without rebuilding the dataset or attaching ad hoc metadata to the kernel
  object.
- **`estimate_rg_pairs(sumstats_tables, ldscore_result, *,
  anchor_index=None, config=None) -> RgResultFamily`** (new). Validates
  configs for all traits up front. Computes `dataset_h2_i =
  build_dataset(trait_i, ldscore_result)` + `hsq_i =
  estimate_h2(dataset_h2_i)` once per trait for `h2_per_trait.tsv`
  (independent of per-pair Hsq inside `RG`). Iterates pairs by all-pairs
  or anchor policy. For each pair `(i, j)`, builds the `RGRegressionDataset`,
  fits it through the same private helper used by `estimate_rg`, and emits a
  row from the fitted `reg.RG` plus dataset metadata. Wraps each per-pair
  build+fit+summary in `try/except Exception` to convert pair-local data or
  numerical failures into NaN rows; `BaseException` subclasses such as
  `KeyboardInterrupt` and `SystemExit` are not caught.
- **`run_rg_from_args(args)`** (rewritten). Resolves
  `--sumstats-sources`, validates arity, loads each `SumstatsTable`,
  applies trait-name collision disambiguation, resolves
  `--anchor-trait-file` to a list index, and calls `estimate_rg_pairs`.
  The return value is always the full `RgResultFamily`. When `--output-dir`
  is supplied, it preflights `rg.tsv`, `rg_full.tsv`, `h2_per_trait.tsv`,
  optional `diagnostics/pairs/`, and `diagnostics/rg.log` before opening `workflow_logging`, then
  passes the result to `RgDirectoryWriter`. When `--output-dir` is omitted,
  it writes no files and creates no log. It never prints directly; CLI stdout
  is handled by `src/ldsc/cli.py`. `diagnostics/rg.log` remains owned by the workflow
  logging context rather than the writer.
- **`RgResultFamily`** (new dataclass; mirrors `PartitionedH2BatchResult` at
  `regression_runner.py:145-157`). Fields: `rg` (the concise 8-column
  table), `rg_full` (the full diagnostic table), `h2_per_trait` (one row per
  input trait), and `per_pair_metadata` (ordered records for the optional
  `diagnostics/pairs/` tree).
- **Small private helpers**: `_validate_rg_inputs(traits,
  ldscore_result)` (the up-front config-compat loop, mirroring
  `build_rg_dataset:282-288`); `_iter_rg_pairs(n, anchor_index)` (yields
  `(i, j)` tuples in input order: `(0,1)`, `(0,2)`, `(1,2)` for all-pairs,
  or `(anchor, other)` in input order for anchor mode);
  `_fit_rg_dataset(dataset, config)` (contains the current array extraction
  and `reg.RG(...)` call from `estimate_rg`);
  `_summarize_rg_pair(rg_result, dataset, trait_a, trait_b, pair_kind)` (builds
  one diagnostic row from a fitted `reg.RG` plus its actual dataset metadata);
  `_nan_rg_row(trait_a, trait_b, pair_kind, error_message)` (builds a NaN
  row for failed pairs); `_format_pair_error(exc)` (returns
  `<ExceptionClass>: <message>`); `_coerce_rg_numeric(value, field_name)` (converts
  kernel `NA` strings to a failed-row condition while preserving finite and
  non-finite numeric values).

### Pre-refactor surfaces dropped (no backward-compatibility shims)

The package is unpublished. The following pre-refactor code paths exist
*only* to support the old CLI shape and are removed wholesale rather
than aliased:

- `add_rg_arguments` (`:884`) — the old `--sumstats-{1,2}-file`,
  `--trait-name-{1,2}`, and `--intercept-h2 nargs=2` registrations are
  deleted. No deprecated stubs.
- `_add_common_regression_arguments(parser, include_h2_intercept: bool)`
  (`:1046-1060`) — the boolean parameter is dropped; every regression
  subcommand registers `--intercept-h2` as a single scalar via the
  shared helper.
- `_select_intercept(value, index, ...)` (`:851-859`) — `index` and the
  `isinstance(value, list)` branch are removed. The helper becomes
  `_select_intercept(value, use_intercept, default_when_disabled)`.
- `RegressionConfig.intercept_h2: float | list[float] | None`
  (`config.py:629`) — collapses to `float | None`. Field docstring at
  `config.py:614` updated.
- `estimate_h2:403` — `not isinstance(config.intercept_h2, list)` guard
  drops out (the list form no longer exists).
- The kernel and the `.ldscore` / `.M` / `.M_5_50` / `w_ld` legacy LDSC
  reference-panel formats are **kept** unchanged. They are the only
  pre-refactor compatibility surface this package preserves.

`src/ldsc/outputs.py`:

- **`RgOutputConfig`** and **`RgDirectoryWriter`** (new). Modeled on
  `PartitionedH2OutputConfig` / `PartitionedH2DirectoryWriter`
  (`outputs.py:265-410`). The writer owns only scientific data artifacts:
  `[rg.tsv, rg_full.tsv, h2_per_trait.tsv, diagnostics/pairs/]`. `diagnostics/rg.log` remains owned
  by `run_rg_from_args()` and `workflow_logging()`. Per-pair
  subtree uses the same staged-then-`os.replace` pattern as
  `PartitionedH2DirectoryWriter:373-405`.
- **`RG_RESULT_FORMAT = "ldsc.rg_result_family.v1"`** — recorded in each
  `diagnostics/pairs/<idx>_<safe_a>_vs_<safe_b>/metadata.json`.

---

## Computational Strategy

`estimate_rg_pairs` is a pure pairwise loop over the same dataset builder and
kernel fit path as the public `estimate_rg`. There is no caching and no
factor-out of `build_rg_dataset`; each pair pays the full cost of its own
LD-score merge and per-pair preprocessing, exactly as today's single-pair
`ldsc rg` does. The only new private helper is the fit helper that accepts an
already-built `RGRegressionDataset`, so the pair loop can report
`n_snps_used` from the dataset it actually fit.

The single statistical correctness principle is maintained automatically:
every `reg.RG` call must see `hsq1`, `hsq2`, and `gencov` fit on the
**same** SNPs in the **same** order so the jackknife `RatioJackknife` at
`_kernel/regression.py:735-745` combines element-aligned delete-block
vectors. `build_rg_dataset` already builds an `RGRegressionDataset.merged`
on the per-pair intersection, and the kernel `reg.RG` runs all three
regressions on that one frame, so the alignment holds by construction.

For N=2, the unified `ldsc rg --sumstats-sources a b` path builds the same
dataset and uses the same private fit helper as public `estimate_rg`, so its
row in `rg.tsv` is the row a single-pair run would emit. Equivalence is
structural, not just numerical — there is no separate preprocessing surface
that could drift.

### Why no caching

A merge-cache design that hoisted `pd.merge(trait_i, ldscore_frame)` out
of the per-pair loop and stored it per trait was considered. Estimated
benefit at N<10 with one baseline LD column: ~10% wall-time reduction
(~1-2 minutes off a ~15-minute run for N=9 in coordinate-family mode; less in
rsID-family mode). Estimated cost: ~0.8-1.2 GB of cached state plus ~3.5x the
implementation surface (~230 vs ~85 LOC) and a permanent byte-identity
test against `build_rg_dataset`. The savings sit below the
user-noticeable wall-clock threshold at this scale and the maintenance
liability is permanent. If preprocessing ever becomes the bottleneck on
real workloads, parallelizing the per-pair loop across processes
(independent pairs, shared read-only `LDScoreResult`) is a bigger and
more direct win than caching, and is more naturally added later.

---

## Output Layout

When `--output-dir` is supplied:

```
<output-dir>/
  rg.tsv             # concise, publication-ready (8 cols)
  rg_full.tsv        # comprehensive diagnostic (32 cols)
  h2_per_trait.tsv   # one row per input trait
  diagnostics/
    metadata.json    # provenance only
    rg.log           # from workflow_logging
    pairs/           # only when --write-per-pair-detail
    0001_<safe_trait_a>_vs_<safe_trait_b>/
      rg_full.tsv    # one-row, comprehensive view
      metadata.json
    0002_...
```

`diagnostics/pairs/` index ordering follows the result row order: input-order all-pairs
(`A-B`, `A-C`, `B-C`) or `(anchor, other)` in input order for anchor mode.
Folder names use a filesystem-safe slug derived from the resolved trait name;
the numeric prefix is the stable identity for ordering and collision
avoidance.

When `--output-dir` is omitted, no files are written. The workflow still
constructs the same result family in memory. The Python-facing return value
contains the concise rg table, full diagnostic rg table, per-trait h2 table,
and per-pair metadata. `run_rg_from_args` returns that object without printing.
For CLI invocations only, `src/ldsc/cli.py` prints the concise
`rg.tsv`-schema table to stdout, so users can pipe quick runs without creating
a directory. The help text should still steer real analyses toward
`--output-dir` because diagnostic files and logs are much easier to inspect in
a result directory.

All table writes, including stdout, render missing values with literal `NaN`
(`to_csv(..., sep="\t", index=False, na_rep="NaN")`).

---

## `rg.tsv` Schema (Concise, Publication-Ready)

8 columns, one row per pair:

| Column | Meaning |
|---|---|
| `trait_1` | Resolved trait name of the left source. |
| `trait_2` | Resolved trait name of the right source. |
| `n_snps_used` | SNPs in the per-pair intersection used by `reg.RG`. |
| `rg` | Genetic correlation point estimate. |
| `rg_se` | Block-jackknife SE of `rg`. |
| `p` | Two-sided normal p-value from `rg / rg_se`. |
| `p_fdr_bh` | Benjamini–Hochberg adjusted p-value over the rows of this file. |
| `note` | Empty for successful pairs; for failed pairs, a short pointer such as `Failed; see rg_full.tsv error column; use --output-dir for diagnostics/rg.log`. |

This file is what users will paste into manuscripts. It mirrors what most
genetic-correlation papers report verbatim (e.g. "rg = 0.68, SE = 0.04,
P_FDR = 4.4e-9").

---

## `rg_full.tsv` Schema (Comprehensive Diagnostic)

32 columns, one row per pair, in this order:

| Group | Columns |
|---|---|
| Identifiers + size | `trait_1`, `trait_2`, `n_snps_used` |
| Headline rg | `rg`, `rg_se`, `z`, `p`, `p_fdr_bh`, `p_bonferroni` |
| Per-pair h2 | `h2_1`, `h2_1_se`, `h2_2`, `h2_2_se` |
| Genetic covariance | `gencov`, `gencov_se` |
| Intercepts | `intercept_h2_1`, `intercept_h2_1_se`, `intercept_h2_2`, `intercept_h2_2_se`, `intercept_gencov`, `intercept_gencov_se` |
| Confounding ratio | `ratio_1`, `ratio_1_se`, `ratio_2`, `ratio_2_se` |
| Inflation | `lambda_gc_1`, `lambda_gc_2`, `mean_chisq_1`, `mean_chisq_2` |
| Provenance | `pair_kind` ∈ {`all_pairs`, `anchor`} |
| Failure status | `status` ∈ {`ok`, `failed`}, `error` |

Notes:
- `p_bonferroni = min(1.0, p * n_valid_pairs)` where
  `n_valid_pairs` is the number of finite p-values being adjusted.
- `p_fdr_bh` is identical to the value in `rg.tsv` for the same row.
- Benjamini-Hochberg is implemented locally with NumPy. Sort finite p-values
  ascending, compute `p_sorted * m / rank`, apply reverse cumulative minima
  to enforce monotonic adjusted values, clip to `[0, 1]`, and scatter back to
  original row order. NaN or non-finite p-values are ignored during adjustment
  and remain NaN in both adjusted-p columns.
- `h2_1` / `h2_2` here are fit on the **per-pair intersection** (so they
  match the Hsq used inside `reg.RG`). They will differ from the
  corresponding rows in `h2_per_trait.tsv`, which are fit on the full
  trait-vs-LD-ref intersection — this is intentional and documented.
- `ratio_t = (intercept_h2_t − 1) / (mean_chisq_t − 1)`. Available
  directly from `RG.hsq{1,2}.ratio` and `.ratio_se`.
- `intercept_gencov_se` is sourced from `RG.gencov.intercept_se`. Together
  with `intercept_gencov` it forms the cross-trait sample-overlap
  diagnostic (Bulik-Sullivan et al. 2015, equation 4 and surrounding
  discussion).
- `status="failed"` rows carry NaN-filled statistics, `n_snps_used=0`, and
  the caught exception message in `error`, formatted as
  `<ExceptionClass>: <message>`.
- Numeric rg estimates are preserved even when they fall outside the usual
  interpretive range. Non-numeric kernel `NA` values are converted to failed
  NaN rows with an explanatory `error`.

---

## `h2_per_trait.tsv` Schema

One row per resolved sumstats source. Columns reuse `summarize_total_h2`
output verbatim:

| Column | Meaning |
|---|---|
| `trait_name` | Resolved trait name (post-collision-disambiguation). |
| `n_snps` | SNPs used in the per-trait h2 fit (trait ∩ LD-ref). |
| `total_h2` | Trait-level heritability point estimate. |
| `total_h2_se` | Block-jackknife SE. |
| `intercept` | Hsq intercept. |
| `intercept_se` | SE of the intercept. |
| `mean_chisq` | Mean χ² across used SNPs. |
| `lambda_gc` | Genomic inflation factor. |
| `ratio` | `(intercept − 1) / (mean_chisq − 1)`. |
| `ratio_se` | SE of the ratio. |

This table is what users would have gotten from running `ldsc h2`
separately on each input. Computing it inside `ldsc rg` saves a manual
loop and keeps the output directory self-contained.

---

## Per-Pair Detail Tree

`diagnostics/pairs/` is an optional per-pair diagnostic tree. It is not needed for the
main scientific result because the root `rg.tsv` and `rg_full.tsv` already
contain every pair. It is created only when `--write-per-pair-detail` is set,
for users who want one small directory per pair with a one-row full diagnostic
table and reproducibility metadata.

When `--write-per-pair-detail` is set:

```
diagnostics/pairs/<idx:04d>_<safe_trait_a>_vs_<safe_trait_b>/
    rg_full.tsv       # exactly one row (the same one as in the root rg_full.tsv)
    metadata.json     # reproducibility record
```

Per-pair `metadata.json` schema:

| Field | Meaning |
|---|---|
| `schema_version` | Current metadata schema version. |
| `artifact_type` | `rg_pair_result`. |
| `trait_1`, `trait_2` | Resolved trait names. |
| `source_1`, `source_2` | Absolute paths of the input sumstats files. |
| `n_snps_used` | After all merges and filters. |
| `n_blocks_used` | `min(n_snps_used, config.n_blocks)`. |
| `intercept_h2_policy` | `"free"` or `"fixed:<value>"`. |
| `intercept_gencov_policy` | `"free"` or `"fixed:<value>"`. |
| `chisq_max_used`, `two_step_cutoff_used` | Effective filter values. |
| `pair_kind` | `"all_pairs"` or `"anchor"`. |
| `status` | `"ok"` or `"failed"`. |
| `error` | Set to the exception message if the pair failed; absent otherwise. |

`safe_trait_a` and `safe_trait_b` are filesystem-safe slugs derived from the
resolved trait names. Use the same slugging policy as partitioned-h2
per-query directories: Unicode normalize to ASCII, lowercase, replace runs of
non `[a-z0-9._-]` characters with `_`, trim leading/trailing punctuation, and
fall back to `trait` if the slug would be empty. The numeric prefix preserves
the authoritative pair order.

The concise `rg.tsv` is **not** duplicated inside `diagnostics/pairs/<...>/` — users
who open the per-pair directory are already in diagnostic mode.

---

## Error Handling and Edge Cases

- **Glob resolves to <2 distinct files** → workflow raises `ValueError`
  with a message pointing at `--sumstats-sources` before any compute.
- **`--anchor-trait-file` does not match any resolved source** (compared
  by absolute path AND resolved trait_name) → hard error before compute.
- **Anchor file matched both directly and via another glob** → dedupe so
  the anchor never self-pairs.
- **Trait-name collisions** across resolved sources → disambiguate
  by `@<parent_dir>`, then by `@<sha1[:8] of full path>`. Each
  disambiguation is INFO-logged.
- **Trait with zero LD-ref overlap** → hard fail at the per-trait h2
  fitting step (matches `estimate_h2`'s current behavior). Traits that
  cannot be analyzed are not silently excluded from the run.
- **`--write-per-pair-detail` without `--output-dir`** → reject before
  loading inputs. This flag requests a filesystem tree, so silently ignoring
  it would be poor UX.
- **Per-pair failure** — any ordinary `Exception` during the pair's
  build+fit+summary phase → caught in `estimate_rg_pairs`; emit a row with
  `n_snps_used=0`, NaN-filled stats, `status="failed"`, an `error` message
  in `rg_full.tsv`, a concise `note` in `rg.tsv`, log WARNING with traceback
  when `diagnostics/rg.log` exists, and continue. Do not catch `BaseException`
  subclasses such as `KeyboardInterrupt` and `SystemExit`.
- **Kernel `NA` values** → if `reg.RG` returns non-numeric `NA` for rg, rg SE,
  z, or p (for example h2 out of bounds), treat the pair as failed and emit a
  NaN row with an explanatory `error`. Preserve numeric out-of-range rg values
  as numeric estimates.
- **`--no-intercept` plus `--intercept-h2`/`--intercept-gencov`** → reject
  before compute with a clear `ValueError`. The current code silently lets
  `--no-intercept` force defaults over fixed intercept values; the new rg
  redesign intentionally replaces that ambiguous precedence with an error.

---

## Testing

`tests/test_regression_workflow.py` (TDD; see plan file for the failing-test-first ordering):

1. Glob/arity at the workflow boundary.
2. `add_rg_arguments` rewrite — `--sumstats-sources` accepts a glob and
   multiple tokens. (No tests for the removed flags; argparse handles
   unknown arguments natively.)
3. Intercept conflict validation (`--no-intercept` with fixed-intercept
   flags rejects before compute).
4. Trait-name collision disambiguation (`@<parent>` → `@<sha>` fallback).
5. `estimate_rg_pairs` skeleton (all-pairs N=3 produces 3 rows).
6. **N=2 equivalence.** The row that the unified path emits for `(t1, t2)`
   matches what `RegressionRunner.estimate_rg(t1, t2, ldscore)` produces
   on the same inputs. Equivalence is structural — the unified path calls
   the same private dataset-fit helper as `estimate_rg`, so the same
   `reg.RG` instance is being summarized.
7. Anchor mode (resolution by path or trait_name; dedup; clean error if
   the anchor is not in sources).
8. Per-pair zero overlap → NaN row, `status="failed"`, `error`, non-empty
   concise `note`, pair iteration continues.
8b. Arbitrary per-pair `Exception` → NaN row, formatted `error`, WARNING with
    traceback, later pairs still run.
8c. Kernel `NA` values → failed NaN row; numeric out-of-range rg values are
    preserved.
9. `h2_per_trait.tsv` matches running `ldsc h2` separately on each input.
10. Multiple-testing columns: local NumPy Benjamini-Hochberg implementation
   matches hand-computed expected values, ignores NaN/non-finite p-values,
   leaves adjusted values as NaN for failed pairs, and keeps `p_fdr_bh`
   identical between `rg.tsv` and `rg_full.tsv` for the same row.
11. Concise vs full schema parity: the eight concise columns in `rg.tsv`
    are byte-identical to the same columns sliced from `rg_full.tsv`.
12. `RgDirectoryWriter` + per-pair detail tree (staged-then-`os.replace`;
    no partial subtree on forced failure).
12b. `--write-per-pair-detail` without `--output-dir` rejects before loading
     inputs.
13. Preflight + overwrite covering the full
    workflow-owned family: `diagnostics/rg.log` plus data-owned
    `[rg.tsv, rg_full.tsv, h2_per_trait.tsv, diagnostics/pairs/]`.
14. No-output-dir mode still returns the full `RgResultFamily` to Python
    callers, `run_rg_from_args` does not print, CLI invocations print only the
    concise table, and no files or log are created.
14b. All file/stdout table rendering uses literal `NaN` for missing values.
15. Slow-marker cross-validation against legacy LDSC `rg` output.

Existing test fixtures `make_ldscore_result`, `make_sumstats_table`,
`write_ldscore_dir` (in `tests/test_regression_workflow.py`) are reused
verbatim.

---

## Out of Scope

- **Per-trait LD-merge caching** (the merge-cache design considered
  during planning). Estimated ~10% wall-time gain at N<10 with one
  baseline LD column, ~0.8–1.2 GB extra cached state, ~3.5×
  implementation surface. Below the user-noticeable threshold for the
  stated scale; the plan-mode review explicitly rejected it. Revisit
  only with profiling evidence from real N≥20 workloads.
- **Hsq caching across pairs.** Statistically incompatible with the
  per-pair SNP intersection (the kernel's `RatioJackknife` requires
  `hsq1`, `hsq2`, and `gencov` delete-block vectors to come from the
  same SNP partition). Would require a `--shared-snp-set` mode that
  trades coverage for speed.
- **Per-pair-loop parallelization across processes.** Each pair is
  independent and the `LDScoreResult` is read-only, so a process pool
  scales near-linearly in cores — a much bigger win than caching for
  large-N runs. Deferred to a future PR; needs separate design for
  worker lifecycle, log aggregation, and result ordering.
- **A shared `IndexedSubtreeWriter`** extracted from
  `PartitionedH2DirectoryWriter` + `RgDirectoryWriter`. Two call sites
  is the right threshold; revisit after this PR lands.
- **Liability-scale conversion.** Kernel already exposes
  `gencov_obs_to_liab` and `h2_obs_to_liab`; CLI exposure
  (`--samp-prev`, `--pop-prev`) is a separate feature. `rg` itself is
  invariant to scale conversion, so its concise table is unaffected.
- **Per-pair full delete-block jackknife arrays.** Useful for users who
  want to redo the SE computation themselves, but storage cost is large
  and the audience is small. Add later if requested.
