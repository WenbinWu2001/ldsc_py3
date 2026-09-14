
# Lessons

Last updated on: 2026-09-14

## Refactoring must preserve attribution outside removable banners

- Summary/root cause/correction: Parser docstring replacement and munging-banner removal dropped upstream notices during refactoring; preserve them as source notices, carry credit through extracted helpers, and check that LICENSE and NOTICE ship in both distribution formats.

## Streamed summary assembly must preserve its schema

- Summary/root cause/correction: Reconstructing a batch summary from per-source pandas Series preserved repeated row labels and inferred plain integer types; reset the shared row index and explicitly retain nullable count and diagnostic dtypes before publishing or comparing summaries.

## Wrapper retirement must account for unmigrated callers

- Summary/root cause/correction: Mechanical wrapper removal during a staged migration left an active caller and an empty import block; search every call site before retirement, replace still-needed forwarding calls directly, and parse edited modules before running the affected suite.
- Recurrence: Removing the last imported symbol during schema-helper consolidation again left an empty parenthesized import; remove its entire import statement and parse every edited module immediately, including import-only edits.

## Parser defaults must not be obtained by parsing invented required paths
- Summary/root cause/correction: Python wrappers for LD-score and reference-panel workflows parsed fake required paths (`placeholder`/`out`) to obtain argparse defaults, which let omitted user arguments become unintended filesystem targets; collect action defaults directly and validate required wrapper arguments before constructing the namespace.

## Cleanup-capable preflights must declare every artifact the current run will write
- Summary: A direct partitioned LD-score overwrite successfully wrote
  `ldscore.overlap.parquet`, then the workflow's post-success stale-artifact
  cleanup deleted that newly written file while `metadata.json` still pointed
  to it.
- Root cause: The direct workflow performs an early preflight before expensive
  computation and the canonical writer performs a second exact preflight. The
  early preflight's `owned_paths` included overlap, but its `produced_paths`
  prediction omitted overlap. With overwrite enabled, the pre-existing overlap
  was therefore captured as stale; path identity meant the later cleanup
  unlinked the replacement written at the same location.
- Correction: Declare artifact ownership once per directory writer and share it with early preflight, including chromosome drop reports; discard early stale predictions and let the final writer reconcile actual production. Cover write/overwrite/reload transitions, checking metadata-listed files, regenerated reports, obsolete overlap removal, and unrelated-file preservation.

## Regression-weight LD scores belong in the shared PLINK projection pass
- Summary: Direct PLINK LD-score calculation and gene-index construction projected baseline/query annotations first, then reset the genotype cursor and recomputed the same correlation blocks for the one-column regression-SNP mask.
- Root cause: The regression-weight output was treated as a separate result family instead of one more annotation projection, obscuring that `ldScoreVarBlocks` recomputes genotype correlations on every call.
- Correction: Append the binary regression mask after all ordinary annotation columns, call `ldScoreVarBlocks` once, then split partitioned scores from `w_ld`; keep the separately reset gene-atom batch calls because those bound index-builder memory.

## ldsc:schema_version is a shared identity contract, not a file-layout version
- Summary: Bumping `ldsc:schema_version` to 2 to mark the new index R2 layout broke `build_reader`'s `validate_identity_artifact_metadata`, which pins `schema_version == SCHEMA_VERSION` (=1) across all artifact types.
- Root cause: `schema_version` is the cross-artifact provenance-contract version (sumstats/ldscore/ref-panel share it), enforced equal to the package constant. It is not a per-format layout marker.
- Correction: Keep `schema_version=1`; identify the index layout structurally (presence of `IDX_1/IDX_2/SIGN` columns + `ldsc:n_snps` + `ldsc:sidecar_identity_sha256`). The build->read parity test missed this because it uses `compute_chrom_from_parquet`, not `build_reader`; the autofill test (build_reader path) caught it.

## Inferred projection builds must reach coordinate-based sibling features
- Summary/root cause/correction: rsID gene-list runs inferred `gene_catalog_build` but default region exclusions consulted only build-independent identity metadata, so reuse the concrete projection build for named regions and cover the combined path rather than testing normalization and exclusion resolution separately.

## scipy CSR @ dense accumulates in the operand dtype — use float64 for LD-score scatter
- Summary: Replacing `np.add.at` with a scipy.sparse CSR `U @ annot` for the parquet LD-score scatter lost ~3e-3 precision (above the int16 quantization floor) when `U.data`/`annot` were float32.
- Root cause: scipy's CSR·dense SpMM sums each row in the operand dtype; float32 operands accumulate hundreds of within-window terms in float32. `np.add.at` had hidden this by accumulating into a float64 `cor_sum` even from float32 products.
- Correction: build `U.data` and cast `annot` to float64 before the SpMM (`cor_sum = annot64 + U@annot64 + U.T@annot64`). Float64 SpMM is exact vs a float64 reference and still ~75x faster than `np.add.at`; the only cost is 8-byte CSR data (12*K bytes/chunk instead of 8*K).

## Manual dataclass reconstruction silently drops newly-added fields
- **Summary:** A new `LDScoreResult.overlap` field was written to disk but missing
  from the result returned by `run_ldscore`.
- **Root cause:** `_replace_result_output_paths` rebuilt `LDScoreResult(...)` by
  listing every field by hand, so it dropped `overlap` (added later). The writer
  ran before this rebuild, so the sidecar existed but the returned object lost it.
- **Correction:** Use `dataclasses.replace(result, ...)` instead of manual
  field-by-field reconstruction so future fields carry over automatically. An
  end-to-end test (write to disk + assert the returned result) caught what every
  unit test missed.

## MAF `>` carried an implicit monomorphic-SNP guard in the PLINK reader
- **Summary:** Flipping the genotype MAF filter from strict `>` to inclusive
  `>=` (threshold-boundary harmonization) silently broke monomorphic-SNP
  exclusion and changed golden LD outputs.
- **Root cause:** `plink_bed.py` used `np.minimum(f, 1 - f) > mafMin` with the
  default `mafMin=0`. The strict `> 0` did double duty: it applied the user MAF
  floor *and* dropped monomorphic (folded MAF == 0, zero-variance) SNPs. Pure
  `>= 0` keeps those invariant SNPs, so polymorphic count `m` and the normalized
  genotype matrix changed (golden `m`: 4 → 6).
- **Correction:** Separate the two concerns — keep `maf > 0` as an explicit
  monomorphic guard and apply the user floor inclusively: `maf > 0 and
  maf >= mafMin`. Other MAF sites did not need this (their thresholds are
  skipped when unset, or default to 0.05 which is always > 0).
- **Takeaway:** Before flipping a comparison operator, check whether the
  boundary value (here, 0) encodes a second invariant beyond the nominal
  threshold.

## Key/row alignment is not genomic ordering — verify order at the boundary
- **Summary:** While auditing CHR:POS sort dependencies, I claimed legacy `.bed`
  ldscore mode was an unguarded gap (raw `.bim` order flowing into the
  forward-only `get_block_lefts` window) and that the regression merge would put
  sumstats into LD-score order. Both were wrong on mechanism.
- **Root cause:** I reasoned about *upstream* order instead of tracing it. (1)
  `compute_chrom_from_plink` discards `.bim` order: `keep_snps` is built from the
  already-sorted annotation metadata, `PlinkBEDFile.__filter_snps_maf__` preserves
  `keep_snps` order, and the panel merge is `sort=False` — so `geno_meta` is sorted
  by construction. (2) `pd.merge(left, right, how="inner", sort=False)` aligns
  *values* by key but keeps the *left* frame's row order; a key match never adopts
  the right frame's order.
- **Correction:** For the window precondition, add a boundary tripwire
  (`validate_window_positions_sorted`, raises `LDSCInternalError`) rather than
  trusting upstream sorts. For jackknife block order (Gap C), make the
  genomically-sorted LD-score frame the *left* side of the merge so the result
  inherits its order — in every identifier mode, since the order originates from
  the LD-score artifact, not the sumstats.
- **Takeaway:** "Aligned" ≠ "ordered." When a downstream step needs genomic order
  (two-pointer windows, contiguous jackknife blocks), assert it at the boundary or
  establish it by construction at the point of use; do not infer it from an
  upstream sort you did not trace.

## Squash-merging `restructure` into `main`: apply the diff, not `git merge --squash`
- **Summary:** A plain `git merge --squash restructure` onto `main` produced
  spurious conflicts and resurrected files that `restructure` had deleted.
- **Root cause:** `main` is a chain of independent squash commits, so it shares
  no real ancestry with `restructure`. The 3-way merge base is stale, so git
  treats `restructure`'s deletions as `main`-side additions and re-adds them
  (e.g. `misc/generate_pgc_top50_munge_commands.py`).
- **Correction:** Reset `main` to `origin/main`, then make its tree exactly equal
  the authoritative `restructure` tip via
  `git diff --binary origin/main restructure | git apply --index`, verify
  `git diff restructure` is empty, and commit one `Squash merge restructure`.
  Full runbook in `docs/release.md`.
- **Takeaway:** When two branches share content but not history (squash/rebase
  workflows), do not rely on 3-way merge; reconstruct the target tree
  deterministically from the authoritative side and assert tree-equality.

## Missing metadata must be explicit in reusable whitespace-parsed tables
- **Summary/root cause/correction:** `ldsc annotate` wrote empty `CM` fields that regex-whitespace readers collapsed, shifting every query annotation left and filling the final column with NaN; write `CM=NA`, require numeric non-missing annotation values at the shared loader seam, and cover the writer-to-LD-score-parser round trip with distinguishable columns.

## Adapter correctness must reach the production numerical path

- Summary/root cause/correction: Raw R² metadata was resolved by the panel adapter but bypassed by the chromosome workflow, silently omitting bias correction; make `RefPanel.prepare_chromosome` the shared production/test path for aligned rows, windows, and reader policy, then assert known ordinary and regression-weight scores through the public calculator; keep physical BIM/BED indices through sorting and close chromosome readers on success or failure.

## Numerical error policies belong around the operations that need them

- Summary/root cause/correction: Import-time `np.seterr` calls made behavior depend on import order, while eager `np.where` division crashed zero-variance quantile tests under strict settings; preserve caller policy, scope required estimator exceptions with `np.errstate`, and use masked division for undefined quantile tests.

## Read-only CLI modes must bypass the entire failure-marker lifecycle

- Summary/root cause/correction: The overwrite marker wrapper treated argparse's successful help exit as a failed run; exclude help and inference-only parsing checks from the marker scope and verify absent directories and existing markers remain untouched on success and failure.

## Munging summaries must use parser accounting

- Summary/root cause/correction: Munging discarded chunk QC counters and reconstructed input counts from physical lines while returning provenance through mutated parser attributes; return a typed result with parser row counts, exclusive stage removals, and provenance, and test conservation on raw files containing blank lines and overlapping rejection reasons, including bad-allele rows deferred by keep-list filtering to final identity cleanup.

## Test runners that share pybedtools temporary files must run sequentially

- Summary/root cause/correction: Concurrent pytest and unittest runs intermittently lost BED overlap files because `pybedtools.cleanup(remove_all=True)` deletes other processes' files in the shared temporary directory; run these compatibility checks sequentially or isolate their temporary roots.

## Empty batch members can change pandas boolean-mask dtype

- Summary/root cause/correction: Concatenating empty and populated gene lists coerced an internal boolean column to object, so bitwise inversion produced integer indices instead of a boolean mask; use an explicit boolean comparison and test mixed empty/nonempty batches through the real workflow.

## Import rewrites need compilation, not only AST parsing

- Summary/root cause/correction: Mechanical import edits left empty import groups and later moved future imports below ordinary imports; AST parsing did not catch the latter, so compile every edited module before running behavioral tests and preserve future-import placement.


## A saved workflow result must detach from its construction owner

- Summary/root cause/correction: Standalone annotation finished writing but retained private preparation shards through the returned bundle, so successful CLI runs left scratch behind; return persistent query descriptors with original baseline dependencies, defer any later read preparation, close the construction owner before return, and verify both immediate scratch cleanup and saved-query reads after original gene inputs are removed.

## Float32 reduction order can invalidate aggregate metadata

- Summary/root cause/correction: A completion audit found batch-dependent quantitative annotation counts after column-major inputs became detached row-major reads and row tiles were reduced separately in float32; reduce normalized annotation values in float64 and verify signed/high-offset counts against independent `math.fsum` values across batch widths and C/F layouts without loosening tolerances. LD-score equivalence on binary fixtures does not establish quantitative count equivalence. See the [completion review](docs/audits/annotation-memory/completion-review.md).

## Optional coordinates must not shift frequency columns in TSV

- Summary/root cause/correction: Frequency round-trip testing exposed blank CHR/POS fields collapsing under the whitespace sumstats reader and shifting Z/N/FRQ; serialize missing fields as `NA` and verify emitted Parquet and gzip TSV preserve selected frequency, including values above 0.5 and sub-three-decimal precision.

## Worker-initialization tests must restore package logging

- Summary/root cause/correction: Calling the pool initializer in the pytest process left the LDSC logger at WARNING and suppressed later INFO capture; scope its logging changes with `caplog.at_level(..., logger="LDSC")` and verify the formerly failing test order.

## Packaged HM3 restriction includes allele identity

- Summary/root cause/correction: Documentation described the packaged HM3 keep-list as allele-free after inspecting the coordinate-only inference view; the full curated map contains A1/A2. Inspect the actual restriction resource and reader, document allele-aware matching, and validate synthetic munging fixtures with compatible reference alleles or explicit base identity.

## Tutorial correlation plots need finite jackknife uncertainty

- Summary/root cause/correction: The 20-SNP rg tutorial produced finite rg with NaN SE because one delete-block heritability product was negative; increase the synthetic Z magnitudes to keep these toy fits positive, assert finite rg/SE before plotting, and execute every tutorial code cell without weakening the plotting validation.

## Diagnostic scratch must follow execution ownership

- Summary/root cause/correction: Prebuilt query batches reused the long-lived input workspace for per-batch drop reports, so completed diagnostics accumulated despite releasing query tables; give every execution batch its own scratch owner, including batches requiring no query projection, and assert earlier diagnostic paths are gone before the next batch.

## NumPy raw staging needs contiguous row tiles

- Summary/root cause/correction: Replacing numeric pickles with raw staging initially slowed wide annotation scans because `ndarray.tofile()` wrote pandas' column-major blocks element by element; convert each bounded tile with `np.ascontiguousarray()` before writing and benchmark scan time separately from downstream reads.

## PLINK filename parsing must be shared across workflows

- Summary/root cause/correction: Dotted chromosome prefixes failed in gene-index construction while direct LD scores worked because filename-based discovery stripped `.22` with `Path.stem` and a separate preflight had already switched to BIM contents; centralize PLINK discovery, trio validation, and chromosome assignment, pass resolved mappings to workers, and test equivalent token forms through both numerical workflows.

## Failure markers must not block their own retry

- Summary/root cause/correction: A failed gene-index overwrite created a marker-only output directory that the next run rejected as unrecognized; both preflight and publication used a diagnostics predicate that omitted `RUN_FAILED.txt`. Recognize the regular root marker alongside owned diagnostics, preserve rejection of unrelated contents, and test repeated failure followed by successful publication and marker removal.

## Output ownership must cover normalization and recovery

- Summary/root cause/correction: A package-wide retry audit found raw-path logs/audits/scratch and markers diverging from normalized result roots, CLI markers choosing the first repeated destination, and index recovery rejecting its own marker while diagnostic scans overlooked empty directories and symlinks; reuse normalized destinations, follow argparse's final value, share index preflight, inspect every diagnostic entry, and test the existing owned-backup recovery with a retained failure marker.

## Worker counts need one CPU-budget policy

- Summary/root cause/correction: The thread-option audit confirmed that gene-index construction used machine-wide CPU counts while both LD-score modes used CPU affinity, and only indexed-query API calls rejected non-integer thread values; share validation and worker resolution, then select inline execution from the resolved count. Implemented in `_parallelism.py` and verified through real builds and both scoring workflows under restricted affinity, including effective-one executor avoidance. See `docs/audits/2026-09-14-thread-option-consistency.md`.
