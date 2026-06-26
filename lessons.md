
## ldsc:schema_version is a shared identity contract, not a file-layout version
- Summary: Bumping `ldsc:schema_version` to 2 to mark the new index R2 layout broke `build_reader`'s `validate_identity_artifact_metadata`, which pins `schema_version == SCHEMA_VERSION` (=1) across all artifact types.
- Root cause: `schema_version` is the cross-artifact provenance-contract version (sumstats/ldscore/ref-panel share it), enforced equal to the package constant. It is not a per-format layout marker.
- Correction: Keep `schema_version=1`; identify the index layout structurally (presence of `IDX_1/IDX_2/SIGN` columns + `ldsc:n_snps` + `ldsc:sidecar_identity_sha256`). The build->read parity test missed this because it uses `compute_chrom_from_parquet`, not `build_reader`; the autofill test (build_reader path) caught it.

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
