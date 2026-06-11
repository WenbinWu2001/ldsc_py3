
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
