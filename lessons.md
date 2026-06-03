
## ldsc:schema_version is a shared identity contract, not a file-layout version
- Summary: Bumping `ldsc:schema_version` to 2 to mark the new index R2 layout broke `build_reader`'s `validate_identity_artifact_metadata`, which pins `schema_version == SCHEMA_VERSION` (=1) across all artifact types.
- Root cause: `schema_version` is the cross-artifact provenance-contract version (sumstats/ldscore/ref-panel share it), enforced equal to the package constant. It is not a per-format layout marker.
- Correction: Keep `schema_version=1`; identify the index layout structurally (presence of `IDX_1/IDX_2/SIGN` columns + `ldsc:n_snps` + `ldsc:sidecar_identity_sha256`). The build->read parity test missed this because it uses `compute_chrom_from_parquet`, not `build_reader`; the autofill test (build_reader path) caught it.
