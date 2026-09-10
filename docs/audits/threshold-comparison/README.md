# Threshold comparison audit archive

Last updated on: 2026-09-10

The completed root [plan](task_plan.md), [progress](progress.md), and [findings](findings.md) are preserved here as historical evidence. Their last source revision was `e8a5182` (2026-06-12); relocation does not reassert their old strict-threshold conclusions.

Current retained-SNP/common-SNP MAF thresholds are inclusive. PLINK additionally excludes monomorphic SNPs (`maf > 0`) independently of the configured floor; see `src/ldsc/_kernel/plink_bed.py::__filter_snps_maf__`, `src/ldsc/_kernel/ref_panel.py::_apply_maf_filter`, and [the current SNP-universe contract](../../current/ldscore-snp-universe-contract.md). Explicit legacy conversion retains its documented strict LDSC2 common-count rule; see [legacy conversion](../../current/legacy-ldscore-conversion.md).

The former `src/ldsc/hm3_reference.py` builder mentioned in the historical findings now lives at `tools/hm3/build_hm3_chr_pos_reference.py::_filter_reference_candidates`. The public packaged resource loader is `ldsc.hm3`. Current regression threshold and window contracts live in [regression configuration](../../current/regression-configuration.md) and [reference-window behavior](../../current/ld-window-parquet-r2-sidecar-behavior.md). [Lessons](../../../lessons.md) preserve the monomorphic-SNP regression and its correction.
