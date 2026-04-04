# Plan: Current Status and Next Steps for `ldscore_new.py`

## Goal
- Replace the old LD-score estimation workflow with a new implementation that:
  - supports both PLINK reference panels and parquet `R2` tables,
  - keeps LDSC-style outputs compatible with downstream use,
  - is callable both from Python and from the command line,
  - remains simple enough to maintain.

## Current Status

### Implemented
- `ldscore_new.py` now consumes SNP-level annotation files only. BED-to-SNP projection remains outside this module.
- Two annotation groups are supported:
  - `--baseline-annot` / `--baseline-annot-chr`
  - `--query-annot` / `--query-annot-chr`
- Two reference-panel modes are implemented in the same script:
  - PLINK mode
  - parquet `R2` mode
- The module now supports both:
  - Python function call via `run_ldscore(...)`
  - CLI execution via `python ldsc_py3_Jerry/ldsc_new.py ...`
- `w_ld` is explicitly generated and written as a separate output:
  - `<out>.w.l2.ldscore.gz`
- Aggregated multi-chromosome output is the default.
- Per-chromosome output is supported with `--per-chr-output`.
- Header documentation was rewritten and includes example usage.

### Implemented Public API
- Python:
  - `ldsc_py3_Jerry.ldscore.ldscore_new.run_ldscore(...)`
- CLI:
  - `python ldsc_py3_Jerry/ldsc_new.py ...`
- BED utility has been aligned to the same naming pattern:
  - Python: `ldsc_py3_Jerry.utils.run_bed_to_annot.run_bed_to_annot(...)`
  - CLI: `python ldsc_py3_Jerry/utils/run_bed_to_annot.py ...`

## Current Interface

### Annotation Inputs
- Annotation files are SNP-level `.annot` / `.annot.gz` style tables.
- Each file must contain the first four columns:
  - `CHR BP SNP CM`
- Each file may contain one or more annotation columns after those metadata columns.
- Files from the same chromosome are column-bound after validating that SNP rows match exactly.

### Reference Panel Inputs
- PLINK:
  - `--bfile`
  - `--bfile-chr`
- parquet `R2`:
  - `--r2-table`
  - `--r2-table-chr`

### Matching and Metadata
- Supported SNP identifier modes:
  - `chr_pos`
  - `rsID`
- No allele matching is implemented in v1.
- parquet `chr_pos` mode requires explicit `--genome-build {hg19,hg38}`.
- Optional metadata inputs:
  - `--frqfile`
  - `--frqfile-chr`
- Optional regression SNP input:
  - `--regression-snps`

### Output Files
- Reference LD scores:
  - `<out>.l2.ldscore.gz`
- Regression-weight LD scores:
  - `<out>.w.l2.ldscore.gz`
- Annotation counts:
  - `<out>.l2.M`
- Common-SNP annotation counts, when MAF is available:
  - `<out>.l2.M_5_50`
- Annotation group manifest:
  - `<out>.annotation_groups.tsv`

## What the New Code Actually Does

### Shared Workflow
1. Resolve chromosomes from the supplied annotation inputs.
2. Read baseline/query annotation files chromosome by chromosome.
3. Normalize SNP identity by the chosen identifier mode.
4. Require all annotation files for a chromosome to have the same retained SNP rows.
5. Build one dense annotation matrix per chromosome.
6. Dispatch to either the parquet backend or the PLINK backend.
7. Write LDSC-style outputs, aggregated by default.

### parquet Backend
1. Read only the needed parquet columns.
2. Derive SNP keys from either `rsID` or chromosome-position columns.
3. Intersect the annotation SNP universe with SNPs present in the parquet table.
4. Build a sparse symmetric `R2` matrix.
5. Add the diagonal explicitly as `r_jj^2 = 1`.
6. If `--r2-bias-mode raw`, apply the LDSC unbiased correction:
   - `r2_unbiased = r2_raw - (1 - r2_raw) / (n - 2)`
7. Compute:
   - partitioned/reference LD scores from the annotation matrix
   - `w_ld` from a one-column regression SNP mask

### PLINK Backend
1. Align annotation SNPs to the BIM table.
2. Reuse the legacy LD-score kernel via the old PLINK genotype backend.
3. Compute:
   - partitioned/reference LD scores
   - `w_ld` using the same LD-score kernel with a one-column regression mask
4. Write outputs in the new output layout.

## Important Clarifications Locked by Current Work

### `w_ld`
- In the new code, `w_ld` is treated as a separate one-column LD-score output.
- It is not one weight per annotation.
- This matches how the old regression loader expects `--w-ld`:
  - one SNP column
  - one LD-weight column

### Regression SNP Meaning
- `--regression-snps` defines the SNP set used to build `w_ld`.
- If omitted, the retained reference SNP set is used.
- This matches the intended LDSC meaning:
  - `w_ld` is precomputed from a chosen regression SNP universe,
  - later downstream regression may merge/subset rows,
  - but LD scores themselves are not recomputed after that merge.

### Missing MAF
- If MAF is unavailable:
  - `MAF` is written as `NA`
  - `.M_5_50` is skipped

### Missing CM
- `--ld-wind-cm` requires non-missing CM for retained SNPs.
- `--ld-wind-kb` and `--ld-wind-snps` do not.

## Differences Between the Original Plan and the Current Code

### Completed from the Original Plan
- Header rewrite
- Separate baseline/query annotation inputs
- parquet `R2` support
- explicit `r2` bias mode handling
- explicit regression SNP support
- separate `w_ld` output
- aggregated default output
- true Python function API for the new module

### Updated Design Decisions
- `run_ldscore(...)` now exists as a keyword-argument Python API.
- The CLI-facing execution body lives in `run_ldscore_from_args(...)`.
- The BED utility module was renamed from `bed_to_annot.py` to `run_bed_to_annot.py` so the module name matches the callable/script name.

### Not Yet Implemented
- Full end-to-end replacement for the old `ldsc.py` regression pipeline.
- Numerical equivalence validation against the old PLINK-based outputs.
- A new `sumstats` / `h2` / `rg` pipeline based on the new interface.
- Allele-aware SNP reconciliation.
- Explicit support for more flexible R2 schemas beyond the currently assumed parquet columns.

## Remaining Work for Functional Parity

### High Priority
- Validate that PLINK mode in `ldscore_new.py` reproduces old `ldsc.py --l2` outputs to numerical tolerance.
- Validate parquet mode against equivalent PLINK-derived LD scores.
- Add tests for:
  - diagonal inclusion,
  - raw vs unbiased parquet `R2`,
  - per-chromosome vs aggregated output,
  - `--regression-snps` changing only `w_ld`

### Medium Priority
- Decide whether `ldsc_new.py` will remain only an LD-score entrypoint or eventually become a full replacement for old `ldsc.py`.
- Refactor downstream regression code if the long-term goal is a unified function-first API for:
  - LD score estimation
  - heritability regression
  - genetic correlation

### Low Priority
- Add a small compatibility shim for any old imports that still reference `utils.bed_to_annot`.
- Consider optional richer output validation or manifests once parity is stable.

## Recommended Near-Term Roadmap
1. Finish numerical validation of `ldscore_new.py` against old PLINK-based LDSC.
2. Lock the parquet `R2` file contract with one documented example schema.
3. Add regression-facing tests confirming the produced `.l2.ldscore`, `.w.l2.ldscore`, `.M`, and `.M_5_50` files are readable by the existing downstream code.
4. Only after that, decide whether to refactor the old `sumstats.py` / `ldsc.py` workflow into the same function-first structure.

## Acceptance Criteria for the Next Milestone
- `run_ldscore(...)` can be called directly from Python with keyword arguments.
- `python ldsc_py3_Jerry/ldsc_new.py ...` runs the same workflow as the Python API.
- PLINK mode and parquet mode both produce:
  - reference LD scores,
  - `w_ld`,
  - `.M`,
  - `.M_5_50` when possible
- Outputs are accepted by the existing downstream regression readers.
- The remaining gap to the old pipeline is clearly limited to downstream regression refactoring and validation, not LD-score generation itself.
