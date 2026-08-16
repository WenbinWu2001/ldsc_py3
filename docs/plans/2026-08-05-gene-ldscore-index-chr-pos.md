# Explicit chr_pos gene LD-score index implementation plan

Last updated on: 2026-08-16

Status: implemented; focused and chromosome-22 validation complete

## Goal and success signal

Extend `ldsc build-gene-ldscore-index` to require an explicit base SNP identity (`rsid` or `chr_pos`) and explicit `hg19`, while preserving direct/index scientific equivalence and the broad-reference-versus-regression-row universe split. Success requires exact direct/index float32 LD scores, rows, ordering, counts, controls, and downstream result tables in both modes; overlap values alone use `rtol=0`, `atol=1e-7`.

The historical governing contract has been superseded by `docs/specs/2026-08-15-defensive-gene-list-input-design.md`. The existing rsID implementation and 2026-08-03 audit were the baseline, not evidence for the new coordinate mode.

## Context and constraints

- Work on branch `restructure`; preserve unrelated changes and do not commit without an explicit request.
- Builder choices are exactly `rsid` and `chr_pos`; both require explicit `--genome-build hg19`. There is no default, `auto`, inference, liftover, hg38, or allele-aware builder mode.
- Direct PLINK-backed LD-score calculation remains the scientific oracle. Indexed Stage 2 inherits mode/build from the index and accepts no live override.
- Mutable baseline and PLINK sources use vectorized drop-all cleanup by active effective key. Restriction keys remain set-like. Immutable artifacts fail on duplicate identities.
- Only regression/output rows are persisted. The broad retained PLINK universe supplies LD contributions, counts, and overlaps.
- Canonical sorting must retain raw BIM/BED column indices and gather genotypes with the final ordered raw-index vector. Sorted DataFrame row numbers must never become BED indices.
- Reuse `ldsc._kernel.snp_identity` and existing coordinate/restriction helpers. Chromosome-level orchestration is allowed; Python per-SNP loops are not. Pause for review if a new per-SNP loop appears necessary.
- The new strict index schema rejects previously built gene indexes. Already published ordinary LD-score directories remain outside that compatibility break.

## Interfaces and invariants

- `GeneLDScoreIndexBuildConfig` requires `genome_build` and `snp_identifier` with no defaults and accepts only `hg19` plus `rsid|chr_pos`.
- Builder omission errors are deliberate and occur before scientific input resolution. Indexed `ldscore` detects live identity/build options by argument presence.
- `chr_pos` effective identity is normalized `(CHR, POS)` with positive 1-based positions. PLINK supplies published `CHR,POS,SNP,A1,A2`; baseline/restriction labels are ignored for coordinate matching.
- Direct and builder cleanup write `diagnostics/dropped_snps/chr<chrom>_dropped.tsv.gz`, including header-only files, and emit summarized warnings.
- Canonical component rows are `CHR SNP POS A1 A2 regression_ld_scores <ordered baseline columns>` and are unique under the declared effective mode.
- Root and components repeat mode/build/index ID. Components store ordered effective-key and published-row metadata digests, which the loader recomputes before assembly.
- Indexed canonical output is self-contained and uses ordinary downstream identity, genome-build, and explicit downgrade behavior.

## Implementation slices

### Slice 1 — Make identity/build explicit at public entry points

**Outcome:** CLI and Python construction cannot enter the builder without explicit approved values; indexed mode cannot accept a live override.

**Likely areas:** `src/ldsc/gene_ldscore_index.py`, `src/ldsc/config.py`, `src/ldsc/ldscore_calculator.py`, CLI dispatch/help, `tests/test_gene_ldscore_index.py`, `tests/test_ldscore_workflow.py`.

**TDD checkpoint:** first add failing tests for each missing/invalid flag, both missing together, required Python fields, accepted `rsid|chr_pos`, rejected allele-aware/auto/hg38 values, and explicit indexed overrides including values equal to the index. Implement argument-provenance tracking without changing ordinary direct-mode defaults.

### Slice 2 — Share vectorized source cleanup, alignment, and BED permutation

**Outcome:** Direct PLINK and builder paths use the same mode-aware drop-all and canonical-alignment seam without per-SNP Python iteration.

**Likely areas:** `src/ldsc/_kernel/snp_identity.py`, `src/ldsc/_kernel/ref_panel.py`, `src/ldsc/_kernel/ldscore.py`, `src/ldsc/ldscore_calculator.py`, `src/ldsc/gene_ldscore_index.py`, identity and workflow tests.

**TDD checkpoint:** add failing base-mode tests for duplicate rsIDs, duplicate coordinates, differing labels at matching coordinates, baseline-only/PLINK-only rows, empty intersections, and diagnostic sidecars. Add an unsorted BIM fixture with distinguishable genotype columns and independently asserted raw BED indices. Then reuse vectorized key/cleanup helpers, carry `_raw_index` through cleanup/intersection/sort, gather BED columns once by the final vector, and verify no new row loop appears in the diff.

### Slice 3 — Extend builder restrictions, semantic identity, and strict artifacts

**Outcome:** The builder constructs `chr_pos` components, writes the new identity-bearing schema, and rejects old/tampered/mismatched indexes before use.

**Likely areas:** `src/ldsc/gene_ldscore_index.py`, `src/ldsc/_kernel/identifiers.py`, `src/ldsc/_kernel/gene_ldscore_index.py`, `src/ldsc/outputs.py`, `tests/test_gene_ldscore_index.py`, `tests/test_gene_ldscore_index_kernel.py`.

**TDD checkpoint:** add failing tests for packaged HM3 hg19 coordinate matching, custom `CHR+POS|hg19_POS` restrictions, ignored coordinate-mode SNP labels, repeated restriction-key collapse, PLINK-authored output metadata, mode-separated `index_id`, required root/component mode/build fields, both row digests, canonical schema/order, and intentional rejection of old or tampered indexes. Preserve the broad contributor universe and fail only when aggregate regression rows are empty.

### Slice 4 — Make indexed output coordinate-aware and self-contained

**Outcome:** Indexed assembly derives configuration and membership exclusively from the validated index and publishes ordinary coordinate-aware LD-score artifacts.

**Likely areas:** `src/ldsc/gene_ldscore_index.py`, `src/ldsc/ldscore_calculator.py`, `src/ldsc/outputs.py`, `src/ldsc/regression_runner.py`, `tests/test_gene_ldscore_index.py`, `tests/test_ldscore_workflow.py`, `tests/test_regression_workflow.py`.

**TDD checkpoint:** add failing tests proving `GlobalConfig`, metadata, and regression membership use `chr_pos+hg19`; PLINK labels are passive during matching but retained in output; baseline/query identity rows agree exactly; output remains usable after the source index is unavailable; h2/rg/partitioned-h2 use ordinary exact-family/build rules; and downgrade remains explicit.

### Slice 5 — Prove direct/index equivalence and finish user-facing documentation

**Outcome:** Focused and empirical evidence covers both modes, and active help/current/wiki/tutorial/troubleshooting surfaces describe only the approved contract.

**Likely areas:** direct/index integration fixtures, `docs/audits/`, `docs/current/`, `docs/wiki/utility-functionalities/build-gene-ldscore-index.md`, `docs/wiki/main-functionalities/ldscore.md`, CLI help, docstrings, troubleshooting.

**Validation checkpoint:** run the deterministic two-mode fixture including non-HM3 and MHC/centromere contributors, controls, focal queries, counts, overlaps, duplicate-mode contrasts, build contradictions, and the BED permutation tripwire. Refresh the chromosome-22 audit for `rsid` and add `chr_pos`; require exact persisted LD scores and downstream partitioned-h2 tables, with only overlap values using `atol=1e-7`. Smoke indexed output through h2 and rg. Document that advanced callers must supply same-build baseline/PLINK inputs and that no liftover or build inference occurs.

## Validation milestones

Run the smallest red test first in each slice, then the relevant focused group:

```bash
pytest -q tests/test_gene_ldscore_index.py tests/test_gene_ldscore_index_kernel.py
pytest -q tests/test_ldscore_workflow.py tests/test_regression_workflow.py
pytest -q tests/test_config_identifiers.py tests/test_snp_identity.py tests/test_annotation.py
```

Before completion:

```bash
pytest -q
python -m unittest discover -s tests -p 'test*.py' -v
ldsc build-gene-ldscore-index --help
ldsc ldscore --help
python -m ldsc --help
git diff --check
```

Record commands, resources, row counts, tolerances, and downstream results for the real chromosome-22 gate. Parser tests or a loadable index alone are not completion evidence.

## Risks and checkpoints

- **BED permutation:** verify early with an unsorted, distinguishable fixture before changing builder numerics.
- **Direct-oracle drift:** share cleanup/alignment primitives; do not maintain two identity implementations.
- **Passive versus semantic labels:** baseline/restriction labels are ignored only in `chr_pos`; PLINK labels and alleles still affect semantic identity because they are published.
- **Artifact break:** strict required metadata intentionally rejects old indexes; errors must say to rebuild.
- **Diagnostics plumbing:** retain structured drop frames through direct and builder orchestration without adding row loops or making diagnostics scientific identity.
- **Numerical equivalence:** investigate any non-overlap discrepancy instead of widening tolerance.
- **Dirty documentation:** layer edits onto existing user changes in current/wiki files and review combined diffs before finalizing.

## Out of scope

- hg38, genome-build inference, or liftover;
- allele-aware gene-index modes;
- R2-backed index construction;
- a package-wide duplicate-policy audit beyond the direct PLINK and gene-index source seams;
- migration/loading compatibility for old gene indexes;
- per-payload numerical checksums or signed artifact authenticity;
- production publication of whole-genome indexes.

## Completion definition

The feature is complete only when explicit CLI/Python behavior, both source cleanup paths, strict artifacts, indexed coordinate output, direct/index scientific equivalence, downstream partitioned-h2 equivalence, h2/rg consumption, active documentation, and the two-mode chromosome-22 audit all pass. Parser-only or unit-only success is insufficient.

## Completion record

All five implementation slices are complete. Focused tests cover the explicit
public contract, mode-aware vectorized drop-all cleanup, raw BIM/BED permutation,
restriction semantics, strict root/component validation and tamper detection,
and exact direct/index operators in both modes. The production-size chromosome-22
gate produced exact direct/index canonical outputs for `rsid` and `chr_pos`;
coordinate output passed `h2` and `rg`, and direct/index `partitioned-h2` tables
were byte-identical. See
[`docs/audits/2026-08-05-gene-ldscore-index-chr-pos-chr22.md`](../audits/2026-08-05-gene-ldscore-index-chr-pos-chr22.md).
The repository-wide gates also passed: 1,211 pytest cases after rerunning five
sandbox-blocked multiprocessing checks with the required OS access, and 962
standard-library unittest cases. Each suite had one expected skip.
