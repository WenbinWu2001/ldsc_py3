# Status

Last updated on: 2026-10-28

Branch `restructure`; checked against HEAD `952703b` on 2026-10-28.
Nothing below is in progress; each open item needs the user's go-ahead.

## Completed: domain subpackage reorganization of `src/ldsc`

- The reorganization landed on `restructure` on 2026-10-28, from `9236b80` (dead-code cleanup) to the documentation commit `docs: describe the domain layout`.
  [Code structure](docs/dev/code-structure.md) and [architecture](docs/dev/architecture.md) describe the layout, and `tests/test_layering.py` enforces its import tiers.
- Follow-up cleanup landed on 2026-10-28 in the `refactor` commits from `cabe537` to `5b3cede`:
  - `window_coordinates` and `get_block_lefts` in `src/ldsc/_refpanel/windows.py` serve both `ldscore` and `build-r2-panel`; each command keeps its own missing-CM error.
  - `_runtime.artifacts.atomic_write_json` is the one atomic JSON writer, and every file it writes ends with a newline; `h2_scale` uses the shared `atomic_write_dataframe`.
  - The gene-index build no longer migrates `<index>.build/` logs or takes the standalone lock of earlier releases.
  - `_core.identifiers.open_text` is the one text-file opener; it decodes UTF-8 and treats `.gz` in any letter case as gzip.
  - Four dead or test-only helpers were removed, `_genome.liftover.load_hm3_curated_map(path)` became `load_hm3_map_file`, and three modules stopped re-exporting names they do not define.
  - The 24 underscore-prefixed names that other domains imported lost the underscore (`_runtime.preflight._readable` became `check_readable`), so an underscore again means module-local.
- Look-alike helpers deliberately kept, because each carries a domain-specific error message or return type: the four `_require_columns`, the two `_read_json`, and the two `get_compression`.
- Remaining follow-up, not started: consider per-domain test folders.

## Annotation SNP grids and LD-score annotation coverage

The earlier handoff for this work was written at `4e2fa3c`, before the implementation; this section replaces it.

### Confirmed decisions

These five decisions are settled; do not reopen them.

1. Annotation definitions (gene lists with their catalog and padding, or BED files) are separate from the annotation SNP grid, which selects the SNPs that receive them.
2. `ldsc annotate` takes exactly one annotation SNP grid: a BIM grid (`--ref-snp-bim-sources`) or a temporarily retained annotation grid (`--ref-snp-annot-sources`).
3. Baseline annotations may supply the annotation grid to `annotate`; downstream they remain the conditioning annotations for LD scores and regression, and choosing a BIM grid does not remove them.
4. The annotation SNP grid defines the candidate SNPs, reference-panel filtering and explicit restrictions select the reference SNPs, and regression SNPs are a separate subset; R² panels take SNP identities and coordinates from their metadata sidecars.
5. Annotation coverage is validated during `ldscore` reference preparation, before LD-score computation, without zero-filling annotation values or silently shrinking the reference SNPs.

### Implemented

- Decisions 1 to 3 landed in `3ef653a` (`feat(annotate): add BIM grids and aligned exports`, 2026-10-20).
  The CLI uses the mutually exclusive `--ref-snp-bim-sources` and `--ref-snp-annot-sources`; the BIM grid reads `.bim`/`.bim.gz` only and creates no baseline values.
  The contract is [docs/reference/annotate.md](docs/reference/annotate.md); later commits `12cb3fd`, `1b928d2` and `e00fca9` added BED-query and reference-panel failures and that reference page.
- Carryover audit findings now resolved:
  successful partial runs print stderr notices (`emit_resolved_only_notice()` and `emit_gene_gate_b_notice()` in `src/ldsc/_annotation/query_outcomes.py`; contract in `annotate.md`, "Command-line notices");
  gene-list naming collisions are reported together with rejected rows (`naming_error` in `_gate_a_error()`);
  [data-flow.md](docs/dev/data-flow.md) describes both BED and gene-list routes.

### Implemented: decision 5 in `ldscore` (2026-10-28)

The user confirmed the design in a grilling session on 2026-10-28 and authorized its implementation the same day.
The contract is in [reference SNPs on the annotation SNP grid](docs/reference/ldscore-snp-universes.md#reference-snps-on-the-annotation-snp-grid).

- The annotation files define the annotation SNP grid; the reference SNPs are the grid SNPs among the reference-panel SNPs that pass the panel's filters.
  Panel SNPs off the grid are excluded without failing: the user's recorded exception to the AGENTS.md rule that unrequested, analysis-changing filtering must fail, because legacy baseline-LD annotations exist only for 1000 Genomes Phase 3 SNPs.
  With the user's approval, that AGENTS.md rule now carries a one-line pointer to the exception.
  No CLI option was added, and `ldscore` does not read the grid kind recorded by `annotate`.
- The exclusion is never silent: `diagnostics/reference_snp_funnel.tsv` (two funnels per chromosome, fixed rows per backend, `NA` for steps that did not run), one log WARNING and a stderr notice after successful CLI runs, and `LDScoreResult.reference_snp_funnel` for in-memory runs.
  PLINK counts the grid step before genotype QC; counts come from batch 1 only.
- A chromosome emptied at the grid step now fails the run, naming every such chromosome (in Gate B for gene-list runs); one emptied by `--ref-panel-snps-file` or genotype QC is skipped with a logged WARNING instead of a Python `UserWarning`.
- `snp_universe_policy.ld_reference_universe.selection` was removed from `metadata.json`; the reference-panel source is in the log's `Inputs:` section.
- Gene LD-score indexes are format version 4: each chromosome stores `reference_snp_funnel.parquet`, which replaces `pre_qc_rows` and `annotation_intersection_removed` (the latter counted every `.bim` row, not the chromosome's); indexed scoring writes the inherited funnel and repeats the WARNING.
- Baseline and prebuilt query files still need identical rows in identical order; that row-mismatch error now explains the shared annotation SNP grid and links a new troubleshooting section.
- The earlier decision-5 example ("a BIM grid A, B, C with a baseline covering A and B") already failed through row alignment before this work, so it was not the missing check; the tests in `tests/test_reference_snp_funnel.py` cover the grid-versus-panel behavior instead.

### Resolved follow-ups (2026-10-28)

- `annotate`'s gene-list preflight error (`_gate_a_error()`) now reuses `gene_list_gate_a_message()`, the `ldscore` wording: per-source rejected/submitted counts, up to ten examples with the omitted count, a suggested audit filter and the `--overwrite` rerun, followed by any naming collision and the troubleshooting link.
  A message with no rejected rows no longer prints an empty `()`.
- rsID-family runs on an R² panel no longer print "genome_build is set but will be ignored": `_artifact_global_config()` drops the build that the panel records instead of passing it to `GlobalConfig`.
  The warning remains for a build the caller supplies.

## Reading map

- [annotate contract](docs/reference/annotate.md) and [LD-score SNP universes](docs/reference/ldscore-snp-universes.md).
- `src/ldsc/_refpanel/adapters.py` (`_align_annotations()`, `PlinkRefPanel._load_genotypes()`, `RefPanel.prepare_chromosome()`), `src/ldsc/_refpanel/plink.py` (per-SNP genotype QC), `src/ldsc/_ldscore/preflight.py` (chromosome-coverage gate) and `src/ldsc/_ldscore/calculator.py` (`_snp_universe_policy()`, batch-1 counts, chromosome skips).
- `src/ldsc/_annotation/annotate.py` (`run_annotate()`, `_gate_a_error()`), `src/ldsc/_annotation/sources.py` (`_aligned_metadata()`) and `src/ldsc/_ldscore/gene_index/build.py` (per-chromosome evidence).
- [Lessons](lessons.md), particularly reference preparation, physical BIM/BED row mapping, and output ownership.
