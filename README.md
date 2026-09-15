# ldsc3_Jerry

Last updated on: 2026-09-15

This repository is the active refactored LDSC package.

Current flag names and retired options are documented in the [LDSC2-to-current CLI flag map](docs/current/legacy-cli-flag-map.md). Munging uses `--input-format` and automatically preserves input frequency as `FRQ`; both `ldscore` and `build-gene-ldscore-index` retain `--threads`. It controls chromosome worker processes for direct `ldscore` and a chromosome thread pool for index building; both default to `1` (sequential).

Both commands and indexed LD scoring require nonzero integer worker requests. Negative requests use process CPU affinity when available, with a machine CPU-count fallback. All requests are capped at chromosome count; an effective count of 1 runs inline. Python booleans and non-integer values are rejected. See the [shared worker policy](docs/current/config-design.md).

Run `ldsc COMMAND --help` for options grouped by task, with required inputs first and descriptions of defaults, dependencies, and mode restrictions. Contributors should follow the [CLI help guidelines](docs/current/cli-help-guidelines.md).

## Authors and maintainer

Authors: **Wenbin Wu, Anthony Abrantes, Brendan Bulik-Sullivan, and Hilary Finucane**.

Maintainer: **Wenbin Wu** ([drkwu999@gmail.com](mailto:drkwu999@gmail.com)). Report package issues in the [issue tracker](https://github.com/WenbinWu2001/ldsc_py3/issues).

Brendan Bulik-Sullivan and Hilary Finucane developed upstream LDSC; Anthony Abrantes contributed the Python 3 port and performance work; Wenbin Wu develops and maintains this refactored package. See [NOTICE](NOTICE) for preserved upstream attribution and dated modifications.

## License

This project continues to use the **GNU General Public License version 3 (GPLv3)**, consistent with upstream LDSC. The complete license is in [LICENSE](LICENSE), which is included in source distributions and wheels together with [NOTICE](NOTICE). Upstream copyright notices remain applicable to derived code; current maintainer information does not replace them.

GPLv3 permits use, modification, and commercial use. Distributing covered derivative software requires compliance with GPL terms, including applicable licensing and corresponding-source obligations. Merely using the package for analysis does not automatically place papers, input data, or ordinary analysis results under GPL. See [LICENSE](LICENSE), sections 0, 2, and 4–6, for the governing terms.

The upstream declarations inspected specify GPLv3 but do not explicitly distinguish `GPL-3.0-only` from `GPL-3.0-or-later`. We retain `GPLv3` in package metadata without assuming either SPDX identifier; the evidence is recorded in [NOTICE](NOTICE). Bundled datasets are assessed independently in [the resource attribution record](src/ldsc/data/ATTRIBUTION.txt), including unresolved source and reuse-term questions.

## Citation

If you use this refactored software, cite **Wenbin Wu, Anthony Abrantes, Brendan Bulik-Sullivan, and Hilary Finucane. ldsc3_Jerry (version 2.0b0)**, with [the repository URL](https://github.com/WenbinWu2001/ldsc_py3) and the exact version or commit used. [CITATION.cff](CITATION.cff) provides the software citation metadata. The installed distribution, import, and command are named `ldsc`; `ldsc3_Jerry` is the current README and citation title.

Also cite the scientific methods used in your analysis:

| Analysis | Method reference |
| --- | --- |
| LDSC, SNP heritability, or regression intercept | [Bulik-Sullivan et al. (2015), basic LD Score regression](https://doi.org/10.1038/ng.3211) |
| Genetic correlation | [Bulik-Sullivan et al. (2015), genetic correlations](https://doi.org/10.1038/ng.3406) |
| Partitioned heritability | [Finucane et al. (2015), functional annotations](https://doi.org/10.1038/ng.3404) |
| Continuous annotations | [Gazal et al. (2017), LD-dependent architecture](https://doi.org/10.1038/ng.3954) |
| Cell-type-specific expression analysis | [Finucane et al. (2018), specifically expressed genes](https://doi.org/10.1038/s41588-018-0081-4) |

These academic credit instructions are separate from the software license. Cite source datasets and any additional methods actually used as appropriate.

## Structure

- `src/ldsc/`: public package surface
- `src/ldsc/_kernel/`: internal compute and file-format modules
- `tests/`: local parity and workflow tests
- `tutorials/`: package-level usage examples
- [Contributor entry](docs/current/code-structure.md): authoritative module navigation and links to current workflow, scientific, and artifact contracts

## Large pathway batches

The annotation and regression workflows are designed for large pathway batches, such as **1,000 pathways in one run, each tested separately against shared baseline categories**. Whole-genome and chromosome-sharded annotation inputs are prepared in bounded chunks sized for the files processed together. Binary columns are detected automatically and stored as packed bits; continuous columns use dense float32. Two input passes validate/classify supplied values, then write retained values directly into private chromosome artifacts, with metadata kept separately for global identity cleanup. Public reads remain float32 in input column order, and no new option is needed. Chromosome working data are released between sequential runs, and parallel work is bounded by `--threads`; annotation preparation itself remains serial. `--query-batch-size` defaults to `1000` for direct/indexed `ldscore` and batch `partitioned-h2`; use a smaller positive value to reduce active query workspace. Final HM3 LD-score tables remain aggregate Parquet files, never public chromosome LD shards. See [annotation format policy](docs/current/annotation-memory-design.md#annotation-format-policy) for classification and measured storage sizes.

Standalone `annotate` accepts BED files or gene lists and writes reusable chromosome annotation shards incrementally. Source-backed Python preparation requires an output directory and explicit bundle ownership; new scratch stays under that destination, with the existing gene-index construction transaction as the exception. See the [developer memory design](docs/current/annotation-memory-design.md) and [pathway batch guide](docs/wiki/main-functionalities/partitioned-h2.md#testing-enrichment-for-a-large-batch-of-pathways).

## Install

Choose the stable branch for routine use, or the `restructure` branch for the
current development version. Both installs use an editable package install so
local source changes are picked up immediately.

### Stable version

Clone the default `main` branch to a directory `ldsc3`:

```bash
git clone https://github.com/WenbinWu2001/ldsc_py3.git ldsc3
cd ldsc3
```

Create the conda environment and install the package with development extras:

```bash
conda env create -f environment.yml
conda activate ldsc3
python -m pip install -e ".[dev]"
```

**For later runs,** activate the environment before launching LDSC commands or scripts:

```bash
conda activate ldsc3
```

To update the package, enter the repository and pull the latest changes from the `main` branch:

```bash
cd ldsc3
git pull origin main
conda activate ldsc3
```

Because the package is installed in editable mode, source-code updates take effect immediately. 
If `environment.yml` or the package dependencies have changed, also update the environment and reinstall the package:

```bash
conda env update -f environment.yml --prune
python -m pip install -e ".[dev]"
```


### Development version

Clone the `restructure` branch to a directory `ldsc3-dev/`:

```bash
git clone --branch restructure https://github.com/WenbinWu2001/ldsc_py3.git ldsc3-dev
cd ldsc3-dev
```

Create a separate development environment and install the package with
development extras:

```bash
conda env create -f environment.yml -n ldsc3-dev
conda activate ldsc3-dev
python -m pip install -e ".[dev]"
```

**For later runs,** activate the environment before launching LDSC commands or scripts:

```bash
conda activate ldsc3-dev
```

To update the development version, enter the repository and pull the latest changes from the `restructure` branch:

```bash
cd ldsc3-dev
git pull origin restructure
conda activate ldsc3-dev
```

Because the package is installed in editable mode, source-code updates take effect immediately.
If `environment.yml` or the package dependencies have changed, also update the development environment and reinstall the package:

```bash
conda env update -f environment.yml -n ldsc3-dev --prune
python -m pip install -e ".[dev]"
```

The package supports Python 3.11 through 3.13. The base install includes
NumPy, pandas, SciPy, PyArrow, and Matplotlib. Optional extras are split by workflow:
`.[plink]` installs `bitarray` for PLINK-backed LD computation, `.[bed]`
installs `pybedtools` for BED projection, and `.[liftover]` installs
`pyliftover` for chain-file liftover in sumstats munging and cross-build
reference-panel output. Matplotlib is installed by default, but plotting is
never automatic and is imported only by figure-producing paths. `.[dev]`
installs the PLINK, BED, and liftover extras plus pytest. BED-based annotation projection also requires the
external `bedtools` executable, which `environment.yml` installs from bioconda.
For non-conda installs, make sure `bedtools` is available on `PATH` before
running BED annotation workflows.

## CLI

```bash
ldsc --help
```

Equivalent development entrypoint:

```bash
python -m ldsc --help
```

Subcommands:

- `ldsc annotate`
- `ldsc build-gene-ldscore-index`
- `ldsc build-r2-panel`
- `ldsc convert-ldsc2-ldscores`
- `ldsc convert-h2-scale`
- `ldsc ldscore`
- `ldsc munge-sumstats`
- `ldsc h2`
- `ldsc plot`
- `ldsc partitioned-h2`
- `ldsc quantile-h2`
- `ldsc query-r2`
- `ldsc rg`

Pass a canonical h2, partitioned-h2, quantile-h2, or rg result root to
`ldsc plot --result-dir RESULT_DIR`. The command selects one
approved plot from result metadata and writes it below `RESULT_DIR/plots/`;
it does not accept an output directory and never runs automatically. Use
`ldsc convert-h2-scale` to recompute liability-scale h2 from a saved observed-scale
h2 result at one population prevalence or over a sensitivity range. See the
[scientist-facing plotting manual](tutorials/plotting-results.md) and the
[developer module contract](docs/current/plotting-module.md).

Both rg plots annotate saved single-trait observed-scale heritability with jackknife SEs: on the all-pairs heatmap diagonal, or in a partner column and anchor subtitle. Missing or unusable heritability values display `failed`.

Reusable `.annot.gz` shards written by `ldsc annotate` keep the legacy
`CHR/BP/SNP/CM` leading layout. Because annotation `CM` is semantically missing,
the writer serializes it as the explicit `NA` token; annotation-value columns
must contain numeric, non-missing values.

`ldsc munge-sumstats` writes `sumstats.parquet` by default, with optional
legacy `sumstats.gz` output through `--output-format tsv.gz` or
`--output-format both`. Supplying `--trait-name BMI` instead names these files `BMI.parquet` and `BMI.sumstats.gz`; unsafe filename characters are sanitized while the metadata label is preserved. Parquet remains the default format. See [output naming](docs/current/munge-sumstats.md#output-artifacts). Package-written artifacts use canonical `SNP`, `CHR`,
`POS`, `Z`, and `N` columns when possible and always include `CHR`/`POS`. Raw
files may provide chromosome and position columns through common aliases such as
`#CHROM`, `CHROM`, `CHR`, `POS`, or `BP`, or explicitly through `--chr` and
`--pos`; pass the raw input as `--raw-sumstats-file`. Leading `##`
metadata/comment lines are skipped before the real header is parsed. The
`sumstats.parquet` is self-describing: its `snp_identifier`, `genome_build`, and
optional `--trait-name` provenance ride in the Parquet footer, so later
regression commands need only that one file -- no `metadata.json` sidecar is
written. The legacy `sumstats.gz` carries no embedded metadata and is
treated as an rsID lookup artifact at regression time. Footerless Parquet is
rejected rather than guessed.
Python run summaries from `munger.build_run_summary()` report parsed input rows, retained rows, exclusive per-stage `drop_counts`, and the sample-size rule actually used. Counts are collected while parsing, excluding headers and blank lines; their totals reconcile with the retained output. Successful CLI runs print the selected restriction, liftover method, and mapping/drop counts to stdout and record the same summary in `diagnostics/sumstats.log` at every log level. See [munging preparation and accounting](docs/current/munge-sumstats.md#preparation-and-run-accounting) and the [liftover behavior table](docs/current/munge-sumstats.md#liftover-rules).

Detailed coordinate and liftover bookkeeping is written to `sumstats.log`. The default
`snp_identifier` is `chr_pos_allele_aware`, which requires usable `A1/A2`; rerun
with `--snp-identifier chr_pos` to use coordinate identity without
allele-aware matching. The legacy `--no-alleles` escape hatch is no longer
accepted; choose the base `rsid` or `chr_pos` identity mode instead. Munging restricts to the packaged HM3 map by default. `--sumstats-snps-file FILE` replaces HM3 with a custom headered keep-list; `--no-snp-restriction` disables keep-list filtering while retaining ordinary QC. These override flags are mutually exclusive. Restriction files may omit
alleles; allele-free restrictions match by
base key. Allele-bearing restrictions, including packaged HM3, in allele-aware modes match by the
effective allele-aware key. Restriction files are identity-only filters:
duplicate restriction keys collapse to one retained key, and non-identity
columns such as `CM` or `MAF` are ignored. These filters are loaded once before
raw chunk parsing and applied while chunks are still streaming, after canonical
columns and coordinate normalization are available. They keep matching rows only
and do not reorder the output. In `chr_pos`-family modes, `SNP` is a label; matching uses
source-build `CHR` and `POS`, plus the allele set in
`chr_pos_allele_aware`, represented internally as effective coordinate-family
keys. Rows with missing or invalid coordinates are dropped and counted at
coordinate match/map stages. The base `chr_pos` mode uses coordinate identity
without allele-aware matching. To convert coordinates after QC and after SNP
restriction, explicitly choose `--output-genome-build`. Matching builds need no liftover. For different builds, packaged HM3 uses automatic quick liftover; an explicit `--liftover-chain-file` overrides it. Custom-list or unrestricted cross-build runs require a chain file. An unresolved source build stops the run. The former HM3 restriction and quick-liftover switches are removed. Liftover is invalid in `rsid`-family modes because
positions are not the row identity there.

`ldsc build-r2-panel` computes pairwise, bias-adjusted R² and writes chromosome Parquet tables plus required SNP metadata sidecars. It replaces `build-ref-panel` without an alias. The four Parquet columns are `IDX_1`, `IDX_2`, `R2`, and `SIGN_R`; `SIGN_R` is true when Pearson r is nonnegative in the sidecar allele orientation. Panels using the former `SIGN` column are unsupported and must be rebuilt. `query-r2` returns `r2`, nullable `sign_r` (+1/-1 in query allele orientation), `r`, and `status`, with no old-name aliases. See the [R² format contract](docs/current/parquet-r2-format-and-read-pipeline.md#2-parquet-format-specification).

`ldsc build-r2-panel` keeps a separate source-build contract for PLINK input:
provide or infer `--source-genome-build`, and a matching chain file emits the
opposite build. Deliberate reference-universe restriction uses an explicit
`--ref-panel-snps-file`; the builder has no HM3-only restriction or quick-liftover
mode. Chain-file liftover is invalid when the active
`snp_identifier` is in the `rsid` family; omit liftover for source-build-only rsID panels. In
`chr_pos`-family modes,
duplicate source or target coordinate groups are dropped by default
(fixed `drop-all` policy, with no CLI switch), with details in the log and dropped-SNP sidecars under `diagnostics/`.

Artifact-writing workflows write completed per-run logs under their output
directories. During `build-gene-ldscore-index`, the open log temporarily lives
under hidden `.<index-name>.build-state/` so it is never part of an atomic index
replacement; after its handler closes, a successful log moves into the index's
`diagnostics/`. `munge-sumstats` keeps the historical `sumstats.log` name; other
commands use `annotate.log`, `ldscore.log`, `build-r2-panel.log`,
`build-gene-ldscore-index.log`,
`h2.log`, `partitioned-h2.log`, `quantile-h2.log`, `rg.log`, `plot.log`, or
`convert-h2-scale.log`. Concrete single-chromosome
`build-r2-panel` runs use `build-r2-panel.chr<chrom>.log` so parallel
per-chromosome jobs can share an output directory without sharing one log file.
Logs are audit artifacts, so Python result objects and `output_paths` mappings
only list scientific data outputs.

An authorized overwrite that fails leaves `RUN_FAILED.txt` in the applicable
result root, or in the nested `plots/` or liability-scale conversion root for
those utilities. Concrete single-chromosome reference-panel attempts use
`RUN_FAILED.chr<chrom>.txt`. The marker supplements the ordinary failed log;
it does not roll back artifacts or change the existing workflow action order.
A successful retry removes its applicable marker.

Help and `munge-sumstats --infer-only` create no output directories, logs, artifacts, or failure markers, even with `--overwrite` or a failed inference check. Existing output contents, including failure markers, remain untouched.

`ldsc ldscore` supports ordinary unpartitioned LD-score generation without
baseline annotations:

```bash
ldsc ldscore --output-dir ldscores --plink-prefix panel --ld-wind-snps 10
ldsc h2 --sumstats-file trait/sumstats.parquet --ldscore-dir ldscores --output-dir h2_out
```

When no baseline and no query annotations are supplied, the workflow writes a
synthetic all-ones baseline column named exactly `base` in
`ldscore.baseline.parquet`.
Query annotation inputs still require explicit `--baseline-annot-sources`.
Regression diagnostics report the SNP population and jackknife block count actually used for fitting. The h2 LD-score regression plot summarizes that same fitted population; see [regression configuration](docs/current/regression-configuration.md).

For large query scans, `ldsc partitioned-h2 --continue-on-query-error` skips queries whose model preparation, regression, or result calculation raises an exception and publishes successful fits. Every attempted query is recorded in `diagnostics/query_status.tsv`; failed queries have no scientific results. The log includes query/stage tracebacks and available jackknife block diagnostics. The default remains strict: any query failure prevents publication. Shared input/output failures and scans with no successful queries still fail. See [query failure handling](docs/current/partitioned-h2-results.md#query-failures-and-continuation).

Use this synthetic `base` directory for `ldsc h2` or `ldsc rg`. A baseline-only
directory is also accepted by `ldsc partitioned-h2` in its functional-category
regime, although a single all-ones `base` column is a degenerate one-category
fit rather than a meaningful partitioned analysis.

Native LD-score runs classify fitted columns as binary or quantitative for interpretation only. Partitioned-h2 preserves numerical legacy enrichment summaries for quantitative columns but logs that those weighted values do not have the ordinary binary-category interpretation; coefficients remain interpretable. `ldsc quantile-h2` consumes one baseline-only or per-query fitted model plus resupplied annotation/reference sources and reports joint-model heritability by target quantile and standardized `tau_star`. Resupplied sources are checked for alignment and agreement with stored counts, annotation sums, and available overlap cross-products; passing these checks does not establish the original annotation value at every SNP. See [the technical contract](docs/current/continuous-annotation-quantile-h2.md) and [the concise workflow](docs/wiki/continuous-annotation-partitioned-ldsc.md).
LD-score generation runs queries in sequential batches (`--query-batch-size`, default 1000). A single batch writes `ldscore.query.parquet`; multiple batches write `ldscore.query.batch00001.parquet` and subsequent ordinals. Every file spans the computed chromosomes with one row group per chromosome. Root `metadata.json.query_batches` records filenames, ordered query columns, and row groups; one shared baseline and overlap are written. Older LD-score directories without this manifest must be regenerated. Direct mode repeats reference work for each batch. Direct and indexed modes support `--threads`, capped at chromosome count.

Writing Python workflows return `LDScoreSource`; use `result.read_queries(["pathway_a", "pathway_b"])` for explicit selections across files. These reads have no cache or width limit; the caller owns their RAM cost. For a small single-batch calculation with zero filesystem writes, prepare `AnnotationBundle.from_frames(...)` and call `LDScoreCalculator.run(..., output_config=None)` with a prepared reference adapter. An output directory is required for multiple batches. See the [memory contracts](docs/current/annotation-memory-design.md).

Gene-list LD scores may also use an explicitly installed exact profile:

```bash
ldsc ldscore \
  --gene-ldscore-index-dir indexes/1000G_EUR_baseline_100kb \
  --query-annot-gene-list-sources gene_lists/immune.txt \
  --output-dir gene_ldscores
```

Build complete indexes offline with `ldsc build-gene-ldscore-index`. Construction requires an explicit one-based `--gene-coordinate-file`, defaults to unpadded gene intervals (`--padding-bp 0`), and requires both `--genome-build hg19` and an explicit `--snp-identifier rsid|chr_pos`; there is no packaged gene catalog, identity/build default, inference, hg38, liftover, or allele-aware index mode. Indexed mode is explicit and fail-closed: it validates and inherits the complete index's identity/build and stored padding, rejects live identity/build or baseline/reference/window overrides, and requires callers to remove `--padding-bp` entirely rather than pass zero. It writes the same self-contained canonical LD-score directory. During the offline build, each completed chromosome is atomically persisted in a hidden run-specific stage, allowing its large in-memory payload to be released while later chromosomes continue. A `Finished chromosome N` log line reports that private durability boundary; the public destination remains absent, empty, or at its prior complete version until every chromosome and shared metadata pass reload validation. These private shards are not resumable checkpoints and never support chromosome append or incremental index updates. See
[the exact gene-index guide](docs/current/gene-ldscore-index.md) and its
[mathematical algorithm](docs/current/gene-ldscore-index-mathematics.md).
Task-oriented walkthroughs cover [building the index](docs/wiki/utility-functionalities/build-gene-ldscore-index.md)
and [using it for gene-list LD scores](docs/wiki/main-functionalities/ldscore-from-gene-list.md).

Gene-index metadata retains the semantic `index_id` and one ordered published-row fingerprint. New indexes omit the redundant effective-row fingerprint and require an updated reader; otherwise-valid existing indexes remain readable. See the [fingerprint and reader contract](docs/current/gene-ldscore-index.md#artifact-and-identity-contract).

Gene-list resolution is strict by default in `annotate`, `ldscore`, and `quantile-h2`. Add the value-free flag `--allow-unresolved-genes` to deliberately use the resolved subset with omission diagnostics. It replaces `--gene-list-resolution-policy resolved-only`; remove `--gene-list-resolution-policy strict` from existing commands. Python policy arguments and metadata retain `strict` / `resolved-only`.

Direct gene-list, BED, and prebuilt-query runs require matching validated baseline/reference chromosome sets. `@` requires autosomes 1–22; globs select actual matches, and file contents determine scope. Every selected focal/control gene must lie within scope under both resolution policies: incomplete coverage fails the whole batch without truncation or automatic pathway skipping. Chromosome-only runs are supported when inputs align. The effective scope is logged and written to `diagnostics/chromosome_scope.json`; input failures, pathway counts, and affected genes appear in diagnostics. See [coverage and repair rules](docs/current/gene-list-diagnostics-and-repair.md#chromosome-scope-and-pathway-coverage), including the trade-off that globs can hide a consistently missing chromosome.

## LDSC2 compatibility boundary

Backward compatibility is deliberately asymmetric and limited to artifacts
that are costly or impractical for users to recreate:

- Regression accepts genuine LDSC2 `.sumstats` and `.sumstats.gz` text files
  automatically. Their `SNP` values are treated as rsID lookup keys and are
  projected onto the canonical LDSC3 LD-score panel. `A1` and `A2` are required;
  allele orientation is validated and `Z` is negated when a swap is needed.
  `FRQ` is optional, is never imputed from the panel, and does not control SNP
  retention. Footerless Parquet is not treated as a legacy artifact.
- LDSC2 LD-score fragments are never accepted directly by regression. Run
  `ldsc convert-ldsc2-ldscores` explicitly with a complete reference directory
  and a complete regression-weight directory. The converter supports a
  one-column unpartitioned suite or a complete baseline partitioned suite; the
  latter also requires a frequency directory. Query/cell-type suites and thin
  annotations are intentionally unsupported.
- Reference and weight suites both use `<prefix><chrom>.l2.ldscore(.gz)` filenames, such as `weights.1.l2.ldscore.gz`, across chromosomes 1-22. The `.w.l2.ldscore(.gz)` suffix is unsupported; discovery explains this spelling when conversion fails and records affected files in its diagnostics. Private `_kernel` code no longer writes legacy sumstats, LD-score, count, or annotation-group files; optional `.sumstats.gz` export belongs to the public munging workflow.
- Converted suites remain allele-unaware (`rsid` by default, or `chr_pos`).
  `.l2.M_5_50` is required and retains the fixed strict LDSC2 common-frequency
  rule. There is no converter threshold flag. Missing `.l2.M` is tolerated only
  under the documented count policy; requesting unavailable all-SNP counts
  later is an error.
- Annotation inputs may omit both `A1` and `A2`, including in allele-aware workflows; annotation membership then matches by base SNP identity, and LD-score calculation uses reference-panel alleles for allele-aware identity. A file containing only one of the two allele columns is rejected. This annotation rule does not relax the allele requirement for legacy sumstats.
- LDSC3 does not promise that its outputs can be fed back into LDSC2. The
  `BP` header written in text `.annot.gz` files is a narrow interoperability
  convenience, and LDSC3 accepts either `BP` or `POS` when reading them.

See the complete policies for
[legacy sumstats](docs/current/legacy-sumstats-compatibility.md) and
[legacy LD-score conversion](docs/current/legacy-ldscore-conversion.md).

For commands, standard filenames, complete suite examples, and practical naming repairs, follow the [legacy LD-score conversion user guide](tutorials/convert-legacy-ldscores.md). The [conversion utility wiki](docs/wiki/utility-functionalities/convert-ldsc2-ldscores.md) provides the concise standard naming and flag reference.

## Python API

```python
from ldsc import (
    AnnotationBuilder,
    ChrPosBuildInference,
    GeneLDScoreIndexBuildConfig,
    LDScoreCalculator,
    ReferencePanelBuilder,
    RegressionRunner,
    SumstatsMunger,
    build_gene_ldscore_index,
    infer_chr_pos_build,
    load_gene_ldscore_index,
    load_sumstats,
    resolve_chr_pos_table,
)
```

`LDScoreCalculator` delegates reference filtering, annotation alignment, and reader lifetime to `RefPanel.prepare_chromosome()`. PLINK and parquet calculations use the same prepared-state contract, also used by gene-index construction and numerical tests. See the [reference preparation boundary](docs/current/architecture.md) for contributor guidance.

Genome-build inference for `chr_pos` tables is public through the Python API:
use `infer_chr_pos_build()` when you only need the `hg19`/`hg38` decision, and
use `resolve_chr_pos_table()` when you also want 0-based inputs converted to
canonical 1-based coordinates. The CLI exposes this behavior inside workflows
with `--genome-build auto`, including during `munge-sumstats`; there is no
standalone inference subcommand.

## Input Path Tokens

Public workflow APIs accept normalized string tokens for inputs:

- exact paths
- Python glob patterns such as `annotations/*.annot.gz`
- explicit chromosome suites using `@`, for example `baseline.@`
- PLINK exact prefixes/member paths, plain suite stems such as `1000G.EUR.QC.`, globs, and `@` patterns

PLINK inputs share one discovery and validation implementation across `ldscore`, `build-gene-ldscore-index`, and `build-r2-panel`. Chromosomes come from BIM contents; a numeric suffix such as `.22` is preserved as part of the prefix. Selected trios must be complete, and different trios cannot claim the same chromosome. See [PLINK resolution](docs/current/path-specification.md#plink-prefix-resolution).

Inputs are resolved before the internal kernel runs. Public outputs use fixed
filenames inside the selected `output_dir`.

## Input validation stages

Input validation is separate from output-overwrite protection. Direct LD scores validate reference inputs before annotation staging; gene-index and R2 builders check companion and auxiliary declarations together. Index consumers inspect every chromosome's component paths and small headers before loading operators. Independent defects at each gate are reported together with machine-readable repair diagnostics. Alignment, identity cleanup, and post-filter support are later content-dependent checks; a successful path gate does not guarantee a valid scientific dataset. See [staged validation and progress](docs/current/workflow-logging.md#staged-input-validation) and [repair guidance](docs/troubleshooting.md#staged-input-validation).

## Output Collision Policy

For workflows that accept `--output-dir` or `output_dir`, the value is always a directory:

- missing output directories are created and logged at INFO
- existing directories are reused
- known output files fail the run before writing starts
- reruns that intentionally replace known files must pass `--overwrite` on the
  CLI or `overwrite=True` in Python
- per-run log files are preflighted with the scientific outputs, so a collision
  fails before a new log is opened

Gene-index retries accept an output directory containing only `RUN_FAILED.txt` and recognized legacy build diagnostics. The marker remains until successful publication; unrelated contents still block replacement. See [gene-index retry guidance](docs/troubleshooting.md#build-gene-ldscore-index-output-directory-is-nonempty-but-invalid).

Failure markers, logs, audits, scratch, and result files use the same normalized destination, including `~` and environment-variable expansion. A root failure marker does not count as a scientific output collision. See the [output-directory retry audit](docs/audits/2026-09-14-output-directory-retries.md) for coverage across all 13 commands.

The overwrite flag applies only to the fixed files owned by that workflow. It
does not remove unrelated files and never cleans a whole directory.
The derived `plot` and `convert-h2-scale` commands instead use the fixed nested destinations described above.

Directory writers share their artifact declarations with workflow preflight. Conditional files, including chromosome drop reports, are checked before the workflow log opens; successful overwrites reconcile the files actually produced. See [Output-Family Preflight](docs/current/workflow-logging.md#output-family-preflight) for ownership and cleanup rules.

At the default INFO level, workflow logs mark validation, staging, computation, and publication phases with elapsed time, plus periodic counts and the current object at long-running chunk checkpoints. Annotation preparation retains its input-reading, identity-checking, and retained-count completion milestones. CM/MAF notices appear once per input file read, and intentional gene exclusions are summarized on one line per gene set. See [annotation preparation logging](docs/current/workflow-logging.md#annotation-preparation) for examples and interpretation.

## Verification

Use the editable development installation described above and run these commands from the repository root. Tests do not inject `src` into Python paths. Run pytest and unittest sequentially because pybedtools cleanup can affect another runner’s temporary files. Goldens are immutable expectations; restore missing NPZs from version control rather than regenerating them during tests.

```bash
pytest
```

The suite also remains compatible with the standard-library runner while the
project transitions fully to pytest:

```bash
python -m unittest discover -s tests -p 'test*.py' -v
```
