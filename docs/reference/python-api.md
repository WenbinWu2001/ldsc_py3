# Python API reference

Last updated on: 2026-10-27

The supported Python interface is the set of names exported by the `ldsc` package (`ldsc.__all__`, 84 names on 2026-10-27).
Modules whose names start with an underscore, including `ldsc._kernel`, are private even when a name defined there is re-exported.
This page lists every exported name with the first paragraph of its docstring, taken from the installed package; `help(ldsc.<name>)` shows the full docstring, and a future documentation build will generate this page.
Task-oriented examples are in [Python usage](../user-guide/python-usage.md), and the command-line equivalents in the [command-line reference](cli.md).

## Workflow classes

| Name | Summary |
| --- | --- |
| `AnnotationBuilder` | Prepare chromosome-backed annotations under explicit resource ownership. |
| `LDScoreCalculator` | Orchestrate chromosome-wise LD-score calculation. |
| `LegacyLDScoreConverter` | Convert a reusable LDSC2 reference/weight suite into LDSC3 format. |
| `ReferencePanelBuilder` | Build standard parquet reference-panel artifacts from PLINK input. |
| `RegressionRunner` | Assemble regression datasets and fit in children with 1 BLAS thread each. |
| `SumstatsMunger` | Run legacy-compatible munging through a typed workflow interface. |

## Functions

| Name | Summary |
| --- | --- |
| `assign_legacy_quantiles` | Assign finite target values with the LDSC2 continuous-annotation rule. |
| `build_gene_ldscore_index` | Build and transactionally publish one complete gene LD-score index. |
| `compute_quantile_h2` | Project one fitted joint model onto target-defined quantiles. |
| `compute_standardized_coefficients` | Compute paper-compatible fixed-scale standardized coefficients. |
| `convert_h2_scale` | Convert a saved observed-scale SNP-heritability estimate. |
| `convert_ldsc2_ldscores` | Convert one selected LDSC2 LD-score suite to a canonical LDSC3 directory. |
| `infer_chr_pos_build` | Infer the genome build and coordinate basis of a ``CHR``/``POS`` table. |
| `infer_raw_sumstats` | Infer raw summary-statistics format and minimal parser hints. |
| `load_fitted_partitioned_model` | Load one baseline-only or per-query LDSC3 partitioned-h2 fitted model. |
| `load_gene_ldscore_index` | Load one complete public autosomes 1--22 gene LD-score index. |
| `load_hm3_curated_map` | Load the packaged curated dual-build HM3 map. |
| `load_ldscore_from_dir` | Open a canonical LD-score directory with query values read on demand. |
| `load_sumstats` | Load one curated LDSC sumstats artifact into a ``SumstatsTable``. |
| `plot_result` | Create the approved plot for a canonical LDSC result. |
| `query_r2` | Open a panel, query ``pairs`` once, and return the annotated DataFrame. |
| `resolve_chr_pos_table` | Infer and normalize the genome build for a ``CHR``/``POS`` table. |
| `resolve_genome_build` | Resolve a genome-build hint to a concrete build or ``None``. |
| `run_annotate` | Write BED or focal gene-list annotations and return a live shard handle. |
| `run_build_ref_panel` | Run the reference-panel builder with CLI-equivalent keyword arguments. |
| `run_h2` | Estimate SNP heritability and save the ``h2`` result directory. |
| `run_ldscore` | Run LD-score calculation from Python using public CLI-style names. |
| `run_partitioned_h2` | Estimate partitioned heritability and save the ``partitioned-h2`` result directory. |
| `run_quantile_h2` | Project a saved partitioned-h2 model onto quantiles and save the result directory. |
| `run_rg` | Estimate genetic correlation and save the ``rg`` result directory. |
| `unbiased_r2_to_pearson_r` | Convert adjusted (unbiased) R² to a signed Pearson correlation r. |

## Configuration

| Name | Summary |
| --- | --- |
| `AnnotationBuildConfig` | Source and projection configuration for chromosome-backed annotations. |
| `GeneLDScoreIndexBuildConfig` | Scientific and resource settings for exact gene-index construction. |
| `GlobalConfig` | Shared configuration used across the refactored workflows. |
| `H2OutputConfig` | Directory-oriented output config for unpartitioned h2 summaries. |
| `LDScoreConfig` | Configuration for chromosome-wise LD-score calculation. |
| `LDScoreOutputConfig` | Directory-oriented output config for canonical LD-score results. |
| `MungeConfig` | Configuration for legacy-compatible summary-statistics munging. |
| `PartitionedH2OutputConfig` | Directory-oriented output config for partitioned-h2 regression summaries. |
| `QuantileH2OutputConfig` | Output directory and replacement policy for ``quantile-h2`` results. |
| `RefPanelConfig` | Configuration for choosing and parameterizing a reference-panel backend. |
| `ReferencePanelBuildConfig` | Configuration for building standard parquet reference panels from PLINK. |
| `RegressionConfig` | Configuration for heritability and genetic-correlation regressions. |
| `RgOutputConfig` | Directory-oriented output config for genetic-correlation summaries. |

## Global configuration

| Name | Summary |
| --- | --- |
| `get_global_config` | Return the process-wide default configuration used by Python workflows. |
| `reset_global_config` | Restore the process-wide default configuration and return it. |
| `set_global_config` | Replace the process-wide default configuration used by Python workflows. |
| `validate_config_compatibility` | Raise when ``a`` and ``b`` disagree on critical LDSC workflow settings. |

## Results and data objects

| Name | Summary |
| --- | --- |
| `AnnotationBundle` | Complete annotation dataset represented by on-demand chromosome shards. |
| `ChrPosBuildInference` | Resolved build and coordinate-basis decision for one ``chr_pos`` table. |
| `ChromLDScoreResult` | Normalized per-chromosome LD-score output in split baseline/query form. |
| `H2ScaleConversionArtifact` | Files written by one liability-scale conversion. |
| `LDScoreResult` | Materialized LD scores for one cross-chromosome calculation batch. |
| `LDScoreSource` | Saved LD scores with resident baseline values and selective query reads. |
| `LegacyLDScoreConversionResult` | Summary of one successful LDSC2 LD-score conversion. |
| `MungeRunSummary` | Compact summary of one munging run. |
| `PartitionedH2BatchResult` | Compact batch summary and valid persistent artifacts for every fitted model. |
| `PlotArtifact` | One saved LDSC figure and its live Matplotlib objects. |
| `PreparedRegressionInputs` | Shared trait alignment before model-dependent column/count filtering. |
| `QuantileAssignment` | Legacy-compatible quantile membership and realized boundaries. |
| `QuantileH2Result` | Compact result family returned by `run_quantile_h2()`. |
| `RawSumstatsInference` | Header-level inference report for one raw summary-statistics file. |
| `RefPanelLoader` | Instantiate the backend requested by a `RefPanelConfig`. |
| `ReferencePanelBuildResult` | Summary of files produced by one reference-panel build run. |
| `RegressionDataset` | Merged regression-ready dataset built from sumstats and LD-score tables. |
| `RgResultFamily` | Result family returned by multi-trait genetic-correlation workflows. |
| `SumstatsTable` | In-memory LDSC-ready summary-statistics table. |

## Reference panels

| Name | Summary |
| --- | --- |
| `ParquetR2RefPanel` | Canonical index-format parquet R2 reference-panel adapter. |
| `PlinkRefPanel` | PLINK-backed reference-panel adapter. |
| `R2Panel` | Handle over a build-r2-panel index-format R² panel for pair queries. |
| `RefPanel` | Abstract chromosome-scoped reference-panel interface. |

## Output writers

| Name | Summary |
| --- | --- |
| `H2DirectoryWriter` | Write the unpartitioned h2 summary, exact fitted-bin diagnostic, and metadata. |
| `LDScoreDirectoryWriter` | Write canonical LD-score result directories. |
| `PartitionedH2DirectoryWriter` | Write aggregate and optional per-query partitioned-h2 outputs. |
| `QuantileH2DirectoryWriter` | Write post-fit quantile summaries, coefficients, and SNP diagnostics. |
| `RgDirectoryWriter` | Write concise, full, per-trait, and optional per-pair rg outputs. |

## Errors

| Name | Summary |
| --- | --- |
| `ConfigMismatchError` | Raised when two result objects carry incompatible config snapshots. |
| `LDSCConfigError` | Raised when configuration values or provenance are incompatible. |
| `LDSCDependencyError` | Raised when an optional dependency is required for the requested workflow. |
| `LDSCError` | Base class for package-specific failures. |
| `LDSCInputError` | Raised when input paths, file schemas, or file contents are invalid. |
| `LDSCInternalError` | Raised when the package reaches an unexpected internal state. |
| `LDSCUsageError` | Raised when a public command or API call uses an unsupported option mix. |
| `LDSCUserError` | Raised when user input, CLI arguments, or configuration are invalid. |

## Configuration fields

Field names of each configuration dataclass; types and defaults are in each class docstring.

| Class | Fields |
| --- | --- |
| `AnnotationBuildConfig` | `baseline_annot_sources`, `query_annot_sources`, `query_annot_bed_sources`, `query_annot_gene_list_sources`, `gene_coordinate_file`, `control_gene_list_file`, `gene_list_resolution_policy`, `gene_exclude_regions`, `padding_bp` |
| `GeneLDScoreIndexBuildConfig` | `baseline_annot_sources`, `plink_prefix`, `output_dir`, `gene_coordinate_file`, `genome_build`, `snp_identifier`, `padding_bp`, `gene_exclude_regions`, `ld_wind_cm`, `ld_wind_kb`, `ld_wind_snps`, `maf_min`, `common_maf_min`, `keep_indivs_file`, `regr_snps_file`, `regr_snps_exclude_regions`, `snp_batch_size`, `atom_batch_size`, `workers` |
| `GlobalConfig` | `snp_identifier`, `genome_build`, `log_level` |
| `H2OutputConfig` | `output_dir`, `overwrite` |
| `LDScoreConfig` | `ld_wind_snps`, `ld_wind_kb`, `ld_wind_cm`, `regr_snps_file`, `snp_batch_size`, `query_batch_size`, `common_maf_min`, `whole_chromosome_ok`, `export_ref_metadata`, `workers` |
| `LDScoreOutputConfig` | `output_dir`, `overwrite`, `parquet_compression` |
| `MungeConfig` | `output_dir`, `raw_sumstats_file`, `compression`, `output_format`, `trait_name`, `column_hints`, `N`, `N_cas`, `N_con`, `info_min`, `maf_min`, `n_min`, `nstudy_min`, `chunk_size`, `sumstats_snps_file`, `no_snp_restriction`, `source_genome_build`, `output_genome_build`, `liftover_chain_file`, `signed_sumstats_spec`, `ignore_columns`, `info_list_columns`, `sumstats_format`, `a1_inc`, `overwrite` |
| `PartitionedH2OutputConfig` | `output_dir`, `overwrite`, `write_per_query_results` |
| `QuantileH2OutputConfig` | `output_dir`, `overwrite` |
| `RefPanelConfig` | `backend`, `plink_prefix`, `r2_dir`, `chromosomes`, `maf_min`, `keep_indivs_file`, `ref_panel_snps_file`, `genetic_map_hg19_sources`, `genetic_map_hg38_sources` |
| `ReferencePanelBuildConfig` | `plink_prefix`, `output_dir`, `source_genome_build`, `genetic_map_hg19_sources`, `genetic_map_hg38_sources`, `liftover_chain_hg19_to_hg38_file`, `liftover_chain_hg38_to_hg19_file`, `ld_wind_snps`, `ld_wind_kb`, `ld_wind_cm`, `maf_min`, `ref_panel_snps_file`, `keep_indivs_file`, `snp_batch_size`, `min_r2`, `overwrite` |
| `RegressionConfig` | `n_blocks`, `use_common_counts`, `use_intercept`, `intercept_h2`, `intercept_gencov`, `two_step_cutoff`, `chisq_max`, `samp_prev`, `pop_prev`, `allow_identity_downgrade` |
| `RgOutputConfig` | `output_dir`, `overwrite`, `write_per_pair_detail` |
