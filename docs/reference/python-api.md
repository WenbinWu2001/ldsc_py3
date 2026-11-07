# Python API reference

Last updated on: 2026-11-07

The supported Python interface is the set of names exported by the `ldsc` package (`ldsc.__all__`).
Modules whose names start with an underscore, including the domain packages such as `ldsc._ldscore`, are private even when a name defined there is re-exported.
The documentation build generates one page per exported name from its docstring, and `help(ldsc.<name>)` shows the same docstring for the installed version.
`tests/package/test_docs_api_reference.py` checks that the lists below name every exported name exactly once.
Task-oriented examples are in [Python usage](../user-guide/python-usage.md), and the command-line equivalents in the [command-line reference](cli.md).

```{eval-rst}
.. currentmodule:: ldsc
```

## Workflow classes

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   AnnotationBuilder
   LDScoreCalculator
   R2MatrixBuilder
   RegressionRunner
   SumstatsMunger
```

## Functions

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   build_gene_ldscore_index
   compute_quantile_h2
   compute_standardized_coefficients
   convert_h2_scale
   convert_legacy_ldscores
   infer_chr_pos_build
   infer_raw_sumstats
   load_fitted_partitioned_model
   load_gene_ldscore_index
   load_hm3_curated_map
   load_ldscore_from_dir
   load_sumstats
   plot_result
   query_r2_matrix
   resolve_chr_pos_table
   resolve_genome_build
   run_annotate
   run_build_r2_matrix
   run_h2
   run_ldscore
   run_partitioned_h2
   run_quantile_h2
   run_rg
   unbiased_r2_to_pearson_r
```

## Configuration

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   AnnotationBuildConfig
   GeneLDScoreIndexBuildConfig
   GlobalConfig
   H2OutputConfig
   LDScoreConfig
   LDScoreOutputConfig
   MungeConfig
   PartitionedH2OutputConfig
   QuantileH2OutputConfig
   RefPanelConfig
   R2MatrixBuildConfig
   RegressionConfig
   RgOutputConfig
```

## Global configuration

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   get_global_config
   reset_global_config
   set_global_config
```

## Results and data objects

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   AnnotationBundle
   AnnotationResult
   ChrPosBuildInference
   H2Result
   H2ScaleConversionArtifact
   LDScoreResult
   LDScoreSource
   LegacyLDScoreConversionResult
   MungeRunSummary
   PartitionedH2BatchResult
   PlotArtifact
   QuantileH2Result
   RawSumstatsInference
   RefPanelLoader
   R2MatrixBuildResult
   RegressionDataset
   RegressionFit
   RgResultFamily
   SumstatsTable
```

## Reference panels

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   R2MatrixRefPanel
   PlinkRefPanel
   R2Matrix
   RefPanel
```

## Output writers

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   H2DirectoryWriter
   LDScoreDirectoryWriter
   PartitionedH2DirectoryWriter
   QuantileH2DirectoryWriter
   RgDirectoryWriter
```

## Errors

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   ConfigMismatchError
   LDSCConfigError
   LDSCDependencyError
   LDSCError
   LDSCInputError
   LDSCInternalError
   LDSCUsageError
   LDSCUserError
```
