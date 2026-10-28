"""``ldsc`` package public surface.

Core functionality:
    Re-export the LDSC3 workflow classes, the R2 matrix builder,
    reference-panel loaders, configuration dataclasses, and convenience
    helpers from a single import location.

Overview
--------
This module is the top-level public API for LDSC3. Import from
``ldsc`` when you want the supported package surface; the private domain
packages (``ldsc._core``, ``ldsc._ldscore`` and the others) hold the
implementation. The exports here mirror the main user
workflows: annotation building, R2 matrix building, exact gene
LD-score index construction, LD-score calculation, summary-statistics munging,
regression, output writing, genome-build inference for ``chr_pos`` inputs, and
lazy access to result plotting and liability-scale post-processing.

Design Notes
------------
- Only stable public objects should be re-exported here.
- Every module whose package name starts with an underscore is private and
  not part of the user-facing API, even when a name it defines is exported here.
- Plotting exports are resolved lazily so importing ``ldsc`` does not import
  Matplotlib until a figure-producing path is used.

Example
-------
>>> from ldsc import GlobalConfig, R2MatrixBuilder
>>> GlobalConfig().snp_identifier
'chr_pos_allele_aware'
>>> isinstance(R2MatrixBuilder(), R2MatrixBuilder)
True
"""

from importlib import import_module

from ._annotation.builder import (
    AnnotationBuilder,
    AnnotationBundle,
    run_annotate,
)
from .config import (
    AnnotationBuildConfig,
    ConfigMismatchError,
    GlobalConfig,
    GeneLDScoreIndexBuildConfig,
    LDScoreConfig,
    MungeConfig,
    R2MatrixBuildConfig,
    RefPanelConfig,
    RegressionConfig,
    get_global_config,
    reset_global_config,
    set_global_config,
    validate_config_compatibility,
)
from .errors import (
    LDSCConfigError,
    LDSCDependencyError,
    LDSCError,
    LDSCInputError,
    LDSCInternalError,
    LDSCUsageError,
    LDSCUserError,
)
from ._genome.build_inference import (
    ChrPosBuildInference,
    infer_chr_pos_build,
    resolve_genome_build,
    resolve_chr_pos_table,
)
from ._genome.hm3 import load_hm3_curated_map
from ._ldscore.calculator import ChromLDScoreResult, LDScoreCalculator, LDScoreResult, run_ldscore
from ._ldscore.gene_index.build import build_gene_ldscore_index
from ._ldscore.gene_index.scoring import load_gene_ldscore_index
from ._regression.outputs import (
    H2DirectoryWriter,
    H2OutputConfig,
    PartitionedH2DirectoryWriter,
    PartitionedH2OutputConfig,
    QuantileH2DirectoryWriter,
    QuantileH2OutputConfig,
    RgDirectoryWriter,
    RgOutputConfig,
)
from ._ldscore.outputs import LDScoreDirectoryWriter, LDScoreOutputConfig
from ._refpanel.builder import R2MatrixBuildResult, R2MatrixBuilder, run_build_r2_matrix
from ._refpanel.adapters import R2MatrixRefPanel, PlinkRefPanel, RefPanel, RefPanelLoader
from ._refpanel.query import R2Matrix, query_r2_matrix, unbiased_r2_to_pearson_r

_LAZY_EXPORTS = {
    "MungeRunSummary": ("._sumstats.munger", "MungeRunSummary"),
    "RawSumstatsInference": ("._sumstats.munger", "RawSumstatsInference"),
    "SumstatsMunger": ("._sumstats.munger", "SumstatsMunger"),
    "SumstatsTable": ("._sumstats.munger", "SumstatsTable"),
    "infer_raw_sumstats": ("._sumstats.munger", "infer_raw_sumstats"),
    "load_sumstats": ("._sumstats.munger", "load_sumstats"),
    "load_ldscore_from_dir": ("._ldscore.source", "load_ldscore_from_dir"),
    "LDScoreSource": ("._ldscore.source", "LDScoreSource"),
    "LegacyLDScoreConverter": ("._ldscore.legacy_converter", "LegacyLDScoreConverter"),
    "LegacyLDScoreConversionResult": ("._ldscore.legacy_converter", "LegacyLDScoreConversionResult"),
    "convert_legacy_ldscores": ("._ldscore.legacy_converter", "convert_legacy_ldscores"),
    "RegressionDataset": ("._regression.runner", "RegressionDataset"),
    "PreparedRegressionInputs": ("._regression.runner", "PreparedRegressionInputs"),
    "PartitionedH2BatchResult": ("._regression.runner", "PartitionedH2BatchResult"),
    "RegressionRunner": ("._regression.runner", "RegressionRunner"),
    "RgResultFamily": ("._regression.runner", "RgResultFamily"),
    "run_h2": ("._regression.runner", "run_h2"),
    "run_partitioned_h2": ("._regression.runner", "run_partitioned_h2"),
    "run_rg": ("._regression.runner", "run_rg"),
    "QuantileAssignment": ("._regression.quantile_h2", "QuantileAssignment"),
    "QuantileH2Result": ("._regression.quantile_h2", "QuantileH2Result"),
    "assign_legacy_quantiles": ("._regression.quantile_h2", "assign_legacy_quantiles"),
    "compute_quantile_h2": ("._regression.quantile_h2", "compute_quantile_h2"),
    "compute_standardized_coefficients": ("._regression.quantile_h2", "compute_standardized_coefficients"),
    "load_fitted_partitioned_model": ("._regression.quantile_h2", "load_fitted_partitioned_model"),
    "run_quantile_h2": ("._regression.quantile_h2", "run_quantile_h2"),
    "H2ScaleConversionArtifact": ("._regression.h2_scale", "H2ScaleConversionArtifact"),
    "convert_h2_scale": ("._regression.h2_scale", "convert_h2_scale"),
    "PlotArtifact": ("._plotting.dispatch", "PlotArtifact"),
    "plot_result": ("._plotting.dispatch", "plot_result"),
}


def __getattr__(name: str):
    """Lazily import optional heavy public workflows on first attribute access."""
    if name not in _LAZY_EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr_name = _LAZY_EXPORTS[name]
    module = import_module(module_name, __name__)
    value = getattr(module, attr_name)
    globals()[name] = value
    return value

__all__ = [
    "AnnotationBuildConfig",
    "AnnotationBuilder",
    "AnnotationBundle",
    "ChromLDScoreResult",
    "ChrPosBuildInference",
    "ConfigMismatchError",
    "GlobalConfig",
    "GeneLDScoreIndexBuildConfig",
    "get_global_config",
    "infer_chr_pos_build",
    "H2DirectoryWriter",
    "H2OutputConfig",
    "H2ScaleConversionArtifact",
    "LDScoreCalculator",
    "LDScoreConfig",
    "LDScoreDirectoryWriter",
    "LDScoreOutputConfig",
    "LDScoreResult",
    "LDScoreSource",
    "LegacyLDScoreConverter",
    "LegacyLDScoreConversionResult",
    "LDSCConfigError",
    "LDSCDependencyError",
    "LDSCError",
    "LDSCInputError",
    "LDSCInternalError",
    "LDSCUsageError",
    "LDSCUserError",
    "load_ldscore_from_dir",
    "load_hm3_curated_map",
    "load_gene_ldscore_index",
    "MungeConfig",
    "MungeRunSummary",
    "RawSumstatsInference",
    "R2MatrixRefPanel",
    "R2Matrix",
    "query_r2_matrix",
    "unbiased_r2_to_pearson_r",
    "PartitionedH2DirectoryWriter",
    "PartitionedH2OutputConfig",
    "QuantileAssignment",
    "QuantileH2DirectoryWriter",
    "QuantileH2OutputConfig",
    "QuantileH2Result",
    "PlinkRefPanel",
    "PlotArtifact",
    "R2MatrixBuildConfig",
    "R2MatrixBuildResult",
    "R2MatrixBuilder",
    "RefPanel",
    "RefPanelConfig",
    "RefPanelLoader",
    "RegressionConfig",
    "RegressionDataset",
    "PreparedRegressionInputs",
    "PartitionedH2BatchResult",
    "RegressionRunner",
    "RgDirectoryWriter",
    "RgOutputConfig",
    "RgResultFamily",
    "SumstatsMunger",
    "SumstatsTable",
    "load_sumstats",
    "plot_result",
    "convert_legacy_ldscores",
    "convert_h2_scale",
    "assign_legacy_quantiles",
    "compute_quantile_h2",
    "compute_standardized_coefficients",
    "infer_raw_sumstats",
    "run_build_r2_matrix",
    "build_gene_ldscore_index",
    "run_annotate",
    "run_ldscore",
    "run_h2",
    "run_partitioned_h2",
    "run_rg",
    "run_quantile_h2",
    "load_fitted_partitioned_model",
    "resolve_chr_pos_table",
    "reset_global_config",
    "resolve_genome_build",
    "set_global_config",
    "validate_config_compatibility",
]
