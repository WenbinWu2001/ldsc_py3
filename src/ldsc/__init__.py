"""ldsc package public surface.

Core functionality:
    Re-export the refactored LDSC workflow classes, reference-panel builders
    and loaders, configuration dataclasses, and convenience helpers from a
    single import location.

Overview
--------
This module is the top-level public API for the refactored package. Import from
``ldsc`` when you want the supported package surface rather than the internal
``ldsc._kernel`` implementation modules. The exports here mirror the main user
workflows: annotation building, parquet reference-panel building, exact gene
LD-score index construction, LD-score calculation, summary-statistics munging,
regression, output writing, genome-build inference for ``chr_pos`` inputs, and
lazy access to result plotting and liability-scale post-processing.

Design Notes
------------
- Only stable public objects should be re-exported here.
- Internal compute kernels remain under ``ldsc._kernel`` and are not part of
  the intended user-facing API.
- Plotting exports are resolved lazily so importing ``ldsc`` does not import
  Matplotlib until a figure-producing path is used.

Example
-------
>>> from ldsc import GlobalConfig, ReferencePanelBuilder
>>> GlobalConfig().snp_identifier
'chr_pos_allele_aware'
>>> isinstance(ReferencePanelBuilder(), ReferencePanelBuilder)
True
"""

from importlib import import_module

from .annotation_builder import (
    AnnotationBuilder,
    AnnotationBundle,
    run_bed_to_annot,
)
from .config import (
    AnnotationBuildConfig,
    ConfigMismatchError,
    GlobalConfig,
    GeneLDScoreIndexBuildConfig,
    LDScoreConfig,
    MungeConfig,
    ReferencePanelBuildConfig,
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
from .genome_build_inference import (
    ChrPosBuildInference,
    infer_chr_pos_build,
    resolve_genome_build,
    resolve_chr_pos_table,
)
from .hm3 import load_hm3_curated_map
from .ldscore_calculator import ChromLDScoreResult, LDScoreCalculator, LDScoreResult, run_ldscore
from .gene_ldscore_index import build_gene_ldscore_index, load_gene_ldscore_index
from .outputs import (
    H2DirectoryWriter,
    H2OutputConfig,
    LDScoreDirectoryWriter,
    LDScoreOutputConfig,
    PartitionedH2DirectoryWriter,
    PartitionedH2OutputConfig,
    QuantileH2DirectoryWriter,
    QuantileH2OutputConfig,
    RgDirectoryWriter,
    RgOutputConfig,
)
from .ref_panel_builder import ReferencePanelBuildResult, ReferencePanelBuilder, run_build_ref_panel
from ._kernel.ref_panel import ParquetR2RefPanel, PlinkRefPanel, RefPanel, RefPanelLoader
from .r2_query import R2Panel, query_r2, unbiased_r2_to_pearson_r

_LAZY_EXPORTS = {
    "MungeRunSummary": (".sumstats_munger", "MungeRunSummary"),
    "RawSumstatsInference": (".sumstats_munger", "RawSumstatsInference"),
    "SumstatsMunger": (".sumstats_munger", "SumstatsMunger"),
    "SumstatsTable": (".sumstats_munger", "SumstatsTable"),
    "infer_raw_sumstats": (".sumstats_munger", "infer_raw_sumstats"),
    "load_sumstats": (".sumstats_munger", "load_sumstats"),
    "load_ldscore_from_dir": (".regression_runner", "load_ldscore_from_dir"),
    "LegacyLDScoreConverter": (".legacy_ldscore_converter", "LegacyLDScoreConverter"),
    "LegacyLDScoreConversionResult": (".legacy_ldscore_converter", "LegacyLDScoreConversionResult"),
    "convert_ldsc2_ldscores": (".legacy_ldscore_converter", "convert_ldsc2_ldscores"),
    "RegressionDataset": (".regression_runner", "RegressionDataset"),
    "RegressionRunner": (".regression_runner", "RegressionRunner"),
    "RgResultFamily": (".regression_runner", "RgResultFamily"),
    "QuantileAssignment": (".quantile_h2", "QuantileAssignment"),
    "QuantileH2Result": (".quantile_h2", "QuantileH2Result"),
    "assign_legacy_quantiles": (".quantile_h2", "assign_legacy_quantiles"),
    "compute_quantile_h2": (".quantile_h2", "compute_quantile_h2"),
    "compute_standardized_coefficients": (".quantile_h2", "compute_standardized_coefficients"),
    "load_fitted_partitioned_model": (".quantile_h2", "load_fitted_partitioned_model"),
    "run_quantile_h2_from_args": (".quantile_h2", "run_quantile_h2_from_args"),
    "H2ScaleConversionArtifact": (".h2_scale", "H2ScaleConversionArtifact"),
    "convert_h2_scale": (".h2_scale", "convert_h2_scale"),
    "PlotArtifact": (".plotting", "PlotArtifact"),
    "plot_result": (".plotting", "plot_result"),
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
    "ParquetR2RefPanel",
    "R2Panel",
    "query_r2",
    "unbiased_r2_to_pearson_r",
    "PartitionedH2DirectoryWriter",
    "PartitionedH2OutputConfig",
    "QuantileAssignment",
    "QuantileH2DirectoryWriter",
    "QuantileH2OutputConfig",
    "QuantileH2Result",
    "PlinkRefPanel",
    "PlotArtifact",
    "ReferencePanelBuildConfig",
    "ReferencePanelBuildResult",
    "ReferencePanelBuilder",
    "RefPanel",
    "RefPanelConfig",
    "RefPanelLoader",
    "RegressionConfig",
    "RegressionDataset",
    "RegressionRunner",
    "RgDirectoryWriter",
    "RgOutputConfig",
    "RgResultFamily",
    "SumstatsMunger",
    "SumstatsTable",
    "load_sumstats",
    "plot_result",
    "convert_ldsc2_ldscores",
    "convert_h2_scale",
    "assign_legacy_quantiles",
    "compute_quantile_h2",
    "compute_standardized_coefficients",
    "infer_raw_sumstats",
    "run_build_ref_panel",
    "build_gene_ldscore_index",
    "run_bed_to_annot",
    "run_ldscore",
    "load_fitted_partitioned_model",
    "run_quantile_h2_from_args",
    "resolve_chr_pos_table",
    "reset_global_config",
    "resolve_genome_build",
    "set_global_config",
    "validate_config_compatibility",
]
