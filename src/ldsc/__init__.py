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
workflows: annotation building, parquet reference-panel building, LD-score
calculation, summary-statistics munging, regression, output writing, and
genome-build inference for ``chr_pos`` inputs.

Design Notes
------------
- Only stable public objects should be re-exported here.
- Internal compute kernels remain under ``ldsc._kernel`` and are not part of
  the intended user-facing API.

Example
-------
>>> from ldsc import GlobalConfig, ReferencePanelBuilder
>>> GlobalConfig(snp_identifier="rsid").snp_identifier
'rsid'
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
from .ldscore_calculator import ChromLDScoreResult, LDScoreCalculator, LDScoreResult, run_ldscore
from .outputs import (
    LDScoreDirectoryWriter,
    LDScoreOutputConfig,
    PartitionedH2DirectoryWriter,
    PartitionedH2OutputConfig,
    RgDirectoryWriter,
    RgOutputConfig,
)
from .ref_panel_builder import ReferencePanelBuildResult, ReferencePanelBuilder, run_build_ref_panel
from ._kernel.ref_panel import ParquetR2RefPanel, PlinkRefPanel, RefPanel, RefPanelLoader

_LAZY_EXPORTS = {
    "MungeRunSummary": (".sumstats_munger", "MungeRunSummary"),
    "SumstatsMunger": (".sumstats_munger", "SumstatsMunger"),
    "SumstatsTable": (".sumstats_munger", "SumstatsTable"),
    "load_sumstats": (".sumstats_munger", "load_sumstats"),
    "load_ldscore_from_dir": (".regression_runner", "load_ldscore_from_dir"),
    "RegressionDataset": (".regression_runner", "RegressionDataset"),
    "RegressionRunner": (".regression_runner", "RegressionRunner"),
    "RgResultFamily": (".regression_runner", "RgResultFamily"),
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
    "get_global_config",
    "infer_chr_pos_build",
    "LDScoreCalculator",
    "LDScoreConfig",
    "LDScoreDirectoryWriter",
    "LDScoreOutputConfig",
    "LDScoreResult",
    "LDSCConfigError",
    "LDSCDependencyError",
    "LDSCError",
    "LDSCInputError",
    "LDSCInternalError",
    "LDSCUsageError",
    "LDSCUserError",
    "load_ldscore_from_dir",
    "MungeConfig",
    "MungeRunSummary",
    "ParquetR2RefPanel",
    "PartitionedH2DirectoryWriter",
    "PartitionedH2OutputConfig",
    "PlinkRefPanel",
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
    "run_build_ref_panel",
    "run_bed_to_annot",
    "run_ldscore",
    "resolve_chr_pos_table",
    "reset_global_config",
    "resolve_genome_build",
    "set_global_config",
    "validate_config_compatibility",
]
