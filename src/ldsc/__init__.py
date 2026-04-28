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
calculation, summary-statistics munging, regression, and output writing.

Design Notes
------------
- Only stable public objects should be re-exported here.
- Internal compute kernels remain under ``ldsc._kernel`` and are not part of
  the intended user-facing API.

Example
-------
>>> from ldsc import GlobalConfig, ReferencePanelBuilder
>>> GlobalConfig().snp_identifier
'chr_pos'
>>> isinstance(ReferencePanelBuilder(), ReferencePanelBuilder)
True
"""

from importlib import import_module

from .annotation_builder import (
    AnnotationBuilder,
    AnnotationBundle,
    AnnotationSourceSpec,
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
from .ldscore_calculator import ChromLDScoreResult, LDScoreCalculator, LDScoreResult, run_ldscore
from .outputs import ArtifactConfig, ArtifactProducer, LDScoreDirectoryWriter, LDScoreOutputConfig, OutputManager, OutputSpec, PostProcessor, ResultFormatter, ResultWriter, RunSummary
from .ref_panel_builder import ReferencePanelBuildResult, ReferencePanelBuilder, run_build_ref_panel
from ._kernel.ref_panel import ParquetR2RefPanel, PlinkRefPanel, RefPanel, RefPanelLoader, RefPanelSpec

_LAZY_EXPORTS = {
    "MungeRunSummary": (".sumstats_munger", "MungeRunSummary"),
    "RawSumstatsSpec": (".sumstats_munger", "RawSumstatsSpec"),
    "SumstatsMunger": (".sumstats_munger", "SumstatsMunger"),
    "SumstatsTable": (".sumstats_munger", "SumstatsTable"),
    "load_sumstats": (".sumstats_munger", "load_sumstats"),
    "load_ldscore_from_dir": (".regression_runner", "load_ldscore_from_dir"),
    "RegressionDataset": (".regression_runner", "RegressionDataset"),
    "RegressionRunner": (".regression_runner", "RegressionRunner"),
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
    "AnnotationSourceSpec",
    "ArtifactConfig",
    "ArtifactProducer",
    "ChromLDScoreResult",
    "ConfigMismatchError",
    "GlobalConfig",
    "get_global_config",
    "LDScoreCalculator",
    "LDScoreConfig",
    "LDScoreDirectoryWriter",
    "LDScoreOutputConfig",
    "LDScoreResult",
    "load_ldscore_from_dir",
    "MungeConfig",
    "MungeRunSummary",
    "OutputManager",
    "OutputSpec",
    "ParquetR2RefPanel",
    "PlinkRefPanel",
    "PostProcessor",
    "RawSumstatsSpec",
    "ReferencePanelBuildConfig",
    "ReferencePanelBuildResult",
    "ReferencePanelBuilder",
    "RefPanel",
    "RefPanelConfig",
    "RefPanelLoader",
    "RefPanelSpec",
    "RegressionConfig",
    "RegressionDataset",
    "RegressionRunner",
    "ResultFormatter",
    "ResultWriter",
    "RunSummary",
    "SumstatsMunger",
    "SumstatsTable",
    "load_sumstats",
    "run_build_ref_panel",
    "run_bed_to_annot",
    "run_ldscore",
    "reset_global_config",
    "set_global_config",
    "validate_config_compatibility",
]
