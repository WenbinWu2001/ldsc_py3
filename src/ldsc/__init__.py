"""ldsc package public surface.

Core functionality:
    Re-export the refactored LDSC workflow classes, configuration dataclasses,
    reference-panel abstractions, and convenience helpers from a single import
    location.

Overview
--------
This module is the top-level public API for the refactored package. Import from
``ldsc`` when you want the supported package surface rather than the internal
``ldsc._kernel`` implementation modules. The exports here mirror the main user
workflows: annotation building, LD-score calculation, summary-statistics
munging, regression, and output writing.

Design Notes
------------
- Only stable public objects should be re-exported here.
- Internal compute kernels remain under ``ldsc._kernel`` and are not part of
  the intended user-facing API.

Example
-------
>>> from ldsc import CommonConfig, LDScoreCalculator, RegressionRunner
>>> CommonConfig().snp_identifier
'chr_pos'
"""

from importlib import import_module

from .annotation_builder import (
    AnnotationBuilder,
    AnnotationBundle,
    AnnotationSourceSpec,
    gene_set_to_bed,
    make_annot_files,
    run_bed_to_annot,
)
from .config import (
    AnnotationBuildConfig,
    CommonConfig,
    LDScoreConfig,
    MungeConfig,
    ReferencePanelBuildConfig,
    RefPanelConfig,
    RegressionConfig,
)
from .ldscore_calculator import ChromLDScoreResult, LDScoreCalculator, LDScoreResult, run_ldscore
from .outputs import ArtifactConfig, ArtifactProducer, OutputManager, OutputSpec, PostProcessor, ResultFormatter, ResultWriter, RunSummary
from .ref_panel_builder import ReferencePanelBuildResult, ReferencePanelBuilder, run_build_ref_panel
from ._kernel.ref_panel import ParquetR2RefPanel, PlinkRefPanel, RefPanel, RefPanelLoader, RefPanelSpec

_LAZY_EXPORTS = {
    "MungeRunSummary": (".sumstats_munger", "MungeRunSummary"),
    "RawSumstatsSpec": (".sumstats_munger", "RawSumstatsSpec"),
    "SumstatsMunger": (".sumstats_munger", "SumstatsMunger"),
    "SumstatsTable": (".sumstats_munger", "SumstatsTable"),
    "load_sumstats": (".sumstats_munger", "load_sumstats"),
    "RegressionDataset": (".regression_runner", "RegressionDataset"),
    "RegressionRunner": (".regression_runner", "RegressionRunner"),
}


def __getattr__(name: str):
    """Lazily import public workflows with optional heavy dependencies."""
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
    "CommonConfig",
    "LDScoreCalculator",
    "LDScoreConfig",
    "LDScoreResult",
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
    "gene_set_to_bed",
    "load_sumstats",
    "make_annot_files",
    "run_build_ref_panel",
    "run_bed_to_annot",
    "run_ldscore",
]
