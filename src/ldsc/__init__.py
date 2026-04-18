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
    RefPanelConfig,
    RegressionConfig,
)
from .ldscore_calculator import ChromLDScoreResult, LDScoreCalculator, LDScoreResult, run_ldscore
from .outputs import ArtifactConfig, ArtifactProducer, OutputManager, OutputSpec, PostProcessor, ResultFormatter, ResultWriter, RunSummary
from .regression_runner import RegressionDataset, RegressionRunner
from .sumstats_munger import MungeRunSummary, RawSumstatsSpec, SumstatsMunger, SumstatsTable
from ._kernel.ref_panel import ParquetR2RefPanel, PlinkRefPanel, RefPanel, RefPanelLoader, RefPanelSpec

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
    "make_annot_files",
    "run_bed_to_annot",
    "run_ldscore",
]
