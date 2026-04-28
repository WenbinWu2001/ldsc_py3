# Config Design Implementation Status

This note records the implementation outcome for the config-provenance design
described in [config-design.md](config-design.md).

## Implemented Behavior

- `GlobalConfig` remains immutable and is still the package-wide defaulting
  mechanism for notebook-style workflows.
- `AnnotationBundle`, `ChromLDScoreResult`, `LDScoreResult`, `SumstatsTable`,
  and `RegressionDataset` now carry `config_snapshot`.
- `RegressionRunner.build_dataset()` validates `SumstatsTable` and
  `LDScoreResult` snapshots with `validate_config_compatibility()` when both
  snapshots are known.
- `LDScoreCalculator.run()` validates `AnnotationBundle.config_snapshot`
  against the active runtime `GlobalConfig`.
- `LDScoreCalculator.run()` validates `RefPanelConfig.genome_build` against the
  active runtime `GlobalConfig.genome_build`.
- LD-score aggregation validates that per-chromosome
  `ChromLDScoreResult.config_snapshot` values agree before constructing the
  aggregate `LDScoreResult`.
- `ConfigMismatchError` and `validate_config_compatibility()` are re-exported
  from `ldsc`.

## Deliberate Compatibility Decisions

- Legacy objects with `config_snapshot=None` are still accepted at merge
  points. Strict compatibility checks run only when both sides carry snapshots.
- `RegressionDataset.config_snapshot` is propagated from the input
  `LDScoreResult` when present. The implementation does not synthesize a new
  snapshot from the current runner config when upstream provenance is missing.
- `load_sumstats()` cannot recover original munge-time config from disk, so it
  warns and leaves `config_snapshot=None`.
- `load_ldscore_from_dir()` warns and leaves `config_snapshot=None` when a
  legacy or malformed manifest lacks usable config provenance.

## Tests Covering This Design

- `tests/test_global_config_registry.py`
  Covers `ConfigMismatchError`, `validate_config_compatibility()`, and warning
  behavior for advisory fields.
- `tests/test_regression_workflow.py`
  Covers regression-time mismatch detection and the legacy `None` snapshot
  compatibility path.
- `tests/test_ldscore_workflow.py`
  Covers LD-score runtime mismatch detection for `AnnotationBundle` vs active
  `GlobalConfig`.
- `tests/test_sumstats_munger.py`
  Covers `config_snapshot` propagation during munging and the
  `load_sumstats()` provenance warning.
