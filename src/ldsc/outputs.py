"""Flexible artifact-oriented output layer for workflow results."""

from __future__ import annotations

import gzip
import json
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass, field, is_dataclass
from os import PathLike
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .config import _normalize_optional_path, _normalize_required_path
from .path_resolution import ensure_output_directory, ensure_output_parent_directory


ArtifactLayout = str
CompressionMode = str


@dataclass(frozen=True)
class OutputSpec:
    """Configuration for artifact emission and on-disk layout."""
    out_prefix: str | PathLike[str]
    output_dir: str | PathLike[str] | None = None
    artifact_layout: ArtifactLayout = "flat"
    write_ldscore: bool = True
    write_w_ld: bool = True
    write_counts: bool = True
    write_annotation_manifest: bool = True
    write_per_chrom: bool = False
    aggregate_across_chromosomes: bool = True
    compression: CompressionMode = "gzip"
    overwrite: bool = False
    log_path: str | PathLike[str] | None = None
    write_summary_json: bool = True
    write_summary_tsv: bool = True
    write_run_metadata: bool = True
    enabled_artifacts: list[str] | None = None

    def __post_init__(self) -> None:
        """Normalize output paths and validate layout and compression choices."""
        if self.artifact_layout not in {"flat", "by_chrom", "run_dir"}:
            raise ValueError("artifact_layout must be 'flat', 'by_chrom', or 'run_dir'.")
        if self.compression not in {"gzip", "none"}:
            raise ValueError("compression must be 'gzip' or 'none'.")
        object.__setattr__(self, "out_prefix", _normalize_required_path(self.out_prefix))
        object.__setattr__(self, "output_dir", _normalize_optional_path(self.output_dir))
        object.__setattr__(self, "log_path", _normalize_optional_path(self.log_path))
        if not self.out_prefix:
            raise ValueError("out_prefix is required.")


@dataclass(frozen=True)
class ArtifactConfig:
    """Optional advanced settings for artifact-specific producers."""
    options: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class RunSummary:
    """Compact run metadata derived from a workflow result."""
    n_reference_snps: int
    n_regression_snps: int
    chromosomes_processed: list[str]
    count_artifacts_available: list[str]
    output_paths: dict[str, str] = field(default_factory=dict)
    config_snapshot: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class Artifact:
    """In-memory artifact payload prior to serialization."""
    name: str
    relative_path: str
    payload: Any
    format: str


class ArtifactProducer(ABC):
    """Extension interface for output producers."""
    name: str

    @abstractmethod
    def supports(self, result: Any) -> bool:
        """Return ``True`` if the producer can emit artifacts for ``result``."""
        raise NotImplementedError

    @abstractmethod
    def build(
        self,
        result: Any,
        run_summary: RunSummary,
        output_spec: OutputSpec,
        artifact_config: ArtifactConfig | None = None,
    ) -> list[Artifact]:
        """Build artifacts from ``result`` and ``run_summary``."""
        raise NotImplementedError


class ResultFormatter:
    """Derive clean summary tables and manifests from raw results."""

    def build_run_summary(
        self,
        result: Any,
        output_spec: OutputSpec,
        config_snapshot: dict[str, Any] | None = None,
    ) -> RunSummary:
        """Summarize row counts, chromosomes, and count-vector keys."""
        ref_meta = getattr(result, "reference_metadata", None)
        reg_meta = getattr(result, "regression_metadata", None)
        chromosomes = getattr(result, "chromosome_results", None) or []
        chromosome_names = [str(getattr(chrom_result, "chrom", "")) for chrom_result in chromosomes]
        count_keys = list(getattr(result, "snp_count_totals", {}).keys())
        return RunSummary(
            n_reference_snps=0 if ref_meta is None else len(ref_meta),
            n_regression_snps=0 if reg_meta is None else len(reg_meta),
            chromosomes_processed=[chrom for chrom in chromosome_names if chrom],
            count_artifacts_available=count_keys,
            config_snapshot=config_snapshot or {},
        )

    def build_annotation_manifest(self, result: Any) -> pd.DataFrame:
        """Build a baseline/query manifest from result column names."""
        baseline = list(getattr(result, "baseline_columns", []))
        query = list(getattr(result, "query_columns", []))
        rows = (
            [{"column": name, "group": "baseline"} for name in baseline]
            + [{"column": name, "group": "query"} for name in query]
        )
        return pd.DataFrame(rows, columns=["column", "group"])

    def build_summary_table(self, run_summary: RunSummary) -> pd.DataFrame:
        """Convert a :class:`RunSummary` into a one-row table."""
        row = {
            "n_reference_snps": run_summary.n_reference_snps,
            "n_regression_snps": run_summary.n_regression_snps,
            "chromosomes_processed": ",".join(run_summary.chromosomes_processed),
            "count_artifacts_available": ",".join(run_summary.count_artifacts_available),
        }
        return pd.DataFrame([row])


class ResultWriter:
    """Serialize artifact payloads to disk."""

    def write(self, artifact: Artifact, root: Path, overwrite: bool = False) -> str:
        """Write ``artifact`` below ``root`` and return the written path."""
        path = root / artifact.relative_path
        path = ensure_output_parent_directory(path, label="artifact output path")
        if path.exists() and not overwrite:
            raise FileExistsError(f"Refusing to overwrite existing artifact: {path}")

        if artifact.format == "dataframe":
            self._write_dataframe(artifact.payload, path)
        elif artifact.format == "json":
            self._write_json(artifact.payload, path)
        elif artifact.format == "text":
            path.write_text(str(artifact.payload), encoding="utf-8")
        else:
            raise ValueError(f"Unsupported artifact format: {artifact.format}")
        return str(path)

    def _write_dataframe(self, df: pd.DataFrame, path: Path) -> None:
        """Write a tabular artifact as TSV, optionally gzip-compressed."""
        if path.suffix == ".gz":
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                df.to_csv(handle, sep="\t", index=False, na_rep="NA")
            return
        df.to_csv(path, sep="\t", index=False, na_rep="NA")

    def _write_json(self, payload: Any, path: Path) -> None:
        """Write one JSON artifact, optionally gzip-compressed."""
        serializable = _to_serializable(payload)
        if path.suffix == ".gz":
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                json.dump(serializable, handle, indent=2, sort_keys=True)
            return
        path.write_text(json.dumps(serializable, indent=2, sort_keys=True), encoding="utf-8")


class LDScoreTableProducer(ArtifactProducer):
    """Producer for reference LD-score tables."""
    name = "ldscore"

    def supports(self, result: Any) -> bool:
        """Return ``True`` when aggregated LD-score tables are available."""
        return getattr(result, "reference_metadata", None) is not None and getattr(result, "ld_scores", None) is not None

    def build(self, result: Any, run_summary: RunSummary, output_spec: OutputSpec, artifact_config: ArtifactConfig | None = None) -> list[Artifact]:
        """Build LD-score artifacts, optionally filtered to regression SNP rows."""
        artifacts: list[Artifact] = []
        filter_to_regression_snps = _regression_output_enabled(result)
        if output_spec.aggregate_across_chromosomes:
            table = pd.concat([result.reference_metadata.reset_index(drop=True), result.ld_scores.reset_index(drop=True)], axis=1)
            table = _filter_ldscore_table(
                table,
                None if not filter_to_regression_snps else getattr(result, "regression_snps", None),
            )
            artifacts.append(Artifact(self.name, _artifact_filename(output_spec, ".l2.ldscore"), table, "dataframe"))
        if output_spec.write_per_chrom:
            for chrom_result in getattr(result, "chromosome_results", []):
                table = pd.concat(
                    [chrom_result.reference_metadata.reset_index(drop=True), chrom_result.ld_scores.reset_index(drop=True)],
                    axis=1,
                )
                table = _filter_ldscore_table(
                    table,
                    None if not filter_to_regression_snps else getattr(chrom_result, "regression_snps", None),
                    raise_on_empty=False,
                )
                if table.empty:
                    continue
                artifacts.append(
                    Artifact(
                        f"{self.name}.chrom_{chrom_result.chrom}",
                        _artifact_filename(output_spec, ".l2.ldscore", chrom=str(chrom_result.chrom)),
                        table,
                        "dataframe",
                    )
                )
        if filter_to_regression_snps and not artifacts:
            raise ValueError("After filtering to regression SNPs, no SNPs remain.")
        return artifacts


class WeightLDProducer(ArtifactProducer):
    """Producer for regression-weight LD-score tables."""
    name = "w_ld"

    def supports(self, result: Any) -> bool:
        """Return ``True`` when regression-weight tables are available."""
        return getattr(result, "regression_metadata", None) is not None and getattr(result, "w_ld", None) is not None

    def build(self, result: Any, run_summary: RunSummary, output_spec: OutputSpec, artifact_config: ArtifactConfig | None = None) -> list[Artifact]:
        """Build aggregate and optional per-chromosome weight-table artifacts."""
        artifacts: list[Artifact] = []
        if output_spec.aggregate_across_chromosomes:
            table = pd.concat([result.regression_metadata.reset_index(drop=True), result.w_ld.reset_index(drop=True)], axis=1)
            artifacts.append(Artifact(self.name, _artifact_filename(output_spec, ".w.l2.ldscore"), table, "dataframe"))
        if output_spec.write_per_chrom:
            for chrom_result in getattr(result, "chromosome_results", []):
                if getattr(chrom_result, "regression_metadata", None) is None or getattr(chrom_result, "w_ld", None) is None:
                    continue
                table = pd.concat(
                    [chrom_result.regression_metadata.reset_index(drop=True), chrom_result.w_ld.reset_index(drop=True)],
                    axis=1,
                )
                artifacts.append(
                    Artifact(
                        f"{self.name}.chrom_{chrom_result.chrom}",
                        _artifact_filename(output_spec, ".w.l2.ldscore", chrom=str(chrom_result.chrom)),
                        table,
                        "dataframe",
                    )
                )
        return artifacts


class CountProducer(ArtifactProducer):
    """Producer for count-vector artifacts such as ``.M`` and ``.M_5_50``."""
    name = "counts"

    def supports(self, result: Any) -> bool:
        """Return ``True`` when named SNP-count vectors are present."""
        return bool(getattr(result, "snp_count_totals", {}))

    def build(self, result: Any, run_summary: RunSummary, output_spec: OutputSpec, artifact_config: ArtifactConfig | None = None) -> list[Artifact]:
        """Build text artifacts for each named SNP-count vector."""
        artifacts: list[Artifact] = []
        count_map = getattr(result, "snp_count_totals", {})
        for key, values in count_map.items():
            suffix = _count_suffix(key)
            payload = "\t".join(str(value) for value in np.asarray(values).tolist())
            artifacts.append(Artifact(f"{self.name}.{key}", _artifact_filename(output_spec, suffix, compressed=False), payload, "text"))
        return artifacts


class AnnotationManifestProducer(ArtifactProducer):
    """Producer for the baseline/query annotation manifest."""
    name = "annotation_manifest"

    def __init__(self, formatter: ResultFormatter):
        """Initialize the producer with the shared result formatter."""
        self._formatter = formatter

    def supports(self, result: Any) -> bool:
        """Return ``True`` when baseline or query annotation columns exist."""
        return bool(getattr(result, "baseline_columns", None) or getattr(result, "query_columns", None))

    def build(self, result: Any, run_summary: RunSummary, output_spec: OutputSpec, artifact_config: ArtifactConfig | None = None) -> list[Artifact]:
        """Build the annotation-group manifest as one TSV artifact."""
        manifest = self._formatter.build_annotation_manifest(result)
        return [Artifact(self.name, _artifact_filename(output_spec, ".annotation_groups.tsv", compressed=False), manifest, "dataframe")]


class SummaryTSVProducer(ArtifactProducer):
    """Producer for a compact TSV run summary."""
    name = "summary_tsv"

    def __init__(self, formatter: ResultFormatter):
        """Initialize the producer with the shared summary formatter."""
        self._formatter = formatter

    def supports(self, result: Any) -> bool:
        """Always emit a TSV summary for supported run results."""
        return True

    def build(self, result: Any, run_summary: RunSummary, output_spec: OutputSpec, artifact_config: ArtifactConfig | None = None) -> list[Artifact]:
        """Build the one-row TSV summary artifact."""
        table = self._formatter.build_summary_table(run_summary)
        return [Artifact(self.name, _artifact_filename(output_spec, ".summary.tsv", compressed=False), table, "dataframe")]


class SummaryJSONProducer(ArtifactProducer):
    """Producer for a JSON-serialized run summary."""
    name = "summary_json"

    def supports(self, result: Any) -> bool:
        """Always emit a JSON summary for supported run results."""
        return True

    def build(self, result: Any, run_summary: RunSummary, output_spec: OutputSpec, artifact_config: ArtifactConfig | None = None) -> list[Artifact]:
        """Build the JSON summary artifact from ``run_summary``."""
        return [Artifact(self.name, _artifact_filename(output_spec, ".summary.json", compressed=False), run_summary, "json")]


class RunMetadataProducer(ArtifactProducer):
    """Producer for JSON run metadata."""
    name = "run_metadata"

    def supports(self, result: Any) -> bool:
        """Always emit run metadata for supported run results."""
        return True

    def build(self, result: Any, run_summary: RunSummary, output_spec: OutputSpec, artifact_config: ArtifactConfig | None = None) -> list[Artifact]:
        """Build the JSON metadata artifact describing the current run."""
        metadata = {
            "result_type": type(result).__name__,
            "config_snapshot": run_summary.config_snapshot,
            "chromosomes_processed": run_summary.chromosomes_processed,
        }
        return [Artifact(self.name, _artifact_filename(output_spec, ".run_metadata.json", compressed=False), metadata, "json")]


class OutputManager:
    """Coordinate output formatting, artifact production, and serialization."""

    def __init__(
        self,
        formatter: ResultFormatter | None = None,
        writer: ResultWriter | None = None,
        producers: Iterable[ArtifactProducer] | None = None,
    ) -> None:
        """Initialize the manager and register the configured producers."""
        self.formatter = formatter or ResultFormatter()
        self.writer = writer or ResultWriter()
        self._registry: dict[str, ArtifactProducer] = {}
        for producer in producers or _default_producers(self.formatter):
            self.register_producer(producer)

    def register_producer(self, producer: ArtifactProducer) -> None:
        """Register or replace an artifact producer by name."""
        self._registry[producer.name] = producer

    def available_artifacts(self) -> list[str]:
        """Return registered artifact names."""
        return sorted(self._registry.keys())

    def resolve_enabled_artifacts(self, output_spec: OutputSpec) -> list[str]:
        """Resolve enabled artifact names implied by ``output_spec``."""
        if output_spec.enabled_artifacts is not None:
            unknown = sorted(set(output_spec.enabled_artifacts) - set(self._registry.keys()))
            if unknown:
                raise ValueError(f"Unknown artifacts requested: {unknown}. Available: {self.available_artifacts()}")
            return list(output_spec.enabled_artifacts)

        enabled: list[str] = []
        if output_spec.write_ldscore:
            enabled.append("ldscore")
        if output_spec.write_w_ld:
            enabled.append("w_ld")
        if output_spec.write_counts:
            enabled.append("counts")
        if output_spec.write_annotation_manifest:
            enabled.append("annotation_manifest")
        if output_spec.write_summary_tsv:
            enabled.append("summary_tsv")
        if output_spec.write_summary_json:
            enabled.append("summary_json")
        if output_spec.write_run_metadata:
            enabled.append("run_metadata")
        return enabled

    def build_run_summary(
        self,
        result: Any,
        output_spec: OutputSpec,
        config_snapshot: dict[str, Any] | None = None,
    ) -> RunSummary:
        """Delegate run-summary construction to :class:`ResultFormatter`."""
        return self.formatter.build_run_summary(result, output_spec, config_snapshot=config_snapshot)

    def write_outputs(
        self,
        result: Any,
        output_spec: OutputSpec,
        artifact_config: ArtifactConfig | None = None,
        config_snapshot: dict[str, Any] | None = None,
    ) -> RunSummary:
        """Build and write all enabled artifacts for ``result``."""
        root = ensure_output_directory(_output_root(output_spec), label="output directory")
        run_summary = self.build_run_summary(result, output_spec, config_snapshot=config_snapshot)
        output_paths: dict[str, str] = {}

        for name in self.resolve_enabled_artifacts(output_spec):
            producer = self._registry[name]
            if not producer.supports(result):
                continue
            for artifact in producer.build(result, run_summary, output_spec, artifact_config=artifact_config):
                output_paths[artifact.name] = self.writer.write(artifact, root, overwrite=output_spec.overwrite)

        return RunSummary(
            n_reference_snps=run_summary.n_reference_snps,
            n_regression_snps=run_summary.n_regression_snps,
            chromosomes_processed=run_summary.chromosomes_processed,
            count_artifacts_available=run_summary.count_artifacts_available,
            output_paths=output_paths,
            config_snapshot=run_summary.config_snapshot,
        )


class PostProcessor:
    """Extension hook for future summary, plotting, and report producers."""

    def __init__(self, output_manager: OutputManager | None = None) -> None:
        """Attach a target output manager for future producer registrations."""
        self.output_manager = output_manager or OutputManager()

    def register_producer(self, producer: ArtifactProducer) -> None:
        """Register an additional producer on the attached output manager."""
        self.output_manager.register_producer(producer)


def _default_producers(formatter: ResultFormatter) -> list[ArtifactProducer]:
    """Return the built-in artifact producers used by the default output manager."""
    return [
        LDScoreTableProducer(),
        WeightLDProducer(),
        CountProducer(),
        AnnotationManifestProducer(formatter),
        SummaryTSVProducer(formatter),
        SummaryJSONProducer(),
        RunMetadataProducer(),
    ]


def _output_root(output_spec: OutputSpec) -> Path:
    """Resolve the filesystem root for the current artifact layout."""
    if output_spec.output_dir is None:
        return Path(".")
    if output_spec.artifact_layout == "run_dir":
        return Path(output_spec.output_dir) / Path(output_spec.out_prefix).name
    return Path(output_spec.output_dir)


def _artifact_filename(output_spec: OutputSpec, suffix: str, chrom: str | None = None, compressed: bool | None = None) -> str:
    """Build one artifact-relative filename from the configured layout rules."""
    compressed = output_spec.compression == "gzip" if compressed is None else compressed
    prefix = Path(output_spec.out_prefix).name
    if chrom is None or output_spec.artifact_layout == "flat":
        stem = f"{prefix}{suffix}"
    else:
        if output_spec.artifact_layout == "by_chrom":
            stem = f"chr{chrom}/{prefix}{suffix}"
        else:
            stem = f"{prefix}.chr{chrom}{suffix}"
    return stem + (".gz" if compressed else "")


def _count_suffix(key: str) -> str:
    """Map an internal count key to the legacy-compatible output suffix."""
    if key == "all_reference_snp_counts":
        return ".l2.M"
    if key == "common_reference_snp_counts_maf_gt_0_05":
        return ".l2.M_5_50"
    return f".{key}.tsv"


def _regression_output_enabled(result: Any) -> bool:
    """Return ``True`` when `regression_snps_path` should filter written LD-score rows."""
    config_snapshot = getattr(result, "config_snapshot", None)
    if config_snapshot is None:
        return False
    return getattr(config_snapshot, "regression_snps_path", None) is not None


def _filter_ldscore_table(
    table: pd.DataFrame,
    regression_snps: set[str] | None,
    *,
    raise_on_empty: bool = True,
) -> pd.DataFrame:
    """Filter one LD-score table to the authoritative regression SNP row set."""
    if regression_snps is None:
        return table
    filtered = table.loc[table["SNP"].astype(str).isin(regression_snps), :].reset_index(drop=True)
    if filtered.empty and raise_on_empty:
        raise ValueError("After filtering to regression SNPs, no SNPs remain.")
    return filtered


def _to_serializable(value: Any) -> Any:
    """Recursively convert workflow objects into JSON-serializable structures."""
    if is_dataclass(value):
        return _to_serializable(asdict(value))
    if isinstance(value, pd.DataFrame):
        return value.to_dict(orient="records")
    if isinstance(value, pd.Series):
        return value.to_list()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _to_serializable(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_to_serializable(item) for item in value]
    return value
