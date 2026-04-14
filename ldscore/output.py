"""Flexible artifact-oriented output layer for workflow results."""

from __future__ import annotations

import gzip
import json
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass, field, is_dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


ArtifactLayout = str
CompressionMode = str


@dataclass(frozen=True)
class OutputSpec:
    out_prefix: str
    output_dir: str | None = None
    artifact_layout: ArtifactLayout = "flat"
    write_ldscore: bool = True
    write_w_ld: bool = True
    write_counts: bool = True
    write_annotation_manifest: bool = True
    write_per_chrom: bool = False
    aggregate_across_chromosomes: bool = True
    compression: CompressionMode = "gzip"
    overwrite: bool = False
    log_path: str | None = None
    write_summary_json: bool = True
    write_summary_tsv: bool = True
    write_run_metadata: bool = True
    enabled_artifacts: list[str] | None = None

    def __post_init__(self) -> None:
        if self.artifact_layout not in {"flat", "by_chrom", "run_dir"}:
            raise ValueError("artifact_layout must be 'flat', 'by_chrom', or 'run_dir'.")
        if self.compression not in {"gzip", "none"}:
            raise ValueError("compression must be 'gzip' or 'none'.")
        if not self.out_prefix:
            raise ValueError("out_prefix is required.")


@dataclass(frozen=True)
class ArtifactConfig:
    options: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class RunSummary:
    n_reference_snps: int
    n_regression_snps: int
    chromosomes_processed: list[str]
    count_artifacts_available: list[str]
    output_paths: dict[str, str] = field(default_factory=dict)
    config_snapshot: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class Artifact:
    name: str
    relative_path: str
    payload: Any
    format: str


class ArtifactProducer(ABC):
    name: str

    @abstractmethod
    def supports(self, result: Any) -> bool:
        raise NotImplementedError

    @abstractmethod
    def build(
        self,
        result: Any,
        run_summary: RunSummary,
        output_spec: OutputSpec,
        artifact_config: ArtifactConfig | None = None,
    ) -> list[Artifact]:
        raise NotImplementedError


class ResultFormatter:
    """Derive clean summary tables and manifests from raw results."""

    def build_run_summary(
        self,
        result: Any,
        output_spec: OutputSpec,
        config_snapshot: dict[str, Any] | None = None,
    ) -> RunSummary:
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
        baseline = list(getattr(result, "baseline_columns", []))
        query = list(getattr(result, "query_columns", []))
        rows = (
            [{"column": name, "group": "baseline"} for name in baseline]
            + [{"column": name, "group": "query"} for name in query]
        )
        return pd.DataFrame(rows, columns=["column", "group"])

    def build_summary_table(self, run_summary: RunSummary) -> pd.DataFrame:
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
        path = root / artifact.relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
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
        if path.suffix == ".gz":
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                df.to_csv(handle, sep="\t", index=False, na_rep="NA")
            return
        df.to_csv(path, sep="\t", index=False, na_rep="NA")

    def _write_json(self, payload: Any, path: Path) -> None:
        serializable = _to_serializable(payload)
        if path.suffix == ".gz":
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                json.dump(serializable, handle, indent=2, sort_keys=True)
            return
        path.write_text(json.dumps(serializable, indent=2, sort_keys=True), encoding="utf-8")


class LDScoreTableProducer(ArtifactProducer):
    name = "ldscore"

    def supports(self, result: Any) -> bool:
        return getattr(result, "reference_metadata", None) is not None and getattr(result, "ld_scores", None) is not None

    def build(self, result: Any, run_summary: RunSummary, output_spec: OutputSpec, artifact_config: ArtifactConfig | None = None) -> list[Artifact]:
        artifacts: list[Artifact] = []
        if output_spec.aggregate_across_chromosomes:
            table = pd.concat([result.reference_metadata.reset_index(drop=True), result.ld_scores.reset_index(drop=True)], axis=1)
            artifacts.append(Artifact(self.name, _artifact_filename(output_spec, ".l2.ldscore"), table, "dataframe"))
        if output_spec.write_per_chrom:
            for chrom_result in getattr(result, "chromosome_results", []):
                table = pd.concat(
                    [chrom_result.reference_metadata.reset_index(drop=True), chrom_result.ld_scores.reset_index(drop=True)],
                    axis=1,
                )
                artifacts.append(
                    Artifact(
                        f"{self.name}.chrom_{chrom_result.chrom}",
                        _artifact_filename(output_spec, ".l2.ldscore", chrom=str(chrom_result.chrom)),
                        table,
                        "dataframe",
                    )
                )
        return artifacts


class WeightLDProducer(ArtifactProducer):
    name = "w_ld"

    def supports(self, result: Any) -> bool:
        return getattr(result, "regression_metadata", None) is not None and getattr(result, "w_ld", None) is not None

    def build(self, result: Any, run_summary: RunSummary, output_spec: OutputSpec, artifact_config: ArtifactConfig | None = None) -> list[Artifact]:
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
    name = "counts"

    def supports(self, result: Any) -> bool:
        return bool(getattr(result, "snp_count_totals", {}))

    def build(self, result: Any, run_summary: RunSummary, output_spec: OutputSpec, artifact_config: ArtifactConfig | None = None) -> list[Artifact]:
        artifacts: list[Artifact] = []
        count_map = getattr(result, "snp_count_totals", {})
        for key, values in count_map.items():
            suffix = _count_suffix(key)
            frame = pd.DataFrame(np.atleast_2d(values))
            artifacts.append(Artifact(f"{self.name}.{key}", _artifact_filename(output_spec, suffix), frame, "dataframe"))
        return artifacts


class AnnotationManifestProducer(ArtifactProducer):
    name = "annotation_manifest"

    def __init__(self, formatter: ResultFormatter):
        self._formatter = formatter

    def supports(self, result: Any) -> bool:
        return bool(getattr(result, "baseline_columns", None) or getattr(result, "query_columns", None))

    def build(self, result: Any, run_summary: RunSummary, output_spec: OutputSpec, artifact_config: ArtifactConfig | None = None) -> list[Artifact]:
        manifest = self._formatter.build_annotation_manifest(result)
        return [Artifact(self.name, _artifact_filename(output_spec, ".annotation_manifest.tsv", compressed=False), manifest, "dataframe")]


class SummaryTSVProducer(ArtifactProducer):
    name = "summary_tsv"

    def __init__(self, formatter: ResultFormatter):
        self._formatter = formatter

    def supports(self, result: Any) -> bool:
        return True

    def build(self, result: Any, run_summary: RunSummary, output_spec: OutputSpec, artifact_config: ArtifactConfig | None = None) -> list[Artifact]:
        table = self._formatter.build_summary_table(run_summary)
        return [Artifact(self.name, _artifact_filename(output_spec, ".summary.tsv", compressed=False), table, "dataframe")]


class SummaryJSONProducer(ArtifactProducer):
    name = "summary_json"

    def supports(self, result: Any) -> bool:
        return True

    def build(self, result: Any, run_summary: RunSummary, output_spec: OutputSpec, artifact_config: ArtifactConfig | None = None) -> list[Artifact]:
        return [Artifact(self.name, _artifact_filename(output_spec, ".summary.json", compressed=False), run_summary, "json")]


class RunMetadataProducer(ArtifactProducer):
    name = "run_metadata"

    def supports(self, result: Any) -> bool:
        return True

    def build(self, result: Any, run_summary: RunSummary, output_spec: OutputSpec, artifact_config: ArtifactConfig | None = None) -> list[Artifact]:
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
        self.formatter = formatter or ResultFormatter()
        self.writer = writer or ResultWriter()
        self._registry: dict[str, ArtifactProducer] = {}
        for producer in producers or _default_producers(self.formatter):
            self.register_producer(producer)

    def register_producer(self, producer: ArtifactProducer) -> None:
        self._registry[producer.name] = producer

    def available_artifacts(self) -> list[str]:
        return sorted(self._registry.keys())

    def resolve_enabled_artifacts(self, output_spec: OutputSpec) -> list[str]:
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
        return self.formatter.build_run_summary(result, output_spec, config_snapshot=config_snapshot)

    def write_outputs(
        self,
        result: Any,
        output_spec: OutputSpec,
        artifact_config: ArtifactConfig | None = None,
        config_snapshot: dict[str, Any] | None = None,
    ) -> RunSummary:
        root = _output_root(output_spec)
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
        self.output_manager = output_manager or OutputManager()

    def register_producer(self, producer: ArtifactProducer) -> None:
        self.output_manager.register_producer(producer)


def _default_producers(formatter: ResultFormatter) -> list[ArtifactProducer]:
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
    if output_spec.output_dir is None:
        return Path(".")
    if output_spec.artifact_layout == "run_dir":
        return Path(output_spec.output_dir) / Path(output_spec.out_prefix).name
    return Path(output_spec.output_dir)


def _artifact_filename(output_spec: OutputSpec, suffix: str, chrom: str | None = None, compressed: bool | None = None) -> str:
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
    if key == "all_reference_snp_counts":
        return ".l2.M"
    if key == "common_reference_snp_counts_maf_gt_0_05":
        return ".l2.M_5_50"
    return f".{key}.tsv"


def _to_serializable(value: Any) -> Any:
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
