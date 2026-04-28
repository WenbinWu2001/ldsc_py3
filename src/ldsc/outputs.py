"""outputs.py

Core functionality:
    Build and serialize workflow artifacts such as canonical LD-score result
    directories, auxiliary tables, count vectors, and run summaries.

Overview
--------
This module owns output writing for the refactored package. The current public
LD-score workflow uses ``LDScoreDirectoryWriter`` to write a canonical result
directory containing ``manifest.json``, ``baseline.parquet``, and optional
``query.parquet``. ``OutputManager`` remains a generic artifact writer for older
or lower-level result shapes that still expose merged ``ldscore_table`` and
``snp_count_totals`` attributes.

Design Notes
------------
- Canonical LD-score writes use fixed filenames inside a directory and do not
  emit public `.M`, `.M_5_50`, `.l2.ldscore.gz`, or `.w.l2.ldscore.gz` files.
- ``LDScoreTableProducer`` mirrors older merged-table result shapes and is not
  used by the public ``ldsc ldscore`` command.
- Additional artifacts can be registered through ``ArtifactProducer``.
"""

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
LDSCORE_RESULT_FORMAT = "ldsc.ldscore_result.v1"


@dataclass(frozen=True)
class OutputSpec:
    """Configuration for artifact emission and on-disk layout.

    This generic artifact spec supports older merged-table result shapes that
    write one ``.l2.ldscore.gz`` file per chromosome, aggregate count vectors,
    and summary metadata. The public LD-score CLI now uses
    ``LDScoreOutputConfig`` and ``LDScoreDirectoryWriter`` instead. When
    ``overwrite`` is ``False``, batched writes fail before the first file is
    emitted if any destination path already exists.
    """
    out_prefix: str | PathLike[str]
    output_dir: str | PathLike[str] | None = None
    artifact_layout: ArtifactLayout = "flat"
    write_ldscore: bool = True
    write_w_ld: bool = False
    write_counts: bool = True
    write_annotation_manifest: bool = True
    write_per_chrom: bool = True
    aggregate_across_chromosomes: bool = False
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
class LDScoreOutputConfig:
    """Directory-oriented output config for canonical LD-score results.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives ``manifest.json``, ``baseline.parquet``, and
        optional ``query.parquet``.
    overwrite : bool, optional
        If ``True``, replace existing canonical LD-score files in
        ``output_dir``. Default is ``False``.
    parquet_compression : {"snappy", "gzip", "brotli", "zstd", "none", None}, optional
        Compression codec used for parquet tables. Default is ``"snappy"``.
    """
    output_dir: str | PathLike[str]
    overwrite: bool = False
    parquet_compression: str = "snappy"

    def __post_init__(self) -> None:
        """Normalize the output directory and validate parquet compression."""
        object.__setattr__(self, "output_dir", _normalize_required_path(self.output_dir))
        allowed = {"snappy", "gzip", "brotli", "zstd", "none", None}
        if self.parquet_compression not in allowed:
            raise ValueError("parquet_compression must be one of 'snappy', 'gzip', 'brotli', 'zstd', 'none', or None.")


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


class LDScoreDirectoryWriter:
    """Write the canonical LD-score result directory."""

    def write(self, result: Any, output_config: LDScoreOutputConfig) -> dict[str, str]:
        """Write ``manifest.json``, ``baseline.parquet``, and optional ``query.parquet``."""
        output_dir = ensure_output_directory(output_config.output_dir, label="LD-score output directory")
        baseline_table = getattr(result, "baseline_table", None)
        query_table = getattr(result, "query_table", None)
        if baseline_table is None:
            raise ValueError("LDScoreResult is missing baseline_table.")
        self._validate_tables(result)

        paths = {
            "manifest": output_dir / "manifest.json",
            "baseline": output_dir / "baseline.parquet",
        }
        if query_table is not None:
            paths["query"] = output_dir / "query.parquet"
        existing = [path for path in paths.values() if path.exists()]
        if existing and not output_config.overwrite:
            raise FileExistsError(_format_existing_artifact_message(existing))

        compression = None if output_config.parquet_compression in {None, "none"} else output_config.parquet_compression
        baseline_table.to_parquet(paths["baseline"], index=False, compression=compression)
        if query_table is not None:
            query_table.to_parquet(paths["query"], index=False, compression=compression)
        manifest = self.build_manifest(result, files={name: path.name for name, path in paths.items() if name != "manifest"})
        paths["manifest"].write_text(json.dumps(_to_serializable(manifest), indent=2, sort_keys=True), encoding="utf-8")
        return {name: str(path) for name, path in paths.items()}

    def build_manifest(self, result: Any, files: dict[str, str]) -> dict[str, Any]:
        """Build the JSON manifest payload for one LD-score result."""
        baseline_table = getattr(result, "baseline_table")
        query_table = getattr(result, "query_table", None)
        config_snapshot = getattr(result, "config_snapshot", None)
        chromosomes = baseline_table["CHR"].astype(str).drop_duplicates().tolist()
        return {
            "format": LDSCORE_RESULT_FORMAT,
            "files": dict(files),
            "snp_identifier": getattr(config_snapshot, "snp_identifier", None),
            "genome_build": getattr(config_snapshot, "genome_build", None),
            "chromosomes": chromosomes,
            "baseline_columns": list(getattr(result, "baseline_columns", [])),
            "query_columns": list(getattr(result, "query_columns", [])),
            "counts": list(getattr(result, "count_records", [])),
            "config_snapshot": config_snapshot,
            "n_baseline_rows": int(len(baseline_table)),
            "n_query_rows": 0 if query_table is None else int(len(query_table)),
        }

    def _validate_tables(self, result: Any) -> None:
        """Validate baseline/query table shape before any files are written."""
        baseline_table = getattr(result, "baseline_table")
        query_table = getattr(result, "query_table", None)
        baseline_columns = list(getattr(result, "baseline_columns", []))
        query_columns = list(getattr(result, "query_columns", []))
        required_baseline = ["CHR", "SNP", "BP", "regr_weight", *baseline_columns]
        missing = [column for column in required_baseline if column not in baseline_table.columns]
        if missing:
            raise ValueError(f"baseline_table is missing required columns: {missing}")
        if query_columns and query_table is None:
            raise ValueError("query_table is required when query_columns are present.")
        if not query_columns and query_table is not None:
            raise ValueError("query_table was provided but query_columns is empty.")
        if query_table is None:
            return
        required_query = ["CHR", "SNP", "BP", *query_columns]
        missing = [column for column in required_query if column not in query_table.columns]
        if missing:
            raise ValueError(f"query_table is missing required columns: {missing}")
        baseline_keys = baseline_table.loc[:, ["CHR", "SNP", "BP"]].reset_index(drop=True)
        query_keys = query_table.loc[:, ["CHR", "SNP", "BP"]].reset_index(drop=True)
        if not baseline_keys.equals(query_keys):
            raise ValueError("query rows must match baseline rows on CHR/SNP/BP.")


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
        ldscore_table = getattr(result, "ldscore_table", None)
        ld_reference_snps = getattr(result, "ld_reference_snps", frozenset())
        ld_regression_snps = getattr(result, "ld_regression_snps", frozenset())
        chromosomes = getattr(result, "chromosome_results", None) or []
        chromosome_names = [str(getattr(chrom_result, "chrom", "")) for chrom_result in chromosomes]
        count_keys = list(getattr(result, "snp_count_totals", {}).keys())
        n_rows = 0 if ldscore_table is None else len(ldscore_table)
        return RunSummary(
            n_reference_snps=len(ld_reference_snps) if ld_reference_snps else n_rows,
            n_regression_snps=len(ld_regression_snps) if ld_regression_snps else n_rows,
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

    def resolve_path(self, artifact: Artifact, root: Path) -> Path:
        """Resolve one artifact path and ensure its parent directory exists."""
        path = root / artifact.relative_path
        return ensure_output_parent_directory(path, label="artifact output path")

    def prepare_writes(self, artifacts: Iterable[Artifact], root: Path, overwrite: bool = False) -> dict[str, Path]:
        """Resolve and validate a batch of artifact destinations.

        ``artifacts`` is typically the full output set for one run. This
        helper ensures every destination path is unique, creates any missing
        parent directories, and raises before the first payload is written when
        existing files would be reused while ``overwrite`` is ``False``.
        """
        prepared: dict[str, Path] = {}
        artifact_names_by_path: dict[Path, str] = {}
        duplicate_targets: list[str] = []
        existing_paths: list[Path] = []

        for artifact in artifacts:
            path = self.resolve_path(artifact, root)
            other_name = artifact_names_by_path.get(path)
            if other_name is not None:
                duplicate_targets.append(f"{path} ({other_name}, {artifact.name})")
            artifact_names_by_path[path] = artifact.name
            prepared[artifact.name] = path
            if path.exists():
                existing_paths.append(path)

        if duplicate_targets:
            joined = ", ".join(duplicate_targets)
            raise ValueError(f"Multiple artifacts resolve to the same output path: {joined}")
        if existing_paths and not overwrite:
            raise FileExistsError(_format_existing_artifact_message(existing_paths))
        return prepared

    def write(self, artifact: Artifact, root: Path, overwrite: bool = False) -> str:
        """Write one artifact below ``root`` and return the written path.

        This method validates only the destination for ``artifact``. Callers
        writing multiple artifacts should preflight the full batch with
        :meth:`prepare_writes` to avoid partial output directories.
        """
        path = self.resolve_path(artifact, root)
        if path.exists() and not overwrite:
            raise FileExistsError(_format_existing_artifact_message([path]))

        return self.write_prepared(artifact, path)

    def write_prepared(self, artifact: Artifact, path: Path) -> str:
        """Write ``artifact`` to a destination already checked for conflicts."""
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
    """Producer for merged LD-score tables with embedded regression weights."""
    name = "ldscore"

    def supports(self, result: Any) -> bool:
        """Return ``True`` when normalized LD-score tables are available."""
        return getattr(result, "ldscore_table", None) is not None

    def build(self, result: Any, run_summary: RunSummary, output_spec: OutputSpec, artifact_config: ArtifactConfig | None = None) -> list[Artifact]:
        """Build per-chromosome `.l2.ldscore.gz` artifacts."""
        artifacts: list[Artifact] = []
        chromosome_results = getattr(result, "chromosome_results", []) or []
        if chromosome_results:
            for chrom_result in chromosome_results:
                if chrom_result.ldscore_table.empty:
                    continue
                artifacts.append(
                    Artifact(
                        f"{self.name}.chrom_{chrom_result.chrom}",
                        _artifact_filename(output_spec, ".l2.ldscore", chrom=str(chrom_result.chrom)),
                        chrom_result.ldscore_table.reset_index(drop=True),
                        "dataframe",
                    )
                )
        elif getattr(result, "ldscore_table", None) is not None:
            artifacts.append(
                Artifact(
                    self.name,
                    _artifact_filename(output_spec, ".l2.ldscore"),
                    result.ldscore_table.reset_index(drop=True),
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
        """Build and write all enabled artifacts for ``result``.

        The manager first materializes the enabled artifact set, then asks the
        writer to validate every destination path. This makes overwrite
        failures fail-fast at the batch level instead of after a partial write.
        """
        root = ensure_output_directory(_output_root(output_spec), label="output directory")
        run_summary = self.build_run_summary(result, output_spec, config_snapshot=config_snapshot)
        output_paths: dict[str, str] = {}
        artifacts = self._build_artifacts(result, run_summary, output_spec, artifact_config=artifact_config)
        prepared_paths = self.writer.prepare_writes(artifacts, root, overwrite=output_spec.overwrite)

        for artifact in artifacts:
            output_paths[artifact.name] = self.writer.write_prepared(artifact, prepared_paths[artifact.name])

        return RunSummary(
            n_reference_snps=run_summary.n_reference_snps,
            n_regression_snps=run_summary.n_regression_snps,
            chromosomes_processed=run_summary.chromosomes_processed,
            count_artifacts_available=run_summary.count_artifacts_available,
            output_paths=output_paths,
            config_snapshot=run_summary.config_snapshot,
        )

    def _build_artifacts(
        self,
        result: Any,
        run_summary: RunSummary,
        output_spec: OutputSpec,
        artifact_config: ArtifactConfig | None = None,
    ) -> list[Artifact]:
        """Build artifacts from all enabled producers without writing them yet."""
        artifacts: list[Artifact] = []

        for name in self.resolve_enabled_artifacts(output_spec):
            producer = self._registry[name]
            if not producer.supports(result):
                continue
            artifacts.extend(producer.build(result, run_summary, output_spec, artifact_config=artifact_config))
        return artifacts


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
    if chrom is None:
        stem = f"{prefix}{suffix}"
    elif output_spec.artifact_layout == "by_chrom":
        stem = f"chr{chrom}/{prefix}.{chrom}{suffix}"
    else:
        stem = f"{prefix}.{chrom}{suffix}"
    return stem + (".gz" if compressed else "")


def _count_suffix(key: str) -> str:
    """Map an internal count key to the legacy-compatible output suffix."""
    if key == "all_reference_snp_counts":
        return ".l2.M"
    if key == "common_reference_snp_counts_maf_gt_0_05":
        return ".l2.M_5_50"
    return f".{key}.tsv"


def _format_existing_artifact_message(paths: list[Path]) -> str:
    """Build the shared overwrite-guard message for one or more destinations."""
    preview = ", ".join(str(path) for path in paths[:3])
    if len(paths) > 3:
        preview = f"{preview}, ... ({len(paths) - 3} more)"
    label = "artifact" if len(paths) == 1 else "artifacts"
    return f"Refusing to overwrite existing {label}: {preview}. Pass overwrite=True to replace them."


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
