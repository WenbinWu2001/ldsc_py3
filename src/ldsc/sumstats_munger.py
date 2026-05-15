"""Workflow wrapper for LDSC summary-statistics munging.

Core functionality:
    Expose a typed, package-level interface to the legacy-compatible munging
    implementation and provide a public loader for curated Parquet or legacy
    ``.sumstats.gz`` artifacts.

Overview
--------
This module converts the historical munging script behavior into explicit
Python dataclasses and a small service object. The public workflow boundary
accepts path-like inputs, normalizes them once, and then passes primitive
values into the internal kernel so numerical behavior stays aligned with
established LDSC outputs. The workflow layer owns CLI orchestration, output
preflight, the fixed ``sumstats.log`` file, metadata sidecars, and result
objects; the kernel keeps the legacy-compatible parsing and filtering
primitives. Run summaries and metadata expose curated data artifacts only;
the log is an audit file and is not included in ``output_paths``.
Summary-statistics metadata sidecars stay thin for downstream compatibility;
coordinate and liftover provenance is written as readable workflow-log text.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field, replace
import json
import logging
from os import PathLike
from pathlib import Path
from typing import Any

import pandas as pd

from .chromosome_inference import chrom_sort_key, normalize_chromosome_series
from ._coordinates import (
    coordinate_missing_mask,
    positive_int_position_series,
)
from .column_inference import (
    INTERNAL_SUMSTATS_ARTIFACT_SPEC_MAP,
    normalize_genome_build,
    normalize_snp_identifier_mode,
    resolve_optional_column,
    resolve_required_column,
)
from .config import GlobalConfig, MungeConfig, _normalize_trait_name, get_global_config
from .errors import LDSCDependencyError
from .hm3 import packaged_hm3_curated_map_path
from .path_resolution import (
    ensure_output_directory,
    preflight_output_artifact_family,
    remove_output_artifacts,
    resolve_scalar_path,
)
from ._logging import log_inputs, log_outputs, workflow_logging
from ._kernel.snp_identity import (
    IDENTITY_DROP_COLUMNS,
    REGENERATE_ARTIFACT_MESSAGE,
    clean_identity_artifact_table,
    coerce_identity_drop_frame,
    effective_merge_key_series,
    identity_artifact_metadata,
    identity_mode_family,
    is_allele_aware_mode,
    validate_identity_artifact_metadata,
)
from ._kernel.liftover import SumstatsLiftoverRequest, default_liftover_metadata
from ._kernel import sumstats_munger as kernel_munge


LOGGER = logging.getLogger("LDSC.sumstats_munger")
parser = kernel_munge.parser
null_values = kernel_munge.null_values
default_cnames = kernel_munge.default_cnames
describe_cname = kernel_munge.describe_cname
numeric_cols = kernel_munge.numeric_cols
read_header = kernel_munge.read_header
get_cname_map = kernel_munge.get_cname_map
get_compression = kernel_munge.get_compression
clean_header = kernel_munge.clean_header
filter_pvals = kernel_munge.filter_pvals
filter_info = kernel_munge.filter_info
filter_frq = kernel_munge.filter_frq
filter_alleles = kernel_munge.filter_alleles
parse_dat = kernel_munge.parse_dat
process_n = kernel_munge.process_n
p_to_z = kernel_munge.p_to_z
check_median = kernel_munge.check_median
parse_flag_cnames = kernel_munge.parse_flag_cnames
munge_sumstats = kernel_munge.munge_sumstats

_SUMSTATS_OUTPUT_FORMATS = {"parquet", "tsv.gz", "both"}
_RAW_SUMSTATS_FORMATS = {"auto", "plain", "daner-old", "daner-new", "pgc-vcf"}
_SUMSTATS_PARQUET_COMPRESSION = "snappy"


@dataclass(frozen=True)
class RawSumstatsInference:
    """Header-level inference report for one raw summary-statistics file.

    The report is produced by ``--infer-only`` and by the normal munging
    workflow before the kernel is called. It records only safe, format-aware
    decisions: detected raw format, column hints that can be applied without
    changing statistical meaning, INFO-list handling, missing required fields,
    and command-line repair suggestions.
    """

    detected_format: str
    column_hints: dict[str, str] = field(default_factory=dict)
    signed_sumstats_spec: str | None = None
    ignore_columns: tuple[str, ...] = ()
    info_list_columns: tuple[str, ...] = ()
    missing_fields: tuple[str, ...] = ()
    suggested_args: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()

    @property
    def runnable(self) -> bool:
        """Return whether inferred configuration contains all required fields."""
        return len(self.missing_fields) == 0


@dataclass(frozen=True)
class SumstatsTable:
    """In-memory LDSC-ready summary-statistics table.

    Parameters
    ----------
    data : pandas.DataFrame
        Munged table containing at least ``SNP``, ``Z``, and ``N``. Current
        package-written artifacts also include canonical ``CHR`` and ``POS``
        columns, filled with missing values when the raw input had no
        coordinates.
    has_alleles : bool
        Whether allele columns are expected to be present and validated.
    source_path : str or None
        Original data source path, if known.
    trait_name : str or None
        Optional trait label for regression summaries.
    provenance : dict, optional
        Lightweight run metadata retained for debugging and output summaries.
        Default is an empty dict.
    config_snapshot : GlobalConfig or None, optional
        Shared configuration captured when the table was produced in-process or
        recovered from a neighboring ``sumstats.metadata.json`` sidecar. Legacy
        disk artifacts without that sidecar use ``None`` because their original
        munge-time configuration is not recoverable from the table alone.
    """
    data: pd.DataFrame
    has_alleles: bool
    source_path: str | None
    trait_name: str | None
    provenance: dict[str, Any] = field(default_factory=dict)
    config_snapshot: GlobalConfig | None = None

    def validate(self) -> None:
        """Validate the minimum LDSC-ready table contract."""
        required = {"SNP", "Z", "N"}
        missing = required - set(self.data.columns)
        if missing:
            raise ValueError(f"SumstatsTable is missing required columns: {sorted(missing)}")
        if self.has_alleles and not {"A1", "A2"}.issubset(self.data.columns):
            raise ValueError("SumstatsTable.has_alleles=True requires A1 and A2 columns.")

    def snp_identifiers(self) -> pd.Series:
        """Return the active SNP identifiers for this table."""
        mode = _sumstats_table_identifier_mode(self.config_snapshot)
        return effective_merge_key_series(
            self.data,
            mode,
            context=f"sumstats identifiers for {self.source_path or self.trait_name or 'sumstats'}",
        ).astype("string")

    def subset_to(self, snps: set[str] | list[str]) -> "SumstatsTable":
        """Return a copy restricted to ``snps`` while preserving metadata."""
        keep = self.snp_identifiers().isin(set(snps))
        return SumstatsTable(
            data=self.data.loc[keep].reset_index(drop=True),
            has_alleles=self.has_alleles,
            source_path=self.source_path,
            trait_name=self.trait_name,
            provenance=dict(self.provenance),
            config_snapshot=self.config_snapshot,
        )

    def align_to_metadata(self, metadata: pd.DataFrame) -> "SumstatsTable":
        """Inner-join the table to ``metadata`` using the active identifier mode."""
        mode = _sumstats_table_identifier_mode(self.config_snapshot)
        left = pd.DataFrame(
            {
                "_ldsc_sumstats_key": effective_merge_key_series(
                    metadata,
                    mode,
                    context="sumstats metadata alignment left table",
                )
            }
        ).dropna(subset=["_ldsc_sumstats_key"])
        right = self.data.copy()
        right["_ldsc_sumstats_key"] = effective_merge_key_series(
            right,
            mode,
            context=f"sumstats metadata alignment for {self.source_path or self.trait_name or 'sumstats'}",
        )
        right = right.dropna(subset=["_ldsc_sumstats_key"])
        merged = pd.merge(left, right, how="inner", on="_ldsc_sumstats_key", sort=False).drop(
            columns=["_ldsc_sumstats_key"]
        )
        return SumstatsTable(
            data=merged.reset_index(drop=True),
            has_alleles=self.has_alleles,
            source_path=self.source_path,
            trait_name=self.trait_name,
            provenance=dict(self.provenance),
            config_snapshot=self.config_snapshot,
        )

    def summary(self) -> dict[str, Any]:
        """Summarize the retained rows and provenance fields."""
        return {
            "n_rows": len(self.data),
            "has_alleles": self.has_alleles,
            "trait_name": self.trait_name,
            "source_path": self.source_path,
        }


@dataclass(frozen=True)
class MungeRunSummary:
    """Compact summary of one munging run.

    ``output_paths`` records curated sumstats data artifacts and the metadata
    sidecar. It intentionally excludes ``sumstats.log`` so Python result
    contracts stay aligned with other workflow modules.
    """
    n_input_rows: int
    n_retained_rows: int
    drop_counts: dict[str, int]
    inferred_columns: dict[str, str]
    used_n_rule: str
    output_paths: dict[str, str]


def _sumstats_table_identifier_mode(config_snapshot: GlobalConfig | None) -> str:
    """Return the mode used for table-local identity helpers."""
    if config_snapshot is None:
        return "chr_pos"
    return normalize_snp_identifier_mode(config_snapshot.snp_identifier)


def load_sumstats(path: str | PathLike[str], trait_name: str | None = None) -> SumstatsTable:
    """Load one curated LDSC sumstats artifact into a ``SumstatsTable``.

    Parameters
    ----------
    path : str or os.PathLike[str]
        Path token for the curated summary-statistics artifact. This may be a
        literal path or an exact-one glob pattern. Resolution happens at the
        workflow layer before suffix inference. ``.parquet`` files are read
        with :func:`pandas.read_parquet`; ``.sumstats.gz`` and ``.sumstats``
        files are read as whitespace-delimited text. Other suffixes raise a
        clear ``ValueError``.
    trait_name : str or None, optional
        Optional trait label propagated into downstream regression summaries.
        When supplied, this value overrides any metadata sidecar label. When
        omitted, the loader uses ``sumstats.metadata.json`` ``trait_name`` when
        present, then falls back to the resolved filename. Default is ``None``.

    Returns
    -------
    SumstatsTable
        Validated in-memory table with canonical LDSC columns such as ``SNP``,
        ``CHR``, ``POS``, ``N``, and ``Z`` when present in the artifact. When a
        ``*.metadata.json`` sidecar is present next to the artifact, the
        returned table also recovers its munge-time ``GlobalConfig`` snapshot.

    Raises
    ------
    ValueError
        If ``path`` does not resolve to exactly one file or if the curated
        artifact is missing required LDSC columns.

    Notes
    -----
    Format inference is suffix-based after exact-one path resolution:
    ``.parquet`` uses :func:`pandas.read_parquet`, ``.sumstats.gz`` uses a
    gzip-compressed whitespace reader, and ``.sumstats`` uses a plain
    whitespace reader. Other suffixes are refused so callers do not
    accidentally parse CSV or raw GWAS input as curated sumstats.

    Sidecar lookup strips one recognized artifact suffix, so
    ``trait.parquet``, ``trait.sumstats.gz``, and ``trait.sumstats`` all look
    for ``trait.metadata.json``.

    Current package-written artifacts must have a neighboring metadata sidecar
    with the minimal identity provenance contract. Missing or old sidecars are
    rejected so callers regenerate artifacts instead of loading unknown
    identity semantics.
    """
    resolved = resolve_scalar_path(path, label="munged sumstats")
    df = _read_curated_sumstats_artifact(resolved)
    resolved_columns = _resolve_curated_sumstats_columns(list(df.columns), context=resolved)
    df = df.loc[:, list(resolved_columns.values())].rename(
        columns={actual: canonical for canonical, actual in resolved_columns.items()}
    )
    metadata_path = _sumstats_metadata_path(resolved)
    metadata = _read_sumstats_metadata(metadata_path)
    if metadata is None:
        raise ValueError(REGENERATE_ARTIFACT_MESSAGE)
    config_snapshot = _global_config_from_sumstats_metadata(metadata)
    mode = normalize_snp_identifier_mode(config_snapshot.snp_identifier)
    if is_allele_aware_mode(mode) and not {"A1", "A2"}.issubset(df.columns):
        raise ValueError(
            f"Munged sumstats artifact is malformed: snp_identifier='{mode}' requires A1/A2 columns. "
            "Regenerate it with the current LDSC package."
        )
    cleanup = clean_identity_artifact_table(
        df,
        mode,
        context="loaded munged sumstats artifact",
        stage="load_sumstats_validation",
        logger=None,
    )
    if not cleanup.dropped.empty:
        reasons = ", ".join(
            f"{reason}={int(count)}"
            for reason, count in cleanup.dropped["reason"].value_counts(sort=False).items()
        )
        raise ValueError(
            "Munged sumstats artifact is malformed: duplicate or invalid SNP identity rows were found "
            f"({reasons}). Regenerate it with the current LDSC package."
        )
    effective_keys = effective_merge_key_series(
        cleanup.cleaned,
        mode,
        context="loaded munged sumstats artifact",
    )
    if bool(effective_keys.isna().any()):
        raise ValueError(
            "Munged sumstats artifact is malformed: missing or invalid SNP identity rows were found. "
            "Regenerate it with the current LDSC package."
        )
    table = SumstatsTable(
        data=df.reset_index(drop=True),
        has_alleles={"A1", "A2"}.issubset(df.columns),
        source_path=resolved,
        trait_name=_resolve_sumstats_trait_name(trait_name, metadata, resolved),
        provenance={
            **({"metadata_path": str(metadata_path), "metadata": metadata} if metadata is not None else {}),
        },
        config_snapshot=config_snapshot,
    )
    table.validate()
    return table


class SumstatsMunger:
    """Run legacy-compatible munging through a typed workflow interface."""

    def __init__(self) -> None:
        """Initialize the workflow wrapper and clear any cached run summary."""
        self._last_summary: MungeRunSummary | None = None

    def run(
        self,
        raw_sumstats_config: MungeConfig,
        munge_config: MungeConfig | None = None,
        global_config: GlobalConfig | None = None,
    ) -> SumstatsTable:
        """Munge one raw summary-statistics file into LDSC-ready form.

        Parameters
        ----------
        raw_sumstats_config : MungeConfig
            Munging config with raw file path and optional column hints. When
            ``munge_config`` is omitted, this object also supplies output and
            QC settings. Common plain-text, old DANER, new DANER, and PGC
            VCF-style headers are inferred by default before the kernel runs;
            explicit column hints still take priority.
        munge_config : MungeConfig
            Munging thresholds, output directory, and curated output format. The
            workflow writes fixed files named ``sumstats.parquet`` and/or
            ``sumstats.sumstats.gz`` plus ``sumstats.log`` and
            ``sumstats.metadata.json`` inside ``munge_config.output_dir``. Any
            existing owned ``sumstats.*`` artifact is refused before the kernel
            runs unless ``munge_config.overwrite`` is true; successful
            overwrites remove stale sibling formats that the current run did
            not produce. If ``munge_config.sumstats_snps_file`` is supplied, it
            is treated as a headered keep-list, loaded once before raw chunk
            parsing, and applied to each parsed chunk after munging QC and
            coordinate normalization. Allele-free keep-lists match by base key,
            and allele-bearing keep-lists in allele-aware modes match by the
            effective allele-aware key; keep-lists do not reorder output rows.
            In coordinate-family modes, keep-list filtering uses source-build
            coordinates before any optional output liftover.
            ``sumstats_format="auto"`` is the default; use
            ``sumstats_format="plain"``, ``"daner-old"``, ``"daner-new"``, or
            ``"pgc-vcf"`` only when overriding auto-detection.
        global_config : GlobalConfig or None, optional
            Shared configuration snapshot to attach to the returned
            ``SumstatsTable``. When omitted, the current package-global
            configuration is captured. Default is ``None``.

        Returns
        -------
        SumstatsTable
            Validated, in-memory table suitable for the regression workflow.
            The table includes canonical ``CHR`` and ``POS`` columns, preserves
            the active or inferred ``GlobalConfig`` snapshot, and writes the
            same provenance into ``sumstats.metadata.json`` so downstream
            regression can detect incompatible LD-score results after reload.
            Optional ``sumstats_snps_file`` filtering is reflected in both the
            returned table and the written curated artifact(s).
            Output paths for the corresponding disk artifacts are available
            through :meth:`build_run_summary`; the workflow log is written to
            ``sumstats.log`` but is not included in that result mapping.
        """
        if munge_config is None:
            munge_config = raw_sumstats_config
        elif isinstance(munge_config, GlobalConfig) and global_config is None:
            global_config = munge_config
            munge_config = raw_sumstats_config
        if raw_sumstats_config.raw_sumstats_file is None:
            raise ValueError("MungeConfig.raw_sumstats_file is required to run summary-statistics munging.")
        if munge_config.output_dir is None:
            raise ValueError("MungeConfig.output_dir is required to run summary-statistics munging.")

        config_snapshot = global_config or get_global_config()
        liftover_request = _liftover_request_from_config(munge_config)
        _validate_liftover_request_before_io(config_snapshot, liftover_request)
        source_path = resolve_scalar_path(raw_sumstats_config.raw_sumstats_file, label="raw sumstats")
        raw_sumstats_config, munge_config, inference = _apply_raw_sumstats_inference(
            source_path, raw_sumstats_config, munge_config
        )
        output_dir = ensure_output_directory(munge_config.output_dir, label="output directory")
        fixed_output_stem = str(output_dir / "sumstats")
        metadata_path = fixed_output_stem + ".metadata.json"
        output_files = _sumstats_output_files(fixed_output_stem, munge_config.output_format)
        log_path = fixed_output_stem + ".log"
        dropped_snps_path = output_dir / "dropped_snps" / "dropped.tsv.gz"
        sumstats_snps_path = _sumstats_snps_file_from_config(munge_config)
        sumstats_snps_label = "none" if sumstats_snps_path is None else str(sumstats_snps_path)
        produced_paths = [*output_files.values(), metadata_path, log_path, dropped_snps_path]
        owned_paths = [
            *_sumstats_output_files(fixed_output_stem, "both").values(),
            metadata_path,
            log_path,
            dropped_snps_path,
        ]
        # See path_resolution.preflight_output_artifact_family for the
        # owned_paths/produced_paths split. The dropped-SNP sidecar is in both
        # lists because it is always written, including clean header-only runs.
        stale_paths = preflight_output_artifact_family(
            produced_paths,
            owned_paths,
            overwrite=munge_config.overwrite,
            label="munged output artifact",
        )
        args = self._build_args(raw_sumstats_config, munge_config, config_snapshot, liftover_request)
        with workflow_logging("munge-sumstats", log_path, log_level=config_snapshot.log_level):
            log_inputs(
                raw_sumstats_file=source_path,
                output_dir=str(output_dir),
                output_format=munge_config.output_format,
                sumstats_snps_file=sumstats_snps_label,
                use_hm3_snps=munge_config.use_hm3_snps,
                hm3_map_file=packaged_hm3_curated_map_path() if munge_config.use_hm3_snps else "none",
                target_genome_build=liftover_request.target_build or "none",
                liftover_method=liftover_request.method or "none",
            )
            LOGGER.info(
                f"Munging summary statistics from '{source_path}' into '{output_dir}' "
                f"with snp_identifier='{config_snapshot.snp_identifier}', "
                f"genome_build='{config_snapshot.genome_build}', "
                f"sumstats_snps_file='{sumstats_snps_label}', "
                f"use_hm3_snps='{munge_config.use_hm3_snps}', "
                f"target_genome_build='{liftover_request.target_build}', "
                f"liftover_method='{liftover_request.method}'."
            )
            data = kernel_munge.munge_sumstats(args, p=False)
            primary_sumstats_file, parquet_row_groups = _write_sumstats_outputs(
                data,
                output_files=output_files,
                output_format=munge_config.output_format,
            )
            coordinate_metadata = dict(getattr(args, "_coordinate_metadata", {}))
            drop_frame = _coerce_sumstats_dropped_snps_frame(
                coordinate_metadata.pop("liftover_drop_frame", None)
            )
            identity_drop_frame = _coerce_sumstats_dropped_snps_frame(
                getattr(args, "_identity_drop_frame", None),
                default_stage=None,
            )
            drop_frame = _coerce_sumstats_dropped_snps_frame(
                pd.concat([drop_frame, identity_drop_frame], ignore_index=True),
                default_stage=None,
            )
            table_config_snapshot = _effective_sumstats_config(config_snapshot, coordinate_metadata)
            _log_sumstats_provenance(
                coordinate_metadata=coordinate_metadata,
                output_format=munge_config.output_format,
                output_files=output_files,
                parquet_compression=(
                    _SUMSTATS_PARQUET_COMPRESSION
                    if "parquet" in output_files
                    else None
                ),
                parquet_row_groups=parquet_row_groups,
            )
            # Always write the audit sidecar so consumers can distinguish a
            # clean run from a missing or stale dropped-SNP artifact.
            _write_sumstats_dropped_snps_sidecar(drop_frame, dropped_snps_path)
            _log_sumstats_dropped_snps_summary(drop_frame, dropped_snps_path)
            _write_sumstats_metadata(
                metadata_path,
                config_snapshot=table_config_snapshot,
                trait_name=raw_sumstats_config.trait_name,
            )
            table = SumstatsTable(
                data=data.reset_index(drop=True),
                has_alleles={"A1", "A2"}.issubset(data.columns),
                source_path=source_path,
                trait_name=raw_sumstats_config.trait_name,
                provenance={
                    "raw_sumstats_file": source_path,
                    "output_dir": str(output_dir),
                    "output_format": munge_config.output_format,
                    "output_files": dict(output_files),
                    "column_hints": dict(raw_sumstats_config.column_hints),
                    "metadata_path": metadata_path,
                    "coordinate_provenance": coordinate_metadata,
                    "metadata": coordinate_metadata,
                },
                config_snapshot=table_config_snapshot,
            )
            table.validate()
            run_output_paths = {
                **({"sumstats_parquet": output_files["parquet"]} if "parquet" in output_files else {}),
                **({"sumstats_gz": output_files["tsv.gz"]} if "tsv.gz" in output_files else {}),
                "metadata_json": metadata_path,
                "dropped_snps_tsv_gz": str(dropped_snps_path),
            }
            self._last_summary = MungeRunSummary(
                n_input_rows=_count_data_rows(source_path),
                n_retained_rows=len(table.data),
                drop_counts={},
                inferred_columns={**dict(raw_sumstats_config.column_hints), "format": inference.detected_format},
                used_n_rule=_infer_used_n_rule(args),
                output_paths=run_output_paths,
            )
            log_outputs(**run_output_paths)
            LOGGER.info(
                f"Munged {self._last_summary.n_input_rows} input rows to {self._last_summary.n_retained_rows} retained rows; "
                f"wrote '{primary_sumstats_file}'."
            )
            remove_output_artifacts(stale_paths)
        return table

    def write_output(
        self,
        sumstats: SumstatsTable,
        output_dir: str | PathLike[str],
        output_format: str = "parquet",
        overwrite: bool = False,
    ) -> str:
        """Write an in-memory sumstats table to fixed curated artifact names.

        Parameters
        ----------
        sumstats : SumstatsTable
            Validated LDSC-ready table to persist. Columns are curated to the
            package output order ``SNP, CHR, POS, A1, A2, Z, N, FRQ`` where
            present; missing ``CHR`` or ``POS`` columns are materialized as
            missing values.
        output_dir : str or os.PathLike[str]
            Destination directory for fixed ``sumstats`` artifacts.
        output_format : {"parquet", "tsv.gz", "both"}, optional
            Disk format to write. ``"parquet"`` is the default and returns
            ``sumstats.parquet``. ``"tsv.gz"`` writes legacy
            ``sumstats.sumstats.gz``. ``"both"`` writes both and returns the
            Parquet path.
        overwrite : bool, optional
            If ``True``, replace current fixed output artifacts and remove
            stale sibling formats after a successful write. If ``False``, any
            existing owned ``sumstats.*`` artifact is refused before writing
            starts. Default is ``False``.

        Returns
        -------
        str
            Primary sumstats artifact path. The primary path is Parquet for
            ``"parquet"`` and ``"both"``, and gzip TSV for ``"tsv.gz"``.

        Notes
        -----
        This helper follows the public workflow naming policy but does not
        create ``sumstats.log`` because no raw munging kernel is run. A config
        snapshot is required so the neighboring ``sumstats.metadata.json``
        sidecar can be written for downstream provenance recovery.
        """
        if sumstats.config_snapshot is None:
            raise ValueError("SumstatsTable.config_snapshot is required to write current sumstats artifacts.")
        output_format = _normalize_output_format(output_format)
        output_root = ensure_output_directory(output_dir, label="output directory")
        fixed_output_stem = str(output_root / "sumstats")
        metadata_path = fixed_output_stem + ".metadata.json"
        output_files = _sumstats_output_files(fixed_output_stem, output_format)
        produced_paths = list(output_files.values())
        produced_paths.append(metadata_path)
        stale_paths = preflight_output_artifact_family(
            produced_paths,
            [*_sumstats_output_files(fixed_output_stem, "both").values(), metadata_path],
            overwrite=overwrite,
            label="munged output artifact",
        )
        primary_sumstats_file, _parquet_row_groups = _write_sumstats_outputs(
            sumstats.data,
            output_files=output_files,
            output_format=output_format,
        )
        _write_sumstats_metadata(
            metadata_path,
            config_snapshot=sumstats.config_snapshot,
            trait_name=sumstats.trait_name,
        )
        remove_output_artifacts(stale_paths)
        return primary_sumstats_file

    def build_run_summary(self, _sumstats: SumstatsTable | None = None) -> MungeRunSummary:
        """Return the summary captured from the most recent call to :meth:`run`."""
        if self._last_summary is None:
            raise ValueError("No munging run has been executed yet.")
        return self._last_summary

    def _build_args(
        self,
        raw_sumstats_config: MungeConfig,
        munge_config: MungeConfig,
        global_config: GlobalConfig,
        liftover_request: SumstatsLiftoverRequest | None = None,
    ) -> argparse.Namespace:
        """Translate dataclass configuration into the legacy parser namespace."""
        args = parser.parse_args("")
        args.sumstats = resolve_scalar_path(raw_sumstats_config.raw_sumstats_file, label="raw sumstats")
        output_dir = ensure_output_directory(munge_config.output_dir, label="output directory")
        args.out = str(output_dir / "sumstats")
        args.N = munge_config.N
        args.N_cas = munge_config.N_cas
        args.N_con = munge_config.N_con
        args.info_min = munge_config.info_min
        args.maf_min = munge_config.maf_min
        args.n_min = munge_config.n_min
        args.nstudy_min = munge_config.nstudy_min
        args.chunksize = munge_config.chunk_size
        args.sumstats_snps = _resolve_sumstats_snps_path(munge_config)
        args.signed_sumstats = munge_config.signed_sumstats_spec
        args.ignore = ",".join(munge_config.ignore_columns) if munge_config.ignore_columns else None
        args.info_list = ",".join(munge_config.info_list_columns) if munge_config.info_list_columns else None
        args.no_alleles = munge_config.no_alleles
        args.a1_inc = munge_config.a1_inc
        args.keep_maf = munge_config.keep_maf
        args.daner_old = munge_config.daner_old
        args.daner_new = munge_config.daner_new
        args.snp_identifier = global_config.snp_identifier
        args.genome_build = global_config.genome_build
        args._liftover_request = liftover_request or _liftover_request_from_config(munge_config)
        for key, value in raw_sumstats_config.column_hints.items():
            attr = _COLUMN_HINT_ATTRS.get(key)
            if attr is not None:
                setattr(args, attr, value)
        if munge_config.info_list_columns:
            args.info_list = ",".join(munge_config.info_list_columns)
        return args


def run_munge_sumstats_from_args(args: argparse.Namespace) -> SumstatsTable | RawSumstatsInference:
    """Run summary-statistics munging from parsed CLI arguments.

    The CLI path normalizes argparse values into the same ``MungeConfig``
    objects used by the Python API, then delegates to :class:`SumstatsMunger`
    so path resolution, output preflight, metadata, and result construction
    stay in one workflow path.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments from :func:`build_parser`. The namespace must include
        ``raw_sumstats_file`` plus any legacy-compatible munging options copied
        from the kernel parser. ``output_dir`` is required unless
        ``infer_only`` is true.

    Returns
    -------
    SumstatsTable or RawSumstatsInference
        Validated in-memory table produced by :meth:`SumstatsMunger.run`, or a
        header-level inference report when ``infer_only`` is true.

    Raises
    ------
    ValueError
        If required config fields are missing or incompatible with the chosen
        SNP identifier mode.
    FileExistsError
        If fixed output artifacts already exist and ``args.overwrite`` is
        false.
    """
    raw_config, munge_config = _munge_configs_from_args(args)
    _validate_sumstats_format_flags(munge_config)
    if getattr(args, "infer_only", False):
        source_path = resolve_scalar_path(raw_config.raw_sumstats_file, label="raw sumstats")
        inference = infer_raw_sumstats(source_path, raw_config, munge_config)
        print(_render_inference_report(inference, source_path))
        return inference
    return SumstatsMunger().run(raw_config, munge_config, _resolve_main_global_config(args))


def main(argv: list[str] | None = None) -> SumstatsTable | RawSumstatsInference:
    """CLI entry point: parse arguments and delegate to workflow orchestration."""
    return run_munge_sumstats_from_args(build_parser().parse_args(argv))


def _resolve_main_global_config(args: argparse.Namespace) -> GlobalConfig:
    mode = normalize_snp_identifier_mode(getattr(args, "snp_identifier", "chr_pos_allele_aware"))
    if identity_mode_family(mode) == "rsid":
        config = GlobalConfig(
            snp_identifier=mode,
            genome_build=normalize_genome_build(getattr(args, "genome_build", None)),
            log_level=getattr(args, "log_level", "INFO"),
        )
        args.genome_build = config.genome_build
        return config
    genome_build = normalize_genome_build(getattr(args, "genome_build", None))
    if genome_build is None:
        raise ValueError(
            "genome_build is required for chr_pos-family snp_identifier modes. "
            "Pass --genome-build auto, --genome-build hg19, or --genome-build hg38."
        )
    args.genome_build = genome_build
    return GlobalConfig(snp_identifier=mode, genome_build=genome_build, log_level=getattr(args, "log_level", "INFO"))


def _sumstats_snps_file_from_config(config: MungeConfig) -> str | None:
    """Return the explicit or packaged SNP keep-list path for one munger run."""
    if config.use_hm3_snps:
        return packaged_hm3_curated_map_path()
    return config.sumstats_snps_file


def _resolve_sumstats_snps_path(config: MungeConfig) -> str | None:
    """Resolve the effective sumstats SNP restriction path."""
    path = _sumstats_snps_file_from_config(config)
    if path is None:
        return None
    label = "packaged HM3 SNP map" if config.use_hm3_snps else "sumstats SNPs file"
    return resolve_scalar_path(path, label=label)


def _liftover_request_from_config(config: MungeConfig) -> SumstatsLiftoverRequest:
    """Build the kernel liftover request from public munger config."""
    hm3_map_file = packaged_hm3_curated_map_path() if config.use_hm3_quick_liftover else None
    return SumstatsLiftoverRequest(
        target_build=config.target_genome_build,
        liftover_chain_file=config.liftover_chain_file,
        use_hm3_quick_liftover=config.use_hm3_quick_liftover,
        hm3_map_file=hm3_map_file,
    )


def _validate_liftover_request_before_io(config: GlobalConfig, request: SumstatsLiftoverRequest) -> None:
    """Reject liftover requests that can be proven invalid before input IO."""
    if request.requested and identity_mode_family(config.snp_identifier) != "chr_pos":
        raise ValueError("Summary-statistics liftover is only valid for chr_pos-family snp_identifier modes.")
    source_build = normalize_genome_build(config.genome_build)
    if source_build not in {"hg19", "hg38"} or request.target_build is None:
        return
    if request.target_build == source_build and request.method is not None:
        raise ValueError("A liftover method was specified, but target_genome_build equals the source genome build.")
    if request.target_build != source_build and request.method is None:
        raise ValueError("target_genome_build differs from the source genome build, but no liftover method was specified.")


def _munge_configs_from_args(args: argparse.Namespace) -> tuple[MungeConfig, MungeConfig]:
    """Convert parsed CLI arguments into raw-input and run configuration."""
    raw_config = MungeConfig(
        raw_sumstats_file=args.raw_sumstats_file,
        trait_name=getattr(args, "trait_name", None),
        column_hints=_column_hints_from_args(args),
    )
    munge_config = MungeConfig(
        output_dir=args.output_dir,
        N=getattr(args, "N", None),
        N_cas=getattr(args, "N_cas", None),
        N_con=getattr(args, "N_con", None),
        info_min=args.info_min,
        maf_min=args.maf_min,
        n_min=getattr(args, "n_min", None),
        nstudy_min=getattr(args, "nstudy_min", None),
        chunk_size=args.chunksize,
        output_format=args.output_format,
        sumstats_snps_file=getattr(args, "sumstats_snps_file", None),
        use_hm3_snps=getattr(args, "use_hm3_snps", False),
        target_genome_build=getattr(args, "target_genome_build", None),
        liftover_chain_file=getattr(args, "liftover_chain_file", None),
        use_hm3_quick_liftover=getattr(args, "use_hm3_quick_liftover", False),
        signed_sumstats_spec=getattr(args, "signed_sumstats", None),
        ignore_columns=_ignore_columns_from_args(args),
        info_list_columns=_info_list_columns_from_args(args),
        sumstats_format=getattr(args, "sumstats_format", "auto"),
        no_alleles=args.no_alleles,
        a1_inc=args.a1_inc,
        keep_maf=args.keep_maf,
        daner_old=args.daner_old,
        daner_new=args.daner_new,
        overwrite=getattr(args, "overwrite", False),
    )
    return raw_config, munge_config


def _column_hints_from_args(args: argparse.Namespace) -> dict[str, str]:
    """Collect explicit raw-column hints from parsed CLI arguments."""
    hints: dict[str, str] = {}
    for key in _COLUMN_HINT_ARG_KEYS:
        value = getattr(args, _COLUMN_HINT_ATTRS[key], None)
        if value is not None:
            hints[key] = value
    return hints


def _ignore_columns_from_args(args: argparse.Namespace) -> tuple[str, ...]:
    """Return normalized ``--ignore`` column tokens from parsed CLI arguments."""
    ignore = getattr(args, "ignore", None)
    if not ignore:
        return ()
    return tuple(token.strip() for token in ignore.split(",") if token.strip())


def _info_list_columns_from_args(args: argparse.Namespace) -> tuple[str, ...]:
    """Return normalized ``--info-list`` column tokens from parsed CLI arguments."""
    info_list = getattr(args, "info_list", None)
    if not info_list:
        return ()
    return tuple(token.strip() for token in info_list.split(",") if token.strip())


def kernel_parser():
    """Expose the public munging parser for CLI action cloning."""
    return build_parser()


def build_parser() -> argparse.ArgumentParser:
    """Build the public summary-statistics munging parser."""
    public = argparse.ArgumentParser(description=getattr(parser, "description", None), allow_abbrev=False)
    public.add_argument("--raw-sumstats-file", required=True, help="Raw summary-statistics file path.")
    public.add_argument("--output-dir", default=None, help="Output directory for munged sumstats and logs.")
    public.add_argument(
        "--format",
        dest="sumstats_format",
        choices=sorted(_RAW_SUMSTATS_FORMATS),
        default="auto",
        help="Raw summary-statistics format profile. Default is auto.",
    )
    public.add_argument(
        "--infer-only",
        action="store_true",
        default=False,
        help="Inspect the raw file, print inferred columns and a suggested command, and write no outputs.",
    )
    public.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Replace sumstats output artifacts and remove stale owned siblings.",
    )
    public.add_argument("--sumstats-snps-file", default=None, help="Optional SNP keep-list for munged summary statistics.")
    public.add_argument(
        "--use-hm3-snps",
        action="store_true",
        default=False,
        help="Restrict munged summary statistics to the packaged curated HM3 SNP map.",
    )
    public.add_argument(
        "--target-genome-build",
        default=None,
        choices=("hg19", "hg37", "GRCh37", "hg38", "GRCh38"),
        help="Optional target genome build for chr_pos-family output coordinates.",
    )
    public.add_argument(
        "--liftover-chain-file",
        default=None,
        help="Optional chain file used to liftover chr_pos-family sumstats coordinates to --target-genome-build.",
    )
    public.add_argument(
        "--use-hm3-quick-liftover",
        action="store_true",
        default=False,
        help="Use the packaged curated dual-build HM3 map for coordinate-only quick liftover.",
    )
    public.add_argument("--trait-name", default=None, help="Optional biological trait label stored in sumstats metadata.")
    public.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"), help="Logging verbosity.")
    public.add_argument(
        "--output-format",
        choices=sorted(_SUMSTATS_OUTPUT_FORMATS),
        default="parquet",
        help="Curated sumstats output format. Default is parquet.",
    )
    for action in parser._actions:
        if action.dest in {"help", "sumstats", "out"}:
            continue
        option_strings = list(action.option_strings)
        if not option_strings:
            continue
        kwargs = {
            "dest": action.dest,
            "default": action.default,
            "required": action.required,
            "help": action.help,
        }
        if getattr(action, "choices", None) is not None:
            kwargs["choices"] = action.choices
        if getattr(action, "type", None) is not None:
            kwargs["type"] = action.type
        if getattr(action, "nargs", None) is not None:
            kwargs["nargs"] = action.nargs
        if action.const is not None:
            kwargs["const"] = action.const
        if action.__class__.__name__ == "_StoreTrueAction":
            public.add_argument(*option_strings, action="store_true", default=action.default, help=action.help)
        elif action.__class__.__name__ == "_StoreFalseAction":
            public.add_argument(*option_strings, action="store_false", default=action.default, help=action.help)
        else:
            public.add_argument(*option_strings, **kwargs)
    return public


def _count_data_rows(path: str) -> int:
    """Count raw data rows after leading ``##`` metadata lines and the header."""
    openfunc, _compression = kernel_munge.get_compression(path)
    skiprows = kernel_munge.count_leading_sumstats_comment_lines(path)
    with openfunc(path) as handle:
        count = sum(1 for _ in handle)
    return max(count - skiprows - 1, 0)


def _infer_used_n_rule(args: argparse.Namespace) -> str:
    """Summarize which sample-size rule the current munging run will apply."""
    if getattr(args, "N", None) is not None:
        return "fixed_N"
    if getattr(args, "N_cas", None) is not None and getattr(args, "N_con", None) is not None:
        return "fixed_case_control_N"
    return "input_columns"


def infer_raw_sumstats(
    raw_sumstats_file: str | PathLike[str],
    raw_config: MungeConfig | None = None,
    munge_config: MungeConfig | None = None,
) -> RawSumstatsInference:
    """Infer raw summary-statistics format and minimal parser hints.

    The inference pass reads the header and first data row only. It detects
    plain text, old DANER, new DANER, and PGC VCF-style inputs; reports missing
    required fields; identifies numeric/NA comma-separated INFO lists; and
    returns exact CLI hints for cases that still need user confirmation. It
    intentionally does not treat ``NEFF`` as total sample size ``N``.
    """
    raw_config = raw_config or MungeConfig(raw_sumstats_file=raw_sumstats_file)
    munge_config = munge_config or MungeConfig()
    path = resolve_scalar_path(raw_sumstats_file, label="raw sumstats")
    file_cnames = read_header(path)
    clean_headers = [clean_header(column) for column in file_cnames]
    clean_set = set(clean_headers)
    sample = _read_first_data_row(path)
    requested_format = munge_config.sumstats_format
    detected_format = _detect_sumstats_format(path, clean_set, requested_format)
    column_hints: dict[str, str] = {}
    info_list_columns: list[str] = list(munge_config.info_list_columns)
    ignore_columns: list[str] = list(munge_config.ignore_columns)
    suggested_args: list[str] = []
    notes: list[str] = []

    if detected_format == "daner-old":
        suggested_args.extend(["--format", "daner-old"])
    elif detected_format == "daner-new":
        suggested_args.extend(["--format", "daner-new"])
        frq_u_column = _first_column_with_clean_prefix(file_cnames, "FRQ_U_")
        if frq_u_column is not None:
            column_hints.setdefault("frq", frq_u_column)
    elif detected_format == "pgc-vcf":
        suggested_args.extend(["--format", "pgc-vcf"])
        if "REF" in clean_set and "ALT" in clean_set and not {"A1", "A2", "EA", "NEA"} & clean_set:
            column_hints.setdefault("a1", _original_column(file_cnames, "REF"))
            column_hints.setdefault("a2", _original_column(file_cnames, "ALT"))

    for column in file_cnames:
        if default_cnames.get(clean_header(column)) == "INFO" and _sample_value_is_info_list(sample.get(column)):
            info_list_columns.append(column)
            notes.append(f"{column} appears to contain comma-separated per-study INFO values.")

    translated = _inferred_targets(file_cnames, {**column_hints, **raw_config.column_hints})
    old_daner_n = detected_format == "daner-old"
    has_case_control_n = {"N_CAS", "N_CON"}.issubset(translated)
    has_total_n = "N" in translated or munge_config.N is not None or (
        munge_config.N_cas is not None and munge_config.N_con is not None
    )
    missing: list[str] = []
    if not (has_total_n or has_case_control_n or old_daner_n):
        missing.append("N")
        if "NEFF" in clean_set:
            notes.append("NEFF is not treated as N automatically; pass --N-col NEFF only if appropriate.")

    has_signed = bool(translated & set(null_values)) or munge_config.signed_sumstats_spec is not None
    signed_sumstats_spec = None
    if not has_signed and not munge_config.a1_inc:
        missing.append("signed statistic")
        signed_column = _likely_signed_sumstat_column(file_cnames)
        if signed_column is not None:
            signed_sumstats_spec = f"{signed_column},0"
            suggested_args.extend(["--signed-sumstats", signed_sumstats_spec])
            notes.append(
                f"{signed_column} may be a signed effect column; pass --signed-sumstats {signed_sumstats_spec} "
                "if its null value is 0 and it is oriented relative to A1."
            )

    has_alleles = {"A1", "A2"}.issubset(translated) or munge_config.no_alleles
    if not has_alleles:
        missing.append("A1/A2")
        if "REF" in clean_set and "ALT" in clean_set:
            notes.append(
                "REF and ALT are present; pass --a1 REF --a2 ALT only if the signed statistic is relative to REF."
            )

    return RawSumstatsInference(
        detected_format=detected_format,
        column_hints=column_hints,
        signed_sumstats_spec=signed_sumstats_spec,
        ignore_columns=tuple(dict.fromkeys(ignore_columns)),
        info_list_columns=tuple(dict.fromkeys(info_list_columns)),
        missing_fields=tuple(missing),
        suggested_args=tuple(suggested_args),
        notes=tuple(notes),
    )


def _read_first_data_row(path: str) -> dict[str, str]:
    """Return the first data row keyed by raw header names."""
    file_cnames = read_header(path)
    openfunc, _compression = get_compression(path)
    skiprows = kernel_munge.count_leading_sumstats_comment_lines(path)
    with openfunc(path) as handle:
        for _idx in range(skiprows + 1):
            handle.readline()
        for line in handle:
            if line.strip():
                values = line.split()
                return {column: values[idx] for idx, column in enumerate(file_cnames) if idx < len(values)}
    return {}


def _detect_sumstats_format(path: str, clean_headers: set[str], requested_format: str) -> str:
    if requested_format != "auto":
        return requested_format
    has_old_daner_n = any(column.startswith("FRQ_A_") for column in clean_headers) and any(
        column.startswith("FRQ_U_") for column in clean_headers
    )
    if has_old_daner_n:
        return "daner-old"
    if kernel_munge.count_leading_sumstats_comment_lines(path) > 0 and "#CHROM" in clean_headers:
        return "pgc-vcf"
    if {"NCA", "NCO"}.issubset(clean_headers) or {"NCAS", "NCON"}.issubset(clean_headers):
        return "daner-new"
    return "plain"


def _original_column(file_cnames: list[str], cleaned: str) -> str:
    for column in file_cnames:
        if clean_header(column) == cleaned:
            return column
    return cleaned


def _first_column_with_clean_prefix(file_cnames: list[str], prefix: str) -> str | None:
    for column in file_cnames:
        if clean_header(column).startswith(prefix):
            return column
    return None


def _likely_signed_sumstat_column(file_cnames: list[str]) -> str | None:
    likely = {"EFFECT_SIZE", "EFFECTSIZE", "LOGOR", "LOG_OR", "BETA_HAT"}
    for column in file_cnames:
        if clean_header(column) in likely:
            return column
    return None


def _sample_value_is_info_list(value: str | None) -> bool:
    return value is not None and "," in str(value)


def _inferred_targets(file_cnames: list[str], column_hints: dict[str, str]) -> set[str]:
    targets = {default_cnames[clean_header(column)] for column in file_cnames if clean_header(column) in default_cnames}
    for key, value in column_hints.items():
        attr = _COLUMN_HINT_ATTRS.get(key)
        if value is None or attr is None:
            continue
        target = {
            "snp": "SNP",
            "chr": "CHR",
            "pos": "POS",
            "N_col": "N",
            "N_cas_col": "N_CAS",
            "N_con_col": "N_CON",
            "a1": "A1",
            "a2": "A2",
            "p": "P",
            "frq": "FRQ",
            "info": "INFO",
            "nstudy": "NSTUDY",
        }.get(attr)
        if target is not None:
            targets.add(target)
    return targets


def _apply_raw_sumstats_inference(
    source_path: str,
    raw_config: MungeConfig,
    munge_config: MungeConfig,
) -> tuple[MungeConfig, MungeConfig, RawSumstatsInference]:
    """Return configs updated with auto-inferred raw-format hints."""
    _validate_sumstats_format_flags(munge_config)
    inference = infer_raw_sumstats(source_path, raw_config, munge_config)
    column_hints = {**inference.column_hints, **raw_config.column_hints}
    daner_old = munge_config.daner_old or inference.detected_format == "daner-old"
    daner_new = munge_config.daner_new or munge_config.sumstats_format == "daner-new"
    raw_config = replace(raw_config, column_hints=column_hints)
    munge_config = replace(
        munge_config,
        daner_old=daner_old,
        daner_new=daner_new,
        info_list_columns=tuple(dict.fromkeys([*munge_config.info_list_columns, *inference.info_list_columns])),
        ignore_columns=tuple(dict.fromkeys([*munge_config.ignore_columns, *inference.ignore_columns])),
    )
    return raw_config, munge_config, inference


def _validate_sumstats_format_flags(munge_config: MungeConfig) -> None:
    fmt = munge_config.sumstats_format
    if fmt == "daner-old" and munge_config.daner_new:
        raise ValueError("--format daner-old conflicts with --daner-new.")
    if fmt == "daner-new" and munge_config.daner_old:
        raise ValueError("--format daner-new conflicts with --daner-old.")
    if fmt in {"plain", "pgc-vcf"} and (munge_config.daner_old or munge_config.daner_new):
        raise ValueError(f"--format {fmt} conflicts with DANER-specific flags.")


def _render_inference_report(inference: RawSumstatsInference, raw_sumstats_file: str) -> str:
    lines = [
        f"Raw sumstats file: {raw_sumstats_file}",
        f"Detected format: {inference.detected_format}",
        f"Runnable: {'yes' if inference.runnable else 'no'}",
        "Missing fields: " + (", ".join(inference.missing_fields) if inference.missing_fields else "none"),
    ]
    if inference.column_hints:
        lines.append("Column hints: " + ", ".join(f"{key}={value}" for key, value in sorted(inference.column_hints.items())))
    if inference.signed_sumstats_spec is not None:
        lines.append(f"Signed statistic hint: {inference.signed_sumstats_spec}")
    if inference.info_list_columns:
        lines.append("INFO list columns: " + ", ".join(inference.info_list_columns))
    if inference.notes:
        lines.extend(f"Note: {note}" for note in inference.notes)
    command = ["ldsc", "munge-sumstats", "--raw-sumstats-file", raw_sumstats_file, "--output-dir", "<out>"]
    command.extend(inference.suggested_args)
    lines.append("Suggested command: " + " ".join(command))
    return "\n".join(lines)


def _normalize_output_format(output_format: str) -> str:
    """Return a validated curated sumstats output-format token."""
    if output_format not in _SUMSTATS_OUTPUT_FORMATS:
        raise ValueError("output_format must be one of 'parquet', 'tsv.gz', or 'both'.")
    return output_format


def _sumstats_output_files(fixed_output_stem: str, output_format: str) -> dict[str, str]:
    """Return selected sumstats artifact paths keyed by format token."""
    output_format = _normalize_output_format(output_format)
    paths: dict[str, str] = {}
    if output_format in {"parquet", "both"}:
        paths["parquet"] = fixed_output_stem + ".parquet"
    if output_format in {"tsv.gz", "both"}:
        paths["tsv.gz"] = fixed_output_stem + ".sumstats.gz"
    return paths


def _prepare_curated_sumstats_frame(data: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with fixed sumstats output columns and coordinate placeholders."""
    frame = data.copy()
    if "CHR" not in frame.columns:
        frame["CHR"] = pd.NA
    if "POS" not in frame.columns:
        frame["POS"] = pd.NA
    columns = [col for col in ("SNP", "CHR", "POS", "A1", "A2", "Z", "N", "FRQ") if col in frame.columns]
    return frame.loc[:, columns]


def _prepare_sumstats_parquet_frame(data: pd.DataFrame) -> pd.DataFrame:
    """Return a precision-preserving frame sorted for chromosome row groups.

    Complete-coordinate rows are normalized and ordered by chromosome sort key,
    ``POS``, and original row order. Rows without a complete ``CHR``/``POS``
    pair are kept, preserve original relative order, and sort after all
    complete-coordinate rows so the writer can place them in one final
    missing-coordinate row group.
    """
    frame = _prepare_curated_sumstats_frame(data)
    frame["_ldsc_original_order"] = range(len(frame))
    chr_missing = coordinate_missing_mask(frame["CHR"])
    pos_missing = coordinate_missing_mask(frame["POS"])
    pos_numeric = pd.to_numeric(frame["POS"], errors="coerce")
    invalid_pos = (~pos_missing) & pos_numeric.isna()
    if invalid_pos.any():
        bad_value = frame.loc[invalid_pos, "POS"].iloc[0]
        raise ValueError(f"POS values must be numeric base-pair positions; got {bad_value!r}.")

    complete = ~(chr_missing | pos_missing)
    if complete.any():
        complete_pos = positive_int_position_series(
            frame.loc[complete, "POS"],
            context="sumstats parquet output",
            label="POS",
        )
        frame.loc[complete, "CHR"] = normalize_chromosome_series(
            frame.loc[complete, "CHR"],
            context="sumstats parquet output",
        ).astype(object)
        frame.loc[complete, "POS"] = complete_pos.astype("int64")

    frame["_ldsc_missing_coordinate"] = ~complete
    frame["_ldsc_pos_sort"] = pos_numeric.where(complete, pd.NA)
    frame["_ldsc_chrom_rank"] = 10_000
    if complete.any():
        unique_chroms = pd.unique(frame.loc[complete, "CHR"])
        rank_map = {chrom: chrom_sort_key(chrom)[1] for chrom in unique_chroms}
        frame.loc[complete, "_ldsc_chrom_rank"] = frame.loc[complete, "CHR"].map(rank_map).astype("int64")

    frame = frame.sort_values(
        by=["_ldsc_missing_coordinate", "_ldsc_chrom_rank", "_ldsc_pos_sort", "_ldsc_original_order"],
        kind="mergesort",
    ).reset_index(drop=True)
    return frame.drop(
        columns=["_ldsc_original_order", "_ldsc_missing_coordinate", "_ldsc_pos_sort", "_ldsc_chrom_rank"]
    )


def _write_sumstats_tsv_gz(data: pd.DataFrame, path: str) -> None:
    """Write the legacy-compatible gzip TSV artifact with rounded floats."""
    frame = _prepare_curated_sumstats_frame(data)
    frame.to_csv(path, sep="\t", index=False, float_format="%.3f", compression="gzip")


def _write_sumstats_parquet(data: pd.DataFrame, path: str) -> list[dict[str, Any]]:
    """Write snappy-compressed Parquet and return row-group metadata.

    The Parquet payload keeps the munger's numeric precision. One row group is
    emitted per normalized chromosome among complete-coordinate rows. If rows
    without complete coordinates exist, they are emitted as the last row group
    with ``chrom`` recorded as ``None``.
    """
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise LDSCDependencyError("Writing sumstats parquet artifacts requires pyarrow.") from exc

    frame = _prepare_sumstats_parquet_frame(data)
    schema = pa.Schema.from_pandas(frame, preserve_index=False)
    row_groups: list[dict[str, Any]] = []
    offset = 0
    with pq.ParquetWriter(path, schema, compression=_SUMSTATS_PARQUET_COMPRESSION) as writer:
        if frame.empty:
            writer.write_table(pa.Table.from_pandas(frame, schema=schema, preserve_index=False))
            return row_groups
        complete = ~(coordinate_missing_mask(frame["CHR"]) | coordinate_missing_mask(frame["POS"]))
        for chrom, chrom_df in frame.loc[complete].groupby("CHR", sort=False):
            writer.write_table(pa.Table.from_pandas(chrom_df, schema=schema, preserve_index=False))
            row_groups.append(
                {
                    "chrom": str(chrom),
                    "row_group_index": len(row_groups),
                    "row_offset": offset,
                    "n_rows": len(chrom_df),
                }
            )
            offset += len(chrom_df)
        missing_df = frame.loc[~complete]
        if not missing_df.empty:
            writer.write_table(pa.Table.from_pandas(missing_df, schema=schema, preserve_index=False))
            row_groups.append(
                {
                    "chrom": None,
                    "row_group_index": len(row_groups),
                    "row_offset": offset,
                    "n_rows": len(missing_df),
                }
            )
    return row_groups


def _write_sumstats_outputs(
    data: pd.DataFrame,
    *,
    output_files: dict[str, str],
    output_format: str,
) -> tuple[str, list[dict[str, Any]]]:
    """Write selected artifacts and return the primary path plus Parquet row groups."""
    output_format = _normalize_output_format(output_format)
    parquet_row_groups: list[dict[str, Any]] = []
    if "tsv.gz" in output_files:
        _write_sumstats_tsv_gz(data, output_files["tsv.gz"])
    if "parquet" in output_files:
        parquet_row_groups = _write_sumstats_parquet(data, output_files["parquet"])
    primary = output_files["parquet"] if output_format in {"parquet", "both"} else output_files["tsv.gz"]
    return primary, parquet_row_groups


def _empty_sumstats_dropped_snps_frame() -> pd.DataFrame:
    """Return the canonical empty sumstats dropped-SNP sidecar frame."""
    return pd.DataFrame(
        {
            "CHR": pd.Series(dtype="string"),
            "SNP": pd.Series(dtype="string"),
            "source_pos": pd.Series(dtype="Int64"),
            "target_pos": pd.Series(dtype="Int64"),
            "reason": pd.Series(dtype="string"),
            "base_key": pd.Series(dtype="string"),
            "identity_key": pd.Series(dtype="string"),
            "allele_set": pd.Series(dtype="string"),
            "stage": pd.Series(dtype="string"),
        },
        columns=IDENTITY_DROP_COLUMNS,
    )


def _coerce_sumstats_dropped_snps_frame(
    frame: pd.DataFrame | None,
    *,
    default_stage: str | None = "liftover",
) -> pd.DataFrame:
    """Return dropped-SNP rows with the unified nullable sidecar schema."""
    if frame is None or frame.empty:
        return _empty_sumstats_dropped_snps_frame()
    output = coerce_identity_drop_frame(frame)
    if default_stage is not None:
        output["stage"] = output["stage"].fillna(default_stage)
    output["CHR"] = output["CHR"].astype("string")
    output["SNP"] = output["SNP"].astype("string")
    output["source_pos"] = pd.to_numeric(output["source_pos"], errors="coerce").astype("Int64")
    output["target_pos"] = pd.to_numeric(output["target_pos"], errors="coerce").astype("Int64")
    output["reason"] = output["reason"].astype("string")
    output["base_key"] = output["base_key"].astype("string")
    output["identity_key"] = output["identity_key"].astype("string")
    output["allele_set"] = output["allele_set"].astype("string")
    output["stage"] = output["stage"].astype("string")
    return output.reset_index(drop=True)


def _write_sumstats_dropped_snps_sidecar(drop_frame: pd.DataFrame, path: Path) -> None:
    """Write the always-owned dropped-SNP audit sidecar, even when header-only."""
    path.parent.mkdir(parents=True, exist_ok=True)
    drop_frame.to_csv(path, sep="\t", index=False, compression="gzip", na_rep="")


def _log_sumstats_dropped_snps_summary(drop_frame: pd.DataFrame, path: Path) -> None:
    """Log a count-only summary of the sumstats dropped-SNP sidecar."""
    if drop_frame.empty:
        LOGGER.info(f"No SNPs dropped during liftover or identity cleanup stages; audit sidecar at '{path}'.")
        return
    counts = drop_frame["reason"].value_counts(sort=False)
    count_text = ", ".join(f"{reason}={int(count)}" for reason, count in counts.items())
    LOGGER.info(
        f"Summary-statistics liftover/identity cleanup drops: {len(drop_frame)} SNPs "
        f"({count_text}); audit sidecar at '{path}'."
    )


def _read_curated_sumstats_artifact(path: str) -> pd.DataFrame:
    """Read one curated sumstats artifact according to the public suffix policy."""
    token = str(path)
    if token.endswith(".parquet"):
        return pd.read_parquet(path)
    if token.endswith(".sumstats.gz"):
        return pd.read_csv(path, sep=r"\s+", compression="gzip")
    if token.endswith(".sumstats"):
        return pd.read_csv(path, sep=r"\s+", compression="infer")
    raise ValueError(
        "Unsupported munged sumstats format. Expected a path ending in "
        "'.parquet', '.sumstats.gz', or '.sumstats'."
    )


def _resolve_curated_sumstats_columns(columns: list[str], *, context: str) -> dict[str, str]:
    """Resolve canonical internal sumstats artifact columns from ``columns``."""
    resolved = {
        canonical: resolve_required_column(columns, INTERNAL_SUMSTATS_ARTIFACT_SPEC_MAP[canonical], context=context)
        for canonical in ("SNP", "N", "Z")
    }
    for canonical in ("CHR", "POS", "A1", "A2", "FRQ"):
        actual = resolve_optional_column(columns, INTERNAL_SUMSTATS_ARTIFACT_SPEC_MAP[canonical], context=context)
        if actual is not None:
            resolved[canonical] = actual
    return {canonical: resolved[canonical] for canonical in ("SNP", "CHR", "POS", "N", "Z", "A1", "A2", "FRQ") if canonical in resolved}


def _sumstats_metadata_path(path: str | PathLike[str]) -> Path:
    """Return ``<stem>.metadata.json`` for a recognized sumstats artifact suffix."""
    token = str(path)
    if token.endswith(".parquet"):
        return Path(token[: -len(".parquet")] + ".metadata.json")
    if token.endswith(".sumstats.gz"):
        return Path(token[: -len(".sumstats.gz")] + ".metadata.json")
    if token.endswith(".sumstats"):
        return Path(token[: -len(".sumstats")] + ".metadata.json")
    return Path(token).with_suffix(Path(token).suffix + ".metadata.json")


def _read_sumstats_metadata(path: Path) -> dict[str, Any] | None:
    """Read a sumstats sidecar metadata file when it exists."""
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(REGENERATE_ARTIFACT_MESSAGE) from exc


def _resolve_sumstats_trait_name(
    explicit_trait_name: str | None,
    metadata: dict[str, Any] | None,
    resolved_path: str,
) -> str:
    """Resolve the public trait label using CLI/API, sidecar, then filename."""
    explicit = _normalize_trait_name(explicit_trait_name)
    if explicit is not None:
        return explicit
    if isinstance(metadata, dict) and "trait_name" in metadata:
        metadata_trait = _normalize_trait_name(metadata.get("trait_name"))
        if metadata_trait is not None:
            return metadata_trait
    return Path(resolved_path).name


def _global_config_from_sumstats_metadata(metadata: dict[str, Any] | None) -> GlobalConfig | None:
    """Recreate a GlobalConfig snapshot from a sumstats sidecar."""
    if not isinstance(metadata, dict):
        raise ValueError(REGENERATE_ARTIFACT_MESSAGE)
    if "genome_build" not in metadata:
        raise ValueError(REGENERATE_ARTIFACT_MESSAGE)
    try:
        mode = validate_identity_artifact_metadata(metadata, expected_artifact_type="sumstats")
        return GlobalConfig(
            snp_identifier=mode,
            genome_build=metadata.get("genome_build"),
            log_level="INFO",
        )
    except Exception as exc:
        if isinstance(exc, ValueError) and str(exc) == REGENERATE_ARTIFACT_MESSAGE:
            raise
        raise ValueError("Sumstats metadata sidecar has invalid identity provenance.") from exc


def _effective_sumstats_config(config: GlobalConfig, coordinate_metadata: dict[str, Any]) -> GlobalConfig:
    """Return the config snapshot implied by coordinate provenance."""
    genome_build = coordinate_metadata.get("genome_build") or config.genome_build
    return GlobalConfig(
        snp_identifier=normalize_snp_identifier_mode(coordinate_metadata.get("snp_identifier", config.snp_identifier)),
        genome_build=normalize_genome_build(genome_build),
        log_level=config.log_level,
        fail_on_missing_metadata=config.fail_on_missing_metadata,
    )


def _write_sumstats_metadata(
    path: str | PathLike[str],
    *,
    config_snapshot: GlobalConfig,
    trait_name: str | None,
) -> None:
    """Write the thin metadata sidecar used for downstream compatibility."""
    payload = {
        **identity_artifact_metadata(
            artifact_type="sumstats",
            snp_identifier=config_snapshot.snp_identifier,
            genome_build=config_snapshot.genome_build,
        ),
        "trait_name": trait_name,
    }
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _log_sumstats_provenance(
    *,
    coordinate_metadata: dict[str, Any],
    output_format: str,
    output_files: dict[str, str],
    parquet_compression: str | None,
    parquet_row_groups: list[dict[str, Any]],
) -> None:
    """Log detailed provenance that is intentionally excluded from the sidecar."""
    coordinate_provenance = dict(coordinate_metadata)
    liftover = coordinate_provenance.pop(
        "liftover",
        default_liftover_metadata(
            source_build=coordinate_provenance.get("genome_build"),
            snp_identifier=coordinate_provenance.get("snp_identifier", "chr_pos_allele_aware"),
        ),
    )
    LOGGER.info(
        "Summary-statistics output bookkeeping: output_format=%s; output_files=%s; "
        "parquet_compression=%s; parquet_row_groups=%d",
        output_format,
        _format_log_mapping(dict(output_files)),
        parquet_compression or "none",
        len(parquet_row_groups),
    )
    LOGGER.info(
        "Summary-statistics coordinate provenance: %s",
        _format_log_mapping(coordinate_provenance),
    )
    if isinstance(liftover, dict):
        LOGGER.info("Summary-statistics liftover report: %s", _format_log_mapping(liftover))
        if liftover.get("method") == "hm3_curated":
            LOGGER.info(
                "Summary-statistics HM3 liftover provenance: method=%s; hm3_map_file=%s",
                liftover.get("method"),
                liftover.get("hm3_map_file"),
            )


def _format_log_mapping(values: dict[str, Any]) -> str:
    """Format a small mapping as stable human-readable key/value text."""
    if not values:
        return "none"
    return "; ".join(f"{key}={_format_log_value(values[key])}" for key in sorted(values))


def _format_log_value(value: Any) -> str:
    """Format one provenance value without emitting JSON payloads."""
    if isinstance(value, dict):
        return "{" + ", ".join(f"{key}={_format_log_value(value[key])}" for key in sorted(value)) + "}"
    if isinstance(value, list):
        return "[" + ", ".join(_format_log_value(item) for item in value) + "]"
    if value is None:
        return "none"
    return str(value)


_COLUMN_HINT_ATTRS = {
    "snp": "snp",
    "chr": "chr",
    "pos": "pos",
    "N": "N_col",
    "N_col": "N_col",
    "N_cas": "N_cas_col",
    "N_cas_col": "N_cas_col",
    "N_con": "N_con_col",
    "N_con_col": "N_con_col",
    "a1": "a1",
    "a2": "a2",
    "p": "p",
    "frq": "frq",
    "signed_sumstats": "signed_sumstats",
    "info": "info",
    "info_list": "info_list",
    "nstudy": "nstudy",
}
_COLUMN_HINT_ARG_KEYS = (
    "snp",
    "chr",
    "pos",
    "N_col",
    "N_cas_col",
    "N_con_col",
    "a1",
    "a2",
    "p",
    "frq",
    "info",
    "info_list",
    "nstudy",
)
