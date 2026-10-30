"""Workflow wrapper for LDSC3 summary-statistics munging.

Core functionality:
    Expose a typed, package-level interface to the legacy-compatible munging
    implementation and provide a public loader for munged Parquet or legacy
    ``.sumstats.gz`` artifacts.

Overview
--------
This module converts the historical munging script behavior into explicit
Python dataclasses and a small service object. The public workflow boundary
accepts path-like inputs and resolves the schema, sample-size choices, source
build and keep-list through ``_sumstats.munge_input.prepare_munge_input``. It passes
a ``ResolvedMungeInput`` to the kernel and consumes an explicit ``MungeResult``.
Packaged HM3 restriction is the default; a custom list replaces it, and
``no_snp_restriction`` disables keep-list filtering. For different coordinate
builds, packaged HM3 uses quick liftover unless an explicit chain overrides it;
custom-list and unrestricted runs require a chain.
The workflow layer owns CLI orchestration, output
preflight, the self-describing Parquet footer, diagnostics, and
result objects; the kernel keeps the legacy-compatible parsing and filtering
primitives. Coordinate-family runs separate raw source-build interpretation
from final output-build compatibility: ``source_genome_build="auto"`` infers
the raw ``CHR``/``POS`` build, while ``output_genome_build`` is the required
build recorded in the parquet footer after any liftover. rsID-family runs reject
build and liftover fields and record ``genome_build=None``. Run summaries
expose munged data artifacts only; ``diagnostics/sumstats.log`` is an audit
file and is not included in ``output_paths``. The Parquet footer
carries only the thin downstream-identity payload; source-build and liftover
provenance is written as readable workflow-log text. Successful CLI runs print
the resolved method and mapping/drop counts to stdout and record the same
summary in the log at every log level. Python API calls record it only in the
log. Trait labels determine sanitized data filenames; without a label, use
``sumstats.parquet`` and optional ``sumstats.gz``.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field, replace
import logging
from os import PathLike
from pathlib import Path
import re
import shlex
from typing import Any

import pandas as pd

from .._runtime.cli_help import CLIHelpFormatter, SCALAR_PATH_HELP
from .._runtime.workflow_logging import LOG_LEVEL_HELP
from .._core.chromosomes import chrom_sort_key, normalize_chromosome_series
from .._core.coordinates import (
    coordinate_missing_mask,
    positive_int_position_series,
)
from .._core.columns import (
    INTERNAL_SUMSTATS_ARTIFACT_SPEC_MAP,
    normalize_genome_build,
    resolve_optional_column,
    resolve_required_column,
)
from ..config import GlobalConfig, MungeConfig, normalize_trait_name, get_global_config
from ..errors import (
    LDSCConfigError,
    LDSCDependencyError,
    LDSCInputError,
    LDSCUsageError,
    LDSCUserError,
)
from .._genome.hm3 import packaged_hm3_curated_map_path
from .._runtime.paths import (
    ensure_output_directory,
    preflight_output_artifact_family,
    remove_output_artifacts,
    resolve_scalar_path,
)
from .._runtime.workflow_logging import log_inputs, log_outputs, log_summary, materializing_overwrite_guard, workflow_logging
from .._core.snp_identity import (
    IDENTITY_DROP_COLUMNS,
    clean_identity_artifact_table,
    coerce_identity_drop_frame,
    effective_merge_key_series,
    identity_artifact_metadata,
    identity_mode_family,
    is_allele_aware_mode,
    validate_identity_artifact_metadata,
    normalize_snp_identifier_mode,
)
from .._genome.liftover import SumstatsLiftoverRequest, default_liftover_metadata
from . import munge_kernel as kernel_munge
from . import munge_input


LOGGER = logging.getLogger("LDSC.sumstats_munger")
null_values = munge_input.null_values
default_cnames = munge_input.default_cnames
read_header = munge_input.read_header
get_cname_map = munge_input.get_cname_map
get_compression = munge_input.get_compression
clean_header = munge_input.clean_header

_SUMSTATS_OUTPUT_FORMATS = {"parquet", "tsv.gz", "both"}
_RAW_SUMSTATS_FORMATS = {"auto", "plain", "daner-old", "daner-new"}
_SUMSTATS_PARQUET_COMPRESSION = "snappy"
_REQUIRED_COLUMN_TROUBLESHOOTING = (
    "docs/user-guide/troubleshooting.md#munge-sumstats-could-not-map-a-required-column"
)
_CURATED_ARTIFACT_TROUBLESHOOTING = (
    "docs/user-guide/troubleshooting.md#munge-sumstats-munged-summary-statistics-are-malformed-or-outdated"
)


@dataclass(frozen=True)
class RawSumstatsInference:
    """Header-level inference report for one raw summary-statistics file.

    The report is produced by ``--infer-only`` and by the normal munging
    workflow before the kernel is called. It records only safe, format-aware
    decisions: detected raw format, column hints that can be applied without
    changing statistical meaning, INFO-list handling, missing required fields,
    and command-line repair suggestions. Suggested commands should prefer
    explicit, resolved, reproducible values. ``auto`` should remain the
    user-friendly input default, not the diagnostic output.
    """

    detected_format: str
    column_hints: dict[str, str] = field(default_factory=dict)
    signed_sumstats_spec: str | None = None
    ignore_columns: tuple[str, ...] = ()
    info_list_columns: tuple[str, ...] = ()
    missing_fields: tuple[str, ...] = ()
    suggested_args: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()
    source_genome_build: str | None = None
    output_genome_build: str | None = None
    liftover_required: bool | None = None
    liftover_method: str | None = None

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
        Assigned trait label for regression summaries, or ``None`` when no
        label was supplied or stored. Regression outputs then label the trait
        by the ``source_path`` filename.
    provenance : dict, optional
        Lightweight run metadata retained for debugging and output summaries.
        Tables from :meth:`SumstatsMunger.run` record ``qc_thresholds``, the
        effective ``info_min`` and ``maf_min`` (``None`` for a filter not
        applied), which :meth:`SumstatsMunger.write_output` writes to the Parquet
        footer. Default is an empty dict.
    config_snapshot : GlobalConfig or None, optional
        Shared configuration captured when the table was produced in-process or
        recovered from the ``sumstats.parquet`` footer. It is ``None`` for a
        legacy ``.sumstats`` or ``.sumstats.gz`` input, whose rsIDs are
        projected onto the regression SNPs in the LD-score directory downstream. Footerless
        Parquet is rejected. Coordinate-family munged artifacts store the final output
        genome build here; rsID-family artifacts store ``genome_build=None``
        because their merge identity is independent of coordinate build.
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
            raise LDSCInputError(
                f"munge-sumstats could not map required column(s) {sorted(missing)} "
                f"from the input header for '{self.source_path or self.trait_name or 'sumstats'}'. "
                "Most likely the file uses unrecognized column names. Pass explicit hints "
                f"(e.g. --snp, --a1) or rename the columns. Other causes & fixes: "
                f"{_REQUIRED_COLUMN_TROUBLESHOOTING}"
            )
        if self.has_alleles and not {"A1", "A2"}.issubset(self.data.columns):
            raise LDSCInputError(
                f"munge-sumstats expected A1/A2 allele columns in "
                f"'{self.source_path or self.trait_name or 'sumstats'}' because this table is allele-aware. "
                "Most likely the artifact was produced without allele columns. Re-run munge-sumstats "
                "with A1/A2 column hints or use a base --snp-identifier mode (rsid or chr_pos)."
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

    ``output_paths`` records munged summary-statistics data artifacts and the dropped-SNP
    audit sidecar. It intentionally excludes ``diagnostics/sumstats.log`` so Python result
    contracts stay aligned with other workflow modules.

    ``n_input_rows`` counts parsed records, excluding headers, leading metadata
    and blank lines. ``drop_counts`` records exclusive removals in stage order:
    ``NA``, ``coordinates``, ``INFO``, ``FRQ``, ``P``, ``sumstats_snps``, ``N``,
    ``NSTUDY``, ``liftover``, and ``identity``. Their sum equals input minus
    retained rows. Coordinate reason counts in provenance may overlap, whereas
    these stage totals count each removed row once. ``used_n_rule`` identifies
    the strategy actually used: ``input_columns``, ``fixed_N``, or
    ``fixed_case_control_N``. ``qc_thresholds`` records the effective ``INFO``
    and ``FRQ`` (minor-allele frequency) thresholds; ``None`` means the filter
    did not run because no such column was selected.
    """
    n_input_rows: int
    n_retained_rows: int
    drop_counts: dict[str, int]
    inferred_columns: dict[str, str]
    used_n_rule: str
    output_paths: dict[str, str]
    qc_thresholds: dict[str, float | None]


def load_sumstats(path: str | PathLike[str], trait_name: str | None = None) -> SumstatsTable:
    """Load one munged summary-statistics artifact into a ``SumstatsTable``.

    Parameters
    ----------
    path : str or os.PathLike[str]
        Path token for the munged summary-statistics artifact. This may be a
        literal path or an exact-one glob pattern. Resolution happens at the
        workflow layer before suffix inference. ``.parquet`` files are read
        with :func:`pandas.read_parquet`; ``.sumstats.gz``, ``sumstats.gz``, and ``.sumstats``
        files are read as whitespace-delimited text. Other suffixes raise a
        clear ``ValueError``.
    trait_name : str or None, optional
        Optional trait label propagated into downstream regression summaries.
        When supplied, this value replaces any Parquet footer ``trait_name``
        for the returned table; the artifact on disk is not modified. When
        omitted, the loader uses the parquet footer ``trait_name`` when
        present; otherwise the returned ``trait_name`` is ``None`` and
        regression outputs label the trait by its filename. Surrounding
        whitespace is trimmed; a blank label is rejected. Default is ``None``.

    Returns
    -------
    SumstatsTable
        Validated in-memory table with canonical LDSC3 columns such as ``SNP``,
        ``CHR``, ``POS``, ``N``, and ``Z`` when present in the artifact. When the
        Parquet footer carries identity metadata, the returned table
        also recovers its munge-time ``GlobalConfig`` snapshot. Legacy text
        inputs have ``config_snapshot=None`` and are marked for panel projection.

    Raises
    ------
    ValueError
        If ``path`` does not resolve to exactly one file or if the munged
        artifact is missing required LDSC3 columns.

    Notes
    -----
    Format inference is suffix-based after exact-one path resolution:
    ``.parquet`` uses :func:`pandas.read_parquet`, ``.sumstats.gz`` (or the
    no-label filename ``sumstats.gz``) uses a
    gzip-compressed whitespace reader, and ``.sumstats`` uses a plain
    whitespace reader. Other suffixes are refused so callers do not
    accidentally parse CSV or raw GWAS input as munged summary statistics.

    Identity metadata is read from the parquet footer (discrete ``ldsc:*``
    keys). No ``metadata.json`` sidecar is consulted. Legacy ``.sumstats`` and
    ``.sumstats.gz`` inputs load with ``config_snapshot=None`` and explicit
    legacy LDSC provenance. A Parquet artifact without the required footer is rejected.
    """
    resolved = resolve_scalar_path(path, label="munged summary statistics")
    df = _read_curated_sumstats_artifact(resolved)
    try:
        resolved_columns = _resolve_curated_sumstats_columns(list(df.columns), context=resolved)
    except ValueError as exc:
        raise LDSCInputError(
            f"Cannot load munged summary statistics at '{resolved}': could not map the required SNP, N, and Z "
            "columns from its header. Most likely this is not a current munged summary-statistics file. "
            "Re-run munge-sumstats from the raw GWAS file. Other causes & fixes: "
            f"{_REQUIRED_COLUMN_TROUBLESHOOTING}"
        ) from exc
    df = df.loc[:, list(resolved_columns.values())].rename(
        columns={actual: canonical for canonical, actual in resolved_columns.items()}
    )
    metadata = _read_sumstats_parquet_footer(resolved)
    if metadata is None:
        if str(resolved).lower().endswith(".parquet"):
            raise LDSCInputError(
                _curated_sumstats_artifact_message(
                    resolved,
                    "the parquet file is missing required LDSC3 identity footer metadata",
                    "Footerless parquet is not a legacy compatibility artifact. Re-run `ldsc munge-sumstats` "
                    "from the raw GWAS input to create a self-describing sumstats.parquet file.",
                )
            )
        LOGGER.info(
            f"Recognized legacy munged summary statistics '{resolved}'. SNP values will be used as rsID "
            "lookup keys and projected onto the regression SNPs in the LD-score directory during regression."
        )
        table = SumstatsTable(
            data=df.reset_index(drop=True),
            has_alleles={"A1", "A2"}.issubset(df.columns),
            source_path=resolved,
            trait_name=_resolve_sumstats_trait_name(trait_name, None, resolved),
            provenance={"source_format": "legacy_sumstats", "legacy_ldsc": True},
            config_snapshot=None,
        )
        table.validate()
        return table
    config_snapshot = _global_config_from_sumstats_metadata(metadata, artifact_path=resolved)
    mode = normalize_snp_identifier_mode(config_snapshot.snp_identifier)
    if is_allele_aware_mode(mode) and not {"A1", "A2"}.issubset(df.columns):
        raise LDSCInputError(
            _curated_sumstats_artifact_message(
                resolved,
                f"snp_identifier='{mode}' requires A1/A2 columns, but the artifact does not contain both.",
                "Most likely the artifact was produced by an older munging run or one in a base mode.",
            )
        )
    cleanup = clean_identity_artifact_table(
        df,
        mode,
        context="loaded munged summary-statistics artifact",
        stage="load_sumstats_validation",
        logger=None,
    )
    if not cleanup.dropped.empty:
        reasons = ", ".join(
            f"{reason}={int(count)}"
            for reason, count in cleanup.dropped["reason"].value_counts(sort=False).items()
        )
        raise LDSCInputError(
            _curated_sumstats_artifact_message(
                resolved,
                f"duplicate or invalid SNP identity rows were found ({reasons}).",
                "Most likely the artifact was hand-edited or produced before identity cleanup was enforced.",
            )
        )
    effective_keys = effective_merge_key_series(
        cleanup.cleaned,
        mode,
        context="loaded munged summary-statistics artifact",
    )
    if bool(effective_keys.isna().any()):
        raise LDSCInputError(
            _curated_sumstats_artifact_message(
                resolved,
                "missing or invalid SNP identity rows were found after identity cleanup.",
                "Most likely the artifact contains incomplete CHR/POS or allele identity columns.",
            )
        )
    table = SumstatsTable(
        data=df.reset_index(drop=True),
        has_alleles={"A1", "A2"}.issubset(df.columns),
        source_path=resolved,
        trait_name=_resolve_sumstats_trait_name(trait_name, metadata, resolved),
        provenance={"metadata": metadata},
        config_snapshot=config_snapshot,
    )
    table.validate()
    return table


class SumstatsMunger:
    """Run legacy-compatible munging through a typed workflow interface."""

    def __init__(self) -> None:
        """Initialize the workflow wrapper and clear any cached run summary."""
        self._last_summary: MungeRunSummary | None = None

    @materializing_overwrite_guard(
        lambda self, raw_sumstats_config, munge_config=None, global_config=None: (
            (
                (raw_sumstats_config if munge_config is None or isinstance(munge_config, GlobalConfig) else munge_config).output_dir,
                (raw_sumstats_config if munge_config is None or isinstance(munge_config, GlobalConfig) else munge_config).overwrite,
                "RUN_FAILED.txt",
            )
            if (raw_sumstats_config if munge_config is None or isinstance(munge_config, GlobalConfig) else munge_config).output_dir
            else None
        ),
        command="SumstatsMunger.run(...)",
    )
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
            QC settings. Common plain-text inputs, including VCF-style headers,
            old DANER, and new DANER are handled before the kernel runs;
            explicit column hints still take priority.
        munge_config : MungeConfig or None, optional
            Munging thresholds, output directory, and munged output format. The
            workflow writes ``<trait>.parquet`` (self-describing via its footer)
            and/or ``<trait>.sumstats.gz``, using the filesystem-safe trait label
            from ``raw_sumstats_config``. Without a label, filenames are
            ``sumstats.parquet`` and ``sumstats.gz``. Diagnostics live under
            ``diagnostics/``. Any existing owned artifact is refused before the kernel
            runs unless ``munge_config.overwrite`` is true; successful
            overwrites remove stale sibling formats that the current run did
            not produce. Packaged HM3 is the default keep-list; a custom
            ``sumstats_snps_file`` replaces it, while ``no_snp_restriction=True``
            disables keep-list filtering. These overrides are mutually exclusive.
            The effective headered keep-list is loaded once before raw chunk
            parsing, and applied to each parsed chunk after munging QC and
            coordinate normalization. Allele-free keep-lists match by base key,
            and allele-bearing keep-lists in allele-aware modes match by the
            effective allele-aware key. Duplicate keep-list keys collapse to one
            retained key, and non-identity columns such as ``CM`` or ``MAF`` are
            ignored; keep-lists do not reorder output rows.
            In coordinate-family modes, keep-list filtering uses source-build
            coordinates before output liftover. The output build must be chosen
            explicitly. Different builds use automatic quick liftover from
            bundled HM3 metadata only under packaged restriction; an explicit
            chain replaces that method and is required for custom/unrestricted
            cross-build runs. Matching builds ignore any chain; unresolved
            source builds stop and request an explicit source build.
            ``sumstats_format="auto"`` is the default; use
            ``sumstats_format="plain"``, ``"daner-old"``, or ``"daner-new"``
            only when overriding auto-detection.
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
            same identity settings into the Parquet footer so downstream
            regression can detect incompatible LD-score results after reload.
            The selected SNP restriction and ordinary QC are reflected in both the
            returned table and the written munged artifact(s).
            Output paths for the corresponding disk artifacts are available
            through :meth:`build_run_summary`; the workflow log is written to
            ``diagnostics/sumstats.log`` but is not included in that result
            mapping. Method selection and mapping/drop counts are always recorded
            in the log; this Python method does not print a console summary.
        """
        if munge_config is None:
            munge_config = raw_sumstats_config
        elif isinstance(munge_config, GlobalConfig) and global_config is None:
            global_config = munge_config
            munge_config = raw_sumstats_config
        if raw_sumstats_config.raw_sumstats_file is None:
            raise LDSCUserError(
                "No input summary statistics given to munge-sumstats. Most likely --raw-sumstats-file was omitted. "
                "Pass --raw-sumstats-file <file>."
            )
        if munge_config.output_dir is None:
            raise LDSCUserError(
                "No output directory given to munge-sumstats. Most likely --output-dir or output_dir was omitted. "
                "Pass --output-dir <directory> on the CLI or set MungeConfig.output_dir."
            )

        config_snapshot = _source_global_config_for_munge(munge_config, global_config or get_global_config())
        _validate_munge_build_contract(munge_config, config_snapshot)
        liftover_request = _liftover_request_from_config(munge_config, config_snapshot.genome_build)
        _validate_liftover_request_before_io(config_snapshot, liftover_request)
        output_dir = ensure_output_directory(munge_config.output_dir, label="output directory")
        diagnostics_dir = output_dir / "diagnostics"
        output_files = _sumstats_output_files(output_dir, munge_config.output_format, raw_sumstats_config.trait_name)
        log_path = str(diagnostics_dir / "sumstats.log")
        dropped_snps_path = diagnostics_dir / "dropped_snps" / "dropped.tsv.gz"
        sumstats_snps_path = _sumstats_snps_file_from_config(munge_config)
        sumstats_snps_label = "none" if sumstats_snps_path is None else str(sumstats_snps_path)
        produced_paths = [*output_files.values(), log_path, dropped_snps_path]
        owned_paths = [
            *_sumstats_output_files(output_dir, "both", raw_sumstats_config.trait_name).values(),
            log_path,
            dropped_snps_path,
        ]
        # See _runtime.paths.preflight_output_artifact_family for the
        # owned_paths/produced_paths split. The dropped-SNP sidecar is in both
        # lists because it is always written, including clean header-only runs.
        stale_paths = preflight_output_artifact_family(
            produced_paths,
            owned_paths,
            overwrite=munge_config.overwrite,
            label="munged output artifact",
        )
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        with workflow_logging("munge-sumstats", log_path, log_level=config_snapshot.log_level):
            from .._runtime.preflight import inspect_declared_inputs
            inspect_declared_inputs(files=[('raw summary statistics', raw_sumstats_config.raw_sumstats_file),
                ('SNP restriction', sumstats_snps_path), ('liftover chain', liftover_request.liftover_chain_file)],
                issues_path=diagnostics_dir/'input_issues.tsv')
            source_path = resolve_scalar_path(raw_sumstats_config.raw_sumstats_file, label="raw sumstats")
            raw_sumstats_config, munge_config, inference = _apply_raw_sumstats_inference(
                source_path, raw_sumstats_config, munge_config
            )
            restriction_path = _resolve_sumstats_snps_path(munge_config)
            log_inputs(
                raw_sumstats_file=source_path,
                output_dir=str(output_dir),
                output_format=munge_config.output_format,
                sumstats_snps_file=sumstats_snps_label,
                snp_restriction=("hm3" if _uses_packaged_hm3(munge_config) else "none" if munge_config.no_snp_restriction else "custom"),
                hm3_map_file=packaged_hm3_curated_map_path() if _uses_packaged_hm3(munge_config) else "none",
                output_genome_build=liftover_request.target_build or "none",
                liftover_method=("pending source inference" if config_snapshot.genome_build == "auto" else liftover_request.method or "none"),
            )
            LOGGER.info(
                f"Munging summary statistics from '{source_path}' into '{output_dir}' "
                f"with snp_identifier='{config_snapshot.snp_identifier}', "
                f"genome_build='{config_snapshot.genome_build}', "
                f"sumstats_snps_file='{sumstats_snps_label}', "
                f"packaged_hm3='{_uses_packaged_hm3(munge_config)}', "
                f"output_genome_build='{liftover_request.target_build}'."
            )
            request = munge_input.prepare_munge_input(
                source_path, raw_sumstats_config, munge_config, config_snapshot,
                liftover_request, restriction_path,
            )
            resolved_liftover = _liftover_request_from_config(munge_config, request.genome_build)
            if resolved_liftover.liftover_chain_file is not None:
                resolved_liftover = replace(
                    resolved_liftover,
                    liftover_chain_file=resolve_scalar_path(resolved_liftover.liftover_chain_file, label="liftover chain file"),
                )
            _validate_liftover_request_before_io(
                replace(config_snapshot, genome_build=request.genome_build), resolved_liftover,
            )
            request = replace(request, liftover_request=resolved_liftover)
            method_label = {
                "hm3_curated": "automatic HM3 quick liftover (package-bundled reference HM3 metadata)",
                "chain_file": "explicit chain file (HM3 quick liftover disabled)",
                None: "none; no liftover required",
            }[resolved_liftover.method]
            LOGGER.info("Selected liftover method: %s; source=%s; output=%s",
                        method_label, request.genome_build, resolved_liftover.target_build)
            result = kernel_munge.munge_sumstats(request)
            data = result.data
            coordinate_metadata = result.coordinate_metadata
            drop_frame = _coerce_sumstats_dropped_snps_frame(
                result.liftover_drop_frame
            )
            identity_drop_frame = _coerce_sumstats_dropped_snps_frame(
                result.identity_drop_frame,
                default_stage=None,
            )
            drop_frame = _coerce_sumstats_dropped_snps_frame(
                pd.concat([drop_frame, identity_drop_frame], ignore_index=True),
                default_stage=None,
            )
            table_config_snapshot = _effective_sumstats_config(config_snapshot, coordinate_metadata)
            qc_thresholds = {"info_min": request.qc.info_min, "maf_min": request.qc.maf_min}
            primary_sumstats_file, parquet_row_groups = _write_sumstats_outputs(
                data,
                output_files=output_files,
                output_format=munge_config.output_format,
                footer_metadata=_sumstats_footer_metadata(
                    table_config_snapshot, raw_sumstats_config.trait_name, qc_thresholds,
                ),
            )
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
                    "coordinate_provenance": coordinate_metadata,
                    "metadata": coordinate_metadata,
                    "qc_thresholds": qc_thresholds,
                },
                config_snapshot=table_config_snapshot,
            )
            table.validate()
            run_output_paths = {
                **({"sumstats_parquet": output_files["parquet"]} if "parquet" in output_files else {}),
                **({"sumstats_gz": output_files["tsv.gz"]} if "tsv.gz" in output_files else {}),
                "dropped_snps_tsv_gz": str(dropped_snps_path),
            }
            self._last_summary = MungeRunSummary(
                n_input_rows=result.n_input_rows,
                n_retained_rows=len(table.data),
                drop_counts=result.drop_counts,
                inferred_columns={**dict(raw_sumstats_config.column_hints), "detected_format": inference.detected_format},
                used_n_rule=result.used_n_rule,
                output_paths=run_output_paths,
                qc_thresholds={"INFO": request.qc.info_min, "FRQ": request.qc.maf_min},
            )
            log_summary(_render_munge_summary(self._last_summary, coordinate_metadata, munge_config))
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
        """Write an in-memory summary-statistics table using its trait label for filenames.

        Parameters
        ----------
        sumstats : SumstatsTable
            Validated LDSC-ready table to persist. Columns are curated to the
            package output order ``SNP, CHR, POS, A1, A2, Z, N, FRQ`` where
            present; missing ``CHR`` or ``POS`` columns are materialized as
            missing values.
        output_dir : str or os.PathLike[str]
            Destination directory for munged summary-statistics artifacts.
        output_format : {"parquet", "tsv.gz", "both"}, optional
            Disk format to write. ``"parquet"`` is the default and returns
            ``<trait>.parquet``. ``"tsv.gz"`` writes ``<trait>.sumstats.gz``.
            Without a trait label, names are ``sumstats.parquet`` and
            ``sumstats.gz``. ``"both"`` writes both and returns the Parquet path.
        overwrite : bool, optional
            If ``True``, replace artifacts for the resolved filename stem and remove
            stale sibling formats after a successful write. If ``False``, any
            existing artifact in that filename family is refused before writing
            starts. Default is ``False``.

        Returns
        -------
        str
            Primary summary-statistics artifact path. The primary path is Parquet for
            ``"parquet"`` and ``"both"``, and gzip TSV for ``"tsv.gz"``.

        Notes
        -----
        This helper follows the public workflow naming policy but does not
        change the trait label stored in metadata. Filename sanitization preserves
        letter case, Unicode word characters, dots, and hyphens; replaces other
        character runs with underscores; strips edge dots, underscores, and
        hyphens; and uses ``trait`` if an explicit label leaves an empty stem.
        Different labels can resolve to the same filename and then follow the
        ordinary overwrite policy. Artifacts for other stems are preserved.
        This helper does not create ``diagnostics/sumstats.log`` because no raw
        munging kernel is run. A config snapshot is required so the Parquet footer can record
        identity provenance for downstream recovery. The footer records the QC
        thresholds only when ``sumstats.provenance`` holds ``qc_thresholds``, as
        tables returned by :meth:`run` do; a table built by hand or loaded from
        disk does not, so the ``ldsc:qc_*`` keys are omitted.
        """
        if sumstats.config_snapshot is None:
            raise LDSCInputError(
                "Cannot write current munge-sumstats artifacts from this SumstatsTable because "
                "config_snapshot is missing. Most likely the table was constructed manually or loaded "
                "from a legacy artifact. Re-run munge-sumstats from the raw GWAS file or attach a "
                "GlobalConfig snapshot before writing."
            )
        output_format = _normalize_output_format(output_format)
        output_root = ensure_output_directory(output_dir, label="output directory")
        output_files = _sumstats_output_files(output_root, output_format, sumstats.trait_name)
        produced_paths = list(output_files.values())
        stale_paths = preflight_output_artifact_family(
            produced_paths,
            list(_sumstats_output_files(output_root, "both", sumstats.trait_name).values()),
            overwrite=overwrite,
            label="munged output artifact",
        )
        primary_sumstats_file, _parquet_row_groups = _write_sumstats_outputs(
            sumstats.data,
            output_files=output_files,
            output_format=output_format,
            footer_metadata=_sumstats_footer_metadata(
                sumstats.config_snapshot, sumstats.trait_name, sumstats.provenance.get("qc_thresholds"),
            ),
        )
        remove_output_artifacts(stale_paths)
        return primary_sumstats_file

    def build_run_summary(self) -> MungeRunSummary:
        """Return the summary captured from the most recent call to :meth:`run`."""
        if self._last_summary is None:
            raise LDSCUsageError(
                "Cannot build a munge-sumstats run summary before a run has completed. "
                "Most likely build_run_summary() was called on a fresh SumstatsMunger. "
                "Call SumstatsMunger.run(...) first, then request the summary."
            )
        return self._last_summary


def run_munge_sumstats_from_args(args: argparse.Namespace) -> SumstatsTable | RawSumstatsInference:
    """Run summary-statistics munging from parsed CLI arguments.

    The CLI path normalizes argparse values into the same ``MungeConfig``
    objects used by the Python API, then delegates to :class:`SumstatsMunger`
    so path resolution, output preflight, metadata, and result construction
    stay in one workflow path. Successful materializing runs print the resolved
    liftover method and mapping/drop summary to stdout at every log level;
    direct Python calls to ``SumstatsMunger.run`` remain console-quiet.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments from :func:`build_parser`. The namespace must include
        ``raw_sumstats_file`` and the public munging options. ``output_dir`` is required for every CLI run,
        including ``infer_only`` runs (which still write no artifacts).

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
        If artifacts for the resolved names already exist and ``args.overwrite`` is
        false.
    """
    if not str(getattr(args, "raw_sumstats_file", "") or "").strip():
        raise LDSCUsageError(
            "munge-sumstats requires a nonempty --raw-sumstats-file, including with --infer-only. "
            "Pass --raw-sumstats-file <file>; if using a shell variable, make sure it contains the input path."
        )
    if not str(getattr(args, "output_dir", "") or "").strip():
        raise LDSCUsageError("munge-sumstats requires a nonempty --output-dir. Pass --output-dir <directory>, including with --infer-only; dry runs write nothing.")
    raw_config, munge_config = _munge_configs_from_args(args)
    if getattr(args, "infer_only", False):
        global_config = _resolve_main_global_config(args)
        source_path = resolve_scalar_path(raw_config.raw_sumstats_file, label="raw sumstats")
        inference = infer_raw_sumstats(source_path, raw_config, munge_config, global_config)
        inference = _apply_build_inference_report(source_path, raw_config, munge_config, global_config, inference)
        inference = _apply_missing_column_report(source_path, raw_config, munge_config, inference)
        inference = _preserve_inference_options(inference, args)
        print(_render_inference_report(inference, source_path, args.output_dir))
        return inference
    munger = SumstatsMunger()
    table = munger.run(raw_config, munge_config, _resolve_main_global_config(args))
    print(_render_munge_summary(munger.build_run_summary(), table.provenance["coordinate_provenance"], munge_config))
    return table


def _render_munge_summary(summary: MungeRunSummary, metadata: dict[str, Any], config: MungeConfig) -> str:
    """Render resolved liftover outcomes and exclusive whole-run drop counts."""
    restriction = "packaged HM3 (default)" if _uses_packaged_hm3(config) else "none" if config.no_snp_restriction else "custom list"
    liftover = metadata["liftover"]
    lines = ["Munge-sumstats summary:", f"  SNP restriction: {restriction}"]
    if liftover["applied"]:
        method = (
            "HM3 quick liftover (automatic; package-bundled HM3 metadata)"
            if liftover["method"] == "hm3_curated"
            else "chain file (explicit; HM3 quick liftover disabled)"
        )
        lines.extend([
            f"  Liftover: {method}; {liftover['source_build']} -> {liftover['target_build']}",
            f"  Mapping: {liftover['n_input']} input; {liftover['n_lifted']} mapped and retained; {liftover['n_dropped']} dropped",
            "  Liftover drop reasons: " + "; ".join(
                f"{label}={liftover[key]}" for label, key in (
                    ("missing coordinates", "n_missing_chr_pos_dropped"),
                    ("duplicate source", "n_duplicate_source_dropped"),
                    ("unmapped", "n_unmapped"),
                    ("cross-chromosome", "n_cross_chrom"),
                    ("duplicate target", "n_duplicate_target_dropped"),
                )
            ),
        ])
    else:
        if identity_mode_family(metadata["snp_identifier"]) == "rsid":
            reason = "rsID identity; liftover not applicable"
        else:
            reason = "source and output builds match"
            if config.liftover_chain_file is not None:
                reason += "; supplied chain ignored"
        lines.extend([f"  Liftover: none ({reason})", "  Mapping: not performed; 0 liftover drops"])
    lines.extend([
        f"  Rows: {summary.n_input_rows} input; {summary.n_retained_rows} retained; {summary.n_input_rows - summary.n_retained_rows} dropped",
        "  Drop counts by stage: " + "; ".join(f"{stage}={count}" for stage, count in summary.drop_counts.items()),
    ])
    skipped = [label for field, label in (("INFO", "INFO (no INFO column)"), ("FRQ", "MAF (no frequency column)"))
               if summary.qc_thresholds[field] is None]
    if skipped:
        lines.append("  QC filters not applied: " + "; ".join(skipped))
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> SumstatsTable | RawSumstatsInference:
    """CLI entry point: parse arguments and delegate to workflow orchestration."""
    return run_munge_sumstats_from_args(build_parser().parse_args(argv))


def _resolve_main_global_config(args: argparse.Namespace) -> GlobalConfig:
    mode = normalize_snp_identifier_mode(getattr(args, "snp_identifier", "chr_pos_allele_aware"))
    source_build = normalize_genome_build(getattr(args, "source_genome_build", "auto"))
    output_build = normalize_genome_build(getattr(args, "output_genome_build", None))
    if identity_mode_family(mode) == "rsid":
        if source_build not in {None, "auto"}:
            raise LDSCUsageError(
                "munge-sumstats cannot use `--source-genome-build` with an rsID-family SNP identifier. "
                "Most likely a coordinate-build option was copied from a chr_pos command into an rsID-based run. "
                "Drop the build flag or switch to a chr_pos identifier."
            )
        if output_build is not None:
            raise LDSCUsageError(
                "munge-sumstats cannot use `--output-genome-build` with an rsID-family SNP identifier. "
                "Most likely a coordinate-liftover option was copied from a chr_pos command into an rsID-based run. "
                "Drop the output build flag or switch to a chr_pos identifier."
            )
        if getattr(args, "liftover_chain_file", None) is not None:
            raise LDSCUsageError(
                "munge-sumstats cannot run summary-statistics liftover with an rsID-family SNP identifier. "
                "Most likely a liftover option was supplied while SNP identity is based on rsIDs rather than coordinates. "
                "Drop the liftover option or switch to a chr_pos identifier."
            )
        return GlobalConfig(snp_identifier=mode, log_level=getattr(args, "log_level", "INFO"))
    if source_build is None:
        raise LDSCUsageError(
            "munge-sumstats needs a source genome build for chr_pos-family --snp-identifier modes. "
            "Most likely --source-genome-build was omitted or set to an empty value. "
            "Pass --source-genome-build auto, hg19, or hg38."
        )
    if output_build is None:
        raise LDSCUsageError(
            "munge-sumstats needs an output genome build for chr_pos-family --snp-identifier modes. "
            "Most likely --output-genome-build was omitted. Pass --output-genome-build hg19 or hg38."
        )
    return GlobalConfig(snp_identifier=mode, genome_build=source_build, log_level=getattr(args, "log_level", "INFO"))


def _source_global_config_for_munge(config: MungeConfig, base: GlobalConfig) -> GlobalConfig:
    """Return the source-build config used while reading raw summary statistics."""
    mode = normalize_snp_identifier_mode(base.snp_identifier)
    if identity_mode_family(mode) == "rsid":
        return GlobalConfig(
            snp_identifier=mode,
            log_level=base.log_level,
        )
    return GlobalConfig(
        snp_identifier=mode,
        genome_build=config.source_genome_build,
        log_level=base.log_level,
    )


def _validate_munge_build_contract(config: MungeConfig, source_config: GlobalConfig) -> None:
    """Validate source/output build options against the active SNP identifier mode."""
    if identity_mode_family(source_config.snp_identifier) == "rsid":
        if config.source_genome_build != "auto":
            raise LDSCUsageError(
                "MungeConfig.source_genome_build cannot be a concrete build in rsID-family snp_identifier modes. "
                "Most likely coordinate-build configuration was reused for an rsID-based run. "
                "Set source_genome_build='auto' or switch to a chr_pos identifier."
            )
        if config.output_genome_build is not None:
            raise LDSCUsageError(
                "MungeConfig.output_genome_build cannot be set in rsID-family snp_identifier modes. "
                "Most likely coordinate-liftover configuration was reused for an rsID-based run. "
                "Set output_genome_build=None or switch to a chr_pos identifier."
            )
        if config.liftover_chain_file is not None:
            raise LDSCUsageError(
                "MungeConfig cannot request summary-statistics liftover in rsID-family snp_identifier modes. "
                "Most likely coordinate-liftover configuration was reused while SNP identity is based on rsIDs. "
                "Remove liftover options or switch to a chr_pos identifier."
            )
        return
    if config.output_genome_build is None:
        raise LDSCUsageError(
            "MungeConfig.output_genome_build is required for chr_pos-family snp_identifier modes. "
            "Most likely output_genome_build was omitted. Set output_genome_build='hg19' or 'hg38'."
        )


def _uses_packaged_hm3(config: MungeConfig) -> bool:
    """Return whether the default packaged restriction is active."""
    return not config.no_snp_restriction and config.sumstats_snps_file is None


def _sumstats_snps_file_from_config(config: MungeConfig) -> str | None:
    """Return the explicit or packaged SNP keep-list path for one munger run."""
    if _uses_packaged_hm3(config):
        return packaged_hm3_curated_map_path()
    return config.sumstats_snps_file


def _resolve_sumstats_snps_path(config: MungeConfig) -> str | None:
    """Resolve the effective summary-statistics SNP restriction path."""
    path = _sumstats_snps_file_from_config(config)
    if path is None:
        return None
    label = "packaged HM3 SNP map" if _uses_packaged_hm3(config) else "summary-statistics SNPs file"
    return resolve_scalar_path(path, label=label)


def _liftover_request_from_config(config: MungeConfig, source_build: str | None) -> SumstatsLiftoverRequest:
    """Resolve the mapping method only when concrete coordinate builds differ."""
    target = config.output_genome_build
    if source_build not in {"hg19", "hg38"} or target is None or source_build == target:
        return SumstatsLiftoverRequest(target_build=target)
    quick = _uses_packaged_hm3(config) and config.liftover_chain_file is None
    return SumstatsLiftoverRequest(
        target_build=target,
        liftover_chain_file=config.liftover_chain_file,
        use_hm3_quick_liftover=quick,
        hm3_map_file=packaged_hm3_curated_map_path() if quick else None,
    )


def _validate_liftover_request_before_io(config: GlobalConfig, request: SumstatsLiftoverRequest) -> None:
    """Reject liftover requests that can be proven invalid before input IO."""
    if request.requested and identity_mode_family(config.snp_identifier) != "chr_pos":
        raise LDSCUsageError(
            "Summary-statistics liftover applies only to chr_pos-family snp_identifier modes. "
            "Most likely liftover was requested with an rsID-based identifier. "
            "Drop the liftover option or switch to a chr_pos identifier."
        )
    source_build = normalize_genome_build(config.genome_build)
    if source_build not in {"hg19", "hg38"} or request.target_build is None:
        return
    if request.target_build != source_build and request.method is None:
        raise LDSCUsageError(
            "munge-sumstats cannot convert source_genome_build to output_genome_build because no liftover "
            "method was specified. Most likely the requested builds differ. Add --liftover-chain-file "
            "<chain.over.chain.gz> for custom-list or unrestricted SNPs; packaged HM3 uses automatic quick liftover."
        )


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
        N_min=getattr(args, "N_min", None),
        nstudy_min=getattr(args, "nstudy_min", None),
        chunk_size=args.chunksize,
        output_format=args.output_format,
        sumstats_snps_file=getattr(args, "sumstats_snps_file", None),
        no_snp_restriction=getattr(args, "no_snp_restriction", False),
        source_genome_build=getattr(args, "source_genome_build", "auto"),
        output_genome_build=getattr(args, "output_genome_build", None),
        liftover_chain_file=getattr(args, "liftover_chain_file", None),
        signed_sumstats_spec=getattr(args, "signed_sumstats", None),
        ignore_columns=_ignore_columns_from_args(args),
        info_list_columns=_info_list_columns_from_args(args),
        sumstats_format=getattr(args, "sumstats_format", "auto"),
        a1_inc=args.a1_inc,
        overwrite=getattr(args, "overwrite", False),
    )
    return raw_config, munge_config


def _column_hints_from_args(args: argparse.Namespace) -> dict[str, str]:
    """Collect explicit raw-column hints from parsed CLI arguments."""
    hints: dict[str, str] = {}
    for key in _COLUMN_HINT_ARG_KEYS:
        value = getattr(args, key, None)
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


def build_parser() -> argparse.ArgumentParser:
    """Build the public summary-statistics munging parser.

    Coordinate-family runs infer the raw source build by default and require an
    explicit output build for downstream-compatible artifacts.
    """
    public = argparse.ArgumentParser(allow_abbrev=False)
    public.prog = 'ldsc munge-sumstats'
    public.formatter_class = CLIHelpFormatter
    public.description = (
        'Prepare GWAS summary statistics for LDSC, restricting to packaged HapMap3 SNPs by default and preserving '
        'input frequency as FRQ when available. Run with --infer-only first to preview the detected columns and '
        'a suggested command without writing files.'
    )
    public.usage = (
        '%(prog)s [--snp-identifier MODE] [--source-genome-build BUILD]\n'
        '    --output-genome-build BUILD --raw-sumstats-file FILE [--trait-name TRAIT_NAME]\n'
        '    --output-dir DIR [options]'
    )
    identity = public.add_argument_group('Step 1: SNP identity and genome build')
    raw = public.add_argument_group('Step 2: Raw summary statistics')
    sample = public.add_argument_group('Step 3: Sample size (only if not detected from the file)', description='Input sample-size columns take precedence over constant fallbacks.')
    output = public.add_argument_group('Step 4: Output')
    filters = public.add_argument_group('Optional: SNP selection and quality filters')
    columns = public.add_argument_group('Optional: column overrides', description='Columns are detected automatically. Supply overrides for unrecognized or ambiguous headers.')
    computing = public.add_argument_group('Computing resources')
    logs = public.add_argument_group('Logging')

    identity.add_argument(
        '--snp-identifier', default='chr_pos_allele_aware', choices=('rsid', 'rsid_allele_aware', 'chr_pos', 'chr_pos_allele_aware'),
        help=(
            'Match SNPs by rsID (rsid) or chromosome and position (chr_pos); allele-aware variants also use '
            'A1/A2. Default: chr_pos_allele_aware; usable alleles are required. Coordinate-family modes also '
            'require --output-genome-build.'
        ),
    )
    identity.add_argument(
        '--source-genome-build', default='auto', choices=('auto', 'hg19', 'hg37', 'GRCh37', 'hg38', 'GRCh38'),
        help=(
            'Genome build of source coordinates. Default: auto, check raw chromosomes and coordinate keep lists; '
            'stop with a grouped report if evidence conflicts or is insufficient. Concrete builds are '
            'checked against declarations and cannot '
            'be used with rsid or rsid_allele_aware identity.'
        ),
    )
    identity.add_argument(
        '--output-genome-build', default=None, choices=('hg19', 'hg37', 'GRCh37', 'hg38', 'GRCh38'),
        help=(
            'Required output build for chr_pos and chr_pos_allele_aware identity. If different from the '
            'source, packaged HM3 uses automatic quick liftover; custom or unrestricted SNPs require '
            '--liftover-chain-file, which overrides quick liftover when supplied. Matching builds need no '
            'liftover. Cannot be used with rsid or rsid_allele_aware identity; no default.'
        ),
    )
    identity.add_argument(
        '--liftover-chain-file', default=None, metavar='FILE',
        help=(
            'Chain file in the source-to-output direction for liftover to --output-genome-build. '
            'Requires chr_pos or chr_pos_allele_aware --snp-identifier and an explicit --output-genome-build. '
            'If omitted, the default packaged HM3 restriction automatically uses quick liftover with '
            'package-bundled reference HM3 metadata when builds differ; supply this flag to disable quick '
            'liftover and use the chain instead. Required when builds differ with --sumstats-snps-file or '
            '--no-snp-restriction; ignored when builds match. '
            + SCALAR_PATH_HELP
        ),
    )

    raw.add_argument(
        '--raw-sumstats-file', required=True, metavar='FILE',
        help=(
            'Required nonempty raw GWAS summary-statistics file path, including with --infer-only. Columns are detected from their headers unless '
            'overridden. '
            + SCALAR_PATH_HELP
        ),
    )
    raw.add_argument(
        '--trait-name', default=None,
        help=(
            (
            'Trait label stored in Parquet metadata and used for data filenames after replacing unsafe '
            'characters with underscores. Default: no explicit label; use sumstats.parquet or sumstats.gz.'
        )
        ),
    )
    raw.add_argument(
        '--input-format', dest='sumstats_format', choices=sorted(_RAW_SUMSTATS_FORMATS), default='auto',
        help=(
            'Input layout: plain reads ordinary tables, daner-old reads counts from frequency headers, and '
            'daner-new reads case/control count columns. Default: auto, detect the layout from headers.'
        ),
    )

    sample.add_argument(
        '--N-col', default=None, type=str, metavar='COLUMN',
        help=(
            'Input per-SNP sample-size column; matched case-insensitively. Overrides inferred case/control '
            'columns; cannot be combined with --N-cas-col or --N-con-col. If omitted, detect sample-size '
            'columns from headers.'
        ),
    )
    sample.add_argument(
        '--N-cas-col', default=None, type=str, metavar='COLUMN',
        help=(
            'Input per-SNP case-count column; matched case-insensitively. Requires --N-con-col and cannot '
            'be combined with --N-col. The pair overrides inferred N; if omitted, detect sample-size '
            'columns from headers.'
        ),
    )
    sample.add_argument(
        '--N-con-col', default=None, type=str, metavar='COLUMN',
        help=(
            'Input per-SNP control-count column; matched case-insensitively. Requires --N-cas-col and '
            'cannot be combined with --N-col. The pair overrides inferred N; if omitted, detect sample-size '
            'columns from headers.'
        ),
    )
    sample.add_argument(
        '--N', default=None, type=float, metavar='VALUE',
        help=(
            'Constant sample-size fallback when per-SNP N or paired case/control columns are absent. Input '
            'columns take precedence; --N takes precedence over --N-cas/--N-con. If omitted, use input '
            'columns or both case/control constants.'
        ),
    )
    sample.add_argument(
        '--N-cas', default=None, type=float, metavar='VALUE',
        help=(
            'Constant case-count fallback. Requires --N-con when used; the constants are added only when '
            'input sample-size columns and --N are absent. Default: no case-count constant.'
        ),
    )
    sample.add_argument(
        '--N-con', default=None, type=float, metavar='VALUE',
        help=(
            'Constant control-count fallback. Requires --N-cas when used; the constants are added only when '
            'input sample-size columns and --N are absent. Default: no control-count constant.'
        ),
    )

    output.add_argument(
        '--output-dir', required=True, metavar='DIR',
        help=(
            'Required destination for munged summary statistics and diagnostics. Still required with '
            '--infer-only, which creates no files or directories.'
        ),
    )
    output.add_argument(
        '--overwrite', action='store_true', default=False,
        help=(
            "Replace this command's existing output files and remove obsolete outputs from an earlier run. "
            'Default: off; stop if output files already exist. Has no effect with --infer-only.'
        ),
    )
    output.add_argument(
        '--output-format', choices=sorted(_SUMSTATS_OUTPUT_FORMATS), default='parquet',
        help=(
            'Write <trait>.parquet, <trait>.sumstats.gz, or both; without --trait-name, use '
            'sumstats.parquet or sumstats.gz. Default: parquet; output choice does not '
            'change SNP filtering or frequency preservation.'
        ),
    )
    output.add_argument(
        '--infer-only', action='store_true', default=False,
        help=(
            'Inspect headers and check selected columns for all-missing values in bounded chunks, then print a suggested command without '
            'writing files. Requires --raw-sumstats-file and --output-dir. Default: off; perform the full '
            'conversion.'
        ),
    )

    restriction = filters.add_mutually_exclusive_group()
    restriction.add_argument(
        '--sumstats-snps-file', default=None, metavar='FILE',
        help=(
            'Replace the default packaged HapMap3 restriction with this headered identity table, matched '
            'in the source genome build. Mutually exclusive with --no-snp-restriction. '
            'Requires --liftover-chain-file when source and output builds differ, even for a custom SNP keep-list. '
            "Same exact-one '*' pattern rules as --raw-sumstats-file; '@' is not expanded."
        ),
    )
    restriction.add_argument(
        '--no-snp-restriction', action='store_true', default=False,
        help=(
            'Disable the default packaged HapMap3 keep-list restriction. Ordinary QC still applies. '
            'Mutually exclusive with --sumstats-snps-file; requires --liftover-chain-file when source and '
            'output builds differ. Default: off; restrict to HapMap3.'
        ),
    )
    filters.add_argument(
        '--N-min', default=None, type=float, metavar='VALUE',
        help=(
            'Minimum per-SNP sample size to retain. If omitted or zero, use the 90th percentile of N '
            'divided by 1.5. Constant sample sizes are assigned after this filter and are not filtered.'
        ),
    )
    filters.add_argument(
        '--nstudy-min', default=None, type=float, metavar='VALUE',
        help=(
            'Minimum number of contributing studies to retain. If omitted or zero, keep the maximum '
            'observed count. Applies only when NSTUDY is present and per-SNP sample sizes are absent.'
        ),
    )
    filters.add_argument(
        '--info-min', default=None, type=float, metavar='VALUE',
        help=(
            'Keep SNPs with INFO at or above this value. Default: 0.9 when an INFO column is selected; '
            'without one the default filter is skipped, and an explicit value is an error. Use --info or '
            '--info-list to select a nonstandard column.'
        ),
    )
    filters.add_argument(
        '--maf-min', default=None, type=float, metavar='VALUE',
        help=(
            'Keep SNPs whose minor-allele frequency min(FRQ, 1-FRQ) is at or above this value. Default: '
            '0.01 when a frequency column is selected; without one the default filter is skipped, and an '
            'explicit value is an error. Output FRQ preserves the original frequency.'
        ),
    )

    columns.add_argument(
        '--snp', default=None, type=str, metavar='COLUMN',
        help=(
            'Input column containing SNP identifiers; matched case-insensitively. If omitted, detect from '
            'headers. SNP identifiers are required for rsid and rsid_allele_aware identity.'
        ),
    )
    columns.add_argument(
        '--chr', default=None, type=str, metavar='COLUMN',
        help=(
            'Input chromosome column; matched case-insensitively. If omitted, detect from headers. Required '
            'in the input for coordinate-family modes; pair with a position column.'
        ),
    )
    columns.add_argument(
        '--pos', default=None, type=str, metavar='COLUMN',
        help=(
            'Input base-pair position column; matched case-insensitively. If omitted, detect from headers. '
            'Required in the input for coordinate-family modes; pair with a chromosome column.'
        ),
    )
    columns.add_argument(
        '--a1', default=None, type=str, metavar='COLUMN',
        help=(
            'Input effect-allele column, relative to which the signed statistic is defined; matched '
            'case-insensitively. If omitted, detect from headers. Allele-aware identity requires both A1 '
            'and A2 in the input.'
        ),
    )
    columns.add_argument(
        '--a2', default=None, type=str, metavar='COLUMN',
        help=(
            'Input other-allele column; matched case-insensitively. If omitted, detect from headers. '
            'Allele-aware identity requires both A1 and A2 in the input.'
        ),
    )
    columns.add_argument(
        '--p', default=None, type=str, metavar='COLUMN',
        help=(
            'Input p-value column used to calculate Z; matched case-insensitively. If omitted, detect from '
            'headers. Input p-values must be in (0, 1].'
        ),
    )
    columns.add_argument(
        '--signed-sumstats', default=None, type=str, metavar='COLUMN,NULL',
        help=(
            'Signed-statistic column and null value, separated by a comma; for example Z,0 or OR,1. Column '
            'names are case-insensitive and effects are relative to A1. If omitted, detect from headers; '
            '--a1-inc bypasses signed-statistic use.'
        ),
    )
    columns.add_argument(
        '--a1-inc', default=False, action='store_true',
        help=(
            'Assert that A1 increases the trait for every SNP and derive positive Z from p-values, '
            'bypassing signed statistics. Default: off; use the detected or --signed-sumstats column to '
            'determine direction.'
        ),
    )
    columns.add_argument(
        '--frq', default=None, type=str, metavar='COLUMN',
        help=(
            'Input allele-frequency or MAF column; matched case-insensitively. If omitted, detect from '
            'headers; if absent, omit FRQ and its filter. Selected values are preserved as FRQ without '
            'conversion to MAF.'
        ),
    )
    columns.add_argument(
        '--info', default=None, type=str, metavar='COLUMN',
        help=(
            'Input scalar INFO column; matched case-insensitively. If omitted, detect from headers; if '
            'absent, skip INFO filtering. Select one INFO source, using --info-list for comma-separated '
            'values within a cell.'
        ),
    )
    columns.add_argument(
        '--info-list', default=None, type=str, metavar='COLUMNS',
        help=(
            (
            'Comma-separated names of input columns whose cells contain comma-separated INFO values, such '
            'as 0.852,0.113,NA. Average nonmissing values for --info-min filtering; names are '
            'case-insensitive. If omitted, detect supported INFO-list columns; competing INFO sources are '
            'rejected.'
        )
        ),
    )
    columns.add_argument(
        '--nstudy', default=None, type=str, metavar='COLUMN',
        help=(
            'Input column counting studies contributing to each SNP; matched case-insensitively. If '
            'omitted, detect from headers. --nstudy-min filtering applies only when per-SNP sample sizes '
            'are absent.'
        ),
    )
    columns.add_argument(
        '--ignore', default=None, type=str, metavar='COLUMNS',
        help=(
            'Comma-separated input column names to exclude from detection and reading; matched '
            'case-insensitively. Default: ignore no columns. Do not also select an ignored column with an '
            'explicit column flag.'
        ),
    )

    computing.add_argument(
        '--chunksize', default=1000000, type=int, metavar='N',
        help=(
            'Read this many input rows per chunk. Larger chunks use more memory. Default: 1,000,000.'
        ),
    )

    logs.add_argument(
        '--log-level', default='INFO', choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'),
        help=(LOG_LEVEL_HELP + ' The selected liftover method and mapping/drop summary are always recorded '
              'in diagnostics/sumstats.log and printed to stdout after a successful run.'),
    )
    return public


def infer_raw_sumstats(
    raw_sumstats_file: str | PathLike[str],
    raw_config: MungeConfig | None = None,
    munge_config: MungeConfig | None = None,
    global_config: GlobalConfig | None = None,
) -> RawSumstatsInference:
    """Infer raw summary-statistics format and minimal parser hints.

    The inference pass reads the header and first data row only. It detects
    plain text, old DANER, and new DANER inputs; reports missing required
    fields; identifies numeric/NA comma-separated INFO lists; and returns exact
    CLI hints for cases that still need user confirmation. VCF-style headers
    with leading ``##`` metadata are treated as plain input. It intentionally
    does not treat ``NEFF`` as total sample size ``N``. Missing
    ``A1/A2`` is reported only when the resolved ``global_config`` uses an
    allele-aware SNP identifier mode; base modes ignore alleles and can munge
    raw summary statistics without allele columns.
    """
    raw_config = raw_config or MungeConfig(raw_sumstats_file=raw_sumstats_file)
    munge_config = munge_config or MungeConfig()
    global_config = global_config or get_global_config()
    path = resolve_scalar_path(raw_sumstats_file, label="raw summary statistics")
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

    _append_suggested_option(suggested_args, "--input-format", detected_format)
    if detected_format == "plain":
        if "REF" in clean_set and "ALT" in clean_set and not {"A1", "A2", "EA", "NEA"} & clean_set:
            column_hints.setdefault("a1", _original_column(file_cnames, "REF"))
            column_hints.setdefault("a2", _original_column(file_cnames, "ALT"))
    elif detected_format == "daner-new":
        frq_u_column = _first_column_with_clean_prefix(file_cnames, "FRQ_U_")
        if frq_u_column is not None:
            column_hints.setdefault("frq", frq_u_column)

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

    requires_alleles = is_allele_aware_mode(global_config.snp_identifier)
    has_alleles = {"A1", "A2"}.issubset(translated)
    if requires_alleles and not has_alleles:
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
    skiprows = munge_input.count_leading_sumstats_comment_lines(path)
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
    targets.update(munge_input.COLUMN_HINT_TARGETS[key] for key, value in column_hints.items()
                   if value is not None and key in munge_input.COLUMN_HINT_TARGETS)
    return targets


def _apply_raw_sumstats_inference(
    source_path: str,
    raw_config: MungeConfig,
    munge_config: MungeConfig,
) -> tuple[MungeConfig, MungeConfig, RawSumstatsInference]:
    """Return configs updated with auto-inferred raw-format hints."""
    inference = infer_raw_sumstats(source_path, raw_config, munge_config)
    column_hints = {**inference.column_hints, **raw_config.column_hints}
    raw_config = replace(raw_config, column_hints=column_hints)
    munge_config = replace(
        munge_config,
        sumstats_format=inference.detected_format,
        info_list_columns=tuple(dict.fromkeys([*munge_config.info_list_columns, *inference.info_list_columns])),
        ignore_columns=tuple(dict.fromkeys([*munge_config.ignore_columns, *inference.ignore_columns])),
    )
    return raw_config, munge_config, inference


def _apply_build_inference_report(
    source_path: str,
    raw_config: MungeConfig,
    munge_config: MungeConfig,
    global_config: GlobalConfig,
    inference: RawSumstatsInference,
) -> RawSumstatsInference:
    """Augment ``--infer-only`` output with source/output build validation."""
    missing = list(inference.missing_fields)
    suggested = list(inference.suggested_args)
    notes = list(inference.notes)
    _append_suggested_option(suggested, "--snp-identifier", global_config.snp_identifier)
    if munge_config.no_snp_restriction:
        _append_suggested_flag(suggested, "--no-snp-restriction")
        notes.append("SNP restriction: none; ordinary QC still applies.")
    elif munge_config.sumstats_snps_file is not None:
        _append_suggested_option(suggested, "--sumstats-snps-file", munge_config.sumstats_snps_file)
        notes.append(f"SNP restriction: custom list ({munge_config.sumstats_snps_file}); replaces packaged HM3.")
    else:
        notes.append("SNP restriction: packaged HM3 (default).")
    if identity_mode_family(global_config.snp_identifier) == "rsid":
        return replace(inference, suggested_args=tuple(suggested), notes=tuple(notes))
    source_hint = normalize_genome_build(munge_config.source_genome_build)
    output_build = normalize_genome_build(munge_config.output_genome_build)
    resolved_source = source_hint
    if output_build in {"hg19", "hg38"}:
        _append_suggested_option(suggested, "--output-genome-build", output_build)
    try:
        hints = {**inference.column_hints, **raw_config.column_hints}
        flags = {clean_header(value): target.upper() for target, value in hints.items() if target in {"chr", "pos"}}
        mapping = get_cname_map(flags, default_cnames, munge_config.ignore_columns)
        translation = {column: mapping[clean_header(column)] for column in read_header(source_path)
                       if clean_header(column) in mapping}
        _, compression = get_compression(source_path)
        metadata = munge_input._resolve_source_coordinates(
            source_path, translation, replace(global_config, genome_build=source_hint), compression,
            munge_input.count_leading_sumstats_comment_lines(source_path),
            restriction_path=munge_config.sumstats_snps_file if not munge_config.no_snp_restriction else None)
        resolved_source = metadata['genome_build']
    except Exception as exc:
        resolved_source = None
        missing.append("source_genome_build")
        notes.append(f"Unable to resolve a consistent source genome build. {exc}")
    if resolved_source in {"hg19", "hg38"}:
        _append_suggested_option(suggested, "--source-genome-build", resolved_source)
    builds_resolved = resolved_source in {"hg19", "hg38"} and output_build in {"hg19", "hg38"}
    liftover_required = resolved_source != output_build if builds_resolved else None
    request = _liftover_request_from_config(munge_config, resolved_source)
    liftover_method = "none"
    if builds_resolved and not liftover_required and munge_config.liftover_chain_file is not None:
        notes.append("Source and output genome builds match; the supplied liftover method will be ignored.")
    if munge_config.liftover_chain_file is not None:
        _append_suggested_option(suggested, "--liftover-chain-file", munge_config.liftover_chain_file)
    if liftover_required and request.method is None:
        missing.append("liftover_method")
        liftover_method = "missing; chain file required"
        notes.append("Source and output genome builds differ; custom-list or unrestricted SNPs require a chain file.")
        notes.append(f"Add --liftover-chain-file <{_expected_chain_label(resolved_source, output_build)}>.")
    elif liftover_required and request.use_hm3_quick_liftover:
        liftover_method = "hm3 quick"
        notes.append("Using automatic HM3 quick liftover from package-bundled reference HM3 metadata; an explicit chain file overrides it.")
    elif liftover_required and request.liftover_chain_file is not None:
        liftover_method = "chain file"
        notes.append("Using the explicit chain file; automatic HM3 quick liftover is disabled.")
        try:
            resolve_scalar_path(request.liftover_chain_file, label="liftover chain file")
        except LDSCInputError as exc:
            missing.append("liftover_chain_file")
            notes.append(
                f"{exc} For this infer-only report, pass one existing chain file for the expected direction "
                f"({_expected_chain_label(resolved_source, output_build)})."
            )
        notes.append(f"Expected chain direction: {resolved_source} -> {output_build}.")
    missing_fields = tuple(dict.fromkeys(missing))
    return replace(
        inference,
        missing_fields=missing_fields,
        suggested_args=tuple(suggested),
        notes=tuple(notes),
        source_genome_build=resolved_source,
        output_genome_build=output_build,
        liftover_required=liftover_required,
        liftover_method=liftover_method,
    )


def _expected_chain_label(source_build: str, output_build: str) -> str:
    """Return a readable source-to-output chain label for infer-only guidance."""
    return f"{source_build}To{output_build[0].upper()}{output_build[1:]}.over.chain"


def _append_suggested_option(args: list[str], flag: str, value: object) -> None:
    """Append a CLI option/value pair to inference suggestions once."""
    if flag not in args:
        args.extend([flag, str(value)])


def _append_suggested_flag(args: list[str], flag: str) -> None:
    """Append a boolean CLI flag to inference suggestions once."""
    if flag not in args:
        args.append(flag)




def _apply_missing_column_report(
    source_path: str,
    raw_config: MungeConfig,
    config: MungeConfig,
    inference: RawSumstatsInference,
) -> RawSumstatsInference:
    """Scan raw columns in bounded chunks and diagnose unusable selected inputs."""
    headers = read_header(source_path)
    config = replace(config, info_list_columns=inference.info_list_columns)
    columns = {munge_input.COLUMN_HINT_TARGETS[key]: value
               for key, value in {**inference.column_hints, **raw_config.column_hints}.items()
               if key in munge_input.COLUMN_HINT_TARGETS}
    flags, _ = munge_input._column_hints(columns, config)
    aliases = {key: target for key, target in default_cnames.items()
               if not (config.signed_sumstats_spec or config.a1_inc) or target not in null_values}
    mapping = munge_input.get_cname_map(flags, aliases, inference.ignore_columns)
    if inference.detected_format == "daner-old":
        frq = _first_column_with_clean_prefix(headers, "FRQ_U_")
        if frq is not None:
            mapping[clean_header(frq)] = "FRQ"
    selected = {column: mapping[clean_header(column)] for column in headers if clean_header(column) in mapping}
    if columns.get("N") is not None or columns.get("N_CAS") is not None:
        munge_input._resolve_sample_size_column_strategy(columns, selected, source_path)
    targets = set(selected.values())
    _, conflicts = munge_input.resolve_qc_thresholds(config, targets)
    if conflicts:
        inference = replace(
            inference,
            notes=(*inference.notes, *(message for _, message in conflicts)),
            missing_fields=(*inference.missing_fields, *(f"column required by {flag}" for flag, _ in conflicts)),
        )
    if "N" in targets or {"N_CAS", "N_CON"}.issubset(targets):
        selected = {column: target for column, target in selected.items() if target != "NSTUDY"}

    empty = set(selected)
    if not empty:
        return inference
    info_lists = {clean_header(column) for column in inference.info_list_columns}
    _, compression = get_compression(source_path)
    with pd.read_csv(
        source_path, sep=r"\s+", compression=compression, na_values=[".", "NA"],
        skiprows=munge_input.count_leading_sumstats_comment_lines(source_path),
        chunksize=config.chunk_size, dtype=str, usecols=list(selected),
    ) as reader:
        for chunk in reader:
            for column in tuple(empty):
                values = chunk[column]
                if clean_header(column) in info_lists:
                    values = values.map(lambda value: kernel_munge._mean_info_list_value(value, column))
                if values.notna().any():
                    empty.remove(column)
            if not empty:
                break
    if not empty:
        return inference
    names = [column for column in headers if column in empty]
    notes = [*inference.notes, "All-missing columns: " + ", ".join(names)]
    notes.append("Selected columns contain no nonmissing values: " + ", ".join(names)
                 + ". Supply populated columns with column overrides or repair the input; use --ignore only for optional columns and provide any required replacement.")
    return replace(inference, notes=tuple(notes), missing_fields=(
        *inference.missing_fields, *(f"all-missing column: {column}" for column in names)))


def _preserve_inference_options(inference: RawSumstatsInference, args: argparse.Namespace) -> RawSumstatsInference:
    """Retain effective CLI options alongside resolved inference suggestions."""
    suggested = list(inference.suggested_args)
    for action in build_parser()._actions:
        if not action.option_strings or action.dest in {"help", "infer_only", "raw_sumstats_file", "output_dir"}:
            continue
        value = getattr(args, action.dest, action.default)
        if value == action.default:
            value = inference.column_hints.get(action.dest, value)
        if value is None or value == action.default:
            continue
        flag = action.option_strings[0]
        if action.nargs == 0:
            if value:
                _append_suggested_flag(suggested, flag)
        else:
            _append_suggested_option(suggested, flag, value)
    return replace(inference, suggested_args=tuple(suggested))


def _render_inference_report(
    inference: RawSumstatsInference,
    raw_sumstats_file: str,
    output_dir: str | Path,
) -> str:
    lines = [
        f"Raw summary-statistics file: {raw_sumstats_file}",
        f"Detected format: {inference.detected_format}",
    ]
    if inference.column_hints:
        lines.append("Column hints: " + ", ".join(f"{key}={value}" for key, value in sorted(inference.column_hints.items())))
    if inference.signed_sumstats_spec is not None:
        lines.append(f"Signed statistic hint: {inference.signed_sumstats_spec}")
    if inference.info_list_columns:
        lines.append("INFO list columns: " + ", ".join(inference.info_list_columns))
    if inference.source_genome_build is not None:
        lines.append(f"Source genome build: {inference.source_genome_build}")
    if inference.output_genome_build is not None:
        lines.append(f"Output genome build: {inference.output_genome_build}")
    if inference.liftover_required is not None:
        method = inference.liftover_method or "unknown"
        lines.append(f"Liftover required: {'yes' if inference.liftover_required else 'no'} (method: {method})")
    if inference.notes:
        lines.extend(f"Note: {note}" for note in inference.notes)
    command = ["ldsc", "munge-sumstats", *_in_help_order(
        ["--raw-sumstats-file", raw_sumstats_file, "--output-dir", str(output_dir), *inference.suggested_args])]
    lines.append("Next step:")
    lines.append(f"  Runnable: {'yes' if inference.runnable else 'no'}")
    lines.append("  Missing fields: " + (", ".join(inference.missing_fields) if inference.missing_fields else "none"))
    lines.append("  Suggested command:")
    lines.append(_format_shell_command(command, base_indent="    ", option_indent="      "))
    return "\n".join(lines)


def _in_help_order(args: list[str]) -> list[str]:
    """Reorder option and value tokens to the order ``munge-sumstats --help`` lists them."""
    parser = build_parser()
    help_actions = [action for group in parser._action_groups for action in group._group_actions]
    actions = {flag: action for action in help_actions for flag in action.option_strings}
    options = []
    idx = 0
    while idx < len(args):
        width = 1 if actions[args[idx]].nargs == 0 else 2
        options.append(args[idx:idx + width])
        idx += width
    options.sort(key=lambda option: help_actions.index(actions[option[0]]))
    return [token for option in options for token in option]


def _format_shell_command(command: list[str], *, base_indent: str = "  ", option_indent: str = "    ") -> str:
    """Return a copy-pasteable multi-line shell command."""
    if len(command) <= 2:
        return base_indent + " ".join(shlex.quote(part) for part in command)
    lines = [f"{base_indent}{shlex.quote(command[0])} {shlex.quote(command[1])} \\"]
    idx = 2
    while idx < len(command):
        token = command[idx]
        if token.startswith("--") and idx + 1 < len(command) and not command[idx + 1].startswith("--"):
            rendered = f"{shlex.quote(token)} {shlex.quote(command[idx + 1])}"
            idx += 2
        else:
            rendered = shlex.quote(token)
            idx += 1
        suffix = " \\" if idx < len(command) else ""
        lines.append(f"{option_indent}{rendered}{suffix}")
    return "\n".join(lines)


def _normalize_output_format(output_format: str) -> str:
    """Return a validated munged summary-statistics output-format token."""
    if output_format not in _SUMSTATS_OUTPUT_FORMATS:
        raise LDSCUsageError(
            f"munge-sumstats does not support output_format={output_format!r}. "
            "Most likely the output format was misspelled. Use one of 'parquet', 'tsv.gz', or 'both'."
        )
    return output_format


def _sumstats_output_files(
    output_dir: Path, output_format: str, trait_name: str | None = None
) -> dict[str, str]:
    """Resolve trait-named artifacts without allowing labels to create paths."""
    output_format = _normalize_output_format(output_format)
    stem = "sumstats" if trait_name is None else re.sub(r"[^\w.-]+", "_", trait_name).strip("._-") or "trait"
    paths: dict[str, str] = {}
    if output_format in {"parquet", "both"}:
        paths["parquet"] = str(output_dir / f"{stem}.parquet")
    if output_format in {"tsv.gz", "both"}:
        name = "sumstats.gz" if trait_name is None else f"{stem}.sumstats.gz"
        paths["tsv.gz"] = str(output_dir / name)
    return paths


def _prepare_curated_sumstats_frame(data: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with fixed summary-statistics output columns and coordinate placeholders."""
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
        raise LDSCInputError(
            f"munge-sumstats could not write Parquet output because POS contains a non-numeric value {bad_value!r}. "
            "Most likely a coordinate column was mapped incorrectly. Pass the correct --pos column or fix the POS values."
        )

    complete = ~(chr_missing | pos_missing)
    if complete.any():
        complete_pos = positive_int_position_series(
            frame.loc[complete, "POS"],
            context="summary-statistics Parquet output",
            label="POS",
        )
        frame.loc[complete, "CHR"] = normalize_chromosome_series(
            frame.loc[complete, "CHR"],
            context="summary-statistics Parquet output",
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
    """Write gzip TSV with explicit missing fields and unrounded frequency."""
    frame = _prepare_curated_sumstats_frame(data)
    if "FRQ" in frame:
        frame["FRQ"] = frame["FRQ"].astype("string")
    frame.to_csv(path, sep="\t", index=False, float_format="%.3f", compression="gzip", na_rep="NA")


def _sumstats_footer_metadata(
    config_snapshot: GlobalConfig, trait_name: str | None, qc_thresholds: dict[str, float | None] | None = None,
) -> dict[bytes, bytes]:
    """Return discrete ``ldsc:*`` Parquet footer keys for a munged summary-statistics artifact.

    Mirrors the reference-panel footer convention in ``_refpanel/adapters.py``. A
    ``None`` ``genome_build`` (rsID-family) or ``trait_name`` is encoded as an
    empty string and decoded back to ``None`` on read. ``qc_thresholds`` (the run's
    effective ``info_min`` and ``maf_min``) adds ``ldsc:qc_info_min`` and
    ``ldsc:qc_maf_min``; an empty value means that filter was not applied. Both
    keys are omitted when the thresholds are unknown, and readers ignore them.
    """
    identity = identity_artifact_metadata(
        artifact_type="sumstats",
        snp_identifier=config_snapshot.snp_identifier,
        genome_build=config_snapshot.genome_build,
    )
    footer = {
        f"ldsc:{key}".encode("utf-8"): ("" if value is None else str(value)).encode("utf-8")
        for key, value in identity.items()
    }
    footer[b"ldsc:trait_name"] = ("" if trait_name is None else str(trait_name)).encode("utf-8")
    if qc_thresholds is not None:
        for key, value in qc_thresholds.items():
            footer[f"ldsc:qc_{key}".encode("utf-8")] = ("" if value is None else repr(float(value))).encode("utf-8")
    return footer


def _write_sumstats_parquet(
    data: pd.DataFrame, path: str, footer_metadata: dict[bytes, bytes] | None = None
) -> list[dict[str, Any]]:
    """Write snappy-compressed Parquet and return row-group metadata.

    The Parquet payload keeps the munger's numeric precision. One row group is
    emitted per normalized chromosome among complete-coordinate rows. If rows
    without complete coordinates exist, they are emitted as the last row group
    with ``chrom`` recorded as ``None``. When ``footer_metadata`` is supplied,
    its ``ldsc:*`` identity keys are merged into the Parquet schema footer so the
    file is self-describing downstream.
    """
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise LDSCDependencyError(
            "munge-sumstats needs the 'pyarrow' package to write Parquet output, but it is not installed. "
            "Most likely Parquet output was requested in an environment missing pyarrow. "
            "Install it (pip install pyarrow) or choose --output-format tsv.gz."
        ) from exc

    frame = _prepare_sumstats_parquet_frame(data)
    schema = pa.Schema.from_pandas(frame, preserve_index=False)
    if footer_metadata:
        schema = schema.with_metadata({**(schema.metadata or {}), **footer_metadata})
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
    footer_metadata: dict[bytes, bytes] | None = None,
) -> tuple[str, list[dict[str, Any]]]:
    """Write selected artifacts and return the primary path plus Parquet row groups.

    The ``.sumstats.gz`` TSV artifact carries no embedded metadata; only the
    Parquet artifact receives the self-describing ``footer_metadata``.
    """
    output_format = _normalize_output_format(output_format)
    parquet_row_groups: list[dict[str, Any]] = []
    if "tsv.gz" in output_files:
        _write_sumstats_tsv_gz(data, output_files["tsv.gz"])
    if "parquet" in output_files:
        parquet_row_groups = _write_sumstats_parquet(data, output_files["parquet"], footer_metadata)
    primary = output_files["parquet"] if output_format in {"parquet", "both"} else output_files["tsv.gz"]
    return primary, parquet_row_groups


def _empty_sumstats_dropped_snps_frame() -> pd.DataFrame:
    """Return the canonical empty summary-statistics dropped-SNP sidecar frame."""
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
    """Log a count-only summary of the summary-statistics dropped-SNP sidecar."""
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
    """Read one munged summary-statistics artifact according to the public suffix policy."""
    token = str(path)
    if token.endswith(".parquet"):
        return pd.read_parquet(path)
    if token.endswith(".sumstats.gz") or Path(token).name == "sumstats.gz":
        return pd.read_csv(path, sep=r"\s+", compression="gzip")
    if token.endswith(".sumstats"):
        return pd.read_csv(path, sep=r"\s+", compression="infer")
    raise LDSCInputError(
        f"Cannot load munged summary statistics at '{path}': unsupported file suffix. "
        "Most likely this is a raw GWAS file or CSV, not a munged summary-statistics file. "
        "Use a path ending in '.parquet', '.sumstats.gz', or '.sumstats', the default 'sumstats.gz', "
        "or run munge-sumstats first."
    )


def _resolve_curated_sumstats_columns(columns: list[str], *, context: str) -> dict[str, str]:
    """Resolve canonical internal summary-statistics artifact columns from ``columns``."""
    resolved = {
        canonical: resolve_required_column(columns, INTERNAL_SUMSTATS_ARTIFACT_SPEC_MAP[canonical], context=context)
        for canonical in ("SNP", "N", "Z")
    }
    for canonical in ("CHR", "POS", "A1", "A2", "FRQ"):
        actual = resolve_optional_column(columns, INTERNAL_SUMSTATS_ARTIFACT_SPEC_MAP[canonical], context=context)
        if actual is not None:
            resolved[canonical] = actual
    return {canonical: resolved[canonical] for canonical in ("SNP", "CHR", "POS", "N", "Z", "A1", "A2", "FRQ") if canonical in resolved}


def _read_sumstats_parquet_footer(path: str) -> dict[str, Any] | None:
    """Return embedded Parquet identity metadata, or ``None`` when absent.

    Text summary statistics and footerless Parquet both return ``None`` here; the public
    loader distinguishes them and rejects the footerless Parquet case. A footer
    without the ``ldsc:snp_identifier`` or ``ldsc:genome_build`` key is rejected.
    Empty ``genome_build`` and ``trait_name`` values decode to ``None``.
    """
    if not str(path).endswith(".parquet"):
        return None
    import pyarrow.parquet as pq

    raw = pq.read_schema(path).metadata or {}
    if b"ldsc:artifact_type" not in raw:
        return None
    missing = [key.decode("utf-8") for key in (b"ldsc:snp_identifier", b"ldsc:genome_build") if key not in raw]
    if missing:
        raise LDSCInputError(
            _curated_sumstats_artifact_message(
                path,
                f"its parquet footer lacks the {' and '.join(missing)} identity metadata.",
                "Most likely it was written by an older LDSC3 version or edited after munging.",
            )
        )
    genome_build = raw[b"ldsc:genome_build"].decode("utf-8") or None
    trait_name = raw.get(b"ldsc:trait_name", b"").decode("utf-8") or None
    return {
        "artifact_type": raw[b"ldsc:artifact_type"].decode("utf-8"),
        "snp_identifier": raw[b"ldsc:snp_identifier"].decode("utf-8"),
        "genome_build": genome_build,
        "trait_name": trait_name,
    }


def _curated_sumstats_artifact_message(path: str | PathLike[str], detail: str, likely: str) -> str:
    """Return the standard malformed munged-artifact diagnostic."""
    return (
        f"Cannot load munged summary statistics at '{path}': {detail} "
        f"{likely} Re-run munge-sumstats from the raw GWAS file. Other causes & fixes: "
        f"{_CURATED_ARTIFACT_TROUBLESHOOTING}"
    )


def _resolve_sumstats_trait_name(
    explicit_trait_name: str | None,
    metadata: dict[str, Any] | None,
    resolved_path: str,
) -> str | None:
    """Resolve an assigned trait label from the CLI/API argument, then the parquet footer.

    Returns ``None`` when neither supplies a label. Regression outputs apply the
    filename fallback themselves, so ``SumstatsMunger.write_output`` never writes
    a derived label back into a footer.
    """
    try:
        explicit = normalize_trait_name(explicit_trait_name)
    except LDSCConfigError as exc:
        raise LDSCConfigError(
            f"Cannot load munged summary statistics at '{resolved_path}': trait_name argument is blank. "
            "Most likely whitespace was passed as the trait label. "
            "Use a non-empty trait name or omit trait_name to derive it from the artifact."
        ) from exc
    if explicit is not None:
        return explicit
    if isinstance(metadata, dict) and "trait_name" in metadata:
        try:
            metadata_trait = normalize_trait_name(metadata.get("trait_name"))
        except LDSCConfigError as exc:
            raise LDSCInputError(
                _curated_sumstats_artifact_message(
                    resolved_path,
                    "its embedded parquet footer has a blank trait_name field.",
                    "Most likely the artifact was hand-edited or written by an outdated workflow.",
                )
            ) from exc
        return metadata_trait
    return None


def _global_config_from_sumstats_metadata(
    metadata: dict[str, Any],
    *,
    artifact_path: str | PathLike[str],
) -> GlobalConfig:
    """Recreate a GlobalConfig snapshot from a summary-statistics Parquet footer payload."""
    try:
        mode = validate_identity_artifact_metadata(metadata, expected_artifact_type="sumstats")
        return GlobalConfig(
            snp_identifier=mode,
            genome_build=metadata.get("genome_build"),
            log_level="INFO",
        )
    except Exception as exc:
        raise LDSCInputError(
            _curated_sumstats_artifact_message(
                artifact_path,
                "its embedded parquet footer has invalid identity provenance.",
                "Most likely it predates the current schema or was edited after munging.",
            )
        ) from exc


def _effective_sumstats_config(config: GlobalConfig, coordinate_metadata: dict[str, Any]) -> GlobalConfig:
    """Return the config snapshot implied by coordinate provenance."""
    genome_build = coordinate_metadata.get("genome_build") or config.genome_build
    return GlobalConfig(
        snp_identifier=normalize_snp_identifier_mode(coordinate_metadata.get("snp_identifier", config.snp_identifier)),
        genome_build=normalize_genome_build(genome_build),
        log_level=config.log_level,
    )


def _log_sumstats_provenance(
    *,
    coordinate_metadata: dict[str, Any],
    output_format: str,
    output_files: dict[str, str],
    parquet_compression: str | None,
    parquet_row_groups: list[dict[str, Any]],
) -> None:
    """Log detailed provenance that is intentionally excluded from the parquet footer."""
    coordinate_provenance = dict(coordinate_metadata)
    liftover = coordinate_provenance.pop(
        "liftover",
        default_liftover_metadata(
            source_build=coordinate_provenance.get("genome_build"),
            snp_identifier=coordinate_provenance.get("snp_identifier", "chr_pos_allele_aware"),
        ),
    )
    LOGGER.info(
        f"Summary-statistics output bookkeeping: output_format={output_format}; "
        f"output_files={_format_log_mapping(dict(output_files))}; "
        f"parquet_compression={parquet_compression or 'none'}; "
        f"parquet_row_groups={len(parquet_row_groups)}"
    )
    LOGGER.info(
        f"Summary-statistics coordinate provenance: {_format_log_mapping(coordinate_provenance)}"
    )
    if isinstance(liftover, dict):
        LOGGER.info(f"Summary-statistics liftover report: {_format_log_mapping(liftover)}")
        if liftover.get("method") == "hm3_curated":
            LOGGER.info(
                f"Summary-statistics HM3 liftover provenance: method={liftover.get('method')}; "
                f"hm3_map_file={liftover.get('hm3_map_file')}"
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
