"""Workflow wrapper for LDSC3 summary-statistics munging.

Core functionality:
    Expose a typed, package-level interface to the legacy-compatible munging
    implementation, its CLI parser and the ``--infer-only`` report. The munged
    artifact itself (``SumstatsTable``, its loader and the Parquet footer
    codec) lives in ``ldsc._sumstats.table``.

Overview
--------
This module converts the historical munging script behavior into explicit
Python dataclasses and a small service object. The public workflow boundary
accepts path-like inputs and resolves the run's munging plan (raw header,
keep-list path, allele orientation, source build, liftover method) with
``_sumstats.munge_input.resolve_munge_plan``, which ``--infer-only`` and
:func:`infer_raw_sumstats` report from as well;
``munge_input.prepare_munge_input`` raises the plan's problems and loads the
keep-list. It passes a ``ResolvedMungeInput`` to the kernel and consumes an explicit ``MungeResult``.
Packaged HM3 restriction is the default; a custom list replaces it, and
``no_snp_restriction`` disables keep-list filtering. For different coordinate
builds, packaged HM3 uses quick liftover unless an explicit chain overrides it;
custom-list and unrestricted runs require a chain.
The workflow layer owns CLI orchestration, output
preflight, diagnostics, and result objects; the kernel keeps the
legacy-compatible parsing and filtering primitives. Coordinate-family runs separate raw source-build interpretation
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
import shlex
from typing import Any

import pandas as pd

from .._runtime.cli_help import CLIHelpFormatter, SCALAR_PATH_HELP, output_directory_argument
from .._core.columns import normalize_genome_build
from ..config import GlobalConfig, MungeConfig, get_global_config
from ..errors import (
    LDSCInputError,
    LDSCUsageError,
    LDSCUserError,
)
from .._runtime.paths import ensure_output_directory
from .._runtime.workflow_logging import (
    LOG_LEVEL_HELP,
    args_output_scope,
    log_inputs,
    log_outputs,
    log_summary,
    materializing_overwrite_guard,
    workflow_logging,
)
from .._core.snp_identity import (
    concat_identity_drop_frames,
    identity_mode_family,
    normalize_snp_identifier_mode,
)
from . import munge_kernel as kernel_munge
from . import munge_input
from .outputs import (
    INPUT_ISSUES_FILE,
    SUMSTATS_OUTPUT_FORMATS,
    SUMSTATS_PARQUET_COMPRESSION,
    SumstatsDirectoryWriter,
    curated_sumstats_frame,
)
from .table import SumstatsTable, encode_sumstats_footer
from .._runtime.preflight import inspect_declared_inputs


LOGGER = logging.getLogger("LDSC.sumstats_munger")

_RAW_SUMSTATS_FORMATS = {"auto", "plain", "daner-old", "daner-new"}


@dataclass(frozen=True)
class RawSumstatsInference:
    """Preview report for one raw summary-statistics file.

    The report is produced by ``--infer-only`` and :func:`infer_raw_sumstats`
    from the munging plan the run itself resolves
    (``munge_input.resolve_munge_plan``), so ``missing_fields`` names the
    problems the run would stop on before its munging pass (header, keep-list
    path, ``INC_ALLELE``/``DEC_ALLELE`` orientation, source build and
    liftover, from the first stage that finds any), plus selected columns
    without any value. It records only safe,
    format-aware decisions: detected raw format, column hints that can be
    applied without changing statistical meaning, INFO-list handling, missing
    or conflicting fields, and command-line repair suggestions. Suggested commands should prefer
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


class SumstatsMunger:
    """Run legacy-compatible munging through a typed workflow interface."""

    def __init__(self) -> None:
        """Initialize the workflow wrapper and clear any cached run summary."""
        self._last_summary: MungeRunSummary | None = None

    @materializing_overwrite_guard(
        lambda self, config, global_config=None: (config.output_dir, config.overwrite) if config.output_dir else None,
        command="SumstatsMunger.run(...)",
    )
    def run(self, config: MungeConfig, *, global_config: GlobalConfig | None = None) -> SumstatsTable:
        """Munge one raw summary-statistics file into LDSC-ready form.

        Parameters
        ----------
        config : MungeConfig
            The whole munging request: the raw file path, trait label and
            column hints, the QC thresholds, the output directory and the
            munged output format. Common plain-text inputs, including
            VCF-style headers, old DANER, and new DANER are handled before the
            kernel runs; explicit column hints still take priority. The
            workflow writes ``<trait>.parquet`` (self-describing via its footer)
            and/or ``<trait>.sumstats.gz``, using the filesystem-safe trait
            label. Without a label, filenames are
            ``sumstats.parquet`` and ``sumstats.gz``. Diagnostics live under
            ``diagnostics/``. Any existing owned artifact is refused before the kernel
            runs unless ``config.overwrite`` is true; successful
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
        global_config : GlobalConfig or None, keyword-only, optional
            Shared configuration snapshot to attach to the returned
            ``SumstatsTable``. When omitted, the current package-global
            configuration is captured. Default is ``None``.

        Returns
        -------
        SumstatsTable
            Validated, in-memory table suitable for the regression workflow.
            It holds the same rows, in the same order and with the same
            ``CHR`` labels, as the written Parquet and ``.sumstats.gz`` files:
            the columns ``SNP, CHR, POS, A1, A2, Z, N, FRQ`` that exist,
            complete-coordinate rows in chromosome and position order, then
            rows without complete coordinates in input order.
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
        munge_input.require_raw_sumstats_file(config)
        if config.output_dir is None:
            raise LDSCUserError(
                "No output directory given to munge-sumstats. Most likely --output-dir or output_dir was omitted. "
                "Pass --output-dir <directory> on the CLI or set MungeConfig.output_dir."
            )

        config_snapshot = munge_input.resolve_source_config(config, global_config or get_global_config())
        # An explicit source build fixes the liftover method now: stop before any input IO.
        requested_liftover = munge_input.liftover_request(config, config_snapshot.genome_build)
        early_issue = munge_input.liftover_issue(requested_liftover, config_snapshot.genome_build)
        if early_issue is not None:
            early_issue.raise_error()
        output_dir = ensure_output_directory(config.output_dir, label="output directory")
        writer = SumstatsDirectoryWriter()
        family = writer.artifact_family(
            output_dir, trait_name=config.trait_name, output_format=config.output_format, run=True,
        )
        output_files = {name: str(family.paths[name]) for name in ("parquet", "tsv.gz") if name in family.paths}
        log_path = str(family.paths["log"])
        dropped_snps_path = family.paths["dropped_snps"]
        input_issues_path = output_dir / INPUT_ISSUES_FILE
        sumstats_snps_path = munge_input.snp_restriction_file(config)
        sumstats_snps_label = "none" if sumstats_snps_path is None else str(sumstats_snps_path)
        # The writer checks the family again when it publishes, then removes
        # stale owned siblings.
        family.preflight(overwrite=config.overwrite)
        family.paths["log"].parent.mkdir(parents=True, exist_ok=True)
        with workflow_logging("munge-sumstats", log_path, log_level=config_snapshot.log_level):
            inspect_declared_inputs(files=[('raw summary statistics', config.raw_sumstats_file),
                ('SNP restriction', sumstats_snps_path), ('liftover chain', requested_liftover.liftover_chain_file)],
                issues_path=input_issues_path)
            plan = munge_input.resolve_munge_plan(config, config_snapshot)
            source_path = plan.schema.source_path
            column_hints = {**plan.schema.inferred_hints, **config.column_hints}
            log_inputs(
                raw_sumstats_file=source_path,
                output_dir=str(output_dir),
                output_format=config.output_format,
                sumstats_snps_file=sumstats_snps_label,
                snp_restriction=("hm3" if munge_input.uses_packaged_hm3(config) else "none" if config.no_snp_restriction else "custom"),
                hm3_map_file=sumstats_snps_label if munge_input.uses_packaged_hm3(config) else "none",
                output_genome_build=config.output_genome_build or "none",
                liftover_method=("not resolved" if plan.coordinates is None else plan.liftover.method or "none"),
            )
            LOGGER.info(
                f"Munging summary statistics from '{source_path}' into '{output_dir}' "
                f"with snp_identifier='{config_snapshot.snp_identifier}', "
                f"genome_build='{config_snapshot.genome_build}', "
                f"sumstats_snps_file='{sumstats_snps_label}', "
                f"packaged_hm3='{munge_input.uses_packaged_hm3(config)}', "
                f"output_genome_build='{config.output_genome_build}'."
            )
            request = munge_input.prepare_munge_input(plan)
            resolved_liftover = request.liftover_request
            method_label = {
                "hm3_curated": "automatic HM3 quick liftover (package-bundled reference HM3 metadata)",
                "chain_file": "explicit chain file (HM3 quick liftover disabled)",
                None: "none; no liftover required",
            }[resolved_liftover.method]
            LOGGER.info("Selected liftover method: %s; source=%s; output=%s",
                        method_label, request.genome_build, resolved_liftover.target_build)
            result = kernel_munge.munge_sumstats(request)
            # One frame, curated once, is the returned table and both data files.
            data = curated_sumstats_frame(result.data)
            coordinate_metadata = result.coordinate_metadata
            drop_frame = concat_identity_drop_frames([result.liftover_drop_frame, result.identity_drop_frame])
            table_config_snapshot = _effective_sumstats_config(config_snapshot, coordinate_metadata)
            qc_thresholds = {"info_min": request.qc.info_min, "maf_min": request.qc.maf_min}
            # Always write the audit sidecar so consumers can distinguish a
            # clean run from a missing or stale dropped-SNP artifact.
            written = writer.write(
                data,
                output_dir,
                trait_name=config.trait_name,
                output_format=config.output_format,
                footer_metadata=encode_sumstats_footer(
                    table_config_snapshot, config.trait_name, qc_thresholds,
                ),
                dropped_snps=drop_frame,
                overwrite=config.overwrite,
                preflighted=True,
            )
            _log_sumstats_provenance(
                coordinate_metadata=coordinate_metadata,
                output_format=config.output_format,
                output_files=output_files,
                parquet_compression=(
                    SUMSTATS_PARQUET_COMPRESSION
                    if "parquet" in output_files
                    else None
                ),
                parquet_row_groups=written.parquet_row_groups,
            )
            _log_sumstats_dropped_snps_summary(drop_frame, dropped_snps_path)
            table = SumstatsTable(
                data=data,
                has_alleles={"A1", "A2"}.issubset(data.columns),
                source_path=source_path,
                trait_name=config.trait_name,
                provenance={
                    "raw_sumstats_file": source_path,
                    "output_dir": str(output_dir),
                    "output_format": config.output_format,
                    "output_files": dict(output_files),
                    "column_hints": dict(column_hints),
                    "coordinate_provenance": coordinate_metadata,
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
                inferred_columns={**column_hints, "detected_format": plan.schema.detected_format},
                used_n_rule=result.used_n_rule,
                output_paths=run_output_paths,
                qc_thresholds={"INFO": request.qc.info_min, "FRQ": request.qc.maf_min},
            )
            log_summary(_render_munge_summary(self._last_summary, coordinate_metadata, config, table.has_alleles))
            log_outputs(**run_output_paths)
            LOGGER.info(
                f"Munged {self._last_summary.n_input_rows} input rows to {self._last_summary.n_retained_rows} retained rows; "
                f"wrote '{written.primary}'."
            )
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
            Validated LDSC-ready table to persist. Its rows are written in the
            form ``curated_sumstats_frame`` gives munged rows: columns in the
            package output order ``SNP, CHR, POS, A1, A2, Z, N, FRQ`` where
            present, missing ``CHR`` or ``POS`` columns materialized as
            missing values, canonical ``CHR`` labels and genomic row order. An
            unsupported ``CHR`` label or invalid ``POS`` on a complete-coordinate
            row raises ``LDSCInputError`` before any file is written.
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
        return SumstatsDirectoryWriter().write(
            curated_sumstats_frame(sumstats.data),
            output_dir,
            trait_name=sumstats.trait_name,
            output_format=output_format,
            footer_metadata=encode_sumstats_footer(
                sumstats.config_snapshot, sumstats.trait_name, sumstats.provenance.get("qc_thresholds"),
            ),
            overwrite=overwrite,
        ).primary

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
    the Python API takes, then delegates to :class:`SumstatsMunger`
    so path resolution, output preflight, metadata, and result construction
    stay in one workflow path. Successful materializing runs print the resolved
    liftover method and mapping/drop summary to stdout at every log level;
    direct Python calls to ``SumstatsMunger.run`` remain console-quiet.
    A failed authorized overwrite, including an invalid option combination,
    leaves a failure marker; ``infer_only`` runs and the printed summary of a
    completed run are outside the marker's scope.

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
    if not getattr(args, "infer_only", False):
        munger, table, config = _munge_from_args(args)
        print(_render_munge_summary(
            munger.build_run_summary(), table.provenance["coordinate_provenance"], config, table.has_alleles,
        ))
        return table
    plan, inference = _infer(_checked_munge_config(args), _resolve_main_global_config(args))
    inference = _preserve_inference_options(inference, args)
    print(_render_inference_report(inference, plan.schema.source_path, args.output_dir))
    return inference


@materializing_overwrite_guard(args_output_scope, command="run_munge_sumstats_from_args(...)")
def _munge_from_args(args: argparse.Namespace) -> tuple[SumstatsMunger, SumstatsTable, MungeConfig]:
    """Run a materializing munging request; the caller prints its summary once the run is complete."""
    config = _checked_munge_config(args)
    munger = SumstatsMunger()
    table = munger.run(config, global_config=_resolve_main_global_config(args))
    return munger, table, config


def _checked_munge_config(args: argparse.Namespace) -> MungeConfig:
    """Return the munging config of parsed arguments whose input and output paths are nonempty."""
    if not str(getattr(args, "raw_sumstats_file", "") or "").strip():
        raise LDSCUsageError(
            "munge-sumstats requires a nonempty --raw-sumstats-file, including with --infer-only. "
            "Pass --raw-sumstats-file <file>; if using a shell variable, make sure it contains the input path."
        )
    if not str(getattr(args, "output_dir", "") or "").strip():
        raise LDSCUsageError("munge-sumstats requires a nonempty --output-dir. Pass --output-dir <directory>, including with --infer-only; dry runs write nothing.")
    return _munge_config_from_args(args)


def _render_munge_summary(
    summary: MungeRunSummary, metadata: dict[str, Any], config: MungeConfig, has_alleles: bool,
) -> str:
    """Render resolved liftover outcomes, exclusive whole-run drop counts and a missing-allele warning."""
    restriction = "packaged HM3 (default)" if munge_input.uses_packaged_hm3(config) else "none" if config.no_snp_restriction else "custom list"
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
    if not has_alleles:
        lines.append("  No A1/A2 columns written: rg cannot use this output.")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> SumstatsTable | RawSumstatsInference:
    """CLI entry point: parse arguments and delegate to workflow orchestration."""
    return run_munge_sumstats_from_args(build_parser().parse_args(argv))


def _resolve_main_global_config(args: argparse.Namespace) -> GlobalConfig:
    """Return the SNP identifier mode and log level of parsed arguments.

    Build and liftover options stay in ``MungeConfig``;
    ``munge_input.resolve_source_config`` checks them against the mode for
    the run and ``--infer-only`` alike.
    """
    return GlobalConfig(
        snp_identifier=normalize_snp_identifier_mode(getattr(args, "snp_identifier", "chr_pos_allele_aware")),
        log_level=getattr(args, "log_level", "INFO"),
    )


def _munge_config_from_args(args: argparse.Namespace) -> MungeConfig:
    """Convert parsed CLI arguments into the munging configuration."""
    return MungeConfig(
        raw_sumstats_file=args.raw_sumstats_file,
        trait_name=getattr(args, "trait_name", None),
        column_hints=_column_hints_from_args(args),
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
        '--output-dir', required=True, metavar='DIR', type=output_directory_argument,
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
        '--output-format', choices=sorted(SUMSTATS_OUTPUT_FORMATS), default='parquet',
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
            'Input column of A1, the effect allele (the allele the signed statistic refers to); matched '
            'case-insensitively. If omitted, only A1, EA, EFFECT_ALLELE or INC_ALLELE is detected; REF/ALT and '
            'ALLELE0/1/2 never are. Allele-aware identity requires both A1 and A2 in the input.'
        ),
    )
    columns.add_argument(
        '--a2', default=None, type=str, metavar='COLUMN',
        help=(
            'Input column of A2, the other allele; matched case-insensitively. If omitted, only A2, NEA, '
            'OTHER_ALLELE, NON_EFFECT_ALLELE or DEC_ALLELE is detected. Allele-aware identity requires both A1 and '
            'A2 in the input.'
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


def infer_raw_sumstats(config: MungeConfig, *, global_config: GlobalConfig | None = None) -> RawSumstatsInference:
    """Preview a munging run: detected format, parser hints and every problem that would stop it.

    This is the Python form of ``munge-sumstats --infer-only`` and returns
    the same report. ``config`` is the ``MungeConfig`` the run would take;
    its ``raw_sumstats_file`` names the input and every other option is read
    from it as the run reads it. The preview resolves the munging plan the
    run resolves (``munge_input.resolve_munge_plan``): the raw header and
    first data row and the keep-list path; the signed statistic, read in
    chunks, when ``INC_ALLELE``/``DEC_ALLELE`` columns are read as A1/A2; in
    coordinate-family modes the source build, sampling coordinate rows when
    it is ``auto``, and the liftover method; then it scans the selected
    columns in bounded chunks for columns without any value. The plan stops
    at the first of these stages that finds a problem, as the run does. The
    preview writes nothing, logs nothing and does not load the keep-list.

    The report lists the problems of that stage: every header problem that
    would stop the run (missing or conflicting columns, ambiguous sample-size
    or signed-statistic columns, unusable old DANER count columns) and a
    keep-list path that does not resolve; allele columns whose names the
    signed statistic contradicts; an unresolved source build; or a missing
    liftover method or chain file. It also lists all-missing columns. It
    detects plain text, old DANER, and new DANER inputs; identifies
    numeric/NA comma-separated INFO lists; and returns exact CLI hints for
    cases that still need user confirmation. VCF-style headers with leading
    ``##`` metadata are treated as plain input. It intentionally does not
    treat ``NEFF`` as total sample size ``N``. Missing ``A1/A2`` is reported
    only when the resolved ``global_config`` uses an allele-aware SNP
    identifier mode; base modes ignore alleles and can munge raw summary
    statistics without allele columns.

    Raises
    ------
    LDSCUsageError
        If build or liftover options do not fit the SNP identifier mode, or a
        coordinate-family mode has no output build, as the run raises.
    """
    return _infer(config, global_config or get_global_config(), logger=None)[1]


def _infer(
    config: MungeConfig, global_config: GlobalConfig, *, logger=LOGGER,
) -> tuple[munge_input.MungePlan, RawSumstatsInference]:
    """Resolve the run's munging plan for ``config`` and report it, as ``--infer-only`` and :func:`infer_raw_sumstats` do.

    ``logger`` receives the source-build check's report; the CLI keeps it, so
    its console shows a failed check as the run's does, and the Python preview
    passes ``None``.
    """
    plan = munge_input.resolve_munge_plan(
        config, munge_input.resolve_source_config(config, global_config), logger=logger,
    )
    inference = _apply_build_inference_report(plan, _inference_from_schema(plan.schema))
    return plan, _apply_missing_column_report(plan.schema, inference)


def _inference_from_schema(schema: munge_input.RawSchema) -> RawSumstatsInference:
    """Report a raw schema: detected format, inferred hints, and each issue as a missing field noted with the run's message."""
    suggested_args = ["--input-format", schema.detected_format]
    notes = [*schema.notes, *(issue.message for issue in schema.issues)]
    signed_sumstats_spec = None
    if "signed statistic" in {issue.label for issue in schema.issues} and schema.config.signed_sumstats_spec is None:
        signed_column = munge_input.likely_signed_sumstat_column(schema.headers)
        if signed_column is not None:
            signed_sumstats_spec = f"{signed_column},0"
            suggested_args.extend(["--signed-sumstats", signed_sumstats_spec])
    return RawSumstatsInference(
        detected_format=schema.detected_format,
        column_hints=dict(schema.inferred_hints),
        signed_sumstats_spec=signed_sumstats_spec,
        ignore_columns=tuple(dict.fromkeys(schema.config.ignore_columns)),
        info_list_columns=schema.config.info_list_columns,
        missing_fields=tuple(dict.fromkeys(issue.label for issue in schema.issues)),
        suggested_args=tuple(suggested_args),
        notes=tuple(notes),
    )


def _apply_build_inference_report(plan: munge_input.MungePlan, inference: RawSumstatsInference) -> RawSumstatsInference:
    """Add the plan's SNP restriction, source and output builds and liftover method to a preview."""
    config, mode = plan.config, plan.source_config.snp_identifier
    missing = list(inference.missing_fields)
    suggested = list(inference.suggested_args)
    notes = list(inference.notes)
    _append_suggested_option(suggested, "--snp-identifier", mode)
    if config.no_snp_restriction:
        _append_suggested_flag(suggested, "--no-snp-restriction")
        notes.append("SNP restriction: none; ordinary QC still applies.")
    elif config.sumstats_snps_file is not None:
        _append_suggested_option(suggested, "--sumstats-snps-file", config.sumstats_snps_file)
        notes.append(f"SNP restriction: custom list ({config.sumstats_snps_file}); replaces packaged HM3.")
    else:
        notes.append("SNP restriction: packaged HM3 (default).")
    issue = None if plan.issue is None else plan.issue.label
    if issue in {"sumstats_snps_file", "A1/A2 orientation"}:
        missing.append(issue)
        notes.append(plan.issue.message)
    if identity_mode_family(mode) == "rsid":
        return replace(inference, missing_fields=tuple(dict.fromkeys(missing)), suggested_args=tuple(suggested),
                       notes=tuple(notes))
    output_build = config.output_genome_build
    resolved_source = plan.source_build
    if output_build in {"hg19", "hg38"}:
        _append_suggested_option(suggested, "--output-genome-build", output_build)
    if issue == "source_genome_build":
        missing.append("source_genome_build")
        notes.append(f"Unable to resolve a consistent source genome build. {plan.issue.message}")
    if resolved_source in {"hg19", "hg38"}:
        _append_suggested_option(suggested, "--source-genome-build", resolved_source)
    builds_resolved = resolved_source in {"hg19", "hg38"} and output_build in {"hg19", "hg38"}
    liftover_required = resolved_source != output_build if builds_resolved else None
    request = plan.liftover
    liftover_method = "none"
    if builds_resolved and not liftover_required and config.liftover_chain_file is not None:
        notes.append("Source and output genome builds match; the supplied liftover method will be ignored.")
    if config.liftover_chain_file is not None:
        _append_suggested_option(suggested, "--liftover-chain-file", config.liftover_chain_file)
    if issue == "liftover_method":
        missing.append("liftover_method")
        liftover_method = "missing; chain file required"
        notes.append("Source and output genome builds differ; custom-list or unrestricted SNPs require a chain file.")
        notes.append(f"Add --liftover-chain-file <{_expected_chain_label(resolved_source, output_build)}>.")
    elif request.method == "hm3_curated":
        liftover_method = "hm3 quick"
        notes.append("Using automatic HM3 quick liftover from package-bundled reference HM3 metadata; an explicit chain file overrides it.")
    elif request.method == "chain_file":
        liftover_method = "chain file"
        notes.append("Using the explicit chain file; automatic HM3 quick liftover is disabled.")
        if issue == "liftover_chain_file":
            missing.append("liftover_chain_file")
            notes.append(
                f"{plan.issue.message} For this infer-only report, pass one existing chain file for the "
                f"expected direction ({_expected_chain_label(resolved_source, output_build)})."
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
    schema: munge_input.RawSchema,
    inference: RawSumstatsInference,
) -> RawSumstatsInference:
    """Scan the selected raw columns in bounded chunks and diagnose those that hold no value."""
    names = munge_input.all_missing_columns(schema)
    if not names:
        return inference
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
    liftover = coordinate_provenance.pop("liftover")
    LOGGER.info(
        f"Summary-statistics output bookkeeping: output_format={output_format}; "
        f"output_files={_format_log_mapping(dict(output_files))}; "
        f"parquet_compression={parquet_compression or 'none'}; "
        f"parquet_row_groups={len(parquet_row_groups)}"
    )
    LOGGER.info(
        f"Summary-statistics coordinate provenance: {_format_log_mapping(coordinate_provenance)}"
    )
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
