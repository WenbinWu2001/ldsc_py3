"""Explicitly convert selected legacy LD-score suites to LDSC3 directories.

The converter is the public import boundary for legacy LD-score suites:
reference and weight ``.l2.ldscore(.gz)`` families, counts, and, for baseline
partitioned suites, annotations and frequencies. Regression consumes the
resulting canonical :class:`ldsc._ldscore.calculator.LDScoreResult` directory.
Standalone annotation workflows retain their separate ``.annot(.gz)`` readers.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import gzip
import logging
from pathlib import Path
import re

import numpy as np
import pandas as pd

from .._runtime.serial_calculation import serial_calculation
from .._runtime.cli_help import CLIHelpFormatter
from .._runtime.workflow_logging import LOG_LEVEL_HELP
from .._runtime.workflow_logging import log_inputs, log_outputs, materializing_overwrite_guard, workflow_logging
from .._core.columns import normalize_genome_build
from .._core.snp_identity import normalize_snp_identifier_mode
from ..config import GlobalConfig, get_global_config
from ..errors import LDSCInputError, LDSCUsageError
from .calculator import LDScoreResult
from .overlap_table import LDScoreOverlap
from .outputs import LDScoreDirectoryWriter, LDScoreOutputConfig
from .._runtime.paths import ensure_output_directory, normalize_path_token


LOGGER = logging.getLogger("LDSC.legacy_ldscore_converter")
_CHROMS = tuple(range(1, 23))
_LDSCORE_RE = re.compile(r"^(?P<prefix>.*?)(?P<chrom>[1-9]|1[0-9]|2[0-2])\.l2\.ldscore(?P<gz>\.gz)?$")
_LDSCORE_FILENAME_GUIDANCE = (
    "Use '<prefix><chrom>.l2.ldscore(.gz)', e.g. '1.l2.ldscore.gz' or 'weights.1.l2.ldscore.gz', "
    "with one constant prefix (possibly empty) and unpadded chromosomes 1-22. "
    "Organize one complete family directly in the supplied directory; put different releases in separate directories. "
    "Restore missing shards or rename misnamed copies to the standard filenames. "
    "See docs/user-guide/troubleshooting.md#convert-legacy-ldscores."
)
_WEIGHT_FILENAME_GUIDANCE = (
    "Found unsupported '.w.l2.ldscore(.gz)' filenames. Rename copies in a separate input directory, "
    "e.g. 1.w.l2.ldscore.gz -> weights.1.l2.ldscore.gz (likewise without .gz); "
    "the converter does not rename source files. " + _LDSCORE_FILENAME_GUIDANCE
)
_COUNT_FILENAME_GUIDANCE = (
    "Place plain-text '<reference-prefix><chrom>.l2.M_5_50' files in the reference directory "
    "with the same prefix as the reference LD scores. Restore the matching counts from the same source release "
    "or rename a misnamed copy to the expected filename; .l2.M is not a substitute."
)
_FREQUENCY_RE = re.compile(r"^(?P<prefix>.*?)(?P<chrom>[1-9]|1[0-9]|2[0-2])\.frq(?P<gz>\.gz)?$")
_ISSUE_COLUMNS = [
    "severity",
    "source_role",
    "file",
    "chromosome",
    "SNP",
    "annotation",
    "reason",
    "observed",
    "expected",
    "details",
]
_METADATA_COLUMNS = {"CHR", "SNP", "BP", "POS", "CM", "MAF", "GENOME_BUILD"}
_COUNT_CONFIG = {
    "common_reference_snp_maf_min": 0.05,
    "common_reference_snp_maf_operator": ">",
    "common_reference_snp_semantics": "legacy_ldsc",
}


@dataclass(frozen=True)
class LegacyLDScoreConversionResult:
    """Summary of one successful legacy LDSC LD-score conversion.

    Parameters
    ----------
    profile : {"unpartitioned", "baseline_partitioned"}
        Structurally detected suite profile.
    n_rows : int
        Regression SNPs retained after the reference/weight rsID inner join.
    baseline_columns : tuple of str
        Canonical LD-score columns written to the baseline parquet table.
    output_paths : dict of str to str
        Files written in the canonical output directory.
    """

    profile: str
    n_rows: int
    baseline_columns: tuple[str, ...]
    output_paths: dict[str, str]


class LegacyLDScoreConverter:
    """Convert a reusable legacy LDSC reference/weight suite into LDSC3 format.

    Source directories are treated as immutable. Matching across legacy
    families uses rsID; output coordinates come from the reference LD-score
    table. The converter accepts only allele-unaware ``rsid`` and ``chr_pos``
    SNP identifier modes and fixes legacy common-SNP count semantics at ``0.05``.
    """

    def convert(
        self,
        *,
        legacy_reference_dir: str | Path,
        legacy_weight_dir: str | Path,
        output_dir: str | Path,
        legacy_frequency_dir: str | Path | None = None,
        snp_identifier: str = "rsid",
        genome_build: str | None = "auto",
        overwrite: bool = False,
        log_level: str | None = None,
    ) -> LegacyLDScoreConversionResult:
        """Convert one complete autosomal legacy suite.

        Conversion runs in one worker with 1 BLAS thread, including
        annotation overlap reconstruction for baseline-partitioned suites.

        Parameters
        ----------
        legacy_reference_dir, legacy_weight_dir : path-like
            Directories containing exactly one coherent autosomal LD-score
            family for the reference and regression-weight roles.
        output_dir : path-like
            New canonical LDSC3 LD-score directory.
        legacy_frequency_dir : path-like, optional
            Required only for baseline partitioned conversion.
        snp_identifier : {"rsid", "chr_pos"}, optional
            Allele-unaware output identity. Default is ``"rsid"``.
        genome_build : {"auto", "hg19", "hg38"}, optional
            Build provenance. ``chr_pos`` with ``"auto"`` requires decisive
            inference from reference coordinates. Default is ``"auto"``.
        overwrite : bool, optional
            Replace converter-owned outputs when true. Source files are never
            changed. Default is false.
        log_level : {"DEBUG", "INFO", "WARNING", "ERROR"}, optional
            Workflow log threshold. Default is ``None``: use the registered
            ``GlobalConfig.log_level`` (``"INFO"`` unless changed).

        Returns
        -------
        LegacyLDScoreConversionResult
            Profile, retained row count, columns, and written paths.

        Raises
        ------
        LDSCInputError
            If family discovery, legacy schemas, counts, coordinates, or
            cross-family invariants are invalid.
        LDSCUsageError
            If an allele-aware output mode is requested.
        """
        # The conversion child starts with the package-default GlobalConfig,
        # so the session's log level is resolved here, in the caller.
        return self._convert(
            legacy_reference_dir=legacy_reference_dir,
            legacy_weight_dir=legacy_weight_dir,
            output_dir=output_dir,
            legacy_frequency_dir=legacy_frequency_dir,
            snp_identifier=snp_identifier,
            genome_build=genome_build,
            overwrite=overwrite,
            log_level=get_global_config().log_level if log_level is None else log_level,
        )

    @serial_calculation
    @materializing_overwrite_guard(
        lambda self, **kwargs: (
            (kwargs["output_dir"], kwargs.get("overwrite", False))
            if kwargs.get("output_dir")
            else None
        ),
        command="LegacyLDScoreConverter.convert(...)",
    )
    def _convert(
        self,
        *,
        legacy_reference_dir: str | Path,
        legacy_weight_dir: str | Path,
        output_dir: str | Path,
        legacy_frequency_dir: str | Path | None,
        snp_identifier: str,
        genome_build: str | None,
        overwrite: bool,
        log_level: str,
    ) -> LegacyLDScoreConversionResult:
        """Run :meth:`convert` with a resolved log level in the worker."""
        mode = normalize_snp_identifier_mode(snp_identifier)
        if mode not in {"rsid", "chr_pos"}:
            raise LDSCUsageError(
                "convert-legacy-ldscores accepts only allele-unaware --snp-identifier rsid or chr_pos; "
                "legacy LD-score suites do not contain authoritative alleles."
            )
        build_hint = normalize_genome_build(genome_build)
        reference_root = _existing_directory(legacy_reference_dir, "legacy reference directory")
        weight_root = _existing_directory(legacy_weight_dir, "legacy weight directory")
        frequency_root = (
            None
            if legacy_frequency_dir is None
            else _existing_directory(legacy_frequency_dir, "legacy frequency directory")
        )
        output_root = ensure_output_directory(output_dir, label="converted LD-score output directory")
        diagnostics = output_root / "diagnostics"
        diagnostics.mkdir(parents=True, exist_ok=True)
        issue_path = diagnostics / "conversion_issues.tsv.gz"
        log_path = diagnostics / "convert-legacy-ldscores.log"
        LDScoreDirectoryWriter.artifact_family(output_root).preflight(
            overwrite=overwrite, additional_paths=[issue_path, log_path],
        )
        issues: list[dict[str, object]] = []
        with workflow_logging("convert-legacy-ldscores", str(log_path), log_level=log_level):
            log_inputs(
                legacy_reference_dir=str(reference_root),
                legacy_weight_dir=str(weight_root),
                legacy_frequency_dir=None if frequency_root is None else str(frequency_root),
            )
            try:
                result = self._convert_validated(
                    reference_root=reference_root,
                    weight_root=weight_root,
                    frequency_root=frequency_root,
                    output_root=output_root,
                    mode=mode,
                    build_hint=build_hint,
                    overwrite=overwrite,
                    issues=issues,
                )
            except Exception as exc:
                if isinstance(exc, LDSCInputError):
                    issues.append(
                        _issue(
                            "error",
                            "converter",
                            reason="conversion_error",
                            details=str(exc),
                        )
                    )
                _write_issues(issue_path, issues)
                raise
            _write_issues(issue_path, issues)
            output_paths = {**result.output_paths, "issues": str(issue_path), "log": str(log_path)}
            log_outputs(**output_paths)
            return LegacyLDScoreConversionResult(
                profile=result.profile,
                n_rows=result.n_rows,
                baseline_columns=result.baseline_columns,
                output_paths=output_paths,
            )

    def _convert_validated(
        self,
        *,
        reference_root: Path,
        weight_root: Path,
        frequency_root: Path | None,
        output_root: Path,
        mode: str,
        build_hint: str | None,
        overwrite: bool,
        issues: list[dict[str, object]],
    ) -> LegacyLDScoreConversionResult:
        reference_family = _discover_ldscore_family(reference_root, "reference", issues)
        weight_family = _discover_ldscore_family(weight_root, "weight", issues)
        reference_prefix = reference_family.prefix
        annotation_family = _discover_associated_family(
            reference_root,
            reference_prefix,
            suffixes=(".annot.gz", ".annot"),
            issues=issues,
        )
        present_annotations = sorted(annotation_family)
        if present_annotations and present_annotations != list(_CHROMS):
            missing = sorted(set(_CHROMS) - set(present_annotations))
            raise LDSCInputError(
                "Legacy reference suite contains a partial annotation family; baseline conversion requires "
                f"all chromosomes 1-22 (missing {missing}) and cannot fall back to unpartitioned conversion. "
                "Place '<reference-prefix><chrom>.annot(.gz)' directly in the reference directory with the "
                f"same prefix as the selected LD scores ('{reference_prefix}'). "
                "Restore missing shards or rename misnamed copies from the matching source suite."
            )

        from .._genome.build_policy import GenomeBuildCheck
        build_check = GenomeBuildCheck(build_hint, context="convert-legacy-ldscores", logger=LOGGER)
        for role, paths in (("reference", reference_family.paths), ("weight", weight_family.paths),
                            ("annotation", annotation_family)):
            for path in paths.values():
                build_check.table(f"legacy {role}", path, infer_unlabeled=(role == "reference" and mode == "chr_pos"))
        if mode == "chr_pos" or build_check.explicit or any(not entry.unchecked for entry in build_check.evidence):
            build_check.resolve()
        first_reference = _read_ldscore_table(reference_family.paths[1], "reference", 1)
        reference_columns = _scientific_columns(first_reference)
        if not reference_columns:
            raise LDSCInputError(
                f"Legacy reference LD-score table '{reference_family.paths[1]}' has no numeric scientific LD-score column. "
                "Supply the original .l2.ldscore(.gz) reference tables containing LD scores, rather than SNP metadata or weight-only files."
            )
        if present_annotations:
            profile = "baseline_partitioned"
            if frequency_root is None:
                raise LDSCInputError(
                    "Baseline partitioned legacy LDSC conversion requires --legacy-frequency-dir so common-SNP counts and "
                    "overlap can be reconstructed with the fixed legacy 0.05 threshold. "
                    "Pass a directory containing '<prefix><chrom>.frq(.gz)' for chromosomes 1-22."
                )
            return self._convert_baseline_partitioned(
                reference_root=reference_root,
                weight_root=weight_root,
                frequency_root=frequency_root,
                output_root=output_root,
                reference_family=reference_family,
                weight_family=weight_family,
                annotation_family=annotation_family,
                reference_columns=reference_columns,
                mode=mode,
                build_hint=build_hint,
                build_check=build_check,
                overwrite=overwrite,
                issues=issues,
            )
        if len(reference_columns) != 1:
            raise LDSCInputError(
                "Unpartitioned legacy reference suite must have exactly one scientific LD-score column when no "
                f"complete annotation family is present; found {reference_columns}."
            )
        profile = "unpartitioned"
        frames: list[pd.DataFrame] = []
        common_total = 0.0
        all_total = 0.0
        all_available = True
        selected_files: list[Path] = []
        intersections: dict[str, dict[str, int]] = {}
        count_origins: dict[str, dict[str, str]] = {}
        seen_reference: set[str] = set()
        seen_weight: set[str] = set()
        for chrom in _CHROMS:
            reference_path = reference_family.paths[chrom]
            weight_path = weight_family.paths[chrom]
            reference = _read_ldscore_table(reference_path, "reference", chrom)
            weight = _read_ldscore_table(weight_path, "weight", chrom)
            _require_columns(reference, {"CHR", "SNP", "BP"}, "reference", reference_path, chrom)
            _require_columns(weight, {"CHR", "SNP", "BP"}, "weight", weight_path, chrom)
            ref_science = _scientific_columns(reference)
            weight_science = _scientific_columns(weight)
            if ref_science != reference_columns:
                raise LDSCInputError(
                    f"Reference LD-score scientific columns differ on chromosome {chrom}: "
                    f"expected {reference_columns}, observed {ref_science} in '{reference_path}'."
                )
            if len(weight_science) != 1:
                raise LDSCInputError(
                    f"Weight LD-score table '{weight_path}' on chromosome {chrom} must contain exactly one "
                    f"scientific value column; found {weight_science}."
                )
            _validate_ldscore_rows(reference, ref_science, "reference", reference_path, chrom)
            _validate_ldscore_rows(weight, weight_science, "weight", weight_path, chrom)
            _validate_global_source_ids(reference["SNP"], seen_reference, "reference", reference_path, chrom)
            _validate_global_source_ids(weight["SNP"], seen_weight, "weight", weight_path, chrom)
            merged = _join_reference_weight(reference, weight, ref_science, weight_science[0], chrom, issues)
            frames.append(merged.rename(columns={ref_science[0]: "base"}))
            intersections[str(chrom)] = {
                "reference_rows": len(reference),
                "weight_rows": len(weight),
                "retained_rows": len(merged),
            }
            common_path = reference_root / f"{reference_prefix}{chrom}.l2.M_5_50"
            if not common_path.exists():
                raise LDSCInputError(
                    f"Required .l2.M_5_50 count file is missing for chromosome {chrom}: '{common_path}'. "
                    + _COUNT_FILENAME_GUIDANCE
                )
            common_total += _read_count_vector(common_path, 1, chrom, "M_5_50")[0]
            all_path = reference_root / f"{reference_prefix}{chrom}.l2.M"
            if all_path.exists():
                all_total += _read_count_vector(all_path, 1, chrom, "M")[0]
                all_origin = "legacy_unvalidated"
                selected_files.append(all_path)
            else:
                all_available = False
                all_origin = "missing"
            count_origins[str(chrom)] = {"base_all": all_origin, "base_common": "legacy_unvalidated"}
            selected_files.extend((reference_path, weight_path, common_path))

        baseline = pd.concat(frames, ignore_index=True)
        frames.clear()
        if baseline.empty:
            raise LDSCInputError(
                "Legacy reference and weight suites have an empty genome-wide rsID intersection. "
                "Check --legacy-reference-dir and --legacy-weight-dir: use matching suites from the same "
                "reference SNP population and confirm their SNP columns use the same rsID identifiers."
            )
        _validate_global_output_identity(baseline, mode)
        baseline, effective_build, coordinate_provenance = _resolve_output_coordinates(
            baseline, mode, build_hint, context="legacy reference LD-score suite", build_check=build_check
        )
        baseline = _sort_baseline(baseline)
        count_records = [
            {
                "group": "baseline",
                "column": "base",
                "all_reference_snp_count": float(all_total) if all_available else None,
                "common_reference_snp_count": float(common_total),
            }
        ]
        provenance = _build_provenance(
            profile=profile,
            reference_root=reference_root,
            weight_root=weight_root,
            frequency_root=None,
            reference_family=reference_family,
            weight_family=weight_family,
            frequency_prefix=None,
            selected_files=selected_files,
            intersections=intersections,
            count_origins=count_origins,
            coordinate_provenance=coordinate_provenance,
            issues=issues,
        )
        ldscore = LDScoreResult(
            baseline_table=baseline,
            query_table=None,
            count_records=count_records,
            baseline_columns=["base"],
            query_columns=[],
            ld_reference_snps=frozenset(),
            ld_regression_snps=frozenset(),  # not written; the published directory is reloaded for use
            chromosome_results=[],
            count_config=dict(_COUNT_CONFIG),
            config_snapshot=GlobalConfig(snp_identifier=mode, genome_build=effective_build),
            legacy_ldscore_conversion=provenance,
        )
        paths = LDScoreDirectoryWriter().write(
            ldscore,
            LDScoreOutputConfig(output_dir=output_root, overwrite=overwrite),
        )
        return LegacyLDScoreConversionResult(profile, len(baseline), ("base",), paths)

    def _convert_baseline_partitioned(
        self,
        *,
        reference_root: Path,
        weight_root: Path,
        frequency_root: Path,
        output_root: Path,
        reference_family: "_DiscoveredFamily",
        weight_family: "_DiscoveredFamily",
        annotation_family: dict[int, Path],
        reference_columns: list[str],
        mode: str,
        build_hint: str | None,
        build_check,
        overwrite: bool,
        issues: list[dict[str, object]],
    ) -> LegacyLDScoreConversionResult:
        """Convert a self-contained full baseline suite chromosome by chromosome."""
        frequency_family = _discover_frequency_family(frequency_root, issues)
        annotation_columns: list[str] | None = None
        column_mapping: dict[str, str] | None = None
        frames: list[pd.DataFrame] = []
        all_counts = np.zeros(len(reference_columns), dtype=float)
        common_counts = np.zeros(len(reference_columns), dtype=float)
        overlap_all = np.zeros((len(reference_columns), len(reference_columns)), dtype=float)
        overlap_common = np.zeros_like(overlap_all)
        n_all = 0
        n_common = 0
        selected_files: list[Path] = []
        intersections: dict[str, dict[str, int]] = {}
        count_origins: dict[str, dict[str, str]] = {}
        seen_reference: set[str] = set()
        seen_weight: set[str] = set()
        seen_annotations: set[str] = set()
        seen_frequency: set[str] = set()
        for chrom in _CHROMS:
            reference_path = reference_family.paths[chrom]
            weight_path = weight_family.paths[chrom]
            annotation_path = annotation_family[chrom]
            frequency_path = frequency_family.paths[chrom]
            reference = _read_ldscore_table(reference_path, "reference", chrom)
            weight = _read_ldscore_table(weight_path, "weight", chrom)
            annotation = _read_whitespace_table(annotation_path, "annotation", chrom)
            frequency = _read_whitespace_table(frequency_path, "frequency", chrom)
            _require_columns(reference, {"CHR", "SNP", "BP"}, "reference", reference_path, chrom)
            _require_columns(weight, {"CHR", "SNP", "BP"}, "weight", weight_path, chrom)
            _require_columns(annotation, {"CHR", "SNP", "BP", "CM"}, "annotation", annotation_path, chrom)
            ref_science = _scientific_columns(reference)
            weight_science = _scientific_columns(weight)
            annot_science = _scientific_columns(annotation)
            if ref_science != reference_columns:
                raise LDSCInputError(
                    f"Reference LD-score scientific columns differ on chromosome {chrom}: expected "
                    f"{reference_columns}, observed {ref_science} in '{reference_path}'."
                )
            if len(weight_science) != 1:
                raise LDSCInputError(
                    f"Weight LD-score table '{weight_path}' on chromosome {chrom} must contain exactly one "
                    f"scientific value column; found {weight_science}."
                )
            if annotation_columns is None:
                annotation_columns = annot_science
                column_mapping = _resolve_annotation_mapping(reference_columns, annotation_columns)
            elif annot_science != annotation_columns:
                raise LDSCInputError(
                    f"Baseline annotation columns differ on chromosome {chrom}: expected {annotation_columns}, "
                    f"observed {annot_science} in '{annotation_path}'."
                )
            assert column_mapping is not None
            _validate_ldscore_rows(reference, ref_science, "reference", reference_path, chrom)
            _validate_ldscore_rows(weight, weight_science, "weight", weight_path, chrom)
            _validate_annotation_rows(annotation, annotation_columns, annotation_path, chrom)
            _validate_global_source_ids(reference["SNP"], seen_reference, "reference", reference_path, chrom)
            _validate_global_source_ids(weight["SNP"], seen_weight, "weight", weight_path, chrom)
            _validate_global_source_ids(annotation["SNP"], seen_annotations, "annotation", annotation_path, chrom)
            _audit_coordinate_disagreements(reference, annotation, chrom, annotation_path, issues)

            freq_column = _frequency_column(frequency, frequency_path, chrom)
            _validate_frequency_rows(frequency, freq_column, frequency_path, chrom)
            _validate_global_source_ids(
                frequency["SNP"], seen_frequency, "frequency", frequency_path, chrom
            )
            frequency_payload = frequency.loc[:, ["SNP", freq_column]]
            # Only SNP and the frequency are read from the join, so join on the SNP column alone.
            annotation_with_frequency = pd.merge(
                annotation.loc[:, ["SNP"]],
                frequency_payload,
                how="left",
                on="SNP",
                sort=False,
                validate="one_to_one",
            )
            missing_frequency = annotation_with_frequency[freq_column].isna()
            if missing_frequency.any():
                sample = annotation_with_frequency.loc[missing_frequency, "SNP"].astype(str).head(10).tolist()
                for snp in annotation_with_frequency.loc[missing_frequency, "SNP"].astype(str):
                    issues.append(
                        _issue(
                            "error",
                            "frequency",
                            file=frequency_path,
                            chromosome=chrom,
                            snp=snp,
                            reason="missing_frequency",
                        )
                    )
                raise LDSCInputError(
                    f"Legacy frequency file '{frequency_path}' is missing {int(missing_frequency.sum())} "
                    f"baseline annotation SNP(s) on chromosome {chrom}; cause SNPs: {sample}."
                )
            annotation_ids = set(annotation["SNP"].astype(str))
            for snp in sorted(set(frequency["SNP"].astype(str)) - annotation_ids):
                issues.append(
                    _issue(
                        "warning",
                        "frequency",
                        file=frequency_path,
                        chromosome=chrom,
                        snp=snp,
                        reason="extra_frequency_rsid",
                    )
                )

            # Fill the float64 matrix column by column in the column-major layout that
            # DataFrame.to_numpy() returns, so sums and products see the same memory order.
            matrix = np.empty((len(annotation), len(reference_columns)), dtype=float, order="F")
            for index, column in enumerate(reference_columns):
                matrix[:, index] = annotation[column_mapping[column]].to_numpy(dtype=float)
            frequency_values = annotation_with_frequency[freq_column].to_numpy(dtype=float)
            if freq_column.upper() == "MAF":
                common_mask = frequency_values > 0.05
            else:
                common_mask = (frequency_values > 0.05) & (frequency_values < 0.95)
            common_matrix = matrix[common_mask]
            reconstructed_all = matrix.sum(axis=0)
            reconstructed_common = common_matrix.sum(axis=0)
            overlap_all += matrix.T @ matrix
            overlap_common += common_matrix.T @ common_matrix
            n_all += len(matrix)
            n_common += int(common_mask.sum())

            common_path = reference_root / f"{reference_family.prefix}{chrom}.l2.M_5_50"
            if not common_path.exists():
                raise LDSCInputError(
                    f"Required .l2.M_5_50 count file is missing for chromosome {chrom}: '{common_path}'. "
                    + _COUNT_FILENAME_GUIDANCE
                )
            legacy_common = _read_count_vector(common_path, len(reference_columns), chrom, "M_5_50")
            _validate_reconstructed_counts(
                legacy_common,
                reconstructed_common,
                matrix,
                reference_columns,
                chrom,
                "M_5_50",
                common_path,
                issues,
            )
            common_counts += legacy_common
            all_path = reference_root / f"{reference_family.prefix}{chrom}.l2.M"
            if all_path.exists():
                legacy_all = _read_count_vector(all_path, len(reference_columns), chrom, "M")
                _validate_reconstructed_counts(
                    legacy_all,
                    reconstructed_all,
                    matrix,
                    reference_columns,
                    chrom,
                    "M",
                    all_path,
                    issues,
                )
                all_counts += legacy_all
                all_origins = {column: "legacy_validated" for column in reference_columns}
                selected_files.append(all_path)
            else:
                all_counts += reconstructed_all
                all_origins = {column: "reconstructed" for column in reference_columns}
            count_origins[str(chrom)] = {
                **{f"{column}_all": origin for column, origin in all_origins.items()},
                **{f"{column}_common": "legacy_validated" for column in reference_columns},
            }

            merged = _join_reference_weight(reference, weight, reference_columns, weight_science[0], chrom, issues)
            frames.append(merged)
            intersections[str(chrom)] = {
                "reference_rows": len(reference),
                "weight_rows": len(weight),
                "retained_rows": len(merged),
            }
            selected_files.extend(
                (reference_path, weight_path, annotation_path, frequency_path, common_path)
            )

        baseline = pd.concat(frames, ignore_index=True)
        frames.clear()
        if baseline.empty:
            raise LDSCInputError(
                "Legacy reference and weight suites have an empty genome-wide rsID intersection. "
                "Check --legacy-reference-dir and --legacy-weight-dir: use matching suites from the same "
                "reference SNP population and confirm their SNP columns use the same rsID identifiers."
            )
        _validate_global_output_identity(baseline, mode)
        baseline, effective_build, coordinate_provenance = _resolve_output_coordinates(
            baseline, mode, build_hint, context="legacy baseline reference LD-score suite", build_check=build_check
        )
        baseline = _sort_baseline(baseline)
        count_records = [
            {
                "group": "baseline",
                "column": column,
                "all_reference_snp_count": float(all_counts[index]),
                "common_reference_snp_count": float(common_counts[index]),
            }
            for index, column in enumerate(reference_columns)
        ]
        overlap = LDScoreOverlap(
            baseline_block_all=pd.DataFrame(overlap_all, index=reference_columns, columns=reference_columns),
            baseline_block_common=pd.DataFrame(
                overlap_common, index=reference_columns, columns=reference_columns
            ),
            query_diagonal_all=pd.Series(dtype=float),
            query_diagonal_common=pd.Series(dtype=float),
            total_all_reference_snps=float(n_all),
            total_common_reference_snps=float(n_common),
        )
        provenance = _build_provenance(
            profile="baseline_partitioned",
            reference_root=reference_root,
            weight_root=weight_root,
            frequency_root=frequency_root,
            reference_family=reference_family,
            weight_family=weight_family,
            frequency_prefix=frequency_family.prefix,
            selected_files=selected_files,
            intersections=intersections,
            count_origins=count_origins,
            coordinate_provenance=coordinate_provenance,
            issues=issues,
        )
        ldscore = LDScoreResult(
            baseline_table=baseline,
            query_table=None,
            count_records=count_records,
            baseline_columns=list(reference_columns),
            query_columns=[],
            ld_reference_snps=frozenset(seen_annotations),
            ld_regression_snps=frozenset(),  # not written; the published directory is reloaded for use
            chromosome_results=[],
            count_config=dict(_COUNT_CONFIG),
            config_snapshot=GlobalConfig(snp_identifier=mode, genome_build=effective_build),
            overlap=overlap,
            legacy_ldscore_conversion=provenance,
        )
        paths = LDScoreDirectoryWriter().write(
            ldscore,
            LDScoreOutputConfig(output_dir=output_root, overwrite=overwrite),
        )
        return LegacyLDScoreConversionResult(
            "baseline_partitioned", len(baseline), tuple(reference_columns), paths
        )


@dataclass(frozen=True)
class _DiscoveredFamily:
    prefix: str
    paths: dict[int, Path]


def convert_legacy_ldscores(
    *,
    legacy_reference_dir: str | Path,
    legacy_weight_dir: str | Path,
    output_dir: str | Path,
    legacy_frequency_dir: str | Path | None = None,
    snp_identifier: str = "rsid",
    genome_build: str | None = "auto",
    overwrite: bool = False,
    log_level: str | None = None,
) -> LegacyLDScoreConversionResult:
    """Convert one selected legacy LD-score suite to a canonical LDSC3 directory.

    Parameters
    ----------
    legacy_reference_dir, legacy_weight_dir : path-like
        Complete chromosomes 1-22 legacy reference and weight directories.
    output_dir : path-like
        Destination for the canonical LDSC3 LD-score directory.
    legacy_frequency_dir : path-like, optional
        Complete frequency family required by baseline partitioned conversion.
    snp_identifier : {"rsid", "chr_pos"}, optional
        Allele-unaware output identity. Default is ``"rsid"``.
    genome_build : {"auto", "hg19", "hg38"}, optional
        Reference coordinate interpretation. Default is ``"auto"``.
    overwrite : bool, optional
        Replace converter-owned destination files. Default is false.
    log_level : {"DEBUG", "INFO", "WARNING", "ERROR"}, optional
        Converter log threshold. Default is ``None``: use the registered
        ``GlobalConfig.log_level`` (``"INFO"`` unless changed).

    Returns
    -------
    LegacyLDScoreConversionResult
        Detected profile, retained rows, columns, and output paths.

    Notes
    -----
    This function mirrors :meth:`LegacyLDScoreConverter.convert` and
    intentionally has no configurable common-MAF threshold.
    """
    return LegacyLDScoreConverter().convert(
        legacy_reference_dir=legacy_reference_dir,
        legacy_weight_dir=legacy_weight_dir,
        output_dir=output_dir,
        legacy_frequency_dir=legacy_frequency_dir,
        snp_identifier=snp_identifier,
        genome_build=genome_build,
        overwrite=overwrite,
        log_level=log_level,
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the standalone converter parser used by the unified CLI.

    Returns
    -------
    argparse.ArgumentParser
        Parser containing the closed converter argument surface.
    """
    parser = argparse.ArgumentParser(prog="ldsc convert-legacy-ldscores", allow_abbrev=False)
    parser.prog = 'ldsc convert-legacy-ldscores'
    parser.formatter_class = CLIHelpFormatter
    parser.description = 'Convert complete legacy reference and regression-weight suites into one LDSC3 LD-score directory.'
    legacy = parser.add_argument_group('Step 1: Legacy LD-score files')
    identity = parser.add_argument_group('Step 2: Output SNP identity and genome build')
    output = parser.add_argument_group('Step 3: Output')
    logs = parser.add_argument_group('Logging')

    legacy.add_argument(
        '--legacy-reference-dir', required=True, metavar='DIR',
        help=(
            'Required directory containing a complete chromosome 1-22 legacy reference LD-score suite and '
            '.l2.M_5_50 common-SNP count files. Baseline partitioned suites also require full annotations and '
            '--legacy-frequency-dir.'
        ),
    )
    legacy.add_argument(
        '--legacy-weight-dir', required=True, metavar='DIR',
        help=(
            'Required directory containing chromosome 1-22 single-column regression-weight LD scores. May '
            'equal --legacy-reference-dir for an unpartitioned suite.'
        ),
    )
    legacy.add_argument(
        '--legacy-frequency-dir', metavar='DIR',
        help=(
            'Directory containing chromosome frequency files used to reconstruct common-SNP counts and '
            'overlap. Required for baseline partitioned --legacy-reference-dir input; optional for '
            'unpartitioned input, with no frequency suite read when omitted.'
        ),
    )

    identity.add_argument(
        '--snp-identifier', choices=('rsid', 'chr_pos'), default='rsid',
        help=(
            'Output SNP matching rule: rsid uses SNP identifiers; chr_pos uses chromosome and position. '
            'Default: rsid. Both ignore alleles; chr_pos requires a resolved genome build.'
        ),
    )
    identity.add_argument(
        '--genome-build', choices=('auto', 'hg19', 'hg38'), default='auto',
        help=(
            'Genome build of legacy reference coordinates; no liftover is performed. Default: '
            'auto, check every retained reference chromosome for chr_pos identity. A concrete build is '
            'required to agree with input declarations; rsid identity needs no coordinate-build inference.'
        ),
    )

    output.add_argument(
        '--output-dir', required=True, metavar='DIR',
        help=(
            'Required destination for the converted LD-score tables, counts, and diagnostics. Source files '
            'are not changed.'
        ),
    )
    output.add_argument(
        '--overwrite', action='store_true',
        help=(
            "Replace this command's existing output files and remove obsolete outputs from an earlier run. "
            'Default: off; stop if output files already exist.'
        ),
    )

    logs.add_argument(
        '--log-level', choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'), default='INFO',
        help=LOG_LEVEL_HELP,
    )
    return parser


def main(argv: list[str] | None = None) -> LegacyLDScoreConversionResult:
    """Parse converter arguments and run one conversion.

    Parameters
    ----------
    argv : list of str or None, optional
        Explicit arguments, or ``None`` to read ``sys.argv``.

    Returns
    -------
    LegacyLDScoreConversionResult
        Completed conversion summary.
    """
    args = build_parser().parse_args(argv)
    return convert_legacy_ldscores(
        legacy_reference_dir=args.legacy_reference_dir,
        legacy_weight_dir=args.legacy_weight_dir,
        legacy_frequency_dir=args.legacy_frequency_dir,
        output_dir=args.output_dir,
        snp_identifier=args.snp_identifier,
        genome_build=args.genome_build,
        overwrite=args.overwrite,
        log_level=args.log_level,
    )


def _existing_directory(path: str | Path, label: str) -> Path:
    root = Path(normalize_path_token(path))
    if not root.is_dir():
        raise LDSCInputError(
            f"Cannot use {label} '{root}': path is not an existing directory. "
            "Pass the directory containing the chromosome files, not a file, filename prefix, or glob."
        )
    return root


def _discover_ldscore_family(root: Path, role: str, issues: list[dict[str, object]]) -> _DiscoveredFamily:
    groups: dict[str, dict[int, list[Path]]] = {}
    unsupported_weight_paths: list[Path] = []
    for path in sorted(root.iterdir(), key=lambda item: item.name):
        match = _LDSCORE_RE.match(path.name)
        if match is None:
            if path.is_file() and ".l2.ldscore" in path.name:
                weight_match = _LDSCORE_RE.match(path.name.replace(".w.l2.ldscore", ".l2.ldscore"))
                if weight_match is not None:
                    unsupported_weight_paths.append(path)
                    issues.append(_issue(
                        "warning", role, file=path, chromosome=int(weight_match.group("chrom")),
                        reason="discarded_unsupported_weight_filename", details=_WEIGHT_FILENAME_GUIDANCE,
                    ))
                else:
                    issues.append(_issue(
                        "warning", role, file=path, reason="discarded_non_family_file",
                        details=_LDSCORE_FILENAME_GUIDANCE,
                    ))
            continue
        groups.setdefault(match.group("prefix"), {}).setdefault(int(match.group("chrom")), []).append(path)
    coherent = {prefix: shards for prefix, shards in groups.items() if set(shards) == set(_CHROMS)}
    if len(coherent) != 1:
        raise LDSCInputError(
            f"Legacy {role} directory '{root}' must contain exactly one coherent chromosome 1-22 "
            f".l2.ldscore family; found prefixes {sorted(coherent)}."
            + " " + (_WEIGHT_FILENAME_GUIDANCE if unsupported_weight_paths else _LDSCORE_FILENAME_GUIDANCE)
        )
    if unsupported_weight_paths:
        LOGGER.warning(
            "Ignoring %d unsupported weight filename(s) in '%s'; using the complete recognized family. %s",
            len(unsupported_weight_paths), root, _WEIGHT_FILENAME_GUIDANCE,
        )
    prefix, shards = next(iter(coherent.items()))
    for discarded_prefix, discarded_shards in groups.items():
        if discarded_prefix == prefix:
            continue
        for representations in discarded_shards.values():
            for path in representations:
                issues.append(
                    _issue("warning", role, file=path, reason="discarded_nonselected_family_file")
                )
    selected: dict[int, Path] = {}
    for chrom, representations in shards.items():
        if len(representations) == 1:
            selected[chrom] = representations[0]
            continue
        if len(representations) != 2 or {path.suffix for path in representations} != {".gz", ".ldscore"}:
            raise LDSCInputError(
                f"Legacy {role} family has ambiguous representations for chromosome {chrom}: "
                f"{[str(path) for path in representations]}. "
                "Keep one correctly named representation per chromosome from the intended source release."
            )
        contents = [_decompressed_bytes(path) for path in representations]
        if contents[0] != contents[1]:
            raise LDSCInputError(
                f"Legacy {role} family has conflicting plain/gzip representations for chromosome {chrom}: "
                f"{[str(path) for path in representations]}. "
                "Select the correct source release and organize one representation per chromosome in a separate "
                "directory; retaining both requires identical decompressed contents."
            )
        selected[chrom] = next(path for path in representations if path.suffix == ".gz")
        for path in representations:
            if path != selected[chrom]:
                issues.append(
                    _issue("warning", role, file=path, chromosome=chrom, reason="discarded_identical_plain_duplicate")
                )
    return _DiscoveredFamily(prefix=prefix, paths=selected)


def _discover_associated_family(
    root: Path,
    prefix: str,
    *,
    suffixes: tuple[str, ...],
    issues: list[dict[str, object]],
) -> dict[int, Path]:
    found: dict[int, Path] = {}
    for chrom in _CHROMS:
        matches = [root / f"{prefix}{chrom}{suffix}" for suffix in suffixes]
        existing = [path for path in matches if path.exists()]
        if len(existing) == 2:
            if _decompressed_bytes(existing[0]) != _decompressed_bytes(existing[1]):
                raise LDSCInputError(
                    f"Legacy suite has conflicting plain/gzip representations for chromosome {chrom}: {existing}. "
                    "Use one annotation copy matching the reference release and prefix; retaining both "
                    "requires identical decompressed contents."
                )
            found[chrom] = next(path for path in existing if path.suffix == ".gz")
            for path in existing:
                if path != found[chrom]:
                    issues.append(
                        _issue(
                            "warning",
                            "annotation",
                            file=path,
                            chromosome=chrom,
                            reason="discarded_identical_plain_duplicate",
                        )
                    )
        elif existing:
            found[chrom] = existing[0]
    return found


def _discover_frequency_family(
    root: Path, issues: list[dict[str, object]]
) -> _DiscoveredFamily:
    groups: dict[str, dict[int, list[Path]]] = {}
    for path in sorted(root.iterdir(), key=lambda item: item.name):
        match = _FREQUENCY_RE.match(path.name)
        if match is None:
            continue
        groups.setdefault(match.group("prefix"), {}).setdefault(int(match.group("chrom")), []).append(path)
    coherent = {prefix: shards for prefix, shards in groups.items() if set(shards) == set(_CHROMS)}
    if len(coherent) != 1:
        raise LDSCInputError(
            f"Legacy frequency directory '{root}' must contain exactly one coherent chromosome 1-22 .frq family; "
            f"found prefixes {sorted(coherent)}. Use '<prefix><chrom>.frq(.gz)', e.g. '1000G.EUR.QC.1.frq.gz', "
            "with one constant prefix (possibly empty) and unpadded chromosomes 1-22. "
            "Restore missing shards or rename misnamed copies; place one complete family directly in this "
            "directory and keep different releases in separate directories. "
            "See docs/user-guide/troubleshooting.md#convert-legacy-ldscores."
        )
    prefix, shards = next(iter(coherent.items()))
    for discarded_prefix, discarded_shards in groups.items():
        if discarded_prefix == prefix:
            continue
        for representations in discarded_shards.values():
            for path in representations:
                issues.append(
                    _issue("warning", "frequency", file=path, reason="discarded_nonselected_family_file")
                )
    selected: dict[int, Path] = {}
    for chrom, representations in shards.items():
        if len(representations) == 1:
            selected[chrom] = representations[0]
            continue
        if len(representations) == 2 and _decompressed_bytes(representations[0]) == _decompressed_bytes(
            representations[1]
        ):
            selected[chrom] = next(path for path in representations if path.suffix == ".gz")
            for path in representations:
                if path != selected[chrom]:
                    issues.append(
                        _issue(
                            "warning",
                            "frequency",
                            file=path,
                            chromosome=chrom,
                            reason="discarded_identical_plain_duplicate",
                        )
                    )
            continue
        raise LDSCInputError(
            f"Legacy frequency family has ambiguous or conflicting representations for chromosome {chrom}: "
            f"{[str(path) for path in representations]}. "
            "Keep one '<prefix><chrom>.frq(.gz)' copy from the matching source release per chromosome; "
            "retaining plain and gzip copies requires identical decompressed contents."
        )
    return _DiscoveredFamily(prefix, selected)


def _read_ldscore_table(path: Path, role: str, chrom: int) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep=r"\s+", compression="infer", comment="#")
    except Exception as exc:
        raise LDSCInputError(
            f"Could not parse legacy {role} LD-score table '{path}' for chromosome {chrom}: {exc}. "
            "Check the named file's readability, compression integrity, and whitespace-delimited header; "
            "restore the original table from the selected legacy suite. "
            "See docs/user-guide/troubleshooting.md#conversion-rejected-the-legacy-suite."
        ) from exc


def _read_whitespace_table(path: Path, role: str, chrom: int) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep=r"\s+", compression="infer", comment="#")
    except Exception as exc:
        raise LDSCInputError(
            f"Could not parse legacy {role} table '{path}' for chromosome {chrom}: {exc}. "
            "Check the named file's readability, compression integrity, and whitespace-delimited header; "
            "restore the original table from the selected legacy suite. "
            "See docs/user-guide/troubleshooting.md#conversion-rejected-the-legacy-suite."
        ) from exc


def _scientific_columns(frame: pd.DataFrame) -> list[str]:
    return [column for column in frame.columns if str(column).upper() not in _METADATA_COLUMNS]


def _require_columns(frame: pd.DataFrame, required: set[str], role: str, path: Path, chrom: int) -> None:
    missing = sorted(required - set(frame.columns))
    if missing:
        raise LDSCInputError(
            f"Legacy {role} table '{path}' on chromosome {chrom} is missing required columns {missing}. "
            "Restore the original table with its header from the selected legacy suite, or export it again "
            "from the source dataset; do not add placeholder scientific columns."
        )


def _validate_ldscore_rows(
    frame: pd.DataFrame,
    scientific: list[str],
    role: str,
    path: Path,
    chrom: int,
) -> None:
    duplicate = frame["SNP"].astype("string").duplicated(keep=False)
    if duplicate.any():
        sample = frame.loc[duplicate, "SNP"].astype(str).head(10).tolist()
        raise LDSCInputError(
            f"Legacy {role} table '{path}' on chromosome {chrom} has duplicate rsIDs; example SNPs: {sample}. "
            "Use the original table with unique rsIDs from the selected legacy suite; check whether rows "
            "from different releases or chromosomes were combined."
        )
    for column in scientific:
        values = pd.to_numeric(frame[column], errors="coerce")
        bad = ~np.isfinite(values.to_numpy(dtype=float))
        if bad.any():
            sample = frame.loc[bad, "SNP"].astype(str).head(10).tolist()
            raise LDSCInputError(
                f"Legacy {role} table '{path}' on chromosome {chrom} has non-finite values in '{column}'; "
                f"example SNPs: {sample}. Restore the original finite LD-score values from the selected "
                "legacy suite, or regenerate the source scores; do not replace missing values with zero."
            )
        frame[column] = values.astype(float)
    frame["BP"] = pd.to_numeric(frame["BP"], errors="coerce")
    bad_bp = frame["BP"].isna() | (frame["BP"] < 0)
    if bad_bp.any():
        sample = frame.loc[bad_bp, "SNP"].astype(str).head(10).tolist()
        raise LDSCInputError(
            f"Legacy {role} table '{path}' on chromosome {chrom} has invalid BP values; example SNPs: {sample}. "
            "Restore numeric, nonnegative BP coordinates from the original matching reference table; "
            "do not infer a different genome build by renaming files."
        )


def _resolve_annotation_mapping(reference_columns: list[str], annotation_columns: list[str]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for column in reference_columns:
        candidates: list[str] = []
        if column in annotation_columns:
            candidates.append(column)
        if column.endswith("L2") and column[:-2] in annotation_columns:
            candidates.append(column[:-2])
        candidates = list(dict.fromkeys(candidates))
        if len(candidates) != 1:
            raise LDSCInputError(
                f"Cannot match legacy LD-score column '{column}' to exactly one full annotation column; "
                f"candidates are {candidates}, annotation columns are {annotation_columns}. "
                "Supply the full annotation tables used to calculate this LD-score suite, "
                "with their original corresponding headers."
            )
        mapping[column] = candidates[0]
    if len(set(mapping.values())) != len(annotation_columns) or set(mapping.values()) != set(annotation_columns):
        raise LDSCInputError(
            "Legacy baseline LD-score and annotation columns do not form a complete one-to-one mapping: "
            f"LD-score columns {reference_columns}, annotation columns {annotation_columns}, mapping {mapping}. "
            "Supply full LD-score and annotation suites from the same reference release with corresponding columns."
        )
    return mapping


def _validate_annotation_rows(frame: pd.DataFrame, columns: list[str], path: Path, chrom: int) -> None:
    duplicate = frame["SNP"].astype("string").duplicated(keep=False)
    if duplicate.any():
        sample = frame.loc[duplicate, "SNP"].astype(str).head(10).tolist()
        raise LDSCInputError(
            f"Legacy annotation table '{path}' on chromosome {chrom} has duplicate rsIDs; example SNPs: {sample}. "
            "Restore the original annotation table with unique rsIDs that matches the selected reference LD scores."
        )
    for column in columns:
        values = pd.to_numeric(frame[column], errors="coerce")
        bad = ~np.isfinite(values.to_numpy(dtype=float))
        if bad.any():
            sample = frame.loc[bad, "SNP"].astype(str).head(10).tolist()
            raise LDSCInputError(
                f"Legacy annotation table '{path}' on chromosome {chrom} has non-finite values in '{column}'; "
                f"example SNPs: {sample}. Restore the original finite annotation values used for the selected "
                "LD scores; do not replace unknown annotation values with zero."
            )
        frame[column] = values.astype(float)


def _validate_global_source_ids(
    values: pd.Series,
    seen: set[str],
    role: str,
    path: Path,
    chrom: int,
) -> None:
    current = set(values.astype(str))
    repeated = sorted(current & seen)
    if repeated:
        raise LDSCInputError(
            f"Legacy {role} family repeats rsIDs across chromosome shards at '{path}' chromosome {chrom}; "
            f"example SNPs: {repeated[:10]}. Restore a coherent chromosome suite from one reference "
            "release and check whether a chromosome file was copied under the wrong name."
        )
    seen.update(current)


def _audit_coordinate_disagreements(
    reference: pd.DataFrame,
    annotation: pd.DataFrame,
    chrom: int,
    annotation_path: Path,
    issues: list[dict[str, object]],
) -> None:
    compared = pd.merge(
        reference.loc[:, ["SNP", "CHR", "BP"]],
        annotation.loc[:, ["SNP", "CHR", "BP"]].rename(columns={"CHR": "annotation_CHR", "BP": "annotation_BP"}),
        how="inner",
        on="SNP",
        sort=False,
    )
    mismatch = (compared["CHR"].astype(str) != compared["annotation_CHR"].astype(str)) | (
        pd.to_numeric(compared["BP"]) != pd.to_numeric(compared["annotation_BP"])
    )
    for _, row in compared.loc[mismatch].iterrows():
        issues.append(
            _issue(
                "warning",
                "annotation",
                file=annotation_path,
                chromosome=chrom,
                snp=row["SNP"],
                reason="coordinate_disagreement",
                observed=f"{row['annotation_CHR']}:{row['annotation_BP']}",
                expected=f"{row['CHR']}:{row['BP']}",
            )
        )


def _frequency_column(frame: pd.DataFrame, path: Path, chrom: int) -> str:
    candidates = [column for column in frame.columns if str(column).upper() in {"MAF", "FRQ"}]
    if len(candidates) != 1:
        raise LDSCInputError(
            f"Legacy frequency table '{path}' on chromosome {chrom} must contain exactly one MAF or FRQ column; "
            f"found {candidates}. Supply the original frequency table for this reference chromosome "
            "with one clearly identified frequency column; do not invent frequencies."
        )
    if "SNP" not in frame.columns:
        raise LDSCInputError(
            f"Legacy frequency table '{path}' on chromosome {chrom} is missing required SNP column. "
            "Supply the original frequency table with rsID identifiers matching the selected reference LD scores."
        )
    return candidates[0]


def _validate_frequency_rows(frame: pd.DataFrame, column: str, path: Path, chrom: int) -> None:
    duplicate = frame["SNP"].astype("string").duplicated(keep=False)
    if duplicate.any():
        sample = frame.loc[duplicate, "SNP"].astype(str).head(10).tolist()
        raise LDSCInputError(
            f"Legacy frequency table '{path}' on chromosome {chrom} has duplicate rsIDs; example SNPs: {sample}. "
            "Restore the original per-chromosome frequency table with unique rsIDs from this reference release."
        )
    values = pd.to_numeric(frame[column], errors="coerce")
    upper = 0.5 if column.upper() == "MAF" else 1.0
    bad = ~np.isfinite(values.to_numpy(dtype=float)) | (values < 0.0) | (values > upper)
    if bad.any():
        sample = frame.loc[bad, "SNP"].astype(str).head(10).tolist()
        raise LDSCInputError(
            f"Legacy frequency table '{path}' on chromosome {chrom} has invalid {column} values; "
            f"allowed range is [0, {upper}], cause SNPs: {sample}."
        )
    frame[column] = values.astype(float)


def _validate_reconstructed_counts(
    legacy: np.ndarray,
    reconstructed: np.ndarray,
    annotation_matrix: np.ndarray,
    columns: list[str],
    chrom: int,
    kind: str,
    path: Path,
    issues: list[dict[str, object]],
) -> None:
    mismatches: list[str] = []
    for index, column in enumerate(columns):
        values = annotation_matrix[:, index]
        integer_annotation = bool(np.all(values == np.round(values)))
        agrees = (
            legacy[index] == reconstructed[index]
            if integer_annotation
            else np.isclose(legacy[index], reconstructed[index], rtol=1e-8, atol=1e-6)
        )
        if agrees:
            continue
        mismatches.append(column)
        issues.append(
            _issue(
                "error",
                "count",
                file=path,
                chromosome=chrom,
                annotation=column,
                reason="count_mismatch",
                observed=float(legacy[index]),
                expected=float(reconstructed[index]),
                details=f"{kind} disagrees with annotation/frequency reconstruction",
            )
        )
    if mismatches:
        first = mismatches[0]
        index = columns.index(first)
        raise LDSCInputError(
            f"Legacy {kind} count conflict on chromosome {chrom} for annotation {first}: source "
            f"{legacy[index]} vs reconstructed {reconstructed[index]}. The suite may mix releases or filters."
        )


def _join_reference_weight(
    reference: pd.DataFrame,
    weight: pd.DataFrame,
    reference_columns: list[str],
    weight_column: str,
    chrom: int,
    issues: list[dict[str, object]],
) -> pd.DataFrame:
    ref = reference.loc[:, ["CHR", "SNP", "BP", *reference_columns]]
    wgt = weight.loc[:, ["CHR", "SNP", "BP", weight_column]].rename(
        columns={"CHR": "weight_CHR", "BP": "weight_BP", weight_column: "regression_ld_scores"}, copy=False
    )
    ref_ids = set(ref["SNP"].astype(str))
    weight_ids = set(wgt["SNP"].astype(str))
    for snp in sorted(ref_ids - weight_ids):
        issues.append(_issue("warning", "reference", chromosome=chrom, snp=snp, reason="reference_only_rsid"))
    for snp in sorted(weight_ids - ref_ids):
        issues.append(_issue("warning", "weight", chromosome=chrom, snp=snp, reason="weight_only_rsid"))
    merged = pd.merge(ref, wgt, how="inner", on="SNP", sort=False)
    disagreement = (merged["CHR"].astype(str) != merged["weight_CHR"].astype(str)) | (
        pd.to_numeric(merged["BP"]) != pd.to_numeric(merged["weight_BP"])
    )
    for _, row in merged.loc[disagreement].iterrows():
        issues.append(
            _issue(
                "warning",
                "weight",
                chromosome=chrom,
                snp=row["SNP"],
                reason="coordinate_disagreement",
                observed=f"{row['weight_CHR']}:{row['weight_BP']}",
                expected=f"{row['CHR']}:{row['BP']}",
            )
        )
    return merged.loc[:, ["CHR", "BP", "SNP", "regression_ld_scores", *reference_columns]].rename(
        columns={"BP": "POS"}, copy=False
    )


def _read_count_vector(path: Path, expected: int, chrom: int, kind: str) -> np.ndarray:
    try:
        values = np.asarray(path.read_text(encoding="utf-8").split(), dtype=float)
    except Exception as exc:
        raise LDSCInputError(f"Could not parse legacy {kind} count file '{path}' on chromosome {chrom}: {exc}") from exc
    if values.size != expected or not np.isfinite(values).all() or (values < 0).any():
        raise LDSCInputError(
            f"Legacy {kind} count file '{path}' on chromosome {chrom} must contain {expected} finite, "
            f"nonnegative value(s); observed {values.tolist()}."
        )
    return values


def _validate_global_output_identity(frame: pd.DataFrame, mode: str) -> None:
    duplicate_rsid = frame["SNP"].astype("string").duplicated(keep=False)
    if duplicate_rsid.any():
        sample = frame.loc[duplicate_rsid, "SNP"].astype(str).head(10).tolist()
        raise LDSCInputError(
            f"Converted suite has duplicate genome-wide rsIDs; cause SNPs: {sample}."
            " Supply a coherent complete legacy suite without repeated SNP identities; restore original "
            "matching shards rather than mixing source releases."
        )
    if mode == "chr_pos":
        keys = frame["CHR"].astype(str) + ":" + frame["POS"].astype(str)
        duplicate_pos = keys.duplicated(keep=False)
        if duplicate_pos.any():
            sample = frame.loc[duplicate_pos, "SNP"].astype(str).head(10).tolist()
            raise LDSCInputError(
                f"Converted chr_pos suite has duplicate coordinates; cause SNPs: {sample}."
                " Supply a coherent complete legacy suite without repeated SNP identities; restore original "
                "matching shards rather than mixing source releases."
            )


def _resolve_output_coordinates(
    frame: pd.DataFrame,
    mode: str,
    build_hint: str | None,
    *,
    context: str,
    build_check=None,
) -> tuple[pd.DataFrame, str | None, dict[str, object]]:
    """Resolve every retained reference chromosome, or honor a concrete build.

    Auto coordinate inference may normalize zero-based positions. A concrete
    declaration assumes one-based coordinates and performs no conversion.
    rsID output carries no coordinate build and requires no inference.
    """
    output = frame.copy(deep=False)  # the 0-based branch below replaces POS in this frame only
    from .._genome.build_policy import GenomeBuildCheck
    check = build_check or GenomeBuildCheck(build_hint, context='convert-legacy-ldscores', logger=LOGGER)
    if build_check is None and mode != "rsid":
        check.frames('--legacy-reference-dir', context, [output])
    if build_hint in {"hg19", "hg38"}:
        check.resolve()
        effective_build = None if mode == "rsid" else build_hint
        return output, effective_build, {
            "requested_build": build_hint, "effective_build": effective_build,
            "inferred_build": None, "coordinate_basis": "declared_1-based", "inference_error": None,
        }
    if mode == "rsid":
        # Reference coordinates are descriptive for rsID joins.
        return output, None, {
            "requested_build": build_hint, "effective_build": None,
            "inferred_build": None, "coordinate_basis": "unknown", "inference_error": None,
        }
    bases = {inference.coordinate_basis for _, inference in check.inferences}
    if len(bases) > 1:
        check.unresolved('--legacy-reference-dir', context, 'inconsistent coordinate basis; use one-based positions')
    effective_build = check.resolve()
    basis = next(iter(bases)) if bases else "declared_1-based"
    if basis == "0-based":
        output["POS"] = pd.to_numeric(output["POS"]).astype("int64") + 1
    return output, effective_build, {
        "requested_build": build_hint, "effective_build": effective_build,
        "inferred_build": effective_build if check.inferences else None, "coordinate_basis": basis, "inference_error": None,
    }


def _sort_baseline(frame: pd.DataFrame) -> pd.DataFrame:
    """Return rows ordered by numeric chromosome, POS and SNP, with string CHR and int64 POS.

    Only the three sort keys are sorted; the full table is copied once by the
    resulting row permutation.
    """
    chrom = frame["CHR"].astype(str)
    pos = pd.to_numeric(frame["POS"], errors="raise").astype("int64")
    keys = pd.DataFrame({"_chrom_order": pd.to_numeric(chrom, errors="coerce").to_numpy(), "POS": pos.to_numpy(),
                         "SNP": frame["SNP"].to_numpy()})
    order = keys.sort_values(["_chrom_order", "POS", "SNP"], kind="mergesort").index.to_numpy()
    ordered = frame.take(order)
    ordered["CHR"] = chrom.to_numpy()[order]
    ordered["POS"] = pos.to_numpy()[order]
    ordered.index = pd.RangeIndex(len(ordered))
    return ordered


def _build_provenance(
    *,
    profile: str,
    reference_root: Path,
    weight_root: Path,
    frequency_root: Path | None,
    reference_family: _DiscoveredFamily,
    weight_family: _DiscoveredFamily,
    frequency_prefix: str | None,
    selected_files: list[Path],
    intersections: dict[str, dict[str, int]],
    count_origins: dict[str, dict[str, str]],
    coordinate_provenance: dict[str, object],
    issues: list[dict[str, object]],
) -> dict[str, object]:
    unique_files = sorted(set(selected_files), key=lambda path: str(path))
    ignored_files = sorted(
        {
            str(row["file"])
            for row in issues
            if row.get("file") and str(row.get("reason", "")).startswith("discarded_")
        }
    )
    return {
        "profile": profile,
        "converter_version": "1",
        "source_directories": {
            "reference": str(reference_root),
            "weight": str(weight_root),
            "frequency": None if frequency_root is None else str(frequency_root),
        },
        "selected_prefixes": {
            "reference": reference_family.prefix,
            "weight": weight_family.prefix,
            **({"frequency": frequency_prefix} if frequency_prefix is not None else {}),
        },
        "selected_files": [str(path) for path in unique_files],
        "ignored_files": ignored_files,
        "intersection_counts": intersections,
        "count_origins": count_origins,
        "common_frequency_rule": dict(_COUNT_CONFIG),
        "coordinate_provenance": coordinate_provenance,
        "diagnostics": {
            "log": "diagnostics/convert-legacy-ldscores.log",
            "issues": "diagnostics/conversion_issues.tsv.gz",
        },
    }


def _issue(
    severity: str,
    source_role: str,
    *,
    file: str | Path | None = None,
    chromosome: int | str | None = None,
    snp: object = None,
    annotation: object = None,
    reason: str,
    observed: object = None,
    expected: object = None,
    details: object = None,
) -> dict[str, object]:
    return {
        "severity": severity,
        "source_role": source_role,
        "file": None if file is None else str(file),
        "chromosome": chromosome,
        "SNP": snp,
        "annotation": annotation,
        "reason": reason,
        "observed": observed,
        "expected": expected,
        "details": details,
    }


def _write_issues(path: Path, issues: list[dict[str, object]]) -> None:
    pd.DataFrame(issues, columns=_ISSUE_COLUMNS).to_csv(
        path,
        sep="\t",
        index=False,
        na_rep="",
        compression={"method": "gzip", "mtime": 0},
    )


def _decompressed_bytes(path: Path) -> bytes:
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as handle:
            return handle.read()
    return path.read_bytes()


__all__ = ["LegacyLDScoreConverter", "LegacyLDScoreConversionResult", "convert_legacy_ldscores"]
