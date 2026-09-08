"""Explicitly convert selected LDSC2 LD-score suites to LDSC3 directories.

The converter is the only public boundary that reads legacy ``.l2.ldscore``,
``.M``, ``.M_5_50``, ``.annot``, and frequency families. Regression never
parses those fragments directly; successful conversion writes an ordinary
canonical :class:`ldsc.ldscore_calculator.LDScoreResult` directory.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import gzip
import hashlib
import logging
from pathlib import Path
import re
from typing import Any

import numpy as np
import pandas as pd

from ._logging import log_inputs, log_outputs, materializing_overwrite_guard, workflow_logging
from .column_inference import normalize_genome_build, normalize_snp_identifier_mode
from .config import GlobalConfig
from .errors import LDSCInputError, LDSCUsageError
from .genome_build_inference import infer_chr_pos_build
from .ldscore_calculator import LDScoreResult
from .overlap_matrix import LDScoreOverlap
from .outputs import LDScoreDirectoryWriter, LDScoreOutputConfig
from .path_resolution import ensure_output_directory, ensure_output_paths_available, normalize_path_token


LOGGER = logging.getLogger("LDSC.legacy_ldscore_converter")
_CHROMS = tuple(range(1, 23))
_LDSCORE_RE = re.compile(r"^(?P<prefix>.*?)(?P<chrom>[1-9]|1[0-9]|2[0-2])\.l2\.ldscore(?P<gz>\.gz)?$")
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
_METADATA_COLUMNS = {"CHR", "SNP", "BP", "POS", "CM", "MAF"}
_COUNT_CONFIG = {
    "common_reference_snp_maf_min": 0.05,
    "common_reference_snp_maf_operator": ">",
    "common_reference_snp_semantics": "legacy_ldsc2",
}


@dataclass(frozen=True)
class LegacyLDScoreConversionResult:
    """Summary of one successful LDSC2 LD-score conversion.

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
    """Convert a reusable LDSC2 reference/weight suite into LDSC3 format.

    Source directories are treated as immutable. Matching across legacy
    families uses rsID; output coordinates come from the reference LD-score
    table. The converter accepts only allele-unaware ``rsid`` and ``chr_pos``
    identity modes and fixes legacy common-count semantics at ``0.05``.
    """

    @materializing_overwrite_guard(
        lambda self, **kwargs: (
            (kwargs["output_dir"], kwargs.get("overwrite", False), "RUN_FAILED.txt")
            if kwargs.get("output_dir")
            else None
        ),
        command="LegacyLDScoreConverter.convert(...)",
    )
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
        log_level: str = "INFO",
    ) -> LegacyLDScoreConversionResult:
        """Convert one complete autosomal legacy suite.

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
            Workflow log threshold. Default is ``"INFO"``.

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
        mode = normalize_snp_identifier_mode(snp_identifier)
        if mode not in {"rsid", "chr_pos"}:
            raise LDSCUsageError(
                "convert-ldsc2-ldscores accepts only allele-unaware --snp-identifier rsid or chr_pos; "
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
        log_path = diagnostics / "convert-ldsc2-ldscores.log"
        owned = [
            output_root / "metadata.json",
            output_root / "ldscore.baseline.parquet",
            output_root / "ldscore.query.parquet",
            output_root / "ldscore.overlap.parquet",
            issue_path,
            log_path,
        ]
        ensure_output_paths_available(owned, overwrite=overwrite, label="converted LD-score output artifact")
        issues: list[dict[str, object]] = []
        with workflow_logging("convert-ldsc2-ldscores", str(log_path), log_level=log_level):
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
                f"all chromosomes 1-22 (missing {missing}) and cannot fall back to unpartitioned conversion."
            )

        first_reference = _read_ldscore_table(reference_family.paths[1], "reference", 1)
        reference_columns = _scientific_columns(first_reference)
        if not reference_columns:
            raise LDSCInputError("Legacy reference LD-score table has no numeric scientific LD-score column.")
        if present_annotations:
            profile = "baseline_partitioned"
            if frequency_root is None:
                raise LDSCInputError(
                    "Baseline partitioned LDSC2 conversion requires --legacy-frequency-dir so common counts and "
                    "overlap can be reconstructed with the fixed legacy 0.05 threshold."
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
                overwrite=overwrite,
                issues=issues,
            )
        if len(reference_columns) != 1:
            raise LDSCInputError(
                "Unpartitioned LDSC2 reference suite must have exactly one scientific LD-score column when no "
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
                    f"Required .l2.M_5_50 count file is missing for chromosome {chrom}: '{common_path}'."
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
        if baseline.empty:
            raise LDSCInputError("Legacy reference and weight suites have an empty genome-wide rsID intersection.")
        _validate_global_output_identity(baseline, mode)
        baseline, effective_build, coordinate_provenance = _resolve_output_coordinates(
            baseline, mode, build_hint, context="legacy reference LD-score suite"
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
            ld_regression_snps=frozenset(baseline["SNP"].astype(str)),
            chromosome_results=[],
            count_config=dict(_COUNT_CONFIG),
            config_snapshot=GlobalConfig(snp_identifier=mode, genome_build=effective_build),
            legacy_ldsc2_import=provenance,
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
            frequency_payload = frequency.loc[:, ["SNP", freq_column]].copy()
            annotation_with_frequency = pd.merge(
                annotation,
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

            matrix = annotation.loc[:, [column_mapping[column] for column in reference_columns]].to_numpy(
                dtype=float
            )
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
                    f"Required .l2.M_5_50 count file is missing for chromosome {chrom}: '{common_path}'."
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
        if baseline.empty:
            raise LDSCInputError("Legacy reference and weight suites have an empty genome-wide rsID intersection.")
        _validate_global_output_identity(baseline, mode)
        baseline, effective_build, coordinate_provenance = _resolve_output_coordinates(
            baseline, mode, build_hint, context="legacy baseline reference LD-score suite"
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
            ld_regression_snps=frozenset(baseline["SNP"].astype(str)),
            chromosome_results=[],
            count_config=dict(_COUNT_CONFIG),
            config_snapshot=GlobalConfig(snp_identifier=mode, genome_build=effective_build),
            overlap=overlap,
            legacy_ldsc2_import=provenance,
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


def convert_ldsc2_ldscores(
    *,
    legacy_reference_dir: str | Path,
    legacy_weight_dir: str | Path,
    output_dir: str | Path,
    legacy_frequency_dir: str | Path | None = None,
    snp_identifier: str = "rsid",
    genome_build: str | None = "auto",
    overwrite: bool = False,
    log_level: str = "INFO",
) -> LegacyLDScoreConversionResult:
    """Convert one selected LDSC2 LD-score suite to a canonical LDSC3 directory.

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
        Converter log threshold. Default is ``"INFO"``.

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
    parser = argparse.ArgumentParser(prog="ldsc convert-ldsc2-ldscores")
    parser.add_argument("--legacy-reference-dir", required=True)
    parser.add_argument("--legacy-weight-dir", required=True)
    parser.add_argument("--legacy-frequency-dir")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--snp-identifier", choices=("rsid", "chr_pos"), default="rsid")
    parser.add_argument("--genome-build", choices=("auto", "hg19", "hg38"), default="auto")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"), default="INFO")
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
    return convert_ldsc2_ldscores(
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
        raise LDSCInputError(f"Cannot use {label} '{root}': path is not an existing directory.")
    return root


def _discover_ldscore_family(root: Path, role: str, issues: list[dict[str, object]]) -> _DiscoveredFamily:
    groups: dict[str, dict[int, list[Path]]] = {}
    for path in sorted(root.iterdir(), key=lambda item: item.name):
        match = _LDSCORE_RE.match(path.name)
        if match is None:
            if path.is_file() and ".l2.ldscore" in path.name:
                issues.append(_issue("warning", role, file=path, reason="discarded_non_family_file"))
            continue
        groups.setdefault(match.group("prefix"), {}).setdefault(int(match.group("chrom")), []).append(path)
    coherent = {prefix: shards for prefix, shards in groups.items() if set(shards) == set(_CHROMS)}
    if len(coherent) != 1:
        raise LDSCInputError(
            f"Legacy {role} directory '{root}' must contain exactly one coherent chromosome 1-22 "
            f".l2.ldscore family; found prefixes {sorted(coherent)}."
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
                f"{[str(path) for path in representations]}."
            )
        contents = [_decompressed_bytes(path) for path in representations]
        if contents[0] != contents[1]:
            raise LDSCInputError(
                f"Legacy {role} family has conflicting plain/gzip representations for chromosome {chrom}: "
                f"{[str(path) for path in representations]}."
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
                    f"Legacy suite has conflicting plain/gzip representations for chromosome {chrom}: {existing}."
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
            f"found prefixes {sorted(coherent)}."
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
            f"{[str(path) for path in representations]}."
        )
    return _DiscoveredFamily(prefix, selected)


def _read_ldscore_table(path: Path, role: str, chrom: int) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep=r"\s+", compression="infer")
    except Exception as exc:
        raise LDSCInputError(
            f"Could not parse legacy {role} LD-score table '{path}' for chromosome {chrom}: {exc}"
        ) from exc


def _read_whitespace_table(path: Path, role: str, chrom: int) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep=r"\s+", compression="infer")
    except Exception as exc:
        raise LDSCInputError(f"Could not parse legacy {role} table '{path}' for chromosome {chrom}: {exc}") from exc


def _scientific_columns(frame: pd.DataFrame) -> list[str]:
    return [column for column in frame.columns if str(column).upper() not in _METADATA_COLUMNS]


def _require_columns(frame: pd.DataFrame, required: set[str], role: str, path: Path, chrom: int) -> None:
    missing = sorted(required - set(frame.columns))
    if missing:
        raise LDSCInputError(
            f"Legacy {role} table '{path}' on chromosome {chrom} is missing required columns {missing}."
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
            f"Legacy {role} table '{path}' on chromosome {chrom} has duplicate rsIDs; cause SNPs: {sample}."
        )
    for column in scientific:
        values = pd.to_numeric(frame[column], errors="coerce")
        bad = ~np.isfinite(values.to_numpy(dtype=float))
        if bad.any():
            sample = frame.loc[bad, "SNP"].astype(str).head(10).tolist()
            raise LDSCInputError(
                f"Legacy {role} table '{path}' on chromosome {chrom} has non-finite values in '{column}'; "
                f"cause SNPs: {sample}."
            )
        frame[column] = values.astype(float)
    frame["BP"] = pd.to_numeric(frame["BP"], errors="coerce")
    bad_bp = frame["BP"].isna() | (frame["BP"] < 0)
    if bad_bp.any():
        sample = frame.loc[bad_bp, "SNP"].astype(str).head(10).tolist()
        raise LDSCInputError(
            f"Legacy {role} table '{path}' on chromosome {chrom} has invalid BP values; cause SNPs: {sample}."
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
                f"Cannot map legacy LD-score column '{column}' bijectively to a full annotation column; "
                f"candidates are {candidates}, annotation columns are {annotation_columns}."
            )
        mapping[column] = candidates[0]
    if len(set(mapping.values())) != len(annotation_columns) or set(mapping.values()) != set(annotation_columns):
        raise LDSCInputError(
            "Legacy baseline LD-score and annotation columns do not form a complete one-to-one mapping: "
            f"LD-score columns {reference_columns}, annotation columns {annotation_columns}, mapping {mapping}."
        )
    return mapping


def _validate_annotation_rows(frame: pd.DataFrame, columns: list[str], path: Path, chrom: int) -> None:
    duplicate = frame["SNP"].astype("string").duplicated(keep=False)
    if duplicate.any():
        sample = frame.loc[duplicate, "SNP"].astype(str).head(10).tolist()
        raise LDSCInputError(
            f"Legacy annotation table '{path}' on chromosome {chrom} has duplicate rsIDs; cause SNPs: {sample}."
        )
    for column in columns:
        values = pd.to_numeric(frame[column], errors="coerce")
        bad = ~np.isfinite(values.to_numpy(dtype=float))
        if bad.any():
            sample = frame.loc[bad, "SNP"].astype(str).head(10).tolist()
            raise LDSCInputError(
                f"Legacy annotation table '{path}' on chromosome {chrom} has non-finite values in '{column}'; "
                f"cause SNPs: {sample}."
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
            f"cause SNPs: {repeated[:10]}."
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
            f"found {candidates}."
        )
    if "SNP" not in frame.columns:
        raise LDSCInputError(f"Legacy frequency table '{path}' on chromosome {chrom} is missing required SNP column.")
    return candidates[0]


def _validate_frequency_rows(frame: pd.DataFrame, column: str, path: Path, chrom: int) -> None:
    duplicate = frame["SNP"].astype("string").duplicated(keep=False)
    if duplicate.any():
        sample = frame.loc[duplicate, "SNP"].astype(str).head(10).tolist()
        raise LDSCInputError(
            f"Legacy frequency table '{path}' on chromosome {chrom} has duplicate rsIDs; cause SNPs: {sample}."
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
    ref = reference.loc[:, ["CHR", "SNP", "BP", *reference_columns]].copy()
    wgt = weight.loc[:, ["CHR", "SNP", "BP", weight_column]].copy().rename(
        columns={"CHR": "weight_CHR", "BP": "weight_BP", weight_column: "regression_ld_scores"}
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
        columns={"BP": "POS"}
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
        raise LDSCInputError(f"Converted suite has duplicate genome-wide rsIDs; cause SNPs: {sample}.")
    if mode == "chr_pos":
        keys = frame["CHR"].astype(str) + ":" + frame["POS"].astype(str)
        duplicate_pos = keys.duplicated(keep=False)
        if duplicate_pos.any():
            sample = frame.loc[duplicate_pos, "SNP"].astype(str).head(10).tolist()
            raise LDSCInputError(f"Converted chr_pos suite has duplicate coordinates; cause SNPs: {sample}.")


def _resolve_output_coordinates(
    frame: pd.DataFrame,
    mode: str,
    build_hint: str | None,
    *,
    context: str,
) -> tuple[pd.DataFrame, str | None, dict[str, object]]:
    output = frame.copy()
    inference = None
    inference_error = None
    try:
        inference = infer_chr_pos_build(output.loc[:, ["CHR", "POS"]], context=context)
    except (ValueError, LDSCInputError) as exc:
        inference_error = str(exc)

    if mode == "rsid":
        if build_hint in {"hg19", "hg38"} and inference is not None and inference.genome_build != build_hint:
            raise LDSCInputError(
                f"Declared reference genome build {build_hint} conflicts with decisive {inference.genome_build} "
                f"coordinate evidence in {context}."
            )
        if inference is None:
            LOGGER.warning(
                "Could not infer reference build for rsID conversion; conversion continues because rsID is "
                "the canonical identity. Details: %s",
                inference_error,
            )
        return output, None, {
            "requested_build": build_hint,
            "effective_build": None,
            "inferred_build": None if inference is None else inference.genome_build,
            "coordinate_basis": "unknown" if inference is None else inference.coordinate_basis,
            "inference_error": inference_error,
        }

    if build_hint in {None, "auto"}:
        if inference is None:
            raise LDSCInputError(
                f"Could not infer the reference genome build for chr_pos conversion of {context}. "
                "Pass --genome-build hg19 or --genome-build hg38 explicitly."
            )
        effective_build = inference.genome_build
    else:
        effective_build = build_hint
        if inference is not None and inference.genome_build != effective_build:
            raise LDSCInputError(
                f"Declared reference genome build {effective_build} conflicts with decisive "
                f"{inference.genome_build} coordinate evidence in {context}."
            )
    basis = "declared_1-based" if inference is None else inference.coordinate_basis
    if basis == "0-based":
        output["POS"] = pd.to_numeric(output["POS"]).astype("int64") + 1
    return output, effective_build, {
        "requested_build": build_hint,
        "effective_build": effective_build,
        "inferred_build": None if inference is None else inference.genome_build,
        "coordinate_basis": basis,
        "inference_error": inference_error,
    }


def _sort_baseline(frame: pd.DataFrame) -> pd.DataFrame:
    ordered = frame.copy()
    ordered["CHR"] = ordered["CHR"].astype(str)
    ordered["POS"] = pd.to_numeric(ordered["POS"], errors="raise").astype("int64")
    ordered["_chrom_order"] = pd.to_numeric(ordered["CHR"], errors="coerce")
    ordered = ordered.sort_values(["_chrom_order", "POS", "SNP"], kind="mergesort").drop(columns="_chrom_order")
    return ordered.reset_index(drop=True)


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
        "source_sha256": {str(path): _sha256(path) for path in unique_files},
        "intersection_counts": intersections,
        "count_origins": count_origins,
        "common_frequency_rule": dict(_COUNT_CONFIG),
        "coordinate_provenance": coordinate_provenance,
        "diagnostics": {
            "log": "diagnostics/convert-ldsc2-ldscores.log",
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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


__all__ = ["LegacyLDScoreConverter", "LegacyLDScoreConversionResult", "convert_ldsc2_ldscores"]
