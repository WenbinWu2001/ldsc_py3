"""Validate gene-coordinate catalogs and resolve gene-list annotations in bulk.

The catalog supplied by the user, or embedded in a gene LD-score index, is the
only gene-resolution authority. Catalog coordinates remain one-based in public
diagnostics and are converted once to zero-based half-open intervals for SNP
projection. This module is pure workflow logic: it does not log or write files.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, replace
import glob
import gzip
from io import StringIO
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

from ._kernel.regions import load_preset_intervals
from .errors import LDSCInputError
from .path_resolution import normalize_path_token


CATALOG_REQUIRED_COLUMNS = ("gene_id", "gene_name", "chrom", "start", "end", "genome_build")
CATALOG_ISSUE_COLUMNS = (
    "source",
    "catalog_line",
    "gene_id",
    "gene_name",
    "chrom",
    "start",
    "end",
    "genome_build",
    "field",
    "reason",
    "related_catalog_lines",
    "observed_value",
    "details",
    "repair",
)
AUDIT_COLUMNS = (
    "argument",
    "input_role",
    "query",
    "source",
    "source_ordinal",
    "line",
    "input_gene",
    "match_type",
    "canonical_gene_id",
    "catalog_lines",
    "disposition",
    "reason",
    "chrom",
    "start",
    "end",
    "details",
)
SUMMARY_COLUMNS = (
    "argument",
    "input_role",
    "query",
    "source",
    "source_ordinal",
    "source_status",
    "source_reasons",
    "resolution_policy",
    "nonblank_input_rows",
    "uniquely_resolved_rows",
    "rejected_rows",
    "unique_resolved_genes",
    "duplicate_rows",
    "excluded_genes",
    "zero_support_genes",
    "genes_with_snp_support",
    "resolution_fraction",
)
RESOLUTION_POLICIES = ("strict", "resolved-only")
RESOLVED_ONLY_REASONS = frozenset(
    {
        "unmatched_identifier",
        "ambiguous_gene_name",
        "identifier_namespace_conflict",
        "catalog_duplicate_id",
        "catalog_invalid_coordinates",
        "outside_supported_chromosome",
    }
)
_AUTOSOMES = frozenset(str(value) for value in range(1, 23))
_BUILD_ALIASES = {
    "hg19": "hg19",
    "hg37": "hg19",
    "grch37": "hg19",
    "hg38": "hg38",
    "grch38": "hg38",
}


class GeneCatalogValidationError(LDSCInputError):
    """Catalog validation failure carrying the complete repairable issue table."""

    def __init__(self, message: str, issues: pd.DataFrame):
        super().__init__(message)
        self.issues = issues.loc[:, CATALOG_ISSUE_COLUMNS].copy()


@dataclass(frozen=True)
class GeneCatalog:
    """Normalized one-build gene-coordinate resolution authority."""

    frame: pd.DataFrame
    source: str
    genome_build: str
    issues: pd.DataFrame

    @classmethod
    def load(cls, path: str | Path, *, require_canonical: bool = False) -> "GeneCatalog":
        """Load a required TSV/TSV.GZ catalog and validate its full contents."""
        source_path = Path(path)
        source = source_path.name
        try:
            payload = source_path.read_bytes()
        except OSError as exc:
            issues = _single_catalog_issue(
                source,
                reason="catalog_unreadable",
                field="file",
                observed_value=str(exc).replace(str(source_path), source),
                repair="Make the catalog readable and rerun.",
            )
            raise GeneCatalogValidationError(
                f"Could not read gene-coordinate catalog '{source}': {exc}", issues
            ) from exc
        try:
            raw = gzip.decompress(payload) if source.lower().endswith(".gz") else payload
        except (OSError, EOFError) as exc:
            issues = _single_catalog_issue(
                source,
                reason="catalog_unparseable",
                field="file",
                observed_value="invalid gzip stream",
                repair="Regenerate the .tsv.gz catalog as a valid gzip-compressed TSV.",
            )
            raise GeneCatalogValidationError(
                f"Gene-coordinate catalog '{source}' is not a readable gzip TSV: {exc}", issues
            ) from exc
        try:
            text = raw.decode("utf-8")
        except UnicodeDecodeError as exc:
            issues = _single_catalog_issue(
                source,
                reason="catalog_unparseable",
                field="file",
                observed_value="invalid UTF-8",
                repair="Encode the catalog as UTF-8 and rerun.",
            )
            raise GeneCatalogValidationError(
                f"Gene-coordinate catalog '{source}' is not valid UTF-8: {exc}", issues
            ) from exc

        frame, parse_issues = _parse_catalog_text(text, source)
        if frame is None:
            reasons = ", ".join(sorted(set(parse_issues["reason"].astype(str))))
            if parse_issues["reason"].eq("missing_required_column").any():
                missing = ", ".join(parse_issues.loc[parse_issues["reason"].eq("missing_required_column"), "field"])
                message = f"Gene-coordinate catalog '{source}' is missing required columns: {missing}."
            else:
                message = (
                    f"Gene-coordinate catalog '{source}' is structurally unusable ({reasons}); "
                    "repair the catalog before resolving gene lists."
                )
            raise GeneCatalogValidationError(
                message,
                parse_issues,
            )
        normalized, validation_issues = _validate_catalog_rows(frame, source)
        issues = _sort_catalog_issues(pd.concat([parse_issues, validation_issues], ignore_index=True))
        structural = issues["reason"].isin(
            {
                "catalog_unreadable",
                "catalog_unparseable",
                "malformed_catalog_row",
                "missing_required_column",
                "missing_genome_build",
                "conflicting_genome_build",
                "unsupported_genome_build",
                "validation_incomplete",
            }
        )
        if structural.any() or (require_canonical and not issues.empty):
            scope = "canonical" if require_canonical else "structurally valid"
            raise GeneCatalogValidationError(
                f"Gene-coordinate catalog '{source}' is not {scope}: found {len(issues)} defect(s).",
                issues,
            )
        genome_build = str(normalized["genome_build"].dropna().iloc[0])
        return cls(
            frame=normalized.reset_index(drop=True),
            source=source,
            genome_build=genome_build,
            issues=issues.reset_index(drop=True),
        )

    @classmethod
    def from_embedded_frame(cls, frame: pd.DataFrame, *, source: str = "gene_catalog.parquet") -> "GeneCatalog":
        """Validate an embedded one-based catalog through the canonical path."""
        required = set(CATALOG_REQUIRED_COLUMNS) | {"catalog_line"}
        missing = sorted(required - set(frame.columns))
        if missing:
            issues = pd.DataFrame(
                [
                    _catalog_issue_dict(
                        source,
                        catalog_line=1,
                        field=column,
                        reason="missing_required_column",
                        observed_value="",
                        repair=f"Rebuild the index with a catalog containing '{column}'.",
                    )
                    for column in missing
                ],
                columns=CATALOG_ISSUE_COLUMNS,
            )
            raise GeneCatalogValidationError(
                f"Embedded gene catalog is missing required columns: {missing}.", issues
            )
        if frame.empty:
            issues = _single_catalog_issue(
                source,
                reason="catalog_unparseable",
                field="file",
                observed_value="empty embedded catalog",
                repair="Rebuild the index from a nonempty canonical coordinate catalog.",
            )
            raise GeneCatalogValidationError(
                f"Embedded gene catalog '{source}' is empty.", issues
            )
        catalog_lines = pd.to_numeric(frame["catalog_line"], errors="coerce")
        invalid_lines = (
            catalog_lines.isna()
            | catalog_lines.lt(2)
            | catalog_lines.mod(1).ne(0)
            | catalog_lines.duplicated(keep=False)
        )
        if invalid_lines.any():
            issue_rows = []
            for index in frame.index[invalid_lines]:
                value = frame.at[index, "catalog_line"]
                issue_rows.append(
                    _catalog_issue_dict(
                        source,
                        catalog_line=(value if pd.notna(value) else pd.NA),
                        row=frame.loc[index],
                        field="catalog_line",
                        reason="catalog_unparseable",
                        observed_value=("" if pd.isna(value) else str(value)),
                        details=(
                            "embedded catalog_line values must be unique integer physical "
                            "data lines greater than or equal to 2"
                        ),
                        repair="Rebuild the index from the original coordinate catalog.",
                    )
                )
            issues = _sort_catalog_issues(
                pd.DataFrame(issue_rows, columns=CATALOG_ISSUE_COLUMNS)
            )
            raise GeneCatalogValidationError(
                f"Embedded gene catalog '{source}' has invalid catalog_line values.",
                issues,
            )
        prepared = frame.loc[:, CATALOG_REQUIRED_COLUMNS].copy()
        prepared["catalog_line"] = catalog_lines.astype("Int64")
        normalized, issues = _validate_catalog_rows(prepared, source, preserve_catalog_lines=True)
        issues = _sort_catalog_issues(issues)
        if not issues.empty:
            raise GeneCatalogValidationError(
                f"Embedded gene catalog '{source}' is not canonical: found {len(issues)} defect(s).",
                issues,
            )
        return cls(
            frame=normalized.reset_index(drop=True),
            source=source,
            genome_build=str(normalized["genome_build"].iloc[0]),
            issues=issues,
        )


@dataclass(frozen=True)
class GeneSourceSelection:
    """Unique usable genes selected for one focal or control source."""

    input_role: str
    argument: str
    query: str
    source: str
    source_ordinal: int
    canonical_gene_ids: tuple[str, ...]
    catalog_indices: tuple[int, ...]
    intervals: tuple[tuple[str, int, int], ...]


@dataclass(frozen=True)
class GeneListBatchResolution:
    """Complete Gate A result for all focal sources and the optional control."""

    audit: pd.DataFrame
    summary: pd.DataFrame
    selections: tuple[GeneSourceSelection, ...]
    resolution_policy: str
    has_fatal_gate_a_issues: bool

    def selection(self, input_role: str, source_ordinal: int) -> GeneSourceSelection:
        """Return one declared source selection by stable role and ordinal."""
        for selection in self.selections:
            if selection.input_role == input_role and selection.source_ordinal == source_ordinal:
                return selection
        raise KeyError((input_role, source_ordinal))

    def with_snp_support(self, support_by_catalog_index: pd.Series | dict[int, int]) -> "GeneListBatchResolution":
        """Return Gate B audit/summary fields updated from per-catalog-row SNP counts."""
        support = pd.Series(support_by_catalog_index, dtype="Int64")
        audit = self.audit.copy()
        selection_map = pd.DataFrame(
            {
                "input_role": [item.input_role for item in self.selections],
                "source_ordinal": [item.source_ordinal for item in self.selections],
                "canonical_gene_id": [item.canonical_gene_ids for item in self.selections],
                "_catalog_index": [item.catalog_indices for item in self.selections],
            }
        ).explode(["canonical_gene_id", "_catalog_index"], ignore_index=True)
        selection_map["_support"] = selection_map["_catalog_index"].map(support)
        counts = audit.merge(
            selection_map.drop(columns="_catalog_index"),
            on=["input_role", "source_ordinal", "canonical_gene_id"],
            how="left",
            sort=False,
        )["_support"].astype("Int64")
        unsupported = audit["disposition"].eq("retained") & counts.fillna(0).eq(0)
        audit.loc[unsupported, "disposition"] = "unsupported"
        audit.loc[unsupported, "reason"] = "zero_reference_snp_support"
        audit.loc[unsupported, "details"] = "No retained reference-panel SNP overlaps this gene interval."
        summary = _summarize_sources(audit, self.summary, self.resolution_policy, support_evaluated=True)
        return replace(self, audit=audit.loc[:, AUDIT_COLUMNS], summary=summary)


def gene_list_query_name(path: str | Path) -> str:
    """Derive a stable query name from one gene-list basename."""
    name = Path(path).name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for suffix in (".txt", ".tsv", ".list"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return name


def _expand_focal_gene_list_sources(
    tokens: Sequence[str | Path],
) -> tuple[str, ...]:
    """Expand focal globs deterministically while retaining unresolved tokens for Gate A."""
    expanded: list[str] = []
    for raw_token in tokens:
        token = normalize_path_token(raw_token)
        if Path(token).is_file():
            matches = [token]
        elif glob.has_magic(token):
            matches = sorted(path for path in glob.glob(token) if Path(path).is_file())
            if not matches:
                matches = [token]
        else:
            matches = [token]
        expanded.extend(matches)
    return tuple(expanded)


def resolve_gene_lists(
    focal_paths: Sequence[str | Path],
    catalog: GeneCatalog,
    *,
    control_path: str | Path | None = None,
    resolution_policy: str = "strict",
    gene_exclude_regions: str = "none",
    index_chromosome_coverage: Iterable[str] | None = None,
) -> GeneListBatchResolution:
    """Resolve all focal/control sources together using vectorized table operations."""
    if resolution_policy not in RESOLUTION_POLICIES:
        raise LDSCInputError(
            f"Unsupported gene-list resolution policy {resolution_policy!r}; use 'strict' or 'resolved-only'."
        )
    if gene_exclude_regions not in {"none", "mhc"}:
        raise LDSCInputError("gene_exclude_regions must be 'none' or 'mhc'.")

    focal_paths = _expand_focal_gene_list_sources(focal_paths)
    if control_path is not None:
        control_path = normalize_path_token(control_path)

    declarations = [
        {
            "argument": "--query-annot-gene-list-sources",
            "input_role": "focal",
            "query": gene_list_query_name(path),
            "source": Path(path).name,
            "source_path": str(path),
            "source_ordinal": ordinal,
        }
        for ordinal, path in enumerate(focal_paths, start=1)
    ]
    if control_path is not None:
        declarations.append(
            {
                "argument": "--control-gene-list-file",
                "input_role": "control",
                "query": "gene_control",
                "source": Path(control_path).name,
                "source_path": str(control_path),
                "source_ordinal": 0,
            }
        )
    declaration_frame = pd.DataFrame(declarations)
    if declaration_frame.empty:
        return GeneListBatchResolution(
            audit=pd.DataFrame(columns=AUDIT_COLUMNS),
            summary=pd.DataFrame(columns=SUMMARY_COLUMNS),
            selections=(),
            resolution_policy=resolution_policy,
            has_fatal_gate_a_issues=False,
        )

    source_frames: list[pd.DataFrame] = []
    source_errors: dict[tuple[str, int], str] = {}
    for declaration in declarations:
        source_frame, source_reason = _read_gene_list_source(declaration)
        if source_frame is not None:
            source_frames.append(source_frame)
        if source_reason:
            source_errors[(declaration["input_role"], declaration["source_ordinal"])] = source_reason
    rows = pd.concat(source_frames, ignore_index=True) if source_frames else pd.DataFrame()
    if rows.empty:
        rows = pd.DataFrame(
            columns=[
                "argument",
                "input_role",
                "query",
                "source",
                "source_ordinal",
                "line",
                "input_gene",
                "_malformed",
            ]
        )

    duplicate_names = declaration_frame["query"].duplicated(keep=False) & declaration_frame["input_role"].eq("focal")
    for row in declaration_frame.loc[duplicate_names].itertuples(index=False):
        source_errors[(row.input_role, row.source_ordinal)] = _combine_reasons(
            source_errors.get((row.input_role, row.source_ordinal), ""), "duplicate_query_name"
        )

    audit = _resolve_rows(rows, catalog, gene_exclude_regions, index_chromosome_coverage)
    source_seed = declaration_frame.drop(columns="source_path").copy()
    source_seed["source_status"] = [
        "error" if source_errors.get((role, ordinal)) else "ok"
        for role, ordinal in zip(source_seed["input_role"], source_seed["source_ordinal"], strict=True)
    ]
    source_seed["source_reasons"] = [
        source_errors.get((role, ordinal), "")
        for role, ordinal in zip(source_seed["input_role"], source_seed["source_ordinal"], strict=True)
    ]
    summary = _summarize_sources(audit, source_seed, resolution_policy, support_evaluated=False)
    selections = _build_selections(audit, catalog, declaration_frame)

    source_fatal = bool(source_errors)
    malformed_fatal = bool(audit["reason"].eq("malformed_input").any())
    rejected_reasons = set(audit.loc[audit["disposition"].eq("rejected"), "reason"].astype(str))
    policy_fatal = bool(rejected_reasons) and (
        resolution_policy == "strict" or not rejected_reasons.issubset(RESOLVED_ONLY_REASONS)
    )
    empty_control = False
    if control_path is not None:
        control_summary = summary[summary["input_role"].eq("control")]
        if not control_summary.empty and control_summary.iloc[0]["source_status"] == "ok":
            selected_control = next(item for item in selections if item.input_role == "control")
            empty_control = len(selected_control.canonical_gene_ids) == 0
    return GeneListBatchResolution(
        audit=audit.loc[:, AUDIT_COLUMNS].reset_index(drop=True),
        summary=summary.loc[:, SUMMARY_COLUMNS].reset_index(drop=True),
        selections=selections,
        resolution_policy=resolution_policy,
        has_fatal_gate_a_issues=source_fatal or malformed_fatal or policy_fatal or empty_control,
    )


def select_index_eligible_gene_indices(
    catalog: GeneCatalog,
    *,
    gene_exclude_regions: str = "none",
) -> tuple[int, ...]:
    """Return canonical catalog rows not removed by the index exclusion policy."""
    if gene_exclude_regions not in {"none", "mhc"}:
        raise LDSCInputError("gene_exclude_regions must be 'none' or 'mhc'.")
    eligible = catalog.frame["_catalog_valid"].copy()
    if gene_exclude_regions == "mhc":
        eligible &= ~_mhc_overlap_mask(catalog.frame, catalog.genome_build)
    return tuple(eligible.index[eligible].astype(int))


def _parse_catalog_text(text: str, source: str) -> tuple[pd.DataFrame | None, pd.DataFrame]:
    physical = [
        (line_number, line)
        for line_number, line in enumerate(text.splitlines(), start=1)
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not physical:
        issues = _single_catalog_issue(
            source,
            reason="catalog_unparseable",
            field="file",
            observed_value="empty catalog",
            repair="Provide a headered TSV with gene-coordinate rows.",
        )
        return None, issues
    header_line, header_text = physical[0]
    try:
        header = [value.strip() for value in next(csv.reader([header_text], delimiter="\t"))]
    except csv.Error as exc:
        issues = _single_catalog_issue(
            source,
            catalog_line=header_line,
            reason="catalog_unparseable",
            field="header",
            observed_value=header_text,
            details=str(exc),
            repair="Regenerate a valid tab-delimited header.",
        )
        return None, issues
    duplicate_headers = sorted(
        {value for value in header if header.count(value) > 1}
    )
    if duplicate_headers:
        issues = _single_catalog_issue(
            source,
            catalog_line=header_line,
            reason="catalog_unparseable",
            field="header",
            observed_value=",".join(duplicate_headers),
            details="column names must be unique",
            repair="Regenerate the header with each required column exactly once.",
        )
        return None, issues
    missing = [column for column in CATALOG_REQUIRED_COLUMNS if column not in header]
    if missing:
        issue_rows = [
                _catalog_issue_dict(
                    source,
                    catalog_line=header_line,
                    field=column,
                    reason="missing_required_column",
                    observed_value="",
                    repair=f"Add the required '{column}' column without changing coordinate semantics.",
                )
                for column in missing
            ]
        issue_rows.append(
            _catalog_issue_dict(
                source,
                catalog_line=header_line,
                field="validation",
                reason="validation_incomplete",
                observed_value=",".join(missing),
                details="row and global canonicality checks were not run because required columns are absent",
                repair="Add the missing columns, then rerun catalog validation to reveal any remaining defects.",
            )
        )
        issues = pd.DataFrame(issue_rows, columns=CATALOG_ISSUE_COLUMNS)
        return None, issues
    if len(physical) == 1:
        issues = _single_catalog_issue(
            source,
            catalog_line=header_line,
            reason="catalog_unparseable",
            field="file",
            observed_value="header without data rows",
            repair="Provide at least one gene-coordinate row below the header.",
        )
        return None, issues
    records: list[dict[str, Any]] = []
    parse_issues: list[dict[str, Any]] = []
    for line_number, raw_line in physical[1:]:
        try:
            values = next(csv.reader([raw_line], delimiter="\t"))
        except csv.Error as exc:
            values = []
            parse_issues.append(
                _catalog_issue_dict(
                    source,
                    catalog_line=line_number,
                    field="row",
                    reason="malformed_catalog_row",
                    observed_value=raw_line,
                    details=str(exc),
                    repair="Regenerate this row as valid tab-delimited text.",
                )
            )
        if len(values) != len(header):
            parse_issues.append(
                _catalog_issue_dict(
                    source,
                    catalog_line=line_number,
                    field="row",
                    reason="malformed_catalog_row",
                    observed_value=raw_line,
                    details=f"expected {len(header)} fields, found {len(values)}",
                    repair="Make this row contain exactly one value for every header column.",
                )
            )
            values = (values + [""] * len(header))[: len(header)]
        record = dict(zip(header, values, strict=True))
        record["catalog_line"] = line_number
        records.append(record)
    frame = pd.DataFrame(records, columns=[*header, "catalog_line"])
    issues = pd.DataFrame(parse_issues, columns=CATALOG_ISSUE_COLUMNS)
    return frame, issues


def _validate_catalog_rows(
    frame: pd.DataFrame,
    source: str,
    *,
    preserve_catalog_lines: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    normalized = frame.copy()
    if "catalog_line" not in normalized:
        normalized["catalog_line"] = np.arange(2, len(normalized) + 2, dtype=int)
    if not preserve_catalog_lines:
        normalized["catalog_line"] = pd.to_numeric(normalized["catalog_line"], errors="coerce").astype("Int64")
    for column in CATALOG_REQUIRED_COLUMNS:
        normalized[column] = normalized[column].fillna("").astype(str).str.strip()
    raw = normalized.loc[:, CATALOG_REQUIRED_COLUMNS].copy()
    normalized["gene_id"] = raw["gene_id"]
    normalized["gene_name"] = raw["gene_name"]
    normalized["genome_build"] = raw["genome_build"].str.lower().map(_BUILD_ALIASES)
    chrom_text = raw["chrom"].str.lower().str.replace(r"^chr", "", regex=True)
    chrom_numeric = pd.to_numeric(chrom_text.where(chrom_text.str.fullmatch(r"[0-9]+")), errors="coerce")
    normalized["chrom"] = chrom_numeric.astype("Int64").astype("string")
    normalized["start"] = pd.to_numeric(
        raw["start"].where(raw["start"].str.fullmatch(r"[0-9]+")), errors="coerce"
    ).astype("Int64")
    normalized["end"] = pd.to_numeric(
        raw["end"].where(raw["end"].str.fullmatch(r"[0-9]+")), errors="coerce"
    ).astype("Int64")
    normalized["start0"] = normalized["start"] - 1

    issues: list[dict[str, Any]] = []

    def add_mask(mask: pd.Series, field: str, reason: str, repair: str, details: str = "") -> None:
        for index in normalized.index[mask.fillna(False)]:
            issues.append(
                _catalog_issue_dict(
                    source,
                    catalog_line=normalized.at[index, "catalog_line"],
                    row=raw.loc[index],
                    field=field,
                    reason=reason,
                    observed_value=raw.at[index, field] if field in raw.columns else "",
                    details=details,
                    repair=repair,
                )
            )

    add_mask(raw["gene_id"].eq(""), "gene_id", "missing_gene_id", "Restore the authoritative gene ID or remove the row upstream.")
    add_mask(raw["genome_build"].eq(""), "genome_build", "missing_genome_build", "Declare hg19/GRCh37 or hg38/GRCh38 for every row.")
    add_mask(
        raw["genome_build"].ne("") & normalized["genome_build"].isna(),
        "genome_build",
        "unsupported_genome_build",
        "Use hg19/GRCh37 or hg38/GRCh38 coordinates.",
    )
    valid_builds = normalized["genome_build"].dropna().unique()
    if len(valid_builds) > 1:
        add_mask(
            normalized["genome_build"].notna(),
            "genome_build",
            "conflicting_genome_build",
            "Regenerate the catalog from one internally consistent genome build.",
            details="catalog contains both hg19 and hg38 declarations",
        )
    invalid_chrom = raw["chrom"].ne("") & chrom_numeric.isna()
    add_mask(invalid_chrom, "chrom", "invalid_chromosome", "Use an autosomal chromosome label from 1 through 22.")
    outside_chrom = chrom_numeric.notna() & ~normalized["chrom"].isin(_AUTOSOMES)
    add_mask(outside_chrom, "chrom", "outside_supported_chromosome", "Remove non-autosomal rows for this workflow.")
    invalid_start = normalized["start"].isna() | normalized["start"].lt(1)
    add_mask(invalid_start, "start", "invalid_start", "Use an integer one-based start greater than or equal to 1.")
    invalid_end = normalized["end"].isna() | normalized["end"].lt(1)
    add_mask(invalid_end, "end", "invalid_end", "Use an integer one-based inclusive end greater than or equal to 1.")
    end_before = normalized["start"].notna() & normalized["end"].notna() & normalized["end"].lt(normalized["start"])
    add_mask(end_before, "end", "end_before_start", "Correct the interval so end is greater than or equal to start.")

    _append_duplicate_issues(issues, normalized, raw, source, "gene_id", "duplicate_gene_id", allow_blank=False)
    _append_duplicate_issues(issues, normalized, raw, source, "gene_name", "duplicate_gene_name", allow_blank=True)
    _append_namespace_issues(issues, normalized, raw, source)

    issue_frame = pd.DataFrame(issues, columns=CATALOG_ISSUE_COLUMNS)
    duplicate_ids = normalized["gene_id"].ne("") & normalized["gene_id"].duplicated(keep=False)
    row_issue_reasons = issue_frame[issue_frame["catalog_line"].notna()].groupby("catalog_line", sort=False)["reason"].agg(
        lambda values: ",".join(dict.fromkeys(values.astype(str)))
    ) if not issue_frame.empty else pd.Series(dtype=str)
    normalized["_row_issue_reasons"] = normalized["catalog_line"].map(row_issue_reasons).fillna("")
    normalized["_duplicate_gene_id"] = duplicate_ids.to_numpy()
    duplicate_id_lines = (
        normalized.loc[duplicate_ids]
        .groupby("gene_id", sort=False)["catalog_line"]
        .agg(lambda lines: ",".join(str(int(value)) for value in sorted(lines.dropna().astype(int))))
    )
    normalized["_duplicate_gene_id_lines"] = normalized["gene_id"].map(duplicate_id_lines).fillna("")
    normalized["_outside_supported_chromosome"] = (invalid_chrom | outside_chrom).to_numpy()
    coordinate_problem = invalid_chrom | outside_chrom | invalid_start | invalid_end | end_before
    normalized["_catalog_valid"] = (
        raw["gene_id"].ne("")
        & normalized["genome_build"].notna()
        & ~duplicate_ids
        & ~coordinate_problem
    )
    return normalized, _sort_catalog_issues(issue_frame)


def _append_duplicate_issues(
    issues: list[dict[str, Any]],
    normalized: pd.DataFrame,
    raw: pd.DataFrame,
    source: str,
    field: str,
    reason: str,
    *,
    allow_blank: bool,
) -> None:
    values = normalized[field]
    duplicate = values.duplicated(keep=False) & (values.ne("") if allow_blank else pd.Series(True, index=values.index))
    if not duplicate.any():
        return
    related = (
        normalized.loc[duplicate].groupby(field, sort=False)["catalog_line"].agg(
            lambda lines: ",".join(str(int(value)) for value in sorted(lines.dropna().astype(int)))
        )
    )
    for index in normalized.index[duplicate]:
        issues.append(
            _catalog_issue_dict(
                source,
                catalog_line=normalized.at[index, "catalog_line"],
                row=raw.loc[index],
                field=field,
                reason=reason,
                related_catalog_lines=related[normalized.at[index, field]],
                observed_value=raw.at[index, field],
                details=f"{field} must identify at most one catalog row",
                repair=f"Regenerate the catalog so {field} is unique.",
            )
        )


def _append_namespace_issues(
    issues: list[dict[str, Any]],
    normalized: pd.DataFrame,
    raw: pd.DataFrame,
    source: str,
) -> None:
    ids = normalized.loc[normalized["gene_id"].ne(""), ["gene_id", "catalog_line"]].rename(
        columns={"gene_id": "token", "catalog_line": "id_line"}
    )
    names = normalized.loc[normalized["gene_name"].ne(""), ["gene_name", "catalog_line"]].rename(
        columns={"gene_name": "token", "catalog_line": "name_line"}
    )
    conflicts = ids.merge(names, on="token", how="inner")
    conflicts = conflicts[conflicts["id_line"] != conflicts["name_line"]]
    if conflicts.empty:
        return
    for token, group in conflicts.groupby("token", sort=True):
        lines = sorted(set(group["id_line"].astype(int)) | set(group["name_line"].astype(int)))
        related = ",".join(str(value) for value in lines)
        for line in lines:
            index = normalized.index[normalized["catalog_line"].astype("Int64").eq(line)][0]
            issues.append(
                _catalog_issue_dict(
                    source,
                    catalog_line=line,
                    row=raw.loc[index],
                    field="gene_id,gene_name",
                    reason="identifier_namespace_conflict",
                    related_catalog_lines=related,
                    observed_value=token,
                    details="token is a gene ID on one row and a gene name on another",
                    repair="Use disjoint gene-ID and gene-name namespaces in the catalog transformation.",
                )
            )


def _read_gene_list_source(declaration: dict[str, Any]) -> tuple[pd.DataFrame | None, str]:
    path = Path(declaration["source_path"])
    if any(token in str(path) for token in ("*", "?", "[", "]")):
        return None, "unreadable_gene_list"
    try:
        payload = path.read_bytes()
    except OSError:
        return None, "unreadable_gene_list"
    try:
        raw = gzip.decompress(payload) if path.name.lower().endswith(".gz") else payload
    except (OSError, EOFError):
        return None, "invalid_gzip"
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        return None, "invalid_utf8"
    lines = pd.Series(text.splitlines(), dtype="string")
    if lines.empty:
        return pd.DataFrame(columns=[*AUDIT_COLUMNS[:7], "_malformed"]), ""
    stripped = lines.str.strip()
    keep = stripped.ne("")
    frame = pd.DataFrame(
        {
            "line": np.arange(1, len(lines) + 1, dtype=int)[keep.to_numpy()],
            "input_gene": stripped[keep].to_numpy(),
            "_malformed": lines[keep].str.contains("\t", regex=False).to_numpy(),
        }
    )
    for field in ("argument", "input_role", "query", "source", "source_ordinal"):
        frame[field] = declaration[field]
    return frame[["argument", "input_role", "query", "source", "source_ordinal", "line", "input_gene", "_malformed"]], ""


def _resolve_rows(
    rows: pd.DataFrame,
    catalog: GeneCatalog,
    gene_exclude_regions: str,
    index_chromosome_coverage: Iterable[str] | None,
) -> pd.DataFrame:
    if rows.empty:
        return pd.DataFrame(columns=AUDIT_COLUMNS)
    working = rows.copy().reset_index(drop=True)
    working["_row_id"] = np.arange(len(working), dtype=int)
    candidates_base = catalog.frame.reset_index().rename(columns={"index": "_catalog_index"})
    candidate_columns = [
        "_catalog_index",
        "gene_id",
        "gene_name",
        "catalog_line",
        "chrom",
        "start",
        "end",
        "start0",
        "_catalog_valid",
        "_duplicate_gene_id",
        "_duplicate_gene_id_lines",
        "_outside_supported_chromosome",
        "_row_issue_reasons",
    ]
    resolvable = working.loc[~working["_malformed"], ["_row_id", "input_gene"]]
    id_matches = resolvable.merge(
        candidates_base[candidate_columns], left_on="input_gene", right_on="gene_id", how="inner", sort=False
    )
    id_matches["_via"] = "gene_id"
    name_matches = resolvable.merge(
        candidates_base[candidate_columns], left_on="input_gene", right_on="gene_name", how="inner", sort=False
    )
    name_matches["_via"] = "gene_name"
    candidates = pd.concat([id_matches, name_matches], ignore_index=True)
    candidates = candidates.drop_duplicates(["_row_id", "_catalog_index", "_via"])

    audit = working[["argument", "input_role", "query", "source", "source_ordinal", "line", "input_gene"]].copy()
    audit["match_type"] = "unmatched"
    audit["canonical_gene_id"] = pd.NA
    audit["catalog_lines"] = pd.NA
    audit["disposition"] = "rejected"
    audit["reason"] = "unmatched_identifier"
    audit["chrom"] = pd.NA
    audit["start"] = pd.NA
    audit["end"] = pd.NA
    audit["details"] = pd.NA
    malformed = working["_malformed"].to_numpy()
    audit.loc[malformed, ["match_type", "reason", "details"]] = [
        "unmatched",
        "malformed_input",
        "expected exactly one identifier field; remove tab-separated extra fields",
    ]

    if not candidates.empty:
        unique_candidates = candidates.sort_values(["_row_id", "catalog_line", "_via"]).drop_duplicates(
            ["_row_id", "_catalog_index"]
        )
        via_counts = candidates.groupby(["_row_id", "_via"], sort=False).size().unstack(fill_value=0)
        groups = unique_candidates.groupby("_row_id", sort=False)
        stats = groups.agg(
            candidate_count=("_catalog_index", "nunique"),
            duplicate_id=("_duplicate_gene_id", "max"),
            duplicate_id_lines=("_duplicate_gene_id_lines", "first"),
            catalog_valid=("_catalog_valid", "all"),
            outside_supported=("_outside_supported_chromosome", "max"),
        )
        stats["id_count"] = via_counts.get("gene_id", pd.Series(0, index=stats.index)).reindex(stats.index, fill_value=0)
        stats["name_count"] = via_counts.get("gene_name", pd.Series(0, index=stats.index)).reindex(stats.index, fill_value=0)
        catalog_lines = (
            unique_candidates.assign(_line=unique_candidates["catalog_line"].astype("Int64").astype(str))
            .groupby("_row_id", sort=False)["_line"]
            .agg(",".join)
        )
        conflict_ids = (
            unique_candidates[["_row_id", "gene_id"]]
            .drop_duplicates()
            .sort_values(["_row_id", "gene_id"])
            .groupby("_row_id", sort=False)["gene_id"]
            .agg(",".join)
        )
        stats["catalog_lines"] = catalog_lines
        stats["conflict_ids"] = conflict_ids
        stats["match_type"] = np.where(stats["candidate_count"].gt(1), "ambiguous", np.where(stats["id_count"].gt(0), "gene_id", "gene_name"))
        stats["reason"] = ""
        duplicate_id = stats["duplicate_id"].astype(bool) & (
            stats["id_count"].gt(0) | stats["candidate_count"].eq(1)
        )
        no_id_ambiguous = stats["candidate_count"].gt(1) & stats["id_count"].eq(0)
        namespace = (
            stats["candidate_count"].gt(1)
            & stats["id_count"].gt(0)
            & ~duplicate_id
        )
        stats.loc[duplicate_id, "reason"] = "catalog_duplicate_id"
        stats.loc[no_id_ambiguous, "reason"] = "ambiguous_gene_name"
        stats.loc[namespace, "reason"] = "identifier_namespace_conflict"
        outside = stats["candidate_count"].eq(1) & stats["reason"].eq("") & stats["outside_supported"].astype(bool)
        stats.loc[outside, "reason"] = "outside_supported_chromosome"
        invalid = stats["candidate_count"].eq(1) & stats["reason"].eq("") & ~stats["catalog_valid"].astype(bool)
        stats.loc[invalid, "reason"] = "catalog_invalid_coordinates"

        coverage = None if index_chromosome_coverage is None else {str(value) for value in index_chromosome_coverage}
        unique_rows = unique_candidates.drop_duplicates("_row_id").set_index("_row_id")
        if coverage is not None:
            outside_coverage = stats["candidate_count"].eq(1) & stats["reason"].eq("") & ~unique_rows["chrom"].astype(str).isin(coverage)
            stats.loc[outside_coverage, "reason"] = "outside_index_chromosome_coverage"
        matched_ids = stats.index
        audit.loc[matched_ids, "match_type"] = stats["match_type"]
        audit.loc[matched_ids, "catalog_lines"] = stats["catalog_lines"]
        audit.loc[matched_ids, "reason"] = stats["reason"]
        conflicts = stats["candidate_count"].gt(1)
        audit.loc[stats.index[conflicts], "details"] = "conflicting gene IDs: " + stats.loc[conflicts, "conflict_ids"]
        duplicate_id_conflicts = stats["reason"].eq("catalog_duplicate_id")
        audit.loc[stats.index[duplicate_id_conflicts], "catalog_lines"] = stats.loc[
            duplicate_id_conflicts, "duplicate_id_lines"
        ]
        audit.loc[stats.index[duplicate_id_conflicts], "details"] = (
            "duplicate gene ID on catalog lines "
            + stats.loc[duplicate_id_conflicts, "duplicate_id_lines"]
        )
        unique = stats["candidate_count"].eq(1)
        unique_index = stats.index[unique]
        audit.loc[unique_index, "canonical_gene_id"] = unique_rows.loc[unique_index, "gene_id"].to_numpy()
        audit.loc[unique_index, "chrom"] = unique_rows.loc[unique_index, "chrom"].to_numpy()
        audit.loc[unique_index, "start"] = unique_rows.loc[unique_index, "start"].to_numpy()
        audit.loc[unique_index, "end"] = unique_rows.loc[unique_index, "end"].to_numpy()
        invalid_unique = unique & stats["reason"].ne("") & ~duplicate_id_conflicts
        audit.loc[stats.index[invalid_unique], "details"] = unique_rows.loc[stats.index[invalid_unique], "_row_issue_reasons"].to_numpy()

    usable = audit["reason"].eq("") & audit["canonical_gene_id"].notna()
    audit.loc[usable, "disposition"] = "retained"
    if gene_exclude_regions == "mhc" and usable.any():
        mhc_mask = _mhc_overlap_mask(audit.loc[usable], catalog.genome_build)
        excluded_indices = audit.loc[usable].index[mhc_mask.to_numpy()]
        audit.loc[excluded_indices, "disposition"] = "excluded"
        audit.loc[excluded_indices, "reason"] = "excluded_gene_region"
        audit.loc[excluded_indices, "details"] = "Unpadded gene interval overlaps the MHC exclusion."

    resolved = audit["canonical_gene_id"].notna() & ~audit["disposition"].eq("rejected")
    duplicate_rank = audit.loc[resolved].groupby(
        ["input_role", "source_ordinal", "canonical_gene_id"], sort=False
    ).cumcount()
    duplicate_indices = duplicate_rank.index[duplicate_rank.gt(0)]
    if len(duplicate_indices):
        first_lines = audit.loc[resolved].groupby(
            ["input_role", "source_ordinal", "canonical_gene_id"], sort=False
        )["line"].transform("first")
        audit.loc[duplicate_indices, "disposition"] = "duplicate"
        audit.loc[duplicate_indices, "reason"] = "duplicate_canonical_gene"
        audit.loc[duplicate_indices, "details"] = (
            "first occurrence at line " + first_lines.loc[duplicate_indices].astype(str)
        )
    audit["_role_order"] = audit["input_role"].map({"focal": 0, "control": 1})
    audit = audit.sort_values(["_role_order", "source_ordinal", "line"], kind="stable").drop(columns="_role_order")
    return audit.loc[:, AUDIT_COLUMNS].reset_index(drop=True)


def _build_selections(
    audit: pd.DataFrame,
    catalog: GeneCatalog,
    declarations: pd.DataFrame,
) -> tuple[GeneSourceSelection, ...]:
    catalog_lookup = catalog.frame.reset_index().set_index("gene_id")
    selections: list[GeneSourceSelection] = []
    ordered = declarations.assign(_role_order=declarations["input_role"].map({"focal": 0, "control": 1})).sort_values(
        ["_role_order", "source_ordinal"], kind="stable"
    )
    for declaration in ordered.itertuples(index=False):
        selected = audit[
            (audit["input_role"] == declaration.input_role)
            & (audit["source_ordinal"] == declaration.source_ordinal)
            & audit["disposition"].eq("retained")
        ]
        gene_ids = tuple(selected["canonical_gene_id"].astype(str))
        if gene_ids:
            rows = catalog_lookup.loc[list(gene_ids)]
            catalog_indices = tuple(rows["index"].astype(int))
            intervals = tuple(
                zip(rows["chrom"].astype(str), rows["start0"].astype(int), rows["end"].astype(int), strict=True)
            )
        else:
            catalog_indices = ()
            intervals = ()
        selections.append(
            GeneSourceSelection(
                input_role=declaration.input_role,
                argument=declaration.argument,
                query=declaration.query,
                source=declaration.source,
                source_ordinal=int(declaration.source_ordinal),
                canonical_gene_ids=gene_ids,
                catalog_indices=catalog_indices,
                intervals=intervals,
            )
        )
    return tuple(selections)


def _summarize_sources(
    audit: pd.DataFrame,
    source_seed: pd.DataFrame,
    resolution_policy: str,
    *,
    support_evaluated: bool,
) -> pd.DataFrame:
    seed_columns = ["argument", "input_role", "query", "source", "source_ordinal", "source_status", "source_reasons"]
    summary = source_seed.loc[:, seed_columns].copy()
    summary["resolution_policy"] = resolution_policy
    key = ["input_role", "source_ordinal"]
    if audit.empty:
        counts = pd.DataFrame(
            columns=[
                *key,
                "nonblank_input_rows",
                "uniquely_resolved_rows",
                "rejected_rows",
                "unique_resolved_genes",
                "duplicate_rows",
                "excluded_genes",
                "zero_support_genes",
                "genes_with_snp_support",
            ]
        )
    else:
        enriched = audit.assign(
            _uniquely_resolved=audit["canonical_gene_id"].notna() & ~audit["disposition"].eq("rejected"),
            _rejected=audit["disposition"].eq("rejected"),
            _duplicate=audit["disposition"].eq("duplicate"),
            _excluded=audit["disposition"].eq("excluded"),
            _unsupported=audit["disposition"].eq("unsupported"),
            _supported=audit["disposition"].eq("retained"),
        )
        enriched["_resolved_gene_id"] = enriched["canonical_gene_id"].where(
            enriched["_uniquely_resolved"]
        )
        counts = enriched.groupby(key, sort=False).agg(
            nonblank_input_rows=("line", "size"),
            uniquely_resolved_rows=("_uniquely_resolved", "sum"),
            rejected_rows=("_rejected", "sum"),
            unique_resolved_genes=("_resolved_gene_id", "nunique"),
            duplicate_rows=("_duplicate", "sum"),
            excluded_genes=("_excluded", "sum"),
            zero_support_genes=("_unsupported", "sum"),
            genes_with_snp_support=("_supported", "sum"),
        ).reset_index()
    summary = summary.merge(counts, on=key, how="left", sort=False)
    readable = summary["source_status"].eq("ok")
    count_columns = [
        "nonblank_input_rows",
        "uniquely_resolved_rows",
        "rejected_rows",
        "unique_resolved_genes",
        "duplicate_rows",
        "excluded_genes",
    ]
    for column in count_columns:
        summary[column] = pd.to_numeric(summary[column], errors="coerce").astype("Int64")
        summary.loc[readable & summary[column].isna(), column] = 0
        summary.loc[~readable, column] = pd.NA
    if support_evaluated:
        for column in ("zero_support_genes", "genes_with_snp_support"):
            summary[column] = pd.to_numeric(summary[column], errors="coerce").astype("Int64")
            summary.loc[readable & summary[column].isna(), column] = 0
            summary.loc[~readable, column] = pd.NA
    else:
        summary["zero_support_genes"] = pd.NA
        summary["genes_with_snp_support"] = pd.NA
    summary["resolution_fraction"] = pd.NA
    nonempty = readable & summary["nonblank_input_rows"].gt(0).fillna(False)
    summary.loc[nonempty, "resolution_fraction"] = (
        summary.loc[nonempty, "uniquely_resolved_rows"].astype(float)
        / summary.loc[nonempty, "nonblank_input_rows"].astype(float)
    )
    summary["_role_order"] = summary["input_role"].map({"focal": 0, "control": 1})
    return summary.sort_values(["_role_order", "source_ordinal"], kind="stable").drop(columns="_role_order").loc[:, SUMMARY_COLUMNS]


def _mhc_overlap_mask(frame: pd.DataFrame, genome_build: str) -> pd.Series:
    regions = load_preset_intervals(("mhc",), genome_build).intervals
    result = pd.Series(False, index=frame.index)
    for chrom, intervals in regions.items():
        same_chrom = frame["chrom"].astype(str).eq(str(chrom))
        if not same_chrom.any():
            continue
        start0 = pd.to_numeric(frame.loc[same_chrom, "start"], errors="coerce") - 1
        end = pd.to_numeric(frame.loc[same_chrom, "end"], errors="coerce")
        overlap = np.zeros(int(same_chrom.sum()), dtype=bool)
        for region_start, region_end in intervals:
            overlap |= (start0.to_numpy() < int(region_end)) & (
                end.to_numpy() > int(region_start)
            )
        result.loc[same_chrom] = overlap
    return result


def _single_catalog_issue(
    source: str,
    *,
    reason: str,
    field: str,
    observed_value: str,
    repair: str,
    catalog_line: int | None = None,
    details: str = "",
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            _catalog_issue_dict(
                source,
                catalog_line=catalog_line,
                field=field,
                reason=reason,
                observed_value=observed_value,
                details=details,
                repair=repair,
            ),
            _catalog_issue_dict(
                source,
                field="validation",
                reason="validation_incomplete",
                observed_value=reason,
                details="row and global canonicality checks could not run against this catalog",
                repair="Repair the structural catalog problem, then rerun validation to reveal any remaining defects.",
            ),
        ],
        columns=CATALOG_ISSUE_COLUMNS,
    )


def _catalog_issue_dict(
    source: str,
    *,
    catalog_line: Any = pd.NA,
    row: pd.Series | None = None,
    field: str,
    reason: str,
    related_catalog_lines: str = "",
    observed_value: str = "",
    details: str = "",
    repair: str,
) -> dict[str, Any]:
    row = pd.Series(dtype=object) if row is None else row
    return {
        "source": source,
        "catalog_line": catalog_line,
        "gene_id": row.get("gene_id", ""),
        "gene_name": row.get("gene_name", ""),
        "chrom": row.get("chrom", ""),
        "start": row.get("start", ""),
        "end": row.get("end", ""),
        "genome_build": row.get("genome_build", ""),
        "field": field,
        "reason": reason,
        "related_catalog_lines": related_catalog_lines,
        "observed_value": observed_value,
        "details": details,
        "repair": repair,
    }


def _sort_catalog_issues(issues: pd.DataFrame) -> pd.DataFrame:
    if issues.empty:
        return pd.DataFrame(columns=CATALOG_ISSUE_COLUMNS)
    ordered = issues.loc[:, CATALOG_ISSUE_COLUMNS].copy()
    ordered["_line"] = pd.to_numeric(ordered["catalog_line"], errors="coerce").fillna(-1)
    return ordered.sort_values(["_line", "field", "reason"], kind="stable").drop(columns="_line").reset_index(drop=True)


def _combine_reasons(left: str, right: str) -> str:
    return ",".join(sorted(set(filter(None, [*left.split(","), right]))))
