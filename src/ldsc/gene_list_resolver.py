"""Resolve gene-list query annotations against the packaged gene catalog.

This workflow-layer module owns catalog validation and identifier resolution. It
does not log or write diagnostics; LD-score orchestration classifies returned
records and owns user-visible artifacts.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
import gzip
import hashlib
from importlib import resources
from io import BytesIO
from pathlib import Path
import re
from typing import Any

import pandas as pd

from .chromosome_inference import normalize_chromosome
from .errors import LDSCInputError


CATALOG_RESOURCE = "protein_coding_genes.tsv.gz"
CATALOG_RELEASE = "GENCODE v49"
CATALOG_COLUMNS = (
    "ensgid",
    "gene_name",
    "hg38_chr",
    "hg38_start0",
    "hg38_end",
    "hg38_strand",
    "hg19_chr",
    "hg19_start0",
    "hg19_end",
    "hg19_strand",
)
UNRESOLVED_REASONS = (
    "unmatched_identifier",
    "invalid_identifier",
    "build_missing",
    "ambiguous_identifier",
    "malformed_input",
)
_VERSIONED_ENSEMBL = re.compile(r"^(ENSG[0-9]+)\.([0-9]+)$")
_ENSEMBL_LIKE = re.compile(r"^ENSG[0-9]+(?:\..*)?$")


@dataclass(frozen=True)
class GeneUnresolvedRecord:
    """One unresolved or structurally invalid gene-list row."""

    query: str
    source: str
    line: int
    input_gene: str
    reason: str
    canonical_ensembl_id: str | None = None
    details: str | None = None

    def as_dict(self) -> dict[str, Any]:
        """Return the fixed unresolved-audit row schema."""
        return {
            "query": self.query,
            "source": self.source,
            "line": self.line,
            "input_gene": self.input_gene,
            "reason": self.reason,
            "canonical_ensembl_id": self.canonical_ensembl_id,
            "details": self.details,
        }


@dataclass(frozen=True)
class GeneListResolution:
    """Resolved in-memory representation of one gene-list query source."""

    source_ordinal: int
    query: str
    source: str
    source_path: str
    input_sha256: str | None
    status: str
    reason: str
    details: str | None
    canonical_ensembl_ids: tuple[str, ...]
    catalog_indices: tuple[int, ...]
    intervals: tuple[tuple[str, int, int], ...]
    counts: dict[str, int]
    unresolved: tuple[GeneUnresolvedRecord, ...]

    def provenance(self) -> dict[str, Any]:
        """Return compact, path-safe metadata for a successful or attempted query."""
        return {
            "source_ordinal": self.source_ordinal,
            "query": self.query,
            "source": self.source,
            "input_sha256": self.input_sha256,
            "counts": dict(self.counts),
        }


@dataclass(frozen=True)
class GeneCatalog:
    """Validated protein-coding gene catalog and exact identifier maps."""

    frame: pd.DataFrame
    resource: str
    release: str
    content_sha256: str
    ensembl_to_index: dict[str, int]
    gene_name_to_indices: dict[str, tuple[int, ...]]

    @classmethod
    def load(cls, path: str | Path | None = None) -> "GeneCatalog":
        """Load and validate the packaged catalog or an explicit test catalog."""
        if path is None:
            resource = resources.files("ldsc").joinpath("data", CATALOG_RESOURCE)
            compressed = resource.read_bytes()
            resource_name = CATALOG_RESOURCE
        else:
            source_path = Path(path)
            try:
                compressed = source_path.read_bytes()
            except OSError as exc:
                raise LDSCInputError(f"Could not read protein-coding gene catalog '{source_path}': {exc}") from exc
            resource_name = source_path.name
        try:
            raw = gzip.decompress(compressed)
        except (OSError, EOFError) as exc:
            raise LDSCInputError(
                f"Protein-coding gene catalog '{resource_name}' is not a readable gzip TSV: {exc}"
            ) from exc
        try:
            frame = pd.read_csv(BytesIO(raw), sep="\t", dtype=str, keep_default_na=True)
        except Exception as exc:
            raise LDSCInputError(
                f"Protein-coding gene catalog '{resource_name}' could not be parsed as TSV: {exc}"
            ) from exc
        _validate_catalog_frame(frame, resource_name)
        frame = _normalize_catalog_frame(frame)
        ensembl_to_index = {str(value): int(index) for index, value in frame["ensgid"].items()}
        name_groups: dict[str, list[int]] = defaultdict(list)
        for index, value in frame["gene_name"].items():
            name_groups[str(value)].append(int(index))
        return cls(
            frame=frame,
            resource=resource_name,
            release=CATALOG_RELEASE,
            content_sha256=hashlib.sha256(raw).hexdigest(),
            ensembl_to_index=ensembl_to_index,
            gene_name_to_indices={name: tuple(indices) for name, indices in name_groups.items()},
        )

    def provenance(self, genome_build: str) -> dict[str, str]:
        """Return concise catalog provenance for LD-score root metadata."""
        return {
            "resource": self.resource,
            "release": self.release,
            "genome_build": genome_build,
            "content_sha256": self.content_sha256,
        }


def gene_list_query_name(path: str | Path) -> str:
    """Derive a query name from one gene-list basename."""
    name = Path(path).name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for suffix in (".txt", ".tsv", ".list"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return name


def resolve_gene_list(
    path: str | Path,
    catalog: GeneCatalog,
    *,
    genome_build: str,
    source_ordinal: int = 1,
) -> GeneListResolution:
    """Resolve one gene list to unique canonical genes and build intervals."""
    if genome_build not in {"hg19", "hg38"}:
        raise LDSCInputError(
            f"Gene-list resolution requires a concrete catalog build, got {genome_build!r}. "
            "Pass or infer hg19/hg38 before resolving gene lists."
        )
    source_path = Path(path)
    query = gene_list_query_name(source_path)
    source = source_path.name
    try:
        payload = source_path.read_bytes()
        text = gzip.decompress(payload).decode("utf-8") if source.lower().endswith(".gz") else payload.decode("utf-8")
    except (OSError, EOFError, UnicodeError) as exc:
        return _unreadable_resolution(source_path, query, source_ordinal, exc)

    unresolved: list[GeneUnresolvedRecord] = []
    normalized_tokens: list[str] = []
    matched_rows = 0
    resolved_rows: list[tuple[str, int, str]] = []
    problem_counts: Counter[str] = Counter()
    blank_rows = 0
    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        stripped = raw_line.strip()
        if not stripped:
            blank_rows += 1
            continue
        if len(raw_line.split("\t")) != 1:
            problem_counts["malformed_input"] += 1
            normalized_tokens.append(stripped)
            unresolved.append(
                GeneUnresolvedRecord(
                    query,
                    source,
                    line_number,
                    stripped,
                    "malformed_input",
                    details="expected one tab-delimited field",
                )
            )
            continue

        token = stripped
        version_match = _VERSIONED_ENSEMBL.fullmatch(token)
        normalized = version_match.group(1) if version_match else token
        normalized_tokens.append(normalized)
        if _ENSEMBL_LIKE.fullmatch(token) and "." in token and version_match is None:
            problem_counts["invalid_identifier"] += 1
            unresolved.append(
                GeneUnresolvedRecord(
                    query,
                    source,
                    line_number,
                    token,
                    "invalid_identifier",
                    details="Ensembl version suffix must be numeric",
                )
            )
            continue

        id_index = catalog.ensembl_to_index.get(normalized)
        name_indices = catalog.gene_name_to_indices.get(token, ())
        candidates = set(name_indices)
        if id_index is not None:
            candidates.add(id_index)
        if len(candidates) > 1:
            conflict_ids = sorted(catalog.frame.loc[list(candidates), "ensgid"].astype(str))
            problem_counts["ambiguous_identifier"] += 1
            unresolved.append(
                GeneUnresolvedRecord(
                    query,
                    source,
                    line_number,
                    token,
                    "ambiguous_identifier",
                    details="conflicting Ensembl IDs: " + ", ".join(conflict_ids),
                )
            )
            continue
        if not candidates:
            problem_counts["unmatched_identifier"] += 1
            unresolved.append(
                GeneUnresolvedRecord(query, source, line_number, token, "unmatched_identifier")
            )
            continue

        index = next(iter(candidates))
        canonical_id = str(catalog.frame.at[index, "ensgid"])
        matched_rows += 1
        chrom = catalog.frame.at[index, f"{genome_build}_chr"]
        start = catalog.frame.at[index, f"{genome_build}_start0"]
        end = catalog.frame.at[index, f"{genome_build}_end"]
        if pd.isna(chrom):
            problem_counts["build_missing"] += 1
            unresolved.append(
                GeneUnresolvedRecord(
                    query,
                    source,
                    line_number,
                    token,
                    "build_missing",
                    canonical_ensembl_id=canonical_id,
                    details=f"no {genome_build} coordinates in catalog",
                )
            )
            continue
        resolved_rows.append((canonical_id, int(index), normalized))

    canonical_ids: list[str] = []
    catalog_indices: list[int] = []
    seen_ids: set[str] = set()
    unique_matched_tokens: set[str] = set()
    for canonical_id, index, normalized in resolved_rows:
        unique_matched_tokens.add(normalized)
        if canonical_id in seen_ids:
            continue
        seen_ids.add(canonical_id)
        canonical_ids.append(canonical_id)
        catalog_indices.append(index)

    intervals = tuple(
        (
            str(catalog.frame.at[index, f"{genome_build}_chr"]),
            int(catalog.frame.at[index, f"{genome_build}_start0"]),
            int(catalog.frame.at[index, f"{genome_build}_end"]),
        )
        for index in catalog_indices
    )
    token_counts = Counter(normalized_tokens)
    counts = {
        "nonblank_input_rows": len(normalized_tokens),
        "unique_normalized_input_tokens": len(token_counts),
        "repeated_token_rows": sum(count - 1 for count in token_counts.values()),
        "matched_input_rows": matched_rows,
        "unique_resolved_canonical_genes": len(canonical_ids),
        "alias_collapsed_rows": max(0, len(unique_matched_tokens) - len(canonical_ids)),
        "blank_rows": blank_rows,
        **{reason: int(problem_counts[reason]) for reason in UNRESOLVED_REASONS},
    }
    if problem_counts["malformed_input"]:
        status, reason = "skipped", "malformed_input"
    elif problem_counts["ambiguous_identifier"]:
        status, reason = "skipped", "ambiguous_identifier"
    elif not normalized_tokens:
        status, reason = "skipped", "empty_input"
    elif not canonical_ids:
        status, reason = "skipped", "fully_unresolved"
    elif any(problem_counts[item] for item in ("unmatched_identifier", "invalid_identifier", "build_missing")):
        status, reason = "warning", "partial_resolution"
    else:
        status, reason = "ok", ""
    details = "diagnostics/gene_list_unresolved.tsv.gz" if unresolved else None
    return GeneListResolution(
        source_ordinal=source_ordinal,
        query=query,
        source=source,
        source_path=str(source_path),
        input_sha256=hashlib.sha256(payload).hexdigest(),
        status=status,
        reason=reason,
        details=details,
        canonical_ensembl_ids=tuple(canonical_ids),
        catalog_indices=tuple(catalog_indices),
        intervals=intervals,
        counts=counts,
        unresolved=tuple(unresolved),
    )


def _unreadable_resolution(
    path: Path,
    query: str,
    source_ordinal: int,
    error: Exception,
) -> GeneListResolution:
    """Return the query-local unreadable-source outcome."""
    counts = {
        "nonblank_input_rows": 0,
        "unique_normalized_input_tokens": 0,
        "repeated_token_rows": 0,
        "matched_input_rows": 0,
        "unique_resolved_canonical_genes": 0,
        "alias_collapsed_rows": 0,
        "blank_rows": 0,
        **{reason: 0 for reason in UNRESOLVED_REASONS},
    }
    return GeneListResolution(
        source_ordinal=source_ordinal,
        query=query,
        source=path.name,
        source_path=str(path),
        input_sha256=None,
        status="skipped",
        reason="unreadable_source",
        details=str(error).replace(str(path), path.name),
        canonical_ensembl_ids=(),
        catalog_indices=(),
        intervals=(),
        counts=counts,
        unresolved=(),
    )


def _validate_catalog_frame(frame: pd.DataFrame, resource: str) -> None:
    """Validate catalog columns and identity fields before normalization."""
    if tuple(frame.columns) != CATALOG_COLUMNS:
        raise LDSCInputError(
            f"Protein-coding gene catalog '{resource}' has incompatible columns. Expected "
            f"{list(CATALOG_COLUMNS)}, found {list(frame.columns)}."
        )
    if frame["ensgid"].isna().any() or frame["ensgid"].astype(str).str.strip().eq("").any():
        raise LDSCInputError(f"Protein-coding gene catalog '{resource}' contains an empty Ensembl gene ID.")
    if frame["ensgid"].duplicated().any():
        duplicate = str(frame.loc[frame["ensgid"].duplicated(keep=False), "ensgid"].iloc[0])
        raise LDSCInputError(
            f"Protein-coding gene catalog '{resource}' contains duplicate Ensembl gene ID {duplicate!r}."
        )
    if frame["gene_name"].isna().any() or frame["gene_name"].astype(str).str.strip().eq("").any():
        raise LDSCInputError(f"Protein-coding gene catalog '{resource}' contains an empty gene name.")
    for build in ("hg19", "hg38"):
        coordinate_columns = [f"{build}_chr", f"{build}_start0", f"{build}_end"]
        missing = frame[coordinate_columns].isna()
        partial = missing.any(axis=1) & ~missing.all(axis=1)
        if partial.any():
            row = int(partial[partial].index[0]) + 2
            raise LDSCInputError(
                f"Protein-coding gene catalog '{resource}' coordinate triplet for {build} is partially missing at TSV line {row}."
            )


def _normalize_catalog_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Normalize valid chromosome and coordinate fields to stable types."""
    normalized = frame.copy()
    normalized["ensgid"] = normalized["ensgid"].astype(str).str.strip()
    normalized["gene_name"] = normalized["gene_name"].astype(str).str.strip()
    for build in ("hg19", "hg38"):
        chr_column = f"{build}_chr"
        start_column = f"{build}_start0"
        end_column = f"{build}_end"
        present = normalized[chr_column].notna()
        for index in normalized.index[present]:
            try:
                chrom = normalize_chromosome(normalized.at[index, chr_column])
                raw_start = float(normalized.at[index, start_column])
                raw_end = float(normalized.at[index, end_column])
                if not raw_start.is_integer() or not raw_end.is_integer():
                    raise ValueError("coordinates must be integers")
                start = int(raw_start)
                end = int(raw_end)
            except (TypeError, ValueError) as exc:
                raise LDSCInputError(
                    f"Protein-coding gene catalog has invalid {build} coordinates for "
                    f"{normalized.at[index, 'ensgid']!r}."
                ) from exc
            if chrom not in {str(value) for value in range(1, 23)} or start < 0 or end <= start:
                raise LDSCInputError(
                    f"Protein-coding gene catalog has invalid {build} interval for "
                    f"{normalized.at[index, 'ensgid']!r}: chromosome={chrom!r}, start0={start}, end={end}."
                )
            normalized.at[index, chr_column] = chrom
            normalized.at[index, start_column] = start
            normalized.at[index, end_column] = end
        normalized[start_column] = pd.to_numeric(normalized[start_column], errors="coerce").astype("Int64")
        normalized[end_column] = pd.to_numeric(normalized[end_column], errors="coerce").astype("Int64")
    return normalized.reset_index(drop=True)
