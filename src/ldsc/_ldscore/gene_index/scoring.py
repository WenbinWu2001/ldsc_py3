"""Loading a gene LD-score index and scoring gene lists with it.

``load_gene_ldscore_index`` opens and verifies a published index;
``run_indexed_ldscore`` resolves query and control gene lists against the
index's catalog and writes their exact LD scores without re-reading the
reference panel.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Sequence

import numpy as np
import pandas as pd

from ..._annotation.gene_lists import GeneCatalog, GeneSourceSelection, select_index_eligible_gene_indices
from ..._annotation.storage import AnnotationWorkspace
from ..._core.snp_identity import normalize_snp_identifier_mode
from ...config import GlobalConfig
from ..outputs import common_count_config
from .artifact import (
    IndexChromosomeData,
    _INDEX_FORMAT_VERSION,
    _catalog_identity_records,
    _load_index_chromosome,
    _read_json,
    calculate_index_id,
    index_global_config,
    inspect_index_paths,
)
from .atoms import ChromosomeAtomModel, assemble_atom_selector
from ..._runtime.parallelism import validate_workers
from ..._runtime.paths import normalize_path_token
from ..._refpanel.funnel import concat_funnels
from ..._runtime.progress import advance, report_phase
from ...errors import LDSCInputError

LOGGER = logging.getLogger("LDSC.gene_ldscore_index")


@dataclass(frozen=True)
class LoadedGeneLDScoreIndex:
    """One validated, self-contained gene LD-score index.

    Attributes
    ----------
    index_id : str
        Content-bound scientific identity for the complete index.
    chromosomes : tuple of str
        Ordered immutable chromosome coverage.
    index_identity : dict
        Scientific settings and canonical input-content identities used to
        calculate ``index_id``.
    gene_catalog : pandas.DataFrame
        Embedded catalog with unpadded zero-based half-open gene coordinates,
        inclusion flags, and chromosome-local row ordering.
    chromosome_paths : dict of str to pathlib.Path
        Validated component locations. No chromosome payload is cached.
    gene_support : pandas.Series
        Compact retained-reference support counts collected during validation.
    reference_snp_funnel : pandas.DataFrame
        Reference-SNP funnel of every chromosome, recorded when the index was
        built; indexed scoring reports it rather than recomputing it.
    """

    index_id: str
    chromosomes: tuple[str, ...]
    snp_identifier: str
    genome_build: str
    index_identity: dict
    gene_catalog: pd.DataFrame
    chromosome_paths: dict[str, Path]
    baseline_columns: tuple[str, ...]
    gene_support: pd.Series
    reference_snp_funnel: pd.DataFrame

    @property
    def global_config(self) -> GlobalConfig:
        """Identity settings the index binds; indexed output records them."""
        return index_global_config(self.snp_identifier, self.genome_build)

    @property
    def count_config(self) -> dict:
        """Common-SNP count metadata for the ``common_maf_min`` the index used."""
        return common_count_config(self.index_identity["common_maf_min"])

    def load_chromosome(self, chrom: str) -> IndexChromosomeData:
        """Load one owned chromosome payload; callers release it after use.

        The immutable component is checked again when opened. Indexed scoring
        reloads it for each active query batch and releases it before the next
        chromosome task; operators are not cached between tasks.
        """
        path = self.chromosome_paths[str(chrom)]
        return _load_index_chromosome(
            path.parent.parent, str(chrom), index_id=self.index_id,
            snp_identifier=self.snp_identifier, genome_build=self.genome_build,
            expected_gene_rows=int(self.gene_catalog.chrom.astype(str).eq(str(chrom)).sum()),
        )


def load_gene_ldscore_index(
    index_dir: str | Path,
) -> LoadedGeneLDScoreIndex:
    """Load one complete public autosomes 1--22 gene LD-score index."""
    return _load_gene_ldscore_index(index_dir, _allow_partial_for_tests=False)


@report_phase(LOGGER, 'validation', 'index content and compatibility')
def _load_gene_ldscore_index(
    index_dir: str | Path,
    *,
    _allow_partial_for_tests: bool,
) -> LoadedGeneLDScoreIndex:
    """Load and fully validate one explicit gene LD-score index directory.

    Parameters
    ----------
    index_dir : path-like
        Directory containing one complete format-version-4 index covering
        autosomes 1 through 22. Earlier formats require a complete rebuild;
        no profile discovery or compatibility migration is performed.

    Returns
    -------
    LoadedGeneLDScoreIndex
        Index components after validating semantic identity,
        chromosome coverage, catalog and row ordering, shapes, dtypes, finite
        values, and CSR invariants.

    Raises
    ------
    LDSCInputError
        If any required component is missing, corrupt, structurally invalid,
        or semantically incompatible.
    """
    index_path = Path(index_dir)
    root = _read_json(index_path / "metadata.json", "root")
    if root.get("format_version") != _INDEX_FORMAT_VERSION:
        raise LDSCInputError(
            f"Unsupported gene LD-score index format version {root.get('format_version')!r} "
            f"(this package reads version {_INDEX_FORMAT_VERSION}). Rebuild the complete index with the "
            "current package; earlier indexes cannot be reused or migrated. Version 4 stores each "
            "chromosome's reference-SNP funnel, which indexed scoring reports; version 3 made PLINK LD "
            "windows exact, and earlier indexes' LD scores differ slightly from current direct ldscore runs."
        )
    if root.get("artifact_type") != "gene_ldscore_index":
        raise LDSCInputError(
            "Gene LD-score index root has the wrong artifact_type."
            " Restore the complete saved index from one compatible build, or regenerate it with ldsc "
            "build-gene-ldscore-index. See "
            "docs/user-guide/troubleshooting.md#ldscore-an-explicit-gene-ld-score-index-is-missing-corrupt-or-incompatible."
        )
    index_identity = root.get("index_identity")
    index_id = root.get("index_id")
    if not isinstance(index_identity, dict) or index_id != calculate_index_id(index_identity):
        raise LDSCInputError(
            "Gene LD-score index semantic identity is invalid."
            " Restore the complete saved index from one compatible build, or regenerate it with ldsc "
            "build-gene-ldscore-index. See "
            "docs/user-guide/troubleshooting.md#ldscore-an-explicit-gene-ld-score-index-is-missing-corrupt-or-incompatible."
        )
    try:
        snp_identifier = normalize_snp_identifier_mode(str(index_identity.get("snp_identifier")))
    except Exception as exc:
        raise LDSCInputError(
            "Gene LD-score index SNP identity metadata is invalid."
            " Restore the complete saved index from one compatible build, or regenerate it with ldsc "
            "build-gene-ldscore-index. See "
            "docs/user-guide/troubleshooting.md#ldscore-an-explicit-gene-ld-score-index-is-missing-corrupt-or-incompatible."
        ) from exc
    genome_build = index_identity.get("genome_build")
    if genome_build not in {"hg19", "hg38"}:
        raise LDSCInputError(
            "Gene LD-score index genome-build metadata must be explicit hg19 or hg38."
            " Restore the complete saved index from one compatible build, or regenerate it with ldsc "
            "build-gene-ldscore-index. See "
            "docs/user-guide/troubleshooting.md#ldscore-an-explicit-gene-ld-score-index-is-missing-corrupt-or-incompatible."
        )
    if root.get("snp_identifier") != snp_identifier or root.get("genome_build") != genome_build:
        raise LDSCInputError(
            "Gene LD-score index root identity metadata disagrees with index_identity."
            " Restore the complete saved index from one compatible build, or regenerate it with ldsc "
            "build-gene-ldscore-index. See "
            "docs/user-guide/troubleshooting.md#ldscore-an-explicit-gene-ld-score-index-is-missing-corrupt-or-incompatible."
        )
    chromosomes = tuple(str(chrom) for chrom in root.get("chromosomes", ()))
    expected_coverage = tuple(str(chrom) for chrom in index_identity.get("chromosomes", ()))
    if not chromosomes or chromosomes != expected_coverage:
        raise LDSCInputError(
            "Gene LD-score index chromosome coverage is invalid."
            " Restore the complete saved index from one compatible build, or regenerate it with ldsc "
            "build-gene-ldscore-index. See "
            "docs/user-guide/troubleshooting.md#ldscore-an-explicit-gene-ld-score-index-is-missing-corrupt-or-incompatible."
        )
    full_coverage = tuple(str(value) for value in range(1, 23))
    if not _allow_partial_for_tests and chromosomes != full_coverage:
        raise LDSCInputError(
            "Gene LD-score index must cover autosomes 1 through 22; partial production indexes are unsupported."
            " Restore the complete saved index from one compatible build, or regenerate it with ldsc "
            "build-gene-ldscore-index. See "
            "docs/user-guide/troubleshooting.md#ldscore-an-explicit-gene-ld-score-index-is-missing-corrupt-or-incompatible."
        )
    advance(0, total=len(chromosomes))
    inspect_index_paths(index_path, metadata=root)
    catalog_path = index_path / "gene_catalog.parquet"
    catalog = pd.read_parquet(catalog_path)
    required_catalog = [
        "gene_index", "gene_id", "gene_name", "chrom", "start", "end", "genome_build",
        "catalog_line", "source", "included", "exclusion_reason", "chromosome_gene_row",
    ]
    if list(catalog.columns) != required_catalog:
        raise LDSCInputError(
            "Gene LD-score index catalog schema is invalid."
            " Restore the complete saved index from one compatible build, or regenerate it with ldsc "
            "build-gene-ldscore-index. See "
            "docs/user-guide/troubleshooting.md#ldscore-an-explicit-gene-ld-score-index-is-missing-corrupt-or-incompatible."
        )
    expected_gene_indices = np.arange(len(catalog), dtype=np.int64)
    if not np.array_equal(catalog["gene_index"].to_numpy(), expected_gene_indices):
        raise LDSCInputError(
            "Gene LD-score index catalog gene_index ordering is invalid."
            " Restore the complete saved index from one compatible build, or regenerate it with ldsc "
            "build-gene-ldscore-index. See "
            "docs/user-guide/troubleshooting.md#ldscore-an-explicit-gene-ld-score-index-is-missing-corrupt-or-incompatible."
        )
    catalog_authority = GeneCatalog.from_embedded_frame(catalog)
    if index_identity.get("catalog") != _catalog_identity_records(catalog):
        raise LDSCInputError(
            "Gene LD-score index embedded catalog disagrees with its semantic identity."
            " Restore the complete saved index from one compatible build, or regenerate it with ldsc "
            "build-gene-ldscore-index. See "
            "docs/user-guide/troubleshooting.md#ldscore-an-explicit-gene-ld-score-index-is-missing-corrupt-or-incompatible."
        )
    if catalog["included"].dtype != np.bool_:
        raise LDSCInputError(
            "Gene LD-score index catalog included flags must use Boolean dtype."
            " Restore the complete saved index from one compatible build, or regenerate it with ldsc "
            "build-gene-ldscore-index. See "
            "docs/user-guide/troubleshooting.md#ldscore-an-explicit-gene-ld-score-index-is-missing-corrupt-or-incompatible."
        )
    source_values = catalog["source"].fillna("").astype(str)
    if (
        source_values.eq("").any()
        or source_values.nunique() != 1
        or source_values.str.contains(r"[/\\]", regex=True).any()
    ):
        raise LDSCInputError(
            "Gene LD-score index catalog source must be one nonempty source basename."
            " Restore the complete saved index from one compatible build, or regenerate it with ldsc "
            "build-gene-ldscore-index. See "
            "docs/user-guide/troubleshooting.md#ldscore-an-explicit-gene-ld-score-index-is-missing-corrupt-or-incompatible."
        )
    gene_policy = index_identity.get("gene_exclude_regions")
    if gene_policy not in {"none", "mhc"}:
        raise LDSCInputError(
            "Gene LD-score index gene-exclusion policy is invalid."
            " Restore the complete saved index from one compatible build, or regenerate it with ldsc "
            "build-gene-ldscore-index. See "
            "docs/user-guide/troubleshooting.md#ldscore-an-explicit-gene-ld-score-index-is-missing-corrupt-or-incompatible."
        )
    eligible = select_index_eligible_gene_indices(
        catalog_authority,
        gene_exclude_regions=str(gene_policy),
    )
    expected_included = catalog.index.isin(eligible)
    expected_reasons = np.where(expected_included, "", "excluded_gene_region")
    if not np.array_equal(catalog["included"].to_numpy(dtype=bool), expected_included):
        raise LDSCInputError(
            "Gene LD-score index catalog inclusion policy disagrees with its coordinates and metadata."
            " Restore the complete saved index from one compatible build, or regenerate it with ldsc "
            "build-gene-ldscore-index. See "
            "docs/user-guide/troubleshooting.md#ldscore-an-explicit-gene-ld-score-index-is-missing-corrupt-or-incompatible."
        )
    if not np.array_equal(catalog["exclusion_reason"].fillna("").astype(str), expected_reasons):
        raise LDSCInputError(
            "Gene LD-score index catalog exclusion reasons disagree with its inclusion policy."
            " Restore the complete saved index from one compatible build, or regenerate it with ldsc "
            "build-gene-ldscore-index. See "
            "docs/user-guide/troubleshooting.md#ldscore-an-explicit-gene-ld-score-index-is-missing-corrupt-or-incompatible."
        )
    chromosome_paths = {}
    baseline_columns = None
    support = pd.Series(pd.NA, index=catalog["gene_index"].astype(int), dtype="Int64")
    funnels = []
    issues = []
    for chrom in chromosomes:
        advance(0, object=f'chromosome {chrom}', force=True)
        record = None
        try:
            chromosome_catalog = catalog.loc[catalog["chrom"].astype(str) == str(chrom)]
            expected_chromosome_rows = np.arange(len(chromosome_catalog), dtype=np.int64)
            if not np.array_equal(
                chromosome_catalog["chromosome_gene_row"].to_numpy(), expected_chromosome_rows
            ):
                raise LDSCInputError(
                    f"Gene LD-score index catalog chromosome_gene_row ordering is invalid for chromosome {chrom}."
                    " Restore the complete saved index from one compatible build, or regenerate it with ldsc "
                    "build-gene-ldscore-index. See "
                    "docs/user-guide/troubleshooting.md#ldscore-an-explicit-gene-ld-score-index-is-missing-corrupt-or-incompatible."
                )
            record = _load_index_chromosome(
                index_path,
                chrom,
                index_id=index_id,
                snp_identifier=snp_identifier,
                genome_build=genome_build,
                expected_gene_rows=len(chromosome_catalog),
            )
            columns = tuple(_baseline_columns_from_rows(record.baseline_rows))
            if baseline_columns is not None and columns != baseline_columns:
                raise LDSCInputError(
                    "Gene LD-score index baseline columns disagree across chromosomes."
                    " Restore the complete saved index from one compatible build, or regenerate it with ldsc "
                    "build-gene-ldscore-index. See "
                    "docs/user-guide/troubleshooting.md#ldscore-an-explicit-gene-ld-score-index-is-missing-corrupt-or-incompatible."
                )
            baseline_columns = columns
            counts = record.atom_model.gene_to_atom.astype(np.int64) @ record.atom_statistics.atom_count_all
            support.loc[chromosome_catalog["gene_index"].astype(int)] = np.asarray(counts).reshape(-1)
            funnels.append(record.reference_snp_funnel)
            chromosome_paths[chrom] = index_path / "chromosomes" / f"chr{chrom}"
            advance(object=f'chromosome {chrom}')
        except (OSError, ValueError, EOFError, LDSCInputError) as exc:
            issues.append({"input_role": "index", "source": str(index_path / "chromosomes" / f"chr{chrom}"),
                           "chrom": chrom, "reason": "invalid_required_input", "details": str(exc),
                           "repair": "Restore or rebuild the complete immutable index."})
        finally:
            record = None
    if issues:
        error = LDSCInputError("Gene LD-score index chromosome preflight failed: " + "; ".join(
            f"chromosome {row['chrom']} {row['source']}: {row['details']}" for row in issues
        ))
        error.input_issues = pd.DataFrame(issues)
        raise error
    return LoadedGeneLDScoreIndex(
        index_id=index_id,
        chromosomes=chromosomes,
        snp_identifier=snp_identifier,
        genome_build=genome_build,
        index_identity=index_identity,
        gene_catalog=catalog,
        chromosome_paths=chromosome_paths,
        baseline_columns=baseline_columns,
        gene_support=support,
        reference_snp_funnel=concat_funnels(funnels),
    )


def run_indexed_ldscore(
    index_dir: str | Path,
    *,
    query_gene_list_sources: Sequence[str | Path],
    control_gene_list_file: str | Path | None = None,
    gene_list_resolution_policy: str = "strict",
    query_batch_size: int = 100,
    workers: int = 1,
    output_dir: str | Path,
    overwrite: bool = False,
    _allow_partial_for_tests: bool = False,
):
    """Assemble exact gene-list LD scores from one explicit index directory.

    Parameters
    ----------
    index_dir : path-like
        Exact validated single-index directory.
    query_gene_list_sources : sequence of path-like
        One-column gene-list files. Exact authoritative gene IDs and
        case-sensitive gene names are resolved against the embedded catalog in
        source order.
    control_gene_list_file : path-like or None, optional
        Optional fixed-control gene-list file. When omitted, no
        ``gene_control`` annotation is added.
    gene_list_resolution_policy : {"strict", "resolved-only"}, optional
        Default ``"strict"`` stops on rejected identifiers after auditing all
        query gene lists and the optional control. ``"resolved-only"`` omits
        allowlisted rejected rows with diagnostics; invalid indexes, incomplete
        chromosome coverage, and an unusable requested control remain fatal.
        This Python parameter retains its string values. In ``ldsc ldscore``,
        omit ``--allow-unresolved-genes`` for strict resolution or supply that
        value-free flag for ``"resolved-only"``.
    query_batch_size : int, optional
        Positive maximum query columns per operator multiplication, default 100.
        Bounds outstanding chromosome fragments and saved query-file width.
        One batch is consumed before the next is computed; the parent writes
        one chromosome at a time without a genome-wide query array.
        Chromosome operators are reloaded for each batch in a persistent pool.
        Each pathway remains a separate annotation and regression model.
    workers : int, optional
        Maximum concurrent workers, default 1, including on HPC.
        Use joblib.cpu_count() as the detected available CPU budget.
        Positive requests above that budget fail before capping by task
        count; negative -k uses max(1, budget + 1 - k). Zero, booleans,
        and non-integers are invalid. Even one worker runs in its own process.
        Every worker uses 1 BLAS thread in serial and parallel execution.
        Numerical OpenMP and Arrow CPU/I/O limits are also one; unsupported
        numerical limits are logged as unverified and computation continues.
        This does not detect idle CPUs, guarantee scheduler allocation
        compliance, or ensure sufficient memory.
    output_dir : path-like
        Destination for the canonical self-contained LD-score directory.
        User and environment tokens are expanded for both scratch and results.
    overwrite : bool, optional
        Replace workflow-owned output artifacts. Default is ``False``.

    Returns
    -------
    LDScoreSource
        Shared baseline/control values, counts, overlap, statuses, provenance,
        and saved paths. ``read_queries(names)`` reads explicit selections
        across the ordered query batch files; no query values are cached.

    Raises
    ------
    LDSCInputError
        If the index is invalid, a control is unusable, or every requested
        query is skipped. The all-skipped case writes diagnostics but no
        scientific LD-score tables.
    ValueError
        Query batch width is not a positive integer or the worker request is
        not a nonzero integer.
    FileExistsError
        An owned output artifact exists and overwrite is disabled.

    Notes
    -----
    Overlapping, nested, duplicated, and alias-selected genes are combined by
    Boolean union. A persistent worker pool processes one active query batch;
    each chromosome task reloads its operator and releases it after projection.
    Private fragments preserve float64 through global query-status checks.
    Assembly reads one chromosome at a time in index order, narrows to the
    public float32 format, and removes consumed fragments before the next batch.
    Completed batches remain private until the complete run succeeds.
    One query file is named ``ldscore.query.parquet``; multiple files use
    numbered names in the required ``query_batches`` manifest. This API writes
    outputs even for a single batch; batch resumption is not supported.
    """
    if isinstance(query_batch_size, bool) or not isinstance(query_batch_size, int) or query_batch_size < 1:
        raise ValueError("query_batch_size must be a positive integer.")
    validate_workers(workers)
    output_dir = normalize_path_token(output_dir)
    from ..._runtime.preflight import inspect_declared_inputs
    from ..outputs import LDScoreDirectoryWriter
    LDScoreDirectoryWriter.artifact_family(output_dir).preflight(overwrite=overwrite)
    try:
        inspect_declared_inputs(files=[('query genes', query_gene_list_sources), ('control genes', control_gene_list_file)],
            checks=[('index', index_dir, lambda: inspect_index_paths(index_dir))],
            issues_path=Path(output_dir)/'diagnostics/input_issues.tsv')
    except LDSCInputError as exc:
        diagnostic = SimpleNamespace(input_issues=exc.input_issues, chromosome_scope={
            'selection': 'validated_immutable_index', 'chromosomes': None,
            'analysis_chromosomes': [], 'validation_status': 'failed'})
        LDScoreDirectoryWriter._write_query_diagnostic_files(diagnostic, {
            'input_issues': Path(output_dir)/'diagnostics/input_issues.tsv',
            'chromosome_scope': Path(output_dir)/'diagnostics/chromosome_scope.json'})
        raise
    with AnnotationWorkspace(output_dir) as workspace:
        return _run_indexed_ldscore_in_workspace(
            index_dir, query_gene_list_sources=query_gene_list_sources,
            control_gene_list_file=control_gene_list_file,
            gene_list_resolution_policy=gene_list_resolution_policy,
            query_batch_size=query_batch_size, workers=workers, output_dir=output_dir,
            overwrite=overwrite, workspace=workspace,
            _allow_partial_for_tests=_allow_partial_for_tests,
        )


def _run_indexed_ldscore_in_workspace(
    index_dir, *, query_gene_list_sources, control_gene_list_file,
    gene_list_resolution_policy, query_batch_size, workers, output_dir, overwrite,
    workspace, _allow_partial_for_tests,
):
    """Own staged gene selections until the canonical result has been written."""
    from ..._annotation.gene_query_storage import resolve_gene_lists_staged
    from ..._annotation.query_outcomes import (
    assess_gene_coverage,
    gene_query_statuses,
    gene_control_errors,
    gene_viability_errors,
)
    from ..outputs import LDScoreDirectoryWriter, LDScoreOutputConfig

    try:
        index = _load_gene_ldscore_index(index_dir, _allow_partial_for_tests=_allow_partial_for_tests)
    except (OSError, ValueError, LDSCInputError) as exc:
        from ..preflight import ISSUE_COLUMNS

        issues = getattr(exc, "input_issues", pd.DataFrame([{
            "input_role": "index", "source": str(index_dir), "chrom": "", "reason": "invalid_required_input",
            "details": str(exc), "repair": "Restore or rebuild the complete immutable autosomes 1-22 index.",
        }], columns=ISSUE_COLUMNS))
        diagnostic = SimpleNamespace(input_issues=issues, chromosome_scope={
            "selection": "validated_immutable_index", "chromosomes": None,
            "analysis_chromosomes": [], "validation_status": "failed",
        })
        LDScoreDirectoryWriter().write_query_diagnostics(diagnostic, LDScoreOutputConfig(output_dir=output_dir, overwrite=overwrite))
        raise
    projection_build = str(index.index_identity.get("projection_build"))
    gene_policy = str(index.index_identity.get("gene_exclude_regions", "none"))
    catalog = GeneCatalog.from_embedded_frame(index.gene_catalog)
    if catalog.genome_build != projection_build:
        raise LDSCInputError("Embedded gene coordinate catalog build disagrees with the index projection build.")
    batch = resolve_gene_lists_staged(
        query_gene_list_sources,
        catalog,
        workspace,
        control_path=control_gene_list_file,
        resolution_policy=gene_list_resolution_policy,
        gene_exclude_regions=gene_policy,
    )
    from ..._annotation.query_outcomes import (
    gene_list_gate_a_message,
    log_gene_list_rejections,
    log_gene_list_snp_support,
    log_query_annotation_statuses,
)
    if batch.has_fatal_gate_a_issues:
        LDScoreDirectoryWriter().write_gene_list_preflight(
            batch,
            LDScoreOutputConfig(output_dir=output_dir, overwrite=overwrite),
        )
        log_gene_list_rejections(batch)
        raise LDSCInputError(gene_list_gate_a_message(batch))
    batch, coverage_errors = assess_gene_coverage(batch, index.chromosomes)
    scope = {"chromosomes": list(index.chromosomes), "analysis_chromosomes": list(index.chromosomes),
             "selection": "validated_immutable_index", "index_id": index.index_id,
             "validation_status": "passed"}
    LOGGER.info("Chromosomes resolved and entering the analysis: %s.", ", ".join(index.chromosomes))
    coverage_errors.extend(gene_control_errors(batch))
    if coverage_errors:
        diagnostic = SimpleNamespace(query_statuses=gene_query_statuses(batch), gene_list_batch=batch,
                                     chromosome_scope={**scope, "validation_status": "failed", "analysis_chromosomes": []})
        LDScoreDirectoryWriter().write_query_diagnostics(diagnostic, LDScoreOutputConfig(output_dir=output_dir, overwrite=overwrite))
        raise LDSCInputError("; ".join(coverage_errors))
    batch = batch.with_snp_support(index.gene_support)
    statuses = gene_query_statuses(batch)
    usable = [
        declaration
        for declaration, status in zip(
            (item for item in batch.declarations if item['input_role'] == "focal"),
            statuses,
            strict=True,
        )
        if status.status in {"ok", "warning"}
    ]
    control_resolution = batch.selection("control", 0) if control_gene_list_file is not None else None
    control_errors = gene_viability_errors(batch, statuses)
    if control_errors:
        diagnostic_result = SimpleNamespace(query_statuses=statuses, gene_list_batch=batch, chromosome_scope=scope)
        LDScoreDirectoryWriter().write_query_diagnostics(
            diagnostic_result,
            LDScoreOutputConfig(output_dir=output_dir, overwrite=overwrite),
        )
        log_gene_list_rejections(batch)
        log_gene_list_snp_support(batch)
        log_query_annotation_statuses(statuses)
        raise LDSCInputError("; ".join(control_errors))
    if any(item["query"] == "gene_control" for item in usable):
        raise LDSCInputError(
            "Query name 'gene_control' collides with the reserved fixed control column. "
            "Rename the query gene-list file so its query name differs from gene_control."
        )

    if "gene_control" in index.baseline_columns:
        raise LDSCInputError(
            "Index supplied baseline columns contain reserved name 'gene_control'. "
            "Rebuild the index with a distinct baseline annotation name; gene_control is reserved for the requested gene control."
        )
    from .batches import indexed_results

    from ..diagnostics import ldscore_logs, log_reference_snp_funnel_warning
    with ldscore_logs(output_dir, index.chromosomes, overwrite=overwrite,
                      log_level=logging.getLevelName(LOGGER.getEffectiveLevel())) as log_paths:
        results = indexed_results(index, batch, usable, control_resolution, query_batch_size,
                                  workspace.path, workers, statuses, scope, log_paths)
        result = LDScoreDirectoryWriter().write_batches(
            results, LDScoreOutputConfig(output_dir=output_dir, overwrite=overwrite),
            query_status_order=[status.query for status in statuses])
        log_reference_snp_funnel_warning(index.reference_snp_funnel, result.output_paths["reference_snp_funnel"])
        log_gene_list_rejections(result.gene_list_batch)
        log_gene_list_snp_support(result.gene_list_batch)
        log_query_annotation_statuses(result.query_statuses)
        return result


def _baseline_columns_from_rows(frame: pd.DataFrame) -> list[str]:
    """Return ordered supplied baseline score columns from one common shard."""
    reserved = {"CHR", "SNP", "POS", "A1", "A2", "regression_ld_scores"}
    return [column for column in frame.columns if column not in reserved]


def _selector_for_resolution(
    embedded_catalog: pd.DataFrame,
    atom_model: ChromosomeAtomModel,
    chrom: str,
    resolution: GeneSourceSelection,
) -> np.ndarray:
    """Translate global embedded-catalog indices into one chromosome atom selector."""
    selected = np.asarray(resolution.catalog_indices, dtype=np.int64)
    rows = embedded_catalog.loc[
        embedded_catalog["gene_index"].isin(selected)
        & (embedded_catalog["chrom"].astype(str) == str(chrom))
        & embedded_catalog["included"].astype(bool),
        "chromosome_gene_row",
    ].astype(int)
    return assemble_atom_selector(atom_model.gene_to_atom, rows.to_numpy(dtype=np.int64))
