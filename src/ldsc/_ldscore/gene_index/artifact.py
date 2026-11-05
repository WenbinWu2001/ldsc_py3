"""The on-disk format of a gene LD-score index: writing and reading it.

Per-chromosome index payloads and their metadata, the index identity hash
(``index_id``) and the published-row hash, staged writes and their
finalization, and the readers: ``load_gene_ldscore_index`` verifies a
complete published index, and ``LoadedGeneLDScoreIndex.chromosome_loader``
reloads one chromosome for indexed scoring.
See ``docs/reference/gene-ldscore-index.md`` for the format contract.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import shutil
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import sparse

from ..._annotation.gene_lists import GeneCatalog, select_index_eligible_gene_indices
from ..._core.chromosomes import AUTOSOMES
from ..._core.snp_identity import identity_cleanup_keys, identity_mode_family, normalize_snp_identifier_mode
from ..._refpanel.funnel import FUNNEL_COLUMNS, concat_funnels, funnel_from_records, funnel_records
from ..._runtime.progress import advance, report_phase
from .atoms import AtomStatistics, ChromosomeAtomModel, validate_ldscore_operator
from ...config import GlobalConfig
from ...errors import LDSCInputError, LDSCInternalError
from ..._runtime.preflight import InputGate, check_readable, input_issue

LOGGER = logging.getLogger("LDSC.gene_ldscore_index")


# Version 6: baseline counts and overlaps accumulate in float64, as in direct
# ldscore, and root metadata stores the annotation_types classification of the
# baseline columns and the snp_universe_policy that indexed output records.
# Version 5: row A1 is the minor
# allele in the retained reference individuals, not the .bim A1; index_id
# hashes only the inputs, so the version marks the change. Version 4: each
# chromosome stores its reference-SNP funnel, which indexed scoring reports.
# Version 3 made PLINK LD windows exact; version 2 operators used legacy
# LDSC's chunk-rounded windows, which depended on --snp-batch-size.
_INDEX_FORMAT_VERSION = 6
# Chromosomes every published index covers. Callers read it through this
# module when they run, so small test indexes can narrow it.
INDEX_CHROMOSOMES = AUTOSOMES


def _semantic_sha256(payload: dict) -> str:
    """Hash one scientific identity with canonical JSON normalization."""
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def calculate_index_id(index_identity: dict) -> str:
    """Return the immutable semantic identity of one complete index."""
    return _semantic_sha256(index_identity)


def index_global_config(snp_identifier: str, genome_build: str) -> GlobalConfig:
    """Return the identity settings an index binds: its literal mode and build.

    The index always records its explicit hg19/hg38 build. rsID-family
    identity carries no coordinate build, so the build is dropped here instead
    of triggering GlobalConfig's warning about an ignored build.
    """
    if identity_mode_family(snp_identifier) == "rsid":
        return GlobalConfig(snp_identifier=snp_identifier)
    return GlobalConfig(snp_identifier=snp_identifier, genome_build=genome_build)


def _remove_gene_index_transaction(path: Path) -> bool:
    """Best-effort removal of builder-owned staging data after a build error.

    A cleanup error is logged and must not replace the primary build exception.
    """
    if not path.exists():
        return True
    try:
        shutil.rmtree(path)
    except Exception as exc:
        LOGGER.warning(
            "Gene LD-score index staging cleanup also failed after an earlier build or "
            f"publication error ({type(exc).__name__}: {exc}). Retained builder-owned path: "
            f"{path}. The earlier exception remains the build failure."
        )
        return False
    return True


def _directory_payload_bytes(path: Path) -> int:
    """Return the byte size of one scientific artifact component."""
    if not path.exists():
        return 0
    return sum(file.stat().st_size for file in path.rglob("*") if file.is_file())


def _scientific_payload_bytes(index_path: Path) -> int:
    """Return index bytes excluding human-readable and JSON diagnostics."""
    total = 0
    for file in index_path.rglob("*"):
        if not file.is_file() or "diagnostics" in file.parts:
            continue
        total += file.stat().st_size
    return total


_HASH_CHUNK_ROWS = 65536


def _canonical_frame_sha256(
    frame: pd.DataFrame,
    *,
    sort_by: Sequence[str] = (),
) -> str:
    """Hash a normalized table without path, compression, or text-layout effects.

    The digest is the SHA-256 of the table as tab-separated text, header
    included, with one canonical token per cell: ``.17g`` for floats,
    ``true``/``false`` for Booleans, ``str`` otherwise and ``<NA>`` for missing
    values. NumPy-typed columns are tokenized and hashed in bounded row
    chunks, so their complete text is never held at once. Extension-typed
    columns are tokenized whole, because pandas may convert their elements
    by the whole column (a nullable integer column with missing values
    formats as floats).
    """
    columns = [str(column) for column in frame.columns]
    order = None
    if sort_by:
        keys = frame.loc[:, list(sort_by)]
        keys.index = pd.RangeIndex(len(keys))
        order = keys.sort_values(list(sort_by), kind="mergesort").index.to_numpy()
    whole = {}
    for position in range(len(columns)):
        if not isinstance(frame.dtypes.iloc[position], np.dtype):
            values = frame.iloc[:, position]
            whole[position] = _canonical_tokens(values if order is None else values.take(order))
    digest = hashlib.sha256()
    digest.update(pd.DataFrame(columns=columns).to_csv(index=False, sep="\t", lineterminator="\n").encode("utf-8"))
    for start in range(0, len(frame), _HASH_CHUNK_ROWS):
        rows = slice(start, start + _HASH_CHUNK_ROWS)
        chunk = frame.iloc[rows] if order is None else frame.take(order[rows])
        tokens = pd.DataFrame({position: whole[position][rows] if position in whole
                               else _canonical_tokens(chunk.iloc[:, position])
                               for position in range(len(columns))}, index=pd.RangeIndex(len(chunk)))
        text = tokens.to_csv(index=False, header=False, sep="\t", lineterminator="\n")
        digest.update(text.encode("utf-8"))
    return digest.hexdigest()


def _canonical_tokens(values: pd.Series) -> np.ndarray:
    """Return the canonical hash token of each value of one column, as an object array."""
    if pd.api.types.is_float_dtype(values):
        array = values.to_numpy()
        if isinstance(values.dtype, np.dtype) and array.dtype.itemsize in (2, 4, 8):
            # Format each distinct bit pattern once: equal patterns are equal
            # tokens, while 0.0 and -0.0 keep their own.
            codes, patterns = pd.factorize(array.view(f"i{array.dtype.itemsize}"))
            formatted = np.array([*("<NA>" if np.isnan(value) else format(float(value), ".17g")
                                    for value in patterns.view(array.dtype))], dtype=object)
            return formatted[codes] if len(codes) else np.array([], dtype=object)
        tokens = values.map(lambda value: "<NA>" if pd.isna(value) else format(float(value), ".17g"))
    elif pd.api.types.is_bool_dtype(values):
        tokens = values.astype(bool).map({True: "true", False: "false"})
    else:
        tokens = values.map(lambda value: "<NA>" if pd.isna(value) else str(value))
    return tokens.to_numpy(dtype=object)


def _catalog_identity_records(frame: pd.DataFrame) -> list[dict[str, object]]:
    """Return JSON-stable normalized catalog rows for the index identity."""
    columns = [
        "gene_id",
        "gene_name",
        "chrom",
        "start",
        "end",
        "genome_build",
        "catalog_line",
    ]
    normalized = frame.loc[:, columns].copy()
    normalized["gene_name"] = normalized["gene_name"].fillna("").astype(str)
    for column in ("gene_id", "chrom", "genome_build"):
        normalized[column] = normalized[column].astype(str)
    for column in ("start", "end", "catalog_line"):
        normalized[column] = pd.to_numeric(normalized[column], errors="raise").astype(int)
    return normalized.sort_values(
        ["gene_id", "catalog_line"], kind="stable"
    ).to_dict(orient="records")


@dataclass(frozen=True)
class IndexChromosomeData:
    """Validated in-memory scientific payload for one index chromosome.

    ``reference_snp_funnel`` counts the SNPs each reference-preparation step
    removed; writing a payload requires it. ``annotation_types`` (each
    baseline column's ``binary`` or ``quantitative`` classification on this
    chromosome) and ``regression_selected_snp_count`` (regression SNPs before
    the region subtraction) are build results the coordinator combines into
    root metadata; the payload does not store them.
    """

    baseline_rows: pd.DataFrame
    baseline_count_all: np.ndarray
    baseline_count_common: np.ndarray
    baseline_overlap_all: np.ndarray
    baseline_overlap_common: np.ndarray
    total_reference_snps_all: int
    total_reference_snps_common: int
    atom_model: ChromosomeAtomModel
    operator: sparse.csr_matrix
    atom_statistics: AtomStatistics
    reference_rows_before_genotype_qc: int = 0
    genotype_qc_removed: int = 0
    maf_removed: int = 0
    cm_source: str = "bim_cm"
    reference_snp_funnel: pd.DataFrame | None = None
    annotation_types: dict[str, str] = field(default_factory=dict)
    regression_selected_snp_count: int = 0


@dataclass(frozen=True)
class StagedIndexChromosome:
    """Compact worker result returned after one chromosome shard is durable.

    Attributes
    ----------
    chromosome : str
        Canonical chromosome label for the staged shard.
    evidence : dict
        Scalar identity inputs and build diagnostics used by the coordinator.
        This mapping contains no chromosome payload tables or sparse matrices.
    component_metadata : dict
        Dimensions, ordered baseline columns, and sparse-format declarations
        that receive the complete ``index_id`` during coordinator finalization.
    """

    chromosome: str
    evidence: dict
    component_metadata: dict


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _read_json(path: Path, label: str) -> dict:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise LDSCInputError(
            f"Gene LD-score index {label} metadata is missing or unreadable: {path}"
            " Check the named metadata file and restore it from the same complete index build; regenerate the "
            "index if the original files are unavailable."
        ) from exc
    if not isinstance(payload, dict):
        raise LDSCInputError(
            f"Gene LD-score index {label} metadata must be a JSON object."
            " Check the named metadata file and restore it from the same complete index build; regenerate the "
            "index if the original files are unavailable."
        )
    return payload


def require_current_index_format(root: dict) -> None:
    """Reject root metadata of an index written in an earlier format, asking for a rebuild.

    Both the loader and the declaration preflight of indexed ``ldscore`` call
    this before any component check: components of an earlier index carry
    the same earlier version, which a component check would report as a
    mismatched restore instead.
    """
    if root.get("format_version") != _INDEX_FORMAT_VERSION:
        raise LDSCInputError(
            f"Unsupported gene LD-score index format version {root.get('format_version')!r} "
            f"(this package reads version {_INDEX_FORMAT_VERSION}). Rebuild the complete index with the "
            "current package; earlier indexes cannot be reused or migrated. Version 6 sums baseline "
            "annotation counts in float64, as direct ldscore does, and stores the annotation_types and "
            "snp_universe_policy that indexed ldscore output records; version 5 gives each SNP's "
            "minor allele as A1, as other PLINK-derived files do; version 4 stores each "
            "chromosome's reference-SNP funnel, which indexed scoring reports; version 3 made PLINK LD "
            "windows exact, and earlier indexes' LD scores differ slightly from current direct ldscore runs."
        )


def _validate_metadata_identity(metadata: dict, *, index_id: str) -> None:
    if (metadata.get("artifact_type") != "gene_ldscore_index" or metadata.get("index_id") != index_id
            or metadata.get("format_version") != _INDEX_FORMAT_VERSION):
        raise LDSCInputError(
            "Gene LD-score index component identity is invalid."
            " Restore components from the same complete index build as its root metadata."
        )


def _cast_index_baseline_rows(frame: pd.DataFrame) -> pd.DataFrame:
    """Match the LD-score directory writer's persisted floating dtypes."""
    float64_columns = [column for column in frame.columns if frame[column].dtype == np.float64]
    return frame.astype({column: np.float32 for column in float64_columns}) if float64_columns else frame


def _published_row_metadata_sha256(frame: pd.DataFrame) -> str:
    """Bind ordered CHR/SNP/POS/A1/A2 metadata for all four SNP identifier modes."""
    required = ["CHR", "SNP", "POS", "A1", "A2"]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise LDSCInputError(
            "New gene LD-score indexes require PLINK-authored A1/A2 and canonical row "
            f"metadata; missing columns: {missing}."
        )
    return _canonical_frame_sha256(frame.loc[:, required])


def _index_chromosome_metadata(
    chrom: str,
    record: IndexChromosomeData,
    *,
    snp_identifier: str,
    genome_build: str,
) -> dict:
    """Return compact component metadata without the not-yet-known index ID."""
    baseline_columns = [
        column
        for column in record.baseline_rows.columns
        if column not in {"CHR", "SNP", "POS", "A1", "A2", "regression_ld_scores"}
    ]
    published_digest = _published_row_metadata_sha256(
        _cast_index_baseline_rows(record.baseline_rows)
    )
    return {
        "artifact_type": "gene_ldscore_index",
        "format_version": _INDEX_FORMAT_VERSION,
        "chromosome": str(chrom),
        "snp_identifier": snp_identifier,
        "genome_build": genome_build,
        "published_row_metadata_sha256": published_digest,
        "n_rows": len(record.baseline_rows),
        "n_baseline": len(baseline_columns),
        "baseline_columns": baseline_columns,
        "n_genes": record.atom_model.gene_to_atom.shape[0],
        "n_atoms": record.atom_model.n_atoms,
        "gene_to_atom_format": "csr_bool_int32",
        "ldscore_operator_format": "csr_float64_int32",
    }


def _write_index_chromosome_payload(
    chrom_path: Path,
    *,
    chrom: str,
    record: IndexChromosomeData,
) -> None:
    """Write one chromosome's large scientific payloads without shared metadata."""
    if record.reference_snp_funnel is None:
        raise LDSCInternalError(
            f"Gene LD-score index chromosome {chrom} has no reference-SNP funnel to store. "
            "Report the command, package version, and traceback."
        )
    _cast_index_baseline_rows(record.baseline_rows).to_parquet(
        chrom_path / "baseline_rows.parquet", index=False
    )
    record.reference_snp_funnel.to_parquet(chrom_path / "reference_snp_funnel.parquet", index=False)
    np.savez(
        chrom_path / "baseline_statistics.npz",
        baseline_count_all=np.asarray(record.baseline_count_all, dtype=np.float64),
        baseline_count_common=np.asarray(record.baseline_count_common, dtype=np.float64),
        baseline_overlap_all=np.asarray(record.baseline_overlap_all, dtype=np.float64),
        baseline_overlap_common=np.asarray(record.baseline_overlap_common, dtype=np.float64),
        total_reference_snps_all=np.asarray(record.total_reference_snps_all, dtype=np.int64),
        total_reference_snps_common=np.asarray(record.total_reference_snps_common, dtype=np.int64),
    )
    atoms = pd.DataFrame(
        {
            "atom_id": np.arange(record.atom_model.n_atoms, dtype=np.int64),
            "CHR": str(chrom),
            "start0": record.atom_model.starts,
            "end": record.atom_model.ends,
        }
    )
    atoms.to_parquet(chrom_path / "atoms.parquet", index=False)
    gene_to_atom = record.atom_model.gene_to_atom.astype(bool).tocsr()
    gene_to_atom.sort_indices()
    operator = record.operator.astype(np.float64, copy=False).tocsr()
    operator.sort_indices()
    sparse.save_npz(chrom_path / "gene_to_atom.npz", gene_to_atom, compressed=True)
    sparse.save_npz(chrom_path / "ldscore_operator.npz", operator, compressed=True)
    np.savez(
        chrom_path / "atom_statistics.npz",
        atom_count_all=np.asarray(record.atom_statistics.atom_count_all, dtype=np.int64),
        atom_count_common=np.asarray(record.atom_statistics.atom_count_common, dtype=np.int64),
        baseline_atom_overlap_all=np.asarray(record.atom_statistics.baseline_atom_overlap_all, dtype=np.float64),
        baseline_atom_overlap_common=np.asarray(record.atom_statistics.baseline_atom_overlap_common, dtype=np.float64),
    )


def _stage_index_chromosome(
    staged_index: Path,
    chrom: str,
    record: IndexChromosomeData,
) -> Path:
    """Atomically install one payload-only shard in the private run stage.

    All payload files are first closed in a chromosome-private temporary
    directory. Renaming that directory to ``chrN`` is the durability boundary
    observed by the coordinator's subsequent ``Finished chromosome N`` log record.
    Shared and component identity metadata are intentionally deferred to the
    coordinator.
    """
    chromosomes_path = staged_index / "chromosomes"
    chromosomes_path.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".chr{chrom}.tmp-", dir=chromosomes_path)
    )
    destination = chromosomes_path / f"chr{chrom}"
    try:
        _write_index_chromosome_payload(temporary, chrom=chrom, record=record)
        os.replace(temporary, destination)
    except Exception:
        _remove_gene_index_transaction(temporary)
        raise
    return destination


def _finalize_staged_index_artifact(
    index_path: Path,
    *,
    index_id: str,
    index_identity: dict,
    gene_catalog: pd.DataFrame,
    chromosomes: dict[str, StagedIndexChromosome],
    annotation_types: dict[str, str],
    snp_universe_policy: dict,
    diagnostic_payload: dict | None = None,
) -> None:
    """Finalize coordinator-owned metadata for a complete staged index.

    ``chromosomes`` is already in canonical requested order. This function
    writes root identity, component identity, the gene coordinate catalog, and successful
    diagnostics only after every payload-only shard is durable.
    ``annotation_types`` (one classification per baseline column) and
    ``snp_universe_policy`` are root metadata outside ``index_identity``;
    indexed ``ldscore`` output records them.
    """
    _write_json(
        index_path / "metadata.json",
        {
            "artifact_type": "gene_ldscore_index",
            "format_version": _INDEX_FORMAT_VERSION,
            "index_id": index_id,
            "snp_identifier": index_identity["snp_identifier"],
            "genome_build": index_identity["genome_build"],
            "index_identity": index_identity,
            "chromosomes": list(chromosomes),
            "annotation_types": dict(annotation_types),
            "snp_universe_policy": dict(snp_universe_policy),
        },
    )
    gene_catalog.to_parquet(index_path / "gene_catalog.parquet", index=False)
    for chrom, result in chromosomes.items():
        _write_json(
            index_path / "chromosomes" / f"chr{chrom}" / "metadata.json",
            {**result.component_metadata, "index_id": index_id},
        )
    if diagnostic_payload is not None:
        diagnostic_payload["payload_bytes"] = _scientific_payload_bytes(index_path)
        chromosome_diagnostics = diagnostic_payload.get("chromosomes", {})
        for chrom in chromosomes:
            values = chromosome_diagnostics.get(chrom)
            if isinstance(values, dict):
                values["payload_bytes"] = _directory_payload_bytes(
                    index_path / "chromosomes" / f"chr{chrom}"
                )
        _write_json(
            index_path / "diagnostics" / "build-gene-ldscore-index.json",
            diagnostic_payload,
        )


def inspect_index_paths(directory, *, metadata=None):
    """Inspect root declarations and all component paths without loading arrays."""
    root = Path(directory)
    if metadata is None:
        metadata = _read_json(root/'metadata.json', 'root')
    require_current_index_format(metadata)
    chromosomes = metadata.get('chromosomes', [])
    if not isinstance(chromosomes, (list, tuple)) or not chromosomes or any(str(chrom) not in AUTOSOMES for chrom in chromosomes):
        raise LDSCInputError(
            'Gene LD-score index chromosome coverage is invalid.'
            " Restore a complete gene LD-score index covering chromosomes 1 through 22, or rebuild it with "
            "ldsc build-gene-ldscore-index. Do not combine components from different builds."
        )
    gate = InputGate('Index component declarations')
    paths = [('gene_catalog.parquet', '')]
    files = ('metadata.json', 'baseline_rows.parquet', 'baseline_statistics.npz', 'atoms.parquet',
             'gene_to_atom.npz', 'ldscore_operator.npz', 'atom_statistics.npz', 'reference_snp_funnel.parquet')
    paths.extend((f'chromosomes/chr{chrom}/{name}', str(chrom)) for chrom in chromosomes for name in files)
    for name, chrom in paths:
        path = root/name
        if not path.is_file():
            gate.issues.append(input_issue('index', path, chrom, 'missing_required_input', 'Required index component is missing.'))
        elif name.endswith('.parquet'):
            gate.check('index', path, lambda path=path: check_readable(path), chrom=chrom)
        elif chrom and name.endswith('metadata.json'):
            def component_header():
                component = _read_json(path, f'chromosome {chrom}')
                _validate_chromosome_metadata(component, chrom=chrom, index_id=metadata.get('index_id'),
                    snp_identifier=metadata.get('snp_identifier'), genome_build=metadata.get('genome_build'))
            gate.check('index', path, component_header, chrom=chrom)
    gate.finish()
    return metadata


def _load_npz_members(path: Path, required: set[str]) -> dict[str, np.ndarray]:
    try:
        with np.load(path, allow_pickle=False) as archive:
            if set(archive.files) != required:
                raise LDSCInputError(
                    f"Gene LD-score index payload {path.name} has required members {sorted(required)}, "
                    f"found {sorted(archive.files)}."
                    " Restore this saved component from the same index build, or rebuild the complete index; "
                    "do not replace array members independently."
                )
            return {name: archive[name] for name in archive.files}
    except (OSError, ValueError) as exc:
        raise LDSCInputError(
            f"Gene LD-score index payload is unreadable: {path}"
            " Restore this saved component from the same index build, or rebuild the complete index; do not "
            "replace array members independently."
        ) from exc


def _validate_chromosome_metadata(component_meta, *, chrom, index_id, snp_identifier, genome_build):
    """Check small component metadata before opening any chromosome payload."""
    _validate_metadata_identity(component_meta, index_id=index_id)
    if (
        component_meta.get("snp_identifier") != snp_identifier
        or component_meta.get("genome_build") != genome_build
    ):
        raise LDSCInputError(
            "Gene LD-score index component identity metadata disagrees with the root."
            " Restore chromosome metadata from the same index build as the root metadata, or rebuild the "
            "complete index."
        )
    if component_meta.get("chromosome") != str(chrom):
        raise LDSCInputError(
            "Gene LD-score index component chromosome identity is invalid."
            " Restore chromosome metadata from the same index build as the root metadata, or rebuild the "
            "complete index."
        )


def _in_position_snp_order(frame: pd.DataFrame) -> bool:
    """Return whether rows are sorted by POS, then SNP, so a stable sort would not move them.

    Missing SNP labels sort last among equal positions, as in ``sort_values``.
    """
    step = np.diff(frame["POS"].to_numpy())
    if (step < 0).any():
        return False
    left = np.flatnonzero(step == 0)
    if not len(left):
        return True
    right = left + 1
    labels = frame["SNP"]
    missing = labels.isna().to_numpy()
    values = labels.to_numpy()
    ordered = missing[right].copy()
    both = ~missing[left] & ~missing[right]
    ordered[both] = values[left[both]] <= values[right[both]]
    return bool(ordered.all())


def _load_index_chromosome(
    index_path: Path,
    chrom: str,
    *,
    index_id: str,
    snp_identifier: str,
    genome_build: str,
    expected_gene_rows: int,
) -> IndexChromosomeData:
    component_path = index_path / "chromosomes" / f"chr{chrom}"
    component_meta = _read_json(component_path / "metadata.json", f"chromosome {chrom}")
    _validate_chromosome_metadata(component_meta, chrom=chrom, index_id=index_id,
                                 snp_identifier=snp_identifier, genome_build=genome_build)
    baseline_rows = pd.read_parquet(component_path / "baseline_rows.parquet")
    if len(baseline_rows) != int(component_meta.get("n_rows", -1)):
        raise LDSCInputError(
            "Gene LD-score index baseline row count disagrees with component metadata."
            " Restore this chromosome component together with its matching root index; do not mix or hand-edit "
            "index components. Rebuild the complete index if the original valid files are unavailable."
        )
    baseline_columns = list(component_meta.get("baseline_columns", ()))
    allele_columns = [column for column in ("A1", "A2") if column in baseline_rows.columns]
    if allele_columns != ["A1", "A2"]:
        raise LDSCInputError(
            "Gene LD-score index baseline rows require both PLINK allele columns A1 and A2."
            " Restore this chromosome component together with its matching root index; do not mix or hand-edit "
            "index components. Rebuild the complete index if the original valid files are unavailable."
        )
    required_baseline = [
        "CHR", "SNP", "POS", *allele_columns, "regression_ld_scores", *baseline_columns
    ]
    if list(baseline_rows.columns) != required_baseline:
        raise LDSCInputError(
            "Gene LD-score index baseline rows have an invalid schema."
            " Restore this chromosome component together with its matching root index; do not mix or hand-edit "
            "index components. Rebuild the complete index if the original valid files are unavailable."
        )
    if any(baseline_rows[column].dtype != np.float32 for column in ["regression_ld_scores", *baseline_columns]):
        raise LDSCInputError(
            "Gene LD-score index baseline LD-score columns must use float32 dtype."
            " Restore this chromosome component together with its matching root index; do not mix or hand-edit "
            "index components. Rebuild the complete index if the original valid files are unavailable."
        )
    if baseline_rows["CHR"].astype(str).ne(str(chrom)).any():
        raise LDSCInputError(
            "Gene LD-score index baseline rows contain the wrong chromosome."
            " Restore this chromosome component together with its matching root index; do not mix or hand-edit "
            "index components. Rebuild the complete index if the original valid files are unavailable."
        )
    if baseline_rows["POS"].dtype != np.int64:
        raise LDSCInputError(
            "Gene LD-score index baseline POS must use int64 dtype."
            " Restore this chromosome component together with its matching root index; do not mix or hand-edit "
            "index components. Rebuild the complete index if the original valid files are unavailable."
        )
    _, identity_drops = identity_cleanup_keys(baseline_rows, snp_identifier, context=f"index chromosome {chrom}",
                                              stage="index_validation")
    if not identity_drops.empty:
        raise LDSCInputError(
            "Gene LD-score index baseline rows contain duplicate effective identities or invalid/ambiguous allele groups."
            " Restore this chromosome component together with its matching root index; do not mix or hand-edit "
            "index components. Rebuild the complete index if the original valid files are unavailable."
        )
    if not _in_position_snp_order(baseline_rows):
        raise LDSCInputError(
            "Gene LD-score index baseline rows are not in canonical genomic order."
            " Restore this chromosome component together with its matching root index; do not mix or hand-edit "
            "index components. Rebuild the complete index if the original valid files are unavailable."
        )
    published_digest = _published_row_metadata_sha256(baseline_rows)
    if component_meta.get("published_row_metadata_sha256") != published_digest:
        raise LDSCInputError(
            "Gene LD-score index published row metadata digest is invalid."
            " Restore this chromosome component together with its matching root index; do not mix or hand-edit "
            "index components. Rebuild the complete index if the original valid files are unavailable."
        )
    stats = _load_npz_members(
        component_path / "baseline_statistics.npz",
        {
            "baseline_count_all", "baseline_count_common", "baseline_overlap_all",
            "baseline_overlap_common", "total_reference_snps_all", "total_reference_snps_common",
        },
    )
    b = int(component_meta.get("n_baseline", -1))
    for name in ("baseline_count_all", "baseline_count_common"):
        if stats[name].shape != (b,) or stats[name].dtype != np.float64 or not np.isfinite(stats[name]).all():
            raise LDSCInputError(
                f"Gene LD-score index {name} has invalid shape, dtype, or values."
                " Restore this chromosome component together with its matching root index; do not mix or "
                "hand-edit index components. Rebuild the complete index if the original valid files are "
                "unavailable."
            )
    for name in ("baseline_overlap_all", "baseline_overlap_common"):
        if stats[name].shape != (b, b) or stats[name].dtype != np.float64 or not np.isfinite(stats[name]).all():
            raise LDSCInputError(
                f"Gene LD-score index {name} has invalid shape, dtype, or values."
                " Restore this chromosome component together with its matching root index; do not mix or "
                "hand-edit index components. Rebuild the complete index if the original valid files are "
                "unavailable."
            )
    atoms = pd.read_parquet(component_path / "atoms.parquet")
    if list(atoms.columns) != ["atom_id", "CHR", "start0", "end"]:
        raise LDSCInputError(
            "Gene LD-score index atoms.parquet schema is invalid."
            " Restore this chromosome component together with its matching root index; do not mix or hand-edit "
            "index components. Rebuild the complete index if the original valid files are unavailable."
        )
    if any(atoms[column].dtype != np.int64 for column in ("atom_id", "start0", "end")):
        raise LDSCInputError(
            "Gene LD-score index atom identifiers and coordinates must use int64 dtype."
            " Restore this chromosome component together with its matching root index; do not mix or hand-edit "
            "index components. Rebuild the complete index if the original valid files are unavailable."
        )
    if atoms["CHR"].astype(str).ne(str(chrom)).any():
        raise LDSCInputError(
            "Gene LD-score index atoms contain the wrong chromosome."
            " Restore this chromosome component together with its matching root index; do not mix or hand-edit "
            "index components. Rebuild the complete index if the original valid files are unavailable."
        )
    if not np.array_equal(atoms["atom_id"].to_numpy(), np.arange(len(atoms), dtype=np.int64)):
        raise LDSCInputError(
            "Gene LD-score index atom_id ordering is invalid."
            " Restore this chromosome component together with its matching root index; do not mix or hand-edit "
            "index components. Rebuild the complete index if the original valid files are unavailable."
        )
    try:
        gene_to_atom = sparse.load_npz(component_path / "gene_to_atom.npz").tocsr()
        operator = sparse.load_npz(component_path / "ldscore_operator.npz").tocsr()
    except (OSError, ValueError) as exc:
        raise LDSCInputError(
            "Gene LD-score index sparse payload is unreadable."
            " Restore this chromosome component together with its matching root index; do not mix or hand-edit "
            "index components. Rebuild the complete index if the original valid files are unavailable."
        ) from exc
    model = ChromosomeAtomModel(
        chromosome=str(chrom),
        starts=atoms["start0"].to_numpy(dtype=np.int64),
        ends=atoms["end"].to_numpy(dtype=np.int64),
        gene_to_atom=gene_to_atom,
    )
    n_atoms = int(component_meta.get("n_atoms", -1))
    n_rows = int(component_meta.get("n_rows", -1))
    if (
        model.n_atoms != n_atoms
        or model.gene_to_atom.shape[0] != int(component_meta.get("n_genes", -1))
        or model.gene_to_atom.shape[0] != expected_gene_rows
    ):
        raise LDSCInputError(
            "Gene LD-score index atom dimensions disagree with metadata."
            " Restore this chromosome component together with its matching root index; do not mix or hand-edit "
            "index components. Rebuild the complete index if the original valid files are unavailable."
        )
    if n_rows != len(baseline_rows):
        raise LDSCInputError(
            "Gene LD-score index operator row count disagrees with baseline rows."
            " Restore this chromosome component together with its matching root index; do not mix or hand-edit "
            "index components. Rebuild the complete index if the original valid files are unavailable."
        )
    try:
        model.validate()
        validate_ldscore_operator(operator, n_rows=n_rows, n_atoms=n_atoms)
    except LDSCInternalError as exc:
        raise LDSCInputError(
            f"Gene LD-score index sparse structure is invalid: {exc}"
            " Restore this chromosome component together with its matching root index; do not mix or hand-edit "
            "index components. Rebuild the complete index if the original valid files are unavailable."
        ) from exc
    atom_stats = _load_npz_members(
        component_path / "atom_statistics.npz",
        {"atom_count_all", "atom_count_common", "baseline_atom_overlap_all", "baseline_atom_overlap_common"},
    )
    for name in ("atom_count_all", "atom_count_common"):
        if atom_stats[name].shape != (n_atoms,) or atom_stats[name].dtype != np.int64:
            raise LDSCInputError(
                f"Gene LD-score index {name} has invalid shape or dtype."
                " Restore this chromosome component together with its matching root index; do not mix or "
                "hand-edit index components. Rebuild the complete index if the original valid files are "
                "unavailable."
            )
    for name in ("baseline_atom_overlap_all", "baseline_atom_overlap_common"):
        if atom_stats[name].shape != (b, n_atoms) or atom_stats[name].dtype != np.float64 or not np.isfinite(atom_stats[name]).all():
            raise LDSCInputError(
                f"Gene LD-score index {name} has invalid shape, dtype, or values."
                " Restore this chromosome component together with its matching root index; do not mix or "
                "hand-edit index components. Rebuild the complete index if the original valid files are "
                "unavailable."
            )
    for name in ("total_reference_snps_all", "total_reference_snps_common"):
        if stats[name].shape != () or stats[name].dtype != np.int64:
            raise LDSCInputError(
                f"Gene LD-score index {name} must be a scalar int64 value."
                " Restore this chromosome component together with its matching root index; do not mix or "
                "hand-edit index components. Rebuild the complete index if the original valid files are "
                "unavailable."
            )
    funnel = pd.read_parquet(component_path / "reference_snp_funnel.parquet")
    if list(funnel.columns) != FUNNEL_COLUMNS or set(funnel["chrom"].astype(str)) != {str(chrom)}:
        raise LDSCInputError(
            "Gene LD-score index reference_snp_funnel.parquet schema or chromosome is invalid."
            " Restore this chromosome component together with its matching root index; do not mix or "
            "hand-edit index components. Rebuild the complete index if the original valid files are unavailable."
        )
    return IndexChromosomeData(
        baseline_rows=baseline_rows,
        baseline_count_all=stats["baseline_count_all"],
        baseline_count_common=stats["baseline_count_common"],
        baseline_overlap_all=stats["baseline_overlap_all"],
        baseline_overlap_common=stats["baseline_overlap_common"],
        total_reference_snps_all=int(stats["total_reference_snps_all"]),
        total_reference_snps_common=int(stats["total_reference_snps_common"]),
        atom_model=model,
        operator=operator,
        atom_statistics=AtomStatistics(
            atom_stats["atom_count_all"],
            atom_stats["atom_count_common"],
            atom_stats["baseline_atom_overlap_all"],
            atom_stats["baseline_atom_overlap_common"],
        ),
        reference_snp_funnel=funnel_from_records(funnel_records(funnel)),
    )


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
    annotation_types : dict of str to {"binary", "quantitative"}
        Classification of each baseline column, made at build time as direct
        ``ldscore`` makes it, in baseline-column order.
    snp_universe_policy : dict
        Reference-SNP and regression-SNP record of the build, in the form
        direct ``ldscore`` writes to ``metadata.json``.
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
    annotation_types: dict[str, str]
    snp_universe_policy: dict

    @property
    def global_config(self) -> GlobalConfig:
        """Identity settings the index binds; indexed output records them."""
        return index_global_config(self.snp_identifier, self.genome_build)

    def chromosome_loader(self) -> "_ChromosomeLoader":
        """Return the small descriptor that chromosome tasks reload components from."""
        labels = self.gene_catalog.chrom.astype(str)
        return _ChromosomeLoader(
            self.index_id, self.snp_identifier, self.genome_build, len(self.baseline_columns),
            {str(chrom): path.parent.parent for chrom, path in self.chromosome_paths.items()},
            {str(chrom): int(labels.eq(str(chrom)).sum()) for chrom in self.chromosome_paths},
        )


@dataclass(frozen=True)
class _ChromosomeLoader:
    """Reloads chromosome components without the index catalog, so tasks stay small."""

    index_id: str
    snp_identifier: str
    genome_build: str
    n_baseline_columns: int
    index_paths: dict[str, Path]
    expected_gene_rows: dict[str, int]

    def load(self, chrom: str) -> IndexChromosomeData:
        """Load one owned chromosome payload; callers release it after use.

        The immutable component is checked again when opened. Indexed scoring
        reloads it for each active query batch and releases it before the next
        chromosome task; operators are not cached between tasks.
        """
        return _load_index_chromosome(
            self.index_paths[str(chrom)], str(chrom), index_id=self.index_id,
            snp_identifier=self.snp_identifier, genome_build=self.genome_build,
            expected_gene_rows=self.expected_gene_rows[str(chrom)],
        )


@report_phase(LOGGER, 'validation', 'index content and compatibility')
def load_gene_ldscore_index(
    index_dir: str | Path,
) -> LoadedGeneLDScoreIndex:
    """Load and fully validate one explicit gene LD-score index directory.

    Parameters
    ----------
    index_dir : path-like
        Directory containing one complete format-version-6 index covering
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
    require_current_index_format(root)
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
    if chromosomes != INDEX_CHROMOSOMES:
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
    annotation_types = root.get("annotation_types")
    snp_universe_policy = root.get("snp_universe_policy")
    if (not isinstance(annotation_types, dict) or set(annotation_types) != set(baseline_columns)
            or not set(annotation_types.values()) <= {"binary", "quantitative"}
            or not isinstance(snp_universe_policy, dict)):
        raise LDSCInputError(
            "Gene LD-score index root annotation_types or snp_universe_policy is missing or invalid."
            " Restore the complete saved index from one compatible build, or regenerate it with ldsc "
            "build-gene-ldscore-index. See "
            "docs/user-guide/troubleshooting.md#ldscore-an-explicit-gene-ld-score-index-is-missing-corrupt-or-incompatible."
        )
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
        annotation_types={name: annotation_types[name] for name in baseline_columns},
        snp_universe_policy=snp_universe_policy,
    )


def _baseline_columns_from_rows(frame: pd.DataFrame) -> list[str]:
    """Return ordered supplied baseline score columns from one common shard."""
    reserved = {"CHR", "SNP", "POS", "A1", "A2", "regression_ld_scores"}
    return [column for column in frame.columns if column not in reserved]
