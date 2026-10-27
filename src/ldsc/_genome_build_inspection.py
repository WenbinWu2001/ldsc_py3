"""Read-only artifact inspection for ``ldsc infer-genome-build``.

Reports coordinate provenance independently of the SNP matching mode. This
helper does not resolve workflow configuration or verify declared coordinates.
"""

from dataclasses import dataclass, field
import gzip
import json
from pathlib import Path

import pandas as pd
from pyarrow import ArrowNotImplementedError
import pyarrow.parquet as pq

from ._coordinates import normalize_chr_pos_frame
from ._genome_build_policy import _SAMPLE_ROWS, _table_frames, _text_options
from ._result_files import declared_result_file
from .column_inference import (
    RESTRICTION_CHRPOS_SPEC_MAP, RESTRICTION_HG19_POS_SPEC, RESTRICTION_HG38_POS_SPEC,
    normalize_genome_build, resolve_optional_column,
)
from .errors import LDSCInputError
from .genome_build_inference import infer_chr_pos_build
from .path_resolution import normalize_path_token


_RSID_NOTICE = "This artifact’s workflow uses rsIDs for SNP matching."
_NO_EVIDENCE = "no build metadata or usable SNP coordinates"
_TROUBLESHOOTING = "docs/user-guide/troubleshooting.md#infer-genome-build-input-cannot-be-inspected"
_JSON_TYPES = {
    "annotation_projection", "ldscore", "ref_panel", "gene_ldscore_index",
    "h2_result", "partitioned_h2_result", "partitioned_h2_query_result", "rg_result",
    "rg_pair_result", "quantile_h2_result", "query_r2_result", "plot_result", "h2_scale_conversion_result",
}
_BUILD_KEYS = ("genome_build", "projection_genome_build", "index_genome_build")


@dataclass(frozen=True)
class _Evidence:
    source: str
    build: str | None
    detail: str


@dataclass
class _Inspection:
    # Each group describes one coordinate system. Named groups are alternatives.
    groups: dict[str, list[_Evidence]] = field(default_factory=dict)
    modes: set[str] = field(default_factory=set)

    def add(self, source, build=None, detail=_NO_EVIDENCE, *, group=""):
        self.groups.setdefault(group, []).append(_Evidence(str(source), build, detail))

    def metadata(self, metadata, source, *, keys=("genome_build", "sorted_by_build")):
        """Record declarations without scoring or assigning coordinates."""
        mode = metadata.get("snp_identifier")
        if mode:
            self.modes.add(str(mode))
        evidence = []
        for key in keys:
            value = metadata.get(key)
            build = normalize_genome_build(value) if value else None
            if build in {"hg19", "hg38"}:
                evidence.append(_Evidence(f"{source} [{key}]", build, "declared in artifact metadata"))
        return evidence

    def merge(self, other, *, group=None):
        for label, entries in other.groups.items():
            self.groups.setdefault(label if group is None else group, []).extend(entries)
        self.modes.update(other.modes)

    def render(self):
        """Summarize complete evidence, distinguishing uncertainty from conflicts."""
        lines, conflicts = [], []
        for group, evidence in sorted(self.groups.items()):
            builds = {entry.build for entry in evidence if entry.build}
            if len(builds) > 1:
                conflicts.extend(f"{entry.source}: {entry.build}" for entry in evidence if entry.build)
                continue
            unresolved = [entry for entry in evidence if entry.build is None]
            if unresolved or not builds:
                build = "unknown"
                detail = unresolved[0].detail if unresolved else _NO_EVIDENCE
                if builds:
                    detail = f"some chromosomes or artifact members remain unresolved: {unresolved[0].source}"
            else:
                build = next(iter(builds))
                detail = "; ".join(dict.fromkeys(entry.detail for entry in evidence))
            label = f" [{group}]" if group else ""
            lines.append(f"Genome build{label}: {build} ({detail}).")
        if conflicts:
            raise LDSCInputError(
                "Conflicting genome-build evidence for the same coordinate set:\n"
                + "\n".join(dict.fromkeys(conflicts))
                + "\nCheck the listed declarations and coordinates; supply a consistent artifact."
            )
        if self.modes and self.modes <= {"rsid", "rsid_allele_aware"}:
            lines.append(_RSID_NOTICE)
        return "\n".join(lines)


def inspect_genome_build(path: str) -> str:
    """Return a terminal report for one exact artifact path, without writes.

    Parameters
    ----------
    path : str
        Exact supported file or recognized artifact directory. Patterns and
        chromosome suites are not expanded.

    Returns
    -------
    str
        Declared or inferred hg19/hg38, labeled coordinate alternatives, or
        ``unknown`` for insufficient evidence, followed by an rsID notice only
        when matching mode is recorded. No coordinate conversion is performed.

    Raises
    ------
    LDSCUserError
        An input is missing, unreadable, malformed, unsupported, or contributes
        contradictory build evidence. Workflow configuration is never changed.
    """
    path = Path(normalize_path_token(path))
    try:
        if path.is_dir():
            return _inspect_directory(path.resolve()).render()
        if not path.is_file():
            raise LDSCInputError(f"Input is not a readable artifact: {path}. Supply one exact file or artifact directory.")
        if path.suffix.lower() == ".json":
            return _inspect_metadata(_read_metadata(path), path).render()
        return _inspect_table(path).render()
    except (OSError, ValueError, EOFError, ArrowNotImplementedError) as exc:
        detail = str(exc).rstrip(".")
        raise LDSCInputError(
            f"Could not inspect {path}: {detail}. Check the file format and readability. "
            f"See {_TROUBLESHOOTING}."
        ) from exc


def _inspect_table(path: Path, *, inherited=None) -> _Inspection:
    result = _Inspection()
    suffix = path.name.removesuffix(".gz").rsplit(".", 1)[-1].lower()
    parquet = suffix == "parquet"
    if parquet:
        schema = pq.read_schema(path)
        metadata = {key.decode().removeprefix("ldsc:"): value.decode()
                    for key, value in (schema.metadata or {}).items() if key.startswith(b"ldsc:")}
        columns = schema.names
    else:
        opener = gzip.open if path.name.endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8") as stream:
            probe = stream.read(4096)
        if any(ord(char) < 32 and char not in "\t\r\n" for char in probe):
            raise LDSCInputError(f"Unsupported binary input: {path}. Supply a BIM, text table, Parquet, or LDSC metadata file.")
        options, metadata = _text_options(path)
        if suffix in {"bim", "bed"}:
            columns = []
        else:
            columns = list(pd.read_csv(path, nrows=0, **options).columns)
    declared = result.metadata(metadata, path)
    owner = inherited if inherited is not None else _owning_metadata(path)
    if owner:
        declared.extend(result.metadata(owner, f"{path} owning metadata", keys=_BUILD_KEYS))
    sidecar = path.with_name(path.name.removesuffix("_r2.parquet") + "_meta.tsv.gz")
    if path.name.endswith("_r2.parquet") and sidecar.is_file():
        sidecar = _member(path.resolve().parent, sidecar.name)
        declared.extend(result.metadata(_text_options(sidecar)[1], sidecar))
        if not declared:
            result.merge(_inspect_table(sidecar, inherited=owner))
            return result
    if declared:
        result.groups[""] = declared
        return result
    if suffix == "bed":
        result.add(path, detail="BED intervals do not establish a genome build without metadata")
        return result
    if suffix == "bim":
        with pd.read_csv(path, sep=r"\s+", header=None, chunksize=50000) as reader:
            def frames():
                for frame in reader:
                    if frame.shape[1] != 6:
                        raise LDSCInputError(f"Malformed BIM {path}: expected six columns.")
                    yield frame.loc[:, [0, 3]].rename(columns={0: "CHR", 3: "POS"})
            _infer_frames(result, path, frames())
        return result
    build_column = next((col for col in columns if str(col).lower() == "genome_build"), None)
    if build_column:
        for frame in _table_frames(path, [build_column], parquet):
            for value in frame[build_column].dropna().unique():
                result.groups.setdefault("", []).extend(result.metadata({"genome_build": value}, path))
        if result.groups.get(""):
            return result
    alternatives = [(build, resolve_optional_column(columns, spec, context=str(path))) for build, spec in (
        ("hg19", RESTRICTION_HG19_POS_SPEC), ("hg38", RESTRICTION_HG38_POS_SPEC),
    )]
    alternatives = [(build, column) for build, column in alternatives if column is not None]
    if alternatives:
        for build, column in alternatives:
            result.add(path, build, f"declared by coordinate column {column}",
                       group=build if len(alternatives) > 1 else "")
        return result
    chr_col, pos_col = [resolve_optional_column(columns, RESTRICTION_CHRPOS_SPEC_MAP[key], context=str(path))
                        for key in ("CHR", "POS")]
    if chr_col is None or pos_col is None:
        result.add(path)
    else:
        frames = (frame.rename(columns={chr_col: "CHR", pos_col: "POS"})
                  for frame in _table_frames(path, [chr_col, pos_col], parquet))
        _infer_frames(result, path, frames)
    return result


def _read_metadata(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if (not isinstance(payload, dict) or not isinstance(payload.get("artifact_type"), str)
            or payload["artifact_type"] not in _JSON_TYPES):
        raise LDSCInputError(f"Unrecognized LDSC metadata: {path}. Supply metadata from a supported LDSC artifact.")
    if "files" in payload and not isinstance(payload["files"], dict):
        raise LDSCInputError(f"Malformed files manifest in {path}: expected an object.")
    if payload["artifact_type"] == "gene_ldscore_index" and not isinstance(payload.get("chromosomes", []), list):
        raise LDSCInputError(f"Malformed chromosomes in {path}: expected a list.")
    return payload


def _inspect_metadata(payload: dict, source: Path) -> _Inspection:
    result = _Inspection()
    if payload["artifact_type"] == "ref_panel":
        result.metadata(payload, source, keys=())
        builds = payload.get("emitted_genome_builds", [])
        if not isinstance(builds, list):
            raise LDSCInputError(f"Malformed emitted_genome_builds in {source}: expected a list.")
        for value in builds:
            build = normalize_genome_build(value)
            if build in {"hg19", "hg38"}:
                result.add(source, build, "declared in artifact metadata", group=build)
    else:
        declared = result.metadata(payload, source, keys=_BUILD_KEYS)
        if payload["artifact_type"] == "gene_ldscore_index" and "index_identity" in payload:
            identity = payload["index_identity"]
            if not isinstance(identity, dict):
                raise LDSCInputError(f"Malformed index_identity in {source}: expected an object.")
            declared.extend(result.metadata(identity, f"{source} index_identity"))
        if declared:
            result.groups[""] = declared
    if not result.groups:
        result.add(source)
    return result


def _file_tokens(payload: dict):
    """Yield only coordinate-bearing members of established artifact manifests."""
    kind = payload["artifact_type"]
    for key, values in payload.get("files", {}).items():
        selected = (
            (kind == "annotation_projection" and key in {"baseline_annotations", "query_annotations"})
            or (kind == "ldscore" and (key in {"baseline", "query"} or key.startswith("query_batch")))
            or (kind == "ref_panel" and key in {"r2_hg19", "r2_hg38", "meta_hg19", "meta_hg38"})
        )
        if selected:
            for token in values if isinstance(values, list) else [values]:
                yield key, token


def _member(root: Path, token) -> Path:
    return declared_result_file(
        root.resolve(), {"files": {"member": token}}, "member", context="infer-genome-build",
    )


def _owning_metadata(path: Path) -> dict | None:
    """Use neighboring metadata only when its manifest owns this exact file."""
    root = path.resolve().parent
    for candidate in (root / "metadata.json", root / "diagnostics" / "metadata.json"):
        if not candidate.is_file():
            continue
        try:
            payload = _read_metadata(_member(root, str(candidate.relative_to(root))))
        except (OSError, ValueError, LDSCInputError):
            continue
        tokens = [token for values in payload.get("files", {}).values()
                  for token in (values if isinstance(values, list) else [values])]
        if payload["artifact_type"] == "gene_ldscore_index":
            tokens.append("gene_catalog.parquet")
        manifest_root = root
        nested_folder = {"partitioned_h2_query_result": "query_annotations", "rg_pair_result": "pairs"}.get(payload["artifact_type"])
        if root.parent.name == nested_folder and root.parent.parent.name == "diagnostics":
            # Nested result writers retain paths relative to the aggregate root.
            manifest_root = root.parents[2]
        for token in tokens:
            if isinstance(token, str) and not Path(token).is_absolute() and (manifest_root / token).resolve() == path.resolve():
                return payload
    return None


def _inspect_directory(root: Path) -> _Inspection:
    """Read a known layout without recursive discovery or upstream traversal."""
    for relative in ("diagnostics", "hg19", "hg38"):
        child = root / relative
        if child.exists() and not child.resolve().is_relative_to(root):
            raise LDSCInputError(f"Artifact directory {child} points outside {root}. Supply a self-contained artifact.")
    metadata_paths = [root / "metadata.json"] if (root / "metadata.json").is_file() else []
    if not metadata_paths and (root / "diagnostics" / "metadata.json").is_file():
        metadata_paths = [root / "diagnostics" / "metadata.json"]
    if not metadata_paths:
        metadata_paths = sorted((root / "diagnostics").glob("metadata.chr*.json"))
    if metadata_paths:
        result = _Inspection()
        for metadata_path in metadata_paths:
            payload = _read_metadata(_member(root, str(metadata_path.relative_to(root))))
            declared = _inspect_metadata(payload, metadata_path)
            kind = payload["artifact_type"]
            members = list(_file_tokens(payload))
            if kind == "gene_ldscore_index":
                result.merge(declared)
                for chrom in payload.get("chromosomes", []):
                    member = _member(root, f"chromosomes/chr{chrom}/metadata.json")
                    result.merge(_inspect_metadata(_read_metadata(member), member))
            elif members:
                if kind == "ref_panel":
                    # GlobalConfig and source_genome_build do not describe emitted coordinates.
                    for key, token in members:
                        build = key.rsplit("_", 1)[-1]
                        inherited = {"snp_identifier": payload.get("snp_identifier")}
                        if build in declared.groups:
                            inherited["genome_build"] = build
                        result.merge(_inspect_table(_member(root, token), inherited=inherited), group=build)
                    for group in declared.groups:
                        if group and group not in result.groups:
                            raise LDSCInputError(f"Declared {group} output has no coordinate members in {metadata_path}.")
                else:
                    for _, token in members:
                        result.merge(_inspect_table(_member(root, token), inherited=payload))
            else:
                result.merge(declared)
        return result

    # Panel build children are also valid standalone inputs. Directory labels
    # organize alternatives; only file evidence establishes the actual build.
    panel_groups = [("", sorted(root.glob("chr*_r2.parquet")))]
    panel_groups.extend((build, sorted((root / build).glob("chr*_r2.parquet"))) for build in ("hg19", "hg38"))
    result = _Inspection()
    for group, members in panel_groups:
        for member in members:
            result.merge(_inspect_table(_member(root, str(member.relative_to(root)))), group=group)
    if result.groups:
        return result

    # Sumstats directories have no JSON manifest and may use a trait filename.
    for member in sorted(root.glob("*.parquet")):
        candidate = _member(root, member.name)
        try:
            metadata = pq.read_schema(candidate).metadata or {}
        except (OSError, ValueError, ArrowNotImplementedError):
            # Without a manifest, an unreadable neighbor cannot establish ownership.
            continue
        if metadata.get(b"ldsc:artifact_type") == b"sumstats":
            result.merge(_inspect_table(_member(root, member.name)))
    if not result.groups:
        legacy = sorted(root.glob("*.sumstats.gz"))
        if (root / "sumstats.gz").is_file():
            legacy.append(root / "sumstats.gz")
        for member in legacy:
            result.merge(_inspect_table(_member(root, member.name)))
    if result.groups:
        return result
    raise LDSCInputError(f"Unrecognized LDSC artifact directory: {root}. Supply an exact supported file or a recognized artifact directory.")


def _infer_frames(result: _Inspection, path: Path, frames) -> None:
    """Score bounded SNP samples from every chromosome, using existing thresholds.

    Coordinates are streamed to discover all chromosomes; only the first 5,000
    usable rows per chromosome are retained. Reader errors propagate as input
    failures, while insufficient HM3 evidence becomes an unknown result.
    """
    samples, counts = {}, {}
    for frame in frames:
        frame, _ = normalize_chr_pos_frame(
            frame, context=str(path), coordinate_policy="drop", min_position=0, logger=None,
        )
        for chrom, group in frame.groupby("CHR", sort=False):
            room = _SAMPLE_ROWS - counts.get(chrom, 0)
            if room > 0:
                sample = group.head(room)
                samples.setdefault(chrom, []).append(sample)
                counts[chrom] = counts.get(chrom, 0) + len(sample)
    if not samples:
        result.add(path)
    for chrom, parts in samples.items():
        source = f"{path} [chr{chrom}]"
        try:
            inference = infer_chr_pos_build(pd.concat(parts, ignore_index=True), context=source)
        except LDSCInputError:
            result.add(source, detail="insufficient or inconclusive SNP-coordinate evidence")
        else:
            result.add(source, inference.genome_build,
                       f"inferred from SNP coordinates; {inference.coordinate_basis}")
