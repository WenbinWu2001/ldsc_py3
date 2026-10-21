"""Collect workflow build evidence before making one consistency decision.

Concrete user declarations are authoritative. Otherwise each selected input
contributes declared or sampled evidence, and unresolved coordinate inputs
remain errors even when other inputs agree. Numerical readers consume the
resolved build; they do not repeat this workflow-level decision.
"""

from dataclasses import dataclass
import json
import gzip
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from ._coordinates import normalize_chr_pos_frame
from .column_inference import (
    normalize_genome_build, resolve_optional_column, resolve_required_column,
    RESTRICTION_CHRPOS_SPEC_MAP, RESTRICTION_HG19_POS_SPEC, RESTRICTION_HG38_POS_SPEC,
)
from .errors import LDSCInputError, LDSCUserError
from .genome_build_inference import infer_chr_pos_build
from .path_resolution import resolve_file_group, resolve_plink_prefix_group, split_cli_path_tokens

_INPUT_ERRORS = (OSError, ValueError, KeyError, EOFError, LDSCUserError)
_BUILDS = frozenset(("hg19", "hg38"))
_SAMPLE_ROWS = 5000


@dataclass(frozen=True)
class BuildEvidence:
    """One input's declared, inferred, unresolved, or unchecked build evidence."""

    role: str
    source: str
    builds: frozenset[str]
    detail: str
    unchecked: bool = False


class GenomeBuildCheck:
    """Inspect all supplied evidence, then resolve one workflow coordinate build.

    Parameters
    ----------
    requested : {'auto', 'hg19', 'hg38'} or None
        A concrete build bypasses inference and compatibility decisions.
        ``None`` behaves as auto here; command-specific omission rules belong
        to the calling workflow.
    context : str
        Workflow name used in actionable errors.
    option : str or None
        Existing CLI override option, or None when no override is available.
    logger : logging.Logger or None
        Receives the complete grouped report on failure.
    """

    def __init__(self, requested, *, context, option="--genome-build", logger=None):
        self.requested = normalize_genome_build(requested)
        self.context, self.option, self.logger = context, option, logger
        self.evidence = []
        self.inferences = []
        self._generic_fallbacks = []

    @property
    def explicit(self):
        return self.requested in _BUILDS

    def declared(self, role, source, build, detail="metadata declaration"):
        """Record a declaration without re-inferring its coordinates."""
        builds = frozenset(normalize_genome_build(value) for value in
                           (build if isinstance(build, (set, frozenset, list, tuple)) else [build]))
        self.evidence.append(BuildEvidence(role, str(source), builds & _BUILDS, detail))

    def unresolved(self, role, source, reason):
        # The grouped decision owns command-specific repair advice.
        reason = str(reason).split("Pass --genome-build")[0].strip()
        self.evidence.append(BuildEvidence(role, str(source), frozenset(), reason))

    def unchecked(self, role, sources, reason):
        for source in split_cli_path_tokens(str(sources) if isinstance(sources, Path) else sources):
            self.evidence.append(BuildEvidence(role, str(source), frozenset(), reason, True))

    def frames(self, role, source, frames):
        """Sample at most 5,000 rows per chromosome, retaining every chromosome."""
        if self.explicit:
            return
        samples, counts = {}, {}
        try:
            for frame in frames:
                missing = {"CHR", "POS"} - set(frame.columns)
                if missing:
                    raise LDSCInputError(f"missing required coordinate columns: {', '.join(sorted(missing))}; supply a coordinate-bearing input")
                frame, _ = normalize_chr_pos_frame(
                    frame.loc[:, ["CHR", "POS"]], context=str(source),
                    coordinate_policy="drop", min_position=0, logger=None,
                )
                for chrom, group in frame.groupby("CHR", sort=False):
                    room = _SAMPLE_ROWS - counts.get(chrom, 0)
                    if room > 0:
                        sample = group.head(room)
                        samples.setdefault(chrom, []).append(sample)
                        counts[chrom] = counts.get(chrom, 0) + len(sample)
            if not samples:
                self.unresolved(role, source, "no coordinate rows")
            for chrom, parts in samples.items():
                label = f"{source} [chr{chrom}]"
                try:
                    inference = infer_chr_pos_build(pd.concat(parts, ignore_index=True), context=label)
                    self.inferences.append((str(source), inference))
                    self.declared(role, label, inference.genome_build,
                                  f"coordinate inference; {inference.inspected_snp_count} informative HM3 matches; {inference.coordinate_basis}")
                except _INPUT_ERRORS as exc:
                    self.unresolved(role, label, exc)
        except _INPUT_ERRORS as exc:
            self.unresolved(role, source, exc)

    def paths(self, role, tokens, *, bim=False, restriction_mode=None):
        """Resolve exact paths, globs, and @ suites without short-circuiting siblings."""
        if self.explicit:
            return
        if restriction_mode and restriction_mode.startswith("rsid"):
            self.unchecked(role, tokens, "rsID filtering; no coordinate build required")
            return
        seen = set()
        for token in split_cli_path_tokens(str(tokens) if isinstance(tokens, Path) else tokens):
            try:
                paths = resolve_file_group([token], allow_chromosome_suite=True, label=role)
            except _INPUT_ERRORS as exc:
                self.unresolved(role, token, exc)
                continue
            for path in paths:
                if path not in seen:
                    seen.add(path)
                    self.table(role, path, bim=bim)

    def plink(self, role, tokens):
        if self.explicit:
            return
        for token in split_cli_path_tokens(str(tokens) if isinstance(tokens, Path) else tokens):
            try:
                prefixes = resolve_plink_prefix_group([token], allow_chromosome_suite=True)
                for prefix in prefixes:
                    self.table(role, prefix + ".bim", bim=True)
            except _INPUT_ERRORS as exc:
                self.unresolved(role, token, exc)

    def table(self, role, path, *, bim=False):
        """Trust build declarations, otherwise infer from bounded coordinate samples."""
        if self.explicit:
            return
        try:
            if bim:
                with pd.read_csv(path, sep=r"\s+", header=None, usecols=[0, 3],
                                 names=["CHR", "POS"], chunksize=50000) as reader:
                    self.frames(role, path, reader)
                return
            parquet = str(path).endswith(".parquet")
            if parquet:
                self._r2_sidecar_declaration(role, path)
                schema = pq.read_schema(path)
                metadata = schema.metadata or {}
                declarations = [(key, metadata[key].decode()) for key in
                                (b"ldsc:genome_build", b"ldsc:sorted_by_build") if metadata.get(key)]
                if declarations:
                    for key, build in declarations:
                        self.declared(role, path, build, f"Parquet {key.decode()} declaration")
                    return
                columns = schema.names
            else:
                options, text_metadata = _text_options(path)
                build = text_metadata.get("genome_build")
                if build in _BUILDS:
                    self.declared(role, path, build, "text artifact build metadata")
                    return
                columns = list(pd.read_csv(path, nrows=0, **options).columns)
                # Only trust an annotation sidecar that explicitly owns this file.
                sidecar = Path(path).parent / "diagnostics" / "metadata.json"
                if sidecar.is_file():
                    payload = json.loads(sidecar.read_text())
                    owned = payload.get("files", {})
                    members = owned.get("baseline_annotations", []) + owned.get("query_annotations", [])
                    if payload.get("artifact_type") == "annotation_projection" and Path(path).name in members:
                        build = payload.get("projection_genome_build") or payload.get("genome_build")
                        if build in _BUILDS:
                            self.declared(role, path, build, "annotation artifact metadata")
                            return
            declared_column = next((col for col in columns if str(col).lower() == "genome_build"), None)
            if declared_column:
                values = set()
                for frame in _table_frames(path, [declared_column], parquet):
                    values.update(frame[declared_column].dropna().astype(str))
                if not values:
                    self.unresolved(role, path, "empty genome_build declaration")
                for value in sorted(values):
                    self.declared(role, path, value, "genome_build column declaration")
                return
            alternatives = {build for build, spec in (("hg19", RESTRICTION_HG19_POS_SPEC),
                                                       ("hg38", RESTRICTION_HG38_POS_SPEC))
                            if resolve_optional_column(columns, spec, context=str(path)) is not None}
            if alternatives:
                self.declared(role, path, alternatives, "build-specific position columns")
                specific_columns = {resolve_optional_column(columns, spec, context=str(path))
                                    for spec in (RESTRICTION_HG19_POS_SPEC, RESTRICTION_HG38_POS_SPEC)}
                generic = resolve_optional_column([col for col in columns if col not in specific_columns],
                                                  RESTRICTION_CHRPOS_SPEC_MAP["POS"], context=str(path))
                if generic is not None and len(alternatives) == 1:
                    self._generic_fallbacks.append((self.evidence[-1], path, columns, generic, parquet))
                return
            chr_col = resolve_required_column(columns, RESTRICTION_CHRPOS_SPEC_MAP["CHR"], context=str(path))
            pos_col = resolve_required_column(columns, RESTRICTION_CHRPOS_SPEC_MAP["POS"], context=str(path))
            self.frames(role, path, (frame.rename(columns={chr_col: "CHR", pos_col: "POS"})
                                    for frame in _table_frames(path, [chr_col, pos_col], parquet)))
        except _INPUT_ERRORS as exc:
            self.unresolved(role, path, exc)

    def _r2_sidecar_declaration(self, role, path):
        """Compare a bound R2 sidecar's declaration without re-inferring it.

        The owning Parquet declaration establishes the artifact build when
        the sidecar has no separate declaration. Runtime binding checks still
        validate that the sidecar describes the same SNP order and identities.
        """
        path = Path(path)
        if not path.name.endswith("_r2.parquet"):
            return
        sidecar = path.with_name(path.name.removesuffix("_r2.parquet") + "_meta.tsv.gz")
        if sidecar.is_file():
            try:
                build = _text_options(sidecar)[1].get("genome_build")
                if build:
                    self.declared(role, sidecar, build, "R2 sidecar build metadata")
            except _INPUT_ERRORS as exc:
                self.unresolved(role, sidecar, exc)

    def resolve(self):
        """Return the shared build or raise one grouped, actionable input error."""
        if self.explicit:
            if self.logger:
                self.logger.info("Using authoritative %s %s.", self.option or "genome build", self.requested)
            return self.requested
        # Do not infer an unused generic column. Select its fallback only when
        # independent evidence establishes a build without a matching column.
        deferred = {entry for entry, *_ in self._generic_fallbacks}
        known = {next(iter(entry.builds)) for entry in self.evidence
                 if entry not in deferred and not entry.unchecked and len(entry.builds) == 1}
        if len(known) == 1:
            target = next(iter(known))
            for entry, path, columns, pos_col, parquet in self._generic_fallbacks:
                if target not in entry.builds:
                    self.evidence.remove(entry)
                    try:
                        chr_col = resolve_required_column(columns, RESTRICTION_CHRPOS_SPEC_MAP["CHR"], context=str(path))
                        self.frames(entry.role, path, (frame.rename(columns={chr_col: "CHR", pos_col: "POS"})
                            for frame in _table_frames(path, [chr_col, pos_col], parquet)))
                    except _INPUT_ERRORS as exc:
                        self.unresolved(entry.role, path, exc)
        self._generic_fallbacks.clear()
        checked = [entry for entry in self.evidence if not entry.unchecked]
        candidates = set(_BUILDS)
        for entry in checked:
            candidates.intersection_update(entry.builds)
        if checked and len(candidates) == 1:
            build = next(iter(candidates))
            if self.logger:
                self.logger.info("Resolved %s genome build: %s.", self.context, build)
                self.logger.debug(self.report(build))
            return build
        message = f"Genome-build consistency check failed for {self.context}.\n{self.report()}\n"
        known = {next(iter(e.builds)) for e in checked if len(e.builds) == 1}
        if len(known) > 1:
            message += "Conflicting genome-build evidence: inputs identify both hg19 and hg38. "
        if any(not entry.builds for entry in checked) or not checked or not known:
            message += "The genome build could not be inferred or resolved for every required coordinate input. "
        message += "Check the listed inputs and ensure their coordinates use a consistent genome build. "
        if self.option:
            message += (f"To make a known build authoritative, explicitly specify {self.option} hg19 or "
                        f"{self.option} hg38. This declaration does not convert coordinates. ")
        else:
            message += "Regenerate or supply matching upstream artifacts. "
        message += "See docs/troubleshooting.md#common-genome-build-consistency-check-failed."
        if self.logger:
            self.logger.error(message)
        error = LDSCInputError(message)
        error.build_evidence = tuple(self.evidence)
        raise error

    def report(self, resolved=None):
        groups = {key: [] for key in ("hg19", "hg38", "unresolved", "unchecked")}
        for entry in self.evidence:
            group = "unchecked" if entry.unchecked else (
                next(iter(entry.builds)) if len(entry.builds) == 1 else
                resolved if resolved in entry.builds else "unresolved")
            groups[group].append(f"    {entry.role}: {entry.source} ({entry.detail})")
        return "Build evidence across inputs:\n" + "\n".join(
            f"  {group}:\n" + ("\n".join(rows) if rows else "    none") for group, rows in groups.items())


def _table_frames(path, columns, parquet=False):
    if parquet:
        with pq.ParquetFile(path) as reader:
            for batch in reader.iter_batches(batch_size=50000, columns=columns):
                yield batch.to_pandas()
    else:
        with pd.read_csv(path, usecols=columns, chunksize=50000, **_text_options(path)[0]) as reader:
            yield from reader


def _text_options(path):
    """Recognize supported delimiters and leading artifact metadata comments."""
    metadata = {}
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as stream:
        for index, line in enumerate(stream):
            stripped = line.strip()
            if stripped.startswith('# ldsc:'):
                key, _, value = stripped[len('# ldsc:'):].partition('=')
                metadata[key] = value
            if not stripped or (stripped.startswith('#') and not stripped.upper().startswith('#CHROM')):
                continue
            separator = '\t' if '\t' in line else ',' if ',' in line else r'\s+'
            return {'sep': separator, 'skiprows': index}, metadata
    return {'sep': r'\s+'}, metadata
