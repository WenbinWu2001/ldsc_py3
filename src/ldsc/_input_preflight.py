"""Bounded input declaration gates, separate from output ownership checks.

Callers authorize diagnostic destinations through their existing output gate.
This module inspects paths and small headers; scientific content validators
remain in their workflow owners. Results retain concrete paths and widths.
"""

from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
import gzip
import logging

import pandas as pd

from ._annotation_preflight import INPUT_ISSUE_COLUMNS, input_issue, resolve_annotation_inputs
from .errors import LDSCInputError, LDSCUserError
from .path_resolution import _discover_plink_members, resolve_file_group, split_cli_path_tokens, normalize_path_tokens, resolve_exact_file, resolve_scalar_path
from ._progress import PhaseProgress

LOGGER = logging.getLogger('LDSC.preflight')
INPUT_ERRORS = (OSError, EOFError, ValueError, UnicodeError, LDSCUserError)


@dataclass
class InputGate:
    """Collect independent failures and retain one machine-readable repair table."""

    name: str
    issues_path: Path | None = None
    issues: list = field(default_factory=list)

    def check(self, role, source, operation, *, chrom=''):
        """Run an independent check, retaining its value or recording its error."""
        try:
            return operation()
        except INPUT_ERRORS as exc:
            nested = getattr(exc, 'input_issues', None)
            self.issues.extend(nested.to_dict('records') if nested is not None else
                               [dict(input_issue(role, source, chrom, 'invalid_required_input', exc),
                                     repair='Correct the listed path, required companions, or schema; preserve the requested scope. See docs/troubleshooting.md#staged-input-validation.')])
            return None

    def finish(self):
        """Replace this gate's repair status, without changing science files."""
        if not self.issues:
            if self.issues_path is not None:
                self.issues_path.unlink(missing_ok=True)
            return
        frame = pd.DataFrame(self.issues, columns=INPUT_ISSUE_COLUMNS)
        location = ''
        if self.issues_path is not None:
            self.issues_path.parent.mkdir(parents=True, exist_ok=True)
            frame.to_csv(self.issues_path, sep='\t', index=False)
            location = f' Complete diagnostics: {self.issues_path}.'
        details = '; '.join(f"{r['input_role']} {r['source']} {r['chrom']}: {r['details']}" for r in self.issues)
        error = LDSCInputError(f'{self.name} failed: {details}.{location} Repair the listed inputs and rerun. '
                               'Other causes & fixes: docs/troubleshooting.md#staged-input-validation')
        error.input_issues = frame
        raise error


def _readable(path):
    with (gzip.open(path, 'rb') if str(path).endswith('.gz') else open(path, 'rb')) as stream:
        stream.read(1)
    if str(path).endswith('.parquet'):
        import pyarrow.parquet as pq
        pq.read_schema(path)
    return str(path)


def _inspect_file_header(role, path):
    """Inspect bounded generic headers using the existing schema vocabulary."""
    _readable(path)
    if role in {'genetic_map_hg19_sources', 'genetic_map_hg38_sources'}:
        from .column_inference import resolve_required_column, CHR_COLUMN_SPEC, POS_COLUMN_SPEC
        from ._kernel.ref_panel_builder import GENETIC_MAP_CM_SPEC
        columns = pd.read_csv(path, sep=r'\s+', nrows=0).columns
        for spec in (CHR_COLUMN_SPEC, POS_COLUMN_SPEC, GENETIC_MAP_CM_SPEC):
            resolve_required_column(columns, spec, context=str(path))
    return str(path)


def inspect_declared_inputs(*, baseline=(), query=(), files=(), scalar_files=(), plink=None, chromosomes=tuple(map(str, range(1, 23))),
                            require_complete_suite=False, mode='rsid', issues_path=None, strict_annotations=True, initial_issues=(), plink_issues_path=None, checks=()):
    """Check every declared object before large scans and return resolved paths.

    ``files`` contains (role, tokens) pairs using the existing group resolver;
    ``scalar_files`` uses the existing exact-one resolver instead.
    Annotation ``@`` completeness remains a workflow-specific policy. Reading
    one annotation row verifies headers and records width without a full scan.
    PLINK inspection here checks companion presence only; BIM/FAM and BED-size
    validation is the next content gate. No numeric reference data are loaded.
    """
    gate = InputGate('Input preflight (declarations and headers)', issues_path, list(initial_issues))
    resolved, declarations, widths, columns = {}, {}, {}, {}
    with PhaseProgress(LOGGER, 'validation', 'input declarations and headers') as progress:
        for role, tokens in (('baseline', baseline), ('query', query)):
            if strict_annotations:
                paths, declared, issues = resolve_annotation_inputs(tokens, role=role)
                declarations.update(declared)
                gate.issues.extend(issues)
            else:
                paths = []
                for token in split_cli_path_tokens(normalize_path_tokens(tokens)):
                    matches = gate.check(role, token, lambda token=token: resolve_file_group([token], allow_chromosome_suite=True))
                    paths.extend(matches or [])
            resolved[role] = paths
            for path in paths:
                def header():
                    from ._annotation_parsing import normalize_annotation_chunk
                    frame = pd.read_csv(path, sep=r'\s+', nrows=1)
                    if frame.empty:
                        raise LDSCInputError(f"Annotation file '{path}' has no SNP rows. Supply a nonempty annotation grid.")
                    metadata, values = normalize_annotation_chunk(frame, path, mode, log_ignored_metadata=False)
                    member = declarations.get(path)
                    if member and set(metadata.CHR.astype(str)) != {member}:
                        raise LDSCInputError(f'@ member for chromosome {member} starts with a different chromosome.')
                    columns[str(path)] = tuple(values.columns)
                    return len(frame.columns)
                width = gate.check(role, path, header, chrom=declarations.get(path, ''))
                if width is not None:
                    widths[str(path)] = width
                progress.advance(object=str(path))
            declared_headers = [(path, columns[path]) for path in paths if path in declarations and path in columns]
            if declared_headers:
                expected = declared_headers[0][1]
                for path, names in declared_headers[1:]:
                    if names != expected:
                        gate.issues.append(input_issue(role, path, declarations[path], 'invalid_required_input',
                            'Annotation columns differ across declared chromosomes. Use the same annotation header for every shard.'))
        for role, tokens, scalar in [*((role, tokens, False) for role, tokens in files),
                                     *((role, tokens, True) for role, tokens in scalar_files)]:
            paths = []
            for token in (normalize_path_tokens(tokens) if scalar else split_cli_path_tokens(normalize_path_tokens(tokens))):
                if role in {'control genes', 'control gene-list file', 'control_gene_list_file'}:
                    match = gate.check('control gene-list file', token,
                        lambda token=token: resolve_exact_file(token, label='control gene-list file'))
                    matches = [match] if match is not None else []
                elif scalar:
                    match = gate.check(role, token, lambda token=token: resolve_scalar_path(token, label=role))
                    matches = [match] if match is not None else []
                else:
                    matches = gate.check(role, token, lambda token=token: resolve_file_group([token], label=role, allow_chromosome_suite=True))
                for path in matches or []:
                    gate.check(role, path, lambda path=path: _inspect_file_header(role, path))
                    paths.append(path)
                    progress.advance(object=str(path))
            resolved[role] = paths
        if plink is not None:
            members, issues = _discover_plink_members(plink, chromosomes=chromosomes,
                allow_chromosome_suite=True, require_complete_suite=require_complete_suite)
            gate.issues.extend(issues)
            if issues and plink_issues_path is not None:
                plink_issues_path.parent.mkdir(parents=True, exist_ok=True)
                pd.DataFrame(issues, columns=INPUT_ISSUE_COLUMNS).to_csv(plink_issues_path, sep='\t', index=False)
            resolved['plink_members'] = members
            progress.advance(object='PLINK companion files')
        for role, source, operation in checks:
            resolved[role] = gate.check(role, source, operation)
        gate.finish()
    return SimpleNamespace(baseline=resolved['baseline'], query=resolved['query'],
                           declarations=declarations, widths=widths, columns=columns, files=resolved)


def inspect_r2_headers(directory, genome_build=None):
    """Check every selected R2/sidecar path and Parquet footer before pair scans."""
    import pyarrow.parquet as pq
    from ._kernel.ref_panel import _resolve_r2_build_dir

    root = _resolve_r2_build_dir(directory, genome_build)
    members = sorted({path.name.split('_')[0][3:] for pattern in ('chr*_r2.parquet', 'chr*_meta.tsv.gz') for path in root.glob(pattern)})
    gate = InputGate('R2 declarations and headers')
    if not members:
        gate.issues.append(input_issue('reference', root, '', 'missing_required_input', 'No R2/sidecar chromosome artifacts found.'))
    for chrom in members:
        pair = root / f'chr{chrom}_r2.parquet'
        sidecar = root / f'chr{chrom}_meta.tsv.gz'
        gate.check('reference', sidecar, lambda: _readable(sidecar), chrom=chrom)
        gate.check('reference', pair, lambda: pq.read_schema(pair), chrom=chrom)
    gate.finish()
    return root


def inspect_index_paths(directory, *, metadata=None):
    """Inspect root declarations and all component paths without loading arrays."""
    root = Path(directory)
    if metadata is None:
        from .gene_ldscore_index import _read_json
        metadata = _read_json(root/'metadata.json', 'root')
    chromosomes = metadata.get('chromosomes', [])
    if not isinstance(chromosomes, (list, tuple)) or not chromosomes or any(str(chrom) not in tuple(map(str, range(1, 23))) for chrom in chromosomes):
        raise LDSCInputError('Gene LD-score index chromosome coverage is invalid.')
    gate = InputGate('Index component declarations')
    paths = [('gene_catalog.parquet', '')]
    files = ('metadata.json', 'baseline_rows.parquet', 'baseline_statistics.npz', 'atoms.parquet',
             'gene_to_atom.npz', 'ldscore_operator.npz', 'atom_statistics.npz')
    paths.extend((f'chromosomes/chr{chrom}/{name}', str(chrom)) for chrom in chromosomes for name in files)
    for name, chrom in paths:
        path = root/name
        if not path.is_file():
            gate.issues.append(input_issue('index', path, chrom, 'missing_required_input', 'Required index component is missing.'))
        elif name.endswith('.parquet'):
            gate.check('index', path, lambda path=path: _readable(path), chrom=chrom)
        elif chrom and name.endswith('metadata.json'):
            from .gene_ldscore_index import _read_json, _validate_chromosome_metadata
            def component_header():
                component = _read_json(path, f'chromosome {chrom}')
                _validate_chromosome_metadata(component, chrom=chrom, index_id=metadata.get('index_id'),
                    snp_identifier=metadata.get('snp_identifier'), genome_build=metadata.get('genome_build'))
            gate.check('index', path, component_header, chrom=chrom)
    gate.finish()
    return metadata


def inspect_artifact_paths(directory, *, metadata_file="metadata.json"):
    """Check files explicitly listed by a saved artifact's metadata manifest."""
    import json
    root = Path(directory)
    with (root/metadata_file).open() as stream:
        metadata = json.load(stream)
    if not isinstance(metadata, dict) or not isinstance(metadata.get('files', {}), dict):
        raise LDSCInputError(f'Saved artifact metadata must contain an object with a files mapping: {root/metadata_file}')
    gate = InputGate('Saved artifact declarations')
    for role, member in metadata.get('files', {}).items():
        for name in ([member] if isinstance(member, str) else member if isinstance(member, list) else []):
            if isinstance(name, str):
                path = root/name
                gate.check(role, path, lambda path=path: _readable(path))
    gate.finish()
    return metadata
