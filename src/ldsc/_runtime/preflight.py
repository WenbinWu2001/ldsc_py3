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

from ..errors import LDSCInputError, LDSCUserError
from .paths import _discover_plink_members, resolve_file_group, split_cli_path_tokens, normalize_path_tokens, resolve_exact_file, resolve_scalar_path
from .progress import PhaseProgress

LOGGER = logging.getLogger('LDSC.preflight')
INPUT_ERRORS = (OSError, EOFError, ValueError, UnicodeError, LDSCUserError)
INPUT_ISSUE_COLUMNS = ['input_role', 'source', 'chrom', 'reason', 'details', 'repair']


def input_issue(role, path, chrom, reason, details):
    return dict(input_role=role, source=str(path), chrom=chrom, reason=reason, details=str(details),
                repair='Supply valid annotation inputs for the declared scope; inspect exact paths, glob matches, and @ declarations.')


@dataclass
class InputGate:
    """Collect independent input failures with diagnostic rows and original causes."""

    name: str
    issues_path: Path | None = None
    issues: list = field(default_factory=list)
    causes: list[Exception] = field(default_factory=list, repr=False)
    troubleshooting: str = 'docs/user-guide/troubleshooting.md#staged-input-validation'

    def check(self, role, source, operation, *, chrom=''):
        """Retain a check's value or record its typed error, repair, and original cause."""
        try:
            return operation()
        except INPUT_ERRORS as exc:
            self.causes.append(exc)
            nested = getattr(exc, 'input_issues', None)
            if nested is not None:
                self.issues.extend(nested.to_dict('records'))
            else:
                cause = exc.__cause__ or exc
                reason = 'invalid_required_input'
                details = f'{type(exc).__name__}: {exc}'
                repair = 'Correct the listed path, required companions, or schema; preserve the requested scope.'
                if isinstance(cause, FileNotFoundError):
                    reason = 'missing_required_input'
                    kind = 'FileNotFoundError' + (f' (errno {cause.errno})' if cause.errno is not None else '')
                    details = (f'{kind}: path resolution found no matching file on this node.'
                               if isinstance(exc, LDSCInputError) else
                               f'{kind}: the named input is missing or not visible on this node.')
                    repair = 'Verify the filename, location, and visibility from the compute node.'
                    if role == 'raw summary statistics':
                        if isinstance(exc, LDSCInputError):
                            details += ' Summary-statistic reading has not started.'
                        repair = 'Set --raw-sumstats-file to an existing file. ' + repair
                elif isinstance(cause, OSError) and cause.errno is not None:
                    reason = 'unreadable_required_input'
                    details = f'{type(cause).__name__} (errno {cause.errno}): {cause.strerror}'
                    repair = 'Check the named input and its access permissions from the compute node.'
                self.issues.append(dict(input_issue(role, source, chrom, reason, details), repair=repair))
            return None

    def finish(self):
        """Publish repair status and raise with all captured original causes."""
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
        displayed = self.issues[:10] if location else self.issues
        details = '; '.join(
            f"{r['input_role']} {r['source']}"
            f"{(' (chromosome ' + str(r['chrom']) + ')') if r['chrom'] else ''}: "
            f"{r['details'].rstrip('.')}. Repair: {r['repair'].rstrip('.')}"
            for r in displayed
        )
        count = len(self.issues)
        omitted = count - len(displayed)
        if omitted:
            details += f'. {omitted} more input issues are listed in the complete repair table'
        error = LDSCInputError(f"{self.name} failed: {count} input issue{'s' if count != 1 else ''}. "
                               f'{details}.{location} '
                               f'Other causes & fixes: {self.troubleshooting}')
        error.input_issues = frame
        if len(self.causes) == 1:
            raise error from self.causes[0]
        if self.causes:
            raise error from ExceptionGroup('Input validation failures', self.causes)
        raise error


def check_readable(path):
    with (gzip.open(path, 'rb') if str(path).endswith('.gz') else open(path, 'rb')) as stream:
        stream.read(1)
    if str(path).endswith('.parquet'):
        import pyarrow.parquet as pq
        pq.read_schema(path)
    return str(path)


def _inspect_file_header(role, path, header_checks):
    """Check that a declared file is readable, then run its role's header check, if any."""
    check_readable(path)
    check = header_checks.get(role)
    if check is not None:
        check(path)
    return str(path)


def inspect_declared_inputs(*, files=(), scalar_files=(), plink=None, chromosomes=tuple(map(str, range(1, 23))),
                            require_complete_suite=False, issues_path=None, initial_issues=(), plink_issues_path=None,
                            checks=(), declaration_check=None, header_checks=None):
    """Check every declared object before large scans and return resolved paths.

    ``files`` contains (role, tokens) pairs using the existing group resolver;
    ``scalar_files`` uses the existing exact-one resolver instead.
    ``declaration_check(gate, progress, resolved)`` runs first so a domain can
    record its own declarations, such as annotation headers, on the same gate;
    it fills ``resolved`` and returns extra result fields. ``header_checks``
    maps a file role to a header check that runs once the file is readable.
    PLINK inspection here checks companion presence only; BIM/FAM and BED-size
    validation is the next content gate. No numeric reference data are loaded.
    """
    gate = InputGate('Input preflight (declarations and headers)', issues_path, list(initial_issues))
    resolved = {'baseline': [], 'query': []}
    extra = dict(declarations={}, widths={}, columns={})
    header_checks = header_checks or {}
    with PhaseProgress(LOGGER, 'validation', 'input declarations and headers') as progress:
        if declaration_check is not None:
            extra.update(declaration_check(gate, progress, resolved))
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
                    gate.check(role, path, lambda path=path: _inspect_file_header(role, path, header_checks))
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
        if gate.issues:
            progress.outcomes['failed_inputs'] = len({
                (issue['input_role'], issue['source'], issue['chrom']) for issue in gate.issues
            })
        gate.finish()
    return SimpleNamespace(baseline=resolved['baseline'], query=resolved['query'], files=resolved, **extra)


def inspect_artifact_paths(directory, *, metadata_file="metadata.json"):
    """Check files explicitly listed by a saved artifact's metadata manifest."""
    import json
    root = Path(directory)
    with (root/metadata_file).open() as stream:
        metadata = json.load(stream)
    if not isinstance(metadata, dict) or not isinstance(metadata.get('files', {}), dict):
        raise LDSCInputError(
            f'Saved artifact metadata must contain an object with a files mapping: {root/metadata_file}'
            " Restore this saved artifact directory with its original metadata and data from the same run, or "
            "regenerate the upstream artifact."
        )
    gate = InputGate('Saved artifact declarations')
    for role, member in metadata.get('files', {}).items():
        for name in ([member] if isinstance(member, str) else member if isinstance(member, list) else []):
            if isinstance(name, str):
                path = root/name
                gate.check(role, path, lambda path=path: check_readable(path))
    gate.finish()
    return metadata
