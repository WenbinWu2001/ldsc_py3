"""Resolve annotation declarations while collecting every missing input."""

from .path_resolution import ANNOTATION_SUFFIXES, resolve_file_group, split_cli_path_tokens, normalize_path_token
from .errors import LDSCUserError

AUTOSOMES = tuple(map(str, range(1, 23)))
INPUT_ISSUE_COLUMNS = ['input_role', 'source', 'chrom', 'reason', 'details', 'repair']


def input_issue(role, path, chrom, reason, details):
    return dict(input_role=role, source=str(path), chrom=chrom, reason=reason, details=str(details),
                repair='Supply valid annotation inputs for the declared scope; inspect exact paths, glob matches, and @ declarations.')


def resolve_annotation_inputs(tokens, *, role='baseline'):
    """Return concrete paths, declared chromosomes, and accumulated path issues."""
    paths, declared, issues = [], {}, []
    for raw in split_cli_path_tokens(tokens):
        token = normalize_path_token(raw)
        members = [(token.replace('@', chrom), chrom) for chrom in AUTOSOMES] if '@' in token else [(token, '')]
        for path, chrom in members:
            try:
                matches = resolve_file_group([path], suffixes=ANNOTATION_SUFFIXES, label=f'{role} annotation')
                for match in matches:
                    if match not in paths:
                        paths.append(match)
                    if chrom:
                        declared[match] = chrom
            except (OSError, LDSCUserError) as exc:
                issues.append(input_issue(role, path, chrom, 'missing_required_input', exc))
    return paths, declared, issues
