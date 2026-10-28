"""Resolve annotation declarations while collecting every missing input."""

import pandas as pd

from .parsing import normalize_annotation_chunk
from .._runtime.preflight import input_issue, inspect_declared_inputs
from .._runtime.paths import (
    ANNOTATION_SUFFIXES, normalize_path_token, normalize_path_tokens, resolve_file_group, split_cli_path_tokens,
)
from ..errors import LDSCInputError, LDSCUserError

AUTOSOMES = tuple(map(str, range(1, 23)))


def resolve_annotation_inputs(tokens, *, role='baseline', suffixes=ANNOTATION_SUFFIXES):
    """Return concrete paths, declared chromosomes, and accumulated path issues."""
    paths, declared, issues = [], {}, []
    for raw in split_cli_path_tokens(tokens):
        token = normalize_path_token(raw)
        members = [(token.replace('@', chrom), chrom) for chrom in AUTOSOMES] if '@' in token else [(token, '')]
        for path, chrom in members:
            try:
                matches = resolve_file_group([path], suffixes=suffixes, label=f'{role} annotation')
                for match in matches:
                    if match not in paths:
                        paths.append(match)
                    if chrom:
                        declared[match] = chrom
            except (OSError, LDSCUserError) as exc:
                issues.append(input_issue(role, path, chrom, 'missing_required_input', exc))
    return paths, declared, issues


def inspect_annotation_inputs(*, baseline=(), query=(), mode='rsid', strict_annotations=True, **kwargs):
    """Run the declared-input gate with baseline and query annotation checks first.

    Annotation ``@`` completeness remains a workflow-specific policy. Reading
    one annotation row verifies headers and records width without a full scan.
    The other keyword arguments go to ``inspect_declared_inputs``, whose result
    gains the annotation ``declarations``, ``widths`` and ``columns``.
    """
    def check(gate, progress, resolved):
        declarations, widths, columns = {}, {}, {}
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
                    frame = pd.read_csv(path, sep=r'\s+', nrows=1)
                    if frame.empty:
                        raise LDSCInputError(f"Annotation file '{path}' has no SNP rows. Supply a nonempty annotation grid.")
                    metadata, values = normalize_annotation_chunk(frame, path, mode, log_ignored_metadata=False)
                    member = declarations.get(path)
                    if member and set(metadata.CHR.astype(str)) != {member}:
                        raise LDSCInputError(
                            f"Annotation file '{path}' selected by @ for chromosome {member} starts with "
                            f"chromosome {metadata.CHR.iloc[0]}. "
                            "Restore the matching chromosome file or correct the path template; "
                            "renaming a file does not change its SNP coordinates."
                        )
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
        return dict(declarations=declarations, widths=widths, columns=columns)

    return inspect_declared_inputs(declaration_check=check, **kwargs)
