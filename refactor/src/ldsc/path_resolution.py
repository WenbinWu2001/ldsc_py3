"""Workflow-layer path-token normalization and resolution helpers.

Core functionality:
    Normalize user-facing path tokens and resolve them into concrete inputs
    before they cross into the primitive-only kernel layer.

Overview
--------
This module defines the public token language for filesystem inputs in the
refactored package. Callers may pass exact paths, standard Python glob
patterns, legacy chromosome-suite placeholders using ``@``, or bare prefixes
that imply chromosome substitution in suite-capable contexts. Resolution is a
workflow-layer concern; the execution kernels should only ever receive concrete
string paths or prefixes.
"""

from __future__ import annotations

import glob
import os
from collections.abc import Sequence
from os import PathLike
from typing import Optional, Union


InputPathToken = Union[str, PathLike[str]]
InputPathCollection = Union[InputPathToken, Sequence[InputPathToken]]

ANNOTATION_SUFFIXES = ("", ".annot.gz", ".annot", ".txt.gz", ".txt", ".tsv.gz", ".tsv")
PARQUET_SUFFIXES = ("", ".parquet")
FREQUENCY_SUFFIXES = ("", ".frq.gz", ".frq", ".txt.gz", ".txt", ".tsv.gz", ".tsv")
STANDARD_CHROMOSOMES = tuple([str(i) for i in range(1, 23)] + ["X", "Y", "MT", "M"])


def normalize_path_token(path: InputPathToken) -> str:
    """Normalize one public path token without resolving it on disk."""
    return os.path.expanduser(os.path.expandvars(str(path)))


def normalize_optional_path_token(path: Optional[InputPathToken]) -> Optional[str]:
    """Normalize an optional path token while preserving ``None``."""
    if path is None:
        return None
    return normalize_path_token(path)


def normalize_path_tokens(values: Optional[InputPathCollection]) -> tuple[str, ...]:
    """Normalize one or many public path tokens into a tuple of strings."""
    if values is None:
        return ()
    if isinstance(values, (str, os.PathLike)):
        return (normalize_path_token(values),)
    return tuple(normalize_path_token(value) for value in values)


def split_cli_path_tokens(value: Optional[Union[str, Sequence[str]]]) -> list[str]:
    """Split CLI-style comma-delimited tokens into normalized path strings."""
    if value is None:
        return []
    if isinstance(value, str):
        raw_values = value.split(",")
    else:
        raw_values = []
        for item in value:
            raw_values.extend(str(item).split(","))
    tokens: list[str] = []
    for raw in raw_values:
        token = normalize_path_token(raw).strip()
        if token:
            tokens.append(token)
    return tokens


def substitute_chromosome(token: str, chrom: str) -> str:
    """Apply the legacy ``@`` chromosome substitution rule to ``token``."""
    if "@" not in token:
        token += "@"
    return token.replace("@", str(chrom))


def resolve_scalar_path(
    token: InputPathToken,
    *,
    suffixes: Sequence[str] = (),
    label: str = "input",
) -> str:
    """Resolve a token that must correspond to exactly one concrete path."""
    matches = _resolve_direct_token(normalize_path_token(token), suffixes=suffixes)
    if len(matches) != 1:
        raise ValueError(f"{label} must resolve to exactly one path; got {len(matches)} matches for {token!r}.")
    return matches[0]


def resolve_file_group(
    tokens: InputPathCollection | None,
    *,
    suffixes: Sequence[str] = (),
    label: str = "input",
    allow_chromosome_suite: bool = False,
    chromosomes: Sequence[str] = STANDARD_CHROMOSOMES,
) -> list[str]:
    """Resolve one or many tokens into a deterministic concrete file group."""
    normalized = normalize_path_tokens(tokens)
    resolved: list[str] = []
    for token in normalized:
        try:
            matches = _resolve_direct_token(token, suffixes=suffixes)
        except FileNotFoundError:
            if not allow_chromosome_suite:
                raise FileNotFoundError(f"Could not resolve {label} token: {token}") from None
            matches = _resolve_suite_token(token, suffixes=suffixes, chromosomes=chromosomes)
            if not matches:
                raise FileNotFoundError(f"Could not resolve {label} token: {token}") from None
        resolved.extend(matches)
    resolved = _dedupe_preserving_order(resolved)
    if not resolved:
        raise FileNotFoundError(f"No {label} files were resolved from the supplied tokens.")
    return resolved


def resolve_chromosome_group(
    tokens: InputPathCollection | None,
    *,
    chrom: str,
    suffixes: Sequence[str] = (),
    label: str = "input",
    required: bool = True,
) -> list[str]:
    """Resolve tokens for one chromosome, supporting legacy suite placeholders."""
    normalized = normalize_path_tokens(tokens)
    resolved: list[str] = []
    for token in normalized:
        try:
            matches = _resolve_direct_token(token, suffixes=suffixes)
        except FileNotFoundError:
            chrom_token = substitute_chromosome(token, chrom)
            try:
                matches = _resolve_direct_token(chrom_token, suffixes=suffixes)
            except FileNotFoundError:
                if required:
                    raise FileNotFoundError(
                        f"Could not resolve {label} token {token!r} for chromosome {chrom}."
                    ) from None
                continue
        resolved.extend(matches)
    resolved = _dedupe_preserving_order(resolved)
    if required and not resolved:
        raise FileNotFoundError(f"No {label} files were resolved for chromosome {chrom}.")
    return resolved


def resolve_plink_prefix(
    token: InputPathToken,
    *,
    chrom: Optional[str] = None,
) -> str:
    """Resolve a PLINK prefix token to exactly one concrete prefix."""
    matches = resolve_plink_prefix_group((token,), chrom=chrom, allow_chromosome_suite=(chrom is None))
    if len(matches) != 1:
        raise ValueError(f"PLINK prefix must resolve to exactly one path; got {len(matches)} matches for {token!r}.")
    return matches[0]


def resolve_plink_prefix_group(
    tokens: InputPathCollection | None,
    *,
    chrom: Optional[str] = None,
    allow_chromosome_suite: bool = False,
    chromosomes: Sequence[str] = STANDARD_CHROMOSOMES,
) -> list[str]:
    """Resolve one or many PLINK prefix tokens into concrete prefixes."""
    normalized = normalize_path_tokens(tokens)
    resolved: list[str] = []
    for token in normalized:
        try:
            matches = _resolve_direct_plink_token(token)
        except FileNotFoundError:
            if chrom is not None:
                chrom_token = substitute_chromosome(token, chrom)
                try:
                    matches = _resolve_direct_plink_token(chrom_token)
                except FileNotFoundError:
                    raise FileNotFoundError(
                        f"Could not resolve PLINK prefix token {token!r} for chromosome {chrom}."
                    ) from None
            elif allow_chromosome_suite:
                matches = _resolve_plink_suite_token(token, chromosomes=chromosomes)
                if not matches:
                    raise FileNotFoundError(f"Could not resolve PLINK prefix token: {token}") from None
            else:
                raise FileNotFoundError(f"Could not resolve PLINK prefix token: {token}") from None
        resolved.extend(matches)
    resolved = _dedupe_preserving_order(resolved)
    if not resolved:
        raise FileNotFoundError("No PLINK prefixes were resolved from the supplied tokens.")
    return resolved


def _resolve_direct_token(token: str, *, suffixes: Sequence[str]) -> list[str]:
    if os.path.exists(token):
        return [token]

    matches = sorted({normalize_path_token(path) for path in glob.glob(token)})
    if matches:
        return matches

    for suffix in suffixes:
        candidate = token + suffix
        if os.path.exists(candidate):
            return [candidate]

    raise FileNotFoundError(token)


def _resolve_suite_token(token: str, *, suffixes: Sequence[str], chromosomes: Sequence[str]) -> list[str]:
    resolved: list[str] = []
    for chrom in chromosomes:
        chrom_token = substitute_chromosome(token, chrom)
        try:
            resolved.extend(_resolve_direct_token(chrom_token, suffixes=suffixes))
        except FileNotFoundError:
            continue
    return _dedupe_preserving_order(resolved)


def _resolve_direct_plink_token(token: str) -> list[str]:
    suffixes = (".bed", ".bim", ".fam")
    for suffix in suffixes:
        if token.endswith(suffix) and os.path.exists(token):
            return [token[: -len(suffix)]]
    for suffix in suffixes:
        if os.path.exists(token + suffix):
            return [token]

    patterns: list[str] = []
    if glob.has_magic(token):
        if any(token.endswith(suffix) for suffix in suffixes):
            patterns.append(token)
        else:
            patterns.extend(token + suffix for suffix in suffixes)
    matches = sorted(
        {
            normalize_path_token(path[: -len(suffix)])
            for pattern in patterns
            for path in glob.glob(pattern)
            for suffix in suffixes
            if path.endswith(suffix)
        }
    )
    if matches:
        return matches

    raise FileNotFoundError(token)


def _resolve_plink_suite_token(token: str, *, chromosomes: Sequence[str]) -> list[str]:
    resolved: list[str] = []
    for chrom in chromosomes:
        try:
            resolved.extend(_resolve_direct_plink_token(substitute_chromosome(token, chrom)))
        except FileNotFoundError:
            continue
    return _dedupe_preserving_order(resolved)


def _dedupe_preserving_order(values: Sequence[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered
