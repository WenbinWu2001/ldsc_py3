"""Workflow-layer path-token normalization and resolution helpers.

Core functionality:
    Normalize user-facing path tokens and resolve them into concrete inputs
    before they cross into the primitive-only kernel layer.

Overview
--------
This module defines the public token language for filesystem inputs in the
refactored package. Callers may pass exact paths, standard Python glob
patterns, or explicit chromosome-suite placeholders using ``@``. Resolution is
a workflow-layer concern; the execution kernels should only ever receive
concrete string paths or prefixes.

Output helpers in this module intentionally use a stricter policy than input
resolution: output directories are literal destinations, missing directories
are created, existing directories are reused, and fixed output files are
preflighted before workflow code writes. Existing files raise by default unless
the caller passes an explicit overwrite flag. Workflows with fixed artifact
families can also preflight the complete owned family and remove stale siblings
after a successful overwrite.
"""

from __future__ import annotations

import glob
import os
import re
import shutil
import warnings
from collections.abc import Iterable, Sequence
from os import PathLike
from pathlib import Path
from typing import Optional, Union

from .chromosome_inference import STANDARD_CHROMOSOMES

InputPathToken = Union[str, PathLike[str]]
InputPathCollection = Union[InputPathToken, Sequence[InputPathToken]]

ANNOTATION_SUFFIXES = ("", ".annot.gz", ".annot", ".txt.gz", ".txt", ".tsv.gz", ".tsv")
PARQUET_SUFFIXES = ("", ".parquet")
FREQUENCY_SUFFIXES = ("", ".frq.gz", ".frq", ".txt.gz", ".txt", ".tsv.gz", ".tsv")


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
        raise FileNotFoundError(token)
    return token.replace("@", str(chrom))


def resolve_scalar_path(
    token: InputPathToken,
    *,
    suffixes: Sequence[str] = (),
    label: str = "input",
) -> str:
    """Resolve a token that must correspond to exactly one concrete path.

    Parameters
    ----------
    token : str or os.PathLike[str]
        Exact path or exact-one glob token.
    suffixes : sequence of str, optional
        Accepted suffixes for this input family. The current strict-token
        contract does not append these suffixes during resolution; they are
        documented here only so callers can pass the same suffix family used by
        the group resolvers. Default is ``()``.
    label : str, optional
        Public argument label used in error messages. Default is ``"input"``.

    Returns
    -------
    str
        The one resolved filesystem path.

    Raises
    ------
    FileNotFoundError
        If ``token`` resolves to no paths.
    ValueError
        If ``token`` resolves to more than one path.
    """
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
    """Resolve one or many path tokens into a deterministic concrete file group.

    This function is the workflow-layer entry point for group-style inputs such
    as annotation files, parquet R2 tables, BED files, and frequency metadata.
    Each token may be an exact path, a standard Python glob, or, when
    ``allow_chromosome_suite`` is true, an explicit ``@`` chromosome-suite
    token. All matches across all tokens are concatenated, deduplicated while
    preserving first-seen order, and returned as one concrete file list.

    Parameters
    ----------
    tokens : path token or sequence of path tokens, optional
        Public user-facing tokens to resolve. ``None`` means no tokens were
        supplied.
    suffixes : sequence of str, optional
        Accepted suffix family for this input class. Resolution remains strict:
        the function does not guess or append missing suffixes. Default is
        ``()``.
    label : str, optional
        Public argument label used in error messages. Default is ``"input"``.
    allow_chromosome_suite : bool, optional
        Whether explicit ``@`` suite tokens are allowed. Default is ``False``.
    chromosomes : sequence of str, optional
        Chromosome labels used when expanding ``@`` suite tokens. Default is
        ``STANDARD_CHROMOSOMES``.

    Returns
    -------
    list of str
        Concrete file paths in deterministic order.

    Raises
    ------
    FileNotFoundError
        If any token resolves to nothing, or if the full group resolves to no
        files.
    """
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
    """Resolve group-style tokens for one chromosome.

    The resolver first tries to interpret each token directly. If a direct
    token expands to multiple files, the function attempts filename-based
    chromosome filtering via :func:`filter_paths_for_chromosome`. If that
    filtering finds at least one match, only those chromosome-matching files
    are retained for this chromosome pass. Otherwise the original matches are
    kept, which allows multi-chromosome files to participate and be filtered by
    ``CHR`` later in the file parser. When a token cannot be resolved directly,
    an explicit ``@`` substitution is attempted.

    Parameters
    ----------
    tokens : path token or sequence of path tokens, optional
        Public user-facing tokens to resolve.
    chrom : str
        Chromosome label for the current processing pass.
    suffixes : sequence of str, optional
        Accepted suffix family for this input class. Resolution remains strict:
        the function does not guess or append missing suffixes. Default is
        ``()``.
    label : str, optional
        Public argument label used in error messages. Default is ``"input"``.
    required : bool, optional
        Whether missing matches should raise immediately. Default is ``True``.

    Returns
    -------
    list of str
        Concrete file paths to use for chromosome ``chrom``.

    Raises
    ------
    FileNotFoundError
        If ``required`` is true and no files can be resolved for ``chrom``.
    """
    normalized = normalize_path_tokens(tokens)
    resolved: list[str] = []
    for token in normalized:
        try:
            matches = _resolve_direct_token(token, suffixes=suffixes)
            if len(matches) > 1:
                chrom_matches = filter_paths_for_chromosome(matches, chrom)
                if chrom_matches:
                    matches = chrom_matches
        except FileNotFoundError:
            try:
                chrom_token = substitute_chromosome(token, chrom)
            except FileNotFoundError:
                if required:
                    raise FileNotFoundError(
                        f"Could not resolve {label} token {token!r} for chromosome {chrom}."
                    ) from None
                continue
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
    """Resolve a PLINK prefix token to exactly one concrete prefix.

    Parameters
    ----------
    token : str or os.PathLike[str]
        Exact PLINK prefix, exact `.bed/.bim/.fam` path, exact-one glob, or an
        explicit ``@`` chromosome-suite token when ``chrom`` is supplied.
    chrom : str or None, optional
        Chromosome label for per-chromosome resolution. Default is ``None``.

    Returns
    -------
    str
        One resolved PLINK prefix without the file extension.

    Raises
    ------
    FileNotFoundError
        If ``token`` cannot be resolved.
    ValueError
        If ``token`` resolves to more than one prefix.
    """
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
    """Resolve one or many PLINK prefix tokens into concrete prefixes.

    This is the prefix-level analogue of :func:`resolve_file_group`. Each token
    is resolved against complete PLINK trios rather than individual files. When
    ``chrom`` is supplied and a direct token expands to multiple prefixes, the
    function first tries filename-based chromosome filtering and otherwise keeps
    the unmatched prefixes as-is.
    """
    normalized = normalize_path_tokens(tokens)
    resolved: list[str] = []
    for token in normalized:
        try:
            matches = _resolve_direct_plink_token(token)
            if chrom is not None and len(matches) > 1:
                chrom_matches = filter_paths_for_chromosome(matches, chrom)
                if chrom_matches:
                    matches = chrom_matches
        except FileNotFoundError:
            if chrom is not None:
                try:
                    chrom_token = substitute_chromosome(token, chrom)
                except FileNotFoundError:
                    raise FileNotFoundError(
                        f"Could not resolve PLINK prefix token {token!r} for chromosome {chrom}."
                    ) from None
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
    """Resolve one token through direct paths and standard globs."""
    if os.path.exists(token):
        return [token]

    matches = sorted({normalize_path_token(path) for path in glob.glob(token)})
    if matches:
        return matches

    raise FileNotFoundError(token)


def _resolve_suite_token(token: str, *, suffixes: Sequence[str], chromosomes: Sequence[str]) -> list[str]:
    """Resolve a chromosome-suite token by substituting each chromosome label."""
    resolved: list[str] = []
    for chrom in chromosomes:
        chrom_token = substitute_chromosome(token, chrom)
        try:
            resolved.extend(_resolve_direct_token(chrom_token, suffixes=suffixes))
        except FileNotFoundError:
            continue
    return _dedupe_preserving_order(resolved)


def _resolve_direct_plink_token(token: str) -> list[str]:
    """Resolve one PLINK token to a concrete prefix with `.bed/.bim/.fam` files."""
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
    """Resolve a PLINK chromosome-suite token across the supplied chromosomes."""
    resolved: list[str] = []
    for chrom in chromosomes:
        try:
            resolved.extend(_resolve_direct_plink_token(substitute_chromosome(token, chrom)))
        except FileNotFoundError:
            continue
    return _dedupe_preserving_order(resolved)


def _dedupe_preserving_order(values: Sequence[str]) -> list[str]:
    """Return ``values`` with duplicates removed while keeping first-seen order."""
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def filter_paths_for_chromosome(paths: Sequence[str], chrom: str) -> list[str]:
    """Return the subset of ``paths`` whose basename cleanly encodes ``chrom``.

    Filenames are matched heuristically against patterns such as ``.1.``,
    ``_1_``, ``-1-``, ``chr1``, ``.X.``, or ``chrX`` after stripping common
    compression and annotation suffix layers from the basename.
    """
    filtered = [path for path in paths if _path_matches_chromosome(path, chrom)]
    return _dedupe_preserving_order(filtered)


def ensure_output_directory(
    path: str | PathLike[str],
    *,
    label: str = "output directory",
) -> Path:
    """Create ``path`` when missing and warn once for the caller."""
    output_dir = Path(normalize_path_token(path))
    if output_dir.exists():
        if not output_dir.is_dir():
            raise NotADirectoryError(f"{label} exists but is not a directory: {output_dir}")
        return output_dir
    warnings.warn(
        f"{label} does not exist and is created: {output_dir}",
        UserWarning,
        stacklevel=2,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def ensure_output_paths_available(
    paths: Iterable[str | PathLike[str]],
    *,
    overwrite: bool = False,
    label: str = "output artifact",
) -> list[Path]:
    """Validate that fixed output files can be written without silent overwrite.

    Parameters
    ----------
    paths : iterable of str or os.PathLike[str]
        Deterministic output file paths that the active workflow may write.
        Parent directories are not created by this helper; call
        :func:`ensure_output_directory` or :func:`ensure_output_parent_directory`
        where appropriate.
    overwrite : bool, optional
        If ``False``, any existing path raises before the caller writes
        anything. If ``True``, existing paths are allowed but are not deleted by
        this helper. Default is ``False``.
    label : str, optional
        Human-readable artifact label used in error messages. Default is
        ``"output artifact"``.

    Returns
    -------
    list of pathlib.Path
        Normalized output paths in the same order as ``paths``.

    Raises
    ------
    FileExistsError
        If one or more paths already exist and ``overwrite`` is ``False``.
    """
    normalized_paths = [Path(normalize_path_token(path)) for path in paths]
    if overwrite:
        return normalized_paths
    existing = [path for path in normalized_paths if path.exists()]
    if existing:
        raise FileExistsError(_format_output_collision_message(existing, label=label))
    return normalized_paths


def preflight_output_artifact_family(
    produced_paths: Iterable[str | PathLike[str]],
    owned_paths: Iterable[str | PathLike[str]],
    *,
    overwrite: bool = False,
    label: str = "output artifact",
) -> list[Path]:
    """Preflight an owned output artifact family and return stale siblings.

    Parameters
    ----------
    produced_paths : iterable of str or os.PathLike
        Paths the **current invocation** will actually write. This subset
        depends on the run's flags (e.g., output format, conditional sidecars).
        Paths here are exempt from the overwrite-protection check because they
        are about to be replaced by the current run.
    owned_paths : iterable of str or os.PathLike
        Every path the workflow could write under any configuration — its
        stable "territory" inside the output directory. This is the union
        across every supported flag combination, including conditional
        sidecars that fire only on some runs. ``produced_paths`` is always a
        subset of ``owned_paths``.
    overwrite : bool, optional
        If ``False``, refuse to start the run when any owned path already
        exists on disk. If ``True``, allow the run to overwrite produced paths
        and return stale-sibling paths for post-success cleanup. Default is
        ``False``.

    Returns
    -------
    list of pathlib.Path
        With ``overwrite=False``: always ``[]`` (the function either passed
        the existence check or raised). With ``overwrite=True``: paths that
        are owned by the workflow, exist on disk, and are **not** in
        ``produced_paths`` — i.e., leftovers from a previous run with
        different flags that should be removed by
        :func:`remove_output_artifacts` after the current run succeeds.

    Raises
    ------
    FileExistsError
        If ``overwrite=False`` and at least one owned path already exists on
        disk.

    Notes
    -----
    Why two lists instead of one — the ``owned_paths`` / ``produced_paths``
    split is the core mechanism for safe stale cleanup across runs with
    different flag combinations. A single list cannot express both "the
    workflow controls this path" and "the current run will write this path,"
    and conflating them produces silent bugs:

    - Listing only what the run *produces* loses stale-cleanup ability: a
      file written by an earlier run with a different flag combination would
      survive forever, because the current run does not know to claim it.
    - Listing the entire family as *produced* makes the preflight refuse to
      overwrite siblings the current run is not actually touching, and
      removes stale-cleanup ability entirely (every owned path is "produced"
      and therefore exempt from the cleanup return list at the bottom of
      this function).
    - Conditionally written artifacts (sidecars or alternate formats that fire
      only on some runs)
      must always be in ``owned_paths`` so a stale copy from a prior run is
      detected and either blocks the run (no-overwrite) or gets cleaned up
      (overwrite). They must be in ``produced_paths`` only when the current
      run will actually write them, so a clean re-run that no longer needs
      the sidecar correctly identifies the old file as stale and removes it.
      Always-written audit sidecars belong in both lists on every invocation,
      including clean runs that write only a header row.

    Concrete example. Run A writes ``output_format='both'`` so both
    ``sumstats.parquet`` and ``sumstats.sumstats.gz`` exist. Run B reuses
    the directory with ``output_format='parquet' --overwrite`` and produces
    only ``sumstats.parquet``. The leftover ``sumstats.sumstats.gz`` is in
    ``owned_paths`` (workflow controls it) but not in ``produced_paths``
    (this run is not writing it), so it is returned here and unlinked by
    :func:`remove_output_artifacts`.

    For always-written audit sidecars such as the current
    ``dropped_snps/dropped.tsv.gz`` sumstats liftover sidecar, include the path
    in both ``owned_paths`` and ``produced_paths`` unconditionally. A clean run
    still produces a header-only file, so the prior sidecar is replaced in
    place rather than treated as a stale conditional artifact.
    """
    produced = _dedupe_paths(Path(normalize_path_token(path)) for path in produced_paths)
    owned = _dedupe_paths(Path(normalize_path_token(path)) for path in owned_paths)
    if not overwrite:
        # No-overwrite contract: any owned path on disk is a hard collision,
        # whether or not the current run intended to produce it. This is what
        # makes a conditional sidecar from a prior run still block a fresh run
        # without --overwrite.
        existing = [path for path in owned if path.exists()]
        if existing:
            raise FileExistsError(_format_output_collision_message(existing, label=label))
        return []
    # Overwrite contract: paths the current run will produce get replaced in
    # place by the writers; everything else the workflow owns that exists on
    # disk is stale (left over from a previous run with different flags) and
    # is returned here for the caller to feed into remove_output_artifacts().
    # Membership test uses produced_set so a path declared as "produced" is
    # never reported as stale, even if both lists contain it.
    produced_set = set(produced)
    return [path for path in owned if path.exists() and path not in produced_set]


def remove_output_artifacts(paths: Iterable[str | PathLike[str]]) -> None:
    """Remove stale workflow-owned output artifacts after a successful write."""
    for path in _dedupe_paths(Path(normalize_path_token(path)) for path in paths):
        if not path.exists():
            continue
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)
        else:
            path.unlink()


def ensure_output_parent_directory(
    path: str | PathLike[str],
    *,
    label: str = "output path",
) -> Path:
    """Create the parent directory for ``path`` when missing and warn."""
    output_path = Path(normalize_path_token(path))
    parent = output_path.parent
    if str(parent) in {"", "."}:
        return output_path
    ensure_output_directory(parent, label=f"{label} parent directory")
    return output_path


def _format_output_collision_message(paths: Sequence[Path], *, label: str) -> str:
    """Build a compact, shared message for existing output file collisions."""
    preview = ", ".join(str(path) for path in paths[:3])
    if len(paths) > 3:
        preview = f"{preview}, ... ({len(paths) - 3} more)"
    noun = label if len(paths) == 1 else f"{label}s"
    return (
        f"Refusing to overwrite existing {noun}: {preview}. "
        "Pass --overwrite on the CLI or overwrite=True in Python to replace them."
    )


def _dedupe_paths(paths: Iterable[Path]) -> list[Path]:
    """Return paths without duplicates while preserving first-seen order."""
    seen: set[Path] = set()
    ordered: list[Path] = []
    for path in paths:
        if path in seen:
            continue
        seen.add(path)
        ordered.append(path)
    return ordered


def _path_matches_chromosome(path: str, chrom: str) -> bool:
    """Heuristically match chromosome-coded filenames such as ``chr1`` or ``.1.``."""
    name = Path(path).stem
    if name.endswith(".annot"):
        name = Path(name).stem
    patterns = (
        rf"(^|[._-]){re.escape(str(chrom))}([._-]|$)",
        rf"(^|[._-])chr{re.escape(str(chrom))}([._-]|$)",
    )
    return any(re.search(pattern, name, flags=re.IGNORECASE) for pattern in patterns)
