"""Munged summary-statistics artifacts: the in-memory table, the loader and the footer codec.

Core functionality:
    Define ``SumstatsTable``, read one munged Parquet or legacy
    ``.sumstats(.gz)`` artifact back into it with :func:`load_sumstats`, and
    encode and decode the ``ldsc:*`` identity footer of the Parquet file.

Overview
--------
This module is the contract between munging and its readers. Regression and
every other reader of a munged artifact import from here instead of from the
munging workflow in ``ldsc._sumstats.munger``. The footer carries only the
thin downstream-identity payload (artifact type, SNP identifier mode, genome
build and trait name, plus the run's QC thresholds); the loader rebuilds the
munge-time ``GlobalConfig`` snapshot from it and rejects footerless Parquet,
while legacy text inputs load with ``config_snapshot=None``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import logging
from os import PathLike
from pathlib import Path
from typing import Any, NamedTuple

import pandas as pd

from .._core.columns import (
    INTERNAL_SUMSTATS_ARTIFACT_SPEC_MAP,
    resolve_optional_column,
    resolve_required_column,
)
from .._core.snp_identity import (
    identity_artifact_metadata,
    identity_cleanup_keys,
    identity_mode_family,
    is_allele_aware_mode,
    normalize_snp_identifier_mode,
    validate_identity_artifact_metadata,
)
from .._runtime.paths import resolve_scalar_path
from ..config import GlobalConfig, normalize_trait_name
from ..errors import LDSCConfigError, LDSCInputError, LDSCUserError


LOGGER = logging.getLogger("LDSC.sumstats_table")

_REQUIRED_COLUMN_TROUBLESHOOTING = (
    "docs/user-guide/troubleshooting.md#munge-sumstats-could-not-map-a-required-column"
)
_CURATED_ARTIFACT_TROUBLESHOOTING = (
    "docs/user-guide/troubleshooting.md#munge-sumstats-munged-summary-statistics-are-malformed-or-outdated"
)


@dataclass(frozen=True)
class SumstatsTable:
    """In-memory LDSC-ready summary-statistics table.

    Parameters
    ----------
    data : pandas.DataFrame
        Munged table containing at least ``SNP``, ``Z``, and ``N``. Current
        package-written artifacts also include canonical ``CHR`` and ``POS``
        columns, filled with missing values when the raw input had no
        coordinates.
    has_alleles : bool
        Whether allele columns are expected to be present and validated.
    source_path : str or None
        Original data source path, if known.
    trait_name : str or None
        Assigned trait label for regression summaries, or ``None`` when no
        label was supplied or stored. Regression outputs then label the trait
        by the ``source_path`` filename.
    provenance : dict, optional
        Lightweight run metadata retained for debugging and output summaries.
        Tables from :meth:`SumstatsMunger.run` record ``coordinate_provenance``
        (source build, coordinate basis and liftover counts) and
        ``qc_thresholds``, the effective ``info_min`` and ``maf_min`` (``None``
        for a filter not applied), which :meth:`SumstatsMunger.write_output`
        writes to the Parquet footer. Tables loaded from Parquet record the
        decoded footer as ``metadata``, and legacy text tables set
        ``legacy_ldsc``. Default is an empty dict.
    config_snapshot : GlobalConfig or None, optional
        Shared configuration captured when the table was produced in-process or
        recovered from the ``sumstats.parquet`` footer. It is ``None`` for a
        legacy ``.sumstats`` or ``.sumstats.gz`` input, whose rsIDs are
        projected onto the regression SNPs in the LD-score directory downstream. Footerless
        Parquet is rejected. Coordinate-family munged artifacts store the final output
        genome build here; rsID-family artifacts store ``genome_build=None``
        because their merge identity is independent of coordinate build.
    """
    data: pd.DataFrame
    has_alleles: bool
    source_path: str | None
    trait_name: str | None
    provenance: dict[str, Any] = field(default_factory=dict)
    config_snapshot: GlobalConfig | None = None

    def validate(self) -> None:
        """Validate the minimum LDSC-ready table contract."""
        required = {"SNP", "Z", "N"}
        missing = required - set(self.data.columns)
        if missing:
            raise LDSCInputError(
                f"munge-sumstats could not map required column(s) {sorted(missing)} "
                f"from the input header for '{self.source_path or self.trait_name or 'sumstats'}'. "
                "Most likely the file uses unrecognized column names. Pass explicit hints "
                f"(e.g. --snp, --a1) or rename the columns. Other causes & fixes: "
                f"{_REQUIRED_COLUMN_TROUBLESHOOTING}"
            )
        if self.has_alleles and not {"A1", "A2"}.issubset(self.data.columns):
            raise LDSCInputError(
                f"munge-sumstats expected A1/A2 allele columns in "
                f"'{self.source_path or self.trait_name or 'sumstats'}' because this table is allele-aware. "
                "Most likely the artifact was produced without allele columns. Re-run munge-sumstats "
                "with --a1 <effect-allele column> --a2 <other column> or use a base --snp-identifier mode "
                "(rsid or chr_pos)."
            )

    def summary(self) -> dict[str, Any]:
        """Summarize the retained rows and provenance fields."""
        return {
            "n_rows": len(self.data),
            "has_alleles": self.has_alleles,
            "trait_name": self.trait_name,
            "source_path": self.source_path,
        }


def load_sumstats(path: str | PathLike[str], trait_name: str | None = None) -> SumstatsTable:
    """Load one munged summary-statistics artifact into a ``SumstatsTable``.

    Parameters
    ----------
    path : str or os.PathLike[str]
        Path token for the munged summary-statistics artifact. This may be a
        literal path or an exact-one glob pattern. Resolution happens at the
        workflow layer before suffix inference. ``.parquet`` files are read
        with :func:`pandas.read_parquet`; ``.sumstats.gz``, ``sumstats.gz``, and ``.sumstats``
        files are read as whitespace-delimited text. Other suffixes raise a
        clear ``ValueError``.
    trait_name : str or None, optional
        Optional trait label propagated into downstream regression summaries.
        When supplied, this value replaces any Parquet footer ``trait_name``
        for the returned table; the artifact on disk is not modified. When
        omitted, the loader uses the parquet footer ``trait_name`` when
        present; otherwise the returned ``trait_name`` is ``None`` and
        regression outputs label the trait by its filename. Surrounding
        whitespace is trimmed; a blank label is rejected. Default is ``None``.

    Returns
    -------
    SumstatsTable
        Validated in-memory table with canonical LDSC3 columns such as ``SNP``,
        ``CHR``, ``POS``, ``N``, and ``Z`` when present in the artifact. When the
        Parquet footer carries identity metadata, the returned table
        also recovers its munge-time ``GlobalConfig`` snapshot. Legacy text
        inputs have ``config_snapshot=None`` and are marked for panel projection.

    Raises
    ------
    ValueError
        If ``path`` does not resolve to exactly one file or if the munged
        artifact is missing required LDSC3 columns.

    Notes
    -----
    Format inference is suffix-based after exact-one path resolution:
    ``.parquet`` uses :func:`pandas.read_parquet`, ``.sumstats.gz`` (or the
    no-label filename ``sumstats.gz``) uses a
    gzip-compressed whitespace reader, and ``.sumstats`` uses a plain
    whitespace reader. Other suffixes are refused so callers do not
    accidentally parse CSV or raw GWAS input as munged summary statistics.

    Identity metadata is read from the parquet footer (discrete ``ldsc:*``
    keys). No ``metadata.json`` sidecar is consulted. Legacy ``.sumstats`` and
    ``.sumstats.gz`` inputs load with ``config_snapshot=None`` and explicit
    legacy LDSC provenance. A Parquet artifact without the required footer is rejected.
    """
    resolved = resolve_scalar_path(path, label="munged summary statistics")
    df = _read_curated_sumstats_artifact(resolved)
    try:
        resolved_columns = _resolve_curated_sumstats_columns(list(df.columns), context=resolved)
    except ValueError as exc:
        raise LDSCInputError(
            f"Cannot load munged summary statistics at '{resolved}': could not map the required SNP, N, and Z "
            "columns from its header. Most likely this is not a current munged summary-statistics file. "
            "Re-run munge-sumstats from the raw GWAS file. Other causes & fixes: "
            f"{_REQUIRED_COLUMN_TROUBLESHOOTING}"
        ) from exc
    df = df.loc[:, list(resolved_columns.values())].rename(
        columns={actual: canonical for canonical, actual in resolved_columns.items()}, copy=False
    )
    if not (isinstance(df.index, pd.RangeIndex) and df.index.start == 0 and df.index.step == 1):
        df = df.reset_index(drop=True)
    footer = _decode_sumstats_footer(resolved)
    if footer is None:
        if str(resolved).lower().endswith(".parquet"):
            raise LDSCInputError(
                _curated_sumstats_artifact_message(
                    resolved,
                    "the parquet file is missing required LDSC3 identity footer metadata",
                    "Footerless parquet is not a legacy compatibility artifact. Re-run `ldsc munge-sumstats` "
                    "from the raw GWAS input to create a self-describing sumstats.parquet file.",
                )
            )
        LOGGER.info(
            f"Recognized legacy munged summary statistics '{resolved}'. SNP values will be used as rsID "
            "lookup keys and projected onto the regression SNPs in the LD-score directory during regression."
        )
        table = SumstatsTable(
            data=df,
            has_alleles={"A1", "A2"}.issubset(df.columns),
            source_path=resolved,
            trait_name=_resolve_sumstats_trait_name(trait_name, None, resolved),
            provenance={"source_format": "legacy_sumstats", "legacy_ldsc": True},
            config_snapshot=None,
        )
        table.validate()
        return table
    config_snapshot = footer.config
    mode = normalize_snp_identifier_mode(config_snapshot.snp_identifier)
    if is_allele_aware_mode(mode) and not {"A1", "A2"}.issubset(df.columns):
        raise LDSCInputError(
            _curated_sumstats_artifact_message(
                resolved,
                f"snp_identifier='{mode}' requires A1/A2 columns, but the artifact does not contain both.",
                "Most likely the artifact was produced by an older munging run or one in a base mode.",
            )
        )
    # One keying pass gives both the identity drops and the kept rows' keys.
    effective_keys, dropped = identity_cleanup_keys(
        df,
        mode,
        context="loaded munged summary-statistics artifact",
        stage="load_sumstats_validation",
    )
    if not dropped.empty:
        reasons = ", ".join(
            f"{reason}={int(count)}"
            for reason, count in dropped["reason"].value_counts(sort=False).items()
        )
        raise LDSCInputError(
            _curated_sumstats_artifact_message(
                resolved,
                f"duplicate or invalid SNP identity rows were found ({reasons}).",
                "Most likely the artifact was hand-edited or produced before identity cleanup was enforced.",
            )
        )
    if bool(pd.isna(effective_keys).any()):
        raise LDSCInputError(
            _curated_sumstats_artifact_message(
                resolved,
                "missing or invalid SNP identity rows were found after identity cleanup.",
                "Most likely the artifact contains incomplete CHR/POS or allele identity columns.",
            )
        )
    table = SumstatsTable(
        data=df,
        has_alleles={"A1", "A2"}.issubset(df.columns),
        source_path=resolved,
        trait_name=_resolve_sumstats_trait_name(trait_name, footer.trait_name, resolved),
        provenance={"metadata": footer.fields},
        config_snapshot=config_snapshot,
    )
    table.validate()
    return table


def _read_curated_sumstats_artifact(path: str) -> pd.DataFrame:
    """Read one munged summary-statistics artifact according to the public suffix policy."""
    token = str(path)
    if token.endswith(".parquet"):
        return pd.read_parquet(path)
    if token.endswith(".sumstats.gz") or Path(token).name == "sumstats.gz":
        return pd.read_csv(path, sep=r"\s+", compression="gzip")
    if token.endswith(".sumstats"):
        return pd.read_csv(path, sep=r"\s+", compression="infer")
    raise LDSCInputError(
        f"Cannot load munged summary statistics at '{path}': unsupported file suffix. "
        "Most likely this is a raw GWAS file or CSV, not a munged summary-statistics file. "
        "Use a path ending in '.parquet', '.sumstats.gz', or '.sumstats', the default 'sumstats.gz', "
        "or run munge-sumstats first."
    )


def _resolve_curated_sumstats_columns(columns: list[str], *, context: str) -> dict[str, str]:
    """Resolve canonical internal summary-statistics artifact columns from ``columns``."""
    resolved = {
        canonical: resolve_required_column(columns, INTERNAL_SUMSTATS_ARTIFACT_SPEC_MAP[canonical], context=context)
        for canonical in ("SNP", "N", "Z")
    }
    for canonical in ("CHR", "POS", "A1", "A2", "FRQ"):
        actual = resolve_optional_column(columns, INTERNAL_SUMSTATS_ARTIFACT_SPEC_MAP[canonical], context=context)
        if actual is not None:
            resolved[canonical] = actual
    return {canonical: resolved[canonical] for canonical in ("SNP", "CHR", "POS", "N", "Z", "A1", "A2", "FRQ") if canonical in resolved}


def encode_sumstats_footer(
    config_snapshot: GlobalConfig, trait_name: str | None, qc_thresholds: dict[str, float | None] | None = None,
) -> dict[bytes, bytes]:
    """Return discrete ``ldsc:*`` Parquet footer keys for a munged summary-statistics artifact.

    Mirrors the reference-panel footer convention in ``_refpanel/adapters.py``. A
    ``None`` ``genome_build`` (rsID-family) or ``trait_name`` is encoded as an
    empty string and decoded back to ``None`` by :func:`_decode_sumstats_footer`. ``qc_thresholds`` (the run's
    effective ``info_min`` and ``maf_min``) adds ``ldsc:qc_info_min`` and
    ``ldsc:qc_maf_min``; an empty value means that filter was not applied. Both
    keys are omitted when the thresholds are unknown, and readers ignore them.
    """
    identity = identity_artifact_metadata(
        artifact_type="sumstats",
        snp_identifier=config_snapshot.snp_identifier,
        genome_build=config_snapshot.genome_build,
    )
    footer = {
        f"ldsc:{key}".encode("utf-8"): ("" if value is None else str(value)).encode("utf-8")
        for key, value in identity.items()
    }
    footer[b"ldsc:trait_name"] = ("" if trait_name is None else str(trait_name)).encode("utf-8")
    if qc_thresholds is not None:
        for key, value in qc_thresholds.items():
            footer[f"ldsc:qc_{key}".encode("utf-8")] = ("" if value is None else repr(float(value))).encode("utf-8")
    return footer


class _SumstatsFooter(NamedTuple):
    """Decoded identity footer of a munged summary-statistics Parquet file."""

    config: GlobalConfig
    trait_name: str | None
    fields: dict[str, Any]


_SUMSTATS_FOOTER_KEYS = ("artifact_type", "snp_identifier", "genome_build", "trait_name")


def _decode_sumstats_footer(path: str) -> _SumstatsFooter | None:
    """Decode the identity footer of a munged summary-statistics Parquet file.

    Returns ``None`` for text summary statistics and for Parquet without an
    ``ldsc:artifact_type`` key; the public loader distinguishes them and rejects
    footerless Parquet. A footer is accepted only in a form
    :func:`encode_sumstats_footer` writes: all four keys present, artifact type
    ``sumstats``, a valid SNP identifier mode, an empty genome build in
    rsID-family modes and hg19 or hg38 in chr_pos-family modes, and a trait name
    that is empty or not blank. The ``ldsc:qc_*`` keys are not read.
    """
    if not str(path).endswith(".parquet"):
        return None
    import pyarrow.parquet as pq

    raw = pq.read_schema(path).metadata or {}
    if b"ldsc:artifact_type" not in raw:
        return None
    missing = [f"ldsc:{key}" for key in _SUMSTATS_FOOTER_KEYS if f"ldsc:{key}".encode("utf-8") not in raw]
    if missing:
        raise LDSCInputError(
            _curated_sumstats_artifact_message(
                path,
                f"its parquet footer lacks the {' and '.join(missing)} key{'s' if len(missing) > 1 else ''}.",
                "Most likely it was written by an older LDSC3 version or edited after munging.",
            )
        )
    fields = {key: raw[f"ldsc:{key}".encode("utf-8")].decode("utf-8") or None for key in _SUMSTATS_FOOTER_KEYS}
    try:
        mode = validate_identity_artifact_metadata(fields, expected_artifact_type="sumstats")
    except LDSCUserError as exc:
        raise LDSCInputError(
            _curated_sumstats_artifact_message(
                path,
                "its embedded parquet footer has invalid identity provenance.",
                "Most likely it predates the current schema or was edited after munging.",
            )
        ) from exc
    build = fields["genome_build"]
    if identity_mode_family(mode) == "rsid":
        valid_build, expected = build is None, "an empty genome_build"
    else:
        valid_build, expected = build in {"hg19", "hg38"}, "genome_build 'hg19' or 'hg38'"
    if not valid_build:
        raise LDSCInputError(
            _curated_sumstats_artifact_message(
                path,
                f"its parquet footer records snp_identifier={mode!r} with genome_build={build!r}, "
                f"but munge-sumstats writes {expected} for that mode.",
                "Most likely the footer was edited after munging or written by another tool.",
            )
        )
    try:
        trait_name = normalize_trait_name(fields["trait_name"])
    except LDSCConfigError as exc:
        raise LDSCInputError(
            _curated_sumstats_artifact_message(
                path,
                "its embedded parquet footer has a blank trait_name field.",
                "Most likely the artifact was hand-edited or written by an outdated workflow.",
            )
        ) from exc
    config = GlobalConfig(snp_identifier=mode, genome_build=build, log_level="INFO")
    return _SumstatsFooter(config, trait_name, {**fields, "snp_identifier": mode, "trait_name": trait_name})


def _curated_sumstats_artifact_message(path: str | PathLike[str], detail: str, likely: str) -> str:
    """Return the standard malformed munged-artifact diagnostic."""
    return (
        f"Cannot load munged summary statistics at '{path}': {detail} "
        f"{likely} Re-run munge-sumstats from the raw GWAS file. Other causes & fixes: "
        f"{_CURATED_ARTIFACT_TROUBLESHOOTING}"
    )


def _resolve_sumstats_trait_name(
    explicit_trait_name: str | None,
    footer_trait_name: str | None,
    resolved_path: str,
) -> str | None:
    """Resolve an assigned trait label from the CLI/API argument, then the decoded parquet footer.

    Returns ``None`` when neither supplies a label. Regression outputs apply the
    filename fallback themselves, so ``SumstatsMunger.write_output`` never writes
    a derived label back into a footer.
    """
    try:
        explicit = normalize_trait_name(explicit_trait_name)
    except LDSCConfigError as exc:
        raise LDSCConfigError(
            f"Cannot load munged summary statistics at '{resolved_path}': trait_name argument is blank. "
            "Most likely whitespace was passed as the trait label. "
            "Use a non-empty trait name or omit trait_name to derive it from the artifact."
        ) from exc
    return explicit if explicit is not None else footer_trait_name
