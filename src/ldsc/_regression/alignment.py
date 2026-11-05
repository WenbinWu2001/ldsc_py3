"""Align summary statistics with LD scores into regression datasets.

Each trait's summary statistics are prepared once per run
(``_prepare_trait_sumstats``): a legacy ``.sumstats(.gz)`` table is projected
by rsID onto the regression SNPs of the LD-score directory, rows whose ``Z`` or
``N`` is missing or infinite are dropped, and keyed rows are cached per SNP
identifier mode. The baseline LD scores are checked and keyed the same way
(``_PreparedLDScores``). This module resolves the effective SNP identifier mode
of each merge, checks the genome build of coordinate-family inputs, merges a
trait with the LD scores (``_prepare_h2_inputs``) or a trait pair with them,
harmonizing trait 2's alleles to trait 1 (``_rg_dataset``), and selects one
model's LD-score columns and SNP-count vectors (``_select_model_columns``).
Fitting and reporting belong to ``fitting``.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import logging
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

from .._core.columns import infer_chr_pos_columns
from .._core.coordinates import CHR_POS_KEY_COLUMN, build_chr_pos_key_frame
from .._core.row_alignment import assert_same_snp_rows, first_identity_mismatch
from .._core.snp_identity import (
    BASE_COMPLEMENT,
    FLIP_ALLELES,
    MATCH_ALLELES,
    clean_identity_artifact_table,
    effective_merge_key_series,
    identity_base_mode,
    identity_mode_family,
    is_allele_aware_mode,
    normalize_snp_identifier_mode,
    resolve_regression_identity_mode,
)
from .._genome.build_policy import GenomeBuildCheck
from .._ldscore.results import LDScoreResult
from .._ldscore.outputs import REGRESSION_LD_SCORE_COLUMN
from .._ldscore.source import LDScoreSource
from .._sumstats.table import SumstatsTable
from ..config import GlobalConfig, RegressionConfig
from ..errors import LDSCInputError, LDSCUsageError
from .outputs import LEGACY_SUMSTATS_DROP_COLUMNS
from .results import RegressionDataset


COMMON_COUNT_KEY = "common_reference_snp_counts"
ALL_COUNT_KEY = "all_reference_snp_counts"
_NONFINITE_LD_DOC = "docs/user-guide/troubleshooting.md#regression-ld-scores-or-regression-weights-are-missing-or-infinite"
REGRESSION_IDENTITY_KEY_COLUMN = "_ldsc_regression_identity_key"
_REGRESSION_NO_OVERLAP_DOC = (
    "docs/user-guide/troubleshooting.md#regression-no-snps-remain-after-merging-regression-inputs"
)
_RG_ALLELES_DOC = "docs/user-guide/troubleshooting.md#rg-a-trait-has-no-a1a2-columns"
LOGGER = logging.getLogger("LDSC.regression_runner")


@dataclass(frozen=True)
class _PreparedRegressionInputs:
    """Shared trait alignment before model-dependent column/count filtering.

    ``dataset`` carries baseline LD values and aligned traits in genomic order.
    ``source_rows`` maps these rows back to the original LD-score table, so
    query batches reuse exactly the same join and allele orientation.
    """

    dataset: RegressionDataset
    source_rows: np.ndarray
    ldscore: LDScoreResult | LDScoreSource
    ldscore_identity_mode: str
    counts_by_column: dict[str, dict]
    ld_variances: dict = field(default_factory=dict)


@dataclass(frozen=True)
class _PartitionedModelInputs:
    """Aligned model state without input readers or the full query source."""

    dataset: RegressionDataset
    counts_by_column: dict[str, dict]
    ld_variances: dict = field(default_factory=dict)


@dataclass(frozen=True)
class _PreparedTraitSumstats:
    """Prepared trait summary statistics: one trait, prepared once for every model that uses it.

    For a legacy ``.sumstats(.gz)`` input, ``table`` is its projection onto the
    LD-score panel and ``legacy_drops`` records the rows that projection dropped;
    otherwise ``table`` is the input and ``legacy_drops`` is ``None``. Either way,
    rows whose ``Z`` or ``N`` is missing or infinite are removed and counted in
    ``dropped_nonfinite_statistic_rows``; for a legacy input they are also
    ``legacy_drops`` rows with reason ``nonfinite_statistic``. Rows keyed
    for a merge are built on first use and cached per effective SNP identifier
    mode and downgrade cleanup, so a trait shared by many rg pairs is keyed at
    most twice: once for merges in the LD-score directory's mode and once for
    downgraded merges, which run in its base mode after cleanup.
    """

    table: SumstatsTable
    legacy_drops: pd.DataFrame | None = None
    dropped_nonfinite_statistic_rows: int = 0
    _keyed: dict[tuple[str, bool], tuple[pd.DataFrame, int]] = field(default_factory=dict, repr=False, compare=False)

    @property
    def context(self) -> str:
        """Name the trait in messages."""
        return self.table.source_path or self.table.trait_name or "sumstats"

    def keyed_rows(self, mode: str, *, cleaned: bool) -> tuple[pd.DataFrame, int]:
        """Return rows keyed for a ``mode`` merge and the rows a downgrade cleanup dropped."""
        if (mode, cleaned) not in self._keyed:
            self._keyed[mode, cleaned] = _key_regression_rows(self.table.data, mode, cleaned=cleaned, context=self.context)
        return self._keyed[mode, cleaned]


@dataclass(frozen=True)
class _PreparedLDScores:
    """Baseline LD scores keyed for merges, cached like ``_PreparedTraitSumstats``.

    Keyed frames carry ``_ldsc_source_row``, each row's position in the
    baseline table, so query batches reuse exactly the same join.
    """

    source: LDScoreSource | LDScoreResult
    _keyed: dict[tuple[str, bool], tuple[pd.DataFrame, int]] = field(default_factory=dict, repr=False, compare=False)

    def __post_init__(self):
        """Stop before any fit if a baseline LD score or regression weight is missing or infinite.

        The whole table is checked once per run: every row is a regression
        SNP, so a bad value means the directory itself is damaged, even where
        no trait matches that SNP.
        """
        table = self.source.baseline_table
        problems = _nonfinite_columns(table, [*self.source.baseline_columns, REGRESSION_LD_SCORE_COLUMN], table["SNP"])
        if problems:
            metadata = self.source.output_paths.get("metadata")
            where = f"LD-score directory '{Path(metadata).parent}'" if metadata else "the supplied LD scores"
            raise LDSCInputError(
                f"Regression cannot use {where}: missing or infinite values in {'; '.join(problems)}. "
                "Every LD score and regression weight must be finite; dropping these SNPs would change the "
                "regression without notice. Most likely the directory was edited or only partly written after it "
                "was created. Regenerate it with `ldsc ldscore` or `ldsc convert-legacy-ldscores`, or restore the "
                f"original files; do not replace missing values with zero. Other causes & fixes: {_NONFINITE_LD_DOC}"
            )

    def keyed_rows(self, mode: str, *, cleaned: bool) -> tuple[pd.DataFrame, int]:
        """Return baseline rows keyed for a ``mode`` merge and the rows a downgrade cleanup dropped."""
        if (mode, cleaned) not in self._keyed:
            frame = _assemble_regression_ldscore_table(self.source)
            frame["_ldsc_source_row"] = np.arange(len(frame), dtype=np.int64)
            self._keyed[mode, cleaned] = _key_regression_rows(frame, mode, cleaned=cleaned, context="LD-score table")
        return self._keyed[mode, cleaned]


@dataclass(frozen=True)
class _PreparedH2Dataset:
    """One trait's h2 dataset and the provenance of its trait preparation.

    ``legacy_drops`` and ``dropped_nonfinite_statistic_rows`` come from the
    trait's ``_PreparedTraitSumstats``; ``h2`` results record both.
    """

    dataset: RegressionDataset
    legacy_drops: pd.DataFrame | None
    dropped_nonfinite_statistic_rows: int


def _prepare_h2_dataset(
    sumstats_table: SumstatsTable,
    ldscore_result: LDScoreSource | LDScoreResult,
    config: RegressionConfig,
    global_config: GlobalConfig,
    *,
    query_columns: Sequence[str] | None = None,
) -> _PreparedH2Dataset:
    """Prepare one trait and build its h2 dataset; ``ldsc h2`` and ``RegressionRunner.build_dataset`` use it.

    ``query_columns`` adds query LD scores to the baseline model.
    ``global_config`` supplies the SNP identifier mode of an input without
    recorded provenance. The trait's keyed rows are released on return; the
    dataset keeps only the aligned rows.
    """
    trait = _prepare_trait_sumstats(sumstats_table, ldscore_result)
    dataset = _h2_dataset(trait, _PreparedLDScores(ldscore_result), config, global_config, query_columns=query_columns)
    return _PreparedH2Dataset(dataset, trait.legacy_drops, trait.dropped_nonfinite_statistic_rows)


def _h2_dataset(
    trait: _PreparedTraitSumstats,
    ldscores: _PreparedLDScores,
    config: RegressionConfig,
    global_config: GlobalConfig,
    *,
    build_checked: bool = False,
    query_columns: Sequence[str] | None = None,
) -> RegressionDataset:
    """Build one model's h2 dataset from prepared trait summary statistics.

    ``rg`` fits each trait's own SNP heritability with it, reusing the trait
    it prepared once for the run; ``build_checked`` is as in
    ``_prepare_h2_inputs``.
    """
    prepared = _prepare_h2_inputs(trait, ldscores, config, global_config, build_checked=build_checked)
    return _dataset_from_prepared(prepared, config, query_columns=query_columns)


def _prepare_rg_dataset(
    sumstats_table_1: SumstatsTable,
    sumstats_table_2: SumstatsTable,
    ldscore_result: LDScoreSource | LDScoreResult,
    config: RegressionConfig,
    global_config: GlobalConfig,
) -> RegressionDataset:
    """Prepare both traits and build their pair dataset; see ``RegressionRunner.build_rg_dataset``.

    Both traits must carry ``A1`` and ``A2``. The rg pair run prepares each
    trait once and calls ``_rg_dataset`` for every pair instead.
    """
    _require_rg_alleles([sumstats_table_1, sumstats_table_2])
    return _rg_dataset(
        _prepare_trait_sumstats(sumstats_table_1, ldscore_result),
        _prepare_trait_sumstats(sumstats_table_2, ldscore_result),
        _PreparedLDScores(ldscore_result),
        config,
        global_config,
    )


def _prepare_h2_inputs(
    trait: _PreparedTraitSumstats,
    ldscores: _PreparedLDScores,
    config: RegressionConfig,
    global_config: GlobalConfig,
    *,
    build_checked: bool = False,
) -> _PreparedRegressionInputs:
    """Merge prepared trait summary statistics with the baseline LD scores.

    This alignment keeps every baseline column, so several models can share
    it; ``_dataset_from_prepared`` selects each model's columns and counts.
    ``build_checked`` skips the genome-build check when the caller already
    checked this trait and LD-score directory together. ``global_config``
    supplies the SNP identifier mode of an input without recorded provenance.
    """
    sumstats_table, ldscore_result = trait.table, ldscores.source
    weight_column = REGRESSION_LD_SCORE_COLUMN
    ref_ld_columns = list(ldscore_result.baseline_columns)
    identity = _resolve_h2_identity(sumstats_table, ldscore_result, global_config, config)
    identifier_mode = identity.effective_mode
    if not build_checked:
        _regression_genome_build(sumstats_table, ldscore_result, identifier_mode, context=trait.context)
    ldscore_mode = _ldscore_identity_mode(ldscore_result, global_config, fallback_mode=identifier_mode)
    sumstats_keyed, sumstats_dropped = trait.keyed_rows(identifier_mode, cleaned=identity.downgrade_applied)
    ldscore_keyed, ldscore_dropped = ldscores.keyed_rows(identifier_mode, cleaned=identity.downgrade_applied)
    if identity.downgrade_applied:
        sumstats_mode = _sumstats_identity_mode(sumstats_table, fallback_mode=ldscore_mode)
        LOGGER.warning(
            f"Identity downgrade enabled: LD-score mode {ldscore_mode}, summary-statistics mode {sumstats_mode}; "
            f"running regression with effective snp_identifier={identifier_mode!r}. "
            f"Dropped {sumstats_dropped + ldscore_dropped} duplicate effective-key rows before merge."
        )
    key = _merge_key_column(identifier_mode)
    payload = {name: ldscore_keyed[name] for name in [key, *ref_ld_columns, weight_column, "_ldsc_source_row"]}
    # LD-score is the left frame so the inner merge inherits its genomic
    # order; the block jackknife needs genomically contiguous rows for valid
    # SEs, and the sumstats order alone (legacy .sumstats.gz is unsorted) does
    # not guarantee that. Values still align by key regardless of side.
    # The payload borrows the keyed columns; the merge makes the one copy.
    # The fit uses chi^2 = Z^2, so Z keeps the summary statistics' A1 orientation.
    merged = pd.merge(pd.DataFrame(payload, copy=False), sumstats_keyed, how="inner", on=key, sort=False)
    if merged.empty:
        raise LDSCInputError(
            f"h2 regression retained no overlapping {identifier_mode} SNPs after merging summary statistics "
            f"'{trait.context}' with {len(ldscore_keyed)} LD-score rows. Most likely the summary statistics and "
            "LD-score directory use different SNP identifiers or genome builds. Regenerate both inputs "
            "with the same `--snp-identifier` and `--genome-build`, or use `--allow-identity-downgrade` "
            f"for same-family allele-aware/base mixes. Other causes & fixes: {_REGRESSION_NO_OVERLAP_DOC}"
        )

    source_rows = _add_panel_coordinates(merged, ldscore_result)
    dataset = RegressionDataset(
        merged=merged,
        ref_ld_columns=ref_ld_columns,
        weight_column=weight_column,
        reference_snp_count_totals={},
        count_key_used_for_regression="",
        retained_ld_columns=ref_ld_columns,
        dropped_zero_variance_ld_columns=[],
        trait_names=[_trait_label(sumstats_table)],
        config_snapshot=ldscore_result.config_snapshot,
        effective_snp_identifier=identifier_mode,
        identity_downgrade_applied=identity.downgrade_applied,
        ldscore_overlap=ldscore_result.overlap,
    )
    dataset.validate()
    # Baseline variances do not depend on the query, so every model reuses them.
    return _PreparedRegressionInputs(dataset, source_rows, ldscore_result, ldscore_mode,
                                    {str(record["column"]): record for record in ldscore_result.count_records},
                                    dict(_ld_column_variances(merged, ref_ld_columns)))


def _dataset_from_prepared(prepared, config, *, query_columns=None) -> RegressionDataset:
    """Select one complete genome-wide model and apply its count/variance rules.

    Query identities are validated when this fit selects query columns.
    Shared preparation never decides model-dependent filters or weights.
    The returned frame borrows aligned columns from ``prepared`` and the
    selected query batch. Treat these inputs as read-only during fitting.
    """
    selected = list(query_columns or [])
    values = _read_query_batch(prepared, selected) if selected else None
    return _dataset_from_query_values(prepared, selected, values, config)


def _read_query_batch(prepared, columns):
    """Validate the batch's stored rows once and select the shared trait rows.

    Query rows whose identity columns equal the baseline's exactly need no
    identifier-key comparison; others are compared by key as before.
    """
    frame = prepared.ldscore.read_queries(columns)
    baseline = prepared.ldscore.baseline_table
    identity = [name for name in ("CHR", "SNP", "POS", "A1", "A2") if name in frame.columns]
    if not set(identity) <= set(baseline.columns) or first_identity_mismatch(baseline, frame, identity) is not None:
        assert_same_snp_rows(
            baseline, frame,
            context="query rows must match baseline rows on CHR/SNP/POS",
            snp_identifier=prepared.ldscore_identity_mode,
        )
    rows = prepared.source_rows
    return pd.DataFrame({column: frame[column].to_numpy()[rows] for column in columns}, copy=False)


def _rg_dataset(
    trait_1: _PreparedTraitSumstats,
    trait_2: _PreparedTraitSumstats,
    ldscores: _PreparedLDScores,
    config: RegressionConfig,
    global_config: GlobalConfig,
    *,
    build_checked: bool = False,
) -> RegressionDataset:
    """Build one pair's dataset from prepared trait summary statistics; see ``build_rg_dataset``.

    ``build_checked`` skips the genome-build check when the caller already
    checked both traits and the LD-score directory together. ``global_config``
    is the SNP identifier fallback, as in ``_prepare_h2_inputs``.
    """
    sumstats_table_1, sumstats_table_2, ldscore_result = trait_1.table, trait_2.table, ldscores.source
    identity = _resolve_rg_identity(sumstats_table_1, sumstats_table_2, ldscore_result, global_config, config)
    identifier_mode = identity.effective_mode
    if not build_checked:
        _regression_genome_build([sumstats_table_1, sumstats_table_2], ldscore_result,
                                 identifier_mode, context="rg")

    weight_column = REGRESSION_LD_SCORE_COLUMN
    ref_ld_columns = list(ldscore_result.baseline_columns)
    left_keyed, left_dropped = trait_1.keyed_rows(identifier_mode, cleaned=identity.downgrade_applied)
    right_keyed, right_dropped = trait_2.keyed_rows(identifier_mode, cleaned=identity.downgrade_applied)
    ldscore_keyed, ldscore_dropped = ldscores.keyed_rows(identifier_mode, cleaned=identity.downgrade_applied)
    if identity.downgrade_applied:
        ldscore_mode = _ldscore_identity_mode(ldscore_result, global_config, fallback_mode=identifier_mode)
        trait_1_mode = _sumstats_identity_mode(sumstats_table_1, fallback_mode=ldscore_mode)
        trait_2_mode = _sumstats_identity_mode(sumstats_table_2, fallback_mode=ldscore_mode)
        LOGGER.warning(
            f"Identity downgrade enabled: LD-score mode {ldscore_mode}, trait 1 mode {trait_1_mode}, "
            f"trait 2 mode {trait_2_mode}; running rg with effective snp_identifier={identifier_mode!r}. "
            f"Dropped {left_dropped + right_dropped + ldscore_dropped} duplicate effective-key rows before merge."
        )

    key = _merge_key_column(identifier_mode)
    # The merge inputs borrow the keyed columns; each merge makes one copy.
    left = left_keyed.rename(columns={"N": "N1", "Z": "Z1", "FRQ": "FRQ1"}, copy=False)
    right = right_keyed.rename(columns={"A1": "A1x", "A2": "A2x", "N": "N2", "Z": "Z2", "FRQ": "FRQ2"}, copy=False)
    right_payload = [column for column in ["A1x", "A2x", "N2", "Z2", "FRQ2"] if column in right.columns]
    ld_payload = {name: ldscore_keyed[name] for name in [key, *ref_ld_columns, weight_column, "_ldsc_source_row"]}
    # LD-score leads the first merge so the result (and the final merge below)
    # inherits its genomic order for valid jackknife block SEs. See build_dataset.
    left_with_ld = pd.merge(pd.DataFrame(ld_payload, copy=False), left, how="inner", on=key, sort=False)
    merged = pd.merge(left_with_ld, pd.DataFrame({name: right[name] for name in [key, *right_payload]}, copy=False),
                      how="inner", on=key, sort=False)
    del left_with_ld
    _add_panel_coordinates(merged, ldscore_result)
    # Z and N are finite after trait preparation, and LD scores and weights
    # after _PreparedLDScores, so no row needs dropping here.
    if merged.empty:
        raise LDSCInputError(
            f"rg regression retained no overlapping {identifier_mode} SNPs after merging both summary-statistics "
            f"tables with {len(ldscore_keyed)} LD-score rows. Most likely one trait was munged with a "
            "different SNP identifier mode, genome build, or allele convention than the LD-score directory. "
            "Regenerate all inputs with matching identity settings, or use `--allow-identity-downgrade` "
            f"for same-family allele-aware/base mixes. Other causes & fixes: {_REGRESSION_NO_OVERLAP_DOC}"
        )

    compatible, flip = _allele_pair_flags([merged["A1"], merged["A2"], merged["A1x"], merged["A2x"]])
    orientation_unknown = ~compatible & _same_strand_ambiguous_pair(merged)
    drops = (int(orientation_unknown.sum()), int((~compatible & ~orientation_unknown).sum()))
    if drops[0] or drops[1]:
        LOGGER.warning(
            f"rg dropped {drops[0]} A/T or C/G SNPs whose strand cannot be aligned between "
            f"'{_trait_label(sumstats_table_1)}' and '{_trait_label(sumstats_table_2)}', and {drops[1]} SNPs whose "
            "two traits list different allele pairs."
        )
        merged = merged.take(np.flatnonzero(compatible))
        merged.index = pd.RangeIndex(len(merged))
    flip = flip[compatible]
    if merged.empty:
        error = LDSCInputError(
            "rg regression retained no allele-compatible SNPs after harmonizing the two summary-statistics tables "
            f"({drops[0]} A/T or C/G SNPs whose strand cannot be aligned, {drops[1]} SNPs whose two traits list "
            "different allele pairs). Most likely the two traits' A1/A2 columns (A1 is the effect allele, the "
            "allele the signed statistic refers to) hold different allele pairs for the shared SNPs, for example "
            "because one trait was munged from columns that are not its alleles. Re-run `ldsc munge-sumstats` for "
            "both traits with --a1 <effect-allele column> --a2 <other column> and the same SNP identifier mode."
        )
        error.allele_alignment_drops = drops
        raise error
    # The integer signs keep the legacy result dtype (float32 Z widens to float64).
    merged["Z2"] = merged["Z2"].to_numpy() * np.where(flip, -1, 1)
    if "FRQ2" in merged.columns:
        flip_index = pd.Series(flip, index=merged.index)
        valid_frq = pd.to_numeric(merged["FRQ2"], errors="coerce").between(0.0, 1.0, inclusive="both")
        transform = flip_index & valid_frq
        merged.loc[transform, "FRQ2"] = 1.0 - pd.to_numeric(
            merged.loc[transform, "FRQ2"], errors="coerce"
        )

    # Model columns are chosen on the final pair SNP set, after allele harmonization.
    dataset = RegressionDataset(
        ref_ld_columns=ref_ld_columns,
        weight_column=weight_column,
        **_select_model_columns(merged, ref_ld_columns, ldscore_result.count_records,
                                config.use_common_counts, command="rg"),
        trait_names=[_trait_label(sumstats_table_1), _trait_label(sumstats_table_2)],
        config_snapshot=ldscore_result.config_snapshot,
        effective_snp_identifier=identifier_mode,
        identity_downgrade_applied=identity.downgrade_applied,
        ldscore_overlap=ldscore_result.overlap,
        dropped_orientation_unknown_rows=drops[0],
        dropped_incompatible_allele_rows=drops[1],
    )
    dataset.validate()
    return dataset


def _dataset_from_query_values(prepared, selected, query_values, config):
    """Borrow aligned columns for one model without copying the baseline.

    Preparation and query loading already provide matching RangeIndexes.
    Borrowed columns must remain unmodified; numerical fitting allocates its
    own design array after the model's chi-square filter has selected rows.
    """
    ref_ld_columns = list(prepared.dataset.ref_ld_columns) + list(selected)
    merged = prepared.dataset.merged
    # A query's damage fails only its own model, inside the query's error handling.
    problems = _nonfinite_columns(query_values, selected, merged["SNP"]) if selected else []
    if problems:
        raise LDSCInputError(
            f"Regression cannot use query annotation {'; '.join(problems)}: its LD scores are missing or "
            "infinite on regression SNPs, and every LD score must be finite. Most likely the LD-score directory "
            "was edited or only partly written after it was created. Regenerate it with `ldsc ldscore`, or "
            f"restore the original query files. Other causes & fixes: {_NONFINITE_LD_DOC}"
        )
    if selected:
        weight_position = merged.columns.get_loc(prepared.dataset.weight_column)
        order = [*merged.columns[:weight_position], *selected, *merged.columns[weight_position:]]
        merged = pd.DataFrame({name: (query_values if name in selected else merged)[name].to_numpy(copy=False)
                               for name in order}, copy=False)
    dataset = replace(prepared.dataset, ref_ld_columns=ref_ld_columns, **_select_model_columns(
        merged, ref_ld_columns, prepared.counts_by_column, config.use_common_counts, command="h2",
        known_variances=prepared.ld_variances,
    ))
    dataset.validate()
    return dataset


def _regression_genome_build(
    sumstats_tables: SumstatsTable | Sequence[SumstatsTable],
    ldscore_result: LDScoreSource | LDScoreResult,
    identifier_mode: str,
    *,
    context: str,
) -> None:
    """Check all coordinate-family artifacts before fitting any model.

    Valid snapshots establish builds without coordinate re-inference.
    In-memory objects lacking a concrete snapshot use bounded coordinate
    samples. rsID joins need no coordinate build. Regression has no build
    override; failures direct users to matching upstream artifacts.
    """
    if identity_mode_family(identifier_mode) != "chr_pos":
        return

    check = GenomeBuildCheck('auto', context='regression', option=None, logger=LOGGER)
    tables = sumstats_tables if isinstance(sumstats_tables, (list, tuple)) else [sumstats_tables]
    role = '--sumstats-sources' if isinstance(sumstats_tables, (list, tuple)) else '--sumstats-file'
    for table in tables:
        source = table.source_path or table.trait_name or context
        build = getattr(table.config_snapshot, 'genome_build', None)
        if _has_legacy_sumstats_source(table):
            check.unchecked(role, [source], 'legacy rsID projection inherits panel coordinates')
        elif build in {'hg19', 'hg38'}:
            check.declared(role, source, build)
        else:
            check.frames(role, source, [table.data])
    panel_source = str(ldscore_result.output_paths.get('metadata', 'LD-score directory'))
    build = getattr(ldscore_result.config_snapshot, 'genome_build', None)
    if build in {'hg19', 'hg38'}:
        check.declared('--ldscore-dir', panel_source, build)
    else:
        check.frames('--ldscore-dir', panel_source, [ldscore_result.baseline_table])
    check.resolve()


def _snapshot_identity_mode(snapshot: GlobalConfig | None) -> str | None:
    """Return a normalized SNP identifier mode from ``snapshot`` if present."""
    value = getattr(snapshot, "snp_identifier", None)
    if value is None:
        return None
    return normalize_snp_identifier_mode(value)


def _has_legacy_sumstats_source(table: SumstatsTable) -> bool:
    """Return whether ``table`` originated from legacy text."""
    return bool(table.provenance.get("legacy_ldsc", False))


def _is_allele_less_legacy(table: SumstatsTable) -> bool:
    """Return whether ``table`` is legacy text without A1/A2, which projects by rsID under the panel's base mode."""
    return _has_legacy_sumstats_source(table) and not {"A1", "A2"}.intersection(table.data.columns)


def _prepare_trait_sumstats(
    table: SumstatsTable,
    ldscore_result: LDScoreSource | LDScoreResult,
) -> _PreparedTraitSumstats:
    """Prepare one trait's summary statistics once for every model that uses it.

    A legacy ``.sumstats(.gz)`` table is projected onto the LD-score panel under
    the panel's recorded SNP identifier mode; other tables are used as given.
    Rows whose ``Z`` or ``N`` is missing or infinite cannot enter a fit and are
    dropped with a warning instead of failing the run (a recorded exception to
    the fail-on-unrequested-filtering rule, as in legacy LDSC).
    """
    context = table.source_path or table.trait_name or "sumstats"
    legacy_drops = None
    if _has_legacy_sumstats_source(table):
        table, legacy_drops = _project_legacy_sumstats_to_panel(table, ldscore_result)
        dropped = int(legacy_drops["reason"].eq("nonfinite_statistic").sum())
    else:
        finite = _finite_statistics(table.data)
        dropped = int((~finite).sum())
        if dropped:
            kept = table.data.take(np.flatnonzero(finite.to_numpy()))
            kept.index = pd.RangeIndex(len(kept))
            table = replace(table, data=kept)
    if dropped:
        LOGGER.warning(
            "Dropped %d row(s) of summary statistics '%s' whose Z or N is missing or infinite; such rows cannot "
            "enter a regression. The count is recorded as dropped_nonfinite_statistic_rows in the result metadata.",
            dropped, context,
        )
    return _PreparedTraitSumstats(table, legacy_drops, dropped)


def _finite_statistics(frame: pd.DataFrame) -> pd.Series:
    """Mark rows whose ``Z`` and ``N`` are both finite numbers."""
    return (np.isfinite(pd.to_numeric(frame["Z"], errors="coerce"))
            & np.isfinite(pd.to_numeric(frame["N"], errors="coerce")))


def _project_legacy_sumstats_to_panel(
    table: SumstatsTable,
    ldscore_result: LDScoreSource | LDScoreResult,
) -> tuple[SumstatsTable, pd.DataFrame]:
    """Project one legacy munged summary-statistics table onto the regression SNPs in the LD-score directory.

    Legacy ``SNP`` values are treated only as rsID lookup keys. The supplied
    LD-score panel owns output coordinates and identifier mode. Allele-aware
    panels also own allele order, except for A/T and C/G SNPs, whose strand
    against the panel is unknown: they keep their source alleles and ``Z``.
    Allele-unaware panels retain the source alleles so two traits can still
    be harmonized for genetic correlation.

    Parameters
    ----------
    table : SumstatsTable
        Legacy text table marked with ``provenance['legacy_ldsc']``.
    ldscore_result : LDScoreSource or LDScoreResult
        LDSC3 LD-score directory used by the pending regression.

    Returns
    -------
    projected : SumstatsTable
        New working table under the panel's identity configuration. The input
        object and source file are not mutated.
    dropped : pandas.DataFrame
        Stable row-level audit records for rejected source rows.

    Raises
    ------
    LDSCInputError
        If required legacy columns or the LD-score directory's SNP identifier metadata are missing,
        or no source row can be projected safely.
    """
    label = table.source_path or table.trait_name or 'sumstats'
    missing = sorted({"SNP", "Z", "N"} - set(table.data.columns))
    if missing:
        raise LDSCInputError(
            f"Cannot use legacy summary-statistics file '{label}': required columns are missing: {missing}. "
            "Legacy regression compatibility requires SNP, Z and N, and A1 and A2 as a pair when present."
        )
    allele_columns = [column for column in ("A1", "A2") if column in table.data.columns]
    if len(allele_columns) == 1:
        present = allele_columns[0]
        absent = "A2" if present == "A1" else "A1"
        raise LDSCInputError(
            f"Cannot use legacy summary-statistics file '{label}': it has {present} but not {absent}. "
            "A1 (the effect allele, the allele the signed statistic refers to) and A2 (the other allele) must come "
            "as a pair. Re-create the file with both allele columns, or without either for h2 and partitioned-h2."
        )
    has_alleles = len(allele_columns) == 2
    panel_snapshot = ldscore_result.config_snapshot
    panel_mode = _snapshot_identity_mode(panel_snapshot)
    if panel_snapshot is None or panel_mode is None:
        raise LDSCInputError(
            "Legacy summary-statistics projection requires an LD-score directory with recorded "
            "snp_identifier metadata. Regenerate or explicitly convert the LD-score directory first."
        )
    panel_source = ldscore_result.baseline_table
    missing_panel = sorted({"CHR", "POS", "SNP"} - set(panel_source.columns))
    if missing_panel:
        raise LDSCInputError(
            f"Legacy summary-statistics projection cannot use the LD-score panel because it is missing {missing_panel}. "
            "Regenerate the LD-score directory."
        )
    # Without alleles, rows project by rsID alone, under the panel's base mode.
    allele_aware = is_allele_aware_mode(panel_mode) and has_alleles
    projected_snapshot = panel_snapshot if has_alleles else replace(
        panel_snapshot, snp_identifier=identity_base_mode(panel_mode))
    if allele_aware and not {"A1", "A2"}.issubset(panel_source.columns):
        raise LDSCInputError(
            f"Legacy summary-statistics projection requires panel A1/A2 for snp_identifier='{panel_mode}', "
            "but the canonical LD-score table does not contain both allele columns."
        )

    source = table.data.copy(deep=False)  # normalized columns are replaced, not modified in place
    source.index = pd.RangeIndex(len(source))
    source["SNP"] = source["SNP"].astype("string").str.strip()
    if has_alleles:
        source["A1"] = source["A1"].astype("string").str.upper().str.strip()
        source["A2"] = source["A2"].astype("string").str.upper().str.strip()
    source_rsids = set(source["SNP"].dropna().astype(str)) - {""}
    panel_snp = panel_source["SNP"].astype("string").str.strip()
    panel_columns = ["CHR", "POS", "SNP", *(["A1", "A2"] if allele_aware else [])]
    in_source = np.flatnonzero(panel_snp.isin(source_rsids).to_numpy())
    panel = panel_source.loc[:, panel_columns].take(in_source)
    panel.index = pd.RangeIndex(len(panel))
    panel["SNP"] = panel_snp.take(in_source).to_numpy()
    if allele_aware:
        panel["A1"] = panel["A1"].astype("string").str.upper().str.strip()
        panel["A2"] = panel["A2"].astype("string").str.upper().str.strip()

    invalid_frq_count = 0
    if "FRQ" in source.columns:
        raw_frq = source["FRQ"]
        numeric_frq = pd.to_numeric(raw_frq, errors="coerce")
        invalid_frq = raw_frq.notna() & (~numeric_frq.between(0.0, 1.0, inclusive="both"))
        invalid_frq_count = int(invalid_frq.sum())
        source["FRQ"] = numeric_frq.where(numeric_frq.between(0.0, 1.0, inclusive="both"))
        if invalid_frq_count:
            LOGGER.warning(
                "Legacy summary-statistics file '%s' has %d nonnumeric, infinite, or out-of-range FRQ value(s); "
                "they were retained as missing because regression does not consume FRQ.",
                table.source_path or table.trait_name or "sumstats",
                invalid_frq_count,
            )

    candidate_counts = panel["SNP"].value_counts(sort=False, dropna=False)
    source_candidate_count = source["SNP"].map(candidate_counts).fillna(0).astype("int64")
    duplicate_source = source["SNP"].notna() & source["SNP"].duplicated(keep=False)
    valid_bases = {"A", "C", "G", "T"}
    reasons = pd.Series(pd.NA, index=source.index, dtype="string")

    def assign_reason(mask: pd.Series, reason: str) -> None:
        reasons.loc[reasons.isna() & mask.fillna(False)] = reason

    assign_reason(duplicate_source, "duplicate_source_rsid")
    if has_alleles:
        assign_reason(
            source["A1"].isna() | source["A2"].isna() | source["A1"].eq("") | source["A2"].eq(""),
            "missing_allele",
        )
        assign_reason(
            ~source["A1"].isin(valid_bases) | ~source["A2"].isin(valid_bases) | source["A1"].eq(source["A2"]),
            "invalid_allele",
        )
    assign_reason(source_candidate_count.eq(0), "missing_panel_rsid")

    eligible = reasons.isna()
    selected_panel: pd.DataFrame
    flip_by_source = pd.Series(False, index=source.index, dtype=bool)
    if allele_aware:
        source_candidates = source.loc[eligible, ["SNP", "A1", "A2"]]
        source_candidates["_legacy_source_row"] = source_candidates.index
        panel_candidates = panel.rename(
            columns={column: f"_legacy_panel_{column}" for column in ["CHR", "POS", "SNP", "A1", "A2"]}
        )
        pairs = source_candidates.merge(
            panel_candidates,
            how="inner",
            left_on="SNP",
            right_on="_legacy_panel_SNP",
            sort=False,
            validate="many_to_many",
        )
        complement_a1 = pairs["A1"].map(BASE_COMPLEMENT)
        complement_a2 = pairs["A2"].map(BASE_COMPLEMENT)
        direct = pairs["A1"].eq(pairs["_legacy_panel_A1"]) & pairs["A2"].eq(pairs["_legacy_panel_A2"])
        complemented = complement_a1.eq(pairs["_legacy_panel_A1"]) & complement_a2.eq(
            pairs["_legacy_panel_A2"]
        )
        swapped = pairs["A1"].eq(pairs["_legacy_panel_A2"]) & pairs["A2"].eq(
            pairs["_legacy_panel_A1"]
        )
        swapped_complement = complement_a1.eq(pairs["_legacy_panel_A2"]) & complement_a2.eq(
            pairs["_legacy_panel_A1"]
        )
        compatible = direct | complemented | swapped | swapped_complement
        compatible_pairs = pairs.take(np.flatnonzero(compatible.to_numpy()))
        # An A/T or C/G pair is its own strand complement, so it matches both
        # ways round and its orientation is unknown: it is never flipped.
        flip = (swapped | swapped_complement) & complement_a1.ne(pairs["A2"])
        compatible_pairs["_legacy_flip"] = flip.loc[compatible_pairs.index].to_numpy()
        compatible_counts = compatible_pairs.groupby("_legacy_source_row", sort=False).size()
        compatible_count_by_source = source.index.to_series().map(compatible_counts).fillna(0).astype("int64")
        assign_reason(eligible & compatible_count_by_source.eq(0), "incompatible_alleles")
        assign_reason(eligible & compatible_count_by_source.gt(1), "ambiguous_panel_mapping")
        accepted_rows = reasons[reasons.isna()].index
        selected_panel = compatible_pairs.loc[
            compatible_pairs["_legacy_source_row"].isin(accepted_rows)
        ].set_index("_legacy_source_row")
        selected_panel = selected_panel.loc[accepted_rows]
        flip_by_source.loc[accepted_rows] = selected_panel["_legacy_flip"].to_numpy()
    else:
        assign_reason(eligible & source_candidate_count.ne(1), "ambiguous_panel_mapping")
        accepted_rows = reasons[reasons.isna()].index
        source_matches = source.loc[accepted_rows, ["SNP"]]
        source_matches["_legacy_source_row"] = source_matches.index
        unique_panel = panel.loc[panel["SNP"].isin(source_matches["SNP"])]
        selected_panel = source_matches.merge(unique_panel, how="left", on="SNP", sort=False, validate="one_to_one")
        selected_panel = selected_panel.set_index("_legacy_source_row").loc[accepted_rows]
    # Checked last, so a row with an earlier reason keeps it.
    assign_reason(~_finite_statistics(source), "nonfinite_statistic")

    accepted_rows = reasons[reasons.isna()].index
    if accepted_rows.empty:
        counts = reasons.dropna().value_counts().to_dict()
        raise LDSCInputError(
            f"Legacy summary-statistics file '{table.source_path or table.trait_name or 'sumstats'}' retained no SNPs "
            f"after projection onto the regression SNPs in the LD-score directory (drop counts: {counts})."
        )
    projected_frame = source.take(accepted_rows.to_numpy())
    panel_column = (lambda column: f"_legacy_panel_{column}") if allele_aware else (lambda column: column)
    projected_frame["CHR"] = selected_panel.loc[accepted_rows, panel_column("CHR")].to_numpy()
    projected_frame["POS"] = selected_panel.loc[accepted_rows, panel_column("POS")].to_numpy()
    projected_frame["SNP"] = selected_panel.loc[accepted_rows, panel_column("SNP")].to_numpy()
    if allele_aware:
        # A/T and C/G rows keep their source alleles, the orientation their Z refers to.
        oriented = projected_frame["A1"].map(BASE_COMPLEMENT).ne(projected_frame["A2"]).to_numpy()
        for column in ("A1", "A2"):
            panel_alleles = selected_panel.loc[accepted_rows, panel_column(column)].to_numpy()
            projected_frame[column] = np.where(oriented, panel_alleles, projected_frame[column].to_numpy())
        flipped_rows = flip_by_source.loc[accepted_rows]
        projected_frame.loc[flipped_rows, "Z"] = -pd.to_numeric(
            projected_frame.loc[flipped_rows, "Z"], errors="raise"
        )
        if "FRQ" in projected_frame.columns:
            valid_flipped_frq = flipped_rows & projected_frame["FRQ"].notna()
            projected_frame.loc[valid_flipped_frq, "FRQ"] = 1.0 - projected_frame.loc[valid_flipped_frq, "FRQ"]
    projected_frame.index = pd.RangeIndex(len(projected_frame))
    projected_frame["POS"] = pd.to_numeric(projected_frame["POS"], errors="raise").astype("int64")
    dropped_mask = reasons.notna()
    drops = source.reindex(columns=["SNP", "A1", "A2"]).loc[dropped_mask]
    drops.insert(0, "source_path", table.source_path)
    drops.insert(0, "trait_name", table.trait_name or Path(table.source_path or "sumstats").name)
    drops["reason"] = reasons.loc[dropped_mask].to_numpy()
    drops["panel_candidate_count"] = source_candidate_count.loc[dropped_mask].to_numpy()
    drops = drops.loc[:, LEGACY_SUMSTATS_DROP_COLUMNS]
    drops.index = pd.RangeIndex(len(drops))
    drop_counts = drops["reason"].value_counts(sort=False)
    if not drop_counts.empty:
        LOGGER.warning(
            "Legacy summary-statistics file '%s' projection dropped %d row(s): %s.",
            table.source_path or table.trait_name or "sumstats",
            len(drops),
            ", ".join(f"{reason}={int(count)}" for reason, count in drop_counts.items()),
        )
    LOGGER.info(
        "Projected %d legacy summary-statistics row(s) from '%s' onto panel SNP identifier mode %s; dropped %d.",
        len(projected_frame),
        table.source_path or table.trait_name or "sumstats",
        projected_snapshot.snp_identifier,
        len(drops),
    )
    return replace(table, data=projected_frame, has_alleles=has_alleles, config_snapshot=projected_snapshot), drops


def _sumstats_identity_mode(
    sumstats_table: SumstatsTable,
    *,
    fallback_mode: str,
) -> str:
    """Return the mode represented by one summary-statistics input for regression."""
    return _snapshot_identity_mode(sumstats_table.config_snapshot) or normalize_snp_identifier_mode(fallback_mode)


def _ldscore_identity_mode(
    ldscore_result: LDScoreSource | LDScoreResult,
    runner_config: GlobalConfig,
    *,
    fallback_mode: str | None = None,
) -> str:
    """Return the mode represented by one LD-score input for regression."""
    return (
        _snapshot_identity_mode(ldscore_result.config_snapshot)
        or normalize_snp_identifier_mode(fallback_mode or runner_config.snp_identifier)
    )


def _resolve_h2_identity(
    sumstats_table: SumstatsTable,
    ldscore_result: LDScoreSource | LDScoreResult,
    runner_config: GlobalConfig,
    regression_config: RegressionConfig,
):
    """Resolve one summary-statistics table and one LD-score result for h2 merging."""
    ldscore_mode = _ldscore_identity_mode(ldscore_result, runner_config)
    sumstats_mode = _sumstats_identity_mode(sumstats_table, fallback_mode=ldscore_mode)
    # Allele-less legacy text has no allele-aware alternative, so its projection is downgraded without the flag.
    return resolve_regression_identity_mode(
        sumstats_mode,
        ldscore_mode,
        allow_identity_downgrade=regression_config.allow_identity_downgrade or _is_allele_less_legacy(sumstats_table),
    )


def _resolve_rg_identity(
    sumstats_table_1: SumstatsTable,
    sumstats_table_2: SumstatsTable,
    ldscore_result: LDScoreSource | LDScoreResult,
    runner_config: GlobalConfig,
    regression_config: RegressionConfig,
):
    """Resolve trait 1, trait 2, and LD-score SNP identifier modes for rg."""
    ldscore_mode = _ldscore_identity_mode(ldscore_result, runner_config)
    trait_1_mode = _sumstats_identity_mode(sumstats_table_1, fallback_mode=ldscore_mode)
    trait_2_mode = _sumstats_identity_mode(sumstats_table_2, fallback_mode=ldscore_mode)
    pairwise = [
        resolve_regression_identity_mode(
            trait_1_mode,
            ldscore_mode,
            allow_identity_downgrade=regression_config.allow_identity_downgrade,
        ),
        resolve_regression_identity_mode(
            trait_2_mode,
            ldscore_mode,
            allow_identity_downgrade=regression_config.allow_identity_downgrade,
        ),
        resolve_regression_identity_mode(
            trait_1_mode,
            trait_2_mode,
            allow_identity_downgrade=regression_config.allow_identity_downgrade,
        ),
    ]
    if len({trait_1_mode, trait_2_mode, ldscore_mode}) == 1:
        return pairwise[0]
    effective_mode = identity_base_mode(trait_1_mode)
    return replace(pairwise[0], effective_mode=effective_mode, downgrade_applied=True)


def _prepare_regression_identity_table(
    frame: pd.DataFrame,
    mode: str,
    *,
    context: str,
    logger: logging.Logger,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Clean duplicate effective identities before downgraded regression merge."""
    cleanup = clean_identity_artifact_table(
        frame,
        mode,
        context=context,
        stage="regression_identity_downgrade",
        logger=logger,
    )
    return cleanup.cleaned, cleanup.dropped


def _key_regression_rows(frame: pd.DataFrame, mode: str, *, cleaned: bool, context: str) -> tuple[pd.DataFrame, int]:
    """Key one side of a regression merge for ``mode``.

    ``cleaned`` first applies the identity-downgrade cleanup. Allele-aware modes
    add the effective identity key and coordinate-family base modes the
    ``CHR:POS`` key; ``rsid`` merges on ``SNP`` as given. Returns the keyed rows
    and the number of duplicate effective-key rows the cleanup dropped.
    """
    dropped = 0
    if cleaned:
        frame, dropped_rows = _prepare_regression_identity_table(frame, mode, context=context, logger=LOGGER)
        dropped = len(dropped_rows)
    if is_allele_aware_mode(mode):
        frame = _with_effective_identity_key(frame, mode, context=context)
    elif identity_mode_family(mode) != "rsid":
        frame = _with_chr_pos_key(frame, context=context)
    return frame, dropped


def _merge_key_column(mode: str) -> str:
    """Return the column ``_key_regression_rows`` merges on under ``mode``."""
    if is_allele_aware_mode(mode):
        return REGRESSION_IDENTITY_KEY_COLUMN
    return "SNP" if identity_mode_family(mode) == "rsid" else CHR_POS_KEY_COLUMN


def _with_effective_identity_key(frame: pd.DataFrame, mode: str, *, context: str) -> pd.DataFrame:
    """Return the rows with an effective regression identity key, with that key added.

    The result shares the input's columns unless rows are dropped; callers
    only read it (merges copy).
    """
    keyed = frame.copy(deep=False)
    keyed[REGRESSION_IDENTITY_KEY_COLUMN] = effective_merge_key_series(keyed, mode, context=context)
    present = keyed[REGRESSION_IDENTITY_KEY_COLUMN].notna().to_numpy()
    if not present.all():
        keyed = keyed.take(np.flatnonzero(present))
    keyed.index = pd.RangeIndex(len(keyed))
    return keyed


def _with_chr_pos_key(frame: pd.DataFrame, *, context: str) -> pd.DataFrame:
    """Return a copy with a private canonical ``CHR:POS`` merge key."""
    chr_col, pos_col = infer_chr_pos_columns(frame.columns, context=context)
    keyed, _report = build_chr_pos_key_frame(
        frame,
        context=context,
        chr_col=chr_col,
        pos_col=pos_col,
        drop_missing=True,
        logger=LOGGER,
    )
    return keyed


def _ld_column_variances(merged: pd.DataFrame, columns: Sequence[str], known: dict | None = None) -> pd.Series:
    """Return each column's variance, reusing ``known`` values for unchanged columns.

    Each variance is the pandas column variance, which equals the value the
    column receives within ``DataFrame.var()``.
    """
    known = known or {}
    return pd.Series({column: known[column] if column in known else merged[column].var() for column in columns})


def _allele_pair_flags(columns: Sequence[pd.Series]) -> tuple[np.ndarray, np.ndarray]:
    """Return per-row (compatible, flip) flags for concatenated allele pairs.

    The pair string is the four alleles concatenated; it is compatible when
    it is in ``MATCH_ALLELES`` and flips when ``FLIP_ALLELES`` says so. The string is built once per distinct allele
    combination with the same pandas operations a per-row build would use.
    """
    codes, sizes = [], []
    for column in columns:
        column_codes, uniques = pd.factorize(column, use_na_sentinel=True)
        codes.append(column_codes + 1)
        sizes.append(len(uniques) + 1)
    if np.prod(sizes, dtype=object) < 2 ** 62:
        combined = np.zeros(len(codes[0]), dtype=np.int64)
        for column_codes, size in zip(codes, sizes):
            combined = combined * size + column_codes
        _, first, inverse = np.unique(combined, return_index=True, return_inverse=True)
    else:
        _, first, inverse = np.unique(np.stack(codes, axis=1), axis=0, return_index=True, return_inverse=True)
    inverse = inverse.reshape(-1)
    parts = [column.iloc[first].reset_index(drop=True) for column in columns]
    pairs = parts[0] + parts[1] + parts[2] + parts[3]
    compatible = pairs.isin(MATCH_ALLELES).to_numpy(dtype=bool)
    flip = np.zeros(len(pairs), dtype=bool)
    flip[compatible] = pairs[compatible].map(FLIP_ALLELES).astype(bool).to_numpy()
    return compatible[inverse], flip[inverse]


def _same_strand_ambiguous_pair(merged: pd.DataFrame) -> np.ndarray:
    """Mark rg rows where both traits list the same A/T or C/G pair, in either order."""
    a1, a2, a1x, a2x = (merged[column] for column in ("A1", "A2", "A1x", "A2x"))
    same_pair = (a1.eq(a1x) & a2.eq(a2x)) | (a1.eq(a2x) & a2.eq(a1x))
    return (a1.map(BASE_COMPLEMENT).eq(a2) & same_pair).to_numpy(dtype=bool, na_value=False)


def _nonfinite_columns(frame: pd.DataFrame, columns: Sequence[str], snps: pd.Series) -> list[str]:
    """Describe each column with missing or infinite values: its name, row count and up to five example SNPs."""
    problems = []
    for column in columns:
        values = frame[column]
        # Float columns are tested as stored; widening them would not change finiteness.
        values = (values.to_numpy() if isinstance(values.dtype, np.dtype) and values.dtype.kind == "f"
                  else values.to_numpy(dtype=np.float64))
        bad = ~np.isfinite(values)
        if bad.any():
            count = int(bad.sum())
            examples = ", ".join(snps.astype(str).to_numpy()[bad][:5])
            problems.append(f"'{column}' ({count} row{'' if count == 1 else 's'}; for example {examples})")
    return problems


def _select_model_columns(
    merged: pd.DataFrame,
    ref_ld_columns: Sequence[str],
    count_records: Iterable[dict[str, Any]] | dict[str, dict[str, Any]],
    use_common_counts: bool,
    *,
    command: str,
    known_variances: dict | None = None,
) -> dict[str, object]:
    """Drop zero-variance LD-score columns and select the requested count vectors.

    Each command calls this after its own SNP alignment, so variances are
    taken on the rows that enter the fit. A dropped column's values are
    constant there and identify no coefficient; its name is kept for the
    fitted-model metadata. Returns the ``RegressionDataset`` fields
    ``merged``, ``retained_ld_columns``, ``dropped_zero_variance_ld_columns``,
    ``reference_snp_count_totals`` and ``count_key_used_for_regression``.
    ``merged`` borrows the input's columns; callers must not modify them.
    ``known_variances`` supplies variances already computed on these rows.
    """
    variances = _ld_column_variances(merged, ref_ld_columns, known_variances)
    dropped = variances.index[variances == 0].tolist()
    retained = [column for column in ref_ld_columns if column not in dropped]
    if not retained:
        raise LDSCInputError(
            f"{command} regression cannot run because all retained LD-score columns have zero variance on the "
            "merged regression SNPs. Most likely the LD-score directory contains constant annotations on the "
            "retained SNP set. Use LD scores with at least one varying annotation or broaden the regression SNP set."
        )
    if dropped:
        merged = pd.DataFrame({name: merged[name].to_numpy(copy=False) for name in merged if name not in dropped},
                              copy=False)
        LOGGER.warning(
            "Dropped %d zero-variance LD-score column(s) after merging with summary statistics. "
            "These annotations are constant on the regression SNP set and receive no coefficient; "
            "names are retained in fitted-model metadata (dropped_zero_variance_ld_columns).",
            len(dropped),
        )
    count_totals = _count_totals_for_columns(count_records, ref_ld_columns)
    count_key = _select_count_key(count_totals, use_common_counts, ref_ld_columns)
    keep = [index for index, column in enumerate(ref_ld_columns) if column in retained]
    return {
        "merged": merged,
        "retained_ld_columns": retained,
        "dropped_zero_variance_ld_columns": dropped,
        "reference_snp_count_totals": {key: np.asarray(values)[keep] for key, values in count_totals.items()},
        "count_key_used_for_regression": count_key,
    }


def _select_count_key(
    count_totals: dict[str, np.ndarray],
    use_common_counts: bool,
    columns: Sequence[str] = (),
) -> str:
    """Return the requested count vector key; the other count kind is never substituted."""
    if use_common_counts:
        if COMMON_COUNT_KEY in count_totals:
            return COMMON_COUNT_KEY
        raise LDSCInputError(
            "Regression requested common-SNP counts (`--count-kind common`, the default), but the LD-score "
            "directory does not record a common-SNP count for every requested LD-score column "
            f"({_preview_names(columns)}). All-SNP counts are never substituted. Most likely the metadata.json "
            "was written by an older LDSC3 version or edited by hand. Regenerate or reconvert the LD-score "
            "directory, or pass `--count-kind all` (`RegressionConfig(use_common_counts=False)`) to scale "
            "estimates by all-SNP counts."
        )
    if ALL_COUNT_KEY in count_totals:
        return ALL_COUNT_KEY
    raise LDSCInputError(
        "Regression requested all-SNP counts, but all-SNP counts are unavailable for LD-score column(s) "
        f"{list(columns)}. This is expected for an imported unpartitioned legacy suite when one or more "
        "chromosome .M files were missing. Use the default common-SNP counts or convert a suite with complete .M files."
    )


def _preview_names(names: Sequence[str]) -> str:
    """Join the first five names for a message and count the rest."""
    names = list(names)
    preview = ", ".join(names[:5])
    return f"{preview}, ... ({len(names) - 5} more)" if len(names) > 5 else preview


def _validate_partitioned_query_columns(ldscore_result: LDScoreSource | LDScoreResult, query_columns: Sequence[str]) -> list[str]:
    """Validate that requested cell-type query columns exist in the LD-score result.

    The functional regime (no query columns) is decided by the caller; this
    helper is only used on the cell-type path, so it just checks that every
    requested query column is available.
    """
    requested = list(query_columns)
    available = list(ldscore_result.query_columns)
    available_set = set(available)
    missing = [column for column in requested if column not in available_set]
    if missing:
        raise LDSCInputError(
            f"partitioned-h2 could not find requested query annotation columns {missing} in the LD-score directory. "
            f"Most likely the query annotation name was misspelled or the LD-score directory was generated from "
            f"different query inputs. Choose one of the available query annotations: {available}."
        )
    return requested


def _assemble_regression_ldscore_table(ldscore_result: LDScoreSource | LDScoreResult) -> pd.DataFrame:
    """Return the baseline LD-score table as a shallow copy with a RangeIndex.

    Callers add columns to it but never modify the shared values in place.
    """
    table = ldscore_result.baseline_table.copy(deep=False)
    table.index = pd.RangeIndex(len(table))
    return table


def _count_totals_for_columns(count_records, columns: Sequence[str]) -> dict[str, np.ndarray]:
    """Return LDSC count vectors aligned to ``columns`` using metadata records.

    ``count_records`` is a sequence of records or a mapping from column name to record.
    """
    records_by_column = (count_records if isinstance(count_records, dict)
                         else {str(record["column"]): record for record in count_records})
    missing = [column for column in columns if column not in records_by_column]
    if missing:
        raise LDSCInputError(
            f"Regression cannot align LD-score count metadata: count records are missing for columns {missing}. "
            "Most likely the LD-score metadata was written by an older version or edited separately from "
            "the parquet tables. Regenerate the LD-score directory with the current `ldsc ldscore`."
        )
    count_totals: dict[str, np.ndarray] = {}
    all_values = [records_by_column[column].get("all_reference_snp_count") for column in columns]
    if all(value is not None and pd.notna(value) for value in all_values):
        all_counts = np.asarray([float(value) for value in all_values], dtype=np.float64)
        if not np.isfinite(all_counts).all():
            raise LDSCInputError(
                "Regression found non-finite all-SNP count metadata in the LD-score directory."
                " Restore consistent finite count metadata from the same LD-score run, or regenerate the "
                "LD-score directory. Do not replace missing or invalid counts with zero."
            )
        count_totals[ALL_COUNT_KEY] = all_counts
    if all("common_reference_snp_count" in records_by_column[column] for column in columns):
        common_values = [records_by_column[column].get("common_reference_snp_count") for column in columns]
        if any(value is None or pd.isna(value) for value in common_values):
            raise LDSCInputError(
                "Regression found a malformed common-SNP count record with a missing value. "
                "Reconvert or regenerate the LD-score directory; common-SNP counts never fall back silently."
            )
        common_counts = np.asarray([float(value) for value in common_values], dtype=np.float64)
        if not np.isfinite(common_counts).all():
            raise LDSCInputError(
                "Regression found non-finite common-SNP count metadata in the LD-score directory."
                " Restore consistent finite count metadata from the same LD-score run, or regenerate the "
                "LD-score directory. Do not replace missing or invalid counts with zero."
            )
        count_totals[COMMON_COUNT_KEY] = common_counts
    return count_totals


def _require_rg_alleles(tables: Sequence[SumstatsTable]) -> None:
    """Reject traits without ``A1``/``A2``; rg needs them to align signed statistics."""
    missing = [table for table in tables if not {"A1", "A2"}.issubset(table.data.columns)]
    if missing:
        listed = ", ".join(f"{_trait_label(table)} ('{table.source_path or 'in-memory table'}')" for table in missing)
        raise LDSCInputError(
            "rg needs A1 (the effect allele, the allele the signed statistic refers to) and A2 in every summary-"
            "statistics input to align the signed statistics of the two traits, but these inputs have no A1/A2 "
            f"columns: {listed}. Most likely `ldsc munge-sumstats` was not given --a1/--a2 and could not infer "
            "them, for example because the raw allele columns were REF/ALT. Re-run `ldsc munge-sumstats` for them "
            "with --a1 <effect-allele column> --a2 <other column>. Other causes & fixes: "
            f"{_RG_ALLELES_DOC}"
        )


def _trait_label(table: SumstatsTable) -> str:
    """Return the public trait label for a summary-statistics table."""
    if table.trait_name:
        return str(table.trait_name)
    if table.source_path:
        return Path(table.source_path).name
    return "trait"


def _log_effective_regression_identity(
    sumstats_tables: Sequence[SumstatsTable],
    ldscore_result: LDScoreSource | LDScoreResult,
    runner_config: GlobalConfig,
    regression_config: RegressionConfig,
) -> None:
    """Log the SNP identifier mode regression will actually use, and its source.

    The GlobalConfig banner reports ``snp_identifier`` as a CLI-level default.
    Regression instead derives the effective mode from each input artifact's
    recorded provenance and only falls back to that default when an artifact
    omits it. This line names the resolved mode, the per-input sources, and any
    allele-awareness downgrade so the banner default is not mistaken for the
    mode in force. It is informational only: incompatible inputs are left for
    the estimation step to reject with its authoritative error.
    """
    default_mode = normalize_snp_identifier_mode(runner_config.snp_identifier)
    used_default = False

    def describe(label: str, snapshot: GlobalConfig | None) -> tuple[str, str]:
        nonlocal used_default
        recorded = _snapshot_identity_mode(snapshot)
        if recorded is None:
            used_default = True
            return default_mode, f"{label}={default_mode} (no provenance recorded; using GlobalConfig default)"
        return recorded, f"{label}={recorded} (from provenance)"

    modes: list[str] = []
    sources: list[str] = []
    ld_mode, ld_source = describe("LD-score", ldscore_result.config_snapshot)
    modes.append(ld_mode)
    sources.append(ld_source)
    for table in sumstats_tables:
        if _is_allele_less_legacy(table):
            # Without alleles, projection gives a legacy table the base form of the LD-score directory's mode.
            mode = identity_base_mode(ld_mode)
            source = f"{_trait_label(table)}={mode} (legacy projection by rsID, without alleles)"
        elif _has_legacy_sumstats_source(table):
            # Projection gives a legacy table the LD-score directory's mode.
            mode, source = ld_mode, f"{_trait_label(table)}={ld_mode} (legacy projection onto the LD-score panel)"
        else:
            mode, source = describe(_trait_label(table), table.config_snapshot)
        modes.append(mode)
        sources.append(source)

    source_details = "; ".join(sources) if len(sources) <= 10 else f"{len(sources)} input artifacts; see input provenance and per-pair metadata"
    if len(sources) > 10 and LOGGER.isEnabledFor(logging.DEBUG):
        LOGGER.debug("Per-input SNP modes: %s", "; ".join(sources))
    effective = modes[0]
    downgrade = False
    try:
        for mode in modes[1:]:
            compat = resolve_regression_identity_mode(
                effective, mode, allow_identity_downgrade=regression_config.allow_identity_downgrade
                or any(_is_allele_less_legacy(table) for table in sumstats_tables),
            )
            effective = compat.effective_mode
            downgrade = downgrade or compat.downgrade_applied
    except LDSCUsageError:
        # Inputs are incompatible; the estimation step raises the authoritative
        # error. Still surface per-input modes so the conflict is visible early.
        LOGGER.info(
            "Regression SNP identifier could not be resolved: the inputs disagree. "
            "Per-input modes: %s. The GlobalConfig snp_identifier %r shown above is "
            "applied only to inputs that have no recorded provenance.",
            source_details,
            default_mode,
        )
        return

    # The GlobalConfig snp_identifier in the banner above is only a fallback. State
    # plainly whether it was used at all, so it is never mistaken for the mode in force.
    if used_default:
        default_clause = (
            f"The GlobalConfig snp_identifier {default_mode!r} shown above was used only as a "
            f"fallback for inputs without recorded provenance"
        )
    else:
        default_clause = (
            f"The GlobalConfig snp_identifier {default_mode!r} shown above was NOT used; "
            f"this mode comes entirely from input provenance"
        )
    downgrade_clause = (
        " An allele-awareness downgrade was applied to reach this shared base mode."
        if downgrade
        else ""
    )
    LOGGER.info(
        "SNP identifier used for this regression: %s.%s %s. Per-input modes: %s.",
        effective,
        downgrade_clause,
        default_clause,
        source_details,
    )


def _add_panel_coordinates(merged: pd.DataFrame, ldscore_result: LDScoreSource | LDScoreResult) -> np.ndarray:
    """Replace ``_ldsc_source_row`` with the LD-score panel's coordinates and return the source rows.

    Delete-block diagnostics report genomic spans from ``_ldsc_panel_CHR`` and
    ``_ldsc_panel_POS`` even when summary statistics omit coordinates (rsID
    inputs) or carry missing placeholders. These columns never reach outputs.
    """
    source_rows = merged.pop("_ldsc_source_row").to_numpy(dtype=np.int64)
    try:
        coordinate_columns = infer_chr_pos_columns(ldscore_result.baseline_table.columns, context="LD-score diagnostics")
    except LDSCInputError:
        coordinate_columns = ()
    for coordinate, source_column in zip(("CHR", "POS"), coordinate_columns):
        merged[f"_ldsc_panel_{coordinate}"] = ldscore_result.baseline_table[source_column].iloc[source_rows].to_numpy()
    return source_rows
