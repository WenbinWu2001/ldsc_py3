"""Gene LD-score index fixtures written through the production staging writer."""

from __future__ import annotations

from pathlib import Path

from ldsc._ldscore.gene_index import artifact as gene_index_artifact
from ldsc._ldscore.gene_index import build as gene_index_build
from ldsc._ldscore.gene_index.artifact import (
    StagedIndexChromosome,
    _finalize_staged_index_artifact,
    _index_chromosome_metadata,
    _remove_gene_index_transaction,
    _stage_index_chromosome,
    calculate_index_id,
)
from ldsc._ldscore.outputs import snp_universe_policy_record


def narrow_index_chromosomes(monkeypatch, chromosomes) -> None:
    """Let one test build and load indexes that cover only ``chromosomes``.

    Published indexes must cover autosomes 1 to 22; small fixtures cover a
    few chromosomes. The narrowed set governs the build's chromosome list and
    the loader's completeness check until the test ends.
    """
    monkeypatch.setattr(
        gene_index_artifact, "INDEX_CHROMOSOMES", tuple(str(chrom) for chrom in chromosomes)
    )


def write_gene_ldscore_index(index_dir, *, index_identity, gene_catalog, chromosomes) -> Path:
    """Publish in-memory chromosome records as a gene LD-score index.

    The steps are those of ``build-gene-ldscore-index`` after its chromosome
    computations: create the staged index in the build directory, stage each
    chromosome's payload, write the coordinator metadata, publish with reload
    validation, then remove the build directory. ``chromosomes`` maps labels to ``IndexChromosomeData`` in index
    order and must equal the narrowed index chromosome set. Root
    ``annotation_types`` classifies every baseline column as ``binary``;
    ``snp_universe_policy`` records no reference-panel or regression-SNP file
    and no region subtraction. ``index_dir`` must be absent or empty; publication-transaction
    behavior is tested through the build command itself.
    """
    destination = Path(index_dir)
    staged_index = gene_index_build._create_gene_index_stage(destination)
    try:
        staged = {}
        for chrom, record in chromosomes.items():
            component_metadata = _index_chromosome_metadata(
                chrom,
                record,
                snp_identifier=str(index_identity["snp_identifier"]),
                genome_build=str(index_identity["genome_build"]),
            )
            _stage_index_chromosome(staged_index, chrom, record)
            staged[chrom] = StagedIndexChromosome(chrom, {}, component_metadata)
        records = list(chromosomes.values())
        baseline_columns = component_metadata["baseline_columns"]
        written = sum(len(record.baseline_rows) for record in records)
        _finalize_staged_index_artifact(
            staged_index,
            index_id=calculate_index_id(index_identity),
            index_identity=index_identity,
            gene_catalog=gene_catalog,
            chromosomes=staged,
            annotation_types=dict.fromkeys(baseline_columns, "binary"),
            snp_universe_policy=snp_universe_policy_record(
                ref_panel_snps_file=None, regr_snps_file=None, regression_regions=None,
                retained_snp_count=sum(record.total_reference_snps_all for record in records),
                selected_snp_count_before_region_exclusion=written,
                region_excluded_snp_count=0, written_snp_count=written),
        )
        return gene_index_build._publish_staged_gene_ldscore_index(staged_index, destination)
    finally:
        _remove_gene_index_transaction(gene_index_build._gene_index_build_dir(destination))
