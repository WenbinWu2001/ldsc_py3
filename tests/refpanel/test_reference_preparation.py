from tests.annotation_fixtures import build_annotation_fixture, make_annotation_bundle, fixture_metadata, fixture_values, fixture_ids
"""Reference policies and numerical projection share the production preparation seam."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ldsc import AnnotationBundle, GlobalConfig, LDScoreCalculator, LDScoreConfig, RefPanelConfig
from ldsc._ldscore import compute as ldscore
from ldsc._refpanel.prepared import ChromosomeAnnotations
from ldsc._refpanel.r2_format import SortedR2BlockReader
from ldsc._refpanel.adapters import ParquetR2RefPanel
from ldsc._refpanel.r2_format import write_r2_parquet, write_runtime_metadata_sidecar
from ldsc._core.snp_identity import sidecar_identity_sha256
from ldsc.errors import LDSCInputError


# write_panel stores R2=0.5 as round(0.5 * 32767); readers decode stored / 32767 in float32.
STORED_R2 = float(np.float32(round(0.5 * 32767)) / np.float32(32767))


def write_panel(root, *, chrom="1", mode="rsid", ld_window=("snps", 10.0), min_r2=0.0):
    """Write one package-format chromosome with the real sidecar and R2 writers.

    ``ld_window`` is the recorded build window; a run must use its unit.
    ``min_r2`` is recorded as the build-time pair floor.
    """
    root.mkdir(parents=True, exist_ok=True)
    metadata = pd.DataFrame({
        "CHR": [chrom, chrom], "POS": [100, 300], "SNP": ["rs1", "rs3"],
        "A1": ["A", "A"], "A2": ["C", "C"], "CM": [0.1, 0.3], "MAF": [0.04, 0.3],
    })
    write_runtime_metadata_sidecar(metadata, root / f"chr{chrom}_meta.tsv.gz",
                                   genome_build="hg19", snp_identifier=mode)
    write_r2_parquet(
        pair_chunks=[(np.array([0]), np.array([1]), np.array([0.5], dtype=np.float32), np.array([1], dtype=np.int8))],
        path=root / f"chr{chrom}_r2.parquet", genome_build="hg19", n_samples=200, snp_identifier=mode,
        n_snps=2, sidecar_identity_sha256=sidecar_identity_sha256(metadata),
        ld_window_mode=ld_window[0], ld_window_value=ld_window[1], min_r2=min_r2,
    )
    return metadata


def public_bundle(metadata, values=None):
    values = pd.DataFrame({"base": np.ones(len(metadata))}) if values is None else values
    return make_annotation_bundle(metadata=metadata, baseline_annotations=values,
                            query_annotations=pd.DataFrame(index=metadata.index),
                            baseline_columns=list(values.columns), query_columns=[],
                            chromosomes=metadata["CHR"].astype(str).unique().tolist(), source_summary={})


def test_package_r2_is_read_as_unbiased_through_public_calculation(tmp_path):
    # Package panels store already-unbiased R2 quantized as round(R2 * 32767), so
    # each SNP scores itself (1) plus its neighbour's decoded R2, with no
    # second sample-size correction.
    metadata = write_panel(tmp_path)
    config = GlobalConfig(snp_identifier="rsid")
    panel = ParquetR2RefPanel(config, RefPanelConfig(backend="parquet_r2", r2_dir=tmp_path))
    result = LDScoreCalculator().compute_chromosome(
        "1", public_bundle(metadata), panel, LDScoreConfig(ld_wind_snps=1, whole_chromosome_ok=True), config,
    )
    np.testing.assert_allclose(result.baseline_table[["base", "regression_ld_scores"]], 1.0 + STORED_R2, rtol=1e-7)


@pytest.mark.parametrize("alleles_in_annotation", [False, True])
def test_prepared_annotations_share_reference_rows_alleles_cm_maf_and_counts(tmp_path, alleles_in_annotation):
    mode = "chr_pos_allele_aware"
    reference = write_panel(tmp_path, mode=mode, ld_window=("cm", 1.0))
    metadata = pd.DataFrame({"CHR": ["1"] * 3, "SNP": ["rs1", "missing", "rs3"],
                             "POS": [100, 200, 300], "CM": [9.0] * 3, "MAF": [0.9] * 3})
    if alleles_in_annotation:
        metadata = metadata.assign(A1="A", A2="C")
    values = pd.DataFrame({"base": [1.0, 99.0, 3.0], "flag": [True, True, False]})
    config = GlobalConfig(snp_identifier=mode, genome_build="hg19")
    panel = ParquetR2RefPanel(config, RefPanelConfig(backend="parquet_r2", r2_dir=tmp_path))
    bundle = ChromosomeAnnotations(metadata, values, ["base", "flag"], [])
    with panel.prepare_chromosome("1", bundle, LDScoreConfig(ld_wind_cm=0.3, whole_chromosome_ok=True)) as prepared:
        assert prepared.annotations.read().dtype == np.float32
        np.testing.assert_array_equal(prepared.annotations.read(), [[1, 1], [3, 0]])
        pd.testing.assert_frame_equal(prepared.metadata[reference.columns], reference)
        result = ldscore.compute_chromosome("1", prepared, snp_identifier=mode, snp_batch_size=3)
        r2 = STORED_R2
        np.testing.assert_allclose(result.ld_scores, [[1 + 3 * r2, 1], [3 + r2, r2]], rtol=1e-7)
        np.testing.assert_array_equal(result.M, [4, 1])
        np.testing.assert_array_equal(result.M_5_50, [3, 0])
        reader = prepared.reader
    assert reader._pf.closed
    np.testing.assert_array_equal(metadata["CM"], [9.0] * 3)


def test_allele_mismatch_is_rejected_before_reader_creation(tmp_path):
    metadata = write_panel(tmp_path, mode="chr_pos_allele_aware").assign(A2="G")
    config = GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg19")
    panel = ParquetR2RefPanel(config, RefPanelConfig(backend="parquet_r2", r2_dir=tmp_path))
    with pytest.raises(LDSCInputError, match="retained no reference SNPs"):
        LDScoreCalculator().compute_chromosome(
            "1", public_bundle(metadata), panel, LDScoreConfig(ld_wind_snps=1), config,
        )


def test_parquet_preparation_releases_reader_on_calculation_failure(tmp_path, monkeypatch):
    metadata = write_panel(tmp_path)
    config = GlobalConfig(snp_identifier="rsid")
    panel = ParquetR2RefPanel(config, RefPanelConfig(backend="parquet_r2", r2_dir=tmp_path))
    bundle = ChromosomeAnnotations(metadata, pd.DataFrame({"base": [1, 1]}), ["base"], [])

    def fail_decode(self, table, row_group_index):
        raise LDSCInputError("injected decode failure")

    monkeypatch.setattr(SortedR2BlockReader, "_decode_index_batch", fail_decode)
    with pytest.raises(LDSCInputError, match="injected decode failure"):
        with panel.prepare_chromosome("1", bundle, LDScoreConfig(ld_wind_snps=1, whole_chromosome_ok=True)) as prepared:
            reader = prepared.reader
            ldscore.compute_chromosome("1", prepared, snp_identifier="rsid", snp_batch_size=3)
    assert reader._pf.closed


def test_public_calculation_uses_panel_maf_filter_for_aligned_annotations(tmp_path):
    metadata = write_panel(tmp_path)
    config = GlobalConfig(snp_identifier="rsid")
    panel = ParquetR2RefPanel(config, RefPanelConfig(backend="parquet_r2", r2_dir=tmp_path, maf_min=0.1))
    result = LDScoreCalculator().compute_chromosome(
        "1", public_bundle(metadata, pd.DataFrame({"base": [2.0, 3.0]})), panel,
        LDScoreConfig(ld_wind_snps=1, whole_chromosome_ok=True), config,
    )
    assert result.baseline_table["SNP"].tolist() == ["rs3"]
    np.testing.assert_allclose(result.baseline_table[["base", "regression_ld_scores"]], [[3.0, 1.0]])
    assert result.count_records[0]["common_reference_snp_count"] == 3.0


def test_public_run_fails_when_the_grid_misses_a_whole_chromosome(tmp_path):
    absent = write_panel(tmp_path, chrom="1").assign(SNP=["absent1", "absent3"])
    present = write_panel(tmp_path, chrom="22")
    config = GlobalConfig(snp_identifier="rsid")
    panel = ParquetR2RefPanel(config, RefPanelConfig(backend="parquet_r2", r2_dir=tmp_path))
    bundle = public_bundle(pd.concat([absent, present], ignore_index=True))
    with pytest.raises(LDSCInputError, match=r"annotation SNP grid for chromosome 1 \(2 annotation rows"):
        LDScoreCalculator().run(bundle, panel, LDScoreConfig(ld_wind_snps=1, whole_chromosome_ok=True), config)


def test_parquet_panel_with_a_pair_floor_is_rejected_at_preparation(tmp_path):
    # The LDScoreCalculator Python API skips the CLI preflight, so preparation
    # itself must refuse a panel that omitted pairs below --min-r2.
    metadata = write_panel(tmp_path, min_r2=0.05)
    panel = ParquetR2RefPanel(GlobalConfig(snp_identifier="rsid"),
                              RefPanelConfig(backend="parquet_r2", r2_dir=tmp_path))
    with pytest.raises(LDSCInputError, match="omits SNP pairs with unbiased R2 below 0.05"):
        panel.prepare_chromosome(
            "1", ChromosomeAnnotations(metadata, pd.DataFrame({"base": [1, 1]}), ["base"], []),
            LDScoreConfig(ld_wind_snps=1, whole_chromosome_ok=True),
        )


def test_rsid_preparation_of_a_build_tagged_panel_raises_no_genome_build_warning(tmp_path, recwarn):
    # Panels always record their build; rsID-family runs ignore it without
    # the GlobalConfig warning meant for a user-supplied genome_build.
    metadata = write_panel(tmp_path)
    panel = ParquetR2RefPanel(GlobalConfig(snp_identifier="rsid"),
                              RefPanelConfig(backend="parquet_r2", r2_dir=tmp_path))
    with panel.prepare_chromosome(
        "1", ChromosomeAnnotations(metadata, pd.DataFrame({"base": [1, 1]}), ["base"], []),
        LDScoreConfig(ld_wind_snps=1, whole_chromosome_ok=True),
    ):
        pass
    assert not [warning for warning in recwarn if "genome_build is set" in str(warning.message)]


def test_parquet_genetic_maps_are_ignored_at_preparation(tmp_path, caplog):
    metadata = write_panel(tmp_path, ld_window=("cm", 1.0))
    panel = ParquetR2RefPanel(GlobalConfig(snp_identifier="rsid"), RefPanelConfig(
        backend="parquet_r2", r2_dir=tmp_path, genetic_map_hg19_sources="unused-map.tsv",
    ))
    with panel.prepare_chromosome(
        "1", ChromosomeAnnotations(metadata, pd.DataFrame({"base": [1, 1]}), ["base"], []),
        LDScoreConfig(ld_wind_cm=0.3, whole_chromosome_ok=True),
    ) as prepared:
        np.testing.assert_allclose(prepared.metadata["CM"], [0.1, 0.3])
    assert "Ignoring --genetic-map" in caplog.text


def test_plink_preparation_maps_sorted_annotations_to_physical_bed_columns(tmp_path):
    pytest.importorskip("bitarray")
    from tests.ref_panel_helpers import prepare_plink
    from argparse import Namespace
    source = Path(__file__).parents[1] / "fixtures" / "plink" / "plink"
    metadata = pd.read_csv(source.with_suffix(".bim"), sep=r"\s+", header=None,
                           names=["CHR", "SNP", "CM", "POS", "A1", "A2"])
    # Reverse the genomic positions while preserving physical BIM/BED column correspondence.
    metadata["POS"] = metadata["POS"].iloc[::-1].to_numpy()
    prefix = tmp_path / "panel"
    metadata.to_csv(prefix.with_suffix(".bim"), sep="\t", index=False, header=False)
    prefix.with_suffix(".bed").write_bytes(source.with_suffix(".bed").read_bytes())
    prefix.with_suffix(".fam").write_bytes(source.with_suffix(".fam").read_bytes())
    metadata = metadata.sort_values("POS").reset_index(drop=True)
    args = Namespace(bfile=str(prefix), snp_identifier="rsid", ld_wind_snps=10,
                     ld_wind_kb=None, ld_wind_cm=None, yes_really=True)
    with prepare_plink("1", ChromosomeAnnotations(metadata, pd.DataFrame({"base": metadata.POS}), ["base"], []), args) as prepared:
        assert prepared.metadata["SNP"].tolist() == ["rs_7", "rs_6", "rs_5", "rs_4"]
        np.testing.assert_array_equal(prepared.annotations.read()[:, 0], [1, 2, 3, 4])
        golden = np.load(source.parents[1] / "golden" / "reader_golden.npz")
        np.testing.assert_allclose(prepared.reader.nextSNPs(4), golden["decoded"][:, ::-1], rtol=0, atol=0)


def test_reference_preparation_defers_value_reads_and_kernel_returns_only_output_rows(tmp_path):
    from ldsc._refpanel.prepared import ArrayAnnotations

    metadata = write_panel(tmp_path)
    reads = []
    class RecordingSource(ArrayAnnotations):
        def read(self, *, rows=None, columns=None):
            reads.append(tuple(columns))
            return super().read(rows=rows, columns=columns)
    source = RecordingSource(np.array([[1,2,0],[1,0,4]],dtype=np.float32), ('base','q1','q2'))
    config = GlobalConfig(snp_identifier='rsid')
    panel = ParquetR2RefPanel(config, RefPanelConfig(backend='parquet_r2',r2_dir=tmp_path))
    bundle = ChromosomeAnnotations(metadata,source,['base'],['q1','q2'])
    with panel.prepare_chromosome('1',bundle,LDScoreConfig(ld_wind_snps=1,whole_chromosome_ok=True)) as prepared:
        assert not reads
        result = ldscore.compute_chromosome('1',prepared,snp_identifier='rsid',snp_batch_size=1,
            regression_keys={'rs1'},query_batch_size=1)
        assert result.metadata.SNP.tolist() == ['rs1']
        # Quantized R2 is not exactly representable, so compare at float32 output precision.
        np.testing.assert_allclose(result.ld_scores, [[1 + STORED_R2, 2, 4 * STORED_R2]], rtol=1e-7)
        np.testing.assert_array_equal(result.w_ld, [[1]])
        np.testing.assert_array_equal(result.M, [2,2,4])
        np.testing.assert_array_equal(result.M_5_50, [1,0,4])
        assert result.reference_snp_count == 2
        assert all(len(names)<=1 for names in reads)
def test_prepared_chromosome_close_releases_state_once():
    from types import SimpleNamespace
    from ldsc._refpanel.prepared import PreparedChromosome

    calls = []
    reader = SimpleNamespace(close=lambda: calls.append(1))
    prepared = PreparedChromosome("plink", reader, pd.DataFrame({"POS": [1]}),
                                 np.ones((1, 1)), np.array([0]), ["base"], [],
                                 identity_drops=pd.DataFrame({"SNP": ["duplicate"]}))
    prepared.close()
    prepared.close()
    assert calls == [1]
    assert all(getattr(prepared, name) is None for name in ("reader", "metadata", "annotations", "block_left", "identity_drops"))
