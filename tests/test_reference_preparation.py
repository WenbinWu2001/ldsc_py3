"""Reference policies and numerical projection share the production preparation seam."""
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from ldsc import AnnotationBundle, GlobalConfig, LDScoreCalculator, LDScoreConfig, RefPanelConfig
from ldsc._kernel import ldscore
from ldsc._kernel.ref_panel import ParquetR2RefPanel
from ldsc._kernel.ref_panel_builder import write_runtime_metadata_sidecar
from ldsc._kernel.snp_identity import sidecar_identity_sha256
from ldsc.errors import LDSCInputError, LDSCUsageError


def write_panel(root, *, bias="unbiased", n_samples=200, chrom="1", mode="rsid"):
    root.mkdir(parents=True, exist_ok=True)
    metadata = pd.DataFrame({
        "CHR": [chrom, chrom], "POS": [100, 300], "SNP": ["rs1", "rs3"],
        "A1": ["A", "A"], "A2": ["C", "C"], "CM": [0.1, 0.3], "MAF": [0.04, 0.3],
    })
    write_runtime_metadata_sidecar(metadata, root / f"chr{chrom}_meta.tsv.gz",
                                   genome_build="hg19", snp_identifier=mode)
    footer = {
        b"ldsc:schema_version": b"1", b"ldsc:artifact_type": b"ref_panel_r2",
        b"ldsc:snp_identifier": mode.encode(), b"ldsc:genome_build": b"hg19",
        b"ldsc:sorted_by_build": b"hg19", b"ldsc:n_snps": b"2",
        b"ldsc:sidecar_identity_sha256": sidecar_identity_sha256(metadata).encode(),
    }
    if bias is not None:
        footer[b"ldsc:r2_bias"] = bias.encode()
    if n_samples is not None:
        footer[b"ldsc:n_samples"] = str(n_samples).encode()
    table = pa.table({"IDX_1": pa.array([0], type=pa.int32()), "IDX_2": pa.array([1], type=pa.int32()),
                      "R2": pa.array([0.5], type=pa.float32()), "SIGN": [True]})
    pq.write_table(table.replace_schema_metadata(footer), root / f"chr{chrom}_r2.parquet")
    return metadata


def public_bundle(metadata, values=None):
    values = pd.DataFrame({"base": np.ones(len(metadata))}) if values is None else values
    return AnnotationBundle(metadata=metadata, baseline_annotations=values,
                            query_annotations=pd.DataFrame(index=metadata.index),
                            baseline_columns=list(values.columns), query_columns=[],
                            chromosomes=metadata["CHR"].astype(str).unique().tolist(), source_summary={})


@pytest.mark.parametrize("bias,stored_n,explicit_n,expected", [
    ("raw", 200, None, 1.4974747474747474),
    ("raw", 200, 100, 1.4948979591836735),
    ("unbiased", 200, 100, 1.5),
    (None, 200, None, 1.4974747474747474),
    (None, None, None, 1.5),
])
def test_r2_policy_is_applied_through_public_calculation(tmp_path, bias, stored_n, explicit_n, expected):
    metadata = write_panel(tmp_path, bias=bias, n_samples=stored_n)
    config = GlobalConfig(snp_identifier="rsid")
    panel = ParquetR2RefPanel(config, RefPanelConfig(backend="parquet_r2", r2_dir=tmp_path,
                                                   sample_size=explicit_n))
    result = LDScoreCalculator().compute_chromosome(
        "1", public_bundle(metadata), panel, LDScoreConfig(ld_wind_snps=1, whole_chromosome_ok=True), config,
    )
    np.testing.assert_allclose(result.baseline_table[["base", "regression_ld_scores"]], expected, rtol=1e-7)


@pytest.mark.parametrize("alleles_in_annotation", [False, True])
def test_prepared_annotations_share_reference_rows_alleles_cm_maf_and_counts(tmp_path, alleles_in_annotation):
    mode = "chr_pos_allele_aware"
    reference = write_panel(tmp_path, mode=mode)
    metadata = pd.DataFrame({"CHR": ["1"] * 3, "SNP": ["rs1", "missing", "rs3"],
                             "POS": [100, 200, 300], "CM": [9.0] * 3, "MAF": [0.9] * 3})
    if alleles_in_annotation:
        metadata = metadata.assign(A1="A", A2="C")
    values = pd.DataFrame({"base": [1.0, 99.0, 3.0], "flag": [True, True, False]})
    config = GlobalConfig(snp_identifier=mode, genome_build="hg19")
    panel = ParquetR2RefPanel(config, RefPanelConfig(backend="parquet_r2", r2_dir=tmp_path))
    bundle = ldscore.AnnotationBundle(metadata, values, ["base", "flag"], [])
    with panel.prepare_chromosome("1", bundle, LDScoreConfig(ld_wind_cm=0.3, whole_chromosome_ok=True)) as prepared:
        assert prepared.annotations.read().dtype == np.float32
        np.testing.assert_array_equal(prepared.annotations.read(), [[1, 1], [3, 0]])
        pd.testing.assert_frame_equal(prepared.metadata[reference.columns], reference)
        result = ldscore.compute_chromosome("1", prepared, snp_identifier=mode, snp_batch_size=3)
        np.testing.assert_allclose(result.ld_scores, [[2.5, 1], [3.5, 0.5]])
        np.testing.assert_array_equal(result.M, [4, 1])
        np.testing.assert_array_equal(result.M_5_50, [3, 0])
        reader = prepared.reader
    assert reader._pf.closed
    np.testing.assert_array_equal(metadata["CM"], [9.0] * 3)


def test_allele_mismatch_is_rejected_before_reader_creation(tmp_path):
    metadata = write_panel(tmp_path, mode="chr_pos_allele_aware").assign(A2="G")
    config = GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg19")
    panel = ParquetR2RefPanel(config, RefPanelConfig(backend="parquet_r2", r2_dir=tmp_path))
    with pytest.raises(LDSCInputError, match="retained no annotation SNPs"):
        LDScoreCalculator().compute_chromosome(
            "1", public_bundle(metadata), panel, LDScoreConfig(ld_wind_snps=1), config,
        )


def test_parquet_preparation_releases_reader_on_calculation_failure(tmp_path):
    metadata = write_panel(tmp_path, bias="raw", n_samples=None)
    config = GlobalConfig(snp_identifier="rsid")
    panel = ParquetR2RefPanel(config, RefPanelConfig(backend="parquet_r2", r2_dir=tmp_path))
    bundle = ldscore.AnnotationBundle(metadata, pd.DataFrame({"base": [1, 1]}), ["base"], [])
    with pytest.raises(LDSCUsageError, match="sample size"):
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


def test_public_run_skips_empty_intersection_and_computes_other_chromosome(tmp_path):
    absent = write_panel(tmp_path, chrom="1").assign(SNP=["absent1", "absent3"])
    present = write_panel(tmp_path, chrom="22")
    config = GlobalConfig(snp_identifier="rsid")
    panel = ParquetR2RefPanel(config, RefPanelConfig(backend="parquet_r2", r2_dir=tmp_path))
    bundle = public_bundle(pd.concat([absent, present], ignore_index=True))
    with pytest.warns(UserWarning, match="Skipping chromosome 1"):
        result = LDScoreCalculator().run(
            bundle, panel, LDScoreConfig(ld_wind_snps=1, whole_chromosome_ok=True), config,
        )
    assert result.baseline_table["CHR"].tolist() == ["22", "22"]
    np.testing.assert_allclose(result.baseline_table["base"], [1.5, 1.5])


def test_parquet_genetic_maps_are_ignored_at_preparation(tmp_path, caplog):
    metadata = write_panel(tmp_path)
    panel = ParquetR2RefPanel(GlobalConfig(snp_identifier="rsid"), RefPanelConfig(
        backend="parquet_r2", r2_dir=tmp_path, genetic_map_hg19_sources="unused-map.tsv",
    ))
    with panel.prepare_chromosome(
        "1", ldscore.AnnotationBundle(metadata, pd.DataFrame({"base": [1, 1]}), ["base"], []),
        LDScoreConfig(ld_wind_cm=0.3, whole_chromosome_ok=True),
    ) as prepared:
        np.testing.assert_allclose(prepared.metadata["CM"], [0.1, 0.3])
    assert "Ignoring --genetic-map" in caplog.text


def test_plink_preparation_maps_sorted_annotations_to_physical_bed_columns(tmp_path):
    pytest.importorskip("bitarray")
    from tests.ref_panel_helpers import prepare_plink
    from argparse import Namespace
    source = Path(__file__).parent / "fixtures" / "plink" / "plink"
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
    with prepare_plink("1", ldscore.AnnotationBundle(metadata, pd.DataFrame({"base": metadata.POS}), ["base"], []), args) as prepared:
        assert prepared.metadata["SNP"].tolist() == ["rs_7", "rs_6", "rs_5", "rs_4"]
        np.testing.assert_array_equal(prepared.annotations.read()[:, 0], [1, 2, 3, 4])
        golden = np.load(source.parents[1] / "golden" / "reader_golden.npz")
        np.testing.assert_allclose(prepared.reader.nextSNPs(4), golden["decoded"][:, ::-1], rtol=0, atol=0)


def test_reference_preparation_defers_value_reads_and_kernel_returns_only_output_rows(tmp_path):
    from ldsc._kernel.ldscore_projection import ArrayAnnotations

    metadata = write_panel(tmp_path)
    reads = []
    class RecordingSource(ArrayAnnotations):
        def read(self, *, rows=None, columns=None):
            reads.append(tuple(columns))
            return super().read(rows=rows, columns=columns)
    source = RecordingSource(np.array([[1,2,0],[1,0,4]],dtype=np.float32), ('base','q1','q2'))
    config = GlobalConfig(snp_identifier='rsid')
    panel = ParquetR2RefPanel(config, RefPanelConfig(backend='parquet_r2',r2_dir=tmp_path))
    bundle = ldscore.AnnotationBundle(metadata,source,['base'],['q1','q2'])
    with panel.prepare_chromosome('1',bundle,LDScoreConfig(ld_wind_snps=1,whole_chromosome_ok=True)) as prepared:
        assert not reads
        result = ldscore.compute_chromosome('1',prepared,snp_identifier='rsid',snp_batch_size=1,
            regression_keys={'rs1'},query_batch_size=1)
        assert result.metadata.SNP.tolist() == ['rs1']
        np.testing.assert_array_equal(result.ld_scores, [[1.5,2,2]])
        np.testing.assert_array_equal(result.w_ld, [[1]])
        np.testing.assert_array_equal(result.M, [2,2,4])
        np.testing.assert_array_equal(result.M_5_50, [1,0,4])
        assert result.reference_snp_count == 2
        assert all(len(names)<=1 for names in reads)
