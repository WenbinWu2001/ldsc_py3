from __future__ import annotations

from argparse import Namespace
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from ldsc import cli, gene_ldscore_index
from ldsc._kernel import ldscore as kernel_ldscore
from ldsc._kernel.gene_ldscore_index import AtomStatistics, build_disjoint_atoms
from ldsc.errors import LDSCInputError
from ldsc.gene_ldscore_index import (
    _build_embedded_gene_catalog,
    _genetic_map_identity,
    _parse_chromosomes,
    _selected_individual_ids,
    build_plink_index_chromosome,
    calculate_index_id,
    IndexChromosomeData,
    load_gene_ldscore_index,
    publish_gene_ldscore_index,
    run_indexed_ldscore,
    intersect_baseline_plink_by_identifier,
)
from ldsc.gene_list_resolver import GeneCatalog


def _rows(*items: tuple[str, int, str]) -> pd.DataFrame:
    return pd.DataFrame(items, columns=["CHR", "POS", "SNP"])


def _artifact_payload():
    atom_model = build_disjoint_atoms("22", np.array([[0, 10], [5, 15]], dtype=np.int64))
    chromosome = IndexChromosomeData(
        baseline_rows=pd.DataFrame(
            {
                "CHR": ["22", "22"], "SNP": ["rs1", "rs2"], "POS": [1, 6],
                "regression_ld_scores": [1.0, 2.0], "base": [3.0, 4.0],
            }
        ),
        baseline_count_all=np.array([2.0]), baseline_count_common=np.array([1.0]),
        baseline_overlap_all=np.array([[2.0]]), baseline_overlap_common=np.array([[1.0]]),
        total_reference_snps_all=2, total_reference_snps_common=1,
        atom_model=atom_model,
        operator=sparse.csr_matrix(np.array([[1.0, -0.25, 0.0], [0.5, 1.0, 0.1]])),
        atom_statistics=AtomStatistics(
            np.array([1, 1, 0], dtype=np.int64), np.array([1, 0, 0], dtype=np.int64),
            np.array([[1.0, 1.0, 0.0]]), np.array([[1.0, 0.0, 0.0]]),
        ),
        reference_metadata=pd.DataFrame(),
    )
    catalog = pd.DataFrame(
        {
            "gene_index": [0, 1], "canonical_ensembl_id": ["ENSG1", "ENSG2"],
            "gene_name": ["G1", "G2"], "CHR": ["22", "22"], "start0": [0, 5],
            "end": [10, 15], "included": [True, True], "exclusion_reason": [None, None],
            "chromosome_gene_row": [0, 1],
        }
    )
    index_identity = {
        "chromosomes": ["22"], "baseline": "fixture", "plink": "fixture",
        "catalog": {"release": "GENCODE v49", "content_sha256": "fixture"},
        "projection_build": "hg19", "padding_bp": 100000, "gene_exclude_regions": "mhc",
    }
    return chromosome, catalog, index_identity


def test_builder_chromosome_selection_is_canonical_and_unique():
    assert _parse_chromosomes("22,1-3,2") == ("1", "2", "3", "22")
    with pytest.raises(LDSCInputError, match="1 through 22"):
        _parse_chromosomes("23")


def test_selected_individual_identity_uses_fam_row_order():
    prefix = Path(__file__).resolve().parent / "fixtures" / "plink" / "plink"
    fam = kernel_ldscore.legacy_parse.PlinkFAMFile(str(prefix.with_suffix(".fam")))

    assert _selected_individual_ids(fam, [1, 0]) == [
        str(fam.IDList.iloc[1, 0]),
        str(fam.IDList.iloc[0, 0]),
    ]


def test_embedded_catalog_keeps_mhc_exclusions_with_empty_gene_rows():
    frame = pd.DataFrame(
        {
            "ensgid": ["ENSG1", "ENSG2", "ENSG3"],
            "gene_name": ["MHC_GENE", "OTHER6", "CHR22"],
            "hg19_chr": ["6", "6", "22"],
            "hg19_start0": [26_000_000, 40_000_000, 100],
            "hg19_end": [26_001_000, 40_001_000, 200],
            "hg38_chr": [pd.NA, pd.NA, pd.NA],
            "hg38_start0": [pd.NA, pd.NA, pd.NA],
            "hg38_end": [pd.NA, pd.NA, pd.NA],
        }
    )
    catalog = GeneCatalog(
        frame=frame,
        resource="test.tsv.gz",
        release="test",
        content_sha256="abc",
        ensembl_to_index={"ENSG1": 0, "ENSG2": 1, "ENSG3": 2},
        gene_name_to_indices={"MHC_GENE": (0,), "OTHER6": (1,), "CHR22": (2,)},
    )

    embedded = _build_embedded_gene_catalog(
        catalog, genome_build="hg19", gene_exclude_regions="mhc"
    )

    assert embedded["included"].tolist() == [False, True, True]
    assert embedded["exclusion_reason"].tolist() == ["excluded_gene_region", "", ""]
    assert embedded["chromosome_gene_row"].tolist() == [0, 1, 0]


def test_baseline_plink_intersection_uses_identifier_mode_and_plink_order(caplog):
    baseline_metadata = _rows(
        ("22", 999, "rs3"), ("22", 10, "rs1"), ("22", 20, "baseline_only")
    )
    annotations = pd.DataFrame({"base": [3.0, 1.0, 2.0]})
    plink = _rows(("22", 10, "rs1"), ("22", 30, "rs3"), ("22", 40, "plink_only"))

    result = intersect_baseline_plink_by_identifier(
        baseline_metadata,
        annotations,
        plink,
        snp_identifier="rsid",
        chrom="22",
    )

    assert result.metadata[["CHR", "POS", "SNP"]].to_records(index=False).tolist() == [
        ("22", 10, "rs1"),
        ("22", 30, "rs3"),
    ]
    assert result.annotations["base"].tolist() == [1.0, 3.0]
    assert result.diagnostics == {
        "baseline_rows": 3,
        "plink_rows": 3,
        "matched_rows": 2,
        "baseline_only_rows": 1,
        "plink_only_rows": 1,
        "coordinate_discordant_rows": 1,
    }
    assert "coordinate disagreement" in caplog.text


@pytest.mark.parametrize("side", ["baseline", "PLINK BIM"])
def test_baseline_plink_intersection_rejects_duplicate_effective_identifiers(side):
    baseline = _rows(("22", 10, "rs1"), ("22", 20, "rs2"))
    plink = _rows(("22", 10, "rs1"), ("22", 20, "rs2"))
    if side == "baseline":
        baseline = pd.concat([baseline, baseline.iloc[[0]]], ignore_index=True)
    else:
        plink = pd.concat([plink, plink.iloc[[0]]], ignore_index=True)
    with pytest.raises(LDSCInputError, match=f"duplicate {side}.*rsid"):
        intersect_baseline_plink_by_identifier(
            baseline,
            pd.DataFrame({"base": np.ones(len(baseline))}),
            plink,
            snp_identifier="rsid",
            chrom="22",
        )


def test_baseline_plink_intersection_rejects_empty_result():
    with pytest.raises(LDSCInputError, match="intersection is empty"):
        intersect_baseline_plink_by_identifier(
            _rows(("22", 10, "rs1")),
            pd.DataFrame({"base": [1.0]}),
            _rows(("22", 20, "rs2")),
            snp_identifier="rsid",
            chrom="22",
        )


def test_prepared_plink_chromosome_is_reused_by_direct_ldscore():
    prefix = Path(__file__).resolve().parent / "fixtures" / "plink" / "plink"
    bim = pd.read_csv(
        prefix.with_suffix(".bim"),
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "POS", "A1", "A2"],
    )
    metadata = bim[["CHR", "SNP", "CM", "POS"]].copy()
    annotations = pd.DataFrame({"base": metadata["POS"].astype(np.float64)})
    bundle = kernel_ldscore.AnnotationBundle(metadata, annotations, ["base"], [])
    args = Namespace(
        bfile=str(prefix),
        keep=None,
        maf_min=None,
        maf=None,
        ld_wind_snps=10,
        ld_wind_kb=None,
        ld_wind_cm=None,
        yes_really=True,
        snp_batch_size=3,
        common_maf_min=0.05,
        snp_identifier="rsid",
        genetic_map=None,
    )

    prepared = kernel_ldscore.prepare_plink_chromosome("1", bundle, args)
    direct = kernel_ldscore.compute_chrom_from_plink("1", bundle, args, None)

    assert prepared.annotation_matrix.shape == (len(prepared.metadata), 1)
    np.testing.assert_array_equal(
        prepared.annotation_matrix[:, 0],
        prepared.metadata["POS"].to_numpy(dtype=np.float64),
    )
    prepared.geno._currentSNP = 0
    reused_scores = prepared.geno.ldScoreVarBlocks(
        prepared.block_left,
        args.snp_batch_size,
        annot=prepared.annotation_matrix,
    )
    np.testing.assert_array_equal(reused_scores.astype(np.float32), direct.ld_scores)


def test_plink_atom_operator_matches_direct_union_in_atom_batches():
    prefix = Path(__file__).resolve().parent / "fixtures" / "plink" / "plink"
    bim = pd.read_csv(
        prefix.with_suffix(".bim"),
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "POS", "A1", "A2"],
    )
    metadata = bim[["CHR", "SNP", "CM", "POS"]].copy()
    baseline = pd.DataFrame({"base": np.ones(len(metadata), dtype=np.float64)})
    base_bundle = kernel_ldscore.AnnotationBundle(metadata, baseline, ["base"], [])
    args = Namespace(
        bfile=str(prefix), keep=None, maf_min=None, maf=None,
        ld_wind_snps=10, ld_wind_kb=None, ld_wind_cm=None,
        yes_really=True, snp_batch_size=3, common_maf_min=0.05,
        snp_identifier="rsid", genetic_map=None,
    )
    regression_keys = set(metadata["SNP"].astype(str))
    gene_intervals = np.array([[0, 4], [2, 7], [3, 5]], dtype=np.int64)

    indexed_one = build_plink_index_chromosome(
        "1",
        base_bundle,
        args,
        regression_keys=regression_keys,
        regression_regions=None,
        gene_intervals=gene_intervals,
        included=np.array([True, True, True]),
        padding_bp=0,
        atom_batch_size=1,
    )
    indexed_many = build_plink_index_chromosome(
        "1",
        base_bundle,
        args,
        regression_keys=regression_keys,
        regression_regions=None,
        gene_intervals=gene_intervals,
        included=np.array([True, True, True]),
        padding_bp=0,
        atom_batch_size=3,
    )
    z = np.asarray(indexed_one.atom_model.gene_to_atom[[0, 1]].max(axis=0).toarray()).reshape(-1).astype(bool)
    direct_annotation = np.array(
        [any(start <= pos - 1 < end for start, end in gene_intervals[:2]) for pos in metadata["POS"]],
        dtype=np.float64,
    )
    direct_bundle = kernel_ldscore.AnnotationBundle(
        metadata,
        pd.DataFrame({"base": baseline["base"], "query": direct_annotation}),
        ["base"],
        ["query"],
    )
    direct = kernel_ldscore.compute_chrom_from_plink("1", direct_bundle, args, regression_keys)
    persisted = kernel_ldscore.regression_mask_from_keys(
        direct.metadata, regression_keys, "rsid"
    ).astype(bool)

    np.testing.assert_array_equal(indexed_one.operator.toarray(), indexed_many.operator.toarray())
    np.testing.assert_array_equal(
        np.asarray(indexed_one.operator @ z.astype(np.float64)).astype(np.float32),
        direct.ld_scores[persisted, 1],
    )
    assert int(indexed_one.atom_statistics.atom_count_all @ z.astype(np.int64)) == int(direct.M[1])


def test_build_index_command_registers_closed_v1_configuration():
    parser = cli.build_parser()
    args = parser.parse_args(
        [
            "build-gene-ldscore-index",
            "--baseline-annot-sources", "baseline.@.annot.gz",
            "--plink-prefix", "panel.@",
            "--output-dir", "suite",
        ]
    )

    assert args.command == "build-gene-ldscore-index"
    assert args.genome_build == "hg19"
    assert args.snp_identifier == "rsid"
    assert args.padding_bp == 100000
    assert args.gene_exclude_regions == "mhc"
    assert args.ld_wind_cm == 1.0
    assert args.common_maf_min == 0.05
    assert args.threads == 1
    assert args.atom_batch_size > 0
    assert args.regression_snps_file is None
    assert args.exclude_regions == "mhc-and-centromeres"


def test_build_index_command_accepts_custom_regression_snps_and_region_policy():
    args = cli.build_parser().parse_args(
        [
            "build-gene-ldscore-index",
            "--baseline-annot-sources", "baseline.@.annot.gz",
            "--plink-prefix", "panel.@",
            "--output-dir", "index",
            "--regression-snps-file", "custom.snplist",
            "--exclude-regions", "centromeres",
        ]
    )
    assert args.regression_snps_file == "custom.snplist"
    assert args.exclude_regions == "centromeres"


def test_build_index_writes_shared_operational_log_and_chromosome_metrics(tmp_path, monkeypatch):
    chromosome, catalog, index_identity = _artifact_payload()
    identity_rows = chromosome.baseline_rows[["CHR", "POS", "SNP"]].copy()
    public_bundle = SimpleNamespace(
        metadata=identity_rows.assign(CM=[0.1, 0.2]),
        baseline_annotations=pd.DataFrame({"base": [1.0, 1.0]}),
        baseline_columns=["base"],
    )
    fake_catalog = SimpleNamespace(resource="catalog.tsv.gz", release="test", content_sha256="catalog")
    embedded_catalog = catalog.copy()
    index_identity = {
        **index_identity,
        "plink_sources": [
            {"chromosome": "22", "kind": kind, "content_sha256": kind}
            for kind in ("bed", "bim")
        ],
        "chromosomes": ["22"],
        "genome_build": "hg19",
        "snp_identifier": "rsid",
        "selected_individuals": {"source": "all", "selected_count": 3},
        "maf_min": None,
        "common_maf_min": 0.05,
        "ld_window": {"unit": "cm", "value": 1.0},
        "genetic_map": "bim_cm",
        "regression_snps": {"kind": "bundled_hapmap3"},
        "exclude_regions": "mhc-and-centromeres",
    }

    class FakeAnnotationBuilder:
        def __init__(self, _global_config):
            pass

        def run(self, _annotation_spec, *, chrom):
            assert chrom == "22"
            return public_bundle

    monkeypatch.setattr(gene_ldscore_index, "AnnotationBuilder", FakeAnnotationBuilder)
    monkeypatch.setattr(gene_ldscore_index.GeneCatalog, "load", lambda: fake_catalog)
    monkeypatch.setattr(gene_ldscore_index, "_build_embedded_gene_catalog", lambda *args, **kwargs: embedded_catalog)
    monkeypatch.setattr(gene_ldscore_index, "_load_builder_genetic_map", lambda _args: None)
    captured_restriction = {}

    def fake_read_restriction(path, *args, **kwargs):
        captured_restriction["path"] = Path(path)
        return {"rs1"}

    monkeypatch.setattr(gene_ldscore_index, "read_snp_restriction_keys", fake_read_restriction)
    monkeypatch.setattr(gene_ldscore_index.kernel_regions, "load_preset_intervals", lambda *args: None)
    monkeypatch.setattr(gene_ldscore_index, "_read_bim_identity", lambda _prefix: identity_rows.copy())
    monkeypatch.setattr(gene_ldscore_index, "build_plink_index_chromosome", lambda *args, **kwargs: chromosome)
    monkeypatch.setattr(gene_ldscore_index, "_builder_index_identity", lambda *args, **kwargs: index_identity)
    monkeypatch.setattr(gene_ldscore_index.kernel_ldscore, "resolve_bfile_prefix", lambda *args, **kwargs: "fixture")

    args = Namespace(
        baseline_annot_sources="baseline.22.annot.gz",
        plink_prefix="reference/1000G.EUR.QC.22",
        output_dir=str(tmp_path / "suite"),
        genome_build="hg19",
        snp_identifier="rsid",
        padding_bp=100000,
        gene_exclude_regions="mhc",
        ld_wind_cm=1.0,
        maf_min=None,
        common_maf_min=0.05,
        keep_indivs_file=None,
        regression_snps_file="custom-regression.tsv",
        exclude_regions="centromeres",
        genetic_map_hg19_sources=None,
        genetic_map_hg38_sources=None,
        chromosomes="22",
        snp_batch_size=128,
        atom_batch_size=64,
        threads=1,
        overwrite=False,
        log_level="INFO",
    )

    index_dir = gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)

    log_path = index_dir.with_name(f"{index_dir.name}.build") / "build-gene-ldscore-index.log"
    json_path = index_dir / "diagnostics" / "build-gene-ldscore-index.json"
    log_text = log_path.read_text(encoding="utf-8")
    diagnostics = json.loads(json_path.read_text(encoding="utf-8"))

    assert "LDSC build-gene-ldscore-index Started" in log_text
    assert "Inputs:" in log_text
    assert "genome_build" in log_text
    assert "custom" in log_text
    assert "exclude-regions=centromeres" in log_text
    assert "baseline/PLINK identifier intersection" in log_text
    assert "Starting chromosome 22" in log_text
    assert "Finished chromosome 22" in log_text
    assert "protein-coding genes=2" in log_text
    assert "Finished " in log_text
    assert diagnostics["chromosomes"]["22"]["protein_coding_genes"] == 2
    assert diagnostics["chromosomes"]["22"]["operator_nnz"] == chromosome.operator.nnz
    assert f"operator_nnz={chromosome.operator.nnz}" in log_text
    assert captured_restriction["path"] == Path("custom-regression.tsv")
    assert not (index_dir / "diagnostics" / "build-gene-ldscore-index.log").exists()


def test_build_index_failure_keeps_stable_log_without_publishing_index(tmp_path, monkeypatch):
    chromosome, catalog, index_identity = _artifact_payload()
    identity_rows = chromosome.baseline_rows[["CHR", "POS", "SNP"]].copy()
    public_bundle = SimpleNamespace(
        metadata=identity_rows.assign(CM=[0.1, 0.2]),
        baseline_annotations=pd.DataFrame({"base": [1.0, 1.0]}),
        baseline_columns=["base"],
    )
    fake_catalog = SimpleNamespace(resource="catalog.tsv.gz", release="test", content_sha256="catalog")
    embedded_catalog = catalog.copy()

    class FakeAnnotationBuilder:
        def __init__(self, _global_config):
            pass

        def run(self, _annotation_spec, *, chrom):
            return public_bundle

    monkeypatch.setattr(gene_ldscore_index, "AnnotationBuilder", FakeAnnotationBuilder)
    monkeypatch.setattr(gene_ldscore_index.GeneCatalog, "load", lambda: fake_catalog)
    monkeypatch.setattr(gene_ldscore_index, "_build_embedded_gene_catalog", lambda *args, **kwargs: embedded_catalog)
    monkeypatch.setattr(gene_ldscore_index, "_builder_index_identity", lambda *args, **kwargs: index_identity)
    monkeypatch.setattr(gene_ldscore_index, "_load_builder_genetic_map", lambda _args: None)
    monkeypatch.setattr(gene_ldscore_index, "read_snp_restriction_keys", lambda *args, **kwargs: {"rs1"})
    monkeypatch.setattr(gene_ldscore_index.kernel_regions, "load_preset_intervals", lambda *args: None)
    monkeypatch.setattr(gene_ldscore_index, "_read_bim_identity", lambda _prefix: identity_rows.copy())
    monkeypatch.setattr(gene_ldscore_index.kernel_ldscore, "resolve_bfile_prefix", lambda *args, **kwargs: "fixture")

    original_intersection = gene_ldscore_index.intersect_baseline_plink_by_identifier

    def fail_intersection(*args, **kwargs):
        raise LDSCInputError("identifier intersection failure in fixture")

    monkeypatch.setattr(gene_ldscore_index, "intersect_baseline_plink_by_identifier", fail_intersection)
    args = Namespace(
        baseline_annot_sources="baseline.22.annot.gz",
        plink_prefix="reference/1000G.EUR.QC.22",
        output_dir=str(tmp_path / "suite"),
        genome_build="hg19",
        snp_identifier="rsid",
        padding_bp=100000,
        gene_exclude_regions="mhc",
        ld_wind_cm=1.0,
        maf_min=None,
        common_maf_min=0.05,
        keep_indivs_file=None,
        regression_snps_file=None,
        exclude_regions="mhc-and-centromeres",
        genetic_map_hg19_sources=None,
        genetic_map_hg38_sources=None,
        chromosomes="22",
        snp_batch_size=128,
        atom_batch_size=64,
        threads=1,
        overwrite=False,
        log_level="INFO",
    )

    with pytest.raises(LDSCInputError, match="identifier intersection failure"):
        gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)

    index_path = tmp_path / "suite"
    log_path = tmp_path / "suite.build" / "build-gene-ldscore-index.log"
    log_text = log_path.read_text(encoding="utf-8")
    assert "LDSC build-gene-ldscore-index Started" in log_text
    assert "Chromosome 22 failed during baseline/PLINK identifier intersection" in log_text
    assert "LDSCInputError" in log_text
    assert "Failed " in log_text
    assert not index_path.exists()

    index_path.mkdir()
    with pytest.raises(LDSCInputError, match="identifier intersection failure"):
        gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)
    assert list(index_path.iterdir()) == []

    monkeypatch.setattr(gene_ldscore_index, "intersect_baseline_plink_by_identifier", original_intersection)
    monkeypatch.setattr(gene_ldscore_index, "build_plink_index_chromosome", lambda *args, **kwargs: chromosome)
    index_dir = gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)
    assert (index_dir / "metadata.json").exists()
    assert any((tmp_path / "suite.build" / "history").glob("*.log"))


def test_plink_operator_matches_independent_dense_adjusted_r2_reference():
    prefix = Path(__file__).resolve().parent / "fixtures" / "plink" / "plink"
    bim = pd.read_csv(
        prefix.with_suffix(".bim"), sep=r"\s+", header=None,
        names=["CHR", "SNP", "CM", "POS", "A1", "A2"],
    )
    bundle = kernel_ldscore.AnnotationBundle(
        bim[["CHR", "SNP", "CM", "POS"]],
        pd.DataFrame({"base": np.ones(len(bim))}),
        ["base"],
        [],
    )
    args = Namespace(
        bfile=str(prefix), keep=None, maf_min=None, maf=None,
        ld_wind_snps=1, ld_wind_kb=None, ld_wind_cm=None,
        yes_really=True, snp_batch_size=1, common_maf_min=0.05,
        snp_identifier="rsid", genetic_map=None,
    )
    persisted_keys = {"rs_5", "rs_6"}
    indexed = build_plink_index_chromosome(
        "1", bundle, args,
        regression_keys=persisted_keys,
        regression_regions=None,
        gene_intervals=np.array([[4, 5], [5, 6], [6, 7], [7, 8]], dtype=np.int64),
        included=np.ones(4, dtype=bool),
        padding_bp=0,
        atom_batch_size=2,
    )
    prepared = kernel_ldscore.prepare_plink_chromosome("1", bundle, args)
    prepared.geno._currentSNP = 0
    x = prepared.geno.nextSNPs(prepared.geno.m)
    correlation = x.T @ x / prepared.geno.n
    adjusted_r2 = correlation**2 - (1.0 - correlation**2) / (prepared.geno.n - 2)
    within_window = np.abs(np.subtract.outer(np.arange(prepared.geno.m), np.arange(prepared.geno.m))) <= 1
    dense_r = np.where(within_window, adjusted_r2, 0.0)
    persisted = prepared.metadata["SNP"].isin(persisted_keys).to_numpy()

    np.testing.assert_allclose(indexed.operator.toarray(), dense_r[persisted], rtol=0, atol=1e-15)
    np.testing.assert_allclose(np.diag(dense_r), np.ones(prepared.geno.m), rtol=0, atol=1e-15)
    assert indexed.operator.toarray()[0, 2] < 0
    assert indexed.operator.toarray()[0, 0] == pytest.approx(1.0)  # non-regression rs_4 still contributes
    assert indexed.operator.toarray()[0, 3] == 0.0  # outside the one-SNP window


def test_semantic_id_uses_canonical_scientific_identity_only():
    index_identity = {
        "baseline_sources": [{"ordinal": 1, "content_sha256": "abc"}],
        "plink": {"bed_sha256": "bed", "bim_sha256": "bim", "fam_sha256": "fam"},
        "chromosomes": ["22"],
        "ld_window": {"mode": "cm", "value": 1.0},
    }
    index_id = calculate_index_id(index_identity)
    assert index_id == calculate_index_id(dict(reversed(list(index_identity.items()))))
    assert index_id != calculate_index_id({**index_identity, "chromosomes": ["21", "22"]})


def test_explicit_genetic_map_identity_is_content_bound(tmp_path):
    first = tmp_path / "map1.txt"
    second = tmp_path / "map2.txt"
    first.write_text("CHR POS CM\n22 1 0.1\n", encoding="utf-8")
    second.write_text("CHR POS CM\n22 2 0.2\n", encoding="utf-8")

    bim_identity = _genetic_map_identity(Namespace(genetic_map_hg19_sources=None))
    explicit = _genetic_map_identity(
        Namespace(genetic_map_hg19_sources=f"{first},{second}")
    )

    assert bim_identity == "bim_cm"
    assert explicit["kind"] == "explicit_hg19"
    assert len(explicit["content_sha256"]) == 64


def test_canonical_table_identity_ignores_row_order_when_key_sort_is_declared():
    frame = pd.DataFrame(
        {"CHR": ["22", "22"], "POS": [2, 1], "SNP": ["rs2", "rs1"], "base": [0.0, 1.0]}
    )
    reversed_frame = frame.iloc[::-1].reset_index(drop=True)

    assert gene_ldscore_index._canonical_frame_sha256(
        frame, sort_by=("CHR", "POS", "SNP")
    ) == gene_ldscore_index._canonical_frame_sha256(
        reversed_frame, sort_by=("CHR", "POS", "SNP")
    )


def test_index_artifact_round_trip_uses_approved_tree_and_payloads(tmp_path):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = publish_gene_ldscore_index(
        tmp_path / "index",
        index_identity=index_identity,
        gene_catalog=catalog,
        chromosomes={"22": chromosome},
        overwrite=False,
    )
    loaded = load_gene_ldscore_index(index_dir)

    assert index_dir == tmp_path / "index"
    assert loaded.index_id == calculate_index_id(index_identity)
    assert loaded.chromosomes == ("22",)
    np.testing.assert_array_equal(loaded.index_chromosomes["22"].operator.toarray(), chromosome.operator.toarray())
    root_metadata = (index_dir / "metadata.json").read_text(encoding="utf-8")
    assert "schema_version" not in root_metadata
    assert "software_version" not in root_metadata
    assert (index_dir / "chromosomes" / "chr22" / "baseline_statistics.npz").exists()
    assert (index_dir / "chromosomes" / "chr22" / "atom_statistics.npz").exists()
    stored_rows = pd.read_parquet(index_dir / "chromosomes" / "chr22" / "baseline_rows.parquet")
    assert stored_rows["regression_ld_scores"].dtype == np.float32
    assert stored_rows["base"].dtype == np.float32


def test_index_loader_rejects_semantic_id_and_sparse_dtype_corruption(tmp_path):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = publish_gene_ldscore_index(
        tmp_path / "index", index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    metadata_path = index_dir / "metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["index_id"] = "corrupt"
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(LDSCInputError, match="semantic identity"):
        load_gene_ldscore_index(index_dir)

    metadata["index_id"] = calculate_index_id(index_identity)
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    operator_path = index_dir / "chromosomes" / "chr22" / "ldscore_operator.npz"
    sparse.save_npz(operator_path, chromosome.operator.astype(np.float32))
    with pytest.raises(Exception, match="float64"):
        load_gene_ldscore_index(index_dir)


def test_index_loader_rejects_noncanonical_catalog_and_baseline_row_order(tmp_path):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = publish_gene_ldscore_index(
        tmp_path / "index", index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    catalog_path = index_dir / "gene_catalog.parquet"
    reversed_catalog = pd.read_parquet(catalog_path).iloc[::-1].reset_index(drop=True)
    reversed_catalog.to_parquet(catalog_path, index=False)
    with pytest.raises(LDSCInputError, match="gene_index ordering"):
        load_gene_ldscore_index(index_dir)

    catalog.to_parquet(catalog_path, index=False)
    baseline_path = index_dir / "chromosomes" / "chr22" / "baseline_rows.parquet"
    reversed_rows = pd.read_parquet(baseline_path).iloc[::-1].reset_index(drop=True)
    reversed_rows.to_parquet(baseline_path, index=False)
    with pytest.raises(LDSCInputError, match="canonical genomic order"):
        load_gene_ldscore_index(index_dir)


def test_index_loader_rejects_chromosome_coverage_corruption(tmp_path):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = publish_gene_ldscore_index(
        tmp_path / "index", index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    metadata_path = index_dir / "metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["chromosomes"] = []
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

    with pytest.raises(LDSCInputError, match="chromosome coverage"):
        load_gene_ldscore_index(index_dir)


def test_index_publication_requires_overwrite_and_replaces_complete_index(tmp_path):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = tmp_path / "index"
    publish_gene_ldscore_index(
        index_dir, index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    with pytest.raises(FileExistsError, match="--overwrite"):
        publish_gene_ldscore_index(
            index_dir, index_identity=index_identity, gene_catalog=catalog,
            chromosomes={"22": chromosome}, overwrite=False,
        )
    replacement_identity = {**index_identity, "padding_bp": 0}
    publish_gene_ldscore_index(
        index_dir, index_identity=replacement_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=True,
    )
    assert load_gene_ldscore_index(index_dir).index_id == calculate_index_id(replacement_identity)


def test_post_commit_cleanup_failure_warns_but_published_index_succeeds(
    tmp_path, monkeypatch, caplog
):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = tmp_path / "index"
    publish_gene_ldscore_index(
        index_dir, index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    replacement_identity = {**index_identity, "padding_bp": 0}
    real_rmtree = shutil.rmtree

    def fail_transaction_cleanup(path, *args, **kwargs):
        candidate = Path(path)
        if candidate.parent == tmp_path and candidate.name.startswith(".index.stage-"):
            raise OSError(16, "injected device or resource busy", str(candidate))
        return real_rmtree(path, *args, **kwargs)

    with monkeypatch.context() as patcher, caplog.at_level(
        "WARNING", logger="LDSC.gene_ldscore_index"
    ):
        patcher.setattr(gene_ldscore_index.shutil, "rmtree", fail_transaction_cleanup)
        published = publish_gene_ldscore_index(
            index_dir, index_identity=replacement_identity, gene_catalog=catalog,
            chromosomes={"22": chromosome}, overwrite=True,
        )

    assert published == index_dir
    assert load_gene_ldscore_index(index_dir).index_id == calculate_index_id(replacement_identity)
    leftovers = list(tmp_path.glob(".index.stage-*"))
    assert len(leftovers) == 1
    assert str(leftovers[0]) in caplog.text
    assert "published index is valid" in caplog.text

    gene_ldscore_index._recover_gene_index_publication(index_dir)
    assert not leftovers[0].exists()


def test_index_preflight_rejects_nonempty_invalid_directory_even_with_overwrite(tmp_path):
    index_dir = tmp_path / "index"
    index_dir.mkdir()
    (index_dir / "user-file.txt").write_text("not an index\n", encoding="utf-8")

    with pytest.raises(FileExistsError, match="nonempty but invalid"):
        gene_ldscore_index._preflight_gene_index_output(index_dir, overwrite=True)


def test_index_build_lock_rejects_a_second_builder_for_same_absolute_target(tmp_path):
    index_dir = tmp_path / "index"
    with gene_ldscore_index._gene_index_build_lock(index_dir):
        with pytest.raises(LDSCInputError, match="Another build-gene-ldscore-index"):
            with gene_ldscore_index._gene_index_build_lock(index_dir):
                pass


def test_index_publication_recovery_restores_one_valid_owned_backup(tmp_path):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = publish_gene_ldscore_index(
        tmp_path / "index", index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    transaction = tmp_path / ".index.stage-interrupted"
    transaction.mkdir()
    marker = {
        "artifact_type": "gene_ldscore_index_publication",
        "target": str(index_dir.resolve()),
    }
    gene_ldscore_index._write_json(transaction / ".gene-index-publication.json", marker)
    backup = transaction / "index.backup"
    index_dir.rename(backup)
    gene_ldscore_index._write_json(backup / ".gene-index-publication.json", marker)

    gene_ldscore_index._recover_gene_index_publication(index_dir)

    assert load_gene_ldscore_index(index_dir).index_id == calculate_index_id(index_identity)
    assert not transaction.exists()


def test_failed_index_overwrite_rolls_back_without_mutating_valid_index(tmp_path, monkeypatch):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = tmp_path / "index"
    publish_gene_ldscore_index(
        index_dir, index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    original_id = load_gene_ldscore_index(index_dir).index_id

    def fail_write(*args, **kwargs):
        raise RuntimeError("injected publication failure")

    monkeypatch.setattr(gene_ldscore_index, "_write_index_artifact", fail_write)
    with pytest.raises(RuntimeError, match="injected"):
        publish_gene_ldscore_index(
            index_dir, index_identity={**index_identity, "padding_bp": 0}, gene_catalog=catalog,
            chromosomes={"22": chromosome}, overwrite=True,
        )

    assert load_gene_ldscore_index(index_dir).index_id == original_id


def test_indexed_gene_lists_assemble_control_queries_and_canonical_output(tmp_path):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = publish_gene_ldscore_index(
        tmp_path / "index", index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    genes = tmp_path / "focal.txt"
    genes.write_text("G1\n", encoding="utf-8")

    result = run_indexed_ldscore(
        index_dir,
        query_gene_list_sources=(genes,),
        control_gene_list_source="all-protein-coding",
        output_dir=tmp_path / "ldscores",
        overwrite=False,
    )

    assert result.baseline_columns == ["base", "gene_control"]
    assert result.query_columns == ["focal"]
    np.testing.assert_allclose(result.baseline_table["gene_control"], [0.75, 1.6])
    np.testing.assert_allclose(result.query_table["focal"], [0.75, 1.5])
    assert [record["column"] for record in result.count_records] == ["base", "gene_control", "focal"]
    assert result.overlap.baseline_block_all.loc["gene_control", "focal"] == 2
    assert (tmp_path / "ldscores" / "ldscore.baseline.parquet").exists()
    assert (tmp_path / "ldscores" / "ldscore.query.parquet").exists()
    metadata = json.loads((tmp_path / "ldscores" / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["index_id"] == calculate_index_id(index_identity)


def test_ldscore_parser_accepts_explicit_indexed_mode():
    args = cli.build_parser().parse_args(
        [
            "ldscore",
            "--output-dir", "out",
            "--gene-ldscore-index-dir", "gene-index",
            "--query-annot-gene-list-sources", "immune.txt,brain.txt",
            "--control-gene-list-source", "none",
        ]
    )

    assert args.gene_ldscore_index_dir == "gene-index"
    assert args.query_annot_gene_list_sources == "immune.txt,brain.txt"
    assert args.control_gene_list_source == "none"


def test_explicit_indexed_mode_dispatches_without_live_reference(monkeypatch, tmp_path):
    sentinel = object()
    captured = {}

    def fake_run(index_dir, **kwargs):
        captured["index_dir"] = index_dir
        captured.update(kwargs)
        return sentinel

    monkeypatch.setattr(gene_ldscore_index, "run_indexed_ldscore", fake_run)

    result = cli.main(
        [
            "ldscore",
            "--output-dir", str(tmp_path / "out"),
            "--gene-ldscore-index-dir", "gene-index",
            "--query-annot-gene-list-sources", "immune.txt,brain.txt",
            "--overwrite",
        ]
    )

    assert result is sentinel
    assert captured == {
        "index_dir": "gene-index",
        "query_gene_list_sources": ("immune.txt", "brain.txt"),
        "control_gene_list_source": "all-protein-coding",
        "output_dir": str(tmp_path / "out"),
        "overwrite": True,
    }


def test_explicit_indexed_cli_writes_the_canonical_workflow_log(tmp_path):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = publish_gene_ldscore_index(
        tmp_path / "index", index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    genes = tmp_path / "focal.txt"
    genes.write_text("G1\n", encoding="utf-8")
    output_dir = tmp_path / "out"

    cli.main(
        [
            "ldscore", "--output-dir", str(output_dir),
            "--gene-ldscore-index-dir", str(index_dir),
            "--query-annot-gene-list-sources", str(genes),
        ]
    )

    log_path = output_dir / "diagnostics" / "ldscore.log"
    assert log_path.exists()
    assert "gene_ldscore_index" in log_path.read_text(encoding="utf-8")


@pytest.mark.parametrize("forbidden", [
    ("--baseline-annot-sources", "baseline.annot.gz"),
    ("--plink-prefix", "panel"),
    ("--r2-dir", "r2"),
])
def test_explicit_indexed_mode_rejects_live_inputs(forbidden, tmp_path):
    with pytest.raises(LDSCInputError, match="indexed mode"):
        cli.main(
            [
                "ldscore",
                "--output-dir", str(tmp_path / "out"),
                "--gene-ldscore-index-dir", "gene-index",
                "--query-annot-gene-list-sources", "immune.txt",
                *forbidden,
            ]
        )


def test_indexed_all_unresolved_writes_diagnostics_without_scientific_outputs(tmp_path):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = publish_gene_ldscore_index(
        tmp_path / "index", index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    genes = tmp_path / "missing.txt"
    genes.write_text("NOT_A_GENE\n", encoding="utf-8")
    output_dir = tmp_path / "out"

    with pytest.raises(LDSCInputError, match="all 1 requested query annotations were skipped"):
        run_indexed_ldscore(
            index_dir,
            query_gene_list_sources=(genes,),
            control_gene_list_source="none",
            output_dir=output_dir,
        )

    assert (output_dir / "diagnostics" / "query_annotation_status.tsv").exists()
    assert (output_dir / "diagnostics" / "gene_list_unresolved.tsv.gz").exists()
    assert not (output_dir / "metadata.json").exists()
    assert not (output_dir / "ldscore.baseline.parquet").exists()
