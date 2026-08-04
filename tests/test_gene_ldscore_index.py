from __future__ import annotations

from argparse import Namespace
import json
from pathlib import Path

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
    calculate_profile_id,
    calculate_suite_id,
    IndexChromosomeData,
    load_gene_ldscore_index,
    publish_gene_ldscore_index,
    run_indexed_ldscore,
    validate_strict_baseline_plink_identity,
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
    suite_identity = {"chromosomes": ["22"], "baseline": "fixture", "plink": "fixture"}
    profile_identity = {
        "catalog": {"release": "GENCODE v49", "content_sha256": "fixture"},
        "projection_build": "hg19", "padding_bp": 100000, "gene_exclude_regions": "mhc",
    }
    return chromosome, catalog, suite_identity, profile_identity


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


def test_strict_baseline_plink_identity_allows_reordering_only():
    baseline = _rows(("22", 30, "rs3"), ("22", 10, "rs1"), ("22", 20, "rs2"))
    plink = _rows(("22", 10, "rs1"), ("22", 20, "rs2"), ("22", 30, "rs3"))

    actual = validate_strict_baseline_plink_identity(baseline, plink, chrom="22")

    assert actual[["CHR", "POS", "SNP"]].to_records(index=False).tolist() == [
        ("22", 10, "rs1"),
        ("22", 20, "rs2"),
        ("22", 30, "rs3"),
    ]


@pytest.mark.parametrize(
    ("baseline", "match"),
    [
        (_rows(("22", 10, "rs1"), ("22", 20, "rs2")), "missing or extra"),
        (_rows(("22", 10, "rs1"), ("22", 20, "different"), ("22", 30, "rs3")), "conflicting"),
        (_rows(("22", 10, "rs1"), ("22", 10, "rs1"), ("22", 30, "rs3")), "duplicate"),
    ],
)
def test_strict_baseline_plink_identity_rejects_any_pre_qc_difference(baseline, match):
    plink = _rows(("22", 10, "rs1"), ("22", 20, "rs2"), ("22", 30, "rs3"))

    with pytest.raises(LDSCInputError, match=match):
        validate_strict_baseline_plink_identity(baseline, plink, chrom="22")


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


def test_semantic_ids_use_canonical_scientific_identity_only():
    suite_identity = {
        "baseline_sources": [{"ordinal": 1, "content_sha256": "abc"}],
        "plink": {"bed_sha256": "bed", "bim_sha256": "bim", "fam_sha256": "fam"},
        "chromosomes": ["22"],
        "ld_window": {"mode": "cm", "value": 1.0},
    }
    suite_id = calculate_suite_id(suite_identity)
    assert suite_id == calculate_suite_id(dict(reversed(list(suite_identity.items()))))
    assert suite_id != calculate_suite_id({**suite_identity, "chromosomes": ["21", "22"]})

    profile_identity = {
        "catalog": {"release": "GENCODE v49", "content_sha256": "catalog"},
        "projection_build": "hg19",
        "padding_bp": 100000,
        "gene_exclude_regions": "mhc",
    }
    profile_id = calculate_profile_id(suite_id, profile_identity)
    assert profile_id == calculate_profile_id(suite_id, dict(reversed(list(profile_identity.items()))))
    assert profile_id != calculate_profile_id(suite_id, {**profile_identity, "padding_bp": 0})


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
    assert [entry["name"] for entry in explicit["sources"]] == ["map1.txt", "map2.txt"]
    assert all(len(entry["sha256"]) == 64 for entry in explicit["sources"])


def test_index_artifact_round_trip_uses_approved_tree_and_payloads(tmp_path):
    atom_model = build_disjoint_atoms("22", np.array([[0, 10], [5, 15]], dtype=np.int64))
    chromosome = IndexChromosomeData(
        baseline_rows=pd.DataFrame(
            {
                "CHR": ["22", "22"], "SNP": ["rs1", "rs2"], "POS": [1, 6],
                "regression_ld_scores": [1.0, 2.0], "base": [3.0, 4.0],
            }
        ),
        baseline_count_all=np.array([2.0]),
        baseline_count_common=np.array([1.0]),
        baseline_overlap_all=np.array([[2.0]]),
        baseline_overlap_common=np.array([[1.0]]),
        total_reference_snps_all=2,
        total_reference_snps_common=1,
        atom_model=atom_model,
        operator=sparse.csr_matrix(np.array([[1.0, -0.25, 0.0], [0.5, 1.0, 0.1]])),
        atom_statistics=AtomStatistics(
            atom_count_all=np.array([1, 1, 0], dtype=np.int64),
            atom_count_common=np.array([1, 0, 0], dtype=np.int64),
            baseline_atom_overlap_all=np.array([[1.0, 1.0, 0.0]]),
            baseline_atom_overlap_common=np.array([[1.0, 0.0, 0.0]]),
        ),
        reference_metadata=pd.DataFrame(),
    )
    catalog = pd.DataFrame(
        {
            "gene_index": [0, 1],
            "canonical_ensembl_id": ["ENSG1", "ENSG2"],
            "gene_name": ["G1", "G2"],
            "CHR": ["22", "22"],
            "start0": [0, 5],
            "end": [10, 15],
            "included": [True, True],
            "exclusion_reason": [None, None],
            "chromosome_gene_row": [0, 1],
        }
    )
    suite_identity = {"chromosomes": ["22"], "baseline": "fixture", "plink": "fixture"}
    profile_identity = {
        "catalog": {"release": "GENCODE v49", "content_sha256": "fixture"},
        "projection_build": "hg19", "padding_bp": 100000, "gene_exclude_regions": "mhc",
    }

    profile_dir = publish_gene_ldscore_index(
        tmp_path / "suite",
        profile_name="padding-100000-mhc",
        suite_identity=suite_identity,
        profile_identity=profile_identity,
        gene_catalog=catalog,
        chromosomes={"22": chromosome},
        overwrite=False,
    )
    loaded = load_gene_ldscore_index(profile_dir)

    assert profile_dir == tmp_path / "suite" / "profiles" / "padding-100000-mhc"
    assert loaded.suite_id == calculate_suite_id(suite_identity)
    assert loaded.profile_id == calculate_profile_id(loaded.suite_id, profile_identity)
    assert loaded.chromosomes == ("22",)
    np.testing.assert_array_equal(loaded.profile_chromosomes["22"].operator.toarray(), chromosome.operator.toarray())
    root_metadata = (tmp_path / "suite" / "metadata.json").read_text(encoding="utf-8")
    assert "schema_version" not in root_metadata
    assert "software_version" not in root_metadata
    assert (tmp_path / "suite" / "common" / "chr22" / "baseline_statistics.npz").exists()
    assert (profile_dir / "chr22" / "atom_statistics.npz").exists()
    stored_rows = pd.read_parquet(
        tmp_path / "suite" / "common" / "chr22" / "baseline_rows.parquet"
    )
    assert stored_rows["regression_ld_scores"].dtype == np.float32
    assert stored_rows["base"].dtype == np.float32


def test_index_loader_rejects_semantic_id_and_sparse_dtype_corruption(tmp_path):
    chromosome, catalog, suite_identity, profile_identity = _artifact_payload()
    profile_dir = publish_gene_ldscore_index(
        tmp_path / "suite", profile_name="profile", suite_identity=suite_identity,
        profile_identity=profile_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    metadata_path = profile_dir / "metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["profile_id"] = "corrupt"
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(LDSCInputError, match="semantic identity"):
        load_gene_ldscore_index(profile_dir)

    metadata["profile_id"] = calculate_profile_id(calculate_suite_id(suite_identity), profile_identity)
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    operator_path = profile_dir / "chr22" / "ldscore_operator.npz"
    sparse.save_npz(operator_path, chromosome.operator.astype(np.float32))
    with pytest.raises(Exception, match="float64"):
        load_gene_ldscore_index(profile_dir)


def test_index_loader_rejects_noncanonical_catalog_and_baseline_row_order(tmp_path):
    chromosome, catalog, suite_identity, profile_identity = _artifact_payload()
    profile_dir = publish_gene_ldscore_index(
        tmp_path / "suite", profile_name="profile", suite_identity=suite_identity,
        profile_identity=profile_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    catalog_path = profile_dir / "gene_catalog.parquet"
    reversed_catalog = pd.read_parquet(catalog_path).iloc[::-1].reset_index(drop=True)
    reversed_catalog.to_parquet(catalog_path, index=False)
    with pytest.raises(LDSCInputError, match="gene_index ordering"):
        load_gene_ldscore_index(profile_dir)

    catalog.to_parquet(catalog_path, index=False)
    baseline_path = tmp_path / "suite" / "common" / "chr22" / "baseline_rows.parquet"
    reversed_rows = pd.read_parquet(baseline_path).iloc[::-1].reset_index(drop=True)
    reversed_rows.to_parquet(baseline_path, index=False)
    with pytest.raises(LDSCInputError, match="canonical genomic order"):
        load_gene_ldscore_index(profile_dir)


def test_index_loader_rejects_common_coverage_corruption(tmp_path):
    chromosome, catalog, suite_identity, profile_identity = _artifact_payload()
    profile_dir = publish_gene_ldscore_index(
        tmp_path / "suite", profile_name="profile", suite_identity=suite_identity,
        profile_identity=profile_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    common_metadata_path = tmp_path / "suite" / "common" / "metadata.json"
    common_metadata = json.loads(common_metadata_path.read_text(encoding="utf-8"))
    common_metadata["chromosomes"] = []
    common_metadata_path.write_text(json.dumps(common_metadata), encoding="utf-8")

    with pytest.raises(LDSCInputError, match="common chromosome coverage"):
        load_gene_ldscore_index(profile_dir)


def test_profile_publication_reuses_common_and_preserves_siblings_on_targeted_overwrite(tmp_path):
    chromosome, catalog, suite_identity, profile_identity = _artifact_payload()
    suite_dir = tmp_path / "suite"
    first = publish_gene_ldscore_index(
        suite_dir, profile_name="first", suite_identity=suite_identity,
        profile_identity=profile_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    common_bytes = (suite_dir / "common" / "chr22" / "baseline_statistics.npz").read_bytes()
    sibling_identity = {**profile_identity, "padding_bp": 0}
    sibling = publish_gene_ldscore_index(
        suite_dir, profile_name="sibling", suite_identity=suite_identity,
        profile_identity=sibling_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    original_sibling_id = load_gene_ldscore_index(sibling).profile_id
    replacement_identity = {**profile_identity, "gene_exclude_regions": "none"}
    publish_gene_ldscore_index(
        suite_dir, profile_name="first", suite_identity=suite_identity,
        profile_identity=replacement_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=True,
    )

    assert (suite_dir / "common" / "chr22" / "baseline_statistics.npz").read_bytes() == common_bytes
    assert load_gene_ldscore_index(sibling).profile_id == original_sibling_id
    assert load_gene_ldscore_index(first).profile_id == calculate_profile_id(
        calculate_suite_id(suite_identity), replacement_identity
    )


def test_failed_profile_overwrite_rolls_back_without_mutating_valid_suite(tmp_path, monkeypatch):
    chromosome, catalog, suite_identity, profile_identity = _artifact_payload()
    suite_dir = tmp_path / "suite"
    profile_dir = publish_gene_ldscore_index(
        suite_dir, profile_name="profile", suite_identity=suite_identity,
        profile_identity=profile_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    original_id = load_gene_ldscore_index(profile_dir).profile_id

    def fail_write(*args, **kwargs):
        raise RuntimeError("injected publication failure")

    monkeypatch.setattr(gene_ldscore_index, "_write_profile_layer", fail_write)
    with pytest.raises(RuntimeError, match="injected"):
        publish_gene_ldscore_index(
            suite_dir, profile_name="profile", suite_identity=suite_identity,
            profile_identity={**profile_identity, "padding_bp": 0}, gene_catalog=catalog,
            chromosomes={"22": chromosome}, overwrite=True,
        )

    assert load_gene_ldscore_index(profile_dir).profile_id == original_id


def test_indexed_gene_lists_assemble_control_queries_and_canonical_output(tmp_path):
    chromosome, catalog, suite_identity, profile_identity = _artifact_payload()
    profile_dir = publish_gene_ldscore_index(
        tmp_path / "suite", profile_name="profile", suite_identity=suite_identity,
        profile_identity=profile_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    genes = tmp_path / "focal.txt"
    genes.write_text("G1\n", encoding="utf-8")

    result = run_indexed_ldscore(
        profile_dir,
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
    assert metadata["suite_id"] == calculate_suite_id(suite_identity)
    assert metadata["profile_id"] == calculate_profile_id(metadata["suite_id"], profile_identity)


def test_ldscore_parser_accepts_explicit_indexed_mode():
    args = cli.build_parser().parse_args(
        [
            "ldscore",
            "--output-dir", "out",
            "--gene-ldscore-index-dir", "suite/profiles/baseline",
            "--query-annot-gene-list-sources", "immune.txt,brain.txt",
            "--control-gene-list-source", "none",
        ]
    )

    assert args.gene_ldscore_index_dir == "suite/profiles/baseline"
    assert args.query_annot_gene_list_sources == "immune.txt,brain.txt"
    assert args.control_gene_list_source == "none"


def test_explicit_indexed_mode_dispatches_without_live_reference(monkeypatch, tmp_path):
    sentinel = object()
    captured = {}

    def fake_run(profile_dir, **kwargs):
        captured["profile_dir"] = profile_dir
        captured.update(kwargs)
        return sentinel

    monkeypatch.setattr(gene_ldscore_index, "run_indexed_ldscore", fake_run)

    result = cli.main(
        [
            "ldscore",
            "--output-dir", str(tmp_path / "out"),
            "--gene-ldscore-index-dir", "suite/profiles/baseline",
            "--query-annot-gene-list-sources", "immune.txt,brain.txt",
            "--overwrite",
        ]
    )

    assert result is sentinel
    assert captured == {
        "profile_dir": "suite/profiles/baseline",
        "query_gene_list_sources": ("immune.txt", "brain.txt"),
        "control_gene_list_source": "all-protein-coding",
        "output_dir": str(tmp_path / "out"),
        "overwrite": True,
    }


def test_explicit_indexed_cli_writes_the_canonical_workflow_log(tmp_path):
    chromosome, catalog, suite_identity, profile_identity = _artifact_payload()
    profile_dir = publish_gene_ldscore_index(
        tmp_path / "suite", profile_name="profile", suite_identity=suite_identity,
        profile_identity=profile_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    genes = tmp_path / "focal.txt"
    genes.write_text("G1\n", encoding="utf-8")
    output_dir = tmp_path / "out"

    cli.main(
        [
            "ldscore", "--output-dir", str(output_dir),
            "--gene-ldscore-index-dir", str(profile_dir),
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
                "--gene-ldscore-index-dir", "suite/profiles/baseline",
                "--query-annot-gene-list-sources", "immune.txt",
                *forbidden,
            ]
        )


def test_indexed_all_unresolved_writes_diagnostics_without_scientific_outputs(tmp_path):
    chromosome, catalog, suite_identity, profile_identity = _artifact_payload()
    profile_dir = publish_gene_ldscore_index(
        tmp_path / "suite", profile_name="profile", suite_identity=suite_identity,
        profile_identity=profile_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome}, overwrite=False,
    )
    genes = tmp_path / "missing.txt"
    genes.write_text("NOT_A_GENE\n", encoding="utf-8")
    output_dir = tmp_path / "out"

    with pytest.raises(LDSCInputError, match="all 1 requested query annotations were skipped"):
        run_indexed_ldscore(
            profile_dir,
            query_gene_list_sources=(genes,),
            control_gene_list_source="none",
            output_dir=output_dir,
        )

    assert (output_dir / "diagnostics" / "query_annotation_status.tsv").exists()
    assert (output_dir / "diagnostics" / "gene_list_unresolved.tsv.gz").exists()
    assert not (output_dir / "metadata.json").exists()
    assert not (output_dir / "ldscore.baseline.parquet").exists()
