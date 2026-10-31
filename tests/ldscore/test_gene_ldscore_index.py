from __future__ import annotations

from tests.ref_panel_helpers import prepare_plink, compute_plink
from tests.runtime.test_path_resolution import _write_plink_trio

from argparse import Namespace
from dataclasses import replace
import gc
import hashlib
import json
from pathlib import Path
import shutil
import threading
import time
from types import SimpleNamespace
import weakref

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from ldsc import cli
from ldsc._ldscore.gene_index import build as gene_ldscore_index
from ldsc._ldscore import calculator as ldscore_calculator
from ldsc._ldscore.gene_index import scoring as _gene_index_scoring
from ldsc._ldscore.gene_index.artifact import _canonical_frame_sha256, _catalog_identity_records
from ldsc.config import GeneLDScoreIndexBuildConfig, GlobalConfig, reset_global_config, set_global_config
from ldsc._core.plink_files import PlinkFAMFile
from ldsc._ldscore import compute as kernel_ldscore
from ldsc._annotation.overlap import compute_overlap
from ldsc._refpanel.prepared import ChromosomeAnnotations
from ldsc._ldscore.gene_index.atoms import AtomStatistics, build_disjoint_atoms
from ldsc.errors import LDSCInputError, LDSCUsageError
from ldsc._ldscore.gene_index.build import (
    _build_embedded_gene_catalog,
    _genetic_map_identity,
    _selected_individual_ids,
    build_plink_index_chromosome,
    intersect_baseline_plink_by_identifier,
)
from ldsc._ldscore.gene_index.artifact import calculate_index_id, IndexChromosomeData
from ldsc._annotation.gene_lists import GeneCatalog, GeneCatalogValidationError
from ldsc._ldscore.calculator import LDScoreCalculator
from tests.cli_parsing import parse_command
from ldsc._core.chromosomes import AUTOSOMES
from ldsc._ldscore.gene_index.scoring import load_gene_ldscore_index, run_indexed_ldscore
from tests.gene_index_helpers import narrow_index_chromosomes, write_gene_ldscore_index


@pytest.fixture(autouse=True)
def _chromosome_22_indexes(monkeypatch):
    """Index fixtures in this module cover chromosome 22 unless a test narrows otherwise."""
    narrow_index_chromosomes(monkeypatch, ("22",))


def _rows(*items: tuple[str, int, str]) -> pd.DataFrame:
    return pd.DataFrame(items, columns=["CHR", "POS", "SNP"])


def test_builder_defaults_snp_batch_size_to_512():
    assert gene_ldscore_index.build_parser().get_default("snp_batch_size") == 512
    assert GeneLDScoreIndexBuildConfig.__dataclass_fields__["snp_batch_size"].default == 512


def test_builder_parser_requires_catalog():
    parser = gene_ldscore_index.build_parser()

    assert "--gene-coordinate-file" in parser.format_help()
    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "--ld-wind-cm", "1",
                "--baseline-annot-sources", "baseline.@.annot.gz",
                "--plink-prefix", "panel.@",
                "--output-dir", "index",
                "--genome-build", "hg19",
                "--snp-identifier", "rsid",
            ]
        )


def test_builder_catalog_failure_writes_full_issues_before_atom_transaction(tmp_path, monkeypatch):
    catalog = tmp_path / "bad-catalog.tsv"
    catalog.write_text(
        "gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\n"
        "G1\tDUP\t1\t10\t20\thg19\n"
        "G1\tDUP\t2\t40\t30\thg19\n",
        encoding="utf-8",
    )
    baseline = tmp_path / "baseline.annot"
    baseline.write_text("CHR SNP POS base\n1 rs1 10 1\n")
    _write_plink_trio(tmp_path / "reference", ["1"])
    output = tmp_path / "index"
    args = gene_ldscore_index.build_parser().parse_args(
        ["--ld-wind-cm", "1",
            "--baseline-annot-sources", str(baseline),
            "--plink-prefix", str(tmp_path / "reference"),
            "--output-dir", str(output),
            "--gene-coordinate-file", str(catalog),
            "--genome-build", "hg19",
            "--snp-identifier", "rsid",
        ]
    )
    monkeypatch.setattr(
        gene_ldscore_index,
        "_create_gene_index_stage",
        lambda *_args, **_kwargs: pytest.fail("staged index created before catalog validation"),
    )

    with pytest.raises(GeneCatalogValidationError):
        gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)

    issues_path = output / ".build" / "gene_coordinate_catalog_issues.tsv.gz"
    issues = pd.read_csv(issues_path, sep="\t")
    assert set(issues["reason"]) == {"duplicate_gene_id", "duplicate_gene_name", "end_before_start"}
    assert [path.name for path in output.iterdir()] == [".build"]


def test_builder_rejects_conflicting_catalog_build(tmp_path, monkeypatch):
    catalog = tmp_path / "catalog.tsv"
    catalog.write_text(
        "gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\n"
        "G1\tGENE1\t1\t10\t20\thg38\n",
        encoding="utf-8",
    )
    baseline = tmp_path / "baseline.annot"
    baseline.write_text("CHR SNP POS base\n1 rs1 10 1\n")
    _write_plink_trio(tmp_path / "reference", ["1"])
    output = tmp_path / "index"
    args = gene_ldscore_index.build_parser().parse_args(
        ["--ld-wind-cm", "1",
            "--baseline-annot-sources", str(baseline),
            "--plink-prefix", str(tmp_path / "reference"),
            "--output-dir", str(output),
            "--gene-coordinate-file", str(catalog),
            "--genome-build", "hg19",
            "--snp-identifier", "rsid",
        ]
    )
    monkeypatch.setattr(gene_ldscore_index, '_create_gene_index_stage',
                        lambda *_a, **_k: pytest.fail('staged index created before build validation'))
    with pytest.raises(LDSCInputError, match='Conflicting genome-build evidence'):
        gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)
    assert 'hg38' in catalog.read_text()
    assert [path.name for path in output.iterdir()] == [".build"]


def _fixture_funnel(chrom, rows=2):
    """A PLINK reference-SNP funnel whose panel and grid both hold the fixture's ``rows`` SNPs."""
    from ldsc._refpanel.funnel import ReferenceSNPFunnel

    funnel = ReferenceSNPFunnel(chrom, "plink")
    for name, start in (("reference_panel", "panel_rows"), ("annotation_snp_grid", "grid_rows")):
        funnel.record(name, start, rows)
        for step in ("annotation_snp_grid", "genotype_qc"):
            funnel.record(name, step, 0)
    return funnel.frame()


def _artifact_payload():
    atom_model = build_disjoint_atoms("22", np.array([[0, 10], [5, 15]], dtype=np.int64))
    chromosome = IndexChromosomeData(
        baseline_rows=pd.DataFrame(
            {
                "CHR": ["22", "22"], "SNP": ["rs1", "rs2"], "POS": [1, 6],
                "A1": ["A", "C"], "A2": ["G", "T"],
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
        reference_snp_funnel=_fixture_funnel("22"),
    )
    catalog = pd.DataFrame(
        {
            "gene_index": [0, 1], "gene_id": ["ENSG1", "ENSG2"],
            "gene_name": ["G1", "G2"], "chrom": ["22", "22"], "start": [1, 6],
            "end": [10, 15], "genome_build": ["hg19", "hg19"],
            "catalog_line": [2, 3], "source": ["fixture.tsv", "fixture.tsv"],
            "included": [True, True], "exclusion_reason": ["", ""],
            "chromosome_gene_row": [0, 1],
        }
    )
    index_identity = {
        "chromosomes": ["22"], "baseline": "fixture", "plink": "fixture",
        "snp_identifier": "rsid", "genome_build": "hg19", "common_maf_min": 0.05,
        "catalog": _catalog_identity_records(catalog),
        "projection_build": "hg19", "padding_bp": 100000, "gene_exclude_regions": "mhc",
    }
    return chromosome, catalog, index_identity


def _configure_two_chromosome_builder(tmp_path, monkeypatch):
    narrow_index_chromosomes(monkeypatch, ("1", "2"))
    chromosome, catalog, index_identity = _artifact_payload()
    identity_rows = chromosome.baseline_rows[["CHR", "POS", "SNP"]].copy()
    embedded_catalog = pd.concat(
        [
            catalog.assign(
                gene_index=[0, 1],
                gene_id=[f"ENSG{chrom}1", f"ENSG{chrom}2"],
                gene_name=[f"G{chrom}1", f"G{chrom}2"],
                chrom=chrom,
                catalog_line=[2 + 2 * offset, 3 + 2 * offset],
                chromosome_gene_row=[0, 1],
            )
            for offset, chrom in enumerate(("1", "2"))
        ],
        ignore_index=True,
    )
    embedded_catalog["gene_index"] = np.arange(len(embedded_catalog))

    baseline_path = tmp_path / "baseline.annot.gz"
    pd.concat([
        identity_rows.assign(CHR=chrom, SNP=[f"rs{chrom}1", f"rs{chrom}2"], base=1.0)
        for chrom in ("1", "2")
    ], ignore_index=True).to_csv(baseline_path, sep="\t", index=False)
    monkeypatch.setattr(
        gene_ldscore_index.GeneCatalog,
        "load",
        lambda *_args, **_kwargs: SimpleNamespace(
            source="catalog.tsv.gz",
            genome_build="hg19",
        ),
    )
    monkeypatch.setattr(
        gene_ldscore_index,
        "_build_embedded_gene_catalog",
        lambda *args, **kwargs: embedded_catalog,
    )
    monkeypatch.setattr(gene_ldscore_index, "_load_builder_genetic_map", lambda _args: None)
    monkeypatch.setattr(
        gene_ldscore_index, "read_snp_restriction_keys", lambda *args, **kwargs: {"rs1"}
    )
    monkeypatch.setattr(
        gene_ldscore_index.kernel_regions, "load_preset_intervals", lambda *args: None
    )
    monkeypatch.setattr(
        gene_ldscore_index,
        "_read_bim_identity",
        lambda prefix: identity_rows.assign(CHR=str(prefix).rsplit(".", 1)[-1],
            SNP=[f"rs{str(prefix).rsplit('.', 1)[-1]}{i}" for i in (1, 2)]),
    )
    monkeypatch.setattr(
        gene_ldscore_index,
        "_builder_index_identity",
        lambda *args, **kwargs: {
            **index_identity,
            "chromosomes": ["1", "2"],
            "catalog": _catalog_identity_records(embedded_catalog),
        },
    )
    for chrom in ("1", "2"):
        _write_plink_trio(tmp_path / f"reference.{chrom}", [chrom])
    output = tmp_path / "index"
    import gzip
    with gzip.open(tmp_path / "catalog.tsv.gz", "wt") as stream:
        stream.write("catalog content supplied by the mocked loader\n")
    args = Namespace(
        baseline_annot_sources=str(baseline_path),
        plink_prefix=str(tmp_path / "reference.@"),
        output_dir=str(output),
        genome_build="hg19",
        snp_identifier="rsid",
        padding_bp=100000,
        gene_exclude_regions="mhc",
        ld_wind_cm=1.0,
        maf_min=None,
        common_maf_min=0.05,
        keep_indivs_file=None,
        regr_snps_file=None,
        regr_snps_exclude_regions="mhc-and-centromeres",
        genetic_map_hg19_sources=None,
        genetic_map_hg38_sources=None,
        gene_coordinate_file=str(tmp_path / "catalog.tsv.gz"),
        snp_batch_size=128,
        atom_batch_size=64,
        workers=2,
        overwrite=False,
        log_level="INFO",
    )
    return chromosome, output, args


def _build_fixture_index(tmp_path, monkeypatch):
    """Configure the two-chromosome build to publish in-memory chromosome records.

    Returns ``(output, args, index_identity)``. Running
    ``run_build_gene_ldscore_index_from_args(args)`` then takes the production
    path from output preflight through staging, commit and recovery.
    """
    from tests.worker_fixtures import local_task_pool
    monkeypatch.setattr("ldsc._runtime.parallelism.ProcessPoolExecutor", local_task_pool)
    chromosome, output, args = _configure_two_chromosome_builder(tmp_path, monkeypatch)

    def fixture_chromosome(chrom, *_args, **_kwargs):
        return replace(
            chromosome,
            baseline_rows=chromosome.baseline_rows.assign(CHR=chrom),
            atom_model=build_disjoint_atoms(chrom, np.array([[0, 10], [5, 15]], dtype=np.int64)),
            reference_snp_funnel=_fixture_funnel(chrom),
        )

    monkeypatch.setattr(gene_ldscore_index, "build_plink_index_chromosome", fixture_chromosome)
    return output, args, gene_ldscore_index._builder_index_identity()


def test_selected_individual_identity_uses_fam_row_order():
    prefix = Path(__file__).resolve().parents[1] / "fixtures" / "plink" / "plink"
    fam = PlinkFAMFile(str(prefix.with_suffix(".fam")))

    assert _selected_individual_ids(fam, [1, 0]) == [
        str(fam.IDList.iloc[1, 0]),
        str(fam.IDList.iloc[0, 0]),
    ]


def test_embedded_catalog_keeps_mhc_exclusions_with_empty_gene_rows():
    frame = pd.DataFrame(
        {
            "gene_id": ["ENSG1", "ENSG2", "ENSG3"],
            "gene_name": ["MHC_GENE", "OTHER6", "CHR22"],
            "chrom": ["6", "6", "22"],
            "start": [26_000_001, 40_000_001, 101],
            "end": [26_001_000, 40_001_000, 200],
            "genome_build": ["hg19", "hg19", "hg19"],
            "catalog_line": [2, 3, 4],
        }
    )
    catalog = GeneCatalog.from_embedded_frame(frame, source="test.tsv.gz")

    embedded = _build_embedded_gene_catalog(
        catalog, genome_build="hg19", gene_exclude_regions="mhc"
    )

    assert embedded["included"].tolist() == [False, True, True]
    assert embedded["exclusion_reason"].tolist() == ["excluded_gene_region", "", ""]
    assert embedded["chromosome_gene_row"].tolist() == [0, 1, 0]


def test_baseline_plink_intersection_uses_identifier_mode_and_plink_order(caplog):
    baseline_metadata = _rows(
        ("22", 30, "rs3"), ("22", 10, "rs1"), ("22", 20, "baseline_only")
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
        "baseline_duplicate_rows_dropped": 0,
        "plink_bim_duplicate_rows_dropped": 0,
        "matched_rows": 2,
        "baseline_only_rows": 1,
        "plink_only_rows": 1,
        "coordinate_discordant_rows": 0,
        "baseline_identity_rows_dropped": 0,
        "plink_bim_identity_rows_dropped": 0,
        "snp_label_discordant_rows": 0,
    }
    assert "coordinate disagreement" not in caplog.text


@pytest.mark.parametrize("side", ["baseline", "PLINK BIM"])
def test_baseline_plink_intersection_drops_duplicate_effective_identifier_group(side, caplog):
    baseline = _rows(("22", 10, "rs1"), ("22", 20, "rs2"))
    plink = _rows(("22", 10, "rs1"), ("22", 20, "rs2"))
    if side == "baseline":
        baseline = pd.concat([baseline, baseline.iloc[[0]]], ignore_index=True)
    else:
        plink = pd.concat([plink, plink.iloc[[0]]], ignore_index=True)
    result = intersect_baseline_plink_by_identifier(
        baseline,
        pd.DataFrame({"base": np.arange(len(baseline), dtype=float)}),
        plink,
        snp_identifier="rsid",
        chrom="22",
    )
    assert result.metadata["SNP"].tolist() == ["rs2"]
    assert result.diagnostics[f"{side.lower().replace(' ', '_')}_duplicate_rows_dropped"] == 2
    assert "duplicate_identity" in caplog.text


def test_chr_pos_intersection_ignores_baseline_label_and_publishes_plink_label(caplog):
    caplog.set_level("INFO")
    baseline = _rows(("22", 20, "baseline_twenty"), ("22", 10, "baseline_ten"))
    plink = _rows(("22", 10, "plink_ten"), ("22", 20, "plink_twenty"))
    result = intersect_baseline_plink_by_identifier(
        baseline,
        pd.DataFrame({"base": [20.0, 10.0]}),
        plink,
        snp_identifier="chr_pos",
        chrom="22",
    )
    assert result.metadata[["POS", "SNP"]].to_records(index=False).tolist() == [
        (10, "plink_ten"), (20, "plink_twenty")
    ]
    assert result.annotations["base"].tolist() == [10.0, 20.0]
    assert result.diagnostics["snp_label_discordant_rows"] == 2
    assert "PLINK SNP labels are authoritative and will be published" in caplog.text


def test_chr_pos_intersection_drops_all_rows_at_duplicate_coordinate():
    baseline = _rows(("22", 10, "base_a"), ("22", 20, "base_b"))
    plink = _rows(("22", 10, "plink_a"), ("22", 10, "plink_b"), ("22", 20, "plink_c"))
    result = intersect_baseline_plink_by_identifier(
        baseline,
        pd.DataFrame({"base": [1.0, 2.0]}),
        plink,
        snp_identifier="chr_pos",
        chrom="22",
    )
    assert result.metadata["SNP"].tolist() == ["plink_c"]
    assert result.diagnostics["plink_bim_duplicate_rows_dropped"] == 2


def test_chr_pos_intersection_allows_repeated_snp_labels_at_distinct_coordinates():
    baseline = _rows(("22", 10, "baseline-a"), ("22", 20, "baseline-b"))
    plink = _rows(("22", 10, "repeated-label"), ("22", 20, "repeated-label"))
    result = intersect_baseline_plink_by_identifier(
        baseline,
        pd.DataFrame({"base": [1.0, 2.0]}),
        plink,
        snp_identifier="chr_pos",
        chrom="22",
    )

    assert result.metadata["POS"].tolist() == [10, 20]
    assert result.metadata["SNP"].tolist() == ["repeated-label", "repeated-label"]


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
    prefix = Path(__file__).resolve().parents[1] / "fixtures" / "plink" / "plink"
    bim = pd.read_csv(
        prefix.with_suffix(".bim"),
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "POS", "A1", "A2"],
    )
    metadata = bim[["CHR", "SNP", "CM", "POS"]].copy()
    annotations = pd.DataFrame({"base": metadata["POS"].astype(np.float64)})
    bundle = ChromosomeAnnotations(metadata, annotations, ["base"], [])
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

    prepared = prepare_plink("1", bundle, args)
    direct = compute_plink("1", bundle, args, None)

    assert prepared.annotations.read().shape == (len(prepared.metadata), 1)
    np.testing.assert_array_equal(
        prepared.annotations.read()[:, 0],
        prepared.metadata["POS"].to_numpy(dtype=np.float64),
    )
    prepared.reader._currentSNP = 0
    reused_scores = kernel_ldscore.ld_scores_from_genotypes(
        prepared.reader,
        prepared.block_left,
        args.snp_batch_size,
        annot=prepared.annotations.read(),
    )
    np.testing.assert_array_equal(reused_scores.astype(np.float32), direct.ld_scores)


def test_direct_plink_projects_annotations_and_regression_weights_in_one_traversal(monkeypatch):
    prefix = Path(__file__).resolve().parents[1] / "fixtures" / "plink" / "plink"
    bim = pd.read_csv(
        prefix.with_suffix(".bim"),
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "POS", "A1", "A2"],
    )
    metadata = bim[["CHR", "SNP", "CM", "POS"]].copy()
    annotations = pd.DataFrame(
        {
            "base": np.ones(len(metadata), dtype=np.float32),
            "baseline_aux": (metadata["POS"] % 2).to_numpy(dtype=np.float32),
            "query": (metadata["POS"] >= 4).to_numpy(dtype=np.float32),
        }
    )
    bundle = ChromosomeAnnotations(
        metadata,
        annotations,
        ["base", "baseline_aux"],
        ["query"],
    )
    args = Namespace(
        bfile=str(prefix), keep=None, maf_min=None, maf=None,
        ld_wind_snps=10, ld_wind_kb=None, ld_wind_cm=None,
        yes_really=True, snp_batch_size=3, common_maf_min=0.05,
        snp_identifier="rsid", genetic_map=None,
    )

    # The former two-pass calculation is the numerical oracle for this
    # optimization: only traversal reuse may change, never the score values.
    prepared = prepare_plink("1", bundle, args)
    regression_keys = set(prepared.metadata.loc[::2, "SNP"].astype(str))
    regression_mask = kernel_ldscore.regression_mask_from_keys(
        prepared.metadata, regression_keys, args.snp_identifier
    )
    prepared.reader._currentSNP = 0
    expected_ld_scores = kernel_ldscore.ld_scores_from_genotypes(
        prepared.reader,
        prepared.block_left, args.snp_batch_size, annot=prepared.annotations.read()
    )
    assert prepared.reader._currentSNP == prepared.reader.m
    prepared.reader._currentSNP = 0
    expected_w_ld = kernel_ldscore.ld_scores_from_genotypes(
        prepared.reader,
        prepared.block_left, args.snp_batch_size, annot=regression_mask.reshape(-1, 1)
    )
    assert prepared.reader._currentSNP == prepared.reader.m

    calls = []
    original = kernel_ldscore.ld_scores_from_genotypes

    def traced_ldscore(reader, block_left, batch_size, annot=None, **options):
        start_cursor = reader._currentSNP
        scores = original(reader, block_left, batch_size, annot=annot, **options)
        calls.append((start_cursor, reader._currentSNP, np.column_stack([annot.read(), options["weight_mask"]])))
        return scores

    monkeypatch.setattr(kernel_ldscore, "ld_scores_from_genotypes", traced_ldscore)
    result = compute_plink(
        "1", bundle, args, regression_keys
    )

    assert len(calls) == 1
    start_cursor, end_cursor, projected = calls[0]
    assert (start_cursor, end_cursor) == (0, result.reference_snp_count)
    assert projected.dtype == np.float32
    np.testing.assert_array_equal(
        projected,
        np.column_stack([prepared.annotations.read(), regression_mask]),
    )
    assert result.ldscore_columns == ["base", "baseline_aux", "query"]
    assert result.baseline_columns == ["base", "baseline_aux"]
    assert result.query_columns == ["query"]
    assert result.ld_scores.dtype == np.float32
    assert result.w_ld.dtype == np.float32
    np.testing.assert_allclose(
        result.ld_scores,
        expected_ld_scores[regression_mask.astype(bool)].astype(np.float32),
        rtol=0.0,
        atol=np.finfo(np.float32).eps,
    )
    np.testing.assert_allclose(
        result.w_ld,
        expected_w_ld[regression_mask.astype(bool)].astype(np.float32),
        rtol=0.0,
        atol=np.finfo(np.float32).eps,
    )
    pd.testing.assert_frame_equal(result.metadata, prepared.metadata.loc[regression_mask.astype(bool)].reset_index(drop=True))
    expected_frame = pd.DataFrame(
        prepared.annotations.read(),
        columns=["base", "baseline_aux", "query"],
    )
    expected_M, expected_M_5_50 = kernel_ldscore.compute_counts(
        prepared.metadata, expected_frame, common_maf_min=args.common_maf_min
    )
    expected_overlap = compute_overlap(
        prepared.metadata,
        expected_frame,
        n_baseline=2,
        common_maf_min=args.common_maf_min,
    )
    np.testing.assert_array_equal(result.M, expected_M)
    np.testing.assert_array_equal(result.M_5_50, expected_M_5_50)
    np.testing.assert_array_equal(
        result.overlap.baseline_block_all, expected_overlap.baseline_block_all
    )
    np.testing.assert_array_equal(
        result.overlap.baseline_block_common, expected_overlap.baseline_block_common
    )
    np.testing.assert_array_equal(
        result.overlap.query_diagonal_all, expected_overlap.query_diagonal_all
    )
    np.testing.assert_array_equal(
        result.overlap.query_diagonal_common, expected_overlap.query_diagonal_common
    )
    public = LDScoreCalculator()._wrap_chrom_result(result, GlobalConfig(snp_identifier="rsid"))
    assert public.baseline_table.columns.tolist() == [
        "CHR", "SNP", "POS", "A1", "A2", "regression_ld_scores", "base", "baseline_aux"
    ]
    assert public.query_table.columns.tolist() == [
        "CHR", "SNP", "POS", "A1", "A2", "query"
    ]
    assert public.baseline_table["SNP"].tolist() == prepared.metadata.loc[
        regression_mask.astype(bool), "SNP"
    ].tolist()
    assert public.query_table["SNP"].tolist() == public.baseline_table["SNP"].tolist()
    np.testing.assert_array_equal(
        public.baseline_table["regression_ld_scores"].to_numpy(dtype=np.float32),
        result.w_ld[:, 0],
    )
    assert public.reference_snp_count == len(prepared.metadata)
    assert public.regression_selected_snp_count == int(regression_mask.sum())
    assert public.regression_region_removed_snp_count == 0
    assert [(record["group"], record["column"]) for record in public.count_records] == [
        ("baseline", "base"),
        ("baseline", "baseline_aux"),
        ("query", "query"),
    ]


def test_direct_plink_output_writes_header_only_duplicate_audit_when_no_rows_drop(tmp_path, monkeypatch):
    from ldsc import config, run_ldscore

    monkeypatch.setattr(config, "_GLOBAL_CONFIG", GlobalConfig(snp_identifier="rsid"))
    prefix = Path(__file__).resolve().parents[1] / "fixtures" / "plink" / "plink"
    bim = pd.read_csv(prefix.with_suffix(".bim"), sep=r"\s+", header=None,
                      names=["CHR", "SNP", "CM", "POS", "A1", "A2"])
    baseline = tmp_path / "baseline.annot"
    bim.assign(base=1).to_csv(baseline, sep="\t", index=False)
    regression = tmp_path / "regression.tsv"
    bim[["SNP"]].to_csv(regression, sep="\t", index=False)
    result = run_ldscore(output_dir=tmp_path / "out", baseline_annot_sources=str(baseline),
                        plink_prefix=str(prefix), regr_snps_file=str(regression),
                        regr_snps_exclude_regions="none", ld_wind_snps=10,
                        yes_really=True, snp_batch_size=3)
    sidecar = tmp_path / "out/diagnostics/dropped_snps/chr1_dropped.tsv.gz"
    assert sidecar.exists()
    assert pd.read_csv(sidecar, sep="\t").empty
    assert result.output_paths["dropped_snps_chr1"] == str(sidecar)
    assert not list((tmp_path / "out").glob(".ldsc-annotation-*"))


def test_direct_plink_duplicate_group_is_dropped_and_recorded(tmp_path):
    fixture = Path(__file__).resolve().parents[1] / "fixtures" / "plink" / "plink"
    prefix = tmp_path / "duplicate_panel"
    shutil.copyfile(fixture.with_suffix(".bed"), prefix.with_suffix(".bed"))
    shutil.copyfile(fixture.with_suffix(".fam"), prefix.with_suffix(".fam"))
    bim = pd.read_csv(
        fixture.with_suffix(".bim"), sep=r"\s+", header=None,
        names=["CHR", "SNP", "CM", "POS", "A1", "A2"],
    )
    duplicate_key = str(bim.loc[0, "SNP"])
    bim.loc[1, "SNP"] = duplicate_key
    bim.to_csv(prefix.with_suffix(".bim"), sep="\t", header=False, index=False)
    bundle = ChromosomeAnnotations(
        bim[["CHR", "SNP", "CM", "POS"]],
        pd.DataFrame({"base": np.ones(len(bim))}),
        ["base"],
        [],
    )
    args = Namespace(
        bfile=str(prefix), keep=None, maf_min=None, maf=None,
        ld_wind_snps=10, ld_wind_kb=None, ld_wind_cm=None,
        yes_really=True, snp_batch_size=3, common_maf_min=0.05,
        snp_identifier="rsid", genetic_map=None,
    )

    result = compute_plink("1", bundle, args, None)

    assert len(result.identity_drops) == 2
    assert set(result.identity_drops["reason"]) == {"duplicate_identity"}
    assert duplicate_key not in set(result.metadata["SNP"])


@pytest.mark.parametrize("snp_identifier", ["rsid", "chr_pos", "rsid_allele_aware", "chr_pos_allele_aware"])
@pytest.mark.parametrize("unit,value", [("snps", 2), ("kb", .002), ("cm", .15)])
@pytest.mark.parametrize("build", ["hg19", "hg38"])
def test_plink_atom_operator_matches_direct_union_in_atom_batches(tmp_path, snp_identifier, unit, value, build):
    fixture = Path(__file__).resolve().parents[1] / "fixtures" / "plink" / "plink"
    prefix = tmp_path / "panel"
    for suffix in (".bed", ".bim", ".fam"):
        shutil.copyfile(fixture.with_suffix(suffix), prefix.with_suffix(suffix))
    bim = pd.read_csv(
        prefix.with_suffix(".bim"),
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "POS", "A1", "A2"],
    )
    bim = bim.assign(A1="A", A2="G", CM=bim.POS * .1)
    bim.to_csv(prefix.with_suffix(".bim"), sep="\t", index=False, header=False)
    metadata = bim[["CHR", "SNP", "CM", "POS"]].copy()
    baseline = pd.DataFrame({"base": np.ones(len(metadata), dtype=np.float64)})
    base_bundle = ChromosomeAnnotations(metadata, baseline, ["base"], [])
    args = Namespace(
        bfile=str(prefix), keep=None, maf_min=None, maf=None,
        ld_wind_snps=None, ld_wind_kb=None, ld_wind_cm=None,
        yes_really=False, snp_batch_size=3, common_maf_min=0.05,
        snp_identifier=snp_identifier, genetic_map=None,
    )
    setattr(args, f"ld_wind_{unit}", value)
    if snp_identifier.startswith("chr_pos"):
        args.genome_build = build
    regression_keys = (
        set(metadata["SNP"].astype(str))
        if snp_identifier.startswith("rsid")
        else set(metadata["CHR"].astype(str) + ":" + metadata["POS"].astype(str))
    )
    from ldsc._core.snp_identity import RestrictionIdentityKeys, empty_identity_drop_frame
    regression_keys = RestrictionIdentityKeys(regression_keys, "base", empty_identity_drop_frame(), len(metadata), len(regression_keys))
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
    direct_bundle = ChromosomeAnnotations(
        metadata,
        pd.DataFrame({"base": baseline["base"], "query": direct_annotation}),
        ["base"],
        ["query"],
    )
    direct = compute_plink("1", direct_bundle, args, regression_keys)
    persisted = kernel_ldscore.regression_mask_from_keys(
        direct.metadata, regression_keys, snp_identifier
    ).astype(bool)

    np.testing.assert_array_equal(indexed_one.operator.toarray(), indexed_many.operator.toarray())
    np.testing.assert_array_equal(
        np.asarray(indexed_one.operator @ z.astype(np.float64)).astype(np.float32),
        direct.ld_scores[persisted, 1],
    )
    assert int(indexed_one.atom_statistics.atom_count_all @ z.astype(np.int64)) == int(direct.M[1])


def test_plink_index_fuses_baseline_and_regression_weights_but_keeps_atom_batches(monkeypatch):
    prefix = Path(__file__).resolve().parents[1] / "fixtures" / "plink" / "plink"
    bim = pd.read_csv(
        prefix.with_suffix(".bim"),
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "POS", "A1", "A2"],
    )
    metadata = bim[["CHR", "SNP", "CM", "POS"]].copy()
    baseline_frame = pd.DataFrame(
        {
            "base": np.ones(len(metadata), dtype=np.float32),
            "baseline_aux": (metadata["POS"] % 2).to_numpy(dtype=np.float32),
        }
    )
    bundle = ChromosomeAnnotations(
        metadata, baseline_frame, ["base", "baseline_aux"], []
    )
    args = Namespace(
        bfile=str(prefix), keep=None, maf_min=None, maf=None,
        ld_wind_snps=10, ld_wind_kb=None, ld_wind_cm=None,
        yes_really=True, snp_batch_size=3, common_maf_min=0.05,
        snp_identifier="rsid", genetic_map=None,
    )

    prepared = prepare_plink("1", bundle, args)
    regression_keys = set(prepared.metadata.loc[::2, "SNP"].astype(str))
    persisted = kernel_ldscore.regression_mask_from_keys(
        prepared.metadata, regression_keys, args.snp_identifier
    ).astype(bool)
    baseline64 = np.asarray(prepared.annotations.read(), dtype=np.float64)
    prepared.reader._currentSNP = 0
    expected_baseline = kernel_ldscore.ld_scores_from_genotypes(
        prepared.reader,
        prepared.block_left, args.snp_batch_size, annot=baseline64
    )
    prepared.reader._currentSNP = 0
    expected_weights = kernel_ldscore.ld_scores_from_genotypes(
        prepared.reader,
        prepared.block_left, args.snp_batch_size, annot=persisted.reshape(-1, 1)
    ).reshape(-1)

    calls = []
    options = []
    original = kernel_ldscore.ld_scores_from_genotypes

    def traced_ldscore(reader, block_left, batch_size, annot=None, **kwargs):
        start_cursor = reader._currentSNP
        scores = original(reader, block_left, batch_size, annot=annot, **kwargs)
        options.append(kwargs)
        assert len(scores) == int(persisted.sum())
        calls.append((start_cursor, reader._currentSNP, np.asarray(annot).copy()))
        return scores

    monkeypatch.setattr(kernel_ldscore, "ld_scores_from_genotypes", traced_ldscore)
    result = build_plink_index_chromosome(
        "1",
        bundle,
        args,
        regression_keys=regression_keys,
        regression_regions=None,
        gene_intervals=np.array([[0, 4], [2, 7], [3, 5]], dtype=np.int64),
        included=np.array([True, True, True]),
        padding_bp=0,
        atom_batch_size=1,
    )

    # One common pass projects baseline and the separate regression mask;
    # each remaining call is the intentional one-column gene-atom operator batch.
    assert len(calls) == 1 + result.atom_model.n_atoms
    assert all((start, end) == (0, len(result.reference_metadata)) for start, end, _ in calls)
    np.testing.assert_array_equal(
        calls[0][2], baseline64
    )
    np.testing.assert_array_equal(options[0]["weight_mask"], persisted)
    for kwargs in options:
        np.testing.assert_array_equal(kwargs["output_rows"], np.flatnonzero(persisted))
    assert calls[0][2].dtype == np.float64
    assert all(projected.shape == (len(result.reference_metadata), 1) for _, _, projected in calls[1:])
    assert all(projected.dtype == np.bool_ for _, _, projected in calls[1:])

    expected_rows = prepared.metadata.loc[
        persisted, ["CHR", "SNP", "POS", "A1", "A2"]
    ].reset_index(drop=True)
    expected_rows["regression_ld_scores"] = expected_weights[persisted]
    expected_rows["base"] = expected_baseline[persisted, 0]
    expected_rows["baseline_aux"] = expected_baseline[persisted, 1]
    pd.testing.assert_frame_equal(result.baseline_rows, expected_rows)
    pd.testing.assert_frame_equal(result.reference_metadata, prepared.metadata)
    expected_counts = kernel_ldscore.compute_counts(
        prepared.metadata,
        pd.DataFrame(baseline64, columns=bundle.baseline_columns),
        common_maf_min=args.common_maf_min,
    )
    np.testing.assert_array_equal(result.baseline_count_all, expected_counts[0])
    np.testing.assert_array_equal(result.baseline_count_common, expected_counts[1])


def test_build_index_command_requires_explicit_identity_and_build(capsys):
    parser = gene_ldscore_index.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "--baseline-annot-sources", "baseline.@.annot.gz",
                "--plink-prefix", "panel.@",
                "--output-dir", "suite",
            ]
        )
    error = capsys.readouterr().err
    assert "--snp-identifier is required; choose rsid, chr_pos, rsid_allele_aware, or chr_pos_allele_aware." in error
    assert "--genome-build is required; choose hg19 or hg38." in error


def test_build_index_help_exposes_matching_build_genetic_maps():
    help_text = gene_ldscore_index.build_parser().format_help()

    assert "--genetic-map-hg19-sources" in help_text
    assert "--genetic-map-hg38-sources" in help_text


@pytest.mark.parametrize("snp_identifier", ["rsid", "chr_pos"])
def test_build_index_command_registers_explicit_base_identity(snp_identifier):
    args = parse_command(
        [
            "build-gene-ldscore-index", "--ld-wind-cm", "1",
            "--baseline-annot-sources", "baseline.@.annot.gz",
            "--plink-prefix", "panel.@",
            "--output-dir", "suite",
            "--gene-coordinate-file", "genes.tsv",
            "--genome-build", "hg19",
            "--snp-identifier", snp_identifier,
        ]
    )

    assert args.command == "build-gene-ldscore-index"
    assert args.genome_build == "hg19"
    assert args.snp_identifier == snp_identifier
    assert args.padding_bp == 0
    assert args.gene_exclude_regions == "mhc"
    assert args.ld_wind_cm == 1.0
    assert args.common_maf_min == 0.05
    assert args.workers == 1
    assert args.atom_batch_size > 0
    assert args.regr_snps_file is None
    assert args.regr_snps_exclude_regions == "mhc-and-centromeres"


def test_build_index_command_accepts_custom_regression_snps_and_region_policy():
    args = parse_command(
        [
            "build-gene-ldscore-index", "--ld-wind-cm", "1",
            "--baseline-annot-sources", "baseline.@.annot.gz",
            "--plink-prefix", "panel.@",
            "--output-dir", "index",
            "--gene-coordinate-file", "genes.tsv",
            "--genome-build", "hg19",
            "--snp-identifier", "chr_pos",
            "--regr-snps-file", "custom.snplist",
            "--regr-snps-exclude-regions", "centromeres",
        ]
    )
    assert args.regr_snps_file == "custom.snplist"
    assert args.regr_snps_exclude_regions == "centromeres"


def test_build_index_python_config_requires_identity_and_build():
    with pytest.raises(TypeError):
        GeneLDScoreIndexBuildConfig(ld_wind_cm=1.0,
            baseline_annot_sources=("baseline.@.annot.gz",),
            plink_prefix="panel.@",
            output_dir="index",
        )

    config = GeneLDScoreIndexBuildConfig(ld_wind_cm=1.0,
        baseline_annot_sources=("baseline.@.annot.gz",),
        plink_prefix="panel.@",
        output_dir="index",
        gene_coordinate_file="genes.tsv",
        genome_build="hg19",
        snp_identifier="chr_pos",
    )
    assert config.snp_identifier == "chr_pos"
    assert config.padding_bp == 0


def test_build_index_writes_shared_operational_log_and_chromosome_metrics(tmp_path, monkeypatch):
    from tests.worker_fixtures import local_task_pool
    monkeypatch.setattr("ldsc._runtime.parallelism.ProcessPoolExecutor", local_task_pool)
    chromosome, catalog, index_identity = _artifact_payload()
    identity_rows = chromosome.baseline_rows[["CHR", "POS", "SNP"]].copy()
    baseline_path = tmp_path / "baseline.22.annot.gz"
    identity_rows.assign(base=1.0).to_csv(baseline_path, sep="\t", index=False)
    fake_catalog = SimpleNamespace(source="catalog.tsv.gz", genome_build="hg19")
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
        "regr_snps_exclude_regions": "mhc-and-centromeres",
    }

    monkeypatch.setattr(
        gene_ldscore_index.GeneCatalog,
        "load",
        lambda *_args, **_kwargs: fake_catalog,
    )
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
    _write_plink_trio(tmp_path / "1000G.EUR.QC.22", ["22"])

    import gzip
    with gzip.open(tmp_path / "catalog.tsv.gz", "wt") as stream:
        stream.write("catalog content supplied by the mocked loader\n")
    (tmp_path / "custom-regression.tsv").write_text("SNP\nrs1\n")
    args = Namespace(
        baseline_annot_sources=str(baseline_path),
        plink_prefix=str(tmp_path / "1000G.EUR.QC."),
        output_dir=str(tmp_path / "suite"),
        genome_build="hg19",
        snp_identifier="rsid",
        padding_bp=100000,
        gene_exclude_regions="mhc",
        ld_wind_cm=1.0,
        maf_min=None,
        common_maf_min=0.05,
        keep_indivs_file=None,
        regr_snps_file=str(tmp_path / "custom-regression.tsv"),
        regr_snps_exclude_regions="centromeres",
        genetic_map_hg19_sources=None,
        genetic_map_hg38_sources=None,
        gene_coordinate_file=str(tmp_path / "catalog.tsv.gz"),
        snp_batch_size=128,
        atom_batch_size=64,
        workers=1,
        overwrite=False,
        log_level="INFO",
    )

    unrelated_sibling = tmp_path / "suite.build"
    unrelated_sibling.mkdir()
    (unrelated_sibling / "build-gene-ldscore-index.log").write_text(
        "not owned by this build\n", encoding="utf-8"
    )
    index_dir = gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)

    log_path = index_dir / "diagnostics" / "build-gene-ldscore-index.log"
    json_path = index_dir / "diagnostics" / "build-gene-ldscore-index.json"
    log_text = log_path.read_text(encoding="utf-8")
    diagnostics = json.loads(json_path.read_text(encoding="utf-8"))

    assert "LDSC build-gene-ldscore-index Started" in log_text
    assert "Inputs:" in log_text
    assert "genome_build" in log_text
    assert "custom" in log_text
    assert "exclude-regions=centromeres" in log_text
    assert "baseline/PLINK identifier intersection" in log_text
    milestones = ["Reading annotation inputs:", "Checking SNP identities and preparing chromosome annotations.",
                  "Annotation preparation complete:", "Starting chromosome 22"]
    assert [log_text.index(message) for message in milestones] == sorted(log_text.index(message) for message in milestones)
    assert "Starting chromosome 22" in log_text
    assert "Finished chromosome 22" in log_text
    assert "catalog genes=2" in log_text
    assert "Finished " in log_text
    assert diagnostics["chromosomes"]["22"]["catalog_genes"] == 2
    assert diagnostics["chromosomes"]["22"]["operator_nnz"] == chromosome.operator.nnz
    assert f"nonzero SNP-atom LD scores={chromosome.operator.nnz}" in log_text
    assert captured_restriction["path"] == tmp_path / "custom-regression.tsv"
    dropped_sidecar = index_dir / "diagnostics" / "dropped_snps" / "chr22_dropped.tsv.gz"
    assert dropped_sidecar.exists()
    assert pd.read_csv(dropped_sidecar, sep="\t").empty
    assert not (index_dir / ".build").exists()
    assert (unrelated_sibling / "build-gene-ldscore-index.log").read_text(encoding="utf-8") == "not owned by this build\n"


def test_completed_chromosome_is_durably_staged_before_other_workers_finish(
    tmp_path, monkeypatch
):
    from tests.worker_fixtures import local_task_pool
    monkeypatch.setattr("ldsc._runtime.parallelism.ProcessPoolExecutor", local_task_pool)
    chromosome, output, args = _configure_two_chromosome_builder(tmp_path, monkeypatch)
    first_started = threading.Event()
    release_first = threading.Event()
    second_computed = threading.Event()
    second_record_released = threading.Event()
    durable_rename_seen = threading.Event()
    build_error: list[BaseException] = []

    def fake_build(chrom, *_args, **_kwargs):
        if chrom == "1":
            first_started.set()
            assert release_first.wait(timeout=10)
        record = replace(
            chromosome,
            baseline_rows=chromosome.baseline_rows.assign(CHR=chrom),
            atom_model=build_disjoint_atoms(
                chrom, np.array([[0, 10], [5, 15]], dtype=np.int64)
            ),
            reference_snp_funnel=_fixture_funnel(chrom),
        )
        if chrom == "2":
            weakref.finalize(record, second_record_released.set)
            second_computed.set()
        return record

    monkeypatch.setattr(gene_ldscore_index, "build_plink_index_chromosome", fake_build)
    live_log = output / ".build" / "build-gene-ldscore-index.log"
    real_replace = gene_ldscore_index.os.replace

    def observe_shard_rename(source, destination):
        if Path(destination).name == "chr2":
            if live_log.exists():
                assert "Finished chromosome 2" not in live_log.read_text(encoding="utf-8")
            result = real_replace(source, destination)
            assert Path(destination).is_dir()
            durable_rename_seen.set()
            return result
        return real_replace(source, destination)

    monkeypatch.setattr(gene_ldscore_index.os, "replace", observe_shard_rename)

    def run_build():
        try:
            gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)
        except BaseException as exc:
            build_error.append(exc)

    thread = threading.Thread(target=run_build)
    thread.start()
    try:
        assert first_started.wait(timeout=10)
        assert second_computed.wait(timeout=10)
        deadline = time.monotonic() + 10
        staged_shard = None
        while time.monotonic() < deadline:
            candidate = output / ".build" / "stage" / "chromosomes" / "chr2"
            if candidate.is_dir():
                staged_shard = candidate
                break
            time.sleep(0.01)

        assert staged_shard is not None
        assert (staged_shard / "ldscore_operator.npz").exists()
        assert not (output / "metadata.json").exists()
        assert durable_rename_seen.is_set()
        assert "Finished chromosome 2" in live_log.read_text(encoding="utf-8")
        gc.collect()
        assert second_record_released.wait(timeout=10)
    finally:
        release_first.set()
        thread.join(timeout=20)

    assert not thread.is_alive()
    assert not build_error
    assert load_gene_ldscore_index(output).chromosomes == ("1", "2")


def test_build_index_failure_keeps_stable_log_without_publishing_index(tmp_path, monkeypatch):
    from tests.worker_fixtures import local_task_pool
    monkeypatch.setattr("ldsc._runtime.parallelism.ProcessPoolExecutor", local_task_pool)
    chromosome, catalog, index_identity = _artifact_payload()
    identity_rows = chromosome.baseline_rows[["CHR", "POS", "SNP"]].copy()
    baseline_path = tmp_path / "baseline.22.annot.gz"
    identity_rows.assign(base=1.0).to_csv(baseline_path, sep="\t", index=False)
    fake_catalog = SimpleNamespace(source="catalog.tsv.gz", genome_build="hg19")
    embedded_catalog = catalog.copy()

    monkeypatch.setattr(
        gene_ldscore_index.GeneCatalog,
        "load",
        lambda *_args, **_kwargs: fake_catalog,
    )
    monkeypatch.setattr(gene_ldscore_index, "_build_embedded_gene_catalog", lambda *args, **kwargs: embedded_catalog)
    monkeypatch.setattr(gene_ldscore_index, "_builder_index_identity", lambda *args, **kwargs: index_identity)
    monkeypatch.setattr(gene_ldscore_index, "_load_builder_genetic_map", lambda _args: None)
    monkeypatch.setattr(gene_ldscore_index, "read_snp_restriction_keys", lambda *args, **kwargs: {"rs1"})
    monkeypatch.setattr(gene_ldscore_index.kernel_regions, "load_preset_intervals", lambda *args: None)
    monkeypatch.setattr(gene_ldscore_index, "_read_bim_identity", lambda _prefix: identity_rows.copy())
    _write_plink_trio(tmp_path / "1000G.EUR.QC.22", ["22"])

    original_intersection = gene_ldscore_index.intersect_baseline_plink_by_identifier

    def fail_intersection(*args, **kwargs):
        raise LDSCInputError("identifier intersection failure in fixture")

    monkeypatch.setattr(gene_ldscore_index, "intersect_baseline_plink_by_identifier", fail_intersection)
    import gzip
    with gzip.open(tmp_path / "catalog.tsv.gz", "wt") as stream:
        stream.write("catalog content supplied by the mocked loader\n")
    args = Namespace(
        baseline_annot_sources=str(baseline_path),
        plink_prefix=str(tmp_path / "1000G.EUR.QC."),
        output_dir=str(tmp_path / "suite"),
        genome_build="hg19",
        snp_identifier="rsid",
        padding_bp=100000,
        gene_exclude_regions="mhc",
        ld_wind_cm=1.0,
        maf_min=None,
        common_maf_min=0.05,
        keep_indivs_file=None,
        regr_snps_file=None,
        regr_snps_exclude_regions="mhc-and-centromeres",
        genetic_map_hg19_sources=None,
        genetic_map_hg38_sources=None,
        gene_coordinate_file=str(tmp_path / "catalog.tsv.gz"),
        snp_batch_size=128,
        atom_batch_size=64,
        workers=1,
        overwrite=False,
        log_level="INFO",
    )

    with pytest.raises(LDSCInputError, match="identifier intersection failure"):
        gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)

    index_path = tmp_path / "suite"
    build_state = index_path / ".build"
    log_path = build_state / "build-gene-ldscore-index.log"
    log_text = log_path.read_text(encoding="utf-8")
    assert "LDSC build-gene-ldscore-index Started" in log_text
    assert "Baseline/PLINK consistency validation failed" in log_text
    assert (build_state / "baseline_plink_issues.tsv").exists()
    assert "Failed " in log_text
    assert [path.name for path in index_path.iterdir()] == [".build"]

    with pytest.raises(LDSCInputError, match="identifier intersection failure"):
        gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)
    assert [path.name for path in index_path.iterdir()] == [".build"]
    assert any((build_state / "history").glob("*.log"))

    monkeypatch.setattr(gene_ldscore_index, "intersect_baseline_plink_by_identifier", original_intersection)
    monkeypatch.setattr(gene_ldscore_index, "build_plink_index_chromosome", lambda *args, **kwargs: chromosome)
    index_dir = gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)
    assert (index_dir / "metadata.json").exists()
    assert (index_dir / "diagnostics" / "build-gene-ldscore-index.log").exists()
    assert not build_state.exists()


def test_chromosome_failure_cleans_staged_shards_without_publication(tmp_path, monkeypatch):
    from tests.worker_fixtures import local_task_pool
    monkeypatch.setattr("ldsc._runtime.parallelism.ProcessPoolExecutor", local_task_pool)
    chromosome, output, args = _configure_two_chromosome_builder(tmp_path, monkeypatch)
    first_durable = threading.Event()

    def fake_build(chrom, *_args, **_kwargs):
        if chrom == "2":
            assert first_durable.wait(timeout=10)
            raise RuntimeError("injected chromosome failure")
        return replace(
            chromosome,
            baseline_rows=chromosome.baseline_rows.assign(CHR=chrom),
            atom_model=build_disjoint_atoms(
                chrom, np.array([[0, 10], [5, 15]], dtype=np.int64)
            ),
            reference_snp_funnel=_fixture_funnel(chrom),
        )

    monkeypatch.setattr(gene_ldscore_index, "build_plink_index_chromosome", fake_build)
    real_replace = gene_ldscore_index.os.replace

    def observe_first_shard(source, destination):
        result = real_replace(source, destination)
        if Path(destination).name == "chr1":
            first_durable.set()
        return result

    monkeypatch.setattr(gene_ldscore_index.os, "replace", observe_first_shard)

    with pytest.raises(RuntimeError, match="injected chromosome failure"):
        gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)

    assert first_durable.is_set()
    assert not (output / "metadata.json").exists()
    assert not (output / ".build" / "stage").exists()
    failed_log = output / ".build" / "build-gene-ldscore-index.log"
    assert "Finished chromosome 1" in failed_log.read_text(encoding="utf-8")


@pytest.mark.parametrize("snp_batch_size", [1, 3])
def test_plink_operator_matches_independent_dense_adjusted_r2_reference(snp_batch_size):
    prefix = Path(__file__).resolve().parents[1] / "fixtures" / "plink" / "plink"
    bim = pd.read_csv(
        prefix.with_suffix(".bim"), sep=r"\s+", header=None,
        names=["CHR", "SNP", "CM", "POS", "A1", "A2"],
    )
    bundle = ChromosomeAnnotations(
        bim[["CHR", "SNP", "CM", "POS"]],
        pd.DataFrame({"base": np.ones(len(bim))}),
        ["base"],
        [],
    )
    args = Namespace(
        bfile=str(prefix), keep=None, maf_min=None, maf=None,
        ld_wind_snps=1, ld_wind_kb=None, ld_wind_cm=None,
        yes_really=True, snp_batch_size=snp_batch_size, common_maf_min=0.05,
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
    prepared = prepare_plink("1", bundle, args)
    prepared.reader._currentSNP = 0
    x = prepared.reader.nextSNPs(prepared.reader.m)
    correlation = x.T @ x / prepared.reader.n
    adjusted_r2 = correlation**2 - (1.0 - correlation**2) / (prepared.reader.n - 2)
    within_window = np.abs(np.subtract.outer(np.arange(prepared.reader.m), np.arange(prepared.reader.m))) <= 1
    dense_r = np.where(within_window, adjusted_r2, 0.0)
    persisted = prepared.metadata["SNP"].isin(persisted_keys).to_numpy()

    np.testing.assert_allclose(indexed.operator.toarray(), dense_r[persisted], rtol=0, atol=1e-15)
    np.testing.assert_allclose(np.diag(dense_r), np.ones(prepared.reader.m), rtol=0, atol=1e-15)
    assert indexed.operator.toarray()[0, 2] < 0
    assert indexed.operator.toarray()[0, 0] == pytest.approx(1.0)  # non-regression rs_4 still contributes
    assert indexed.operator.toarray()[0, 3] == 0.0  # outside the one-SNP window
    # A batch wider than the window must not widen it, as legacy LDSC chunk rounding did.
    np.testing.assert_allclose(indexed.baseline_rows["base"], dense_r[persisted].sum(axis=1), rtol=0, atol=1e-12)


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


def test_semantic_id_and_baseline_content_are_mode_aware():
    rsid_identity = {
        "chromosomes": ["22"], "snp_identifier": "rsid", "genome_build": "hg19"
    }
    assert calculate_index_id(rsid_identity) != calculate_index_id(
        {**rsid_identity, "snp_identifier": "chr_pos"}
    )
    rows_a = _rows(("22", 10, "label-a"), ("22", 20, "label-b"))
    rows_b = rows_a.assign(SNP=["different-a", "different-b"])
    annotations = pd.DataFrame({"base": [1.0, 0.0]})
    assert gene_ldscore_index._retained_baseline_content_sha256(
        rows_a, annotations, "chr_pos"
    ) == gene_ldscore_index._retained_baseline_content_sha256(
        rows_b, annotations, "chr_pos"
    )
    assert gene_ldscore_index._retained_baseline_content_sha256(
        rows_a, annotations, "rsid"
    ) != gene_ldscore_index._retained_baseline_content_sha256(
        rows_b, annotations, "rsid"
    )


def test_explicit_genetic_map_identity_is_content_bound(tmp_path):
    first = tmp_path / "map1.txt"
    second = tmp_path / "map2.txt"
    first.write_text("CHR POS CM\n22 1 0.1\n", encoding="utf-8")
    second.write_text("CHR POS CM\n22 2 0.2\n", encoding="utf-8")

    bim_identity = _genetic_map_identity(Namespace(ld_wind_cm=1.0, genome_build="hg19", genetic_map_hg19_sources=None))
    explicit = _genetic_map_identity(
        Namespace(ld_wind_cm=1.0, genome_build="hg19", genetic_map_hg19_sources=f"{first},{second}")
    )

    assert bim_identity == "bim_cm"
    assert explicit["kind"] == "explicit_hg19"
    assert len(explicit["content_sha256"]) == 64


def _whole_text_canonical_sha256(frame, sort_by=()):
    """The digest of the complete normalized table text, which chunked hashing must reproduce."""
    normalized = frame.copy()
    if sort_by:
        normalized = normalized.sort_values(list(sort_by), kind="mergesort").reset_index(drop=True)
    normalized.columns = [str(column) for column in normalized.columns]
    for column in normalized.columns:
        if pd.api.types.is_float_dtype(normalized[column]):
            normalized[column] = normalized[column].map(
                lambda value: "<NA>" if pd.isna(value) else format(float(value), ".17g"))
        elif pd.api.types.is_bool_dtype(normalized[column]):
            normalized[column] = normalized[column].astype(bool).map({True: "true", False: "false"})
        else:
            normalized[column] = normalized[column].map(lambda value: "<NA>" if pd.isna(value) else str(value))
    payload = normalized.to_csv(index=False, sep="\t", lineterminator="\n")
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@pytest.mark.parametrize("sort_by", [(), ("CHR", "POS")])
def test_chunked_canonical_digest_equals_the_whole_table_digest(monkeypatch, sort_by):
    from ldsc._ldscore.gene_index import artifact

    monkeypatch.setattr(artifact, "_HASH_CHUNK_ROWS", 3)
    frame = pd.DataFrame({
        "CHR": ["2", "1", "1", "2", "1", "1", "2"],
        "POS": [5, 3, 3, 1, 2, 9, 4],
        "float": [0.0, -0.0, np.nan, np.inf, 1 / 3, 1e-310, 2.5],
        "float32": np.array([1.5, -0.0, np.nan, 3, 3, 3, 0.1], dtype=np.float32),
        "flag": [True, False, True, True, False, False, True],
        "text": ["a", "b\tc", 'q"x', "", None, "a", "z"],
        # One missing value makes pandas format this whole nullable column as floats.
        "count": pd.array([1, 2, 3, 4, 5, 6, None], dtype="Int64"),
    })
    assert _canonical_frame_sha256(frame, sort_by=sort_by) == _whole_text_canonical_sha256(frame, sort_by)


@pytest.mark.parametrize("pos, snp, ordered", [
    ([1, 2, 2, 3], ["c", "a", "b", "a"], True),
    ([1, 2, 2], ["a", "b", "a"], False),
    ([2, 1], ["a", "b"], False),
    ([1, 1, 1], ["a", None, None], True),
    ([1, 1], [None, "a"], False),
    ([], [], True),
])
def test_genomic_order_check_agrees_with_a_stable_sort(pos, snp, ordered):
    from ldsc._ldscore.gene_index.artifact import _in_position_snp_order

    frame = pd.DataFrame({"POS": np.array(pos, dtype=np.int64), "SNP": pd.Series(snp, dtype=object)})
    stable = frame.sort_values(["POS", "SNP"], kind="mergesort").reset_index(drop=True)
    assert frame.equals(stable) is ordered
    assert _in_position_snp_order(frame) is ordered


def test_common_atom_overlap_equals_the_common_row_product():
    from ldsc._ldscore.gene_index.atoms import compute_atom_statistics

    rng = np.random.default_rng(8)
    mapped = rng.integers(-1, 6, 400)
    baseline = rng.normal(size=(400, 3))
    common = rng.random(400) < .6
    statistics = compute_atom_statistics(mapped, 6, baseline, common)
    assigned = mapped >= 0
    membership = sparse.csr_matrix((np.ones(int(assigned.sum())), (np.flatnonzero(assigned), mapped[assigned])),
                                   shape=(400, 6))
    np.testing.assert_array_equal(statistics.baseline_atom_overlap_all, baseline.T @ membership)
    np.testing.assert_array_equal(statistics.baseline_atom_overlap_common, baseline[common].T @ membership[common])


def test_canonical_table_identity_ignores_row_order_when_key_sort_is_declared():
    frame = pd.DataFrame(
        {"CHR": ["22", "22"], "POS": [2, 1], "SNP": ["rs2", "rs1"], "base": [0.0, 1.0]}
    )
    reversed_frame = frame.iloc[::-1].reset_index(drop=True)

    assert _canonical_frame_sha256(
        frame, sort_by=("CHR", "POS", "SNP")
    ) == _canonical_frame_sha256(
        reversed_frame, sort_by=("CHR", "POS", "SNP")
    )


@pytest.mark.parametrize(
    "snp_identifier, expected_index_id",
    [
        ("rsid", "03c3143a86dfa28a23a02f3beccd2de477d4a4f71c2e839209a5353be8802824"),
        ("chr_pos", "d3b2d6ca4a3dda71028970e90ded1e14d297477917eb88cc18234ae57ef8add5"),
    ],
)
def test_index_artifact_round_trip_uses_approved_tree_and_payloads(tmp_path, snp_identifier, expected_index_id):
    chromosome, catalog, index_identity = _artifact_payload()
    index_identity = {**index_identity, "snp_identifier": snp_identifier}
    index_dir = write_gene_ldscore_index(
        tmp_path / "index",
        index_identity=index_identity,
        gene_catalog=catalog,
        chromosomes={"22": chromosome},
    )
    loaded = load_gene_ldscore_index(index_dir)

    assert index_dir == tmp_path / "index"
    assert loaded.index_id == calculate_index_id(index_identity)
    assert loaded.index_id == expected_index_id
    assert loaded.snp_identifier == snp_identifier
    assert loaded.genome_build == "hg19"
    assert loaded.chromosomes == ("22",)
    np.testing.assert_array_equal(loaded.chromosome_loader().load("22").operator.toarray(), chromosome.operator.toarray())
    root_metadata = (index_dir / "metadata.json").read_text(encoding="utf-8")
    assert "schema_version" not in root_metadata
    assert "software_version" not in root_metadata
    component_metadata = json.loads(
        (index_dir / "chromosomes" / "chr22" / "metadata.json").read_text(encoding="utf-8")
    )
    assert "effective_identity_sha256" not in component_metadata
    # Persisted digest from the fixture at the pre-removal revision 9d3fe47.
    assert component_metadata["published_row_metadata_sha256"] == (
        "5d56bf6d59eb23a9df273321249f8e341b1deaac155f59c9ad477903fff62463"
    )
    assert (index_dir / "chromosomes" / "chr22" / "baseline_statistics.npz").exists()
    assert (index_dir / "chromosomes" / "chr22" / "atom_statistics.npz").exists()
    stored_rows = pd.read_parquet(index_dir / "chromosomes" / "chr22" / "baseline_rows.parquet")
    assert stored_rows["regression_ld_scores"].dtype == np.float32
    assert stored_rows["base"].dtype == np.float32
    assert list(stored_rows.columns[:5]) == ["CHR", "SNP", "POS", "A1", "A2"]
    component_metadata["effective_identity_sha256"] = "unused_retired_field"
    (index_dir / "chromosomes" / "chr22" / "metadata.json").write_text(
        json.dumps(component_metadata), encoding="utf-8"
    )
    assert load_gene_ldscore_index(index_dir).index_id == expected_index_id


def test_public_loader_rejects_partial_index_before_gene_resolution(tmp_path, monkeypatch):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = write_gene_ldscore_index(
        tmp_path / "index",
        index_identity=index_identity,
        gene_catalog=catalog,
        chromosomes={"22": chromosome},
    )

    narrow_index_chromosomes(monkeypatch, AUTOSOMES)
    with pytest.raises(LDSCInputError, match="autosomes 1 through 22"):
        load_gene_ldscore_index(index_dir)


def test_chr_pos_index_round_trip_uses_coordinate_identity_and_publishes_plink_labels(tmp_path):
    chromosome, catalog, index_identity = _artifact_payload()
    chr_pos_identity = {**index_identity, "snp_identifier": "chr_pos"}
    index_dir = write_gene_ldscore_index(
        tmp_path / "index", index_identity=chr_pos_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome},
    )

    loaded = load_gene_ldscore_index(index_dir)
    assert loaded.snp_identifier == "chr_pos"
    assert loaded.genome_build == "hg19"
    assert loaded.chromosome_loader().load("22").baseline_rows["SNP"].tolist() == ["rs1", "rs2"]

    genes = tmp_path / "focal.txt"
    genes.write_text("G1\n", encoding="utf-8")
    result = run_indexed_ldscore(
        index_dir, query_gene_list_sources=(genes,), output_dir=tmp_path / "out"
    )
    assert result.config_snapshot.snp_identifier == "chr_pos"
    assert result.config_snapshot.genome_build == "hg19"
    metadata = json.loads((tmp_path / "out" / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["snp_identifier"] == "chr_pos"
    assert metadata["genome_build"] == "hg19"
    assert metadata["index_snp_identifier"] == "chr_pos"
    assert metadata["index_genome_build"] == "hg19"


@pytest.mark.parametrize("snp_identifier", ["rsid", "chr_pos"])
@pytest.mark.parametrize("column, changed_value", [("SNP", "rs_changed"), ("POS", 2), ("A1", "T"), ("A2", "C")])
def test_index_loader_rejects_component_identity_metadata_and_published_row_tampering(
    tmp_path, snp_identifier, column, changed_value
):
    chromosome, catalog, index_identity = _artifact_payload()
    index_identity = {**index_identity, "snp_identifier": snp_identifier}
    index_dir = write_gene_ldscore_index(
        tmp_path / "index", index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome},
    )
    component_metadata_path = index_dir / "chromosomes" / "chr22" / "metadata.json"
    component_metadata = json.loads(component_metadata_path.read_text(encoding="utf-8"))
    component_metadata["snp_identifier"] = "chr_pos" if snp_identifier == "rsid" else "rsid"
    component_metadata_path.write_text(json.dumps(component_metadata), encoding="utf-8")
    with pytest.raises(LDSCInputError, match="component identity"):
        load_gene_ldscore_index(index_dir)

    component_metadata["snp_identifier"] = snp_identifier
    component_metadata.pop("effective_identity_sha256", None)
    component_metadata_path.write_text(json.dumps(component_metadata), encoding="utf-8")
    assert load_gene_ldscore_index(index_dir).index_id == calculate_index_id(index_identity)
    rows_path = index_dir / "chromosomes" / "chr22" / "baseline_rows.parquet"
    rows = pd.read_parquet(rows_path)
    rows.loc[0, column] = changed_value
    rows.to_parquet(rows_path, index=False)
    with pytest.raises(LDSCInputError, match="published row metadata digest"):
        load_gene_ldscore_index(index_dir)


def test_chr_pos_index_loader_rejects_repeated_coordinates_even_when_labels_differ(tmp_path):
    chromosome, catalog, index_identity = _artifact_payload()
    chromosome = replace(
        chromosome,
        baseline_rows=chromosome.baseline_rows.assign(POS=[1, 1], SNP=["label-a", "label-b"]),
    )
    with pytest.raises(LDSCInputError, match="duplicate effective identities"):
        write_gene_ldscore_index(
            tmp_path / "index",
            index_identity={**index_identity, "snp_identifier": "chr_pos"},
            gene_catalog=catalog,
            chromosomes={"22": chromosome},
        )


def test_index_loader_rejects_semantic_id_and_sparse_dtype_corruption(tmp_path):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = write_gene_ldscore_index(
        tmp_path / "index", index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome},
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
    index_dir = write_gene_ldscore_index(
        tmp_path / "index", index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome},
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
    index_dir = write_gene_ldscore_index(
        tmp_path / "index", index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome},
    )
    metadata_path = index_dir / "metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["chromosomes"] = []
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

    with pytest.raises(LDSCInputError, match="chromosome coverage"):
        load_gene_ldscore_index(index_dir)


def test_index_publication_requires_overwrite_and_replaces_complete_index(tmp_path, monkeypatch):
    output, args, identity = _build_fixture_index(tmp_path, monkeypatch)
    gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)
    with pytest.raises(FileExistsError, match="--overwrite"):
        gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)

    replacement_identity = {**identity, "padding_bp": 0}
    monkeypatch.setattr(gene_ldscore_index, "_builder_index_identity", lambda *_a, **_k: replacement_identity)
    args.overwrite = True
    gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)
    assert load_gene_ldscore_index(output).index_id == calculate_index_id(replacement_identity)


def test_post_publication_cleanup_failure_warns_but_published_index_succeeds(tmp_path, monkeypatch):
    output, args, identity = _build_fixture_index(tmp_path, monkeypatch)
    gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)
    replacement_identity = {**identity, "padding_bp": 0}
    monkeypatch.setattr(gene_ldscore_index, "_builder_index_identity", lambda *_a, **_k: replacement_identity)
    args.overwrite = True
    real_rmtree = shutil.rmtree

    def fail_build_directory_cleanup(path, *args, **kwargs):
        if Path(path) == output / ".build":
            raise OSError(16, "injected device or resource busy", str(path))
        return real_rmtree(path, *args, **kwargs)

    with monkeypatch.context() as patcher:
        patcher.setattr(gene_ldscore_index.shutil, "rmtree", fail_build_directory_cleanup)
        with pytest.warns(UserWarning, match="build directory could not be removed"):
            published = gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)

    assert published == output
    assert load_gene_ldscore_index(output).index_id == calculate_index_id(replacement_identity)
    assert (output / ".build").is_dir()
    assert [path.name for path in tmp_path.iterdir() if path.name.startswith(".")] == []

    # The next build discards the leftover build directory and removes it after success.
    gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)
    assert not (output / ".build").exists()


def test_index_build_lock_rejects_a_second_builder_for_same_absolute_target(tmp_path):
    index_dir = tmp_path / "index"
    with gene_ldscore_index._gene_index_build_lock(index_dir):
        with pytest.raises(
            LDSCInputError,
            match=r"index/\.build/build-gene-ldscore-index\.log",
        ):
            with gene_ldscore_index._gene_index_build_lock(index_dir):
                pass
    assert (index_dir / ".build" / "build-gene-ldscore-index.lock").exists()
    assert [path.name for path in tmp_path.iterdir()] == ["index"]


def test_index_build_lock_retries_when_a_finishing_build_removes_the_lock_file(tmp_path, monkeypatch):
    index_dir = tmp_path / "index"
    lock_path = index_dir / ".build" / "build-gene-ldscore-index.lock"
    real_flock = gene_ldscore_index.fcntl.flock
    removed = []

    def flock_after_a_finishing_build(descriptor, operation):
        real_flock(descriptor, operation)
        if not removed and operation & gene_ldscore_index.fcntl.LOCK_EX:
            # A finishing build removed its build directory between this open and this lock.
            shutil.rmtree(lock_path.parent)
            removed.append(True)

    monkeypatch.setattr(gene_ldscore_index.fcntl, "flock", flock_after_a_finishing_build)
    with gene_ldscore_index._gene_index_build_lock(index_dir):
        assert removed
        assert "pid=" in lock_path.read_text(encoding="utf-8")


def test_final_log_move_failure_warns_and_retains_hidden_log(tmp_path, monkeypatch):
    index_dir = tmp_path / "index"
    (index_dir / "diagnostics").mkdir(parents=True)
    build_state = index_dir / ".build"
    build_state.mkdir()
    live_log = build_state / "build-gene-ldscore-index.log"
    live_log.write_text("Finished\n", encoding="utf-8")
    final_log = index_dir / "diagnostics" / live_log.name
    real_replace = gene_ldscore_index.os.replace

    def fail_final_log_move(source, destination):
        if Path(source) == live_log and Path(destination) == final_log:
            raise OSError(16, "injected device or resource busy", str(destination))
        return real_replace(source, destination)

    monkeypatch.setattr(gene_ldscore_index.os, "replace", fail_final_log_move)
    with pytest.warns(UserWarning, match="published index remains valid"):
        authoritative_log = gene_ldscore_index._finalize_gene_index_log(
            live_log, index_dir
        )

    assert authoritative_log == live_log
    assert live_log.exists()
    assert not final_log.exists()


def test_retry_discards_a_partial_stage_instead_of_reusing_it(tmp_path):
    index_dir = tmp_path / "index"
    stale_shard = index_dir / ".build" / "stage" / "chromosomes" / "chr22"
    stale_shard.mkdir(parents=True)
    (stale_shard / "stale-payload.txt").write_text("partial\n", encoding="utf-8")

    gene_ldscore_index._prepare_gene_index_log(index_dir)

    assert not (index_dir / ".build" / "stage").exists()
    assert [path.name for path in tmp_path.iterdir()] == ["index"]


def test_index_loader_rejects_embedded_catalog_policy_and_source_tampering(tmp_path):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = write_gene_ldscore_index(
        tmp_path / "index",
        index_identity=index_identity,
        gene_catalog=catalog,
        chromosomes={"22": chromosome},
    )
    catalog_path = index_dir / "gene_catalog.parquet"
    tampered = pd.read_parquet(catalog_path)
    tampered.loc[0, "included"] = False
    tampered.loc[0, "exclusion_reason"] = "excluded_gene_region"
    tampered.to_parquet(catalog_path, index=False)

    with pytest.raises(LDSCInputError, match="inclusion policy"):
        load_gene_ldscore_index(index_dir)

    tampered.loc[0, "included"] = True
    tampered.loc[0, "exclusion_reason"] = ""
    tampered["source"] = ""
    tampered.to_parquet(catalog_path, index=False)
    with pytest.raises(LDSCInputError, match="source basename"):
        load_gene_ldscore_index(index_dir)


def test_failed_index_overwrite_before_publication_keeps_the_valid_index(tmp_path, monkeypatch):
    output, args, identity = _build_fixture_index(tmp_path, monkeypatch)
    gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)
    original_id = load_gene_ldscore_index(output).index_id
    monkeypatch.setattr(
        gene_ldscore_index, "_builder_index_identity", lambda *_a, **_k: {**identity, "padding_bp": 0}
    )
    args.overwrite = True

    def fail_finalizing_the_stage(*_args, **_kwargs):
        raise OSError(16, "injected staging failure")

    monkeypatch.setattr(gene_ldscore_index, "_finalize_staged_index_artifact", fail_finalizing_the_stage)
    with pytest.raises(OSError, match="injected"):
        gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)

    # Publication had not begun, so the earlier index is still complete and loadable.
    assert load_gene_ldscore_index(output).index_id == original_id
    assert (output / "RUN_FAILED.txt").is_file()


def test_indexed_gene_lists_assemble_control_queries_and_canonical_output(tmp_path, monkeypatch):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = write_gene_ldscore_index(
        tmp_path / "index", index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome},
    )
    genes = tmp_path / "focal.txt"
    control = tmp_path / "control.txt"
    genes.write_text("G1\n", encoding="utf-8")
    control.write_text("G1\nG2\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("LDSC_INDEXED_OUT", str(tmp_path / "ldscores"))

    result = run_indexed_ldscore(
        index_dir,
        query_gene_list_sources=(genes,),
        control_gene_list_file=control,
        output_dir="$LDSC_INDEXED_OUT",
        overwrite=False,
    )

    assert result.baseline_columns == ["base", "gene_control"]
    assert result.query_columns == ["focal"]
    assert not (tmp_path / "$LDSC_INDEXED_OUT").exists()
    np.testing.assert_allclose(result.baseline_table["gene_control"], [0.75, 1.6])
    np.testing.assert_allclose(result.read_queries(["focal"])["focal"], [0.75, 1.5])
    assert [record["column"] for record in result.count_records] == ["base", "gene_control", "focal"]
    assert result.overlap.baseline_block_all.loc["gene_control", "focal"] == 2
    assert (tmp_path / "ldscores" / "ldscore.baseline.parquet").exists()
    assert (tmp_path / "ldscores" / "ldscore.query.parquet").exists()
    metadata = json.loads((tmp_path / "ldscores" / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["index_id"] == calculate_index_id(index_identity)


def test_ldscore_parser_accepts_control_gene_list_file_in_indexed_mode():
    args = parse_command(
        [
            "ldscore",
            "--output-dir", "out",
            "--gene-ldscore-index-dir", "gene-index",
            "--query-annot-gene-list-sources", "immune.txt,brain.txt",
            "--control-gene-list-file", "control.txt",
        ]
    )

    assert args.gene_ldscore_index_dir == "gene-index"
    assert args.query_annot_gene_list_sources == "immune.txt,brain.txt"
    assert args.control_gene_list_file == "control.txt"


def test_explicit_indexed_mode_dispatches_without_live_reference(monkeypatch, tmp_path):
    sentinel = object()
    captured = {}

    def fake_run(index_dir, **kwargs):
        captured["index_dir"] = index_dir
        captured.update(kwargs)
        return sentinel

    monkeypatch.setattr(_gene_index_scoring, "run_indexed_ldscore", fake_run)

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
        "control_gene_list_file": None,
        "gene_list_resolution_policy": "strict",
        "query_batch_size": 100,
        "workers": 1,
        "output_dir": str(tmp_path / "out"),
        "overwrite": True,
    }


def test_explicit_indexed_cli_writes_the_canonical_workflow_log(tmp_path, monkeypatch):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = write_gene_ldscore_index(
        tmp_path / "index", index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome},
    )
    genes = tmp_path / "focal.txt"
    genes.write_text("G1\n", encoding="utf-8")
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    marker = output_dir / "RUN_FAILED.txt"
    marker.write_text("earlier failed overwrite")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("LDSC_INDEXED_OUT", str(output_dir))

    cli.main(
        [
            "ldscore", "--output-dir", "$LDSC_INDEXED_OUT",
            "--gene-ldscore-index-dir", str(index_dir),
            "--query-annot-gene-list-sources", str(genes),
        ]
    )

    log_path = output_dir / "diagnostics" / "ldscore.log"
    assert not marker.exists()
    assert not (tmp_path / "$LDSC_INDEXED_OUT").exists()
    assert log_path.exists()
    assert "gene_ldscore_index" in log_path.read_text(encoding="utf-8")
    assert "gene_control" not in pd.read_parquet(output_dir / "ldscore.baseline.parquet").columns


def test_indexed_control_gene_list_file_must_resolve_to_a_file(tmp_path):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = write_gene_ldscore_index(
        tmp_path / "index", index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome},
    )
    genes = tmp_path / "focal.txt"
    genes.write_text("G1\n", encoding="utf-8")

    with pytest.raises(LDSCInputError, match="control gene-list file"):
        run_indexed_ldscore(
            index_dir,
            query_gene_list_sources=(genes,),
            control_gene_list_file=tmp_path / "missing-control.txt",
            output_dir=tmp_path / "out",
        )

    control = tmp_path / "control.txt"
    control.write_text("G1\n", encoding="utf-8")
    with pytest.raises(LDSCInputError, match="glob patterns"):
        run_indexed_ldscore(
            index_dir,
            query_gene_list_sources=(genes,),
            control_gene_list_file=tmp_path / "*.txt",
            output_dir=tmp_path / "out-glob",
        )


@pytest.mark.parametrize("forbidden", [
    ("--baseline-annot-sources", "baseline.annot.gz"),
    ("--plink-prefix", "panel"),
    ("--r2-matrix-dir", "r2"),
    ("--snp-identifier", "chr_pos_allele_aware"),
    ("--genome-build", "hg19"),
    ("--gene-exclude-regions", "none"),
    ("--common-maf-min", "0.05"),
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


def test_explicit_indexed_mode_rejects_ignored_boolean_live_option(tmp_path):
    with pytest.raises(LDSCInputError, match="--export-ref-metadata"):
        cli.main(
            [
                "ldscore",
                "--output-dir", str(tmp_path / "out"),
                "--gene-ldscore-index-dir", "gene-index",
                "--query-annot-gene-list-sources", "immune.txt",
                "--export-ref-metadata",
            ]
        )


def test_python_indexed_mode_rejects_explicit_identity_even_when_equal_to_ordinary_default(tmp_path):
    with pytest.raises(LDSCUsageError, match="snp_identifier"):
        ldscore_calculator.run_ldscore(
            output_dir=tmp_path / "out",
            gene_ldscore_index_dir=tmp_path / "index",
            query_annot_gene_list_sources=(tmp_path / "genes.txt",),
            snp_identifier="chr_pos_allele_aware",
        )


@pytest.mark.parametrize("registered", [
    GlobalConfig(),
    GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
    GlobalConfig(snp_identifier="rsid_allele_aware"),
])
def test_python_indexed_mode_inherits_index_identity_over_registered_global_config(
    registered, tmp_path, monkeypatch
):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = write_gene_ldscore_index(
        tmp_path / "index", index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome},
    )
    genes = tmp_path / "focal.txt"
    genes.write_text("G1\n", encoding="utf-8")

    set_global_config(registered)
    try:
        result = ldscore_calculator.run_ldscore(
            output_dir=tmp_path / "out",
            gene_ldscore_index_dir=index_dir,
            query_annot_gene_list_sources=str(genes),
        )
    finally:
        reset_global_config()

    assert result.config_snapshot.snp_identifier == "rsid"
    assert result.config_snapshot.genome_build is None
    assert (tmp_path / "out" / "ldscore.query.parquet").exists()


def test_indexed_strict_unresolved_stops_at_gate_a_with_complete_audit(tmp_path):
    chromosome, catalog, index_identity = _artifact_payload()
    index_dir = write_gene_ldscore_index(
        tmp_path / "index", index_identity=index_identity, gene_catalog=catalog,
        chromosomes={"22": chromosome},
    )
    genes = tmp_path / "missing.txt"
    genes.write_text("NOT_A_GENE\n", encoding="utf-8")
    output_dir = tmp_path / "out"

    with pytest.raises(LDSCInputError, match="catalog preflight rejected 1 of 1"):
        run_indexed_ldscore(
            index_dir,
            query_gene_list_sources=(genes,),
            output_dir=output_dir,
        )

    assert not (output_dir / "diagnostics" / "query_annotation_status.tsv").exists()
    audit_path = output_dir / "diagnostics" / "gene_list_audit.tsv.gz"
    assert audit_path.exists()
    audit = pd.read_csv(audit_path, sep="\t")
    assert audit[["input_gene", "disposition", "reason"]].to_records(index=False).tolist() == [
        ("NOT_A_GENE", "rejected", "unmatched_identifier")
    ]
    assert (output_dir / "diagnostics" / "gene_list_resolution_summary.tsv").exists()
    assert not (output_dir / "metadata.json").exists()
    assert not (output_dir / "ldscore.baseline.parquet").exists()


@pytest.mark.parametrize("overwrite", [False, True])
def test_index_preflight_ignores_failure_marker_without_removing_it(tmp_path, overwrite):
    marker = tmp_path / "RUN_FAILED.txt"
    marker.write_text("old failed run\n")
    gene_ldscore_index._preflight_gene_index_output(tmp_path, overwrite=overwrite)
    assert marker.read_text() == "old failed run\n"


@pytest.mark.parametrize("extra", [
    "user-file.txt", "user-directory", "metadata.json", "marker-directory", "marker-symlink",
])
def test_failure_marker_does_not_authorize_replacing_unrecognized_contents(tmp_path, extra):
    output = tmp_path / "index"
    output.mkdir()
    marker = output / "RUN_FAILED.txt"
    if extra == "marker-directory":
        marker.mkdir()
    elif extra == "marker-symlink":
        marker.symlink_to(tmp_path / "missing-marker-target")
    else:
        marker.write_text("old failure\n")
        if extra == "user-directory":
            (output / extra).mkdir()
        else:
            (output / extra).write_text("unrecognized contents\n")
    before = sorted(path.name for path in output.iterdir())
    with pytest.raises(FileExistsError, match="nonempty but invalid"):
        gene_ldscore_index._preflight_gene_index_output(output, overwrite=True)
    assert sorted(path.name for path in output.iterdir()) == before


def test_index_build_accepts_failure_marker_only_destination(tmp_path, monkeypatch):
    output, args, identity = _build_fixture_index(tmp_path, monkeypatch)
    output.mkdir()
    marker = output / "RUN_FAILED.txt"
    marker.write_text("old failure\n")
    args.overwrite = True
    published = gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)
    assert load_gene_ldscore_index(published).index_id == calculate_index_id(identity)
    assert not marker.exists()
