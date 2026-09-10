"""Independent multi-chromosome expectations for streamed quantile inputs."""

import json
from pathlib import Path
import weakref

import numpy as np
import pandas as pd
import pytest

import tests.test_quantile_h2 as fixtures
from ldsc.errors import LDSCInputError
from ldsc.quantile_h2 import run_quantile_h2_from_args


@pytest.fixture
def workflow():
    case = fixtures.QuantileH2WorkflowTest()
    case.setUp()
    for name in ("annotations.tsv", "target.tsv", "ref.tsv"):
        path = case.root / name
        frame = pd.read_csv(path, sep="\t")
        frame["CHR"] = [1, 2, 1, 2]
        frame.to_csv(path, sep="\t", index=False)
    yield case
    case.doCleanups()


def test_exact_float64_targets_across_chromosomes(workflow):
    path = workflow.root / "target.tsv"
    target = pd.read_csv(path, sep="\t")
    # float32 would collapse all four values and produce an empty quantile.
    target["target"] = [1.0, 1.000000001, 1.000000002, 1.000000003]
    target.to_csv(path, sep="\t", index=False)
    result = run_quantile_h2_from_args(workflow.args)
    np.testing.assert_allclose(result.quantile_h2.h2_obs, [.45, .35])
    assert result.quantile_h2.n_snps.tolist() == [3, 1]
    assert Path(result.alignment_issues_path).is_file()
    assert not hasattr(result, "alignment_issues")
    assert not list(Path(workflow.args.output_dir).glob(".ldsc-annotation-*"))


@pytest.mark.parametrize("filename", ["ref.tsv", "target.tsv"])
def test_global_duplicates_report_both_chromosomes(workflow, filename):
    path = workflow.root / filename
    frame = pd.read_csv(path, sep="\t")
    frame.loc[1, "SNP"] = "rs1"
    frame.to_csv(path, sep="\t", index=False)
    with pytest.raises(LDSCInputError, match="duplicate effective SNP identities"):
        run_quantile_h2_from_args(workflow.args)
    issues = pd.read_csv(Path(workflow.args.output_dir) / "diagnostics/snp_alignment_issues.tsv.gz", sep="\t")
    assert issues.SNP.tolist() == ["rs1", "rs1"]
    assert issues.CHR.tolist() == [1, 2]
    assert not list(Path(workflow.args.output_dir).glob(".ldsc-annotation-*"))


def test_excluded_targets_still_contribute_to_standardization(workflow):
    path = workflow.root / "target.tsv"
    target = pd.read_csv(path, sep="\t")
    target.loc[3, "target"] = -999
    target.to_csv(path, sep="\t", index=False)
    workflow.args.target_missing_value = "-999"
    result = run_quantile_h2_from_args(workflow.args)
    assert result.quantile_h2.n_snps.tolist() == [2, 1]
    np.testing.assert_allclose(result.quantile_h2.h2_obs, [.2, .25])
    np.testing.assert_allclose(result.standardized_coefficients.annotation_sd, [0, np.std([0, 1, 2, 3])])
    np.testing.assert_allclose(result.standardized_coefficients.tau_star, [0, .1 * np.std([0, 1, 2, 3]) * 4 / .8])


def test_high_offset_annotation_variance_is_stable(workflow):
    path = workflow.root / "annotations.tsv"
    frame = pd.read_csv(path, sep="\t")
    frame["cont"] = np.array([1e6, 1e6 + .125, 1e6 + .25, 1e6 + .375], dtype=np.float32)
    frame.to_csv(path, sep="\t", index=False)
    metadata_path = workflow.root / "ldscores/metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["counts"][1]["common_reference_snp_count"] = float(frame.cont.sum())
    metadata["files"] = {}
    metadata_path.write_text(json.dumps(metadata))
    result = run_quantile_h2_from_args(workflow.args)
    np.testing.assert_allclose(result.standardized_coefficients.annotation_sd, [0, np.std(frame.cont.to_numpy(dtype=np.float64))], rtol=1e-12)


def test_bounded_reads_and_released_chromosome_arrays(workflow, monkeypatch):
    from ldsc import _quantile_storage
    from ldsc._annotation_bundle import AnnotationBundle

    monkeypatch.setattr(_quantile_storage, "ROW_CHUNK_SIZE", 1)
    original = AnnotationBundle.read
    previous, calls = [], []

    def read(bundle, chrom, *, rows=None, columns=None):
        assert all(value() is None for value in previous)
        assert len(rows) <= 1
        assert columns == ["base", "cont"]
        values = original(bundle, chrom, rows=rows, columns=columns)
        previous.append(weakref.ref(values))
        calls.append(chrom)
        return values

    monkeypatch.setattr(AnnotationBundle, "read", read)
    result = run_quantile_h2_from_args(workflow.args)
    assert calls == ["1", "1", "2", "2"] * 2
    assert all(value() is None for value in previous)
    np.testing.assert_allclose(result.quantile_h2.h2_obs, [.45, .35])


def test_omitted_alleles_require_globally_unique_reference_bases(workflow):
    metadata_path = workflow.root / "ldscores/metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["snp_identifier"] = "rsid_allele_aware"
    metadata_path.write_text(json.dumps(metadata))
    path = workflow.root / "ref.tsv"
    frame = pd.read_csv(path, sep="\t")
    frame["A1"], frame["A2"] = "A", "C"
    frame.loc[1, ["SNP", "A2"]] = ["rs1", "G"]
    frame.to_csv(path, sep="\t", index=False)
    with pytest.raises(LDSCInputError, match="cannot infer omitted alleles.*non-unique"):
        run_quantile_h2_from_args(workflow.args)
    issues = pd.read_csv(Path(workflow.args.output_dir) / "diagnostics/snp_alignment_issues.tsv.gz", sep="\t")
    assert issues.issue.tolist() == ["ambiguous_allele_inference"] * 2
    assert issues.CHR.tolist() == [1, 2]


def test_global_ties_and_sharded_baselines(workflow):
    annotation = pd.read_csv(workflow.root / "annotations.tsv", sep="\t")
    paths = []
    for chrom, frame in annotation.groupby("CHR"):
        path = workflow.root / f"baseline.{chrom}.annot.gz"
        frame.to_csv(path, sep="\t", index=False)
        paths.append(str(path))
    workflow.args.baseline_annot_sources = paths
    path = workflow.root / "target.tsv"
    frame = pd.read_csv(path, sep="\t")
    frame["target"] = [0, 1, 1, 2]
    frame.to_csv(path, sep="\t", index=False)
    result = run_quantile_h2_from_args(workflow.args)
    assert result.quantile_h2.n_snps.tolist() == [3, 1]
    np.testing.assert_allclose(result.quantile_h2.target_value_upper, [1, 2])
    np.testing.assert_allclose(result.quantile_h2.h2_obs, [.45, .35])


@pytest.mark.parametrize("route", ["annot", "bed", "gene"])
def test_single_query_model_reconstructs_each_query_route(workflow, route):
    path = workflow.root / "annotations.tsv"
    frame = pd.read_csv(path, sep="\t")
    frame.drop(columns="cont").to_csv(path, sep="\t", index=False)
    metadata_path = workflow.root / "ldscores/metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata.update(snp_identifier="chr_pos", genome_build="hg38", baseline_columns=["base"], query_columns=["cont"])
    metadata["counts"][1]["common_reference_snp_count"] = 2.0
    metadata_path.write_text(json.dumps(metadata))
    overlap_path = workflow.root / "ldscores/ldscore.overlap.parquet"
    overlap = pd.read_parquet(overlap_path)
    overlap["overlap_all_snps"] = overlap["overlap_common_snps"] = [4., 2., 2., 2.]
    overlap = overlap.loc[~(overlap.row_annotation.eq("cont") & overlap.col_annotation.eq("base"))]
    overlap.to_parquet(overlap_path, index=False)
    model_dir = workflow.root / "single_query"
    model_dir.mkdir()
    fitted = workflow.root / "fitted"
    (model_dir / "metadata.json").write_text((fitted / "diagnostics/metadata.json").read_text())
    (model_dir / "partitioned_h2_full.tsv").write_bytes((fitted / "partitioned_h2.tsv").read_bytes())
    (model_dir / "coefficient_delete_values.parquet").write_bytes((fitted / "diagnostics/coefficient_delete_values.parquet").read_bytes())
    workflow.args.partitioned_h2_result_dir = str(model_dir)
    if route == "annot":
        query = frame.drop(columns="base")
        query["cont"] = [0., 1., 1., 0.]
        path = workflow.root / "query.annot.gz"
        query.to_csv(path, sep="\t", index=False)
        workflow.args.query_annot_sources = [str(path)]
    elif route == "bed":
        path = workflow.root / "cont.bed"
        path.write_text("1\t14\t35\n2\t14\t35\n")
        workflow.args.query_annot_bed_sources = [str(path)]
    else:
        catalog = workflow.root / "genes.tsv"
        catalog.write_text("gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\nG1\tONE\t1\t15\t35\thg38\nG2\tTWO\t2\t15\t35\thg38\n")
        path = workflow.root / "cont.txt"
        path.write_text("G1\nG2\n")
        workflow.args.query_annot_gene_list_sources = [str(path)]
        workflow.args.gene_coordinate_file = str(catalog)
    result = run_quantile_h2_from_args(workflow.args)
    assert result.metadata["selected_model_type"] == "baseline_plus_query"
    np.testing.assert_allclose(result.quantile_h2.h2_obs, [.35, .05])
    assert not list(Path(workflow.args.output_dir).glob(".ldsc-annotation-*"))
