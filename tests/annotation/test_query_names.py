"""Query names that cannot become annotation columns stop every preparation path."""

import pandas as pd
import pytest

from ldsc import AnnotationBuildConfig, AnnotationBuilder, GlobalConfig, run_annotate, run_ldscore, set_global_config
from ldsc.errors import LDSCInputError
from tests.ldscore.test_ldscore_parallelism import _PANEL_CHROMS, _write_index_panel


@pytest.fixture
def grid(tmp_path):
    metadata = pd.concat([pair[0] for pair in _PANEL_CHROMS.values()], ignore_index=True)
    baseline = tmp_path / "baseline.annot"
    metadata.loc[:, ["CHR", "POS", "SNP", "A1", "A2"]].assign(base=1).to_csv(baseline, sep="\t", index=False)
    bim = tmp_path / "grid.bim"
    metadata.assign(CM=0)[["CHR", "SNP", "CM", "POS", "A1", "A2"]].to_csv(bim, sep="\t", header=False, index=False)
    catalog = tmp_path / "genes.tsv"
    catalog.write_text("gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\nG1\tONE\t1\t100\t300\thg19\n")
    return {"baseline": baseline, "bim": bim, "catalog": catalog, "metadata": metadata}


def write_bed(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("chr1\t150\t350\n")
    return path


def collision_issues(output):
    issues = pd.read_csv(output / "diagnostics/input_issues.tsv", sep="\t")
    assert issues.reason.eq("annotation_name_collision").all()
    return issues


def refuse_grid_preparation(monkeypatch, *modules_and_names):
    def unexpected(*args, **kwargs):
        pytest.fail("A query-name failure must stop before BED parsing or grid preparation")
    for module, name in modules_and_names:
        monkeypatch.setattr(module, name, unexpected)


@pytest.mark.parametrize("kind", ["annot", "bim"])
def test_annotate_rejects_every_unusable_bed_name_before_the_grid(tmp_path, grid, monkeypatch, kind):
    from ldsc._annotation import annotate as workflow

    beds = [write_bed(tmp_path / name) for name in ("SNP.bed", "gene_control.bed", "a/query.bed", "b/query.bed", "ok.bed")]
    if kind == "annot":
        beds.append(write_bed(tmp_path / "base.bed"))
    refuse_grid_preparation(monkeypatch, (workflow, "prepare_bed_queries"), (workflow, "prepare_annotation_sources"),
                            (workflow, "prepare_bim_grid"))
    grid_option = {"ref_snp_annot_sources": [grid["baseline"]]} if kind == "annot" else {"ref_snp_bim_sources": [grid["bim"]]}
    output = tmp_path / "out"
    with pytest.raises(LDSCInputError, match="annotation_name_collision|query name") as failure:
        run_annotate(**grid_option, query_annot_bed_sources=beds, genome_build="hg19",
                     global_config=GlobalConfig(snp_identifier="rsid"), output_dir=output)
    message = str(failure.value)
    rejected = [path for path in beds if path.name != "ok.bed"]
    for path in rejected:
        assert str(path) in message
    assert "troubleshooting.md#common-a-bed-query-name-is-not-usable" in message
    assert sorted(collision_issues(output).source) == sorted(map(str, rejected))
    assert not list(output.glob("*.annot.gz"))


def test_bed_names_are_checked_before_bed_contents(tmp_path, grid):
    malformed = tmp_path / "malformed.bed"
    malformed.write_text("chr1\t300\t200\n")
    output = tmp_path / "out"
    with pytest.raises(LDSCInputError) as failure:
        run_annotate(ref_snp_annot_sources=[grid["baseline"]], query_annot_bed_sources=[malformed, write_bed(tmp_path / "SNP.bed")],
                     global_config=GlobalConfig(snp_identifier="rsid"), output_dir=output)
    assert "SNP.bed" in str(failure.value)
    assert "malformed.bed" not in str(failure.value)


def test_ldscore_rejects_reserved_bed_names_before_the_grid(tmp_path, grid, monkeypatch):
    from ldsc._ldscore import direct_annotation

    _write_index_panel(tmp_path / "panel", _PANEL_CHROMS)
    regression = tmp_path / "regression.tsv"
    grid["metadata"][["SNP"]].to_csv(regression, sep="\t", index=False)
    refuse_grid_preparation(monkeypatch, (direct_annotation, "prepare_bed_queries"),
                            (direct_annotation, "prepare_annotation_sources"))
    set_global_config(GlobalConfig(snp_identifier="rsid"))
    output = tmp_path / "ld"
    with pytest.raises(LDSCInputError, match="SNP.bed"):
        run_ldscore(baseline_annot_sources=[str(grid["baseline"])],
                    query_annot_bed_sources=[str(write_bed(tmp_path / "SNP.bed")), str(write_bed(tmp_path / "base.bed"))],
                    r2_matrix_dir=str(tmp_path / "panel"), ld_wind_kb=1, yes_really=True, regr_snps_file=str(regression),
                    regr_snps_exclude_regions="none", output_dir=str(output))
    assert sorted(collision_issues(output).source) == sorted(str(tmp_path / name) for name in ("SNP.bed", "base.bed"))


def test_ldscore_gene_list_name_error_explains_the_rename(tmp_path, grid):
    _write_index_panel(tmp_path / "panel", _PANEL_CHROMS)
    genes = tmp_path / "SNP.txt"
    genes.write_text("G1\n")
    set_global_config(GlobalConfig(snp_identifier="rsid"))
    with pytest.raises(LDSCInputError) as failure:
        run_ldscore(baseline_annot_sources=[str(grid["baseline"])], query_annot_gene_list_sources=[str(genes)],
                    gene_coordinate_file=str(grid["catalog"]), padding_bp=0,
                    r2_matrix_dir=str(tmp_path / "panel"), ld_wind_kb=1, yes_really=True,
                    regr_snps_exclude_regions="none", output_dir=str(tmp_path / "ld"))
    assert "annotation_name_collision" in str(failure.value)
    assert "Rename conflicting query source files" in str(failure.value)


@pytest.mark.parametrize("name", ["SNP.bed", "base.bed", "gene_control.bed"])
def test_annotation_builder_rejects_unusable_bed_names(tmp_path, grid, name):
    spec = AnnotationBuildConfig(baseline_annot_sources=[str(grid["baseline"])],
                                 query_annot_bed_sources=[str(write_bed(tmp_path / name))])
    output = tmp_path / "builder"
    with pytest.raises(LDSCInputError, match=name):
        AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), projection_genome_build="hg19").run(spec, output_dir=output)
    assert collision_issues(output).source.tolist() == [str(tmp_path / name)]


@pytest.mark.parametrize("name", ["SNP", "base"])
def test_annotation_builder_rejects_unusable_gene_list_names(tmp_path, grid, name):
    genes = tmp_path / f"{name}.txt"
    genes.write_text("G1\n")
    spec = AnnotationBuildConfig(baseline_annot_sources=[str(grid["baseline"])], query_annot_gene_list_sources=[str(genes)],
                                 gene_coordinate_file=str(grid["catalog"]), padding_bp=0)
    output = tmp_path / "builder"
    with pytest.raises(LDSCInputError, match="annotation_name_collision") as failure:
        AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), projection_genome_build="hg19").run(spec, output_dir=output)
    assert "Rename conflicting query source files" in str(failure.value)
    summary = pd.read_csv(output / "diagnostics/gene_list_resolution_summary.tsv", sep="\t")
    assert summary.source_reasons.tolist() == ["annotation_name_collision"]
