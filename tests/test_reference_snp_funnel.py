"""Reference-SNP funnels recorded by reference preparation on both panel backends.

Expected counts come from worked examples over the small fixtures: the
two-SNP R2 panel from ``tests/test_reference_preparation.py`` and the 32-SNP
chromosome-22 PLINK subset.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ldsc import GlobalConfig, LDScoreConfig, RefPanelConfig
from ldsc._refpanel.adapters import ParquetR2RefPanel, PlinkRefPanel
from ldsc._refpanel.funnel import annotation_snp_grid_emptied, with_grid_rows
from ldsc._refpanel.prepared import ChromosomeAnnotations
from ldsc.errors import EmptyReferenceSNPs
from tests.test_reference_preparation import write_panel

PLINK_PREFIX = Path(__file__).resolve().parent / "fixtures" / "minimal_external_resources" / "plink" / "hm3_chr22_subset"
WINDOW = LDScoreConfig(ld_wind_snps=1, whole_chromosome_ok=True)


def steps(frame, funnel):
    """Return ``{step: (removed, remaining)}`` with missing counts as ``None``."""
    rows = frame.loc[frame["funnel"].eq(funnel)].sort_values("step_index")
    return {row.step: (None if pd.isna(row.removed) else int(row.removed),
                       None if pd.isna(row.remaining) else int(row.remaining)) for row in rows.itertuples()}


def assert_invariant(frame, retained):
    """Each funnel's start minus its recorded removals equals the retained reference SNPs."""
    for funnel in ("reference_panel", "annotation_snp_grid"):
        rows = frame.loc[frame["funnel"].eq(funnel)].sort_values("step_index")
        start = int(rows["remaining"].iloc[0])
        assert start - int(rows["removed"].fillna(0).sum()) == int(rows["remaining"].iloc[-1]) == retained


def r2_annotations(snps):
    positions = {"rs1": 100, "rs2": 200, "rs3": 300, "absent": 400}
    metadata = pd.DataFrame({"CHR": "1", "SNP": snps, "POS": [positions[snp] for snp in snps]})
    return ChromosomeAnnotations(metadata, pd.DataFrame({"base": np.ones(len(snps))}), ["base"], [])


def test_r2_funnel_counts_grid_rows_missing_from_the_panel(tmp_path):
    write_panel(tmp_path)
    panel = ParquetR2RefPanel(GlobalConfig(snp_identifier="rsid"), RefPanelConfig(backend="parquet_r2", r2_dir=tmp_path))
    with panel.prepare_chromosome("1", r2_annotations(["rs1", "rs2", "rs3"]), WINDOW) as prepared:
        funnel = prepared.reference_snp_funnel
        retained = len(prepared.metadata)
    assert steps(funnel, "reference_panel") == {
        "panel_rows": (None, 2), "ref_panel_snps_file": (None, 2), "maf_floor": (None, 2),
        "allele_aware_identity_cleanup": (None, 2), "annotation_snp_grid": (0, 2),
    }
    assert steps(funnel, "annotation_snp_grid") == {"grid_rows": (None, 3), "annotation_snp_grid": (1, 2)}
    assert set(funnel["backend"]) == {"parquet_r2"} and set(funnel["chrom"]) == {"1"}
    assert_invariant(funnel, retained)


def test_r2_funnel_applies_maf_before_the_grid_and_counts_off_grid_panel_snps(tmp_path):
    write_panel(tmp_path)  # rs1 has MAF 0.04 and rs3 has MAF 0.3
    panel = ParquetR2RefPanel(GlobalConfig(snp_identifier="rsid"),
                              RefPanelConfig(backend="parquet_r2", r2_dir=tmp_path, maf_min=0.01))
    with panel.prepare_chromosome("1", r2_annotations(["rs1"]), WINDOW) as prepared:
        funnel = prepared.reference_snp_funnel
        assert prepared.metadata["SNP"].tolist() == ["rs1"]
    assert steps(funnel, "reference_panel")["maf_floor"] == (0, 2)
    assert steps(funnel, "reference_panel")["annotation_snp_grid"] == (1, 1)
    assert steps(funnel, "annotation_snp_grid") == {"grid_rows": (None, 1), "annotation_snp_grid": (0, 1)}
    assert_invariant(funnel, 1)


def test_r2_empty_chromosome_carries_its_funnel_and_names_the_emptying_step(tmp_path):
    write_panel(tmp_path)
    panel = ParquetR2RefPanel(GlobalConfig(snp_identifier="rsid"), RefPanelConfig(backend="parquet_r2", r2_dir=tmp_path))
    with pytest.raises(EmptyReferenceSNPs) as caught:
        panel.prepare_chromosome("1", r2_annotations(["absent"]), WINDOW)
    funnel = caught.value.reference_snp_funnel
    assert steps(funnel, "reference_panel")["annotation_snp_grid"] == (2, 0)
    assert annotation_snp_grid_emptied(funnel)

    restriction = tmp_path / "keep.tsv"
    restriction.write_text("SNP\nabsent\n")
    restricted = ParquetR2RefPanel(GlobalConfig(snp_identifier="rsid"), RefPanelConfig(
        backend="parquet_r2", r2_dir=tmp_path, ref_panel_snps_file=str(restriction)))
    with pytest.raises(EmptyReferenceSNPs) as caught:
        restricted.prepare_chromosome("1", r2_annotations(["rs1"]), WINDOW)
    funnel = caught.value.reference_snp_funnel
    assert steps(funnel, "reference_panel")["ref_panel_snps_file"] == (2, 0)
    assert not annotation_snp_grid_emptied(funnel)


def plink_annotations(rows, *, alleles):
    bim = pd.read_csv(PLINK_PREFIX.with_suffix(".bim"), sep=r"\s+", header=None,
                      names=["CHR", "SNP", "CM", "POS", "A1", "A2"]).iloc[rows]
    columns = ["CHR", "SNP", "POS", *(["A1", "A2"] if alleles else [])]
    metadata = bim.loc[:, columns].astype({"CHR": str}).reset_index(drop=True)
    return ChromosomeAnnotations(metadata, pd.DataFrame({"base": np.ones(len(metadata))}), ["base"], [])


@pytest.mark.parametrize(("mode", "alleles"), [("rsid", False), ("chr_pos_allele_aware", True)])
def test_plink_funnel_counts_off_grid_panel_snps_before_genotype_qc(mode, alleles):
    pytest.importorskip("bitarray")
    config = GlobalConfig(snp_identifier=mode, **({} if mode == "rsid" else {"genome_build": "hg19"}))
    panel = PlinkRefPanel(config, RefPanelConfig(backend="plink", plink_prefix=str(PLINK_PREFIX)))
    with panel.prepare_chromosome("22", plink_annotations(slice(0, 20), alleles=alleles), WINDOW) as prepared:
        funnel = prepared.reference_snp_funnel
        retained, qc_removed = len(prepared.metadata), prepared.genotype_qc_removed
    reference = steps(funnel, "reference_panel")
    assert reference["panel_rows"] == (None, 32)
    assert reference["identity_cleanup"][0] == (0 if mode == "rsid" else None)
    assert reference["ref_panel_snps_file"] == (None, 32)
    assert reference["allele_aware_identity_cleanup"][0] == (None if mode == "rsid" else 0)
    assert reference["annotation_snp_grid"] == (12, 20)
    assert reference["genotype_qc"] == (qc_removed, 20 - qc_removed)
    assert reference["maf_floor"] == (None, retained)
    assert steps(funnel, "annotation_snp_grid")["grid_rows"] == (None, 20)
    assert_invariant(funnel, retained)


def test_plink_grid_that_misses_every_panel_snp_is_emptied_at_the_grid_step():
    pytest.importorskip("bitarray")
    annotations = plink_annotations(slice(0, 2), alleles=False)
    annotations.metadata["SNP"] = ["absent1", "absent2"]
    panel = PlinkRefPanel(GlobalConfig(snp_identifier="rsid"), RefPanelConfig(backend="plink", plink_prefix=str(PLINK_PREFIX)))
    with pytest.raises(EmptyReferenceSNPs) as caught:
        panel.prepare_chromosome("22", annotations, WINDOW)
    funnel = caught.value.reference_snp_funnel
    assert steps(funnel, "reference_panel")["annotation_snp_grid"] == (32, 0)
    assert steps(funnel, "reference_panel")["genotype_qc"] == (None, 0)
    assert steps(funnel, "annotation_snp_grid") == {
        "grid_rows": (None, 2), "annotation_snp_grid": (2, 0), "genotype_qc": (None, 0), "maf_floor": (None, 0)}
    assert annotation_snp_grid_emptied(funnel)


def test_index_grid_funnel_restarts_before_the_baseline_intersection(tmp_path):
    write_panel(tmp_path)
    panel = ParquetR2RefPanel(GlobalConfig(snp_identifier="rsid"), RefPanelConfig(backend="parquet_r2", r2_dir=tmp_path))
    with panel.prepare_chromosome("1", r2_annotations(["rs1", "rs3"]), WINDOW) as prepared:
        funnel = with_grid_rows(prepared.reference_snp_funnel, 5)
    assert steps(funnel, "annotation_snp_grid") == {"grid_rows": (None, 5), "annotation_snp_grid": (3, 2)}
    assert_invariant(funnel, 2)


# --- LD-score workflow -------------------------------------------------------

def r2_workflow_inputs(root, baseline_snps):
    """Write the three-SNP chromosome-22 R2 panel and a baseline on ``baseline_snps``."""
    from tests.test_ldscore_chromosome_scope import write_r2_inputs

    inputs = write_r2_inputs(root)
    positions = {"rs1": 10, "rs2": 20, "rs3": 30, "rs4": 40}
    pd.DataFrame({"CHR": "22", "SNP": baseline_snps, "POS": [positions[snp] for snp in baseline_snps],
                  "base": 1}).to_csv(inputs.baseline_annot_sources, sep="\t", index=False)
    (root / "regression.txt").write_text("SNP\nrs1\nrs2\nrs3\n")
    return inputs


def workflow_argv(root, inputs, *extra):
    return ["--output-dir", str(root / "out"), "--r2-dir", inputs.r2_dir,
            "--baseline-annot-sources", inputs.baseline_annot_sources,
            "--regr-snps-file", str(root / "regression.txt"), "--regr-snps-exclude-regions", "none",
            "--snp-identifier", "rsid", "--ld-wind-snps", "2", "--yes-really", *extra]


def read_funnel(root):
    from ldsc._refpanel.funnel import read_reference_snp_funnel

    return read_reference_snp_funnel(root / "out" / "diagnostics" / "reference_snp_funnel.tsv")


def test_panel_larger_than_the_grid_succeeds_and_reports_the_excluded_panel_snps(tmp_path, capsys):
    import json
    from ldsc._ldscore.calculator import main

    inputs = r2_workflow_inputs(tmp_path, ["rs1", "rs2"])
    result = main(workflow_argv(tmp_path, inputs))
    assert result.baseline_table["SNP"].tolist() == ["rs1", "rs2"]
    funnel = read_funnel(tmp_path)
    assert steps(funnel, "reference_panel")["annotation_snp_grid"] == (1, 2)
    assert_invariant(funnel, 2)
    log = (tmp_path / "out" / "diagnostics" / "ldscore.log").read_text()
    warning = "ldscore excluded 1 reference-panel SNPs on 1 chromosome that are not on the annotation SNP grid"
    assert any("WARNING" in line and warning in line for line in log.splitlines())
    assert ["reference_panel_source", inputs.r2_dir] in [line.split() for line in log.splitlines()]
    assert f"WARNING: {warning}" in capsys.readouterr().err
    metadata = json.loads((tmp_path / "out" / "metadata.json").read_text())
    universe = metadata["snp_universe_policy"]["ld_reference_universe"]
    assert "selection" not in universe and universe["retained_snp_count"] == 2
    assert "reference_snp_funnel" not in metadata["files"]


def test_grid_larger_than_the_panel_succeeds_without_a_warning(tmp_path, capsys):
    from ldsc._ldscore.calculator import main

    inputs = r2_workflow_inputs(tmp_path, ["rs1", "rs2", "rs3", "rs4"])
    main(workflow_argv(tmp_path, inputs))
    funnel = read_funnel(tmp_path)
    assert steps(funnel, "reference_panel")["annotation_snp_grid"] == (0, 3)
    assert steps(funnel, "annotation_snp_grid") == {"grid_rows": (None, 4), "annotation_snp_grid": (1, 3)}
    assert "excluded" not in (tmp_path / "out" / "diagnostics" / "ldscore.log").read_text()
    assert "excluded" not in capsys.readouterr().err


def test_chromosome_emptied_by_the_panel_restriction_is_skipped_with_a_logged_warning(tmp_path):
    from ldsc._ldscore.calculator import build_parser, run_ldscore_from_args

    inputs = r2_workflow_inputs(tmp_path, ["rs1", "rs2", "rs3"])
    write_panel(Path(inputs.r2_dir), chrom="1")  # chromosome 1 holds rs1 and rs3
    baseline = pd.read_csv(inputs.baseline_annot_sources, sep="\t")
    chrom1 = pd.DataFrame({"CHR": [1], "SNP": ["absent"], "POS": [100], "base": [1]})
    pd.concat([chrom1, baseline], ignore_index=True).to_csv(inputs.baseline_annot_sources, sep="\t", index=False)
    restriction = tmp_path / "keep.txt"
    restriction.write_text("SNP\nrs2\n")
    result = run_ldscore_from_args(build_parser().parse_args(
        workflow_argv(tmp_path, inputs, "--ref-panel-snps-file", str(restriction))))
    assert result.baseline_table["CHR"].tolist() == ["22"]
    funnel = read_funnel(tmp_path)
    chrom1 = funnel.loc[funnel["chrom"].eq("1")]
    assert steps(chrom1, "reference_panel")["ref_panel_snps_file"] == (2, 0)
    assert steps(chrom1, "reference_panel")["annotation_snp_grid"] == (0, 0)
    log = (tmp_path / "out" / "diagnostics" / "ldscore.log").read_text()
    assert any("WARNING" in line and "Skipping chromosome 1" in line for line in log.splitlines())


def test_multi_batch_run_records_each_chromosome_funnel_once(tmp_path):
    from ldsc._ldscore.calculator import build_parser, run_ldscore_from_args

    inputs = r2_workflow_inputs(tmp_path, ["rs1", "rs2"])
    query = tmp_path / "query.22.annot"
    pd.DataFrame({"CHR": "22", "SNP": ["rs1", "rs2"], "POS": [10, 20], "q1": [1, 0], "q2": [0, 1]}).to_csv(
        query, sep="\t", index=False)
    run_ldscore_from_args(build_parser().parse_args(workflow_argv(
        tmp_path, inputs, "--query-annot-sources", str(query), "--query-batch-size", "1")))
    funnel = read_funnel(tmp_path)
    assert len(funnel) == 7
    assert not funnel.duplicated(["funnel", "chrom", "step"]).any()
    assert steps(funnel, "reference_panel")["annotation_snp_grid"] == (1, 2)


def test_in_memory_result_carries_the_funnel(tmp_path):
    from ldsc import LDScoreCalculator
    from tests.test_reference_preparation import public_bundle

    metadata = write_panel(tmp_path)
    config = GlobalConfig(snp_identifier="rsid")
    panel = ParquetR2RefPanel(config, RefPanelConfig(backend="parquet_r2", r2_dir=tmp_path))
    result = LDScoreCalculator().run(public_bundle(metadata.iloc[[0]]), panel, WINDOW, config)
    assert steps(result.reference_snp_funnel, "reference_panel")["annotation_snp_grid"] == (1, 1)
    assert list(result.reference_snp_funnel.columns) == [
        "funnel", "chrom", "backend", "step_index", "step", "removed", "remaining"]


def test_gene_list_run_fails_in_gate_b_when_the_grid_misses_a_whole_chromosome(tmp_path):
    from ldsc._ldscore.calculator import build_parser, run_ldscore_from_args
    from ldsc.errors import LDSCInputError

    inputs = r2_workflow_inputs(tmp_path, ["rs1", "rs2"])
    baseline = pd.read_csv(inputs.baseline_annot_sources, sep="\t").assign(SNP=["absent1", "absent2"])
    baseline.to_csv(inputs.baseline_annot_sources, sep="\t", index=False)
    genes, catalog = tmp_path / "focal.txt", tmp_path / "genes.tsv"
    genes.write_text("G1\n")
    catalog.write_text("gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\nG1\tONE\t22\t10\t10\thg19\n")
    argv = workflow_argv(tmp_path, inputs, "--query-annot-gene-list-sources", str(genes),
                         "--gene-coordinate-file", str(catalog), "--genome-build", "hg19", "--padding-bp", "0")
    with pytest.raises(LDSCInputError, match=r"annotation SNP grid for chromosome 22 \(2 annotation rows, "
                                             r"3 reference-panel SNPs\)"):
        run_ldscore_from_args(build_parser().parse_args(argv))
    assert not (tmp_path / "out" / "ldscore.baseline.parquet").exists()
    # Gate B writes query diagnostics before failing; batch-1 calculation never starts.
    assert (tmp_path / "out" / "diagnostics" / "query_annotation_status.tsv").exists()
    assert not (tmp_path / "out" / "diagnostics" / "chromosome-logs").exists()
