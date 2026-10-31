"""Every file a run touches stays inside its output directory.

An empty output directory would otherwise normalize to the working directory,
where an authorized overwrite can replace or remove the user's own files.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

import ldsc
from ldsc import cli
from ldsc._annotation.storage import AnnotationWorkspace
from ldsc._runtime.paths import ensure_output_directory
from ldsc._runtime.workflow_logging import overwrite_failure_marker
from ldsc.errors import LDSCConfigError, LDSCInputError
from tests.runtime.test_failure_markers import FAILING_CLI_RUNS

PLINK_FIXTURE = Path(__file__).parents[1] / "fixtures/minimal_external_resources/plink/hm3_chr22_subset"


@pytest.fixture
def empty_cwd(tmp_path, monkeypatch):
    """Run from an empty working directory, with inputs kept in a separate folder."""
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    return cwd


@pytest.mark.parametrize("command", FAILING_CLI_RUNS)
@pytest.mark.parametrize("blank", ["", "  "])
def test_cli_rejects_an_empty_destination_before_writing(empty_cwd, command, blank, capsys):
    options, _ = FAILING_CLI_RUNS[command]

    assert cli.run_cli([command, *(token.format(root=blank) for token in options), "--overwrite"]) == 2

    assert list(empty_cwd.iterdir()) == []
    assert "must name a directory" in capsys.readouterr().err


def _legacy_suite(root: Path) -> dict[str, str]:
    for role in ("reference", "weights"):
        (root / role).mkdir(parents=True)
    for chrom in range(1, 23):
        with gzip.open(root / "reference" / f"{chrom}.l2.ldscore.gz", "wt") as handle:
            handle.write(f"CHR SNP BP L2\n{chrom} rs{chrom} {chrom * 100} {1.0 + chrom / 100}\n")
        with gzip.open(root / "weights" / f"weights.{chrom}.l2.ldscore.gz", "wt") as handle:
            handle.write(f"CHR SNP BP L2\n{chrom} rs{chrom} {chrom * 100 + 1} {2.0 + chrom / 100}\n")
        (root / "reference" / f"{chrom}.l2.M_5_50").write_text("10\n")
    return {"legacy_reference_dir": str(root / "reference"), "legacy_weight_dir": str(root / "weights")}


def _python_call(name: str, inputs: Path, blank: str):
    """Return a Python call that would succeed with valid inputs if the empty destination were accepted."""
    inputs.mkdir(exist_ok=True)
    bed = inputs / "query.bed"
    bed.write_text("chr22\t16000000\t17000000\n")
    raw = inputs / "raw.tsv"
    raw.write_text("SNP A1 A2 P BETA N\nrs1 A G 0.05 0.2 1000\nrs2 C T 0.5 -0.1 1000\n")
    calls = {
        "run_annotate": lambda: ldsc.run_annotate(
            output_dir=blank, ref_snp_bim_sources=[f"{PLINK_FIXTURE}.bim"], query_annot_bed_sources=[str(bed)],
            genome_build="hg38",
        ),
        "convert_legacy_ldscores": lambda: ldsc.convert_legacy_ldscores(output_dir=blank, **_legacy_suite(inputs)),
        "SumstatsMunger.run": lambda: ldsc.SumstatsMunger().run(
            ldsc.MungeConfig(raw_sumstats_file=str(raw), output_dir=blank, no_snp_restriction=True),
            global_config=ldsc.GlobalConfig(snp_identifier="rsid"),
        ),
        "AnnotationBuilder.run": lambda: ldsc.AnnotationBuilder(
            ldsc.AnnotationBuildConfig(query_annot_bed_sources=(str(bed),)),
        ).run(output_dir=blank),
        "plot_result": lambda: ldsc.plot_result(blank, overwrite=True),
        "plot_result output_dir": lambda: ldsc.plot_result(str(inputs), output_dir=blank, overwrite=True),
        "convert_h2_scale": lambda: ldsc.convert_h2_scale(blank, samp_prev=0.5, pop_prev=0.01, overwrite=True),
        "convert_h2_scale output_dir": lambda: ldsc.convert_h2_scale(
            str(inputs), samp_prev=0.5, pop_prev=0.01, output_dir=blank, overwrite=True,
        ),
    }
    return calls[name]


@pytest.mark.parametrize("name", [
    "run_annotate", "convert_legacy_ldscores", "SumstatsMunger.run", "AnnotationBuilder.run",
    "plot_result", "plot_result output_dir", "convert_h2_scale", "convert_h2_scale output_dir",
])
def test_python_rejects_an_empty_destination_before_writing(empty_cwd, tmp_path, name):
    call = _python_call(name, tmp_path / "inputs", "")

    with pytest.raises(LDSCConfigError, match="empty"):
        call()

    assert list(empty_cwd.iterdir()) == []


@pytest.mark.parametrize("helper", [
    lambda: ensure_output_directory(""),
    lambda: AnnotationWorkspace(" "),
    lambda: overwrite_failure_marker("", overwrite=True, command="unit").__enter__(),
])
def test_output_helpers_reject_an_empty_destination(empty_cwd, helper):
    with pytest.raises(LDSCConfigError, match="empty"):
        helper()

    assert list(empty_cwd.iterdir()) == []


@pytest.mark.parametrize("command,option", [("plot", "--result-dir"), ("convert-h2-scale", "--h2-result-dir")])
@pytest.mark.parametrize("overwrite", [True, False])
def test_derived_commands_create_nothing_for_a_missing_result_directory(empty_cwd, command, option, overwrite):
    extra = ["--samp-prev", "0.5", "--pop-prev", "0.01"] if command == "convert-h2-scale" else []

    assert cli.run_cli([command, option, "resutls_typo", *extra, *(["--overwrite"] if overwrite else [])]) == 1

    assert list(empty_cwd.iterdir()) == []


def _tree(root: Path) -> set[str]:
    return {path.relative_to(root).as_posix() for path in root.rglob("*")}


@pytest.fixture
def index_inputs(tmp_path):
    from tests.ldscore.test_plink_workflow_resolution import write_inputs

    inputs = tmp_path / "inputs"
    inputs.mkdir()
    write_inputs(inputs)
    return inputs


def _build_index(inputs, index, monkeypatch, *, prefix="1000G.EUR.QC.", overwrite=False):
    from ldsc._ldscore.gene_index import build as gene_index_build
    from tests.ldscore.test_plink_workflow_resolution import index_args

    args = index_args(inputs, str(inputs / prefix), index, monkeypatch)
    args.overwrite = overwrite
    return gene_index_build.run_build_gene_ldscore_index_from_args(args)


def test_gene_index_build_writes_only_inside_its_output_directory(tmp_path, index_inputs, empty_cwd, monkeypatch):
    from ldsc._ldscore.gene_index.scoring import load_gene_ldscore_index

    index = tmp_path / "parent" / "index"
    inputs_before = _tree(index_inputs)

    _build_index(index_inputs, index, monkeypatch)
    _build_index(index_inputs, index, monkeypatch, overwrite=True)

    assert [path.name for path in index.parent.iterdir()] == ["index"]
    assert not (index / ".build").exists()
    assert (index / "diagnostics" / "build-gene-ldscore-index.log").is_file()
    assert _tree(index_inputs) == inputs_before
    assert list(empty_cwd.iterdir()) == []
    load_gene_ldscore_index(index)


@pytest.mark.parametrize("overwrite", [False, True])
def test_failed_gene_index_build_keeps_its_state_inside_and_retries(tmp_path, index_inputs, monkeypatch, overwrite):
    index = tmp_path / "parent" / "index"

    with pytest.raises(LDSCInputError):
        _build_index(index_inputs, index, monkeypatch, prefix="missing-panel", overwrite=overwrite)

    assert [path.name for path in index.parent.iterdir()] == ["index"]
    assert {"build-gene-ldscore-index.log", "plink_input_issues.tsv"} <= {
        path.name for path in (index / ".build").iterdir()}
    assert (index / "RUN_FAILED.txt").is_file() == overwrite

    _build_index(index_inputs, index, monkeypatch)

    assert [path.name for path in index.parent.iterdir()] == ["index"]
    assert not (index / ".build").exists()
    assert not (index / "RUN_FAILED.txt").exists()


def test_an_incomplete_gene_index_is_never_read_and_needs_overwrite(tmp_path, index_inputs, monkeypatch):
    from ldsc._ldscore.gene_index.scoring import load_gene_ldscore_index

    index = tmp_path / "parent" / "index"
    _build_index(index_inputs, index, monkeypatch)
    # An interrupted publication leaves owned components without the root metadata written last,
    # beside the build directory that only a successful build removes.
    (index / "metadata.json").unlink()
    (index / ".build").mkdir()

    with pytest.raises(LDSCInputError):
        load_gene_ldscore_index(index)
    with pytest.raises(FileExistsError, match="incomplete gene LD-score index"):
        _build_index(index_inputs, index, monkeypatch)
    _build_index(index_inputs, index, monkeypatch, overwrite=True)
    load_gene_ldscore_index(index)


def test_a_failure_during_publication_leaves_an_index_no_reader_accepts(tmp_path, index_inputs, monkeypatch):
    from ldsc._ldscore.gene_index import build as gene_index_build
    from ldsc._ldscore.gene_index.scoring import load_gene_ldscore_index

    index = tmp_path / "parent" / "index"
    _build_index(index_inputs, index, monkeypatch)
    calls = []

    def fail_on_published_index(path):
        # The first load of the destination is the overwrite preflight; the second follows publication.
        if Path(path) == index:
            calls.append(path)
            if len(calls) == 2:
                raise LDSCInputError("simulated failure after publication began")
        return load_gene_ldscore_index(path)

    monkeypatch.setattr(gene_index_build, "load_gene_ldscore_index", fail_on_published_index)
    with pytest.raises(LDSCInputError, match="simulated"):
        _build_index(index_inputs, index, monkeypatch, overwrite=True)

    assert [path.name for path in index.parent.iterdir()] == ["index"]
    assert (index / "RUN_FAILED.txt").is_file()
    monkeypatch.undo()
    with pytest.raises(LDSCInputError):
        load_gene_ldscore_index(index)


@pytest.mark.parametrize("mode", ["chr_pos", "rsid"])
def test_ldscore_from_a_saved_annotation_prepares_inside_its_own_output_directory(tmp_path, mode):
    import numpy as np
    import pandas as pd

    from ldsc import GlobalConfig, LDScoreCalculator, LDScoreConfig, RefPanelConfig, run_annotate
    from ldsc._ldscore.outputs import LDScoreOutputConfig
    from ldsc._refpanel.adapters import R2MatrixRefPanel
    from tests.ldscore.test_ldscore_parallelism import _PANEL_CHROMS, _write_index_panel

    config = GlobalConfig(snp_identifier=mode, genome_build="hg19" if mode == "chr_pos" else None)
    _write_index_panel(tmp_path / "panel", _PANEL_CHROMS, snp_identifier=mode)
    panel = R2MatrixRefPanel(config, RefPanelConfig(backend="r2_matrix", r2_matrix_dir=str(tmp_path / "panel")))
    numerical = LDScoreConfig(ld_wind_kb=1, whole_chromosome_ok=True)
    metadata = pd.concat([pair[0] for pair in _PANEL_CHROMS.values()], ignore_index=True)
    baseline = tmp_path / "baseline.annot"
    metadata.loc[:, ["CHR", "POS", "SNP", "A1", "A2"]].assign(base=1).to_csv(baseline, sep="\t", index=False)
    bed = tmp_path / "query.bed"
    bed.write_text("chr1\t150\t350\nchr2\t50\t250\n")
    annotated = tmp_path / "annotated"
    saved = run_annotate(ref_snp_annot_sources=[str(baseline)], query_annot_bed_sources=[str(bed)],
                         genome_build="hg19", global_config=config, output_dir=annotated)
    annotated_before = _tree(annotated)

    result = LDScoreCalculator().run(saved, panel, numerical, config, LDScoreOutputConfig(tmp_path / "ldscores"))

    # The annotate output is this run's input: nothing appears or changes there.
    assert _tree(annotated) == annotated_before
    assert not list((tmp_path / "ldscores").glob(".ldsc-annotation-*"))
    # Preparing in the reading run's own workspace leaves the LD scores unchanged.
    with saved:
        saved.read("1")
        expected = LDScoreCalculator().run(saved, panel, numerical, config, LDScoreOutputConfig(tmp_path / "ldscores2"))
    np.testing.assert_array_equal(result.baseline_table.to_numpy(), expected.baseline_table.to_numpy())
    np.testing.assert_array_equal(result.read_queries(result.query_columns).to_numpy(),
                                  expected.read_queries(expected.query_columns).to_numpy())
