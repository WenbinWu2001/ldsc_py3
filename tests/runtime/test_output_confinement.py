"""Every file a run touches stays inside its output directory.

An empty output directory would otherwise normalize to the working directory,
where an authorized overwrite can replace or remove the user's own files.
"""

from __future__ import annotations

import gzip
import os
from pathlib import Path

import pytest

import ldsc
from ldsc import cli
from ldsc._runtime.scratch import ScratchWorkspace
from ldsc._runtime.paths import ensure_output_directory
from ldsc._runtime.workflow_logging import overwrite_failure_marker
from ldsc.errors import LDSCConfigError, LDSCInputError
from tests.runtime.test_failure_markers import FAILING_CLI_RUNS
from tests.ref_panel_helpers import HM3_CHR22_PLINK

PLINK_FIXTURE = HM3_CHR22_PLINK


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
    lambda: ScratchWorkspace(" "),
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
    from ldsc._ldscore.gene_index.artifact import load_gene_ldscore_index

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
    from ldsc._ldscore.gene_index.artifact import load_gene_ldscore_index

    index = tmp_path / "parent" / "index"
    _build_index(index_inputs, index, monkeypatch)
    # An interrupted publication leaves owned components without the root metadata written last,
    # beside the build directory, holding the lock it ran under, that only a successful build removes.
    (index / "metadata.json").unlink()
    (index / ".build").mkdir()
    (index / ".build" / "build-gene-ldscore-index.lock").write_text("pid=1 host=node target=index\n")

    with pytest.raises(LDSCInputError):
        load_gene_ldscore_index(index)
    with pytest.raises(FileExistsError, match="incomplete gene LD-score index"):
        _build_index(index_inputs, index, monkeypatch)
    _build_index(index_inputs, index, monkeypatch, overwrite=True)
    load_gene_ldscore_index(index)


def test_index_component_names_without_a_build_directory_are_never_replaced(tmp_path, index_inputs, monkeypatch):
    # A user directory whose entries merely carry index component names; no build ever ran here.
    index = tmp_path / "parent" / "index"
    (index / "chromosomes").mkdir(parents=True)
    (index / "diagnostics").mkdir()
    (index / "chromosomes" / "my_notes.txt").write_text("precious\n")
    (index / "diagnostics" / "notes.txt").write_text("notes\n")
    before = _tree(index)

    # Refused with and without --overwrite; a refused attempt leaves no .build/ that would make
    # the directory look like an interrupted index to the next attempt.
    for overwrite in (False, True):
        with pytest.raises(FileExistsError, match="nonempty but invalid"):
            _build_index(index_inputs, index, monkeypatch, overwrite=overwrite)
        assert not (index / ".build").exists()

    # Only the failure marker of the refused --overwrite attempt is added.
    assert _tree(index) == before | {"RUN_FAILED.txt"}
    assert (index / "chromosomes" / "my_notes.txt").read_text() == "precious\n"
    assert (index / "diagnostics" / "notes.txt").read_text() == "notes\n"


@pytest.mark.parametrize("overwrite", [False, True])
@pytest.mark.parametrize("components", [False, True], ids=["build-only", "with-components"])
def test_a_build_directory_without_the_builders_files_is_never_treated_as_owned(
    tmp_path, index_inputs, monkeypatch, components, overwrite,
):
    # A .build/ that holds neither the lock nor the log every build attempt writes there first
    # was not left by the builder, so neither it nor component names beside it are the builder's.
    index = tmp_path / "parent" / "index"
    (index / ".build" / "notes").mkdir(parents=True)
    (index / ".build" / "notes" / "keep.txt").write_text("user notes\n")
    if components:
        (index / "chromosomes").mkdir()
        (index / "chromosomes" / "my_notes.txt").write_text("precious\n")
    before = _snapshot(index)

    with pytest.raises(FileExistsError, match="nonempty but invalid"):
        _build_index(index_inputs, index, monkeypatch, overwrite=overwrite)

    after = _snapshot(index)
    after.pop("RUN_FAILED.txt", None)
    assert after == before


@pytest.mark.parametrize("beside", ["nothing", "complete-index", "component-names"])
def test_an_empty_build_directory_is_reused_but_owns_nothing_beside_it(tmp_path, index_inputs, monkeypatch, beside):
    from ldsc._ldscore.gene_index.artifact import load_gene_ldscore_index

    # On a shared filesystem a successful build can fail to remove its build directory: the lock it
    # still holds open is renamed to a temporary file that disappears later, leaving .build/ empty.
    index = tmp_path / "parent" / "index"
    if beside == "complete-index":
        _build_index(index_inputs, index, monkeypatch)
    elif beside == "component-names":
        (index / "chromosomes").mkdir(parents=True)
        (index / "chromosomes" / "my_notes.txt").write_text("precious\n")
    (index / ".build").mkdir(parents=True)

    if beside == "component-names":
        # An empty .build/ holds no record of an earlier attempt, so the component names stay the user's.
        for overwrite in (False, True):
            with pytest.raises(FileExistsError, match="nonempty but invalid"):
                _build_index(index_inputs, index, monkeypatch, overwrite=overwrite)
        assert (index / "chromosomes" / "my_notes.txt").read_text() == "precious\n"
        return
    _build_index(index_inputs, index, monkeypatch, overwrite=beside == "complete-index")
    assert not (index / ".build").exists()
    load_gene_ldscore_index(index)


@pytest.mark.skipif(os.geteuid() == 0, reason="directory permissions do not bind the superuser")
def test_a_build_directory_a_successful_build_removed_partway_is_removed_by_the_next_build(
    tmp_path, index_inputs, monkeypatch,
):
    from ldsc._ldscore.gene_index.artifact import load_gene_ldscore_index

    index = tmp_path / "parent" / "index"
    with pytest.raises(LDSCInputError):
        _build_index(index_inputs, index, monkeypatch, prefix="missing-panel")
    # An entry the successful build cannot remove (on a cluster, typically "Device or resource busy").
    busy = index / ".build" / "history" / "busy"
    busy.mkdir(parents=True)
    (busy / "held.txt").write_text("held\n")
    busy.chmod(0o500)
    try:
        with pytest.warns(UserWarning, match="build directory could not be removed"):
            _build_index(index_inputs, index, monkeypatch)
        load_gene_ldscore_index(index)
        # Whatever the removal reached, the lock file that marks .build/ as the builder's is left.
        assert (index / ".build" / "build-gene-ldscore-index.lock").is_file()
    finally:
        busy.chmod(0o700)

    # As the warning says, the next build of this index removes the directory.
    _build_index(index_inputs, index, monkeypatch, overwrite=True)
    assert not (index / ".build").exists()
    load_gene_ldscore_index(index)


def _snapshot(root: Path) -> dict[str, bytes | None]:
    """Map every entry under ``root`` to its bytes (``None`` for a directory)."""
    return {path.relative_to(root).as_posix(): path.read_bytes() if path.is_file() else None
            for path in root.rglob("*")}


@pytest.mark.parametrize("existing_index", [False, True], ids=["empty", "complete-index"])
@pytest.mark.parametrize("link", [".build", ".build/build-gene-ldscore-index.lock"])
def test_gene_index_build_never_follows_a_symlink_in_its_build_directory(
    tmp_path, index_inputs, monkeypatch, link, existing_index,
):
    index = tmp_path / "parent" / "index"
    index.mkdir(parents=True)
    if existing_index:
        _build_index(index_inputs, index, monkeypatch)
    elsewhere = tmp_path / "elsewhere"
    (elsewhere / "projectA").mkdir(parents=True)
    (elsewhere / "projectA" / "results.txt").write_text("data\n")
    (elsewhere / "keep.txt").write_text("top\n")
    (index / link).parent.mkdir(exist_ok=True)
    (index / link).symlink_to(elsewhere if link == ".build" else elsewhere / "keep.txt")
    before = _snapshot(elsewhere)

    with pytest.raises(FileExistsError, match=rf"{index / link} is a symlink"):
        _build_index(index_inputs, index, monkeypatch, overwrite=existing_index)

    # Nothing was created in, written to or removed from the link target.
    assert _snapshot(elsewhere) == before
    assert (index / link).is_symlink()


def test_build_directory_cleanup_never_deletes_through_a_symlinked_build_directory(tmp_path, monkeypatch):
    from ldsc._ldscore.gene_index import build as gene_index_build

    # The state a race leaves if .build/ is replaced by a symlink after the build took its lock:
    # the cleanup that follows a successful publication must not empty the link target.
    monkeypatch.setattr(gene_index_build, "set_workflow_log_path", lambda path: None)
    index = tmp_path / "index"
    index.mkdir()
    elsewhere = tmp_path / "elsewhere"
    (elsewhere / "results").mkdir(parents=True)
    (elsewhere / "results" / "table.tsv").write_text("data\n")
    (elsewhere / "thesis.tex").write_text("text\n")
    (elsewhere / "build-gene-ldscore-index.log").write_text("log\n")
    (index / ".build").symlink_to(elsewhere)

    with pytest.warns(UserWarning, match="build directory could not be removed"):
        gene_index_build._finalize_gene_index_log(index / ".build" / "build-gene-ldscore-index.log", index)

    assert (elsewhere / "results" / "table.tsv").read_text() == "data\n"
    assert (elsewhere / "thesis.tex").read_text() == "text\n"
    assert (index / ".build").is_symlink()


def test_gene_index_build_refuses_a_build_entry_that_is_not_a_directory(tmp_path, index_inputs, monkeypatch):
    index = tmp_path / "parent" / "index"
    _build_index(index_inputs, index, monkeypatch)
    (index / ".build").write_text("not the builder's\n")

    with pytest.raises(FileExistsError, match=rf"{index / '.build'} is not a directory"):
        _build_index(index_inputs, index, monkeypatch, overwrite=True)

    assert (index / ".build").read_text() == "not the builder's\n"


def test_a_failure_during_publication_leaves_an_index_no_reader_accepts(tmp_path, index_inputs, monkeypatch):
    from ldsc._ldscore.gene_index import build as gene_index_build
    from ldsc._ldscore.gene_index.artifact import load_gene_ldscore_index

    index = tmp_path / "parent" / "index"
    _build_index(index_inputs, index, monkeypatch)
    real_replace = gene_index_build.os.replace

    def fail_before_the_root_metadata(source, target):
        # Publication writes the root metadata last; fail just before it.
        if Path(target) == index / "metadata.json":
            raise OSError("simulated failure during publication")
        return real_replace(source, target)

    with monkeypatch.context() as patcher:
        patcher.setattr(gene_index_build.os, "replace", fail_before_the_root_metadata)
        with pytest.raises(OSError, match="simulated"):
            _build_index(index_inputs, index, monkeypatch, overwrite=True)

    assert [path.name for path in index.parent.iterdir()] == ["index"]
    assert (index / "RUN_FAILED.txt").is_file()
    with pytest.raises(LDSCInputError):
        load_gene_ldscore_index(index)

    # The build directory the interrupted attempt left marks the remains as the builder's own.
    with pytest.raises(FileExistsError, match="incomplete gene LD-score index"):
        _build_index(index_inputs, index, monkeypatch)
    _build_index(index_inputs, index, monkeypatch, overwrite=True)
    load_gene_ldscore_index(index)


@pytest.mark.parametrize("existing_index", [False, True], ids=["empty", "complete-index"])
def test_a_concurrent_overwrite_refused_by_the_lock_does_not_stop_the_running_build(
    tmp_path, index_inputs, monkeypatch, existing_index,
):
    import threading

    from ldsc._ldscore.gene_index import build as gene_index_build
    from ldsc._ldscore.gene_index.artifact import load_gene_ldscore_index
    from tests.ldscore.test_plink_workflow_resolution import index_args

    index = tmp_path / "parent" / "index"
    if existing_index:
        _build_index(index_inputs, index, monkeypatch)
    duplicate_args = index_args(index_inputs, str(index_inputs / "1000G.EUR.QC."), index, monkeypatch)
    duplicate_args.overwrite = True
    duplicate = {}
    real_stage = gene_index_build._create_gene_index_stage

    def stage_while_a_duplicate_job_starts(destination):
        # The same --overwrite job submitted again while this build holds the lock: it is refused
        # by the lock, and its own failure-marker scope (another thread) writes RUN_FAILED.txt.
        def run_duplicate():
            try:
                gene_index_build.run_build_gene_ldscore_index_from_args(duplicate_args)
            except Exception as exc:  # noqa: BLE001 - the refusal is inspected below
                duplicate["error"] = exc

        thread = threading.Thread(target=run_duplicate)
        thread.start()
        thread.join()
        duplicate["marker"] = (index / "RUN_FAILED.txt").is_file()
        return real_stage(destination)

    monkeypatch.setattr(gene_index_build, "_create_gene_index_stage", stage_while_a_duplicate_job_starts)
    _build_index(index_inputs, index, monkeypatch, overwrite=existing_index)

    assert isinstance(duplicate["error"], LDSCInputError)
    assert "Another build-gene-ldscore-index process is active" in str(duplicate["error"])
    assert duplicate["marker"]
    # The running build still publishes, and its successful marker scope removes the marker.
    load_gene_ldscore_index(index)
    assert not (index / "RUN_FAILED.txt").exists()
    assert not (index / ".build").exists()


@pytest.mark.parametrize("intruder", ["user-file", "another-index"])
@pytest.mark.parametrize("overwrite", [False, True])
def test_publication_refuses_a_destination_that_changed_during_the_build(
    tmp_path, index_inputs, monkeypatch, intruder, overwrite,
):
    import shutil

    from ldsc._ldscore.gene_index import build as gene_index_build
    from ldsc._ldscore.gene_index.artifact import load_gene_ldscore_index

    index = tmp_path / "parent" / "index"
    if intruder == "another-index":
        # A failed attempt's leftovers, so every component of the other index is a new name.
        with pytest.raises(LDSCInputError):
            _build_index(index_inputs, index, monkeypatch, prefix="missing-panel", overwrite=overwrite)
    elif overwrite:
        _build_index(index_inputs, index, monkeypatch)
    other = tmp_path / "other"
    _build_index(index_inputs, other, monkeypatch)
    real_stage = gene_index_build._create_gene_index_stage

    def stage_while_the_destination_changes(destination):
        # The output check ran before the lock; meanwhile a user adds a file, or another
        # builder that held the lock first publishes its index.
        staged = real_stage(destination)
        if intruder == "user-file":
            (index / "results.txt").write_text("written during the build\n")
        else:
            # Index components alone, no other new name, must stop publication.
            for entry in other.iterdir():
                shutil.copytree(entry, index / entry.name) if entry.is_dir() else shutil.copy2(entry, index)
        return staged

    monkeypatch.setattr(gene_index_build, "_create_gene_index_stage", stage_while_the_destination_changes)
    before = None

    def snapshot_before_publication(staged_index, destination, **kwargs):
        nonlocal before
        before = _snapshot(destination)
        return real_publish(staged_index, destination, **kwargs)

    real_publish = gene_index_build._publish_staged_gene_ldscore_index
    monkeypatch.setattr(gene_index_build, "_publish_staged_gene_ldscore_index", snapshot_before_publication)
    with pytest.raises(FileExistsError, match="changed during the build"):
        _build_index(index_inputs, index, monkeypatch, overwrite=overwrite)

    # Nothing outside this run's own build directory was replaced or removed.
    def outside_the_build(snapshot):
        return {name: data for name, data in snapshot.items()
                if not name.startswith(".build") and name != "RUN_FAILED.txt"}

    assert outside_the_build(_snapshot(index)) == outside_the_build(before)
    if intruder == "user-file":
        assert (index / "results.txt").is_file()
    if intruder == "another-index" or overwrite:
        load_gene_ldscore_index(index)


def _saved_annotation(tmp_path, mode):
    """Run annotate once on the two-chromosome R2 test panel; return what LD scoring needs."""
    import pandas as pd

    from ldsc import GlobalConfig, LDScoreConfig, RefPanelConfig, run_annotate
    from ldsc._refpanel.adapters import R2MatrixRefPanel
    from tests.ldscore.test_ldscore_parallelism import _PANEL_CHROMS, _write_index_panel

    config = GlobalConfig(snp_identifier=mode, genome_build="hg19" if mode == "chr_pos" else None)
    _write_index_panel(tmp_path / "panel", _PANEL_CHROMS, snp_identifier=mode)
    panel = R2MatrixRefPanel(config, RefPanelConfig(backend="r2_matrix", r2_matrix_dir=str(tmp_path / "panel")))
    metadata = pd.concat([pair[0] for pair in _PANEL_CHROMS.values()], ignore_index=True)
    baseline = tmp_path / "baseline.annot"
    metadata.loc[:, ["CHR", "POS", "SNP", "A1", "A2"]].assign(base=1).to_csv(baseline, sep="\t", index=False)
    bed = tmp_path / "query.bed"
    bed.write_text("chr1\t150\t350\nchr2\t50\t250\n")
    saved = run_annotate(ref_snp_annot_sources=[str(baseline)], query_annot_bed_sources=[str(bed)],
                         genome_build="hg19", global_config=config, output_dir=tmp_path / "annotated")
    return saved, panel, LDScoreConfig(ld_wind_kb=1, whole_chromosome_ok=True), config


@pytest.mark.parametrize("mode", ["chr_pos", "rsid"])
def test_ldscore_from_a_saved_annotation_prepares_inside_its_own_output_directory(tmp_path, mode):
    import numpy as np

    from ldsc import LDScoreCalculator
    from ldsc._ldscore.outputs import LDScoreOutputConfig

    saved, panel, numerical, config = _saved_annotation(tmp_path, mode)
    annotated = tmp_path / "annotated"
    annotated_before = _tree(annotated)

    result = LDScoreCalculator().run(saved, panel, numerical, config, LDScoreOutputConfig(tmp_path / "ldscores"))

    # The annotate output is this run's input: nothing appears or changes there.
    assert _tree(annotated) == annotated_before
    assert not list((tmp_path / "ldscores").glob(".ldsc-annotation-*"))
    # Preparing in a directory the caller chooses leaves the LD scores unchanged.
    with saved.prepared_in(tmp_path / "scratch") as bundle:
        expected = LDScoreCalculator().run(bundle, panel, numerical, config, LDScoreOutputConfig(tmp_path / "ldscores2"))
    np.testing.assert_array_equal(result.baseline_table.to_numpy(), expected.baseline_table.to_numpy())
    np.testing.assert_array_equal(result.read_queries(result.query_columns).to_numpy(),
                                  expected.read_queries(expected.query_columns).to_numpy())


@pytest.mark.parametrize("mode", ["chr_pos", "rsid"])
def test_a_saved_annotation_is_read_only_after_explicit_preparation(tmp_path, mode):
    from ldsc import AnnotationResult, LDScoreCalculator
    from ldsc.errors import LDSCUsageError

    saved, panel, numerical, config = _saved_annotation(tmp_path, mode)
    annotated = tmp_path / "annotated"
    annotated_before = _tree(annotated)
    assert isinstance(saved, AnnotationResult)
    assert {chrom: shard.n_rows for chrom, shard in saved.shards.items()} == {"1": 4, "2": 3}

    # Without an output directory, LD scoring refuses a saved result in every mode.
    with pytest.raises(LDSCUsageError, match="prepared_in"):
        LDScoreCalculator().run(saved, panel, numerical, config, output_config=None)
    assert _tree(annotated) == annotated_before

    scratch = tmp_path / "scratch"
    with saved.prepared_in(scratch) as bundle:
        assert bundle.read("1", columns=["query"])[:, 0].tolist() == [0, 1, 1, 0]
        assert bundle.read("2", columns=["query"])[:, 0].tolist() == [1, 1, 0]
        in_memory = LDScoreCalculator().run(bundle, panel, numerical, config, output_config=None)
        assert list(in_memory.query_columns) == ["query"]
        assert len(list(scratch.glob(".ldsc-annotation-*"))) == 1
    assert list(scratch.iterdir()) == []
    assert _tree(annotated) == annotated_before
