"""Output ownership and retry checks shared by materializing workflows."""

import pandas as pd
import pytest

from ldsc import cli
from ldsc._logging import overwrite_failure_marker
from ldsc.outputs import (
    AnnotationDirectoryWriter,
    H2DirectoryWriter,
    H2OutputConfig,
    LDScoreDirectoryWriter,
    PartitionedH2DirectoryWriter,
    QuantileH2DirectoryWriter,
    QueryR2DirectoryWriter,
    RgDirectoryWriter,
)
from ldsc.path_resolution import ensure_output_directory, remove_output_artifacts
from ldsc.errors import LDSCInputError


@pytest.mark.parametrize("overwrite", [False, True])
@pytest.mark.parametrize("writer,artifact", [
    (AnnotationDirectoryWriter, "query.22.annot.gz"),
    (H2DirectoryWriter, "h2.tsv"),
    (LDScoreDirectoryWriter, "ldscore.query.batch00001.parquet"),
    (PartitionedH2DirectoryWriter, "partitioned_h2.tsv"),
    (QuantileH2DirectoryWriter, "quantile_h2.tsv"),
    (QueryR2DirectoryWriter, "query_r2.tsv"),
    (RgDirectoryWriter, "rg.tsv"),
])
def test_writer_families_ignore_failure_markers_but_protect_owned_outputs(
    tmp_path, overwrite, writer, artifact
):
    marker = tmp_path / "RUN_FAILED.txt"
    marker.write_text("old failure")
    unrelated = tmp_path / "notes.txt"
    unrelated.write_text("keep")
    family = writer.artifact_family(tmp_path)
    assert family.preflight(overwrite=overwrite) == []
    assert "RUN_FAILED" not in str(family.metadata_files())

    owned = tmp_path / artifact
    owned.write_text("old output")
    family = writer.artifact_family(tmp_path)
    with pytest.raises(FileExistsError, match="overwrite"):
        family.preflight(overwrite=False)
    stale = family.preflight(overwrite=True)
    remove_output_artifacts(stale)
    assert marker.read_text() == "old failure"
    assert unrelated.read_text() == "keep"
    assert owned.exists() == (owned not in stale)


def test_failure_marker_uses_the_normalized_output_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("LDSC_RETRY_OUTPUT", str(tmp_path / "actual"))
    token = "$LDSC_RETRY_OUTPUT/result"
    output = ensure_output_directory(token)
    with pytest.raises(RuntimeError, match="failed"):
        with overwrite_failure_marker(token, overwrite=True, command="test"):
            raise RuntimeError("failed")
    assert (output / "RUN_FAILED.txt").is_file()
    assert not (tmp_path / "$LDSC_RETRY_OUTPUT").exists()
    with overwrite_failure_marker(token, overwrite=True, command="retry"):
        pass
    assert not (output / "RUN_FAILED.txt").exists()


@pytest.mark.parametrize("command,option,nested", [
    ("h2", "--output-dir", ""),
    ("plot", "--result-dir", "plots"),
    ("convert-h2-scale", "--h2-result-dir", "postprocessing/liability-scale"),
])
def test_cli_marker_follows_last_destination_option(tmp_path, command, option, nested):
    first, last = tmp_path / "first", tmp_path / "last"
    # Invalid input fails before computation; argparse chooses the last value.
    assert cli.run_cli([
        command, option, str(first), f"{option}={last}", "--overwrite",
    ]) != 0
    assert not first.exists()
    assert (last / nested / "RUN_FAILED.txt").is_file()


@pytest.mark.parametrize("workflow,nested", [
    ("plot", "plots"), ("scale", "postprocessing/liability-scale"),
])
@pytest.mark.parametrize("custom", [False, True])
def test_derived_workflows_retry_in_normalized_output_scope(tmp_path, monkeypatch, workflow, nested, custom):
    from ldsc.h2_scale import convert_h2_scale
    from ldsc.plotting import plot_result

    monkeypatch.chdir(tmp_path)
    source = tmp_path / "source"
    H2DirectoryWriter().write(
        pd.DataFrame([{"trait_name": "trait", "total_h2_obs": .2, "total_h2_obs_se": .03}]),
        H2OutputConfig(source), metadata={"trait_name": "trait"},
        diagnostic_bins=pd.DataFrame([
            {"bin": i, "n_snps": 10, "ld_score_min": i, "ld_score_max": i + 1,
             "mean_ld_score": i + .5, "mean_chi_square": 1 + .1 * i, "sd_chi_square": .2,
             "mean_sample_size": 1000, "mean_fitted_chi_square": 1 + .1 * i,
             "mean_regression_weight": .5} for i in (1, 2)
        ]),
    )
    monkeypatch.setenv("LDSC_RETRY_SOURCE", str(source))
    monkeypatch.setenv("LDSC_RETRY_DESTINATION", str(tmp_path / "custom"))
    output = tmp_path / "custom" if custom else source / nested
    kwargs = {"overwrite": True}
    if custom:
        kwargs["output_dir"] = "$LDSC_RETRY_DESTINATION"
    if workflow == "scale":
        kwargs.update(samp_prev=.5, pop_prev=.01)
    run = plot_result if workflow == "plot" else convert_h2_scale
    metadata_path = source / "diagnostics/metadata.json"
    original = metadata_path.read_bytes()
    metadata_path.unlink()
    with pytest.raises(LDSCInputError):
        run("$LDSC_RETRY_SOURCE", **kwargs)
    marker = output / "RUN_FAILED.txt"
    assert marker.is_file()
    metadata_path.write_bytes(original)
    unrelated = output / "notes.txt"
    unrelated.write_text("keep")
    artifact = run("$LDSC_RETRY_SOURCE", **kwargs)
    path = artifact.path if workflow == "plot" else artifact.table_path
    assert path.parent == output
    assert path.is_file()
    assert not marker.exists()
    assert unrelated.read_text() == "keep"
    assert metadata_path.read_bytes() == original
    assert not (tmp_path / "$LDSC_RETRY_SOURCE").exists()
    assert not (tmp_path / "$LDSC_RETRY_DESTINATION").exists()
    if workflow == "plot":
        from matplotlib import pyplot
        pyplot.close(artifact.figure)
