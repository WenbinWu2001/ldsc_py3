from __future__ import annotations

import pandas as pd
import pytest

from ldsc.outputs import (
    H2DirectoryWriter,
    H2OutputConfig,
    H2_REGRESSION_BIN_COLUMNS,
    PARTITIONED_H2_COLUMNS,
    PartitionedH2DirectoryWriter,
    PartitionedH2OutputConfig,
)


def _bins():
    return pd.DataFrame([{column: 0 for column in H2_REGRESSION_BIN_COLUMNS}])


def test_h2_no_overwrite_treats_reserved_derived_root_as_collision(tmp_path):
    output_dir = tmp_path / "h2"
    (output_dir / "plots").mkdir(parents=True)
    (output_dir / "plots" / "old.png").write_bytes(b"old")

    with pytest.raises(FileExistsError, match="plots"):
        H2DirectoryWriter().write(
            pd.DataFrame([{"trait_name": "trait"}]),
            H2OutputConfig(output_dir=output_dir),
            metadata={},
            diagnostic_bins=_bins(),
        )


def test_h2_successful_overwrite_removes_default_plots_and_postprocessing(tmp_path):
    output_dir = tmp_path / "h2"
    (output_dir / "plots").mkdir(parents=True)
    (output_dir / "plots" / "old.png").write_bytes(b"old")
    (output_dir / "postprocessing" / "liability-scale").mkdir(parents=True)
    (output_dir / "postprocessing" / "liability-scale" / "old.tsv").write_text("old")

    H2DirectoryWriter().write(
        pd.DataFrame([{"trait_name": "trait"}]),
        H2OutputConfig(output_dir=output_dir, overwrite=True),
        metadata={},
        diagnostic_bins=_bins(),
    )

    assert not (output_dir / "plots").exists()
    assert not (output_dir / "postprocessing").exists()


def test_partitioned_successful_overwrite_removes_default_plots_only(tmp_path):
    output_dir = tmp_path / "partitioned"
    (output_dir / "plots").mkdir(parents=True)
    (output_dir / "plots" / "old.png").write_bytes(b"old")
    unrelated = output_dir / "notes.txt"
    unrelated.write_text("keep")
    summary = pd.DataFrame(
        [{"category": "base", **{column: 0 for column in PARTITIONED_H2_COLUMNS if column != "category"}}]
    )

    PartitionedH2DirectoryWriter().write(
        summary,
        PartitionedH2OutputConfig(output_dir=output_dir, overwrite=True),
        metadata={"analysis_type": "functional_category", "headline_metric": "enrichment"},
    )

    assert not (output_dir / "plots").exists()
    assert unrelated.read_text() == "keep"


@pytest.mark.parametrize("workflow", ["plot", "scale"])
@pytest.mark.parametrize("source", ["absolute", "parent", "symlink", "missing"])
def test_derived_workflows_reject_invalid_declared_inputs_before_writing(tmp_path, workflow, source):
    import json
    from ldsc.errors import LDSCInputError
    from ldsc.h2_scale import convert_h2_scale
    from ldsc.plotting import plot_result

    outside = tmp_path / "outside.tsv"
    outside.write_text("immutable input\n")
    root = tmp_path / "result"
    (root / "diagnostics").mkdir(parents=True)
    if source == "symlink":
        (root / "linked.tsv").symlink_to(outside)
    token = {"absolute": str(outside), "parent": "../outside.tsv",
             "symlink": "linked.tsv", "missing": "missing.tsv"}[source]
    key = "summary" if workflow == "scale" else "ld_score_regression_bins"
    metadata_path = root / "diagnostics" / "metadata.json"
    metadata_path.write_text(json.dumps({"artifact_type": "h2_result", "files": {key: token}}))
    original = metadata_path.read_bytes()
    with pytest.raises(LDSCInputError, match="relative|escapes|missing"):
        if workflow == "scale":
            convert_h2_scale(root, samp_prev=0.5, pop_prev=0.01)
        else:
            plot_result(root)
    assert not (root / "plots").exists()
    assert not (root / "postprocessing").exists()
    assert metadata_path.read_bytes() == original
    assert outside.read_text() == "immutable input\n"
