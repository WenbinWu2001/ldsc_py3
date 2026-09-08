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
