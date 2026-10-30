"""Shared provenance and file ownership in fitted-query publication."""

import json
import shutil

import numpy as np
import pandas as pd
import pytest

from ldsc.errors import LDSCInputError
from ldsc._regression.outputs import (
    PARTITIONED_H2_COLUMNS,
    PartitionedH2DirectoryWriter,
    PartitionedH2FitArtifacts,
    PartitionedH2OutputConfig,
    stage_partitioned_h2_fit,
)
from ldsc._regression.quantile_h2 import load_fitted_partitioned_model


@pytest.fixture
def fitted_query(tmp_path):
    full = pd.DataFrame({column: [1.0, 2.0] for column in PARTITIONED_H2_COLUMNS})
    full["category"] = ["base", "q"]
    full["overlap_annot"] = True
    deletes = pd.DataFrame({"delete_block": [0, 1, 2], "base": [.1, .2, .3], "q": [.4, .5, .6]})
    linked = tmp_path / "ldscores"
    linked.mkdir()
    metadata = {
        "annotation_types": {"base": "binary", "q": "quantitative", "unused": "binary"},
        "ldscore_dir": str(linked), "count_kind": "common", "analysis_type": "cell_type_specific",
        "n_snps": 999,
    }
    fits = stage_partitioned_h2_fit(tmp_path / "private", full, deletes, {"n_snps": 17})
    return full.iloc[[1]].reset_index(drop=True), fits, metadata, full, deletes


def publish(root, fitted_query):
    summary, fits, metadata, _, _ = fitted_query
    PartitionedH2DirectoryWriter().write(
        summary, PartitionedH2OutputConfig(root, write_per_query_results=True), {"q": fits}, metadata,
    )
    return root / "diagnostics/query_annotations/0001_q"


def test_shared_metadata_and_owned_files_survive_result_relocation(tmp_path, fitted_query):
    _, fits, metadata, full, deletes = fitted_query
    before = {path.name: (path.read_bytes(), path.stat().st_ino) for path in (fits.full, fits.delete_values)}
    query = publish(tmp_path / "result", fitted_query)
    payload = json.loads((query / "metadata.json").read_text())
    assert payload["global_metadata"] == "../../metadata.json"
    assert "annotation_types" not in payload
    assert "ldscore_dir" not in payload
    for path in (fits.full, fits.delete_values):
        assert not path.exists()
        assert (query / path.name).read_bytes() == before[path.name][0]
        assert (query / path.name).stat().st_ino == before[path.name][1]
    # The shared reference must follow a copied result tree, independently of cwd.
    shutil.copytree(tmp_path / "result", tmp_path / "relocated")
    shutil.rmtree(tmp_path / "result")
    model = load_fitted_partitioned_model(tmp_path / "relocated/diagnostics/query_annotations/0001_q")
    assert model.metadata["annotation_types"] == metadata["annotation_types"]
    assert model.metadata["n_snps"] == 17
    assert model.metadata["artifact_type"] == "partitioned_h2_query_result"
    np.testing.assert_array_equal(model.tau, full["coefficient"])
    np.testing.assert_array_equal(model.tau_delete, deletes[["base", "q"]])


def test_persistent_fit_files_are_copied_and_remain_reusable(tmp_path, fitted_query):
    summary, owned, metadata, full, deletes = fitted_query
    persistent = PartitionedH2FitArtifacts(owned.full, owned.delete_values, owned.metadata)
    before = {path: path.read_bytes() for path in (persistent.full, persistent.delete_values, persistent.metadata)}
    for name in ("first", "second"):
        query = publish(tmp_path / name, (summary, persistent, metadata, full, deletes))
        for path, content in before.items():
            assert path.read_bytes() == content
        for path in (persistent.full, persistent.delete_values):
            assert (query / path.name).read_bytes() == before[path]
            assert (query / path.name).stat().st_ino != path.stat().st_ino


def test_failed_owned_publication_does_not_expose_partial_query_tree(tmp_path, monkeypatch, fitted_query):
    from ldsc._regression import outputs

    original = outputs.atomic_write_json

    def fail_query(payload, path):
        if payload.get("artifact_type") == "partitioned_h2_query_result":
            raise OSError("injected publication failure")
        return original(payload, path)

    monkeypatch.setattr(outputs, "atomic_write_json", fail_query)
    with pytest.raises(OSError, match="injected publication failure"):
        publish(tmp_path / "result", fitted_query)
    diagnostics = tmp_path / "result/diagnostics"
    assert not (diagnostics / "query_annotations").exists()
    assert not list(diagnostics.glob(".query_annotations.tmp.*"))


@pytest.mark.parametrize("damage", ["missing_parent", "invalid_reference"])
def test_shared_metadata_cannot_be_silently_lost(tmp_path, fitted_query, damage):
    query = publish(tmp_path / "result", fitted_query)
    if damage == "missing_parent":
        (tmp_path / "result/diagnostics/metadata.json").unlink()
    else:
        path = query / "metadata.json"
        payload = json.loads(path.read_text())
        payload["global_metadata"] = "../../../unrelated.json"
        path.write_text(json.dumps(payload))
    with pytest.raises(LDSCInputError, match="metadata"):
        load_fitted_partitioned_model(query)


@pytest.mark.parametrize(
    ("relative_metadata", "artifact_type"),
    [("metadata.json", "rg_pair_result"), ("diagnostics/metadata.json", "h2_result")],
)
def test_fitted_model_loader_names_a_result_of_another_type(tmp_path, relative_metadata, artifact_type):
    metadata_path = tmp_path / "result" / relative_metadata
    metadata_path.parent.mkdir(parents=True)
    metadata_path.write_text(json.dumps({"artifact_type": artifact_type}))
    with pytest.raises(LDSCInputError, match=artifact_type):
        load_fitted_partitioned_model(tmp_path / "result")
