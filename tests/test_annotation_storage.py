"""Independent selected-read and ownership checks for annotation staging."""

import numpy as np
import pandas as pd
import pytest

from ldsc._annotation_storage import AnnotationWorkspace, ColumnStore


def test_selected_reads_preserve_row_order_and_column_values(tmp_path):
    with AnnotationWorkspace(tmp_path) as workspace:
        store = ColumnStore.create(workspace.path / "values.npy", 5, ["base", "a", "b"])
        store.write(0, np.array([[1, 11, 21], [1, 12, 22]], dtype=np.float32))
        store.write(2, np.array([[1, 13, 23], [1, 14, 24], [1, 15, 25]], dtype=np.float32))
        actual = store.read(rows=[4, 0, 3, 0], columns=["b", "a"], max_read_rows=2)
        np.testing.assert_array_equal(actual, [[25, 15], [21, 11], [24, 14], [21, 11]])
        assert actual.dtype == np.float32
        assert not isinstance(actual, np.memmap)
        np.testing.assert_array_equal(np.load(store.path), [[1, 11, 21], [1, 12, 22], [1, 13, 23], [1, 14, 24], [1, 15, 25]])


def test_workspace_closes_only_its_owned_scratch(tmp_path):
    persistent = tmp_path / "query.1.annot.gz"
    persistent.write_bytes(b"persistent output")
    workspace = AnnotationWorkspace(tmp_path)
    private = workspace.path
    assert private.parent == tmp_path
    with pytest.raises(RuntimeError, match="failed consumption"):
        with workspace:
            (private / "staged").write_text("private")
            raise RuntimeError("failed consumption")
    assert not private.exists()
    assert persistent.read_bytes() == b"persistent output"
    workspace.close()
    with pytest.raises(ValueError, match="closed"):
        workspace.require_open()


def test_float64_storage_does_not_collapse_close_target_values(tmp_path):
    store = ColumnStore.create(tmp_path / "target.npy", 2, ["target"], dtype=np.float64)
    store.write(0, np.array([[1.00000001], [1.00000002]], dtype=np.float64))
    np.testing.assert_array_equal(store.read(), [[1.00000001], [1.00000002]])


def test_column_store_empty_selection_and_out_of_bounds(tmp_path):
    store = ColumnStore.create(tmp_path / "values.npy", 2, ["x"])
    store.write(0, np.array([[3], [4]]))
    assert store.read(rows=[], columns=["x"]).shape == (0, 1)
    assert store.read(columns=[]).shape == (2, 0)
    with pytest.raises(IndexError):
        store.read(rows=[2])
    with pytest.raises(KeyError):
        store.read(columns=["missing"])


@pytest.mark.parametrize("mode, expected", [("rsid", ["unique"]), ("chr_pos", ["dup", "unique", "dup"])])
def test_disk_identity_cleanup_is_global_across_chunks_and_chromosomes(tmp_path, mode, expected):
    from ldsc._annotation_identity import DiskIdentityIndex

    chunks = [pd.DataFrame({"CHR": [1, 1], "POS": [10, 20], "SNP": ["dup", "unique"]}),
              pd.DataFrame({"CHR": [2], "POS": [30], "SNP": ["dup"]})]
    with DiskIdentityIndex(tmp_path / "identity.sqlite", mode) as index:
        for chunk in chunks:
            index.add(chunk)
        actual = []
        drops = []
        for chunk in chunks:
            keep, dropped = index.select(chunk)
            actual.extend(chunk.loc[keep, "SNP"])
            drops.extend(dropped["reason"])
        assert actual == expected
        assert drops == (["duplicate_identity", "duplicate_identity"] if mode == "rsid" else [])
        assert index.connection.execute("PRAGMA temp_store").fetchone()[0] == 2
        assert index.connection.execute("PRAGMA cache_size").fetchone()[0] == -8192


def test_disk_identity_allele_policy_excludes_invalid_before_global_clusters(tmp_path):
    from ldsc._annotation_identity import DiskIdentityIndex

    frame = pd.DataFrame({"CHR": [1] * 7, "POS": list(range(1, 8)),
                          "SNP": ["multi", "dup", "valid", "bad", "multi", "dup", "valid"],
                          "A1": ["A", "A", "A", "A", "A", "T", "N"],
                          "A2": ["C", "G", "C", "T", "G", "C", "C"]})
    with DiskIdentityIndex(tmp_path / "identity.sqlite", "rsid_allele_aware") as index:
        index.add(frame.iloc[:4])
        index.add(frame.iloc[4:])
        keep, dropped = index.select(frame)
    assert frame.loc[keep, "POS"].tolist() == [3]
    assert dropped.groupby("reason")["SNP"].apply(list).to_dict() == {
        "multi_allelic_base_key": ["multi", "multi"], "duplicate_identity": ["dup", "dup"],
        "strand_ambiguous_allele": ["bad"], "invalid_allele": ["valid"],
    }


@pytest.mark.parametrize("sharded", [False, True])
def test_source_preparation_aligns_columns_then_cleans_global_duplicates(tmp_path, sharded):
    from ldsc._annotation_sources import prepare_annotation_sources

    meta = pd.DataFrame({"CHR": [1, 1, 2, 2], "SNP": ["repeat", "rs2", "repeat", "rs4"], "POS": [10, 20, 30, 40]})
    baseline = meta.assign(base=1, extra=[11, 12, 13, 14])
    query = meta.assign(first=[21, 22, 23, 24], second=[31, 32, 33, 34])
    paths = []
    for name, frame in [("base", baseline), ("query", query)]:
        selected = []
        for label, chunk in (frame.groupby("CHR") if sharded else [("all", frame)]):
            path = tmp_path / f"{name}.{label}.annot.gz"
            chunk.to_csv(path, sep="\t", index=False)
            selected.append(path)
        paths.append(selected)
    with AnnotationWorkspace(tmp_path / "output") as workspace:
        source = prepare_annotation_sources(workspace, *paths, mode="rsid", chunk_rows=1)
        assert list(source.shards) == ["1", "2"]
        assert source.baseline_columns == ("base", "extra")
        assert source.query_columns == ("first", "second")
        assert source.shards["1"].metadata()["SNP"].tolist() == ["rs2"]
        assert source.shards["2"].metadata()["SNP"].tolist() == ["rs4"]
        np.testing.assert_array_equal(source.shards["1"].read(columns=["second", "extra"]), [[32, 12]])
        np.testing.assert_array_equal(source.shards["2"].read(columns=["second", "extra"]), [[34, 14]])
        assert source.drops.n_rows == 2
        assert [row for frame in source.drops.frames() for row in frame["SNP"]] == ["repeat", "repeat"]


def test_aligned_whole_genome_column_sources_are_not_duplicate_observations(tmp_path):
    from ldsc._annotation_sources import prepare_annotation_sources

    meta = pd.DataFrame({"CHR": [1, 2], "SNP": ["rs1", "rs2"], "POS": [10, 20]})
    paths = []
    for name, values in [("base", [1, 1]), ("extra", [8, 9])]:
        path = tmp_path / f"{name}.annot"
        meta.assign(**{name: values}).to_csv(path, sep="\t", index=False)
        paths.append(path)
    with AnnotationWorkspace(tmp_path / "output") as workspace:
        source = prepare_annotation_sources(workspace, paths, [], mode="rsid", chunk_rows=1)
        assert source.drops.n_rows == 0
        np.testing.assert_array_equal(source.shards["2"].read(), [[1, 9]])


def test_bundle_does_not_cache_chromosome_arrays_and_close_invalidates_reads(tmp_path):
    import gc
    import weakref
    from ldsc._annotation_bundle import AnnotationBundle
    from ldsc._annotation_storage import AnnotationShard

    workspace = AnnotationWorkspace(tmp_path)
    metadata_path = workspace.path / "rows.parquet"
    pd.DataFrame({"CHR": [1], "POS": [10], "SNP": ["rs1"]}).to_parquet(metadata_path)
    store = ColumnStore.create(workspace.path / "values.npy", 1, ["base"])
    store.write(0, np.array([[1]]))
    bundle = AnnotationBundle({"1": AnnotationShard(metadata_path, (store,), 1)}, ["base"], [], workspace)
    bundle.validate()
    values = bundle.read("1")
    reference = weakref.ref(values)
    del values
    gc.collect()
    assert reference() is None
    bundle.close()
    with pytest.raises(ValueError, match="closed"):
        bundle.read("1")
