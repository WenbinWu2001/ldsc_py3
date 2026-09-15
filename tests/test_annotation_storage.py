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


def test_binary_store_packs_snp_bits_and_preserves_partial_writes(tmp_path):
    values = np.array([[1, 0], [0, 1], [1, 0], [0, 1], [0, 1], [0, 1],
                       [0, 1], [1, 0], [1, 0], [0, 1], [1, 0]], dtype=bool)
    store = ColumnStore.create(tmp_path / "binary.npy", 11, ["a", "b"], dtype=bool)
    for first, last in [(0, 3), (3, 9), (9, 11)]:
        store.write(first, values[first:last])
    # SNP zero is the low bit; the unused five bits of the last byte are zero.
    np.testing.assert_array_equal(np.load(store.path), np.array([[133, 122], [5, 2]], dtype=np.uint8))
    assert store.path.stat().st_size - store.offset == 4
    assert store.n_rows == 11 and store.bitorder == "little"
    np.testing.assert_array_equal(store.read(), values)
    np.testing.assert_array_equal(store.read(rows=[9, 7, 3, 7], columns=["b", "a", "b"], max_read_rows=1),
                                  [[1, 0, 1], [0, 1, 0], [1, 0, 1], [0, 1, 0]])
    store.write(7, np.array([[0, 1], [0, 1]]))
    np.testing.assert_array_equal(np.load(store.path), [[5, 250], [4, 3]])
    with pytest.raises(ValueError, match="zero or one"):
        store.write(0, np.array([[.5, 1]]))


@pytest.mark.parametrize("n_rows", [0, 1, 7, 8, 9, 17])
def test_binary_store_padding_empty_reads_and_worker_descriptor(tmp_path, n_rows):
    import pickle

    store = ColumnStore.create(tmp_path / "binary.npy", n_rows, ["zero", "one"], dtype=bool)
    store.write(0, np.ones((n_rows, 1), dtype=bool), columns=["one"])
    store.write(0, np.zeros((n_rows, 1), dtype=bool), columns=["zero"])
    # Worker transfer retains the logical length; NPY itself describes bytes.
    restored = pickle.loads(pickle.dumps(store))
    np.testing.assert_array_equal(restored.read(columns=["one", "zero"]),
                                  np.column_stack([np.ones(n_rows), np.zeros(n_rows)]))
    physical = np.load(store.path)
    assert physical.dtype == np.uint8
    assert physical.shape == ((n_rows + 7) // 8, 2)
    if n_rows:
        assert physical[-1, 1] == (1 << ((n_rows - 1) % 8 + 1)) - 1
    assert restored.read(rows=[]).shape == (0, 2)
    assert restored.read(columns=[]).shape == (n_rows, 0)


def test_binary_selected_reads_decode_bounded_runs_and_detect_truncation(tmp_path, monkeypatch):
    values = (np.arange(4097) % 3 == 0)[:, None]
    store = ColumnStore.create(tmp_path / "binary.npy", len(values), ["x"], dtype=bool)
    store.write(0, values)
    unpack = np.unpackbits
    decoded_bytes = []

    def unpackbits(packed, **kwargs):
        decoded_bytes.append(packed.nbytes)
        return unpack(packed, **kwargs)

    monkeypatch.setattr(np, "unpackbits", unpackbits)
    for rows in ([4096, 7, 8, 9, 7, 0], slice(19, 3, -3), np.arange(len(values)) % 257 == 0):
        np.testing.assert_array_equal(store.read(rows=rows, max_read_rows=7), values[rows])
    assert max(decoded_bytes) <= 2
    with store.path.open("r+b") as stream:
        stream.truncate(store.path.stat().st_size - 1)
    with pytest.raises(OSError, match="Truncated annotation shard"):
        store.read(rows=[4096])


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


def test_frame_spool_keeps_one_file_and_replays_exact_frames(tmp_path):
    from ldsc._annotation_storage import FrameSpool

    spool = FrameSpool(tmp_path / "spool")
    frames = [pd.DataFrame({"position": pd.Series([i], dtype="Int64"),
                            "label": [f"row{i}"]}, index=[0]) for i in range(100)]
    for frame in frames:
        spool.append(frame)
    spool.append(pd.DataFrame())
    assert spool.n_rows == spool.n_parts == 100
    assert len([p for p in spool.path.iterdir() if p.is_file()]) == 1
    for _ in range(2):
        for actual, expected in zip(spool.frames(), frames, strict=True):
            pd.testing.assert_frame_equal(actual, expected)
    spool.append(frames[0])
    assert len(list(spool.frames())) == 101


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


def test_identity_lookup_batches_respect_sqlite_limits_and_global_counts(tmp_path):
    from ldsc._annotation_identity import DiskIdentityIndex
    import sqlite3

    frame = pd.DataFrame({"CHR": ["1"] * 1100, "POS": np.arange(1100) + 1,
                          "SNP": [f"rs'{i}" for i in range(1100)]})
    with DiskIdentityIndex(tmp_path / "identity.sqlite", "rsid") as index:
        index.add(frame)
        index.add(frame.iloc[[1, 1000]].assign(CHR="2"))
        statements = []
        index.connection.setlimit(sqlite3.SQLITE_LIMIT_VARIABLE_NUMBER, 127)
        index.connection.set_trace_callback(statements.append)
        keep, dropped = index.select(frame)
        assert np.flatnonzero(~keep).tolist() == [1, 1000]
        assert dropped.SNP.tolist() == ["rs'1", "rs'1000"]
        assert dropped.reason.tolist() == ["duplicate_identity"] * 2
        lookups = [sql for sql in statements if sql.startswith("SELECT")]
        assert len(lookups) <= 9


def test_identity_add_batches_transactions_and_finishes_before_selection(tmp_path):
    from ldsc._annotation_identity import DiskIdentityIndex

    frame = pd.DataFrame({"CHR": ["1"], "POS": [1], "SNP": ["rs1"]})
    with DiskIdentityIndex(tmp_path / "identity.sqlite", "rsid") as index:
        for _ in range(10):
            index.add(frame)
        assert index.connection.in_transaction
        keep, dropped = index.select(frame)
        assert not index.connection.in_transaction
        assert keep.tolist() == [False]
        assert dropped.reason.tolist() == ["duplicate_identity"]


def test_identity_transaction_limit_spans_chunks_without_losing_duplicates(tmp_path, monkeypatch):
    from ldsc import _annotation_identity as identity

    monkeypatch.setattr(identity, "_TRANSACTION_ROWS", 8)
    frame = pd.DataFrame({"CHR": ["1"] * 11, "POS": np.arange(11) + 1,
                          "SNP": [f"rs{i}" for i in range(11)]})
    with identity.DiskIdentityIndex(tmp_path / "identity.sqlite", "rsid") as index:
        statements = []
        index.connection.set_trace_callback(statements.append)
        index.add(frame)
        assert statements.count("COMMIT") == 1
        index.add(frame.iloc[[0, 2, 3, 4, 5]])
        assert statements.count("COMMIT") == 2
        assert not index.connection.in_transaction
        keep, dropped = index.select(frame)
        assert np.flatnonzero(~keep).tolist() == [0, 2, 3, 4, 5]
        assert dropped.reason.eq("duplicate_identity").all()


def test_mixed_shards_preserve_omitted_alleles_and_global_duplicate_cleanup(tmp_path):
    from ldsc._annotation_sources import prepare_annotation_sources

    first = tmp_path / "first.annot"
    second = tmp_path / "second.annot"
    first.write_text("CHR SNP POS A1 A2 base\n1 unique1 10 A C 1\n1 duplicate 20 A C 2\n")
    second.write_text("CHR SNP POS base\n2 unique2 10 3\n2 duplicate 20 4\n")
    with AnnotationWorkspace(tmp_path / "out") as workspace:
        prepared = prepare_annotation_sources(workspace, [first, second], [], mode="rsid_allele_aware", chunk_rows=1)
        assert prepared.shards["1"].metadata().SNP.tolist() == ["unique1"]
        assert prepared.shards["2"].metadata().SNP.tolist() == ["unique2"]
        assert "A1" not in prepared.shards["2"].metadata()
        np.testing.assert_array_equal(prepared.shards["2"].read(columns=prepared.baseline_columns), [[3]])
        drops = pd.concat(prepared.drops.frames(), ignore_index=True)
        assert drops.reason.tolist() == ["duplicate_identity", "duplicate_identity"]
        assert drops.SNP.tolist() == ["duplicate", "duplicate"]


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
        np.testing.assert_array_equal(source.shards["2"].read(columns=source.baseline_columns + source.query_columns), [[1, 9]])


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
