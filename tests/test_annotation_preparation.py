"""Resource and scientific contracts for source-backed annotation preparation."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ldsc._annotation_sources import prepare_annotation_sources
from ldsc._annotation_storage import AnnotationWorkspace
from ldsc._logging import workflow_logging
from ldsc.errors import LDSCInputError


@pytest.mark.parametrize("chrom, chromosomes, snps", [(None, 2, 2), ("2", 1, 1)])
def test_preparation_log_marks_main_steps_and_retained_logical_rows(tmp_path, capsys, chrom, chromosomes, snps):
    paths = [tmp_path / f"{role}.annot" for role in ("baseline", "query")]
    for path, column in zip(paths, ("base", "query")):
        path.write_text(f"CHR SNP POS {column}\n1 duplicate 1 1\n1 rs1 2 .5\n2 duplicate 1 1\n2 rs2 2 .25\n")
    log_path = tmp_path / "preparation.log"
    with workflow_logging("annotate", log_path, log_level="INFO"), AnnotationWorkspace(tmp_path / "out") as workspace:
        prepared = prepare_annotation_sources(workspace, paths[:1], paths[1:], mode="rsid", chunk_rows=1, chrom=chrom)
        messages = log_path.read_text().splitlines()
        assert "Reading annotation inputs: baseline files=1, query files=1." in messages
        assert "Checking SNP identities and preparing chromosome annotations." in messages
        completed = [line for line in messages if line.startswith("Annotation preparation complete:")]
        assert len(completed) == 1
        assert completed[0].startswith(f"Annotation preparation complete: chromosomes={chromosomes}, retained SNPs={snps}, elapsed=")
        assert messages.index(completed[0]) > messages.index("Checking SNP identities and preparing chromosome annotations.") > messages.index("Reading annotation inputs: baseline files=1, query files=1.")
        np.testing.assert_array_equal(prepared.shards["2"].read(columns=["base", "query"]), [[.25, .25]])
    assert capsys.readouterr() == ("", "")


@pytest.mark.parametrize("rows, reached_identity", [("1 rs1 1 bad\n", False), ("1 rs1 1 1\n2 rs1 2 1\n", True)])
def test_preparation_failure_does_not_log_completion(tmp_path, rows, reached_identity):
    path = tmp_path / "baseline.annot"
    path.write_text("CHR SNP POS base\n" + rows)
    log_path = tmp_path / "preparation.log"
    with pytest.raises(LDSCInputError), workflow_logging("annotate", log_path, log_level="INFO"), AnnotationWorkspace(tmp_path / "out") as workspace:
        prepare_annotation_sources(workspace, [path], [], mode="rsid", chunk_rows=1)
    log = log_path.read_text()
    assert "Phase failed: validation" in log or "Reading annotation inputs:" in log
    assert ("Checking SNP identities" in log) == reached_identity
    assert "Annotation preparation complete:" not in log
    assert "Failed" in log and "Traceback" in log


def test_metadata_notice_once_per_source_on_each_preparation(tmp_path, caplog):
    paths = [tmp_path / f"baseline.{chrom}.annot" for chrom in (1, 2)]
    for chrom, path in enumerate(paths, 1):
        path.write_text(f"CHR SNP POS CM MAF base\n{chrom} rs{chrom}a 1 0 .2 1\n{chrom} rs{chrom}b 2 0 .3 .5\n")
    for run in range(2):
        caplog.clear()
        with caplog.at_level("INFO", logger="LDSC"), AnnotationWorkspace(tmp_path / f"out{run}") as workspace:
            prepared = prepare_annotation_sources(workspace, paths, [], mode="rsid", chunk_rows=1)
            notices = [record.getMessage() for record in caplog.records if "contains CM/MAF" in record.getMessage()]
            assert len(notices) == 2
            assert all(sum(str(path) in message for message in notices) == 1 for path in paths)
            np.testing.assert_array_equal(prepared.shards["1"].read(columns=["base"]), [[1], [.5]])
            assert prepared.shards["1"].metadata().CM.isna().all()


@pytest.mark.parametrize("declared", [False, True])
@pytest.mark.parametrize("chunk_rows", [None, 2, 2.0])
def test_sharded_scan_budget_does_not_shrink_with_chromosome_count(tmp_path, monkeypatch, declared, chunk_rows):
    paths = []
    for chromosome in range(1, 4):
        path = tmp_path / f"baseline.{chromosome}.annot"
        metadata = pd.DataFrame({"CHR": [chromosome] * 3,
                                 "SNP": [f"rs{chromosome}_{i}" for i in range(3)], "POS": [1, 2, 3]})
        values = pd.DataFrame(np.ones((3, 40)), columns=[f"b{i}" for i in range(40)])
        pd.concat([metadata, values], axis=1).to_csv(path, sep="\t", index=False)
        paths.append(path)
    original = pd.read_csv
    requested = []

    def read_csv(*args, **kwargs):
        if "chunksize" in kwargs:
            requested.append(kwargs["chunksize"])
        return original(*args, **kwargs)

    monkeypatch.setattr(pd, "read_csv", read_csv)
    with AnnotationWorkspace(tmp_path / "out") as workspace:
        prepared = prepare_annotation_sources(
            workspace, paths, [], mode="rsid", chunk_rows=chunk_rows,
            declared_chromosomes={str(p): str(i) for i, p in enumerate(paths, 1)} if declared else None,
        )
        assert list(prepared.shards) == ["1", "2", "3"]
        assert requested == [chunk_rows or 48770] * 6


def test_identity_preparation_uses_only_metadata_and_never_stages_numeric_values(tmp_path, monkeypatch):
    from ldsc._annotation_identity import DiskIdentityIndex

    path = tmp_path / "base.annot"
    path.write_text("CHR SNP POS base\n1 rs1 1 1\n2 rs2 2 -0.25\n")
    original_open = Path.open
    original_select = DiskIdentityIndex.select
    indexing = True

    def select(self, frame):
        nonlocal indexing
        indexing = False
        return original_select(self, frame)

    def open_path(self, mode="r", *args, **kwargs):
        assert self.suffix != ".bin"
        if indexing and mode == "rb":
            assert self.suffix != ".npy"
        return original_open(self, mode, *args, **kwargs)

    import pickle
    original_load = pickle.load

    def load(stream, *args, **kwargs):
        frame = original_load(stream, *args, **kwargs)
        if indexing and isinstance(frame, pd.DataFrame):
            assert "base" not in frame.columns
        return frame

    monkeypatch.setattr(DiskIdentityIndex, "select", select)
    monkeypatch.setattr(Path, "open", open_path)
    monkeypatch.setattr(pickle, "load", load)
    with AnnotationWorkspace(tmp_path / "out") as workspace:
        prepared = prepare_annotation_sources(workspace, [path], [], mode="rsid", chunk_rows=1)
        np.testing.assert_array_equal(prepared.shards["2"].read(columns=["base"]), [[-0.25]])


@pytest.mark.parametrize("sharded", [False, True])
def test_automatic_tiles_align_unequal_widths_and_preserve_logical_rows(tmp_path, sharded):
    rows, width = 4400, 500
    metadata = pd.DataFrame({"CHR": np.repeat([2, 1], rows),
                             "SNP": [f"rs{i}" for i in range(rows * 2)],
                             "POS": np.tile(np.arange(rows) + 1, 2)})
    metadata.loc[[0, rows], "SNP"] = "cross_chromosome_duplicate"
    values = (np.arange(rows * 2)[:, None] % 33 - 16) / 16 + np.arange(width) / 2
    query = pd.concat([metadata.assign(A1="A", A2="C"),
                       pd.DataFrame(values, columns=[f"q{i}" for i in range(width)])], axis=1)
    sources = []
    for kind, frame in [("baseline", metadata.assign(base=1)), ("query", query)]:
        paths = []
        for label, part in frame.groupby("CHR", sort=False) if sharded else [("all", frame)]:
            path = tmp_path / f"{kind}.{label}.annot"
            part.to_csv(path, sep="\t", index=False)
            paths.append(path)
        sources.append(paths)
    with AnnotationWorkspace(tmp_path / "out") as workspace:
        prepared = prepare_annotation_sources(workspace, *sources, mode="rsid_allele_aware")
        assert list(prepared.shards) == ["1", "2"]
        assert prepared.baseline_columns == ("base",)
        assert prepared.query_columns == tuple(f"q{i}" for i in range(width))
        for chrom, start in [("2", 0), ("1", rows)]:
            shard = prepared.shards[chrom]
            actual = shard.metadata()
            assert actual.SNP.tolist() == metadata.SNP.iloc[start + 1:start + rows].tolist()
            assert actual.A1.eq("A").all() and actual.A2.eq("C").all()
            np.testing.assert_array_equal(shard.read(columns=["q499", "base", "q0"]),
                                          np.column_stack([values[start + 1:start + rows, -1],
                                                           np.ones(rows - 1), values[start + 1:start + rows, 0]]))
            assert shard.read(rows=[0], columns=["base", "q0"]).dtype == np.float32
            assert np.load(shard.stores[0].path, mmap_mode="r").flags.f_contiguous
        drops = pd.concat(prepared.drops.frames(), ignore_index=True)
        assert drops.CHR.tolist() == (["1", "2"] if sharded else ["2", "1"])
        assert drops.reason.tolist() == ["duplicate_identity"] * 2
        assert not list(workspace.path.glob("source-*"))
        assert not list(workspace.path.glob("selection-*"))
        assert not list(workspace.path.glob("identity.sqlite*"))


@pytest.mark.parametrize("chunk_rows", [1, 3, None])
def test_chromosome_selection_follows_global_conflicts_and_preserves_order(tmp_path, chunk_rows):
    first = tmp_path / "first.annot"
    second = tmp_path / "second.annot"
    first.write_text("CHR SNP POS A1 A2 base\n1 duplicate 10 A C 1\n1 conflict 20 A C 2\n1 retained1 30 A C 3\n")
    second.write_text("CHR SNP POS A1 A2 base\n2 retained2 40 A C -1.25\n2 duplicate 20 G T 5\n2 conflict 30 A G 6\n2 invalid 10 N C 7\n2 retained3 5 A C 0.125\n")
    with AnnotationWorkspace(tmp_path / "out") as workspace:
        prepared = prepare_annotation_sources(workspace, [second, first], [],
                                              mode="rsid_allele_aware", chunk_rows=chunk_rows, chrom="2")
        assert list(prepared.shards) == ["2"]
        assert prepared.scope_chromosomes == ("1", "2")
        assert prepared.shards["2"].metadata().POS.tolist() == [40, 5]
        np.testing.assert_array_equal(prepared.shards["2"].read(columns=["base"]), [[-1.25], [.125]])
        drops = pd.concat(prepared.drops.frames(), ignore_index=True)
        assert drops.reason.tolist() == ["invalid_allele", "multi_allelic_base_key", "multi_allelic_base_key",
                                         "duplicate_identity", "duplicate_identity"]
        assert drops.CHR.tolist() == ["2", "1", "2", "1", "2"]


@pytest.mark.parametrize("defect", ["row-count", "row-order", "numeric"])
def test_failed_preparation_keeps_validation_and_workspace_ownership(tmp_path, defect):
    from ldsc.errors import LDSCInputError

    baseline = tmp_path / "baseline.annot"
    query = tmp_path / "query.annot"
    baseline.write_text("CHR SNP POS base\n1 rs1 1 1\n2 rs2 2 1\n1 rs3 3 1\n2 rs4 4 1\n")
    if defect == "row-count":
        query.write_text("CHR SNP POS q\n1 rs1 1 1\n2 rs2 2 2\n1 rs3 3 3\n")
        message = "row count mismatch"
    elif defect == "row-order":
        query.write_text("CHR SNP POS q\n1 rs1 1 1\n2 rs2 2 2\n2 rs4 4 4\n1 rs3 3 3\n")
        message = "SNP row mismatch"
    else:
        query.write_text("CHR SNP POS q\n1 rs1 1 1\n2 rs2 2 invalid\n1 rs3 3 3\n2 rs4 4 4\n")
        message = "input preflight failed"
    output = tmp_path / "out"
    output.mkdir()
    persistent = output / "existing.txt"
    persistent.write_text("user-owned")
    with pytest.raises(LDSCInputError, match=message):
        with AnnotationWorkspace(output) as workspace:
            prepare_annotation_sources(workspace, [baseline], [query], mode="rsid", chunk_rows=2, chrom="1")
    assert list(output.iterdir()) == [persistent]
    assert persistent.read_text() == "user-owned"


def test_empty_retained_result_preserves_complete_diagnostics(tmp_path):
    from ldsc.errors import LDSCInputError

    path = tmp_path / "duplicate.annot"
    path.write_text("CHR SNP POS base\n1 duplicated 1 1\n2 duplicated 2 2\n")
    with AnnotationWorkspace(tmp_path / "out") as workspace:
        with pytest.raises(LDSCInputError, match="retained no annotation rows") as error:
            prepare_annotation_sources(workspace, [path], [], mode="rsid", chunk_rows=1)
        assert error.value.annotation_drops.n_rows == 2
        assert [snp for frame in error.value.annotation_drops.frames() for snp in frame.SNP] == ["duplicated"] * 2
        assert not list(workspace.path.glob("source-*"))


def test_second_source_pass_failure_closes_handles_and_cleans_owned_scratch(tmp_path, monkeypatch):
    path = tmp_path / "base.annot"
    path.write_text("CHR SNP POS base\n1 rs1 1 1\n2 rs2 2 2\n")
    handles = []
    original = pd.read_csv

    def read_csv(*args, **kwargs):
        reader = original(*args, **kwargs)
        if "chunksize" in kwargs and "usecols" in kwargs:
            handles.extend(reader.handles.created_handles)

            def fail_read(*args, **kwargs):
                raise OSError("interrupted numeric read")

            monkeypatch.setattr(reader, "read", fail_read)
        return reader

    monkeypatch.setattr(pd, "read_csv", read_csv)
    output = tmp_path / "out"
    with pytest.raises(OSError, match="interrupted numeric read"):
        with AnnotationWorkspace(output) as workspace:
            prepare_annotation_sources(workspace, [path], [], mode="rsid", chunk_rows=1)
    assert handles and all(stream.closed for stream in handles)
    assert not list(output.iterdir())
