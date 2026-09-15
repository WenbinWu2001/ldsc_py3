"""Annotation encoding preserves logical identities, order, and numeric values."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ldsc import AnnotationBundle
from ldsc._annotation_sources import prepare_annotation_sources
from ldsc._annotation_storage import AnnotationWorkspace


@pytest.mark.parametrize("sharded", [False, True])
def test_supplied_annotations_use_global_exact_binary_storage(tmp_path, monkeypatch, sharded):
    metadata = pd.DataFrame({"CHR": [2, 1] * 11, "POS": np.repeat(np.arange(11) + 1, 2),
                             "SNP": [f"rs{i}" for i in range(22)]})
    metadata.loc[[0, 21], "SNP"] = "duplicate"
    baseline = metadata.assign(continuous=-np.arange(22) / 4, binary=np.arange(22) % 2,
                               late=1., near_one=1.00000001)
    baseline.loc[20, "late"] = .5
    query = metadata.assign(query_binary=(np.arange(22) % 3 == 0).astype(int),
                            query_continuous=np.arange(22) / 8)
    sources = []
    for role, frame in [("baseline", baseline), ("query", query)]:
        paths = []
        for label, part in frame.groupby("CHR", sort=False) if sharded else [("all", frame)]:
            path = tmp_path / f"{role}.{label}.annot.gz"
            part.to_csv(path, sep="\t", index=False)
            paths.append(path)
        sources.append(paths)
    reads, numeric_spools = [], []
    original_read, original_open = pd.read_csv, Path.open
    def read_csv(path, *args, **kwargs):
        if "chunksize" in kwargs and Path(path) in sources[0] + sources[1]:
            reads.append(Path(path))
        return original_read(path, *args, **kwargs)
    def open_path(path, mode="r", *args, **kwargs):
        if path.suffix == ".bin" and any(flag in mode for flag in "wax+"):
            numeric_spools.append(path)
        return original_open(path, mode, *args, **kwargs)
    monkeypatch.setattr(pd, "read_csv", read_csv)
    monkeypatch.setattr(Path, "open", open_path)
    names = ["continuous", "binary", "late", "near_one", "query_binary", "query_continuous"]
    expected = pd.concat([baseline.iloc[:, 3:], query.iloc[:, 3:]], axis=1).astype(np.float32)
    with AnnotationWorkspace(tmp_path / "out") as workspace:
        prepared = prepare_annotation_sources(workspace, *sources, mode="rsid", chunk_rows=3)
        bundle = AnnotationBundle(prepared.shards, list(prepared.baseline_columns),
                                  list(prepared.query_columns), workspace)
        for chrom in bundle.chromosomes:
            keep = metadata.CHR.eq(int(chrom)) & metadata.SNP.ne("duplicate")
            shard = bundle.shard(chrom)
            types = {name: store.dtype for store in shard.stores for name in store.columns}
            assert types == {name: np.dtype(bool if name in {"binary", "query_binary"} else np.float32) for name in names}
            np.testing.assert_array_equal(bundle.read(chrom), expected.loc[keep])
            selected = ["query_continuous", "binary", "continuous", "query_binary"]
            np.testing.assert_array_equal(bundle.read(chrom, rows=[9, 0, 9], columns=selected),
                                          expected.loc[keep, selected].iloc[[9, 0, 9]])
            assert bundle.read(chrom, columns=["binary", "query_binary"]).dtype == np.float32
            with pytest.raises(TypeError):
                shard.read()
        assert prepared.drops.n_rows == 2
    assert all(reads.count(path) == 2 for paths in sources for path in paths)
    assert not numeric_spools


@pytest.mark.parametrize("kind", ["gene", "bed"])
def test_generated_annotations_are_packed_at_projection(tmp_path, kind):
    from ldsc import AnnotationBuilder, AnnotationBuildConfig, GlobalConfig
    from tests.test_annotate_streaming import inputs

    baseline, genes, bed, catalog = inputs(tmp_path)
    options = dict(query_annot_gene_list_sources=(genes,), gene_coordinate_file=catalog, padding_bp=0) if kind == "gene" else dict(query_annot_bed_sources=(bed,))
    spec = AnnotationBuildConfig(baseline_annot_sources=tuple(baseline), **options)
    with AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), spec, projection_genome_build="hg19").run(output_dir=tmp_path / "out") as bundle:
        np.testing.assert_array_equal(bundle.read("1", columns=["pathway"]), [[0], [1], [1], [0]])
        assert bundle.read("1", columns=["pathway"]).dtype == np.float32
        for chrom, expected in [("1", [[6]]), ("2", [[0]])]:
            store = next(store for store in bundle.shard(chrom).stores if "pathway" in store.columns)
            assert store.dtype == bool and store.bitorder == "little"
            np.testing.assert_array_equal(np.load(store.path), expected)


def test_from_frames_packs_binary_columns_without_writing_and_preserves_public_reads(tmp_path):
    from ldsc import GlobalConfig

    metadata = pd.DataFrame({"CHR": [1] * 11, "POS": [11, *range(1, 11)],
                             "SNP": [f"rs{i}" for i in range(11)]})
    baseline = pd.DataFrame({"near_one": [1.00000001] * 11, "binary": np.arange(11) % 2})
    query = pd.DataFrame({"query_binary": np.arange(11) % 3 == 0, "signed": -np.arange(11) / 4})
    expected = pd.concat([baseline, query], axis=1).astype(np.float32).iloc[[*range(1, 11), 0]]
    with AnnotationBundle.from_frames(metadata, baseline, query, config_snapshot=GlobalConfig(snp_identifier="rsid")) as bundle:
        shard = bundle.shard("1")
        assert shard.binary_values.nbytes == 4
        assert shard.values.nbytes == 88
        np.testing.assert_array_equal(bundle.read("1"), expected)
        np.testing.assert_array_equal(bundle.read("1", rows=[10, 0, 10], columns=["signed", "binary"]),
                                      expected[["signed", "binary"]].iloc[[10, 0, 10]])
        assert bundle.read("1", columns=["binary"]).dtype == np.float32
        assert bundle.read("1", rows=[], columns=[]).shape == (0, 0)
    assert not list(tmp_path.iterdir())
