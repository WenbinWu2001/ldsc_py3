"""Known-value and failure checks for the benchmark-only parallel adapter."""

import argparse
import json
from pathlib import Path
import sqlite3

import numpy as np
import pandas as pd

from annotation_preparation_parallel import ReadOnlyIdentity, parallel_prototype
from ldsc._annotation_identity import DiskIdentityIndex
from ldsc._annotation_sources import prepare_annotation_sources
from ldsc._annotation_storage import AnnotationWorkspace


def check(root):
    root.mkdir(parents=True, exist_ok=False)
    rows = {
        1: [("dup", 10, "A", "C", 1, 11), ("multi", 20, "A", "C", 2, 12),
            ("good1", 30, "A", "C", -1.25, .125), ("bad", 40, "N", "C", 4, 14)],
        2: [("good2", 90, "A", "C", 2.5, 1.5), ("dup", 50, "G", "T", 5, 15),
            ("multi", 60, "A", "G", 6, 16), ("good22", 5, "A", "C", 3.25, -.5)],
        3: [("good3", 10, None, None, 0, 2)],
        4: [("dup", 10, "A", "C", 9, 19)],
    }
    baselines, queries = [], []
    for chrom in [3, 1, 4, 2]:
        frame = pd.DataFrame(rows[chrom], columns=["SNP", "POS", "A1", "A2", "base", "q"]).assign(CHR=chrom)
        base, query = root / f"base.{chrom}.annot", root / f"query.{chrom}.annot"
        frame[["CHR", "SNP", "POS", "base"]].to_csv(base, sep="\t", index=False)
        columns = ["CHR", "SNP", "POS", "q"] if chrom == 3 else ["CHR", "SNP", "POS", "A1", "A2", "q"]
        frame[columns].to_csv(query, sep="\t", index=False)
        baselines.append(base)
        queries.append(query)
    expected = {"1": (["good1"], [[-1.25, .125]]),
                "2": (["good2", "good22"], [[2.5, 1.5], [3.25, -.5]]), "3": (["good3"], [[0, 2]])}
    outcomes = []
    for executor in ["thread", "process"]:
        for workers, selected in [(1, None), (2, None), (4, None), (8, None), (2, "2")]:
            label = f"{executor}-{workers}-chrom-{selected}"
            with AnnotationWorkspace(root / label) as workspace, parallel_prototype(workers, executor) as report:
                result = prepare_annotation_sources(workspace, baselines, queries,
                                                    mode="rsid_allele_aware", chunk_rows=2, chrom=selected)
                assert list(result.shards) == (["2"] if selected else ["1", "2", "3"])
                assert result.scope_chromosomes == ("1", "2", "3", "4")
                for chrom, shard in result.shards.items():
                    assert shard.metadata().SNP.tolist() == expected[chrom][0]
                    np.testing.assert_array_equal(shard.read(), expected[chrom][1])
                    assert shard.read().dtype == np.float32
                drops = pd.concat(result.drops.frames(), ignore_index=True)
                assert drops.SNP.tolist() == ["bad", "multi", "multi", "dup", "dup", "dup"]
                assert drops.CHR.tolist() == ["1", "1", "2", "1", "2", "4"]
                assert drops.reason.tolist() == ["invalid_allele"] + ["multi_allelic_base_key"] * 2 + ["duplicate_identity"] * 3
                assert report["effective_workers"] == min(workers, 4)
            assert not list((root / label).glob(".ldsc-annotation-*"))
            outcomes.append(label)
        events = root / f"{executor}-failure-events"
        events.mkdir()
        owner = AnnotationWorkspace(root / f"{executor}-failure")
        try:
            with owner, parallel_prototype(2, executor, fail_group=0, events=events):
                prepare_annotation_sources(owner, baselines, queries, mode="rsid_allele_aware", chunk_rows=2)
        except OSError as error:
            assert "injected chromosome worker failure" in str(error)
        else:
            raise AssertionError("Worker failure was not propagated")
        assert not owner.path.exists()
        started = {p.name.split("-")[1] for p in events.glob("started-*")}
        finished = {p.name.split("-")[1] for p in events.glob("finished-*")}
        assert started == finished and started
        assert all(p.read_text() == "True" for p in events.glob("finished-*"))
        outcomes.append(f"{executor}-failure-cleanup")
    whole_base, whole_query = root / "whole-base.annot", root / "whole-query.annot"
    pd.concat([pd.read_csv(root / f"base.{c}.annot", sep="\t") for c in [1, 2]]).to_csv(whole_base, sep="\t", index=False)
    pd.concat([pd.read_csv(root / f"query.{c}.annot", sep="\t") for c in [1, 2]]).to_csv(whole_query, sep="\t", index=False)
    with AnnotationWorkspace(root / "whole") as owner, parallel_prototype(8, "process") as report:
        result = prepare_annotation_sources(owner, [whole_base], [whole_query], mode="rsid_allele_aware", chunk_rows=2)
        assert not report["parallel_used"] and report["effective_workers"] == 1
        for chrom, shard in result.shards.items():
            np.testing.assert_array_equal(shard.read(), expected[chrom][1])
    outcomes.append("whole-genome-sequential-fallback")
    database = root / "readonly.sqlite"
    with DiskIdentityIndex(database, "rsid") as writer:
        writer.add(pd.DataFrame({"SNP": ["unique"], "CHR": ["1"], "POS": [1]}))
        writer._finish_transaction()
        with ReadOnlyIdentity(database, "rsid") as reader:
            try:
                reader.connection.execute("DELETE FROM identities")
            except sqlite3.OperationalError as error:
                assert "readonly" in str(error)
            else:
                raise AssertionError("Read-only connection permitted mutation")
    outcomes.append("readonly-sqlite")
    result = {"passed": len(outcomes), "checks": outcomes}
    (root / "checks.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    check(parser.parse_args().output)
