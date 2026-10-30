"""Direct writing streams chromosome scores without changing the fitted universe."""

from dataclasses import replace
from pathlib import Path
import multiprocessing
import weakref
import json

import numpy as np
import pandas as pd
import pytest

from ldsc import AnnotationBundle, LDScoreCalculator
from ldsc._ldscore import calculator as ldscore_calculator
from ldsc._ldscore.outputs import LDScoreOutputConfig
from tests.ldscore.test_ldscore_query_batches import prepared_case
from tests.ldscore.test_gene_query_integration import inputs, direct_args


@pytest.mark.parametrize("width", [1, 2, 10])
@pytest.mark.parametrize("empty_first_chromosome", [False, True])
def test_direct_global_filtering_matches_materialized_output(tmp_path, width, empty_first_chromosome):
    from ldsc._annotation.query_outcomes import QueryAnnotationStatus
    from tests.ldscore.test_ldscore_parallelism import _PANEL_CHROMS, _write_index_panel

    metadata, baseline, _, config, panel, numerical = prepared_case(tmp_path)
    # Zero off-diagonals make each LD score equal its annotation value.
    _write_index_panel(tmp_path / "panel", {chrom: (rows, [(0, 1, 0.)])
                                            for chrom, (rows, _) in _PANEL_CHROMS.items()})
    queries = pd.DataFrame({"constant": np.ones(7), "zero": np.zeros(7),
                            "between": [1., 1., 1., 1., 2., 2., 2.],
                            "within": [0., 1., 0., 1., 0., 1., 0.]})
    restriction = {"rs5", "rs6"} if empty_first_chromosome else None
    statuses = tuple(QueryAnnotationStatus(name, name, "bed", "ok") for name in queries)
    with AnnotationBundle.from_frames(metadata, baseline, queries, config_snapshot=config) as bundle:
        bundle = replace(bundle, query_statuses=statuses)
        calculator = LDScoreCalculator()
        expected = calculator.run(bundle, panel, numerical, config, regression_snps=restriction)
        saved = calculator.run(bundle, panel, replace(numerical, query_batch_size=width), config,
            regression_snps=restriction, output_config=LDScoreOutputConfig(tmp_path / "out"))
    retained = ["within"] if empty_first_chromosome else ["between", "within"]
    assert saved.query_columns == retained
    assert saved.query_statuses == expected.query_statuses
    assert [status.reason for status in saved.query_statuses] == [
        "zero_variance_ld_scores", "zero_annotation_snps",
        "zero_variance_ld_scores" if empty_first_chromosome else "", ""]
    pd.testing.assert_frame_equal(saved.read_queries(retained), expected.query_table, check_exact=True)
    assert saved.count_records == expected.count_records
    assert saved.snp_universe_policy == expected.snp_universe_policy
    assert saved.snp_universe_policy["ld_reference_universe"]["retained_snp_count"] == 7
    assert saved.count_records[-1]["all_reference_snp_count"] == 3
    assert not list((tmp_path / "out").glob(".ldsc-*"))


def test_direct_parent_releases_and_deletes_each_consumed_fragment(tmp_path, monkeypatch):
    from ldsc._ldscore import fragments

    metadata, baseline, queries, config, panel, numerical = prepared_case(tmp_path)
    read = fragments.pd.read_parquet
    live, consumed = [], []

    def measured(path, *args, **kwargs):
        if Path(path).parent.name.startswith(".ldsc-fragments-"):
            assert all(reference() is None for reference in live)
            assert all(not previous.exists() for previous in consumed)
            frame = read(path, *args, **kwargs)
            assert frame.CHR.nunique() == 1
            live.append(weakref.ref(frame))
            consumed.append(Path(path))
            return frame
        return read(path, *args, **kwargs)

    monkeypatch.setattr(fragments.pd, "read_parquet", measured)
    with AnnotationBundle.from_frames(metadata, baseline, queries, config_snapshot=config) as bundle:
        LDScoreCalculator().run(bundle, panel, replace(numerical, query_batch_size=2), config,
                               output_config=LDScoreOutputConfig(tmp_path / "out"))
    assert len(consumed) == 4
    assert all(reference() is None for reference in live)
    assert all(not path.exists() for path in consumed)


def _fail_last_batch(*args, **kwargs):
    if kwargs["batch_label"] == "batch 2/2" and args[0] == "2":
        raise OSError("injected final-batch chromosome failure")
    return ldscore_calculator._compute_one_chromosome(*args, **kwargs)


@pytest.mark.parametrize("failure", ["worker", "alignment", "writer"])
def test_direct_final_batch_failure_removes_private_outputs(tmp_path, monkeypatch, failure):
    from ldsc._ldscore import fragments
    from ldsc.errors import LDSCInputError

    metadata, baseline, queries, config, panel, numerical = prepared_case(tmp_path)
    output = tmp_path / "out"
    original_write = fragments.write_query_fragments
    writes = []

    def write(result, baseline_table, path, compression):
        assert not (output / "metadata.json").exists()
        assert not list(output.glob("ldscore.*.parquet"))
        writes.append(result.query_columns)
        if len(writes) == 2:
            if failure == "writer":
                raise OSError("injected final-batch publication failure")
            if failure == "alignment":
                bad = result.query_fragments[-1]
                frame = pd.read_parquet(bad)
                frame.loc[0, "SNP"] = "wrong-row"
                frame.to_parquet(bad, index=False)
        return original_write(result, baseline_table, path, compression)

    if failure == "worker":
        monkeypatch.setattr(ldscore_calculator, "_compute_one_chromosome", _fail_last_batch)
    monkeypatch.setattr(fragments, "write_query_fragments", write)
    before = {child.pid for child in multiprocessing.active_children()}
    with AnnotationBundle.from_frames(metadata, baseline, queries, config_snapshot=config) as bundle:
        with pytest.raises((OSError, LDSCInputError), match="injected final-batch|SNP row mismatch"):
            LDScoreCalculator().run(bundle, panel, replace(numerical, query_batch_size=2, workers=2), config,
                                   output_config=LDScoreOutputConfig(output))
    assert len(writes) == (1 if failure == "worker" else 2)
    assert {child.pid for child in multiprocessing.active_children()} == before
    assert not (output / "metadata.json").exists()
    assert not list(output.glob("ldscore.*.parquet"))
    assert not list(output.glob(".ldsc-*"))
    assert (output / "diagnostics/chromosome-logs/chr2.log").is_file()


def test_real_plink_gene_batches_preserve_all_metadata_and_diagnostics(tmp_path, monkeypatch, inputs):
    from ldsc._ldscore.calculator import run_ldscore_from_args

    _, annotation_path, _, _ = inputs
    annotations = pd.read_csv(annotation_path, sep="\t")
    pd.concat([annotations, annotations.iloc[[0]]], ignore_index=True).to_csv(
        annotation_path, sep="\t", index=False)
    mixed, second, empty, zero, control = [tmp_path / f"{name}.txt"
        for name in ("mixed", "second", "empty", "zero", "control")]
    for path, text in [(mixed, "G1\nONE\nG0\nUNKNOWN\n"), (second, "G2\n"),
                       (empty, ""), (zero, "G0\n"), (control, "G2\n")]:
        path.write_text(text)

    def run(name):
        args = direct_args(tmp_path, inputs, [mixed, second, empty, zero],
                           policy="resolved-only", control=control, name=name)
        args.query_batch_size = 1
        args.export_ref_metadata = True
        return run_ldscore_from_args(args)

    from contextlib import contextmanager
    from ldsc._ldscore.calculator import _DirectProducer
    streamed_batch = _DirectProducer.open_batch

    @contextmanager
    def materialized(self, *args, **kwargs):
        # Keep the original materialized assembly as the transport comparator.
        with streamed_batch(self, *args, **kwargs) as batch:
            yield replace(batch, fragment_dir=None)

    with monkeypatch.context() as patch:
        patch.setattr(_DirectProducer, "open_batch", materialized)
        expected = run("materialized")
    actual = run("streamed")

    assert actual.query_columns == expected.query_columns == ["mixed", "second"]
    assert actual.query_statuses == expected.query_statuses
    assert len(actual.query_batches) == len(expected.query_batches) == 2
    left, right = tmp_path / "materialized", tmp_path / "streamed"
    assert json.loads((left / "metadata.json").read_text()) == json.loads((right / "metadata.json").read_text())
    paths = {p.relative_to(left) for p in left.rglob("*") if p.is_file() and p.suffix != ".log"}
    assert paths == {p.relative_to(right) for p in right.rglob("*") if p.is_file() and p.suffix != ".log"}
    assert Path("diagnostics/gene_list_audit.tsv.gz") in paths
    assert Path("diagnostics/chromosome_scope.json") in paths
    assert Path("ref_metadata/chr1_meta.tsv.gz") in paths
    assert Path("diagnostics/dropped_snps/chr1_dropped.tsv.gz") in paths
    for path in sorted(paths):
        if path.suffix == ".parquet":
            import pyarrow.parquet as pq
            pd.testing.assert_frame_equal(pd.read_parquet(left / path), pd.read_parquet(right / path), check_exact=True)
            assert pq.read_schema(left / path) == pq.read_schema(right / path)
        elif path.suffix == ".json":
            assert json.loads((left / path).read_text()) == json.loads((right / path).read_text())
        else:
            pd.testing.assert_frame_equal(pd.read_csv(left / path, sep="\t"), pd.read_csv(right / path, sep="\t"),
                                          check_exact=True)
