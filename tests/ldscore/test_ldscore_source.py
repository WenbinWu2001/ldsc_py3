"""Opening and validating LD-score directories with ``load_ldscore_from_dir``."""

import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd
import pytest

from ldsc.config import ConfigMismatchError, GlobalConfig
from ldsc.errors import LDSCInputError, LDSCInternalError
from tests.ldscore_helpers import write_ldscore_dir, write_one_snp_ldscore_dir


class LoadLDScoreFromDirTest(unittest.TestCase):
    def test_load_ldscore_from_dir_reads_metadata_and_parquet_files(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            identity = pd.DataFrame({"CHR": ["1"], "SNP": ["rs1"], "POS": [10]})
            write_ldscore_dir(
                tmpdir, identity.assign(regression_ld_scores=3.0, base=1.0), baseline_columns=["base"],
                query=identity.assign(query=2.0), query_columns=["query"],
                counts={"base": (5.0, 4.0), "query": (6.0, 5.0)}, snp_identifier="rsid",
            )
            result = load_ldscore_from_dir(str(tmpdir))
            query = result.read_queries(["query"])

        self.assertEqual(result.baseline_table.columns.tolist(), ["CHR", "SNP", "POS", "regression_ld_scores", "base"])
        self.assertEqual(query.columns.tolist(), ["CHR", "SNP", "POS", "query"])
        self.assertEqual(result.count_records[0]["column"], "base")
        self.assertEqual(result.baseline_table["SNP"].tolist(), ["rs1"])
        self.assertEqual(result.config_snapshot, GlobalConfig(snp_identifier="rsid"))

    def test_load_ldscore_from_dir_rejects_baseline_without_pos(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            write_one_snp_ldscore_dir(tmpdir, include_query=False)
            baseline = tmpdir / "ldscore.baseline.parquet"
            pd.read_parquet(baseline).drop(columns="POS").to_parquet(baseline, index=False)

            with self.assertRaisesRegex(LDSCInternalError, "baseline_table is missing required columns.*POS"):
                load_ldscore_from_dir(str(tmpdir))

    def test_load_ldscore_from_dir_rejects_missing_snp_identifier(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            write_one_snp_ldscore_dir(tmpdir, include_query=False)
            metadata_path = tmpdir / "metadata.json"
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            del metadata["snp_identifier"]
            metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

            with self.assertRaisesRegex(LDSCInputError, "snp_identifier=.*invalid"):
                load_ldscore_from_dir(str(tmpdir))

    def test_load_ldscore_from_dir_rejects_invalid_minimal_identity(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            write_one_snp_ldscore_dir(tmpdir, include_query=False)
            metadata_path = tmpdir / "metadata.json"
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            metadata["snp_identifier"] = "bad"
            metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

            with self.assertRaisesRegex(LDSCInputError, "snp_identifier=.*invalid"):
                load_ldscore_from_dir(str(tmpdir))

    def test_load_ldscore_from_dir_rejects_allele_aware_baseline_without_alleles(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            write_one_snp_ldscore_dir(tmpdir, include_query=True)
            pd.read_parquet(tmpdir / "ldscore.baseline.parquet").drop(columns=["A1", "A2"]).to_parquet(
                tmpdir / "ldscore.baseline.parquet",
                index=False,
            )
            metadata_path = tmpdir / "metadata.json"
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            metadata["snp_identifier"] = "rsid_allele_aware"
            metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

            with self.assertRaisesRegex(LDSCInputError, "baseline_table.*allele columns are missing"):
                load_ldscore_from_dir(str(tmpdir))

    def test_load_ldscore_from_dir_rejects_allele_aware_query_without_alleles(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            write_one_snp_ldscore_dir(tmpdir, include_query=True)
            pd.read_parquet(tmpdir / "ldscore.query.parquet").drop(columns=["A1", "A2"]).to_parquet(
                tmpdir / "ldscore.query.parquet",
                index=False,
            )
            metadata_path = tmpdir / "metadata.json"
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            metadata["snp_identifier"] = "rsid_allele_aware"
            metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

            with self.assertRaisesRegex(LDSCInputError, "query_table.*allele columns are missing"):
                load_ldscore_from_dir(str(tmpdir))

    def test_load_ldscore_from_dir_rejects_snp_identifier_override_mismatch(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            write_one_snp_ldscore_dir(tmpdir, include_query=False)

            with self.assertRaisesRegex(ConfigMismatchError, "snp_identifier mismatch"):
                load_ldscore_from_dir(str(tmpdir), snp_identifier="chr_pos")


def _edit_metadata(root, edit):
    """Apply ``edit`` to the parsed ``metadata.json`` of the LD-score directory ``root`` and save it."""
    path = root / "metadata.json"
    metadata = json.loads(path.read_text(encoding="utf-8"))
    edit(metadata)
    path.write_text(json.dumps(metadata), encoding="utf-8")


def _quantile_ldscores(tmp_path):
    from tests.quantile_fixtures import write_quantile_ldscores
    return write_quantile_ldscores(tmp_path)


def test_manifest_reads_metadata_and_the_overlap_without_the_tables(tmp_path, monkeypatch):
    from ldsc import load_ldscore_from_dir
    from ldsc._ldscore.source import read_ldscore_manifest

    root = _quantile_ldscores(tmp_path)
    loaded = load_ldscore_from_dir(str(root))
    monkeypatch.setattr(pd, "read_parquet", lambda *args, **kwargs: pytest.fail("the manifest read a table"))
    manifest = read_ldscore_manifest(root)

    assert manifest.snp_identifier == "rsid" and manifest.genome_build is None
    assert manifest.baseline_columns == ["base", "cont"] and manifest.query_columns == []
    assert manifest.count_records == loaded.count_records
    assert manifest.count_config == {"common_reference_snp_maf_min": 0.05, "common_reference_snp_maf_operator": ">="}
    assert manifest.annotation_types == {"base": "binary", "cont": "quantitative"}
    assert manifest.overlap_config["total_common_reference_snps"] == 4.0
    monkeypatch.undo()
    overlap = manifest.load_overlap()
    pd.testing.assert_frame_equal(overlap.baseline_block_common, loaded.overlap.baseline_block_common)
    assert overlap.total_all_reference_snps == loaded.overlap.total_all_reference_snps


@pytest.mark.parametrize("content", ["{not json", "[]"])
def test_unreadable_metadata_is_an_input_error(tmp_path, content):
    from ldsc import load_ldscore_from_dir

    root = write_one_snp_ldscore_dir(tmp_path / "ldscores", include_query=False)
    (root / "metadata.json").write_text(content, encoding="utf-8")

    with pytest.raises(LDSCInputError, match="metadata.json.*JSON object"):
        load_ldscore_from_dir(str(root))


@pytest.mark.parametrize("name", ["ldscore.baseline.parquet", "ldscore.overlap.parquet"])
def test_a_missing_declared_table_is_an_input_error(tmp_path, name):
    from ldsc import load_ldscore_from_dir

    root = _quantile_ldscores(tmp_path)
    (root / name).unlink()

    with pytest.raises(LDSCInputError, match=f"^Could not load LD-score directory .*{name}.*missing.*troubleshooting"):
        load_ldscore_from_dir(str(root))


@pytest.mark.parametrize("token", ["../ldscore.overlap.parquet", "absolute"])
def test_a_declared_table_outside_the_directory_is_rejected(tmp_path, token):
    from ldsc import load_ldscore_from_dir

    root = _quantile_ldscores(tmp_path)
    (root / "ldscore.overlap.parquet").rename(tmp_path / "ldscore.overlap.parquet")
    token = str(tmp_path / "ldscore.overlap.parquet") if token == "absolute" else token
    _edit_metadata(root, lambda metadata: metadata["files"].update(overlap=token))

    with pytest.raises(LDSCInputError, match="^Could not load LD-score directory .*files.overlap`=.* points outside the directory.*troubleshooting"):
        load_ldscore_from_dir(str(root))


def test_declared_files_may_be_symbolic_links_to_shared_storage(tmp_path):
    # Large tables linked from shared storage are inside the directory by name, as before the manifest read.
    from ldsc import load_ldscore_from_dir

    root = _quantile_ldscores(tmp_path)
    expected = load_ldscore_from_dir(str(root))
    shared = tmp_path / "shared"
    (shared / "diagnostics").mkdir(parents=True)
    for name in ("ldscore.baseline.parquet", "ldscore.overlap.parquet"):
        (root / name).rename(shared / name)
        (root / name).symlink_to(shared / name)
    (shared / "diagnostics" / "notes.tsv").write_text("x\n", encoding="utf-8")
    (root / "diagnostics").symlink_to(shared / "diagnostics", target_is_directory=True)
    _edit_metadata(root, lambda metadata: metadata["files"].update(notes="diagnostics/notes.tsv"))

    loaded = load_ldscore_from_dir(str(root))

    pd.testing.assert_frame_equal(loaded.baseline_table, expected.baseline_table)
    pd.testing.assert_frame_equal(loaded.overlap.baseline_block_all, expected.overlap.baseline_block_all)


@pytest.mark.parametrize("name", ["ldscore.query.parquet", "ldscore.overlap.parquet"])
def test_an_unreadable_parquet_footer_is_an_input_error(tmp_path, name):
    from ldsc import load_ldscore_from_dir

    root = (write_one_snp_ldscore_dir(tmp_path / "ldscores", include_query=True)
            if name == "ldscore.query.parquet" else _quantile_ldscores(tmp_path))
    (root / name).write_text("not parquet", encoding="utf-8")

    with pytest.raises(LDSCInputError, match=f"^Could not read LD-score file .*{name}.*troubleshooting"):
        load_ldscore_from_dir(str(root))


@pytest.mark.parametrize("overlap_config", [
    None, {}, [4.0],
    {"total_common_reference_snps": 4.0, "common_maf_min": 0.05, "common_maf_operator": ">="},
    {"total_all_reference_snps": "4", "total_common_reference_snps": 4.0,
     "common_maf_min": 0.05, "common_maf_operator": ">="},
    {"total_all_reference_snps": 4.0, "total_common_reference_snps": "x",
     "common_maf_min": 0.05, "common_maf_operator": ">="},
    # JSON NaN and Infinity decode to floats; a total must be finite and the common total within the all total.
    *({"total_all_reference_snps": total_all, "total_common_reference_snps": total_common,
       "common_maf_min": 0.05, "common_maf_operator": ">="}
      for total_all, total_common in [(float("nan"), 4.0), (float("inf"), 4.0), (0.0, None), (-4.0, None),
                                      (4.0, float("nan")), (4.0, -1.0), (4.0, 5.0)]),
])
def test_a_declared_overlap_needs_a_usable_overlap_config(tmp_path, overlap_config):
    # The overlap file stores cross-products only; its reference-SNP totals live in overlap_config.
    from ldsc import load_ldscore_from_dir

    root = _quantile_ldscores(tmp_path)
    _edit_metadata(root, lambda metadata: metadata.update(overlap_config=overlap_config))

    with pytest.raises(LDSCInputError,
                       match="^Could not load LD-score directory .*overlap_config.*total_all_reference_snps.*troubleshooting"):
        load_ldscore_from_dir(str(root))


@pytest.mark.parametrize("overlap_config", [[4.0], "x", 0, False])
def test_overlap_config_without_an_overlap_file_must_be_an_object_or_null(tmp_path, overlap_config):
    from ldsc import load_ldscore_from_dir

    root = write_one_snp_ldscore_dir(tmp_path / "ldscores", include_query=False)
    _edit_metadata(root, lambda metadata: metadata.update(overlap_config=overlap_config))

    with pytest.raises(LDSCInputError,
                       match="^Could not load LD-score directory .*`overlap_config`.*JSON object or null.*troubleshooting"):
        load_ldscore_from_dir(str(root))


@pytest.mark.parametrize("overlap_config", [None, {"total_all_reference_snps": float("nan")}])
def test_load_overlap_fails_closed_without_usable_totals(tmp_path, overlap_config):
    # A manifest built by hand skips read_ldscore_manifest's checks; the overlap read checks its totals itself.
    from dataclasses import replace
    from ldsc._ldscore.source import read_ldscore_manifest

    manifest = replace(read_ldscore_manifest(_quantile_ldscores(tmp_path)), overlap_config=overlap_config)

    with pytest.raises(LDSCInputError, match="^Could not load LD-score directory .*overlap_config.*total_all_reference_snps"):
        manifest.load_overlap()


@pytest.mark.parametrize("count_config", [None, {}, *({"common_reference_snp_maf_min": threshold,
                                                       "common_reference_snp_maf_operator": ">="}
                                                      for threshold in ["x", float("nan"), float("inf")])])
def test_a_missing_or_unusable_count_config_is_rejected(tmp_path, count_config):
    from ldsc import load_ldscore_from_dir

    root = write_one_snp_ldscore_dir(tmp_path / "ldscores", include_query=False)
    _edit_metadata(root, lambda metadata: metadata.update(count_config=count_config))

    with pytest.raises(LDSCInputError, match="count_config"):
        load_ldscore_from_dir(str(root))


@pytest.mark.parametrize("missing", [None, "index_snp_identifier", "index_genome_build"])
def test_an_index_id_needs_the_rest_of_the_index_provenance(tmp_path, missing):
    from ldsc import load_ldscore_from_dir

    root = write_one_snp_ldscore_dir(tmp_path / "ldscores", include_query=False)
    provenance = {"index_id": "0" * 64, "index_snp_identifier": "rsid", "index_genome_build": "hg19"}
    provenance.pop(missing, None)
    _edit_metadata(root, lambda metadata: metadata.update(provenance))

    if missing is None:
        assert load_ldscore_from_dir(str(root)).index_provenance == provenance
        return
    with pytest.raises(LDSCInputError, match=f"^Could not load LD-score directory .*index_id.*{missing}.*troubleshooting"):
        load_ldscore_from_dir(str(root))
