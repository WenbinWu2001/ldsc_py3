"""Opening and validating LD-score directories with ``load_ldscore_from_dir``."""

import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

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
