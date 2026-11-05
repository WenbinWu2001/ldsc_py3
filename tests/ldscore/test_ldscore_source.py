"""Opening and validating LD-score directories with ``load_ldscore_from_dir``."""

import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from ldsc.config import ConfigMismatchError, GlobalConfig
from ldsc.errors import LDSCInputError, LDSCInternalError


def write_ldscore_dir(root: Path, *, include_query: bool = True) -> Path:
    """Write a one-SNP rsID LD-score directory, with or without one query column."""
    root.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "CHR": ["1"],
            "SNP": ["rs1"],
            "POS": [10],
            "A1": ["A"],
            "A2": ["C"],
            "regression_ld_scores": [2.0],
            "base": [1.0],
        }
    ).to_parquet(root / "ldscore.baseline.parquet", index=False)
    files = {"baseline": "ldscore.baseline.parquet"}
    query_columns = []
    counts = [
        {
            "group": "baseline",
            "column": "base",
            "all_reference_snp_count": 5.0,
            "common_reference_snp_count": 5.0,
        }
    ]
    if include_query:
        pd.DataFrame(
            {
                "CHR": ["1"],
                "SNP": ["rs1"],
                "POS": [10],
                "A1": ["A"],
                "A2": ["C"],
                "query": [2.0],
            }
        ).to_parquet(root / "ldscore.query.parquet", index=False)
        files["query"] = "ldscore.query.parquet"
        query_columns = ["query"]
        counts.append(
            {
                "group": "query",
                "column": "query",
                "all_reference_snp_count": 6.0,
                "common_reference_snp_count": 6.0,
            }
        )
    (root / "metadata.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "artifact_type": "ldscore",
                "files": files,
                "snp_identifier": "rsid",
                "genome_build": "hg38",
                "chromosomes": ["1"],
                "baseline_columns": ["base"],
                "query_columns": query_columns,
                "query_batches": [{"file": "ldscore.query.parquet", "query_columns": ["query"], "row_groups": [{"chrom": "1", "row_group_index": 0, "row_offset": 0, "n_rows": 1}]}] if include_query else [],
                "baseline_row_groups": [{"chrom": "1", "row_group_index": 0, "row_offset": 0, "n_rows": 1}],
                "counts": counts,
            }
        ),
        encoding="utf-8",
    )
    return root


class LoadLDScoreFromDirTest(unittest.TestCase):
    def test_load_ldscore_from_dir_reads_metadata_and_parquet_files(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "POS": [10],
                    "regression_ld_scores": [3.0],
                    "base": [1.0],
                }
            ).to_parquet(tmpdir / "ldscore.baseline.parquet", index=False)
            pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "POS": [10],
                    "query": [2.0],
                }
            ).to_parquet(tmpdir / "ldscore.query.parquet", index=False)
            (tmpdir / "metadata.json").write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "artifact_type": "ldscore",
                        "files": {"baseline": "ldscore.baseline.parquet", "query": "ldscore.query.parquet"},
                        "snp_identifier": "rsid",
                        "genome_build": "hg38",
                        "chromosomes": ["1"],
                        "baseline_columns": ["base"],
                        "query_columns": ["query"],
                        "query_batches": [{"file": "ldscore.query.parquet", "query_columns": ["query"], "row_groups": [{"chrom": "1", "row_group_index": 0, "row_offset": 0, "n_rows": 1}]}],
                        "baseline_row_groups": [{"chrom": "1", "row_group_index": 0, "row_offset": 0, "n_rows": 1}],
                        "counts": [
                            {
                                "group": "baseline",
                                "column": "base",
                                "all_reference_snp_count": 5.0,
                                "common_reference_snp_count": 4.0,
                            },
                            {
                                "group": "query",
                                "column": "query",
                                "all_reference_snp_count": 6.0,
                                "common_reference_snp_count": 5.0,
                            },
                        ],
                    }
                ),
                encoding="utf-8",
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
            pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "regression_ld_scores": [3.0],
                    "base": [1.0],
                }
            ).to_parquet(tmpdir / "ldscore.baseline.parquet", index=False)
            (tmpdir / "metadata.json").write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "artifact_type": "ldscore",
                        "files": {"baseline": "ldscore.baseline.parquet"},
                        "snp_identifier": "rsid",
                        "genome_build": "hg38",
                        "chromosomes": ["1"],
                        "baseline_columns": ["base"],
                        "query_columns": [],
                        "query_batches": [],
                        "baseline_row_groups": [{"chrom": "1", "row_group_index": 0, "row_offset": 0, "n_rows": 1}],
                        "counts": [
                            {
                                "group": "baseline",
                                "column": "base",
                                "all_reference_snp_count": 5.0,
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(LDSCInternalError, "baseline_table is missing required columns.*POS"):
                load_ldscore_from_dir(str(tmpdir))

    def test_load_ldscore_from_dir_rejects_missing_snp_identifier(self):
        from ldsc import load_ldscore_from_dir

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            write_ldscore_dir(tmpdir, include_query=False)
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
            write_ldscore_dir(tmpdir, include_query=False)
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
            write_ldscore_dir(tmpdir, include_query=True)
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
            write_ldscore_dir(tmpdir, include_query=True)
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
            write_ldscore_dir(tmpdir, include_query=False)

            with self.assertRaisesRegex(ConfigMismatchError, "snp_identifier mismatch"):
                load_ldscore_from_dir(str(tmpdir), snp_identifier="chr_pos")
