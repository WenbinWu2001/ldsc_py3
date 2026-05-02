from __future__ import annotations

from dataclasses import replace as dataclass_replace
from pathlib import Path
import json
import sys
import tempfile
import unittest
import warnings

import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc.config import GlobalConfig
from ldsc.ldscore_calculator import LDScoreResult
from ldsc.outputs import (
    LDScoreDirectoryWriter,
    LDScoreOutputConfig,
    PartitionedH2DirectoryWriter,
    PartitionedH2OutputConfig,
)


def make_split_ldscore_result(query: bool = True) -> LDScoreResult:
    baseline_table = pd.DataFrame(
        {
            "CHR": ["1", "1"],
            "SNP": ["rs1", "rs2"],
            "POS": [10, 20],
            "regr_weight": [1.5, 2.5],
            "base": [1.0, 2.0],
        }
    )
    query_table = None
    query_columns: list[str] = []
    if query:
        query_table = pd.DataFrame(
            {
                "CHR": ["1", "1"],
                "SNP": ["rs1", "rs2"],
                "POS": [10, 20],
                "query": [0.5, 0.7],
            }
        )
        query_columns = ["query"]
    count_records = [
        {
            "group": "baseline",
            "column": "base",
            "all_reference_snp_count": 10.0,
            "common_reference_snp_count": 8.0,
        }
    ]
    if query:
        count_records.append(
            {
                "group": "query",
                "column": "query",
                "all_reference_snp_count": 20.0,
                "common_reference_snp_count": 18.0,
            }
        )
    return LDScoreResult(
        baseline_table=baseline_table,
        query_table=query_table,
        count_records=count_records,
        baseline_columns=["base"],
        query_columns=query_columns,
        ld_reference_snps=frozenset(),
        ld_regression_snps=frozenset({"rs1", "rs2"}),
        chromosome_results=[],
        config_snapshot=GlobalConfig(genome_build="hg38", snp_identifier="rsid"),
    )


def make_multi_chrom_result(chromosomes: list[str] | None = None) -> LDScoreResult:
    """Three-chromosome result for row-group layout tests."""
    if chromosomes is None:
        chromosomes = ["1", "2", "22"]
    rows_per_chrom = {"1": 3, "2": 2, "22": 1}
    baseline_rows = []
    query_rows = []
    snp_idx = 1
    for chrom in chromosomes:
        n = rows_per_chrom.get(chrom, 2)
        for i in range(n):
            baseline_rows.append(
                {
                    "CHR": chrom,
                    "SNP": f"rs{snp_idx}",
                    "POS": (i + 1) * 100,
                    "regr_weight": 1.0,
                    "base": float(snp_idx),
                }
            )
            query_rows.append(
                {
                    "CHR": chrom,
                    "SNP": f"rs{snp_idx}",
                    "POS": (i + 1) * 100,
                    "query": float(snp_idx) * 0.5,
                }
            )
            snp_idx += 1
    baseline_table = pd.DataFrame(baseline_rows)
    query_table = pd.DataFrame(query_rows)
    return LDScoreResult(
        baseline_table=baseline_table,
        query_table=query_table,
        count_records=[],
        baseline_columns=["base"],
        query_columns=["query"],
        ld_reference_snps=frozenset(),
        ld_regression_snps=frozenset(),
        chromosome_results=[],
        config_snapshot=GlobalConfig(genome_build="hg38", snp_identifier="rsid"),
    )


class LDScoreDirectoryWriterTest(unittest.TestCase):
    def test_outputs_module_does_not_expose_prefix_based_artifact_pipeline(self):
        import ldsc.outputs as outputs

        removed_names = [
            "ArtifactConfig",
            "ArtifactOutputConfig",
            "ArtifactProducer",
            "OutputManager",
            "PostProcessor",
            "ResultFormatter",
            "ResultWriter",
        ]
        for name in removed_names:
            self.assertFalse(hasattr(outputs, name), f"{name} should not be exposed by ldsc.outputs")

    def test_writes_manifest_baseline_and_query_parquet(self):
        result = make_split_ldscore_result(query=True)
        writer = LDScoreDirectoryWriter()
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "ldscores"
            output_paths = writer.write(result, LDScoreOutputConfig(output_dir=output_dir))

            self.assertEqual(set(output_paths), {"manifest", "baseline", "query"})
            manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["format"], "ldsc.ldscore_result.v1")
            self.assertEqual(manifest["files"], {"baseline": "baseline.parquet", "query": "query.parquet"})
            self.assertEqual(manifest["baseline_columns"], ["base"])
            self.assertEqual(manifest["query_columns"], ["query"])
            self.assertEqual(manifest["counts"][1]["column"], "query")
            self.assertEqual(
                manifest["count_config"],
                {
                    "common_reference_snp_maf_min": 0.05,
                    "common_reference_snp_maf_operator": ">=",
                },
            )

            baseline = pd.read_parquet(output_dir / "baseline.parquet")
            query = pd.read_parquet(output_dir / "query.parquet")
            self.assertEqual(baseline.columns.tolist(), ["CHR", "SNP", "POS", "regr_weight", "base"])
            self.assertEqual(query.columns.tolist(), ["CHR", "SNP", "POS", "query"])

    def test_one_row_group_per_chromosome(self):
        import pyarrow.parquet as pq

        result = make_multi_chrom_result()
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "out"
            LDScoreDirectoryWriter().write(result, LDScoreOutputConfig(output_dir=output_dir))

            for fname in ("baseline.parquet", "query.parquet"):
                pf = pq.ParquetFile(output_dir / fname)
                self.assertEqual(pf.metadata.num_row_groups, 3, f"{fname}: expected 3 row groups")
                for i in range(pf.metadata.num_row_groups):
                    chroms = pf.read_row_group(i)["CHR"].unique().to_pylist()
                    self.assertEqual(len(chroms), 1, f"{fname} row group {i} mixes chromosomes: {chroms}")

    def test_chromosome_read_via_row_group_index_excludes_other_chromosomes(self):
        import pyarrow.parquet as pq

        result = make_multi_chrom_result()
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "out"
            LDScoreDirectoryWriter().write(result, LDScoreOutputConfig(output_dir=output_dir))

            manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
            rg_by_chrom = {e["chrom"]: e["row_group_index"] for e in manifest["baseline_row_groups"]}

            pf = pq.ParquetFile(output_dir / "baseline.parquet")
            df = pf.read_row_group(rg_by_chrom["1"]).to_pandas()
            self.assertTrue((df["CHR"] == "1").all())
            self.assertNotIn("2", df["CHR"].values)
            self.assertNotIn("22", df["CHR"].values)

    def test_manifest_row_group_metadata_is_consistent(self):
        result = make_multi_chrom_result()
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "out"
            LDScoreDirectoryWriter().write(result, LDScoreOutputConfig(output_dir=output_dir))

            manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["row_group_layout"], "one_per_chromosome")

            for field_name, total_field in [
                ("baseline_row_groups", "n_baseline_rows"),
                ("query_row_groups", "n_query_rows"),
            ]:
                entries = manifest[field_name]
                self.assertIsNotNone(entries)
                expected_offset = 0
                for e in entries:
                    self.assertEqual(e["row_offset"], expected_offset)
                    expected_offset += e["n_rows"]
                self.assertEqual(expected_offset, manifest[total_field])

    def test_query_row_groups_null_when_no_query_table(self):
        result = make_split_ldscore_result(query=False)
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "out"
            LDScoreDirectoryWriter().write(result, LDScoreOutputConfig(output_dir=output_dir))

            manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
            self.assertIsNone(manifest["query_row_groups"])
            self.assertIsNotNone(manifest["baseline_row_groups"])
            self.assertEqual(len(manifest["baseline_row_groups"]), 1)

    def test_omits_query_parquet_for_baseline_only_result(self):
        result = make_split_ldscore_result(query=False)
        writer = LDScoreDirectoryWriter()
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "ldscores"
            output_paths = writer.write(result, LDScoreOutputConfig(output_dir=output_dir))

            self.assertEqual(set(output_paths), {"manifest", "baseline"})
            self.assertFalse((output_dir / "query.parquet").exists())
            manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["files"], {"baseline": "baseline.parquet"})
            self.assertEqual(manifest["query_columns"], [])

    def test_rejects_query_table_with_mismatched_row_keys(self):
        result = make_split_ldscore_result(query=True)
        bad_query = result.query_table.copy()
        bad_query.loc[1, "POS"] = 999
        bad_result = dataclass_replace(result, query_table=bad_query)

        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(ValueError, "query rows must match baseline rows"):
                LDScoreDirectoryWriter().write(bad_result, LDScoreOutputConfig(output_dir=tmpdir))

    def test_refuses_overwrite_before_writing_any_file(self):
        result = make_split_ldscore_result(query=True)
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            (output_dir / "baseline.parquet").write_text("existing", encoding="utf-8")

            with self.assertRaisesRegex(FileExistsError, "overwrite=True"):
                LDScoreDirectoryWriter().write(result, LDScoreOutputConfig(output_dir=output_dir))

            self.assertFalse((output_dir / "manifest.json").exists())
            self.assertEqual((output_dir / "baseline.parquet").read_text(encoding="utf-8"), "existing")

    def test_overwrite_true_replaces_existing_canonical_files(self):
        result = make_split_ldscore_result(query=True)
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            (output_dir / "manifest.json").write_text("existing", encoding="utf-8")
            (output_dir / "baseline.parquet").write_text("existing", encoding="utf-8")

            LDScoreDirectoryWriter().write(result, LDScoreOutputConfig(output_dir=output_dir, overwrite=True))

            manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
            baseline = pd.read_parquet(output_dir / "baseline.parquet")
            self.assertEqual(manifest["format"], "ldsc.ldscore_result.v1")
            self.assertEqual(baseline["SNP"].tolist(), ["rs1", "rs2"])

    def test_omits_missing_common_count_values_from_manifest_records(self):
        result = make_split_ldscore_result(query=False)
        result = dataclass_replace(
            result,
            count_records=[
                {
                    "group": "baseline",
                    "column": "base",
                    "all_reference_snp_count": 10.0,
                }
            ],
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "ldscores"
            LDScoreDirectoryWriter().write(result, LDScoreOutputConfig(output_dir=output_dir))

            manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))

        self.assertEqual(manifest["counts"][0]["all_reference_snp_count"], 10.0)
        self.assertNotIn("common_reference_snp_count", manifest["counts"][0])
        self.assertEqual(manifest["count_config"]["common_reference_snp_maf_operator"], ">=")


class PartitionedH2DirectoryWriterTest(unittest.TestCase):
    def make_summary(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "query_annotation": ["IL-6/JAK STAT (Hallmark)", "IL 6 JAK STAT Hallmark"],
                "coefficient": [1.0, 2.0],
                "coefficient_se": [0.1, 0.2],
                "category_h2": [0.3, 0.4],
            }
        )

    def make_category_tables(self) -> dict[str, pd.DataFrame]:
        return {
            "IL-6/JAK STAT (Hallmark)": pd.DataFrame(
                {
                    "query_annotation": ["base", "IL-6/JAK STAT (Hallmark)"],
                    "coefficient": [0.5, 1.0],
                    "category_h2": [0.2, 0.3],
                }
            ),
            "IL 6 JAK STAT Hallmark": pd.DataFrame(
                {
                    "query_annotation": ["base", "IL 6 JAK STAT Hallmark"],
                    "coefficient": [0.6, 2.0],
                    "category_h2": [0.25, 0.4],
                }
            ),
        }

    def test_writes_aggregate_only_by_default(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "partitioned"
            paths = PartitionedH2DirectoryWriter().write(
                self.make_summary(),
                PartitionedH2OutputConfig(output_dir=output_dir),
            )

            self.assertTrue((output_dir / "partitioned_h2.tsv").exists())
            self.assertFalse((output_dir / "query_annotations").exists())
            self.assertEqual(paths, {"summary": str(output_dir / "partitioned_h2.tsv")})

    def test_writes_per_query_tree_with_sanitized_manifest(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "partitioned"
            paths = PartitionedH2DirectoryWriter().write(
                self.make_summary(),
                PartitionedH2OutputConfig(output_dir=output_dir, write_per_query_results=True),
                per_query_category_tables=self.make_category_tables(),
                metadata={"trait_name": "trait", "count_kind": "common", "ldscore_dir": "ldscores"},
            )

            query_root = output_dir / "query_annotations"
            manifest = pd.read_csv(query_root / "manifest.tsv", sep="\t")
            self.assertEqual(manifest["ordinal"].tolist(), [1, 2])
            self.assertEqual(
                manifest["query_annotation"].tolist(),
                ["IL-6/JAK STAT (Hallmark)", "IL 6 JAK STAT Hallmark"],
            )
            self.assertEqual(
                manifest["folder"].tolist(),
                ["0001_il-6_jak_stat_hallmark", "0002_il_6_jak_stat_hallmark"],
            )
            self.assertTrue((query_root / "0001_il-6_jak_stat_hallmark" / "partitioned_h2.tsv").exists())
            self.assertTrue((query_root / "0001_il-6_jak_stat_hallmark" / "model_categories.tsv").exists())
            metadata = json.loads(
                (query_root / "0001_il-6_jak_stat_hallmark" / "metadata.json").read_text(encoding="utf-8")
            )
            categories = pd.read_csv(query_root / "0001_il-6_jak_stat_hallmark" / "model_categories.tsv", sep="\t")
            self.assertEqual(metadata["query_annotation"], "IL-6/JAK STAT (Hallmark)")
            self.assertEqual(metadata["trait_name"], "trait")
            self.assertEqual(categories["query_annotation"].tolist(), ["base", "IL-6/JAK STAT (Hallmark)"])
            self.assertEqual(paths["per_query_root"], str(query_root))

    def test_refuses_existing_outputs_before_writing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "partitioned"
            output_dir.mkdir()
            existing = output_dir / "partitioned_h2.tsv"
            existing.write_text("existing\n", encoding="utf-8")

            with self.assertRaisesRegex(FileExistsError, "overwrite"):
                PartitionedH2DirectoryWriter().write(
                    self.make_summary(),
                    PartitionedH2OutputConfig(output_dir=output_dir, write_per_query_results=True),
                    per_query_category_tables=self.make_category_tables(),
                )

            self.assertEqual(existing.read_text(encoding="utf-8"), "existing\n")
            self.assertFalse((output_dir / "query_annotations").exists())

    def test_overwrite_replaces_existing_per_query_tree(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "partitioned"
            stale = output_dir / "query_annotations" / "stale"
            stale.mkdir(parents=True)
            (stale / "old.txt").write_text("old\n", encoding="utf-8")
            (output_dir / "partitioned_h2.tsv").write_text("old\n", encoding="utf-8")

            PartitionedH2DirectoryWriter().write(
                self.make_summary(),
                PartitionedH2OutputConfig(
                    output_dir=output_dir,
                    overwrite=True,
                    write_per_query_results=True,
                ),
                per_query_category_tables=self.make_category_tables(),
            )

            self.assertFalse(stale.exists())
            self.assertTrue((output_dir / "query_annotations" / "manifest.tsv").exists())
            self.assertIn("coefficient", (output_dir / "partitioned_h2.tsv").read_text(encoding="utf-8"))


class FixedOutputDirectoryTest(unittest.TestCase):
    def test_missing_munge_output_dir_warns_and_is_created(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "nested" / "trait"
            munger = __import__("ldsc.sumstats_munger", fromlist=["SumstatsMunger"]).SumstatsMunger()
            table = __import__("ldsc.sumstats_munger", fromlist=["SumstatsTable"]).SumstatsTable(
                data=pd.DataFrame({"SNP": ["rs1"], "N": [1000.0], "Z": [1.5]}),
                has_alleles=False,
                source_path="source.sumstats.gz",
                trait_name="trait",
            )

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                written = munger.write_output(table, str(output_dir))

            self.assertTrue(output_dir.exists())
            self.assertEqual(len(caught), 1)
            self.assertIn("output directory", str(caught[0].message).lower())
            self.assertIn("created", str(caught[0].message).lower())
            self.assertTrue(Path(written).exists())

    def test_munge_write_output_refuses_existing_file_by_default(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            existing = output_dir / "sumstats.parquet"
            existing.write_text("existing", encoding="utf-8")
            munger = __import__("ldsc.sumstats_munger", fromlist=["SumstatsMunger"]).SumstatsMunger()
            table = __import__("ldsc.sumstats_munger", fromlist=["SumstatsTable"]).SumstatsTable(
                data=pd.DataFrame({"SNP": ["rs1"], "N": [1000.0], "Z": [1.5]}),
                has_alleles=False,
                source_path="source.sumstats.gz",
                trait_name="trait",
            )

            with self.assertRaisesRegex(FileExistsError, "overwrite"):
                munger.write_output(table, output_dir)

            self.assertEqual(existing.read_text(encoding="utf-8"), "existing")

    def test_munge_write_output_allows_existing_file_with_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            existing = output_dir / "sumstats.parquet"
            existing.write_text("existing", encoding="utf-8")
            munger = __import__("ldsc.sumstats_munger", fromlist=["SumstatsMunger"]).SumstatsMunger()
            table = __import__("ldsc.sumstats_munger", fromlist=["SumstatsTable"]).SumstatsTable(
                data=pd.DataFrame({"SNP": ["rs1"], "N": [1000.0], "Z": [1.5]}),
                has_alleles=False,
                source_path="source.sumstats.gz",
                trait_name="trait",
            )

            written = munger.write_output(table, output_dir, overwrite=True)

            self.assertEqual(written, str(existing))
            self.assertGreater(existing.stat().st_size, len("existing"))
