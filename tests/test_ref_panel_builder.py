import gzip
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock
from typing import Optional

import numpy as np
import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc._kernel import ref_panel_builder as kernel_builder
from ldsc.config import CommonConfig, ReferencePanelBuildConfig
from ldsc import ldscore_calculator, ref_panel_builder


_HAS_BITARRAY = importlib.util.find_spec("bitarray") is not None
_HAS_PYARROW = importlib.util.find_spec("pyarrow") is not None
_HAS_PYLIFTOVER = importlib.util.find_spec("pyliftover") is not None


def _find_resources_root() -> Optional[Path]:
    for parent in Path(__file__).resolve().parents:
        candidate = parent / "resources"
        if candidate.exists():
            return candidate
    return None


class GeneticMapParserTest(unittest.TestCase):
    def test_load_genetic_map_accepts_alkes_headers_and_normalizes_chr23_to_x(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "map.txt"
            path.write_text(
                "chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n"
                "23 10 0.0 0.1\n"
                "23 20 0.0 0.2\n",
                encoding="utf-8",
            )

            genetic_map = kernel_builder.load_genetic_map(path)

            self.assertEqual(genetic_map.columns.tolist(), ["CHR", "POS", "CM"])
            self.assertEqual(genetic_map["CHR"].tolist(), ["X", "X"])
            self.assertEqual(genetic_map["POS"].tolist(), [10, 20])
            self.assertEqual(genetic_map["CM"].tolist(), [0.1, 0.2])

    def test_load_genetic_map_accepts_basic_alias_family(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "map.tsv"
            path.write_text(
                "chromosome bp genetic_map_cm\n"
                "1 100 0.0\n"
                "1 200 1.0\n",
                encoding="utf-8",
            )

            genetic_map = kernel_builder.load_genetic_map(path)

            self.assertEqual(genetic_map["CHR"].tolist(), ["1", "1"])
            self.assertEqual(genetic_map["POS"].tolist(), [100, 200])
            self.assertEqual(genetic_map["CM"].tolist(), [0.0, 1.0])

    def test_load_genetic_map_rejects_unsorted_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "map.txt"
            path.write_text(
                "chr position Genetic_Map(cM)\n"
                "1 200 1.0\n"
                "1 100 0.0\n",
                encoding="utf-8",
            )

            with self.assertRaises(ValueError):
                kernel_builder.load_genetic_map(path)

    def test_load_genetic_map_rejects_duplicate_positions(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "map.txt"
            path.write_text(
                "chr position Genetic_Map(cM)\n"
                "1 100 0.0\n"
                "1 100 0.1\n",
                encoding="utf-8",
            )

            with self.assertRaises(ValueError):
                kernel_builder.load_genetic_map(path)


class GeneticMapInterpolationTest(unittest.TestCase):
    def test_interpolate_cm_linearly_and_clamps_endpoints(self):
        genetic_map = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "POS": [100, 200, 400],
                "CM": [0.0, 1.0, 3.0],
            }
        )

        interpolated = kernel_builder.interpolate_genetic_map_cm(
            chrom="1",
            positions=np.array([50, 150, 300, 500], dtype=np.int64),
            genetic_map=genetic_map,
        )

        np.testing.assert_allclose(interpolated, np.array([0.0, 0.5, 2.0, 3.0]))


class LiftOverTranslatorTest(unittest.TestCase):
    def test_map_positions_filters_unmapped_and_cross_chrom_hits(self):
        class _FakeLiftOver:
            def convert_coordinate(self, chrom, pos0):
                table = {
                    ("chr22", 99): [("chr22", 999, "+", 0)],
                    ("chr22", 199): [("chr21", 1999, "+", 0)],
                    ("chr22", 299): [],
                }
                return table[(chrom, pos0)]

        translator = kernel_builder.LiftOverTranslator.__new__(kernel_builder.LiftOverTranslator)
        translator.source_build = "hg38"
        translator.target_build = "hg19"
        translator.chain_path = None
        translator._identity = False
        translator._liftover = _FakeLiftOver()

        result = translator.map_positions("22", np.array([100, 200, 300], dtype=np.int64))

        np.testing.assert_array_equal(result.translated_positions, np.array([1000], dtype=np.int64))
        np.testing.assert_array_equal(result.keep_mask, np.array([True, False, False]))
        self.assertEqual(result.cross_chrom_count, 1)
        self.assertEqual(result.unmapped_count, 1)


class RestrictionModeDetectionTest(unittest.TestCase):
    def test_detect_restriction_identifier_mode_accepts_rsid_lists(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.txt"
            path.write_text("rs1\nrs2\n", encoding="utf-8")

            mode = kernel_builder.detect_restriction_identifier_mode(path)

            self.assertEqual(mode, "rsid")

    def test_detect_restriction_identifier_mode_accepts_chr_pos_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.txt"
            path.write_text("1:100\n1:200\n", encoding="utf-8")

            mode = kernel_builder.detect_restriction_identifier_mode(path)

            self.assertEqual(mode, "chr_pos")


class _SequentialSNPGetter:
    def __init__(self, matrix: np.ndarray):
        self.matrix = matrix
        self.cursor = 0

    def __call__(self, width: int) -> np.ndarray:
        block = self.matrix[:, self.cursor:self.cursor + width]
        self.cursor += width
        return block


class PairwiseEmissionTest(unittest.TestCase):
    def test_iter_pairwise_r2_rows_emits_one_unordered_pair_per_window(self):
        standardized = np.array(
            [
                [-1.0, -1.0, 1.0, 1.0],
                [-1.0, 1.0, -1.0, 1.0],
                [1.0, -1.0, -1.0, 1.0],
                [1.0, 1.0, 1.0, -1.0],
            ]
        )
        block_left = np.array([0, 0, 1, 2], dtype=float)

        rows = kernel_builder.iter_pairwise_r2_rows(
            block_left=block_left,
            chunk_size=2,
            standardized_snp_getter=_SequentialSNPGetter(standardized),
            m=4,
            n=4,
        )

        expected_pairs = {(0, 1), (1, 2), (2, 3)}
        observed_pairs = {(row["i"], row["j"]) for row in rows}
        self.assertEqual(observed_pairs, expected_pairs)
        self.assertTrue(all(row["i"] < row["j"] for row in rows))

        brute_force = {}
        for j in range(1, standardized.shape[1]):
            for i in range(int(block_left[j]), j):
                corr = float(np.dot(standardized[:, i], standardized[:, j]) / standardized.shape[0])
                r2 = corr**2 - (1.0 - corr**2) / (standardized.shape[0] - 2)
                brute_force[(i, j)] = (r2, "+" if corr >= 0 else "-")

        for row in rows:
            expected_r2, expected_sign = brute_force[(row["i"], row["j"])]
            self.assertAlmostEqual(row["R2"], expected_r2)
            self.assertEqual(row["sign"], expected_sign)


class StandardTableFormattingTest(unittest.TestCase):
    def test_build_standard_annotation_table_uses_exact_schema(self):
        metadata = pd.DataFrame(
            {
                "CHR": ["1", "X"],
                "SNP": ["1:100:A:G", "rsX"],
                "POS": [100, 200],
                "MAF": [0.2, 0.3],
                "A1": ["A", "C"],
                "A2": ["G", "T"],
            }
        )

        table = kernel_builder.build_standard_annotation_table(
            metadata=metadata,
            hg19_positions=np.array([100, 150], dtype=np.int64),
            hg38_positions=np.array([110, 250], dtype=np.int64),
        )

        self.assertEqual(
            table.columns.tolist(),
            ["chr", "hg19_pos", "hg38_pos", "hg19_Uniq_ID", "hg38_Uniq_ID", "rsID", "MAF", "REF", "ALT"],
        )
        self.assertEqual(table["hg19_Uniq_ID"].tolist(), ["1:100:A:G", "X:150:C:T"])
        self.assertEqual(table["hg38_Uniq_ID"].tolist(), ["1:110:A:G", "X:250:C:T"])
        self.assertEqual(table["rsID"].tolist(), ["1:100:A:G", "rsX"])

    def test_build_standard_ld_table_uses_exact_schema(self):
        annotation_table = pd.DataFrame(
            {
                "chr": ["1", "1"],
                "hg19_pos": [100, 120],
                "hg38_pos": [110, 130],
                "hg19_Uniq_ID": ["1:100:A:G", "1:120:C:T"],
                "hg38_Uniq_ID": ["1:110:A:G", "1:130:C:T"],
                "rsID": ["rs1", "rs2"],
                "MAF": [0.2, 0.3],
                "REF": ["A", "C"],
                "ALT": ["G", "T"],
            }
        )
        pair_rows = [{"i": 0, "j": 1, "R2": 0.75, "sign": "-"}]

        table = kernel_builder.build_standard_ld_table(pair_rows=pair_rows, annotation_table=annotation_table)

        self.assertEqual(
            table.columns.tolist(),
            [
                "chr",
                "rsID_1",
                "rsID_2",
                "hg38_pos_1",
                "hg38_pos_2",
                "hg19_pos_1",
                "hg19_pos_2",
                "hg38_Uniq_ID_1",
                "hg38_Uniq_ID_2",
                "hg19_Uniq_ID_1",
                "hg19_Uniq_ID_2",
                "R2",
                "Dprime",
                "+/-corr",
            ],
        )
        self.assertEqual(table.loc[0, "rsID_1"], "rs1")
        self.assertEqual(table.loc[0, "rsID_2"], "rs2")
        self.assertTrue(pd.isna(table.loc[0, "Dprime"]))
        self.assertEqual(table.loc[0, "+/-corr"], "-")

    def test_build_runtime_metadata_table_is_build_specific(self):
        metadata = pd.DataFrame(
            {
                "CHR": ["1"],
                "SNP": ["rs1"],
                "MAF": [0.2],
            }
        )

        table = kernel_builder.build_runtime_metadata_table(
            metadata=metadata,
            positions=np.array([120], dtype=np.int64),
            cm_values=np.array([0.5], dtype=float),
        )

        self.assertEqual(table.columns.tolist(), ["CHR", "POS", "SNP", "CM", "MAF"])
        self.assertEqual(table.loc[0, "POS"], 120)
        self.assertEqual(table.loc[0, "CM"], 0.5)


class ReferencePanelBuilderWorkflowTest(unittest.TestCase):
    def _write_dummy_plink_prefix(self, root: Path, stem: str, chrom: str):
        prefix = root / stem
        Path(str(prefix) + ".bed").write_bytes(b"")
        Path(str(prefix) + ".fam").write_text("fam iid 0 0 0 -9\n", encoding="utf-8")
        Path(str(prefix) + ".bim").write_text(
            f"{chrom} rs{chrom} 0.0 100 A G\n",
            encoding="utf-8",
        )
        return prefix

    def _build_config(self, tmpdir: Path) -> ReferencePanelBuildConfig:
        map_hg19 = tmpdir / "hg19.map"
        map_hg38 = tmpdir / "hg38.map"
        map_hg19.write_text("chr position Genetic_Map(cM)\n1 100 0.0\n2 100 0.0\n", encoding="utf-8")
        map_hg38.write_text("chr position Genetic_Map(cM)\n1 100 0.0\n2 100 0.0\n", encoding="utf-8")
        return ReferencePanelBuildConfig(
            panel_label="EUR",
            plink_prefix=tmpdir / "panel.@",
            source_genome_build="hg19",
            genetic_map_hg19_path=map_hg19,
            genetic_map_hg38_path=map_hg38,
            output_dir=tmpdir / "out",
            ld_wind_kb=1.0,
        )

    def test_builder_run_collects_artifact_paths_from_resolved_suite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_dummy_plink_prefix(tmpdir, "panel.1", "1")
            self._write_dummy_plink_prefix(tmpdir, "panel.2", "2")
            config = self._build_config(tmpdir)
            builder = ref_panel_builder.ReferencePanelBuilder(common_config=CommonConfig())

            def fake_build(prefix, chrom, config, build_state):
                out_root = Path(config.output_dir) / "parquet"
                return {
                    "ann": str(out_root / "ann" / f"{config.panel_label}_chr{chrom}_ann.parquet"),
                    "ld": str(out_root / "ld" / f"{config.panel_label}_chr{chrom}_LD.parquet"),
                    "meta_hg19": str(out_root / "meta" / f"{config.panel_label}_chr{chrom}_meta_hg19.tsv.gz"),
                    "meta_hg38": str(out_root / "meta" / f"{config.panel_label}_chr{chrom}_meta_hg38.tsv.gz"),
                }

            with mock.patch.object(
                ref_panel_builder.ReferencePanelBuilder,
                "_build_chromosome",
                side_effect=fake_build,
            ) as patched:
                result = builder.run(config)

            self.assertEqual(result.panel_label, "EUR")
            self.assertEqual(result.chromosomes, ["1", "2"])
            self.assertEqual(patched.call_count, 2)
            self.assertEqual(
                result.output_paths["ann"],
                [
                    str(tmpdir / "out" / "parquet" / "ann" / "EUR_chr1_ann.parquet"),
                    str(tmpdir / "out" / "parquet" / "ann" / "EUR_chr2_ann.parquet"),
                ],
            )
            self.assertEqual(
                result.output_paths["meta_hg38"],
                [
                    str(tmpdir / "out" / "parquet" / "meta" / "EUR_chr1_meta_hg38.tsv.gz"),
                    str(tmpdir / "out" / "parquet" / "meta" / "EUR_chr2_meta_hg38.tsv.gz"),
                ],
            )

    def test_builder_run_rejects_duplicate_chromosome_across_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_dummy_plink_prefix(tmpdir, "panel_a", "1")
            self._write_dummy_plink_prefix(tmpdir, "panel_b", "1")
            map_hg19 = tmpdir / "hg19.map"
            map_hg38 = tmpdir / "hg38.map"
            map_hg19.write_text("chr position Genetic_Map(cM)\n1 100 0.0\n", encoding="utf-8")
            map_hg38.write_text("chr position Genetic_Map(cM)\n1 100 0.0\n", encoding="utf-8")
            config = ReferencePanelBuildConfig(
                panel_label="EUR",
                plink_prefix=str(tmpdir / "panel_*"),
                source_genome_build="hg19",
                genetic_map_hg19_path=map_hg19,
                genetic_map_hg38_path=map_hg38,
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(common_config=CommonConfig())

            with mock.patch.object(
                ref_panel_builder.ReferencePanelBuilder,
                "_build_chromosome",
                return_value={"ann": "ann", "ld": "ld", "meta_hg19": "m19", "meta_hg38": "m38"},
            ):
                with self.assertRaises(ValueError):
                    builder.run(config)


@unittest.skipUnless(
    _HAS_BITARRAY and _HAS_PYARROW and _HAS_PYLIFTOVER,
    "bitarray, pyarrow, and pyliftover are required for the smoke/parity builder test",
)
class ReferencePanelBuilderParityTest(unittest.TestCase):
    def test_example_1kg30x_chr22_matches_direct_plink_ldscore(self):
        resources = _find_resources_root()
        if resources is None:
            self.skipTest("resources directory is not available from this workspace")

        prefix = resources / "example_1kg_30x" / "genomes_30x_chr22"
        if not (Path(str(prefix) + ".bed").exists() and Path(str(prefix) + ".bim").exists() and Path(str(prefix) + ".fam").exists()):
            self.skipTest("example 1KG 30x chr22 PLINK inputs are unavailable")

        map_hg19 = resources / "genetic_maps" / "genetic_map_alkesgroup" / "genetic_map_hg19_withX.txt"
        map_hg38 = resources / "genetic_maps" / "genetic_map_alkesgroup" / "genetic_map_hg38_withX.txt"
        if not (map_hg19.exists() and map_hg38.exists()):
            self.skipTest("required bundled genetic maps are unavailable")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            build_result = ref_panel_builder.run_build_ref_panel(
                bfile=str(prefix),
                panel_label="SMOKE",
                source_genome_build="hg38",
                genetic_map_hg19=str(map_hg19),
                genetic_map_hg38=str(map_hg38),
                out=str(tmpdir / "panel"),
                ld_wind_kb=10.0,
                chunk_size=64,
            )

            self.assertEqual(build_result.chromosomes, ["22"])
            ld_path = build_result.output_paths["ld"][0]
            meta_hg38_path = build_result.output_paths["meta_hg38"][0]
            self.assertTrue(Path(ld_path).exists())
            self.assertTrue(Path(meta_hg38_path).exists())

            with gzip.open(meta_hg38_path, "rt", encoding="utf-8") as handle:
                meta = pd.read_csv(handle, sep="\t")
            baseline = tmpdir / "baseline.annot.gz"
            baseline_df = meta.loc[:, ["CHR", "POS", "SNP", "CM"]].rename(columns={"POS": "BP"})
            baseline_df["base"] = 1.0
            with gzip.open(baseline, "wt", encoding="utf-8") as handle:
                baseline_df.to_csv(handle, sep="\t", index=False)

            direct = ldscore_calculator.run_ldscore(
                out=str(tmpdir / "direct"),
                baseline_annot=str(baseline),
                bfile=str(prefix),
                snp_identifier="rsid",
                ld_wind_kb=10.0,
                chunk_size=64,
            )
            parquet = ldscore_calculator.run_ldscore(
                out=str(tmpdir / "parquet"),
                baseline_annot=str(baseline),
                r2_table=str(ld_path),
                frqfile=str(meta_hg38_path),
                snp_identifier="rsid",
                genome_build="hg38",
                r2_bias_mode="unbiased",
                ld_wind_kb=10.0,
                chunk_size=64,
            )

            self.assertEqual(
                direct.reference_metadata["SNP"].tolist(),
                parquet.reference_metadata["SNP"].tolist(),
            )
            np.testing.assert_allclose(
                direct.ld_scores["base"].to_numpy(dtype=float),
                parquet.ld_scores["base"].to_numpy(dtype=float),
                rtol=1e-5,
                atol=1e-5,
            )
            np.testing.assert_allclose(
                direct.w_ld["L2"].to_numpy(dtype=float),
                parquet.w_ld["L2"].to_numpy(dtype=float),
                rtol=1e-5,
                atol=1e-5,
            )
