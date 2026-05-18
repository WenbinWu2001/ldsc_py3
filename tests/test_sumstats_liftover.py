from __future__ import annotations

import gzip
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc._kernel.liftover import (
    Hm3DualBuildLifter,
    LiftOverMappingResult,
    SumstatsLiftoverRequest,
    apply_sumstats_liftover,
    load_hm3_curated_map,
)


class SumstatsLiftoverTest(unittest.TestCase):
    SIDEcar_DTYPES = {
        "CHR": "string",
        "SNP": "string",
        "source_pos": "Int64",
        "target_pos": "Int64",
        "reason": "string",
    }

    def write_hm3_map(self, root: Path, text: str) -> Path:
        path = root / "hm3.tsv.gz"
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
        return path

    def test_load_hm3_curated_map_accepts_aliases_and_ignores_extra_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_hm3_map(
                Path(tmpdir),
                "chromosome\tGRCh37_bp\tGRCh38_position\trsid\textra\n"
                "chr1\t100\t1000\trs1\tignored\n"
                "1\t200\t2000\trs2\tignored\n",
            )

            loaded = load_hm3_curated_map(path)

            self.assertEqual(loaded.columns.tolist(), ["CHR", "hg19_POS", "hg38_POS", "SNP"])
            self.assertEqual(loaded["CHR"].tolist(), ["1", "1"])
            self.assertEqual(loaded["hg19_POS"].tolist(), [100, 200])
            self.assertEqual(loaded["hg38_POS"].tolist(), [1000, 2000])
            self.assertEqual(loaded["SNP"].tolist(), ["rs1", "rs2"])

    def test_load_hm3_curated_map_rejects_duplicate_source_or_target_coordinates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            duplicate_hg19 = self.write_hm3_map(
                root,
                "CHR\thg19_POS\thg38_POS\tSNP\n"
                "1\t100\t1000\trs1\n"
                "1\t100\t2000\trs2\n",
            )
            with self.assertRaisesRegex(ValueError, "duplicate.*hg19"):
                load_hm3_curated_map(duplicate_hg19)

            duplicate_hg38 = self.write_hm3_map(
                root,
                "CHR\thg19_POS\thg38_POS\tSNP\n"
                "1\t100\t1000\trs1\n"
                "1\t200\t1000\trs2\n",
            )
            with self.assertRaisesRegex(ValueError, "duplicate.*hg38"):
                load_hm3_curated_map(duplicate_hg38)

    def test_hm3_liftover_uses_coordinates_only_and_preserves_snp_labels(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_hm3_map(
                Path(tmpdir),
                "CHR\thg19_POS\thg38_POS\tSNP\n"
                "1\t100\t1000\tmap_rs1\n"
                "1\t200\t2000\tmap_rs2\n",
            )
            frame = pd.DataFrame(
                {
                    "CHR": ["1", "1", "1"],
                    "POS": [100, 200, 999],
                    "SNP": ["raw_a", "raw_b", "raw_missing"],
                    "Z": [1.0, 2.0, 3.0],
                    "N": [100.0, 100.0, 100.0],
                }
            )
            request = SumstatsLiftoverRequest(target_build="hg38", use_hm3_quick_liftover=True, hm3_map_file=path)

            lifted, report, drop_frame = apply_sumstats_liftover(
                frame,
                request,
                source_build="hg19",
                snp_identifier="chr_pos",
            )

            self.assertEqual(lifted["SNP"].tolist(), ["raw_a", "raw_b"])
            self.assertEqual(lifted["POS"].tolist(), [1000, 2000])
            self.assertEqual(report["method"], "hm3_curated")
            self.assertEqual(report["n_input"], 3)
            self.assertEqual(report["n_lifted"], 2)
            self.assertEqual(report["n_unmapped"], 1)
            self.assertEqual(drop_frame["reason"].tolist(), ["unmapped_liftover"])

    def test_hm3_dual_build_lifter_supports_reverse_and_identity(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_hm3_map(
                Path(tmpdir),
                "CHR\thg19_POS\thg38_POS\tSNP\n"
                "1\t100\t1000\trs1\n"
                "1\t200\t2000\trs2\n",
            )
            frame = pd.DataFrame({"CHR": ["1", "1"], "POS": [1000, 9999], "SNP": ["raw_a", "raw_missing"]})

            lifted, unmapped = Hm3DualBuildLifter("hg38", "hg19", path).lift(frame)

            self.assertEqual(lifted["SNP"].tolist(), ["raw_a"])
            self.assertEqual(lifted["POS"].tolist(), [100])
            self.assertEqual(unmapped.tolist(), [1])

            identity, identity_unmapped = Hm3DualBuildLifter("hg38", "hg38", path).lift(frame)
            pd.testing.assert_frame_equal(identity, frame)
            self.assertEqual(identity_unmapped.tolist(), [])

    def test_chain_liftover_counts_cross_chromosome_drops(self):
        frame = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "POS": [100, 200, 300],
                "SNP": ["rs1", "rs2", "rs3"],
                "Z": [1.0, 2.0, 3.0],
                "N": [100.0, 100.0, 100.0],
            }
        )
        request = SumstatsLiftoverRequest(target_build="hg38", liftover_chain_file="chain.over")

        fake_translator = mock.Mock()
        fake_translator.map_positions.return_value = LiftOverMappingResult(
            translated_positions=np.asarray([1000, 3000], dtype=np.int64),
            keep_mask=np.asarray([True, False, True], dtype=bool),
            unmapped_count=0,
            cross_chrom_count=1,
        )
        with mock.patch("ldsc._kernel.liftover.LiftOverTranslator", return_value=fake_translator):
            lifted, report, drop_frame = apply_sumstats_liftover(
                frame,
                request,
                source_build="hg19",
                snp_identifier="chr_pos",
            )

        self.assertEqual(lifted["SNP"].tolist(), ["rs1", "rs3"])
        self.assertEqual(lifted["POS"].tolist(), [1000, 3000])
        self.assertEqual(report["method"], "chain_file")
        self.assertEqual(report["n_cross_chrom"], 1)
        self.assertEqual(report["n_unmapped"], 0)
        self.assertEqual(report["n_dropped"], 1)
        self.assertEqual(drop_frame["reason"].tolist(), ["cross_chromosome_liftover"])

    def test_source_duplicate_coordinate_groups_are_removed_before_mapping(self):
        frame = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "POS": [100, 100, 300],
                "SNP": ["rs1", "rs2", "rs3"],
                "Z": [1.0, 2.0, 3.0],
                "N": [100.0, 100.0, 100.0],
            }
        )
        request = SumstatsLiftoverRequest(target_build="hg38", liftover_chain_file="chain.over")

        fake_translator = mock.Mock()

        def fake_map(chrom, positions):
            self.assertEqual(chrom, "1")
            np.testing.assert_array_equal(positions, np.asarray([300], dtype=np.int64))
            return LiftOverMappingResult(
                translated_positions=np.asarray([3000], dtype=np.int64),
                keep_mask=np.asarray([True], dtype=bool),
                unmapped_count=0,
                cross_chrom_count=0,
            )

        fake_translator.map_positions.side_effect = fake_map
        with mock.patch("ldsc._kernel.liftover.LiftOverTranslator", return_value=fake_translator):
            with self.assertLogs("LDSC.liftover", level="INFO") as caught:
                lifted, report, drop_frame = apply_sumstats_liftover(
                    frame,
                    request,
                    source_build="hg19",
                    snp_identifier="chr_pos",
                )

        self.assertEqual(lifted["SNP"].tolist(), ["rs3"])
        self.assertEqual(lifted["POS"].tolist(), [3000])
        self.assertEqual(report["n_duplicate_source_dropped"], 2)
        self.assertEqual(report["n_duplicate_target_dropped"], 0)
        self.assertIn("source_duplicate", "\n".join(caught.output))
        self.assertNotIn("rs1", "\n".join(caught.output))
        self.assertEqual(drop_frame["reason"].tolist(), ["source_duplicate", "source_duplicate"])

    def test_duplicate_target_coordinate_groups_are_removed_entirely(self):
        frame = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "POS": [100, 200, 300],
                "SNP": ["rs1", "rs2", "rs3"],
                "Z": [1.0, 2.0, 3.0],
                "N": [100.0, 100.0, 100.0],
            }
        )
        request = SumstatsLiftoverRequest(target_build="hg38", liftover_chain_file="chain.over")
        fake_translator = mock.Mock()
        fake_translator.map_positions.return_value = LiftOverMappingResult(
            translated_positions=np.asarray([1000, 1000, 3000], dtype=np.int64),
            keep_mask=np.asarray([True, True, True], dtype=bool),
            unmapped_count=0,
            cross_chrom_count=0,
        )

        with mock.patch("ldsc._kernel.liftover.LiftOverTranslator", return_value=fake_translator):
            lifted, report, drop_frame = apply_sumstats_liftover(
                frame,
                request,
                source_build="hg19",
                snp_identifier="chr_pos",
            )

        self.assertEqual(lifted["SNP"].tolist(), ["rs3"])
        self.assertEqual(lifted["POS"].tolist(), [3000])
        self.assertEqual(report["n_duplicate_source_dropped"], 0)
        self.assertEqual(report["n_duplicate_target_dropped"], 2)
        self.assertEqual(report["n_lifted"], 1)
        self.assertEqual(drop_frame["reason"].tolist(), ["target_collision", "target_collision"])

    def test_liftover_duplicate_coordinate_detection_is_allele_blind(self):
        frame = pd.DataFrame(
            {
                "CHR": ["1", "1", "1", "1", "1"],
                "POS": [100, 100, 200, 300, 400],
                "SNP": ["source_a", "source_b", "target_a", "target_b", "kept"],
                "A1": ["A", "A", "A", "A", "A"],
                "A2": ["C", "G", "C", "G", "C"],
                "Z": [1.0, 2.0, 3.0, 4.0, 5.0],
                "N": [100.0, 100.0, 100.0, 100.0, 100.0],
            }
        )
        request = SumstatsLiftoverRequest(target_build="hg38", liftover_chain_file="chain.over")
        fake_translator = mock.Mock()
        fake_translator.map_positions.return_value = LiftOverMappingResult(
            translated_positions=np.asarray([900, 900, 1200], dtype=np.int64),
            keep_mask=np.asarray([True, True, True], dtype=bool),
            unmapped_count=0,
            cross_chrom_count=0,
        )

        with mock.patch("ldsc._kernel.liftover.LiftOverTranslator", return_value=fake_translator):
            lifted, report, drop_frame = apply_sumstats_liftover(
                frame,
                request,
                source_build="hg19",
                snp_identifier="chr_pos_allele_aware",
            )

        self.assertEqual(lifted["SNP"].tolist(), ["kept"])
        self.assertEqual(report["n_duplicate_source_dropped"], 2)
        self.assertEqual(report["n_duplicate_target_dropped"], 2)
        self.assertEqual(
            drop_frame["reason"].tolist(),
            ["source_duplicate", "source_duplicate", "target_collision", "target_collision"],
        )

    def test_missing_coordinates_are_dropped_before_hm3_liftover(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_hm3_map(
                Path(tmpdir),
                "CHR\thg19_POS\thg38_POS\tSNP\n"
                "1\t100\t1000\trs1\n",
            )
            request = SumstatsLiftoverRequest(target_build="hg38", use_hm3_quick_liftover=True, hm3_map_file=path)
            frame = pd.DataFrame(
                {
                    "CHR": ["1", "1", pd.NA, "1"],
                    "POS": [100, pd.NA, 200, 999],
                    "SNP": ["rs1", "missing_pos", "missing_chr", "unmapped"],
                    "Z": [1.0, 2.0, 3.0, 4.0],
                    "N": [100.0, 100.0, 100.0, 100.0],
                }
            )

            with self.assertLogs("LDSC.liftover", level="INFO") as caught:
                lifted, report, drop_frame = apply_sumstats_liftover(
                    frame,
                    request,
                    source_build="hg19",
                    snp_identifier="chr_pos",
                )

            self.assertEqual(lifted["SNP"].tolist(), ["rs1"])
            self.assertEqual(lifted["POS"].tolist(), [1000])
            self.assertEqual(report["n_input"], 4)
            self.assertEqual(report["n_lifted"], 1)
            self.assertEqual(report["n_missing_chr_pos_dropped"], 2)
            self.assertEqual(report["n_unmapped"], 1)
            self.assertEqual(report["n_dropped"], 3)
            log_text = "\n".join(caught.output)
            self.assertIn("Dropped 2 SNPs with invalid or missing CHR/POS", log_text)
            self.assertIn("missing CHR=1", log_text)
            self.assertIn("missing POS=1", log_text)
            self.assertNotIn("missing_pos", log_text)
            self.assertEqual(
                drop_frame["reason"].tolist(),
                ["missing_coordinate", "missing_coordinate", "unmapped_liftover"],
            )

    def test_chain_liftover_drop_frame_contains_all_sumstats_reasons(self):
        frame = pd.DataFrame(
            {
                "CHR": ["1", pd.NA, "1", "1", "1", "1", "1", "1", "1"],
                "POS": [pd.NA, 150, 100, 100, 200, 300, 400, 500, 600],
                "SNP": [
                    "missing_pos",
                    "missing_chr",
                    "source_dup_a",
                    "source_dup_b",
                    "unmapped",
                    "cross_chrom",
                    "target_a",
                    "target_b",
                    "kept",
                ],
                "Z": [1.0] * 9,
                "N": [100.0] * 9,
            }
        )
        request = SumstatsLiftoverRequest(target_build="hg38", liftover_chain_file="chain.over")
        fake_translator = mock.Mock()
        fake_translator.map_positions.return_value = LiftOverMappingResult(
            translated_positions=np.asarray([900, 900, 1600], dtype=np.int64),
            keep_mask=np.asarray([False, False, True, True, True], dtype=bool),
            unmapped_count=1,
            cross_chrom_count=1,
            unmapped_mask=np.asarray([True, False, False, False, False], dtype=bool),
            cross_chrom_mask=np.asarray([False, True, False, False, False], dtype=bool),
        )

        with mock.patch("ldsc._kernel.liftover.LiftOverTranslator", return_value=fake_translator):
            lifted, report, drop_frame = apply_sumstats_liftover(
                frame,
                request,
                source_build="hg19",
                snp_identifier="chr_pos",
            )

        self.assertEqual(lifted["SNP"].tolist(), ["kept"])
        self.assertEqual(report["n_dropped"], 8)
        self.assertEqual(
            set(drop_frame["reason"]),
            {
                "missing_coordinate",
                "source_duplicate",
                "unmapped_liftover",
                "cross_chromosome_liftover",
                "target_collision",
            },
        )
        self.assertEqual(drop_frame.columns.tolist(), ["CHR", "SNP", "source_pos", "target_pos", "reason"])
        self.assertEqual({column: str(dtype) for column, dtype in drop_frame.dtypes.items()}, self.SIDEcar_DTYPES)

    def test_liftover_errors_after_missing_coordinate_drop_removes_everything(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self.write_hm3_map(
                Path(tmpdir),
                "CHR\thg19_POS\thg38_POS\tSNP\n"
                "1\t100\t1000\trs1\n",
            )
            request = SumstatsLiftoverRequest(target_build="hg38", use_hm3_quick_liftover=True, hm3_map_file=path)
            missing = pd.DataFrame({"CHR": ["1"], "POS": [pd.NA], "SNP": ["rs1"], "Z": [1.0], "N": [100.0]})
            with self.assertRaisesRegex(ValueError, "dropped all rows"):
                apply_sumstats_liftover(missing, request, source_build="hg19", snp_identifier="chr_pos")

            unmapped = pd.DataFrame({"CHR": ["1"], "POS": [999], "SNP": ["rs2"], "Z": [1.0], "N": [100.0]})
            with self.assertRaisesRegex(ValueError, "dropped all rows"):
                apply_sumstats_liftover(unmapped, request, source_build="hg19", snp_identifier="chr_pos")

            duplicate_source = pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "POS": [100, 100],
                    "SNP": ["rs1", "rs2"],
                    "Z": [1.0, 2.0],
                    "N": [100.0, 100.0],
                }
            )
            with self.assertRaisesRegex(ValueError, "dropped all rows"):
                apply_sumstats_liftover(duplicate_source, request, source_build="hg19", snp_identifier="chr_pos")

    def test_noop_liftover_metadata_schema(self):
        frame = pd.DataFrame({"CHR": ["1"], "POS": [100], "SNP": ["rs1"], "Z": [1.0], "N": [100.0]})
        request = SumstatsLiftoverRequest()

        lifted, report, drop_frame = apply_sumstats_liftover(
            frame,
            request,
            source_build="hg38",
            snp_identifier="chr_pos",
        )

        pd.testing.assert_frame_equal(lifted, frame)
        self.assertEqual(
            report,
            {
                "applied": False,
                "source_build": None,
                "target_build": None,
                "method": None,
                "chain_file": None,
                "hm3_map_file": None,
                "n_input": None,
                "n_lifted": None,
                "n_dropped": None,
                "n_missing_chr_pos_dropped": None,
                "n_unmapped": None,
                "n_cross_chrom": None,
                "n_duplicate_source_dropped": None,
                "n_duplicate_target_dropped": None,
            },
        )
        self.assertEqual(drop_frame.columns.tolist(), ["CHR", "SNP", "source_pos", "target_pos", "reason"])
        self.assertEqual(len(drop_frame), 0)
        self.assertEqual({column: str(dtype) for column, dtype in drop_frame.dtypes.items()}, self.SIDEcar_DTYPES)

    def test_noop_liftover_accepts_coordinate_family_mode(self):
        frame = pd.DataFrame({"CHR": ["1"], "POS": [100], "SNP": ["rs1"], "Z": [1.0], "N": [100.0]})
        request = SumstatsLiftoverRequest(target_build="hg38")

        lifted, report, drop_frame = apply_sumstats_liftover(
            frame,
            request,
            source_build="hg38",
            snp_identifier="chr_pos_allele_aware",
        )

        pd.testing.assert_frame_equal(lifted, frame)
        self.assertFalse(report["applied"])
        self.assertEqual(len(drop_frame), 0)

    def test_liftover_rejects_rsid_family_mode(self):
        frame = pd.DataFrame({"CHR": ["1"], "POS": [100], "SNP": ["rs1"], "Z": [1.0], "N": [100.0]})
        request = SumstatsLiftoverRequest(target_build="hg38")

        with self.assertRaisesRegex(ValueError, "chr_pos-family"):
            apply_sumstats_liftover(
                frame,
                request,
                source_build="hg38",
                snp_identifier="rsid_allele_aware",
            )


if __name__ == "__main__":
    unittest.main()
