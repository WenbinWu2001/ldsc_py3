from pathlib import Path
import importlib
import sys
import tempfile
import unittest
from unittest import mock

import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc.config import GlobalConfig, RefPanelConfig
from ldsc._kernel.ref_panel import ParquetR2RefPanel, PlinkRefPanel, RefPanelLoader


FIXTURES = Path(__file__).resolve().parent / "fixtures" / "minimal_external_resources" / "plink"
PLINK_PREFIX = FIXTURES / "hm3_chr22_subset"
PLINK_SNPS = ["22:16406147:A:G", "22:16805059:T:C"]


def _has_module(name: str) -> bool:
    try:
        importlib.import_module(name)
        return True
    except Exception:
        return False


class RefPanelLoaderTest(unittest.TestCase):
    def test_loader_selects_plink_backend(self):
        loader = RefPanelLoader(GlobalConfig(snp_identifier="rsid"), RefPanelConfig())
        spec = RefPanelConfig(backend="plink", plink_path=str(PLINK_PREFIX))
        panel = loader.load(spec)
        self.assertIsInstance(panel, PlinkRefPanel)

    def test_ref_panel_config_accepts_r2_bias_mode(self):
        spec = RefPanelConfig(backend="parquet_r2", r2_bias_mode="unbiased")
        self.assertEqual(spec.r2_bias_mode, "unbiased")

    def test_ref_panel_config_accepts_ref_panel_snps_path(self):
        spec = RefPanelConfig(backend="plink", ref_panel_snps_path="/path/to/snps.txt")
        self.assertEqual(spec.ref_panel_snps_path, "/path/to/snps.txt")

    def test_loader_selects_parquet_backend(self):
        loader = RefPanelLoader(GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"), RefPanelConfig())
        spec = RefPanelConfig(backend="parquet_r2", chromosomes=("1",), metadata_paths=())
        panel = loader.load(spec)
        self.assertIsInstance(panel, ParquetR2RefPanel)


class PlinkRefPanelTest(unittest.TestCase):
    def test_available_chromosomes_and_metadata(self):
        panel = PlinkRefPanel(
            GlobalConfig(snp_identifier="rsid"),
            RefPanelConfig(backend="plink", plink_path=str(PLINK_PREFIX)),
        )
        self.assertEqual(panel.available_chromosomes(), ["22"])
        metadata = panel.load_metadata("22")
        self.assertEqual(list(metadata.columns), ["CHR", "SNP", "CM", "POS"])
        self.assertGreater(len(metadata), 0)

    def test_filter_to_snps_rsid(self):
        panel = PlinkRefPanel(
            GlobalConfig(snp_identifier="rsid"),
            RefPanelConfig(backend="plink", plink_path=str(PLINK_PREFIX)),
        )
        filtered = panel.filter_to_snps("22", set(PLINK_SNPS))
        self.assertEqual(set(filtered["SNP"]), set(PLINK_SNPS))

    def test_duplicate_rsid_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "dup"
            prefix.with_suffix(".bed").write_bytes(b"")
            prefix.with_suffix(".fam").write_text("fam iid 0 0 0 -9\n", encoding="utf-8")
            prefix.with_suffix(".bim").write_text(
                "1 rs_dup 0.0 10 A G\n"
                "1 rs_dup 0.0 20 A G\n",
                encoding="utf-8",
            )
            panel = PlinkRefPanel(
                GlobalConfig(snp_identifier="rsid"),
                RefPanelConfig(backend="plink", plink_path=str(prefix)),
            )
            with self.assertRaises(ValueError):
                panel.load_metadata("1")

    def test_global_restriction_applies(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            restrict = Path(tmpdir) / "restrict.txt"
            restrict.write_text(f"SNP\n{PLINK_SNPS[0]}\n", encoding="utf-8")
            panel = PlinkRefPanel(
                GlobalConfig(snp_identifier="rsid"),
                RefPanelConfig(backend="plink", plink_path=str(PLINK_PREFIX), ref_panel_snps_path=str(restrict)),
            )
            metadata = panel.load_metadata("22")
            self.assertEqual(set(metadata["SNP"]), {PLINK_SNPS[0]})

    @unittest.skipUnless(_has_module("bitarray"), "bitarray is not installed")
    def test_build_reader(self):
        panel = PlinkRefPanel(
            GlobalConfig(snp_identifier="rsid"),
            RefPanelConfig(backend="plink", plink_path=str(PLINK_PREFIX)),
        )
        reader = panel.build_reader("22")
        self.assertIsNotNone(reader)


class ParquetRefPanelTest(unittest.TestCase):
    def test_sidecar_metadata_loading(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            meta = Path(tmpdir) / "meta.tsv"
            meta.write_text(
                "CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n1\t20\trs2\t0.2\t0.8\n2\t30\trs3\t0.3\t0.4\n",
                encoding="utf-8",
            )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", metadata_paths=(str(meta),)),
            )
            self.assertEqual(panel.available_chromosomes(), ["1", "2"])
            metadata = panel.load_metadata("1")
            self.assertEqual(metadata["CHR"].tolist(), ["1", "1"])
            self.assertEqual(metadata["MAF"].round(3).tolist(), [0.2, 0.2])

    def test_sidecar_metadata_loading_normalizes_float_chromosomes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            meta = Path(tmpdir) / "meta.tsv"
            meta.write_text(
                "CHR\tBP\tSNP\tCM\tMAF\n1.0\t10\trs1\t0.1\t0.2\nchr2.0\t20\trs2\t0.2\t0.3\n",
                encoding="utf-8",
            )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", metadata_paths=(str(meta),)),
            )

            self.assertEqual(panel.available_chromosomes(), ["1", "2"])
            metadata = panel.load_metadata("2")
            self.assertEqual(metadata["CHR"].tolist(), ["2"])

    def test_filter_to_snps_chr_pos(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            meta = Path(tmpdir) / "meta.tsv"
            meta.write_text(
                "CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n1\t20\trs2\t0.2\t0.3\n",
                encoding="utf-8",
            )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", metadata_paths=(str(meta),), chromosomes=("1",)),
            )
            filtered = panel.filter_to_snps("1", {"1:20"})
            self.assertEqual(filtered["POS"].tolist(), [20])

    def test_sidecar_metadata_loading_accepts_suite_tokens(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            meta1 = tmpdir / "meta.1.tsv"
            meta2 = tmpdir / "meta.2.tsv"
            meta1.write_text(
                "CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n",
                encoding="utf-8",
            )
            meta2.write_text(
                "CHR\tBP\tSNP\tCM\tMAF\n2\t20\trs2\t0.2\t0.3\n",
                encoding="utf-8",
            )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", metadata_paths=str(tmpdir / "meta.@.tsv")),
            )

            self.assertEqual(panel.available_chromosomes(), ["1", "2"])
            metadata = panel.load_metadata("2")
            self.assertEqual(metadata["POS"].tolist(), [20])

    def test_sidecar_metadata_loading_auto_genome_build_shifts_zero_based_positions_and_logs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            meta = Path(tmpdir) / "meta.tsv"
            rows = ["CHR\tBP\tSNP\tCM\tMAF"]
            positions = [1000 + (idx * 10) for idx in range(250)]
            for idx, pos in enumerate(positions, start=1):
                rows.append(f"1\t{pos - 1}\trs{idx}\t0.1\t0.2")
            meta.write_text("\n".join(rows) + "\n", encoding="utf-8")
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="auto"),
                RefPanelConfig(backend="parquet_r2", metadata_paths=(str(meta),)),
            )

            with mock.patch(
                "ldsc._kernel.ref_panel.load_packaged_reference_table",
                return_value=pd.DataFrame(
                    {
                        "CHR": ["1"] * 250,
                        "hg19_POS": positions,
                        "hg38_POS": [5000 + (idx * 10) for idx in range(250)],
                    }
                ),
            ), self.assertLogs("LDSC.ref_panel", level="INFO") as logs:
                metadata = panel.load_metadata("1")

            self.assertEqual(metadata["POS"].tolist()[:3], [1000, 1010, 1020])
            self.assertTrue(any("0-based" in message and "hg19" in message for message in logs.output))

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_build_reader(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            meta = Path(tmpdir) / "meta.tsv"
            meta.write_text("CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n", encoding="utf-8")
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", metadata_paths=(str(meta),), r2_paths=(str(Path(tmpdir) / "missing.parquet"),)),
            )
            with self.assertRaises(Exception):
                panel.build_reader("1")

    def test_plink_prefix_accepts_suite_token(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            prefix = tmpdir / "panel.1"
            (tmpdir / "panel.1.bed").write_text("", encoding="utf-8")
            (tmpdir / "panel.1.bim").write_text(
                "1 rs1 0.1 10 A G\n",
                encoding="utf-8",
            )
            (tmpdir / "panel.1.fam").write_text(
                "f i 0 0 0 -9\n",
                encoding="utf-8",
            )
            panel = PlinkRefPanel(
                GlobalConfig(snp_identifier="rsid"),
                RefPanelConfig(backend="plink", plink_path=str(tmpdir / "panel.@")),
            )

            self.assertEqual(panel.available_chromosomes(), ["1"])
            metadata = panel.load_metadata("1")
            self.assertEqual(metadata["SNP"].tolist(), ["rs1"])
