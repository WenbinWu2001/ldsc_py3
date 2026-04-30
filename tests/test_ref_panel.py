from pathlib import Path
import gzip
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


def _write_canonical_r2_parquet(path: Path, chrom: str = "1") -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq

    path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(
        {
            "CHR": [chrom, chrom],
            "POS_1": [10, 20],
            "POS_2": [20, 30],
            "R2": pd.Series([0.4, 0.2], dtype="float32"),
            "SNP_1": ["rs1", "rs2"],
            "SNP_2": ["rs2", "rs3"],
        }
    )
    table = pa.Table.from_pandas(frame, preserve_index=False)
    table = table.replace_schema_metadata({b"ldsc:sorted_by_build": b"hg38"})
    pq.write_table(table, path)


class RefPanelLoaderTest(unittest.TestCase):
    def test_loader_selects_plink_backend(self):
        loader = RefPanelLoader(GlobalConfig(snp_identifier="rsid"), RefPanelConfig())
        spec = RefPanelConfig(backend="plink", plink_prefix=str(PLINK_PREFIX))
        panel = loader.load(spec)
        self.assertIsInstance(panel, PlinkRefPanel)

    def test_ref_panel_config_accepts_r2_bias_mode(self):
        spec = RefPanelConfig(backend="parquet_r2", r2_bias_mode="unbiased")
        self.assertEqual(spec.r2_bias_mode, "unbiased")

    def test_ref_panel_config_accepts_ref_panel_snps_file(self):
        spec = RefPanelConfig(backend="plink", ref_panel_snps_file="/path/to/snps.txt")
        self.assertEqual(spec.ref_panel_snps_file, "/path/to/snps.txt")

    def test_ref_panel_config_accepts_ref_panel_dir(self):
        spec = RefPanelConfig(backend="parquet_r2", ref_panel_dir="/path/to/panel/hg38")
        self.assertEqual(spec.ref_panel_dir, "/path/to/panel/hg38")

    def test_loader_selects_parquet_backend(self):
        loader = RefPanelLoader(GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"), RefPanelConfig())
        spec = RefPanelConfig(backend="parquet_r2", chromosomes=("1",), metadata_sources=())
        panel = loader.load(spec)
        self.assertIsInstance(panel, ParquetR2RefPanel)


class PlinkRefPanelTest(unittest.TestCase):
    def test_available_chromosomes_and_metadata(self):
        panel = PlinkRefPanel(
            GlobalConfig(snp_identifier="rsid"),
            RefPanelConfig(backend="plink", plink_prefix=str(PLINK_PREFIX)),
        )
        self.assertEqual(panel.available_chromosomes(), ["22"])
        metadata = panel.load_metadata("22")
        self.assertEqual(metadata.loc[:, ["CHR", "SNP", "CM", "POS"]].columns.tolist(), ["CHR", "SNP", "CM", "POS"])
        self.assertGreater(len(metadata), 0)

    def test_filter_to_snps_rsid(self):
        panel = PlinkRefPanel(
            GlobalConfig(snp_identifier="rsid"),
            RefPanelConfig(backend="plink", plink_prefix=str(PLINK_PREFIX)),
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
                RefPanelConfig(backend="plink", plink_prefix=str(prefix)),
            )
            with self.assertRaises(ValueError):
                panel.load_metadata("1")

    def test_global_restriction_applies(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            restrict = Path(tmpdir) / "restrict.txt"
            restrict.write_text(f"SNP\n{PLINK_SNPS[0]}\n", encoding="utf-8")
            panel = PlinkRefPanel(
                GlobalConfig(snp_identifier="rsid"),
                RefPanelConfig(backend="plink", plink_prefix=str(PLINK_PREFIX), ref_panel_snps_file=str(restrict)),
            )
            metadata = panel.load_metadata("22")
            self.assertEqual(set(metadata["SNP"]), {PLINK_SNPS[0]})

    @unittest.skipUnless(_has_module("bitarray"), "bitarray is not installed")
    def test_plink_metadata_includes_genotype_maf_and_applies_maf_min(self):
        panel = PlinkRefPanel(
            GlobalConfig(snp_identifier="rsid"),
            RefPanelConfig(backend="plink", plink_prefix=str(PLINK_PREFIX), maf_min=0.2),
        )

        metadata = panel.load_metadata("22")

        self.assertIn("MAF", metadata.columns)
        self.assertTrue((metadata["MAF"] > 0.2).all())

    @unittest.skipUnless(_has_module("bitarray"), "bitarray is not installed")
    def test_build_reader(self):
        panel = PlinkRefPanel(
            GlobalConfig(snp_identifier="rsid"),
            RefPanelConfig(backend="plink", plink_prefix=str(PLINK_PREFIX)),
        )
        reader = panel.build_reader("22")
        self.assertIsNotNone(reader)


class ParquetRefPanelTest(unittest.TestCase):
    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_ref_panel_dir_loads_sidecar_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet")
            with gzip.open(build_dir / "chr1_meta.tsv.gz", "wt", encoding="utf-8") as handle:
                handle.write("CHR\tPOS\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n1\t20\trs2\t0.2\t0.3\n")

            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", ref_panel_dir=str(build_dir)),
            )

            self.assertEqual(panel.available_chromosomes(), ["1"])
            metadata = panel.load_metadata("1")
            self.assertEqual(metadata["SNP"].tolist(), ["rs1", "rs2"])
            self.assertEqual(metadata["CM"].round(3).tolist(), [0.1, 0.2])

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_ref_panel_dir_falls_back_to_r2_endpoints_without_sidecar(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet")

            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", ref_panel_dir=str(build_dir)),
            )

            self.assertEqual(panel.available_chromosomes(), ["1"])
            metadata = panel.load_metadata("1")
            self.assertEqual(metadata["CHR"].tolist(), ["1", "1", "1"])
            self.assertEqual(metadata["POS"].tolist(), [10, 20, 30])
            self.assertEqual(metadata["SNP"].tolist(), ["rs1", "rs2", "rs3"])
            self.assertTrue(metadata["CM"].isna().all())
            self.assertNotIn("MAF", metadata.columns)

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_ref_panel_dir_requires_r2_for_requested_chromosome(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            build_dir.mkdir(parents=True)
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", ref_panel_dir=str(build_dir)),
            )

            with self.assertRaisesRegex(FileNotFoundError, "chr1_r2.parquet"):
                panel.load_metadata("1")

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_ref_panel_dir_missing_metadata_hard_fails_when_configured(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet")
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38", fail_on_missing_metadata=True),
                RefPanelConfig(backend="parquet_r2", ref_panel_dir=str(build_dir)),
            )

            with self.assertRaisesRegex(ValueError, "metadata sidecar"):
                panel.load_metadata("1")

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_ref_panel_dir_top_level_root_descends_by_genome_build(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "panel"
            _write_canonical_r2_parquet(root / "hg38" / "chr1_r2.parquet")

            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", ref_panel_dir=str(root)),
            )

            self.assertEqual(panel.available_chromosomes(), ["1"])

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_ref_panel_dir_top_level_root_rejects_ambiguous_build(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "panel"
            _write_canonical_r2_parquet(root / "hg19" / "chr1_r2.parquet")
            _write_canonical_r2_parquet(root / "hg38" / "chr1_r2.parquet")

            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="rsid"),
                RefPanelConfig(backend="parquet_r2", ref_panel_dir=str(root)),
            )

            with self.assertRaisesRegex(ValueError, "ambiguous"):
                panel.available_chromosomes()

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_ref_panel_dir_top_level_root_rejects_build_dirs_even_when_only_one_has_r2(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "panel"
            (root / "hg19").mkdir(parents=True)
            _write_canonical_r2_parquet(root / "hg38" / "chr1_r2.parquet")

            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="rsid"),
                RefPanelConfig(backend="parquet_r2", ref_panel_dir=str(root)),
            )

            with self.assertRaisesRegex(ValueError, "ambiguous"):
                panel.available_chromosomes()

    def test_sidecar_metadata_loading(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            meta = Path(tmpdir) / "meta.tsv"
            meta.write_text(
                "CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n1\t20\trs2\t0.2\t0.8\n2\t30\trs3\t0.3\t0.4\n",
                encoding="utf-8",
            )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", metadata_sources=(str(meta),)),
            )
            self.assertEqual(panel.available_chromosomes(), ["1", "2"])
            metadata = panel.load_metadata("1")
            self.assertEqual(metadata["CHR"].tolist(), ["1", "1"])
            self.assertEqual(metadata["MAF"].round(3).tolist(), [0.2, 0.2])

    def test_sidecar_metadata_applies_maf_min_filter(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            meta = Path(tmpdir) / "meta.tsv"
            meta.write_text(
                "CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n1\t20\trs2\t0.2\t0.3\n1\t30\trs3\t0.3\t0.4\n",
                encoding="utf-8",
            )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", metadata_sources=(str(meta),), maf_min=0.25),
            )

            metadata = panel.load_metadata("1")

            self.assertEqual(metadata["SNP"].tolist(), ["rs2", "rs3"])
            self.assertTrue((metadata["MAF"] > 0.25).all())

    def test_sidecar_metadata_loading_normalizes_float_chromosomes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            meta = Path(tmpdir) / "meta.tsv"
            meta.write_text(
                "CHR\tBP\tSNP\tCM\tMAF\n1.0\t10\trs1\t0.1\t0.2\nchr2.0\t20\trs2\t0.2\t0.3\n",
                encoding="utf-8",
            )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", metadata_sources=(str(meta),)),
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
                RefPanelConfig(backend="parquet_r2", metadata_sources=(str(meta),), chromosomes=("1",)),
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
                RefPanelConfig(backend="parquet_r2", metadata_sources=str(tmpdir / "meta.@.tsv")),
            )

            self.assertEqual(panel.available_chromosomes(), ["1", "2"])
            metadata = panel.load_metadata("2")
            self.assertEqual(metadata["POS"].tolist(), [20])

    def test_sidecar_metadata_loading_rejects_auto_genome_build_at_kernel_boundary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            meta = Path(tmpdir) / "meta.tsv"
            rows = ["CHR\tBP\tSNP\tCM\tMAF"]
            positions = [1000 + (idx * 10) for idx in range(250)]
            for idx, pos in enumerate(positions, start=1):
                rows.append(f"1\t{pos - 1}\trs{idx}\t0.1\t0.2")
            meta.write_text("\n".join(rows) + "\n", encoding="utf-8")
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="auto"),
                RefPanelConfig(backend="parquet_r2", metadata_sources=(str(meta),)),
            )

            with self.assertRaisesRegex(AssertionError, "workflow entry"):
                panel.load_metadata("1")

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_build_reader(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            meta = Path(tmpdir) / "meta.tsv"
            meta.write_text("CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n", encoding="utf-8")
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", metadata_sources=(str(meta),), r2_sources=(str(Path(tmpdir) / "missing.parquet"),)),
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
                RefPanelConfig(backend="plink", plink_prefix=str(tmpdir / "panel.@")),
            )

            self.assertEqual(panel.available_chromosomes(), ["1"])
            metadata = panel.load_metadata("1")
            self.assertEqual(metadata["SNP"].tolist(), ["rs1"])
