from pathlib import Path
import gzip
import importlib
import tempfile
import unittest
from unittest import mock

import pandas as pd
import pytest


from ldsc.config import GlobalConfig, RefPanelConfig
from ldsc.errors import LDSCInputError
from ldsc._refpanel.adapters import (
    ParquetR2RefPanel,
    PlinkRefPanel,
    RefPanelLoader,
    _read_metadata_table,
    _snp_id_series_for_matching,
)
from ldsc._refpanel.build_kernel import build_restriction_mask

import numpy as np


def dict_chunks(rows):
    """Wrap legacy {i,j,R2,sign} dict rows as the chunk iterable write_r2_parquet expects."""
    rows = list(rows)
    if not rows:
        return iter([])
    i = np.array([r["i"] for r in rows], dtype=np.int64)
    j = np.array([r["j"] for r in rows], dtype=np.int64)
    r2 = np.array([r["R2"] for r in rows], dtype=np.float32)
    sign = np.array([1 if r["sign"] == "+" else -1 for r in rows], dtype=np.int8)
    return [(i, j, r2, sign)]


FIXTURES = Path(__file__).resolve().parent / "fixtures" / "minimal_external_resources" / "plink"
PLINK_PREFIX = FIXTURES / "hm3_chr22_subset"
PLINK_SNPS = ["22:16406147:A:G", "22:16805059:T:C"]


def _has_module(name: str) -> bool:
    try:
        importlib.import_module(name)
        return True
    except Exception:
        return False


def _write_canonical_r2_parquet(
    path: Path,
    chrom: str = "1",
    *,
    identity_metadata: bool = True,
    snp_identifier: str = "chr_pos",
    genome_build: str = "hg38",
) -> None:
    """Write a package-format R2 parquet with no pairs; sidecar tests pair it with their own rows."""
    import pyarrow.parquet as pq
    from ldsc._refpanel.r2_format import write_r2_parquet

    write_r2_parquet(
        pair_chunks=[], path=path, genome_build=genome_build, n_samples=10,
        snp_identifier=snp_identifier, n_snps=0, sidecar_identity_sha256="0" * 64,
        ld_window_mode="kb", ld_window_value=1000.0,
    )
    if not identity_metadata:
        table = pq.read_table(path)
        stripped = {key: value for key, value in table.schema.metadata.items()
                    if key not in {b"ldsc:artifact_type", b"ldsc:snp_identifier", b"ldsc:genome_build"}}
        pq.write_table(table.replace_schema_metadata(stripped), path)


def _write_meta_sidecar(
    path: Path,
    rows: str,
    *,
    identity_metadata: bool = True,
    snp_identifier: str = "chr_pos",
    genome_build: str = "hg38",
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        if identity_metadata:
            handle.write("# ldsc:artifact_type=ref_panel_metadata\n")
            handle.write(f"# ldsc:snp_identifier={snp_identifier}\n")
            handle.write(f"# ldsc:genome_build={genome_build}\n")
        handle.write(rows)


class RefPanelLoaderTest(unittest.TestCase):
    def test_loader_selects_plink_backend(self):
        loader = RefPanelLoader(GlobalConfig(snp_identifier="rsid"), RefPanelConfig())
        spec = RefPanelConfig(backend="plink", plink_prefix=str(PLINK_PREFIX))
        panel = loader.load(spec)
        self.assertIsInstance(panel, PlinkRefPanel)

    def test_ref_panel_config_accepts_ref_panel_snps_file(self):
        spec = RefPanelConfig(backend="plink", ref_panel_snps_file="/path/to/snps.txt")
        self.assertEqual(spec.ref_panel_snps_file, "/path/to/snps.txt")

    def test_ref_panel_config_accepts_r2_dir(self):
        spec = RefPanelConfig(backend="parquet_r2", r2_dir="/path/to/panel/hg38")
        self.assertEqual(spec.r2_dir, "/path/to/panel/hg38")

    def test_ref_panel_config_rejects_removed_parquet_source_fields(self):
        for name in ("r2_ref_panel_dir", "ref_panel_dir", "r2_sources", "metadata_sources"):
            with self.subTest(name=name):
                with self.assertRaises(TypeError):
                    RefPanelConfig(backend="parquet_r2", **{name: "/path/to/panel/hg38"})

    def test_loader_selects_parquet_backend(self):
        loader = RefPanelLoader(GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"), RefPanelConfig())
        spec = RefPanelConfig(backend="parquet_r2", chromosomes=("1",))
        panel = loader.load(spec)
        self.assertIsInstance(panel, ParquetR2RefPanel)

    def test_auto_build_resolution_uses_coordinate_family(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            restriction = Path(tmpdir) / "snps.tsv"
            restriction.write_text("CHR\tHG38_POS\n1\t10\n", encoding="utf-8")
            loader = RefPanelLoader(
                GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="auto"),
                RefPanelConfig(),
            )

            resolved = loader._global_config_for_spec(
                RefPanelConfig(backend="plink", ref_panel_snps_file=str(restriction))
            )

            self.assertEqual(resolved.snp_identifier, "chr_pos_allele_aware")
            self.assertEqual(resolved.genome_build, "hg38")


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

    def test_matching_ids_use_rsid_family(self):
        metadata = pd.DataFrame({"SNP": ["rs1", "rs2"], "A1": ["A", "A"], "A2": ["C", "G"]})

        ids = _snp_id_series_for_matching(metadata, "rsid_allele_aware", context="test")

        self.assertEqual(ids.tolist(), ["rs1:A:C", "rs2:A:G"])

    def test_build_restriction_mask_uses_rsid_family(self):
        metadata = pd.DataFrame({"SNP": ["rs1", "rs2"]})

        mask = build_restriction_mask(metadata, {"rs2"}, "rsid_allele_aware")

        self.assertEqual(mask.tolist(), [False, True])

    def test_ref_panel_restriction_uses_allele_aware_identity_keys(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            restrict_path = Path(tmpdir) / "restrict.tsv"
            restrict_path.write_text("SNP\tA1\tA2\nrs1\tA\tC\n", encoding="utf-8")
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="rsid_allele_aware"),
                RefPanelConfig(backend="parquet_r2", r2_dir=tmpdir, ref_panel_snps_file=str(restrict_path)),
            )
            metadata = pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "POS": [10, 10],
                    "SNP": ["rs1", "rs1"],
                    "A1": ["A", "A"],
                    "A2": ["C", "G"],
                }
            )

            filtered = panel._apply_snp_restriction(metadata)

            self.assertEqual(filtered["A2"].tolist(), ["C"])

    def test_plink_ref_panel_restriction_uses_bim_alleles(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            restrict_path = Path(tmpdir) / "restrict.tsv"
            restrict_path.write_text(f"SNP\tA1\tA2\n{PLINK_SNPS[0]}\tG\tA\n", encoding="utf-8")
            panel = PlinkRefPanel(
                GlobalConfig(snp_identifier="rsid_allele_aware"),
                RefPanelConfig(
                    backend="plink",
                    plink_prefix=str(PLINK_PREFIX),
                    ref_panel_snps_file=str(restrict_path),
                ),
            )

            metadata = panel.load_metadata("22")

            self.assertEqual(metadata["SNP"].tolist(), [PLINK_SNPS[0]])
            self.assertEqual(metadata["A1"].tolist(), ["G"])
            self.assertEqual(metadata["A2"].tolist(), ["A"])

    def test_metadata_reader_requires_snp_for_rsid_family(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "meta.tsv.gz"
            _write_meta_sidecar(path, "CHR\tPOS\n1\t10\n", snp_identifier="rsid_allele_aware", genome_build="hg19")

            with self.assertRaisesRegex(LDSCInputError, "no SNP column was found"):
                _read_metadata_table(
                    path,
                    None,
                    GlobalConfig(snp_identifier="rsid_allele_aware"),
                )

    def test_duplicate_rsid_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "dup"
            prefix.with_suffix(".bed").write_bytes(b"\x6c\x1b\x01\x02\x02")
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
            with self.assertRaisesRegex(LDSCInputError, "after duplicate SNP identity cleanup"):
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

    def test_small_panel_needs_override_even_with_declared_keep_list(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            restrict = Path(tmpdir) / "restrict.tsv"
            restrict.write_text(
                "CHR\thg38_POS\n"
                "22\t16406147\n"
                "22\t16805059\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(LDSCInputError, 'unresolved:'):
                RefPanelLoader(GlobalConfig(snp_identifier='chr_pos')).load(
                    RefPanelConfig(backend='plink', plink_prefix=str(PLINK_PREFIX), ref_panel_snps_file=str(restrict)))
            panel = RefPanelLoader(GlobalConfig(snp_identifier="chr_pos", genome_build="hg38")).load(
                RefPanelConfig(
                    backend="plink",
                    plink_prefix=str(PLINK_PREFIX),
                    ref_panel_snps_file=str(restrict),
                )
            )

            metadata = panel.load_metadata("22")

        self.assertEqual(set(metadata["SNP"]), set(PLINK_SNPS))
        self.assertEqual(panel.global_config.genome_build, "hg38")

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
    def test_missing_sidecar_is_a_hard_error(self):
        from ldsc._refpanel import build_kernel as kb
        from ldsc._refpanel import r2_format

        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            build_dir.mkdir(parents=True)
            # index parquet present, but no chr1_meta.tsv.gz sidecar alongside it
            r2_format.write_r2_parquet(
                pair_chunks=dict_chunks([{"i": 0, "j": 1, "R2": 0.5, "sign": "+"}]),
                path=build_dir / "chr1_r2.parquet", genome_build="hg38", n_samples=10,
                snp_identifier="chr_pos", n_snps=2, sidecar_identity_sha256="0" * 64,
                ld_window_mode="kb", ld_window_value=1000.0,
            )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir)),
            )
            with self.assertRaisesRegex(LDSCInputError, "sidecar"):
                panel.load_metadata("1")

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_r2_dir_loads_sidecar_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet")
            _write_meta_sidecar(build_dir / "chr1_meta.tsv.gz", "CHR\tPOS\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n1\t20\trs2\t0.2\t0.3\n")

            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir)),
            )

            self.assertEqual(panel.available_chromosomes(), ["1"])
            metadata = panel.load_metadata("1")
            self.assertEqual(metadata["SNP"].tolist(), ["rs1", "rs2"])
            self.assertEqual(metadata["CM"].round(3).tolist(), [0.1, 0.2])

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_r2_dir_loads_allele_metadata_sidecar_in_allele_aware_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(
                build_dir / "chr1_r2.parquet",
                snp_identifier="chr_pos_allele_aware",
            )
            _write_meta_sidecar(
                build_dir / "chr1_meta.tsv.gz",
                "CHR\tPOS\tSNP\tCM\tMAF\tA1\tA2\n1\t10\trs1\t0.1\t0.2\tA\tC\n1\t20\trs2\t0.2\t0.3\tA\tG\n",
                snp_identifier="chr_pos_allele_aware",
            )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir)),
            )

            metadata = panel.load_metadata("1")

            self.assertEqual(metadata["SNP"].tolist(), ["rs1", "rs2"])
            self.assertEqual(metadata["A1"].tolist(), ["A", "A"])
            self.assertEqual(metadata["A2"].tolist(), ["C", "G"])

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_metadata_sidecar_missing_alleles_is_malformed_in_allele_aware_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(
                build_dir / "chr1_r2.parquet",
                snp_identifier="chr_pos_allele_aware",
            )
            _write_meta_sidecar(
                build_dir / "chr1_meta.tsv.gz",
                "CHR\tPOS\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n",
                snp_identifier="chr_pos_allele_aware",
            )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir)),
            )

            with self.assertRaisesRegex(LDSCInputError, "malformed.*A1/A2"):
                panel.load_metadata("1")

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_package_canonical_r2_rejects_missing_identity_schema_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet", identity_metadata=False)
            _write_meta_sidecar(build_dir / "chr1_meta.tsv.gz", "CHR\tPOS\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n")

            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir)),
            )

            with self.assertRaisesRegex(LDSCInputError, "ldsc:artifact_type is missing"):
                panel.build_reader("1")

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_package_metadata_sidecar_rejects_missing_identity_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet")
            _write_meta_sidecar(
                build_dir / "chr1_meta.tsv.gz",
                "CHR\tPOS\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n",
                identity_metadata=False,
            )

            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir)),
            )

            with self.assertRaisesRegex(
                LDSCInputError,
                "no LDSC3 identity metadata",
            ):
                panel.load_metadata("1")

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_package_metadata_sidecar_rejects_identity_mismatch(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet")
            _write_meta_sidecar(
                build_dir / "chr1_meta.tsv.gz",
                "CHR\tPOS\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n",
                snp_identifier="chr_pos_allele_aware",
            )

            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir)),
            )

            with self.assertRaisesRegex(ValueError, "snp_identifier mismatch"):
                panel.load_metadata("1")

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_package_canonical_r2_rejects_snp_identifier_mismatch(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet", snp_identifier="chr_pos")
            _write_meta_sidecar(build_dir / "chr1_meta.tsv.gz", "CHR\tPOS\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n")

            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="rsid"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir)),
            )

            with self.assertRaisesRegex(ValueError, "snp_identifier mismatch"):
                panel.build_reader("1")

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_concrete_build_rejects_conflicting_r2_declaration(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet", genome_build="hg38")
            _write_meta_sidecar(build_dir / "chr1_meta.tsv.gz", "CHR\tPOS\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n")

            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir)),
            )

            with self.assertRaisesRegex(LDSCInputError, "Conflicting genome-build evidence"):
                panel.build_reader("1")

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_r2_dir_requires_r2_for_requested_chromosome(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            build_dir.mkdir(parents=True)
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir)),
            )

            with self.assertRaisesRegex(LDSCInputError, "chr1_r2.parquet"):
                panel.load_metadata("1")

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_r2_dir_missing_metadata_hard_fails_when_configured(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet")
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir)),
            )

            with self.assertRaisesRegex(LDSCInputError, "metadata sidecar"):
                panel.load_metadata("1")

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_r2_dir_top_level_root_descends_by_genome_build(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "panel"
            _write_canonical_r2_parquet(root / "hg38" / "chr1_r2.parquet")

            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(root)),
            )

            self.assertEqual(panel.available_chromosomes(), ["1"])

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_r2_dir_top_level_root_rejects_ambiguous_build(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "panel"
            _write_canonical_r2_parquet(root / "hg19" / "chr1_r2.parquet")
            _write_canonical_r2_parquet(root / "hg38" / "chr1_r2.parquet")

            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="rsid"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(root)),
            )

            with self.assertRaisesRegex(LDSCInputError, "ambiguous"):
                panel.available_chromosomes()

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_r2_dir_top_level_root_rejects_build_dirs_even_when_only_one_has_r2(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "panel"
            (root / "hg19").mkdir(parents=True)
            _write_canonical_r2_parquet(root / "hg38" / "chr1_r2.parquet")

            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="rsid"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(root)),
            )

            with self.assertRaisesRegex(LDSCInputError, "ambiguous"):
                panel.available_chromosomes()

    def test_sidecar_metadata_loading(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet")
            _write_canonical_r2_parquet(build_dir / "chr2_r2.parquet", chrom="2")
            _write_meta_sidecar(
                build_dir / "chr1_meta.tsv.gz",
                "CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n1\t20\trs2\t0.2\t0.2\n2\t30\trs3\t0.3\t0.4\n",
            )
            _write_meta_sidecar(build_dir / "chr2_meta.tsv.gz", "CHR\tBP\tSNP\tCM\tMAF\n2\t30\trs3\t0.3\t0.4\n")
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir)),
            )
            self.assertEqual(panel.available_chromosomes(), ["1", "2"])
            metadata = panel.load_metadata("1")
            self.assertEqual(metadata["CHR"].tolist(), ["1", "1"])
            self.assertEqual(metadata["MAF"].round(3).tolist(), [0.2, 0.2])

    def test_sidecar_metadata_applies_maf_min_filter(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet")
            _write_meta_sidecar(
                build_dir / "chr1_meta.tsv.gz",
                "CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n1\t20\trs2\t0.2\t0.3\n1\t30\trs3\t0.3\t0.4\n",
            )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir), maf_min=0.25),
            )

            metadata = panel.load_metadata("1")

            self.assertEqual(metadata["SNP"].tolist(), ["rs2", "rs3"])
            self.assertTrue((metadata["MAF"] > 0.25).all())

    def test_sidecar_metadata_maf_min_is_inclusive_at_threshold(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet")
            _write_meta_sidecar(
                build_dir / "chr1_meta.tsv.gz",
                "CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.24\n1\t20\trs2\t0.2\t0.25\n1\t30\trs3\t0.3\t0.4\n",
            )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir), maf_min=0.25),
            )

            metadata = panel.load_metadata("1")

            # MAF == 0.25 is retained (inclusive >=); 0.24 dropped.
            self.assertEqual(metadata["SNP"].tolist(), ["rs2", "rs3"])

    def test_sidecar_metadata_rejects_invalid_maf_values(self):
        # build-r2-panel records MAF = freq(A1) <= 0.5 from genotypes for every
        # sidecar row, so a missing or out-of-range value means the artifact was
        # altered; it must never be read as "below threshold", folded into range,
        # or silently skip a requested --maf-min filter.
        sidecars = {
            "partial": ("CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n1\t20\trs2\t0.2\tNA\n1\t30\trs3\t0.3\t0.4\n",
                        "1 of 3", "rs2"),
            "all": ("CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\tNA\n1\t20\trs2\t0.2\tNA\n",
                    "2 of 2", "rs1"),
            "out_of_range": ("CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0\n1\t20\trs2\t0.2\t0.7\n"
                             "1\t30\trs3\t0.3\t-0.2\n1\t40\trs4\t0.4\tinf\n1\t50\trs5\t0.5\t1.5\n"
                             "1\t60\trs6\t0.6\t0.5\n",
                             "4 of 6", "rs2, rs3, rs4, rs5"),
        }
        for label, (rows, count, example) in sidecars.items():
            for maf_min in (None, 0.25):
                with self.subTest(sidecar=label, maf_min=maf_min), tempfile.TemporaryDirectory() as tmpdir:
                    build_dir = Path(tmpdir) / "panel" / "hg38"
                    _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet")
                    _write_meta_sidecar(build_dir / "chr1_meta.tsv.gz", rows)
                    panel = ParquetR2RefPanel(
                        GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                        RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir), maf_min=maf_min),
                    )

                    with self.assertRaisesRegex(
                        LDSCInputError, f"missing, non-numeric, or out-of-range MAF for {count} SNP rows"
                    ) as ctx:
                        panel.load_metadata("1")
                    self.assertIn(example, str(ctx.exception))

    def test_sidecar_metadata_keeps_boundary_maf_values_unchanged(self):
        # 0 (monomorphic) and 0.5 are the closed ends of the canonical MAF range.
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet")
            _write_meta_sidecar(
                build_dir / "chr1_meta.tsv.gz",
                "CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0\n1\t20\trs2\t0.2\t0.5\n1\t30\trs3\t0.3\t0.25\n",
            )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir)),
            )

            metadata = panel.load_metadata("1")

            self.assertEqual(metadata["MAF"].tolist(), [0.0, 0.5, 0.25])

    def test_sidecar_metadata_reads_back_the_values_the_builder_used(self):
        # build-r2-panel windows on round_to_sidecar_precision(CM); the reader must
        # return exactly those values for CM and MAF.
        from ldsc._refpanel.r2_format import round_to_sidecar_precision, write_runtime_metadata_sidecar

        cm = np.array([12.3456789012345, 123.456789012345, 0.000123456789012])
        maf = np.array([0.0123456789012, 0.333333333333, 1.23456789012e-05])
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet")
            write_runtime_metadata_sidecar(
                pd.DataFrame({"CHR": ["1"] * 3, "POS": [10, 20, 30], "SNP": ["rs1", "rs2", "rs3"],
                              "A1": ["A"] * 3, "A2": ["G"] * 3, "CM": cm, "MAF": maf}),
                build_dir / "chr1_meta.tsv.gz", genome_build="hg38", snp_identifier="chr_pos",
            )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir)),
            )

            metadata = panel.load_metadata("1")

            np.testing.assert_array_equal(metadata["CM"].to_numpy(dtype=float), round_to_sidecar_precision(cm))
            np.testing.assert_array_equal(metadata["MAF"].to_numpy(dtype=float), round_to_sidecar_precision(maf))

    def test_sidecar_metadata_loading_normalizes_float_chromosomes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet")
            _write_canonical_r2_parquet(build_dir / "chr2_r2.parquet", chrom="2")
            _write_meta_sidecar(
                build_dir / "chr1_meta.tsv.gz",
                "CHR\tBP\tSNP\tCM\tMAF\n1.0\t10\trs1\t0.1\t0.2\n",
            )
            _write_meta_sidecar(
                build_dir / "chr2_meta.tsv.gz",
                "CHR\tBP\tSNP\tCM\tMAF\n1.0\t10\trs1\t0.1\t0.2\nchr2.0\t20\trs2\t0.2\t0.3\n",
            )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir)),
            )

            self.assertEqual(panel.available_chromosomes(), ["1", "2"])
            metadata = panel.load_metadata("2")
            self.assertEqual(metadata["CHR"].tolist(), ["2"])

    def test_filter_to_snps_chr_pos(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet")
            _write_meta_sidecar(
                build_dir / "chr1_meta.tsv.gz",
                "CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n1\t20\trs2\t0.2\t0.3\n",
            )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir), chromosomes=("1",)),
            )
            filtered = panel.filter_to_snps("1", {"1:20"})
            self.assertEqual(filtered["POS"].tolist(), [20])

    def test_r2_dir_sidecar_metadata_loading_by_chromosome(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet")
            _write_canonical_r2_parquet(build_dir / "chr2_r2.parquet", chrom="2")
            _write_meta_sidecar(build_dir / "chr1_meta.tsv.gz", "CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n")
            _write_meta_sidecar(build_dir / "chr2_meta.tsv.gz", "CHR\tBP\tSNP\tCM\tMAF\n2\t20\trs2\t0.2\t0.3\n")
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir)),
            )

            self.assertEqual(panel.available_chromosomes(), ["1", "2"])
            metadata = panel.load_metadata("2")
            self.assertEqual(metadata["POS"].tolist(), [20])

    def test_sidecar_metadata_loading_rejects_auto_genome_build_at_kernel_boundary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_canonical_r2_parquet(build_dir / "chr1_r2.parquet")
            rows = ["CHR\tBP\tSNP\tCM\tMAF"]
            positions = [1000 + (idx * 10) for idx in range(250)]
            for idx, pos in enumerate(positions, start=1):
                rows.append(f"1\t{pos - 1}\trs{idx}\t0.1\t0.2")
            _write_meta_sidecar(build_dir / "chr1_meta.tsv.gz", "\n".join(rows) + "\n")
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="auto"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir)),
            )

            with self.assertRaisesRegex(AssertionError, "workflow entry"):
                panel.load_metadata("1")

    @unittest.skipUnless(_has_module("pyarrow"), "pyarrow is not installed")
    def test_build_reader(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_meta_sidecar(build_dir / "chr1_meta.tsv.gz", "CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n")
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir)),
            )
            with self.assertRaises(LDSCInputError):
                panel.build_reader("1")

    def test_plink_prefix_accepts_suite_token(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            prefix = tmpdir / "panel.1"
            (tmpdir / "panel.1.bed").write_bytes(b"\x6c\x1b\x01\x02")
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


_CHR22 = Path(__file__).resolve().parent / "fixtures" / "minimal_external_resources" / "plink" / "hm3_chr22_subset"


def _chr22_available() -> bool:
    return all(Path(str(_CHR22) + ext).exists() for ext in (".bed", ".bim", ".fam"))


def test_ref_panel_has_no_region_exclusion_configuration(tmp_path):
    if not _chr22_available():
        pytest.skip("chr22 PLINK fixture unavailable; run tests/fixtures/generate_minimal_external_resources.py")
    gc = GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38")

    base_panel = PlinkRefPanel(gc, RefPanelConfig(backend="plink", plink_prefix=str(_CHR22)))
    base_meta = base_panel.load_metadata("22")
    # The reference panel always retains this SNP unless explicitly filtered by
    # a reference SNP list.
    target_pos = int(base_meta["POS"].iloc[0])

    # Region pruning is no longer a reference-panel operation.
    bed = tmp_path / "exclude.bed"
    bed.write_text(f"22\t{target_pos - 1}\t{target_pos}\n", encoding="utf-8")

    with pytest.raises(TypeError, match="exclude_regions_bed"):
        RefPanelConfig(backend="plink", plink_prefix=str(_CHR22), exclude_regions_bed=(str(bed),))
    assert target_pos in set(base_meta["POS"])
