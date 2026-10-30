from pathlib import Path
import gzip
import importlib
import tempfile
import unittest

import pandas as pd


from ldsc.config import GlobalConfig, RefPanelConfig
from ldsc.errors import LDSCConfigError, LDSCInputError
from ldsc._refpanel.adapters import (
    R2MatrixRefPanel,
    PlinkRefPanel,
    RefPanelLoader,
    _snp_id_series_for_matching,
)
from ldsc._refpanel.build_kernel import build_restriction_mask
from tests.ref_panel_helpers import write_r2_matrix_chromosome

import numpy as np


FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "minimal_external_resources" / "plink"
PLINK_PREFIX = FIXTURES / "hm3_chr22_subset"
PLINK_SNPS = ["22:16406147:A:G", "22:16805059:T:C"]


def _has_module(name: str) -> bool:
    try:
        importlib.import_module(name)
        return True
    except Exception:
        return False


class RefPanelBackendTest(unittest.TestCase):
    def test_backend_is_required_and_has_no_auto_value(self):
        with self.assertRaises(TypeError):
            RefPanelConfig(plink_prefix=str(PLINK_PREFIX))
        with self.assertRaisesRegex(LDSCConfigError, "backend"):
            RefPanelConfig(backend="auto", plink_prefix=str(PLINK_PREFIX))

    def test_adapters_reject_a_config_for_the_other_backend(self):
        config = GlobalConfig(snp_identifier="chr_pos", genome_build="hg38")
        with self.assertRaisesRegex(LDSCConfigError, "PlinkRefPanel.*backend='r2_matrix'"):
            PlinkRefPanel(config, RefPanelConfig(backend="r2_matrix", r2_matrix_dir="/path/to/panel/hg38"))
        with self.assertRaisesRegex(LDSCConfigError, "R2MatrixRefPanel.*backend='plink'"):
            R2MatrixRefPanel(config, RefPanelConfig(backend="plink", plink_prefix=str(PLINK_PREFIX)))


class RefPanelLoaderTest(unittest.TestCase):
    def test_loader_selects_plink_backend(self):
        loader = RefPanelLoader(GlobalConfig(snp_identifier="rsid"))
        spec = RefPanelConfig(backend="plink", plink_prefix=str(PLINK_PREFIX))
        panel = loader.load(spec)
        self.assertIsInstance(panel, PlinkRefPanel)

    def test_ref_panel_config_accepts_ref_panel_snps_file(self):
        spec = RefPanelConfig(backend="plink", ref_panel_snps_file="/path/to/snps.txt")
        self.assertEqual(spec.ref_panel_snps_file, "/path/to/snps.txt")

    def test_ref_panel_config_accepts_r2_matrix_dir(self):
        spec = RefPanelConfig(backend="r2_matrix", r2_matrix_dir="/path/to/panel/hg38")
        self.assertEqual(spec.r2_matrix_dir, "/path/to/panel/hg38")

    def test_loader_selects_parquet_backend(self):
        loader = RefPanelLoader(GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))
        spec = RefPanelConfig(backend="r2_matrix", chromosomes=("1",))
        panel = loader.load(spec)
        self.assertIsInstance(panel, R2MatrixRefPanel)

    def test_auto_build_resolution_uses_coordinate_family(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            restriction = Path(tmpdir) / "snps.tsv"
            restriction.write_text("CHR\tHG38_POS\n1\t10\n", encoding="utf-8")
            loader = RefPanelLoader(GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="auto"))

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
            panel = R2MatrixRefPanel(
                GlobalConfig(snp_identifier="rsid_allele_aware"),
                RefPanelConfig(backend="r2_matrix", r2_matrix_dir=tmpdir, ref_panel_snps_file=str(restrict_path)),
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


def _r2_panel(chrom="1", positions=(10, 20), maf=None, cm=None, snps=None):
    """SNP rows for one R2 matrix chromosome, with A1/A2 and MAF as build-r2-matrix records them."""
    n = len(positions)
    panel = pd.DataFrame({
        "CHR": [chrom] * n, "POS": list(positions),
        "SNP": list(snps) if snps is not None else [f"rs{i + 1}" for i in range(n)],
        "A1": ["A"] * n, "A2": ["G", "C", "T", "G", "C", "T"][:n],
        "MAF": list(maf) if maf is not None else [0.2] * n,
    })
    if cm is not None:
        panel["CM"] = list(cm)
    return panel


def _write_r2(build_dir, panel=None, *, snp_identifier="chr_pos", genome_build="hg38"):
    """Write one R2 matrix chromosome with the package writers; return (parquet_path, sidecar_path)."""
    return write_r2_matrix_chromosome(build_dir, _r2_panel() if panel is None else panel,
                                      snp_identifier=snp_identifier, genome_build=genome_build)


def _edit_sidecar(path, edit):
    """Rewrite a sidecar's decompressed text."""
    path.write_bytes(gzip.compress(edit(gzip.decompress(path.read_bytes()).decode()).encode()))


def _r2_ref_panel(build_dir, snp_identifier="chr_pos", genome_build="hg38", **spec):
    return R2MatrixRefPanel(
        GlobalConfig(snp_identifier=snp_identifier, genome_build=genome_build),
        RefPanelConfig(backend="r2_matrix", r2_matrix_dir=str(build_dir), **spec),
    )


class ParquetRefPanelTest(unittest.TestCase):
    def test_missing_sidecar_is_a_hard_error(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_r2(build_dir)[1].unlink()
            with self.assertRaisesRegex(LDSCInputError, "the sidecar .*chr1_meta.tsv.gz' is missing"):
                _r2_ref_panel(build_dir).load_metadata("1")

    def test_r2_matrix_dir_loads_sidecar_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_r2(build_dir, _r2_panel(cm=(0.1, 0.2)), snp_identifier="chr_pos_allele_aware")
            panel = _r2_ref_panel(build_dir, snp_identifier="chr_pos_allele_aware")

            self.assertEqual(panel.available_chromosomes(), ["1"])
            metadata = panel.load_metadata("1")
            self.assertEqual(metadata["SNP"].tolist(), ["rs1", "rs2"])
            self.assertEqual(metadata["CM"].round(3).tolist(), [0.1, 0.2])
            self.assertEqual(metadata["A1"].tolist(), ["A", "A"])
            self.assertEqual(metadata["A2"].tolist(), ["G", "C"])

    def test_package_canonical_r2_rejects_missing_identity_schema_metadata(self):
        import pyarrow.parquet as pq

        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            parquet_path, _ = _write_r2(build_dir)
            table = pq.read_table(parquet_path)
            stripped = {key: value for key, value in table.schema.metadata.items()
                        if key not in {b"ldsc:artifact_type", b"ldsc:snp_identifier", b"ldsc:genome_build"}}
            pq.write_table(table.replace_schema_metadata(stripped), parquet_path)

            with self.assertRaisesRegex(LDSCInputError, "ldsc:artifact_type is missing"):
                _r2_ref_panel(build_dir).build_reader("1")

    def test_package_metadata_sidecar_rejects_missing_identity_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _edit_sidecar(_write_r2(build_dir)[1],
                          lambda text: "".join(line for line in text.splitlines(keepends=True) if not line.startswith("#")))

            with self.assertRaisesRegex(LDSCInputError, "records ldsc:artifact_type=None"):
                _r2_ref_panel(build_dir).load_metadata("1")

    def test_package_metadata_sidecar_rejects_identity_differing_from_its_r2_parquet_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _edit_sidecar(_write_r2(build_dir)[1], lambda text: text.replace(
                "# ldsc:snp_identifier=chr_pos\n", "# ldsc:snp_identifier=chr_pos_allele_aware\n"))

            with self.assertRaisesRegex(LDSCInputError, "come from different `build-r2-matrix` runs"):
                _r2_ref_panel(build_dir).load_metadata("1")

    def test_package_canonical_r2_rejects_snp_identifier_mismatch(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_r2(build_dir, snp_identifier="chr_pos")
            with self.assertRaisesRegex(ValueError, "snp_identifier mismatch"):
                _r2_ref_panel(build_dir, snp_identifier="rsid", genome_build=None).build_reader("1")

    def test_concrete_build_rejects_conflicting_r2_declaration(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_r2(build_dir, genome_build="hg38")
            with self.assertRaisesRegex(LDSCInputError, "Conflicting genome-build evidence"):
                _r2_ref_panel(build_dir, genome_build="hg19").build_reader("1")

    def test_r2_matrix_dir_requires_r2_for_requested_chromosome(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            build_dir.mkdir(parents=True)
            with self.assertRaisesRegex(LDSCInputError, "chr1_r2.parquet"):
                _r2_ref_panel(build_dir).load_metadata("1")

    def test_r2_matrix_dir_top_level_root_descends_by_genome_build(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "panel"
            _write_r2(root / "hg38")
            self.assertEqual(_r2_ref_panel(root).available_chromosomes(), ["1"])

    def test_r2_matrix_dir_top_level_root_rejects_ambiguous_build(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "panel"
            _write_r2(root / "hg19", genome_build="hg19")
            _write_r2(root / "hg38")
            with self.assertRaisesRegex(LDSCInputError, "ambiguous"):
                _r2_ref_panel(root, snp_identifier="rsid", genome_build=None).available_chromosomes()

    def test_r2_matrix_dir_top_level_root_rejects_build_dirs_even_when_only_one_has_r2(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "panel"
            (root / "hg19").mkdir(parents=True)
            _write_r2(root / "hg38")
            with self.assertRaisesRegex(LDSCInputError, "ambiguous"):
                _r2_ref_panel(root, snp_identifier="rsid", genome_build=None).available_chromosomes()

    def test_sidecar_metadata_loading_by_chromosome(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_r2(build_dir, _r2_panel(maf=(0.2, 0.25)))
            _write_r2(build_dir, _r2_panel(chrom="2", positions=(30,), maf=(0.4,), snps=("rs3",)))
            panel = _r2_ref_panel(build_dir)

            self.assertEqual(panel.available_chromosomes(), ["1", "2"])
            metadata = panel.load_metadata("1")
            self.assertEqual(metadata["CHR"].tolist(), ["1", "1"])
            self.assertEqual(metadata["MAF"].tolist(), [0.2, 0.25])
            self.assertEqual(panel.load_metadata("2")["POS"].tolist(), [30])

    def test_sidecar_listing_another_chromosome_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_r2(build_dir, pd.concat([_r2_panel(), _r2_panel(chrom="2", positions=(30,))], ignore_index=True))
            with self.assertRaisesRegex(LDSCInputError, r"lists chromosome\(s\) 1, 2, but it belongs to chromosome 1"):
                _r2_ref_panel(build_dir).load_metadata("1")

    def test_sidecar_row_without_a_chromosome_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _edit_sidecar(_write_r2(build_dir)[1], lambda text: text.replace("\n1\t10\t", "\n\t10\t", 1))
            with self.assertRaisesRegex(LDSCInputError, r"lists chromosome\(s\) 1, nan, but it belongs to chromosome 1"):
                _r2_ref_panel(build_dir).load_metadata("1")

    def test_sidecar_metadata_applies_maf_min_filter(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_r2(build_dir, _r2_panel(positions=(10, 20, 30), maf=(0.2, 0.3, 0.4)))

            metadata = _r2_ref_panel(build_dir, maf_min=0.25).load_metadata("1")

            self.assertEqual(metadata["SNP"].tolist(), ["rs2", "rs3"])
            self.assertTrue((metadata["MAF"] > 0.25).all())

    def test_sidecar_metadata_maf_min_is_inclusive_at_threshold(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_r2(build_dir, _r2_panel(positions=(10, 20, 30), maf=(0.24, 0.25, 0.4)))

            metadata = _r2_ref_panel(build_dir, maf_min=0.25).load_metadata("1")

            # MAF == 0.25 is retained (inclusive >=); 0.24 dropped.
            self.assertEqual(metadata["SNP"].tolist(), ["rs2", "rs3"])

    def test_sidecar_metadata_rejects_invalid_maf_values(self):
        # build-r2-matrix records MAF = freq(A1) <= 0.5 from genotypes for every
        # sidecar row, so a missing or out-of-range value means the artifact was
        # altered; it must never be read as "below threshold", folded into range,
        # or silently skip a requested --maf-min filter.
        sidecars = {
            "partial": ((0.2, None, 0.4), "1 of 3", "rs2"),
            "all": ((None, None), "2 of 2", "rs1"),
            "out_of_range": ((0, 0.7, -0.2, float("inf"), 1.5, 0.5), "4 of 6", "rs2, rs3, rs4, rs5"),
        }
        for label, (maf, count, example) in sidecars.items():
            for maf_min in (None, 0.25):
                with self.subTest(sidecar=label, maf_min=maf_min), tempfile.TemporaryDirectory() as tmpdir:
                    build_dir = Path(tmpdir) / "panel" / "hg38"
                    _write_r2(build_dir, _r2_panel(positions=tuple(10 * (i + 1) for i in range(len(maf))), maf=maf))

                    with self.assertRaisesRegex(
                        LDSCInputError, f"missing, non-numeric, or out-of-range MAF for {count} rows"
                    ) as ctx:
                        _r2_ref_panel(build_dir, maf_min=maf_min).load_metadata("1")
                    self.assertIn(example, str(ctx.exception))

    def test_sidecar_metadata_keeps_boundary_maf_values_unchanged(self):
        # 0 (monomorphic) and 0.5 are the closed ends of the canonical MAF range.
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_r2(build_dir, _r2_panel(positions=(10, 20, 30), maf=(0, 0.5, 0.25)))

            metadata = _r2_ref_panel(build_dir).load_metadata("1")

            self.assertEqual(metadata["MAF"].tolist(), [0.0, 0.5, 0.25])

    def test_sidecar_metadata_reads_back_the_values_the_builder_used(self):
        # build-r2-matrix windows on round_to_sidecar_precision(CM); the reader must
        # return exactly those values for CM and MAF.
        from ldsc._refpanel.r2_format import round_to_sidecar_precision

        cm = np.array([12.3456789012345, 123.456789012345, 0.000123456789012])
        maf = np.array([0.0123456789012, 0.333333333333, 1.23456789012e-05])
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_r2(build_dir, _r2_panel(positions=(10, 20, 30), maf=maf, cm=cm))

            metadata = _r2_ref_panel(build_dir).load_metadata("1")

            np.testing.assert_array_equal(metadata["CM"].to_numpy(dtype=float), round_to_sidecar_precision(cm))
            np.testing.assert_array_equal(metadata["MAF"].to_numpy(dtype=float), round_to_sidecar_precision(maf))

    def test_sidecar_metadata_loading_rejects_auto_genome_build_at_kernel_boundary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_r2(build_dir)
            with self.assertRaisesRegex(AssertionError, "workflow entry"):
                _r2_ref_panel(build_dir, genome_build="auto").load_metadata("1")

    def test_build_reader_requires_the_r2_parquet_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "panel" / "hg38"
            _write_r2(build_dir)[0].unlink()
            with self.assertRaisesRegex(LDSCInputError, "the R2 Parquet file .*chr1_r2.parquet' is missing"):
                _r2_ref_panel(build_dir).build_reader("1")

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
