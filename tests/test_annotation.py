from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import pandas as pd

from ldscore.annotation import AnnotationBuilder, AnnotationSourceSpec, run_bed_to_annot
from ldscore.config import AnnotationBuildConfig, CommonConfig


FIXTURES = Path(__file__).resolve().parent / "fixtures" / "legacy"
ANNOT_FIXTURE = FIXTURES / "annot_test" / "test.annot"
RUN_BED_WRAPPER = Path(__file__).resolve().parents[1] / "utils" / "run_bed_to_annot.py"


def _write_annot(path: Path, rows: list[tuple], annotation_columns: dict[str, list[int]]) -> None:
    df = pd.DataFrame(rows, columns=["CHR", "BP", "SNP", "CM"])
    for column, values in annotation_columns.items():
        df[column] = values
    df.to_csv(path, sep="\t", index=False)


class AnnotationBuilderTest(unittest.TestCase):
    def test_run_builds_bundle(self):
        builder = AnnotationBuilder(CommonConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            query = tmpdir / "query.annot"
            rows = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2), ("2", 30, "rs3", 0.3)]
            _write_annot(base, rows, {"base_a": [1, 0, 1]})
            _write_annot(query, rows, {"query_a": [0, 1, 1]})

            bundle = builder.run(
                AnnotationSourceSpec(
                    baseline_annot_paths=(str(base),),
                    query_annot_paths=(str(query),),
                )
            )

            self.assertEqual(bundle.baseline_columns, ["base_a"])
            self.assertEqual(bundle.query_columns, ["query_a"])
            self.assertEqual(bundle.chromosomes, ["1", "2"])
            self.assertEqual(bundle.annotation_matrix(include_query=True).shape, (3, 2))
            self.assertEqual(bundle.reference_snps("rsid"), {"rs1", "rs2", "rs3"})

    def test_row_mismatch_raises(self):
        builder = AnnotationBuilder(CommonConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            query = tmpdir / "query.annot"
            rows = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2)]
            _write_annot(base, rows, {"base_a": [1, 0]})
            _write_annot(query, list(reversed(rows)), {"query_a": [0, 1]})

            with self.assertRaises(ValueError):
                builder.run(
                    AnnotationSourceSpec(
                        baseline_annot_paths=(str(base),),
                        query_annot_paths=(str(query),),
                    )
                )

    def test_global_rsid_restriction_applies(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            query = tmpdir / "query.annot"
            restrict = tmpdir / "restrict.txt"
            rows = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2), ("2", 30, "rs3", 0.3)]
            _write_annot(base, rows, {"base_a": [1, 0, 1]})
            _write_annot(query, rows, {"query_a": [0, 1, 1]})
            restrict.write_text("rs1\nrs3\n", encoding="utf-8")

            builder = AnnotationBuilder(
                CommonConfig(snp_identifier="rsid", global_snp_restriction_path=str(restrict)),
                AnnotationBuildConfig(),
            )
            bundle = builder.run(
                AnnotationSourceSpec(
                    baseline_annot_paths=(str(base),),
                    query_annot_paths=(str(query),),
                )
            )

            self.assertEqual(bundle.reference_snps("rsid"), {"rs1", "rs3"})
            self.assertEqual(len(bundle.metadata), 2)

    def test_global_chr_pos_restriction_applies(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            rows = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2), ("2", 30, "rs3", 0.3)]
            _write_annot(base, rows, {"base_a": [1, 0, 1]})
            restrict = tmpdir / "restrict.tsv"
            restrict.write_text("CHR\tBP\n1\t20\n2\t30\n", encoding="utf-8")

            builder = AnnotationBuilder(
                CommonConfig(snp_identifier="chr_pos", global_snp_restriction_path=str(restrict)),
                AnnotationBuildConfig(),
            )
            bundle = builder.run(AnnotationSourceSpec(baseline_annot_paths=(str(base),)))
            self.assertEqual(bundle.reference_snps("chr_pos"), {"1:20", "2:30"})

    def test_parse_fixture_annotation(self):
        builder = AnnotationBuilder(CommonConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        metadata, annotations = builder.parse_annotation_file(ANNOT_FIXTURE)
        self.assertEqual(list(metadata.columns), ["CHR", "BP", "SNP", "CM"])
        self.assertEqual(list(annotations.columns), ["C1", "C2", "C3"])
        self.assertEqual(len(metadata), 3)


class AnnotationWrapperTest(unittest.TestCase):
    def test_make_annot_wrapper_calls_main(self):
        import make_annot

        with mock.patch.object(make_annot, "main_make_annot", return_value=7) as patched:
            rc = make_annot.main(["--help"])
        patched.assert_called_once()
        self.assertEqual(rc, 7)

    def test_run_bed_wrapper_calls_main(self):
        spec = spec_from_file_location("run_bed_to_annot_wrapper", RUN_BED_WRAPPER)
        assert spec is not None and spec.loader is not None
        module = module_from_spec(spec)
        spec.loader.exec_module(module)
        with mock.patch.object(module, "main_bed_to_annot", return_value=5) as patched:
            rc = module.main(["--help"])
        patched.assert_called_once()
        self.assertEqual(rc, 5)
