from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc import annotation_builder
from ldsc.annotation_builder import AnnotationBuilder, AnnotationSourceSpec, run_bed_to_annot
from ldsc.config import AnnotationBuildConfig, CommonConfig


FIXTURES = Path(__file__).resolve().parent / "fixtures" / "legacy"
ANNOT_FIXTURE = FIXTURES / "annot_test" / "test.annot"


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
        with mock.patch.object(annotation_builder, "main_make_annot", return_value=7) as patched:
            rc = annotation_builder.main_make_annot(["--help"])
        patched.assert_called_once()
        self.assertEqual(rc, 7)

    def test_run_bed_wrapper_calls_main(self):
        with mock.patch.object(annotation_builder, "main_bed_to_annot", return_value=5) as patched:
            rc = annotation_builder.main_bed_to_annot(["--help"])
        patched.assert_called_once()
        self.assertEqual(rc, 5)
