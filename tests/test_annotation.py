from pathlib import Path
import contextlib
import io
import sys
import tempfile
import unittest
from unittest import mock

import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc._kernel import annotation as kernel_annotation
from ldsc import annotation_builder
from ldsc.annotation_builder import AnnotationBuilder, AnnotationBundle, run_bed_to_annot
from ldsc.config import AnnotationBuildConfig, GlobalConfig, reset_global_config, set_global_config


ANNOT_FIXTURE = Path(__file__).resolve().parent / "fixtures" / "annotation" / "test.annot"


def _write_annot(path: Path, rows: list[tuple], annotation_columns: dict[str, list[int]]) -> None:
    df = pd.DataFrame(rows, columns=["CHR", "BP", "SNP", "CM"])
    for column, values in annotation_columns.items():
        df[column] = values
    compression = "gzip" if path.suffix == ".gz" else None
    df.to_csv(path, sep="\t", index=False, compression=compression)


class AnnotationBuilderTest(unittest.TestCase):
    def test_bed_to_annot_parser_genome_build_help_documents_chr_pos_requirement(self):
        stdout = io.StringIO()
        with self.assertRaises(SystemExit):
            with contextlib.redirect_stdout(stdout):
                annotation_builder.parse_bed_to_annot_args(["--help"])

        help_text = stdout.getvalue()
        self.assertIn("Required when", help_text)
        self.assertIn("Not used for rsid-family modes", help_text)

    def test_bed_to_annot_parser_rejects_removed_no_batch_flag(self):
        stderr = io.StringIO()
        with self.assertRaises(SystemExit):
            with contextlib.redirect_stderr(stderr):
                annotation_builder.parse_bed_to_annot_args(
                    [
                        "--query-annot-bed-sources",
                        "query.bed",
                        "--baseline-annot-sources",
                        "baseline.annot.gz",
                        "--output-dir",
                        "out",
                        "--no-batch",
                    ]
                )

        self.assertIn("unrecognized arguments: --no-batch", stderr.getvalue())

    def test_annotation_build_config_has_no_gene_set_paths(self):
        spec = AnnotationBuildConfig(baseline_annot_sources=("baseline.1.annot.gz",))
        self.assertFalse(hasattr(spec, "gene_set_paths"))

    def test_gene_set_functions_not_exported(self):
        for name in ("gene_set_to_bed", "make_annot_files", "main_make_annot", "parse_make_annot_args"):
            self.assertFalse(hasattr(annotation_builder, name), f"{name} should not be exported")

    def test_kernel_annotation_does_not_export_workflow_api(self):
        for name in (
            "AnnotationBuilder",
            "AnnotationBundle",
            "run_bed_to_annot",
            "parse_bed_to_annot_args",
            "main_bed_to_annot",
        ):
            self.assertFalse(hasattr(kernel_annotation, name), f"{name} should live in ldsc.annotation_builder")

    def test_run_builds_bundle(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            query = tmpdir / "query.annot"
            rows = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2), ("2", 30, "rs3", 0.3)]
            _write_annot(base, rows, {"base_a": [1, 0, 1]})
            _write_annot(query, rows, {"query_a": [0, 1, 1]})

            bundle = builder.run(
                AnnotationBuildConfig(
                    baseline_annot_sources=(str(base),),
                    query_annot_sources=(str(query),),
                )
            )

            self.assertEqual(bundle.baseline_columns, ["base_a"])
            self.assertEqual(bundle.query_columns, ["query_a"])
            self.assertEqual(bundle.chromosomes, ["1", "2"])
            self.assertEqual(bundle.annotation_matrix(include_query=True).shape, (3, 2))
            self.assertEqual(bundle.reference_snps("rsid"), {"rs1", "rs2", "rs3"})

    def test_run_uses_annotation_config_from_builder_when_source_is_omitted(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            rows = [("1", 10, "rs1", 0.1)]
            _write_annot(base, rows, {"base_a": [1]})
            builder = AnnotationBuilder(
                GlobalConfig(snp_identifier="rsid"),
                AnnotationBuildConfig(baseline_annot_sources=(str(base),)),
            )

            bundle = builder.run()

            self.assertEqual(bundle.baseline_columns, ["base_a"])
            self.assertEqual(bundle.reference_snps("rsid"), {"rs1"})

    def test_row_mismatch_raises(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            query = tmpdir / "query.annot"
            rows = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2)]
            _write_annot(base, rows, {"base_a": [1, 0]})
            _write_annot(query, list(reversed(rows)), {"query_a": [0, 1]})

            with self.assertRaises(ValueError):
                builder.run(
                    AnnotationBuildConfig(
                        baseline_annot_sources=(str(base),),
                        query_annot_sources=(str(query),),
                    )
                )

    def test_run_in_rsid_mode_keeps_all_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            query = tmpdir / "query.annot"
            rows = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2), ("2", 30, "rs3", 0.3)]
            _write_annot(base, rows, {"base_a": [1, 0, 1]})
            _write_annot(query, rows, {"query_a": [0, 1, 1]})

            builder = AnnotationBuilder(
                GlobalConfig(snp_identifier="rsid"),
                AnnotationBuildConfig(),
            )
            bundle = builder.run(
                AnnotationBuildConfig(
                    baseline_annot_sources=(str(base),),
                    query_annot_sources=(str(query),),
                )
            )

            self.assertEqual(bundle.reference_snps("rsid"), {"rs1", "rs2", "rs3"})
            self.assertEqual(len(bundle.metadata), 3)

    def test_run_in_chr_pos_mode_keeps_all_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            rows = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2), ("2", 30, "rs3", 0.3)]
            _write_annot(base, rows, {"base_a": [1, 0, 1]})

            builder = AnnotationBuilder(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                AnnotationBuildConfig(),
            )
            bundle = builder.run(AnnotationBuildConfig(baseline_annot_sources=(str(base),)))
            self.assertEqual(bundle.reference_snps("chr_pos"), {"1:10", "1:20", "2:30"})

    def test_run_with_bed_paths_returns_bundle_with_binary_query_columns(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            bed = tmpdir / "query.bed"
            rows = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2), ("2", 30, "rs3", 0.3)]
            _write_annot(base, rows, {"base_a": [1, 0, 1]})
            bed.write_text("chr1\t0\t100\tfeature\n", encoding="utf-8")

            class _FakeBedTool:
                def __init__(self, path: str):
                    self.path = path

            fake_pybedtools = mock.Mock()
            fake_pybedtools.BedTool = _FakeBedTool

            with mock.patch.object(kernel_annotation, "_get_pybedtools", return_value=fake_pybedtools), mock.patch.object(
                kernel_annotation,
                "_compute_bed_overlap_mask",
                return_value=[True, False, True],
            ):
                bundle = builder.run(
                    AnnotationBuildConfig(
                        baseline_annot_sources=(str(base),),
                        query_annot_bed_sources=(str(bed),),
                    )
                )

            self.assertEqual(bundle.query_columns, ["query"])
            self.assertEqual(bundle.query_annotations["query"].tolist(), [1.0, 0.0, 1.0])
            self.assertEqual(len(bundle.metadata), len(bundle.query_annotations))

    def test_project_bed_annotations_refuses_existing_output_before_writing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "out"
            output_dir.mkdir()
            existing = output_dir / "query.1.annot.gz"
            existing.write_text("existing\n", encoding="utf-8")
            bundle = AnnotationBundle(
                metadata=pd.DataFrame({"CHR": ["1"], "POS": [10], "SNP": ["rs1"], "CM": [0.1]}),
                baseline_annotations=pd.DataFrame({"base": [1.0]}),
                query_annotations=pd.DataFrame({"query": [1.0]}),
                baseline_columns=["base"],
                query_columns=["query"],
                chromosomes=["1"],
                source_summary={},
                config_snapshot=GlobalConfig(snp_identifier="rsid"),
            )
            builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
            builder._workflow_log_path = output_dir / "annotate.log"

            with mock.patch.object(builder, "run", return_value=bundle):
                with self.assertRaisesRegex(FileExistsError, "overwrite"):
                    builder.project_bed_annotations(
                        query_annot_bed_sources=("query.bed",),
                        baseline_annot_sources=("baseline.1.annot.gz",),
                        output_dir=output_dir,
                    )

            self.assertEqual(existing.read_text(encoding="utf-8"), "existing\n")
            self.assertFalse((output_dir / "annotate.log").exists())

    def test_project_bed_annotations_refuses_stale_query_shard_without_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "out"
            output_dir.mkdir()
            stale = output_dir / "query.2.annot.gz"
            stale.write_text("stale\n", encoding="utf-8")
            bundle = AnnotationBundle(
                metadata=pd.DataFrame({"CHR": ["1"], "POS": [10], "SNP": ["rs1"], "CM": [0.1]}),
                baseline_annotations=pd.DataFrame({"base": [1.0]}),
                query_annotations=pd.DataFrame({"query": [1.0]}),
                baseline_columns=["base"],
                query_columns=["query"],
                chromosomes=["1"],
                source_summary={},
                config_snapshot=GlobalConfig(snp_identifier="rsid"),
            )
            builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())

            with mock.patch.object(builder, "run", return_value=bundle):
                with self.assertRaisesRegex(FileExistsError, "overwrite"):
                    builder.project_bed_annotations(
                        query_annot_bed_sources=("query.bed",),
                        baseline_annot_sources=("baseline.1.annot.gz",),
                        output_dir=output_dir,
                    )

            self.assertEqual(stale.read_text(encoding="utf-8"), "stale\n")
            self.assertFalse((output_dir / "query.1.annot.gz").exists())

    def test_project_bed_annotations_overwrite_removes_stale_query_shard(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "out"
            output_dir.mkdir()
            stale = output_dir / "query.2.annot.gz"
            stale.write_text("stale\n", encoding="utf-8")
            bundle = AnnotationBundle(
                metadata=pd.DataFrame({"CHR": ["1"], "POS": [10], "SNP": ["rs1"], "CM": [0.1]}),
                baseline_annotations=pd.DataFrame({"base": [1.0]}),
                query_annotations=pd.DataFrame({"query": [1.0]}),
                baseline_columns=["base"],
                query_columns=["query"],
                chromosomes=["1"],
                source_summary={},
                config_snapshot=GlobalConfig(snp_identifier="rsid"),
            )
            builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())

            with mock.patch.object(builder, "run", return_value=bundle):
                builder.project_bed_annotations(
                    query_annot_bed_sources=("query.bed",),
                    baseline_annot_sources=("baseline.1.annot.gz",),
                    output_dir=output_dir,
                    overwrite=True,
                )

            self.assertTrue((output_dir / "query.1.annot.gz").exists())
            self.assertFalse(stale.exists())

    def test_run_annotate_from_args_writes_workflow_log(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            bed = tmpdir / "query.bed"
            output_dir = tmpdir / "out"
            rows = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2)]
            _write_annot(base, rows, {"base": [1, 1]})
            bed.write_text("chr1\t0\t100\tfeature\n", encoding="utf-8")
            args = type(
                "Args",
                (),
                {
                    "query_annot_bed_sources": str(bed),
                    "baseline_annot_sources": str(base),
                    "output_dir": str(output_dir),
                    "overwrite": False,
                    "snp_identifier": "rsid",
                    "genome_build": None,
                    "log_level": "INFO",
                },
            )()

            class _FakeBedTool:
                def __init__(self, path: str):
                    self.path = path

            fake_pybedtools = mock.Mock()
            fake_pybedtools.BedTool = _FakeBedTool
            with mock.patch.object(kernel_annotation, "_get_pybedtools", return_value=fake_pybedtools), mock.patch.object(
                kernel_annotation,
                "_compute_bed_overlap_mask",
                return_value=[True, False],
            ):
                bundle = annotation_builder.run_annotate_from_args(args)

            self.assertEqual(bundle.query_columns, ["query"])
            self.assertTrue((output_dir / "query.1.annot.gz").exists())
            self.assertTrue((output_dir / "annotate.log").exists())

    def test_parse_fixture_annotation(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        metadata, annotations = builder.parse_annotation_file(ANNOT_FIXTURE)
        self.assertEqual(list(metadata.columns), ["CHR", "POS", "SNP", "CM"])
        self.assertEqual(list(annotations.columns), ["C1", "C2", "C3"])
        self.assertEqual(len(metadata), 3)

    def test_parse_annotation_file_normalizes_bp_header_to_pos(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "legacy_bp.annot"
            path.write_text(
                "CHR\tBP\tSNP\tCM\tbase_a\n1\t10\trs1\t0.1\t1\n",
                encoding="utf-8",
            )

            metadata, annotations = builder.parse_annotation_file(path)

            self.assertEqual(list(metadata.columns), ["CHR", "POS", "SNP", "CM"])
            self.assertEqual(metadata["POS"].tolist(), [10])
            self.assertEqual(list(annotations.columns), ["base_a"])

    def test_annotation_builder_accepts_auto_genome_build(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="chr_pos", genome_build="auto"), AnnotationBuildConfig())

        self.assertEqual(builder.global_config.genome_build, "auto")

    def test_run_auto_bundles_per_chromosome_files(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base1 = tmpdir / "baseline.1.annot.gz"
            base2 = tmpdir / "baseline.2.annot.gz"
            query1 = tmpdir / "baseline_query.1.annot.gz"
            query2 = tmpdir / "baseline_query.2.annot.gz"
            rows1 = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2)]
            rows2 = [("2", 30, "rs3", 0.3), ("2", 40, "rs4", 0.4)]
            _write_annot(base1, rows1, {"base_a": [1, 0]})
            _write_annot(base2, rows2, {"base_a": [0, 1]})
            _write_annot(query1, rows1, {"query_a": [0, 1]})
            _write_annot(query2, rows2, {"query_a": [1, 0]})

            bundle = builder.run(
                AnnotationBuildConfig(
                    baseline_annot_sources=(str(base1), str(base2)),
                    query_annot_sources=(str(query1), str(query2)),
                )
            )

            self.assertEqual(bundle.chromosomes, ["1", "2"])
            self.assertEqual(bundle.baseline_columns, ["base_a"])
            self.assertEqual(bundle.query_columns, ["query_a"])
            self.assertEqual(bundle.reference_snps("rsid"), {"rs1", "rs2", "rs3", "rs4"})
            self.assertEqual(bundle.annotation_matrix().shape, (4, 2))

    def test_run_accepts_single_glob_tokens_for_annotation_groups(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base1 = tmpdir / "baseline.1.annot.gz"
            base2 = tmpdir / "baseline.2.annot.gz"
            query1 = tmpdir / "query.1.annot.gz"
            query2 = tmpdir / "query.2.annot.gz"
            rows1 = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2)]
            rows2 = [("2", 30, "rs3", 0.3), ("2", 40, "rs4", 0.4)]
            _write_annot(base1, rows1, {"base_a": [1, 0]})
            _write_annot(base2, rows2, {"base_a": [0, 1]})
            _write_annot(query1, rows1, {"query_a": [0, 1]})
            _write_annot(query2, rows2, {"query_a": [1, 0]})

            bundle = builder.run(
                AnnotationBuildConfig(
                    baseline_annot_sources=str(tmpdir / "baseline.*.annot.gz"),
                    query_annot_sources=str(tmpdir / "query.*.annot.gz"),
                )
            )

            self.assertEqual(bundle.chromosomes, ["1", "2"])
            self.assertEqual(bundle.reference_snps("rsid"), {"rs1", "rs2", "rs3", "rs4"})

    def test_run_accepts_legacy_suite_tokens_for_annotation_groups(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base1 = tmpdir / "baseline.1.annot.gz"
            base2 = tmpdir / "baseline.2.annot.gz"
            rows1 = [("1", 10, "rs1", 0.1)]
            rows2 = [("2", 20, "rs2", 0.2)]
            _write_annot(base1, rows1, {"base_a": [1]})
            _write_annot(base2, rows2, {"base_a": [0]})

            bundle = builder.run(
                AnnotationBuildConfig(
                    baseline_annot_sources=str(tmpdir / "baseline.@.annot.gz"),
                )
            )

            self.assertEqual(bundle.chromosomes, ["1", "2"])
            self.assertEqual(bundle.reference_snps("rsid"), {"rs1", "rs2"})

    def test_run_prefers_compressed_suite_shards_over_stale_plain_shards(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base1 = tmpdir / "baseline.1.annot.gz"
            base2 = tmpdir / "baseline.2.annot.gz"
            query1_plain = tmpdir / "query.1.annot"
            query1_gz = tmpdir / "query.1.annot.gz"
            query2_gz = tmpdir / "query.2.annot.gz"
            rows1 = [("1", 10, "rs1", 0.1)]
            rows2 = [("2", 20, "rs2", 0.2)]
            _write_annot(base1, rows1, {"base_a": [1]})
            _write_annot(base2, rows2, {"base_a": [0]})
            _write_annot(query1_plain, rows1, {"query_a": [1]})
            _write_annot(query1_gz, rows1, {"query_a": [1], "query_b": [0]})
            _write_annot(query2_gz, rows2, {"query_a": [0], "query_b": [1]})

            bundle = builder.run(
                AnnotationBuildConfig(
                    baseline_annot_sources=str(tmpdir / "baseline.@.annot.gz"),
                    query_annot_sources=str(tmpdir / "query.@.annot.gz"),
                )
            )

            self.assertEqual(bundle.chromosomes, ["1", "2"])
            self.assertEqual(bundle.query_columns, ["query_a", "query_b"])

    def test_run_auto_bundles_baseline_only_per_chromosome_files(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base1 = tmpdir / "baseline.1.annot.gz"
            base2 = tmpdir / "baseline.2.annot.gz"
            rows1 = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2)]
            rows2 = [("2", 30, "rs3", 0.3)]
            _write_annot(base1, rows1, {"base_a": [1, 0]})
            _write_annot(base2, rows2, {"base_a": [0]})

            bundle = builder.run(
                AnnotationBuildConfig(
                    baseline_annot_sources=(str(base1), str(base2)),
                )
            )

            self.assertEqual(bundle.chromosomes, ["1", "2"])
            self.assertEqual(bundle.query_columns, [])
            self.assertEqual(bundle.reference_snps("rsid"), {"rs1", "rs2", "rs3"})

    def test_run_rejects_mixed_whole_genome_and_per_chromosome_inputs(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base_all = tmpdir / "baseline.annot.gz"
            base1 = tmpdir / "baseline.1.annot.gz"
            rows_all = [("1", 10, "rs1", 0.1), ("2", 20, "rs2", 0.2)]
            rows1 = [("1", 10, "rs1", 0.1)]
            _write_annot(base_all, rows_all, {"base_a": [1, 0]})
            _write_annot(base1, rows1, {"base_b": [1]})

            with self.assertRaisesRegex(ValueError, "mixed|chromosome|shard"):
                builder.run(
                    AnnotationBuildConfig(
                        baseline_annot_sources=(str(base_all), str(base1)),
                    )
                )

    def test_run_rejects_missing_query_chromosome_shard(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base1 = tmpdir / "baseline.1.annot.gz"
            base2 = tmpdir / "baseline.2.annot.gz"
            query1 = tmpdir / "baseline_query.1.annot.gz"
            rows1 = [("1", 10, "rs1", 0.1)]
            rows2 = [("2", 20, "rs2", 0.2)]
            _write_annot(base1, rows1, {"base_a": [1]})
            _write_annot(base2, rows2, {"base_a": [0]})
            _write_annot(query1, rows1, {"query_a": [1]})

            with self.assertRaisesRegex(ValueError, "query|chromosome|shard"):
                builder.run(
                    AnnotationBuildConfig(
                        baseline_annot_sources=(str(base1), str(base2)),
                        query_annot_sources=(str(query1),),
                    )
                )

    def test_run_rejects_unsharded_query_when_baseline_is_sharded(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base1 = tmpdir / "baseline.1.annot.gz"
            base2 = tmpdir / "baseline.2.annot.gz"
            query = tmpdir / "query.annot.gz"
            rows1 = [("1", 10, "rs1", 0.1)]
            rows2 = [("2", 20, "rs2", 0.2)]
            _write_annot(base1, rows1, {"base_a": [1]})
            _write_annot(base2, rows2, {"base_a": [0]})
            _write_annot(query, rows1 + rows2, {"query_a": [1, 0]})

            with self.assertRaisesRegex(ValueError, "query|chromosome|shard"):
                builder.run(
                    AnnotationBuildConfig(
                        baseline_annot_sources=(str(base1), str(base2)),
                        query_annot_sources=(str(query),),
                    )
                )

    def test_run_rejects_ambiguous_duplicate_shards(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base1a = tmpdir / "baseline.1.annot.gz"
            base1b = tmpdir / "baseline_extra.1.annot.gz"
            rows = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2)]
            _write_annot(base1a, rows, {"base_a": [1, 0]})
            _write_annot(base1b, rows, {"base_b": [0, 1]})

            with self.assertRaisesRegex(ValueError, "ambiguous|duplicate|chromosome"):
                builder.run(
                    AnnotationBuildConfig(
                        baseline_annot_sources=(str(base1a), str(base1b)),
                    )
                )

    def test_run_still_raises_for_row_mismatch_within_chromosome_shard(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base1 = tmpdir / "baseline.1.annot.gz"
            base2 = tmpdir / "baseline.2.annot.gz"
            query1 = tmpdir / "baseline_query.1.annot.gz"
            query2 = tmpdir / "baseline_query.2.annot.gz"
            rows1 = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2)]
            rows2 = [("2", 30, "rs3", 0.3)]
            bad_rows1 = [("1", 20, "rs2", 0.2), ("1", 10, "rs1", 0.1)]
            _write_annot(base1, rows1, {"base_a": [1, 0]})
            _write_annot(base2, rows2, {"base_a": [0]})
            _write_annot(query1, bad_rows1, {"query_a": [1, 0]})
            _write_annot(query2, rows2, {"query_a": [1]})

            with self.assertRaisesRegex(ValueError, "Annotation SNP rows do not match"):
                builder.run(
                    AnnotationBuildConfig(
                        baseline_annot_sources=(str(base1), str(base2)),
                        query_annot_sources=(str(query1), str(query2)),
                    )
                )


class AnnotationWrapperTest(unittest.TestCase):
    def test_make_annot_wrapper_is_removed(self):
        self.assertFalse(hasattr(annotation_builder, "main_make_annot"))

    def test_run_bed_to_annot_returns_annotation_bundle(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            baseline = tmpdir / "baseline.annot"
            bed = tmpdir / "query.bed"
            _write_annot(
                baseline,
                [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2)],
                {"base_a": [1, 0]},
            )
            bed.write_text("chr1\t0\t100\tfeature\n", encoding="utf-8")

            class _FakeBedTool:
                def __init__(self, path: str):
                    self.path = path

            fake_pybedtools = mock.Mock()
            fake_pybedtools.BedTool = _FakeBedTool

            set_global_config(GlobalConfig(snp_identifier="chr_pos"))
            try:
                with mock.patch.object(kernel_annotation, "_get_pybedtools", return_value=fake_pybedtools), mock.patch.object(
                    kernel_annotation,
                    "_compute_bed_overlap_mask",
                    return_value=[True, False],
                ):
                    result = run_bed_to_annot(
                        query_annot_bed_sources=[str(bed)],
                        baseline_annot_sources=[str(baseline)],
                    )
            finally:
                reset_global_config()

            self.assertIsInstance(result, AnnotationBundle)
            self.assertEqual(result.query_columns, ["query"])

    def test_main_chr_pos_auto_resolves_build_before_running(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            baseline = tmpdir / "baseline.1.annot"
            bed = tmpdir / "query.bed"
            baseline.write_text("CHR\tBP\tSNP\tCM\tbase\n1\t100\trs1\t0.0\t1\n", encoding="utf-8")
            bed.write_text("chr1\t0\t100\tfeature\n", encoding="utf-8")

            with mock.patch.object(
                annotation_builder,
                "sample_frame_from_chr_pattern",
                return_value=(pd.DataFrame({"CHR": ["1"], "POS": [100]}), str(baseline)),
            ) as patched_sample, mock.patch.object(
                annotation_builder,
                "resolve_genome_build",
                return_value="hg19",
            ) as patched_resolve, mock.patch.object(
                annotation_builder,
                "_run_bed_to_annot_with_global_config",
                return_value=mock.sentinel.bundle,
            ) as patched_run:
                result = annotation_builder.main(
                    [
                        "--query-annot-bed-sources",
                        str(bed),
                        "--baseline-annot-sources",
                        str(tmpdir / "baseline.@.annot"),
                        "--output-dir",
                        str(tmpdir / "out"),
                        "--snp-identifier",
                        "chr_pos",
                        "--genome-build",
                        "auto",
                    ]
                )

            self.assertIs(result, mock.sentinel.bundle)
            self.assertEqual(patched_sample.call_args.kwargs["context"], "annotation inputs")
            self.assertEqual(patched_resolve.call_args.args[0], "auto")
            self.assertEqual(patched_run.call_args.kwargs["global_config"].genome_build, "hg19")

    def test_main_bed_to_annot_public_entrypoint_is_removed(self):
        self.assertFalse(hasattr(annotation_builder, "main_bed_to_annot"))
