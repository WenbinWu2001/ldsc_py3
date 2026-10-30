from tests.annotation_fixtures import build_annotation_fixture, make_annotation_bundle, fixture_metadata, fixture_values, fixture_ids, read_annotation_fixture
from pathlib import Path
import contextlib
import gzip
import io
import json
import tempfile
import unittest
from unittest import mock

import numpy as np
import pandas as pd


from ldsc._annotation import tables as kernel_annotation
from ldsc._annotation.sources import prepare_annotation_sources
from ldsc._annotation.storage import AnnotationWorkspace
from ldsc._annotation import builder as annotation_builder
from ldsc._annotation.builder import AnnotationBuilder, AnnotationBundle, run_annotate
from ldsc.config import AnnotationBuildConfig, GlobalConfig
from ldsc.errors import LDSCInputError
from ldsc._annotation.outputs import AnnotationDirectoryWriter


ANNOT_FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "annotation" / "test.annot"


def _write_annot(path: Path, rows: list[tuple], annotation_columns: dict[str, list[int]]) -> None:
    df = pd.DataFrame(rows, columns=["CHR", "BP", "SNP", "CM"])
    for column, values in annotation_columns.items():
        df[column] = values
    compression = "gzip" if path.suffix == ".gz" else None
    df.to_csv(path, sep="\t", index=False, compression=compression)


def _write_gene_catalog(path: Path, *, genome_build: str = "hg38") -> Path:
    path.write_text(
        "gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\n"
        f"ENSG00000186092\tOR4F5\t1\t65419\t71585\t{genome_build}\n"
        f"ENSG00000206503\tHLA-A\t6\t29940000\t29950000\t{genome_build}\n",
        encoding="utf-8",
    )
    return path


class AnnotationBuilderTest(unittest.TestCase):
    def test_preparation_requires_owned_output_workspace(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "out"
            source = Path(tmpdir) / "baseline.annot"
            _write_annot(source, [("1", 10, "rs1", .1)], {"base": [1]})
            builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"))
            spec = AnnotationBuildConfig(baseline_annot_sources=(source,))
            with self.assertRaises(TypeError):
                builder.run(spec)
            with builder.run(spec, output_dir=output) as bundle:
                self.assertEqual(bundle.workspace.path.parent, output)
                self.assertEqual(bundle.read("1").tolist(), [[1.]])
            self.assertEqual(list(output.iterdir()), [])

    def test_bundle_rejects_duplicate_names_across_baseline_and_query(self):
        bundle = make_annotation_bundle(
            metadata=pd.DataFrame({"CHR": ["1"], "POS": [10], "SNP": ["rs1"], "CM": [np.nan]}),
            baseline_annotations=pd.DataFrame({"shared": [1.0]}),
            query_annotations=pd.DataFrame({"shared": [0.0]}),
            baseline_columns=["shared"],
            query_columns=["shared"],
            chromosomes=["1"],
            source_summary={},
        )

        with self.assertRaisesRegex(LDSCInputError, "globally unique.*shared"):
            bundle.validate()

    def test_annotate_output_round_trips_through_annotation_loader(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            expected = pd.DataFrame(
                {
                    "adipose": [True, False, False],
                    "cerebellum": [False, True, False],
                    "cortex": [False, False, True],
                    "whole_blood": [True, True, True],
                }
            )
            bundle = make_annotation_bundle(
                metadata=pd.DataFrame(
                    {
                        "CHR": ["1", "1", "1"],
                        "POS": [10, 20, 30],
                        "SNP": ["rs1", "rs2", "rs3"],
                        "CM": [float("nan")] * 3,
                    }
                ),
                baseline_annotations=pd.DataFrame({"base": [1.0, 1.0, 1.0]}),
                query_annotations=expected,
                baseline_columns=["base"],
                query_columns=expected.columns.tolist(),
                chromosomes=["1"],
                source_summary={},
                config_snapshot=GlobalConfig(snp_identifier="rsid"),
            )

            AnnotationDirectoryWriter().write(bundle, output_dir, overwrite=False, provenance={}, scope={})
            path = output_dir / "query.1.annot.gz"
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                written_rows = [line.split("\t") for line in handle.read().splitlines()[1:]]
            with AnnotationWorkspace(output_dir / "reload") as workspace:
                prepared = prepare_annotation_sources(workspace, [path], [], mode="rsid")
                loaded = AnnotationBundle(prepared.shards, list(prepared.baseline_columns), [], workspace)
                metadata = fixture_metadata(loaded)
                actual = pd.DataFrame(loaded.read("1"), columns=list(prepared.baseline_columns))

        self.assertEqual(metadata["SNP"].tolist(), ["rs1", "rs2", "rs3"])
        self.assertEqual(actual.columns.tolist(), expected.columns.tolist())
        pd.testing.assert_frame_equal(actual, expected.astype("float32"))
        self.assertFalse(actual.isna().any().any())
        self.assertEqual([row[-4:] for row in written_rows], [["1", "0", "0", "1"], ["0", "1", "0", "1"], ["0", "0", "1", "1"]])

    def test_gene_list_projection_matches_catalog_intervals_and_retains_resolution_status(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            baseline = tmpdir / "baseline.annot"
            genes = tmpdir / "immune_genes.txt"
            catalog = _write_gene_catalog(tmpdir / "catalog.tsv")
            baseline.write_text(
                "CHR\tPOS\tSNP\tCM\tbase\n"
                "1\t65419\trs1\t0\t1\n"
                "1\t71586\trs2\t0\t1\n"
                "1\t70000\trs3\t0\t1\n",
                encoding="utf-8",
            )
            genes.write_text("OR4F5\nNOT_A_GENE\n", encoding="utf-8")

            bundle = build_annotation_fixture(AnnotationBuilder(
                GlobalConfig(snp_identifier="rsid"),
                projection_genome_build="hg38",
            ),
                AnnotationBuildConfig(
                    baseline_annot_sources=(baseline,),
                    query_annot_gene_list_sources=(genes,),
                    gene_coordinate_file=catalog,
                    gene_list_resolution_policy="resolved-only",
                )
            )

        self.assertEqual(bundle.query_columns, ["immune_genes"])
        self.assertEqual(bundle.baseline_columns, ["base"])
        self.assertFalse(any(item.input_role == "control" for item in bundle.gene_list_batch.selections))
        self.assertEqual(fixture_values(bundle, 'query')["immune_genes"].dtype, np.dtype("float32"))
        self.assertEqual(fixture_values(bundle, 'query')["immune_genes"].tolist(), [True, False, True])
        self.assertEqual(len(bundle.query_statuses), 1)
        self.assertEqual(bundle.query_statuses[0].status, "warning")
        self.assertEqual(bundle.query_statuses[0].reason, "partial_gene_resolution")
        self.assertEqual(bundle.gene_list_batch.selection("focal", 1).canonical_gene_ids, ("ENSG00000186092",))
        self.assertEqual(
            pd.concat(bundle.gene_list_batch.audit_frames(), ignore_index=True).loc[pd.concat(bundle.gene_list_batch.audit_frames(), ignore_index=True)["disposition"].eq("rejected"), "input_gene"].tolist(),
            ["NOT_A_GENE"],
        )

    def test_gene_control_is_disabled_by_default(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            baseline = tmpdir / "baseline.annot"
            genes = tmpdir / "genes.txt"
            catalog = _write_gene_catalog(tmpdir / "catalog.tsv")
            baseline.write_text(
                "CHR\tPOS\tSNP\tCM\tbase\n1\t65419\trs1\t0\t1\n",
                encoding="utf-8",
            )
            genes.write_text("OR4F5\n", encoding="utf-8")

            bundle = build_annotation_fixture(AnnotationBuilder(
                GlobalConfig(snp_identifier="rsid"),
                projection_genome_build="hg38",
            ),
                AnnotationBuildConfig(
                    baseline_annot_sources=(baseline,),
                    query_annot_gene_list_sources=(genes,),
                    gene_coordinate_file=catalog,
                )
            )

        self.assertEqual(bundle.baseline_columns, ["base"])
        self.assertFalse(any(item.input_role == "control" for item in bundle.gene_list_batch.selections))

    def test_custom_gene_control_uses_the_same_resolution_and_projection_rules(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            baseline = tmpdir / "baseline.annot"
            genes = tmpdir / "focal.txt"
            control = tmpdir / "control.txt"
            catalog = _write_gene_catalog(tmpdir / "catalog.tsv")
            baseline.write_text(
                "CHR\tPOS\tSNP\tCM\tbase\n"
                "1\t65410\trs1\t0\t1\n"
                "1\t65419\trs2\t0\t1\n"
                "1\t71590\trs3\t0\t1\n",
                encoding="utf-8",
            )
            genes.write_text("OR4F5\n", encoding="utf-8")
            control.write_text("OR4F5\nNOT_A_GENE\n", encoding="utf-8")

            bundle = build_annotation_fixture(AnnotationBuilder(
                GlobalConfig(snp_identifier="rsid"),
                projection_genome_build="hg38",
            ),
                AnnotationBuildConfig(
                    baseline_annot_sources=(baseline,),
                    query_annot_gene_list_sources=(genes,),
                    gene_coordinate_file=catalog,
                    control_gene_list_file=control,
                    gene_list_resolution_policy="resolved-only",
                    padding_bp=10,
                )
            )

        self.assertEqual(fixture_values(bundle, 'baseline')["gene_control"].tolist(), [1.0, 1.0, 1.0])
        control_summary = bundle.gene_list_batch.summary[bundle.gene_list_batch.summary["input_role"].eq("control")].iloc[0]
        self.assertEqual(control_summary["rejected_rows"], 1)

    def test_unusable_custom_gene_control_aborts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            baseline = tmpdir / "baseline.annot"
            genes = tmpdir / "focal.txt"
            control = tmpdir / "control.txt"
            catalog = _write_gene_catalog(tmpdir / "catalog.tsv")
            baseline.write_text(
                "CHR\tPOS\tSNP\tCM\tbase\n1\t65419\trs1\t0\t1\n",
                encoding="utf-8",
            )
            genes.write_text("OR4F5\n", encoding="utf-8")
            control.write_text("NOT_A_GENE\n", encoding="utf-8")

            with self.assertRaisesRegex(LDSCInputError, "control.txt.*unmatched_identifier"):
                bundle = build_annotation_fixture(AnnotationBuilder(
                    GlobalConfig(snp_identifier="rsid"),
                    projection_genome_build="hg38",
                ),
                    AnnotationBuildConfig(
                        baseline_annot_sources=(baseline,),
                        query_annot_gene_list_sources=(genes,),
                        gene_coordinate_file=catalog,
                        control_gene_list_file=control,
                    )
                )

    def test_fatal_batched_source_issues_return_diagnostic_bundle_without_projection(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            baseline = tmpdir / "baseline.annot"
            catalog = _write_gene_catalog(tmpdir / "catalog.tsv")
            left_dir = tmpdir / "left"
            right_dir = tmpdir / "right"
            left_dir.mkdir()
            right_dir.mkdir()
            left = left_dir / "same.txt"
            right = right_dir / "same.tsv"
            baseline.write_text(
                "CHR\tPOS\tSNP\tCM\tbase\n1\t65419\trs1\t0\t1\n",
                encoding="utf-8",
            )
            left.write_text("OR4F5\n", encoding="utf-8")
            right.write_text("OR4F5\n", encoding="utf-8")

            with self.assertRaisesRegex(LDSCInputError, "duplicate_query_name"):
                run_annotate(ref_snp_annot_sources=[baseline], query_annot_gene_list_sources=[left, right],
                             gene_coordinate_file=catalog, padding_bp=0, genome_build="hg38",
                             global_config=GlobalConfig(snp_identifier="rsid"), output_dir=tmpdir / "out")
            summary = pd.read_csv(tmpdir / "out/diagnostics/gene_list_resolution_summary.tsv", sep="\t")
            self.assertTrue(summary.source_reasons.str.contains("duplicate_query_name").all())

    def test_gene_control_rejects_reserved_baseline_column_collision(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            baseline = tmpdir / "baseline.annot"
            genes = tmpdir / "genes.txt"
            control = tmpdir / "control.txt"
            catalog = _write_gene_catalog(tmpdir / "catalog.tsv")
            baseline.write_text(
                "CHR\tPOS\tSNP\tCM\tgene_control\n1\t65419\trs1\t0\t1\n",
                encoding="utf-8",
            )
            genes.write_text("OR4F5\n", encoding="utf-8")
            control.write_text("OR4F5\n", encoding="utf-8")

            with self.assertRaisesRegex(Exception, "gene_control"):
                build_annotation_fixture(AnnotationBuilder(
                    GlobalConfig(snp_identifier="rsid"),
                    projection_genome_build="hg38",
                ),
                    AnnotationBuildConfig(
                        baseline_annot_sources=(baseline,),
                        query_annot_gene_list_sources=(genes,),
                        gene_coordinate_file=catalog,
                        control_gene_list_file=control,
                    )
                )

    def test_bed_query_parse_failures_stop_preparation_and_are_all_reported(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            baseline = tmpdir / "baseline.annot"
            valid = tmpdir / "valid.bed"
            malformed = tmpdir / "malformed.bed"
            unordered = tmpdir / "unordered.bed"
            baseline.write_text("CHR\tPOS\tSNP\tCM\tbase\n1\t10\trs1\t0\t1\n", encoding="utf-8")
            valid.write_text("chr1\t9\t11\n", encoding="utf-8")
            malformed.write_text("chr1\tbad\t11\n", encoding="utf-8")
            unordered.write_text("chr1\t11\t9\n", encoding="utf-8")

            with self.assertRaises(LDSCInputError) as caught:
                build_annotation_fixture(AnnotationBuilder(GlobalConfig(snp_identifier="rsid")),
                    AnnotationBuildConfig(
                        baseline_annot_sources=(baseline,),
                        query_annot_bed_sources=(valid, malformed, unordered),
                    )
                )

        message = str(caught.exception)
        self.assertIn("malformed.bed", message)
        self.assertIn("unordered.bed", message)
        self.assertNotIn("valid.bed", message)
        self.assertEqual(sorted(caught.exception.input_issues.source), ["malformed.bed", "unordered.bed"])

    def test_gene_list_projection_is_identical_to_equivalent_bed_with_padding(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            baseline = tmpdir / "baseline.annot"
            genes = tmpdir / "or4f5.txt"
            catalog = _write_gene_catalog(tmpdir / "catalog.tsv")
            bed = tmpdir / "or4f5.bed"
            baseline.write_text(
                "CHR\tPOS\tSNP\tCM\tbase\n"
                "1\t65410\trs1\t0\t1\n"
                "1\t65419\trs2\t0\t1\n"
                "1\t71590\trs3\t0\t1\n"
                "1\t71600\trs4\t0\t1\n",
                encoding="utf-8",
            )
            genes.write_text("OR4F5\n", encoding="utf-8")
            bed.write_text("chr1\t65418\t71585\n", encoding="utf-8")
            builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), projection_genome_build="hg38")

            gene_bundle = build_annotation_fixture(builder,
                AnnotationBuildConfig(
                    baseline_annot_sources=(baseline,),
                    query_annot_gene_list_sources=(genes,),
                    gene_coordinate_file=catalog,
                    padding_bp=10,
                )
            )
            bed_bundle = build_annotation_fixture(builder,
                AnnotationBuildConfig(
                    baseline_annot_sources=(baseline,),
                    query_annot_bed_sources=(bed,),
                    padding_bp=10,
                )
            )

        self.assertEqual(fixture_values(gene_bundle, 'query')["or4f5"].tolist(), [1.0, 1.0, 1.0, 0.0])
        self.assertEqual(
            fixture_values(gene_bundle, 'query')["or4f5"].tolist(),
            fixture_values(bed_bundle, 'query')["or4f5"].tolist(),
        )

    def test_bed_to_annot_parser_genome_build_help_documents_chr_pos_requirement(self):
        stdout = io.StringIO()
        with self.assertRaises(SystemExit):
            with contextlib.redirect_stdout(stdout):
                annotation_builder.parse_annotate_args(["--help"])

        help_text = " ".join(stdout.getvalue().split())
        self.assertIn("Required with chr_pos or chr_pos_allele_aware", help_text)
        self.assertIn("auto checks every grid shard and the gene coordinate catalog", help_text)
        self.assertIn("Concrete builds must agree with input declarations", help_text)

    def test_bed_to_annot_parser_accepts_padding_bp(self):
        args = annotation_builder.parse_annotate_args(
            [
                "--query-annot-bed-sources",
                "query.bed",
                "--ref-snp-annot-sources",
                "baseline.annot.gz",
                "--output-dir",
                "out",
                "--padding-bp",
                "100000",
                "--snp-identifier",
                "rsid",
            ]
        )

        self.assertEqual(args.padding_bp, 100000)

    def test_annotation_build_config_defaults_to_no_bed_padding(self):
        spec = AnnotationBuildConfig(baseline_annot_sources=("baseline.1.annot.gz",))

        self.assertEqual(spec.padding_bp, 0)

    def test_annotation_build_config_accepts_integer_string_bed_padding(self):
        spec = AnnotationBuildConfig(padding_bp="100000")

        self.assertEqual(spec.padding_bp, 100000)

    def test_annotation_build_config_rejects_negative_bed_padding(self):
        with self.assertRaisesRegex(Exception, "padding_bp"):
            AnnotationBuildConfig(padding_bp=-1)

    def test_annotation_build_config_rejects_non_integer_bed_padding(self):
        for value in (1.2, "10.5", "abc"):
            with self.subTest(value=value):
                with self.assertRaisesRegex(Exception, "padding_bp"):
                    AnnotationBuildConfig(padding_bp=value)

    def test_kernel_annotation_does_not_export_workflow_api(self):
        for name in ("run_annotate", "parse_annotate_args"):
            self.assertFalse(hasattr(kernel_annotation, name), f"{name} should live in ldsc._annotation.builder")

    def test_run_builds_bundle(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            query = tmpdir / "query.annot"
            rows = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2), ("2", 30, "rs3", 0.3)]
            _write_annot(base, rows, {"base_a": [1, 0, 1]})
            _write_annot(query, rows, {"query_a": [0, 1, 1]})

            bundle = build_annotation_fixture(builder,
                AnnotationBuildConfig(
                    baseline_annot_sources=(str(base),),
                    query_annot_sources=(str(query),),
                )
            )

            self.assertEqual(bundle.baseline_columns, ["base_a"])
            self.assertEqual(bundle.query_columns, ["query_a"])
            self.assertEqual(bundle.chromosomes, ["1", "2"])
            self.assertEqual(fixture_values(bundle, include_query=True).shape, (3, 2))
            self.assertEqual(fixture_ids(bundle, "rsid"), {"rs1", "rs2", "rs3"})

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

            bundle = build_annotation_fixture(builder)

            self.assertEqual(bundle.baseline_columns, ["base_a"])
            self.assertEqual(fixture_ids(bundle, "rsid"), {"rs1"})

    def test_row_mismatch_raises(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            query = tmpdir / "query.annot"
            rows = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2)]
            _write_annot(base, rows, {"base_a": [1, 0]})
            _write_annot(query, list(reversed(rows)), {"query_a": [0, 1]})

            with self.assertRaises(LDSCInputError):
                build_annotation_fixture(builder,
                    AnnotationBuildConfig(
                        baseline_annot_sources=(str(base),),
                        query_annot_sources=(str(query),),
                    )
                )

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
            bundle = build_annotation_fixture(builder, AnnotationBuildConfig(baseline_annot_sources=(str(base),)))
            self.assertEqual(fixture_ids(bundle, "chr_pos"), {"1:10", "1:20", "2:30"})

    def test_rsid_mode_drops_all_duplicate_identity_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            rows = [("1", 10, "rs1", 0.1), ("1", 20, "rs1", 0.2), ("1", 30, "rs2", 0.3)]
            _write_annot(base, rows, {"base_a": [1, 0, 1]})

            builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
            bundle = build_annotation_fixture(builder, AnnotationBuildConfig(baseline_annot_sources=(str(base),)))

            self.assertEqual(fixture_metadata(bundle)["SNP"].tolist(), ["rs2"])
            self.assertEqual(fixture_values(bundle, 'baseline')["base_a"].tolist(), [1])

    def test_chr_pos_mode_drops_all_duplicate_identity_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            rows = [("1", 10, "rs1", 0.1), ("1", 10, "rs2", 0.2), ("1", 30, "rs3", 0.3)]
            _write_annot(base, rows, {"base_a": [1, 0, 1]})

            builder = AnnotationBuilder(GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"), AnnotationBuildConfig())
            bundle = build_annotation_fixture(builder, AnnotationBuildConfig(baseline_annot_sources=(str(base),)))

            self.assertEqual(fixture_metadata(bundle)["SNP"].tolist(), ["rs3"])
            self.assertEqual(fixture_ids(bundle, "chr_pos"), {"1:30"})

    def test_run_in_allele_aware_mode_accepts_allele_free_annotations(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            rows = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2)]
            _write_annot(base, rows, {"base_a": [1, 0]})

            builder = AnnotationBuilder(
                GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"),
                AnnotationBuildConfig(),
            )
            bundle = build_annotation_fixture(builder, AnnotationBuildConfig(baseline_annot_sources=(str(base),)))

        self.assertNotIn("A1", fixture_metadata(bundle).columns)
        self.assertEqual(fixture_ids(bundle, "chr_pos_allele_aware"), {"1:10", "1:20"})

    def test_run_in_allele_aware_mode_aligns_allele_free_baseline_and_query(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            query = tmpdir / "query.annot"
            rows = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2)]
            _write_annot(base, rows, {"base_a": [1, 0]})
            _write_annot(query, rows, {"query_a": [0, 1]})

            builder = AnnotationBuilder(
                GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"),
                AnnotationBuildConfig(),
            )
            bundle = build_annotation_fixture(builder,
                AnnotationBuildConfig(
                    baseline_annot_sources=(str(base),),
                    query_annot_sources=(str(query),),
                )
            )

        self.assertEqual(bundle.baseline_columns, ["base_a"])
        self.assertEqual(bundle.query_columns, ["query_a"])
        self.assertNotIn("A1", fixture_metadata(bundle).columns)

    def test_run_in_allele_aware_mode_promotes_later_query_alleles_after_alignment(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            query = tmpdir / "query.annot"
            base.write_text(
                "CHR\tBP\tSNP\tCM\tbase_a\n"
                "1\t10\trs1\t0.1\t1\n"
                "1\t20\trs2\t0.2\t0\n",
                encoding="utf-8",
            )
            query.write_text(
                "CHR\tBP\tSNP\tCM\tA1\tA2\tquery_a\n"
                "1\t10\trs1\t0.1\tA\tG\t0\n"
                "1\t20\trs2\t0.2\tC\tT\t1\n",
                encoding="utf-8",
            )

            builder = AnnotationBuilder(
                GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"),
                AnnotationBuildConfig(),
            )
            bundle = build_annotation_fixture(builder,
                AnnotationBuildConfig(
                    baseline_annot_sources=(str(base),),
                    query_annot_sources=(str(query),),
                )
            )

        self.assertEqual(fixture_metadata(bundle)["A1"].tolist(), ["A", "C"])
        self.assertEqual(fixture_metadata(bundle)["A2"].tolist(), ["G", "T"])
        self.assertEqual(bundle.baseline_columns, ["base_a"])
        self.assertEqual(bundle.query_columns, ["query_a"])

    def test_rsid_allele_aware_mode_drops_multiallelic_base_key_clusters(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            base.write_text(
                "CHR\tBP\tSNP\tCM\tA1\tA2\tbase_a\n"
                "1\t10\trs1\t0.1\tA\tC\t1\n"
                "1\t20\trs1\t0.2\tA\tG\t0\n"
                "1\t30\trs2\t0.3\tA\tC\t1\n",
                encoding="utf-8",
            )

            builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid_allele_aware"), AnnotationBuildConfig())
            bundle = build_annotation_fixture(builder, AnnotationBuildConfig(baseline_annot_sources=(str(base),)))

            self.assertEqual(fixture_metadata(bundle)["SNP"].tolist(), ["rs2"])
            self.assertEqual(fixture_values(bundle, 'baseline')["base_a"].tolist(), [1])

    def test_chr_pos_allele_aware_mode_drops_duplicate_effective_keys(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            base.write_text(
                "CHR\tBP\tSNP\tCM\tA1\tA2\tbase_a\n"
                "1\t10\trs1\t0.1\tA\tC\t1\n"
                "1\t10\trs2\t0.2\tC\tA\t0\n"
                "1\t30\trs3\t0.3\tA\tG\t1\n",
                encoding="utf-8",
            )

            builder = AnnotationBuilder(
                GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"),
                AnnotationBuildConfig(),
            )
            bundle = build_annotation_fixture(builder, AnnotationBuildConfig(baseline_annot_sources=(str(base),)))

            self.assertEqual(fixture_metadata(bundle)["SNP"].tolist(), ["rs3"])
            self.assertEqual(fixture_values(bundle, 'baseline')["base_a"].tolist(), [1])

    def test_allele_aware_mode_projects_bed_on_cleaned_grid(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            bed = tmpdir / "query.bed"
            rows = [("1", 10, "rs1", 0.1), ("1", 10, "rs2", 0.2), ("1", 30, "rs3", 0.3)]
            _write_annot(base, rows, {"base_a": [1, 0, 1]})
            bed.write_text("chr1\t29\t30\n", encoding="utf-8")

            builder = AnnotationBuilder(
                GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"),
                AnnotationBuildConfig(),
            )
            bundle = build_annotation_fixture(builder,
                AnnotationBuildConfig(
                    baseline_annot_sources=(str(base),),
                    query_annot_bed_sources=(str(bed),),
                )
            )

            self.assertEqual(fixture_metadata(bundle)["SNP"].tolist(), ["rs3"])
            self.assertEqual(fixture_values(bundle, 'query')[bed.stem].tolist(), [1.0])

    def test_identity_cleanup_raises_when_all_rows_drop(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            rows = [("1", 10, "rs1", 0.1), ("1", 10, "rs2", 0.2)]
            _write_annot(base, rows, {"base_a": [1, 0]})

            builder = AnnotationBuilder(GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"), AnnotationBuildConfig())
            with self.assertRaisesRegex(LDSCInputError, "retained no annotation rows.*Other causes"):
                build_annotation_fixture(builder, AnnotationBuildConfig(baseline_annot_sources=(str(base),)))

    def test_run_with_bed_paths_returns_bundle_with_binary_query_columns(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            bed = tmpdir / "query.bed"
            rows = [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2), ("2", 30, "rs3", 0.3)]
            _write_annot(base, rows, {"base_a": [1, 0, 1]})
            bed.write_text("chr1\t0\t100\tfeature\n", encoding="utf-8")

            bundle = build_annotation_fixture(builder, AnnotationBuildConfig(
                baseline_annot_sources=(str(base),), query_annot_bed_sources=(str(bed),)))

            self.assertEqual(bundle.query_columns, ["query"])
            self.assertEqual(fixture_values(bundle, 'query')["query"].dtype, np.dtype("float32"))
            self.assertEqual(fixture_values(bundle, 'query')["query"].tolist(), [1.0, 1.0, 0.0])
            self.assertEqual(len(fixture_metadata(bundle)), len(fixture_values(bundle, 'query')))

    def test_project_bed_annotations_refuses_existing_output_before_writing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "out"
            output_dir.mkdir()
            existing = output_dir / "query.1.annot.gz"
            existing.write_text("existing\n", encoding="utf-8")
            bundle = make_annotation_bundle(
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

            base = tmpdir / "baseline.annot"
            bed = tmpdir / "query.bed"
            _write_annot(base, [("1", 10, "rs1", .1)], {"base": [1]})
            bed.write_text("1\t9\t11\n")
            with contextlib.nullcontext():
                with self.assertRaisesRegex(FileExistsError, "overwrite"):
                    run_annotate(global_config=builder.global_config,
                        query_annot_bed_sources=(str(bed),),
                        ref_snp_annot_sources=(str(base),),
                        output_dir=output_dir,
                    )

            self.assertEqual(existing.read_text(encoding="utf-8"), "existing\n")
            self.assertFalse((output_dir / "diagnostics" / "annotate.log").exists())

    def test_project_bed_annotations_refuses_stale_query_shard_without_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "out"
            output_dir.mkdir()
            stale = output_dir / "query.2.annot.gz"
            stale.write_text("stale\n", encoding="utf-8")
            bundle = make_annotation_bundle(
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

            base = tmpdir / "baseline.annot"
            bed = tmpdir / "query.bed"
            _write_annot(base, [("1", 10, "rs1", .1)], {"base": [1]})
            bed.write_text("1\t9\t11\n")
            with contextlib.nullcontext():
                with self.assertRaisesRegex(FileExistsError, "overwrite"):
                    run_annotate(global_config=builder.global_config,
                        query_annot_bed_sources=(str(bed),),
                        ref_snp_annot_sources=(str(base),),
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
            bundle = make_annotation_bundle(
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

            base = tmpdir / "baseline.annot"
            bed = tmpdir / "query.bed"
            _write_annot(base, [("1", 10, "rs1", .1)], {"base": [1]})
            bed.write_text("1\t9\t11\n")
            with contextlib.nullcontext():
                run_annotate(global_config=builder.global_config,
                    query_annot_bed_sources=(str(bed),),
                    ref_snp_annot_sources=(str(base),),
                    output_dir=output_dir,
                    overwrite=True,
                )

            self.assertTrue((output_dir / "query.1.annot.gz").exists())
            self.assertTrue((output_dir / "diagnostics" / "metadata.json").exists())
            self.assertFalse((output_dir / "metadata.json").exists())
            self.assertFalse(stale.exists())

    def test_identity_cleanup_sidecar_is_written_for_annotate_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            bed = tmpdir / "query.bed"
            output_dir = tmpdir / "out"
            rows = [("1", 10, "rs1", 0.1), ("1", 10, "rs2", 0.2), ("1", 30, "rs3", 0.3)]
            _write_annot(base, rows, {"base_a": [1, 0, 1]})
            bed.write_text("chr1\t29\t30\n", encoding="utf-8")

            builder = AnnotationBuilder(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                AnnotationBuildConfig(),
            )
            run_annotate(global_config=builder.global_config,
                query_annot_bed_sources=(str(bed),),
                ref_snp_annot_sources=(str(base),),
                output_dir=output_dir,
                padding_bp=100000,
            )

            sidecar = output_dir / "diagnostics" / "dropped_snps" / "dropped.tsv.gz"
            self.assertTrue(sidecar.exists())
            self.assertFalse((output_dir / "metadata.json").exists())
            metadata = json.loads((output_dir / "diagnostics" / "metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["artifact_type"], "annotation_projection")
            self.assertEqual(metadata["snp_identifier"], "chr_pos")
            self.assertEqual(metadata["genome_build"], "hg38")
            self.assertEqual(metadata["padding_bp"], 100000)
            self.assertEqual(
                metadata["files"],
                {"baseline_annotations": ["baseline.1.annot.gz"], "query_annotations": ["query.1.annot.gz"], "dropped_snps": "diagnostics/dropped_snps/dropped.tsv.gz", "query_annotation_status": "diagnostics/query_annotation_status.tsv", "chromosome_scope": "diagnostics/chromosome_scope.json"},
            )
            dropped = pd.read_csv(sidecar, sep="\t", compression="gzip")
            self.assertEqual(dropped["reason"].tolist(), ["duplicate_identity", "duplicate_identity"])
            self.assertEqual(set(dropped["SNP"]), {"rs1", "rs2"})

    def test_annotate_output_has_explicit_missing_cm_and_no_maf(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            bed = tmpdir / "query.bed"
            output_dir = tmpdir / "out"
            _write_annot(base, [("1", 10, "rs1", 0.1), ("1", 30, "rs3", 0.3)], {"base_a": [1, 1]})
            bed.write_text("chr1\t9\t11\n", encoding="utf-8")
            builder = AnnotationBuilder(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                AnnotationBuildConfig(),
            )
            run_annotate(global_config=builder.global_config,
                query_annot_bed_sources=(str(bed),),
                ref_snp_annot_sources=(str(base),),
                output_dir=output_dir,
            )
            output_path = output_dir / "query.1.annot.gz"
            with gzip.open(output_path, "rt", encoding="utf-8") as handle:
                raw_lines = handle.read().splitlines()
            out = pd.read_csv(output_path, sep="\t", compression="gzip")
            self.assertEqual(out.columns[:4].tolist(), ["CHR", "BP", "SNP", "CM"])
            self.assertNotIn("POS", out.columns)
            self.assertIn("CM", out.columns)  # legacy CHR/BP/SNP/CM layout column present
            self.assertTrue(out["CM"].isna().all())  # explicit token normalizes to NaN
            self.assertTrue(all("\tNA\t" in line for line in raw_lines[1:]))
            self.assertTrue(all("\t\t" not in line for line in raw_lines[1:]))
            self.assertTrue(all(len(line.split("\t")) == len(raw_lines[0].split("\t")) for line in raw_lines[1:]))
            self.assertNotIn("MAF", out.columns)  # MAF is never written

    def test_identity_cleanup_sidecar_is_header_only_when_no_rows_drop(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            bed = tmpdir / "query.bed"
            output_dir = tmpdir / "out"
            rows = [("1", 10, "rs1", 0.1), ("1", 30, "rs2", 0.3)]
            _write_annot(base, rows, {"base_a": [1, 0]})
            bed.write_text("chr1\t29\t30\n", encoding="utf-8")

            builder = AnnotationBuilder(GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"), AnnotationBuildConfig())
            run_annotate(global_config=builder.global_config,
                query_annot_bed_sources=(str(bed),),
                ref_snp_annot_sources=(str(base),),
                output_dir=output_dir,
            )

            sidecar = output_dir / "diagnostics" / "dropped_snps" / "dropped.tsv.gz"
            self.assertTrue(sidecar.exists())
            self.assertTrue((output_dir / "diagnostics" / "metadata.json").exists())
            self.assertFalse((output_dir / "metadata.json").exists())
            dropped = pd.read_csv(sidecar, sep="\t", compression="gzip")
            self.assertEqual(
                list(dropped.columns),
                [
                    "CHR",
                    "SNP",
                    "source_pos",
                    "target_pos",
                    "reason",
                    "base_key",
                    "identity_key",
                    "allele_set",
                    "stage",
                ],
            )
            self.assertEqual(len(dropped), 0)

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
                    "ref_snp_annot_sources": str(base),
                    "output_dir": str(output_dir),
                    "overwrite": False,
                    "snp_identifier": "rsid",
                    "genome_build": None,
                    "log_level": "INFO",
                },
            )()

            bundle = annotation_builder.run_annotate_from_args(args)

            self.assertEqual(bundle.query_columns, ["query"])
            self.assertTrue((output_dir / "query.1.annot.gz").exists())
            self.assertTrue((output_dir / "diagnostics" / "metadata.json").exists())
            self.assertFalse((output_dir / "metadata.json").exists())
            self.assertTrue((output_dir / "diagnostics" / "annotate.log").exists())

    def test_parse_fixture_annotation(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        metadata, annotations = read_annotation_fixture(ANNOT_FIXTURE, builder.global_config.snp_identifier)
        self.assertEqual(list(metadata.columns), ["CHR", "POS", "SNP", "CM"])
        self.assertEqual(list(annotations.columns), ["C1", "C2", "C3"])
        self.assertEqual(len(metadata), 3)

    def test_parse_annotation_file_snp_required_in_rsid_mode(self):
        builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"), AnnotationBuildConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "no_snp.annot"
            path.write_text("CHR\tBP\tcoding\n1\t10\t1\n", encoding="utf-8")
            with self.assertRaises(LDSCInputError):
                read_annotation_fixture(path, builder.global_config.snp_identifier)

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

            bundle = build_annotation_fixture(builder,
                AnnotationBuildConfig(
                    baseline_annot_sources=(str(base1), str(base2)),
                    query_annot_sources=(str(query1), str(query2)),
                )
            )

            self.assertEqual(bundle.chromosomes, ["1", "2"])
            self.assertEqual(bundle.baseline_columns, ["base_a"])
            self.assertEqual(bundle.query_columns, ["query_a"])
            self.assertEqual(fixture_ids(bundle, "rsid"), {"rs1", "rs2", "rs3", "rs4"})
            self.assertEqual(fixture_values(bundle).shape, (4, 2))

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

            bundle = build_annotation_fixture(builder,
                AnnotationBuildConfig(
                    baseline_annot_sources=str(tmpdir / "baseline.*.annot.gz"),
                    query_annot_sources=str(tmpdir / "query.*.annot.gz"),
                )
            )

            self.assertEqual(bundle.chromosomes, ["1", "2"])
            self.assertEqual(fixture_ids(bundle, "rsid"), {"rs1", "rs2", "rs3", "rs4"})

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

            bundle = build_annotation_fixture(builder,
                AnnotationBuildConfig(
                    baseline_annot_sources=str(tmpdir / "baseline.@.annot.gz"),
                )
            )

            self.assertEqual(bundle.chromosomes, ["1", "2"])
            self.assertEqual(fixture_ids(bundle, "rsid"), {"rs1", "rs2"})

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

            bundle = build_annotation_fixture(builder,
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

            bundle = build_annotation_fixture(builder,
                AnnotationBuildConfig(
                    baseline_annot_sources=(str(base1), str(base2)),
                )
            )

            self.assertEqual(bundle.chromosomes, ["1", "2"])
            self.assertEqual(bundle.query_columns, [])
            self.assertEqual(fixture_ids(bundle, "rsid"), {"rs1", "rs2", "rs3"})

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

            with self.assertRaisesRegex(LDSCInputError, "mix whole-genome and chromosome-sharded"):
                build_annotation_fixture(builder,
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

            with self.assertRaisesRegex(LDSCInputError, "Query annotation shards do not match baseline chromosome shards"):
                build_annotation_fixture(builder,
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

            with self.assertRaisesRegex(LDSCInputError, "Query annotation shards.*do not match"):
                build_annotation_fixture(builder,
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

            with self.assertRaisesRegex(LDSCInputError, "Multiple baseline annotation files describe chromosome"):
                build_annotation_fixture(builder,
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

            with self.assertRaisesRegex(LDSCInputError, "Annotation SNP rows do not match"):
                build_annotation_fixture(builder,
                    AnnotationBuildConfig(
                        baseline_annot_sources=(str(base1), str(base2)),
                        query_annot_sources=(str(query1), str(query2)),
                    )
                )


class AnnotationWrapperTest(unittest.TestCase):
    def test_run_annotate_returns_annotation_bundle(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            baseline, bed = root / "baseline.annot", root / "query.bed"
            _write_annot(baseline, [("1", 10, "rs1", 0.1), ("1", 20, "rs2", 0.2)], {"base_a": [1, 0]})
            bed.write_text("chr1\t0\t100\tfeature\n", encoding="utf-8")
            with run_annotate(query_annot_bed_sources=[bed], ref_snp_annot_sources=[baseline],
                              output_dir=root / "out", global_config=GlobalConfig(snp_identifier="rsid")) as result:
                self.assertEqual(result.query_columns, ["query"])
                self.assertEqual(result.read("1", columns=["query"]).ravel().tolist(), [1, 1])
                self.assertEqual(result.source_summary["padding_bp"], 0)
                self.assertTrue(Path(result.output_paths["query_annotations"][0]).is_file())

    def test_main_chr_pos_auto_resolves_build_before_running(self):
        from ldsc._genome.build_inference import load_packaged_reference_table
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            baseline, bed = root / "baseline.2.annot", root / "query.bed"
            reference = load_packaged_reference_table()
            rows = reference.loc[reference.CHR.astype(str).eq('2')].head(500)
            rows[['CHR', 'hg19_POS']].rename(columns={'hg19_POS': 'BP'}).assign(base=1).to_csv(baseline, sep='\t', index=False)
            bed.write_text("chr2\t0\t1000000000\n")
            with annotation_builder.main([
                "--query-annot-bed-sources", str(bed), "--ref-snp-annot-sources", str(root / "baseline.*.annot"),
                "--output-dir", str(root / "out"), "--snp-identifier", "chr_pos", "--genome-build", "auto",
            ]) as result:
                self.assertEqual(result.config_snapshot.genome_build, "hg19")
                self.assertEqual(result.read("2", columns=["query"]).ravel().tolist(), [1] * 500)

