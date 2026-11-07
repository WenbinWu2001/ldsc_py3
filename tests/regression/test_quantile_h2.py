import math
import json
from argparse import Namespace
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd
import pytest

import ldsc
from ldsc.errors import LDSCConfigError, LDSCInputError
from ldsc._regression.quantile_storage import parse_missing_token, target_chunks
from ldsc._regression.quantile_h2 import (
    compute_quantile_h2,
    compute_standardized_coefficients,
    load_fitted_partitioned_model,
    run_quantile_h2_from_args,
)
from tests.ldscore_helpers import edit_json, write_ldscore_dir
from tests.quantile_fixtures import write_quantile_inputs, write_quantile_ldscores


class QuantileH2NumericsTest(unittest.TestCase):
    def test_literal_nan_missing_token_is_excluded_explicitly(self):
        tmpdir = self.tmp_path
        target = tmpdir / "target.tsv"
        target.write_text(
            "CHR\tPOS\tSNP\ttarget\n1\t10\trs1\tNaN\n1\t20\trs2\t1.5\n",
            encoding="utf-8",
        )

        chunk = next(target_chunks([str(target)], "target", parse_missing_token("NaN")))
        values, excluded = chunk.target_value, chunk.target_excluded

        self.assertEqual(excluded.tolist(), [True, False])
        self.assertTrue(math.isnan(values.iloc[0]))
        self.assertEqual(values.iloc[1], 1.5)

    def test_missing_token_is_parsed_once_into_a_number_or_text(self):
        # A finite number matches numerically equal values; any other token matches its stripped text.
        self.assertEqual(parse_missing_token(" -999 "), -999.0)
        self.assertIsInstance(parse_missing_token("-999"), float)
        self.assertEqual(parse_missing_token(" NaN "), "NaN")
        self.assertEqual(parse_missing_token("inf"), "inf")
        self.assertEqual(parse_missing_token("NA"), "NA")
        self.assertIsNone(parse_missing_token(None))

    def test_projection_matches_hand_worked_jackknife_example(self):
        result = compute_quantile_h2(
            annotation_sums=np.array([[2.0, 1.0], [0.0, 3.0]]),
            tau=np.array([0.1, 0.2]),
            tau_delete=np.array([[0.08, 0.22], [0.12, 0.18]]),
            snp_counts=np.array([2, 2]),
            lower=np.array([0.0, 1.0]),
            upper=np.array([1.0, 2.0]),
        )

        np.testing.assert_allclose(result["h2_obs"], [0.2, 0.7])
        np.testing.assert_allclose(result["h2_obs_se"], [0.04, 0.04])
        np.testing.assert_allclose(result["prop_h2"], [2.0 / 9.0, 7.0 / 9.0])
        np.testing.assert_allclose(result["prop_h2_se"], [2.0 / 45.0, 2.0 / 45.0])
        np.testing.assert_allclose(result["enrichment"], [4.0 / 9.0, 14.0 / 9.0])
        np.testing.assert_allclose(result["enrichment_se"], [4.0 / 45.0, 4.0 / 45.0])
        expected_p = math.erfc(6.25 / math.sqrt(2.0))
        np.testing.assert_allclose(result["enrichment_p"], [expected_p, expected_p])

    def test_tau_star_uses_fixed_full_model_scale(self):
        result = compute_standardized_coefficients(
            annotation_names=["a", "b"],
            annotation_types={"a": "binary", "b": "quantitative"},
            annotation_sd=np.array([0.5, 1.0]),
            tau=np.array([0.1, 0.2]),
            tau_se=np.array([0.01, 0.02]),
            total_h2=0.9,
            n_common_snps=4,
        )

        np.testing.assert_allclose(result["tau_star"], [2.0 / 9.0, 8.0 / 9.0])
        np.testing.assert_allclose(result["tau_star_se"], [1.0 / 45.0, 4.0 / 45.0])
        np.testing.assert_allclose(result["tau_star_p"], result["tau_p"])

    def test_zero_variance_contrast_has_missing_p_under_strict_numpy_policy(self):
        with np.errstate(divide="raise", invalid="raise"):
            before = np.geterr().copy()
            result = compute_quantile_h2(
                annotation_sums=np.array([[1.0, 2.0]]),
                tau=np.array([0.125]),
                tau_delete=np.array([[0.125], [0.125]]),
                snp_counts=np.array([1, 1]),
                lower=np.array([0.0, 1.0]),
                upper=np.array([1.0, 2.0]),
            )
            self.assertEqual(np.geterr(), before)

        np.testing.assert_allclose(result["h2_obs"], [0.125, 0.25])
        np.testing.assert_allclose(result["h2_obs_se"], [0.0, 0.0])
        self.assertTrue(result["enrichment_p"].isna().all())

    def test_negative_nonzero_quantile_total_retains_ratios(self):
        result = compute_quantile_h2(
            annotation_sums=np.array([[1.0, 2.0]]),
            tau=np.array([-0.1]),
            tau_delete=np.array([[-0.09], [-0.11]]),
            snp_counts=np.array([1, 1]),
            lower=np.array([0.0, 1.0]),
            upper=np.array([1.0, 2.0]),
        )

        np.testing.assert_allclose(result["prop_h2"], [1.0 / 3.0, 2.0 / 3.0])
        np.testing.assert_allclose(result["enrichment"], [2.0 / 3.0, 4.0 / 3.0])

    def test_nonpositive_full_model_total_hides_only_tau_star(self):
        result = compute_standardized_coefficients(
            annotation_names=["a"],
            annotation_types={"a": "quantitative"},
            annotation_sd=np.array([1.0]),
            tau=np.array([-0.1]),
            tau_se=np.array([0.02]),
            total_h2=-0.1,
            n_common_snps=10,
        )

        self.assertEqual(result.loc[0, "tau"], -0.1)
        self.assertTrue(math.isnan(result.loc[0, "tau_star"]))
        self.assertTrue(math.isnan(result.loc[0, "tau_star_p"]))

    def test_loader_accepts_one_complete_per_query_model(self):
        root = self.tmp_path
        ldscore_dir = root / "ldscores"
        ldscore_dir.mkdir()
        model_dir = root / "query"
        model_dir.mkdir()
        (model_dir / "metadata.json").write_text(
            json.dumps({"artifact_type": "partitioned_h2_query_result", "ldscore_dir": str(ldscore_dir),
                        "retained_ld_columns": ["base", "query"]}),
            encoding="utf-8",
        )
        pd.DataFrame(
            {
                "category": ["base", "query"],
                "coefficient": [0.1, 0.2],
                "coefficient_se": [0.01, 0.02],
            }
        ).to_csv(model_dir / "partitioned_h2_full.tsv", sep="\t", index=False)
        pd.DataFrame(
            {"delete_block": [0, 1], "base": [0.09, 0.11], "query": [0.19, 0.21]}
        ).to_parquet(model_dir / "coefficient_delete_values.parquet", index=False)

        model = load_fitted_partitioned_model(model_dir)

        self.assertEqual(model.model_type, "baseline_plus_query")
        self.assertEqual(model.annotation_names, ["base", "query"])
        self.assertEqual(model.tau_delete.shape, (2, 2))

    def test_loader_rejects_aggregate_multi_query_root(self):
        root = self.tmp_path
        (root / "diagnostics").mkdir()
        (root / "diagnostics" / "metadata.json").write_text(
            json.dumps({"artifact_type": "partitioned_h2_result", "analysis_type": "cell_type_specific"}),
            encoding="utf-8",
        )

        with self.assertRaisesRegex(LDSCInputError, "aggregate multi-query"):
            load_fitted_partitioned_model(root)



class QuantileH2WorkflowTest(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.args = write_quantile_inputs(self.root)

    def test_baseline_only_workflow_writes_complete_result_family(self):
        output = Path(self.args.output_dir)
        output.mkdir()
        marker = output / "RUN_FAILED.txt"
        marker.write_text("earlier failed overwrite")
        result = run_quantile_h2_from_args(self.args)
        self.assertFalse(marker.exists())
        self.assertEqual(result.quantile_h2["n_snps"].tolist(), [3, 1])
        self.assertTrue((output / "quantile_h2.tsv").is_file())
        self.assertTrue((output / "standardized_coefficients.tsv").is_file())
        self.assertTrue((output / "diagnostics" / "snp_alignment_issues.tsv.gz").is_file())
        metadata = json.loads((output / "diagnostics" / "metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(metadata["artifact_type"], "quantile_h2_result")
        self.assertNotIn("verification_level", metadata)
        self.assertNotIn("fingerprint_canonicalization", metadata)
        log = (output / "diagnostics" / "quantile-h2.log").read_text(encoding="utf-8")
        self.assertNotIn("fingerprint", log.lower())
        np.testing.assert_allclose(result.quantile_h2["h2_obs"], [0.45, 0.35])

    def _write_target(self, values):
        rows = "".join(f"1\t{10 * (k + 1)}\trs{k + 1}\t{value}\n" for k, value in enumerate(values))
        (self.root / "target.tsv").write_text("CHR\tPOS\tSNP\ttarget\n" + rows, encoding="utf-8")

    def test_legacy_boundaries_keep_ties_in_lower_quantile(self):
        # Boundary index i is floor(i * (n - 1) / q + 0.5): sorted values 0, 1, 1, 3 give
        # boundaries 0, 1, 3, and both tied values stay in the lower quantile.
        self._write_target([0, 1, 1, 3])

        result = run_quantile_h2_from_args(self.args)

        self.assertEqual(result.quantile_h2["n_snps"].tolist(), [3, 1])
        self.assertEqual(result.quantile_h2["target_value_lower"].tolist(), [0.0, 1.0])
        self.assertEqual(result.quantile_h2["target_value_upper"].tolist(), [1.0, 3.0])

    def test_quantile_count_is_checked_before_anything_is_written(self):
        # The quantile count is part of the typed settings built before the output preflight, so it no
        # longer waits for annotation staging and SNP alignment.
        for value in (1, 2.5):
            with self.subTest(num_quantiles=value):
                self.args.num_quantiles = value
                with self.assertRaisesRegex(LDSCInputError, "requires --num-quantiles to be an integer of at least 2"):
                    run_quantile_h2_from_args(self.args)
                self.assertFalse(Path(self.args.output_dir).exists())

    def test_gene_list_settings_are_checked_before_anything_is_written(self):
        # The resupplied annotation settings become an AnnotationBuildConfig before the output preflight.
        cases = {
            "gene coordinate file without gene lists": {"gene_coordinate_file": str(self.root / "genes.tsv")},
            "two query source types": {
                "query_annot_sources": [str(self.root / "query.annot")],
                "query_annot_bed_sources": [str(self.root / "query.bed")],
            },
        }
        for label, settings in cases.items():
            with self.subTest(label):
                args = Namespace(**{**vars(self.args), **settings})
                with self.assertRaisesRegex(LDSCConfigError, "Could not construct AnnotationBuildConfig"):
                    run_quantile_h2_from_args(args)
                self.assertFalse(Path(args.output_dir).exists())

    def test_empty_quantiles_are_rejected(self):
        self._write_target([1, 1, 1, 1])
        self.args.num_quantiles = 3

        with self.assertRaisesRegex(LDSCInputError, "empty quantile"):
            run_quantile_h2_from_args(self.args)

    def test_unusable_bed_query_name_writes_the_input_issue_table(self):
        bed = self.root / "SNP.bed"
        bed.write_text("chr1\t5\t15\n", encoding="utf-8")
        self.args.query_annot_bed_sources = [str(bed)]

        with self.assertRaisesRegex(LDSCInputError, "BED query name validation failed"):
            run_quantile_h2_from_args(self.args)

        issues = pd.read_csv(Path(self.args.output_dir) / "diagnostics" / "input_issues.tsv", sep="\t")
        self.assertEqual(issues["reason"].tolist(), ["annotation_name_collision"])
        self.assertEqual(Path(issues["source"].iloc[0]).name, "SNP.bed")

    def test_gene_list_preflight_failure_writes_the_gene_list_diagnostics(self):
        catalog = self.root / "genes.tsv"
        catalog.write_text("gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\nG1\tONE\t1\t15\t25\thg38\n",
                           encoding="utf-8")
        genes = self.root / "query.txt"
        genes.write_text("G1\nG9\n", encoding="utf-8")
        self.args.query_annot_gene_list_sources = [str(genes)]
        self.args.gene_coordinate_file = str(catalog)
        # The fit's saved build is the projection build of the resupplied gene lists.
        edit_json(self.root / "ldscores" / "metadata.json", lambda metadata: metadata.update(genome_build="hg38"))

        with self.assertRaisesRegex(LDSCInputError, "Gene-list catalog preflight rejected"):
            run_quantile_h2_from_args(self.args)

        diagnostics = Path(self.args.output_dir) / "diagnostics"
        summary = pd.read_csv(diagnostics / "gene_list_resolution_summary.tsv", sep="\t")
        self.assertEqual(summary["rejected_rows"].tolist(), [1])
        audit = pd.read_csv(diagnostics / "gene_list_audit.tsv.gz", sep="\t")
        self.assertEqual(audit.loc[audit.disposition.eq("rejected"), "input_gene"].tolist(), ["G9"])

        # Both files belong to the output directory, so a corrected overwrite removes them.
        self.args.query_annot_gene_list_sources = self.args.gene_coordinate_file = None
        self.args.overwrite = True
        run_quantile_h2_from_args(self.args)
        self.assertFalse((diagnostics / "gene_list_resolution_summary.tsv").exists())
        self.assertFalse((diagnostics / "gene_list_audit.tsv.gz").exists())

    def test_invalid_gene_catalog_names_its_defects_and_writes_the_catalog_table(self):
        catalog = self.root / "genes.tsv"
        catalog.write_text("gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\nG1\tONE\t1\t15\t25\tGRCh-38\n",
                           encoding="utf-8")
        genes = self.root / "query.txt"
        genes.write_text("G1\n", encoding="utf-8")
        self.args.query_annot_gene_list_sources = [str(genes)]
        self.args.gene_coordinate_file = str(catalog)

        with self.assertRaises(LDSCInputError) as caught:
            run_quantile_h2_from_args(self.args)

        message = str(caught.exception)
        self.assertIn("is not structurally valid: found 1 defect(s): line 2, genome_build: "
                      "unsupported_genome_build, value 'GRCh-38'", message)
        self.assertIn("troubleshooting.md#common-gene-coordinate-catalog-is-not-valid", message)
        diagnostics = Path(self.args.output_dir) / "diagnostics"
        issues = pd.read_csv(diagnostics / "gene_catalog_issues.tsv", sep="\t", keep_default_na=False)
        self.assertEqual(issues["reason"].tolist(), ["unsupported_genome_build"])
        self.assertEqual(issues["source"].tolist(), ["genes.tsv"])

        # The table belongs to the output directory, so a corrected overwrite removes it.
        self.args.query_annot_gene_list_sources = self.args.gene_coordinate_file = None
        self.args.overwrite = True
        result = run_quantile_h2_from_args(self.args)
        self.assertFalse((diagnostics / "gene_catalog_issues.tsv").exists())
        self.assertNotIn("gene_catalog_issues", result.output_paths)

    def test_run_quantile_h2_names_catalog_defects_and_writes_the_catalog_table(self):
        # The Python form of quantile-h2 writes the same files as the command, this table included.
        catalog = self.root / "genes.tsv"
        catalog.write_text("gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\nG1\tONE\t1\t15\t25\tGRCh-38\n",
                           encoding="utf-8")
        genes = self.root / "query.txt"
        genes.write_text("G1\n", encoding="utf-8")
        options = {
            name: value[0] if isinstance(value, list) else value
            for name, value in vars(self.args).items() if name != "log_level"
        }
        options.update(query_annot_gene_list_sources=str(genes), gene_coordinate_file=str(catalog))

        with self.assertRaises(LDSCInputError) as caught:
            ldsc.run_quantile_h2(**options)

        self.assertIn("is not structurally valid: found 1 defect(s): line 2, genome_build: "
                      "unsupported_genome_build, value 'GRCh-38'", str(caught.exception))
        issues = pd.read_csv(Path(self.args.output_dir) / "diagnostics" / "gene_catalog_issues.tsv", sep="\t",
                             keep_default_na=False)
        self.assertEqual(issues["reason"].tolist(), ["unsupported_genome_build"])

    def test_run_quantile_h2_takes_the_command_option_names(self):
        # Single path tokens may be plain strings, as in the documented example.
        options = {
            name: value[0] if isinstance(value, list) else value
            for name, value in vars(self.args).items() if name != "log_level"
        }

        result = ldsc.run_quantile_h2(**options)

        np.testing.assert_allclose(result.quantile_h2["h2_obs"], [0.45, 0.35])
        self.assertTrue((Path(self.args.output_dir) / "quantile_h2.tsv").is_file())
        self.assertTrue((Path(self.args.output_dir) / "diagnostics" / "quantile-h2.log").is_file())

    def test_aggregate_preserving_value_reassignment_is_accepted(self):
        path = self.root / "annotations.tsv"
        annotations = pd.read_csv(path, sep="\t")
        # Reassign values to SNPs while preserving sums and the full Gram matrix.
        annotations["cont"] = [3, 2, 1, 0]
        annotations.to_csv(path, sep="\t", index=False)

        result = run_quantile_h2_from_args(self.args)

        self.assertEqual(result.quantile_h2["n_snps"].tolist(), [3, 1])
        np.testing.assert_allclose(result.quantile_h2["h2_obs"], [0.75, 0.05])
        self.assertAlmostEqual(result.quantile_h2["h2_obs"].sum(), 0.8)

    def test_reordered_source_rows_keep_values_aligned(self):
        for filename in ("annotations.tsv", "target.tsv", "ref.tsv"):
            path = self.root / filename
            frame = pd.read_csv(path, sep="\t")
            frame.iloc[::-1].to_csv(path, sep="\t", index=False)

        result = run_quantile_h2_from_args(self.args)

        np.testing.assert_allclose(result.quantile_h2["h2_obs"], [0.45, 0.35])

    def test_annotation_sum_mismatch_is_rejected(self):
        path = self.root / "annotations.tsv"
        annotations = pd.read_csv(path, sep="\t")
        annotations["cont"] = [0, 1, 2, 4]
        annotations.to_csv(path, sep="\t", index=False)

        with self.assertRaisesRegex(LDSCInputError, "annotation 'cont'.*sum 7.0, expected 6.0"):
            run_quantile_h2_from_args(self.args)

    def test_overlap_mismatch_with_matching_sums_is_rejected(self):
        path = self.root / "annotations.tsv"
        annotations = pd.read_csv(path, sep="\t")
        annotations["cont"] = [0, 0, 3, 3]
        annotations.to_csv(path, sep="\t", index=False)

        with self.assertRaisesRegex(LDSCInputError, "cross-products.*disagree"):
            run_quantile_h2_from_args(self.args)

    def test_both_snp_universe_counts_are_checked(self):
        path = self.root / "ldscores" / "metadata.json"
        original = path.read_text(encoding="utf-8")
        # Both fixture totals are 4; the common total stays within the all total, as the manifest requires.
        for universe, total in (("all", 5), ("common", 3)):
            with self.subTest(universe=universe):
                metadata = json.loads(original)
                metadata["overlap_config"][f"total_{universe}_reference_snps"] = total
                path.write_text(json.dumps(metadata), encoding="utf-8")
                self.args.output_dir = str(self.root / universe)
                with self.assertRaisesRegex(LDSCInputError, f"different {universe}.*count than the fitted"):
                    run_quantile_h2_from_args(self.args)

    def test_duplicate_reference_and_target_identities_are_rejected(self):
        for filename in ("ref.tsv", "target.tsv"):
            with self.subTest(source=filename):
                path = self.root / filename
                original = path.read_text(encoding="utf-8")
                frame = pd.read_csv(path, sep="\t")
                pd.concat([frame, frame.iloc[[0]]]).to_csv(path, sep="\t", index=False)
                self.args.output_dir = str(self.root / filename.replace(".tsv", "_duplicate"))
                with self.assertRaisesRegex(LDSCInputError, "duplicate effective SNP identities"):
                    run_quantile_h2_from_args(self.args)
                path.write_text(original, encoding="utf-8")

    def test_missing_target_row_is_rejected(self):
        path = self.root / "target.tsv"
        target = pd.read_csv(path, sep="\t")
        target.iloc[:-1].to_csv(path, sep="\t", index=False)

        with self.assertRaisesRegex(LDSCInputError, "does not cover the common reference SNPs"):
            run_quantile_h2_from_args(self.args)
        issues = pd.read_csv(
            Path(self.args.output_dir) / "diagnostics" / "snp_alignment_issues.tsv.gz", sep="\t"
        )
        self.assertEqual(issues["SNP"].tolist(), ["rs4"])
        self.assertEqual(issues["issue"].tolist(), ["missing_target_annotation"])

    def test_same_named_fitted_target_values_must_match(self):
        path = self.root / "target.tsv"
        target = pd.read_csv(path, sep="\t").rename(columns={"target": "cont"})
        target["cont"] = [3, 2, 1, 0]
        target.to_csv(path, sep="\t", index=False)
        self.args.target_annot_column = "cont"

        with self.assertRaisesRegex(LDSCInputError, "matches a fitted annotation name but its values differ"):
            run_quantile_h2_from_args(self.args)

    def test_matching_fitted_target_is_accepted(self):
        self.args.target_annot_column = "cont"
        self.args.target_annot_sources = self.args.baseline_annot_sources

        result = run_quantile_h2_from_args(self.args)

        np.testing.assert_allclose(result.quantile_h2["h2_obs"], [0.45, 0.35])

    def test_counts_are_usable_without_an_overlap_artifact(self):
        write_quantile_ldscores(self.root, overlap=False, overwrite=True)

        result = run_quantile_h2_from_args(self.args)

        np.testing.assert_allclose(result.quantile_h2["h2_obs"], [0.45, 0.35])


if __name__ == "__main__":
    unittest.main()


def test_loader_rejects_a_relative_ldscore_directory(tmp_path, monkeypatch):
    # A same-named directory in the working directory must not be read silently.
    (tmp_path / "elsewhere" / "ldscores").mkdir(parents=True)
    model_dir = tmp_path / "fit" / "query"
    (model_dir.parent / "ldscores").mkdir(parents=True)
    model_dir.mkdir()
    (model_dir / "metadata.json").write_text(
        json.dumps({"artifact_type": "partitioned_h2_query_result", "ldscore_dir": "ldscores",
                    "retained_ld_columns": ["base", "query"]}), encoding="utf-8")
    pd.DataFrame({"category": ["base", "query"], "coefficient": [0.1, 0.2], "coefficient_se": [0.01, 0.02]}
                 ).to_csv(model_dir / "partitioned_h2_full.tsv", sep="\t", index=False)
    pd.DataFrame({"delete_block": [0, 1], "base": [0.09, 0.11], "query": [0.19, 0.21]}
                 ).to_parquet(model_dir / "coefficient_delete_values.parquet", index=False)
    monkeypatch.chdir(tmp_path / "elsewhere")
    with pytest.raises(LDSCInputError, match="relative"):
        load_fitted_partitioned_model(model_dir)


def test_loader_requires_retained_ld_columns(tmp_path):
    # partitioned-h2 records the coefficient order whenever it writes the delete table.
    model_dir = tmp_path / "fit" / "query"
    model_dir.mkdir(parents=True)
    (model_dir / "metadata.json").write_text(
        json.dumps({"artifact_type": "partitioned_h2_query_result", "ldscore_dir": str(tmp_path)}), encoding="utf-8")
    pd.DataFrame({"category": ["base", "query"], "coefficient": [0.1, 0.2], "coefficient_se": [0.01, 0.02]}
                 ).to_csv(model_dir / "partitioned_h2_full.tsv", sep="\t", index=False)
    pd.DataFrame({"delete_block": [0, 1], "base": [0.09, 0.11], "query": [0.19, 0.21]}
                 ).to_parquet(model_dir / "coefficient_delete_values.parquet", index=False)
    with pytest.raises(LDSCInputError, match="retained_ld_columns"):
        load_fitted_partitioned_model(model_dir)


@pytest.mark.parametrize("edit, message", [
    (lambda metadata: metadata.update(artifact_type="annotation"), "expected artifact_type='ldscore'"),
    (lambda metadata: metadata.pop("artifact_type"), "lacks `artifact_type`"),
    (lambda metadata: metadata["overlap_config"].update(common_maf_min=0.01, common_maf_operator=">"),
     "count_config and overlap_config disagree"),
    (lambda metadata: metadata.pop("count_config"), "missing a usable `count_config`"),
])
def test_ldscore_directories_the_loader_rejects_stop_before_computation(tmp_path, edit, message):
    # Each edit makes load_ldscore_from_dir reject the directory; quantile-h2 used to accept all four.
    args = write_quantile_inputs(tmp_path)
    edit_json(tmp_path / "ldscores" / "metadata.json", edit)
    # Annotation staging would fail on this value, so the LD-score error shows the check ran before it.
    annotations = Path(args.baseline_annot_sources[0])
    annotations.write_text(annotations.read_text(encoding="utf-8").replace("\t1\t3\n", "\t1\tnot-a-number\n"),
                           encoding="utf-8")

    with pytest.raises(LDSCInputError, match=message):
        run_quantile_h2_from_args(args)


@pytest.mark.parametrize("damage", ["missing overlap", "missing query batch", "corrupt metadata"])
def test_a_damaged_ldscore_directory_fails_in_the_preflight(tmp_path, damage):
    args = write_quantile_inputs(tmp_path)
    ldscore_dir = tmp_path / "ldscores"
    missing = {"missing overlap": ("overlap", "ldscore.overlap.parquet"),
               "missing query batch": ("query", "ldscore.query.parquet")}.get(damage)
    if damage == "missing query batch":
        # Rewrite the LD scores with one query batch beside the fitted baseline annotations.
        baseline = pd.read_parquet(ldscore_dir / "ldscore.baseline.parquet")
        write_ldscore_dir(ldscore_dir, baseline, baseline_columns=["base", "cont"],
                          query=baseline.loc[:, ["CHR", "SNP", "POS", "A1", "A2"]].assign(extra=1.0),
                          query_columns=["extra"], counts={"base": (4.0, 4.0), "cont": (6.0, 6.0), "extra": (4.0, 4.0)},
                          snp_identifier="rsid", overwrite=True)
    if missing is not None:
        (ldscore_dir / missing[1]).unlink()
    else:
        (ldscore_dir / "metadata.json").write_text("{not json", encoding="utf-8")

    with pytest.raises(LDSCInputError, match="Input preflight.*LD scores"):
        run_quantile_h2_from_args(args)

    issues = pd.read_csv(Path(args.output_dir) / "diagnostics" / "input_issues.tsv", sep="\t")
    if missing is not None:
        # The manifest check and the declared-file check each name the missing file.
        role, name = missing
        assert issues["input_role"].tolist() == ["LD scores", role]
        assert issues["source"].tolist() == [str(ldscore_dir), str(ldscore_dir / name)]
        assert issues["details"].iloc[0].startswith("LDSCInputError: Could not ")
        assert name in issues["details"].iloc[0]
        # The directory row reports an inconsistent manifest; the file row reports the missing file.
        assert issues["reason"].tolist() == ["invalid_required_input", "missing_required_input"]
    else:
        # Undecodable metadata is reported once; the declared files cannot be listed.
        assert issues["input_role"].tolist() == ["LD scores"]
        assert issues["source"].tolist() == [str(tmp_path / "ldscores")]
    assert not (Path(args.output_dir) / "quantile_h2.tsv").exists()


def test_the_ldscore_preflight_reports_manifest_and_file_defects_in_one_run(tmp_path):
    args = write_quantile_inputs(tmp_path)
    ldscore_dir = tmp_path / "ldscores"
    (ldscore_dir / "diagnostics").mkdir(exist_ok=True)
    (ldscore_dir / "diagnostics" / "notes.tsv").write_text("x\n", encoding="utf-8")
    edit_json(ldscore_dir / "metadata.json", lambda metadata: (metadata.pop("count_config"),
                                                              metadata["files"].update(notes="diagnostics/notes.tsv")))
    (ldscore_dir / "diagnostics" / "notes.tsv").unlink()

    with pytest.raises(LDSCInputError, match="2 input issues"):
        run_quantile_h2_from_args(args)

    issues = pd.read_csv(Path(args.output_dir) / "diagnostics" / "input_issues.tsv", sep="\t")
    assert issues["input_role"].tolist() == ["LD scores", "notes"]
    assert "count_config" in issues["details"].iloc[0]
    assert issues["source"].iloc[1] == str(ldscore_dir / "diagnostics" / "notes.tsv")


def test_an_ldscore_directory_is_not_read_as_a_fitted_model(tmp_path):
    # A directory holding its own metadata.json is a per-query folder only under diagnostics/query_annotations.
    args = write_quantile_inputs(tmp_path)
    args.partitioned_h2_result_dir = str(tmp_path / "ldscores")

    with pytest.raises(LDSCInputError, match="expected a partitioned-h2 fitted model.*artifact_type 'ldscore'"):
        run_quantile_h2_from_args(args)


@pytest.mark.parametrize("target", ["query", "shared"])
def test_undecodable_fitted_metadata_is_an_input_error(tmp_path, target):
    from ldsc._regression.outputs import (
        PARTITIONED_H2_COLUMNS, PartitionedH2DirectoryWriter, PartitionedH2OutputConfig, stage_partitioned_h2_fit,
    )

    args = write_quantile_inputs(tmp_path)
    fitted = Path(args.partitioned_h2_result_dir)
    full = pd.read_csv(fitted / "partitioned_h2.tsv", sep="\t").reindex(columns=PARTITIONED_H2_COLUMNS)
    deletes = pd.read_parquet(fitted / "diagnostics/coefficient_delete_values.parquet")
    artifacts = stage_partitioned_h2_fit(tmp_path / "private_fit", full, deletes, {"n_snps": 4})
    shared = json.loads((fitted / "diagnostics/metadata.json").read_text(encoding="utf-8"))
    published = tmp_path / "published"
    PartitionedH2DirectoryWriter().write(full.iloc[[1]], PartitionedH2OutputConfig(published, write_per_query_results=True),
                                         {"cont": artifacts}, dict(shared, analysis_type="cell_type_specific"))
    model_dir = published / "diagnostics/query_annotations/0001_cont"
    damaged = model_dir / "metadata.json" if target == "query" else published / "diagnostics/metadata.json"
    damaged.write_text("{not json", encoding="utf-8")

    with pytest.raises(LDSCInputError, match=f"could not read.*{damaged.name}"):
        load_fitted_partitioned_model(model_dir)
