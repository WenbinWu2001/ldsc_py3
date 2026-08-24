import math
import json
from argparse import Namespace
from pathlib import Path
import tempfile
import unittest

import numpy as np

from ldsc.errors import LDSCInputError
from ldsc.quantile_h2 import (
    _read_target_annotation,
    assign_legacy_quantiles,
    compute_quantile_h2,
    compute_standardized_coefficients,
    load_fitted_partitioned_model,
    run_quantile_h2_from_args,
)


class QuantileH2NumericsTest(unittest.TestCase):
    def test_literal_nan_missing_token_is_excluded_explicitly(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            target = Path(tmpdir) / "target.tsv"
            target.write_text(
                "CHR\tPOS\tSNP\ttarget\n1\t10\trs1\tNaN\n1\t20\trs2\t1.5\n",
                encoding="utf-8",
            )

            _, values, excluded, _ = _read_target_annotation([str(target)], "target", "NaN")

            self.assertEqual(excluded.tolist(), [True, False])
            self.assertTrue(math.isnan(values.iloc[0]))
            self.assertEqual(values.iloc[1], 1.5)

    def test_legacy_boundaries_keep_ties_in_lower_quantile(self):
        assignment = assign_legacy_quantiles(np.array([0.0, 1.0, 1.0, 2.0, 3.0]), 2)

        self.assertEqual(assignment.quantile.tolist(), [1, 1, 1, 2, 2])
        self.assertEqual(assignment.lower.tolist(), [0.0, 1.0])
        self.assertEqual(assignment.upper.tolist(), [1.0, 3.0])
        self.assertEqual(assignment.counts.tolist(), [3, 2])

    def test_empty_quantiles_are_rejected(self):
        with self.assertRaisesRegex(LDSCInputError, "empty quantile"):
            assign_legacy_quantiles(np.ones(5), 3)

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
        import pandas as pd

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            ldscore_dir = root / "ldscores"
            ldscore_dir.mkdir()
            model_dir = root / "query"
            model_dir.mkdir()
            (model_dir / "metadata.json").write_text(
                json.dumps({"ldscore_dir": str(ldscore_dir), "retained_ld_columns": ["base", "query"]}),
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
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "diagnostics").mkdir()
            (root / "diagnostics" / "metadata.json").write_text(
                json.dumps({"analysis_type": "cell_type_specific"}),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(LDSCInputError, "aggregate multi-query"):
                load_fitted_partitioned_model(root)

    def test_baseline_only_workflow_writes_complete_result_family(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            ldscore_dir = root / "ldscores"
            ldscore_dir.mkdir()
            (ldscore_dir / "metadata.json").write_text(
                json.dumps(
                    {
                        "artifact_type": "ldscore",
                        "snp_identifier": "rsid",
                        "genome_build": None,
                        "baseline_columns": ["base", "cont"],
                        "query_columns": [],
                        "count_config": {
                            "common_reference_snp_maf_min": 0.05,
                            "common_reference_snp_maf_operator": ">=",
                        },
                        "overlap_config": {"total_common_reference_snps": 4},
                        "counts": [
                            {"column": "base", "common_reference_snp_count": 4.0},
                            {"column": "cont", "common_reference_snp_count": 6.0},
                        ],
                        "annotation_types": {"base": "binary", "cont": "quantitative"},
                    }
                ),
                encoding="utf-8",
            )
            fitted_dir = root / "fitted"
            (fitted_dir / "diagnostics").mkdir(parents=True)
            (fitted_dir / "diagnostics" / "metadata.json").write_text(
                json.dumps(
                    {
                        "artifact_type": "partitioned_h2_result",
                        "analysis_type": "functional_category",
                        "ldscore_dir": str(ldscore_dir),
                        "retained_ld_columns": ["base", "cont"],
                    }
                ),
                encoding="utf-8",
            )
            import pandas as pd

            pd.DataFrame(
                {
                    "category": ["base", "cont"],
                    "coefficient": [0.05, 0.10],
                    "coefficient_se": [0.01, 0.02],
                    "samp_prev": [np.nan, np.nan],
                    "pop_prev": [np.nan, np.nan],
                }
            ).to_csv(fitted_dir / "partitioned_h2.tsv", sep="\t", index=False, na_rep="NaN")
            pd.DataFrame(
                {
                    "delete_block": [0, 1],
                    "base": np.array([0.04, 0.06], dtype=np.float64),
                    "cont": np.array([0.09, 0.11], dtype=np.float64),
                }
            ).to_parquet(fitted_dir / "diagnostics" / "coefficient_delete_values.parquet", index=False)
            annot = root / "annotations.tsv"
            annot.write_text(
                "CHR\tPOS\tSNP\tCM\tbase\tcont\n"
                "1\t10\trs1\t0\t1\t0\n"
                "1\t20\trs2\t0\t1\t1\n"
                "1\t30\trs3\t0\t1\t2\n"
                "1\t40\trs4\t0\t1\t3\n",
                encoding="utf-8",
            )
            target = root / "target.tsv"
            target.write_text(
                "CHR\tPOS\tSNP\ttarget\n"
                "1\t10\trs1\t0\n"
                "1\t20\trs2\t1\n"
                "1\t30\trs3\t2\n"
                "1\t40\trs4\t3\n",
                encoding="utf-8",
            )
            reference = root / "ref.tsv"
            reference.write_text(
                "CHR\tPOS\tSNP\tMAF\n"
                "1\t10\trs1\t0.2\n"
                "1\t20\trs2\t0.2\n"
                "1\t30\trs3\t0.2\n"
                "1\t40\trs4\t0.2\n",
                encoding="utf-8",
            )
            output = root / "quantile"
            result = run_quantile_h2_from_args(
                Namespace(
                    partitioned_h2_result_dir=str(fitted_dir),
                    baseline_annot_sources=[str(annot)],
                    query_annot_sources=None,
                    query_annot_bed_sources=None,
                    query_annot_gene_list_sources=None,
                    gene_coordinate_file=None,
                    control_gene_list_file=None,
                    gene_list_resolution_policy="strict",
                    gene_exclude_regions="none",
                    padding_bp=0,
                    target_annot_sources=[str(target)],
                    target_annotation="target",
                    ref_metadata_sources=[str(reference)],
                    target_missing_value=None,
                    num_quantiles=2,
                    output_dir=str(output),
                    overwrite=False,
                    log_level="INFO",
                )
            )

            self.assertEqual(result.quantile_h2["n_snps"].tolist(), [3, 1])
            self.assertTrue((output / "quantile_h2.tsv").is_file())
            self.assertTrue((output / "standardized_coefficients.tsv").is_file())
            self.assertTrue((output / "diagnostics" / "snp_alignment_issues.tsv.gz").is_file())
            metadata = json.loads((output / "diagnostics" / "metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["artifact_type"], "quantile_h2_result")
            self.assertEqual(metadata["verification_level"], "aggregate_only")


if __name__ == "__main__":
    unittest.main()
