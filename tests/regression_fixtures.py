"""Shared inputs for the regression dataset, fitting and command tests.

``RegressionFixtures`` is a mixin for ``unittest.TestCase`` classes: it resets
the registered ``GlobalConfig``, runs calculations in this process with the
local query-worker seam, and builds small summary-statistics tables, LD-score
results and directories, and kernel results.
"""

import contextlib
from dataclasses import replace
import json
from pathlib import Path
import tempfile
from unittest import mock

import numpy as np
import pandas as pd

from ldsc.config import GlobalConfig, reset_global_config
from ldsc._ldscore.results import LDScoreResult
from ldsc._regression import estimators as reg
from ldsc._sumstats.table import SumstatsTable
from tests.annotation_fixtures import make_annotation_bundle
from tests.cli_parsing import parse_command
from tests.ldscore_helpers import write_one_snp_ldscore_dir


@contextlib.contextmanager
def stub_h2_kernel(estimator):
    """Stub numerical estimation while exercising real selection and diagnostics."""
    if "coef" not in vars(estimator):
        estimator.coef = np.array([0.001])
    if "intercept" not in vars(estimator):
        estimator.intercept = np.array([1.0])
    weights = reg.Hsq.weights
    with mock.patch.object(reg, "Hsq", wraps=reg.Hsq, return_value=estimator) as kernel:
        kernel.weights = weights
        yield kernel


def fake_total_h2_hsq(**overrides):
    """A fitted-estimator stub carrying the total-h2 summary fields; ``overrides`` replace or add fields."""
    return mock.Mock(**{
        "tot": np.array([0.1]), "tot_se": np.array([0.01]), "intercept": np.array([1.0]), "intercept_se": 0.01,
        "mean_chisq": np.array([1.1]), "lambda_gc": np.array([1.0]), "ratio": 0.0, "ratio_se": 0.0, **overrides,
    })


def fake_partitioned_hsq(n_annot=2, **overrides):
    """A fitted partitioned-estimator stub over one (baseline only) or two (baseline + query) annotations."""
    if n_annot == 1:
        fields = {
            "tot": np.array([0.25]), "tot_se": np.array([0.03]),
            "coef": np.array([0.025]), "coef_cov": np.array([[9e-6]]), "coef_se": np.array([0.003]),
            "cat": np.array([0.25]), "cat_se": np.array([0.03]),
            "prop": np.array([1.0]), "prop_cov": np.array([[0.0]]), "prop_se": np.array([0.0]),
            "enrichment": np.array([1.0]),
        }
    else:
        fields = {
            "tot": np.array([0.3]), "tot_se": np.array([0.03]),
            "coef": np.array([0.0, 1.0]), "coef_cov": np.diag([0.01, 0.04]), "coef_se": np.array([0.1, 0.2]),
            "cat": np.array([0.0, 0.3]), "cat_se": np.array([0.01, 0.03]),
            "prop": np.array([0.0, 0.4]), "prop_cov": np.diag([0.0001, 0.0016]), "prop_se": np.array([0.01, 0.04]),
            "enrichment": np.array([0.0, 2.0]),
        }
    return mock.Mock(**{**fields, "n_blocks": 200, "n_annot": n_annot, **overrides})


def parsed_regression_args(command, options):
    """Parse ``options`` (option dest to value) with the ``ldsc <command>`` parser, so the others take its defaults.

    ``None`` and ``False`` leave an option unset, ``True`` passes a flag, and a
    list passes several values. ``output_dir=None`` is set after parsing, for
    checks that a programmatic namespace without an output directory is refused.
    """
    argv = [command]
    for dest, value in options.items():
        if value is None or value is False:
            continue
        flag = "--" + dest.replace("_", "-")
        argv += [flag] if value is True else [flag, *map(str, value)] if isinstance(value, list) else [flag, str(value)]
    if options.get("output_dir") is None:
        argv += ["--output-dir", "unused"]
    args = parse_command(argv)
    if options.get("output_dir") is None:
        args.output_dir = None
    return args


class RegressionFixtures:
    def setUp(self):
        reset_global_config()
        # This class tests model selection/formatting with mocked estimators.
        # Real process behavior is covered in test_partitioned_query_workers.
        from tests.worker_fixtures import calculation_body, LocalQueryWorkers
        self.enterContext(calculation_body())
        self.enterContext(mock.patch("ldsc._regression.partitioned_parallel._QueryWorkers", LocalQueryWorkers))

    def tearDown(self):
        reset_global_config()

    def temporary_output_dir(self):
        return Path(self.enterContext(tempfile.TemporaryDirectory()))

    def write_sumstats_sidecar(
        self,
        path: Path,
        *,
        trait_name: str | None = "trait",
        snp_identifier: str = "rsid",
        genome_build: str | None = None,
    ) -> None:
        path.write_text(
            json.dumps(
                {
                    "artifact_type": "sumstats",
                    "files": {},
                    "trait_name": trait_name,
                    "snp_identifier": snp_identifier,
                    "genome_build": genome_build,
                }
            ),
            encoding="utf-8",
        )

    def write_footer_sumstats_parquet(
        self,
        path: Path,
        *,
        frame: pd.DataFrame | None = None,
        trait_name: str | None = "trait",
        snp_identifier: str = "rsid",
        genome_build: str | None = None,
    ) -> None:
        """Write a self-describing sumstats parquet carrying identity footer metadata."""
        from ldsc._sumstats.outputs import curated_sumstats_frame, write_sumstats_parquet
        from ldsc._sumstats.table import encode_sumstats_footer

        if frame is None:
            frame = pd.DataFrame({"SNP": ["rs1"], "Z": [1.0], "N": [1000.0]})
        footer = encode_sumstats_footer(
            GlobalConfig(snp_identifier=snp_identifier, genome_build=genome_build), trait_name
        )
        write_sumstats_parquet(curated_sumstats_frame(frame), path, footer)

    def make_ldscore_result(self):
        baseline_table = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "SNP": ["rs1", "rs2", "rs3"],
                "POS": [10, 20, 30],
                "regression_ld_scores": [2.0, 2.0, 2.0],
                "base": [1.0, 2.0, 3.0],
            }
        )
        query_table = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "SNP": ["rs1", "rs2", "rs3"],
                "POS": [10, 20, 30],
                "query1": [1.0, 2.0, 3.0],
                "query2": [0.5, 1.5, 2.5],
            }
        )
        return LDScoreResult(
            baseline_table=baseline_table,
            query_table=query_table,
            count_records=[
                {
                    "group": "baseline",
                    "column": "base",
                    "all_reference_snp_count": 10.0,
                    "common_reference_snp_count": 8.0,
                },
                {
                    "group": "query",
                    "column": "query1",
                    "all_reference_snp_count": 20.0,
                    "common_reference_snp_count": 18.0,
                },
                {
                    "group": "query",
                    "column": "query2",
                    "all_reference_snp_count": 30.0,
                    "common_reference_snp_count": 28.0,
                },
            ],
            baseline_columns=["base"],
            query_columns=["query1", "query2"],
            config_snapshot=GlobalConfig(snp_identifier="rsid"),
            overlap=self._overlap_fixture(),
        )

    def _overlap_fixture(self):
        # Consistent overlap: diag(O) equals each annotation's M_annot, and
        # M_tot (total_*_reference_snps) is larger than every category count.
        from ldsc._annotation.overlap import OverlapContribution
        from ldsc._ldscore.overlap_table import LDScoreOverlap

        return LDScoreOverlap.from_contribution(
            OverlapContribution(
                baseline_block_all=np.array([[10.0, 5.0, 8.0]]),   # base x [base, query1, query2]
                baseline_block_common=np.array([[8.0, 4.0, 7.0]]),
                query_diagonal_all=np.array([20.0, 30.0]),          # query1, query2
                query_diagonal_common=np.array([18.0, 28.0]),
                n_all=50, n_common=45,
            ),
            baseline_columns=["base"],
            query_columns=["query1", "query2"],
        )

    def make_sumstats_table(self):
        return SumstatsTable(
            data=pd.DataFrame(
                {
                    "SNP": ["rs1", "rs2", "rs3"],
                    "Z": [2.0, 1.0, 0.5],
                    "N": [1000.0, 1000.0, 1000.0],
                    "A1": ["A", "C", "G"],
                    "A2": ["G", "T", "A"],
                }
            ),
            has_alleles=True,
            source_path="sumstats.gz",
            trait_name="trait",
            provenance={},
            config_snapshot=GlobalConfig(snp_identifier="rsid"),
        )

    def make_sumstats_table_from_frame(self, frame, *, snp_identifier="rsid", trait_name="trait"):
        return SumstatsTable(
            data=frame,
            has_alleles={"A1", "A2"}.issubset(frame.columns),
            source_path=f"{trait_name}.sumstats.gz",
            trait_name=trait_name,
            provenance={},
            config_snapshot=GlobalConfig(snp_identifier=snp_identifier, genome_build="hg38"),
        )

    def make_ldscore_result_from_frame(self, frame, *, snp_identifier="rsid"):
        return LDScoreResult(
            baseline_table=frame,
            query_table=None,
            count_records=[
                {
                    "group": "baseline",
                    "column": "base",
                    "all_reference_snp_count": 10.0,
                    "common_reference_snp_count": 8.0,
                }
            ],
            baseline_columns=["base"],
            query_columns=[],
            config_snapshot=GlobalConfig(snp_identifier=snp_identifier, genome_build="hg38"),
        )

    def make_rg_kernel_result(self, *, rg=0.25, rg_se=0.05, z=5.0, p=0.01, h2_2=0.30):
        hsq1 = mock.Mock(
            tot=np.array([0.20]),
            tot_se=np.array([0.02]),
            intercept=np.array([1.01]),
            intercept_se=0.01,
            mean_chisq=np.array([1.20]),
            lambda_gc=np.array([1.05]),
            ratio=0.05,
            ratio_se=0.01,
        )
        hsq2 = mock.Mock(
            tot=np.array([h2_2]),
            tot_se=np.array([0.03]),
            intercept=np.array([1.02]),
            intercept_se=0.02,
            mean_chisq=np.array([1.30]),
            lambda_gc=np.array([1.06]),
            ratio=0.06,
            ratio_se=0.02,
        )
        gencov = mock.Mock(
            tot=np.array([0.10]),
            tot_se=np.array([0.01]),
            intercept=np.array([0.001]),
            intercept_se=0.0001,
        )
        return mock.Mock(
            rg_ratio=rg,
            rg_se=rg_se,
            z=z,
            p=p,
            hsq1=hsq1,
            hsq2=hsq2,
            gencov=gencov,
            # As the kernel sets it: True when a heritability estimate is at most 0, else None.
            _negative_hsq=True if h2_2 <= 0 else None,
        )

    def make_annotation_bundle(self):
        metadata = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "POS": [10, 20, 30],
                "SNP": ["rs1", "rs2", "rs3"],
                "CM": [0.1, 0.2, 0.3],
            }
        )
        return make_annotation_bundle(
            metadata=metadata,
            baseline_annotations=pd.DataFrame({"base": [1, 1, 1]}),
            query_annotations=pd.DataFrame({"query1": [1, 0, 0], "query2": [0, 1, 0]}),
            baseline_columns=["base"],
            query_columns=["query1", "query2"],
            chromosomes=["1"],
            source_summary={},
        )

    def write_ldscore_dir(self, root: Path, *, include_query: bool = True) -> Path:
        return write_one_snp_ldscore_dir(root, include_query=include_query)

    def allele_free_table(self, trait_name):
        table = self.make_sumstats_table()
        return replace(table, has_alleles=False, trait_name=trait_name, source_path=f"{trait_name}.parquet",
                       data=table.data.drop(columns=["A1", "A2"]))

    def _high_chisq_sumstats(self):
        # Z=[12,10,2] -> chi^2=[144,100,4]; N.max()=100000 -> default cap=100.
        return self.make_sumstats_table_from_frame(
            pd.DataFrame(
                {
                    "SNP": ["rs1", "rs2", "rs3"],
                    "Z": [12.0, 10.0, 2.0],
                    "N": [100000.0, 100000.0, 100000.0],
                    "A1": ["A", "C", "G"],
                    "A2": ["G", "T", "A"],
                }
            )
        )
