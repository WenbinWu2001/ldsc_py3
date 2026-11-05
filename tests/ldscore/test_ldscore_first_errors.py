"""The first error each class of invalid ``ldscore`` input raises, on both routes.

Each case has one defect and is run through ``ldsc ldscore`` (``ldsc.cli.main``)
or ``ldsc.run_ldscore``. The error type and the complete message are pinned, so
a change to where a check runs, which check runs first, or its wording shows
up here. Inputs a check does not need are absent paths: a check that started
reading files earlier would raise a different error. The keyword route cannot
express a missing build, because a registered coordinate-family
``GlobalConfig`` always carries one.
"""

from pathlib import Path

import pytest

import ldsc
from ldsc import GlobalConfig, cli, reset_global_config, set_global_config

PLINK = str(Path(__file__).parents[1] / "fixtures/minimal_external_resources/plink/hm3_chr22_subset")
OUT = "{tmp}/out"
PLINK_RUN = ["--genome-build", "hg38", "--plink-prefix", PLINK, "--ld-wind-kb", "1000"]
INDEX_RUN = ["--gene-ldscore-index-dir", "{tmp}/index", "--query-annot-gene-list-sources", "{tmp}/genes.txt"]

BACKEND = (
    "ldscore could not choose a reference-panel backend. Most likely both R2 matrix input and PLINK input were "
    "supplied, or neither was supplied. Pass exactly one of `--r2-matrix-dir` for R2 matrix mode or "
    "`--plink-prefix` for PLINK mode."
)
KEEP_WITH_R2 = (
    "ldscore cannot apply `--keep-indivs-file` in R2 matrix mode. Most likely `--keep-indivs-file` was combined "
    "with `--r2-matrix-dir`, but individual-level filtering only exists before PLINK genotype LD calculation. "
    "Remove `--keep-indivs-file`, or rerun in PLINK mode with `--plink-prefix`."
)
R2_BUILD_MISSING = (
    "ldscore cannot run R2 matrix mode with a coordinate-family SNP identifier mode without a genome build. "
    "Most likely `--snp-identifier chr_pos` or `--snp-identifier chr_pos_allele_aware` was supplied without "
    "`--genome-build`. Pass `--genome-build hg19`, `--genome-build hg38`, or `--genome-build auto`."
)
BUILD_MISSING = (
    "ldscore cannot resolve a coordinate-family SNP identifier mode without a genome build. Most likely "
    "`--snp-identifier` was set to a coordinate-family mode but `--genome-build` was omitted. Pass "
    "`--genome-build auto`, `--genome-build hg19`, or `--genome-build hg38`."
)
RSID_AUTO = (
    "Could not construct GlobalConfig: genome_build='auto' is not valid for rsid-family snp_identifier='rsid'. "
    "Most likely an auto-build coordinate setting was reused with an rsID-family mode. Omit genome_build or "
    "pass genome_build=None for rsID-family workflows."
)
MAF = (
    "ldscore received invalid `--maf-min={}`. Most likely the minor-allele frequency threshold was entered "
    "outside the valid [0, 0.5] range. Pass a value between 0 and 0.5, or omit the option."
)
COMMON_MAF = (
    "ldscore received invalid `--common-maf-min={}`. Most likely the common-SNP MAF threshold was entered "
    "outside the valid [0, 0.5] range. Pass a value between 0 and 0.5."
)
SNP_BATCH = (
    "ldscore received invalid `--snp-batch-size=0`. Most likely the parquet query batch size was set to zero "
    "or a negative value. Pass a positive integer batch size."
)
QUERY_BATCH = (
    "Cannot configure LD-score query batches: query_batch_size=0 must be a positive integer. Set "
    "--query-batch-size (query_batch_size in Python) to 1 or a larger integer."
)
WORKERS = "workers must be a nonzero integer (1, positive workers, or a negative CPU offset)."
NO_WINDOW = (
    "Cannot configure the LD window: set exactly one of --ld-wind-snps, --ld-wind-kb or --ld-wind-cm "
    "(ld_wind_snps, ld_wind_kb or ld_wind_cm in Python); none was set."
)
TWO_WINDOWS = (
    "Cannot configure the LD window: set exactly one of --ld-wind-snps, --ld-wind-kb or --ld-wind-cm "
    "(ld_wind_snps, ld_wind_kb or ld_wind_cm in Python); --ld-wind-snps and --ld-wind-kb were both set."
)
NEGATIVE_WINDOW = (
    "Cannot configure the LD window: --ld-wind-kb (ld_wind_kb in Python) is -5.0. Use a finite number greater "
    "than 0; --ld-wind-snps takes a whole number of SNPs."
)
REGION_BUILD = (
    "Reference/control content validation failed: 1 input issue. regression regions mhc-and-centromeres: "
    "LDSCUsageError: ldscore cannot select named regression exclusion regions without a concrete genome build. "
    "Pass `--genome-build hg19` or `--genome-build hg38`, or use `--regr-snps-exclude-regions none`. "
    "Repair: Correct the listed path, required companions, or schema; preserve the requested scope. "
    "Complete diagnostics: <tmp>/out/diagnostics/input_issues.tsv. "
    "Other causes & fixes: docs/user-guide/troubleshooting.md#staged-input-validation"
)
LIVE_OVERRIDES = (
    "ldscore indexed mode inherits immutable scientific inputs, SNP identity, and genome build from the index "
    "and cannot accept live overrides: {}. Remove these options, or remove --gene-ldscore-index-dir to run "
    "direct mode."
)
INDEX_WITHOUT_GENES = "ldscore indexed mode requires --query-annot-gene-list-sources."
GENE_OPTIONS_WITHOUT_GENES = (
    "Gene coordinate, exclusion, control-list, and resolution policy options are valid only with gene-list "
    "query annotations. Remove the gene-specific option or use gene-list mode."
)
GENES_WITHOUT_PADDING = (
    "Live gene-list mode requires an explicit --padding-bp value; it must be chosen deliberately. Use "
    "--padding-bp 0 for gene bodies or a positive value for padded intervals."
)
GENES_WITHOUT_CATALOG = (
    "Live gene-list mode requires --gene-coordinate-file. Supply a one-build, one-based gene-coordinate "
    "TSV/TSV.GZ; there is no packaged catalog fallback."
)
PADDING_WITHOUT_INTERVALS = (
    "ldscore accepts `--padding-bp` only with live `--query-annot-bed-sources` or "
    "`--query-annot-gene-list-sources`. It cannot be used with prebuilt annotation queries, without a query, or "
    "in indexed mode. Remove `--padding-bp` from the command (or remove the Python `padding_bp` keyword)."
)
QUERY_WITHOUT_BASELINE = (
    "ldscore cannot run query annotations without baseline annotations. Most likely a prebuilt, BED, or "
    "gene-list query source was supplied without `--baseline-annot-sources`. Pass matching baseline "
    "annotations, or create an explicit all-ones `base` baseline annotation over the same SNP rows as the query "
    "annotations before estimating partitioned heritability."
)

CLI_CASES = {
    # backend choice
    "neither-backend": (["--genome-build", "hg38", "--ld-wind-kb", "1000"], "LDSCUsageError", BACKEND),
    "both-backends": ([*PLINK_RUN, "--r2-matrix-dir", "{tmp}/r2"], "LDSCUsageError", BACKEND),
    # --keep-indivs-file with an R2 matrix
    "keep-with-r2": (["--genome-build", "hg38", "--r2-matrix-dir", "{tmp}/r2", "--keep-indivs-file",
                      "{tmp}/keep.txt", "--ld-wind-kb", "1000"], "LDSCUsageError", KEEP_WITH_R2),
    # build missing
    "r2-build-missing": (["--snp-identifier", "chr_pos", "--r2-matrix-dir", "{tmp}/r2", "--ld-wind-kb", "1000"],
                         "LDSCUsageError", R2_BUILD_MISSING),
    "plink-build-missing": (["--snp-identifier", "chr_pos", "--plink-prefix", PLINK, "--ld-wind-kb", "1000"],
                            "LDSCUsageError", BUILD_MISSING),
    "default-mode-build-missing": (["--plink-prefix", PLINK, "--ld-wind-kb", "1000"], "LDSCUsageError",
                                   BUILD_MISSING),
    "rsid-auto-build": (["--snp-identifier", "rsid", "--genome-build", "auto", "--plink-prefix", PLINK,
                         "--ld-wind-kb", "1000", "--regr-snps-exclude-regions", "none"],
                        "LDSCConfigError", RSID_AUTO),
    # range values
    "maf-min": ([*PLINK_RUN, "--maf-min", "0.7"], "LDSCConfigError", MAF.format(0.7)),
    "common-maf-min": ([*PLINK_RUN, "--common-maf-min", "-0.1"], "LDSCConfigError", COMMON_MAF.format(-0.1)),
    "snp-batch-size": ([*PLINK_RUN, "--snp-batch-size", "0"], "LDSCConfigError", SNP_BATCH),
    "query-batch-size": ([*PLINK_RUN, "--query-batch-size", "0"], "LDSCConfigError", QUERY_BATCH),
    "no-window": (["--genome-build", "hg38", "--plink-prefix", PLINK], "LDSCConfigError", NO_WINDOW),
    # region build
    "region-build": (["--snp-identifier", "rsid", "--plink-prefix", PLINK, "--ld-wind-kb", "1000"],
                     "LDSCInputError", REGION_BUILD),
    # indexed-mode live overrides
    "indexed-genome-build": ([*INDEX_RUN, "--genome-build", "hg19"], "LDSCInputError",
                             LIVE_OVERRIDES.format("--genome-build")),
    "indexed-many-overrides": ([*INDEX_RUN, "--snp-identifier", "chr_pos_allele_aware", "--regr-snps-exclude-regions",
                                "none", "--export-ref-metadata", "--snp-batch-size", "10", "--maf-min", "0.1"],
                               "LDSCInputError", LIVE_OVERRIDES.format(
                                   "--export-ref-metadata, --maf-min, --regr-snps-exclude-regions, "
                                   "--snp-batch-size, --snp-identifier")),
    "indexed-without-genes": (["--gene-ldscore-index-dir", "{tmp}/index"], "LDSCInputError", INDEX_WITHOUT_GENES),
    # gene-list mode option conflicts
    "catalog-without-genes": ([*PLINK_RUN, "--gene-coordinate-file", "{tmp}/genes.tsv"], "LDSCUsageError",
                              GENE_OPTIONS_WITHOUT_GENES),
    "allow-unresolved-without-genes": ([*PLINK_RUN, "--allow-unresolved-genes"], "LDSCUsageError",
                                       GENE_OPTIONS_WITHOUT_GENES),
    "genes-without-padding": ([*PLINK_RUN, "--baseline-annot-sources", "{tmp}/base.annot",
                               "--query-annot-gene-list-sources", "{tmp}/genes.txt", "--gene-coordinate-file",
                               "{tmp}/catalog.tsv"], "LDSCUsageError", GENES_WITHOUT_PADDING),
    "genes-without-catalog": ([*PLINK_RUN, "--baseline-annot-sources", "{tmp}/base.annot",
                               "--query-annot-gene-list-sources", "{tmp}/genes.txt", "--padding-bp", "0"],
                              "LDSCUsageError", GENES_WITHOUT_CATALOG),
    "padding-without-intervals": ([*PLINK_RUN, "--padding-bp", "0"], "LDSCUsageError", PADDING_WITHOUT_INTERVALS),
    "query-without-baseline": ([*PLINK_RUN, "--query-annot-bed-sources", "{tmp}/query.bed"], "LDSCUsageError",
                               QUERY_WITHOUT_BASELINE),
}

PLINK_KWARGS = {"plink_prefix": PLINK, "ld_wind_kb": 1000}
INDEX_KWARGS = {"gene_ldscore_index_dir": "{tmp}/index", "query_annot_gene_list_sources": "{tmp}/genes.txt"}
RSID = GlobalConfig(snp_identifier="rsid")

KWARGS_CASES = {
    # backend choice
    "neither-backend": (None, {"ld_wind_kb": 1000}, "LDSCUsageError", BACKEND),
    # --keep-indivs-file with an R2 matrix
    "keep-with-r2": (None, {"r2_matrix_dir": "{tmp}/r2", "keep_indivs_file": "{tmp}/keep.txt", "ld_wind_kb": 1000},
                     "LDSCUsageError", KEEP_WITH_R2),
    # range values
    "maf-min": (None, {**PLINK_KWARGS, "maf_min": 0.7}, "LDSCConfigError", MAF.format(0.7)),
    "common-maf-min": (None, {**PLINK_KWARGS, "common_maf_min": 0.6}, "LDSCConfigError", COMMON_MAF.format(0.6)),
    "snp-batch-size": (None, {**PLINK_KWARGS, "snp_batch_size": 0}, "LDSCConfigError", SNP_BATCH),
    "query-batch-size": (None, {**PLINK_KWARGS, "query_batch_size": 0}, "LDSCConfigError", QUERY_BATCH),
    "workers": (None, {**PLINK_KWARGS, "workers": 0}, "LDSCConfigError", WORKERS),
    "two-windows": (None, {**PLINK_KWARGS, "ld_wind_snps": 5}, "LDSCConfigError", TWO_WINDOWS),
    "negative-window": (None, {"plink_prefix": PLINK, "ld_wind_kb": -5.0}, "LDSCConfigError", NEGATIVE_WINDOW),
    # region build: a registered rsID-family configuration carries no build
    "region-build": (RSID, PLINK_KWARGS, "LDSCInputError", REGION_BUILD),
    # indexed-mode live overrides; the registered identity is not one
    "indexed-overrides": (None, {**INDEX_KWARGS, "maf_min": 0.1, "common_maf_min": 0.05}, "LDSCInputError",
                          LIVE_OVERRIDES.format("--common-maf-min, --maf-min")),
    "indexed-registered-identity": (RSID, {"gene_ldscore_index_dir": "{tmp}/index"}, "LDSCInputError",
                                    INDEX_WITHOUT_GENES),
    "indexed-workers": (None, {**INDEX_KWARGS, "workers": 0}, "LDSCConfigError", WORKERS),
    # gene-list mode option conflicts
    "policy-without-genes": (None, {**PLINK_KWARGS, "gene_list_resolution_policy": "resolved-only"},
                             "LDSCUsageError", GENE_OPTIONS_WITHOUT_GENES),
    "exclusion-without-genes": (None, {**PLINK_KWARGS, "gene_exclude_regions": "mhc"}, "LDSCUsageError",
                                GENE_OPTIONS_WITHOUT_GENES),
    "genes-without-catalog": (None, {**PLINK_KWARGS, "baseline_annot_sources": "{tmp}/base.annot",
                                     "query_annot_gene_list_sources": "{tmp}/genes.txt", "padding_bp": 0},
                              "LDSCUsageError", GENES_WITHOUT_CATALOG),
    "query-without-baseline": (None, {**PLINK_KWARGS, "query_annot_sources": "{tmp}/query.annot"}, "LDSCUsageError",
                               QUERY_WITHOUT_BASELINE),
}


@pytest.fixture(autouse=True)
def _registry():
    reset_global_config()
    yield
    reset_global_config()


def _fill(value, tmp_path):
    return value.replace("{tmp}", str(tmp_path)) if isinstance(value, str) else value


def _first_error(tmp_path, run):
    with pytest.raises(Exception) as raised:
        run()
    return type(raised.value).__name__, str(raised.value).replace(str(tmp_path), "<tmp>")


@pytest.mark.parametrize(("argv", "error_type", "message"), CLI_CASES.values(), ids=CLI_CASES)
def test_command_line_first_error(tmp_path, argv, error_type, message):
    argv = ["ldscore", *(_fill(token, tmp_path) for token in argv), "--output-dir", _fill(OUT, tmp_path)]
    assert _first_error(tmp_path, lambda: cli.main(argv)) == (error_type, message)


@pytest.mark.parametrize(("registered", "kwargs", "error_type", "message"), KWARGS_CASES.values(), ids=KWARGS_CASES)
def test_keyword_first_error(tmp_path, registered, kwargs, error_type, message):
    if registered is not None:
        set_global_config(registered)
    kwargs = {name: _fill(value, tmp_path) for name, value in {**kwargs, "output_dir": OUT}.items()}
    assert _first_error(tmp_path, lambda: ldsc.run_ldscore(**kwargs)) == (error_type, message)
