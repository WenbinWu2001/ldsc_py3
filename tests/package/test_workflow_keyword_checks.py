"""Keyword workflows reject what the command line would reject, before writing anything."""

import re

import pytest

import ldsc
from ldsc.errors import LDSCUsageError

# Plausible inputs; validation fails before any of them is read.
_INPUTS = {
    "run_ldscore": {"plink_prefix": "panel", "ld_wind_cm": 1.0},
    "run_build_ref_panel": {"plink_prefix": "panel.@", "ld_wind_kb": 1.0},
}
_SHARED = "set_global_config(GlobalConfig("
# Each case would otherwise run with a default or a truthy string in place of the intended value.
_INVALID_KEYWORDS = [
    ("run_ldscore", {"ld_wind_cn": 1.0}, "unknown keyword argument(s) ld_wind_cn"),
    ("run_ldscore", {"bfile": "panel"}, "unknown keyword argument(s) bfile"),
    ("run_ldscore", {"regr_snps_exclude_regions": "MHC"},
     "regr_snps_exclude_regions='MHC' is not one of 'none', 'mhc', 'centromeres', 'mhc-and-centromeres'"),
    ("run_ldscore", {"yes_really": "no"}, "yes_really='no' must be False or True"),
    ("run_ldscore", {"overwrite": "False"}, "overwrite='False' must be False or True"),
    ("run_ldscore", {"gene_list_resolution_policy": "resolved_only"},
     "gene_list_resolution_policy='resolved_only' must be 'strict' or 'resolved-only'"),
    ("run_ldscore", {"snp_identifier": "rsid"}, _SHARED + "snp_identifier="),
    ("run_ldscore", {"genome_build": "hg19"}, _SHARED + "genome_build="),
    ("run_ldscore", {"log_level": "DEBUG"}, _SHARED + "log_level="),
    ("run_build_ref_panel", {"maf": 0.01}, "unknown keyword argument(s) maf"),
    ("run_build_ref_panel", {"source_genome_build": "HG19"}, "source_genome_build='HG19' is not one of"),
    ("run_build_ref_panel", {"overwrite": "False"}, "overwrite='False' must be False or True"),
    ("run_build_ref_panel", {"snp_identifier": "rsid"}, _SHARED + "snp_identifier="),
    ("run_build_ref_panel", {"log_level": "DEBUG"}, _SHARED + "log_level="),
    ("run_build_ref_panel", {"genome_build": "hg19"}, "pass source_genome_build"),
]


@pytest.mark.parametrize(("function", "extra", "message"), _INVALID_KEYWORDS)
def test_keyword_workflows_reject_invalid_keywords_before_writing(tmp_path, function, extra, message):
    output_dir = tmp_path / "out"
    with pytest.raises(LDSCUsageError, match=re.escape(f"{function}() cannot run: ") + ".*" + re.escape(message)):
        getattr(ldsc, function)(**_INPUTS[function], output_dir=str(output_dir), **extra)
    assert not output_dir.exists()


@pytest.mark.parametrize(("function", "required"), [
    ("run_ldscore", ["output_dir"]),
    ("run_build_ref_panel", ["output_dir", "plink_prefix"]),
])
def test_keyword_workflows_list_every_missing_required_keyword(function, required):
    with pytest.raises(LDSCUsageError, match=re.escape("missing required keyword argument(s)")) as raised:
        getattr(ldsc, function)()
    listed = str(raised.value).split("missing required keyword argument(s) ")[1].split(";")[0].rstrip(".").split(", ")
    assert sorted(listed) == required
