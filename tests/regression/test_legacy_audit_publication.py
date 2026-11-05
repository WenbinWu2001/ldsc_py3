"""The regression writers publish the legacy summary-statistics audit with their results."""

import gzip
import json

import pandas as pd

from ldsc._regression import runner as workflow
from ldsc._regression.outputs import (
    LEGACY_SUMSTATS_DROP_COLUMNS, RG_CONCISE_COLUMNS, RG_FULL_COLUMNS, H2DirectoryWriter, H2OutputConfig,
    RgDirectoryWriter, RgOutputConfig,
)
from ldsc._regression.results import RgResultFamily
from ldsc.config import RegressionConfig
from tests.regression.test_nonfinite_statistics import functional_panel
from tests.regression.test_trait_preparation import legacy
from tests.worker_fixtures import calculation_body

AUDIT = "diagnostics/dropped_snps/legacy_sumstats.tsv.gz"
DROPS = pd.DataFrame([["trait", "trait.sumstats.gz", "rs9", "A", "C", "missing_panel_rsid", 0]],
                     columns=LEGACY_SUMSTATS_DROP_COLUMNS)


def files(root):
    return json.loads((root / "diagnostics/metadata.json").read_text())["files"]


def write_h2(root, drops, *, overwrite=False):
    return H2DirectoryWriter().write(pd.DataFrame({"trait_name": ["trait"]}),
                                     H2OutputConfig(output_dir=root, overwrite=overwrite),
                                     metadata={}, diagnostic_bins=None, legacy_sumstats_drops=drops)


def test_h2_writer_publishes_and_lists_the_legacy_audit(tmp_path):
    paths = write_h2(tmp_path, DROPS)

    assert paths["legacy_sumstats_drops"] == str(tmp_path / AUDIT)
    assert files(tmp_path)["legacy_sumstats_drops"] == AUDIT
    assert pd.read_csv(tmp_path / AUDIT, sep="\t").equals(DROPS)
    # The gzip header names the final file and carries no timestamp, so equal drops give equal bytes.
    header = (tmp_path / AUDIT).read_bytes()[:30]
    assert header[4:8] == b"\x00\x00\x00\x00" and b"legacy_sumstats.tsv\x00" in header
    assert [path.name for path in (tmp_path / AUDIT).parent.iterdir()] == ["legacy_sumstats.tsv.gz"]


def test_an_overwrite_without_legacy_input_removes_the_audit_and_its_entry(tmp_path):
    write_h2(tmp_path, DROPS)

    paths = write_h2(tmp_path, None, overwrite=True)

    assert not (tmp_path / AUDIT).exists()
    assert "legacy_sumstats_drops" not in paths
    assert "legacy_sumstats_drops" not in files(tmp_path)


def test_rg_writer_publishes_the_audit_its_result_carries(tmp_path):
    result = RgResultFamily(rg=pd.DataFrame(columns=RG_CONCISE_COLUMNS),
                            rg_full=pd.DataFrame(columns=RG_FULL_COLUMNS),
                            h2_per_trait=pd.DataFrame(), per_pair_metadata=[], legacy_sumstats_drops=DROPS)

    paths = RgDirectoryWriter().write(result, RgOutputConfig(output_dir=tmp_path))

    assert paths["legacy_sumstats_drops"] == str(tmp_path / AUDIT)
    assert files(tmp_path)["legacy_sumstats_drops"] == AUDIT


def test_a_direct_partitioned_batch_run_writes_the_audit_of_its_legacy_input(tmp_path):
    table = legacy("trait", 5, extra_row={"SNP": "absent", "A1": "A", "A2": "C", "Z": 1.0, "N": 10_000.0})
    runner = workflow.RegressionRunner(regression_config=RegressionConfig(num_blocks=10))

    with calculation_body():
        result = runner.estimate_partitioned_h2_batch(table, functional_panel(), output_dir=tmp_path)

    assert result.output_paths["legacy_sumstats_drops"] == str(tmp_path / AUDIT)
    assert files(tmp_path)["legacy_sumstats_drops"] == AUDIT
    with gzip.open(tmp_path / AUDIT, "rt") as handle:
        audit = pd.read_csv(handle, sep="\t")
    assert audit["SNP"].tolist() == ["absent"]
    pd.testing.assert_frame_equal(audit, result.legacy_sumstats_drops.reset_index(drop=True), check_dtype=False)
