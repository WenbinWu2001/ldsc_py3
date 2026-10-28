"""R2 parquet metadata contract shared by the ldscore and query-r2 readers.

A reference-panel R2 parquet is accepted only when it matches what the current
``ldsc build-r2-panel`` writes: int16 ``R2`` with ``ldsc:r2_encoding`` /
``ldsc:r2_scale``, already-unbiased values with their sample size, and complete
binding, identity, and LD-window metadata. Every defect is rejected before any
pair is decoded, and all defects of all files are reported together.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pa = pytest.importorskip("pyarrow")
pq = pytest.importorskip("pyarrow.parquet")

from ldsc.config import GlobalConfig, RefPanelConfig
from ldsc.errors import LDSCInputError
from ldsc._ref_panel_preflight import inspect_r2_headers
from ldsc._kernel import ref_panel_builder as kb
from ldsc._kernel.ref_panel import ParquetR2RefPanel
from ldsc._core.snp_identity import sidecar_identity_sha256
from ldsc.r2_query import R2Panel

ANCHOR = "docs/user-guide/troubleshooting.md#ldscore-r2-panel-is-incompatible"


def _write_panel(root: Path, chrom: str = "1", *, snp_identifier: str = "chr_pos_allele_aware") -> Path:
    """Write one package-format chromosome (sidecar + R2 parquet) with the real writers."""
    panel = pd.DataFrame({
        "CHR": [chrom] * 3, "POS": [100, 120, 140], "SNP": ["rs1", "rs2", "rs3"],
        "A1": ["A", "A", "C"], "A2": ["C", "G", "T"], "CM": [0.0, 0.0, 0.0], "MAF": [0.3, 0.2, 0.4],
    })
    meta = root / f"chr{chrom}_meta.tsv.gz"
    kb.write_runtime_metadata_sidecar(panel, meta, genome_build="hg19", snp_identifier=snp_identifier)
    r2 = meta.with_name(f"chr{chrom}_r2.parquet")
    pairs = [(np.array([0, 0]), np.array([1, 2]), np.array([1.0, 0.25], dtype=np.float32),
              np.array([1, -1], dtype=np.int8))]
    kb.write_r2_parquet(
        pair_chunks=pairs, path=r2, genome_build="hg19", n_samples=100, snp_identifier=snp_identifier,
        n_snps=3, sidecar_identity_sha256=sidecar_identity_sha256(panel),
        ld_window_mode="kb", ld_window_value=1000.0,
    )
    return r2


def _rewrite(path: Path, *, drop=(), set_meta=None, float_r2=False) -> None:
    """Alter one written R2 parquet the way a stripped or foreign artifact would differ."""
    table = pq.read_table(path)
    metadata = dict(table.schema.metadata)
    for key in drop:
        metadata.pop(key)
    metadata.update(set_meta or {})
    if float_r2:
        index = table.schema.get_field_index("R2")
        decoded = pa.array(table.column("R2").to_numpy().astype(np.float32) / np.float32(32767))
        table = table.set_column(index, "R2", decoded)
    pq.write_table(table.replace_schema_metadata(metadata), path)


# Each defect names the metadata key (or column) the error must identify.
DEFECTS = {
    "missing scale": (dict(drop=[b"ldsc:r2_scale"]), "ldsc:r2_scale"),
    "non-writer scale": (dict(set_meta={b"ldsc:r2_scale": b"16383"}), "ldsc:r2_scale"),
    "missing encoding": (dict(drop=[b"ldsc:r2_encoding"]), "ldsc:r2_encoding"),
    "wrong encoding": (dict(set_meta={b"ldsc:r2_encoding": b"float32"}), "ldsc:r2_encoding"),
    "float R2": (dict(float_r2=True), "R2 column"),
    "bias missing, n kept": (dict(drop=[b"ldsc:r2_bias"]), "ldsc:r2_bias"),
    "bias and n missing": (dict(drop=[b"ldsc:r2_bias", b"ldsc:n_samples"]), "ldsc:r2_bias"),
    "raw bias": (dict(set_meta={b"ldsc:r2_bias": b"raw"}), "ldsc:r2_bias"),
    "missing n": (dict(drop=[b"ldsc:n_samples"]), "ldsc:n_samples"),
    "nonpositive n": (dict(set_meta={b"ldsc:n_samples": b"0"}), "ldsc:n_samples"),
    # The unbiased correction divides by n - 2, so fewer than 3 samples cannot be decoded.
    "two samples": (dict(set_meta={b"ldsc:n_samples": b"2"}), "ldsc:n_samples"),
    "missing window": (dict(drop=[b"ldsc:ld_window_mode", b"ldsc:ld_window_value"]), "ldsc:ld_window_mode"),
    "missing pair floor": (dict(drop=[b"ldsc:min_r2"]), "ldsc:min_r2"),
    "non-numeric pair floor": (dict(set_meta={b"ldsc:min_r2": b"none"}), "ldsc:min_r2"),
    "pair floor above one": (dict(set_meta={b"ldsc:min_r2": b"1.5"}), "ldsc:min_r2"),
    "missing binding": (dict(drop=[b"ldsc:sidecar_identity_sha256"]), "ldsc:sidecar_identity_sha256"),
    "missing identity": (dict(drop=[b"ldsc:artifact_type"]), "ldsc:artifact_type"),
}


def _assert_contract_error(error: Exception, path: Path, *keys: str) -> None:
    message = str(error)
    assert str(path) in message
    for key in keys:
        assert key in message
    assert "ldsc build-r2-panel" in message


@pytest.mark.parametrize("defect", sorted(DEFECTS))
def test_ldscore_reader_rejects_contract_defect(tmp_path, defect):
    rewrite, key = DEFECTS[defect]
    r2 = _write_panel(tmp_path / "hg19")
    _rewrite(r2, **rewrite)
    panel = ParquetR2RefPanel(
        GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg19"),
        RefPanelConfig(backend="parquet_r2", r2_dir=str(tmp_path / "hg19")),
    )
    with pytest.raises(LDSCInputError) as failure:
        panel.build_reader("1").close()
    _assert_contract_error(failure.value, r2, key)


@pytest.mark.parametrize("defect", sorted(DEFECTS))
def test_query_r2_panel_rejects_contract_defect_at_open(tmp_path, defect):
    rewrite, key = DEFECTS[defect]
    r2 = _write_panel(tmp_path / "hg19")
    _rewrite(r2, **rewrite)
    with pytest.raises(LDSCInputError) as failure:
        R2Panel.open(tmp_path / "hg19", genome_build="hg19")
    _assert_contract_error(failure.value, r2, key)


def test_one_error_lists_every_defect_of_every_file(tmp_path):
    root = tmp_path / "hg19"
    first, second = _write_panel(root, "1"), _write_panel(root, "2")
    _rewrite(first, drop=[b"ldsc:r2_scale"], set_meta={b"ldsc:r2_bias": b"raw"})
    _rewrite(second, drop=[b"ldsc:n_samples"])
    with pytest.raises(LDSCInputError) as failure:
        inspect_r2_headers(root, "hg19")
    _assert_contract_error(failure.value, first, "ldsc:r2_scale", "ldsc:r2_bias")
    _assert_contract_error(failure.value, second, "ldsc:n_samples")
    with pytest.raises(LDSCInputError) as failure:
        R2Panel.open(root, genome_build="hg19")
    _assert_contract_error(failure.value, first, "ldsc:r2_scale", "ldsc:r2_bias")
    _assert_contract_error(failure.value, second, "ldsc:n_samples")


def test_unreadable_file_is_reported_with_the_other_files_defects(tmp_path):
    root = tmp_path / "hg19"
    first, second = _write_panel(root, "1"), _write_panel(root, "2")
    first.write_bytes(b"not a parquet file")
    _rewrite(second, drop=[b"ldsc:r2_scale"])
    for open_panel in (lambda: inspect_r2_headers(root, "hg19"), lambda: R2Panel.open(root, genome_build="hg19")):
        with pytest.raises(LDSCInputError) as failure:
            open_panel()
        _assert_contract_error(failure.value, first, "cannot be read as parquet")
        _assert_contract_error(failure.value, second, "ldsc:r2_scale")


def test_query_r2_panel_rejects_chromosomes_with_different_sample_sizes(tmp_path):
    # query-r2 converts unbiased R2 back to r with each panel's sample size, so
    # chromosome files from runs with different samples cannot form one panel.
    root = tmp_path / "hg19"
    _write_panel(root, "1")
    second = _write_panel(root, "2")
    _rewrite(second, set_meta={b"ldsc:n_samples": b"5000"})
    with pytest.raises(LDSCInputError, match=r"different sample sizes \(ldsc:n_samples 100, 5000\)"):
        R2Panel.open(root, genome_build="hg19")


def test_writer_output_satisfies_contract_and_decodes_exactly(tmp_path):
    root = tmp_path / "hg19"
    _write_panel(root)
    inspect_r2_headers(root, "hg19")
    panel = ParquetR2RefPanel(
        GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg19"),
        RefPanelConfig(backend="parquet_r2", r2_dir=str(root)),
    )
    reader = panel.build_reader("1")
    try:
        _, _, values = next(reader.iter_all_pairs())
    finally:
        reader.close()
    # Writer stores round(R2 * 32767); the endpoint 1.0 and 0.25 decode within one step.
    np.testing.assert_allclose(np.sort(values), [0.25, 1.0], rtol=0, atol=1 / 32767)
    assert values.dtype == np.float32
    result = R2Panel.open(root, genome_build="hg19").query_pairs(pd.DataFrame({
        "CHR_1": ["1"], "POS_1": [100], "A1_1": ["A"], "A2_1": ["C"],
        "CHR_2": ["1"], "POS_2": [120], "A1_2": ["A"], "A2_2": ["G"],
    }))
    assert result["r2"].tolist() == [1.0]
    assert np.isfinite(result["r"]).all()


@pytest.mark.parametrize("requested,accepted", [
    ("chr_pos_allele_aware", True), ("chr_pos", True), ("rsid", False), ("rsid_allele_aware", False),
])
def test_query_r2_identifier_override_must_be_panel_mode_or_its_base(tmp_path, requested, accepted):
    root = tmp_path / "hg19"
    _write_panel(root)
    if accepted:
        assert R2Panel.open(root, snp_identifier=requested, genome_build="hg19").snp_identifier == requested
        return
    with pytest.raises(LDSCInputError, match="chr_pos_allele_aware"):
        R2Panel.open(root, snp_identifier=requested, genome_build="hg19")


def test_query_r2_rejects_allele_aware_override_on_base_panel(tmp_path):
    root = tmp_path / "hg19"
    _write_panel(root, snp_identifier="chr_pos")
    with pytest.raises(LDSCInputError, match="chr_pos"):
        R2Panel.open(root, snp_identifier="chr_pos_allele_aware", genome_build="hg19")


def test_python_api_run_reports_every_chromosome_defect_before_computing(tmp_path):
    """``LDScoreCalculator.run`` checks every footer and sidecar before any chromosome starts.

    Without the check, chromosome workers stop at the first defect they reach,
    after computing earlier chromosomes, and report only that one.
    """
    import gzip

    from ldsc import AnnotationBundle, LDScoreCalculator
    from ldsc.outputs import LDScoreOutputConfig
    from tests.test_ldscore_query_batches import prepared_case

    metadata, baseline, queries, config, panel, numerical = prepared_case(tmp_path)
    _rewrite(tmp_path / "panel" / "chr1_r2.parquet", set_meta={b"ldsc:n_samples": b"2"})
    sidecar = tmp_path / "panel" / "chr2_meta.tsv.gz"
    lines = gzip.open(sidecar, "rt").read().splitlines()
    lines = [line.rsplit("\t", 1)[0] + "\t0.7" if "\trs6\t" in line else line for line in lines]
    with gzip.open(sidecar, "wt") as handle:
        handle.write("\n".join(lines) + "\n")

    output = tmp_path / "out"
    with AnnotationBundle.from_frames(metadata, baseline, queries, config_snapshot=config) as bundle:
        with pytest.raises(LDSCInputError) as failure:
            LDScoreCalculator().run(bundle, panel, numerical, config, LDScoreOutputConfig(output))
    message = str(failure.value)
    assert "chr1_r2.parquet" in message and "ldsc:n_samples" in message
    assert "chr2_meta.tsv.gz" in message and "out-of-range MAF" in message
    assert not (output / "diagnostics" / "chromosome-logs").exists()
