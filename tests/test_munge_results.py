"""Whole-run accounting through the same munging path used by scientists."""

import bz2
import gzip

import numpy as np
import pandas as pd
import pytest

from ldsc import GlobalConfig, MungeConfig, SumstatsMunger, load_sumstats
from ldsc import _sumstats_input as munge_input
from ldsc._kernel import sumstats_munger as kernel_munge
from ldsc._kernel.liftover import LiftOverMappingResult, SumstatsLiftoverRequest


@pytest.mark.parametrize("suffix,compress", [(".tsv", str.encode), (".gz", lambda text: gzip.compress(text.encode())), (".bz2", lambda text: bz2.compress(text.encode()))])
@pytest.mark.parametrize("chunk_size", [2, 100])
def test_summary_counts_parsed_rows_and_exclusive_filter_stages(tmp_path, suffix, compress, chunk_size):
    raw = tmp_path / ("raw" + suffix)
    raw.write_bytes(compress(
        "## producer metadata\nSNP P BETA N INFO FRQ\n\n"
        "missing .05 NA 1000 1 .2\n"
        "info 0 0 1000 .1 .2\n"
        "maf 0 0 1000 1 0\n"
        "pvalue 0 0 1000 1 .2\n"
        "outside .05 0 1000 1 .2\n"
        "small .05 0 10 1 .2\n"
        "duplicate .05 0 1000 1 .2\n"
        "keep1 .05 -.05 1000 1 .2\n"
        "duplicate .05 0 1000 1 .2\n"
        "keep2 .05 .05 1000 1 .2\n\n"
    ))
    keep = tmp_path / "keep.tsv"
    keep.write_text("SNP\nsmall\nduplicate\nkeep1\nkeep2\n")
    munger = SumstatsMunger()
    table = munger.run(
        MungeConfig(raw_sumstats_file=raw, output_dir=tmp_path / "out", sumstats_snps_file=keep, chunk_size=chunk_size),
        global_config=GlobalConfig(snp_identifier="rsid"),
    )
    summary = munger.build_run_summary()
    assert table.data.SNP.tolist() == ["keep1", "keep2"]
    np.testing.assert_allclose(table.data.Z, [-1.95996398454, 1.95996398454])
    assert summary.n_input_rows == 10
    assert summary.n_retained_rows == 2
    assert summary.drop_counts == {
        "NA": 1, "coordinates": 0, "INFO": 1, "FRQ": 1, "P": 1,
        "sumstats_snps": 1, "N": 1, "NSTUDY": 0, "liftover": 0, "identity": 2,
    }
    assert sum(summary.drop_counts.values()) == summary.n_input_rows - summary.n_retained_rows
    assert summary.used_n_rule == "input_columns"
    assert load_sumstats(summary.output_paths["sumstats_parquet"]).data.SNP.tolist() == ["keep1", "keep2"]


def test_resolved_request_can_run_twice_without_accumulating_state(tmp_path):
    raw = tmp_path / "raw.tsv"
    raw.write_text("SNP P BETA N\nkeep .05 .05 1000\noutside .05 -.05 1000\n")
    keep = tmp_path / "keep.tsv"
    keep.write_text("SNP\nkeep\n")
    request = munge_input.prepare_munge_input(
        str(raw), MungeConfig(no_snp_restriction=True), MungeConfig(no_snp_restriction=True), GlobalConfig(snp_identifier="rsid"),
        SumstatsLiftoverRequest(), str(keep),
    )
    first = kernel_munge.munge_sumstats(request)
    second = kernel_munge.munge_sumstats(request)
    assert first.n_input_rows == second.n_input_rows == 2
    assert first.drop_counts == second.drop_counts
    assert first.drop_counts["sumstats_snps"] == 1
    assert first.coordinate_metadata == second.coordinate_metadata
    pd.testing.assert_frame_equal(first.data, second.data)
    assert first.data.SNP.tolist() == ["keep"]
    assert sorted(path.name for path in tmp_path.iterdir()) == ["keep.tsv", "raw.tsv"]


def test_allele_restriction_counts_deferred_bad_alleles_only_at_identity_cleanup(tmp_path):
    raw = tmp_path / "raw.tsv"
    raw.write_text(
        "SNP A1 A2 P BETA N\n"
        "bad . G .05 0 1000\n"
        "keep A G .05 .05 1000\n"
        "outside A G .05 0 1000\n"
        "mismatch A C .05 0 1000\n"
    )
    keep = tmp_path / "keep.tsv"
    keep.write_text("SNP A1 A2\nbad A G\nkeep A G\nmismatch A G\n")
    munger = SumstatsMunger()
    table = munger.run(
        MungeConfig(raw_sumstats_file=raw, output_dir=tmp_path / "out", sumstats_snps_file=keep, chunk_size=2),
        global_config=GlobalConfig(snp_identifier="rsid_allele_aware"),
    )
    summary = munger.build_run_summary()
    assert table.data.SNP.tolist() == ["keep"]
    assert summary.drop_counts["sumstats_snps"] == 2
    assert summary.drop_counts["identity"] == 1
    assert sum(summary.drop_counts.values()) == 3
    dropped = pd.read_csv(summary.output_paths["dropped_snps_tsv_gz"], sep="\t")
    assert dropped.SNP.tolist() == ["bad"]
    assert dropped.reason.tolist() == ["missing_allele"]


@pytest.mark.parametrize("columns,rows,options,rule,n_value,nstudy_drops", [
    ("", ["", ""], {"N": 500}, "fixed_N", 500, 0),
    ("", ["", ""], {"N_cas": 200, "N_con": 300}, "fixed_case_control_N", 500, 0),
    ("N", ["1000", "1000"], {"N": 500}, "input_columns", 1000, 0),
    ("N_CAS N_CON", ["400 600", "400 600"], {}, "input_columns", 1000, 0),
    ("NSTUDY", ["3", "5"], {"N": 500}, "fixed_N", 500, 1),
])
def test_summary_reports_the_sample_size_rule_actually_used(tmp_path, columns, rows, options, rule, n_value, nstudy_drops):
    raw = tmp_path / "raw.tsv"
    raw.write_text(f"SNP P BETA {columns}\nrs1 .05 -.05 {rows[0]}\nrs2 .05 .05 {rows[1]}\n")
    munger = SumstatsMunger()
    table = munger.run(
        MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw, output_dir=tmp_path / "out", **options),
        global_config=GlobalConfig(snp_identifier="rsid"),
    )
    summary = munger.build_run_summary()
    assert summary.used_n_rule == rule
    assert table.data.N.tolist() == [n_value] * (2 - nstudy_drops)
    assert summary.drop_counts["NSTUDY"] == nstudy_drops
    assert sum(summary.drop_counts.values()) == nstudy_drops


def test_production_build_inference_uses_bounded_raw_evidence_before_qc(tmp_path, monkeypatch):
    from ldsc.genome_build_inference import load_packaged_reference_table

    reference = load_packaged_reference_table().head(300)
    rows = [f"{row.CHR} {int(row.hg19_POS) - 1} evidence{idx} 0 0 1000\n"
            for idx, row in enumerate(reference.itertuples(index=False))]
    rows += [f"1 {900_000_000 + idx} tail{idx} 0 0 1000\n" for idx in range(5000)]
    rows += ["1 950000000 keep1 .05 -.05 1000\n", "1 950000100 keep2 .05 .05 1000\n"]
    raw = tmp_path / "raw.tsv"
    raw.write_text("CHR POS SNP P BETA N\n" + "".join(rows))
    evidence_sizes = []
    resolve = munge_input.resolve_chr_pos_table

    def observe_evidence(frame, **kwargs):
        assert frame.columns.tolist() == ["CHR", "POS"]
        evidence_sizes.append(len(frame))
        return resolve(frame, **kwargs)

    monkeypatch.setattr(munge_input, "resolve_chr_pos_table", observe_evidence)
    munger = SumstatsMunger()
    table = munger.run(
        MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw, output_dir=tmp_path / "out", source_genome_build="auto",
                    output_genome_build="hg19", chunk_size=1000),
        global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="auto"),
    )
    assert len(evidence_sizes) == 1
    assert 300 <= evidence_sizes[0] <= 5000
    assert table.config_snapshot.genome_build == "hg19"
    assert table.data.POS.tolist() == [950000001, 950000101]
    assert table.provenance["coordinate_provenance"]["coordinate_basis"] == "0-based"
    summary = munger.build_run_summary()
    assert summary.n_input_rows == 5302
    assert summary.n_retained_rows == 2
    assert summary.drop_counts["P"] == 5300


def test_coordinate_liftover_and_identity_counts_are_exclusive(tmp_path, monkeypatch):
    from ldsc._kernel import liftover

    raw = tmp_path / "raw.tsv"
    raw.write_text(
        "CHR POS SNP A1 A2 P BETA N INFO\n"
        ". . missing A G .05 0 1000 1\n"
        "1 bad invalid A G .05 0 1000 1\n"
        "1 50 info A G 0 0 1000 .1\n"
        "1 100 duplicate1 A G .05 0 1000 1\n"
        "1 100 duplicate2 A G .05 0 1000 1\n"
        "1 200 unmapped A G .05 0 1000 1\n"
        "1 300 ambiguous A T .05 0 1000 1\n"
        "1 400 keep A G .05 .05 1000 1\n"
    )

    class Translator:
        def __init__(self, **kwargs):
            pass

        def map_positions(self, chrom, positions):
            assert chrom == "1"
            np.testing.assert_array_equal(positions, [200, 300, 400])
            return LiftOverMappingResult(
                translated_positions=np.array([1300, 1400]), keep_mask=np.array([False, True, True]),
                unmapped_count=1, cross_chrom_count=0,
            )

    monkeypatch.setattr(liftover, "LiftOverTranslator", Translator)
    chain = tmp_path / "test.chain"
    chain.write_text("translation supplied by the test double\n")
    munger = SumstatsMunger()
    table = munger.run(
        MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw, output_dir=tmp_path / "out", source_genome_build="hg19",
                    output_genome_build="hg38", liftover_chain_file=chain, chunk_size=2),
        global_config=GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg19"),
    )
    summary = munger.build_run_summary()
    assert table.data.SNP.tolist() == ["keep"]
    assert table.data.POS.tolist() == [1400]
    assert table.config_snapshot.genome_build == "hg38"
    assert summary.n_input_rows == 8
    assert summary.n_retained_rows == 1
    assert {key: value for key, value in summary.drop_counts.items() if value} == {
        "coordinates": 2, "INFO": 1, "liftover": 3, "identity": 1,
    }
    provenance = table.provenance["coordinate_provenance"]
    assert provenance["coordinate_drop_report"]["n_dropped"] == 2
    assert provenance["liftover"]["n_dropped"] == 3
    dropped = pd.read_csv(summary.output_paths["dropped_snps_tsv_gz"], sep="\t")
    assert dropped.reason.value_counts().to_dict() == {
        "source_duplicate": 2, "unmapped_liftover": 1, "strand_ambiguous_allele": 1,
    }
