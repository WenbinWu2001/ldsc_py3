"""Whole-run accounting through the same munging path used by scientists."""

import bz2
import gzip

import numpy as np
import pandas as pd
import pytest

from ldsc import GlobalConfig, LDSCInputError, MungeConfig, SumstatsMunger, load_sumstats
from ldsc._sumstats import munge_input
from ldsc._sumstats import munge_kernel as kernel_munge
from ldsc._genome.liftover import LiftOverMappingResult


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
    marker = tmp_path / "out/RUN_FAILED.txt"
    marker.parent.mkdir()
    marker.write_text("earlier failed overwrite")
    munger = SumstatsMunger()
    table = munger.run(
        MungeConfig(raw_sumstats_file=raw, output_dir=tmp_path / "out", sumstats_snps_file=keep, chunk_size=chunk_size),
        global_config=GlobalConfig(snp_identifier="rsid"),
    )
    summary = munger.build_run_summary()
    assert not marker.exists()
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
    config = MungeConfig(raw_sumstats_file=raw, sumstats_snps_file=keep)
    plan = munge_input.resolve_munge_plan(config, munge_input.resolve_source_config(config, GlobalConfig(snp_identifier="rsid")))
    request = munge_input.prepare_munge_input(plan)
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
    from ldsc._genome.build_inference import load_packaged_reference_table

    reference = load_packaged_reference_table().head(300)
    rows = [f"{row.CHR} {int(row.hg19_POS) - 1} evidence{idx} 0 0 1000\n"
            for idx, row in enumerate(reference.itertuples(index=False))]
    rows += [f"1 {900_000_000 + idx} tail{idx} 0 0 1000\n" for idx in range(5000)]
    rows += ["1 950000000 keep1 .05 -.05 1000\n", "1 950000100 keep2 .05 .05 1000\n"]
    raw = tmp_path / "raw.tsv"
    raw.write_text("CHR POS SNP P BETA N\n" + "".join(rows))
    evidence_sizes = []
    from ldsc._genome import build_policy as policy
    resolve = policy.infer_chr_pos_build

    def observe_evidence(frame, **kwargs):
        assert frame.columns.tolist() == ["CHR", "POS"]
        evidence_sizes.append(len(frame))
        return resolve(frame, **kwargs)

    monkeypatch.setattr(policy, "infer_chr_pos_build", observe_evidence)
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


def test_build_inference_reads_the_raw_file_with_the_munging_read_rules(tmp_path):
    # Audit D11 (i): a tab in the header line made the build check read the file as tab-delimited, so it found no
    # coordinate rows in the space-delimited body that munging reads; the run stopped on an unresolved build.
    from ldsc._genome.build_inference import load_packaged_reference_table

    reference = load_packaged_reference_table()
    positions = reference.loc[reference.CHR.astype(str).eq("22"), "hg19_POS"].head(1500).astype(int)
    body = "".join(f"rs{idx} 22 {pos} A G {(-1) ** idx * 2.0} 0.05 1000\n" for idx, pos in enumerate(positions))
    tables = {}
    for name, header in (("spaces", "SNP CHR POS A1 A2 Z P N\n"), ("tab_header", "SNP\tCHR\tPOS\tA1\tA2\tZ\tP\tN\n")):
        raw = tmp_path / f"{name}.txt"
        raw.write_text(header + body)
        tables[name] = SumstatsMunger().run(
            MungeConfig(raw_sumstats_file=raw, output_dir=tmp_path / name, no_snp_restriction=True,
                        output_genome_build="hg19"),
            global_config=GlobalConfig(snp_identifier="chr_pos"),
        )
    assert tables["tab_header"].config_snapshot.genome_build == "hg19"
    assert tables["tab_header"].provenance["coordinate_provenance"]["genome_build_inferred"] is True
    pd.testing.assert_frame_equal(tables["tab_header"].data, tables["spaces"].data)


def test_coordinate_liftover_and_identity_counts_are_exclusive(tmp_path, monkeypatch):
    from ldsc._sumstats import liftover

    raw = tmp_path / "raw.tsv"
    raw.write_text(
        "CHR POS SNP A1 A2 P BETA N INFO\n"
        ". . missing A G .05 0 1000 1\n"
        "1 bad invalid A G .05 0 1000 1\n"
        "1 50 info A G 0 0 1000 .1\n"
        "1 100 duplicate1 A G .05 0 1000 1\n"
        "1 100 duplicate2 A G .05 0 1000 1\n"
        "1 200 unmapped A G .05 0 1000 1\n"
        "1 300 identical G G .05 0 1000 1\n"
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
        "source_duplicate": 2, "unmapped_liftover": 1, "invalid_allele": 1,
    }


_UNSORTED_ROWS = (
    "SNP CHR POS A1 A2 Z P N\n"
    "rs4 chr2 300 A G 2.0 0.05 1000\n"
    "rs1 chr1 200 A G -2.0 0.05 1000\n"
    "rs3 chr2 100 A G 2.0 0.05 1000\n"
    "rs2 chr1 100 A G -2.0 0.05 1000\n"
    "rs6 chrX 50 A G 0.5 0.05 1000\n"
    "rs5 chr1 150 A G -0.5 0.05 1000\n"
)


@pytest.mark.parametrize("mode", ["rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"])
def test_returned_table_and_both_files_share_one_row_order_and_labels(tmp_path, mode):
    # Audit D8: only the Parquet writer normalized CHR and sorted rows, so the returned table and the .sumstats.gz
    # kept the input order and, in rsID-family modes, the 'chr' prefixes.
    raw = tmp_path / "raw.tsv"
    raw.write_text(_UNSORTED_ROWS)
    builds = {} if mode.startswith("rsid") else {"source_genome_build": "hg19", "output_genome_build": "hg19"}
    table = SumstatsMunger().run(
        MungeConfig(raw_sumstats_file=raw, output_dir=tmp_path / "out", no_snp_restriction=True, output_format="both",
                    **builds),
        global_config=GlobalConfig(snp_identifier=mode),
    )
    # Chromosome order 1 < 2 < X, then position.
    snps, chroms = ["rs2", "rs5", "rs1", "rs3", "rs4", "rs6"], ["1", "1", "1", "2", "2", "X"]
    text = pd.read_csv(tmp_path / "out" / "sumstats.gz", sep="\t", dtype={"CHR": str})
    parquet = load_sumstats(tmp_path / "out" / "sumstats.parquet").data
    for frame in (table.data, text, parquet):
        assert frame.SNP.tolist() == snps
        assert frame.CHR.astype(str).tolist() == chroms
        assert frame.POS.tolist() == [100, 150, 200, 100, 300, 50]
    assert table.data.columns.tolist() == text.columns.tolist() == ["SNP", "CHR", "POS", "A1", "A2", "Z", "N"]


@pytest.mark.parametrize("output_format", ["tsv.gz", "parquet"])
@pytest.mark.parametrize("mode", ["rsid", "chr_pos"])
def test_unsupported_contig_rows_are_dropped_before_any_output_in_every_mode(tmp_path, mode, output_format):
    # Audit D8: in rsID mode an unsupported CHR label passed QC and the .sumstats.gz output but stopped the Parquet
    # writer after all QC. rsID-family modes now drop such rows as coordinate-family modes do, whatever the format.
    raw = tmp_path / "raw.tsv"
    raw.write_text("SNP CHR POS A1 A2 Z P N\n"
                   "rs0 chrUn_gl000220 100 A G -0.01 0.05 1000\n"
                   "rs1 1 100 A G -2.0 0.05 1000\n"
                   "rs2 1 200 A G 0.02 0.05 1000\n"
                   "rs3 2 100 A G 2.0 0.05 1000\n")
    builds = {} if mode == "rsid" else {"source_genome_build": "hg19", "output_genome_build": "hg19"}
    munger = SumstatsMunger()
    table = munger.run(
        MungeConfig(raw_sumstats_file=raw, output_dir=tmp_path / "out", no_snp_restriction=True,
                    output_format=output_format, **builds),
        global_config=GlobalConfig(snp_identifier=mode),
    )
    summary = munger.build_run_summary()
    assert table.data.SNP.tolist() == ["rs1", "rs2", "rs3"]
    assert summary.drop_counts["coordinates"] == 1
    assert sum(summary.drop_counts.values()) == 1
    provenance = table.provenance["coordinate_provenance"]
    # The coordinate stage saw 4 rows and kept 3, in both mode families.
    assert provenance["n_invalid_chr_pos"] == provenance["n_dropped_invalid_chr_pos"] == 1
    assert (provenance["n_rows"], provenance["n_retained_after_chr_pos_policy"]) == (4, 3)
    assert [row["SNP"] for row in provenance["coordinate_drop_report"]["examples"]] == ["rs0"]


@pytest.mark.parametrize("output_format", ["parquet", "tsv.gz"])
@pytest.mark.parametrize("prefix", ["", "chr"])
@pytest.mark.parametrize("mode", ["rsid", "rsid_allele_aware"])
def test_rows_without_complete_coordinates_carry_string_chromosome_labels(tmp_path, mode, prefix, output_format):
    # rsID-family rows without a complete CHR/POS pair kept their parsed CHR (a float when CHR has a missing value,
    # or a 'chr'-prefixed string) next to the canonical labels of complete rows: the Parquet writer failed with an
    # internal error, and the .sumstats.gz wrote '1' and '1.0' in one column.
    raw = tmp_path / "raw.tsv"
    raw.write_text("SNP CHR POS A1 A2 Z P N\n"
                   f"rs1 {prefix}1 100 A G -2.0 0.05 1000\n"
                   f"rs2 {prefix}1 NA A G 2.0 0.05 1000\n"
                   "rs3 NA 100 A G -0.01 0.05 1000\n"
                   f"rs4 {prefix}2 100 A G 0.02 0.05 1000\n")
    out = tmp_path / "out"
    table = SumstatsMunger().run(
        MungeConfig(raw_sumstats_file=raw, output_dir=out, no_snp_restriction=True, output_format=output_format),
        global_config=GlobalConfig(snp_identifier=mode),
    )
    if output_format == "parquet":
        written = load_sumstats(out / "sumstats.parquet").data
    else:
        written = pd.read_csv(out / "sumstats.gz", sep="\t", dtype={"CHR": "string"})
    for frame in (table.data, written):
        assert frame.SNP.tolist() == ["rs1", "rs4", "rs2", "rs3"]
        assert frame.CHR.iloc[:3].tolist() == ["1", "2", "1"]
        assert pd.isna(frame.CHR.iloc[3])
    assert all(isinstance(value, str) for value in table.data.CHR.iloc[:3])


@pytest.mark.parametrize("alleles,reason", [(("1", "2"), "invalid_allele"), (("NA", "NA"), "missing_allele")])
@pytest.mark.parametrize("mode", ["rsid_allele_aware", "chr_pos_allele_aware"])
def test_identity_cleanup_that_drops_every_row_stops_before_writing(tmp_path, mode, alleles, reason):
    # A 0-row artifact is a silent analysis change: stop as the parse step does when it retains no row.
    raw = tmp_path / "raw.tsv"
    raw.write_text("SNP CHR POS A1 A2 Z P N\n" + "".join(
        f"rs{idx} 1 {1000 + idx} {alleles[0]} {alleles[1]} {(-1) ** idx * 2.0} 0.05 1000\n" for idx in range(20)))
    builds = {} if mode.startswith("rsid") else {"source_genome_build": "hg19", "output_genome_build": "hg19"}
    out = tmp_path / "out"
    with pytest.raises(LDSCInputError, match="retained 0 of 20 input rows after SNP identity cleanup") as caught:
        SumstatsMunger().run(
            MungeConfig(raw_sumstats_file=raw, output_dir=out, no_snp_restriction=True, output_format="both", **builds),
            global_config=GlobalConfig(snp_identifier=mode),
        )
    assert "identity=20" in str(caught.value) and f"{reason}=20" in str(caught.value)
    assert sorted(path.relative_to(out).as_posix() for path in out.rglob("*") if path.is_file()) == [
        "diagnostics/sumstats.log"]


def test_sample_size_filtering_that_drops_every_row_stops_before_writing(tmp_path):
    raw = tmp_path / "raw.tsv"
    raw.write_text("SNP A1 A2 Z P N\nrs1 A G 2.0 0.05 1000\nrs2 A G -2.0 0.05 1000\n")
    out = tmp_path / "out"
    with pytest.raises(LDSCInputError, match="retained 0 of 2 input rows after sample-size filtering") as caught:
        SumstatsMunger().run(
            MungeConfig(raw_sumstats_file=raw, output_dir=out, no_snp_restriction=True, N_min=5000),
            global_config=GlobalConfig(snp_identifier="rsid"),
        )
    assert "N=2" in str(caught.value)
    assert sorted(path.relative_to(out).as_posix() for path in out.rglob("*") if path.is_file()) == [
        "diagnostics/sumstats.log"]
