"""Tests for the maintainer tool that rebuilds the packaged HapMap3 map.

Fixtures are small files in the formats of the real inputs: the w_hm3
snplist, dbSNP RefSNP VCFs whose CHROM column holds RefSeq accessions, NCBI
assembly reports mapping those accessions to chromosomes, and the
one-record-per-line refsnp-merged JSON. GRCh37 autosome accessions end in the
versions below and GRCh38 accessions are one version higher, as in the
GCF_000001405.25 and GCF_000001405.40 assembly reports.
"""

from __future__ import annotations

import bz2
import contextlib
import gzip
import io
import json
from pathlib import Path
import runpy
import tempfile
import unittest

import pandas as pd

from ldsc._genome.hm3 import _load_hm3_curated_map_from_path

TOOL = Path(__file__).resolve().parents[2] / "tools" / "hm3" / "build_hm3_curated_map.py"

# RefSeq versions of NC_0000NN for chromosomes 1-22 in GRCh37.p13.
GRCH37_VERSIONS = [10, 11, 11, 11, 9, 11, 13, 10, 11, 10, 9, 11, 10, 8, 9, 9, 10, 9, 9, 10, 8, 10]


def accession(chrom: int, build: str) -> str:
    version = GRCH37_VERSIONS[chrom - 1] + (1 if build == "hg38" else 0)
    return f"NC_{chrom:06d}.{version}"


def write_assembly_report(path: Path, build: str) -> None:
    """Write an assembly report with all autosomes, chrX, and one alt scaffold."""
    header = [
        "Sequence-Name", "Sequence-Role", "Assigned-Molecule", "Assigned-Molecule-Location/Type",
        "GenBank-Accn", "Relationship", "RefSeq-Accn", "Assembly-Unit", "Sequence-Length", "UCSC-style-name",
    ]
    lines = [f"# Assembly name:  {'GRCh37.p13' if build == 'hg19' else 'GRCh38.p14'}", "# " + "\t".join(header)]
    for chrom in range(1, 23):
        lines.append("\t".join([str(chrom), "assembled-molecule", str(chrom), "Chromosome", "CM0", "=",
                                accession(chrom, build), "Primary Assembly", "1000000", f"chr{chrom}"]))
    x_accession = "NC_000023.10" if build == "hg19" else "NC_000023.11"
    lines.append("\t".join(["X", "assembled-molecule", "X", "Chromosome", "CM0", "=", x_accession,
                            "Primary Assembly", "1000000", "chrX"]))
    lines.append("\t".join(["HSCHR6_ALT", "alt-scaffold", "6", "Chromosome", "GL0", "=", "NT_167244.1",
                            "ALT_REF_LOCI_1", "1000", "chr6_alt"]))
    path.write_text("\r\n".join(lines) + "\r\n", encoding="utf-8")


def write_vcf(path: Path, records: list[tuple[str, int, str, str, str]]) -> None:
    """Write a bgzip-like gzip VCF from (CHROM, POS, ID, REF, ALT) records."""
    lines = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for chrom, pos, rsid, ref, alt in records:
        lines.append(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\t.\tRS={rsid[2:]}")
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def write_merged_json(path: Path, merges: dict[str, list[str]]) -> None:
    """Write refsnp-merged records: numeric retired ID -> numeric current IDs."""
    with bz2.open(path, "wt", encoding="utf-8") as handle:
        for old, targets in merges.items():
            record = {
                "refsnp_id": old,
                "last_update_build_id": "157",
                "merged_snapshot_data": {"proxy_time": "2020-01-01T00:00Z", "merged_into": targets},
            }
            handle.write(json.dumps(record) + "\n")


def write_sources(
    tmpdir: Path,
    hm3: list[tuple[str, str, str]],
    hg19: list[tuple[str, int, str, str, str]],
    hg38: list[tuple[str, int, str, str, str]],
    merges: dict[str, list[str]] | None = None,
) -> dict[str, Path]:
    paths = {
        "hm3_snplist": tmpdir / "w_hm3.snplist.gz",
        "dbsnp_vcf_hg19": tmpdir / "GCF_000001405.25.gz",
        "dbsnp_vcf_hg38": tmpdir / "GCF_000001405.40.gz",
        "dbsnp_merged_json": tmpdir / "refsnp-merged.json.bz2",
        "assembly_report_hg19": tmpdir / "grch37_assembly_report.txt",
        "assembly_report_hg38": tmpdir / "grch38_assembly_report.txt",
        "output_path": tmpdir / "out" / "hm3_curated_map.tsv.gz",
        "drop_report_path": tmpdir / "out" / "hm3_curated_map.dropped.tsv.gz",
    }
    pd.DataFrame(hm3, columns=["SNP", "A1", "A2"]).to_csv(paths["hm3_snplist"], sep="\t", index=False)
    write_vcf(paths["dbsnp_vcf_hg19"], hg19)
    write_vcf(paths["dbsnp_vcf_hg38"], hg38)
    write_merged_json(paths["dbsnp_merged_json"], merges or {})
    write_assembly_report(paths["assembly_report_hg19"], "hg19")
    write_assembly_report(paths["assembly_report_hg38"], "hg38")
    return paths


class Hm3CuratedMapBuilderTest(unittest.TestCase):
    def setUp(self):
        self.tool = runpy.run_path(str(TOOL))
        self._tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def test_places_snps_natively_in_both_builds_with_six_columns(self):
        paths = write_sources(
            self.tmpdir,
            hm3=[("rs1", "A", "G"), ("rs2", "C", "T")],
            hg19=[(accession(1, "hg19"), 100, "rs1", "A", "G"), (accession(22, "hg19"), 200, "rs2", "C", "T")],
            hg38=[(accession(1, "hg38"), 1100, "rs1", "A", "G"), (accession(22, "hg38"), 2200, "rs2", "C", "T")],
        )

        result = self.tool["build_hm3_curated_map"](**paths)

        expected = pd.DataFrame(
            {
                "CHR": ["1", "22"],
                "hg19_POS": [100, 200],
                "hg38_POS": [1100, 2200],
                "SNP": ["rs1", "rs2"],
                "A1": ["A", "C"],
                "A2": ["G", "T"],
            }
        )
        pd.testing.assert_frame_equal(result.reset_index(drop=True), expected, check_dtype=False)
        packaged_view = _load_hm3_curated_map_from_path(paths["output_path"], preserve_extra_columns=True)
        self.assertEqual(list(packaged_view.columns), ["CHR", "hg19_POS", "hg38_POS", "SNP", "A1", "A2"])
        pd.testing.assert_frame_equal(packaged_view, expected, check_dtype=False)

    def test_accepts_only_allele_matched_autosomal_placements(self):
        hg38_ok = [(accession(1, "hg38"), 1000 + i, f"rs{i}", "A", "G") for i in (3, 4, 5, 6)]
        paths = write_sources(
            self.tmpdir,
            hm3=[("rs3", "A", "G"), ("rs4", "A", "G"), ("rs5", "A", "G"), ("rs6", "A", "G")],
            hg19=[
                (accession(1, "hg19"), 30, "rs3", "A", "C,G"),  # multi-allelic: A/G is one pair
                (accession(1, "hg19"), 40, "rs4", "T", "C"),  # strand flip of A/G
                (accession(1, "hg19"), 50, "rs5", "A", "C"),  # different alleles
                ("NT_167244.1", 60, "rs6", "A", "G"),  # alt scaffold
                ("NC_000023.10", 60, "rs6", "A", "G"),  # chrX
            ],
            hg38=hg38_ok,
        )

        result = self.tool["build_hm3_curated_map"](**paths)

        self.assertEqual(result["SNP"].tolist(), ["rs3", "rs4"])
        self.assertEqual(result["hg19_POS"].tolist(), [30, 40])
        self.assertEqual(result["A1"].tolist(), ["A", "A"])
        self.assertEqual(result["A2"].tolist(), ["G", "G"])
        drops = read_drop_report(paths)
        self.assertEqual(drops["SNP"].tolist(), ["rs5", "rs6"])
        self.assertEqual(drops["reason"].tolist(), ["no_placement_hg19", "no_placement_hg19"])

    def test_places_retired_rsids_under_their_current_id_but_keeps_w_hm3_label(self):
        paths = write_sources(
            self.tmpdir,
            hm3=[("rs10", "A", "G"), ("rs11", "C", "T"), ("rs14", "A", "C")],
            hg19=[
                (accession(2, "hg19"), 100, "rs20", "A", "G"),
                (accession(2, "hg19"), 200, "rs13", "C", "T"),
                (accession(2, "hg19"), 300, "rs15", "A", "C"),
            ],
            hg38=[
                (accession(2, "hg38"), 1100, "rs20", "A", "G"),
                (accession(2, "hg38"), 1200, "rs13", "C", "T"),
                (accession(2, "hg38"), 1300, "rs15", "A", "C"),
            ],
            merges={"10": ["20"], "11": ["12"], "12": ["13"], "14": ["15", "16"]},
        )

        result = self.tool["build_hm3_curated_map"](**paths)

        self.assertEqual(result["SNP"].tolist(), ["rs10", "rs11"])
        self.assertEqual(result["hg19_POS"].tolist(), [100, 200])
        self.assertEqual(result["hg38_POS"].tolist(), [1100, 1200])
        drops = read_drop_report(paths)
        self.assertEqual(drops["SNP"].tolist(), ["rs14"])
        self.assertEqual(drops["reason"].tolist(), ["ambiguous_merge"])

    def test_skips_snps_that_break_map_identity_rules(self):
        one = lambda pos, rsid, ref, alt, chrom=1, build="hg19": (accession(chrom, build), pos, rsid, ref, alt)
        paths = write_sources(
            self.tmpdir,
            hm3=[
                ("rs30", "A", "G"),  # two accepted hg38 placements
                ("rs31", "A", "G"),  # chr1 in hg19, chr2 in hg38
                ("rs32", "A", "T"),  # strand-ambiguous
                ("rs33", "C", "T"),  # shares hg19 chr1:700 with rs34
                ("rs34", "C", "T"),
                ("rs35", "A", "C"),  # no hg38 placement
                ("rs36", "G", "T"),  # kept
                ("rs40", "A", "G"),  # merged into rs41, which is also listed
                ("rs41", "A", "G"),
            ],
            hg19=[
                one(300, "rs30", "A", "G"), one(310, "rs31", "A", "G"), one(320, "rs32", "A", "T"),
                one(700, "rs33", "C", "T"), one(700, "rs34", "C", "T"), one(350, "rs35", "A", "C"),
                one(360, "rs36", "G", "T"), one(410, "rs41", "A", "G"),
            ],
            hg38=[
                one(1300, "rs30", "A", "G", build="hg38"), one(1900, "rs30", "A", "G", build="hg38"),
                one(1310, "rs31", "A", "G", chrom=2, build="hg38"), one(1320, "rs32", "A", "T", build="hg38"),
                one(1700, "rs33", "C", "T", build="hg38"), one(1701, "rs34", "C", "T", build="hg38"),
                one(1360, "rs36", "G", "T", build="hg38"), one(1410, "rs41", "A", "G", build="hg38"),
            ],
            merges={"40": ["41"]},
        )

        result = self.tool["build_hm3_curated_map"](**paths)

        self.assertEqual(result["SNP"].tolist(), ["rs36"])
        drops = read_drop_report(paths)
        self.assertEqual(
            dict(zip(drops["SNP"], drops["reason"])),
            {
                "rs30": "multiple_placements_hg38",
                "rs31": "chromosome_mismatch",
                "rs32": "strand_ambiguous_allele",
                "rs33": "duplicate_coordinate_hg19",
                "rs34": "duplicate_coordinate_hg19",
                "rs35": "no_placement_hg38",
                "rs40": "duplicate_coordinate_hg19",
                "rs41": "duplicate_coordinate_hg19",
            },
        )
        self.assertEqual(drops.set_index("SNP").at["rs30", "detail"], "1:1300;1:1900")

    def test_rejects_an_assembly_report_that_matches_no_vcf_chromosome(self):
        paths = write_sources(
            self.tmpdir,
            hm3=[("rs1", "A", "G")],
            hg19=[(accession(1, "hg19"), 100, "rs1", "A", "G")],
            hg38=[(accession(1, "hg38"), 1100, "rs1", "A", "G")],
        )
        paths["assembly_report_hg19"] = paths["assembly_report_hg38"]

        with self.assertRaisesRegex(Exception, "no VCF record lies on an autosome"):
            self.tool["build_hm3_curated_map"](**paths)
        self.assertFalse(paths["output_path"].exists())

    def test_command_line_writes_map_and_drop_report_and_summarizes_counts(self):
        paths = write_sources(
            self.tmpdir,
            hm3=[("rs1", "A", "G"), ("rs2", "A", "T")],
            hg19=[(accession(1, "hg19"), 100, "rs1", "A", "G")],
            hg38=[(accession(1, "hg38"), 1100, "rs1", "A", "G")],
        )
        argv = []
        for name, path in paths.items():
            flag = {"output_path": "--output", "drop_report_path": "--drop-report"}.get(name, "--" + name.replace("_", "-"))
            argv += [flag, str(path)]

        with contextlib.redirect_stdout(io.StringIO()) as stdout:
            status = self.tool["main"](argv)

        self.assertEqual(status, 0)
        summary = json.loads(stdout.getvalue())
        self.assertEqual(summary["input_snps"], 2)
        self.assertEqual(summary["retained_snps"], 1)
        self.assertEqual(summary["dropped_by_reason"], {"strand_ambiguous_allele": 1})
        self.assertEqual(len(pd.read_csv(paths["output_path"], sep="\t")), 1)
        self.assertEqual(read_drop_report(paths)["SNP"].tolist(), ["rs2"])


def read_drop_report(paths: dict[str, Path]) -> pd.DataFrame:
    drops = pd.read_csv(paths["drop_report_path"], sep="\t", dtype=str, keep_default_na=False)
    return drops.sort_values("SNP").reset_index(drop=True)


if __name__ == "__main__":
    unittest.main()
