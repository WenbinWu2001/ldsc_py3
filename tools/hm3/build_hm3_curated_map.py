"""One-time curation tool: rebuild the packaged dual-build HapMap3 SNP map.

Run manually from the repo root (requires an editable install of the package):
    python tools/hm3/build_hm3_curated_map.py \
        --hm3-snplist <sources>/w_hm3/w_hm3.snplist.gz \
        --dbsnp-vcf-hg19 <sources>/dbsnp_b157/VCF/GCF_000001405.25.gz \
        --dbsnp-vcf-hg38 <sources>/dbsnp_b157/VCF/GCF_000001405.40.gz \
        --dbsnp-merged-json <sources>/dbsnp_b157/JSON/refsnp-merged.json.bz2 \
        --assembly-report-hg19 <sources>/ncbi_assembly/GCF_000001405.25_GRCh37.p13_assembly_report.txt \
        --assembly-report-hg38 <sources>/ncbi_assembly/GCF_000001405.40_GRCh38.p14_assembly_report.txt \
        --output src/ldsc/data/hm3_curated_map.tsv.gz \
        --drop-report <build-dir>/hm3_curated_map.dropped.tsv.gz

Afterwards rebuild the genome-build inference subset with
tools/hm3/build_hm3_chr_pos_reference.py. Like that builder, this tool is not
shipped with the package.

Sources (all pinned; accessed 2026-10-27)
-----------------------------------------
w_hm3.snplist
    Legacy LDSC HapMap3 SNP list (SNP, A1, A2; 1,217,311 SNPs) distributed by the
    Alkes group. The original https://data.broadinstitute.org/alkesgroup/
    LDSCORE/w_hm3.snplist.bz2 now returns 404; use the gzipped mirror
    https://doi.org/10.5281/zenodo.7773502 (md5
    153ecc2bcfa740afafe656e6a384d769).
dbSNP build 157 RefSNP VCFs, released 2025-01-15
    https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/GCF_000001405.25.gz
    (GRCh37.p13, used for hg19) and GCF_000001405.40.gz (GRCh38.p14, used for
    hg38), each with a published .md5. CHROM holds RefSeq accessions such as
    NC_000001.10, and POS is one-based.
dbSNP build 157 merge history
    https://ftp.ncbi.nih.gov/snp/archive/b157/JSON/refsnp-merged.json.bz2, one
    JSON record per retired rsID with ``merged_snapshot_data.merged_into``.
NCBI assembly reports
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/ for
    GCF_000001405.25_GRCh37.p13 and GCF_000001405.40_GRCh38.p14, verified
    against each directory's md5checksums.txt. They map RefSeq accessions of
    assembled autosomes to chromosome names.

Column derivation
-----------------
SNP, A1, A2
    Copied from w_hm3.snplist, A1/A2 uppercased. SNP keeps the w_hm3 label even
    when dbSNP has merged it into a newer rsID; the newer rsID is used only to
    look up placements.
CHR, hg19_POS, hg38_POS
    The single accepted dbSNP b157 placement of the SNP's current rsID in each
    build, with the accession translated to a chromosome name. No liftover is
    used. CHR must agree between builds.

A placement is accepted when it lies on an autosome's assembled molecule and
one of its single-base REF/ALT pairs reduces to the same allele-set token as
the w_hm3 A1/A2 pair. The token is the package's allele identity
(``ldsc._core.snp_identity``): unordered, with a strand flip treated as the
same pair, so a dbSNP T/C placement matches a w_hm3 A/G SNP.

Skipped SNPs are written to the drop report with one reason each.
"""

from __future__ import annotations

import argparse
import bz2
import gzip
import json
from os import PathLike
from pathlib import Path
import re

import pandas as pd

from ldsc._core.snp_identity import allele_set_series, normalize_allele_set
from ldsc.errors import LDSCInputError

AUTOSOMES = tuple(str(chrom) for chrom in range(1, 23))
BUILDS = ("hg19", "hg38")
OUTPUT_COLUMNS = ["CHR", "hg19_POS", "hg38_POS", "SNP", "A1", "A2"]
DROP_REPORT_COLUMNS = ["SNP", "resolved_SNP", "reason", "detail"]
_REFSNP_ID = re.compile(rb'"refsnp_id"\s*:\s*"(\d+)"')


def build_hm3_curated_map(
    hm3_snplist: str | PathLike[str],
    dbsnp_vcf_hg19: str | PathLike[str],
    dbsnp_vcf_hg38: str | PathLike[str],
    dbsnp_merged_json: str | PathLike[str],
    assembly_report_hg19: str | PathLike[str],
    assembly_report_hg38: str | PathLike[str],
    output_path: str | PathLike[str],
    drop_report_path: str | PathLike[str],
) -> pd.DataFrame:
    """
    Build and write the six-column dual-build HM3 map and its drop report.

    Returns the written map: one row per retained w_hm3 SNP with columns
    ``CHR, hg19_POS, hg38_POS, SNP, A1, A2``, sorted by chromosome and hg19
    position.
    """
    hm3 = _read_hm3_snplist(hm3_snplist)
    _resolve_merges(hm3, dbsnp_merged_json)
    wanted = _wanted_tokens(hm3.loc[hm3["reason"].isna()])
    vcfs = {"hg19": dbsnp_vcf_hg19, "hg38": dbsnp_vcf_hg38}
    reports = {"hg19": assembly_report_hg19, "hg38": assembly_report_hg38}
    placements = {
        build: _scan_dbsnp_vcf(vcfs[build], reports[build], _read_autosome_accessions(reports[build]), wanted)
        for build in BUILDS
    }
    for build in BUILDS:
        _assign_build_placement(hm3, placements[build], build)
    _drop_chromosome_mismatches(hm3)
    _drop_duplicate_coordinates(hm3)

    kept = hm3["reason"].isna()
    final = hm3.loc[kept].assign(CHR=hm3["hg19_CHR"]).loc[:, OUTPUT_COLUMNS]
    final["hg19_POS"] = final["hg19_POS"].astype("int64")
    final["hg38_POS"] = final["hg38_POS"].astype("int64")
    final = final.sort_values(["CHR", "hg19_POS"], key=_genomic_sort_key).reset_index(drop=True)
    _write_tsv(final, output_path)
    _write_tsv(hm3.loc[~kept, DROP_REPORT_COLUMNS].sort_values(["reason", "SNP"]), drop_report_path)
    return final


def _read_hm3_snplist(path: str | PathLike[str]) -> pd.DataFrame:
    """
    Read the w_hm3 SNP list with its allele-set token and drop bookkeeping.

    Pairs the package cannot use as an allele identity (missing, non-ACGT,
    identical, or strand-ambiguous) are dropped under the package's reason.
    """
    frame = pd.read_csv(path, sep=r"\s+", dtype=str)
    frame = frame.loc[:, ["SNP", "A1", "A2"]].copy()
    frame["A1"] = frame["A1"].str.upper()
    frame["A2"] = frame["A2"].str.upper()
    tokens, reasons = allele_set_series(frame, context=f"HM3 SNP list {path}")
    frame["allele_set"] = tokens
    frame["reason"] = reasons.astype("object").where(reasons.notna(), pd.NA)
    frame["detail"] = ""
    invalid = frame["reason"].notna()
    frame.loc[invalid, "detail"] = "A1=" + frame.loc[invalid, "A1"] + ";A2=" + frame.loc[invalid, "A2"]
    return frame


def _resolve_merges(hm3: pd.DataFrame, merged_json: str | PathLike[str]) -> None:
    """
    Set ``resolved_SNP`` to the current dbSNP rsID of each w_hm3 SNP.

    Merges are followed to their end. A record merged into more than one rsID,
    or a merge loop, has no single current identity and is dropped as
    ``ambiguous_merge``.
    """
    merges: dict[str, list[str]] = {}
    queried: set[str] = set()
    pending = {snp[2:] for snp in hm3["SNP"] if snp.startswith("rs") and snp[2:].isdigit()}
    while pending:
        queried |= pending
        found = _read_merge_targets(merged_json, pending)
        merges.update(found)
        pending = {target for targets in found.values() for target in targets} - queried

    resolved: list[str] = []
    for row, snp in hm3["SNP"].items():
        current, seen = snp[2:], set()
        while current in merges:
            targets = merges[current]
            if len(targets) != 1 or current in seen:
                if pd.isna(hm3.at[row, "reason"]):
                    hm3.at[row, "reason"] = "ambiguous_merge"
                    hm3.at[row, "detail"] = "merged_into=" + ",".join(f"rs{target}" for target in targets)
                break
            seen.add(current)
            current = targets[0]
        resolved.append(snp if current == snp[2:] else f"rs{current}")
    hm3["resolved_SNP"] = resolved


def _read_merge_targets(path: str | PathLike[str], ids: set[str]) -> dict[str, list[str]]:
    """Return ``merged_into`` rsID numbers for the requested retired rsID numbers."""
    found: dict[str, list[str]] = {}
    with bz2.open(path, "rb") as handle:
        for line in handle:
            match = _REFSNP_ID.search(line)
            if match is None or match.group(1).decode() not in ids:
                continue
            record = json.loads(line)
            snapshot = record.get("merged_snapshot_data")
            if not snapshot or "merged_into" not in snapshot:
                raise LDSCInputError(
                    f"dbSNP merge record rs{record['refsnp_id']} in '{path}' has no merged_snapshot_data.merged_into. "
                    "Most likely the file is not refsnp-merged.json from a dbSNP RefSNP JSON release. "
                    "Pass the refsnp-merged.json.bz2 file from the same dbSNP build as the VCFs."
                )
            found[str(record["refsnp_id"])] = [str(target) for target in snapshot["merged_into"]]
    return found


def _wanted_tokens(hm3: pd.DataFrame) -> dict[bytes, set[str]]:
    """Map each rsID to look up in dbSNP to the allele-set tokens it must match."""
    wanted: dict[bytes, set[str]] = {}
    for rsid, token in zip(hm3["resolved_SNP"], hm3["allele_set"]):
        wanted.setdefault(str(rsid).encode(), set()).add(str(token))
    return wanted


def _read_autosome_accessions(path: str | PathLike[str]) -> dict[bytes, str]:
    """Map RefSeq accessions of assembled autosomes to chromosome names."""
    header: list[str] | None = None
    accessions: dict[bytes, str] = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith("# Sequence-Name"):
                header = line[2:].split("\t")
                continue
            if line.startswith("#") or not line or header is None:
                continue
            row = dict(zip(header, line.split("\t")))
            if row["Sequence-Role"] == "assembled-molecule" and row["Assigned-Molecule"] in AUTOSOMES:
                accessions[row["RefSeq-Accn"].encode()] = row["Assigned-Molecule"]
    return accessions


def _scan_dbsnp_vcf(
    path: str | PathLike[str],
    assembly_report: str | PathLike[str],
    chromosomes: dict[bytes, str],
    wanted: dict[bytes, set[str]],
) -> pd.DataFrame:
    """
    Return allele-matched autosomal placements of the wanted rsIDs.

    Fails when no record of a wanted rsID lies on an autosome named by the
    assembly report: that pairing would otherwise drop every SNP silently.
    """
    records: list[tuple[str, str, str, int]] = []
    n_on_autosomes = 0
    with gzip.open(path, "rb") as handle:
        for line in handle:
            if line[0] == 35:  # "#"
                continue
            # About a billion records per build; slice only the ID column until it is wanted.
            id_start = line.find(b"\t", line.find(b"\t") + 1) + 1
            ids = line[id_start:line.find(b"\t", id_start)]
            if ids not in wanted and b";" not in ids:
                continue
            chrom_accession, pos, ids, ref, alts = line.split(b"\t", 5)[:5]
            for rsid in ids.split(b";"):
                tokens = wanted.get(rsid)
                if tokens is None:
                    continue
                chrom = chromosomes.get(chrom_accession)
                if chrom is None:
                    continue
                n_on_autosomes += 1
                if len(ref) != 1:
                    continue
                for alt in alts.split(b","):
                    if len(alt) != 1:
                        continue
                    try:
                        token = normalize_allele_set(ref.decode(), alt.decode())
                    except ValueError:
                        continue
                    if token in tokens:
                        records.append((rsid.decode(), token, chrom, int(pos)))
    if wanted and n_on_autosomes == 0:
        raise LDSCInputError(
            f"Cannot place HM3 SNPs with '{path}': no VCF record lies on an autosome named by assembly report "
            f"'{assembly_report}'. Most likely the report belongs to the other genome build or the VCF names "
            "chromosomes differently. Pair GCF_000001405.25 with the GRCh37.p13 report and GCF_000001405.40 "
            "with the GRCh38.p14 report."
        )
    columns = ["resolved_SNP", "allele_set", "CHR", "POS"]
    return pd.DataFrame(records, columns=columns).drop_duplicates()


def _assign_build_placement(hm3: pd.DataFrame, placements: pd.DataFrame, build: str) -> None:
    """Record each active SNP's single accepted placement in ``build`` or its drop reason."""
    hm3[f"{build}_CHR"] = pd.NA
    hm3[f"{build}_POS"] = pd.NA
    active = hm3.index[hm3["reason"].isna()]
    rows = hm3.loc[active, ["resolved_SNP", "allele_set"]].rename_axis("row").reset_index()
    matched = rows.merge(placements, on=["resolved_SNP", "allele_set"], how="inner")
    counts = matched.groupby("row").size().reindex(active, fill_value=0)
    hm3.loc[counts.index[counts == 0], "reason"] = f"no_placement_{build}"

    multiple = matched.loc[matched["row"].map(counts) > 1]
    if not multiple.empty:
        multiple = multiple.assign(_chrom=multiple["CHR"].astype(int)).sort_values(["row", "_chrom", "POS"])
        sites = (multiple["CHR"] + ":" + multiple["POS"].astype(str)).groupby(multiple["row"]).agg(";".join)
        hm3.loc[sites.index, "reason"] = f"multiple_placements_{build}"
        hm3.loc[sites.index, "detail"] = sites

    single = matched.loc[matched["row"].map(counts) == 1].set_index("row")
    hm3.loc[single.index, f"{build}_CHR"] = single["CHR"]
    hm3.loc[single.index, f"{build}_POS"] = single["POS"]


def _drop_chromosome_mismatches(hm3: pd.DataFrame) -> None:
    """Drop SNPs placed on different chromosomes in hg19 and hg38."""
    mismatch = hm3["reason"].isna() & (hm3["hg19_CHR"] != hm3["hg38_CHR"])
    hm3.loc[mismatch, "reason"] = "chromosome_mismatch"
    hm3.loc[mismatch, "detail"] = "hg19=" + hm3.loc[mismatch, "hg19_CHR"] + ";hg38=" + hm3.loc[mismatch, "hg38_CHR"]


def _drop_duplicate_coordinates(hm3: pd.DataFrame) -> None:
    """
    Drop every SNP sharing a (CHR, POS) with another SNP in either build.

    Both builds are checked on the same candidate set, so the result does not
    depend on build order; a SNP duplicated in both builds reports hg19.
    """
    active = hm3["reason"].isna()
    duplicated = {
        build: active & hm3.loc[:, [f"{build}_CHR", f"{build}_POS"]].where(active).duplicated(keep=False)
        for build in BUILDS
    }
    for build in reversed(BUILDS):
        mask = duplicated[build]
        hm3.loc[mask, "reason"] = f"duplicate_coordinate_{build}"
        hm3.loc[mask, "detail"] = hm3.loc[mask, f"{build}_CHR"] + ":" + hm3.loc[mask, f"{build}_POS"].astype(str)


def _genomic_sort_key(column: pd.Series) -> pd.Series:
    return column.astype("int64")


def _write_tsv(frame: pd.DataFrame, path: str | PathLike[str]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, compression="gzip")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Rebuild the packaged six-column HapMap3 map from w_hm3.snplist and dbSNP b157."
    )
    parser.add_argument("--hm3-snplist", required=True, help="w_hm3.snplist(.gz) with SNP, A1, A2.")
    parser.add_argument("--dbsnp-vcf-hg19", required=True, help="dbSNP b157 GCF_000001405.25.gz (GRCh37.p13).")
    parser.add_argument("--dbsnp-vcf-hg38", required=True, help="dbSNP b157 GCF_000001405.40.gz (GRCh38.p14).")
    parser.add_argument("--dbsnp-merged-json", required=True, help="dbSNP b157 refsnp-merged.json.bz2.")
    parser.add_argument("--assembly-report-hg19", required=True, help="GCF_000001405.25_GRCh37.p13 assembly report.")
    parser.add_argument("--assembly-report-hg38", required=True, help="GCF_000001405.40_GRCh38.p14 assembly report.")
    parser.add_argument("--output", required=True, help="Map to write, e.g. src/ldsc/data/hm3_curated_map.tsv.gz.")
    parser.add_argument("--drop-report", required=True, help="Skipped SNPs with one reason each (.tsv.gz).")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Build the map and print input, retained, and per-reason drop counts as JSON."""
    args = build_arg_parser().parse_args(argv)
    final = build_hm3_curated_map(
        hm3_snplist=args.hm3_snplist,
        dbsnp_vcf_hg19=args.dbsnp_vcf_hg19,
        dbsnp_vcf_hg38=args.dbsnp_vcf_hg38,
        dbsnp_merged_json=args.dbsnp_merged_json,
        assembly_report_hg19=args.assembly_report_hg19,
        assembly_report_hg38=args.assembly_report_hg38,
        output_path=args.output,
        drop_report_path=args.drop_report,
    )
    drops = pd.read_csv(args.drop_report, sep="\t", dtype=str)
    summary = {
        "input_snps": len(final) + len(drops),
        "retained_snps": len(final),
        "dropped_by_reason": drops["reason"].value_counts().sort_index().to_dict(),
        "output": str(args.output),
        "drop_report": str(args.drop_report),
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
