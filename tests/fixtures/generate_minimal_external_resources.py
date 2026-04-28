"""Generate tiny external-resource fixtures for reference-panel parity tests.

The generated PLINK files are a SNP subset of the local 1KG 30x chromosome 22
resource. SNPs are selected by intersecting the source ``.bim`` hg38 positions
with the packaged HapMap3 chromosome-position reference.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[2]
WORKSPACE_ROOT = REPO_ROOT.parent
SOURCE_PREFIX = WORKSPACE_ROOT / "resources" / "example_1kg_30x" / "genomes_30x_chr22"
SOURCE_HG19_MAP = (
    WORKSPACE_ROOT
    / "resources"
    / "genetic_maps"
    / "genetic_map_alkesgroup"
    / "genetic_map_hg19_withX.txt"
)
SOURCE_HG38_MAP = (
    WORKSPACE_ROOT
    / "resources"
    / "genetic_maps"
    / "genetic_map_alkesgroup"
    / "genetic_map_hg38_withX.txt"
)
HM3_REFERENCE = REPO_ROOT / "src" / "ldsc" / "data" / "hm3_chr_pos_reference.tsv.gz"
OUTPUT_ROOT = REPO_ROOT / "tests" / "fixtures" / "minimal_external_resources"
OUTPUT_PREFIX = OUTPUT_ROOT / "plink" / "hm3_chr22_subset"
OUTPUT_HG19_MAP = OUTPUT_ROOT / "genetic_maps" / "genetic_map_hg19_chr22_subset.txt"
OUTPUT_HG38_MAP = OUTPUT_ROOT / "genetic_maps" / "genetic_map_hg38_chr22_subset.txt"
OUTPUT_SELECTION = OUTPUT_ROOT / "hm3_chr22_subset_snps.tsv"

CHROM = "22"
N_SNPS = 32


def main() -> None:
    selected = select_hm3_bim_rows(n_snps=N_SNPS)
    write_plink_subset(selected)
    write_map_slice(SOURCE_HG19_MAP, OUTPUT_HG19_MAP, [row["hg19_pos"] for row in selected])
    write_map_slice(SOURCE_HG38_MAP, OUTPUT_HG38_MAP, [row["hg38_pos"] for row in selected])
    write_selection_table(selected)
    print(f"Wrote {len(selected)} SNPs to {OUTPUT_ROOT}")


def load_hm3_chr22_positions() -> dict[int, int]:
    """Return ``hg38_pos -> hg19_pos`` for chromosome 22 HM3 sites."""
    positions: dict[int, int] = {}
    with gzip.open(HM3_REFERENCE, "rt", encoding="utf-8") as handle:
        header = handle.readline().strip().split("\t")
        column_index = {name: idx for idx, name in enumerate(header)}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[column_index["CHR"]] != CHROM:
                continue
            hg19_pos = int(fields[column_index["hg19_POS"]])
            hg38_pos = int(fields[column_index["hg38_POS"]])
            positions[hg38_pos] = hg19_pos
    return positions


def select_hm3_bim_rows(*, n_snps: int) -> list[dict[str, object]]:
    """Select the first ``n_snps`` source BIM rows whose hg38 BP is in HM3."""
    hm3_positions = load_hm3_chr22_positions()
    selected: list[dict[str, object]] = []
    with Path(str(SOURCE_PREFIX) + ".bim").open("rt", encoding="utf-8") as handle:
        for index, line in enumerate(handle):
            fields = line.rstrip("\n").split()
            if len(fields) < 6 or fields[0] != CHROM:
                continue
            hg38_pos = int(fields[3])
            hg19_pos = hm3_positions.get(hg38_pos)
            if hg19_pos is None:
                continue
            selected.append(
                {
                    "source_index": index,
                    "chrom": fields[0],
                    "snp": fields[1],
                    "cm": fields[2],
                    "hg38_pos": hg38_pos,
                    "hg19_pos": hg19_pos,
                    "a1": fields[4],
                    "a2": fields[5],
                    "bim_line": line,
                }
            )
            if len(selected) == n_snps:
                break
    if len(selected) < n_snps:
        raise RuntimeError(f"Only found {len(selected)} HM3-overlapping chr22 SNPs in {SOURCE_PREFIX}.bim")
    return selected


def write_plink_subset(selected: list[dict[str, object]]) -> None:
    """Write a SNP-major PLINK BED/BIM/FAM subset."""
    OUTPUT_PREFIX.parent.mkdir(parents=True, exist_ok=True)
    fam_source = Path(str(SOURCE_PREFIX) + ".fam")
    fam_target = Path(str(OUTPUT_PREFIX) + ".fam")
    fam_text = fam_source.read_text(encoding="utf-8")
    fam_target.write_text(fam_text, encoding="utf-8")
    sample_count = len([line for line in fam_text.splitlines() if line.strip()])
    bytes_per_snp = (sample_count + 3) // 4

    source_bed = Path(str(SOURCE_PREFIX) + ".bed")
    target_bed = Path(str(OUTPUT_PREFIX) + ".bed")
    selected_by_index = {int(row["source_index"]): row for row in selected}
    with source_bed.open("rb") as source, target_bed.open("wb") as target:
        magic = source.read(3)
        if magic != b"\x6c\x1b\x01":
            raise RuntimeError(f"{source_bed} is not a SNP-major PLINK BED file")
        target.write(magic)
        for source_index in sorted(selected_by_index):
            source.seek(3 + source_index * bytes_per_snp)
            target.write(source.read(bytes_per_snp))

    bim_target = Path(str(OUTPUT_PREFIX) + ".bim")
    with bim_target.open("wt", encoding="utf-8") as handle:
        for row in selected:
            handle.write(str(row["bim_line"]))


def write_map_slice(source_path: Path, output_path: Path, positions: Iterable[int]) -> None:
    """Write bracketing map rows around selected positions."""
    target_positions = sorted(set(int(pos) for pos in positions))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not target_positions:
        raise ValueError("positions must not be empty")

    rows: dict[int, str] = {}
    target_idx = 0
    previous_pos: int | None = None
    previous_line: str | None = None
    header = ""
    with source_path.open("rt", encoding="utf-8") as handle:
        header = handle.readline()
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            chrom = fields[0].removeprefix("chr")
            if chrom != CHROM:
                if rows and chrom != CHROM:
                    break
                continue
            pos = int(fields[1])
            while target_idx < len(target_positions) and pos >= target_positions[target_idx]:
                if previous_line is not None and previous_pos is not None:
                    rows[previous_pos] = previous_line
                rows[pos] = line
                target_idx += 1
            previous_pos = pos
            previous_line = line
            if target_idx == len(target_positions):
                # Keep one final right bracket after the last selected position.
                continue
    if previous_line is not None and previous_pos is not None:
        rows[previous_pos] = previous_line
    if not rows:
        raise RuntimeError(f"No chromosome {CHROM} rows were found in {source_path}")

    with output_path.open("wt", encoding="utf-8") as handle:
        handle.write(header)
        for pos in sorted(rows):
            handle.write(rows[pos])


def write_selection_table(selected: list[dict[str, object]]) -> None:
    OUTPUT_SELECTION.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_SELECTION.open("wt", encoding="utf-8") as handle:
        handle.write("CHR\tSNP\tsource_bim_index\thg19_POS\thg38_POS\tA1\tA2\n")
        for row in selected:
            handle.write(
                "{chrom}\t{snp}\t{source_index}\t{hg19_pos}\t{hg38_pos}\t{a1}\t{a2}\n".format(
                    **row
                )
            )


if __name__ == "__main__":
    main()
