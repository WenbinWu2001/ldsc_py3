"""Reference-panel fixtures exercised through the production preparation interface."""
from pathlib import Path

import numpy as np

from ldsc import GlobalConfig, LDScoreConfig, PlinkRefPanel, RefPanelConfig
from ldsc._core.snp_identity import sidecar_identity_sha256
from ldsc._ldscore import compute as ldscore
from ldsc._refpanel.r2_format import (
    build_runtime_metadata_table,
    r2_matrix_file_paths,
    write_r2_parquet,
    write_runtime_metadata_sidecar,
)


def write_r2_matrix_chromosome(build_dir, panel, pairs=(), *, snp_identifier="chr_pos", genome_build="hg19",
                               n_samples=100, ld_window=("kb", 1000.0), r2_min=0.0, row_group_size=50_000):
    """Write one chromosome of an R2 matrix (R2 Parquet file and sidecar) with the package writers.

    ``panel`` has CHR, POS, SNP, A1, A2 and MAF columns, and CM when the sidecar
    should record it; ``pairs`` are ``{"i", "j", "R2", "sign"}`` rows over panel
    rows, with ``sign`` ``"+"`` or ``"-"``. Returns ``(parquet_path, sidecar_path)``.
    """
    pairs = list(pairs)
    chunks = [] if not pairs else [(
        np.array([row["i"] for row in pairs], dtype=np.int64), np.array([row["j"] for row in pairs], dtype=np.int64),
        np.array([row["R2"] for row in pairs], dtype=np.float32),
        np.array([1 if row["sign"] == "+" else -1 for row in pairs], dtype=np.int8),
    )]
    parquet_path, sidecar_path = r2_matrix_file_paths(build_dir, str(panel["CHR"].iloc[0]))
    sidecar = build_runtime_metadata_table(
        metadata=panel, positions=panel["POS"].to_numpy(),
        cm_values=panel["CM"].to_numpy() if "CM" in panel else None,
    )
    write_runtime_metadata_sidecar(sidecar, sidecar_path, genome_build=genome_build, snp_identifier=snp_identifier)
    write_r2_parquet(
        pair_chunks=chunks, path=parquet_path, genome_build=genome_build, n_samples=n_samples,
        snp_identifier=snp_identifier, n_snps=len(sidecar), sidecar_identity_sha256=sidecar_identity_sha256(sidecar),
        ld_window_mode=ld_window[0], ld_window_value=ld_window[1], r2_min=r2_min, row_group_size=row_group_size,
    )
    return parquet_path, sidecar_path


def write_tiny_plink(prefix, metadata):
    """Write a four-sample polymorphic panel on the supplied SNP coordinates."""
    prefix.with_suffix(".bim").write_text("".join(
        f"{row.CHR} {row.SNP} 0 {row.POS} A C\n" for row in metadata.itertuples(index=False)
    ))
    prefix.with_suffix(".fam").write_text("".join(f"F{i} I{i} 0 0 0 -9\n" for i in range(4)))
    prefix.with_suffix(".bed").write_bytes(b"\x6c\x1b\x01" + b"\xf0" * len(metadata))
    return prefix


# Five SNPs on four individuals: (SNP, .bim A1, .bim A2, .bim A1 copies per individual).
# The less frequent allele in these counts is the A1 that panel preparation reports.
ORIENTATION_PANEL = (
    ("rs1", "A", "G", (2, 2, 2, 1)),     # A1 frequency 7/8: the .bim A1 is the major allele
    ("rs2", "C", "T", (0, 1, 0, 0)),     # 1/8: the .bim A1 is already minor
    ("rs3", "A", "C", (2, 2, 0, 0)),     # 1/2: a tie keeps the .bim order
    ("rs4", "G", "T", (2, 2, None, 1)),  # 5/6 among the three individuals with a call
    ("rs5", "C", "A", (2, 2, 1, 0)),     # 5/8 overall, 1/4 among individuals I2 and I3
)


def write_counted_plink(prefix, rows, *, chrom="1"):
    """Write a PLINK trio whose genotypes are counts of the ``.bim`` A1 allele.

    ``rows`` are ``(SNP, A1, A2, counts)`` at positions 10, 20, ...; each count
    is 0, 1 or 2 copies of the ``.bim`` A1, or None for a missing call. The BED
    uses PLINK's SNP-major two-bit codes (00 two A1 copies, 10 heterozygous,
    11 two A2 copies, 01 missing), with the first individual in the lowest
    bits of each byte. Individuals are ``F<i> I<i>``.
    """
    code = {2: 0b00, 1: 0b10, 0: 0b11, None: 0b01}
    n_individuals = len(rows[0][3])
    bed = bytearray(b"\x6c\x1b\x01")
    for *_, counts in rows:
        packed = bytearray((n_individuals + 3) // 4)
        for individual, count in enumerate(counts):
            packed[individual // 4] |= code[count] << 2 * (individual % 4)
        bed += packed
    Path(f"{prefix}.bim").write_text("".join(
        f"{chrom} {snp} 0 {10 * (row + 1)} {a1} {a2}\n" for row, (snp, a1, a2, _) in enumerate(rows)))
    Path(f"{prefix}.fam").write_text("".join(f"F{i} I{i} 0 0 0 -9\n" for i in range(n_individuals)))
    Path(f"{prefix}.bed").write_bytes(bytes(bed))
    return prefix


def prepare_plink(chrom, bundle, args):
    global_options = {"snp_identifier": args.snp_identifier}
    if hasattr(args, "genome_build"):
        global_options["genome_build"] = args.genome_build
    panel = PlinkRefPanel(
        GlobalConfig(**global_options),
        RefPanelConfig(backend="plink", plink_prefix=args.bfile,
                       keep_indivs_file=getattr(args, "keep", None), maf_min=getattr(args, "maf_min", None)),
    )
    config = LDScoreConfig(
        ld_wind_snps=args.ld_wind_snps, ld_wind_kb=args.ld_wind_kb, ld_wind_cm=args.ld_wind_cm,
        whole_chromosome_ok=args.yes_really,
    )
    return panel.prepare_chromosome(chrom, bundle, config, genetic_map=getattr(args, "genetic_map", None))


def compute_plink(chrom, bundle, args, regression_keys, regression_regions=None):
    with prepare_plink(chrom, bundle, args) as prepared:
        return ldscore.compute_chromosome(
            chrom, prepared, snp_identifier=args.snp_identifier, snp_batch_size=args.snp_batch_size,
            common_maf_min=getattr(args, "common_maf_min", 0.05), regression_keys=regression_keys,
            regression_regions=regression_regions,
        )
