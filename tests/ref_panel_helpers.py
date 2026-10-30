"""Reference-panel fixtures exercised through the production preparation interface."""
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
