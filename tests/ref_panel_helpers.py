"""Reference-panel fixtures exercised through the production preparation interface."""
from ldsc import GlobalConfig, LDScoreConfig, PlinkRefPanel, RefPanelConfig
from ldsc._kernel import ldscore


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
