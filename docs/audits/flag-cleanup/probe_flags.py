"""Reproduce flag-audit observations without changing package code.

Run from the repository root with the ldsc3-dev Python environment.
All generated scientific artifacts are confined to a temporary directory.
"""
from contextlib import redirect_stdout, redirect_stderr
from io import StringIO
import json
import tempfile
from pathlib import Path

from ldsc.cli import build_parser
from ldsc.sumstats_munger import main as munge_main
from ldsc.h2_scale import _population_prevalence_grid


def main():
    observations = []
    with tempfile.TemporaryDirectory(prefix="ldsc-flag-audit-") as temporary:
        root = Path(temporary)
        def run(name, header, rows, flags=()):
            raw = root / f"{name}.tsv"
            raw.write_text(header + "\n" + "\n".join(rows) + "\n")
            args = ["--raw-sumstats-file", str(raw), "--output-dir", str(root / name),
                    "--snp-identifier", "rsid", *flags]
            with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
                try:
                    result = munge_main(args)
                    if hasattr(result, "data"):
                        detail = result.data.to_dict(orient="list")
                    else:
                        detail = {"runnable": result.runnable, "suggested_args": result.suggested_args,
                                  "missing_fields": result.missing_fields}
                    observations.append({"case": name, "outcome": "success", "detail": detail})
                except Exception as exc:
                    observations.append({"case": name, "outcome": type(exc).__name__, "detail": str(exc)})
        run("constant_n_with_n_column", "SNP P BETA N", ["rs1 .05 -.05 1000", "rs2 .05 .05 1000"], ["--N", "9999"])
        run("constant_case_control", "SNP P BETA", ["rs1 .05 -.05", "rs2 .05 .05"], ["--N-cas", "400", "--N-con", "600"])
        run("n_min_default", "SNP P BETA N", ["rs1 .05 -.05 550", "rs2 .05 .05 1000"])
        run("n_min_zero", "SNP P BETA N", ["rs1 .05 -.05 550", "rs2 .05 .05 1000"], ["--n-min", "0"])
        run("n_min_constant", "SNP P BETA", ["rs1 .05 -.05", "rs2 .05 .05"], ["--N", "1000", "--n-min", "2000"])
        run("nstudy_with_n", "SNP P BETA N NSTUDY", ["rs1 .05 -.05 1000 1", "rs2 .05 .05 1000 2"], ["--nstudy-min", "2"])
        run("nstudy_zero", "SNP P BETA NSTUDY", ["rs1 .05 -.05 1", "rs2 .05 .05 2"], ["--N", "1000", "--nstudy-min", "0"])
        run("keep_maf", "SNP P BETA N FRQ", ["rs1 .05 -.05 1000 .8", "rs2 .05 .05 1000 .7"], ["--keep-maf"])
        run("a1_inc", "SNP P BETA N", ["rs1 .05 -.05 1000", "rs2 .05 .05 1000"], ["--a1-inc"])
        run("multi_info_list", "SNP P BETA N Q1 Q2", ["rs1 .05 -.05 1000 .99,.99 .99,.99", "rs2 .05 .05 1000 .99,.99 .99,.99"], ["--info-list", "Q1,Q2"])
        run("single_info_list", "SNP P BETA N Q1", ["rs1 .05 -.05 1000 .99,NA", "rs2 .05 .05 1000 .99,.99"], ["--info-list", "Q1"])
        run("daner_auto_without_frequency", "SNP P BETA Nca Nco", ["rs1 .05 -.05 40 60", "rs2 .05 .05 40 60"])
        run("daner_explicit_without_frequency", "SNP P BETA Nca Nco", ["rs1 .05 -.05 40 60", "rs2 .05 .05 40 60"], ["--format", "daner-new"])
        run("daner_explicit_aliases", "SNP P BETA NCAS NCON FRQ_U_60", ["rs1 .05 -.05 40 60 .3", "rs2 .05 .05 40 60 .3"], ["--format", "daner-new"])
        run("infer_drops_explicit_options", "SNP P BETA N", ["rs1 .05 -.05 1000", "rs2 .05 .05 1000"], ["--infer-only", "--n-min", "1", "--keep-maf", "--trait-name", "example"])
        run("infer_missing_p", "SNP BETA N", ["rs1 -.05 1000", "rs2 .05 1000"], ["--infer-only"])
        run("infer_ambiguous_n", "SNP P BETA N NCAS NCON", ["rs1 .05 -.05 1000 40 60", "rs2 .05 .05 1000 40 60"], ["--infer-only"])
        run("info_threshold_without_info", "SNP P BETA N", ["rs1 .05 -.05 1000", "rs2 .05 .05 1000"], ["--info-min", "10"])
        run("maf_threshold_without_frequency", "SNP P BETA N", ["rs1 .05 -.05 1000", "rs2 .05 .05 1000"], ["--maf-min", ".5"])
    parser = build_parser()
    for command, required in (("ldscore", ["--output-dir", "unused"]),
                              ("build-gene-ldscore-index", ["--baseline-annot-sources", "unused", "--plink-prefix", "unused", "--output-dir", "unused", "--gene-coordinate-file", "unused", "--genome-build", "hg19", "--snp-identifier", "rsid"])):
        parsed = parser.parse_args([command, *required, "--exclude-regions", "none"])
        observations.append({"case": command + "_hidden_alias", "outcome": "success", "detail": parsed.regr_snps_exclude_regions})
    values, mode = _population_prevalence_grid(pop_prev=.1, pop_prev_range=None, num_points=-7)
    observations.append({"case": "num_points_exact_prevalence", "outcome": "success", "detail": {"values": values.tolist(), "mode": mode}})
    print(json.dumps(observations, indent=2, default=str))


if __name__ == "__main__":
    main()
