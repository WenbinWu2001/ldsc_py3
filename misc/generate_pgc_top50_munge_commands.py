#!/usr/bin/env python3
"""
Generate minimal munge-sumstats commands for masked PGC top50 extracts.

The script discovers ``*.txt`` files under the masked extraction directory,
creates cleaned copies without the extraction banner line, and emits bash
commands that use the refactored ``ldsc munge-sumstats`` CLI without SNP
restriction flags. Commands add column flags only when the masked header cannot
be inferred by the current munger.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
import shlex
import subprocess
import sys


DEFAULT_INPUT_ROOT = Path(
    "/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/results_local/"
    "2026-05-10_extract_pgc_GWAS_top50_lines_masked"
)
DEFAULT_OUTPUT_ROOT = Path(
    "/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/results_local/"
    "2026-05-11_munge_pgc_GWAS_top50_minimal"
)
DEFAULT_CLEAN_ROOT = DEFAULT_OUTPUT_ROOT / "_clean_inputs"


@dataclass(frozen=True)
class MungePlan:
    """File-specific arguments needed after automatic header inference."""

    label: str
    extra_args: tuple[str, ...]
    note: str = ""
    requires_n: bool = False


def discover_sumstats(root: Path, include_checkpoints: bool = False) -> list[Path]:
    """Return recursively discovered top50 text files in stable order."""

    files = sorted(path for path in root.rglob("*.txt") if path.is_file())
    if include_checkpoints:
        return files
    return [path for path in files if ".ipynb_checkpoints" not in path.parts]


def is_extract_banner(line: str) -> bool:
    """Return True for the one-line banner added by the top50 extractor."""

    stripped = line.strip()
    return stripped.startswith("=====") and stripped.endswith("=====")


def write_clean_input(raw_path: Path, clean_path: Path) -> bool:
    """Copy a raw extract, dropping only the leading extractor banner if present."""

    clean_path.parent.mkdir(parents=True, exist_ok=True)
    lines = raw_path.read_text(encoding="utf-8").splitlines(keepends=True)
    skipped_banner = bool(lines and is_extract_banner(lines[0]))
    if skipped_banner:
        lines = lines[1:]
    clean_path.write_text("".join(lines), encoding="utf-8")
    return skipped_banner


def read_header(path: Path) -> list[str]:
    """Read the first non-comment header line from a cleaned sumstats file."""

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            stripped = line.strip()
            if stripped:
                return stripped.split()
    raise ValueError(f"could not find a header line in {path}")


def plan_for_header(header: list[str]) -> MungePlan:
    """Choose the smallest extra argument set implied by a masked header."""

    upper = {column.upper() for column in header}
    extra_args = ["--info-min", "0", "--maf-min", "0"]
    notes: list[str] = ["relax INFO/MAF QC because these are 50-line masked extracts"]

    has_signed = bool({"BETA", "LOG_ODDS", "OR", "Z"} & upper)
    if not has_signed and "BEAA" in upper:
        extra_args.extend(["--signed-sumstats", "BEAA,0"])
        notes.append("BEAA is privacy-masked BETA")

    if "REF" in upper and "ALA" in upper and not bool({"A1", "A2", "EA", "NEA"} & upper):
        extra_args.extend(["--a1", "REF", "--a2", "ALA"])
        notes.append("REF/ALA are privacy-masked REF/ALT alleles")

    if "IMPINFO" in upper:
        extra_args.extend(["--ignore", "IMPINFO"])
        notes.append("ignore comma-list IMPINFO in masked cross-disorder file")

    has_case_control_n = bool({"NCA", "NCO", "NCAS", "NCON", "N_CAS", "N_CON"} & upper)
    has_total_n = "N" in upper
    has_daner_old_n = any(column.startswith("FRQ_A_") for column in upper) and any(
        column.startswith("FRQ_U_") for column in upper
    )
    if has_daner_old_n and not has_case_control_n and not has_total_n:
        extra_args.append("--daner-old")
        notes.append("infer case/control N from FRQ_A_/FRQ_U_ column names")

    requires_n = not has_case_control_n and not has_total_n and not has_daner_old_n
    if requires_n:
        notes.append("sample size is not inferable; fill in --N before running")

    label = "requires_n" if requires_n else "minimal"
    return MungePlan(label=label, extra_args=tuple(extra_args), note="; ".join(notes), requires_n=requires_n)


def relative_output_stem(path: Path, root: Path) -> Path:
    """Return the input-relative stem used for clean inputs and output dirs."""

    relative = path.relative_to(root)
    return relative.with_suffix("")


def build_munge_command(raw_path: Path, clean_path: Path, output_dir: Path, plan: MungePlan) -> list[str]:
    """Build one ``ldsc munge-sumstats`` command with no SNP restriction flag."""

    command = [
        "ldsc",
        "munge-sumstats",
        "--raw-sumstats-file",
        str(clean_path),
        "--output-dir",
        str(output_dir),
        "--overwrite",
    ]
    command.extend(plan.extra_args)
    if plan.requires_n:
        command.extend(["--N", "<N>"])
    return command


def quote_command(command: list[str]) -> str:
    """Return a shell-quoted command line."""

    return " ".join(
        token if token == "<N>" or token.startswith("$") else shlex.quote(token)
        for token in command
    )


def emit_bash(
    input_root: Path,
    clean_root: Path,
    output_root: Path,
    include_checkpoints: bool = False,
) -> list[str]:
    """Create cleaned inputs and return a bash script as a list of lines."""

    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "# Run inside conda env ldsc3-dev, or set LDSC_CMD='conda run -n ldsc3-dev ldsc'.",
        'LDSC_CMD="${LDSC_CMD:-ldsc}"',
        "",
    ]
    for raw_path in discover_sumstats(input_root, include_checkpoints=include_checkpoints):
        stem = relative_output_stem(raw_path, input_root)
        clean_path = clean_root / stem
        output_dir = output_root / stem
        skipped_banner = write_clean_input(raw_path, clean_path)
        header = read_header(clean_path)
        plan = plan_for_header(header)
        command = build_munge_command(raw_path, clean_path, output_dir, plan)
        command[0] = "$LDSC_CMD"

        lines.append(f"# {raw_path.relative_to(input_root)}")
        if skipped_banner:
            lines.append("# cleaned input drops the leading extraction banner")
        if plan.note:
            lines.append(f"# {plan.note}")
        prefix = "# TODO set sample size: " if plan.requires_n else ""
        lines.append(prefix + quote_command(command))
        lines.append("")
    return lines


def build_report_rows(
    input_root: Path,
    clean_root: Path,
    output_root: Path,
    include_checkpoints: bool = False,
) -> list[dict[str, str]]:
    """Return one audit row per discovered input file and generated command."""

    rows: list[dict[str, str]] = []
    for raw_path in discover_sumstats(input_root, include_checkpoints=include_checkpoints):
        stem = relative_output_stem(raw_path, input_root)
        clean_path = clean_root / stem
        output_dir = output_root / stem
        skipped_banner = write_clean_input(raw_path, clean_path)
        plan = plan_for_header(read_header(clean_path))
        command = build_munge_command(raw_path, clean_path, output_dir, plan)
        command[0] = "$LDSC_CMD"
        rows.append(
            {
                "relative_path": str(raw_path.relative_to(input_root)),
                "status": "requires_n" if plan.requires_n else "runnable",
                "skipped_banner": str(skipped_banner),
                "extra_args": " ".join(plan.extra_args),
                "note": plan.note,
                "command": quote_command(command),
            }
        )
    return rows


def write_report(path: Path, rows: list[dict[str, str]]) -> None:
    """Write the command audit report as a TSV file."""

    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["relative_path", "status", "skipped_banner", "extra_args", "note", "command"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run_commands(input_root: Path, clean_root: Path, output_root: Path, include_checkpoints: bool = False) -> int:
    """Execute runnable commands and skip files whose N must be supplied."""

    failures = 0
    for raw_path in discover_sumstats(input_root, include_checkpoints=include_checkpoints):
        stem = relative_output_stem(raw_path, input_root)
        clean_path = clean_root / stem
        output_dir = output_root / stem
        write_clean_input(raw_path, clean_path)
        plan = plan_for_header(read_header(clean_path))
        rel = raw_path.relative_to(input_root)
        if plan.requires_n:
            print(f"SKIP\t{rel}\trequires --N", file=sys.stderr)
            continue
        command = build_munge_command(raw_path, clean_path, output_dir, plan)
        print(f"RUN\t{rel}", file=sys.stderr)
        result = subprocess.run(command, check=False)
        if result.returncode:
            failures += 1
            print(f"FAIL\t{rel}\treturncode={result.returncode}", file=sys.stderr)
        else:
            print(f"OK\t{rel}", file=sys.stderr)
    return 1 if failures else 0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--clean-root", type=Path, default=DEFAULT_CLEAN_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--include-checkpoints", action="store_true")
    parser.add_argument("--run", action="store_true", help="Run commands that do not require manual --N.")
    parser.add_argument("--bash-out", type=Path, help="Write emitted bash commands to this file.")
    parser.add_argument("--report-out", type=Path, help="Write a TSV audit report for generated commands.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.run:
        return run_commands(args.input_root, args.clean_root, args.output_root, args.include_checkpoints)

    lines = emit_bash(args.input_root, args.clean_root, args.output_root, args.include_checkpoints)
    if args.report_out:
        write_report(
            args.report_out,
            build_report_rows(args.input_root, args.clean_root, args.output_root, args.include_checkpoints),
        )
    text = "\n".join(lines)
    if args.bash_out:
        args.bash_out.parent.mkdir(parents=True, exist_ok=True)
        args.bash_out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
