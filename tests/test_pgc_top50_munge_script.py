from __future__ import annotations

import importlib.util
from pathlib import Path
import sys


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "misc" / "generate_pgc_top50_munge_commands.py"


def load_script():
    spec = importlib.util.spec_from_file_location("generate_pgc_top50_munge_commands", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_discovers_txt_files_recursively_and_skips_checkpoints(tmp_path):
    script = load_script()
    root = tmp_path / "sumstats"
    keep = root / "mdd" / "trait.top50.txt"
    skip = root / "mdd" / ".ipynb_checkpoints" / "trait-checkpoint.top50.txt"
    other = root / "mdd" / "notes.tsv"
    keep.parent.mkdir(parents=True)
    skip.parent.mkdir(parents=True)
    keep.write_text("header\n", encoding="utf-8")
    skip.write_text("header\n", encoding="utf-8")
    other.write_text("header\n", encoding="utf-8")

    assert script.discover_sumstats(root) == [keep]


def test_cleaned_input_skips_only_extract_banner(tmp_path):
    script = load_script()
    raw = tmp_path / "raw.top50.txt"
    clean = tmp_path / "clean.top50.txt"
    raw.write_text("===== trait =====\n##metadata\nCHR SNP BP OR P\n", encoding="utf-8")

    skipped = script.write_clean_input(raw, clean)

    assert skipped is True
    assert clean.read_text(encoding="utf-8") == "##metadata\nCHR SNP BP OR P\n"


def test_command_omits_restriction_snps_and_adds_only_needed_flags(tmp_path):
    script = load_script()
    raw = tmp_path / "trait.top50.txt"
    clean = tmp_path / "clean" / "trait.top50.txt"
    output_dir = tmp_path / "out" / "trait"
    header = ["CHROM", "POS", "ID", "REF", "ALA", "BEAA", "SE", "PVAL", "NCAS", "NCON"]
    plan = script.plan_for_header(header)

    command = script.build_munge_command(raw, clean, output_dir, plan)

    assert "--sumstats-snps-file" not in command
    assert "--use-hm3-snps" not in command
    assert command[:7] == [
        "ldsc",
        "munge-sumstats",
        "--raw-sumstats-file",
        str(clean),
        "--output-dir",
        str(output_dir),
        "--overwrite",
    ]
    assert ["--info-min", "0"] == command[7:9]
    assert ["--maf-min", "0"] == command[9:11]
    assert "--signed-sumstats" in command
    assert "BEAA,0" in command


def test_plan_handles_masked_ref_alt_and_list_info():
    script = load_script()
    header = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALA",
        "BEAA",
        "SE",
        "PVAL",
        "IMPINFO",
        "NCAS",
        "NCON",
    ]

    plan = script.plan_for_header(header)

    assert plan.requires_n is False
    assert "--a1" in plan.extra_args
    assert "--a2" in plan.extra_args
    assert plan.extra_args[-2:] == ("--ignore", "IMPINFO")


def test_neff_column_is_not_treated_as_inferable_n():
    script = load_script()
    header = ["#CHROM", "ID", "POS", "A1", "A2", "FREQ", "NEFF", "Z", "P", "DIRE"]

    plan = script.plan_for_header(header)

    assert plan.requires_n is True


def test_emitted_bash_keeps_ldsc_cmd_expandable(tmp_path):
    script = load_script()
    root = tmp_path / "sumstats"
    raw = root / "trait.top50.txt"
    raw.parent.mkdir(parents=True)
    raw.write_text("===== trait =====\nCHR SNP BP A1 A2 OR P N\n", encoding="utf-8")

    lines = script.emit_bash(root, tmp_path / "clean", tmp_path / "out")

    command_lines = [line for line in lines if "munge-sumstats" in line]
    assert command_lines
    assert command_lines[0].startswith("$LDSC_CMD munge-sumstats")


def test_build_report_rows_marks_manual_n_and_runnable_files(tmp_path):
    script = load_script()
    root = tmp_path / "sumstats"
    runnable = root / "new" / "trait.top50.txt"
    manual_n = root / "old" / "trait.top50.txt"
    runnable.parent.mkdir(parents=True)
    manual_n.parent.mkdir(parents=True)
    runnable.write_text("===== trait =====\nCHR SNP BP A1 A2 OR P N\n", encoding="utf-8")
    manual_n.write_text("===== trait =====\nsnpid hg18chr bp a1 a2 or se pval info\n", encoding="utf-8")

    rows = script.build_report_rows(root, tmp_path / "clean", tmp_path / "out")

    statuses = {row["relative_path"]: row["status"] for row in rows}
    assert statuses["new/trait.top50.txt"] == "runnable"
    assert statuses["old/trait.top50.txt"] == "requires_n"
    assert all("--sumstats-snps-file" not in row["command"] for row in rows)
