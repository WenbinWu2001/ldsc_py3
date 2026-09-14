import os
import warnings
from pathlib import Path
import tempfile
import unittest

import pytest



from ldsc.path_resolution import (
    ensure_output_paths_available,
    ensure_output_directory,
    normalize_path_token,
    normalize_path_tokens,
    resolve_chromosome_group,
    resolve_file_group,
    resolve_plink_prefix,
    resolve_plink_prefix_group,
    resolve_scalar_path,
)
from ldsc.errors import LDSCInputError


class PathResolutionTest(unittest.TestCase):
    def test_requested_directory_creation_is_information(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            destination = Path(tmpdir) / "results"
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                with self.assertLogs("LDSC.path_resolution", level="INFO") as logged:
                    self.assertEqual(ensure_output_directory(destination), destination)
            self.assertTrue(destination.is_dir())
            self.assertEqual(caught, [])
            self.assertIn(str(destination), logged.output[0])

    def test_normalize_path_token_expands_user_and_env(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            os.environ["LDSC_TEST_ROOT"] = tmpdir
            self.assertEqual(
                normalize_path_token("$LDSC_TEST_ROOT/data.tsv"),
                str(Path(tmpdir) / "data.tsv"),
            )
        self.assertEqual(
            normalize_path_token("~/ldsc-test"),
            str(Path.home() / "ldsc-test"),
        )

    def test_normalize_path_tokens_accepts_single_string_token(self):
        self.assertEqual(
            normalize_path_tokens("resources/baseline.*.annot.gz"),
            ("resources/baseline.*.annot.gz",),
        )

    def test_resolve_scalar_path_requires_exactly_one_match(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            match = tmpdir / "trait.raw.tsv"
            match.write_text("x\n", encoding="utf-8")

            resolved = resolve_scalar_path(str(tmpdir / "trait.*.tsv"), label="sumstats")
            self.assertEqual(resolved, str(match))

            second = tmpdir / "trait.extra.tsv"
            second.write_text("x\n", encoding="utf-8")
            with self.assertRaisesRegex(LDSCInputError, "matched 2 files"):
                resolve_scalar_path(str(tmpdir / "trait.*.tsv"), label="sumstats")

    def test_resolve_scalar_path_reports_cause_fix_and_reference(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            (tmpdir / "trait.raw.tsv").write_text("x\n", encoding="utf-8")
            (tmpdir / "trait.extra.tsv").write_text("x\n", encoding="utf-8")

            with self.assertRaisesRegex(
                LDSCInputError,
                "Could not resolve sumstats path.*matched 2 files.*"
                "Most likely the glob is too broad.*"
                "docs/troubleshooting.md#common-input-path-did-not-resolve-to-one-file",
            ):
                resolve_scalar_path(str(tmpdir / "trait.*.tsv"), label="sumstats")

    def test_resolve_file_group_handles_globs_prefixes_and_deduplicates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base1 = tmpdir / "baseline.1.annot.gz"
            base2 = tmpdir / "baseline.2.annot.gz"
            base1.write_text("x\n", encoding="utf-8")
            base2.write_text("x\n", encoding="utf-8")

            resolved = resolve_file_group(
                [str(tmpdir / "baseline.*.annot.gz"), str(base1)],
                label="annotation",
            )

            self.assertEqual(resolved, [str(base1), str(base2)])

    def test_resolve_chromosome_group_supports_explicit_at_tokens_only(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base1 = tmpdir / "baseline.1.annot.gz"
            base1.write_text("x\n", encoding="utf-8")

            with_at = resolve_chromosome_group(
                [str(tmpdir / "baseline.@.annot.gz")],
                chrom="1",
                label="annotation",
            )

            self.assertEqual(with_at, [str(base1)])

    def test_resolve_chromosome_group_rejects_bare_prefix_tokens(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            (tmpdir / "baseline.1.annot.gz").write_text("x\n", encoding="utf-8")

            with self.assertRaisesRegex(LDSCInputError, "annotation path for chromosome 1"):
                resolve_chromosome_group(
                    [str(tmpdir / "baseline.")],
                    chrom="1",
                    label="annotation",
                )

    def test_resolve_file_group_discovers_chromosome_suite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base1 = tmpdir / "baseline.1.annot.gz"
            base2 = tmpdir / "baseline.2.annot.gz"
            base1.write_text("x\n", encoding="utf-8")
            base2.write_text("x\n", encoding="utf-8")

            resolved = resolve_file_group(
                str(tmpdir / "baseline.@.annot.gz"),
                label="annotation",
                allow_chromosome_suite=True,
            )

            self.assertEqual(resolved, [str(base1), str(base2)])

    def test_resolve_chromosome_group_can_ignore_missing_chromosome(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self.assertEqual(
                resolve_chromosome_group(
                    [str(tmpdir / "baseline.@.annot.gz")],
                    chrom="1",
                    label="annotation",
                    required=False,
                ),
                [],
            )

    def test_resolve_plink_prefix_accepts_glob_and_returns_prefix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            prefix = tmpdir / "panel.1"
            (tmpdir / "panel.1.bed").write_text("", encoding="utf-8")
            (tmpdir / "panel.1.bim").write_text("", encoding="utf-8")
            (tmpdir / "panel.1.fam").write_text("", encoding="utf-8")

            resolved = resolve_plink_prefix(str(tmpdir / "panel.*"), chrom=None)
            self.assertEqual(resolved, str(prefix))

            (tmpdir / "panel.2.bed").write_text("", encoding="utf-8")
            (tmpdir / "panel.2.bim").write_text("", encoding="utf-8")
            (tmpdir / "panel.2.fam").write_text("", encoding="utf-8")
            with self.assertRaisesRegex(LDSCInputError, "matched 2 files"):
                resolve_plink_prefix(str(tmpdir / "panel.*"), chrom=None)

    def test_resolve_plink_prefix_group_discovers_chromosome_suite_from_plain_stem(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            for chrom in ("1", "2"):
                for suffix in (".bed", ".bim", ".fam"):
                    (tmpdir / f"panel_chr{chrom}{suffix}").write_bytes({
                        ".bed": b"\x6c\x1b\x01\x00",
                        ".bim": f"{chrom} rs{chrom} 0 10 A G\n".encode(),
                        ".fam": b"F I 0 0 1 -9\n",
                    }[suffix])

            stem = str(tmpdir / "panel_chr")

            self.assertEqual(
                resolve_plink_prefix_group(stem, allow_chromosome_suite=True),
                [str(tmpdir / "panel_chr1"), str(tmpdir / "panel_chr2")],
            )
            self.assertEqual(resolve_plink_prefix(stem, chrom="2"), str(tmpdir / "panel_chr2"))

    def test_resolve_plink_prefix_rejects_incomplete_file_trio(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "panel"
            prefix.with_suffix(".bed").write_text("", encoding="utf-8")
            prefix.with_suffix(".bim").write_text("", encoding="utf-8")

            with self.assertRaisesRegex(LDSCInputError, "complete .bed/.bim/.fam trio"):
                resolve_plink_prefix(str(prefix), chrom=None)

    def test_resolve_scalar_path_rejects_suffix_inference(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            path = tmpdir / "reference.1.parquet"
            path.write_text("x\n", encoding="utf-8")

            with self.assertRaisesRegex(LDSCInputError, str(tmpdir / "reference.1")):
                resolve_scalar_path(
                    str(tmpdir / "reference.1"),
                    label="parquet",
                )

    def test_resolve_scalar_path_missing_reports_cause_fix_and_reference(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            with self.assertRaisesRegex(
                LDSCInputError,
                "Could not resolve parquet path.*matched 0 files.*"
                "Most likely the path is misspelled.*"
                "docs/troubleshooting.md#common-input-path-did-not-resolve-to-one-file",
            ):
                resolve_scalar_path(str(tmpdir / "reference.1"), label="parquet")

    def test_resolve_file_group_requires_full_explicit_suite_pattern(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            gz = tmpdir / "baseline.1.annot.gz"
            gz.write_text("gz\n", encoding="utf-8")

            with self.assertRaisesRegex(LDSCInputError, "annotation path"):
                resolve_file_group(
                    str(tmpdir / "baseline.@"),
                    label="annotation",
                    allow_chromosome_suite=True,
                )

    def test_ensure_output_paths_available_allows_missing_paths(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "results" / "h2.tsv"

            ensure_output_paths_available([path], overwrite=False)

            self.assertFalse(path.exists())

    def test_ensure_output_paths_available_rejects_existing_paths(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "results.tsv"
            path.write_text("existing\n", encoding="utf-8")

            with self.assertRaisesRegex(FileExistsError, "overwrite"):
                ensure_output_paths_available([path], overwrite=False)

    def test_output_collision_message_reports_cause_fix_and_reference(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "results.tsv"
            path.write_text("existing\n", encoding="utf-8")

            with self.assertRaisesRegex(
                FileExistsError,
                "Cannot write output artifact.*already exists.*"
                "Most likely this output directory contains results from an earlier run.*"
                "docs/troubleshooting.md#common-output-artifact-already-exists",
            ):
                ensure_output_paths_available([path], overwrite=False)

    def test_ensure_output_paths_available_allows_existing_paths_with_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "results.tsv"
            path.write_text("existing\n", encoding="utf-8")

            ensure_output_paths_available([path], overwrite=True)

            self.assertEqual(path.read_text(encoding="utf-8"), "existing\n")

    def test_ensure_output_paths_available_error_lists_existing_artifacts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            first = Path(tmpdir) / "a.tsv"
            second = Path(tmpdir) / "b.tsv"
            first.write_text("a\n", encoding="utf-8")
            second.write_text("b\n", encoding="utf-8")

            with self.assertRaises(FileExistsError) as caught:
                ensure_output_paths_available([first, second], overwrite=False, label="summary table")

            message = str(caught.exception)
            self.assertIn(str(first), message)
            self.assertIn(str(second), message)
            self.assertIn("--overwrite", message)
            self.assertIn("overwrite=True", message)


def test_plain_plink_prefix_uses_bim_chromosomes(tmp_path):
    from ldsc.path_resolution import resolve_plink_prefix

    # The filename suggests chr1, but the data contain chr22.
    prefix = tmp_path / "1000G.EUR.QC.1"
    prefix.with_name(prefix.name + ".bed").write_bytes(b"\x6c\x1b\x01\x00")
    prefix.with_name(prefix.name + ".bim").write_text("22 rs1 0 10 A G\n")
    prefix.with_name(prefix.name + ".fam").write_text("F I 0 0 1 -9\n")
    assert resolve_plink_prefix(str(tmp_path / "1000G.EUR.QC."), chrom="22") == str(prefix)


def _write_plink_trio(prefix, chroms):
    Path(str(prefix) + ".bed").write_bytes(b"\x6c\x1b\x01" + b"\x00" * len(chroms))
    Path(str(prefix) + ".bim").write_text("".join(f"{chrom} rs{i} 0 {10+i} A G\n" for i, chrom in enumerate(chroms)))
    Path(str(prefix) + ".fam").write_text("F I 0 0 1 -9\n")


@pytest.mark.parametrize("style", ["1000G.EUR.QC.{chrom}", "panel_chr{chrom}"])
@pytest.mark.parametrize("form", ["plain", "glob", "at", "exact", "bed", "bim", "fam"])
def test_plink_forms_resolve_identical_chromosome_sources(tmp_path, style, form):
    from ldsc.path_resolution import inspect_plink_inputs

    expected = {}
    for chrom in map(str, range(1, 23)):
        prefix = tmp_path / style.format(chrom=chrom)
        _write_plink_trio(prefix, [chrom])
        expected[chrom] = str(prefix)
    if form in {"plain", "glob", "at"}:
        token = str(tmp_path / style.format(chrom={"plain": "", "glob": "*", "at": "@"}[form]))
    else:
        suffix = "" if form == "exact" else "." + form
        token = [prefix + suffix for prefix in expected.values()]
    result = inspect_plink_inputs(token, chromosomes=tuple(map(str, range(1, 23))), require_complete_suite=True)
    result.require_valid(required_chromosomes=tuple(expected))
    assert result.chromosome_prefixes == expected
    assert set(resolve_plink_prefix_group(token, allow_chromosome_suite=True)) == set(expected.values())


def test_plink_exact_prefix_precedes_stem_discovery(tmp_path):
    from ldsc.path_resolution import inspect_plink_inputs
    _write_plink_trio(tmp_path / "panel", ["22"])
    _write_plink_trio(tmp_path / "panel.extra", ["22"])
    result = inspect_plink_inputs(tmp_path / "panel")
    result.require_valid()
    assert result.chromosome_prefixes == {"22": str(tmp_path / "panel")}


def test_plink_joint_file_and_wrong_requested_chromosome(tmp_path):
    from ldsc.path_resolution import inspect_plink_inputs
    prefix = tmp_path / "arbitrary.17"
    _write_plink_trio(prefix, ["chr1", "22"])
    result = inspect_plink_inputs(prefix)
    result.require_valid()
    assert result.chromosome_prefixes == {"1": str(prefix), "22": str(prefix)}
    with pytest.raises(LDSCInputError, match="No validated PLINK trio contains chromosome 17"):
        resolve_plink_prefix(prefix, chrom="17")


def test_plink_reports_all_incomplete_and_conflicting_sources(tmp_path):
    from ldsc.path_resolution import inspect_plink_inputs
    for name in ("panel.a", "panel.b", "panel.c", "panel.d"):
        _write_plink_trio(tmp_path / name, ["22"])
    (tmp_path / "panel.c.bed").unlink()
    (tmp_path / "panel.d.fam").unlink()
    result = inspect_plink_inputs(tmp_path / "panel.")
    assert {row["reason"] for row in result.issues} == {"missing_required_input", "ambiguous_chromosome_input"}
    with pytest.raises(LDSCInputError) as caught:
        result.require_valid()
    assert all(name in str(caught.value) for name in ("panel.a", "panel.b", "panel.c.bed", "panel.d.fam"))
    with pytest.raises(LDSCInputError, match="panel.c.bed"):
        resolve_plink_prefix_group(tmp_path / "panel.", allow_chromosome_suite=True)


def test_plink_at_declares_content_and_workflow_completeness(tmp_path):
    from ldsc.path_resolution import inspect_plink_inputs
    _write_plink_trio(tmp_path / "panel.22", ["22"])
    token = tmp_path / "panel.@"
    inspect_plink_inputs(token).require_valid()
    full = inspect_plink_inputs(token, chromosomes=tuple(map(str, range(1, 23))), require_complete_suite=True)
    assert [row["chrom"] for row in full.issues] == list(map(str, range(1, 22)))
    _write_plink_trio(tmp_path / "panel.22", ["1"])
    with pytest.raises(LDSCInputError, match="contains chromosomes.*1"):
        inspect_plink_inputs(token).require_valid()
