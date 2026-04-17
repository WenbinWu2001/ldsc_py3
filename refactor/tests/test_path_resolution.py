import os
from pathlib import Path
import sys
import tempfile
import unittest


SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc.path_resolution import (
    ANNOTATION_SUFFIXES,
    PARQUET_SUFFIXES,
    normalize_path_token,
    normalize_path_tokens,
    resolve_chromosome_group,
    resolve_file_group,
    resolve_plink_prefix,
    resolve_scalar_path,
)


class PathResolutionTest(unittest.TestCase):
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
            with self.assertRaisesRegex(ValueError, "exactly one"):
                resolve_scalar_path(str(tmpdir / "trait.*.tsv"), label="sumstats")

    def test_resolve_file_group_handles_globs_prefixes_and_deduplicates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base1 = tmpdir / "baseline.1.annot.gz"
            base2 = tmpdir / "baseline.2.annot.gz"
            base1.write_text("x\n", encoding="utf-8")
            base2.write_text("x\n", encoding="utf-8")

            resolved = resolve_file_group(
                [str(tmpdir / "baseline.*.annot.gz"), str(tmpdir / "baseline.1")],
                suffixes=ANNOTATION_SUFFIXES,
                label="annotation",
            )

            self.assertEqual(resolved, [str(base1), str(base2)])

    def test_resolve_chromosome_group_supports_at_and_bare_prefix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base1 = tmpdir / "baseline.1.annot.gz"
            base1.write_text("x\n", encoding="utf-8")

            with_at = resolve_chromosome_group(
                [str(tmpdir / "baseline.@")],
                chrom="1",
                suffixes=ANNOTATION_SUFFIXES,
                label="annotation",
            )
            with_bare_prefix = resolve_chromosome_group(
                [str(tmpdir / "baseline.")],
                chrom="1",
                suffixes=ANNOTATION_SUFFIXES,
                label="annotation",
            )

            self.assertEqual(with_at, [str(base1)])
            self.assertEqual(with_bare_prefix, [str(base1)])

    def test_resolve_file_group_discovers_chromosome_suite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base1 = tmpdir / "baseline.1.annot.gz"
            base2 = tmpdir / "baseline.2.annot.gz"
            base1.write_text("x\n", encoding="utf-8")
            base2.write_text("x\n", encoding="utf-8")

            resolved = resolve_file_group(
                str(tmpdir / "baseline.@"),
                suffixes=ANNOTATION_SUFFIXES,
                label="annotation",
                allow_chromosome_suite=True,
            )

            self.assertEqual(resolved, [str(base1), str(base2)])

    def test_resolve_chromosome_group_can_ignore_missing_chromosome(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self.assertEqual(
                resolve_chromosome_group(
                    [str(tmpdir / "baseline.@")],
                    chrom="1",
                    suffixes=ANNOTATION_SUFFIXES,
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
            with self.assertRaisesRegex(ValueError, "exactly one"):
                resolve_plink_prefix(str(tmpdir / "panel.*"), chrom=None)

    def test_resolve_scalar_path_supports_suffix_inference(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            path = tmpdir / "reference.1.parquet"
            path.write_text("x\n", encoding="utf-8")

            resolved = resolve_scalar_path(
                str(tmpdir / "reference.1"),
                suffixes=PARQUET_SUFFIXES,
                label="parquet",
            )
            self.assertEqual(resolved, str(path))

    def test_resolve_file_group_prefers_compressed_suffixes_over_plain_text(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            plain = tmpdir / "baseline.1.annot"
            gz = tmpdir / "baseline.1.annot.gz"
            plain.write_text("plain\n", encoding="utf-8")
            gz.write_text("gz\n", encoding="utf-8")

            resolved = resolve_file_group(
                str(tmpdir / "baseline.@"),
                suffixes=ANNOTATION_SUFFIXES,
                label="annotation",
                allow_chromosome_suite=True,
            )

            self.assertEqual(resolved, [str(gz)])
