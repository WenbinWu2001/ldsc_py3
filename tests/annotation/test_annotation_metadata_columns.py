"""Annotation files recognize metadata columns by exact name or alias only.

A column recognized as SNP, CM or MAF metadata is removed from the annotation
values, so a header that merely ends in a metadata name (``lowMAF``,
``eQTL_SNP``, ``EUR_MAF``) or in a short SNP alias (``Enhancers`` ends in
``rs``, ``Thyroid`` in ``id``) must stay an annotation column.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ldsc import cli
from ldsc._annotation.parsing import normalize_annotation_chunk
from ldsc._regression.quantile_storage import target_chunks
from ldsc.config import AnnotationBuildConfig, GlobalConfig
from ldsc._annotation.builder import AnnotationBuilder
from ldsc.errors import LDSCInputError
from tests.annotation_fixtures import build_annotation_fixture, fixture_values

PLINK = Path(__file__).resolve().parents[1] / "fixtures" / "minimal_external_resources" / "plink" / "hm3_chr22_subset"
BIM = pd.read_csv(f"{PLINK}.bim", sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "POS", "A1", "A2"])

# Annotation names that end in a metadata name or alias, with one value per row.
LOOKALIKES = {"lowMAF": [0, 1, 1], "eQTL_SNP": [1, 0, 1], "Enhancers": [1, 1, 0], "Thyroid": [0, 0, 1],
              "EUR_MAF": [0.1, 0.2, 0.3], "eQTL_CM": [2, 0, 0]}


def _frame(**identity):
    return pd.DataFrame({"CHR": ["1", "1", "1"], "POS": [100, 200, 300], **identity,
                         "base": [1, 1, 1], **LOOKALIKES})


@pytest.mark.parametrize("mode", ["chr_pos", "chr_pos_allele_aware"])
def test_coordinate_mode_file_without_snp_keeps_every_lookalike_as_an_annotation(mode):
    metadata, annotations = normalize_annotation_chunk(_frame(), "x.annot", mode)

    assert list(annotations.columns) == ["base", *LOOKALIKES]
    assert "SNP" not in metadata.columns
    for name, values in LOOKALIKES.items():
        np.testing.assert_array_equal(annotations[name].to_numpy(), np.asarray(values, dtype=np.float32))


@pytest.mark.parametrize("mode", ["rsid", "rsid_allele_aware", "chr_pos_allele_aware"])
def test_file_with_exact_snp_keeps_every_lookalike_as_an_annotation(mode):
    metadata, annotations = normalize_annotation_chunk(_frame(SNP=["rs1", "rs2", "rs3"]), "x.annot", mode)

    assert list(annotations.columns) == ["base", *LOOKALIKES]
    assert metadata["SNP"].tolist() == ["rs1", "rs2", "rs3"]


@pytest.mark.parametrize(("header", "snp"), [
    (["CHR", "POS", "SNP", "CM", "MAF"], "SNP"),
    (["chr", "bp", "rsid", "centimorgan", "maf"], "rsid"),
    (["#CHROM", "POSITION", "MarkerName", "CMBP", "MAF"], "MarkerName"),
    (["CHR", "POS", "ID", "cM", "MAF"], "ID"),
])
@pytest.mark.parametrize("mode", ["rsid", "chr_pos_allele_aware"])
def test_exact_names_and_aliases_are_still_metadata(header, snp, mode):
    frame = pd.DataFrame([["1", 100, "rs1", 0.5, 0.2, 1]], columns=[*header, "base"])

    metadata, annotations = normalize_annotation_chunk(frame, "x.annot", mode)

    assert list(annotations.columns) == ["base"]
    assert metadata["SNP"].tolist() == ["rs1"]
    assert metadata["POS"].tolist() == [100]


def test_rsid_mode_does_not_take_a_lookalike_as_the_required_snp_column():
    with pytest.raises(LDSCInputError, match="SNP identifier column"):
        normalize_annotation_chunk(_frame(), "x.annot", "rsid")


def test_annotation_builder_keeps_lookalike_annotation_columns(tmp_path):
    path = tmp_path / "baseline.annot"
    _frame().to_csv(path, sep="\t", index=False)
    builder = AnnotationBuilder(GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"), AnnotationBuildConfig())

    bundle = build_annotation_fixture(builder, AnnotationBuildConfig(baseline_annot_sources=(str(path),)))

    assert bundle.baseline_columns == ["base", *LOOKALIKES]
    np.testing.assert_array_equal(fixture_values(bundle, "baseline")["lowMAF"].to_numpy(), [0, 1, 1])


def test_quantile_target_annotation_does_not_take_a_lookalike_as_the_snp_column(tmp_path):
    lookalike = tmp_path / "lookalike.tsv"
    lookalike.write_text("CHR\tPOS\tThyroid\ttarget\n1\t10\t1\t0.5\n", encoding="utf-8")
    alias = tmp_path / "alias.tsv"
    alias.write_text("CHR\tPOS\tRSID\ttarget\n1\t10\trs1\t0.5\n", encoding="utf-8")

    with pytest.raises(LDSCInputError, match="SNP identifier column"):
        next(target_chunks([str(lookalike)], "target", None))
    chunk = next(target_chunks([str(alias)], "target", None))
    assert chunk["SNP"].tolist() == ["rs1"]
    assert chunk["target_value"].tolist() == [0.5]


def _ldscore(output_dir, annotation, mode):
    assert cli.run_cli(["ldscore", "--output-dir", str(output_dir), "--snp-identifier", mode,
                        "--genome-build", "hg38", "--ld-wind-snps", "3", "--plink-prefix", str(PLINK),
                        "--baseline-annot-sources", str(annotation)]) == 0
    return pd.read_parquet(output_dir / "ldscore.baseline.parquet")


@pytest.mark.parametrize("mode", ["chr_pos_allele_aware", "chr_pos"])
def test_ld_scores_include_lookalike_annotations_whatever_their_names(tmp_path, mode):
    """LD scores do not depend on annotation names; ``Thyroid`` is twice ``base``, so its LD score is too."""
    rng = np.random.default_rng(20261106)
    values = {"base": np.ones(len(BIM)), "lowMAF": rng.integers(0, 2, len(BIM)),
              "eQTL_SNP": rng.integers(0, 2, len(BIM)), "Enhancers": rng.integers(0, 2, len(BIM)),
              "Thyroid": np.full(len(BIM), 2.0), "EUR_MAF": rng.uniform(0, 0.5, len(BIM))}
    lookalike = BIM[["CHR", "POS"]].assign(**values)
    neutral_names = {name: f"c{index}" for index, name in enumerate(values)}
    lookalike_path, neutral_path = tmp_path / "lookalike.annot", tmp_path / "neutral.annot"
    lookalike.to_csv(lookalike_path, sep="\t", index=False)
    lookalike.rename(columns=neutral_names).to_csv(neutral_path, sep="\t", index=False)

    observed = _ldscore(tmp_path / "lookalike", lookalike_path, mode)
    reference = _ldscore(tmp_path / "neutral", neutral_path, mode)

    for name in values:
        np.testing.assert_array_equal(observed[name].to_numpy(), reference[neutral_names[name]].to_numpy(), err_msg=name)
    np.testing.assert_allclose(observed["Thyroid"].to_numpy(), 2 * observed["base"].to_numpy(), rtol=1e-6)
