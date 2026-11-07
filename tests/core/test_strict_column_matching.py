"""Every input header matches a field only by its canonical name or a listed alias.

Headers are compared after lowercasing and removing punctuation, so case and
punctuation never matter, but a header that merely ends in a field name or an
alias (``EUR_MAF``, ``eQTL_POS``, ``my_hg19_pos``) matches nothing. When a
required field is missing, the error names such look-alike headers as a hint.
"""
from __future__ import annotations

import pandas as pd
import pytest

from ldsc._annotation.parsing import normalize_annotation_chunk
from ldsc._core.columns import (
    A1_COLUMN_SPEC,
    A2_COLUMN_SPEC,
    CHR_COLUMN_SPEC,
    CM_COLUMN_SPEC,
    MAF_COLUMN_SPEC,
    POS_COLUMN_SPEC,
    RESTRICTION_HG19_POS_SPEC,
    RESTRICTION_HG38_POS_SPEC,
    SNP_COLUMN_SPEC,
    resolve_optional_column,
    resolve_required_column,
    resolve_restriction_chr_pos_columns,
    resolve_restriction_rsid_column,
)
from ldsc._genome.hm3 import _load_hm3_curated_map_from_path
from ldsc._refpanel.genetic_map import GENETIC_MAP_CM_SPEC, load_genetic_map
from ldsc._regression.quantile_storage import _reference_columns
from ldsc.errors import LDSCInputError

SHARED_SPECS = (CHR_COLUMN_SPEC, POS_COLUMN_SPEC, SNP_COLUMN_SPEC, CM_COLUMN_SPEC, MAF_COLUMN_SPEC,
                A1_COLUMN_SPEC, A2_COLUMN_SPEC, RESTRICTION_HG19_POS_SPEC, RESTRICTION_HG38_POS_SPEC,
                GENETIC_MAP_CM_SPEC)


@pytest.mark.parametrize(("spec", "header"), [
    (CHR_COLUMN_SPEC, "eQTL_CHR"),
    (POS_COLUMN_SPEC, "eQTL_POS"),
    (POS_COLUMN_SPEC, "start_bp"),
    (POS_COLUMN_SPEC, "hg19_POS"),
    (SNP_COLUMN_SPEC, "eQTL_SNP"),
    (SNP_COLUMN_SPEC, "Enhancers"),
    (SNP_COLUMN_SPEC, "Thyroid"),
    (CM_COLUMN_SPEC, "eQTL_CM"),
    (MAF_COLUMN_SPEC, "EUR_MAF"),
    (A1_COLUMN_SPEC, "EUR_A1"),
    (A2_COLUMN_SPEC, "EUR_A2"),
    (RESTRICTION_HG19_POS_SPEC, "my_hg19_pos"),
    (RESTRICTION_HG38_POS_SPEC, "liftover_GRCh38_BP"),
    (GENETIC_MAP_CM_SPEC, "eQTL_CM"),
])
def test_a_header_that_only_ends_in_a_name_matches_no_field(spec, header):
    assert resolve_optional_column(["unrelated", header], spec) is None


@pytest.mark.parametrize("spec", SHARED_SPECS, ids=lambda spec: spec.label)
def test_every_listed_name_matches_in_any_case_and_punctuation(spec):
    for name in (spec.canonical, *spec.aliases):
        for header in (name, name.lower(), name.replace("_", "-")):
            assert resolve_optional_column(["unrelated", header], spec) == header


def test_build_specific_position_matches_after_normalization():
    assert resolve_optional_column(["GRCh37-BP"], RESTRICTION_HG19_POS_SPEC) == "GRCh37-BP"
    assert resolve_optional_column(["grch38.position"], RESTRICTION_HG38_POS_SPEC) == "grch38.position"


def test_missing_required_column_names_lookalike_headers():
    with pytest.raises(LDSCInputError, match=r"EUR_MAF.*rename it to MAF") as caught:
        resolve_required_column(["CHR", "POS", "SNP", "EUR_MAF", "lowMAF", "notes"], MAF_COLUMN_SPEC, context="ref.tsv")
    hint = str(caught.value).split("The header has", 1)[1]
    # Only a name after a separator counts: lowMAF ends in MAF inside one word.
    assert "lowMAF" not in hint and "notes" not in hint


def test_missing_required_column_links_its_troubleshooting_section():
    with pytest.raises(LDSCInputError) as caught:
        resolve_required_column(["CHR", "POS"], MAF_COLUMN_SPEC)
    assert str(caught.value).endswith(
        "Other causes & fixes: https://wenbinwu2001.github.io/ldsc_py3/user-guide/troubleshooting.html#common-a-required-column-is-not-recognized"
    )


def test_missing_required_column_without_lookalikes_adds_no_hint():
    with pytest.raises(LDSCInputError) as caught:
        resolve_required_column(["CHR", "POS", "Thyroid"], SNP_COLUMN_SPEC)
    assert "rename it to" not in str(caught.value)


def test_missing_required_column_with_no_similar_header_adds_no_hint():
    with pytest.raises(LDSCInputError) as caught:
        resolve_required_column(["CHR", "POS"], MAF_COLUMN_SPEC)
    assert "rename it to" not in str(caught.value)


def test_annotation_file_needs_a_named_position_column():
    frame = pd.DataFrame({"CHR": ["1"], "eQTL_POS": [100], "base": [1.0]})
    with pytest.raises(LDSCInputError, match=r"eQTL_POS.*rename it to POS"):
        normalize_annotation_chunk(frame, "x.annot", "chr_pos")


def test_restriction_file_needs_a_named_snp_column():
    with pytest.raises(LDSCInputError, match=r"eQTL_SNP.*rename it to SNP"):
        resolve_restriction_rsid_column(["eQTL_SNP", "score"], context="keep.tsv")


def test_restriction_file_prefers_the_runs_build_column_over_generic_position():
    assert resolve_restriction_chr_pos_columns(["CHR", "POS", "hg19_POS"], genome_build="hg19") == ("CHR", "hg19_POS")
    assert resolve_restriction_chr_pos_columns(["CHR", "POS", "hg19_POS"], genome_build="hg38") == ("CHR", "POS")


def test_restriction_file_with_only_another_builds_position_stops():
    with pytest.raises(LDSCInputError, match="missing recognizable POS") as caught:
        resolve_restriction_chr_pos_columns(["CHR", "hg19_POS"], genome_build="hg38", context="keep.tsv")
    # Another build's positions must never be suggested as the run's POS.
    assert "rename it to POS" not in str(caught.value)


def test_restriction_file_ignores_a_position_that_only_ends_in_a_name():
    with pytest.raises(LDSCInputError, match=r"eQTL_POS.*rename it to POS"):
        resolve_restriction_chr_pos_columns(["CHR", "eQTL_POS"], genome_build="hg19", context="keep.tsv")


def test_user_hm3_map_needs_named_columns(tmp_path):
    path = tmp_path / "map.tsv"
    path.write_text("CHR\thg19_POS\thg38_POS\tmy_SNP\tA1\tA2\n1\t100\t200\trs1\tA\tG\n", encoding="utf-8")
    with pytest.raises(LDSCInputError, match=r"my_SNP.*rename it to SNP"):
        _load_hm3_curated_map_from_path(path, preserve_extra_columns=False)


def test_genetic_map_needs_a_named_position_column(tmp_path):
    path = tmp_path / "map.txt"
    path.write_text("chr eQTL_position Genetic_Map(cM)\n1 100 0.0\n", encoding="utf-8")
    with pytest.raises(LDSCInputError, match=r"eQTL_position.*rename it to POS"):
        load_genetic_map(path, command="ldscore")


def test_genetic_map_reads_the_hapmap_header(tmp_path):
    path = tmp_path / "genetic_map_chr1.txt"
    path.write_text("Chromosome Position(bp) Rate(cM/Mb) Map(cM)\n1 100 0.5 0.0\n1 200 0.5 0.1\n", encoding="utf-8")
    loaded = load_genetic_map(path, command="ldscore")
    assert loaded["POS"].tolist() == [100, 200]
    assert loaded["CM"].tolist() == [0.0, 0.1]


def test_quantile_reference_metadata_does_not_read_eur_maf_as_maf():
    with pytest.raises(LDSCInputError, match=r"missing required column\(s\) MAF.*EUR_MAF.*rename it to MAF"):
        _reference_columns(["CHR", "POS", "SNP", "EUR_MAF"], "ref.tsv")
