"""A keep file is a set of individuals: a repeated IID selects its ``.fam`` row once.

``--keep-indivs-file`` selects the ``.fam`` rows whose IID it lists, in ``.fam``
order. A keep file that lists an IID more than once, as concatenated keep lists
do, must select the same individuals, and give the same LD scores, as the file
that lists each IID once, in every PLINK route.

The expected rows come from the ``.fam`` text itself. The fixture is the
32-SNP chromosome-22 panel with 3,202 individuals; the keep list names fifteen
individuals near the start of the ``.fam`` and the last ten, so a repeated IID
can push a selected row past the last individual.

IIDs are text: a numeric-looking IID matches only the same text, so ``01003``
does not match ``1003``, a numeric keep list matches a ``.fam`` that also
holds non-numeric IIDs, and ``NA`` or ``nan`` is an IID, not a missing value.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import pyarrow.parquet as pq

from ldsc.cli import run_cli
from ldsc._core.plink_files import PlinkFAMFile
from ldsc.errors import LDSCInputError
from ldsc._refpanel.plink_source import keep_individual_indices, read_keep_individuals
from tests.ref_panel_helpers import HM3_CHR22_PLINK, read_bim

FIXTURE = HM3_CHR22_PLINK
FAM_IIDS = [line.split()[1] for line in Path(f"{FIXTURE}.fam").read_text(encoding="utf-8").splitlines()]
LISTED = FAM_IIDS[1:16] + FAM_IIDS[-10:]
REPEATED = {
    "early_individual": [FAM_IIDS[1]],
    "last_fam_individual": [FAM_IIDS[-1]],
    "both": [FAM_IIDS[1], FAM_IIDS[-1]],
}
needs_bitarray = pytest.mark.skipif(importlib.util.find_spec("bitarray") is None, reason="bitarray required")


def write_keep(path: Path, iids: list[str]) -> Path:
    path.write_text("\n".join(iids) + "\n", encoding="utf-8")
    return path


def keep_files(tmp_path: Path, repeats: dict[str, list[str]]) -> dict[str, Path]:
    """Write the keep file that lists each IID once and one file per repeat pattern."""
    files = {"once": write_keep(tmp_path / "once.keep", LISTED)}
    for name, extra in repeats.items():
        files[name] = write_keep(tmp_path / f"{name}.keep", LISTED + extra)
    return files


@pytest.mark.parametrize("repeated", REPEATED.values(), ids=REPEATED.keys())
def test_a_repeated_iid_selects_its_fam_row_once(tmp_path, repeated):
    listed = set(LISTED)
    expected = [row for row, iid in enumerate(FAM_IIDS) if iid in listed]
    assert expected == list(range(1, 16)) + list(range(len(FAM_IIDS) - 10, len(FAM_IIDS)))
    fam = PlinkFAMFile(f"{FIXTURE}.fam")

    selected = keep_individual_indices(read_keep_individuals(str(write_keep(tmp_path / "keep.txt", LISTED + repeated))), fam)

    assert selected == expected


SILENT_REPEATS = {name: REPEATED[name] for name in ("early_individual", "last_fam_individual")}


@needs_bitarray
def test_ldscore_with_a_repeated_iid_matches_the_keep_file_listing_it_once(tmp_path):
    tables = {}
    for name, keep in keep_files(tmp_path, SILENT_REPEATS).items():
        out = tmp_path / f"ldscore_{name}"
        assert run_cli(["ldscore", "--plink-prefix", str(FIXTURE), "--genome-build", "hg38", "--ld-wind-kb", "1000",
                        "--keep-indivs-file", str(keep), "--output-dir", str(out)]) == 0
        tables[name] = pd.read_parquet(out / "ldscore.baseline.parquet")

    for name in SILENT_REPEATS:
        pd.testing.assert_frame_equal(tables[name], tables["once"], check_exact=True, obj=name)


@needs_bitarray
def test_build_r2_matrix_with_a_repeated_iid_matches_the_keep_file_listing_it_once(tmp_path):
    builds = {}
    for name, keep in keep_files(tmp_path, SILENT_REPEATS).items():
        out = tmp_path / f"r2_{name}"
        assert run_cli(["build-r2-matrix", "--plink-prefix", str(FIXTURE), "--source-genome-build", "hg38",
                        "--ld-wind-kb", "1000", "--keep-indivs-file", str(keep), "--output-dir", str(out)]) == 0
        [r2_path] = out.rglob("chr22_r2.parquet")
        [meta_path] = out.rglob("chr22_meta.tsv.gz")
        n_samples = pq.ParquetFile(r2_path).schema_arrow.metadata[b"ldsc:n_samples"]
        builds[name] = (pd.read_parquet(r2_path), pd.read_csv(meta_path, sep="\t", comment="#"), n_samples)

    assert builds["once"][2] == str(len(LISTED)).encode()
    for name in SILENT_REPEATS:
        pd.testing.assert_frame_equal(builds[name][0], builds["once"][0], check_exact=True, obj=f"{name} R2")
        pd.testing.assert_frame_equal(builds[name][1], builds["once"][1], check_exact=True, obj=f"{name} sidecar")
        assert builds[name][2] == builds["once"][2]


@needs_bitarray
def test_gene_index_with_a_repeated_iid_matches_the_keep_file_listing_it_once(tmp_path, monkeypatch):
    from ldsc._ldscore.gene_index.build import main as build_index
    from tests.gene_index_helpers import narrow_index_chromosomes

    narrow_index_chromosomes(monkeypatch, ["22"])
    bim = read_bim(FIXTURE)
    bim.loc[:, ["CHR", "SNP", "POS", "CM", "A1", "A2"]].assign(base=1.0).to_csv(
        tmp_path / "baseline.annot", sep="\t", index=False)
    bim.loc[:, ["SNP", "CHR", "POS"]].to_csv(tmp_path / "regression_snps.tsv", sep="\t", index=False)
    (tmp_path / "genes.tsv").write_text(
        "gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\n"
        f"G1\tGene1\t22\t{int(bim.POS.min())}\t{int(bim.POS.median())}\thg38\n")
    indexes = {}
    for name, keep in keep_files(tmp_path, SILENT_REPEATS).items():
        out = tmp_path / f"index_{name}"
        build_index(["--baseline-annot-sources", str(tmp_path / "baseline.annot"), "--plink-prefix", str(FIXTURE),
                     "--gene-coordinate-file", str(tmp_path / "genes.tsv"), "--genome-build", "hg38",
                     "--snp-identifier", "chr_pos_allele_aware", "--ld-wind-kb", "1000", "--padding-bp", "0", "--gene-exclude-regions", "none",
                     "--regr-snps-exclude-regions", "none", "--regr-snps-file", str(tmp_path / "regression_snps.tsv"),
                     "--keep-indivs-file", str(keep), "--output-dir", str(out)])
        metadata = json.loads((out / "metadata.json").read_text(encoding="utf-8"))
        with np.load(out / "chromosomes" / "chr22" / "ldscore_operator.npz") as operator:
            arrays = {key: operator[key] for key in operator.files}
        indexes[name] = (metadata, pd.read_parquet(out / "chromosomes" / "chr22" / "baseline_rows.parquet"), arrays)

    assert indexes["once"][0]["index_identity"]["selected_individuals"]["selected_count"] == len(LISTED)
    for name in SILENT_REPEATS:
        assert indexes[name][0]["index_id"] == indexes["once"][0]["index_id"]
        pd.testing.assert_frame_equal(indexes[name][1], indexes["once"][1], check_exact=True, obj=name)
        for key, values in indexes["once"][2].items():
            np.testing.assert_array_equal(indexes[name][2][key], values, err_msg=f"{name} {key}")


def write_fam(path: Path, iids: list[str]) -> Path:
    path.write_text("".join(f"F{row} {iid} 0 0 0 -9\n" for row, iid in enumerate(iids)), encoding="utf-8")
    return path


NUMERIC_LOOKING = {
    # FAM IIDs, keep-file IIDs. The .fam and keep file are text: an ID matches only its exact text.
    "mixed_fam_numeric_keep": (["1001", "A2", "1003"], ["1001", "1003"]),
    "leading_zero": (["1003", "01003", "1004"], ["01003"]),
    "numeric_fam_mixed_keep": (["1001", "1002"], ["1002", "B7"]),
}


@pytest.mark.parametrize("fam_iids, keep_iids", NUMERIC_LOOKING.values(), ids=NUMERIC_LOOKING.keys())
def test_numeric_looking_iids_match_by_their_exact_text(tmp_path, fam_iids, keep_iids):
    expected = [row for row, iid in enumerate(fam_iids) if iid in set(keep_iids)]
    fam = PlinkFAMFile(str(write_fam(tmp_path / "panel.fam", fam_iids)))

    selected = keep_individual_indices(read_keep_individuals(str(write_keep(tmp_path / "keep.txt", keep_iids))), fam)

    assert selected == expected


def test_a_keep_iid_naming_several_fam_rows_is_rejected(tmp_path):
    # Two .fam rows share IID S1 (different FIDs); matching by IID cannot tell them apart.
    fam_path = write_fam(tmp_path / "panel.fam", ["S1", "S1", "S2", "S3", "S3"])
    fam = PlinkFAMFile(str(fam_path))
    keep = read_keep_individuals(str(write_keep(tmp_path / "keep.txt", ["S1", "S2", "S3"])))

    with pytest.raises(LDSCInputError, match="more than one row") as raised:
        keep_individual_indices(keep, fam)

    message = str(raised.value)
    assert "2 listed IIDs name" in message and "S1, S3" in message and str(fam_path) in message


def test_a_repeated_fam_iid_the_keep_file_does_not_list_is_accepted(tmp_path):
    fam = PlinkFAMFile(str(write_fam(tmp_path / "panel.fam", ["S1", "S1", "S2"])))

    assert keep_individual_indices(read_keep_individuals(str(write_keep(tmp_path / "keep.txt", ["S2"]))), fam) == [2]


NA_TOKENS = {
    # FAM IIDs, keep-file IIDs. Pandas reads these words as missing values by default;
    # as IIDs they are ordinary text and match only the same text.
    "na_listed": (["S0", "NA", "nan", "NULL"], ["NA"]),
    "na_and_null_listed": (["S0", "NA", "nan", "NULL"], ["NULL", "NA", "S0"]),
    "na_not_listed": (["S0", "NA", "nan"], ["S0"]),
}


@pytest.mark.parametrize("fam_iids, keep_iids", NA_TOKENS.values(), ids=NA_TOKENS.keys())
def test_missing_value_words_are_iids_matched_by_their_exact_text(tmp_path, fam_iids, keep_iids):
    expected = [row for row, iid in enumerate(fam_iids) if iid in set(keep_iids)]
    fam = PlinkFAMFile(str(write_fam(tmp_path / "panel.fam", fam_iids)))

    selected = keep_individual_indices(read_keep_individuals(str(write_keep(tmp_path / "keep.txt", keep_iids))), fam)

    assert selected == expected
    assert fam.IDList.iloc[:, 0].tolist() == fam_iids
