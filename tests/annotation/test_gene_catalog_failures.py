"""A structurally invalid gene coordinate catalog names its defects on every route.

``annotate``, ``ldscore`` and ``quantile-h2``, on the command line and through
their output-directory Python functions, write the complete defect table to
``diagnostics/gene_catalog_issues.tsv`` of their output directory; the
in-memory ``AnnotationBuilder`` writes no table and relies on the message and
the exception's ``issues`` frame.
"""

from pathlib import Path
import pickle

import pandas as pd
import pytest

import ldsc
from ldsc import cli
from ldsc._annotation.builder import AnnotationBuilder
from ldsc._annotation.gene_lists import CATALOG_ISSUE_COLUMNS, GeneCatalog, GeneCatalogValidationError
from ldsc.config import AnnotationBuildConfig, GlobalConfig, reset_global_config, set_global_config
from tests.ldscore.test_gene_query_integration import inputs  # noqa: F401  (pytest fixture)

ANCHOR = "docs/user-guide/troubleshooting.md#common-gene-coordinate-catalog-is-not-valid"
HEADER = "gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\n"


def _unsupported_build_catalog(path: Path) -> Path:
    path.write_text(HEADER + "G1\tONE\t1\t5\t6\tGRCh-38\nG2\tTWO\t1\t6\t8\tGRCh-38\n")
    return path


def _assert_names_the_build_defect(message: str, source: str) -> None:
    assert f"Gene coordinate catalog '{source}' is not structurally valid: found 2 defect(s)" in message
    assert "line 2, genome_build: unsupported_genome_build, value 'GRCh-38'" in message
    assert "line 3, genome_build: unsupported_genome_build, value 'GRCh-38'" in message
    assert "Repair: Use hg19/GRCh37 or hg38/GRCh38 coordinates" in message
    assert ANCHOR in message


def test_message_names_each_defect_with_its_line_value_and_repair(tmp_path):
    catalog = _unsupported_build_catalog(tmp_path / "catalog.tsv")

    with pytest.raises(GeneCatalogValidationError) as caught:
        GeneCatalog.load(catalog)

    _assert_names_the_build_defect(str(caught.value), "catalog.tsv")
    assert caught.value.issues["reason"].tolist() == ["unsupported_genome_build"] * 2


def test_message_shows_the_first_defects_and_counts_the_rest(tmp_path):
    catalog = tmp_path / "catalog.tsv"
    catalog.write_text(HEADER + "".join(f"G{i}\tN{i}\t1\t5\n" for i in range(8)))

    with pytest.raises(GeneCatalogValidationError) as caught:
        GeneCatalog.load(catalog)

    message = str(caught.value)
    issues = caught.value.issues
    assert f"found {len(issues)} defect(s)" in message
    assert "line 2, end: invalid_end. Repair:" in message
    assert "line 2, row: malformed_catalog_row (expected 6 fields, found 4), value 'G0\\tN0\\t1\\t5'" in message
    assert f"{len(issues) - 5} more defect(s) not shown" in message
    assert "line 4," not in message
    assert ANCHOR in message


def test_message_states_a_repair_shared_by_consecutive_defects_once(tmp_path):
    catalog = tmp_path / "catalog.tsv"
    catalog.write_text(HEADER + "".join(f"G{i}\tN{i}\t1\t{10 * i + 5}\t{10 * i + 6}\tGRCh-38\n" for i in range(7)))

    with pytest.raises(GeneCatalogValidationError) as caught:
        GeneCatalog.load(catalog)

    message = str(caught.value)
    assert ("line 5, genome_build: unsupported_genome_build, value 'GRCh-38'; "
            "line 6, genome_build: unsupported_genome_build, value 'GRCh-38'. "
            "Repair: Use hg19/GRCh37 or hg38/GRCh38 coordinates; 2 more defect(s) not shown") in message
    assert message.count("Repair:") == 1


@pytest.mark.parametrize(("name", "payload", "reason"), [
    ("catalog.tsv.gz", b"not a gzip stream", "catalog_unparseable"),
    ("catalog.tsv", b"gene_id\xff\n", "catalog_unparseable"),
    ("catalog.tsv", None, "catalog_unreadable"),
], ids=["gzip", "utf8", "unreadable"])
def test_file_level_failures_point_to_the_troubleshooting_section(tmp_path, name, payload, reason):
    catalog = tmp_path / name
    if payload is None:
        catalog.mkdir()  # reading a directory fails with an OSError
    else:
        catalog.write_bytes(payload)

    with pytest.raises(GeneCatalogValidationError) as caught:
        GeneCatalog.load(catalog)

    assert f"gene coordinate catalog '{name}'" in str(caught.value).lower()
    assert str(caught.value).endswith(f"Other causes & fixes: {ANCHOR}")
    assert caught.value.issues["reason"].iloc[0] == reason


def test_error_keeps_its_message_and_table_across_a_process_boundary(tmp_path):
    # quantile-h2 loads the catalog in a worker process and re-raises the error in the parent.
    catalog = _unsupported_build_catalog(tmp_path / "catalog.tsv")
    with pytest.raises(GeneCatalogValidationError) as caught:
        GeneCatalog.load(catalog)

    # The worker adds a note with its traceback before the error crosses back.
    caught.value.add_note("Calculation child traceback: ...")

    restored = pickle.loads(pickle.dumps(caught.value))

    assert str(restored) == str(caught.value)
    pd.testing.assert_frame_equal(restored.issues, caught.value.issues)
    assert restored.__notes__ == ["Calculation child traceback: ..."]


def test_embedded_catalog_message_names_its_defects():
    frame = pd.DataFrame({"gene_id": ["G1"], "gene_name": ["ONE"], "chrom": ["1"], "start": [9], "end": [5],
                          "genome_build": ["hg19"], "catalog_line": [2]})

    with pytest.raises(GeneCatalogValidationError) as caught:
        GeneCatalog.from_embedded_frame(frame)

    message = str(caught.value)
    assert "is not canonical: found 1 defect(s): line 2" in message
    assert caught.value.issues["reason"].iloc[0] in message


def test_ldscore_cli_names_the_defects_and_writes_the_catalog_table(tmp_path, inputs):  # noqa: F811
    prefix, baseline, regression, _ = inputs
    catalog = _unsupported_build_catalog(tmp_path / "catalog.tsv")
    focal = tmp_path / "focal.txt"
    focal.write_text("G1\n")
    output = tmp_path / "out"
    argv = ["ldscore", "--output-dir", str(output), "--baseline-annot-sources", str(baseline),
            "--query-annot-gene-list-sources", str(focal), "--gene-coordinate-file", str(catalog),
            "--padding-bp", "0", "--plink-prefix", str(prefix), "--genome-build", "hg19",
            "--snp-identifier", "rsid", "--ld-wind-cm", "0.2", "--regr-snps-file", str(regression),
            "--regr-snps-exclude-regions", "none"]

    with pytest.raises(GeneCatalogValidationError) as caught:
        cli.main(argv)

    _assert_names_the_build_defect(str(caught.value), "catalog.tsv")
    table = pd.read_csv(output / "diagnostics/gene_catalog_issues.tsv", sep="\t", keep_default_na=False)
    assert tuple(table.columns) == CATALOG_ISSUE_COLUMNS
    assert table["reason"].tolist() == ["unsupported_genome_build"] * 2
    assert table["catalog_line"].tolist() == [2, 3]
    assert not (output / "metadata.json").exists()

    # The table belongs to the LD-score output directory: a corrected overwrite removes it.
    catalog.write_text(catalog.read_text().replace("GRCh-38", "hg19"))
    result = cli.main([*argv, "--overwrite"])
    assert result.query_columns == ["focal"]
    assert not (output / "diagnostics/gene_catalog_issues.tsv").exists()
    assert "gene_catalog_issues" not in result.output_paths


def test_run_ldscore_names_the_defects_and_writes_the_catalog_table(tmp_path, inputs):  # noqa: F811
    # The Python form of ldscore shares the command's output directory and its diagnostics.
    prefix, baseline, regression, _ = inputs
    catalog = _unsupported_build_catalog(tmp_path / "catalog.tsv")
    focal = tmp_path / "focal.txt"
    focal.write_text("G1\n")
    output = tmp_path / "py"
    set_global_config(GlobalConfig(snp_identifier="rsid"))
    try:
        with pytest.raises(GeneCatalogValidationError) as caught:
            ldsc.run_ldscore(output_dir=output, baseline_annot_sources=[str(baseline)],
                             query_annot_gene_list_sources=[str(focal)], gene_coordinate_file=str(catalog),
                             padding_bp=0, plink_prefix=str(prefix), ld_wind_cm=0.2,
                             regr_snps_file=str(regression), regr_snps_exclude_regions="none")
    finally:
        reset_global_config()

    _assert_names_the_build_defect(str(caught.value), "catalog.tsv")
    table = pd.read_csv(output / "diagnostics/gene_catalog_issues.tsv", sep="\t", keep_default_na=False)
    assert tuple(table.columns) == CATALOG_ISSUE_COLUMNS
    assert table["catalog_line"].tolist() == [2, 3]
    assert not (output / "metadata.json").exists()


def test_annotation_builder_names_the_defects_and_writes_no_catalog_table(tmp_path, inputs):  # noqa: F811
    _, baseline, _, _ = inputs
    catalog = _unsupported_build_catalog(tmp_path / "catalog.tsv")
    focal = tmp_path / "focal.txt"
    focal.write_text("G1\n")
    spec = AnnotationBuildConfig(baseline_annot_sources=(str(baseline),), query_annot_gene_list_sources=(str(focal),),
                                 gene_coordinate_file=str(catalog))

    with pytest.raises(GeneCatalogValidationError) as caught:
        AnnotationBuilder(GlobalConfig(snp_identifier="rsid")).run(spec, output_dir=tmp_path / "work")

    _assert_names_the_build_defect(str(caught.value), "catalog.tsv")
    assert not (tmp_path / "work/diagnostics/gene_catalog_issues.tsv").exists()


def test_annotate_cli_names_the_defects_and_writes_the_catalog_table(tmp_path, inputs):  # noqa: F811
    _, baseline, _, _ = inputs
    catalog = _unsupported_build_catalog(tmp_path / "catalog.tsv")
    focal = tmp_path / "focal.txt"
    focal.write_text("G1\n")
    output = tmp_path / "annotated"

    with pytest.raises(GeneCatalogValidationError) as caught:
        cli.main(["annotate", "--ref-snp-annot-sources", str(baseline), "--query-annot-gene-list-sources", str(focal),
                  "--gene-coordinate-file", str(catalog), "--snp-identifier", "rsid", "--genome-build", "hg19",
                  "--padding-bp", "0", "--output-dir", str(output)])

    _assert_names_the_build_defect(str(caught.value), "catalog.tsv")
    table = pd.read_csv(output / "diagnostics/gene_catalog_issues.tsv", sep="\t", keep_default_na=False)
    assert tuple(table.columns) == CATALOG_ISSUE_COLUMNS
    assert table["reason"].tolist() == ["unsupported_genome_build"] * 2
