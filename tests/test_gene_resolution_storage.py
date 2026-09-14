"""Gene resolution remains exact across streamed source boundaries."""

import pandas as pd
import pytest

from ldsc._annotation_storage import AnnotationWorkspace
from ldsc.gene_list_resolver import GeneCatalog


@pytest.mark.parametrize("chunk_rows", [1, 3])
def test_exclusion_log_summarizes_each_source_across_audit_chunks(tmp_path, caplog, chunk_rows):
    from ldsc._gene_query_storage import resolve_gene_lists_staged
    from ldsc.query_annotations import _log_gene_list_rejections

    catalog_path = tmp_path / "catalog.tsv"
    catalog_path.write_text("gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\n"
                           "MHC1\tFIRST\t6\t30000000\t30000100\thg19\n"
                           "MHC2\tSECOND\t6\t31000000\t31000100\thg19\n"
                           "G1\tONE\t1\t10\t20\thg19\n")
    focal, retained, other = [tmp_path / f"{name}.txt" for name in ("focal", "retained", "other")]
    focal.write_text("MHC1\n\nSECOND\nG1\nUNKNOWN\n")
    retained.write_text("G1\n")
    other.write_text("MHC2\n")
    with AnnotationWorkspace(tmp_path / "out") as workspace:
        batch = resolve_gene_lists_staged([focal, retained, other], GeneCatalog.load(catalog_path), workspace,
                                         control_path=focal, gene_exclude_regions="mhc",
                                         resolution_policy="resolved-only", chunk_rows=chunk_rows)
        before = pd.concat(batch.audit_frames(), ignore_index=True)
        with caplog.at_level("INFO", logger="LDSC"):
            _log_gene_list_rejections(batch)
        messages = [record.getMessage() for record in caplog.records]
        exclusions = [message for message in messages if "intentionally excluded" in message]
        assert exclusions == [
            "Genes intentionally excluded by region policy: role=focal source=focal.txt count=2 line:gene=[1:MHC1, 3:SECOND->MHC2]",
            "Genes intentionally excluded by region policy: role=focal source=other.txt count=1 line:gene=[1:MHC2]",
            "Genes intentionally excluded by region policy: role=control source=focal.txt count=2 line:gene=[1:MHC1, 3:SECOND->MHC2]",
        ]
        rejected = [message for message in messages if "Gene-list row rejected:" in message]
        assert len(rejected) == 2 and all("line=5 input_gene='UNKNOWN'" in message for message in rejected)
        pd.testing.assert_frame_equal(pd.concat(batch.audit_frames(), ignore_index=True), before)
        assert batch.selection("focal", 1).canonical_gene_ids == ("G1",)


def test_streamed_resolution_deduplicates_aliases_across_chunks(tmp_path):
    from ldsc._gene_query_storage import resolve_gene_lists_staged

    catalog_path = tmp_path / "catalog.tsv"
    pd.DataFrame([['G1', 'ONE', '1', 101, 110, 'hg19'], ['G2', 'TWO', '2', 201, 210, 'hg19']],
                 columns=['gene_id', 'gene_name', 'chrom', 'start', 'end', 'genome_build']).to_csv(catalog_path, sep='\t', index=False)
    genes = tmp_path / 'focal.txt'
    genes.write_text('ONE\n\nG2\nG1\nUNKNOWN\n')
    empty = tmp_path / 'empty.txt'
    empty.write_text('')
    catalog = GeneCatalog.load(catalog_path)
    with AnnotationWorkspace(tmp_path / 'output') as workspace:
        batch = resolve_gene_lists_staged([genes, empty], catalog, workspace, resolution_policy='resolved-only', chunk_rows=1)
        assert not batch.has_fatal_gate_a_issues
        audit = pd.concat(batch.audit_frames(), ignore_index=True)
        assert audit['line'].tolist() == [1, 3, 4, 5]
        assert audit['disposition'].tolist() == ['retained', 'retained', 'duplicate', 'rejected']
        assert audit.loc[2, 'details'] == 'first occurrence at line 1'
        assert batch.selection('focal', 1).canonical_gene_ids == ('G1', 'G2')
        assert batch.selection('focal', 2).canonical_gene_ids == ()
        assert batch.summary['unique_resolved_genes'].tolist() == [2, 0]
        assert batch.summary['duplicate_rows'].tolist() == [1, 0]
        supported = batch.with_snp_support({0: 2, 1: 0}, support_kind='annotation')
        audit = pd.concat(supported.audit_frames(), ignore_index=True)
        assert audit['reference_snp_count'].isna().all()
        assert audit['annotation_snp_count'].tolist()[:3] == [2, 0, 2]
        assert audit.loc[1, 'reason'] == 'zero_annotation_snp_support'
