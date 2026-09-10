"""Gene resolution remains exact across streamed source boundaries."""

import pandas as pd

from ldsc._annotation_storage import AnnotationWorkspace
from ldsc.gene_list_resolver import GeneCatalog


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
