"""Direct workflow contracts for output-contained, chromosome-owned annotations."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ldsc import run_ldscore
from tests.test_ldscore_chromosome_scope import write_r2_inputs


@pytest.mark.parametrize('batch_size', [1,2,1000])
def test_direct_file_queries_use_shards_and_leave_only_persistent_outputs(tmp_path,monkeypatch,batch_size):
    from ldsc.annotation_builder import AnnotationBuilder
    from ldsc import config
    monkeypatch.setattr(config, "_GLOBAL_CONFIG", config.GlobalConfig(snp_identifier="rsid"))

    args = write_r2_inputs(tmp_path)
    metadata = pd.read_csv(args.baseline_annot_sources,sep='\t')
    query = tmp_path/'queries.annot'
    metadata.drop(columns='base').assign(first=[0,1,0],second=[1,0,2]).to_csv(query,sep='\t',index=False)
    regression = tmp_path/'regression.txt'
    regression.write_text('SNP\nrs1\nrs3\n')
    def eager_path(*args,**kwargs):
        raise AssertionError('Whole-genome annotation preparation was invoked')
    monkeypatch.setattr(AnnotationBuilder,'run',eager_path)
    output = tmp_path/'out'
    result = run_ldscore(output_dir=output,baseline_annot_sources=args.baseline_annot_sources,
        query_annot_sources=str(query),r2_dir=args.r2_dir,
        regr_snps_file=regression,regr_snps_exclude_regions='none',
        ld_wind_snps=2,yes_really=True,query_batch_size=batch_size)
    assert result.query_table.SNP.tolist() == ['rs1','rs3']
    np.testing.assert_allclose(result.query_table[['first','second']],[[.5,1],[.25,2]],atol=2e-5,rtol=0)
    np.testing.assert_allclose(result.baseline_table.regression_ld_scores,[1,1],atol=0,rtol=0)
    assert [row['all_reference_snp_count'] for row in result.count_records] == [3,1,3]
    assert (output/'ldscore.baseline.parquet').is_file()
    assert (output/'ldscore.query.parquet').is_file()
    assert not list(output.glob('.ldsc-annotation-*'))
    pd.testing.assert_frame_equal(pd.read_parquet(output/'ldscore.query.parquet'),result.query_table)


def test_sequential_run_releases_prepared_chromosome_before_loading_next(tmp_path,monkeypatch):
    import weakref
    from ldsc import config
    from ldsc._kernel.ref_panel import ParquetR2RefPanel
    from tests.test_ldscore_parallelism import _PANEL_CHROMS, _write_index_panel

    monkeypatch.setattr(config,'_GLOBAL_CONFIG',config.GlobalConfig(snp_identifier='rsid'))
    panel = tmp_path/'panel'
    _write_index_panel(panel,_PANEL_CHROMS)
    regression = tmp_path/'regression.txt'
    regression.write_text('SNP\n'+'\n'.join(f'rs{i}' for i in range(1,8))+'\n')
    previous = []
    loads = []
    original = ParquetR2RefPanel.prepare_chromosome
    def measured(self,chrom,*args,**kwargs):
        assert all(reference() is None for reference in previous)
        state = original(self,chrom,*args,**kwargs)
        previous[:] = [weakref.ref(state),weakref.ref(state.metadata),weakref.ref(state.annotations)]
        loads.append(chrom)
        return state
    monkeypatch.setattr(ParquetR2RefPanel,'prepare_chromosome',measured)
    run_ldscore(output_dir=str(tmp_path/'out'),r2_dir=str(panel),
        regr_snps_file=str(regression),regr_snps_exclude_regions='none',
        ld_wind_kb=1,yes_really=True,threads=1)
    assert loads == ['1','2']
    assert all(reference() is None for reference in previous)


@pytest.mark.parametrize('batch_size',[0,-1,1.5,True,None])
def test_query_batch_size_requires_a_positive_integer(batch_size):
    from ldsc.config import LDScoreConfig
    from ldsc.errors import LDSCConfigError

    with pytest.raises(LDSCConfigError,match='positive integer'):
        LDScoreConfig(ld_wind_cm=1,query_batch_size=batch_size)
