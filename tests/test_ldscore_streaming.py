"""Direct workflow contracts for output-contained, chromosome-owned annotations."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ldsc.regression_runner import load_ldscore_from_dir
from ldsc import run_ldscore
from tests.test_ldscore_chromosome_scope import write_r2_inputs


@pytest.mark.parametrize('batch_size', [1,2,1000])
def test_direct_file_queries_use_shards_and_leave_only_persistent_outputs(tmp_path,monkeypatch,batch_size):
    from ldsc.annotation_builder import AnnotationBuilder
    from ldsc import ldscore_calculator
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
    original_diagnostics = ldscore_calculator._stage_chromosome_drops
    diagnostic_paths = []
    def stage_diagnostics(*args, **kwargs):
        assert all(not path.exists() for path in diagnostic_paths), "previous batch diagnostics remain in scratch"
        result = original_diagnostics(*args, **kwargs)
        diagnostic_paths.append(result.path)
        return result
    monkeypatch.setattr(ldscore_calculator, '_stage_chromosome_drops', stage_diagnostics)
    output = tmp_path/'out'
    result = run_ldscore(output_dir=output,baseline_annot_sources=args.baseline_annot_sources,
        query_annot_sources=str(query),r2_dir=args.r2_dir,
        regr_snps_file=regression,regr_snps_exclude_regions='none',
        ld_wind_snps=2,yes_really=True,query_batch_size=batch_size)
    assert result.read_queries(result.query_columns).SNP.tolist() == ['rs1','rs3']
    np.testing.assert_allclose(result.read_queries(result.query_columns)[['first','second']],[[.5,1],[.25,2]],atol=2e-5,rtol=0)
    np.testing.assert_allclose(result.baseline_table.regression_ld_scores,[1,1],atol=0,rtol=0)
    assert [row['all_reference_snp_count'] for row in result.count_records] == [3,1,3]
    assert (output/'ldscore.baseline.parquet').is_file()
    assert all(entry['path'].is_file() for entry in result.query_batches)
    assert not list(output.glob('.ldsc-annotation-*'))
    pd.testing.assert_frame_equal(load_ldscore_from_dir(str(output)).read_queries(result.query_columns),result.read_queries(result.query_columns))


def test_direct_results_keep_complete_persistent_drop_records(tmp_path, monkeypatch):
    from ldsc import config
    monkeypatch.setattr(config, "_GLOBAL_CONFIG", config.GlobalConfig(snp_identifier="rsid"))
    args = write_r2_inputs(tmp_path)
    path = Path(args.baseline_annot_sources)
    frame = pd.read_csv(path, sep="\t")
    removed_chromosome = frame.iloc[[0, 0]].assign(CHR=2, SNP="removed", POS=[10, 20])
    pd.concat([frame, frame.iloc[[0]], removed_chromosome], ignore_index=True).to_csv(path, sep="\t", index=False)
    regression = tmp_path / "regression.txt"
    regression.write_text("SNP\nrs1\nrs3\n")
    output = tmp_path / "out"
    result = run_ldscore(output_dir=output, baseline_annot_sources=str(path), r2_dir=args.r2_dir,
                        regr_snps_file=regression, regr_snps_exclude_regions="none",
                        ld_wind_snps=2, yes_really=True)
    assert not list(output.glob(".ldsc-annotation-*"))
    artifact = result.identity_drops_by_chrom["22"]
    assert not isinstance(artifact, pd.DataFrame)
    assert artifact.path.is_file()
    drops = pd.concat(artifact.frames(), ignore_index=True)
    cleaned = drops.loc[drops.stage.eq("annotation_identity_cleanup")]
    assert cleaned.SNP.tolist() == ["rs1", "rs1"]
    assert cleaned.reason.tolist() == ["duplicate_identity", "duplicate_identity"]
    assert pd.concat(result.identity_drops_by_chrom["2"].frames()).SNP.tolist() == ["removed", "removed"]


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
