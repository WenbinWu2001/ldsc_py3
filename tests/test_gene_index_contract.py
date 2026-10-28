"""Public gene-index settings and shared SNP identity contracts."""

from argparse import Namespace
from pathlib import Path
import json
import shutil

import pandas as pd
import pytest

from ldsc import gene_ldscore_index as index
from ldsc._gene_index_scoring import load_gene_ldscore_index
from ldsc.config import GeneLDScoreIndexBuildConfig
from ldsc.errors import LDSCConfigError, LDSCInputError


def _arguments():
    return ["--baseline-annot-sources", "baseline.@.annot.gz", "--plink-prefix", "panel.@",
            "--gene-coordinate-file", "genes.tsv", "--output-dir", "index",
            "--genome-build", "hg38", "--snp-identifier", "chr_pos_allele_aware"]


@pytest.mark.parametrize("windows", [[], ["--ld-wind-cm", "1", "--ld-wind-kb", "100"]])
def test_builder_requires_exactly_one_window(windows):
    # Use the previously supported build/mode so this tests the window alone.
    args = _arguments()
    args[args.index("hg38")] = "hg19"
    args[args.index("chr_pos_allele_aware")] = "chr_pos"
    with pytest.raises(SystemExit):
        index.build_parser().parse_args(args + windows)


@pytest.mark.parametrize("unit,value", [("cm", "1.5"), ("kb", "100"), ("snps", "20")])
@pytest.mark.parametrize("mode", ["rsid", "chr_pos", "rsid_allele_aware", "chr_pos_allele_aware"])
def test_builder_accepts_all_windows_and_identity_modes(unit, value, mode):
    args = _arguments()
    args[-1] = mode
    parsed = index.build_parser().parse_args(args + [f"--ld-wind-{unit}", value])
    assert parsed.genome_build == "hg38"
    assert parsed.snp_identifier == mode
    assert getattr(parsed, f"ld_wind_{unit}") == float(value)
    assert sum(getattr(parsed, f"ld_wind_{name}") is not None for name in ("cm", "kb", "snps")) == 1


@pytest.mark.parametrize("window", [{}, {"ld_wind_cm": 1, "ld_wind_kb": 100}, {"ld_wind_cm": float("nan")}, {"ld_wind_snps": 1.5}])
def test_builder_python_config_rejects_invalid_windows(window):
    with pytest.raises(LDSCConfigError):
        GeneLDScoreIndexBuildConfig(baseline_annot_sources=("baseline",), plink_prefix="panel",
                                   output_dir="index", gene_coordinate_file="genes", genome_build="hg19",
                                   snp_identifier="rsid", **window)


@pytest.mark.parametrize("build,canonical", [("GRCh37", "hg19"), ("GRCh38", "hg38")])
def test_builder_normalizes_explicit_build_aliases(build, canonical):
    config = GeneLDScoreIndexBuildConfig(baseline_annot_sources=("baseline",), plink_prefix="panel",
                                       output_dir="index", gene_coordinate_file="genes", genome_build=build,
                                       snp_identifier="rsid_allele_aware", ld_wind_kb=100)
    assert config.genome_build == canonical


@pytest.mark.parametrize("mode", ["rsid_allele_aware", "chr_pos_allele_aware"])
@pytest.mark.parametrize("with_alleles", [False, True])
def test_builder_allele_matching_uses_shared_cleanup(mode, with_alleles):
    reference = pd.DataFrame({"CHR": ["1"] * 6, "POS": [10, 20, 30, 40, 50, 50],
                              "SNP": ["rs1", "rs2", "rs3", "rs4", "rs5", "rs5"],
                              "A1": ["A", "A", "A", "N", "A", "A"],
                              "A2": ["G", "C", "T", "G", "G", "C"]})
    baseline = reference.iloc[:4].copy()
    if with_alleles:
        baseline["A1"] = ["C", "T", "A", "N"]
        baseline["A2"] = ["T", "G", "T", "G"]
    else:
        baseline = baseline.drop(columns=["A1", "A2"])
    result = index.intersect_baseline_plink_by_identifier(
        baseline, pd.DataFrame({"base": [11, 22, 33, 44]}), reference, snp_identifier=mode, chrom="1")
    assert result.metadata.SNP.tolist() == ["rs1", "rs2"]
    assert result.annotations.base.tolist() == [11, 22]
    assert {"strand_ambiguous_allele", "invalid_allele", "multi_allelic_base_key"}.issubset(set(result.dropped_rows.reason))


@pytest.mark.parametrize("mode", ["rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"])
def test_builder_rejects_contradictory_alleles_with_repair_rows(mode):
    baseline = pd.DataFrame({"CHR": ["1"], "POS": [10], "SNP": ["rs1"], "A1": ["A"], "A2": ["G"]})
    reference = baseline.assign(A2="C")
    with pytest.raises(LDSCInputError, match="contradict") as raised:
        index.intersect_baseline_plink_by_identifier(
            baseline, pd.DataFrame({"base": [1]}), reference, snp_identifier=mode, chrom="1")
    assert raised.value.input_issues.iloc[0].reason == "allele_mismatch"


@pytest.mark.parametrize("mode", ["rsid", "rsid_allele_aware"])
def test_builder_rejects_shared_rsid_with_different_coordinates(mode):
    baseline = pd.DataFrame({"CHR": ["1"], "POS": [10], "SNP": ["rs1"]})
    reference = baseline.assign(POS=11, A1="A", A2="G")
    with pytest.raises(LDSCInputError, match="contradict") as raised:
        index.intersect_baseline_plink_by_identifier(
            baseline, pd.DataFrame({"base": [1]}), reference, snp_identifier=mode, chrom="1")
    assert raised.value.input_issues.iloc[0].reason == "coordinate_mismatch"


def test_unused_genetic_map_is_not_required_for_physical_window():
    assert index._load_builder_genetic_map(Namespace(genome_build="hg38", ld_wind_cm=None)) is None


def test_builder_reports_shared_rsid_on_different_chromosomes(tmp_path, monkeypatch):
    baseline = {chrom: pd.DataFrame({'CHR': [chrom, chrom], 'POS': [10, 20],
                                     'SNP': [f'rs{chrom}', 'moved' if chrom == '1' else 'base2']})
                for chrom in ('1', '2')}
    reference = {chrom: frame.assign(SNP=[f'rs{chrom}', 'moved' if chrom == '2' else 'ref1'])
                 for chrom, frame in baseline.items()}
    sources = Namespace(shards={chrom: Namespace(metadata=lambda chrom=chrom: baseline[chrom]) for chrom in baseline})
    monkeypatch.setattr(index, '_read_bim_identity', lambda prefix: reference[prefix])
    config = Namespace(snp_identifier='rsid', baseline_annot_sources=('baseline',))
    with pytest.raises(LDSCInputError, match='baseline_plink_issues.tsv'):
        index._validate_baseline_plink_inputs(sources, {'1': '1', '2': '2'}, config, ('1', '2'), tmp_path)
    issues = pd.read_csv(tmp_path / 'baseline_plink_issues.tsv', sep='\t')
    assert 'cross_chromosome_coordinate_mismatch' in set(issues.reason)


@pytest.mark.parametrize('build,mode,unit,value', [('hg19', 'rsid_allele_aware', 'cm', .15),
                                                ('hg38', 'chr_pos_allele_aware', 'snps', 2)])
def test_complete_autosomal_index_build_and_reload(tmp_path, build, mode, unit, value):
    fixture = Path(__file__).parent / 'fixtures' / 'plink' / 'plink'
    chromosomes = list(map(str, range(1, 23)))
    variants = []
    for chrom in chromosomes:
        prefix = tmp_path / f'panel.{chrom}'
        for suffix in ('.bed', '.fam'):
            shutil.copyfile(fixture.with_suffix(suffix), prefix.with_suffix(prefix.suffix + suffix))
        bim = pd.read_csv(fixture.with_suffix('.bim'), sep=r'\s+', header=None,
                          names=['CHR', 'SNP', 'CM', 'POS', 'A1', 'A2'])
        bim = bim.assign(CHR=chrom, SNP=[f'rs{chrom}_{i}' for i in range(len(bim))],
                         CM=bim.POS * .1, A1='A', A2='G')
        bim.to_csv(str(prefix) + '.bim', sep='\t', header=False, index=False)
        baseline = bim[['CHR', 'SNP', 'POS']].assign(base=1)
        if build == 'hg38':
            baseline = baseline.assign(A1='A', A2='G')
            baseline = pd.concat([baseline, baseline.iloc[[0]].assign(SNP=f'invalid_{chrom}', POS=999, A1='N')])
        baseline.to_csv(tmp_path / f'baseline.{chrom}.annot', sep='\t', index=False)
        variants.append(bim[['CHR', 'SNP', 'POS']])
    catalog = tmp_path / 'genes.tsv'
    pd.DataFrame({'gene_id': [f'G{c}' for c in chromosomes], 'gene_name': [f'Gene{c}' for c in chromosomes],
                  'chrom': chromosomes, 'start': 1, 'end': 5, 'genome_build': build}).to_csv(catalog, sep='\t', index=False)
    restriction = tmp_path / 'keep.tsv'
    pd.concat(variants).to_csv(restriction, sep='\t', index=False)
    config = GeneLDScoreIndexBuildConfig(
        baseline_annot_sources=(str(tmp_path / 'baseline.@.annot'),), plink_prefix=str(tmp_path / 'panel.@'),
        output_dir=str(tmp_path / 'index'), gene_coordinate_file=str(catalog), genome_build=build,
        snp_identifier=mode, regr_snps_file=str(restriction), regr_snps_exclude_regions='none',
        gene_exclude_regions='none', **{f'ld_wind_{unit}': value})
    path = index.build_gene_ldscore_index(config)
    loaded = load_gene_ldscore_index(path)
    assert loaded.chromosomes == tuple(chromosomes)
    assert loaded.genome_build == build and loaded.snp_identifier == mode
    assert loaded.index_identity['ld_window'] == {'unit': unit, 'value': value}
    if build == 'hg38':
        diagnostics = json.loads((path / 'diagnostics' / 'build-gene-ldscore-index.json').read_text())
        for chrom in chromosomes:
            assert diagnostics['chromosomes'][chrom]['baseline_identity_rows_dropped'] == 1
            dropped = pd.read_csv(path / 'diagnostics' / 'dropped_snps' / f'chr{chrom}_dropped.tsv.gz', sep='\t')
            assert len(dropped) == 1 and dropped.iloc[0].SNP == f'invalid_{chrom}'
    root = json.loads((path / 'metadata.json').read_text())
    assert root['format_version'] == 3
    # Version 2 indexes were built with legacy LDSC chunk-rounded PLINK windows.
    root['format_version'] = 2
    (path / 'metadata.json').write_text(json.dumps(root))
    with pytest.raises(LDSCInputError, match='Rebuild the complete index') as rejected:
        load_gene_ldscore_index(path)
    assert 'exact PLINK LD windows' in str(rejected.value)
    root.pop('format_version')
    (path / 'metadata.json').write_text(json.dumps(root))
    with pytest.raises(LDSCInputError, match='Rebuild the complete index'):
        load_gene_ldscore_index(path)
