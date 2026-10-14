from pathlib import Path
from types import SimpleNamespace
import json
import shutil
import subprocess
import numpy as np
import pandas as pd
import pytest
from KEA import KEA
from test_modules import fixture


def install_mock_renderer(monkeypatch, fail=False):
    calls = []
    monkeypatch.setattr(shutil, 'which', lambda name: 'Rscript-test')
    def run(command, **kwargs):
        calls.append(command)
        if '-e' in command:
            return SimpleNamespace(returncode=0, stdout='', stderr='')
        if fail:
            return SimpleNamespace(returncode=1, stdout='', stderr='renderer failed')
        script = Path(command[2]).read_text(encoding='utf-8')
        assert "method='ward.D2'" in script and 'profile <- p_cap - rowMeans(p_cap)' in script
        assert 'max_domains' not in script and "domain_limit" in script
        folder = Path(command[4])
        for stem in ['main_domains','supplement_positions_frames']:
            for ext in ['pdf','png','svg']:
                (folder/f'{stem}.{ext}').write_text('mock renderer artifact')
        pd.DataFrame({'kmer':['AAA','CCC'],'cluster':['CDS start','CDS body']}).to_csv(
            folder/'cluster_assignments.tsv',sep='\t',index=False)
        return SimpleNamespace(returncode=0, stdout='rendered', stderr='')
    monkeypatch.setattr(subprocess, 'run', run)
    return calls


def positional():
    return pd.DataFrame({'kmer':['AAA','AAA','CCC','CCC'],
                         'bin':['CDS_1','CDS_2']*2,
                         'log2_or':[1.5,0.0,0.3,0.3], 'fdr':[1.0]*4})


def test_kcds_reuses_memory_and_only_emits_new_panels(tmp_path, monkeypatch):
    a,ref,ctrl,gtf,bed=fixture(tmp_path)
    saved=positional()
    a.enrichment=saved.drop(columns='kmer').assign(
        analysis=saved.kmer+'_Enr_VS_Ctrl',type_analysis='kmers')
    old=a.enrichment.copy(deep=True)
    monkeypatch.setattr(a,'EnrichmentAnalysis',lambda *args,**kwargs:pytest.fail('cached analysis recomputed'))
    calls=install_mock_renderer(monkeypatch)
    result=a.K_CDS(gtf,bed,ref_fasta=ref,ctrl_fasta=ctrl,kmers=['AAA','CCC'],tracks=('unipDomain',),show_final=False)
    assert len(calls)==2
    assert calls[-1][-5:]==['0.5','0.1','0.05','2','60']
    pd.testing.assert_frame_equal(a.enrichment,old)
    assert len(result['cluster_assignments'])==2
    assert result['figure'] is None and result['frame_figure'] is None
    assert set(result['final_figures'])=={'main_domains','supplement_positions_frames'}
    assert all(Path(path).is_file() for paths in result['final_figures'].values() for path in paths.values())
    assert not (Path(result['output_dir'])/'heatmap.png').exists()
    assert not (Path(result['output_dir'])/'frame_distribution.png').exists()
    assert json.loads((Path(result['output_dir'])/'settings.json').read_text())['final_plots'] is True


def test_kcds_saved_profile_reuse_and_r_failure_preserves_tables(tmp_path,monkeypatch):
    a,ref,ctrl,gtf,bed=fixture(tmp_path)
    cache=Path(a.dir_out)/'K_CDS';cache.mkdir()
    positional().to_csv(cache/'positional_enrichment.tsv',sep='\t',index=False)
    install_mock_renderer(monkeypatch,fail=True)
    with pytest.raises(RuntimeError,match='renderer failed'):
        a.K_CDS(gtf,bed,ref_fasta=ref,kmers=['AAA','CCC'],tracks=('unipDomain',),show_final=False)
    assert (Path(a.k_cds['output_dir'])/'statistics.tsv').exists()
    assert 'renderer failed' in (Path(a.k_cds['output_dir'])/'heatmap_R.log').read_text()


def test_kcds_explicit_invalid_profile_fails_before_counting(tmp_path,monkeypatch):
    a,ref,ctrl,gtf,bed=fixture(tmp_path)
    install_mock_renderer(monkeypatch)
    with pytest.raises(ValueError,match='incomplete'):
        a.K_CDS(gtf,bed,ref_fasta=ref,kmers=['AAA','CCC'],tracks=('unipDomain',),
                positional_enrichment=positional().iloc[:2],show_final=False)


def test_kcds_missing_r_is_actionable(tmp_path,monkeypatch):
    a,ref,ctrl,gtf,bed=fixture(tmp_path)
    monkeypatch.setattr(shutil,'which',lambda _:None)
    with pytest.raises(RuntimeError,match='Rscript not found'):
        a.K_CDS(gtf,bed,ref_fasta=ref,kmers=['AAA','CCC'],tracks=('unipDomain',))


def test_kcds_computes_only_missing_sevenmer_and_restores_enrichment(tmp_path,monkeypatch):
    a,ref,ctrl,gtf,bed=fixture(tmp_path)
    ref.write_text('>ref.1\nAAAAAAACCCCCCCAAAAAAA\n')
    ctrl.write_text('>ctrl.1\nCCCCCCCAAAAAAACCCCCCC\n')
    gtf.write_text(''.join(f'1\ttest\t{feature}\t{start}\t{start+20}\t.\t+\t0\ttranscript_id "{name}";\n'
                           for name,start in [('ref',1),('ctrl',101)] for feature in ['exon','CDS']))
    a.metatx=pd.DataFrame([(name,lo,hi,label,1,'+') for name in ['ref.1','ctrl.1']
                          for lo,hi,label in [(0,7,'CDS_1'),(7,14,'CDS_2'),(14,21,'3UTR_1')]],
                         columns=['chrom','start','end','name','score','strand'])
    a.enrichment=pd.DataFrame({'analysis':['AAAAAAA_Enr_VS_Ctrl']*3,
                              'bin':['CDS_1','CDS_2','3UTR_1'],'log2_or':[1,0,0],'fdr':[1]*3})
    original=a.enrichment
    computed=[]
    def enrich(features,kmer,*args,**kwargs):
        computed.append(kmer)
        assert len(features)>0 and features['name'].eq(kmer).all()
        a.enrichment=None
        return pd.DataFrame({'bin':['CDS_1','CDS_2','3UTR_1'],'log2_or':[0,0,0],'fdr':[1]*3})
    monkeypatch.setattr(a,'EnrichmentAnalysis',enrich)
    install_mock_renderer(monkeypatch)
    result=a.K_CDS(gtf,bed,ref_fasta=ref,ctrl_fasta=ctrl,kmers=['AAAAAAA','CCCCCCC'],
                  tracks=('unipDomain',),show_final=False)
    assert computed==['CCCCCCC']
    assert a.enrichment is original
    assert set(result['positional_enrichment'].kmer)=={'AAAAAAA','CCCCCCC'}


def test_krbp_original_word_scores_and_unchanged_return(tmp_path):
    a=KEA(str(tmp_path),'unused','unused')
    scores=pd.DataFrame({'RBP_A_HepG2':[1,4,7,20,0,0,0],
                         'RBP_A_K562':[2,np.nan,8,5,0,0,0]},
                        index=['ACGTT','CGTTA','GTTAC','AAAAA','CCCCC','GGGGG','ACACA'])
    original=scores.copy(deep=True)
    result=a.KmersProteinEnrichment(['ACGTTAC','AAAAAAA','TTTTTTT','ACGTTAC'],scores,k_out=5,plot=False)
    table=a.kmer_protein_scores
    assert table.index.tolist()==['ACGUUAC','AAAAAAA','UUUUUUU']
    assert table.columns.tolist()==scores.columns.tolist()
    assert table.loc['ACGUUAC','RBP_A_HepG2']==4
    assert np.isnan(table.loc['ACGUUAC','RBP_A_K562'])
    assert table.loc['AAAAAAA','RBP_A_HepG2']==20
    assert table.loc['UUUUUUU'].isna().all()
    assert {'pvalue','fdr','bonf','zscore_diff'}.issubset(result.columns)
    pd.testing.assert_frame_equal(scores,original)
