from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from KEA import KEA


def fixture(tmp_path):
    ref=tmp_path/'ref.fa'; ctrl=tmp_path/'ctrl.fa'; gtf=tmp_path/'models.gtf'; bed=tmp_path/'bed'; bed.mkdir()
    ref.write_text('>ref.1\nAAAAAACCCCCC\n>minus.1\nAAAAAACCCCCC\n'); ctrl.write_text('>ctrl.1\nCCCAAACCCCCC\n')
    lines=[]
    for name,start,strand in [('ref',1,'+'),('ctrl',101,'+'),('minus',201,'-')]:
        for feature in ['exon','CDS']:
            lines.append(f'1\ttest\t{feature}\t{start}\t{start+11}\t.\t{strand}\t0\ttranscript_id "{name}";\n')
    gtf.write_text(''.join(lines)); (bed/'unipDomain.bed').write_text('chr1\t0\t3\tD\t0\t+\nchr1\t209\t212\tD\t0\t-\n')
    return KEA(str(tmp_path/'out'),'unused','unused'),ref,ctrl,gtf,bed


def test_domain_projection_counts_and_overrides(tmp_path):
    a,ref,ctrl,gtf,bed=fixture(tmp_path)
    result=a.DomainEnrichment(gtf,bed,ref_fasta=ref,ctrl_fasta=ctrl,kmers=['AAA'],plot=False)
    row=result['statistics'].iloc[0]
    assert [row.A,row.B,row.C,row.D]==[6,2,0,1]
    assert np.isinf(row.odds_ratio) and np.isfinite(row.log2_or)
    assert result['qc'].status.eq('used').all()
    assert a.reference_fasta=='unused'
    assert result['occurrences'].query('transcript == "minus" and start == 0').iloc[0].domains==('D',)
    with pytest.raises(FileExistsError):
        a.DomainEnrichment(gtf,bed,kmers=['AAA'],output_dir=result['output_dir'],plot=False)


def test_domain_shuffle_frame_and_heatmap(tmp_path):
    a,ref,ctrl,gtf,bed=fixture(tmp_path)
    kwargs=dict(ref_fasta=ref,kmers=['AAA'],background='shuffle',in_frame=True,seed=42)
    one=a.DomainEnrichment(gtf,bed,plot=True,**kwargs)
    two=a.DomainEnrichment(gtf,bed,plot=False,**kwargs)
    pd.testing.assert_frame_equal(one['occurrences'],two['occurrences'])
    assert one['occurrences'].frame.eq(0).all()
    assert len(one['occurrences'].query('group == "reference"'))==4
    assert len(one['occurrences'].query('group == "control"'))==4
    assert (Path(one['output_dir'])/'heatmap.png').exists()
    import matplotlib.pyplot as plt
    plt.close(one['figure'])


def test_kmap_metrics_and_sankey(tmp_path):
    a=KEA(str(tmp_path),'unused','unused')
    a.kmers_count_table_by_species['input']=pd.DataFrame({'query':[6,3,1,0],'same':[6,3,1,0],'other':[0,1,3,6],'k':[1]*4},index=['A','C','G','T'])
    for metric in ['spearman','pearson','cosine','jsd']:
        result=a.CompareTranscriptKmerProfiles('query',['same','other'],metric=metric)
        assert np.isclose(result.loc['same'].iloc[0],1)
        assert result.loc['same'].iloc[0]>result.loc['other'].iloc[0]
    path=tmp_path/'sankey.png'
    assert a.plot_kmer_rank_alignment('query','other',bin_size=2,save=str(path)) is None
    assert path.exists()
