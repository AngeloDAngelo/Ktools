from pathlib import Path
import pandas as pd
import pytest
from KEA import KEA

RESOURCE=Path(__file__).resolve().parents[1]/'data/mirna/miR_Family_Info.human.tsv'


def test_mir7_orientation_and_human_filter(tmp_path):
    a=KEA(str(tmp_path),'unused','unused')
    r=a.AnnotateMiRNASeeds(kmers=['GTCTTCC','GGAAGAC'],source='manual',plot=False)
    hit=r['matches'].query('mirna_family == "miR-7-5p"')
    assert len(hit)==1 and hit.iloc[0].kmer=='GTCTTCC'
    assert hit.iloc[0].seed_rna=='GGAAGAC'
    assert 'hsa-miR-7-5p' in hit.iloc[0].mirnas
    with pytest.raises(ValueError):a.AnnotateMiRNASeeds(kmers=['AAAA'],source='manual',plot=False)


def test_krs_scores_plot_and_empty(tmp_path):
    a=KEA(str(tmp_path),'unused','unused')
    a.KEA_results={'input':{'KRS':{'reference':'target','signature':['GTCTTCC'],'kmer_ranks':pd.DataFrame({'kmer':['GTCTTCC'],'target_frequency':[0.04],'target_rank_percentile':[0.99],'selected_KRS':[True]})}}}
    result=a.AnnotateMiRNASeeds(source='KRS')
    assert result['matches'].target_rank_percentile.eq(0.99).all()
    assert (Path(result['output_dir'])/'seed_scores.png').exists()
    import matplotlib.pyplot as plt
    for fig in result['figures']:plt.close(fig)
    a.KEA_results['input']['KRS']['signature']=[]
    assert a.AnnotateMiRNASeeds(source='KRS',plot=False)['matches'].empty


def test_kea_intersection_preserves_two_scores(tmp_path):
    a=KEA(str(tmp_path),'unused','unused')
    a.KEA_results={'input':{'delta_median':{'enriched':['GTCTTCC','AAAAAAA'],'all_kmers':pd.DataFrame({'delta_median':[2.,1.]},index=['GTCTTCC','AAAAAAA'])},'stat_log2fc':{'enriched':['GTCTTCC'],'all_kmers':pd.DataFrame({'log2FoldChange':[1.5],'FDR':[.01]},index=['GTCTTCC'])}}}
    result=a.AnnotateMiRNASeeds(source='KEA',plot=False)
    assert result['query'].kmer.tolist()==['GTCTTCC']
    assert result['matches'].delta_median.eq(2).all()
    assert result['matches'].log2FoldChange.eq(1.5).all()


def test_streaming_case_study_matches_krs_definition(tmp_path):
    import importlib.util
    path=Path(__file__).resolve().parents[1]/'examples/run_cdr1as_mirna.py'
    spec=importlib.util.spec_from_file_location('cdr_example',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    sequences={'target':'GTCTTCCGTCTTCCAAAAAAA','b1':'GTCTTCCAAAAAAAAAAAAAA','b2':'AAAAAAAAAAAAAAAAAAAA','b3':'GTCTTCCGTCTTCCAAAAAAA'}
    counts={name:module.counts(seq) for name,seq in sequences.items()}
    table=pd.DataFrame(counts).fillna(0);table['k']=7
    a=KEA(str(tmp_path),'unused','unused');a.kmers_count_table_by_species['input']=table
    a.KRS('target',save=False,verbose=False)
    ranked=a.KEA_results['input']['KRS']['kmer_ranks'].set_index('kmer')
    for word in counts['target']:
        target=counts['target'][word]/sum(counts['target'].values())
        population=[counter[word]/sum(counter.values()) for counter in counts.values() if counter[word]>0]
        lower=sum(value<target for value in population);equal=sum(value==target for value in population)
        streamed=(lower+(equal+1)/2)/len(population)
        assert abs(streamed-ranked.loc[word,'target_rank_percentile'])<1e-12
