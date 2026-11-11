import numpy as np
import pandas as pd
from KEA import KEA

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
