"""Real CDR1as case study: exact KRS ranks against a user-provided FASTA.

Streams background transcripts to avoid a dense transcriptome count matrix.
Only target-present canonical 7-mers need ranks; absent target words cannot
enter KRS. Uses the same nonzero population and average-tie rank as KEA.KRS.
"""
import argparse
from collections import Counter
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from Bio import SeqIO
from KEA import KEA


def counts(seq):
    seq=str(seq).upper().replace('U','T')
    return Counter(seq[i:i+7] for i in range(len(seq)-6) if set(seq[i:i+7])<=set('ACGT'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--background',required=True)
    parser.add_argument('--out',default='results/cdr1as_mirna')
    parser.add_argument('--top-pct',type=float,default=5)
    args=parser.parse_args()
    if not 0<args.top_pct<=100:raise ValueError('top-pct must be in (0,100]')
    root=Path(__file__).resolve().parents[1]
    target_path=root/'examples/data/CDR1as_GRCh38.fa'
    target=next(SeqIO.parse(target_path,'fasta'))
    c=counts(target.seq); words=sorted(c); freq=np.array([c[w]/sum(c.values()) for w in words])
    lower=np.zeros(len(words),dtype=int);equal=np.ones(len(words),dtype=int);n=np.ones(len(words),dtype=int)
    seen={target.id}; used=0; identical=0
    for rec in SeqIO.parse(args.background,'fasta'):
        if rec.id in seen:raise ValueError('Duplicate target/background identifier: '+rec.id)
        seen.add(rec.id)
        if str(rec.seq).upper().replace('U','T')==str(target.seq):
            identical+=1;continue
        d=counts(rec.seq);total=sum(d.values())
        if total==0:continue
        values=np.array([d.get(w,0)/total for w in words]);present=values>0
        n+=present;lower+=present & (values<freq);equal+=present & (values==freq);used+=1
    if not used:raise ValueError('No usable background transcripts')
    ranks=(lower+(equal+1)/2)/n
    table=pd.DataFrame({'kmer':words,'target_frequency':freq,'target_rank_percentile':ranks,'selected_KRS':ranks>=1-args.top_pct/100,'target_count':[c[w] for w in words],'nonzero_population':n})
    out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    table.to_csv(out/'target_kmer_ranks.tsv',sep='\t',index=False)
    analysis=KEA(str(out),str(target_path),args.background)
    analysis.KEA_results={'input':{'KRS':{'reference':target.id,'signature':table.loc[table.selected_KRS,'kmer'].tolist(),'kmer_ranks':table}}}
    result=analysis.AnnotateMiRNASeeds(root/'data/mirna/miR_Family_Info.human.tsv',source='KRS',top_n=20)
    # Also preserve miR-7 even if outside the top 20 plotted annotations.
    mir7=table[table.kmer=='GTCTTCC'];mir7.to_csv(out/'mir7_target_score.tsv',sep='\t',index=False)
    summary={'target':target.id,'target_length':len(target),'background_transcripts_used':used,'identical_sequences_excluded':identical,'background_sha256':hashlib.sha256(Path(args.background).read_bytes()).hexdigest(),'top_pct':args.top_pct,'selected_kmers':int(table.selected_KRS.sum()),'miR7_7mer_m8':mir7.to_dict(orient='records'),'annotation_output':Path(result['output_dir']).name,'circular_junction_windows':'Not counted; sequence treated as one linearized traversal'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
