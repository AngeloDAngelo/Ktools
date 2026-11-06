import os
import sys
import shutil
import shlex
import subprocess
from pathlib import Path
from typing import List, Optional, Union
import multiprocessing as mp
from sklearn.cluster import AgglomerativeClustering
import pandas as pd
pd.set_option("future.no_silent_downcasting", True) # added to stop warning
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
import seaborn as sns
import logomaker
import matplotlib.patches as mpatches
from scipy.stats import ks_2samp, ttest_ind, mannwhitneyu, spearmanr, pearsonr
from statsmodels.stats.multitest import multipletests
from scipy.spatial.distance import cosine, jensenshannon
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pybiomart import Dataset
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.metrics import silhouette_score
from scipy import stats
from statsmodels.stats import multitest
from collections import defaultdict
from scipy.cluster.hierarchy import linkage, leaves_list
from wordcloud import WordCloud
import matplotlib as mpl
import matplotlib.font_manager as fm
from matplotlib_venn import venn2, venn3
import gzip
import re
import tempfile
import warnings
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import json
import re
import tempfile
from pathlib import Path
import textwrap
def LocalizeRangeInTx(table,annotation_file,tx_id_column,start_col,end_col):
    
    COLONNA_ID=tx_id_column#"ID"#"id1"
    COLONNA_START=start_col#"I"#"start1"
    COLONNA_END=end_col#"F"#"end1"

    utr5_col="UTR5"
    cds_col="CDS"
    utr3_col="UTR3"

    ### file info_sequences

    ## info file is located here ##
    ## "/data01/common_files/utilities/ensembl/ens99/annotation/annotation_"+specie_file+"_"+genome+".txt"

    ### file output

    ### CODE ###

    # importa dataframe con coordinate da localizzare
    Df=table
    try:
        del Df["Unnamed: 0"]
    except:
        pass

    print("numero righe file iniziale: ",len(Df))

    # importa file con coordinate delle regioni di ciascun trascritto
    Info=pd.read_table(annotation_file,sep="\t")
    try:
        del Info["Unnamed: 0"]
    except:
        pass

    # parsing del file co coordinate delle regioni
    Info.loc[Info[utr5_col]=="/",[utr5_col,cds_col,utr3_col]]=0

    Info[utr5_col]=Info[utr5_col].astype("int64")
    Info[cds_col]=Info[cds_col].astype("int64")
    Info[utr3_col]=Info[utr3_col].astype("int64")

    # localizzazione
    
    def Localization(Df,Info,ColMergeDf,ColMergeFileRegions,StartCol,EndCol,ColName):
        G4=Df
        Af=Info

        G4[StartCol]=G4[StartCol].astype("int64")
        G4[EndCol]=G4[EndCol].astype("int64")

        G4=pd.merge(G4,Af,left_on=ColMergeDf,right_on=ColMergeFileRegions)
        G4=G4.reset_index()
        del G4["index"]

        G4[ColName+"_localization"]="/"

        G4.loc[(G4["transcript_biotype"]=="protein_coding")&(G4[EndCol]<=G4[utr5_col]),[ColName+"_localization"]]="5UTR"
        G4.loc[(G4["transcript_biotype"]=="protein_coding")&(G4[utr5_col]<G4[StartCol])&(G4[EndCol]<=G4[cds_col]),[ColName+"_localization"]]="CDS"
        G4.loc[(G4["transcript_biotype"]=="protein_coding")&(G4[cds_col]<G4[StartCol])&(G4[EndCol]<=G4[utr3_col]),[ColName+"_localization"]]="3UTR" 
        G4.loc[(G4["transcript_biotype"]=="protein_coding")&(G4[utr5_col]>G4[StartCol])&(G4[EndCol]>=G4[utr5_col]),[ColName+"_localization"]]="5UTR-CDS" 
        G4.loc[(G4["transcript_biotype"]=="protein_coding")&(G4[cds_col]>G4[StartCol])&(G4[EndCol]>=G4[cds_col]),[ColName+"_localization"]]="CDS-3UTR" 
        return G4

    Df=Localization(Df,Info,COLONNA_ID,"ID",COLONNA_START,COLONNA_END,COLONNA_ID)

    # parsing file output

    del Df["5utr"]
    del Df["cds"]
    del Df["3utr"]
    del Df["length"]
    del Df["cdna"]
    del Df["ID"]
    del Df["UTR5"]
    del Df["CDS"]
    del Df["UTR3"]
    
    print("numero righe file finale: ",len(Df))
    
    del Df['ensembl_gene_id']
    del Df['ensembl_transcript_id']
    del Df['external_gene_name']
    del Df['transcript_biotype']
    del Df['gene_biotype']
    
    return Df

def Localize(self,ann_file):
        print()
        print("Localize")
        print("N names:",len(set(self.bed.loc[:,"name"])))
        
        self.region_bed=LocalizeRangeInTx(self.bed,ann_file,"chrom","start","end")
        print()
        print("Summary Ranges localization in transcripts:")
        print(self.region_bed.groupby("chrom_localization")["name"].apply(lambda x: len(set(x)))) 
        print(self.region_bed.head())


def FileByList(file,lista):
    
    with open(file,"w") as f:
        for i in lista:
            f.write(i+"\n")
        f.close()

def ReindexMetagene(res_both, col_name="bin"):

    res_both = res_both.reset_index()
    res_both["region"] = res_both[col_name].astype(str).str.split("_").str[0]
    res_both["n_bin"] = res_both[col_name].astype(str).str.split("_").str[1]
    res_both["region"] = (
        res_both["region"]
        .astype(str)
        .str.upper()
        .str.strip())

    res_both["n_bin"] = pd.to_numeric(res_both["n_bin"], errors="coerce")
    cat = pd.CategoricalDtype(categories=["5UTR", "CDS", "3UTR"], ordered=True)
    res_both["region"] = res_both["region"].astype(cat)
    res_both = res_both.sort_values(["region", "n_bin"], na_position="last")
    return res_both.set_index(col_name)

    
def FisherExactTest_MultiTest(counts_final,ref_group_name,col_a,col_b,comp_group_name,col_c,col_d):
    
    print("tested bins:",len(counts_final.index))
    print("modified fisher")
    
    for i in counts_final.index:

        print(i,end="\r")

        a=counts_final.loc[i,col_a].item()
        b=counts_final.loc[i,col_b].item()
        c=counts_final.loc[i,col_c].item()
        d=counts_final.loc[i,col_d].item()

        OR=stats.fisher_exact([[a,b],[c,d]])[0]
        
        if abs(OR) == np.inf:
            pseudo=0.001
            OR=((a+pseudo)/(b+pseudo))/((c+pseudo)/(d+pseudo))
            
            
        if abs(np.log2(OR)) == np.inf:
            pseudo=0.001
            OR=((a+pseudo)/(b+pseudo))/((c+pseudo)/(d+pseudo))

        pvalue=stats.fisher_exact([[a,b],[c,d]])[1]

        counts_final.loc[i,"pvalue"]=pvalue
        counts_final.loc[i,"or"]=OR

    
    counts_final.loc[:,"freq_1"]=counts_final.loc[:,col_a]/(counts_final.loc[:,col_a]+counts_final.loc[:,col_b])
    counts_final.loc[:,"freq_2"]=counts_final.loc[:,col_c]/(counts_final.loc[:,col_c]+counts_final.loc[:,col_d])
    
    counts_final.loc[:,"fdr"]=multitest.multipletests(list(counts_final.loc[:,"pvalue"]),method="fdr_bh")[1]
    counts_final.loc[:,"bonf"]=counts_final.loc[:,"pvalue"]*len(counts_final)
    counts_final.loc[:,"log2_or"]=np.log2(counts_final.loc[:,"or"].astype(float))
    
    counts_final.loc[:,"analysis"]=ref_group_name+"_VS_"+comp_group_name
    
    counts_final=counts_final[[col_a,col_b,col_c,col_d,"freq_1","freq_2","or","log2_or","pvalue","fdr","bonf","analysis"]]
    counts_final.columns=["A","B","C","D","freq_1","freq_2","or","log2_or","pvalue","fdr","bonf","analysis"]
    
    return counts_final
    
def IntersectBedDf(bed1,bed2,options,dir_temp="./temp_bedtools/"):
    
    os.makedirs(dir_temp,exist_ok=True)
    bed1.loc[:,"start"]=bed1.loc[:,"start"].astype(int)
    bed1.loc[:,"end"]=bed1.loc[:,"end"].astype(int)
    bed2.loc[:,"start"]=bed2.loc[:,"start"].astype(int)
    bed2.loc[:,"end"]=bed2.loc[:,"end"].astype(int)
    
    bed1=bed1.sort_values(["chrom","start","end","strand"])
    bed2=bed2.sort_values(["chrom","start","end","strand"])
    
    file_a=dir_temp+"bed1.bed"
    file_b=dir_temp+"bed2.bed"
    output=dir_temp+"output.txt"
    
    bed1.to_csv(file_a,sep="\t",index=None,header=None)
    bed2.to_csv(file_b,sep="\t",index=None,header=None)
    bedtools = shutil.which("bedtools")
    if bedtools is None:
        raise FileNotFoundError("Install BEDTools and add bedtools to PATH")
    with open(output, "w") as handle:
        subprocess.run([bedtools, "intersect", *shlex.split(options),
                        "-a", file_a, "-b", file_b], stdout=handle, check=True)
    out=pd.read_table(output,header=None)
    if "temp_bedtools" in dir_temp:
        os.system("rm -r "+dir_temp)
    return out
    
def KmersCountsInMetaTx(metatx_bed,kmers_bed_input,type_metatx="tx",k=7):
    
    kmers_bed=kmers_bed_input.copy(deep=True)
    
    ### pointed kmer bed
    kmers_bed.loc[:,"end"]=kmers_bed.loc[:,"start"]+1
    kmer_counts = IntersectBedDf(metatx_bed,kmers_bed,'-wao -s')
    kmer_counts.columns=[str(x)+"_1" for x in metatx_bed.columns]+[str(x)+"_2" for x in kmers_bed.columns]+["coverage"]

    kmer_counts['length_bin'] =  kmer_counts["end_1"]-kmer_counts["start_1"]
    
    if type_metatx == "tx":
        
        kmer_counts.loc[:,"bin_1"]=kmer_counts.loc[:,"name_1"]
        
    #if type_metatx ==" gx":   
    #    kmer_counts.loc[:,"bin"]=kmer_counts.loc[:,"name_1"].apply(lambda x: x.split("_")[4]+"_"+x.split("_")[5])
        
        
    bins = list(set(kmer_counts.loc[:,"bin_1"]))
    
    # define last bin
    utr_3_elements = [item for item in bins if item.startswith('3UTR')]
    numbers = [int(item.split('_')[1]) for item in utr_3_elements]
    max_index = numbers.index(max(numbers))
    result = utr_3_elements[max_index]
    
    print("Kmers length: ",k)
    #print("BIN EXPLORED:",result)
    
    ### formula per escludere i bin finali non lunghi 7
    kmer_counts['tot_sum'] = kmer_counts.apply(
        lambda row: row['length_bin'] if row["bin_1"] != result else row['length_bin'] - (k - 1),
        axis=1
    )
    ### se il bin minore di k diventa 0 tot
    kmer_counts.loc[kmer_counts["tot_sum"]<0,"tot_sum"]=0

    kmer_counts[['tot_sum','coverage']] = kmer_counts[['tot_sum','coverage']].astype(int)
    kmer_counts_tot_grouped = kmer_counts.groupby("bin_1")['tot_sum'].sum()
    kmer_counts_tot_grouped = pd.DataFrame(kmer_counts_tot_grouped)
    
    
    kmer_counts_filtered = kmer_counts[kmer_counts["start_2"] != -1]
    
    kmer_counts_grouped = kmer_counts_filtered.groupby("bin_1")['coverage'].sum()
    kmer_counts_cov_grouped = pd.DataFrame(kmer_counts_grouped)
    
    kmer_counts = pd.merge(kmer_counts_tot_grouped,kmer_counts_cov_grouped,left_index=True,right_index=True,how='outer').fillna(0).reset_index()
    kmer_counts[['tot_sum','coverage']] = kmer_counts[['tot_sum','coverage']].astype(int)
    
    #print(kmer_counts.loc[kmer_counts["bin_1"]=="3UTR_150","tot_sum"])

    # Rename the columns for clarity
    kmer_counts = kmer_counts.rename(columns={
        'coverage': 'yes'
    })
    kmer_counts['no'] = kmer_counts['tot_sum'] - kmer_counts['yes']
    kmer_counts = kmer_counts[['bin_1','yes','no']].set_index('bin_1')
    
    #print("ERRORE yes:",kmer_counts.loc[kmer_counts["yes"]<0])
    #print("ERRORE no:",kmer_counts.loc[kmer_counts["no"]<0])
    
    return kmer_counts

def DfScoresByBed(bed_df,file_score,type_extraction="mean",cores=1):

    bed_df.loc[:,"file"]=file_score
    bed_df.loc[:,"type_extraction"]=type_extraction
    
    if cores > 1:
        print("used cores",cores)
        with mp.Pool(cores) as pool:
            bed_df.loc[:,'score'] = pool.starmap(RetriveScores, zip(bed_df['chrom'].to_list(),
                                                                    bed_df['start'].to_list(),
                                                                    bed_df['end'].to_list(),
                                                                    bed_df['file'].to_list(),
                                                                    bed_df['type_extraction'].to_list()))
    else:
        print("no multicore",cores)
        bed_df.loc[:,"score"]=bed_df.apply(lambda row: RetriveScores(row["chrom"],row["start"],row["end"],row["file"],row["type_extraction"]),axis=1)
    
    del bed_df["file"]
    del bed_df["type_extraction"]
    return bed_df
    
def MapperTxToGenome_Df(bed_df,python_bin,mapper_specie_release,output_dir,chunk_size="none",file_name="bed_to_map.bed"):
    mapper_dir = os.environ.get("KTOOLS_MAPPER_DIR")
    script_bin = os.environ.get("KTOOLS_MAPPER_SCRIPT")
    if not mapper_dir or not script_bin:
        raise ValueError("Set KTOOLS_MAPPER_DIR and KTOOLS_MAPPER_SCRIPT for genome mapping")
    mapper_dir = os.path.join(mapper_dir, "")
    
    if os.path.exists(output_dir)==False:
        os.makedirs(output_dir)
        
    file_bed=output_dir+file_name
    
    
    
    if chunk_size == "none":
        
        bed_df.to_csv(file_bed,sep="\t",index=None,header=None)

        #print("mappers available:\n")


        lista_mapper=os.listdir(mapper_dir)
        #for i in lista_mapper:
        #    print(i)

        # mapper file
        mapper_file=mapper_dir+"parsed_genomic_"+mapper_specie_release+".txt"

        # launching

        print(python_bin+" "+script_bin+" "+file_bed+" "+mapper_file+" "+output_dir)
        os.system(python_bin+" "+script_bin+" "+file_bed+" "+mapper_file+" "+output_dir)

        out_mapper=pd.read_table(output_dir+file_name.replace(".bed","")+".mapped.bed",header=None)
        out_mapper.columns=["chrom","start","end","name","score","strand"]
    
    else:
        
        chunks_list=Chunks(len(bed_df),chunk_size)
        print()
        print(chunks_list)
        print()
        
        count=0
        
        for chunk_i in chunks_list:
            count+=1
            
            print("chunked processing:",chunk_i,count/len(chunks_list))
            bed_df.loc[chunk_i[0]:chunk_i[1]].to_csv(file_bed,sep="\t",index=None,header=None)

            #print("mappers available:\n")


            lista_mapper=os.listdir(mapper_dir)
            #for i in lista_mapper:
            #    print(i)

            # mapper file
            mapper_file=mapper_dir+"parsed_genomic_"+mapper_specie_release+".txt"

            # launching


            os.system(python_bin+" "+script_bin+" "+file_bed+" "+mapper_file+" "+output_dir)

            out_mapper_i=pd.read_table(output_dir+file_name.replace(".bed","")+".mapped.bed",header=None)
            out_mapper_i.columns=["chrom","start","end","name","score","strand"]
            
            if chunk_i == chunks_list[0]:
        
                out_mapper=out_mapper_i
            
            else:
                
                out_mapper=out_mapper.append(out_mapper_i)
            
        
    #out_mapper.loc[:,"bed_type"]=out_mapper.loc[:,"name"].apply(lambda x: x.split("_")[-1])
    #out_mapper.loc[:,"name"]=out_mapper.loc[:,"name"].apply(lambda x: x.split("__")[0])
    out_mapper=out_mapper.reset_index()
    del out_mapper["index"]
    return out_mapper



def ChunksByNBins(lunghezza,n_bins,add_start=0):
    
    # devo generare TOT bins
    step=lunghezza/n_bins
    a=0
    all_chunks=[[0,step]]
    for i in range(0,n_bins-1):

        a=a+step

        all_chunks+=[[a,a+step]]
    
    chunks_df=pd.DataFrame(all_chunks)
    chunks_df.columns=["start","end"]
    
    # la coordinata puntiforme è localizzata nel trascritto mentre
    # i bin sono fatti per regione. Aggiungendo la grandezza delle porzioni 
    # che stanno prima permette di localizzare i bin nel trascritto e di
    # risalire alla posizione della coordinata
    
    chunks_df.loc[:,"start"]=chunks_df.loc[:,"start"]+add_start
    chunks_df.loc[:,"end"]=chunks_df.loc[:,"end"]+add_start
    chunks_df.loc[list(chunks_df.index)[-1],"end"]=lunghezza+add_start
    chunks_df.index=chunks_df.index+1
    return chunks_df
    
def RegionBinner(tx_i,region_name,size_region,n_bins,add_start=0):

    bed_range_tx=ChunksByNBins(size_region,n_bins,add_start=add_start)
    bed_range_tx.loc[:,"start"]=bed_range_tx.loc[:,"start"].apply(lambda x: round(x))
    bed_range_tx.loc[:,"end"]=bed_range_tx.loc[:,"end"].apply(lambda x: round(x))
    bed_range_tx=bed_range_tx.reset_index()
    bed_range_tx.loc[:,"name"]=region_name+"_"+bed_range_tx.loc[:,"index"].astype(str)
    bed_range_tx.loc[:,"chr"]=tx_i
    bed_range_tx=bed_range_tx[["chr","start","end","name"]]
    
    bed_range_tx=bed_range_tx.loc[bed_range_tx["start"]!=bed_range_tx["end"]]
    return bed_range_tx

def RegionSizeByTxList(tx_list,annfile,type_proportion="mean"):

    tx_df=pd.DataFrame(tx_list).drop_duplicates()
    tx_df.columns=["tx"]

    ann_df=pd.read_table(annfile)

    print("input transcripts:",len(tx_df))
    tx_info=pd.merge(tx_df,ann_df,left_on="tx",right_on="ensembl_transcript_id")
    tx_info=tx_info.loc[tx_info["transcript_biotype"]=="protein_coding"]
    print("input transcripts identified as protein coding transcripts:",len(tx_info))

    #tx_info.loc[:,"5utr"].astype(int).plot.box(showfliers=False)
    #tx_info.loc[:,"cds"].astype(int).plot.box(showfliers=False)
    #tx_info.loc[:,"3utr"].astype(int).plot.box(showfliers=False)

    ### in questo caso prendo in considerazione la grandezza delle regioni per
    ### poi costruire un modello porporzionato
    tx_info[["5utr","cds","3utr","UTR5","CDS","UTR3"]] = tx_info[["5utr","cds","3utr","UTR5","CDS","UTR3"]].astype(int)
    print(tx_info.loc[:,["5utr","cds","3utr"]].plot.box(showfliers=False))

    if type_proportion == "mean":
    
        param_5utr=tx_info.loc[:,"5utr"].astype(int).mean()
        param_cds=tx_info.loc[:,"cds"].astype(int).mean()
        param_3utr=tx_info.loc[:,"3utr"].astype(int).mean()
        
    if type_proportion == "median":
        
        param_5utr=tx_info.loc[:,"5utr"].astype(int).median()
        param_cds=tx_info.loc[:,"cds"].astype(int).median()
        param_3utr=tx_info.loc[:,"3utr"].astype(int).median()

    total_param_regions=param_5utr+param_cds+param_3utr

    frac_5utr=(param_5utr/total_param_regions)*100
    frac_cds=(param_cds/total_param_regions)*100
    frac_3utr=(param_3utr/total_param_regions)*100
    
    print()
    print("fraction 5'UTR:",frac_5utr)
    print("fraction CDS:",frac_cds)
    print("fraction 3'UTR:",frac_3utr)

    return [frac_5utr,frac_cds,frac_3utr]

def FindAll(string,substring):
    res = [i for i in range(len(string)) if string.startswith(substring, i)]
    return res

def kmersInTranscripts(tx_ids: list, ann_df, kmers_list: list):

    if isinstance(ann_df, str):
        ann_df = pd.read_csv(ann_df, sep='\t')

    # Filter for relevant transcripts
    tx_ann_df = ann_df.loc[ann_df["ensembl_transcript_id"].isin(tx_ids)].copy()

    if len(kmers_list) == 0:
        raise ValueError("kmers_list is empty")

    k = len(kmers_list[0])

    # Scan each kmer in each transcript using FindAll
    for i, kmer in enumerate(kmers_list):
        print(f"Scanning kmer {i+1}/{len(kmers_list)}: {kmer}", end="\r")
        tx_ann_df[kmer] = tx_ann_df["cdna"].apply(lambda x: FindAll(x, kmer))

    # Melt to long format
    tx_ann_df = pd.melt(
        tx_ann_df,
        id_vars=['ensembl_transcript_id'],
        value_vars=kmers_list,
        var_name='name',
        value_name='start')

    # Keep only positions with hits
    tx_ann_df = tx_ann_df[tx_ann_df['start'].apply(lambda x: len(x) > 0)]
    tx_ann_df = tx_ann_df.explode("start")

    # Generate BED-like columns
    tx_ann_df["end"] = tx_ann_df["start"] + k
    tx_ann_df["strand"] = "+"
    tx_ann_df["score"] = 1
    tx_ann_df = tx_ann_df.rename(columns={"ensembl_transcript_id": "chrom"})
    tx_ann_df = tx_ann_df[["chrom", "start", "end", "name", "score", "strand"]]

    # Reset index
    tx_ann_df = tx_ann_df.reset_index(drop=True)

    print(f"\nN transcripts scanned: {len(set(tx_ids))}")
    print(f"N kmers found: {len(tx_ann_df)}")
    print(f"N kmer types: {len(set(kmers_list))}")

    return tx_ann_df


def NtFrequencyByDf(df_nt,mode="dna"):

        positions=[]
        A_ratio=[]
        T_ratio=[]
        C_ratio=[]
        G_ratio=[]

        for pos in df_nt.columns:

            tot=len("".join(df_nt.loc[:,pos]))

            if tot > 0:

                A_count="".join(df_nt.loc[:,pos]).count("A")
                if mode =="dna":
                    T_count="".join(df_nt.loc[:,pos]).count("T")
                if mode == "rna":
                    T_count="".join(df_nt.loc[:,pos]).count("U")
                C_count="".join(df_nt.loc[:,pos]).count("C")
                G_count="".join(df_nt.loc[:,pos]).count("G")


                positions+=[pos]
                A_ratio+=[A_count/tot]
                T_ratio+=[T_count/tot]
                C_ratio+=[C_count/tot]
                G_ratio+=[G_count/tot]

        perc_nt=pd.DataFrame([A_ratio,T_ratio,C_ratio,G_ratio]).T
        
        if mode == "dna":
        
            perc_nt.columns=["A","T","C","G"]
        if mode == "rna":
            perc_nt.columns=["A","U","C","G"]
        

        return perc_nt
    
def SeqLogoByList(lista_kmers,mode = "dna"):
        
        nt_df=lista_kmers
        nt_df=pd.DataFrame(nt_df)
        nt_df.iloc[:,0] = nt_df.iloc[:,0].str.rstrip()
        nt_df=nt_df.iloc[:,0].apply(lambda x: pd.Series(list(x)))

        perc_nt=NtFrequencyByDf(nt_df,mode)

        crp_logo=logomaker.Logo(perc_nt,
                            shade_below=.5,
                            fade_below=.5,
                            font_name='Arial Rounded MT Bold')


        # style using Logo methods
        crp_logo.style_spines(visible=False)
        crp_logo.style_spines(spines=['left', 'bottom'], visible=True)
        crp_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

        # style using Axes methods
        crp_logo.ax.set_ylabel("Nt freqeuncy in 7mers", labelpad=-1)
        crp_logo.ax.xaxis.set_ticks_position('none')
        crp_logo.ax.xaxis.set_tick_params(pad=-1)
        
        print(crp_logo)  
        
def get_fasta_ids(fasta_path):
    return [record.id for record in SeqIO.parse(fasta_path, "fasta")]

def kmers_counter(seq, k, all_kmers_df=None, tx_id=None):
    if all_kmers_df is not None and not isinstance(all_kmers_df, pd.DataFrame):
        all_kmers_df = pd.DataFrame(all_kmers_df)

    kmers = [seq[i:i+k] for i in range(0, len(seq)-k)]
    uniq, counts = np.unique(kmers, return_counts=True)
    conte_kmers = pd.DataFrame({0: uniq, 1: counts})

    if all_kmers_df is not None:
        conte_kmers = (
            pd.merge(all_kmers_df, conte_kmers, on=0, how="left")
            .fillna(0)
        )

    conte_kmers = conte_kmers.set_index(0)

    if tx_id is not None:
        conte_kmers = conte_kmers.rename(columns={1: tx_id})

    return conte_kmers
    
def process_kmer_counts(kmers_df, label):
    kmers_df['sum'] = kmers_df.sum(axis=1)
    kmers_df = kmers_df[['sum']].reset_index()
    kmers_df.columns = ['kmer', 'count']
    df_repeated = kmers_df.loc[kmers_df.index.repeat(kmers_df['count'])].reset_index(drop=True)
    df_repeated = df_repeated.drop(columns=['count'])
    df_repeated['label'] = label
    return df_repeated

def run_stat_test(x, y, test):
    # Convert inputs to numeric arrays (float)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    # Skip empty arrays
    if x.size == 0 or y.size == 0:
        return np.nan, np.nan

    # Skip constant arrays (all values identical)
    if np.all(x == x[0]) and np.all(y == y[0]):
        return np.nan, np.nan

    if test == "ks":
        stat, pval = ks_2samp(x, y, alternative="two-sided")
    elif test == "mwu":
        stat, pval = mannwhitneyu(x, y, alternative="two-sided")
    elif test == "ttest":
        stat, pval = ttest_ind(x, y, equal_var=False, nan_policy="omit")
    else:
        raise ValueError(f"Unknown stat test: {test}")

    return stat, pval



class KEA():
    
    def __init__(self, dir_out, ref_fasta, ctrl_fasta):
        
        print("Motifs object")
        
        # Kmer-related attributes
        self.kmer_list = None
        self.kmers_count_table_by_species = {}  # dict to store species-specific kmer tables
        self.kmers_zscore = None
        self.kmer_bed = None
        self.enrichment = None
        self.score_matrix = None
        self.kmers_frequencies = None
        self.kmer_cluster_map = None
        
        # Directories
        os.makedirs(dir_out, exist_ok=True)
        self.dir_out = dir_out
        
        # Fasta inputs
        self.reference_fasta = ref_fasta
        self.control_fasta = ctrl_fasta
        
        # Combined FASTA
        self.combined_fasta = None  # default combined fasta
        self.combined_fasta_by_species = {}  # dict for species-specific combined FASTAs
        
        # Results
        self.KEA_results = None
        from importlib.resources import files
        self.seed_family_file = str(files("ktools_resources").joinpath("miR_Family_Info.txt"))
        
    def get_kmer_table_for_species(self, species):
        """
        Returns kmer table, reference IDs, control IDs for a given species.

        Parameters
        ----------
        species : str
            Species key used in combined FASTA / kmer table.

        Returns
        -------
        kmer_table : pd.DataFrame
        ref_group : list of str
        comp_group : list of str
        """
        # Species-specific kmer table
        if species not in self.kmers_count_table_by_species:
            raise ValueError(f"Kmer table for species '{species}' not found. Run KmersCountsTable first.")

        kmer_table = self.kmers_count_table_by_species[species]

        # Reference and control transcript IDs
        if species == "input" or species is None:
            ref_group = get_fasta_ids(self.reference_fasta)
            comp_group = get_fasta_ids(self.control_fasta)
        else:
            # Ortholog species: assume control_fasta contains ortholog transcripts
            ref_group = get_fasta_ids(self.reference_fasta)
            comp_group = get_fasta_ids(self.combined_fasta_by_species[species])

        return kmer_table, ref_group, comp_group
    
    
    
    def CreateCombinedFasta(self, species=None):
        output_dir = self.dir_out
        if not hasattr(self, "combined_fasta_by_species"):
            self.combined_fasta_by_species = {}

        # Output file path
        output_fasta = os.path.join(output_dir, 'combined_sequences.fa' if species is None else f'combined_sequences_{species}.fa')

        # Determine input FASTAs
        if species is None:
            ref_fasta = self.reference_fasta
            ctrl_fasta = self.control_fasta
        else:
            # Look up species-specific FASTAs created by FetchOrthologsAndSaveFasta
            ref_key = f"{species}_ref"
            ctrl_key = f"{species}_ctrl"
            ref_fasta = self.combined_fasta_by_species.get(ref_key)
            ctrl_fasta = self.combined_fasta_by_species.get(ctrl_key)

            if ref_fasta is None or ctrl_fasta is None:
                raise ValueError(f"Species-specific FASTAs not found for {species}. Run FetchOrthologsAndSaveFasta first.")

        # Combine sequences
        with open(output_fasta, 'w') as outfile:
            for fasta in [ref_fasta, ctrl_fasta]:
                with open(fasta, 'r') as infile:
                    for line in infile:
                        outfile.write(line)

        # Store path
        if species is None:
            self.combined_fasta = output_fasta
        else:
            self.combined_fasta_by_species[species] = output_fasta

    def KmersCountsTable(self, k, cores=5, subset_by_subseq=None, species=None,
                         *, cache=True, force=False):
        """Build the complete k-mer count table, automatically saving/reusing it.

        By default write a compressed pandas table and JSON manifest inside
        self.dir_out (the directory supplied when constructing KEA). Reuse
        requires matching FASTA content, k, species, subsequence filter,
        counter implementation and numpy/pandas versions. Core count does
        not affect the cache key. force=True recounts; cache=False disables
        disk caching. Counting semantics and the existing None return are
        unchanged. The table remains in kmers_count_table_by_species.
        """
        import hashlib
        import inspect
        import json
        import pickle
        import tempfile

        species_key = "input" if species is None else species
        if species is None:
            INPUT = self.combined_fasta
        else:
            INPUT = self.combined_fasta_by_species.get(species)
        if INPUT is None:
            existing = Path(self.dir_out) / (
                'combined_sequences.fa' if species is None else f'combined_sequences_{species}.fa')
            if existing.is_file():
                INPUT = str(existing)
                if species is None:
                    self.combined_fasta = INPUT
                else:
                    self.combined_fasta_by_species[species] = INPUT
            else:
                raise ValueError(f"Combined FASTA for species '{species_key}' not found. Run CreateCombinedFasta first.")
        input_path = Path(INPUT)
        if not input_path.is_file():
            raise FileNotFoundError(input_path)

        # Preserve the original order and duplicate-filter semantics.
        filters = list(subset_by_subseq) if subset_by_subseq else None
        all_kmers = list(map("".join, product("ATCG", repeat=k)))
        if filters:
            all_kmers = [word for subseq in filters for word in all_kmers if subseq in word]

        def digest(path):
            sha = hashlib.sha256()
            with Path(path).open('rb') as stream:
                for block in iter(lambda: stream.read(1024*1024), b''):
                    sha.update(block)
            return sha.hexdigest()

        if not hasattr(self, 'kmers_count_table_by_species'):
            self.kmers_count_table_by_species = {}
        if not hasattr(self, 'kmers_count_cache'):
            self.kmers_count_cache = {}
        identity = cache_path = metadata_path = None
        if cache:
            identity = dict(format_version=1, k=k, species=species_key,
                            subset_by_subseq=filters, input_sha256=digest(input_path),
                            counter_sha256=hashlib.sha256(inspect.getsource(kmers_counter).encode()).hexdigest(),
                            numpy_version=np.__version__, pandas_version=pd.__version__)
            slug = ''.join(c if c.isalnum() or c in '._-' else '_' for c in str(species_key))
            if slug != str(species_key):
                slug += '_'+hashlib.sha256(str(species_key).encode()).hexdigest()[:12]
            suffix = '' if filters is None else '_subset_'+hashlib.sha256(
                json.dumps(filters, ensure_ascii=True).encode()).hexdigest()[:16]
            cache_path = Path(self.dir_out)/f'kmers_counts_{slug}_k{k}{suffix}.pkl.gz'
            metadata_path = cache_path.with_suffix(cache_path.suffix+'.json')
            if not force and cache_path.is_file() and metadata_path.is_file():
                try:
                    metadata = json.loads(metadata_path.read_text())
                    if metadata['identity'] == identity:
                        if digest(cache_path) != metadata['table_sha256']:
                            raise ValueError('cache checksum mismatch')
                        table = pd.read_pickle(cache_path, compression='gzip')
                        if (not isinstance(table, pd.DataFrame)
                                or table.index.tolist() != all_kmers
                                or table.columns.tolist() != metadata['columns']
                                or list(table.shape) != metadata['shape']
                                or not table['k'].eq(k).all()):
                            raise ValueError('cache schema mismatch')
                        if digest(input_path) != identity['input_sha256']:
                            raise RuntimeError('Input FASTA changed while loading the cache.')
                        self.kmers_count_table_by_species[species_key] = table
                        self.kmers_count_cache[species_key] = dict(
                            state='reused', table=str(cache_path), metadata=str(metadata_path))
                        print(f"[KEA] Loaded cached k-mer table: {cache_path}")
                        return
                    print('[KEA] K-mer cache differs from the current inputs or parameters; recounting.')
                except (OSError, ValueError, KeyError, TypeError, EOFError,
                        pickle.UnpicklingError, ImportError, AttributeError) as exc:
                    print(f'[KEA] Invalid k-mer cache; recounting ({type(exc).__name__}).')

        print(f"[KEA] Counting kmers (k={k}) for species '{species_key}' | N kmers: {len(all_kmers)}")
        records = list(SeqIO.parse(INPUT, "fasta"))
        input_id = [rec.id for rec in records]
        input_seq = [str(rec.seq) for rec in records]
        all_k = [k] * len(input_seq)
        all_kmers_list = [all_kmers] * len(input_seq)
        with mp.Pool(cores) as pool:
            kmers_list = pool.starmap(kmers_counter, zip(input_seq, all_k, all_kmers_list, input_id))
        kmer_counts = pd.concat(kmers_list, axis=1)
        kmer_counts["k"] = k

        if cache:
            if digest(input_path) != identity['input_sha256']:
                raise RuntimeError('Input FASTA changed during counting; no cache was saved.')
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            table_temp = metadata_temp = None
            try:
                with tempfile.NamedTemporaryFile(dir=cache_path.parent,
                        prefix='.'+cache_path.name+'.', suffix='.tmp', delete=False) as stream:
                    table_temp = Path(stream.name)
                kmer_counts.to_pickle(table_temp, compression={'method':'gzip', 'compresslevel':1})
                metadata = dict(identity=identity, input_fasta=str(input_path.resolve()),
                                table_sha256=digest(table_temp), shape=list(kmer_counts.shape),
                                columns=kmer_counts.columns.tolist())
                with tempfile.NamedTemporaryFile(mode='w', dir=cache_path.parent,
                        prefix='.'+metadata_path.name+'.', suffix='.tmp', delete=False) as stream:
                    metadata_temp = Path(stream.name)
                    json.dump(metadata, stream, indent=2)
                    stream.write('\n')
                os.replace(table_temp, cache_path)
                os.replace(metadata_temp, metadata_path)
            finally:
                for temporary in [table_temp, metadata_temp]:
                    if temporary is not None:
                        temporary.unlink(missing_ok=True)
            self.kmers_count_cache[species_key] = dict(
                state='computed_saved', table=str(cache_path), metadata=str(metadata_path))
            print(f'[KEA] Saved k-mer table: {cache_path}')
        else:
            self.kmers_count_cache[species_key] = dict(state='computed_cache_disabled', table=None, metadata=None)
        self.kmers_count_table_by_species[species_key] = kmer_counts

    def KRS(
        self,
        reference,
        species=None,
        top_pct=5,
        store=True,
        save=True,
        save_full_rank_table=False,
        verbose=True):
        """
        K-mer Rank Signature (KRS)
    
        Parameters
        ----------
        reference : str
            Column name corresponding to the target transcript.
    
        species : str, optional
            Species key. If None, "input" is used.
    
        top_pct : float, default=5
            Select k-mers for which the target transcript falls within
            the specified upper percentile. For example, top_pct=5
            selects target percentile ranks >= 0.95.
    
        store : bool, default=True
            Store results in self.KEA_results.
    
        save : bool, default=True
            Save the signature and compact ranking table to disk.
    
        save_full_rank_table : bool, default=False
            Save the complete k-mer-by-transcript percentile-rank table.
            This table can be extremely large for transcriptome-scale
            backgrounds.
    
        verbose : bool, default=True
            Print progress and output paths.
    
        Returns
        -------
        list
            K-mers belonging to the target transcript signature.
        """
    
        if not 0 < top_pct <= 100:
            raise ValueError("top_pct must be greater than 0 and <= 100.")
    
        species_key = species if species is not None else "input"
    
        # --------------------------------------------------
        # Load k-mer count table
        # --------------------------------------------------
        if species_key not in self.kmers_count_table_by_species:
            raise ValueError(
                f"K-mer table for species '{species_key}' not found. "
                "Run KmersCountsTable first."
            )
    
        kmers_count_table = (
            self.kmers_count_table_by_species[species_key]
            .copy()
            .drop(columns=["k"], errors="ignore")
        )
    
        if reference not in kmers_count_table.columns:
            raise ValueError(
                f"Reference '{reference}' not found in the k-mer table."
            )
    
        if verbose:
            print(
                f"[KEA] KRS | species={species_key} | "
                f"reference={reference}"
            )
    
        # --------------------------------------------------
        # Calculate per-transcript k-mer frequencies
        # --------------------------------------------------
        transcript_totals = kmers_count_table.sum(axis=0)
    
        # Avoid division by zero for empty transcripts
        transcript_totals = transcript_totals.replace(0, np.nan)
        kmers_freq_table = kmers_count_table.div(transcript_totals,axis=1)
    
        # --------------------------------------------------
        # Calculate percentile ranks
        # --------------------------------------------------
        # Transcripts become rows and k-mers become columns
        frequencies_by_transcript = kmers_freq_table.T
    
        # Exclude absent k-mers from percentile calculation
        frequencies_by_transcript = (frequencies_by_transcript.replace(0, np.nan))
    
        # For each k-mer, rank its frequency across transcripts
        ranks_by_transcript = frequencies_by_transcript.rank(pct=True,axis=0)
    
        # Return to k-mers as rows and transcripts as columns
        kmers_rankpct_table = ranks_by_transcript.T
    
        # --------------------------------------------------
        # Extract the KRS signature
        # --------------------------------------------------
        threshold = 1 - (top_pct / 100)
        target_rank = kmers_rankpct_table[reference]
        target_frequency = kmers_freq_table[reference]
        selected_mask = target_rank >= threshold
        kmers_signature = (target_rank[selected_mask].sort_values(ascending=False).index.tolist())
    
        # --------------------------------------------------
        # Create a compact output table
        # --------------------------------------------------
        kmer_ranks = pd.DataFrame({
            "kmer": kmers_rankpct_table.index,
            "target_frequency": target_frequency,
            "target_rank_percentile": target_rank,
            "selected_KRS": selected_mask
        })
    
        kmer_ranks = (
            kmer_ranks
            .sort_values(
                by=[
                    "selected_KRS",
                    "target_rank_percentile",
                    "target_frequency"
                ],
                ascending=[False, False, False]
            )
            .reset_index(drop=True)
        )
    
        if verbose:
            print(
                f"[KEA] Extracted {len(kmers_signature)} k-mers "
                f"with target percentile rank >= {threshold:.3f}"
            )
    
        # --------------------------------------------------
        # Save results
        # --------------------------------------------------
        output_dir = None
    
        if save:
            # Make the reference name safe for use as a folder name
            safe_reference = "".join(
                character
                if character.isalnum() or character in "._-"
                else "_"
                for character in str(reference)
            )
    
            output_dir = os.path.join(
                self.dir_out,
                "ExtractKmers",
                species_key,
                "KRS",
                safe_reference
            )
    
            os.makedirs(output_dir, exist_ok=True)
    
            # One k-mer per line
            FileByList(
                os.path.join(output_dir, "signature.txt"),
                kmers_signature
            )
    
            # Compact table with target frequency and percentile rank
            kmer_ranks.to_csv(
                os.path.join(output_dir, "kmer_ranks.tsv"),
                sep="\t",
                index=False
            )
    
            # Optional complete rank matrix
            if save_full_rank_table:
                kmers_rankpct_table.to_csv(
                    os.path.join(
                        output_dir,
                        "full_rank_table.tsv"
                    ),
                    sep="\t",
                    index=True,
                    index_label="kmer"
                )
    
            if verbose:
                print(f"[KEA] KRS results saved in: {output_dir}")
    
        # --------------------------------------------------
        # Store results in the object
        # --------------------------------------------------
        if store:
            if (
                not hasattr(self, "KEA_results") or
                self.KEA_results is None
            ):
                self.KEA_results = {}
    
            self.KEA_results.setdefault(species_key, {})
    
            self.KEA_results[species_key]["KRS"] = {
                "reference": reference,
                "top_pct": top_pct,
                "threshold": threshold,
                "signature": kmers_signature,
                "kmer_ranks": kmer_ranks,
                "rank_table": kmers_rankpct_table,
                "output_dir": output_dir
            }
    
        return kmers_signature

    def KmerSlidingWindowProfile(
    self,
    kmers,
    tx_id,
    ann_df,
    fasta_file=None, 
    kmer_groups=None,
    window_size=50,
    step=25,
    strand="+",
    n_clusters="auto",
    cluster_range=(2, 10),
    convert_to_rna=True,
    plot_line=True,
    plot_heatmap=True,
    plot_cluster_heatmap=True,
    plot_silhouette=True,
    save_prefix=None,
    verbose=True,
    improvement_threshold=0.10):

    
        # ------------------------------------------------------------------
        # 0. Load annotation OR FASTA (ONLY NEW PART)
        # ------------------------------------------------------------------
        if fasta_file is not None:
    
            if verbose:
                print(f"[KEA] Loading sequence from FASTA: {fasta_file}")
    
            seq_list = []
            current_id = None
            current_seq = []
    
            with open(fasta_file) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        if current_id == tx_id and current_seq:
                            seq_list = ["".join(current_seq)]
                            break
                        current_id = line[1:].split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line)
    
                if current_id == tx_id and current_seq:
                    seq_list = ["".join(current_seq)]
    
            if len(seq_list) == 0:
                raise ValueError(f"Transcript '{tx_id}' not found in FASTA.")
    
        else:
            if isinstance(ann_df, str):
                ann_df = pd.read_table(ann_df, sep="\t")
    
            id_col = "ID" if "ID" in ann_df.columns else "ensembl_transcript_id"
            tx_row = ann_df[ann_df[id_col] == tx_id]
            if len(tx_row) == 0:
                raise ValueError(f"Transcript '{tx_id}' not found in ann_df column '{id_col}'.")
    
            seq_list = list(tx_row["cdna"])
    
        if verbose:
            print(f"[KEA] KmerSlidingWindowProfile | tx={tx_id} | "
                  f"windows: size={window_size}, step={step} | "
                  f"kmers={len(kmers) if kmer_groups is None else 'grouped'} | convert_to_rna={convert_to_rna}")
    
        # ------------------------------------------------------------------
        # 1. Sliding window BED-like dataframe
        # ------------------------------------------------------------------
        records = []
        for seq in seq_list:
            for i in range(0, len(seq) - window_size + 1, step):
                records.append([tx_id, i, i + window_size, seq[i:i + window_size], 0, strand])
    
        df = pd.DataFrame(records, columns=["chrom", "start", "end", "name", "score", "strand"])
        if verbose:
            print(f"[KEA] Generated {len(df)} windows.")
    
        # ------------------------------------------------------------------
        # 2. T -> U conversion
        # ------------------------------------------------------------------
        if kmer_groups is not None:
            working_kmers = [
                k.replace("T", "U") if convert_to_rna else k
                for group in kmer_groups
                for k in group
            ]
            group_sizes = [len(g) for g in kmer_groups]
            group_boundaries = np.cumsum(group_sizes)
        
        else:
            working_kmers = (
                [k.replace("T", "U") for k in kmers]
                if convert_to_rna
                else list(kmers)
            )
        
        df["seq"] = df["name"].str.upper()
        if convert_to_rna:
            df["seq"] = df["seq"].str.replace("T", "U")
    
        # ------------------------------------------------------------------
        # 3. Vectorised kmer counting
        # ------------------------------------------------------------------
        def _count_kmers_in_windows(seqs, kmers_list):
            results_counts, results_freqs = [], []
            total_counts, combined_freqs = [], []
    
            for seq in seqs:
                seq_arr = np.array(list(seq))
                counts, freqs = {}, {}
                total_count = 0
                total_positions = 0
    
                for kmer in kmers_list:
                    k = len(kmer)
                    positions = len(seq) - k + 1
                    if positions <= 0:
                        counts[kmer] = 0
                        freqs[kmer] = 0.0
                        continue
                    kmer_arr = np.array(list(kmer))
                    windows_view = np.lib.stride_tricks.sliding_window_view(seq_arr, k)
                    count = int(np.all(windows_view == kmer_arr, axis=1).sum())
                    counts[kmer] = count
                    freqs[kmer] = count / positions
                    total_count += count
                    total_positions += positions
    
                combined = total_count / total_positions if total_positions > 0 else 0.0
                results_counts.append(counts)
                results_freqs.append(freqs)
                total_counts.append(total_count)
                combined_freqs.append(combined)
    
            return (pd.DataFrame(results_counts),
                pd.DataFrame(results_freqs),
                total_counts,
                combined_freqs)
    
        counts_df, freqs_df, total_counts, combined_freqs = _count_kmers_in_windows(df["seq"].tolist(), working_kmers)
    
        df = pd.concat([df, counts_df.add_prefix("count_"), freqs_df.add_prefix("freq_")], axis=1)
        df["total_kmer_count"] = total_counts
        df["combined_freq"] = combined_freqs
    
        # ------------------------------------------------------------------
        # 4. Kmer x bin frequency matrix
        # ------------------------------------------------------------------
        df["bin"] = df["chrom"] + "_" + df["start"].astype(str) + "_" + df["end"].astype(str)
        df_indexed = df.set_index("bin")
        freq_cols = [c for c in df_indexed.columns if c.startswith("freq_") and "combined" not in c]
        df_matrix = df_indexed[freq_cols].T
        df_matrix.index = df_matrix.index.str.replace("freq_", "")
    
        # ------------------------------------------------------------------
        # 5. Hierarchical clustering
        # ------------------------------------------------------------------
        Z = linkage(df_matrix.values, method="ward", metric="euclidean")
    
        # ------------------------------------------------------------------
        # 6. Cluster number selection
        # ------------------------------------------------------------------
        if n_clusters == "auto":
        
            k_min, k_max = cluster_range
            k_max = min(k_max, len(df_matrix) - 1)
            k_min = max(k_min, 2)
        
            if k_min >= k_max:
                best_k = 2
                sil_scores = {}
            else:
                sil_scores = {}
                best_k = None
                best_score = -1
        
                for k in range(k_min, k_max + 1):
                    labels = fcluster(Z, k, criterion="maxclust")
                    if len(np.unique(labels)) < 2:
                        continue
        
                    score = silhouette_score(df_matrix.values, labels)
                    sil_scores[k] = score
        
                    if score > best_score * (1 + improvement_threshold):
                        best_score = score
                        best_k = k
        
                if best_k is None:
                    best_k = max(sil_scores, key=sil_scores.get)
        
            if plot_silhouette and sil_scores:
                plt.figure(figsize=(6,3))
                plt.plot(list(sil_scores.keys()), list(sil_scores.values()), marker="o")
                plt.axvline(best_k, linestyle="--")
                plt.title(f"Silhouette — {tx_id}")
                plt.show()
        else:
            best_k = int(n_clusters)
            sil_scores = {}

        # ------------------------------------------------------------------
        # 7. Final clusters
        # ------------------------------------------------------------------
        row_cluster_labels = fcluster(Z, best_k, criterion="maxclust")
        kmer_clusters = pd.Series(row_cluster_labels, index=df_matrix.index, name="cluster")
    
        # ------------------------------------------------------------------
        # 8. Clustered heatmap
        # ------------------------------------------------------------------
        cg = sns.clustermap(
            df_matrix,
            cmap="Reds",
            row_cluster=True,
            col_cluster=False,
            linewidths=0.5,
            figsize=(14, 12)
        )

        # NEW: draw separation lines between kmer groups
        if kmer_groups is not None:
            row_order = cg.dendrogram_row.reordered_ind
            ordered_kmers = df_matrix.index[row_order]

            current_pos = 0
            boundaries = np.cumsum(group_sizes)

            for b in boundaries[:-1]:
                cg.ax_heatmap.axhline(b, color="black", linewidth=2)

        plt.suptitle(f"K-mer frequency heatmap — {tx_id} (k={best_k})", y=1.01)
    
        if plot_heatmap:
            if save_prefix:
                cg.savefig(f"{save_prefix}_heatmap.pdf", bbox_inches="tight")
            plt.show()
        else:
            plt.close()
    
        # ------------------------------------------------------------------
        # 9. Line plot (unchanged)
        # ------------------------------------------------------------------
        if plot_line:
            fig, ax = plt.subplots(figsize=(12, 3))
            ax.plot(range(len(df)), df["combined_freq"].values, linewidth=1.2)
            plt.show()
            plt.close()
    
        # ------------------------------------------------------------------
        # 10. Cluster heatmap (unchanged)
        # ------------------------------------------------------------------
        if plot_cluster_heatmap:
            df_combined=df_matrix.merge(kmer_clusters,left_index=True,right_index=True)
            df_long=(df_combined.reset_index().melt(id_vars=["index","cluster"],var_name="bin",value_name="freq").rename(columns={"index":"kmer"}))
            df_cluster_mean=df_long.groupby(["cluster","bin"])["freq"].mean().reset_index()
            df_cluster_mean["bin"]=pd.Categorical(df_cluster_mean["bin"],categories=df_matrix.columns,ordered=True)
            heatmap_data=(df_cluster_mean.pivot(index="cluster",columns="bin",values="freq")[df_matrix.columns])
            heatmap_data=heatmap_data.apply(lambda x:(x-x.mean())/(x.std()+1e-9),axis=1)
        
            fig_width=max(10,heatmap_data.shape[1]*0.25); fig_height=max(4,heatmap_data.shape[0]*0.75)
            fig,ax=plt.subplots(figsize=(fig_width,fig_height))
            sns.heatmap(heatmap_data,cmap="RdBu_r",center=0,linewidths=0.2,linecolor="lightgray",cbar_kws={"label":"Z-scored mean k-mer frequency"},ax=ax)
            ax.set_xlabel("Window bin"); ax.set_ylabel("Cluster"); ax.set_title(f"Clustered k-mer dynamics — {tx_id} (k={best_k})")
            ax.tick_params(axis="x",labelrotation=90,labelsize=7); ax.tick_params(axis="y",labelsize=8)
            plt.tight_layout()
            if save_prefix: fig.savefig(f"{save_prefix}_cluster_heatmap.pdf",bbox_inches="tight")
            plt.show(); plt.close()
    
        return df, kmer_clusters, sil_scores
    
    def ExtractKmers(
    self,
    k_selected,
    species=None,
    extraction_methods=("delta_rank_percentile", "delta_median", "stat_log2fc"),
    n_iter=1000,
    sign_emp_pval=0.0025,
    log2fc_threshold=1.0,
    fdr_alpha=0.05,
    stat_test="mwu",
    correction_method="bonferroni",
    median=True):
        """Reuse complete results when only selection thresholds change."""
        import hashlib
        import json

        species_key = species if species is not None else "input"
        self.KEA_results = self.KEA_results or {}
        state = self.KEA_results.setdefault(species_key, {})
        methods = tuple(extraction_methods)
        for method in methods:
            if method not in {"delta_rank_percentile", "delta_median", "stat_log2fc"}:
                raise ValueError(f"Unknown extraction method: {method}")
        if "stat_log2fc" in methods and correction_method not in {"bonferroni", "fdr", "PValue"}:
            raise ValueError(f"Unknown correction_method: {correction_method}")

        if species_key == "input":
            kmer_table = self.kmers_count_table_by_species["input"]
            ref_group = get_fasta_ids(self.reference_fasta)
            comp_group = get_fasta_ids(self.control_fasta)
        else:
            if species_key not in self.kmers_count_table_by_species:
                raise ValueError(f"Kmer table for species '{species_key}' not found. Run KmersCountsTable first.")
            kmer_table = self.kmers_count_table_by_species[species_key]
            ref_key, ctrl_key = f"{species_key}_ref", f"{species_key}_ctrl"
            if ref_key not in self.combined_fasta_by_species or ctrl_key not in self.combined_fasta_by_species:
                raise ValueError(f"Mapped FASTAs for species '{species_key}' not found. Run FetchOrthologsAndSaveFastaOneToOne first.")
            ref_group = get_fasta_ids(self.combined_fasta_by_species[ref_key])
            comp_group = get_fasta_ids(self.combined_fasta_by_species[ctrl_key])

        ref_group = [x for x in ref_group if x in kmer_table.columns]
        comp_group = [x for x in comp_group if x in kmer_table.columns]
        if not ref_group or not comp_group:
            raise ValueError(f"No matching transcripts found for species '{species_key}' in kmer table.")
        kmer_matrix = kmer_table.loc[kmer_table["k"] == k_selected].drop(columns="k").copy()
        base_out = Path(self.dir_out) / "ExtractKmers" / species_key
        base_out.mkdir(parents=True, exist_ok=True)

        # Thresholds and correction choice only filter results; exclude them
        # from the identity of the calculation. Track inputs to avoid stale reuse.
        context = {
            "version": 1, "species": species_key, "k": int(k_selected),
            "n_iter": int(n_iter), "stat_test": stat_test, "median": bool(median),
            "reference": ref_group, "control": comp_group,
            "columns": list(kmer_matrix.columns),
            "dtypes": [str(dtype) for dtype in kmer_matrix.dtypes],
            "counts_sha256": hashlib.sha256(
                pd.util.hash_pandas_object(kmer_matrix, index=True).values.tobytes()
            ).hexdigest(),
        }
        cached = {}
        for method in methods:
            required = {"freq_in_group_1", "freq_in_group_2"}
            required.update({"log2FoldChange", "PValue", "bonf", "FDR"}
                            if method == "stat_log2fc" else {method})
            previous = state.get(method, {})
            table = previous.get("all_kmers")
            if (previous.get("_calculation") == context
                    and isinstance(table, pd.DataFrame)
                    and required.issubset(table.columns)):
                cached[method] = table
                continue
            # Reuse a compatible saved calculation even after restarting Python.
            for metadata_path in sorted(base_out.glob(f"{method}_*/calculation.json")):
                try:
                    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
                    if metadata.get("context") != context:
                        continue
                    table_path = metadata_path.parent / "all_kmers.tsv"
                    data = table_path.read_bytes()
                    if hashlib.sha256(data).hexdigest() != metadata.get("table_sha256"):
                        continue
                    table = pd.read_csv(table_path, sep="\t", index_col=0,
                                        float_precision="round_trip")
                    if required.issubset(table.columns):
                        cached[method] = table
                        break
                except (OSError, ValueError, TypeError, AttributeError):
                    continue

        # No frequencies, control draws or statistical tests on a full cache hit.
        if any(method not in cached for method in methods):
            perc_ref = kmer_matrix[ref_group].sum(axis=1) / kmer_matrix[ref_group].sum().sum()
            if len(comp_group) >= 2 * len(ref_group):
                freq_ctrl_iter = []
                for _ in range(n_iter):
                    sampled = np.random.choice(comp_group, len(ref_group), replace=False)
                    subset = kmer_matrix[sampled]
                    freq_ctrl_iter.append(subset.sum(axis=1) / subset.sum().sum())
                perc_ctrl = pd.concat(freq_ctrl_iter, axis=1).mean(axis=1)
            else:
                perc_ctrl = kmer_matrix[comp_group].sum(axis=1) / kmer_matrix[comp_group].sum().sum()

        for method in methods:
            print(f"[KEA] {species_key} | {method}" +
                  (" | reusing all_kmers" if method in cached else ""))
            if method in cached:
                all_kmers = cached[method]
            elif method in {"delta_rank_percentile", "delta_median"}:
                if method == "delta_rank_percentile":
                    score = (perc_ref.rank(method="first") / len(perc_ref) -
                             perc_ctrl.rank(method="first") / len(perc_ctrl))
                else:
                    q1_r, med_r, q3_r = np.percentile(perc_ref, [25, 50, 75])
                    q1_c, med_c, q3_c = np.percentile(perc_ctrl, [25, 50, 75])
                    diff_ref = np.where(perc_ref < med_r,
                                        (perc_ref - med_r) / (med_r - q1_r),
                                        (perc_ref - med_r) / (q3_r - med_r))
                    diff_ctrl = np.where(perc_ctrl < med_c,
                                         (perc_ctrl - med_c) / (med_c - q1_c),
                                         (perc_ctrl - med_c) / (q3_c - med_c))
                    score = diff_ref - diff_ctrl
                all_kmers = pd.DataFrame({"freq_in_group_1": perc_ref,
                                          "freq_in_group_2": perc_ctrl,
                                          method: score}).sort_values(method, ascending=False)
            else:
                testing = kmer_matrix.T
                testing = testing.div(testing.sum(axis=1), axis=0).fillna(0)
                test_ref, test_ctrl = testing.loc[ref_group], testing.loc[comp_group]
                rows = []
                for kmer in kmer_matrix.index:
                    if kmer not in test_ref.columns:
                        continue
                    x, y = test_ref[kmer].values, test_ctrl[kmer].values
                    if len(x) == 0 or len(y) == 0:
                        continue
                    if np.all(x == x[0]) and np.all(y == y[0]):
                        continue
                    stat, pval = run_stat_test(x, y, stat_test)
                    aggregate = np.median if median else np.mean
                    label = "Median" if median else "Mean"
                    rows.append({"kmer": kmer,
                                 "log2FoldChange": np.log2((aggregate(x)+1e-9)/(aggregate(y)+1e-9)),
                                 "stat": stat, "PValue": pval,
                                 f"{label}_group1": aggregate(x),
                                 f"{label}_group2": aggregate(y)})
                all_kmers = pd.DataFrame(rows).set_index("kmer")
                _, all_kmers["bonf"], _, _ = multipletests(all_kmers["PValue"], method="bonferroni")
                _, all_kmers["FDR"], _, _ = multipletests(all_kmers["PValue"], method="fdr_bh")
                all_kmers["freq_in_group_1"] = perc_ref.loc[all_kmers.index]
                all_kmers["freq_in_group_2"] = perc_ctrl.loc[all_kmers.index]

            if method in {"delta_rank_percentile", "delta_median"}:
                top_n = max(1, int(len(all_kmers) * sign_emp_pval))
                enriched = all_kmers.head(top_n).index.tolist()
                depleted = all_kmers.tail(top_n).index.tolist()
                subdir = f"{method}_emp{sign_emp_pval}"
            else:
                p_column = {"bonferroni": "bonf", "fdr": "FDR", "PValue": "PValue"}[correction_method]
                signif = all_kmers[p_column] <= fdr_alpha
                enriched = all_kmers[(all_kmers["log2FoldChange"] >= log2fc_threshold) & signif].index.tolist()
                depleted = all_kmers[(all_kmers["log2FoldChange"] <= -log2fc_threshold) & signif].index.tolist()
                subdir = f"{method}_log2fc{log2fc_threshold}_{correction_method}_alpha{fdr_alpha}"

            outdir = base_out / subdir
            outdir.mkdir(parents=True, exist_ok=True)
            FileByList(str(outdir / "enriched.txt"), enriched)
            FileByList(str(outdir / "depleted.txt"), depleted)
            table_path = outdir / "all_kmers.tsv"
            all_kmers.to_csv(table_path, sep="\t")
            metadata = {"context": context,
                        "table_sha256": hashlib.sha256(table_path.read_bytes()).hexdigest()}
            (outdir / "calculation.json").write_text(
                json.dumps(metadata, indent=2), encoding="utf-8")
            state[method] = {"enriched": enriched, "depleted": depleted,
                             "all_kmers": all_kmers, "_calculation": context}

        return self.KEA_results

        
    def KmerCorrelationWithScore(self, txids_and_score: pd.DataFrame, species=None, method="spearman", sign_emp_pval=0.05):
       
        species_key = species if species is not None else "input"
        print(f"[KEA] KmerCorrelationWithScore for species '{species_key}' using {method}")

        # Get kmers table and transcript IDs for this species
        kmer_table, ref_group, comp_group = self.get_kmer_table_for_species(species_key)

        # Ensure txids_and_score has correct columns
        txids_and_score = txids_and_score.copy()
        if txids_and_score.shape[1] != 2:
            raise ValueError("txids_and_score must have exactly 2 columns: transcript_id, score")
        txids_and_score.columns = ["ensembl_transcript_id", "score"]

        # Keep only transcripts present in the kmer table
        common_tx = [tx for tx in txids_and_score["ensembl_transcript_id"] if tx in kmer_table.columns]
        if len(common_tx) == 0:
            raise ValueError(f"No matching transcripts found for species '{species_key}' in kmer table.")
        txids_and_score = txids_and_score[txids_and_score["ensembl_transcript_id"].isin(common_tx)]

        # Normalize kmer counts by transcript
        fract_table = kmer_table[common_tx].div(kmer_table[common_tx].sum(axis=0), axis=1).T
        fract_table = pd.merge(fract_table, txids_and_score, left_index=True, right_on="ensembl_transcript_id").set_index("ensembl_transcript_id")

        kmers = [col for col in fract_table.columns if col != "score"]
        corr_values = []

        print(f"Computing correlation for {len(kmers)} kmers...", end="\r")
        for i, kmer in enumerate(kmers):
            corr_i = fract_table[[kmer, "score"]].corr(method=method)
            corr_values.append(corr_i.loc["score", kmer])
            if (i+1) % 1000 == 0:
                print(f"{i+1}/{len(kmers)} kmers processed...", end="\r")

        kmers_corr = pd.DataFrame({
            "kmer": kmers,
            "corr": corr_values
        }).sort_values("corr", ascending=False)
        kmers_corr["rank"] = range(len(kmers_corr))
        kmers_corr["rank_percentage"] = kmers_corr["rank"] / len(kmers_corr)

        # Classify enriched/depleted kmers
        kmers_corr["kmers_groups"] = "none"
        kmers_corr.loc[kmers_corr["rank_percentage"] < sign_emp_pval, "kmers_groups"] = "enriched"
        kmers_corr.loc[kmers_corr["rank_percentage"] > (1 - sign_emp_pval), "kmers_groups"] = "depleted"
        kmers_corr["color"] = "lightgrey"
        kmers_corr.loc[kmers_corr["kmers_groups"] == "enriched", "color"] = "red"
        kmers_corr.loc[kmers_corr["kmers_groups"] == "depleted", "color"] = "blue"
        kmers_corr["method"] = method

        # Save results
        subdir_out = os.path.join(self.dir_out, f"KmerCorrelationWithScore/{species_key}/")
        os.makedirs(subdir_out,exist_ok=True)
        kmers_corr.to_csv(os.path.join(subdir_out, f"corr_analysis.{method}.txt"), sep="\t", index=False)

        enriched_kmers = kmers_corr.loc[kmers_corr["kmers_groups"] == "enriched", "kmer"].tolist()
        depleted_kmers = kmers_corr.loc[kmers_corr["kmers_groups"] == "depleted", "kmer"].tolist()

        FileByList(os.path.join(subdir_out, f"corr_analysis.{method}.enriched.txt"), enriched_kmers)
        FileByList(os.path.join(subdir_out, f"corr_analysis.{method}.depleted.txt"), depleted_kmers)

        # -----------------------------
        # Store results in KEA_results dictionary under species
        # -----------------------------
        if not hasattr(self, "KEA_results"):
            self.KEA_results = {}
        if species_key not in self.KEA_results:
            self.KEA_results[species_key] = {}
        self.KEA_results[species_key]["correlation"] = kmers_corr

        print(f"[KEA] Correlation completed. Enriched: {len(enriched_kmers)}, Depleted: {len(depleted_kmers)}")
        return kmers_corr

    
    def LogoByList(self, lista_seq, nome_logo="LogoKmers", rna=False):
    
        if len(lista_seq) == 0:
            print("LogoByList: no valid kmers after cleaning")
            return
    
        clean_kmers = lista_seq.copy()
    
        if rna:
            clean_kmers = [k.replace("T", "U") for k in clean_kmers]
    
        mode = "rna" if rna else "dna"
        SeqLogoByList(clean_kmers, mode)
    
        fig = plt.gcf()
        fig.savefig(os.path.join(self.dir_out, f"{nome_logo}.pdf"))
        plt.close()


    def CompareTranscriptKmerProfiles(
    self,
    target_tx,
    other_txs,
    target_species=None,
    other_species=None,
    metric="spearman",
    pseudocount=1e-9,
    k=None,
    fasta_input=False,
    fasta_species_key="_external_db"):
        """
        Compare transcript k-mer profiles using ALL kmers in the matrix.
        Frequencies are computed globally (per transcript).

        Parameters
        ----------
        target_tx : str
            Transcript ID to use as query (must exist in target species kmer table).
        other_txs : list of str or str (path to fasta)
            Either a list of transcript IDs already in the kmer table,
            or a path to a FASTA file if fasta_input=True.
        target_species : str, optional
            Species key for the target transcript.
        other_species : str, optional
            Species key for the other transcripts (ignored if fasta_input=True).
        metric : str
            One of: spearman, pearson, cosine, jsd.
        pseudocount : float
            Added after normalization to avoid log(0) / division issues.
        k : int, optional
            K-mer size to use when fasta_input=True.
            If None, inferred from the target kmer table index.
        fasta_input : bool
            If True, other_txs is treated as a path to a FASTA file.
            KmersCountsTable will be run on it automatically.
        fasta_species_key : str
            Key under which the fasta kmer table will be stored in
            kmers_count_table_by_species. Defaults to '_fasta_input_tmp'.
            Set a meaningful name to reuse the table across multiple calls.
        """

        target_key = target_species if target_species is not None else "input"
        tgt_table = self.kmers_count_table_by_species.get(target_key)

        if tgt_table is None:
            raise ValueError(f"Missing kmer table for target species '{target_key}'.")
        if target_tx not in tgt_table.columns:
            raise ValueError(f"{target_tx} not found in '{target_key}'.")

        tgt_table = tgt_table.drop(columns=["k"], errors="ignore")

        # ----------------------------------------------------------------
        # BRANCH A — FASTA file input
        # ----------------------------------------------------------------
        if fasta_input:
            if not os.path.isfile(other_txs):
                raise ValueError(f"fasta_input=True but file not found: {other_txs}")

            # infer k from target table if not provided
            if k is None:
                sample_kmer = tgt_table.index[0]
                k = len(sample_kmer)
                print(f"[KEA] Inferred k={k} from target kmer table.")

            # check if table already exists under this key (reuse if so)
            if fasta_species_key in self.kmers_count_table_by_species:
                print(f"[KEA] Reusing existing kmer table for key '{fasta_species_key}'.")
            else:
                print(f"[KEA] Running KmersCountsTable on FASTA: {other_txs}")
                # register fasta path so KmersCountsTable can find it
                self.combined_fasta_by_species[fasta_species_key] = other_txs
                self.KmersCountsTable(k=k, species=fasta_species_key)
                print(f"[KEA] Kmer table stored under key '{fasta_species_key}'.")

            oth_table = self.kmers_count_table_by_species[fasta_species_key].drop(columns=["k"], errors="ignore")
            tx_list = list(oth_table.columns)
            print(f"[KEA] Comparing against {len(tx_list)} sequences from FASTA.")

        # ----------------------------------------------------------------
        # BRANCH B — list of transcript IDs
        # ----------------------------------------------------------------
        else:
            other_key = other_species if other_species is not None else target_key
            oth_table = self.kmers_count_table_by_species.get(other_key)

            if oth_table is None:
                raise ValueError(f"Missing kmer table for other species '{other_key}'.")

            oth_table = oth_table.drop(columns=["k"], errors="ignore")

            missing = [tx for tx in other_txs if tx not in oth_table.columns]
            if missing:
                raise ValueError(f"Transcripts not found in '{other_key}': {missing}")

            tx_list = list(other_txs)

        # ----------------------------------------------------------------
        # NORMALIZATION
        # ----------------------------------------------------------------
        tgt_freq = tgt_table.div(tgt_table.sum(axis=0), axis=1).fillna(0) + pseudocount
        oth_freq = oth_table.div(oth_table.sum(axis=0), axis=1).fillna(0) + pseudocount

        # ----------------------------------------------------------------
        # COMMON KMER SPACE
        # ----------------------------------------------------------------
        common_kmers = tgt_freq.index.intersection(oth_freq.index)
        if len(common_kmers) == 0:
            raise ValueError("No common kmers between target table and other sequences.")

        v_target = tgt_freq.loc[common_kmers, target_tx].values

        # ----------------------------------------------------------------
        # SIMILARITY
        # ----------------------------------------------------------------
        scores = {}
        for tx in tx_list:
            v = oth_freq.loc[common_kmers, tx].values

            if v.sum() == 0 or v_target.sum() == 0:
                scores[tx] = 0.0
                continue

            try:
                if metric == "spearman":
                    score, _ = spearmanr(v_target, v)
                elif metric == "pearson":
                    score, _ = pearsonr(v_target, v)
                elif metric == "cosine":
                    score = 1 - cosine(v_target, v)
                elif metric == "jsd":
                    score = 1 - jensenshannon(v_target, v)
                else:
                    raise ValueError(f"Unknown metric: {metric}")

                if np.isnan(score):
                    score = 0.0

            except Exception:
                score = 0.0

            scores[tx] = float(score)

        return (pd.DataFrame.from_dict(scores, orient="index", columns=[f"{metric}_similarity"]).sort_values(by=f"{metric}_similarity", ascending=False))
        


    def plot_kmer_rank_alignment(
    self,
    ref_tx,
    ctrl_tx,
    ref_species="input",
    ctrl_species="input",
    bg_species="input",
    kmers=None,
    exclude_ctrl_from_bg=True,
    top_n=None,
    min_freq=None,
    bin_size=100,           # kmers per rank bin
    alpha_band=0.4,         # transparency of flow bands
    color_by="gc",          # "gc" | "ref_rank" | "category"
    kmer_categories=None,   # optional pd.Series(index=kmer, values=category)
                              # used when color_by="category"
    figsize=(10, 10),
    save=None):
        
        """
        Sankey-style kmer rank alignment plot.
    
        Kmers are binned into rank blocks of `bin_size`. Each block is drawn
        as a filled band connecting the reference rank position to the control
        and background rank positions. Band width = number of kmers in the block.
        Band color reflects GC content, reference rank, or kmer category.
    
        Parameters
        ----------
        bin_size : int
            Number of kmers per rank bin. Smaller = finer resolution.
        alpha_band : float
            Transparency of the Bezier bands.
        color_by : str
            "gc"       — mean GC content of kmers in the bin (coolwarm palette)
            "ref_rank" — position in reference ranking (viridis)
            "category" — kmer_categories Series must be provided
        kmer_categories : pd.Series, optional
            Index = kmer strings, values = category labels (str).
            Required when color_by="category".
        """
        import matplotlib.patches as mpatches
        from matplotlib.path import Path
        import matplotlib.patches as patches
    
        # ------------------------------------------------------------------
        # 1. Load + normalize tables
        # ------------------------------------------------------------------
        ref_df  = self.kmers_count_table_by_species[ref_species].drop(columns=["k"], errors="ignore")
        ctrl_df = self.kmers_count_table_by_species[ctrl_species].drop(columns=["k"], errors="ignore")
        bg_df   = self.kmers_count_table_by_species[bg_species].drop(columns=["k"], errors="ignore")
    
        if ref_tx  not in ref_df.columns:  raise ValueError(f"{ref_tx} not in {ref_species}")
        if ctrl_tx not in ctrl_df.columns: raise ValueError(f"{ctrl_tx} not in {ctrl_species}")
    
        ref_df  = ref_df.div(ref_df.sum(axis=0),   axis=1)
        ctrl_df = ctrl_df.div(ctrl_df.sum(axis=0), axis=1)
        bg_df   = bg_df.div(bg_df.sum(axis=0),     axis=1)
    
        if exclude_ctrl_from_bg and ctrl_tx in bg_df.columns:
            bg_df = bg_df.drop(columns=[ctrl_tx])
    
        # ------------------------------------------------------------------
        # 2. Common kmer space + optional filters
        # ------------------------------------------------------------------
        common_kmers = ref_df.index.intersection(ctrl_df.index).intersection(bg_df.index)
    
        kmers_use = common_kmers if kmers is None else list(set(kmers).intersection(common_kmers))
        if len(kmers_use) == 0:
            raise ValueError("No overlapping kmers between input and datasets")
    
        ref  = ref_df.loc[kmers_use, ref_tx]
        ctrl = ctrl_df.loc[kmers_use, ctrl_tx]
        bg   = bg_df.loc[kmers_use].mean(axis=1)
    
        if min_freq is not None:
            keep = ref > min_freq
            ref, ctrl, bg = ref[keep], ctrl[keep], bg[keep]
            kmers_use = ref.index
    
        if top_n is not None:
            top_kmers = ref.sort_values(ascending=False).head(top_n).index
            ref, ctrl, bg = ref.loc[top_kmers], ctrl.loc[top_kmers], bg.loc[top_kmers]
            kmers_use = top_kmers
    
        # ------------------------------------------------------------------
        # 3. Compute normalised ranks (0-1, lower = higher ranked)
        # ------------------------------------------------------------------
        N = len(ref)
        r_ref  = ref.rank( ascending=False, method="first") / N
        r_ctrl = ctrl.rank(ascending=False, method="first") / N
        r_bg   = bg.rank(  ascending=False, method="first") / N
    
        # ------------------------------------------------------------------
        # 4. Sort by reference rank and build bins
        # ------------------------------------------------------------------
        rank_df = pd.DataFrame({
            "kmer":   r_ref.index,
            "r_ref":  r_ref.values,
            "r_ctrl": r_ctrl.loc[r_ref.index].values,
            "r_bg":   r_bg.loc[r_ref.index].values,
        }).sort_values("r_ref").reset_index(drop=True)
    
        # Bin index for each kmer
        rank_df["bin"] = rank_df.index // bin_size
    
        # ------------------------------------------------------------------
        # 5. Aggregate bins
        # ------------------------------------------------------------------
        bin_df = rank_df.groupby("bin").agg(
            r_ref_mean  = ("r_ref",  "mean"),
            r_ref_min   = ("r_ref",  "min"),
            r_ref_max   = ("r_ref",  "max"),
            r_ctrl_mean = ("r_ctrl", "mean"),
            r_ctrl_min  = ("r_ctrl", "min"),
            r_ctrl_max  = ("r_ctrl", "max"),
            r_bg_mean   = ("r_bg",   "mean"),
            r_bg_min    = ("r_bg",   "min"),
            r_bg_max    = ("r_bg",   "max"),
            kmers       = ("kmer",   list),
            n           = ("kmer",   "count"),
        ).reset_index()
    
        # Width of each band proportional to bin size (in data units)
        band_half = (bin_df["r_ref_max"] - bin_df["r_ref_min"]) / 2
    
        # ------------------------------------------------------------------
        # 6. Colour mapping
        # ------------------------------------------------------------------
        def _gc(kmer):
            return (kmer.count("G") + kmer.count("C")) / len(kmer)
    
        if color_by == "gc":
            cmap = plt.cm.coolwarm
            bin_df["color_val"] = bin_df["kmers"].apply(lambda ks: np.mean([_gc(k) for k in ks]))
            norm = plt.Normalize(0, 1)
            colors = [cmap(norm(v)) for v in bin_df["color_val"]]
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar_label = "Mean GC content"
    
        elif color_by == "ref_rank":
            cmap = plt.cm.viridis_r
            norm = plt.Normalize(0, 1)
            colors = [cmap(norm(v)) for v in bin_df["r_ref_mean"]]
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar_label = "Reference rank (0=top)"
    
        elif color_by == "category":
            if kmer_categories is None:
                raise ValueError("color_by='category' requires kmer_categories Series")
            # Map each bin to its most common category
            cat_palette = {
                cat: plt.cm.tab10(i)
                for i, cat in enumerate(kmer_categories.unique())
            }
            def _dominant_cat(ks):
                cats = [kmer_categories.get(k, "Unknown") for k in ks]
                return max(set(cats), key=cats.count)
            bin_df["dominant_cat"] = bin_df["kmers"].apply(_dominant_cat)
            colors = [cat_palette.get(c, "grey") for c in bin_df["dominant_cat"]]
            sm = None
            cbar_label = None
    
        else:
            raise ValueError(f"Unknown color_by='{color_by}'")
    
        # ------------------------------------------------------------------
        # 7. Draw Sankey bands using cubic Bezier paths
        # ------------------------------------------------------------------
        def _bezier_band(ax, x_left, x_right,
                         y_left_center, y_left_half,
                         y_right_center, y_right_half,
                         color, alpha):
            """
            Draw a filled Bezier band between two columns.
            The band has vertical extent [center-half, center+half] on each side.
            """
            top_left   = y_left_center  - y_left_half
            bot_left   = y_left_center  + y_left_half
            top_right  = y_right_center - y_right_half
            bot_right  = y_right_center + y_right_half
    
            # Control points: pull horizontally toward the midpoint
            cx = (x_left + x_right) / 2
    
            verts = [
                # top edge: left -> right (Bezier)
                (x_left,  top_left),
                (cx,      top_left),
                (cx,      top_right),
                (x_right, top_right),
                # bottom edge: right -> left (Bezier, reversed)
                (x_right, bot_right),
                (cx,      bot_right),
                (cx,      bot_left),
                (x_left,  bot_left),
                (x_left,  top_left),  # close
            ]
            codes = [
                Path.MOVETO,
                Path.CURVE4, Path.CURVE4, Path.CURVE4,   # top Bezier (needs 3 pts after MOVETO)
                Path.LINETO,
                Path.CURVE4, Path.CURVE4, Path.CURVE4,   # bottom Bezier reversed
                Path.CLOSEPOLY,
            ]
            path = Path(verts, codes)
            patch = patches.PathPatch(path, facecolor=color, edgecolor="none", alpha=alpha)
            ax.add_patch(patch)
    
        fig, ax = plt.subplots(figsize=figsize)
    
        x_ctrl, x_ref, x_bg = 0, 1, 2   # column positions
    
        for i, row in bin_df.iterrows():
            c     = colors[i]
            half  = max(row["n"] / (2 * N), 0.002)   # minimum visible band width
    
            # ref -> ctrl band
            _bezier_band(
                ax,
                x_left=x_ctrl,  x_right=x_ref,
                y_left_center=row["r_ctrl_mean"],  y_left_half=half,
                y_right_center=row["r_ref_mean"],  y_right_half=half,
                color=c, alpha=alpha_band
            )
            # ref -> bg band
            _bezier_band(
                ax,
                x_left=x_ref,  x_right=x_bg,
                y_left_center=row["r_ref_mean"],  y_left_half=half,
                y_right_center=row["r_bg_mean"],  y_right_half=half,
                color=c, alpha=alpha_band
            )
    
        # Column spine lines
        for x in [x_ctrl, x_ref, x_bg]:
            ax.axvline(x, color="black", linewidth=1.5, zorder=5)
    
        # ------------------------------------------------------------------
        # 8. Similarity metrics as subtitle
        # ------------------------------------------------------------------
        delta_ctrl = np.abs(r_ref - r_ctrl).mean()
        delta_bg   = np.abs(r_ref - r_bg).mean()
        sim_ctrl   = 1 - delta_ctrl
        sim_bg     = 1 - delta_bg
    
        ax.set_xticks([x_ctrl, x_ref, x_bg])
        ax.set_xticklabels(["Control", "Query", "Background"], fontsize=13)
        ax.set_xlim(-0.4, 2.4)
        ax.set_ylim(0, 1)
        ax.invert_yaxis()
        ax.set_ylabel("Normalised rank  (0 = most frequent)", fontsize=11)
        ax.set_title(
            f"K-mer Rank Alignment  —  {ref_tx}\n"
            f"n={N} kmers  |  bins of {bin_size}  |  "
            f"sim(ctrl)={sim_ctrl:.3f}   sim(bg)={sim_bg:.3f}",
            fontsize=11
        )
        ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
        ax.yaxis.set_ticks_position("left")
    
        # Colorbar or legend
        if sm is not None:
            cbar = fig.colorbar(sm, ax=ax, shrink=0.4, pad=0.02)
            cbar.set_label(cbar_label, fontsize=9)
        elif color_by == "category":
            legend_handles = [
                mpatches.Patch(color=cat_palette[c], label=c)
                for c in kmer_categories.unique()
            ]
            ax.legend(handles=legend_handles, fontsize=8,
                      loc="lower right", title="Category")
    
        plt.tight_layout()
    
        if save:
            fig.savefig(save, bbox_inches="tight", dpi=150)
            print(f"[KEA] Plot saved → {save}")
    
        plt.show()
        plt.close()
    
        print(f"Similarity (query vs control):    {sim_ctrl:.4f}")
        print(f"Similarity (query vs background): {sim_bg:.4f}")

    
 
    
    def FetchOrthologsAndSaveFastaOneToOne(self, ref_species, target_species, out_prefix="homologs", biomart_host="http://www.ensembl.org"):

        # --------------------------------------------------
        # 1. Input transcript IDs
        # --------------------------------------------------
        ref_tx = get_fasta_ids(self.reference_fasta)
        ctrl_tx = get_fasta_ids(self.control_fasta)
        all_tx = list(set(ref_tx + ctrl_tx))

        # --------------------------------------------------
        # 2. BioMart datasets
        # --------------------------------------------------
        ref_ds = Dataset(name=f"{ref_species}_gene_ensembl", host=biomart_host)
        tgt_ds = Dataset(name=f"{target_species}_gene_ensembl", host=biomart_host)

        # --------------------------------------------------
        # 3. Query orthology (GENE LEVEL)
        # --------------------------------------------------
        attrs = [
            "ensembl_transcript_id",
            "ensembl_gene_id",
            f"{target_species}_homolog_ensembl_gene",
            f"{target_species}_homolog_orthology_type",
            f"{target_species}_homolog_orthology_confidence",
            f"{target_species}_homolog_perc_id"
        ]
        df = ref_ds.query(attributes=attrs)
        df = df[df["Transcript stable ID"].isin(all_tx)]
        df.columns = ["ref_tx", "ref_gene", "target_gene", "orthology_type", "orthology_confidence", "perc_id"]

        # label ref / ctrl
        df["label"] = "ctrl"
        df.loc[df["ref_tx"].isin(ref_tx), "label"] = "ref"

        # --------------------------------------------------
        # 4. Pick BEST TARGET GENE per reference transcript
        # --------------------------------------------------
        df_best = (
            df.sort_values(["ref_tx", "orthology_confidence", "perc_id"], ascending=[True, False, False])
            .drop_duplicates(subset="ref_tx", keep="first")
        )

        # Save mapping (optional)
        out_dir = os.path.join(self.dir_out, target_species)
        os.makebio.Dirs(out_dir, exist_ok=True)
        mapping_file = os.path.join(out_dir, f"{out_prefix}_mapping_{target_species}.csv")
        df_best.to_csv(mapping_file, index=False)
        print(f"[KEA] Mapping saved to {mapping_file}")

        # --------------------------------------------------
        # 5. Fetch TARGET transcripts (choose longest CDS)
        # --------------------------------------------------
        target_genes = df_best["target_gene"].dropna().unique()
        tgt_df = tgt_ds.query(attributes=[
            "ensembl_transcript_id",
            "ensembl_gene_id",
            "transcript_length",
            "transcript_biotype",
            "cdna"
        ])
        tgt_df = tgt_df[
            (tgt_df["Gene stable ID"].isin(target_genes)) &
            (tgt_df["Transcript type"] == "protein_coding")
        ]

        # Longest transcript per TARGET gene
        tgt_longest = tgt_df.loc[
            tgt_df.groupby("Gene stable ID")["Transcript length (including UTRs and CDS)"].idxmax()
        ]

        gene_to_seq = {
            row["Gene stable ID"]: (row["Transcript stable ID"], row["cDNA sequences"])
            for _, row in tgt_longest.iterrows() if row["cDNA sequences"]
        }

        # --------------------------------------------------
        # 6. Write separate FASTA files for ref / ctrl
        # --------------------------------------------------
        if not hasattr(self, "combined_fasta_by_species"):
            self.combined_fasta_by_species = {}

        for label in ["ref", "ctrl"]:
            records = []
            df_subset = df_best[df_best["label"] == label]

            for _, row in df_subset.iterrows():
                tgt_gene = row["target_gene"]
                if tgt_gene not in gene_to_seq:
                    continue
                tgt_tx, seq = gene_to_seq[tgt_gene]
                records.append(SeqRecord(Seq(seq), id=tgt_tx, description=""))

            out_fasta = os.path.join(out_dir, f"{label}_{target_species}.fa")
            with open(out_fasta, "w") as f:
                for rec in records:
                    f.write(f">{rec.id}\n{str(rec.seq)}\n")
            print(f"[KEA] Saved {len(records)} sequences to {out_fasta}")


            self.combined_fasta_by_species[f"{target_species}_{label}"] = out_fasta

    def IntersectEnrichedKmersAndPlotVenn(self, save_plot=None):

        if not hasattr(self, "KEA_results") or not self.KEA_results:
            raise ValueError("KEA_results is empty. Run ExtractKmers first.")

        species_method_intersections = {}

        # Step 1: Intersection across methods for each species
        for species, methods_dict in self.KEA_results.items():
            method_sets = []
            for method, result in methods_dict.items():
                enriched = set(result.get("enriched", []))
                if enriched:
                    method_sets.append(enriched)
            if method_sets:
                species_method_intersections[species] = set.intersection(*method_sets)
            else:
                species_method_intersections[species] = set()

        # Step 2: Intersection across species
        if species_method_intersections:
            final_intersection = set.intersection(*species_method_intersections.values())
        else:
            final_intersection = set()

        # Step 3: Venn diagram for 2–3 species
        species_list = list(species_method_intersections.keys())
        n_species = len(species_list)

        if n_species == 2:
            s1, s2 = species_list
            venn2(
                [species_method_intersections[s1], species_method_intersections[s2]],
                set_labels=[s1, s2]
            )
            plt.title("Venn Diagram of Enriched Kmers")
        elif n_species == 3:
            s1, s2, s3 = species_list
            venn3(
                [species_method_intersections[s1], species_method_intersections[s2], species_method_intersections[s3]],
                set_labels=[s1, s2, s3]
            )
            plt.title("Venn Diagram of Enriched Kmers")
        else:
            print(f"Venn diagram not supported for {n_species} species. Only 2 or 3 species can be plotted.")

        if n_species in [2,3]:
            if save_plot:
                plt.savefig(save_plot, bbox_inches="tight")
                plt.show()
                plt.close()
                print(f"Venn diagram saved to {save_plot}")
            else:
                plt.show()
                plt.close()

        return species_method_intersections, final_intersection

    def KmersProteinEnrichment(
    self,
    kmers,
    df_zscore,
    # --- kmer preprocessing ---
    k_out=None,
    deconvolute_unique=True,
    convert_to_rna=True,
    # --- analysis mode ---
    mode="continuous",
    test="mwu",
    # --- binary mode ---
    threshold=None,
    top_n=None,
    # --- multiple testing ---
    fdr_alpha=0.05,
    # --- significance / plotting ---
    plot=True,
    p_col="fdr",
    p_thresh=0.05,
    fc_thresh=1.0,
    highlight=None,
    highlight_top_n=10,
    # --- output ---
    save=None,
    # --- background ---
    bg_kmers=None,
    # --- plotting ---
    plot_type="lollipop",
    wordcloud_score="cliff_delta",
    wordcloud_top_n=10,
    wordcloud_font_path=None):

            # ----------------------------------------------------------------
            # NEW DEFAULT OUTPUT: original k-mer × protein/profile PEKA scores
            #
            # Retrieve the additional DataFrame with:
            #     self.kmer_protein_scores
            # e.g.:
            #     synaptosome.kmer_protein_scores
            #
            # Each 7-mer contributes its three overlapping pentamers.
            # Repeated pentamers are counted at each position.
            # A missing/nonfinite score in any pentamer gives NaN for that cell.
            # Protein/profile names remain exactly those supplied in df_zscore.
            # The original return value and enrichment analysis remain unchanged.
            # ----------------------------------------------------------------
            original_kmers = list(kmers)
    
            score_lookup = df_zscore.copy()
            score_lookup.index = pd.Index([
                str(k).strip().upper().replace("T", "U")
                for k in score_lookup.index
            ])
            score_lookup = score_lookup.apply(
                pd.to_numeric, errors="coerce"
            )
            score_lookup = score_lookup.replace(
                [np.inf, -np.inf], np.nan
            )
    
            if not score_lookup.index.is_unique:
                score_lookup = score_lookup.groupby(
                    level=0, sort=False
                ).mean()
    
            row_labels = []
            score_rows = []
            seen_labels = set()
    
            for original_kmer in original_kmers:
                sequence = str(original_kmer).strip().upper()
                label = (
                    sequence.replace("T", "U")
                    if convert_to_rna
                    else sequence
                )
    
                if label in seen_labels:
                    continue
    
                seen_labels.add(label)
                row_labels.append(label)
                sequence = sequence.replace("T", "U")
    
                pentamers = [
                    sequence[i:i + 5]
                    for i in range(len(sequence) - 4)
                ]
    
                if pentamers:
                    values = score_lookup.reindex(
                        pentamers
                    ).to_numpy(dtype=float)
                    score_rows.append(np.mean(values, axis=0))
                else:
                    score_rows.append(
                        np.full(len(df_zscore.columns), np.nan)
                    )
    
            self.kmer_protein_scores = pd.DataFrame(
                score_rows,
                index=pd.Index(row_labels, name="kmer"),
                columns=df_zscore.columns.copy(),
                dtype=float
            )
    
            # ----------------------------------------------------------------
            # 1. DECONVOLUTE KMERS (optional)
            # ----------------------------------------------------------------
            if k_out is not None:
                k_in = len(kmers[0])
                if k_out > k_in:
                    raise ValueError(
                        f"k_out ({k_out}) must be <= k_in ({k_in})"
                    )
    
                print(
                    f"[KEA] Deconvoluting {len(kmers)} kmers: "
                    f"k={k_in} -> k={k_out}"
                )
    
                def _deconv(lst):
                    out = set() if deconvolute_unique else []
                    for kmer in lst:
                        if len(kmer) != k_in:
                            raise ValueError(
                                f"Kmer '{kmer}' has length {len(kmer)}, "
                                f"expected {k_in}"
                            )
                        for i in range(k_in - k_out + 1):
                            sub = kmer[i:i + k_out]
                            if deconvolute_unique:
                                out.add(sub)
                            else:
                                out.append(sub)
                    return list(out)
    
                kmers = _deconv(kmers)
    
                if bg_kmers is not None:
                    print(
                        f"[KEA] Deconvoluting bg_kmers: "
                        f"{len(bg_kmers)} kmers"
                    )
                    bg_kmers = _deconv(bg_kmers)
                    print(
                        f"[KEA] After deconvolution: {len(kmers)} kmers | "
                        f"bg={len(bg_kmers) if bg_kmers is not None else 'default'}"
                    )
    
            # ----------------------------------------------------------------
            # 2. CONVERT T -> U (optional)
            # ----------------------------------------------------------------
            if convert_to_rna:
                df_zscore = df_zscore.copy()
                df_zscore.index = df_zscore.index.str.replace("T", "U")
                kmers = [k.replace("T", "U") for k in kmers]
                if bg_kmers is not None:
                    bg_kmers = [k.replace("T", "U") for k in bg_kmers]
    
            # ----------------------------------------------------------------
            # 3. VALIDATE KMERS AGAINST DF INDEX
            # ----------------------------------------------------------------
            index_kmers = set(df_zscore.index)
            kmers_valid = list(set(kmers).intersection(index_kmers))
            n_missing = len(set(kmers)) - len(kmers_valid)
    
            print(
                f"[KEA] Kmers of interest : {len(set(kmers))} total | "
                f"{len(kmers_valid)} found in df | {n_missing} missing"
            )
    
            if len(kmers_valid) == 0:
                raise ValueError(
                    "No valid kmers found in df_zscore index."
                )
    
            # ----------------------------------------------------------------
            # VALIDATE bg_kmers
            # ----------------------------------------------------------------
            if bg_kmers is not None:
                bg_kmers_valid = list(
                    set(bg_kmers).intersection(index_kmers)
                )
                print(
                    f"[KEA] Custom bg_kmers: {len(bg_kmers)} total | "
                    f"{len(bg_kmers_valid)} found"
                )
                if len(bg_kmers_valid) == 0:
                    raise ValueError(
                        "None of the provided bg_kmers are present "
                        "in df_zscore index."
                    )
            else:
                bg_kmers_valid = None
    
            # ----------------------------------------------------------------
            # 4. RUN ENRICHMENT
            # ----------------------------------------------------------------
            if mode == "continuous":
    
                if bg_kmers_valid is None:
                    background_kmers = list(
                        index_kmers - set(kmers_valid)
                    )
                else:
                    background_kmers = bg_kmers_valid
    
                MAX_BG_RATIO = 10
                if (
                    len(background_kmers)
                    > MAX_BG_RATIO * len(kmers_valid)
                ):
                    rng = np.random.default_rng(42)
                    background_kmers = list(
                        rng.choice(
                            background_kmers,
                            size=MAX_BG_RATIO * len(kmers_valid),
                            replace=False
                        )
                    )
    
                def cliffs_delta(x, y):
                    import numpy as np
    
                    x = np.asarray(x)
                    y = np.asarray(y)
                    n = len(x)
                    m = len(y)
                    greater = sum((xi > y).sum() for xi in x)
                    less = sum((xi < y).sum() for xi in x)
                    return (greater - less) / (n * m)
    
                results = []
    
                for col in df_zscore.columns:
                    g1 = (
                        df_zscore.loc[kmers_valid, col]
                        .dropna().values
                    )
                    g2 = (
                        df_zscore.loc[background_kmers, col]
                        .dropna().values
                    )
    
                    if len(g1) < 2 or len(g2) < 2:
                        continue
    
                    mean1, mean2 = np.mean(g1), np.mean(g2)
                    med1, med2 = np.median(g1), np.median(g2)
                    sd1, sd2 = np.std(g1), np.std(g2)
    
                    if test == "mwu":
                        stat, pval = mannwhitneyu(
                            g1, g2, alternative="greater"
                        )
                        effect = (
                            1 - (2 * stat) / (len(g1) * len(g2))
                        )
                    elif test == "ks":
                        stat, pval = ks_2samp(
                            g1, g2, alternative="greater"
                        )
                        effect = stat
                    else:
                        raise ValueError(
                            f"Unknown test '{test}'."
                        )
    
                    zscore_diff = cliffs_delta(g1, g2)
    
                    results.append({
                        "protein": col,
                        "mean_group1": mean1,
                        "mean_group2": mean2,
                        "median_group1": med1,
                        "median_group2": med2,
                        "sd1": sd1,
                        "sd2": sd2,
                        "effect_size": effect,
                        "pvalue": pval,
                        "zscore_diff": zscore_diff,
                        "n1": len(g1),
                        "n2": len(g2)
                    })
    
            # ----------------------------------------------------------------
            # 5. MULTIPLE TESTING CORRECTION
            # ----------------------------------------------------------------
            results = pd.DataFrame(results)
    
            if len(results) == 0:
                print("[KEA] Warning: no results generated.")
                return results
    
            results["fdr"] = multipletests(
                results["pvalue"], method="fdr_bh"
            )[1]
            results["bonf"] = multipletests(
                results["pvalue"], method="bonferroni"
            )[1]
            results["significant_fdr"] = (
                results["fdr"] < fdr_alpha
            )
            results["significant_bonf"] = (
                results["bonf"] < fdr_alpha
            )
            results = (
                results.sort_values("fdr")
                .reset_index(drop=True)
            )
    
            n_sig_fdr = results["significant_fdr"].sum()
            n_sig_bonf = results["significant_bonf"].sum()
            print(
                f"[KEA] Significant proteins — "
                f"FDR<{fdr_alpha}: {n_sig_fdr} | "
                f"Bonferroni<{fdr_alpha}: {n_sig_bonf}"
            )
    
            # ----------------------------------------------------------------
            # 6. VISUALIZATION
            # ----------------------------------------------------------------
            if plot:
                df_plot = results.copy()
                df_plot["neglog10p"] = -np.log10(
                    df_plot["pvalue"] + 1e-300
                )
                sig_mask = (
                    (df_plot[p_col] < p_thresh)
                    & (np.abs(df_plot["zscore_diff"]) > fc_thresh)
                )
    
                hl_mask = pd.Series(False, index=df_plot.index)
                if highlight is not None:
                    if isinstance(highlight, str):
                        highlight = [highlight]
                    hl_mask = df_plot["protein"].str.contains(
                        "|".join(highlight),
                        case=False,
                        na=False
                    )
    
                hl_top_mask = pd.Series(
                    False, index=df_plot.index
                )
                if highlight_top_n is not None:
                    sig = df_plot[sig_mask]
                    top_idx = sig.nlargest(
                        highlight_top_n, "zscore_diff"
                    ).index
                    hl_top_mask = df_plot.index.isin(top_idx)
    
                # ------------------------------------------------------------
                # VOLCANO PLOT
                # ------------------------------------------------------------
                if plot_type == "volcano":
                    fig, ax = plt.subplots(figsize=(9, 7))
                    ax.scatter(
                        df_plot.loc[~sig_mask, "zscore_diff"],
                        df_plot.loc[~sig_mask, "neglog10p"],
                        color="lightgrey",
                        alpha=0.5,
                        s=25,
                        label="n.s."
                    )
                    ax.scatter(
                        df_plot.loc[
                            sig_mask & ~hl_mask, "zscore_diff"
                        ],
                        df_plot.loc[
                            sig_mask & ~hl_mask, "neglog10p"
                        ],
                        color="tomato",
                        alpha=0.8,
                        s=35,
                        label="significant"
                    )
                    ax.scatter(
                        df_plot.loc[hl_mask, "zscore_diff"],
                        df_plot.loc[hl_mask, "neglog10p"],
                        color="royalblue",
                        edgecolor="black",
                        s=110,
                        label="highlight"
                    )
                    for _, row in df_plot.loc[
                        hl_mask | hl_top_mask
                    ].iterrows():
                        ax.text(
                            row["zscore_diff"],
                            row["neglog10p"],
                            row["protein"],
                            fontsize=8,
                            ha="left",
                            va="bottom"
                        )
    
                    ax.axhline(
                        -np.log10(p_thresh),
                        linestyle="--",
                        color="grey",
                        linewidth=0.8
                    )
                    ax.axvline(
                        fc_thresh,
                        linestyle="--",
                        color="grey",
                        linewidth=0.8
                    )
                    ax.axvline(
                        -fc_thresh,
                        linestyle="--",
                        color="grey",
                        linewidth=0.8
                    )
                    ax.set_xlabel(
                        "Cliff’s delta (group1 vs background)"
                    )
                    ax.set_ylabel("-log10(pvalue)")
                    ax.set_title("Volcano plot (Cliff’s delta)")
                    ax.legend(fontsize=8)
    
                    if save is not None:
                        filenames = (
                            save
                            if isinstance(save, (list, tuple))
                            else [save]
                        )
                        for filename in filenames:
                            fig.savefig(
                                filename,
                                dpi=600,
                                bbox_inches="tight",
                                facecolor="white"
                            )
    
                    plt.tight_layout()
                    plt.show()
                    plt.close()
    
                # ------------------------------------------------------------
                # LOLLIPOP PLOT
                # ------------------------------------------------------------
                elif plot_type == "lollipop":
                    sig = df_plot[sig_mask].sort_values(
                        "zscore_diff"
                    )
                    fig, ax = plt.subplots(
                        figsize=(9, max(4, len(sig) * 0.35))
                    )
                    ax.hlines(
                        sig["protein"],
                        0,
                        sig["zscore_diff"],
                        color="grey",
                        alpha=0.6
                    )
                    ax.scatter(
                        sig["zscore_diff"],
                        sig["protein"],
                        c="tomato",
                        s=40
                    )
                    sig_hl_mask = sig["protein"].isin(
                        df_plot.loc[hl_mask, "protein"]
                    )
                    ax.scatter(
                        sig.loc[sig_hl_mask, "zscore_diff"],
                        sig.loc[sig_hl_mask, "protein"],
                        c="royalblue",
                        edgecolor="black",
                        s=110
                    )
                    ax.set_xlabel("Cliff’s delta")
                    ax.set_ylabel("Protein")
                    ax.set_title(
                        "Lollipop plot (significant proteins)"
                    )
    
                    if save is not None:
                        filenames = (
                            save
                            if isinstance(save, (list, tuple))
                            else [save]
                        )
                        for filename in filenames:
                            fig.savefig(
                                filename,
                                dpi=600,
                                bbox_inches="tight",
                                facecolor="white"
                            )
    
                    plt.tight_layout()
                    plt.show()
    
                # ------------------------------------------------------------
                # BARPLOT
                # ------------------------------------------------------------
                elif plot_type == "barplot":
                    df_bar = df_plot.copy().sort_values(
                        "zscore_diff"
                    )
                    colors = [
                        "tomato"
                        if (
                            row[p_col] < p_thresh
                            and abs(row["zscore_diff"]) > fc_thresh
                        )
                        else "lightgrey"
                        for _, row in df_bar.iterrows()
                    ]
                    fig, ax = plt.subplots(
                        figsize=(max(9, len(df_bar) * 0.25), 6)
                    )
                    ax.bar(
                        df_bar["protein"],
                        df_bar["zscore_diff"],
                        color=colors
                    )
                    ax.set_xlim(-0.5, len(df_bar) - 0.5)
                    ax.axhline(
                        0, color="black", linewidth=0.8
                    )
                    ax.set_ylabel("Cliff’s delta")
                    ax.set_xlabel("Protein")
                    ax.set_title(
                        "Cliff’s delta — All proteins (barplot)"
                    )
                    plt.xticks(rotation=90)
                    plt.tight_layout()
    
                    if save is not None:
                        filenames = (
                            save
                            if isinstance(save, (list, tuple))
                            else [save]
                        )
                        for filename in filenames:
                            fig.savefig(
                                filename,
                                dpi=600,
                                bbox_inches="tight",
                                facecolor="white"
                            )
    
                    plt.show()
    
                # ------------------------------------------------------------
                # WORDCLOUD
                # ------------------------------------------------------------
                elif plot_type == "wordcloud":
                    df_wc = df_plot.loc[
                        (df_plot[p_col] < p_thresh)
                        & (df_plot["zscore_diff"] >= fc_thresh)
                    ].copy()
                    df_wc = df_wc.dropna(
                        subset=["protein", "zscore_diff"]
                    )
                    df_wc = df_wc.loc[
                        np.isfinite(df_wc["zscore_diff"])
                    ].copy()
                    df_wc["protein"] = (
                        df_wc["protein"].astype(str).str.strip()
                    )
                    df_wc = df_wc.loc[
                        df_wc["protein"].ne("")
                    ]
                    df_wc = df_wc.drop_duplicates(
                        subset="protein", keep="first"
                    )
    
                    if wordcloud_top_n is not None:
                        if (
                            not isinstance(
                                wordcloud_top_n, (int, np.integer)
                            )
                            or wordcloud_top_n < 1
                        ):
                            raise ValueError(
                                "wordcloud_top_n must be a "
                                "positive integer or None."
                            )
                        df_wc = df_wc.head(wordcloud_top_n)
    
                    if df_wc.empty:
                        raise ValueError(
                            "No proteins passed the selected thresholds."
                        )
    
                    if wordcloud_score == "cliff_delta":
                        df_wc["wordcloud_score"] = (
                            df_wc["zscore_diff"]
                        )
                        score_label = "Cliff’s delta"
                        title_score = "Cliff’s delta"
                    elif wordcloud_score == "median_ratio":
                        denominator = df_wc["median_group2"] + 1
                        df_wc["wordcloud_score"] = (
                            df_wc["median_group1"]
                            / denominator.mask(
                                np.isclose(denominator, 0)
                            )
                        )
                        score_label = (
                            "Median group 1 / (median group 2 + 1)"
                        )
                        title_score = "median-ratio score"
                    else:
                        raise ValueError(
                            "wordcloud_score must be "
                            "'cliff_delta' or 'median_ratio'."
                        )
    
                    invalid = (
                        ~np.isfinite(df_wc["wordcloud_score"])
                        | (df_wc["wordcloud_score"] <= 0)
                    )
                    if invalid.any():
                        names = df_wc.loc[
                            invalid, "protein"
                        ].tolist()
                        raise ValueError(
                            "Selected proteins have undefined or "
                            f"non-positive plotting scores: {names}"
                        )
    
                    print(
                        f"\n[KEA] WordCloud selection: "
                        f"top {len(df_wc)} by Cliff's delta"
                    )
                    print(
                        df_wc[
                            ["protein", "zscore_diff", "wordcloud_score"]
                        ].to_string(index=False)
                    )
    
                    score_by_protein = dict(zip(
                        df_wc["protein"],
                        df_wc["wordcloud_score"]
                    ))
                    word_dict = score_by_protein.copy()
                    score_min = df_wc["wordcloud_score"].min()
                    score_max = df_wc["wordcloud_score"].max()
    
                    if score_min == score_max:
                        padding = max(
                            abs(score_min) * 0.01, 0.01
                        )
                        norm = mpl.colors.Normalize(
                            vmin=score_min - padding,
                            vmax=score_max + padding
                        )
                    else:
                        norm = mpl.colors.Normalize(
                            vmin=score_min, vmax=score_max
                        )
    
                    cmap = (
                        mpl.colors.LinearSegmentedColormap.from_list(
                            "rbp_score_greys",
                            ["#BDBDBD", "#737373", "#000000"]
                        )
                    )
    
                    def score_color_func(word, **kwargs):
                        return mpl.colors.to_hex(
                            cmap(norm(score_by_protein[word]))
                        )
    
                    if wordcloud_font_path is None:
                        available_fonts = {
                            font.name
                            for font in fm.fontManager.ttflist
                        }
                        font_name = next(
                            (
                                name
                                for name in [
                                    "Helvetica",
                                    "Liberation Sans",
                                    "Arial",
                                    "DejaVu Sans"
                                ]
                                if name in available_fonts
                            ),
                            "DejaVu Sans"
                        )
                        font_path = fm.findfont(
                            fm.FontProperties(family=font_name)
                        )
                    else:
                        font_path = str(wordcloud_font_path)
                        fm.fontManager.addfont(font_path)
                        font_name = fm.FontProperties(
                            fname=font_path
                        ).get_name()
    
                    canvas_width, canvas_height = 1600, 700
                    wc = WordCloud(
                        width=canvas_width,
                        height=canvas_height,
                        background_color="white",
                        color_func=score_color_func,
                        font_path=font_path,
                        min_font_size=12,
                        max_font_size=None,
                        relative_scaling="auto",
                        prefer_horizontal=1.0,
                        random_state=42,
                        collocations=False,
                        max_words=len(word_dict),
                        repeat=False,
                        margin=6
                    ).generate_from_frequencies(word_dict)
    
                    if len(wc.layout_) != len(word_dict):
                        raise ValueError(
                            "Not all proteins fit. Increase canvas_width "
                            "or canvas_height, or reduce min_font_size."
                        )
    
                    with plt.rc_context({
                        "font.family": font_name,
                        "font.size": 12,
                        "axes.titleweight": "bold",
                        "axes.labelweight": "bold"
                    }):
                        fig, ax = plt.subplots(
                            figsize=(15, 6), dpi=150
                        )
                        ax.imshow(
                            wc.to_array(), interpolation="bilinear"
                        )
                        ax.axis("off")
                        ax.set_title(
                            f"RBP enrichment by {title_score}",
                            fontsize=16,
                            fontweight="bold",
                            pad=15
                        )
    
                        scalar_mappable = mpl.cm.ScalarMappable(
                            norm=norm, cmap=cmap
                        )
                        scalar_mappable.set_array([])
                        cbar = fig.colorbar(
                            scalar_mappable,
                            ax=ax,
                            fraction=0.025,
                            pad=0.025,
                            shrink=0.80
                        )
                        cbar.set_label(
                            score_label,
                            fontsize=12,
                            fontweight="bold",
                            labelpad=12
                        )
    
                        if score_min == score_max:
                            cbar.set_ticks([score_min])
    
                        for label in cbar.ax.get_yticklabels():
                            label.set_fontsize(10)
                            label.set_fontweight("bold")
    
                        cbar.outline.set_linewidth(0.8)
                        fig.tight_layout()
    
                        if save is not None:
                            filenames = (
                                save
                                if isinstance(save, (list, tuple))
                                else [save]
                            )
                            for filename in filenames:
                                fig.savefig(
                                    filename,
                                    dpi=600,
                                    bbox_inches="tight",
                                    facecolor="white"
                                )
    
                        plt.show()
                        plt.close(fig)
    
            return results

    def ImportMetaTx(self,file_name,header="yes"):
        if header == "yes":
            self.metatx=pd.read_table(file_name)
        else:
            self.metatx=pd.read_table(file_name,header=None)
            self.metatx.columns=["chrom","start","end","name","score","strand"]
    
    def ImportMetaTx_Mapped(self,file_name,header="yes"):
        
        if header == "yes":
            self.metatx_mapped=pd.read_table(file_name)
            
        else:
            self.metatx_mapped=pd.read_table(file_name,header=None)
            self.metatx_mapped.columns=["chrom","start","end","name","score","strand"]
        
    def ImportMetaTx_MappedScore(self,file_name,header="yes"):
        
        if header == "yes":
            self.metatx_mapped_score=pd.read_table(file_name)
        else:
            self.metatx_mapped_score=pd.read_table(file_name,header=None)
            self.metatx_mapped_score.columns=["chrom","start","end","name","score","strand"]
        
    def ImportMetaTx_Enrichment(self,file_name,header="yes"):
        
        if header=="yes":
            self.enrichment=pd.read_table(file_name)
        else:
            self.enrichment=pd.read_table(file_name,header=None)

    
    def GenerateMetaTx(self,all_tx,annfile,region_sizes="by_tx"):
            
        all_tx_interacting=get_fasta_ids(all_tx)
        ann_df=pd.read_table(annfile)

        if region_sizes == "by_tx":
            ### calcolo il numero di bin per regione in base ai tx
            region_sizes=RegionSizeByTxList(all_tx_interacting,annfile)
            region_sizes=[round(x)*3 for x in region_sizes]

        else:
            region_sizes=region_sizes

        ##### TRANSCRIPTS BINNER

        tx_ann_df=ann_df.loc[ann_df["ensembl_transcript_id"].isin(all_tx_interacting)]
        tx_ann_df_pc=tx_ann_df.loc[tx_ann_df["transcript_biotype"]=="protein_coding"]
        all_tx=list(set(tx_ann_df_pc.loc[:,"ensembl_transcript_id"]))
        regions=["5UTR","CDS","3UTR"]

        count_i=0
        for tx_i in all_tx:
            count_i+=1

            print(count_i/len(all_tx),end="\r")
            for region_name in regions:

                ### 5UTR ###
                if region_name == "5UTR":
                    #print(tx_i)
                    add_start=0# nulla da aggiungere allo start
                    size_region=int(tx_ann_df_pc.loc[tx_ann_df_pc["ensembl_transcript_id"]==tx_i,"5utr"].item())
                    n_bins=region_sizes[0]

                    if size_region >= 0:
                        # se è l'inizio genero il df finale
                        if [tx_i,region_name] == [all_tx[0],regions[0]]:

                            binned_tx_bed=RegionBinner(tx_i,region_name,size_region,n_bins,add_start)
                        else:
                            # altrimenti impilo a lui
                            binned_i=RegionBinner(tx_i,region_name,size_region,n_bins,add_start)
                            binned_tx_bed = pd.concat([binned_tx_bed, binned_i], ignore_index=True)


                ### CDS ###
                if region_name == "CDS":

                    add_start=int(tx_ann_df_pc.loc[tx_ann_df_pc["ensembl_transcript_id"]==tx_i,"UTR5"].item())
                    size_region=int(tx_ann_df_pc.loc[tx_ann_df_pc["ensembl_transcript_id"]==tx_i,"cds"].item())
                    n_bins=region_sizes[1]

                    if size_region >= 0:

                        binned_i=RegionBinner(tx_i,region_name,size_region,n_bins,add_start)
                        binned_tx_bed = pd.concat([binned_tx_bed, binned_i], ignore_index=True)


                ### 3UTR ###
                if region_name == "3UTR":

                    add_start=int(tx_ann_df_pc.loc[tx_ann_df_pc["ensembl_transcript_id"]==tx_i,"CDS"].item())
                    size_region=int(tx_ann_df_pc.loc[tx_ann_df_pc["ensembl_transcript_id"]==tx_i,"3utr"].item())
                    n_bins=region_sizes[2]

                    if size_region >= 0:

                        binned_i=RegionBinner(tx_i,region_name,size_region,n_bins,add_start)
                        binned_tx_bed = pd.concat([binned_tx_bed, binned_i], ignore_index=True)


        print("region sizes - i bin sono basati su questi numeri: "+str(region_sizes)+"\n")
        print("n protein coding transcripts: "+str(len(all_tx))+"\n")
        print("n protein coding in binned bed: "+str(len(set(binned_tx_bed.loc[:,"chr"])))+"\n")


        binned_tx_bed.loc[:,"score"]=1
        binned_tx_bed.loc[:,"strand"]="+"

        binned_tx_bed.columns=["chrom","start","end","name","score","strand"]
        binned_tx_bed[['start','end']] = binned_tx_bed[['start','end']].astype(int)
        self.metatx=binned_tx_bed
        
        return binned_tx_bed

    
    def MapperTxToGenome(self,python_bin,map_to_genome,dir_temp="./mapper_bed/"):
    
        os.makedirs(dir_temp, exist_ok=True)                
        self.metatx_mapped=MapperTxToGenome_Df(self.metatx,python_bin,map_to_genome,dir_temp)
        #del self.bed_mapped["bed_type"]
        os.system("rm -r "+dir_temp)
        self.metatx_mapped_score = None
        return self.metatx_mapped

    def ScoreByBigWig(self,bigwig,name_score,cores=5):
        
        if self.metatx_mapped_score is None:
            mapped_bed=self.metatx_mapped
            
        else:
            mapped_bed=self.metatx_mapped_score
    
        mapped_bed=DfScoresByBed(mapped_bed,bigwig,cores=cores)
    
        mapped_bed=mapped_bed.fillna(0)
    
        mapped_bed=mapped_bed.sort_values(["chrom","start","end","strand"])
        mapped_bed.loc[:,name_score]=mapped_bed.loc[:,"score"]
        
        self.metatx_mapped_score=mapped_bed
    
        return self.metatx_mapped_score 

    def EnrichmentAnalysis(self,bed_feature,name_bed,
                               ref_group,ref_group_name,
                               comp_group,comp_group_name,
                               type_analysis="kmers",
                               metatx_type="tx"):
            
            print()
            print("Enrichment Analysis:")
            print("BedInfo: ",name_bed)
            print("RefGroup: ",ref_group_name)
            print("ComparedGroup: ",comp_group_name)
            print("ATTENTION: type analysis ",type_analysis)
            print("ATTENTION: metatx type ",metatx_type)
            print()
            
            if metatx_type == "tx":
                bed=self.metatx.copy(deep=True)
                bed.loc[:,"group"]=None
                bed.loc[bed["chrom"].isin(ref_group),"group"]= ref_group_name
                bed.loc[bed["chrom"].isin(comp_group),"group"]= comp_group_name
                
                print(bed.groupby("group").size())
                
            if metatx_type == "gx":
                try:
                    bed=self.metatx_mapped.copy(deep=True)
                    bed.loc[:,"bin"]=bed.loc[:,"name"].apply(lambda x: x.split("_")[4]+"_"+x.split("_")[5])
                    bed.loc[:,"tx"]=bed.loc[:,"name"].apply(lambda x: x.split("_")[1])
                    bed.loc[:,"group"]=None
                    bed.loc[bed["tx"].isin(ref_group),"group"]= ref_group_name
                    bed.loc[bed["tx"].isin(comp_group),"group"]= comp_group_name
                except:
                    print("NO GENOME MAPPING")
    
            ### Types:
            # tx - tx with or without features
            # nucleotides
            # n feature bed1 nfeature bed2
            
            if type_analysis=="tx":
                if metatx_type == "gx":
                    ref_group_counts=CountsTxWithFeature(bed.loc[bed["group"]==ref_group_name],bed_feature,"tx","bin")
                    comp_group_counts=CountsTxWithFeature(bed.loc[bed["group"]==comp_group_name],bed_feature,"tx","bin")
                    
                else:
                    ref_group_counts=CountsTxWithFeature(bed.loc[bed["group"]==ref_group_name],bed_feature)
                    comp_group_counts=CountsTxWithFeature(bed.loc[bed["group"]==comp_group_name],bed_feature)
            
            if type_analysis=="kmers":
                ref_group_counts=KmersCountsInMetaTx(bed.loc[bed["group"]==ref_group_name],bed_feature,type_metatx=metatx_type)
                comp_group_counts=KmersCountsInMetaTx(bed.loc[bed["group"]==comp_group_name],bed_feature,type_metatx=metatx_type)
            
            ref_group_counts=ref_group_counts.reset_index()
            comp_group_counts=comp_group_counts.reset_index()
            
            ref_group_counts.columns=["bin",f"yes_{ref_group_name}",f"no_{ref_group_name}"]
            comp_group_counts.columns=["bin",f"yes_{comp_group_name}",f"no_{comp_group_name}"]
    
            counts_table=pd.merge(ref_group_counts,comp_group_counts,left_on="bin",right_on="bin")
    
            counts_table=counts_table.set_index("bin")
            counts_final=FisherExactTest_MultiTest(counts_table,ref_group_name,f"yes_{ref_group_name}",f"no_{ref_group_name}"
                                                              ,comp_group_name,f"yes_{comp_group_name}",f"no_{comp_group_name}")
            counts_final=counts_final.reset_index()
    
            counts_final.loc[:,"region"]=counts_final.loc[:,"bin"].apply(lambda x: x.split("_")[0])
            categoria=pd.CategoricalDtype(["5UTR","CDS","3UTR"],ordered=True)
            counts_final.loc[:,"region"]=counts_final.loc[:,"region"].astype(categoria)
            counts_final.loc[:,"number"]=counts_final.loc[:,"bin"].apply(lambda x: x.split("_")[1])
            counts_final.loc[:,"number"]=counts_final.loc[:,"number"].astype(int)
            counts_final=counts_final.sort_values(["region","number"])
            
            counts_final.loc[:,"analysis"]=name_bed+"_"+counts_final.loc[:,"analysis"]
            
            counts_final=counts_final.reset_index()
            del counts_final["index"]
            del counts_final["number"]
            del counts_final["region"]
            #counts_final=counts_final.set_index("bin")
            
            counts_final.loc[:,"type_analysis"]=type_analysis
    
            if self.enrichment is None:
    
                self.enrichment = counts_final
                
    
            else:
    
                self.enrichment = pd.concat([self.enrichment, counts_final], ignore_index=True)

                
            self.enrichment=self.enrichment.reset_index()
            del self.enrichment["index"]
    
            return counts_final

    def ScoreMatrix_Enrichment(self, col_sign="fdr", sign=0.05, col_score="log2_or", selected_analysis="all"):
    
        if self.enrichment is None:
            print("No Kmers Enrichment Analysis Available")
            return None
    
        matrix = self.enrichment.copy(deep=True)
        if selected_analysis != "all":
            matrix = matrix.loc[matrix["analysis"].isin(selected_analysis)]
    
        matrix["sign_enrichment"] = np.nan
        matrix.loc[matrix[col_sign] < sign, "sign_enrichment"] = matrix.loc[matrix[col_sign] < sign, col_score]
    
        matrix = matrix.pivot_table(index="bin", columns="analysis", values="sign_enrichment", dropna=False)
        matrix = matrix.reset_index()
        matrix = ReindexMetagene(matrix, "bin")
    
        # Rimuovi colonne temporanee
        matrix = matrix.drop(columns=["region", "n_bin", "index"], errors="ignore")
    
        if "freq" in col_score:
            matrix.columns = [col_score + "_" + str(x) for x in matrix.columns]
    
        if self.score_matrix is None:
            self.score_matrix = matrix
        else:
            self.score_matrix = pd.merge(self.score_matrix, matrix, left_index=True, right_index=True, how="outer")
            self.score_matrix = self.score_matrix.loc[matrix.index]
    
        return matrix


    def KmerBinFrequencyMatrix(self, metatx_bed, kmers_list, metatx_type="tx"):
        """
        Versione efficiente: UNA sola intersect per tutti i kmer.
        """
    
        # 1. Intersect unico
        kmers_bed = self.kmer_bed.copy()
        kmers_bed["end"] = kmers_bed["start"] + 1
    
        kmer_counts = IntersectBedDf(metatx_bed, kmers_bed, '-wao -s')
        kmer_counts.columns = (
            [c+"_1" for c in metatx_bed.columns] +
            [c+"_2" for c in kmers_bed.columns] +
            ["coverage"]
        )
    
        # 2. Calcolo tot_sum per ogni bin (come in KmersCountsInMetaTx)
        k = len(kmers_list[0])
        kmer_counts["length_bin"] = kmer_counts["end_1"] - kmer_counts["start_1"]
    
        # trova ultimo bin 3UTR
        bins = list(set(kmer_counts["name_1"]))
        utr3 = [b for b in bins if b.startswith("3UTR")]
        last_bin = max(utr3, key=lambda x: int(x.split("_")[1]))
    
        kmer_counts["tot_sum"] = kmer_counts.apply(
            lambda row: row["length_bin"] if row["name_1"] != last_bin else row["length_bin"] - (k - 1),
            axis=1
        )
        kmer_counts.loc[kmer_counts["tot_sum"] < 0, "tot_sum"] = 0
    
        # 3. Prepara struttura per output
        freq_matrix = pd.DataFrame(index=sorted(bins))
    
        # 4. Per ogni kmer → groupby in memoria (velocissimo)
        for kmer in kmers_list:
            df_k = kmer_counts[kmer_counts["name_2"] == kmer]
    
            yes = df_k.groupby("name_1")["coverage"].sum()
            tot = kmer_counts.groupby("name_1")["tot_sum"].sum()
    
            freq = (yes / tot).fillna(0)
            freq_matrix[kmer] = freq
    
        # 5. Ordina i bin biologicamente
        freq_matrix = freq_matrix.reset_index().rename(columns={"index": "bin"})
        freq_matrix = ReindexMetagene(freq_matrix, "bin")
        del freq_matrix['n_bin']
        del freq_matrix['region'] 
        del freq_matrix['index']
        self.kmers_frequencies = freq_matrix

        return freq_matrix

    
    def PlotKmerMatrixFreq(self, n_clusters, figsize=(50, 25)):
        def sort_key(idx):
            region, num = idx.split("_")
            region_order = {"5UTR": 0, "CDS": 1, "3UTR": 2}
            return (region_order[region], int(num))
            
        df_final = self.kmers_frequencies
        df_final = df_final.loc[sorted(df_final.index, key=sort_key)]
        # --------------------------------------------------
        # 2) Z-SCORE ROW-WISE
        # --------------------------------------------------
        df_z = df_final.sub(df_final.mean(axis=1), axis=0).div(df_final.std(axis=1), axis=0)
        df_z = df_z.fillna(0)
        # --------------------------------------------------
        # 3) HIERARCHICAL CLUSTERING OF KMERS (BINS)
        # --------------------------------------------------
        agg = AgglomerativeClustering(n_clusters=n_clusters, metric="euclidean", linkage="average")
        bin_clusters = agg.fit_predict(df_z.T)
        # --------------------------------------------------
        # 4) CLUSTER COLORS (ALIGNED WITH KMERS)
        # --------------------------------------------------
        cluster_palette = sns.color_palette("tab10", n_clusters)
        row_colors = pd.Series(bin_clusters,index=df_z.columns).map(dict(enumerate(cluster_palette)))
        # --------------------------------------------------
        # 5) CLUSTERMAP WITH FIXED Z-SCORE COLOR SCALE
        # --------------------------------------------------
        sns.set(context="notebook", style="white")
        g = sns.clustermap(
            df_z.T,
            row_cluster=True,
            col_cluster=False,
            row_colors=row_colors,
            cmap="RdBu_r",     # blue -> white -> red
            center=0,
            vmin=-3,
            vmax=3,
            xticklabels=True,
            yticklabels=True,
            figsize=figsize)
        
        # --------------------------------------------------
        # 6) CLUSTER LEGEND
        # --------------------------------------------------
        for i, color in enumerate(cluster_palette):
            g.ax_row_dendrogram.bar(0, 0, color=color, label=f"Cluster {i}", linewidth=0)
        
        g.ax_row_dendrogram.legend(title="Bin clusters",loc="center",ncol=1)
        #plt.savefig("",dpi=300)
        plt.show()

        kmer_cluster_map = pd.DataFrame({
        "kmer": df_z.columns,
        "cluster": bin_clusters})
        self.kmer_cluster_map = kmer_cluster_map
        return kmer_cluster_map

 

    def AnnotateMiRNASeeds(self, kmers=None, source="KEA", reference=None, species=None, species_id=9606, plot=True, top_n=20, output_dir=None):

        """Annotate RNA-oriented 7-mers against miRNA positions 2-8.
    
        Exact target 7mer-m8 matches only, using reverse-complemented Seed+m8.
        KEA defaults to the enriched delta_median/stat_log2fc intersection;
        KRS defaults to its selected signature and target percentile scores.
        source='manual' accepts explicit kmers without analysis scores.
        No new miRNA-level hypothesis test or combined KEA score is computed.
        """
    
        if source not in {"KEA", "KRS", "manual"}:
            raise ValueError("source must be KEA, KRS or manual")
        if not isinstance(top_n, (int, np.integer)) or top_n < 1:
            raise ValueError("top_n must be a positive integer")
        if output_dir is not None and Path(output_dir).exists():
            raise FileExistsError("output_dir must be a new directory")
        key = species or "input"
        family = pd.read_csv(self.seed_family_file, sep="\t")
        required = {"miR family", "Seed+m8", "Species ID", "MiRBase ID"}
        if not required.issubset(family.columns):
            raise ValueError(f"Family file requires columns: {sorted(required)}")
        family = family.loc[pd.to_numeric(family["Species ID"], errors="coerce").eq(int(species_id))].copy()
        if family.empty:
            raise ValueError("No family records for selected species_id")
        family["seed_rna"] = family["Seed+m8"].astype(str).str.upper().str.replace("T", "U", regex=False)
        if not family.seed_rna.map(lambda s: bool(re.fullmatch("[ACGU]{7}", s))).all():
            raise ValueError("Seed+m8 must contain canonical 7-nt seeds")
        if "Mature sequence" in family:
            mature = family["Mature sequence"].astype(str).str.upper().str.replace("T", "U", regex=False)
            if not mature.str[1:8].eq(family.seed_rna).all():
                raise ValueError("Seed+m8 does not match mature miRNA positions 2-8")
    
        family["target_kmer"] = family.seed_rna.map(lambda s: str(Seq(s.replace("U", "T")).reverse_complement()))
        family["target_rna"] = family.target_kmer.str.replace("T", "U", regex=False)
    
        def names(series):
            return ";".join(sorted(set(series.dropna().astype(str))))
        aggregations = {"mirnas": ("MiRBase ID", names)}
    
        if "MiRBase Accession" in family:
            aggregations["accessions"] = ("MiRBase Accession", names)
        annotations = family.groupby(["miR family", "seed_rna", "target_kmer", "target_rna"], as_index=False).agg(**aggregations)
        annotations = annotations.rename(columns={"miR family": "mirna_family"})
        scores = None
        state = (self.KEA_results or {}).get(key, {})
    
        if source == "KEA":
            if not {"delta_median", "stat_log2fc"}.issubset(state):
                raise ValueError("Run both delta_median and stat_log2fc before KEA annotation")
            if kmers is None:
                kmers = sorted(set(state["delta_median"]["enriched"]) & set(state["stat_log2fc"]["enriched"]))
            delta = state["delta_median"]["all_kmers"]
            stat = state["stat_log2fc"]["all_kmers"]
            scores = delta[["delta_median"]].join(stat[[c for c in ["log2FoldChange", "PValue", "FDR", "bonf"] if c in stat]], how="outer")
        elif source == "KRS":
            krs = state.get("KRS")
            if not isinstance(krs, dict):
                raise ValueError("Run KRS before KRS annotation")
            if reference is not None and reference != krs["reference"]:
                raise ValueError("Stored KRS result belongs to a different reference; rerun KRS")
            reference = krs["reference"]
            if kmers is None:
                kmers = krs["signature"]
            scores = krs["kmer_ranks"].set_index("kmer")[["target_frequency", "target_rank_percentile", "selected_KRS"]]
        elif kmers is None:
            raise ValueError("Supply kmers for source=manual")
        if isinstance(kmers, str):
            kmers = [kmers]
        words = list(dict.fromkeys(str(k).upper().replace("U", "T") for k in kmers))
        if any(not re.fullmatch("[ACGT]{7}", k) for k in words):
            raise ValueError("Query words must be canonical 7-mers; no automatic decomposition")
        query = pd.DataFrame({"kmer": pd.Series(words, dtype=str)})
        if scores is not None:
            scores = scores.copy()
            scores.index = scores.index.astype(str).str.upper().str.replace("U", "T", regex=False)
            if not scores.index.is_unique:
                raise ValueError("Analysis scores must have unique normalized k-mer rows")
            query = query.merge(scores, left_on="kmer", right_index=True, how="left", validate="one_to_one")
        matches = query.merge(annotations, left_on="kmer", right_on="target_kmer", how="inner", validate="one_to_many")
        query["matched_family_count"] = query.kmer.map(matches.groupby("kmer").mirna_family.nunique()).fillna(0).astype(int)
        figures = []
        if plot and len(matches) and source != "manual":
            import matplotlib.pyplot as plt
            metrics = ["delta_median", "log2FoldChange"] if source == "KEA" else ["target_rank_percentile"]
            sorting = "log2FoldChange" if source == "KEA" else "target_rank_percentile"
            shown = matches.sort_values(sorting, ascending=False, na_position="last").head(top_n).copy()
            shown["label"] = shown.mirna_family + " | " + shown.kmer
            fig, axes = plt.subplots(1, len(metrics), figsize=(7*len(metrics), max(3, 0.35*len(shown)+1)), squeeze=False)
            for ax, metric in zip(axes[0], metrics):
                ax.barh(range(len(shown)), shown[metric], color="#4477AA")
                ax.set_yticks(range(len(shown)), shown.label)
                ax.invert_yaxis()
                ax.set_xlabel(metric)
                ax.set_title("Matched k-mer score (not family-level inference)")
                if source == "KRS":
                    ax.set_xlim(0, 1)
            fig.tight_layout()
            figures.append(fig)
        base = Path(self.dir_out)
        if output_dir is None:
            base.mkdir(parents=True, exist_ok=True)
            destination = Path(tempfile.mkdtemp(prefix="MiRNASeeds_", dir=base))
        else:
            destination = Path(output_dir)
            destination.mkdir(parents=True, exist_ok=False)
        query.to_csv(destination/"query_kmers.tsv", sep="\t", index=False)
        matches.to_csv(destination/"seed_matches.tsv", sep="\t", index=False)
        annotations.to_csv(destination/"family_seed_map.tsv", sep="\t", index=False)
        for fig in figures:
            fig.savefig(destination/"seed_scores.png", dpi=180, bbox_inches="tight")
            fig.savefig(destination/"seed_scores.pdf", bbox_inches="tight")
        settings = dict(family_file=str(self.seed_family_file), species_id=int(species_id), source=source,
                        reference=reference, species=key, query_kmers=words,
                        site_type="7mer-m8", seed_orientation="miRNA positions 2-8; reverse complement to target",
                        plot=plot, top_n=top_n)
        (destination/"settings.json").write_text(json.dumps(settings, indent=2))
        result = dict(matches=matches, query=query, seed_map=annotations, figures=figures,
                      settings=settings, output_dir=str(destination))
        self.mirna_seed_annotation = result
        return result

    def K_CDS(
    self, gtf_file=None, uniprot_dir=None,
    kmers=None, method='stat_log2fc', species=None,
    in_frame=False, background='shuffle', tracks=("unipLocCytopl","unipLocExtra","unipLocSignal","unipLocTransMemb",),
    alpha=0.05, fdr_scope='global', split_blocks=True,
    seed=42, output_dir=None, plot=True, max_domains=60,
    ref_fasta=None, ctrl_fasta=None, *,
    positional_enrichment=None, bed_feature=None,
    recompute_positional=False, final_plots=None, rscript='Rscript',
    frame_threshold=0.5, start_fraction=0.10, color_cap=2,
    show_final=True):

        """
    
            Domain/frame tables and statistics are preserved. Only the final R panels are plotted.
            final_plots defaults to plot. R, ComplexHeatmap and circlize are required
            only for final panels. Pass final_plots=False for the original behaviour.
            positional_enrichment accepts a TSV path or a DataFrame. When omitted,
            search the last result directory, dir_out/K_CDS and dir_out for the TSV.
            A supplied cache is validated and reused; it is never silently recomputed.
            Missing positional data are computed automatically from existing BED and
            metatranscript inputs, or from FASTA/GTF if those inputs are unavailable.
            Previously computed per-k-mer EnrichmentAnalysis results are reused.
            recompute_positional=True explicitly requests this positional calculation.
            Frame/domain tests still compare reference with the chosen background.
            """

        if gtf_file is None or uniprot_dir is None:
            raise ValueError('Supply gtf_file and uniprot_dir for matching GRCh38 resources')
        if background not in {'control', 'shuffle'}:
            raise ValueError("background must be 'control' or 'shuffle'")
        if fdr_scope not in {'global', 'kmer'} or not 0 < alpha < 1:
            raise ValueError('Invalid fdr_scope or alpha')
        if max_domains is not None and max_domains < 1:
            raise ValueError('max_domains must be positive or None')
        if output_dir is not None and Path(output_dir).exists():
            raise FileExistsError(f'Refusing to overwrite existing directory: {output_dir}')
        key = species or 'input'
        if kmers is None:
            try:
                kmers = self.KEA_results[key][method]['enriched']
            except (KeyError, TypeError):
                raise ValueError('Run ExtractKmers first or supply kmers explicitly') from None
        if isinstance(kmers, str):
            kmers = [kmers]
        kmers = list(dict.fromkeys((str(k).upper().replace('U', 'T') for k in kmers)))
        if not kmers or any((not re.fullmatch('[ACGT]+', k) for k in kmers)):
            raise ValueError('Supply at least one DNA/RNA k-mer without ambiguous bases')
        if key != 'input':
            raise ValueError('K-CDS supports the human input KEA object; species must be None or input')
        if isinstance(tracks, str):
            tracks = (tracks,)
        allowed = {'unipDomain', 'unipInterest', 'unipStruct', 'unipLocTransMemb', 'unipLocExtra', 'unipLocSignal', 'unipLocCytopl', 'unipRepeat'}
        if not tracks or not set(tracks) <= allowed:
            raise ValueError('Unsupported track; use a UniProt feature track')

        final_plots = bool(plot) if final_plots is None else bool(final_plots)
        positional = None
        positional_source = None
        r_executable = None
        if final_plots:
            if len(kmers) < 2:
                raise ValueError('At least two k-mers are required for the final clusters')
            if not np.isfinite(color_cap) or color_cap <= 0:
                raise ValueError('color_cap must be finite and positive')
            if not np.isfinite(frame_threshold) or not 0 <= frame_threshold <= 1:
                raise ValueError('frame_threshold must be between 0 and 1')
            if not np.isfinite(start_fraction) or not 0 < start_fraction < 1:
                raise ValueError('start_fraction must be between 0 and 1')
            r_executable = shutil.which(str(rscript))
            if r_executable is None:
                raise RuntimeError('Rscript not found. Install R or set final_plots=False.')
            probe = subprocess.run(
                [r_executable, '--vanilla', '-e',
                 'if (!all(vapply(c("ComplexHeatmap", "circlize"), requireNamespace, logical(1), quietly=TRUE))) stop("Install ComplexHeatmap and circlize")'],
                capture_output=True, text=True
            )
            if probe.returncode:
                raise RuntimeError('R dependency check failed:\n' + probe.stdout + probe.stderr)
            source = positional_enrichment
            if source is None and not recompute_positional:
                # Only per-k-mer analyses are reusable: a pooled signature is not.
                existing = getattr(self, 'enrichment', None)
                memory_parts = []
                if isinstance(existing, pd.DataFrame) and {'analysis','bin','log2_or'}.issubset(existing.columns):
                    for kmer in kmers:
                        part = existing.loc[existing['analysis'].eq(kmer + '_Enr_VS_Ctrl')].copy()
                        if 'type_analysis' in part:
                            part = part.loc[part['type_analysis'].eq('kmers')]
                        if part.empty:
                            continue
                        # EnrichmentAnalysis appends results on repeated calls.
                        part = part.drop_duplicates('bin', keep='last')
                        cds_part = part.loc[part['bin'].astype(str).str.fullmatch(r'CDS_\d+')]
                        if len(cds_part) < 2 or not np.isfinite(pd.to_numeric(cds_part['log2_or'], errors='coerce')).all():
                            continue
                        part['kmer'] = kmer
                        if 'fdr' not in part:
                            part['fdr'] = np.nan
                        memory_parts.append(part[['kmer','bin','log2_or','fdr']])
                if memory_parts:
                    source = pd.concat(memory_parts, ignore_index=True)
                    positional_source = 'self.enrichment: per-k-mer Enr/Ctrl results'
                else:
                    previous = getattr(self, 'k_cds', None)
                    if isinstance(previous, dict) and isinstance(previous.get('positional_enrichment'), pd.DataFrame):
                        source = previous['positional_enrichment']
                    candidates = []
                    if isinstance(previous, dict) and previous.get('output_dir'):
                        candidates.append(Path(previous['output_dir']) / 'positional_enrichment.tsv')
                    candidates.extend([
                        Path(self.dir_out) / 'K_CDS' / 'positional_enrichment.tsv',
                        Path(self.dir_out) / 'positional_enrichment.tsv'
                    ])
                    if source is None:
                        source = next((path for path in candidates if path.is_file()), None)
            if not recompute_positional and source is not None:
                if isinstance(source, pd.DataFrame):
                    positional = source.copy()
                    positional_source = positional_source or 'supplied DataFrame'
                else:
                    cache = Path(source).expanduser()
                    if not cache.is_file():
                        raise FileNotFoundError(cache)
                    positional = pd.read_csv(cache, sep='\t')
                    positional_source = str(cache.resolve())
                print('Reusing positional enrichment:', positional_source)

        def txid(s):
            return re.sub('\\.\\d+$', '', str(s))

        def fasta(path):
            result = {}
            op = gzip.open if str(path).endswith('.gz') else open
            with op(path, 'rt') as handle:
                for rec in SeqIO.parse(handle, 'fasta'):
                    name = txid(rec.id)
                    if name in result:
                        raise ValueError(f'Duplicate normalized FASTA transcript ID: {name}')
                    result[name] = str(rec.seq).upper().replace('U', 'T')
            return result
        reference_path = ref_fasta if ref_fasta is not None else self.reference_fasta
        control_path = ctrl_fasta if ctrl_fasta is not None else self.control_fasta
        ref = fasta(reference_path)
        ctrl = fasta(control_path) if background == 'control' else {}
        if set(ref) & set(ctrl):
            raise ValueError('Reference and control transcript sets must be disjoint')
        seqs = {**ref, **ctrl}
        if final_plots:
            # Determine which per-k-mer profiles are missing/incomplete.
            reusable = {}
            if positional is not None:
                required = {'kmer','bin','log2_or'}
                if not required.issubset(positional.columns):
                    raise ValueError('Positional input requires kmer, bin and log2_or')
                positional = positional.copy()
                positional['kmer'] = positional['kmer'].astype(str).str.upper().str.replace('U','T',regex=False)
                if positional.duplicated(['kmer','bin']).any():
                    raise ValueError('Duplicate positional k-mer/bin pairs')
                for kmer in kmers:
                    part = positional.loc[positional['kmer'].eq(kmer)].copy()
                    cds_part = part.loc[part['bin'].astype(str).str.fullmatch(r'CDS_\d+')]
                    if len(cds_part) >= 2 and np.isfinite(pd.to_numeric(cds_part['log2_or'],errors='coerce')).all():
                        if 'fdr' not in part:
                            part['fdr'] = np.nan
                        reusable[kmer] = part[['kmer','bin','log2_or','fdr']]
            # Require a shared, complete set of CDS bins.
            expected_bins = set()
            for part in reusable.values():
                expected_bins.update(part.loc[part['bin'].astype(str).str.fullmatch(r'CDS_\d+'),'bin'])
            meta = getattr(self,'metatx',None)
            if isinstance(meta,pd.DataFrame) and 'name' in meta:
                expected_bins.update(meta.loc[meta['name'].astype(str).str.fullmatch(r'CDS_\d+'),'name'])
            reusable = {k:part for k,part in reusable.items()
                        if set(part.loc[part['bin'].astype(str).str.fullmatch(r'CDS_\d+'),'bin']) == expected_bins}
            missing = [k for k in kmers if k not in reusable]
            if positional_enrichment is not None and not recompute_positional and missing:
                raise ValueError('Explicit positional input is incomplete for: ' + ', '.join(missing))
            if missing:
                if any(len(kmer) != 7 for kmer in missing):
                    raise ValueError('Automatic EnrichmentAnalysis positional counting uses k=7. Supply positional_enrichment for other word lengths.')
                # Load full transcripts from both populations, even with shuffle.
                full_sequences = {}
                group_ids = []
                for path in [reference_path,control_path]:
                    ids = []
                    op = gzip.open if str(path).endswith('.gz') else open
                    with op(path,'rt') as handle:
                        for rec in SeqIO.parse(handle,'fasta'):
                            if rec.id in full_sequences:
                                raise ValueError('Duplicate/disjoint FASTA ID problem: ' + rec.id)
                            ids.append(rec.id)
                            full_sequences[rec.id] = str(rec.seq).upper().replace('U','T')
                    group_ids.append(ids)
                ref_ids,ctrl_ids = group_ids
                # Prefer the same metatranscript already used by the notebook.
                if meta is None:
                    candidates = list(Path(self.dir_out).glob('metatx*.bed'))
                    valid_meta = []
                    for path in candidates:
                        table_meta = pd.read_csv(path,sep='\t')
                        if {'chrom','start','end','name','score','strand'}.issubset(table_meta.columns):
                            if table_meta['name'].astype(str).str.fullmatch(r'CDS_\d+').any():
                                valid_meta.append((path,table_meta))
                    if len(valid_meta) == 1:
                        meta = valid_meta[0][1]
                        print('Using saved metatranscript:',valid_meta[0][0])
                if meta is None:
                    # No previous metatranscript: build relative-position bins
                    # from the matching GTF, using mean regional proportions.
                    normalized_ids = {txid(name):name for name in full_sequences}
                    if len(normalized_ids) != len(full_sequences):
                        raise ValueError('Duplicate normalized FASTA IDs')
                    annotations = defaultdict(lambda:{'exon':[],'CDS':[]})
                    op = gzip.open if str(gtf_file).endswith('.gz') else open
                    with op(gtf_file,'rt') as handle:
                        for line in handle:
                            if line.startswith('#'):
                                continue
                            fields = line.rstrip('\n').split('\t')
                            if len(fields)!=9 or fields[2] not in {'exon','CDS'}:
                                continue
                            hit = re.search(r'transcript_id "([^"]+)"',fields[8])
                            if not hit or txid(hit[1]) not in normalized_ids:
                                continue
                            model = annotations[normalized_ids[txid(hit[1])]]
                            model['strand'] = fields[6]
                            model[fields[2]].append((int(fields[3])-1,int(fields[4])))
                    lengths = {}
                    for name,sequence in full_sequences.items():
                        model = annotations.get(name)
                        if model is None or not model['exon'] or not model['CDS']:
                            continue
                        exons = sorted(set(model['exon']),reverse=model['strand']=='-')
                        offsets,offset = [],0
                        for lo,hi in exons:
                            offsets.append((lo,hi,offset)); offset += hi-lo
                        if offset != len(sequence):
                            continue
                        intervals = []
                        for lo,hi in set(model['CDS']):
                            for left,right,off in offsets:
                                if left<=lo<hi<=right:
                                    start = off+(lo-left if model['strand']=='+' else right-hi)
                                    intervals.append((start,start+hi-lo)); break
                        intervals = sorted(intervals)
                        if len(intervals)!=len(set(model['CDS'])) or any(right!=left for (_,right),(left,_) in zip(intervals,intervals[1:])):
                            continue
                        cs,ce = intervals[0][0],intervals[-1][1]
                        lengths[name] = (cs,ce-cs,len(sequence)-ce)
                    if not lengths:
                        raise ValueError('Cannot build metatranscript from the supplied FASTA/GTF')
                    means = np.mean(list(lengths.values()),axis=0)
                    n_bins = np.maximum(1,np.rint(100*means/means.sum()).astype(int)*3)
                    meta_rows = []
                    for name,region_lengths in lengths.items():
                        offset = 0
                        for label,size,count in zip(['5UTR','CDS','3UTR'],region_lengths,n_bins):
                            for j in range(int(count)):
                                lo = int(round(offset+j*size/count))
                                hi = int(round(offset+(j+1)*size/count))
                                if lo<hi:
                                    meta_rows.append((name,lo,hi,label+'_'+str(j+1),1,'+'))
                            offset += size
                    meta = pd.DataFrame(meta_rows,columns=['chrom','start','end','name','score','strand'])
                    print('Built metatranscript from FASTA/GTF:',n_bins.tolist())
                self.metatx = meta
                target_bins = set(meta.loc[meta['name'].astype(str).str.fullmatch(r'CDS_\d+'),'name'])
                # Partial cached profiles must use the same binning as newly
                # computed ones; recompute incompatible profiles as needed.
                reusable = {k:part for k,part in reusable.items()
                            if set(part.loc[part['bin'].astype(str).str.fullmatch(r'CDS_\d+'),'bin']) == target_bins}
                missing = [k for k in kmers if k not in reusable]
                # Use existing occurrence coordinates, or scan overlapping hits.
                if bed_feature is not None:
                    features = bed_feature.copy() if isinstance(bed_feature,pd.DataFrame) else pd.read_csv(bed_feature,sep='\t')
                elif (Path(self.dir_out)/'kmers_bed.bed').is_file():
                    features = pd.read_csv(Path(self.dir_out)/'kmers_bed.bed',sep='\t')
                else:
                    feature_rows = []
                    for name,sequence in full_sequences.items():
                        for kmer in missing:
                            pos = sequence.find(kmer)
                            while pos>=0:
                                feature_rows.append((name,pos,pos+len(kmer),kmer,1,'+'))
                                pos = sequence.find(kmer,pos+1)
                    features = pd.DataFrame(feature_rows,columns=['chrom','start','end','name','score','strand'])
                if not {'chrom','start','end','name','score','strand'}.issubset(features.columns):
                    raise ValueError('Invalid k-mer BED table')
                features['name'] = features['name'].astype(str).str.upper().str.replace('U','T',regex=False)
                old_enrichment = getattr(self,'enrichment',None)
                try:
                    for kmer in missing:
                        print('Computing missing positional profile:',kmer)
                        piece = self.EnrichmentAnalysis(
                            features.loc[features['name'].eq(kmer)].copy(), kmer,
                            ref_ids,'Enr',ctrl_ids,'Ctrl',
                            type_analysis='kmers',metatx_type='tx'
                        ).copy()
                        piece['kmer'] = kmer
                        reusable[kmer] = piece[['kmer','bin','log2_or','fdr']]
                finally:
                    self.enrichment = old_enrichment
                positional_source = (positional_source or '') + '; computed missing profiles with EnrichmentAnalysis'
            positional = pd.concat([reusable[k] for k in kmers],ignore_index=True)
            required = {'kmer', 'bin', 'log2_or'}
            if not required.issubset(positional.columns):
                raise ValueError('Positional columns missing: ' + str(sorted(required - set(positional.columns))))
            if positional[['kmer', 'bin']].isna().any().any():
                raise ValueError('Missing positional k-mer/bin identifiers')
            positional['kmer'] = positional['kmer'].astype(str).str.upper().str.replace('U', 'T', regex=False)
            positional['bin'] = positional['bin'].astype(str)
            positional = positional.loc[positional['kmer'].isin(kmers)].copy()
            if positional.duplicated(['kmer', 'bin']).any():
                raise ValueError('Duplicate positional k-mer/bin pairs')
            if set(positional['kmer']) != set(kmers):
                raise ValueError('Positional data do not contain all selected k-mers')
            positional['log2_or'] = pd.to_numeric(positional['log2_or'], errors='raise')
            cds_rows = positional[positional['bin'].str.fullmatch(r'CDS_\d+', na=False)]
            cds_bins = sorted(cds_rows['bin'].unique(), key=lambda b: int(b.split('_')[1]))
            positional_matrix = cds_rows.pivot(index='kmer', columns='bin', values='log2_or').reindex(index=kmers, columns=cds_bins)
            if len(cds_bins) < 2 or not np.isfinite(positional_matrix.to_numpy(dtype=float)).all():
                raise ValueError('Need at least two CDS bins and a complete finite positional matrix')

        models = defaultdict(lambda: {'exon': [], 'CDS': []})
        op = gzip.open if str(gtf_file).endswith('.gz') else open
        with op(gtf_file, 'rt') as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) != 9 or f[2] not in {'exon', 'CDS'}:
                    continue
                match = re.search('transcript_id "([^"]+)"', f[8])
                if not match or txid(match[1]) not in seqs:
                    continue
                name = txid(match[1])
                chrom = f[0] if f[0].startswith('chr') else 'chr' + f[0]
                if chrom == 'chrMT':
                    chrom = 'chrM'
                m = models[name]
                if 'chrom' in m and (m['chrom'], m['strand']) != (chrom, f[6]):
                    raise ValueError(f'Conflicting GTF loci for {name}')
                m.update(chrom=chrom, strand=f[6])
                m[f[2]].append((int(f[3]) - 1, int(f[4]), f[7]))
        bins = defaultdict(list)
        bin_size = 100000
        for track in dict.fromkeys(tracks):
            with open(Path(uniprot_dir) / (track + '.bed')) as handle:
                for line in handle:
                    if not line.strip() or line.startswith(('#', 'track', 'browser')):
                        continue
                    f = line.rstrip('\n').split('\t')
                    a, b = (int(f[1]), int(f[2]))
                    label = f[26] if track in {'unipDomain', 'unipInterest'} and len(f) > 26 else f[3]
                    if track == 'unipInterest' and label != 'Disordered':
                        continue
                    if track == 'unipLocSignal' and label != 'Signal peptide':
                        continue
                    label = label.strip()
                    if not label:
                        continue
                    if len(set(tracks)) > 1:
                        label = track + ': ' + label
                    blocks = [(a, b)]
                    if split_blocks and len(f) >= 12:
                        sizes = [int(x) for x in f[10].rstrip(',').split(',')]
                        starts = [int(x) for x in f[11].rstrip(',').split(',')]
                        if len(sizes) != int(f[9]) or len(starts) != len(sizes):
                            raise ValueError('Malformed UniProt BED12 block count')
                        blocks = [(a + s, a + s + n) for s, n in zip(starts, sizes)]
                    for lo, hi in blocks:
                        if not a <= lo < hi <= b:
                            raise ValueError('Invalid UniProt BED block coordinates')
                        for bucket in range(lo // bin_size, (hi - 1) // bin_size + 1):
                            bins[f[0], f[5], bucket].append((lo, hi, label))
        rows, qc, frame_rows = ([], [], [])
        frame_rng = np.random.default_rng(seed)
        rng = np.random.default_rng(seed)
        for name, seq in seqs.items():
            m = models.get(name)
            reason = None
            if m is None or not m['CDS']:
                reason = 'missing_GTF_or_non_coding'
            elif not m['exon']:
                reason = 'missing_exons'
            if reason:
                qc.append((name, reason))
                continue
            exons = sorted(set(((a, b) for a, b, _ in m['exon'])), reverse=m['strand'] == '-')
            genomic_order = sorted(exons)
            if any((b > c for (a, b), (c, d) in zip(genomic_order, genomic_order[1:]))):
                raise ValueError(f'Overlapping GTF exons for {name}')
            offsets, offset = ([], 0)
            for a, b in exons:
                offsets.append((a, b, offset))
                offset += b - a
            if offset != len(seq):
                qc.append((name, 'FASTA_GTF_length_mismatch'))
                continue
            cds = []
            for a, b, phase in m['CDS']:
                contained = False
                for lo, hi, off in offsets:
                    if lo <= a < b <= hi:
                        start = off + (a - lo if m['strand'] == '+' else hi - b)
                        cds.append((start, start + b - a, int(phase)))
                        contained = True
                        break
                if not contained:
                    raise ValueError(f'CDS is not contained in an exon: {name}')
            cds = sorted(set(cds))
            if any((b != c for (a, b, _), (c, d, _) in zip(cds, cds[1:]))):
                qc.append((name, 'non_contiguous_CDS'))
                continue
            cs, ce, phase0 = (cds[0][0], cds[-1][1], cds[0][2])
            anchor = cs + phase0
            if any((p != (anchor - a) % 3 for a, b, p in cds)):
                qc.append((name, 'inconsistent_CDS_phase'))
                continue
            projected = set()
            for lo, hi, off in offsets:
                candidates = set()
                for bucket in range(lo // bin_size, (hi - 1) // bin_size + 1):
                    candidates.update(bins.get((m['chrom'], m['strand'], bucket), ()))
                for a, b, label in candidates:
                    left, right = (max(lo, a), min(hi, b))
                    if left < right:
                        start = off + (left - lo if m['strand'] == '+' else hi - right)
                        projected.add((start, start + right - left, label))

            def record(k, pos, group):
                domains = tuple(sorted({label for a, b, label in projected if a < pos + len(k) and pos < b}))
                rows.append((name, k, pos, pos + len(k), (pos - anchor) % 3, group, domains))
            for k in kmers:
                pos = seq.find(k, cs, ce)
                while pos >= 0:
                    frame_rows.append((k, (pos - anchor) % 3, 'coding' if name in ref else 'control'))
                    if name in ref:
                        shuffled_frame_pos = cs + int(frame_rng.integers(ce - len(k) - cs + 1))
                        frame_rows.append((k, (shuffled_frame_pos - anchor) % 3, 'shuffle'))
                    if not in_frame or (pos >= anchor and (pos - anchor) % 3 == 0):
                        record(k, pos, 'reference' if name in ref else 'control')
                        if background == 'shuffle':
                            first = anchor if in_frame else cs
                            step = 3 if in_frame else 1
                            n_positions = (ce - len(k) - first) // step + 1
                            shuffled = first + step * int(rng.integers(n_positions))
                            record(k, shuffled, 'control')
                    pos = seq.find(k, pos + 1, ce)
            qc.append((name, 'used'))
        occurrences = pd.DataFrame(rows, columns=['transcript', 'kmer', 'start', 'end', 'frame', 'group', 'domains'])
        frame_data = pd.DataFrame(frame_rows, columns=['kmer', 'frame', 'group'])
        frame_index = pd.MultiIndex.from_product([kmers, ['coding', 'shuffle', 'control'], range(3)], names=['kmer', 'group', 'frame'])
        frame_counts = frame_data.groupby(['kmer', 'group', 'frame']).size().reindex(frame_index, fill_value=0).rename('count').reset_index()
        frame_totals = frame_counts.groupby(['kmer', 'group'])['count'].transform('sum')
        frame_counts['percent'] = 100 * frame_counts['count'] / frame_totals.replace(0, np.nan)
        quality = pd.DataFrame(qc, columns=['transcript', 'status'])
        if len(quality) and (quality.status != 'used').any():
            warnings.warn("Some transcripts excluded; inspect result['qc'] for reasons", stacklevel=2)
        if occurrences.empty:
            raise ValueError('No eligible CDS occurrences. Check IDs, full-length FASTAs, GTF and frame filter. ' + str(quality.status.value_counts().to_dict()))
        totals = Counter(zip(occurrences.kmer, occurrences.group))
        hits = Counter()
        for row in occurrences.itertuples(index=False):
            hits.update(((row.kmer, row.group, domain) for domain in row.domains))
        domains = sorted({domain for _, _, domain in hits})
        if not domains:
            raise ValueError('No same-strand UniProt overlaps; check assembly, coordinates and selected tracks')
        stats = []
        for k in kmers:
            n1, n0 = (totals[k, 'reference'], totals[k, 'control'])
            for domain in domains:
                a, c = (hits[k, 'reference', domain], hits[k, 'control', domain])
                b, d = (n1 - a, n0 - c)
                testable = n1 > 0 and n0 > 0 and (0 < a + c < n1 + n0)
                odds, p = fisher_exact([[a, b], [c, d]], alternative='two-sided') if testable else (np.nan, np.nan)
                delta = 0.5 if 0 in (a, b, c, d) else 0.0
                logor = np.log2((a + delta) * (d + delta) / ((b + delta) * (c + delta))) if testable else np.nan
                stats.append((k, domain, a, b, c, d, odds, logor, p, testable))
        table = pd.DataFrame(stats, columns=['kmer', 'domain', 'A', 'B', 'C', 'D', 'odds_ratio', 'log2_or', 'pvalue', 'testable'])
        table['fdr'] = np.nan
        groups = [table.index] if fdr_scope == 'global' else table.groupby('kmer').groups.values()
        for idx in groups:
            valid = table.loc[idx, 'pvalue'].dropna().index
            if len(valid):
                table.loc[valid, 'fdr'] = multipletests(table.loc[valid, 'pvalue'], method='fdr_bh')[1]
        table['significant'] = table.fdr.le(alpha)
        table['enriched'] = table.odds_ratio.gt(1)
        table['significant_enrichment'] = table.significant & table.enriched
        matrix = table.pivot(index='kmer', columns='domain', values='log2_or').reindex(kmers)
        qmatrix = table.pivot(index='kmer', columns='domain', values='fdr').reindex(kmers)
        figure = None
        frame_figure = None
        plotted = list(matrix.columns)
        if output_dir is None:
            base = Path(self.dir_out)
            base.mkdir(parents=True, exist_ok=True)
            destination = Path(tempfile.mkdtemp(prefix='K_CDS_', dir=base))
        else:
            destination = Path(output_dir)
            destination.mkdir(parents=True, exist_ok=False)
        frame_counts.to_csv(destination / 'frame_distribution.tsv', sep='\t', index=False)
        table.to_csv(destination / 'statistics.tsv', sep='\t', index=False)
        occurrences.to_csv(destination / 'occurrences.tsv', sep='\t', index=False)
        quality.to_csv(destination / 'qc.tsv', sep='\t', index=False)
        matrix.to_csv(destination / 'log2_or.tsv', sep='\t')
        qmatrix.to_csv(destination / 'fdr.tsv', sep='\t')
        import json
        settings = dict(gtf_file=str(gtf_file), uniprot_dir=str(uniprot_dir), kmers=kmers, method=method, in_frame=in_frame, background=background, tracks=list(tracks), alpha=alpha, fdr_scope=fdr_scope, split_blocks=split_blocks, seed=seed, reference_fasta=str(reference_path), control_fasta=str(control_path), plotted_domains=plotted, species=key, plot=plot, max_domains=max_domains)
        (destination / 'settings.json').write_text(json.dumps(settings, indent=2))
        result = dict(statistics=table, occurrences=occurrences, qc=quality, log2_or=matrix, fdr=qmatrix, figure=figure, frame_distribution=frame_counts, frame_figure=frame_figure, output_dir=str(destination), settings=settings)
        result['final_figures'] = {}
        result['cluster_assignments'] = None
        result['positional_enrichment'] = positional
        settings.update(final_plots=final_plots, positional_source=positional_source,
                        frame_threshold=frame_threshold, start_fraction=start_fraction,
                        color_cap=color_cap, recompute_positional=recompute_positional)
        if final_plots:
            positional.to_csv(destination / 'positional_enrichment.tsv', sep='\t', index=False)
            figures_dir = destination / 'figures'
            figures_dir.mkdir()
            r_code = r"""
        #!/usr/bin/env Rscript
        # Usage: Rscript heatmap_script.R K_CDS_DIRECTORY [OUTPUT_DIRECTORY]
        # Dependencies: BiocManager::install('ComplexHeatmap'); install.packages('circlize')
        suppressPackageStartupMessages({
          library(ComplexHeatmap)
          library(circlize)
          library(grid)
        })
        args <- commandArgs(trailingOnly=TRUE)
        if (!length(args)) stop('Usage: Rscript heatmap_script.R K_CDS_DIRECTORY [OUTPUT_DIRECTORY]')
        input <- args[1]
        out <- if (length(args)>1) args[2] else file.path(input, 'R_figures')
        dir.create(out, recursive=TRUE, showWarnings=FALSE)
        # Adjustable analysis/display parameters.
        argument_number <- function(i, default) {
          value <- if (length(args)>=i) suppressWarnings(as.numeric(args[i])) else default
          if (length(value)!=1L || !is.finite(value)) stop('Invalid numeric argument ', i)
          value
        }
        threshold <- argument_number(3, 0.5)     # Inclusive max frame fraction.
        start_fraction <- argument_number(4, 0.10)
        alpha <- argument_number(5, 0.05)
        CAP <- argument_number(6, 2)
        if (threshold<0 || threshold>1 || start_fraction<=0 || start_fraction>=1 || alpha<=0 || alpha>=1 || CAP<=0)
          stop('Invalid threshold, start fraction, alpha or colour cap')
        frame_colors <- c(First='lightgray', Second='white', Third='black')
        col_fun <- colorRamp2(c(-CAP, -CAP/2, 0, CAP/2, CAP),
                              c('#2166AC', '#92C5DE', '#FFFFFF', '#F4A582', '#B2182B'))
        read_tsv <- function(name, required) {
          path <- file.path(input, name)
          if (!file.exists(path)) stop('Missing input: ', path)
          x <- read.delim(path, check.names=FALSE, stringsAsFactors=FALSE)
          if (!all(required %in% names(x))) stop(name, ': required columns: ', paste(required, collapse=', '))
          x
        }
        frames <- read_tsv('frame_distribution.tsv', c('kmer','group','frame','count'))
        stats <- read_tsv('statistics.tsv', c('kmer','domain','log2_or','fdr'))
        pos <- read_tsv('positional_enrichment.tsv', c('kmer','bin','log2_or'))
        if (anyDuplicated(stats[c('kmer','domain')])) stop('Duplicate kmer/domain pairs')
        if (anyDuplicated(pos[c('kmer','bin')])) stop('Duplicate kmer/bin pairs')
        if (anyDuplicated(frames[c('kmer','group','frame')])) stop('Duplicate frame rows')
        if (anyNA(frames$count) || any(frames$count<0) || any(!frames$frame %in% 0:2)) stop('Invalid frame counts or frame labels')
        kmers <- unique(stats$kmer)
        if (length(kmers)<2) stop('At least two k-mers are required for two positional clusters')
        if (!all(kmers %in% pos$kmer)) stop('Positional data missing for some k-mers')
        bins <- unique(pos$bin[grepl('^CDS_[0-9]+$', pos$bin)])
        bins <- bins[order(as.integer(sub('CDS_', '', bins)))]
        if (length(bins)<2) stop('Need at least two CDS bins named CDS_0, CDS_1, ...')
        make_matrix <- function(d, rows, cols, column, column_id) {
          m <- matrix(NA_real_, length(rows), length(cols), dimnames=list(rows,cols))
          ri <- match(d$kmer, rows); ci <- match(d[[column_id]], cols)
          ok <- !is.na(ri) & !is.na(ci)
          m[cbind(ri[ok],ci[ok])] <- as.numeric(d[[column]][ok])
          m
        }
        p <- make_matrix(pos, kmers, bins, 'log2_or', 'bin')
        # Missing estimates must not silently become neutral enrichment for clustering.
        if (any(!is.finite(p))) stop('CDS positional matrix contains missing/nonfinite estimates. Supply complete, unmasked log2_or values; do not use the pooled or FDR-masked matrix.')
        # Ward groups are computed from capped, row-centred positional profiles.
        # The displayed heatmap retains capped enrichment values, without centring.
        p_cap <- pmin(pmax(p, -CAP), CAP)
        ns <- max(1L, min(ncol(p)-1L, ceiling(ncol(p)*start_fraction)))
        contrast_raw <- rowMeans(p[,seq_len(ns),drop=FALSE]) - rowMeans(p[,-seq_len(ns),drop=FALSE])
        contrast <- rowMeans(p_cap[,seq_len(ns),drop=FALSE]) - rowMeans(p_cap[,-seq_len(ns),drop=FALSE])
        profile <- p_cap - rowMeans(p_cap)
        hc <- hclust(dist(profile, method='euclidean'), method='ward.D2')
        cl <- cutree(hc, k=2)
        start_cluster <- as.integer(names(which.max(tapply(contrast, cl, mean))))
        group <- ifelse(cl == start_cluster, 'CDS start', 'CDS body')
        frame_matrix <- function(g) {
          d <- frames[frames$group==g,,drop=FALSE]
          m <- make_matrix(d, kmers, as.character(0:2), 'count', 'frame')
          if (anyNA(m)) stop('Missing frame rows for group ', g, '; include all three frames, even zero counts')
          total <- rowSums(m)
          m <- m / ifelse(total>0, total, NA_real_)
          colnames(m) <- c('First','Second','Third')
          m
        }
        cds <- frame_matrix('coding'); shuffle <- frame_matrix('shuffle')
        score <- apply(cds, 1, function(x) if (all(is.na(x))) NA_real_ else max(x))
        passes <- !is.na(score) & score >= threshold
        winner <- apply(cds, 1, function(x) if (all(is.na(x))) NA_character_ else paste(colnames(cds)[which(x==max(x))],collapse='/'))
        # Shared row order: CDS start first, CDS body second; decreasing In Frame Index.
        # All row-dependent vectors and matrices must follow exactly the same order.
        in_frame_index <- 100 * score
        ord <- order(factor(group, levels=c('CDS start','CDS body')), -in_frame_index, kmers)
        kmers <- kmers[ord]
        rna <- gsub('T', 'U', kmers)
        p <- p_cap[ord,,drop=FALSE]
        cds <- cds[ord,,drop=FALSE]; shuffle <- shuffle[ord,,drop=FALSE]
        score <- score[ord]; in_frame_index <- in_frame_index[ord]
        passes <- passes[ord]; winner <- winner[ord]
        contrast <- contrast[ord]; contrast_raw <- contrast_raw[ord]
        split <- factor(group[ord], levels=c('CDS start','CDS body'))
        message('Ward on capped, row-centred profiles: ', sum(split == 'CDS body'),
                ' CDS body; ', sum(split == 'CDS start'), ' CDS start')

        domains <- unique(stats$domain)
        domain_limit <- if (args[7]=='all') Inf else as.integer(args[7])
        if (length(domains) > domain_limit) {
          best_fdr <- vapply(domains, function(domain) {
            values <- stats$fdr[stats$domain==domain]
            if (all(is.na(values))) Inf else min(values, na.rm=TRUE)
          }, numeric(1))
          domains <- domains[order(best_fdr, seq_along(domains))][seq_len(domain_limit)]
        }
        d <- make_matrix(stats, kmers, domains, 'log2_or', 'domain')
        q <- make_matrix(stats, kmers, domains, 'fdr', 'domain')
        d[!is.finite(d)] <- NA_real_
        d[!passes,] <- NA_real_
        write.table(data.frame(kmer=kmers, cluster=as.character(split), score=score,
                              passes=passes, preferred_frame=winner,
                              start_minus_body=contrast,
                              start_minus_body_raw=contrast_raw,
                              in_frame_index=in_frame_index, kmer_rna=rna),
                    file.path(out,'cluster_assignments.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
        # Shared ordering is explicit: neither heatmap can independently reorder rows.
        base_hm <- function(mat, name, colors=col_fun,
                            row_title_gp=gpar(fontsize=8,fontface='bold'),
                            column_names_gp=gpar(fontsize=8), ...) Heatmap(mat, name=name, col=colors,
          cluster_rows=FALSE, cluster_columns=FALSE, cluster_row_slices=FALSE,
          row_order=seq_len(nrow(mat)), row_split=split, row_gap=unit(2,'mm'),
          row_title_rot=90, row_title_gp=row_title_gp,
          row_names_gp=gpar(fontsize=8), column_names_gp=column_names_gp,
          column_names_rot=90, border=TRUE, ...)
        flag <- function() rowAnnotation(
          `Frame preference`=anno_simple(ifelse(passes,'Pass','Fail'),
            col=c(Pass='black',Fail='white'), border=TRUE, simple_anno_size=unit(3,'mm')),
          annotation_name_rot=90)
        # Compact main panel matching the supplied reference; scientific calculations above are unchanged.
        cell_mm <- 4.2
        main_row_mm <- 2.0
        main_gap_mm <- 2 * (length(unique(split)) - 1L)
        main_col_fun <- colorRamp2(c(-CAP, -CAP/2, 0, CAP/2, CAP),
                                   c('#9B63E6', '#C7A2F7', '#FFFFFF', '#FFAB94', '#EF6548'))
        main_ns_color <- '#EEEEEE'
        main_not_frame_color <- '#D5D5D5'
        main_left <- rowAnnotation(
          `k-mer`=anno_text(rna, just='right', location=unit(1,'npc'),
            gp=gpar(fontsize=5.1), width=max_text_width(rna, gp=gpar(fontsize=5.1)) + unit(1,'mm')),
          `In Frame Index`=anno_lines(in_frame_index, ylim=c(30,100), extend=0,
            add_points=TRUE, pt_gp=gpar(col='black',fill='black'),
            size=unit(0.25,'mm'), pch=16, gp=gpar(col='black',lwd=0.45),
            border=TRUE, width=unit(17,'mm'),
            axis_param=list(at=c(50,100),labels=c('50','100'), side='bottom',
                            labels_rot=90, gp=gpar(fontsize=5))),
          `Frame strip`=anno_simple(ifelse(passes,'Pass','Fail'),
            col=c(Pass='black',Fail='white'), border=TRUE, simple_anno_size=unit(1.8,'mm')),
          gap=unit(c(1.0,0.8),'mm'), show_legend=FALSE,
          show_annotation_name=c(`k-mer`=FALSE,`In Frame Index`=TRUE,`Frame strip`=FALSE),
          annotation_name_side='top', annotation_name_rot=0,
          annotation_name_gp=gpar(fontsize=6))
        main <- base_hm(d, 'log2 OR', colors=main_col_fun,
          na_col=main_not_frame_color, column_title=NULL, show_row_names=FALSE,
          row_title_gp=gpar(fontsize=7,fontface='plain'),
          column_names_gp=gpar(fontsize=6),
          column_labels=sub('^[^:]+: *','',colnames(d)),
          left_annotation=main_left,
          width=unit(ncol(d)*cell_mm,'mm'),
          height=unit(nrow(d)*main_row_mm + main_gap_mm,'mm'),
          heatmap_legend_param=list(title='log2 OR',at=c(-CAP,0,CAP),
            legend_height=unit(17,'mm'),grid_width=unit(2,'mm'),
            title_gp=gpar(fontsize=6),labels_gp=gpar(fontsize=5)),
          rect_gp=gpar(col='#D8D8D8',lwd=0.25),
          cell_fun=function(j,i,x,y,w,h,fill) {
            if (!passes[i]) {
              grid.rect(x,y,w,h,gp=gpar(fill=main_not_frame_color,col='#D8D8D8',lwd=0.25))
            } else if (!is.finite(d[i,j]) || !is.finite(q[i,j])) {
              grid.rect(x,y,w,h,gp=gpar(fill=main_ns_color,col='#D8D8D8',lwd=0.25))
              grid.text('NA',x,y,gp=gpar(fontsize=3.8,col='#777777'))
            } else if (q[i,j]>alpha) {
              grid.rect(x,y,w,h,gp=gpar(fill=main_ns_color,col='#D8D8D8',lwd=0.25))
              grid.text('NS',x,y,gp=gpar(fontsize=3.8,col='#777777'))
            }
          })
        main_leg_cell <- Legend(title='Cell',
          labels=c(paste0('FDR <= ',format(alpha,trim=TRUE)), 'NS: not significant', 'not in frame'),
          legend_gp=gpar(fill=c('#FFAB94',main_ns_color,main_not_frame_color),col=NA),
          grid_height=unit(2,'mm'), grid_width=unit(2,'mm'),
          title_gp=gpar(fontsize=6),labels_gp=gpar(fontsize=5))
        main_leg_frame <- Legend(title='Frame strip',labels=c('in frame','not in frame'),
          legend_gp=gpar(fill=c('black','white'),col='black',lwd=0.4),
          grid_height=unit(2,'mm'),grid_width=unit(2,'mm'),
          title_gp=gpar(fontsize=6),labels_gp=gpar(fontsize=5))
        main_leg_threshold <- Legend(title='Dashed line',type='lines',background='white',
          labels=paste0('threshold ',format(100*threshold,trim=TRUE),'%'),
          legend_gp=gpar(col='#999999',lty=2,lwd=0.7),
          title_gp=gpar(fontsize=6),labels_gp=gpar(fontsize=5))
        # Supplementary: aligned positional profiles, frame proportions and frame index.
        # Means use the same capped, uncentred values displayed in the positional heatmap.
        supp_row_mm <- 2.0
        cluster_colors <- c(`CDS body`='#69A99B', `CDS start`='#E39A65')
        cluster_means <- sapply(names(cluster_colors), function(g) colMeans(p[split==g,,drop=FALSE]))
        write.table(data.frame(bin=colnames(p), cluster_means, check.names=FALSE),
          file.path(out,'cluster_mean_profiles.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
        profile_limits <- range(c(0, cluster_means), finite=TRUE)
        if (diff(profile_limits)==0) profile_limits <- profile_limits + c(-0.1,0.1)
        profile_ticks <- pretty(profile_limits, n=4)
        profile_limits <- range(profile_limits, profile_ticks)
        supp_top <- HeatmapAnnotation(
          `Mean log2 OR`=anno_lines(cluster_means, add_points=FALSE,
            gp=gpar(col=unname(cluster_colors),lwd=0.65),
            ylim=profile_limits, extend=0, border=TRUE, height=unit(15,'mm'),
            axis_param=list(at=profile_ticks,labels=format(profile_ticks,trim=TRUE),
                            side='left',gp=gpar(fontsize=5))),
          annotation_name_side='left',annotation_name_rot=90,
          annotation_name_gp=gpar(fontsize=5.5))
        bar <- function(m) {
          m[is.na(m)] <- 0
          anno_barplot(m*100, beside=FALSE,
            gp=gpar(fill=frame_colors,col='grey50',lwd=0.3),
            border=TRUE,bar_width=0.8,ylim=c(0,100),extend=0,width=unit(19,'mm'),
            axis_param=list(at=c(0,50,100),labels=c('0','50','100'),
                            side='top',labels_rot=0,gp=gpar(fontsize=5)))
        }
        anno <- rowAnnotation(`CDS (%)`=bar(cds),`Shuffle (%)`=bar(shuffle),
          `In Frame Index`=anno_lines(in_frame_index,ylim=c(30,100),extend=0,
            add_points=TRUE,pt_gp=gpar(col='black',fill='black'),
            size=unit(0.25,'mm'),pch=16,gp=gpar(col='black',lwd=0.45),
            border=TRUE,width=unit(19,'mm'),
            axis_param=list(at=c(40,60,80,100),labels=c('40','60','80','100'),
                            side='top',labels_rot=0,gp=gpar(fontsize=5))),
          annotation_name_rot=0,annotation_name_side='top',
          annotation_name_gp=gpar(fontsize=6),gap=unit(3,'mm'))
        supp_left <- rowAnnotation(
          `k-mer`=anno_text(rna,just='right',location=unit(1,'npc'),
            gp=gpar(fontsize=5.1),
            width=max_text_width(rna,gp=gpar(fontsize=5.1)) + unit(1,'mm')),
          show_annotation_name=FALSE)
        supp <- base_hm(p,'Position log2 OR',colors=main_col_fun,
          column_title="CDS positional enrichment (5' -> 3')",
          column_title_gp=gpar(fontsize=7,fontface='plain'),
          show_column_names=FALSE,show_row_names=FALSE,
          row_title_gp=gpar(fontsize=7,fontface='plain'),
          top_annotation=supp_top,left_annotation=supp_left,width=unit(62,'mm'),
          height=unit(nrow(p)*supp_row_mm + main_gap_mm,'mm'),
          heatmap_legend_param=list(title='Position log2 OR',at=c(-CAP,-CAP/2,0,CAP/2,CAP),
            legend_height=unit(23,'mm'),grid_width=unit(2,'mm'),
            title_gp=gpar(fontsize=6),labels_gp=gpar(fontsize=5)),
          rect_gp=gpar(col=NA)) + anno
        # Explicit labels distinguish absent frame hits from an observed zero fraction.
        missing_labels <- ifelse(is.na(score),'No CDS hits',
                                 ifelse(is.na(shuffle[,1]),'No shuffle hits',''))
        if (any(nzchar(missing_labels))) {
          supp <- supp + rowAnnotation(
            ` `=anno_text(missing_labels,gp=gpar(fontsize=5)),show_annotation_name=FALSE)
        }
        leg_cluster <- Legend(title='Cluster mean',labels=names(cluster_colors),type='lines',
          background='white',legend_gp=gpar(col=unname(cluster_colors),lwd=0.65),
          grid_height=unit(2,'mm'),grid_width=unit(4,'mm'),
          title_gp=gpar(fontsize=6),labels_gp=gpar(fontsize=5))
        leg_frame <- Legend(title='Codon start position',labels=c('First','Second','Third'),
          legend_gp=gpar(fill=frame_colors,col='grey45',lwd=0.3),
          grid_height=unit(2,'mm'),grid_width=unit(2,'mm'),
          title_gp=gpar(fontsize=6),labels_gp=gpar(fontsize=5))
        height <- max(6, length(kmers)*0.17+2.8)
        save_plot <- function(obj, name, width, legends, plot_height=height) {
          for (ext in c('pdf','png','svg')) {
            path <- file.path(out,paste0(name,'.',ext))
            if (ext=='pdf') pdf(path,width=width,height=plot_height,onefile=FALSE)
            else if (ext=='svg') svg(path,width=width,height=plot_height,family='Arial')
            else png(path,width=width,height=plot_height,units='in',res=300,family='Arial')
            tryCatch({
              draw(obj, merge_legends=TRUE, heatmap_legend_side='right',
                   annotation_legend_side='right', align_heatmap_legend=if (name=='main_domains') 'heatmap_top' else 'heatmap_center',
                   align_annotation_legend=if (name=='main_domains') 'heatmap_top' else 'heatmap_center', annotation_legend_list=legends,
                   padding=unit(c(2,2,2,2),'mm'))
              if (name=='main_domains') {
                # anno_lines uses exactly [30,100] because extend=0; draw one threshold per slice.
                for (slice in seq_along(levels(droplevels(split)))) {
                  decorate_annotation('In Frame Index', {
                    x <- (100*threshold-30)/70
                    if (x>=0 && x<=1) grid.segments(x0=unit(x,'npc'),x1=unit(x,'npc'),
                      y0=unit(0,'npc'),y1=unit(1,'npc'),gp=gpar(col='#999999',lty=2,lwd=0.6))
                  }, slice=slice)
                }
              }
            }, finally=dev.off())
          }
        }
        save_plot(main,'main_domains',max(3.6, (ncol(d)*cell_mm + 75)/25.4),
                  list(main_leg_cell,main_leg_frame,main_leg_threshold),
                  plot_height=max(4.5,(nrow(d)*main_row_mm + main_gap_mm + 37)/25.4))
        save_plot(supp,'supplement_positions_frames',7.4,list(leg_cluster,leg_frame),
                  plot_height=max(4.5,(nrow(p)*supp_row_mm + main_gap_mm + 40)/25.4))
        writeLines(c(capture.output(sessionInfo()),
          paste('Threshold (inclusive):',threshold), paste('Start fraction:',start_fraction),
          paste('Colour scale and clustering cap: +/-',CAP),
          'Classification: Euclidean distance, Ward.D2, capped CDS log2 OR centred on each row mean; cutree(k=2).',
          'The cluster with the higher mean first-10%-minus-rest capped contrast is CDS start.',
          'Within-group ordering: decreasing In Frame Index (100 * maximum CDS frame fraction), then k-mer.',
          'Positional heatmap shows capped, uncentred log2 OR; domain values use the same colour scale.',
          'cluster_assignments.tsv: start_minus_body is capped; start_minus_body_raw is uncapped.',
          'Frame colours: First=lightgray, Second=white, Third=black; displayed k-mers use U.',
          'Main panel: CDS start above CDS body; RNA row labels, In Frame Index, threshold line, frame strip and domain heatmap.',
          'Main colors: purple-white-salmon; NS cells light grey; below-frame cells grey; stars and N.C. labels removed from main.',
          'Main labels: In Frame Index is 100 * maximum CDS frame fraction, not exclusively the first reading frame.',
          'Main row height: 2 mm; compact legends aligned at the top; short domain column labels.',
          'Supplement: capped uncentred positional heatmap; mean profiles of both clusters above; CDS and shuffle frame fractions and In Frame Index at right.',
          'Supplement: 2 mm rows; purple-white-salmon colours; cluster_mean_profiles.tsv contains the displayed mean profiles.',
          'Frame-threshold failure is not proof of noncoding.'),
          file.path(out,'plot_settings.txt'))
        message('Saved figures and assignments in ', normalizePath(out))
        """
            script_path = destination / 'heatmap_script.R'
            script_path.write_text(textwrap.dedent(r_code).lstrip(), encoding='utf-8')
            command = [r_executable, '--vanilla', str(script_path.resolve()),
                       str(destination.resolve()), str(figures_dir.resolve()),
                       str(frame_threshold), str(start_fraction), str(alpha), str(color_cap),
                       "all" if max_domains is None else str(max_domains)]
            completed = subprocess.run(command, capture_output=True, text=True)
            (destination / 'heatmap_R.log').write_text(completed.stdout + completed.stderr, encoding='utf-8')
            if completed.returncode:
                self.k_cds = result
                self.domain_enrichment = result
                raise RuntimeError('Final R panels failed; existing K-CDS outputs are saved in ' + str(destination) + '\n' + completed.stderr)
            names = ['main_domains', 'supplement_positions_frames']
            for name in names:
                paths = {ext: str(figures_dir / (name + '.' + ext)) for ext in ['pdf', 'png', 'svg']}
                if not all(Path(path).is_file() for path in paths.values()):
                    raise RuntimeError('R did not generate all final figure files: ' + name)
                result['final_figures'][name] = paths
            result['cluster_assignments'] = pd.read_csv(figures_dir / 'cluster_assignments.tsv', sep='\t')
            result['figures_dir'] = str(figures_dir)
            print(result['cluster_assignments']['cluster'].value_counts().to_string())
            if show_final:
                try:
                    from IPython import get_ipython
                    from IPython.display import display, Image
                    if get_ipython() is not None:
                        for name in names:
                            display(Image(filename=result['final_figures'][name]['png']))
                except ImportError:
                    pass
        (destination / 'settings.json').write_text(json.dumps(settings, indent=2), encoding='utf-8')
        self.k_cds = result
        self.domain_enrichment = result
        return result

    def cooccurrence(self, kmers=None, *, pairs=None, species='input',
                     fasta_file=None, n_resamples=1000, n_background=1000000000,
                     window_width=20, step=1, min_distance=8,
                     min_enrichment=2.0, pvalue_threshold=0.05,
                     min_support_fraction=0.10, fdr_threshold=None,
                     fdr_method='fdr_bh', masking=False,
                     tantan_executable=None, tantan_options=None,
                     seed=42, output_dir=None, plot=True, verbose=True,
                     fit_tail=True, threads=1, keep_windows=True, resume_dir=None, jobs=1):
        """WORK IN PROGRESS: enrichment of nonoverlapping motif spacings.

        jobs controls concurrent pair processes (Linux fork); threads controls
        simulation threads PER process, so jobs*threads is the requested maximum.
        Pair computations and FDR table compression both use process workers.
        Completed, hash-verified branches and pairs are reused on resume.

        v1.2 retains the empirical scan p-value unchanged and additionally fits
        a data-parameterized generalized Pareto upper tail to the global maxima.
        Inspect tail_fit.json diagnostics; fitted p-values are exploratory.
        Both empirical and fitted pair FDR are reported, with explicit empirical
        fallback for unavailable fits; only empirical p/FDR gates selection.
        threads parallelizes independent batches with unchanged per-batch seeds.
        keep_windows=False stores compressed per-pair tables without pooling
        all windows into RAM. resume_dir reuses hash-verified completed pairs.


        All unordered distinct motif pairs use every Cartesian forward-strand
        occurrence distance within each reference transcript. Distances >=
        min_distance (default 8, i.e. strictly >7 nt) enter sliding 20-nt windows.
        The expected distribution is computed exactly from valid coordinate-pair
        geometry, weighted by nA*nB BEFORE the distance filter, then conditioned
        globally on the same distance filter. n_background is only a multinomial
        sample of this distribution for display, not the inferential sample size.

        Each MC replicate jointly samples nA+nB distinct valid starts within each
        transcript, assigns A/B labels, rebuilds all Cartesian distances and
        normalizes AFTER filtering. Its maximum excess across ALL windows
        calibrates each observed window's scan-adjusted one-sided p-value.
        This is a conditional uniform-site model; local repeat clustering is
        NOT preserved. Real motif coordinates are allowed in randomizations.

        Apply enrichment >= min_enrichment, pvalue_scan < pvalue_threshold and
        support >= ceil(min_support_fraction * all_reference_transcripts) to
        EACH original window, THEN merge overlapping passing windows. A TX
        counts once per window/merged interval. Report the strongest constituent
        window p-value explicitly; no new region-level test is computed.
        fdr_threshold=None disables correction between pairs. Otherwise compute
        pair_FDR by fdr_method ('fdr_bh' or conservative 'fdr_by'), separately
        within each masking branch, using the global max-scan p-value of EVERY
        unordered pair of the declared motif set, before effect/support filtering.
        An incomplete pairs subset is rejected; at least two pairs are required.
        Untestable pairs count as p=1 internally, but retain NaN in exported FDR.
        pair_FDR is a PAIR-level gate, not an interval-level FDR. All local
        enrichment/pvalue/support criteria remain required. No FDR is inferred
        across the masked and unmasked branches.

        masking=True runs paired unmasked/masked branches on the SAME input
        motifs using tantan hard masking (-x N, standard DNA parameters unless
        explicitly changed). All masked windows are excluded from observed hits
        and random valid coordinates; lengths/IDs/denominator are unchanged.
        Temporary masked FASTA is removed. Inputs and other KEA methods are
        untouched. Return a dict of tables, plots, metadata and output directory.
        """
        from pathlib import Path
        import hashlib, json, math, tempfile, subprocess, shutil, time
        from itertools import combinations
        import numpy as np
        import pandas as pd
        from Bio import SeqIO
        from scipy.signal import correlate
        from numba import njit

        from concurrent.futures import ThreadPoolExecutor

        VERSION = '1.3.0-spacing-jobs-wip'
        if not isinstance(masking, (bool, np.bool_)):
            raise TypeError('masking must be True or False')
        if not isinstance(fit_tail, (bool, np.bool_)) or not isinstance(keep_windows, (bool, np.bool_)):
            raise TypeError('fit_tail and keep_windows must be boolean')
        if isinstance(threads,bool) or int(threads)!=threads or threads<1:
            raise ValueError('threads must be a positive integer')
        threads=int(threads)
        if isinstance(jobs,(bool,np.bool_)) or int(jobs)!=jobs or jobs<1:
            raise ValueError('jobs must be a positive integer')
        jobs=int(jobs)
        if jobs>1:
            import multiprocessing as multiprocessing_local
            if 'fork' not in multiprocessing_local.get_all_start_methods():
                raise ValueError('jobs>1 currently requires Linux/POSIX fork; use jobs=1 on this platform')
        for name, value in [('n_resamples', n_resamples), ('n_background', n_background),
                            ('window_width', window_width), ('step', step), ('min_distance', min_distance)]:
            if isinstance(value, bool) or int(value) != value or value <= 0:
                raise ValueError(name + ' must be a positive integer')
        n_resamples, n_background = int(n_resamples), int(n_background)
        window_width, step, min_distance = int(window_width), int(step), int(min_distance)
        if not (0 < pvalue_threshold < 1 and 0 <= min_support_fraction <= 1 and min_enrichment >= 1):
            raise ValueError('Invalid p-value, support or enrichment threshold')
        if fdr_threshold is not None and not (0 < fdr_threshold < 1):
            raise ValueError('fdr_threshold must be None or a number strictly between 0 and 1')
        if fdr_method not in {'fdr_bh','fdr_by'}:
            raise ValueError("fdr_method must be 'fdr_bh' or 'fdr_by'")
        if not isinstance(seed, (int, np.integer)) or seed < 0:
            raise ValueError('seed must be a nonnegative integer')
        source = Path(fasta_file or self.reference_fasta).resolve()
        if not source.is_file():
            raise FileNotFoundError(source)
        source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
        if kmers is None:
            results = (getattr(self, 'KEA_results', None) or {}).get(species, {})
            methods = ['delta_median', 'stat_log2fc']
            if any(m not in results for m in methods):
                raise ValueError('Pass the KEA-enriched kmers explicitly, or load delta_median and stat_log2fc results first.')
            kmers = set.intersection(*(set(results[m]['enriched']) for m in methods))
        words = sorted({str(w).upper().replace('U', 'T') for w in kmers})
        if len(words) < 2 or any(not w or set(w) - set('ACGT') for w in words):
            raise ValueError('At least two distinct ACGT k-mers are required')
        lengths = {len(w) for w in words}
        if len(lengths) != 1:
            raise ValueError('All k-mers must have the same length')
        k = lengths.pop()
        if min_distance < k:
            raise ValueError('min_distance must be >= k to exclude overlapping motif windows')
        selected_pairs = list(combinations(words, 2)) if pairs is None else sorted({tuple(sorted(str(w).upper().replace('U','T') for w in pair)) for pair in pairs})
        if not selected_pairs or any(len(pair) != 2 or pair[0] == pair[1] or not set(pair).issubset(words) for pair in selected_pairs):
            raise ValueError('pairs must contain distinct unordered pairs from kmers')
        declared_pair_family = list(combinations(words,2))
        if fdr_threshold is not None:
            if len(declared_pair_family) < 2:
                raise ValueError('FDR requires multiple declared motif pairs; use fdr_threshold=None for a single-pair toy')
            if set(selected_pairs) != set(declared_pair_family):
                raise ValueError('FDR requires ALL unordered pairs of the declared kmers; remove the pairs subset or disable FDR')
        records = list(SeqIO.parse(str(source), 'fasta'))
        ids = [r.id for r in records]
        if not ids or len(set(ids)) != len(ids):
            raise ValueError('Reference FASTA must have unique IDs and at least one transcript')
        sequences = [str(r.seq).upper().replace('U','T') for r in records]
        if any(set(s)-set('ACGTNRYWSKMBDHVX-') for s in sequences):
            raise ValueError('Unsupported non-DNA symbols in reference')
        n_tx = len(ids)
        support_required = int(math.ceil(min_support_fraction * n_tx))
        maxd = max(map(len, sequences))-k
        if maxd-min_distance+1 < window_width:
            raise ValueError('Reference too short for the requested distance windows')
        starts = np.arange(min_distance, maxd-window_width+2, step, dtype=np.int64)
        distance = np.arange(maxd+1, dtype=np.int64)
        base = Path(output_dir).resolve() if output_dir is not None else Path(self.dir_out).resolve()/'Cooccurrence'/species/'spacing_v1'
        base.mkdir(parents=True, exist_ok=True)
        if resume_dir is None:
            destination = Path(tempfile.mkdtemp(prefix='run_', dir=str(base)))
        else:
            destination = Path(resume_dir).resolve()
            if not destination.is_dir() or not destination.is_relative_to(base):
                raise ValueError('resume_dir must be an existing run within output_dir')
        started = time.time()
        settings = dict(version=VERSION, state='running', reference_fasta=str(source), reference_sha256=source_sha,
                        kmers=words, pairs=[list(p) for p in selected_pairs], n_transcripts=n_tx, k=k,
                        n_resamples=n_resamples, n_background=n_background, window_width=window_width, step=step,
                        min_distance=min_distance, min_enrichment=float(min_enrichment), pvalue_threshold=float(pvalue_threshold),
                        min_support_fraction=float(min_support_fraction), required_support_transcripts=support_required,
                        masking=bool(masking), seed=int(seed), fdr_computed=False,
                        fit_tail=bool(fit_tail), threads=threads, keep_windows=bool(keep_windows), jobs=jobs,
                        primary_pvalue='empirical scan-adjusted pvalue_scan; fitting never changes primary selection',
                        fitted_pvalue='GPD upper tail of the SAME 1000 global maxima; empirical body; diagnostic failures or finite-endpoint extrapolation are NaN',
                        fitted_fdr='secondary pair_FDR_fitted uses fitted pair p-values where available and empirical fallback otherwise; never gates primary selection',
                        fdr_threshold=None if fdr_threshold is None else float(fdr_threshold),
                        fdr_method=fdr_method, fdr_level='pair; separately within each masking branch',
                        fdr_pair_family_size=len(declared_pair_family),
                        fdr_family='all unordered pairs of declared input motifs, before enrichment/support filtering; untestable use p=1 internally and exported NaN',
                        weighting='nA*nB per TX before distance filtering; normalize globally after filtering',
                        null='joint uniform distinct valid starts within each TX; counts preserved; local clustering not preserved; real A/B coordinates allowed',
                        pvalue='one-sided excess vs fixed exact expected; maximum across ALL windows in each complete MC replicate; plus-one correction',
                        merge='overlapping passing windows merged only AFTER all three per-window filters; bookended windows remain separate',
                        region_pvalue='best_window_pvalue_scan is a constituent-window p-value, not a test of the merged interval',
                        masked_fasta_policy='temporary hard-masked FASTA only; original input never overwritten')
        def write_json(path, value):
            temporary=path.with_name(path.name+'.writing')
            temporary.write_text(json.dumps(value,indent=2),encoding='utf-8')
            temporary.replace(path)
        import inspect, sys, scipy, numba
        settings['software_versions'] = dict(python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__, numba=numba.__version__)
        settings['method_source_sha256'] = hashlib.sha256(inspect.getsource(type(self).cooccurrence).encode()).hexdigest()
        if resume_dir is not None:
            previous=json.loads((destination/'settings.json').read_text())
            ignored={'state','threads','jobs','software_versions','fdr_computed','version','method_source_sha256'}
            # v1.2 hash pinned after AST/RNG equivalence tests.
            compatible_previous={'2d0041c6f8b33f20f2653bd58528415ca818226b26ceba204540d9ebeb0dc84e'}
            if previous.get('method_source_sha256') not in compatible_previous|{settings['method_source_sha256']}:
                raise ValueError('Resume method source is not a verified compatible version')
            comparable=[key for key in settings if key not in ignored]
            if any(previous.get(key)!=settings[key] for key in comparable):
                raise ValueError('Resume inputs, parameters or method source differ; existing run preserved')
            settings['resume_history']=previous.get('resume_history',[])+[
                dict(previous_version=previous.get('version'),previous_method_source_sha256=previous.get('method_source_sha256'),
                     current_version=VERSION,jobs=jobs,threads=threads,resumed_unix=time.time())]
        write_json(destination/'settings.json', settings)
        window_columns = ['kmer_a','kmer_b','start_nt','end_nt','observed_count','observed_pct','expected_pct',
                          'excess_percentage_points','enrichment','support_transcripts','support_fraction',
                          'pvalue_scan','pvalue_scan_fitted','passes_enrichment','passes_pvalue','passes_support','selected']
        region_columns = ['region','kmer_a','kmer_b','start_nt','end_nt','width_nt','n_windows',
                          'observed_pct','expected_pct','enrichment','support_transcripts','support_fraction',
                          'best_window_start_nt','best_window_end_nt','best_window_pvalue_scan','best_window_pvalue_scan_fitted']
        pair_rows = []
        branch_results = {}

        def fit_maximum_tail(maxima, scores):
            """Semiparametric upper-tail fit; empirical body, no zero extrapolations.

            GPD shape/scale are fitted from the data. Three tail thresholds
            are screened; choose the lowest acceptable threshold (most data).
            The KS discrepancy is a diagnostic, NOT a calibrated GOF p-value.
            Failure/finite-endpoint extrapolation yields NaN, never p=0.
            """
            from scipy.stats import genpareto
            import warnings
            scores = np.asarray(scores, dtype=float)
            fitted = np.full(scores.shape, np.nan)
            meta = dict(status='disabled' if not fit_tail else 'insufficient_simulations',
                        family='generalized Pareto upper tail; empirical body',
                        threshold_selection='lowest quantile passing descriptive KS screen, with >=50 exceedances and >=20 distinct values',
                        ks_screen='D <= 1.36/sqrt(n_tail); heuristic diagnostic only, NOT a GOF p-value',
                        candidates=[], n_simulations=len(maxima))
            if not fit_tail or len(maxima) < 250:
                return fitted, meta
            empirical = (1+len(maxima)-np.searchsorted(np.sort(maxima),scores-1e-14,side='left'))/(len(maxima)+1)
            chosen = None
            for quantile in (0.80,0.90,0.95):
                threshold = float(np.quantile(maxima,quantile))
                tail = np.asarray(maxima[maxima>threshold]-threshold,dtype=float)
                candidate = dict(quantile=quantile,threshold=threshold,n_tail=len(tail),accepted=False)
                meta['candidates'].append(candidate)
                if len(tail)<50 or len(np.unique(tail))<20:
                    candidate['reason']='too few exceedances or distinct values'
                    continue
                unit = float(tail.max())
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter('error',RuntimeWarning)
                        shape,_,scale = genpareto.fit(tail/unit,floc=0)
                    if not np.isfinite([shape,scale]).all() or scale<=0:
                        raise ValueError('nonfinite parameters')
                    cdf = genpareto.cdf(np.sort(tail)/unit,shape,loc=0,scale=scale)
                    ranks = np.arange(1,len(tail)+1)/len(tail)
                    discrepancy = float(max(np.max(ranks-cdf),np.max(cdf-(ranks-1/len(tail)))))
                    limit = float(1.36/np.sqrt(len(tail)))
                    accepted = bool(np.isfinite(discrepancy) and discrepancy<=limit)
                    candidate.update(shape=float(shape),scale=float(scale*unit),
                                     ks_distance=discrepancy,ks_screen_limit=limit,accepted=accepted)
                    if accepted and chosen is None:
                        chosen=(threshold,len(tail),float(shape),float(scale*unit),quantile,discrepancy)
                except Exception as error:
                    candidate['reason']=f'{type(error).__name__}: {error}'
            if chosen is None:
                meta['status']='no_acceptable_tail_fit'
                return fitted,meta
            threshold,n_tail,shape,scale,quantile,discrepancy=chosen
            fitted[:]=empirical
            use_tail=scores>threshold
            logp=np.log(n_tail/len(maxima))+genpareto.logsf(scores[use_tail]-threshold,shape,loc=0,scale=scale)
            finite=np.isfinite(logp) & (logp>=np.log(np.finfo(float).tiny))
            tail_p=np.full(logp.shape,np.nan)
            tail_p[finite]=np.exp(logp[finite])
            fitted[use_tail]=tail_p
            meta.update(status='ok',threshold=threshold,threshold_quantile=quantile,
                        n_tail=n_tail,shape=shape,scale=scale,ks_distance=discrepancy,
                        tail_probability=n_tail/len(maxima),maximum_simulated=float(np.max(maxima)),
                        finite_upper_endpoint=float(threshold-scale/shape) if shape<0 else None,
                        windows_in_tail=int(use_tail.sum()),
                        windows_beyond_simulated_maximum=int((scores>np.max(maxima)).sum()),
                        windows_fit_unavailable=int(np.isnan(fitted).sum()),
                        unavailable_policy='NaN for unsupported or underflowing tail estimates; empirical p-value retained',
                        interpretation='exploratory model estimate; tail screening does not establish p-value calibration')
            return fitted,meta

        def window_path(pair_dir):
            return pair_dir/('windows.tsv' if keep_windows else 'windows.tsv.gz')

        def pair_checkpoint(pair_dir):
            names=[window_path(pair_dir).name,'distributions.npz','transcript_counts.tsv',
                   'enriched_regions.tsv','selected_regions.tsv','support_transcripts.tsv',
                   'tail_fit.json','summary.json']
            info={name:dict(bytes=(pair_dir/name).stat().st_size,
                           sha256=hashlib.sha256((pair_dir/name).read_bytes()).hexdigest()) for name in names}
            write_json(pair_dir/'pair_complete.json',dict(outputs=info))

        def completed_pair(pair_dir):
            try:
                checkpoint=json.loads((pair_dir/'pair_complete.json').read_text())
                return bool(checkpoint['outputs']) and all(
                    (pair_dir/name).is_file() and (pair_dir/name).stat().st_size==info['bytes']
                    and hashlib.sha256((pair_dir/name).read_bytes()).hexdigest()==info['sha256']
                    for name,info in checkpoint['outputs'].items())
            except (OSError,ValueError,KeyError):
                return False

        def parallel_results(action, items):
            """Linux fork workers inherit read-only sequence/mapping caches.

            Only pair identifiers and small summaries cross queues. Workers
            never update global summaries/FDR. Seeds remain pair/batch based.
            """
            if jobs == 1:
                for item in items:
                    yield item,action(item)
                return
            import multiprocessing as multiprocessing_local
            import queue, traceback
            context=multiprocessing_local.get_context('fork')
            count=min(jobs,len(items))
            tasks=context.Queue(maxsize=max(2,count*2))
            results=context.Queue(maxsize=max(2,count*2))
            def worker():
                while True:
                    item=tasks.get()
                    if item is None:
                        return
                    try:
                        answer=action(item)
                        results.put((True,item,answer))
                    except BaseException:
                        results.put((False,item,traceback.format_exc()))
                        return
            workers=[context.Process(target=worker) for _ in range(count)]
            try:
                for process in workers:
                    process.start()
                pending=iter(items)
                for _ in workers:
                    tasks.put(next(pending))
                for _ in range(len(items)):
                    while True:
                        try:
                            success,item,answer=results.get(timeout=1)
                            break
                        except queue.Empty:
                            dead=[p for p in workers if p.exitcode is not None]
                            if dead:
                                raise RuntimeError('Cooccurrence worker exited unexpectedly: '+str([(p.pid,p.exitcode) for p in dead]))
                    if not success:
                        raise RuntimeError('Cooccurrence worker failed for '+str(item)+'\n'+answer)
                    yield item,answer
                    try:
                        tasks.put(next(pending))
                    except StopIteration:
                        pass
                for _ in workers:
                    tasks.put(None)
                for process in workers:
                    process.join(timeout=10)
            finally:
                for process in workers:
                    if process.is_alive():
                        process.terminate()
                for process in workers:
                    process.join(timeout=5)
                    if process.is_alive():
                        process.kill();process.join(timeout=5)
                tasks.cancel_join_thread();results.cancel_join_thread()
                tasks.close();results.close()

        def find_all(s, word):
            hits = []
            p = s.find(word)
            while p >= 0:
                hits.append(p)
                p = s.find(word, p+1)
            return np.asarray(hits, dtype=np.int64)

        def merge_intervals(windows, mask, a, b, unique_distances):
            groups = []
            for j in np.flatnonzero(mask):
                lo, stop = int(starts[j]), int(starts[j]+window_width)
                if groups and lo < groups[-1]['stop']:
                    groups[-1]['stop'] = max(groups[-1]['stop'], stop)
                    groups[-1]['members'].append(int(j))
                else:
                    groups.append(dict(start=lo, stop=stop, members=[int(j)]))
            merged, membership = [], []
            for i, group in enumerate(groups, 1):
                lo, stop = group['start'], group['stop']
                members = windows.iloc[group['members']]
                best = members.sort_values(['excess_percentage_points','start_nt'], ascending=[False,True]).iloc[0]
                support_ids = []
                for t, ds in unique_distances.items():
                    j = np.searchsorted(ds, lo)
                    if j < len(ds) and ds[j] < stop:
                        support_ids.append(ids[t])
                probability = float(observed[lo:stop].sum()/observed.sum())
                exp = float(expected[lo:stop].sum())
                name = f'R{i:03d}'
                merged.append(dict(region=name,kmer_a=a,kmer_b=b,start_nt=lo,end_nt=stop-1,width_nt=stop-lo,
                                   n_windows=len(members),observed_pct=100*probability,expected_pct=100*exp,
                                   enrichment=probability/exp if exp else np.nan,support_transcripts=len(support_ids),
                                   support_fraction=len(support_ids)/n_tx,best_window_start_nt=int(best.start_nt),
                                   best_window_end_nt=int(best.end_nt),best_window_pvalue_scan=float(best.pvalue_scan),
                                   best_window_pvalue_scan_fitted=float(best.pvalue_scan_fitted)))
                membership.extend(dict(region=name,kmer_a=a,kmer_b=b,transcript_id=tx) for tx in support_ids)
            return pd.DataFrame(merged, columns=region_columns), pd.DataFrame(membership, columns=['region','kmer_a','kmer_b','transcript_id'])

        @njit(nogil=True)
        def draw_batch(flat, offsets, na, nb, n_dist, first, width, stride, expected_mass, B, rng_seed):
            np.random.seed(rng_seed)
            max_sites = 0
            max_draw = 0
            for t in range(len(na)):
                max_sites = max(max_sites, offsets[t+1]-offsets[t])
                max_draw = max(max_draw, na[t]+nb[t])
            pool = np.arange(max_sites, dtype=np.int64)
            swaps = np.empty(max_draw, dtype=np.int64)
            coordinates = np.empty(max_draw, dtype=np.int64)
            maxima = np.empty(B, dtype=np.float64)
            totals = np.empty(B, dtype=np.int64)
            for trial in range(B):
                hist = np.zeros(n_dist, dtype=np.int64)
                for t in range(len(na)):
                    off = offsets[t]
                    m = offsets[t+1]-off
                    ns = na[t]+nb[t]
                    for i in range(ns):
                        j = np.random.randint(i,m)
                        swaps[i] = j
                        temp = pool[i]; pool[i] = pool[j]; pool[j] = temp
                        coordinates[i] = flat[off+pool[i]]
                    for i in range(na[t]):
                        for j in range(na[t],ns):
                            delta = abs(coordinates[i]-coordinates[j])
                            if delta >= first:
                                hist[delta] += 1
                    for i in range(ns-1,-1,-1):
                        j = swaps[i]
                        temp = pool[i]; pool[i] = pool[j]; pool[j] = temp
                total = hist.sum()
                if total == 0:
                    raise ValueError('A random replica has no retained distances; change the model or input instead of silently dropping it.')
                totals[trial] = total
                prefix = np.empty(n_dist+1, dtype=np.int64)
                prefix[0] = 0
                for i in range(n_dist):
                    prefix[i+1] = prefix[i]+hist[i]
                maximum = -1.0
                for q in range(len(expected_mass)):
                    start = first+q*stride
                    score = (prefix[start+width]-prefix[start])/total-expected_mass[q]
                    maximum = max(maximum, score)
                maxima[trial] = maximum
            for i in range(max_sites):
                if pool[i] != i:
                    raise ValueError('Sampling pool restoration failed')
            return maxima, totals

        def make_plots(pair_dir, a, b, windows, candidates, selected_regions, summary):
            import matplotlib.pyplot as plt
            from matplotlib.patches import Rectangle
            paths = {}
            centers = starts+(window_width-1)/2
            for phase in [1,2,3]:
                fig = plt.figure(figsize=(10,5.8))
                grid = fig.add_gridspec(2,1,height_ratios=[8,.7],hspace=.08)
                ax = fig.add_subplot(grid[0]); band = fig.add_subplot(grid[1],sharex=ax)
                ax.set_xscale('log'); band.set_xscale('log')
                if len(windows):
                    ax.fill_between(centers,0,windows.expected_pct,color='#BDBDBD',alpha=.65,label='Expected: matched TX weights')
                    ax.plot(centers,windows.expected_pct,color='#888888',lw=.9)
                    ax.plot(centers,windows.observed_pct,color='#D62728',lw=.9,label='Observed A/B distances')
                    if phase == 3:
                        ok = windows.selected.to_numpy(bool)
                        ax.scatter(centers[ok],windows.loc[ok,'observed_pct'],s=3,color='#1565C0',linewidths=0,label='All filters passed',zorder=4)
                    regions = candidates if phase == 2 else selected_regions if phase == 3 else pd.DataFrame()
                    for r in regions.itertuples():
                        band.add_patch(Rectangle((r.start_nt,.12),r.width_nt,.76,fill=False,edgecolor='#1565C0',lw=1.1,linestyle=(0,(4,3))))
                    if phase == 3 and selected_regions.empty:
                        ax.text(.98,.69,'No window passes all three filters',transform=ax.transAxes,ha='right',color='#1565C0',fontsize=10)
                    ax.set_ylim(bottom=0)
                else:
                    ax.text(.5,.5,'No observed distances after filtering\nNo p-value can be computed',transform=ax.transAxes,ha='center',va='center')
                    ax.set_ylim(0,1)
                titles = {1:'Phase 1: observed vs exact expected',2:f'Phase 2: enrichment >= {min_enrichment:g}x; overlapping windows merged',3:f'Phase 3: enrichment + scan p-value + TX support; {len(selected_regions)} merged intervals'}
                ax.set(title=f'{a} / {b} | {summary["analysis"]}\n{titles[phase]}',ylabel=f'Distances in {window_width}-nt window (%)',xlim=(min_distance,maxd))
                ax.title.set_fontsize(11)
                ax.tick_params(axis='x',labelbottom=False)
                ax.spines[['top','right']].set_visible(False)
                if len(windows): ax.legend(loc='upper right',frameon=False,fontsize=9)
                band.set(ylim=(0,1),yticks=[],xlabel='Start-start distance (nt, log scale; window centers)')
                band.set_ylabel('Merge',rotation=0,labelpad=25,fontsize=9)
                band.spines[['top','right','left']].set_visible(False)
                fig.subplots_adjust(left=.11,right=.98,top=.84,bottom=.19)
                foot = f'{window_width}-nt windows, step {step}; distance >= {min_distance} nt. Same resolution for curves and tests.'
                if phase == 3:
                    foot += f'\n{n_resamples} replicates; p < {pvalue_threshold:g}; support >= {support_required}/{n_tx} TX.'
                    if fdr_threshold is None:
                        foot += ' No FDR.'
                    else:
                        q = summary.get('pair_FDR', None)
                        foot += f' Pair FDR < {fdr_threshold:g}; ' + (f'q={q:.4g}' if q is not None else 'not testable')
                fig.text(.11,.045,foot,fontsize=8)
                for ext in ['png','pdf']:
                    path = pair_dir/f'phase_{phase}.{ext}'
                    fig.savefig(path,dpi=180)
                    paths[f'phase_{phase}_{ext}'] = str(path)
                plt.close(fig)
            return paths

        observed=expected=None
        masking_stats = []
        try:
            branches = [('unmasked',sequences)]
            if masking:
                exe = str(tantan_executable or shutil.which('tantan') or '')
                if not exe:
                    local = Path(self.dir_out)/'_runtime/tantan51/bin/tantan'
                    if local.is_file(): exe = str(local)
                found = shutil.which(exe) if exe else None
                if not found:
                    raise FileNotFoundError('tantan executable not found; pass tantan_executable explicitly')
                exe = str(Path(found).resolve())
                options = list(tantan_options or [])
                if any(not isinstance(x,str) for x in options) or any(x in {'-x','-f','-p'} for x in options):
                    raise ValueError('tantan_options must be DNA model options only; hard masking/output mode are managed internally')
                version = subprocess.run([exe,'--version'],capture_output=True,text=True,check=True)
                settings['tantan'] = dict(executable=exe,version=(version.stdout+version.stderr).strip(),
                                          executable_sha256=hashlib.sha256(Path(exe).read_bytes()).hexdigest(),options=options,
                                          command=[exe,'-x','N']+options+[str(source)])
                with tempfile.TemporaryDirectory(prefix='kea_spacing_tantan_') as tmp:
                    masked_file = Path(tmp)/'reference.masked.fa'
                    with masked_file.open('w') as stream:
                        done = subprocess.run(settings['tantan']['command'],stdout=stream,stderr=subprocess.PIPE,text=True,check=True)
                    (destination/'tantan.log').write_text(done.stderr)
                    masked_records = list(SeqIO.parse(str(masked_file),'fasta'))
                    if [r.id for r in masked_records] != ids:
                        raise ValueError('tantan changed FASTA IDs/order')
                    masked_sequences = [str(r.seq).upper().replace('U','T') for r in masked_records]
                    mask_intervals = []
                    for tx,s,m in zip(ids,sequences,masked_sequences):
                        if len(s) != len(m): raise ValueError('tantan changed sequence lengths')
                        before = np.frombuffer(s.encode(),np.uint8)
                        after = np.frombuffer(m.encode(),np.uint8)
                        changed = before != after
                        if np.any(changed & (after != ord('N'))): raise ValueError('Unexpected tantan base change')
                        padded = np.r_[False,changed,False].astype(np.int8)
                        edges = np.diff(padded)
                        mask_intervals.extend((tx,int(a),int(b)) for a,b in zip(np.flatnonzero(edges==1),np.flatnonzero(edges==-1)))
                        masking_stats.append(dict(transcript_id=tx,length_nt=len(s),masked_nt=int(changed.sum()),masked_fraction=float(changed.mean()) if len(s) else 0.0))
                    pd.DataFrame(mask_intervals,columns=['transcript_id','start','end']).to_csv(destination/'tantan_masks.bed',sep='\t',header=False,index=False)
                    pd.DataFrame(masking_stats).to_csv(destination/'masking_stats.tsv',sep='\t',index=False)
                    settings['tantan']['masked_fasta_sha256'] = hashlib.sha256(masked_file.read_bytes()).hexdigest()
                branches.append(('masked',masked_sequences))
                write_json(destination/'settings.json',settings)
            for label, seqs in branches:
                branch_dir = destination/label
                branch_dir.mkdir(exist_ok=True)
                # A previously finalized branch is verified and reused wholesale;
                # never rewrite its windows just because worker count changed.
                cached=None
                if resume_dir is not None and (destination/'pair_summary.tsv').is_file():
                    saved_table=pd.read_csv(destination/'pair_summary.tsv',sep='\t',low_memory=False)
                    candidate=saved_table.loc[saved_table['analysis'].eq(label)]
                    eligible=len(candidate)==len(selected_pairs) and (branch_dir/'selected_regions.tsv').is_file() and (branch_dir/'support_transcripts.tsv').is_file()
                    if eligible and fdr_threshold is not None:
                        eligible='fdr_family_size' in candidate and candidate['fdr_family_size'].eq(len(declared_pair_family)).all()
                    if eligible:
                        def verify_saved(pair):
                            p=branch_dir/f'{pair[0]}__{pair[1]}'
                            return json.loads((p/'summary.json').read_text()) if completed_pair(p) else None
                        verified={}
                        for pair,summary in parallel_results(verify_saved,selected_pairs):
                            if summary is not None:verified[tuple(pair)]=summary
                        if len(verified)==len(selected_pairs):
                            cached=[verified[tuple(pair)] for pair in selected_pairs]
                if cached is not None:
                    pair_rows.extend(cached)
                    saved_windows=[pd.read_csv(window_path(branch_dir/f'{a}__{b}'),sep='\t') for a,b in selected_pairs] if keep_windows else []
                    branch_results[label]=dict(windows=pd.concat(saved_windows,ignore_index=True) if saved_windows else pd.DataFrame(columns=window_columns),
                                              windows_in_memory=bool(keep_windows),window_files=[str(window_path(branch_dir/f'{a}__{b}')) for a,b in selected_pairs],
                                              regions=pd.read_csv(branch_dir/'selected_regions.tsv',sep='\t'),
                                              support_transcripts=pd.read_csv(branch_dir/'support_transcripts.tsv',sep='\t'),
                                              plots={f"{r['kmer_a']}|{r['kmer_b']}":r.get('plots',{}) for r in cached},output_dir=str(branch_dir))
                    if verbose:print(f'[KEA spacing] {label}: reused {len(cached)} verified pairs and completed FDR',flush=True)
                    continue
                occurrences = {word:[find_all(s,word) for s in seqs] for word in words}
                valid_sites = []
                for s in seqs:
                    if len(s)<k:
                        valid_sites.append(np.empty(0,dtype=np.int32));continue
                    invalid = ~np.isin(np.frombuffer(s.encode(),np.uint8),np.frombuffer(b'ACGT',np.uint8))
                    cs = np.r_[0,np.cumsum(invalid,dtype=np.int64)]
                    valid_sites.append(np.flatnonzero((cs[k:]-cs[:-k])==0).astype(np.int32))
                geometry_cache = {}
                branch_windows = []; branch_regions = []; branch_support = []; plots = {}
                pair_indices={tuple(pair):j for j,pair in enumerate(selected_pairs)}
                def compute_pair(pair):
                    nonlocal observed,expected
                    a,b=pair
                    pair_number=pair_indices[tuple(pair)]+1
                    if verbose: print(f'[KEA spacing] {label} {pair_number}/{len(selected_pairs)}: {a}/{b}',flush=True)
                    pair_dir=branch_dir/f'{a}__{b}'
                    pair_dir.mkdir(exist_ok=True)
                    if resume_dir is not None and completed_pair(pair_dir):
                        return json.loads((pair_dir/'summary.json').read_text())
                    pair_started = time.time()
                    raw_pairs = sum(len(x)*len(y) for x,y in zip(occurrences[a],occurrences[b]))
                    observed = np.zeros(maxd+1,dtype=np.int64)
                    unique_distances = {}
                    support = np.zeros(len(starts),dtype=np.int64)
                    expected = np.zeros(maxd+1)
                    active = []; flat = []; offsets = [0]; na = []; nb = []
                    tx_rows = []
                    for t,(aa,bb) in enumerate(zip(occurrences[a],occurrences[b])):
                        tx_rows.append(dict(transcript_id=ids[t],length_nt=len(seqs[t]),n_a=len(aa),n_b=len(bb),raw_pairs=len(aa)*len(bb),valid_random_starts=len(valid_sites[t])))
                        if not len(aa) or not len(bb): continue
                        hist = np.zeros(maxd+1,dtype=np.int64)
                        for pos in range(0,len(aa),max(1,1000000//len(bb))):
                            ds = np.abs(aa[pos:pos+max(1,1000000//len(bb)),None]-bb[None,:]).ravel()
                            hist += np.bincount(ds,minlength=maxd+1)
                        if int(hist.sum()) != len(aa)*len(bb): raise RuntimeError('Observed Cartesian pair count mismatch')
                        hist[:min_distance] = 0
                        observed += hist
                        unique_distances[t] = np.flatnonzero(hist)
                        cumulative = np.r_[0,np.cumsum(hist)]
                        support += (cumulative[starts+window_width]-cumulative[starts])>0
                        sites = valid_sites[t]
                        if len(sites)<len(aa)+len(bb): raise RuntimeError('Not enough distinct valid random coordinates')
                        if t not in geometry_cache:
                            n = len(seqs[t])-k+1
                            marks = np.zeros(n)
                            marks[sites] = 1
                            corr = np.rint(correlate(marks,marks,mode='full',method='fft')).astype(np.int64)
                            geom = corr[n:]
                            total = len(sites)*(len(sites)-1)//2
                            if np.any(geom<0) or int(geom.sum())!=total: raise RuntimeError('Coordinate-pair geometry mismatch')
                            geometry_cache[t] = (geom.astype(np.int32),total)
                        geom,total = geometry_cache[t]
                        expected[1:len(geom)+1] += (len(aa)*len(bb)/raw_pairs)*(geom/total)
                        active.append(t);flat.append(sites);offsets.append(offsets[-1]+len(sites));na.append(len(aa));nb.append(len(bb))
                    pd.DataFrame(tx_rows).to_csv(pair_dir/'transcript_counts.tsv',sep='\t',index=False)
                    status = 'complete' if observed.sum() else 'no_observed_distances'
                    summary = dict(analysis=label,kmer_a=a,kmer_b=b,status=status,n_transcripts=n_tx,
                                   n_transcripts_with_both=len(active),raw_observed_pairs=int(raw_pairs),
                                   retained_observed_distances=int(observed.sum()),required_support_transcripts=support_required)
                    windows = pd.DataFrame(columns=window_columns)
                    candidates = pd.DataFrame(columns=region_columns)
                    merged = pd.DataFrame(columns=region_columns)
                    membership = pd.DataFrame(columns=['region','kmer_a','kmer_b','transcript_id'])
                    maxima = np.empty(0); background = np.zeros(maxd+1,dtype=np.int64)
                    expected_mass = np.zeros(len(starts)); retained = np.empty(0,dtype=np.int64)
                    fitted_pvalue=np.full(len(starts),np.nan)
                    tail_fit=dict(status='no_observed_distances')
                    if observed.sum():
                        expected[:min_distance] = 0
                        if expected.sum() <= 0: raise RuntimeError('No valid random distances after filtering')
                        expected /= expected.sum()
                        pair_seed = int.from_bytes(hashlib.sha256(f'{seed}|{label}|{a}|{b}'.encode()).digest()[:4],'little')
                        rng = np.random.default_rng(pair_seed)
                        background = rng.multinomial(n_background,expected)
                        if int(background.sum()) != n_background: raise RuntimeError('Background sample total mismatch')
                        ce = np.r_[0,np.cumsum(expected)]
                        co = np.r_[0,np.cumsum(observed)]
                        expected_mass = ce[starts+window_width]-ce[starts]
                        observed_count = co[starts+window_width]-co[starts]
                        observed_mass = observed_count/observed.sum()
                        excess = observed_mass-expected_mass
                        maxima_parts=[];retained_parts=[]
                        flat_array=np.concatenate(flat)
                        offset_array=np.asarray(offsets,np.int64)
                        na_array=np.asarray(na,np.int64);nb_array=np.asarray(nb,np.int64)
                        def simulate_batch(first):
                            return draw_batch(flat_array,offset_array,na_array,nb_array,maxd+1,min_distance,
                                              window_width,step,expected_mass,min(100,n_resamples-first),
                                              (pair_seed+first+1)%(2**32))
                        batch_starts=list(range(0,n_resamples,100))
                        with ThreadPoolExecutor(max_workers=threads) as executor:
                            for first,(mx,totals) in zip(batch_starts,executor.map(simulate_batch,batch_starts)):
                                maxima_parts.append(mx);retained_parts.append(totals)
                                if verbose: print(f'  Monte Carlo {min(first+100,n_resamples)}/{n_resamples}',flush=True)
                        maxima=np.concatenate(maxima_parts);retained=np.concatenate(retained_parts)
                        pvalue=(1+n_resamples-np.searchsorted(np.sort(maxima),excess-1e-14,side='left'))/(n_resamples+1)
                        fitted_pvalue,tail_fit=fit_maximum_tail(maxima,excess)
                        enrichment=np.divide(observed_mass,expected_mass,out=np.full_like(observed_mass,np.nan),where=expected_mass>0)
                        pass_enrichment=(enrichment>=min_enrichment)&(excess>0)
                        pass_p=pvalue<pvalue_threshold
                        pass_support=support>=support_required
                        selected=pass_enrichment&pass_p&pass_support
                        windows=pd.DataFrame(dict(kmer_a=a,kmer_b=b,start_nt=starts,end_nt=starts+window_width-1,
                                                  observed_count=observed_count,observed_pct=100*observed_mass,expected_pct=100*expected_mass,
                                                  excess_percentage_points=100*excess,enrichment=enrichment,
                                                  support_transcripts=support,support_fraction=support/n_tx,pvalue_scan=pvalue,pvalue_scan_fitted=fitted_pvalue,
                                                  passes_enrichment=pass_enrichment,passes_pvalue=pass_p,passes_support=pass_support,selected=selected),columns=window_columns)
                        candidates,_=merge_intervals(windows,pass_enrichment,a,b,unique_distances)
                        merged,membership=merge_intervals(windows,selected,a,b,unique_distances)
                        summary.update(windows_tested=len(starts),windows_enrichment=int(pass_enrichment.sum()),windows_pvalue=int(pass_p.sum()),
                                       windows_support=int(pass_support.sum()),windows_selected=int(selected.sum()),merged_intervals=len(merged),
                                       max_window_support_transcripts=int(support.max()),max_window_support_fraction=float(support.max()/n_tx),
                                       minimum_scan_pvalue=float(pvalue.min()),pair_scan_pvalue=float((1+(maxima>=excess.max()-1e-14).sum())/(n_resamples+1)),
                                       pair_scan_pvalue_fitted=float(fitted_pvalue[np.argmax(excess)]) if np.isfinite(fitted_pvalue[np.argmax(excess)]) else None)
                    else:
                        summary.update(windows_tested=0,windows_enrichment=0,windows_pvalue=0,windows_support=0,windows_selected=0,
                                       merged_intervals=0,max_window_support_transcripts=0,max_window_support_fraction=0.,minimum_scan_pvalue=None,pair_scan_pvalue=None,pair_scan_pvalue_fitted=None)
                    summary.update(tail_fit_status=tail_fit['status'],tail_fit_threshold=tail_fit.get('threshold'),
                                   tail_fit_shape=tail_fit.get('shape'),tail_fit_scale=tail_fit.get('scale'),
                                   tail_fit_n=tail_fit.get('n_tail'),tail_fit_ks_distance=tail_fit.get('ks_distance'))
                    write_json(pair_dir/'tail_fit.json',tail_fit)
                    windows.to_csv(window_path(pair_dir),sep='\t',index=False)
                    candidates.to_csv(pair_dir/'enriched_regions.tsv',sep='\t',index=False)
                    merged.to_csv(pair_dir/'selected_regions.tsv',sep='\t',index=False)
                    membership.to_csv(pair_dir/'support_transcripts.tsv',sep='\t',index=False)
                    archive=dict(observed_count=observed,expected_probability=expected,
                                 null_max_excess=maxima,null_retained_distances=retained)
                    if keep_windows:
                        archive.update(distance_nt=distance,background_count=background,
                                       window_starts=starts,expected_window_mass=expected_mass)
                    np.savez_compressed(pair_dir/'distributions.npz',**archive)
                    summary['elapsed_seconds']=round(time.time()-pair_started,3)
                    pair_plots=make_plots(pair_dir,a,b,windows,candidates,merged,summary) if plot and fdr_threshold is None else {}
                    summary['plots']=pair_plots;plots[f'{a}|{b}']=pair_plots
                    write_json(pair_dir/'summary.json',summary)
                    pair_checkpoint(pair_dir)
                    return summary
                summaries={}
                for pair,summary in parallel_results(compute_pair,selected_pairs):
                    summaries[tuple(pair)]=summary
                    write_json(destination/'progress.json',dict(state='running',phase='simulations',branch=label,
                               branch_completed=len(summaries),branch_total=len(selected_pairs),jobs=jobs,threads=threads,
                               completed_pair_branches=len(pair_rows)+len(summaries),total_pair_branches=len(selected_pairs)*(2 if masking else 1),
                               last_pair=list(pair),updated_unix=time.time()))
                    if len(summaries)%10==0 or len(summaries)==len(selected_pairs):
                        pd.DataFrame([{key:value for key,value in r.items() if key!='plots'} for r in pair_rows+list(summaries.values())]).to_csv(destination/'pair_summary.tsv',sep='\t',index=False)
                canonical=[summaries[tuple(pair)] for pair in selected_pairs]
                pair_rows.extend(canonical)
                branch_windows=[pd.read_csv(window_path(branch_dir/f'{a}__{b}'),sep='\t') for a,b in selected_pairs] if keep_windows else [None]*len(selected_pairs)
                branch_regions=[pd.read_csv(branch_dir/f'{a}__{b}'/'selected_regions.tsv',sep='\t') for a,b in selected_pairs]
                branch_support=[pd.read_csv(branch_dir/f'{a}__{b}'/'support_transcripts.tsv',sep='\t') for a,b in selected_pairs]
                if fdr_threshold is not None:
                    from statsmodels.stats.multitest import multipletests
                    rows_this_branch = [r for r in pair_rows if r['analysis'] == label]
                    if len(rows_this_branch) != len(declared_pair_family):
                        raise RuntimeError('Incomplete pair family: FDR cannot be finalized')
                    valid = np.array([r['pair_scan_pvalue'] is not None for r in rows_this_branch],dtype=bool)
                    probabilities = np.array([r['pair_scan_pvalue'] if good else 1.0 for r,good in zip(rows_this_branch,valid)],dtype=float)
                    adjusted = multipletests(probabilities,method=fdr_method)[1]
                    fitted_family=np.array([r.get('pair_scan_pvalue_fitted') if r.get('pair_scan_pvalue_fitted') is not None else p for r,p in zip(rows_this_branch,probabilities)],dtype=float)
                    adjusted_fitted=multipletests(fitted_family,method=fdr_method)[1]
                    write_json(destination/'progress.json',dict(state='running',phase='finalizing_fdr',branch=label,
                               completed_pair_branches=len(pair_rows),total_pair_branches=len(selected_pairs)*(2 if masking else 1),updated_unix=time.time()))
                    def finalize_pair(item):
                        j,summary,good,q=item
                        a,b=summary['kmer_a'],summary['kmer_b']
                        pair_dir=branch_dir/f'{a}__{b}'
                        qvalue=float(q) if good else None
                        pass_fdr=bool(good and q < fdr_threshold)
                        summary.update(pair_FDR=qvalue,passes_pair_fdr=pass_fdr,fdr_method=fdr_method,fdr_family_size=len(rows_this_branch),
                                       pair_FDR_fitted=float(adjusted_fitted[j]) if good else None,
                                       fitted_fdr_used_empirical_fallback=bool(good and summary.get('pair_scan_pvalue_fitted') is None))
                        windows=pd.read_csv(window_path(pair_dir),sep='\t')
                        windows['pair_FDR']=np.nan if qvalue is None else qvalue
                        windows['pair_FDR_fitted']=summary['pair_FDR_fitted'] if good else np.nan
                        windows['passes_pair_fdr']=pass_fdr
                        windows['selected_pre_fdr']=windows['selected'].astype(bool)
                        windows['selected']=windows['selected_pre_fdr'] & pass_fdr
                        merged=pd.read_csv(pair_dir/'selected_regions.tsv',sep='\t')
                        membership=pd.read_csv(pair_dir/'support_transcripts.tsv',sep='\t')
                        if not pass_fdr:
                            merged=merged.iloc[:0].copy()
                            membership=membership.iloc[:0].copy()
                        merged['pair_FDR']=np.nan if qvalue is None else qvalue
                        merged['pair_FDR_fitted']=summary['pair_FDR_fitted'] if good else np.nan
                        summary.update(windows_selected=int(windows.selected.sum()),merged_intervals=len(merged))
                        candidates=pd.read_csv(pair_dir/'enriched_regions.tsv',sep='\t')
                        pair_plots=make_plots(pair_dir,a,b,windows,candidates,merged,summary) if plot else {}
                        summary['plots']=pair_plots
                        plots[f'{a}|{b}']=pair_plots
                        windows.to_csv(window_path(pair_dir),sep='\t',index=False)
                        merged.to_csv(pair_dir/'selected_regions.tsv',sep='\t',index=False)
                        membership.to_csv(pair_dir/'support_transcripts.tsv',sep='\t',index=False)
                        write_json(pair_dir/'summary.json',summary)
                        pair_checkpoint(pair_dir)
                        return summary
                    items=[(j,summary,bool(good),float(q)) for j,(summary,good,q) in enumerate(zip(rows_this_branch,valid,adjusted))]
                    for item,updated in parallel_results(finalize_pair,items):
                        j=item[0]
                        if updated is not rows_this_branch[j]:
                            rows_this_branch[j].clear();rows_this_branch[j].update(updated)
                    branch_windows=[pd.read_csv(window_path(branch_dir/f'{a}__{b}'),sep='\t') for a,b in selected_pairs] if keep_windows else [None]*len(selected_pairs)
                    branch_regions=[pd.read_csv(branch_dir/f'{a}__{b}'/'selected_regions.tsv',sep='\t') for a,b in selected_pairs]
                    branch_support=[pd.read_csv(branch_dir/f'{a}__{b}'/'support_transcripts.tsv',sep='\t') for a,b in selected_pairs]
                plots={f"{r['kmer_a']}|{r['kmer_b']}":r.get('plots',{}) for r in pair_rows if r['analysis']==label}
                nonempty_windows=[frame for frame in branch_windows if frame is not None and not frame.empty]
                nonempty_regions=[frame for frame in branch_regions if not frame.empty]
                nonempty_support=[frame for frame in branch_support if not frame.empty]
                all_windows=pd.concat(nonempty_windows,ignore_index=True) if nonempty_windows else pd.DataFrame(columns=window_columns)
                all_regions=pd.concat(nonempty_regions,ignore_index=True) if nonempty_regions else branch_regions[0].copy()
                all_support=pd.concat(nonempty_support,ignore_index=True) if nonempty_support else branch_support[0].copy()
                all_regions.to_csv(branch_dir/'selected_regions.tsv',sep='\t',index=False)
                all_support.to_csv(branch_dir/'support_transcripts.tsv',sep='\t',index=False)
                branch_results[label]=dict(windows=all_windows,windows_in_memory=bool(keep_windows),
                                           window_files=[str(window_path(branch_dir/f'{a}__{b}')) for a,b in selected_pairs],
                                           regions=all_regions,support_transcripts=all_support,plots=plots,output_dir=str(branch_dir))
                pd.DataFrame([{key:value for key,value in r.items() if key!='plots'} for r in pair_rows]).to_csv(destination/'pair_summary.tsv',sep='\t',index=False)
            if hashlib.sha256(source.read_bytes()).hexdigest()!=source_sha:
                raise RuntimeError('Original reference FASTA changed during analysis')
            settings.update(state='complete',elapsed_seconds=round(time.time()-started,3),input_verified_unchanged=True, fdr_computed=fdr_threshold is not None)
            write_json(destination/'settings.json',settings)
            pair_summary=pd.DataFrame([{key:value for key,value in r.items() if key!='plots'} for r in pair_rows])
            result=dict(pair_summary=pair_summary,branches=branch_results,settings=settings,output_dir=str(destination))
            self.cooccurrence_results=result
            return result
        except BaseException as error:
            settings.update(state='interrupted' if isinstance(error,KeyboardInterrupt) else 'failed',error=f'{type(error).__name__}: {error}')
            write_json(destination/'settings.json',settings)
            raise



