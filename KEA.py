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
    print("/home/angelo/data01/packages/bedtools2/bin/bedtools intersect "+options+" -a "+file_a+" -b "+file_b+" > "+output)
    os.system("/home/angelo/data01/packages/bedtools2/bin/bedtools intersect "+options+" -a "+file_a+" -b "+file_b+" > "+output)
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
    
    if os.path.exists(output_dir)==False:
        os.makedirs(output_dir)
        
    file_bed=output_dir+file_name
    
    
    
    if chunk_size == "none":
        
        bed_df.to_csv(file_bed,sep="\t",index=None,header=None)

        #print("mappers available:\n")
        mapper_dir="/data01/common_data_files/utilities/mapper_files/"

        lista_mapper=os.listdir(mapper_dir)
        #for i in lista_mapper:
        #    print(i)

        # mapper file
        mapper_file=mapper_dir+"parsed_genomic_"+mapper_specie_release+".txt"

        # launching
        script_bin="/data01/common_scripts/utilities/range_tools/mapper_to_genome.py"
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
            mapper_dir="/data01/common_data_files/utilities/mapper_files/"

            lista_mapper=os.listdir(mapper_dir)
            #for i in lista_mapper:
            #    print(i)

            # mapper file
            mapper_file=mapper_dir+"parsed_genomic_"+mapper_specie_release+".txt"

            # launching
            python_bin="/home/angelo/anaconda3/bin/python"
            script_bin="/data01/common_scripts/utilities/range_tools/mapper_to_genome.py"
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

    kmers = [seq[i:i+k] for i in range(0, len(seq)-k+1)]
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

    def KmersCountsTable(self, k, cores=5, subset_by_subseq=None, species=None):

        if species is None:
            INPUT = self.combined_fasta
            species_key = "input"
        else:
            INPUT = self.combined_fasta_by_species.get(species)
            species_key = species
            if INPUT is None:
                raise ValueError(f"Combined FASTA for species '{species}' not found. Run CreateCombinedFasta first.")

        # -------------------------------
        # Generate all kmers
        # -------------------------------
        prodotto = list(product("ATCG", repeat=k))
        all_kmers = list(map("".join, prodotto))

        if subset_by_subseq:
            selected = []
            for subseq in subset_by_subseq:
                selected += [x for x in all_kmers if subseq in x]
            all_kmers = selected

        print(f"[KEA] Counting kmers (k={k}) for species '{species_key}' | N kmers: {len(all_kmers)}")
        all_kmers_df = pd.DataFrame(all_kmers)

        # -------------------------------
        # Read FASTA as dataframe
        # -------------------------------
        # input_file = pd.read_table(INPUT, header=None, dtype=str)
        # input_seq = input_file[input_file[0].str.count(">") == 0].reset_index(drop=True)
        # input_id = input_file[input_file[0].str.count(">") == 1].reset_index(drop=True)

        records = list(SeqIO.parse(INPUT, "fasta"))
        input_id = [rec.id for rec in records]
        input_seq = [str(rec.seq) for rec in records]
        

        

        all_k = [k] * len(input_seq)
        all_kmers_list = [all_kmers] * len(input_seq)
        #all_ids = [x.replace(">", "") for x in input_id[0].tolist()]
        all_ids = input_id

        # -------------------------------
        # Count kmers using multiprocessing
        # -------------------------------
        with mp.Pool(cores) as pool:
            kmers_list = pool.starmap(kmers_counter,zip(input_seq, all_k, all_kmers_list, all_ids)) #input_seq[0].tolist()


        kmer_counts = pd.concat(kmers_list, axis=1)
        kmer_counts["k"] = k

        # -------------------------------
        # Store in species-specific dictionary
        # -------------------------------
        if not hasattr(self, "kmers_count_table_by_species"):
            self.kmers_count_table_by_species = {}

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
       
        species_key = species if species is not None else "input"
        self.KEA_results = self.KEA_results or {}
        self.KEA_results.setdefault(species_key, {})

        # ------------------------------
        # Load kmer table and groups
        # ------------------------------
        if species_key == "input":
            kmer_table = self.kmers_count_table_by_species["input"]
            ref_group = get_fasta_ids(self.reference_fasta)
            comp_group = get_fasta_ids(self.control_fasta)
        else:
            # Species-specific kmer table and mapped FASTAs
            if species_key not in self.kmers_count_table_by_species:
                raise ValueError(f"Kmer table for species '{species_key}' not found. Run KmersCountsTable first.")

            kmer_table = self.kmers_count_table_by_species[species_key]

            ref_key = f"{species_key}_ref"
            ctrl_key = f"{species_key}_ctrl"

            if ref_key not in self.combined_fasta_by_species or ctrl_key not in self.combined_fasta_by_species:
                raise ValueError(f"Mapped FASTAs for species '{species_key}' not found. Run FetchOrthologsAndSaveFastaOneToOne first.")

            ref_group = get_fasta_ids(self.combined_fasta_by_species[ref_key])
            comp_group = get_fasta_ids(self.combined_fasta_by_species[ctrl_key])

        # Only keep transcript IDs that exist in kmer_table
        ref_group = [x for x in ref_group if x in kmer_table.columns]
        comp_group = [x for x in comp_group if x in kmer_table.columns]

        if len(ref_group) == 0 or len(comp_group) == 0:
            raise ValueError(f"No matching transcripts found for species '{species_key}' in kmer table.")

        # ------------------------------
        # Subset by k
        # ------------------------------
        kmer_matrix = kmer_table.loc[kmer_table["k"] == k_selected].copy()
        del kmer_matrix["k"]

        # ------------------------------
        # Compute frequencies
        # ------------------------------
        perc_ref = kmer_matrix[ref_group].sum(axis=1) / kmer_matrix[ref_group].sum().sum()
        perc_ctrl_raw = kmer_matrix[comp_group].sum(axis=1) / kmer_matrix[comp_group].sum().sum()

        # Automatic control iteration if control group is much larger
        use_iteration = len(comp_group) >= 2 * len(ref_group)

        if use_iteration:
            ref_size = len(ref_group)
        
            freq_ctrl_iter = []
        
            for _ in range(n_iter):
                sampled = np.random.choice(comp_group, ref_size, replace=False)
                subset = kmer_matrix[sampled]
        
                # SAME normalization as perc_ref
                freq = subset.sum(axis=1) / subset.sum().sum()
        
                freq_ctrl_iter.append(freq)
        
            perc_ctrl = pd.concat(freq_ctrl_iter, axis=1).mean(axis=1)
        
        else:
            perc_ctrl = perc_ctrl_raw

        # ------------------------------
        # Base output directory
        # ------------------------------
        base_out = os.path.join(self.dir_out, "ExtractKmers", species_key)
        os.makedirs(base_out, exist_ok=True)

        # ------------------------------
        # Loop over extraction methods
        # ------------------------------
        for method in extraction_methods:
            print(f"[KEA] {species_key} | {method}")

            # --------------------------
            # DELTA METHODS
            # --------------------------
            if method in {"delta_rank_percentile", "delta_median"}:

                if method == "delta_rank_percentile":
                    score = (
                        perc_ref.rank(method="first") / len(perc_ref) -
                        perc_ctrl.rank(method="first") / len(perc_ctrl)
                    )
                    score_col = "delta_rank_percentile"

                else:  # delta_median
                   
                    q1_r, med_r, q3_r = np.percentile(perc_ref, [25, 50, 75])
                    q1_c, med_c, q3_c = np.percentile(perc_ctrl, [25, 50, 75])

                    diff_ref = np.where(
                        perc_ref < med_r,
                        (perc_ref - med_r) / (med_r - q1_r),
                        (perc_ref - med_r) / (q3_r - med_r)
                    )

                    diff_ctrl = np.where(
                        perc_ctrl < med_c,
                        (perc_ctrl - med_c) / (med_c - q1_c),
                        (perc_ctrl - med_c) / (q3_c - med_c)
                    )

                    score = diff_ref - diff_ctrl
                    score_col = "delta_median"

                all_kmers = pd.DataFrame({
                    "freq_in_group_1": perc_ref,
                    "freq_in_group_2": perc_ctrl,
                    score_col: score
                }).sort_values(score_col, ascending=False)

                top_n = max(1, int(len(all_kmers) * sign_emp_pval))
                enriched = all_kmers.head(top_n).index.tolist()
                depleted = all_kmers.tail(top_n).index.tolist()

                subdir = f"{method}_emp{sign_emp_pval}"

            # --------------------------
            # STAT + LOG2FC
            # --------------------------
            elif method == "stat_log2fc":

                testing = kmer_matrix.T
                testing = testing.div(testing.sum(axis=1), axis=0).fillna(0)

                test_ref = testing.loc[ref_group]
                test_ctrl = testing.loc[comp_group]

                rows = []
                for kmer in kmer_matrix.index:
                    if kmer not in test_ref.columns:
                        continue

                    x, y = test_ref[kmer].values, test_ctrl[kmer].values
                    #print(y)
                    if len(x) == 0 or len(y) == 0:
                        continue
                    if np.all(x == x[0]) and np.all(y == y[0]):
                        continue

                    stat, pval = run_stat_test(x, y, stat_test)
                    if median:
                        log2fc = np.log2((np.median(x)+1e-9)/(np.median(y)+1e-9))  # or mean
                        rows.append({
                            "kmer": kmer,
                            "log2FoldChange": log2fc,
                            "stat": stat,
                            "PValue": pval,
                            "Median_group1": np.median(x),
                            "Median_group2": np.median(y)
                        })
                    else:
                        log2fc = np.log2((np.mean(x)+1e-9)/(np.mean(y)+1e-9)) 
                        rows.append({
                            "kmer": kmer,
                            "log2FoldChange": log2fc,
                            "stat": stat,
                            "PValue": pval,
                            "Mean_group1": np.mean(x),
                            "Mean_group2": np.mean(y)
                        })

                all_kmers = pd.DataFrame(rows).set_index("kmer")

                _, all_kmers["bonf"], _, _ = multipletests(all_kmers["PValue"], method="bonferroni")
                _, all_kmers["FDR"], _, _ = multipletests(all_kmers["PValue"], method="fdr_bh")

                all_kmers["freq_in_group_1"] = perc_ref.loc[all_kmers.index]
                all_kmers["freq_in_group_2"] = perc_ctrl.loc[all_kmers.index]

                if correction_method == "bonferroni":
                    signif = all_kmers["bonf"] <= fdr_alpha
                
                elif correction_method == "fdr":
                    signif = all_kmers["FDR"] <= fdr_alpha
                
                elif correction_method == "PValue":
                    signif = all_kmers["PValue"] <= fdr_alpha   # usa p-value grezzo
                
                else:
                    raise ValueError(f"Unknown correction_method: {correction_method}")

                enriched = all_kmers[(all_kmers["log2FoldChange"] >= log2fc_threshold) & signif].index.tolist()
                depleted = all_kmers[(all_kmers["log2FoldChange"] <= -log2fc_threshold) & signif].index.tolist()

                subdir = f"{method}_log2fc{log2fc_threshold}_{correction_method}_alpha{fdr_alpha}"

            else:
                raise ValueError(f"Unknown extraction method: {method}")

            # --------------------------
            # SAVE RESULTS
            # --------------------------
            outdir = os.path.join(base_out, subdir)
            os.makedirs(outdir, exist_ok=True)

            FileByList(os.path.join(outdir, "enriched.txt"), enriched)
            FileByList(os.path.join(outdir, "depleted.txt"), depleted)
            all_kmers.to_csv(os.path.join(outdir, "all_kmers.tsv"), sep="\t")

            self.KEA_results[species_key][method] = {
                "enriched": enriched,
                "depleted": depleted,
                "all_kmers": all_kmers
            }

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


    def DomainEnrichment(
        self,
        gtf_file=None,
        uniprot_dir=None,
        kmers=None, method="stat_log2fc", species=None,
        in_frame=False, background="control", tracks=("unipDomain",),
        alpha=0.05, fdr_scope="global", split_blocks=True,
        seed=42, output_dir=None, plot=True, max_domains=60,
        ref_fasta=None, ctrl_fasta=None,
    ):
        """CDS k-mer / UniProt-domain enrichment (GRCh38 inputs required).

        Supply explicit paths to a matching GRCh38 GTF and UniProt BED directory.
        FASTA inputs must be FULL spliced transcripts in transcript orientation,
        with transcript IDs as the first header token, matching the GTF release.
        Default: same k-mer occurrences in reference versus control transcripts.
        background='shuffle': one seeded random CDS position per reference hit,
        preserving transcript, length and (when in_frame=True) codon-start frame.
        Fisher tests on occurrence counts are exploratory: overlapping occurrences
        and repeated hits within transcripts are not independent biological samples.
        Domain assignment is genomic, same-strand, >=1 nucleotide overlap; it does
        not establish an isoform-specific protein-domain annotation.

        in_frame=False (default) keeps ALL CDS occurrences, regardless of frame.
        in_frame=True means START at a codon boundary; k need not be divisible by 3.
        GTF CDS excludes the stop codon. Partial CDS phases are respected.
        All observed domains across both groups are tested for every selected k-mer.
        FDR is BH across the whole matrix by default. Untestable cells remain NaN.
        Heatmap: corrected log2 OR; * means FDR <= alpha, grey means untestable.
        Raw OR (including 0/inf) and raw p-values are retained in the table.
        Output always goes to a NEW directory; existing directories are refused.
        Returns dict and sets self.domain_enrichment.
        """
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

        if gtf_file is None or uniprot_dir is None:
            raise ValueError("Supply gtf_file and uniprot_dir for matching GRCh38 resources")
        if background not in {"control", "shuffle"}:
            raise ValueError("background must be 'control' or 'shuffle'")
        if fdr_scope not in {"global", "kmer"} or not 0 < alpha < 1:
            raise ValueError("Invalid fdr_scope or alpha")
        if max_domains is not None and max_domains < 1:
            raise ValueError("max_domains must be positive or None")
        if output_dir is not None and Path(output_dir).exists():
            raise FileExistsError(f"Refusing to overwrite existing directory: {output_dir}")
        key = species or "input"
        if kmers is None:
            try:
                kmers = self.KEA_results[key][method]["enriched"]
            except (KeyError, TypeError):
                raise ValueError("Run ExtractKmers first or supply kmers explicitly") from None
        if isinstance(kmers, str):
            kmers = [kmers]
        kmers = list(dict.fromkeys(str(k).upper().replace("U", "T") for k in kmers))
        if not kmers or any(not re.fullmatch("[ACGT]+", k) for k in kmers):
            raise ValueError("Supply at least one DNA/RNA k-mer without ambiguous bases")
        if key != "input":
            raise ValueError("DomainEnrichment supports the human input KEA object; species must be None or input")
        if isinstance(tracks, str):
            tracks = (tracks,)
        allowed = {"unipDomain", "unipInterest", "unipStruct", "unipLocTransMemb",
                   "unipLocExtra", "unipLocSignal", "unipLocCytopl", "unipRepeat"}
        if not tracks or not set(tracks) <= allowed:
            raise ValueError("Unsupported track; use a UniProt feature track")

        def txid(s):
            return re.sub(r"\.\d+$", "", str(s))

        def fasta(path):
            result = {}
            op = gzip.open if str(path).endswith(".gz") else open
            with op(path, "rt") as handle:
                for rec in SeqIO.parse(handle, "fasta"):
                    name = txid(rec.id)
                    if name in result:
                        raise ValueError(f"Duplicate normalized FASTA transcript ID: {name}")
                    result[name] = str(rec.seq).upper().replace("U", "T")
            return result

        reference_path = ref_fasta if ref_fasta is not None else self.reference_fasta
        control_path = ctrl_fasta if ctrl_fasta is not None else self.control_fasta
        ref = fasta(reference_path)
        ctrl = fasta(control_path) if background == "control" else {}
        if set(ref) & set(ctrl):
            raise ValueError("Reference and control transcript sets must be disjoint")
        seqs = {**ref, **ctrl}
        models = defaultdict(lambda: {"exon": [], "CDS": []})
        op = gzip.open if str(gtf_file).endswith(".gz") else open
        with op(gtf_file, "rt") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) != 9 or f[2] not in {"exon", "CDS"}:
                    continue
                match = re.search(r'transcript_id "([^"]+)"', f[8])
                if not match or txid(match[1]) not in seqs:
                    continue
                name = txid(match[1])
                chrom = f[0] if f[0].startswith("chr") else "chr" + f[0]
                if chrom == "chrMT":
                    chrom = "chrM"
                m = models[name]
                if "chrom" in m and (m["chrom"], m["strand"]) != (chrom, f[6]):
                    raise ValueError(f"Conflicting GTF loci for {name}")
                m.update(chrom=chrom, strand=f[6])
                m[f[2]].append((int(f[3])-1, int(f[4]), f[7]))

        # Spatial bins avoid a genome-wide interval scan for every occurrence.
        bins = defaultdict(list)
        bin_size = 100000
        for track in dict.fromkeys(tracks):
            with open(Path(uniprot_dir) / (track + ".bed")) as handle:
                for line in handle:
                    if not line.strip() or line.startswith(("#", "track", "browser")):
                        continue
                    f = line.rstrip("\n").split("\t")
                    a, b = int(f[1]), int(f[2])
                    label = f[26] if track in {"unipDomain", "unipInterest"} and len(f) > 26 else f[3]
                    if track == "unipInterest" and label != "Disordered":
                        continue
                    if track == "unipLocSignal" and label != "Signal peptide":
                        continue
                    label = label.strip()
                    if not label:
                        continue
                    if len(set(tracks)) > 1:
                        label = track + ": " + label
                    blocks = [(a, b)]
                    if split_blocks and len(f) >= 12:
                        sizes = [int(x) for x in f[10].rstrip(",").split(",")]
                        starts = [int(x) for x in f[11].rstrip(",").split(",")]
                        if len(sizes) != int(f[9]) or len(starts) != len(sizes):
                            raise ValueError("Malformed UniProt BED12 block count")
                        blocks = [(a+s, a+s+n) for s, n in zip(starts, sizes)]
                    for lo, hi in blocks:
                        if not a <= lo < hi <= b:
                            raise ValueError("Invalid UniProt BED block coordinates")
                        for bucket in range(lo//bin_size, (hi-1)//bin_size+1):
                            bins[(f[0], f[5], bucket)].append((lo, hi, label))

        rows, qc = [], []
        rng = np.random.default_rng(seed)
        for name, seq in seqs.items():
            m = models.get(name)
            reason = None
            if m is None or not m["CDS"]:
                reason = "missing_GTF_or_non_coding"
            elif not m["exon"]:
                reason = "missing_exons"
            if reason:
                qc.append((name, reason))
                continue
            exons = sorted(set((a, b) for a, b, _ in m["exon"]), reverse=m["strand"] == "-")
            genomic_order = sorted(exons)
            if any(b > c for (a, b), (c, d) in zip(genomic_order, genomic_order[1:])):
                raise ValueError(f"Overlapping GTF exons for {name}")
            offsets, offset = [], 0
            for a, b in exons:
                offsets.append((a, b, offset))
                offset += b-a
            if offset != len(seq):
                qc.append((name, "FASTA_GTF_length_mismatch"))
                continue
            cds = []
            for a, b, phase in m["CDS"]:
                contained = False
                for lo, hi, off in offsets:
                    if lo <= a < b <= hi:
                        start = off + (a-lo if m["strand"] == "+" else hi-b)
                        cds.append((start, start+b-a, int(phase)))
                        contained = True
                        break
                if not contained:
                    raise ValueError(f"CDS is not contained in an exon: {name}")
            cds = sorted(set(cds))
            if any(b != c for (a, b, _), (c, d, _) in zip(cds, cds[1:])):
                qc.append((name, "non_contiguous_CDS"))
                continue
            cs, ce, phase0 = cds[0][0], cds[-1][1], cds[0][2]
            anchor = cs + phase0
            if any(p != (anchor-a) % 3 for a, b, p in cds):
                qc.append((name, "inconsistent_CDS_phase"))
                continue
            # Project genomic domain blocks onto this spliced transcript.
            projected = set()
            for lo, hi, off in offsets:
                candidates = set()
                for bucket in range(lo//bin_size, (hi-1)//bin_size+1):
                    candidates.update(bins.get((m["chrom"], m["strand"], bucket), ()))
                for a, b, label in candidates:
                    left, right = max(lo, a), min(hi, b)
                    if left < right:
                        start = off + (left-lo if m["strand"] == "+" else hi-right)
                        projected.add((start, start+right-left, label))

            def record(k, pos, group):
                domains = tuple(sorted({label for a, b, label in projected if a < pos+len(k) and pos < b}))
                rows.append((name, k, pos, pos+len(k), (pos-anchor) % 3, group, domains))

            for k in kmers:
                pos = seq.find(k, cs, ce)
                while pos >= 0:
                    if not in_frame or (pos >= anchor and (pos-anchor) % 3 == 0):
                        record(k, pos, "reference" if name in ref else "control")
                        if background == "shuffle":
                            first = anchor if in_frame else cs
                            step = 3 if in_frame else 1
                            n_positions = (ce-len(k)-first)//step+1
                            shuffled = first + step*int(rng.integers(n_positions))
                            record(k, shuffled, "control")
                    pos = seq.find(k, pos+1, ce)
            qc.append((name, "used"))
        occurrences = pd.DataFrame(rows, columns=["transcript", "kmer", "start", "end", "frame", "group", "domains"])
        quality = pd.DataFrame(qc, columns=["transcript", "status"])
        if len(quality) and (quality.status != "used").any():
            warnings.warn("Some transcripts excluded; inspect result['qc'] for reasons", stacklevel=2)
        if occurrences.empty:
            raise ValueError("No eligible CDS occurrences. Check IDs, full-length FASTAs, GTF and frame filter. " + str(quality.status.value_counts().to_dict()))
        totals = Counter(zip(occurrences.kmer, occurrences.group))
        hits = Counter()
        for row in occurrences.itertuples(index=False):
            hits.update((row.kmer, row.group, domain) for domain in row.domains)
        domains = sorted({domain for _, _, domain in hits})
        if not domains:
            raise ValueError("No same-strand UniProt overlaps; check assembly, coordinates and selected tracks")
        stats = []
        for k in kmers:
            n1, n0 = totals[k, "reference"], totals[k, "control"]
            for domain in domains:
                a, c = hits[k, "reference", domain], hits[k, "control", domain]
                b, d = n1-a, n0-c
                testable = n1 > 0 and n0 > 0 and 0 < a+c < n1+n0
                odds, p = fisher_exact([[a, b], [c, d]], alternative="two-sided") if testable else (np.nan, np.nan)
                # Continuity correction is for visualization ONLY, never Fisher counts.
                delta = 0.5 if 0 in (a, b, c, d) else 0.0
                logor = np.log2(((a+delta)*(d+delta))/((b+delta)*(c+delta))) if testable else np.nan
                stats.append((k, domain, a, b, c, d, odds, logor, p, testable))
        table = pd.DataFrame(stats, columns=["kmer", "domain", "A", "B", "C", "D", "odds_ratio", "log2_or", "pvalue", "testable"])
        table["fdr"] = np.nan
        groups = [table.index] if fdr_scope == "global" else table.groupby("kmer").groups.values()
        for idx in groups:
            valid = table.loc[idx, "pvalue"].dropna().index
            if len(valid):
                table.loc[valid, "fdr"] = multipletests(table.loc[valid, "pvalue"], method="fdr_bh")[1]
        table["significant"] = table.fdr.le(alpha)
        table["enriched"] = table.odds_ratio.gt(1)
        table["significant_enrichment"] = table.significant & table.enriched
        matrix = table.pivot(index="kmer", columns="domain", values="log2_or").reindex(kmers)
        qmatrix = table.pivot(index="kmer", columns="domain", values="fdr").reindex(kmers)
        figure = None
        plotted = list(matrix.columns)
        if plot:
            import matplotlib.pyplot as plt
            import seaborn as sns
            if max_domains is not None and len(plotted) > max_domains:
                plotted = list(qmatrix.min().sort_values(kind="stable").head(max_domains).index)
            shown = matrix[plotted]
            marks = qmatrix[plotted].map(lambda v: "*" if v <= alpha else "") if hasattr(qmatrix, "map") else qmatrix[plotted].applymap(lambda v: "*" if v <= alpha else "")
            figure, ax = plt.subplots(figsize=(max(8, len(plotted)*0.48), max(3, len(kmers)*0.32)))
            ax.set_facecolor("#dddddd")
            sns.heatmap(shown, mask=shown.isna(), cmap="RdBu_r", center=0,
                        annot=marks, fmt="", ax=ax, xticklabels=True, yticklabels=True,
                        cbar_kws={"label": "log2 OR (0.5 correction if zero cells)"})
            ax.set(xlabel="UniProt domain", ylabel="k-mer", title=f"CDS domains: reference vs {background}; * BH FDR <= {alpha}")
            figure.tight_layout()
        if output_dir is None:
            base = Path(self.dir_out)
            base.mkdir(parents=True, exist_ok=True)
            destination = Path(tempfile.mkdtemp(prefix="DomainEnrichment_", dir=base))
        else:
            destination = Path(output_dir)
            destination.mkdir(parents=True, exist_ok=False)
        table.to_csv(destination / "statistics.tsv", sep="\t", index=False)
        occurrences.to_csv(destination / "occurrences.tsv", sep="\t", index=False)
        quality.to_csv(destination / "qc.tsv", sep="\t", index=False)
        matrix.to_csv(destination / "log2_or.tsv", sep="\t")
        qmatrix.to_csv(destination / "fdr.tsv", sep="\t")
        if figure is not None:
            figure.savefig(destination / "heatmap.png", dpi=180, bbox_inches="tight")
            figure.savefig(destination / "heatmap.pdf", bbox_inches="tight")
        import json
        settings = dict(gtf_file=str(gtf_file), uniprot_dir=str(uniprot_dir), kmers=kmers,
                        method=method, in_frame=in_frame, background=background, tracks=list(tracks),
                        alpha=alpha, fdr_scope=fdr_scope, split_blocks=split_blocks, seed=seed,
                        reference_fasta=str(reference_path), control_fasta=str(control_path),
                        plotted_domains=plotted, species=key, plot=plot, max_domains=max_domains)
        (destination / "settings.json").write_text(json.dumps(settings, indent=2))
        result = dict(statistics=table, occurrences=occurrences, qc=quality,
                      log2_or=matrix, fdr=qmatrix, figure=figure, output_dir=str(destination), settings=settings)
        self.domain_enrichment = result
        return result


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
        # 1. DECONVOLUTE KMERS (optional)
        # ----------------------------------------------------------------
        if k_out is not None:
            k_in = len(kmers[0])
            if k_out > k_in:
                raise ValueError(f"k_out ({k_out}) must be <= k_in ({k_in})")
    
            
            print(f"[KEA] Deconvoluting {len(kmers)} kmers: k={k_in} -> k={k_out}")
    
            def _deconv(lst):
                out = set() if deconvolute_unique else []
                for kmer in lst:
                    if len(kmer) != k_in:
                        raise ValueError(f"Kmer '{kmer}' has length {len(kmer)}, expected {k_in}")
                    for i in range(k_in - k_out + 1):
                        sub = kmer[i:i + k_out]
                        if deconvolute_unique:
                            out.add(sub)
                        else:
                            out.append(sub)
                return list(out)
    
            kmers = _deconv(kmers)
    
            if bg_kmers is not None:
                
                print(f"[KEA] Deconvoluting bg_kmers: {len(bg_kmers)} kmers")
                bg_kmers = _deconv(bg_kmers)
                print(f"[KEA] After deconvolution: {len(kmers)} kmers | " f"bg={len(bg_kmers) if bg_kmers is not None else 'default'}")
    
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
        kmers_valid = sorted(set(kmers).intersection(index_kmers))
        n_missing = len(set(kmers)) - len(kmers_valid)
    
       
        print(f"[KEA] Kmers of interest : {len(set(kmers))} total | "f"{len(kmers_valid)} found in df | {n_missing} missing")
    
        if len(kmers_valid) == 0:
            raise ValueError("No valid kmers found in df_zscore index.")
    
        # ----------------------------------------------------------------
        # NEW: VALIDATE bg_kmers
        # ----------------------------------------------------------------
        if bg_kmers is not None:
            bg_kmers_valid = sorted(set(bg_kmers).intersection(index_kmers))
            print(f"[KEA] Custom bg_kmers: {len(bg_kmers)} total | {len(bg_kmers_valid)} found")
            if len(bg_kmers_valid) == 0:
                raise ValueError("None of the provided bg_kmers are present in df_zscore index.")
        else:
            bg_kmers_valid = None
    
        # ----------------------------------------------------------------
        # 4. RUN ENRICHMENT (FIXED)
        # ----------------------------------------------------------------
        if mode == "continuous":
    
            if bg_kmers_valid is None:
                background_kmers = sorted(index_kmers - set(kmers_valid))
            else:
                background_kmers = bg_kmers_valid
    
            MAX_BG_RATIO = 10
            if len(background_kmers) > MAX_BG_RATIO * len(kmers_valid):
                rng = np.random.default_rng(42)
                background_kmers = list(
                    rng.choice(background_kmers,
                               size=MAX_BG_RATIO * len(kmers_valid),
                               replace=False)
                )
    
            # Cliff’s delta
            def cliffs_delta(x, y):
                import numpy as np
                x = np.asarray(x)
                y = np.asarray(y)
                n = len(x)
                m = len(y)
                greater = sum((xi > y).sum() for xi in x)
                less    = sum((xi < y).sum() for xi in x)
                return (greater - less) / (n * m)
    
            results = []
    
            for col in df_zscore.columns:
                g1 = df_zscore.loc[kmers_valid, col].dropna().values
                g2 = df_zscore.loc[background_kmers, col].dropna().values
    
                if len(g1) < 2 or len(g2) < 2:
                    continue
    
                mean1, mean2 = np.mean(g1), np.mean(g2)
                med1,  med2  = np.median(g1), np.median(g2)
                sd1,   sd2   = np.std(g1), np.std(g2)
    
                if test == "mwu":
                    stat, pval = mannwhitneyu(g1, g2, alternative="greater")
                    effect = 1 - (2 * stat) / (len(g1) * len(g2))
                elif test == "ks":
                    stat, pval = ks_2samp(g1, g2, alternative="greater")
                    effect = stat
                else:
                    raise ValueError(f"Unknown test '{test}'.")
    
                # REPLACE Δ‑median WITH CLIFF’S DELTA
                zscore_diff = cliffs_delta(g1, g2)
    
                results.append({
                    "protein":        col,
                    "mean_group1":    mean1,
                    "mean_group2":    mean2,
                    "median_group1":  med1,
                    "median_group2":  med2,
                    "sd1":            sd1,
                    "sd2":            sd2,
                    "effect_size":    effect,
                    "pvalue":         pval,
                    "zscore_diff":    zscore_diff,
                    "n1":             len(g1),
                    "n2":             len(g2)
                })
    
        # ----------------------------------------------------------------
        # 5. MULTIPLE TESTING CORRECTION
        # ----------------------------------------------------------------
        results = pd.DataFrame(results)
    
        if len(results) == 0:
            print("[KEA] Warning: no results generated.")
            return results
    
        results["fdr"]  = multipletests(results["pvalue"], method="fdr_bh")[1]
        results["bonf"] = multipletests(results["pvalue"], method="bonferroni")[1]
        results["significant_fdr"]  = results["fdr"]  < fdr_alpha
        results["significant_bonf"] = results["bonf"] < fdr_alpha
        results = results.sort_values("fdr").reset_index(drop=True)
    
        n_sig_fdr  = results["significant_fdr"].sum()
        n_sig_bonf = results["significant_bonf"].sum()
        print(f"[KEA] Significant proteins — FDR<{fdr_alpha}: {n_sig_fdr} | "f"Bonferroni<{fdr_alpha}: {n_sig_bonf}")
    

        # ----------------------------------------------------------------
        # 6. VISUALIZATION (VOLCANO + OPTIONAL ALTERNATIVES)
        # ----------------------------------------------------------------
        if plot:
            df_plot = results.copy(); df_plot["neglog10p"] = -np.log10(df_plot["pvalue"] + 1e-300)
            sig_mask = (df_plot[p_col] < p_thresh) & (np.abs(df_plot["zscore_diff"]) > fc_thresh)
            hl_mask = pd.Series(False, index=df_plot.index)
            if highlight is not None:
                if isinstance(highlight, str): highlight = [highlight]
                hl_mask = df_plot["protein"].str.contains("|".join(highlight), case=False, na=False)
            hl_top_mask = pd.Series(False, index=df_plot.index)
            if highlight_top_n is not None:
                sig = df_plot[sig_mask]; top_idx = sig.nlargest(highlight_top_n, "zscore_diff").index
                hl_top_mask = df_plot.index.isin(top_idx)
        
            # ================================================================  # 6A. VOLCANO PLOT  # ================================================================
            if plot_type == "volcano":
                fig, ax = plt.subplots(figsize=(9, 7))
                ax.scatter(df_plot.loc[~sig_mask, "zscore_diff"], df_plot.loc[~sig_mask, "neglog10p"], color="lightgrey", alpha=0.5, s=25, label="n.s.")
                ax.scatter(df_plot.loc[sig_mask & ~hl_mask, "zscore_diff"], df_plot.loc[sig_mask & ~hl_mask, "neglog10p"], color="tomato", alpha=0.8, s=35, label="significant")
                ax.scatter(df_plot.loc[hl_mask, "zscore_diff"], df_plot.loc[hl_mask, "neglog10p"], color="royalblue", edgecolor="black", s=110, label="highlight")
                for _, row in df_plot.loc[hl_mask | hl_top_mask].iterrows(): ax.text(row["zscore_diff"], row["neglog10p"], row["protein"], fontsize=8, ha="left", va="bottom")
                ax.axhline(-np.log10(p_thresh), linestyle="--", color="grey", linewidth=0.8); ax.axvline(fc_thresh, linestyle="--", color="grey", linewidth=0.8); ax.axvline(-fc_thresh, linestyle="--", color="grey", linewidth=0.8)
                ax.set_xlabel("Cliff’s delta (group1 vs background)"); ax.set_ylabel("-log10(pvalue)"); ax.set_title("Volcano plot (Cliff’s delta)"); ax.legend(fontsize=8)
                if save is not None:
                        filenames = save if isinstance(save, (list, tuple)) else [save]
                        for filename in filenames:
                            fig.savefig(filename, dpi=600, bbox_inches="tight", facecolor="white")
                plt.tight_layout(); plt.show(); plt.close()
        
            # ================================================================  # 6D. LOLLIPOP PLOT  # ================================================================
            elif plot_type == "lollipop":
                sig = df_plot[sig_mask].sort_values("zscore_diff")
                fig, ax = plt.subplots(figsize=(9, max(4, len(sig) * 0.35)))
                ax.hlines(sig["protein"], 0, sig["zscore_diff"], color="grey", alpha=0.6)
                ax.scatter(sig["zscore_diff"], sig["protein"], c="tomato", s=40)
                sig_hl_mask = sig["protein"].isin(df_plot.loc[hl_mask, "protein"])
                ax.scatter(sig.loc[sig_hl_mask, "zscore_diff"], sig.loc[sig_hl_mask, "protein"], c="royalblue", edgecolor="black", s=110)
                # NO LABELS HERE — axis already shows protein names
                ax.set_xlabel("Cliff’s delta"); ax.set_ylabel("Protein"); ax.set_title("Lollipop plot (significant proteins)")
                if save is not None:
                        filenames = save if isinstance(save, (list, tuple)) else [save]
                        for filename in filenames:
                            fig.savefig(filename, dpi=600, bbox_inches="tight", facecolor="white")
                plt.tight_layout(); plt.show()

            elif plot_type == "barplot":
                df_bar = df_plot.copy().sort_values("zscore_diff")
                colors = ["tomato" if ((row[p_col] < p_thresh) and (abs(row["zscore_diff"]) > fc_thresh)) else "lightgrey" for _, row in df_bar.iterrows()]
                fig, ax = plt.subplots(figsize=(max(9, len(df_bar)*0.25), 6))
                ax.bar(df_bar["protein"], df_bar["zscore_diff"], color=colors)
                ax.set_xlim(-0.5, len(df_bar)-0.5)
                ax.axhline(0, color="black", linewidth=0.8)
                ax.set_ylabel("Cliff’s delta")
                ax.set_xlabel("Protein")
                ax.set_title("Cliff’s delta — All proteins (barplot)")
                plt.xticks(rotation=90)
                plt.tight_layout()
                if save is not None:
                        filenames = save if isinstance(save, (list, tuple)) else [save]
                        for filename in filenames:
                            fig.savefig(filename, dpi=600, bbox_inches="tight", facecolor="white")
                plt.show()

            elif plot_type == "wordcloud":
                # Select proteins ONLY by significance and Cliff's delta
                df_wc = df_plot.loc[(df_plot[p_col] < p_thresh) & (df_plot["zscore_diff"] >= fc_thresh)].copy()
                df_wc = df_wc.dropna(subset=["protein", "zscore_diff"])
                df_wc = df_wc.loc[np.isfinite(df_wc["zscore_diff"])].copy()
                df_wc["protein"] = df_wc["protein"].astype(str).str.strip()
                df_wc = df_wc.loc[df_wc["protein"].ne("")]
                df_wc = df_wc.drop_duplicates(subset="protein", keep="first")

                if wordcloud_top_n is not None:
                    if not isinstance(wordcloud_top_n, (int, np.integer)) or wordcloud_top_n < 1:
                        raise ValueError("wordcloud_top_n must be a positive integer or None.")
                    df_wc = df_wc.head(wordcloud_top_n)

                if df_wc.empty:
                    raise ValueError("No proteins passed the selected thresholds.")

                # Selection is now fixed: calculate plotting scores afterwards
                if wordcloud_score == "cliff_delta":
                    df_wc["wordcloud_score"] = df_wc["zscore_diff"]
                    score_label = "Cliff’s delta"
                    title_score = "Cliff’s delta"
                elif wordcloud_score == "median_ratio":
                    denominator = df_wc["median_group2"] + 1
                    df_wc["wordcloud_score"] = df_wc["median_group1"] / denominator.mask(np.isclose(denominator, 0))
                    score_label = "Median group 1 / (median group 2 + 1)"
                    title_score = "median-ratio score"
                else:
                    raise ValueError("wordcloud_score must be 'cliff_delta' or 'median_ratio'.")

                # Do not remove or replace selected proteins based on plotting scores
                invalid = ~np.isfinite(df_wc["wordcloud_score"]) | (df_wc["wordcloud_score"] <= 0)
                if invalid.any():
                    names = df_wc.loc[invalid, "protein"].tolist()
                    raise ValueError(f"Selected proteins have undefined or non-positive plotting scores: {names}")

                # Audit the exact proteins passed to WordCloud
                print(f"\n[KEA] WordCloud selection: top {len(df_wc)} by Cliff's delta")
                print(df_wc[["protein", "zscore_diff", "wordcloud_score"]].to_string(index=False))

                score_by_protein = dict(zip(df_wc["protein"], df_wc["wordcloud_score"]))
                word_dict = score_by_protein.copy()
                score_min, score_max = df_wc["wordcloud_score"].min(), df_wc["wordcloud_score"].max()

                # Colour scale: higher scores are darker
                if score_min == score_max:
                    padding = max(abs(score_min) * 0.01, 0.01)
                    norm = mpl.colors.Normalize(vmin=score_min - padding, vmax=score_max + padding)
                else:
                    norm = mpl.colors.Normalize(vmin=score_min, vmax=score_max)

                cmap = mpl.colors.LinearSegmentedColormap.from_list("rbp_score_greys", ["#BDBDBD", "#737373", "#000000"])

                def score_color_func(word, **kwargs):
                    return mpl.colors.to_hex(cmap(norm(score_by_protein[word])))

                # Choose an installed font
                if wordcloud_font_path is None:
                    available_fonts = {font.name for font in fm.fontManager.ttflist}
                    font_name = next((name for name in ["Helvetica", "Liberation Sans", "Arial", "DejaVu Sans"] if name in available_fonts), "DejaVu Sans")
                    font_path = fm.findfont(fm.FontProperties(family=font_name))
                else:
                    font_path = str(wordcloud_font_path)
                    fm.fontManager.addfont(font_path)
                    font_name = fm.FontProperties(fname=font_path).get_name()

                # Automatic starting size and relative scaling
                canvas_width, canvas_height = 1600, 700
                wc = WordCloud(
                    width=canvas_width, height=canvas_height, background_color="white",
                    color_func=score_color_func, font_path=font_path,
                    min_font_size=12, max_font_size=None, relative_scaling="auto",
                    prefer_horizontal=1.0, random_state=42, collocations=False,
                    max_words=len(word_dict), repeat=False, margin=6
                ).generate_from_frequencies(word_dict)

                if len(wc.layout_) != len(word_dict):
                    raise ValueError("Not all proteins fit. Increase canvas_width or canvas_height, or reduce min_font_size.")

                # Plot
                with plt.rc_context({"font.family": font_name, "font.size": 12, "axes.titleweight": "bold", "axes.labelweight": "bold"}):
                    fig, ax = plt.subplots(figsize=(15, 6), dpi=150)
                    ax.imshow(wc.to_array(), interpolation="bilinear")
                    ax.axis("off")
                    ax.set_title(f"RBP enrichment by {title_score}", fontsize=16, fontweight="bold", pad=15)

                    # Colour legend
                    scalar_mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
                    scalar_mappable.set_array([])
                    cbar = fig.colorbar(scalar_mappable, ax=ax, fraction=0.025, pad=0.025, shrink=0.80)
                    cbar.set_label(score_label, fontsize=12, fontweight="bold", labelpad=12)

                    if score_min == score_max:
                        cbar.set_ticks([score_min])

                    for label in cbar.ax.get_yticklabels():
                        label.set_fontsize(10)
                        label.set_fontweight("bold")

                    cbar.outline.set_linewidth(0.8)
                    fig.tight_layout()

                    if save is not None:
                        filenames = save if isinstance(save, (list, tuple)) else [save]
                        for filename in filenames:
                            fig.savefig(filename, dpi=600, bbox_inches="tight", facecolor="white")

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

    # def ContaminationAnalysis(
    #     self,
    #     bed_cds_genomic: pd.DataFrame,
    #     uniprot_bed_dir: str,
    #     uniprot_types: list,
    #     analysis_name: str = "UniProtZipcodes"):
        
    #     """
    #     Contamination analysis: enrichment of each k-mer in UniProt protein regions.
    #     Clean version: 
    #     - NO shuffle
    #     - NO synaptic
    #     - NO transcript logic
    #     - NO mapping (bed_cds_genomic is already genomic)
    #     - Uses str.contains to match kmers inside the 'name' column
    #     """
    #     out_dir = self.dir_out + '/contamination_analysis/'
    #     os.makedirs(out_dir, exist_ok=True)
    #     # -----------------------------
    #     # 1. Load UniProt regions
    #     # -----------------------------
    #     def load_uniprot_track(track_type):
    #         tm = pd.read_table(f"{uniprot_bed_dir}/unipLocTransMemb.bed", header=None).iloc[:, :6]
    #         tm.columns = ["chrom","start","end","name","score","strand"]
    
    #         em = pd.read_table(f"{uniprot_bed_dir}/unipLocExtra.bed", header=None).iloc[:, :6]
    #         em.columns = ["chrom","start","end","name","score","strand"]
    
    #         sig = pd.read_table(f"{uniprot_bed_dir}/unipLocSignal.bed", header=None).iloc[:, :6]
    #         sig.columns = ["chrom","start","end","name","score","strand"]
    #         sig = sig[sig["name"] == "Signal peptide"]
    
    #         cyto = pd.read_table(f"{uniprot_bed_dir}/unipLocCytopl.bed", header=None).iloc[:, :6]
    #         cyto.columns = ["chrom","start","end","name","score","strand"]
    
    #         struct = pd.read_table(f"{uniprot_bed_dir}/unipStruct.bed", header=None).iloc[:, :6]
    #         struct.columns = ["chrom","start","end","name","score","strand"]
    
    #         dom = pd.read_table(f"{uniprot_bed_dir}/unipDomain.bed", header=None)
    #         dom = dom.iloc[:, [0,1,2,26,4,5]]
    #         dom.columns = ["chrom","start","end","name","score","strand"]
    
    #         inter = pd.read_table(f"{uniprot_bed_dir}/unipInterest.bed", header=None)
    #         inter = inter.iloc[:, [0,1,2,26,4,5]]
    #         inter.columns = ["chrom","start","end","name","score","strand"]
    
    #         if track_type == "TransExtracellSignalCyto":
    #             return pd.concat([tm, em, sig, cyto], ignore_index=True)
    
    #         if track_type == "DomainsDisordered":
    #             return pd.concat([dom, inter[inter["name"] == "Disordered"]], ignore_index=True)
    
    #         if track_type == "Struct":
    #             return struct.copy()
    
    #         raise ValueError(f"Unknown UniProt type: {track_type}")
    
    #     # merge all requested UniProt tracks
    #     all_re = pd.concat([load_uniprot_track(t) for t in uniprot_types], ignore_index=True)
    
    #     # -----------------------------
    #     # 2. Extract unique kmers
    #     # -----------------------------
    #     zipcodes = sorted(list(set(bed_cds_genomic["name"])))
    
    #     # -----------------------------
    #     # 3. For each kmer: intersect vs UniProt
    #     # -----------------------------
    #     results = []
    
    #     for kmer_family in zipcodes:
    
    #         # select occurrences of this kmer
    #         study_kmer = bed_cds_genomic[ bed_cds_genomic["name"].str.contains(kmer_family) ]
    
    #         # select all other kmers
    #         others = bed_cds_genomic[ ~bed_cds_genomic["name"].str.contains(kmer_family) ]
    
    #         # intersect
    #         inter_kmer = IntersectBedDf(study_kmer, all_re, "-s -wao")
    #         inter_others = IntersectBedDf(others, all_re, "-s -wao")
    
    #         # region name is column 9 (0-based)
    #         region_col = inter_kmer.columns[9]
    
    #         # count occurrences
    #         A_counts = inter_kmer.groupby(region_col)[inter_kmer.columns[3]].nunique()
    #         C_counts = inter_others.groupby(region_col)[inter_others.columns[3]].nunique()
    
    #         # total occurrences
    #         tot_kmer = len(study_kmer)
    #         tot_others = len(others)
    
    #         # build Fisher table per region
    #         for region in A_counts.index:
    
    #             if region == ".":
    #                 continue
    
    #             A = A_counts.get(region, 0)
    #             B = tot_kmer - A
    
    #             C = C_counts.get(region, 0)
    #             D = tot_others - C
    
    #             df_tmp = pd.DataFrame({
    #                 "A":[A], "B":[B], "C":[C], "D":[D]
    #             })
    
    #             stat = FisherExactTest_MultiTest(
    #                 df_tmp,
    #                 "KMER", "A","B",
    #                 "OTHERS","C","D"
    #             )
    
    #             stat["zipcode_family"] = kmer_family
    #             stat["region"] = region
    
    #             results.append(stat)
    
    #     # -----------------------------
    #     # 4. Combine results
    #     # -----------------------------
    #     stat_final = pd.concat(results, ignore_index=True)
    
    #     # log2_or_sign
    #     stat_final["log2_or_sign"] = 0
    #     stat_final.loc[stat_final["fdr"] < 0.05, "log2_or_sign"] = stat_final["log2_or"]
    
    #     # -----------------------------
    #     # 5. Pivot matrix
    #     # -----------------------------
    #     pivot = stat_final.pivot_table(
    #         index="zipcode_family",
    #         columns="region",
    #         values="log2_or_sign",
    #         dropna=True
    #     ).fillna(0)
    
    #     # save
    #     pivot.to_csv(f"{out_dir}/score_matrix.{analysis_name}.txt", sep="\t")
    
    #     return pivot, stat_final
