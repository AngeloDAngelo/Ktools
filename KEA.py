import os
import multiprocessing as mp
from pathlib import Path
from itertools import product
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
import logomaker
from matplotlib_venn import venn2, venn3
from wordcloud import WordCloud
from Bio import SeqIO
from scipy.stats import ks_2samp, ttest_ind, mannwhitneyu, fisher_exact
from statsmodels.stats.multitest import multipletests

pd.set_option("future.no_silent_downcasting", True)

def FileByList(file,lista):
    
    with open(file,"w") as f:
        for i in lista:
            f.write(i+"\n")
        f.close()


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


class KEA:
    def __init__(self, dir_out, ref_fasta, ctrl_fasta):
        os.makedirs(dir_out, exist_ok=True)
        self.dir_out = dir_out
        self.reference_fasta = ref_fasta
        self.control_fasta = ctrl_fasta
        self.combined_fasta = None
        self.combined_fasta_by_species = {}
        self.kmers_count_table_by_species = {}
        self.KEA_results = None
        self.kmer_protein_scores = None

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
            # Look up registered species-specific FASTAs
            ref_key = f"{species}_ref"
            ctrl_key = f"{species}_ctrl"
            ref_fasta = self.combined_fasta_by_species.get(ref_key)
            ctrl_fasta = self.combined_fasta_by_species.get(ctrl_key)

            if ref_fasta is None or ctrl_fasta is None:
                raise ValueError(f"Species-specific FASTAs not found for {species}. Register the species-specific reference and control FASTAs first.")

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
                raise ValueError(f"Mapped FASTAs for species '{species_key}' not found. Register the species-specific reference and control FASTAs first.")
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


