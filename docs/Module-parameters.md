# Module parameters

[API signatures](API-reference.md) · [Methodology](Methodology.md) · [Tutorial](Tutorial.md)

All three modules use `from KEA import KEA`. Construct `KEA(dir_out, ref_fasta, ctrl_fasta)` with an output directory and the two FASTA paths, then call `CreateCombinedFasta()` and `KmersCountsTable(...)`. FASTA IDs must be unique across inputs; use uppercase A/C/G/T sequences. For the ordinary two-file workflow leave `species=None` (internal key `input`).

## Shared counting

| Parameter | Meaning and choice |
| :--- | :--- |
| `dir_out` | Directory for generated results; use a separate directory for each analysis. |
| `ref_fasta`, `ctrl_fasta` | Reference and control RNA FASTA files. For KRS, include the target among the reference sequences and choose a biologically appropriate comparison population. |
| `k` | K-mer length, required for counting. The universe contains 4^k types; memory grows rapidly with k. |
| `cores` | Number of counting workers; default 5. Start with 1 for the small demo. Run within a guarded Python script. |
| `subset_by_subseq` | Optional sequence-subset filter for advanced counting; leave unset for the tutorial. |
| `species` | Optional species key for the advanced multi-species workflow; leave unset for ordinary inputs. |

## KEA: `ExtractKmers`

| Parameter | Meaning and choice |
| :--- | :--- |
| `k_selected` | Select the previously counted k-mer length. |
| `extraction_methods` | Request `('delta_median', 'stat_log2fc')` for the integrated workflow, then intersect their enriched or depleted lists explicitly. `delta_rank_percentile` is also available. |
| `n_iter` | Default 1000 control-subset draws when control size is at least twice reference size. Set NumPy's random seed before extraction for reproducibility. |
| `sign_emp_pval` | Tail fraction for population-level selection; default 0.0025. Despite the name, it is not a calibrated P value. Reported benchmarks use 0.01. |
| `log2fc_threshold` | Absolute log2 fold-change cutoff for transcript-level selection; default 1.0. |
| `fdr_alpha` | Significance cutoff applied using the selected correction; default 0.05. |
| `stat_test` | `mwu` (default), `ks` or `ttest`; methodology explains their differences. |
| `correction_method` | `bonferroni` (default), `fdr` or `PValue` for the raw P value. |
| `median` | True (default) uses transcript medians for fold change; False uses means. |

The returned dictionary is indexed by species, then extraction method. Inspect `enriched`, `depleted` and `all_kmers`. More permissive thresholds increase signature size and can change downstream RBP associations.

## KRS: `KRS`

| Parameter | Meaning and choice |
| :--- | :--- |
| `reference` | Exact target transcript FASTA identifier, required. |
| `top_pct` | Upper percentage of percentile ranks to retain; default 5 means target rank ≥0.95. It is not a P-value threshold. |
| `store` | True (default) retains results on the analysis object. |
| `save` | True (default) writes the target signature and rank output. |
| `save_full_rank_table` | False (default); enable to save the complete population rank table, which can be large. |
| `verbose` | True (default) prints progress. |

Returns a list of k-mers. A rare element observed only in the target can attain rank 1; inspect prevalence as well as ranks.

## K-RBP: `KmersProteinEnrichment`

| Parameter | Meaning and choice |
| :--- | :--- |
| `kmers` | KEA or KRS signature list, required. |
| `df_zscore` | Numeric pandas DataFrame with k-mers as index and RBP profiles as columns, required. See [PEKA matrix format](../data/peka/README.md). |
| `k_out` | Default None leaves query length unchanged. Set to matrix k-mer length (e.g. 5) to decompose longer query elements. Cannot exceed query length. |
| `deconvolute_unique` | True (default) deduplicates decomposed elements. |
| `convert_to_rna` | True (default) converts T to U; match the matrix alphabet. |
| `bg_kmers` | Default None uses matrix rows outside the query. A custom list must have valid rows and be disjoint from the query. |
| `mode`, `test` | Use `continuous` and `mwu` (defaults) for the documented workflow. Historical binary options `threshold` and `top_n` are not implemented. |
| `fdr_alpha` | Default 0.05 controls the reported significance summary. Apply your chosen adjusted-P and effect filters explicitly to the returned table. |
| `plot` | True (default) produces a figure; False is useful for scripts or batch analyses. |
| `p_col`, `p_thresh`, `fc_thresh` | Plot significance/effect filters; defaults `fdr`, 0.05 and 1.0. Set these explicitly; the default effect threshold is at the upper bound of Cliff's delta. |
| `plot_type` | `lollipop` (default), `volcano`, `barplot` or `wordcloud`; see implementation for plot-specific options. |
| `highlight`, `highlight_top_n` | Optional highlighted profiles and maximum number of highlighted hits (default 10). |
| `save` | Optional figure output path. Save the result DataFrame separately with `to_csv`. |
| `wordcloud_score`, `wordcloud_top_n`, `wordcloud_font_path` | Word-cloud score selection (default `cliff_delta`), number of entries (10) and optional font file. |

Returns all evaluated profiles with P values and corrections. The column `zscore_diff` contains Cliff's delta, not a difference of Z scores. Positive values indicate higher signature scores. Keep profile/cell-line identity when reporting results.

## Choosing parameters: exact behavior and worked examples

### Counting: `k`, `cores` and `subset_by_subseq`

`k` is the length of the words counted, not a significance setting. The class enumerates all 4^k DNA words and stores one count column per transcript plus a `k` metadata column. Each call replaces the count table for that species; calling k=5 after k=7 does not accumulate both lengths. Therefore `k_selected` must match the currently stored table. At k=3 there are 64 rows, at k=5 there are 1,024 and at k=7 there are 16,384. A length-100 sequence contributes 94 overlapping starts for k=7 when every base is canonical. Longer k increases sequence specificity but also increases sparsity and memory requirements. This distribution corrects the original counter's omission of the final valid start; see the changelog.

`cores` is passed directly to `multiprocessing.Pool`. Increasing it changes parallel execution, not thresholds or the intended statistical definition. Even `cores=1` creates a worker process. Put analysis inside `if __name__ == "__main__":` and begin with one worker before increasing it for large inputs.

`subset_by_subseq` filters the enumerated universe to words containing one of the supplied substrings. For example, `subset_by_subseq=["AT"]` keeps words containing AT. It does not select transcripts or crop their sequences. The implementation concatenates matches without deduplicating words matching multiple substrings. It also normalizes downstream frequencies by the retained counts rather than by every possible word; restricted-universe analyses therefore have a different denominator. Leave it unset for the standard KEA/KRS workflow.

### KEA: why there are two selection methods

The pooled frequency of a word is its total count across a group divided by the total count of all retained words in that group. Longer transcripts can contribute more counted windows. `delta_median` centers each group's pooled frequency distribution around its median, divides by the appropriate lower or upper quartile distance, and subtracts control from reference. It then ranks the resulting scores. A positive score means the word is more prominent relative to the reference distribution's center than it is in the control distribution. It is not a log fold change.

`stat_log2fc` first divides each transcript's word counts by that transcript's total. The statistical test then compares those per-transcript frequencies between the two groups. Consequently this branch asks whether the difference is supported across transcripts, rather than only in the pooled counts. The integrated signature is the intersection of the two branches' enriched lists (or separately their depleted lists). `ExtractKmers` does not compute that intersection automatically.

**`extraction_methods`** controls which branches execute. The default includes `delta_rank_percentile` as well as the two integrated branches. Specify `("delta_median", "stat_log2fc")` when you intend the documented two-branch analysis. The rank branch instead subtracts each word's normalized frequency rank between populations; its tied ranks use `method="first"`, so row order can affect ties.

**`sign_emp_pval`** controls only the two delta branches. For N rows, each tail contains `max(1, int(N * sign_emp_pval))` words. With the full k=7 universe, 0.01 selects 163 enriched and 163 depleted words before intersection; 0.0025 selects 40 in each tail. Increasing it admits more words without changing the scores. It does not control the Mann–Whitney P values or provide a false-positive probability. Even zero would retain one word because of the minimum; use a positive fraction below 0.5 to obtain distinct upper/lower tails and record it explicitly. Zero quartile distances can produce nonfinite delta scores: inspect the score table before interpreting the selected tails.

**`n_iter`** controls how the pooled control frequency is estimated when the control has at least twice as many transcripts as the reference. If there are 24 reference and 48 control transcripts, each iteration selects 24 distinct control transcripts; selections can recur across iterations. The code averages the normalized pooled frequencies across `n_iter` draws. If there are 24 reference and 47 control transcripts, it uses the entire control and this parameter has no effect. This sampling does not replace the full control group in `stat_log2fc`: that test still uses all matching control transcript columns. Larger iteration counts can stabilize the estimated control frequencies but cost time and memory; the implementation retains the intermediate frequency vectors. Use a positive integer and set `np.random.seed(42)` before extraction to repeat the sampling.

**`median`** chooses the summary used for the fold-change filter, not the input to the statistical test. With True, `log2FC = log2((median(reference) + 1e-9)/(median(control) + 1e-9))`; with False the means replace the medians. Means are more influenced by a minority of high-frequency transcripts. If both medians are zero, the median-based fold change is zero even when a smaller subset differs. The 1e-9 pseudocount prevents division by zero and can strongly influence ratios near zero.

**`log2fc_threshold`** applies signed fold-change filters in addition to significance. At 1.0, enriched words require log2FC ≥1 and depleted words require log2FC ≤−1, corresponding to at least a twofold pseudocount-adjusted ratio in either direction. At 0.5 the ratio is approximately 1.414. This parameter does not alter the test P values, the delta tails or the KRS ranks. Use a nonnegative value; smaller values permit weaker frequency differences.

**`stat_test`** accepts `mwu`, `ks` or `ttest`. The shared helper uses a two-sided Mann–Whitney U test, a two-sided Kolmogorov–Smirnov test or Welch's unequal-variance t-test, respectively. These assess different properties of the frequency distributions. The current KEA implementation skips a word whenever both groups are internally constant, even if their constant values differ. Therefore absence from `all_kmers.tsv` can reflect this skip rather than a nonsignificant result.

**`correction_method` and `fdr_alpha`** act together. `bonferroni` filters `bonf`, `fdr` filters the BH-adjusted `FDR` column, and `PValue` filters the unadjusted `PValue` column. KEA uses `<= fdr_alpha`. Both corrected columns are calculated regardless of the chosen filter, over the words that were actually tested. For example, a word with log2FC=1.2, raw P=0.001, BH P=0.02 and Bonferroni P=0.08 passes the fold-change threshold 1.0 and significance cutoff 0.05 with `fdr`, but fails with `bonferroni`. Despite its name, `fdr_alpha` is not necessarily an FDR cutoff: its meaning follows `correction_method`.

```python
selected = analysis.ExtractKmers(
    k_selected=7,
    extraction_methods=("delta_median", "stat_log2fc"),
    n_iter=1000,
    sign_emp_pval=0.01,
    log2fc_threshold=1.0,
    stat_test="mwu",
    correction_method="bonferroni",
    fdr_alpha=0.05,
    median=True,
)["input"]
signature = sorted(set(selected["delta_median"]["enriched"])
                   & set(selected["stat_log2fc"]["enriched"]))
```

These are explicit analysis settings, not a guarantee of optimal settings for every dataset. Keep an empty intersection empty. The legacy Venn helper omits empty sets and is unsuitable for computing this intersection.

### KRS: `reference` and `top_pct`

`reference` is the exact target transcript column name, derived from the first whitespace-delimited token of its FASTA header. The target must already be in the count table; KRS does not load a new target sequence. KRS uses all transcript columns in the selected species table as the population, including both input FASTAs and the target itself. It does not treat the control FASTA alone as a separate ranking population.

For each word, KRS normalizes counts within each transcript, replaces zero frequencies with missing values, and computes percentile ranks across the remaining transcripts. Equal frequencies receive average ranks. Suppose a word's positive frequencies are 0.01, 0.02, 0.02 and 0.08. Their percentile ranks are 0.25, 0.625, 0.625 and 1.0. A target with frequency 0.08 is selected with `top_pct=5`, but a target with 0.02 is not. Transcripts where that word is absent do not enter this ranking denominator.

`top_pct=5` sets a rank cutoff of 0.95; 10 sets 0.90; 1 sets 0.99. Increasing `top_pct` lowers the cutoff and can add words. This is evaluated independently for each word, so the signature need not contain 5% of the word universe. A word unique to the target gets rank 1 even if its absolute count is small. A word absent from the target has a missing rank and is not selected, including when `top_pct=100`. The method validates `0 < top_pct <= 100`.

`save=False` suppresses disk outputs, while `store=False` suppresses storage in `analysis.KEA_results`; neither changes the returned signature. `save_full_rank_table=True` writes the complete matrix only when `save=True`. With `store=True`, the full rank table is still retained in memory regardless of the disk-save setting. The compact `kmer_ranks.tsv` contains `target_frequency`, `target_rank_percentile` and `selected_KRS` for auditing each word.

### K-RBP: query representation and background

`df_zscore` must already contain the scores to compare. This method does not compute PEKA scores or standardize the columns into Z scores. Each column is tested independently as a profile; the column name becomes `protein` in the output. Preserve RBP and experimental context in that name. Supply numeric values, unique row labels and a consistent alphabet/word length.

`k_out` decomposes every query word into overlapping words of that length. For example, `ATCGATC` with `k_out=5` yields `ATCGA`, `TCGAT` and `CGATC`. With RNA conversion these become `AUCGA`, `UCGAU` and `CGAUC`. Use this when a 7-mer signature is compared against a 5-mer matrix. With `k_out=None`, the original words must match matrix rows directly. The method infers the input length from the first query word and validates equal input lengths during decomposition. A custom `bg_kmers` list is decomposed with that same inferred input length, so when `k_out` is supplied its elements must have the query's original length, not merely the output length.

`deconvolute_unique=False` retains repetitions at the decomposition stage, but the subsequent matrix intersection uses `set(kmers)`. Thus the actual statistical query remains unique: this option does not weight repeated words more heavily in continuous-mode testing.

`convert_to_rna=True` replaces uppercase T with U in both matrix row labels and query/background words. It does not uppercase strings, validate ambiguous letters or reverse complement sequences. Ensure conversion does not create duplicate row labels (for example by supplying both DNA and RNA spellings of the same word).

With `bg_kmers=None`, background is every matrix row outside the valid query. With a custom list, only matching matrix rows are kept; overlaps with the query are not removed automatically. An empty matching query or custom background raises an error. If the background contains more than ten times as many words as the query, the code samples exactly 10× query-size background words without replacement using a fixed seed of 42. Missing scores are then dropped independently for each profile. A profile is skipped if either group has fewer than two observations. The reported `n1` and `n2` are therefore profile-specific valid observation counts after sampling and missing-value removal.

### K-RBP: tests, returned statistics and visualization

Use `mode="continuous"`. Other modes do not initialize the result list in this implementation; the signature's historical `threshold` and `top_n` options are unused and do not implement a binary enrichment test.

With `test="mwu"`, the one-sided alternative is that query scores are higher than background scores. This differs from KEA's two-sided test. `test="ks"` calls SciPy's KS test with `alternative="greater"`; that alternative concerns the ordering of the cumulative distribution functions and must not be described as equivalent to the MWU higher-score hypothesis. Use MWU for the documented RBP score-enrichment question.

The returned `zscore_diff` is Cliff's delta, `(query>background pairs − query<background pairs)/(n1*n2)`. It is in [−1,1]. A value 0.6 means an excess of 60 percentage points of favorable pairwise comparisons over unfavorable ones; it does not mean the mean Z score increased by 0.6. In the MWU branch, the separate `effect_size` field uses the opposite sign. Use `zscore_diff` for the positive-effect filter.

`fdr_alpha` sets `significant_fdr` and `significant_bonf` using strict `<`, and controls the printed significance counts. It does not remove rows from the returned DataFrame. BH and Bonferroni corrections cover the evaluated profiles, not every original column when some profiles were skipped. Save all results and apply a documented positive-effect filter explicitly:

```python
rbp = analysis.KmersProteinEnrichment(
    signature, matrix, k_out=5, convert_to_rna=True,
    mode="continuous", test="mwu", fdr_alpha=0.05, plot=False,
)
# Illustrative effect cutoff: choose and report it for your analysis.
candidates = rbp[(rbp["fdr"] < 0.05) & (rbp["zscore_diff"] > 0.3)]
rbp.to_csv("results/rbp_all_profiles.tsv", sep="\t", index=False)
candidates.to_csv("results/rbp_candidates.tsv", sep="\t", index=False)
```

`p_col` chooses the significance column for coloring/selecting plotted profiles (`fdr`, `bonf` or `pvalue`); `p_thresh` sets its plotting cutoff. `fc_thresh` is an effect-size threshold, despite the historical fold-change name. Volcano/lollipop/barplot significance uses `abs(zscore_diff) > fc_thresh`, whereas word-cloud selection uses positive `zscore_diff >= fc_thresh`. At the default `fc_thresh=1.0`, the strict comparison cannot select a significant point because delta is bounded by 1. Set a smaller, explicitly reported effect threshold such as 0.3 for an illustrative plot. These settings change visualization, not the underlying tests or returned rows.

The volcano y-axis always uses the raw `pvalue`, even if colors use `p_col="fdr"`; its horizontal line is at `-log10(p_thresh)` on that raw-P scale. Use the tabular adjusted values to audit the chosen significance rule.

`plot_type="lollipop"` shows plot-selected profiles; `barplot` shows all profiles with significance-dependent colors; `volcano` shows effect versus raw-P evidence; `wordcloud` shows selected positive effects. `highlight` matches case-insensitive regular-expression substrings of profile names and changes highlighting only. `highlight_top_n` labels the largest positive deltas among plot-selected profiles in the volcano plot. `save` accepts a filename or list of filenames for figures; it never saves the statistical table.

For the word cloud, `wordcloud_top_n` limits entries after significance/effect filtering. The current table is sorted by FDR and `.head(...)` is applied before plotting weights: the retained subset is therefore the first passing rows in FDR order, despite the console message saying “by Cliff's delta”. `wordcloud_score="cliff_delta"` uses delta for word size; `median_ratio` uses `median_group1/(median_group2+1)`, which is not a fold enrichment for arbitrary signed scores. Undefined or nonpositive plotting weights raise an error. `wordcloud_font_path` selects a local font file. Neither weight option changes the tested scores or P values.

## Additional analysis parameters

### `KmerCorrelationWithScore`

Pass `txids_and_score` as exactly two columns: transcript identifiers first, numeric scores second. Only identifiers found in the count table are retained. `method="spearman"` correlates normalized word frequencies with score ranks; pandas also accepts `pearson` and `kendall`. `sign_emp_pval` selects the upper/lower fractions of the correlation ranking, not significant correlations. Here selection is based on `rank_percentage = zero_based_rank/N`, so its boundaries differ from KEA's fixed tail counts. Pass `species="input"` explicitly for the ordinary workflow. This analysis does not compute adjusted hypothesis-test P values.

### `KmerSlidingWindowProfile`

`tx_id` selects one transcript from `fasta_file`, or from `ann_df` using `ID`/`ensembl_transcript_id` and `cdna`. With FASTA input, the required positional `ann_df` argument can be None. `window_size` is the window length in nucleotides; `step` is the distance between successive starts. At size 50 and step 25, adjacent windows overlap by 25 nt. Only complete windows are emitted; an incomplete terminal segment is excluded. Use positive integers and a window at least as long as the words of interest.

`kmers` defines the elements profiled. `kmer_groups` supplies grouped lists and supersedes the ordinary list for plotting/group boundaries. Frequencies are overlapping counts divided by `window_size − word_length + 1`. The combined frequency divides the sum of all query-word counts by the sum of those per-word denominators; for equal-length words this is their average frequency, not the fraction of positions matching any signature element. `strand` labels output intervals and does not reverse complement the sequence. `convert_to_rna` changes T to U in sequence and words.

`n_clusters` can be a fixed cluster count or `"auto"`; `cluster_range` bounds the automatic search. `improvement_threshold` controls the relative silhouette-improvement rule. The clustering and plotting flags are advanced options; this distribution's installation/demo checks do not establish biological validity or broad robustness of positional clustering. `save_prefix` saves the ordinary and cluster heatmaps as PDFs; the implementation does not save the line or silhouette figures with this prefix. The ordinary clustermap is constructed even when `plot_heatmap=False`, which closes it after construction. Meanwhile, `plot_line`, `plot_heatmap`, `plot_cluster_heatmap` and `plot_silhouette` toggle their corresponding figures. The method returns the window table, cluster assignments and silhouette-score results.

## Kmap and CDS domain analysis

- [Kmap: similarity metrics, parameters and Sankey visualization](Kmap.md)
- [CDS domain enrichment: inputs, frame, backgrounds and Fisher tests](Domain-enrichment.md)
