# Module parameters

[API signatures](API-reference) · [Methodology](Methodology) · [Tutorial](Tutorial)

All modules use `from KEA import KEA`. Construct `KEA(dir_out, ref_fasta, ctrl_fasta)` with an output directory and the two FASTA paths, then call `CreateCombinedFasta()` and `KmersCountsTable(...)`. FASTA IDs must be unique across inputs; use uppercase A/C/G/T sequences. For the ordinary two-file workflow leave `species=None` (internal key `input`).

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
| `sign_emp_pval` | Tail fraction for population-level selection; default 0.0025. Despite the name, it is not a calibrated P value. |
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
| `df_zscore` | Numeric pandas DataFrame with k-mers as index and RBP profiles as columns, required. See [PEKA matrix format](https://github.com/AngeloDAngelo/Ktools/blob/main/data/peka/README.md). |
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


### Original k-mer PEKA score table

Every `KmersProteinEnrichment` call additionally creates `analysis.kmer_protein_scores`, without changing parameters, enrichment tests, plots or the returned profile-level DataFrame. Rows are the original supplied words, before deconvolution; columns are exactly the supplied matrix profiles, including cell-line identifiers. Original order is retained and duplicate normalized words are represented once. Row labels use U when `convert_to_rna=True`.

For a 7-mer, each cell is the arithmetic mean of its three overlapping 5-mer scores. Repeated pentamers count at each position, regardless of `deconvolute_unique`. Missing, nonnumeric or nonfinite component scores make that word/profile cell NaN: an average of fewer than three scores is not reported. Lookup accepts DNA/RNA notation after uppercase normalization. Duplicate normalized matrix rows are averaged within each profile for this additional table only. Words shorter than five have unavailable scores; longer words use all overlapping pentamers. A matrix of another word length cannot supply these pentamer means and will yield unavailable cells.

These are summaries of supplied PEKA scores, **not measured binding affinities or new significance tests**. Profile columns are not merged by protein name. The added table is independent of background sampling and plot significance filters.

```python
results = analysis.KmersProteinEnrichment(
    signature, peka_matrix, k_out=5, plot=False,
)
results.to_csv("rbp_enrichment.tsv", sep="\t", index=False)
analysis.kmer_protein_scores.to_csv("original_kmer_peka_scores.tsv", sep="\t")
```

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

## K-map: transcript similarity and Sankey visualization

K-map compares a target RNA's complete k-mer frequency profile with candidate transcript profiles. The Python method is **`CompareTranscriptKmerProfiles`**. It is independent of KEA/KRS signature extraction and returns ranked similarity scores, not significance tests.

### Compare profiles

| Parameter | Default | Behavior |
| --- | --- | --- |
| `target_tx` | Required | Exact target transcript ID in the target count table. |
| `other_txs` | Required | List of comparison IDs, or a FASTA path when `fasta_input=True`. For one ID pass `["id"]`, not a bare string. |
| `target_species` | None | Count-table key, default `input`. |
| `other_species` | None | Candidate count-table key; default same as target, ignored in FASTA mode. |
| `metric` | `spearman` | `spearman`, `pearson`, `cosine` or `jsd`, defined below. |
| `pseudocount` | 1e-9 | Added to every frequency after per-transcript normalization. It is not an extra raw count. |
| `k` | None | Used only for external FASTA counting; inferred from the first target-table word if omitted. Match the target word length. |
| `fasta_input` | False | True automatically counts the supplied FASTA, using the counting method's default five workers. Use a guarded script. |
| `fasta_species_key` | `_external_db` | Key for caching the external table. If already present it is reused without checking the new file or k. Use a new key for a changed FASTA/word length. |

Counts are normalized within each transcript over its complete table. Comparison uses the intersection of row labels between target/candidate tables; frequencies are not renormalized to that intersection except implicitly by the Jensen–Shannon distance calculation. Consistent full universes are preferred.

- Spearman: correlation of frequency ranks; tied frequencies receive average ranks.
- Pearson: linear correlation of frequency values.
- Cosine: `1 - cosine_distance`, measuring vector direction.
- `jsd`: `1 - scipy.spatial.distance.jensenshannon(...)`, the complement of the square-root Jensen–Shannon divergence, using SciPy's default logarithm base. It is not `1 - raw_divergence`.

Correlation values can be negative. These metrics have different scales and should not be treated as calibrated probabilities. Undefined scores or exceptions in the metric calculation are converted to zero by the inherited implementation; even an unsupported metric can silently yield zeros. Validate metric names and inspect empty/constant profiles. Positive pseudocounts can make empty count profiles appear uniform rather than empty.

The return is a DataFrame indexed by candidate ID, with `<metric>_similarity`, sorted descending. It does not save itself. Including the query among candidates permits a self-comparison.

```python
similarity = analysis.CompareTranscriptKmerProfiles(
    target_tx="reference_006",
    other_txs=["reference_001", "control_001"],
    metric="spearman",
)
similarity.to_csv("results/kmap_similarity.tsv", sep="\t")

## Alternative external candidate FASTA:
## similarity = analysis.CompareTranscriptKmerProfiles(
##     "reference_006", "candidate_transcripts.fa",
##     metric="cosine", fasta_input=True, k=3,
##     fasta_species_key="candidate_set_k3",
## )
```

### Sankey-style rank alignment

`plot_kmer_rank_alignment` visualizes a chosen pair and a background population. It does not choose the candidate automatically; select a candidate ID from the comparison results.

| Parameter | Default | Behavior |
| --- | --- | --- |
| `ref_tx`, `ctrl_tx` | Required | Query and comparison transcript IDs. Plot labels are Query and Control. |
| `ref_species`, `ctrl_species`, `bg_species` | `input` | Count-table keys for query, comparison and background. |
| `kmers` | None | All common words, or an explicit subset. Subset matching uses set intersection; tied ranks may therefore depend on set-derived ordering. |
| `exclude_ctrl_from_bg` | True | Drop the comparison transcript from the background. The query is not automatically excluded. Use a separate background table when both should be excluded. |
| `top_n` | None | Keep the N most frequent query words after the minimum-frequency filter. |
| `min_freq` | None | Keep query frequencies strictly greater than this value; frequency units, not counts. |
| `bin_size` | 100 | Number of words per band, after sorting by query rank. Use a positive integer; smaller bins show more detail. One-word bins are supported. |
| `alpha_band` | 0.4 | Band opacity; choose between 0 and 1. |
| `color_by` | `gc` | `gc`, `ref_rank` or `category`. |
| `kmer_categories` | None | pandas Series indexed by word, required for category coloring. Missing labels become Unknown/grey. Dominant category colors each bin; category ties can depend on ordering. |
| `figsize` | `(10,10)` | Matplotlib figure dimensions in inches. |
| `save` | None | Figure filename; extension selects format, e.g. PNG/PDF. Parent directory must exist. |

Each transcript is normalized separately. Background frequency is the unweighted mean of normalized transcript frequencies, not pooled counts. After filters, each column is ranked among the retained words only, using descending ranks with `method="first"`. Rank positions are 1/N through 1, so smaller positions mean higher frequency. This differs from KRS, which ranks one word across transcripts.

Bands group consecutive query-ranked words and connect their mean ranks in Control, Query and Background. Half-width is `max(bin_word_count/(2*N), 0.002)`, so a minimum visible width applies. This is a binned alignment, not independent flow for every individual word.

GC coloring uses mean GC fraction per band; reference-rank coloring uses its mean query rank; category coloring uses its most frequent category. The inherited tab10 palette has limited category capacity.

The title's similarity is `1 - mean(abs(query_rank - comparison_rank))`, calculated for the plotted words. It is not the selected K-map correlation/cosine/JSD score and is not a P value. Filtering changes its interpretation.

```python
analysis.plot_kmer_rank_alignment(
    ref_tx="reference_006",
    ctrl_tx="control_001",
    ref_species="input",
    ctrl_species="input",
    bg_species="input",
    top_n=32,
    bin_size=4,
    alpha_band=0.4,
    color_by="gc",
    save="results/kmap_rank_alignment.png",
)
```

The method saves if requested, displays and closes the figure, prints the two rank-alignment similarities and returns None. Ensure a nonempty word set and nonempty background remain after filtering. For external candidates, use their cached table key as `ctrl_species`.

## K-CDS: CDS occurrences and UniProt features

`analysis.K_CDS(...)` performs the occurrence-level UniProt association analysis and, by default, directly creates the two final panels: domain enrichment and CDS positional enrichment with coding/shuffle frame bars and In Frame Index. No subsequent notebook cells or external R script preparation are needed. Occurrence-level Fisher tests remain exploratory; genomic overlaps do not establish isoform-specific protein function.

### Required resources and software

Supply a matching GRCh38 GTF with exon/CDS features and numeric phases, a directory of the selected UniProt BED tracks, and full spliced reference/control transcript FASTAs in transcript orientation. UTRs must be retained. GTF/FASTA IDs are matched after stripping terminal version suffixes for the CDS projection; duplicate normalized IDs are rejected. GTF and FASTA inputs may be gzip-compressed. UniProt inputs are uncompressed genomic BED files. GTF coordinates are converted to zero-based half-open coordinates and chromosome names prefixed with `chr`; `chrMT` becomes `chrM`. BED names are used unchanged.

The final panels require `Rscript` on PATH, R packages **ComplexHeatmap** and **circlize**, and a graphics-capable R installation. Install these separately from the Python package:

```r
install.packages(c("BiocManager", "circlize"))
BiocManager::install("ComplexHeatmap")
```

BEDTools on PATH is needed only when missing positional profiles must be computed through `EnrichmentAnalysis`. It is not used when complete cached profiles are available. The embedded R script is written to the output directory and run by Python. Dependency failures stop before the expensive positional analysis. `plot=False` disables rendering by default and requires neither R nor positional inputs; domain/frame tables are still calculated.

### Parameters

| Parameter | Default | Behavior |
| --- | --- | --- |
| `gtf_file`, `uniprot_dir` | None | Required matching GTF and directory of `<track>.bed` files. |
| `ref_fasta`, `ctrl_fasta` | None | Call-specific overrides; otherwise use the constructor's paths. The control FASTA is ignored by shuffle-domain counting, but is needed if reference-versus-control positional enrichment must be computed. |
| `kmers` | None | List or single canonical DNA/RNA word; uppercase, U→T and deduplicate. None uses `KEA_results[species][method]["enriched"]`, not an automatic method intersection. Final panels require at least two words. |
| `method`, `species` | `stat_log2fc`, None | Stored extraction-result key and input species. Only None/`input` species is supported. |
| `background` | `shuffle` | Domain tests compare reference occurrences with randomized intervals. `control` instead uses actual control-transcript occurrences. Positional profiles separately compare Enr with Ctrl. |
| `tracks` | `("unipLocCytopl", "unipLocExtra", "unipLocSignal", "unipLocTransMemb")` | Feature files to load. Labels receive track prefixes when more than one track is selected. |
| `in_frame` | False | Retain all CDS-contained occurrences; True restricts domain-analysis occurrences and their shuffled backgrounds to codon-boundary starts, respecting the first CDS phase. Frame distributions always include all three positions. |
| `alpha`, `fdr_scope` | 0.05, `global` | Inclusive BH cutoff; `global` corrects across testable word/feature cells, `kmer` separately within each word. |
| `split_blocks` | True | Project BED12 blocks rather than entire outer intervals. |
| `seed` | 42 | Independent seeded NumPy generators for domain shuffle and the unrestricted frame-distribution shuffle. |
| `output_dir` | None | Create a unique `K_CDS_*` directory under `dir_out`; an explicit destination must not exist. |
| `plot` | True | Controls final-panel generation unless `final_plots` is explicitly set. The old Matplotlib domain and CODING/SHUFFLE figures are removed. |
| `final_plots` | None | None follows `plot`; False disables R rendering and positional preparation. True can render panels even with `plot=False`. |
| `positional_enrichment` | None | Optional DataFrame or TSV with `kmer`, `bin`, `log2_or`; `fdr` is retained when present. A supplied input must contain complete finite profiles for selected words. |
| `bed_feature` | None | Optional occurrence BED DataFrame or headered TSV used only to calculate missing positional profiles. |
| `recompute_positional` | False | True bypasses positional caches and recalculates profiles; False reuses them and calculates missing/incompatible profiles only. |
| `rscript` | `Rscript` | Executable name or path for the R renderer. |
| `frame_threshold` | 0.5 | Inclusive maximum coding-frame fraction required to display a domain row. Below-threshold or unavailable rows are marked `N.C.`; this does not classify their origin as noncoding. Must be in [0,1]. |
| `start_fraction` | 0.10 | Leading fraction of CDS bins used to name the clusters; count is rounded upward and bounded to leave at least one remaining bin. Must be strictly between 0 and 1. |
| `color_cap` | 2 | Positive finite clipping bound, applied before centering/clustering; also sets the shared color scale. This affects cluster membership, not just colors. |
| `max_domains` | 60 | Maximum displayed domain columns, ranked by minimum FDR, with stable original-order ties; None shows all. Full statistics and matrices are retained. |
| `show_final` | True | Display both final PNGs in IPython/Jupyter; outside IPython files are still saved. |

### Positional input reuse and automatic calculation

With `recompute_positional=False`, the method first looks for per-word `EnrichmentAnalysis` entries in `self.enrichment`, identified by exact `<kmer>_Enr_VS_Ctrl` analysis labels. Repeated stored entries use the last row per bin. Pooled signature analyses such as `ZipCodes_Enr_VS_Ctrl` cannot provide individual word profiles. If no reusable per-word entries exist, the previous `self.k_cds` positional table/output directory and then `dir_out/K_CDS/positional_enrichment.tsv` and `dir_out/positional_enrichment.tsv` are considered. An explicit `positional_enrichment` overrides discovery.

If profiles are missing, the method loads both full-transcript FASTAs and calls `EnrichmentAnalysis` separately for missing words. It uses `self.metatx` when available; otherwise it looks for one compatible `metatx*.bed` file under `dir_out`, or reconstructs a metatranscript from the matching FASTA/GTF using mean regional length proportions. Reconstructed regional bin counts are three times rounded percentage proportions, with a minimum of one bin per region; zero-width transcript bins are omitted. The reconstructed coordinates can differ from a separately supplied annotation table: retain the original metatranscript/cache to reproduce an earlier analysis exactly.

Occurrence coordinates come from `bed_feature`, then a saved `dir_out/kmers_bed.bed`, or an automatic scan for overlapping exact matches in both FASTA populations. Automatic positional counting uses the existing KEA helper's **k=7** convention; for other word lengths, supply complete positional estimates. All computed and cached profiles must share the same CDS bin set; incompatible partial caches are recalculated. `self.enrichment` is restored after internal calls, and the assembled profiles are stored in the returned result and saved for reuse. Original FASTA IDs are used for transcript-coordinate enrichment and must match `self.metatx`/BED IDs.

### Clusters, ordering and panels

CDS bins are sorted numerically. Complete unmasked finite log2 OR estimates are clipped to ±`color_cap`, then each row is centered on its own mean. Euclidean distance and Ward.D2 are applied to these profiles and the tree is cut into two groups. The cluster with the larger mean contrast between the first `start_fraction` of bins and the remaining bins is called **CDS start**; the other is **CDS body**. The contrast uses clipped, uncentered values. These are relative labels, not independent significance tests, and body rows may still have positive initial contrast. The numbers of words are data-dependent and are not fixed at 22/52.

Both panels use body first, then start, with decreasing **In Frame Index** inside each group and DNA sequence as a tie-breaker. In Frame Index is 100 times the largest coding-frame fraction. The positional heatmap displays clipped values without centering. Its grey/white/black bars represent first/second/third starting codon positions for coding and shuffle groups. Missing frame groups have unavailable fractions and explicit missing-data labels; finite groups sum to 100%. K-mer labels use U. Domain rows below `frame_threshold` are marked `N.C.`; stars mark FDR ≤`alpha` in either direction.

Frame 0/1/2 is `(start - anchor) % 3`, with anchor equal to CDS start plus first GTF phase. Domain shuffle relocates intervals; it does not shuffle nucleotides or require randomized intervals to contain the word. Random starts can repeat. The frame-distribution background draws one unrestricted random start per reference CDS occurrence, independently of domain filtering and background mode. No frame significance test is performed.

### Supported tracks

`unipDomain`, `unipInterest`, `unipStruct`, `unipLocTransMemb`, `unipLocExtra`, `unipLocSignal`, `unipLocCytopl`, `unipRepeat`. For `unipDomain` and `unipInterest`, use BED column 27 when present, otherwise column 4. `unipInterest` retains only `Disordered`; `unipLocSignal` only `Signal peptide`. Other tracks use column 4.

### Statistical interpretation

For each selected word and each feature label observed in either group:

| | Overlaps this feature | Does not overlap this feature |
| --- | --- | --- |
| Reference occurrences | A | B |
| Control/shuffled occurrences | C | D |

A same occurrence can contribute to several different feature tests if its interval overlaps multiple labels; repeated blocks of the same label count only once for that occurrence.

A two-sided Fisher exact test compares these occurrence counts. Raw odds ratio is `(A*D)/(B*C)`. Tests require nonzero group totals and at least one overlapping and one non-overlapping occurrence across groups; otherwise the cell remains untestable with NaN statistics.

Raw odds ratios, including zero/infinity, and P values are preserved. If a contingency-table cell is zero, 0.5 is added to all four cells only for the displayed `log2_or`; Fisher counts are unchanged.

Positive log2 odds ratio indicates enrichment; negative values indicate depletion. `significant_enrichment` requires adjusted significance and raw odds ratio >1. Heatmap stars indicate significant association in either direction, not enrichment only. Grey cells are untestable.

Overlapping words and repeated hits in a transcript are dependent observations. These P values are exploratory; do not interpret them as transcript-level replicated inference. Genomic overlap does not establish that a feature belongs to a particular protein isoform.

### Example

```python
result = analysis.K_CDS(
    gtf_file="resources/annotation.gtf",
    uniprot_dir="data/uniprot",
    kmers=signature,
    background="shuffle",
    plot=True,
)
# Both final panels are already saved and displayed in a notebook.
print(result["figures_dir"])
print(result["cluster_assignments"])
print(result["final_figures"]["main_domains"]["svg"])
```

To reuse an explicitly chosen positional export, add `positional_enrichment="previous_run/positional_enrichment.tsv"`. To compute missing profiles, the reference/control FASTAs and BEDTools must be available. No cells need to run after `K_CDS`.

For annotation/statistics only:

```python
result = analysis.K_CDS(
    gtf_file="resources/annotation.gtf",
    uniprot_dir="data/uniprot",
    kmers=signature,
    tracks=("unipDomain",),
    background="control",
    plot=False,
)
```

### Outputs and quality control

The returned dictionary is stored as both `analysis.k_cds` and `analysis.domain_enrichment`.

- `statistics`, `occurrences`, `qc`, `log2_or`, `fdr`, `frame_distribution`, `settings`, `output_dir`: original tables and metadata, saved as their matching TSV/JSON files.
- `positional_enrichment`: assembled positional DataFrame; also `positional_enrichment.tsv` when final panels are enabled.
- `cluster_assignments`: per-word DataFrame with group, maximum frame fraction, preferred frame, capped/raw initial contrasts, In Frame Index and RNA label. None without final panels.
- `final_figures`: dictionaries of `pdf`, `png`, `svg` paths for `main_domains` and `supplement_positions_frames`; empty without final panels.
- `figures_dir`: final-panel directory, present when rendering succeeds. Contains `cluster_assignments.tsv`, both panels in all three formats, and `plot_settings.txt` with the R session information.
- `heatmap_script.R` and `heatmap_R.log`: generated renderer and its stdout/stderr in the run directory.
- `figure` and `frame_figure`: compatibility keys set to None. The old `heatmap.png/pdf` and `frame_distribution.png/pdf` are no longer generated.

QC exclusions include missing/noncoding models, missing exons, sequence/GTF length mismatch, noncontiguous CDS and inconsistent phase. Conflicting loci, overlapping exons, invalid BED blocks and CDS outside an exon raise errors. Empty eligible occurrence sets or absent same-strand UniProt overlaps also raise errors. GTF CDS excludes the stop codon. On an R rendering failure the existing domain/frame tables remain saved, with the R log for diagnosis.


## K-miR: miRNA family seed annotation

`AnnotateMiRNASeeds` merges RNA sequence signatures with family seed annotations. It reports candidate seed-compatible miRNAs, not a new miRNA-binding probability or a family-level enrichment test.

#### Orientation and site definition

The supplied `Seed+m8` is positions 2-8 of the mature miRNA. All 2,606 human rows were checked against their mature sequences. [TargetScan defines families using these positions](https://www.targetscan.org/docs/seed.html).

To match a transcript-oriented query word, reverse complement the seed. Example: miR-7-5p seed `GGAAGAC` gives target `GUCUUCC` in RNA, or `GTCTTCC` in DNA.

This module annotates exact **7mer-m8** sites only. It does not include wobble, mismatches, 6mer, 7mer-A1, 8mer/context classification, accessibility, conservation scoring or miRNA expression. An 8mer includes a matching 7mer-m8 segment, but this join does not classify it as an 8mer. See [TargetScan site definitions](https://www.targetscan.org/docs/7mer.html).

### Parameters

```python
AnnotateMiRNASeeds(self, kmers=None, source="KEA",
    reference=None, species=None, species_id=9606,
    plot=True, top_n=20, output_dir=None)
```

| Parameter | Behavior |
| --- | --- |
| `analysis.seed_family_file` | Path to the tab-separated family table. Defaults to the complete table installed with K-tools; set this attribute to use a custom file. Required columns: `miR family`, `Seed+m8`, `Species ID`, `MiRBase ID`. If `Mature sequence` is present, positions 2–8 are validated. |
| `species_id` | Taxonomy ID for row filtering, default human 9606. The bundled human subset cannot annotate other species. |
| `source` | `KEA`: score join from stored delta-median and stat_log2fc tables; `KRS`: score join from stored KRS target results; `manual`: sequence annotation without analysis scores. |
| `kmers` | Explicit canonical 7-mers. If omitted in KEA, use the intersection of enriched delta-median and stat_log2fc lists. If omitted in KRS, use the stored selected signature. Required for manual mode. DNA/RNA alphabets are accepted and normalized to DNA. Empty signatures return empty matches. Other lengths raise an error rather than silently decomposing or transferring parent scores. |
| `reference` | Optional validation of the stored KRS target ID. To analyze another target, rerun KRS first. |
| `species` | Count/results key, default `input`; distinct from taxonomy `species_id`. |
| `plot` | Draw/save scored annotation plots, default True. Manual mode has no analysis scores and does not plot. |
| `top_n` | Positive integer, default 20; limit plotted matched word/family rows, not returned annotations. KEA sorts by log2FoldChange; KRS by target percentile. |
| `output_dir` | New output directory; None creates a unique `MiRNASeeds_*` subdirectory of the analysis output. Existing directories are refused. |

### Analysis scores and interpretation

KEA has an integrated **selection** but no defined scalar combined score. The method retains both `delta_median` and `log2FoldChange` and plots them in separate panels for matched words from the integrated signature. Original `PValue`, `FDR` and `bonf` are retained when available. They are k-mer test values, not new miRNA-family P values.

KRS retains `target_frequency`, `target_rank_percentile` and `selected_KRS`. The plot shows target percentile, not a calibrated P value or fold enrichment.

Family rows sharing a target word are preserved. Multiple mature miRNAs within a family/seed are grouped into semicolon-separated IDs/accessions. This prevents treating each mature ID as an independent discovered site; family annotations sharing a word still share the same evidence. Plots label each family and word explicitly and do not sum scores across miRNAs.

Unmatched words remain in `query_kmers.tsv`; the annotation does not discard them silently. Supplied explicit words may be outside the selected signature; their scores are joined where available, and missing score values remain missing.

### Launch examples

The bundled table is selected automatically. To use another table, set `analysis.seed_family_file = "path/to/miR_Family_Info.txt"` before calling the method.

```python
# After k=7 counting and both KEA extraction branches:
mirna = analysis.AnnotateMiRNASeeds(
    kmers=signature,
    source="KEA",
    plot=True,
)

# After KRS(reference="my_target", ...):
mirna = analysis.AnnotateMiRNASeeds(
    source="KRS",
    reference="my_target",
    plot=True,
)

# Orientation check, independent of any analysis:
mirna = analysis.AnnotateMiRNASeeds(
    kmers=["GTCTTCC"], source="manual", plot=False,
)
print(mirna["matches"])
```

### Outputs

Returns `matches`, `query`, `seed_map`, `figures`, `settings`, `output_dir`; also stores `analysis.mirna_seed_annotation`.

- `seed_matches.tsv`: matched words, family/seed, mature miRNA IDs and original scores.
- `query_kmers.tsv`: all unique supplied words and matched-family count.
- `family_seed_map.tsv`: filtered family-to-target-word mapping.
- `settings.json`: inputs, site type, selected words and parameters.
- `seed_scores.png` / `seed_scores.pdf` when scored matches are plotted.

Figure objects remain available for the caller to display or close. No extra runtime dependencies are needed.

## K-CDS figures and BED compatibility

See [K-CDS](K-CDS) for the final panels, frame interpretation and compatible BED inputs.
