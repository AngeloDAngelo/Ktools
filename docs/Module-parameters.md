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

## Kmap: transcript similarity and Sankey visualization

Kmap compares a target RNA's complete k-mer frequency profile with candidate transcript profiles. The Python method is **`CompareTranscriptKmerProfiles`**. It is independent of KEA/KRS signature extraction and returns ranked similarity scores, not significance tests.

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

The title's similarity is `1 - mean(abs(query_rank - comparison_rank))`, calculated for the plotted words. It is not the selected Kmap correlation/cosine/JSD score and is not a P value. Filtering changes its interpretation.

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


## DomainEnrichment: CDS occurrences and UniProt features

`DomainEnrichment` asks whether occurrences of selected sequence words overlap genomic UniProt domain/feature annotations more often in reference transcripts than in a chosen background. This is an exploratory occurrence-level association, not a test of independent biological replicates or isoform-specific domain function.

The repository reserves [data/uniprot](../data/uniprot/README.md) for BED resources uploaded by the authors. Pass that folder as `uniprot_dir`; the files are not downloaded automatically.

### Required resources

Supply `gtf_file` and `uniprot_dir` explicitly; there are no bundled annotation datasets or laboratory-path defaults.

- Full spliced reference/control transcript FASTAs, in transcript orientation, including UTRs. CDS-only FASTAs will generally fail the length check.
- A matching GRCh38 GTF containing exon and CDS features and numeric CDS phases. Use the same transcript release as the FASTAs.
- A directory containing the selected UniProt genomic BED files, for example `unipDomain.bed`. Coordinate assembly and chromosome names must agree with the GTF projection.

FASTA and GTF transcript IDs have terminal version suffixes stripped. Duplicate IDs after stripping are rejected; control/reference sets must remain disjoint. FASTAs and GTF can be gzip-compressed. UniProt BED files are read as uncompressed text.

GTF coordinates are converted from 1-based inclusive to 0-based half-open. GTF chromosomes are prefixed with `chr` if needed; `chrMT` becomes `chrM`. BED chromosome names are used as supplied. Domains are assigned on the same strand with at least one nucleotide overlap.

### Parameters

| Parameter | Default | Exact behavior |
| --- | --- | --- |
| `ref_fasta`, `ctrl_fasta` | None | Optional FASTA overrides for this call; None uses paths from the KEA constructor. The object is not modified. Control is ignored in shuffle mode. |
| `gtf_file` | None | Required path to the matching GTF; missing path argument raises an error. |
| `uniprot_dir` | None | Required directory of `<track>.bed` files. |
| `kmers` | None | Explicit word list or a single word. If omitted, use `KEA_results[species][method]["enriched"]`; this is one method's list, not an automatically intersected KEA signature. Words are uppercased, U becomes T, duplicates are removed, and ambiguous words are rejected. Mixed lengths are allowed. |
| `method` | `stat_log2fc` | Extraction-result key used only when `kmers=None`. Pass an explicit signature to analyze KEA intersections or KRS results. |
| `species` | None | Only None or `input` is supported. Cross-species keys are rejected. |
| `in_frame` | False | False retains all CDS-contained occurrences. True requires the occurrence start at a codon boundary determined from the first CDS phase. Word length need not be divisible by three. |
| `background` | `control` | `control` finds the same words in control transcript CDSs. `shuffle` draws one random eligible CDS position per reference occurrence, within the same transcript and of the same length. |
| `tracks` | `("unipDomain",)` | Selected UniProt feature files, listed below. With multiple tracks, labels receive track prefixes to preserve their origin. |
| `alpha` | 0.05 | BH-adjusted significance threshold, inclusive `fdr <= alpha`; must be between 0 and 1. |
| `fdr_scope` | `global` | `global`: BH across every testable word/feature cell. `kmer`: separate BH correction within each word. These define different testing families. |
| `split_blocks` | True | For BED12, use its blocks instead of the complete outer interval, avoiding assignment through intervening gaps. False uses the outer interval. |
| `seed` | 42 | NumPy generator seed for shuffled positions. It has no statistical effect in control mode. |
| `output_dir` | None | None creates a new `DomainEnrichment_*` directory under `dir_out`. An explicit destination must not already exist, even if empty. |
| `plot` | True | Generate and save the heatmap; False still calculates and saves all tables. |
| `max_domains` | 60 | Maximum heatmap columns, prioritized by each domain's smallest FDR. None shows all. This does not reduce the tested domains or saved matrices. |

#### Supported tracks

`unipDomain`, `unipInterest`, `unipStruct`, `unipLocTransMemb`, `unipLocExtra`, `unipLocSignal`, `unipLocCytopl`, `unipRepeat`.

For `unipDomain` and `unipInterest`, labels use BED column 27 when available, otherwise column 4. `unipInterest` retains only `Disordered`; `unipLocSignal` retains only `Signal peptide`. Other tracks use column 4. A generic BED with a different label convention must be adapted before use.

#### Frame and shuffled background

The codon anchor is the start of the spliced CDS plus its first phase. Frame is `(occurrence_start - anchor) % 3`. Frame filtering does not translate the word and does not test amino-acid motifs.

Shuffle mode relocates intervals, not sequence words: the random interval need not contain the original k-mer. Positions are sampled independently and can repeat or overlap real hits. With frame filtering, shuffled starts also respect the codon anchor. Control FASTA is not read in shuffle mode.

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
result = analysis.DomainEnrichment(
    gtf_file="resources/annotation.gtf",
    uniprot_dir="data/uniprot",
    ref_fasta="reference_full_transcripts.fa",
    ctrl_fasta="control_full_transcripts.fa",
    kmers=signature,
    background="control",
    in_frame=False,
    tracks=("unipDomain",),
    alpha=0.05,
    fdr_scope="global",
    split_blocks=True,
    plot=True,
    max_domains=60,
)
print(result["output_dir"])
print(result["statistics"])
print(result["qc"]["status"].value_counts())
```

For the position-matched null, change `background="shuffle"` and specify `seed=42`. If you pass `output_dir`, choose a new directory for each call.

### Outputs and quality control

The returned dictionary is also stored as `analysis.domain_enrichment`.

- `statistics.tsv`: A/B/C/D, raw odds ratio, corrected display log2 OR, raw P, BH FDR and selection flags.
- `occurrences.tsv`: transcript, word, spliced start/end, frame, group and overlapping labels. Shuffle-control rows are random intervals, not observed sequence matches.
- `qc.tsv`: transcript status.
- `log2_or.tsv` and `fdr.tsv`: full word-by-feature matrices.
- `settings.json`: input paths, selected words and analysis/plot settings.
- `heatmap.png` and `heatmap.pdf` when `plot=True`.

QC exclusions include missing/noncoding GTF models, missing exons, FASTA/GTF length mismatch, noncontiguous spliced CDS and inconsistent phase. Overlapping exons, conflicting loci, invalid BED blocks or CDS outside an exon raise errors. GTF CDS does not include the stop codon.

An empty eligible occurrence set or absence of same-strand overlaps raises an error. Inspect assembly, IDs, full-transcript sequences and tracks before changing thresholds.

The returned `figure` can be displayed or closed by the caller. Plot generation saves figures but does not call `plt.show()`.
