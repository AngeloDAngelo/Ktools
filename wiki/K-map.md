# K-map: transcript k-mer profile comparison

[Documentation](Documentation) · [Tutorial](Tutorial) · [Module parameters](Module-parameters)

K-map compares a target RNA's complete k-mer frequency profile with candidate transcript profiles. The Python method is **`CompareTranscriptKmerProfiles`**. It is independent of KEA/KRS signature extraction and returns ranked similarity scores, not significance tests.

![K-map: compare the target RNA k-mer frequency profile with matching and low-similarity profiles](https://github.com/AngeloDAngelo/Ktools/raw/refs/heads/main/docs/assets/kmap-workflow.png)

## Compare profiles

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

# Alternative external candidate FASTA:
# similarity = analysis.CompareTranscriptKmerProfiles(
#     "reference_006", "candidate_transcripts.fa",
#     metric="cosine", fasta_input=True, k=3,
#     fasta_species_key="candidate_set_k3",
# )
```

## Sankey-style rank alignment

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
