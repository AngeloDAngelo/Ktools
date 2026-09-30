# Methodology

[Documentation](Documentation.md) · [Tutorial](Tutorial.md) · [Benchmarks](Benchmarks.md)

## Counting and normalization

For a transcript of length L and k-mer length k, enumerate overlapping starts from 0 through L−k inclusive. There are L−k+1 possible windows when L≥k. Count each canonical k-mer without reverse-complement collapsing: RNA direction matters. Let Cᵤₜ denote the count of k-mer u in transcript t.

Per-transcript frequency is fᵤₜ = Cᵤₜ / Σᵥ Cᵥₜ. Pooled population frequency is Fᵤᴳ = Σₜ∈ᴳ Cᵤₜ / ΣᵥΣₜ∈ᴳ Cᵥₜ. The first normalizes each transcript independently; the second weights contributions by the number of counted windows.

The supplied source omitted the final window. This development snapshot corrects that boundary. Manuscript-reported benchmark values have not been recomputed with the correction; see [source changes](../CHANGELOG.md).

## KEA group-level selection

When the control population has at least twice the reference population's transcript count, sample reference-sized control subsets without replacement and average their pooled frequencies (`n_iter=1000`). Otherwise use the full control population. Sampling uses NumPy's global RNG; set `np.random.seed` before extraction and record the seed.

For each group, calculate the 25th percentile Q1, median M and 75th percentile Q3 across k-mer frequencies. Transform a frequency F as:

- (F−M)/(M−Q1) when F<M;
- (F−M)/(Q3−M) otherwise.

The **delta-median** score is the transformed reference frequency minus the transformed control frequency. Select the first and last `max(1, floor(number_of_kmers × sign_emp_pval))` rows after sorting. At k=7 and fraction 0.01, each tail has 163 entries before intersection. The parameter name `sign_emp_pval` denotes a selected tail fraction here, **not a calibrated P value**. Ties can make a boundary selection arbitrary. Zero quartile spreads make this scaling undefined; inspect outputs and avoid interpreting nonfinite scores.

The optional **delta-rank-percentile** selection subtracts frequency ranks normalized by the number of k-mer types. Its implementation uses `rank(method="first")`, so tied ranks depend on row order.

## KEA transcript-level selection

Compare fᵤₜ distributions with a two-sided Mann–Whitney U test (`stat_test="mwu"`). The implementation skips a k-mer if both groups have constant frequencies, even when the two constants differ. It also supports a two-sample KS test and Welch's t-test; these are not the manuscript's default procedure.

With `median=True`, log₂ fold change is:

`log2((median(reference frequencies) + 1e-9) / (median(control frequencies) + 1e-9))`.

A zero median in both groups produces log₂FC=0 even if a subset of transcripts carries a real difference. Selection combines the signed fold-change threshold with raw (`PValue`), BH-adjusted (`FDR`, selected with `correction_method="fdr"`) or Bonferroni-adjusted (`bonf`, selected with `correction_method="bonferroni"`) significance. Both corrected columns are stored in the table. Correction is across the tested k-mers; skipped k-mers are excluded.

The integrated signature is the set intersection of **delta-median** and **stat_log2fc** selections. Intersect the explicit result lists, including empty lists. The legacy Venn convenience function drops empty method sets and should not be used to compute this two-method signature.

## Continuous transcript scores

`KmerCorrelationWithScore` correlates normalized k-mer frequencies with a numeric transcript score, using Spearman by default. It selects tails of the correlation ranking and saves the coefficients. Its tail selection is not a multiple-testing-adjusted significance test. The supplied implementation requires `species="input"` for the ordinary combined FASTA table.

## KRS ranks

Normalize counts within each transcript. For each k-mer, exclude transcripts with zero frequency, rank remaining values with average ranks for ties, and divide by the number of nonzero observations. Include the target in this population. Retain target ranks ≥ `1 − top_pct/100`.

A k-mer present only in the target attains rank 1; this need not indicate many occurrences. Background choice, population size and prevalence therefore matter. No hypothesis-test P value is produced. The default `top_pct=5` corresponds to rank ≥0.95. The percentile threshold is an analysis parameter and should be reported explicitly.

## K-RBP associations

Decompose longer query k-mers into overlapping k-mers of the matrix length (`k_out=5` for the PEKA matrices described in the manuscript). Deduplicate, convert T to U when requested, and intersect with matrix rows. The default background is the matrix universe excluding the query. A user-defined background should also be disjoint; the current code does not remove overlaps automatically.

If the background exceeds ten times the query size, subsample it without replacement with a generator seeded at 42. This snapshot sorts set-derived lists before sampling for reproducibility across Python processes. Drop missing scores independently per profile and skip profiles with fewer than two query or background observations.

For each RBP profile, the one-sided Mann–Whitney U alternative is that query scores exceed background scores. Cliff's delta is:

`δ = [number of query > background pairs − number of query < background pairs] / (n_query × n_background)`.

It lies between −1 and 1. Positive values indicate larger query scores. **`zscore_diff` is the returned column containing Cliff's delta**; its historical name does not mean Z score. The separate `effect_size` column from the MWU calculation has the opposite sign. Use `zscore_diff` when applying manuscript Cliff's-delta thresholds.

BH and Bonferroni corrections are computed across tested RBP profiles. Profile-level results retain cell-line identity; do not silently merge them into protein-level tests. Apply positive effect-size and adjusted-P filters explicitly before biological interpretation.

## Positional analysis

The sliding-window method counts overlapping occurrences within complete windows, normalizes by possible starts, and can cluster k-mer spatial profiles. Window size, step and number of clusters are configurable. This positional analysis does not establish functional domains without additional evidence.

