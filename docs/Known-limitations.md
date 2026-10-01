# Known limitations

[Documentation](Documentation.md) · [Methodology](Methodology.md)

The FASTA tutorial covers a tested subset of the supplied implementation. This repository does not claim that every inherited helper has been validated.

- **Scientific version matching:** counting now includes the final valid window. Stable K-RBP background ordering may also change selected background members. Paper results need rechecking.
- **Sparse inputs:** delta-median scaling divides by quartile distances without a zero-spread guard. Inspect for NaN/Inf. All-skipped transcript tests can cause an empty-DataFrame error.
- **Input handling:** alphabet normalization, duplicate FASTA IDs, empty groups and short/ambiguous sequences need explicit validation by the caller. The demo checks its inputs.
- **Constant distributions:** the test routine skips a k-mer when both groups are constant, even if constants differ.
- **Intersections:** `IntersectEnrichedKmersAndPlotVenn` omits empty method selections and includes all populated methods. It is not a safe substitute for the explicit two-method intersection shown in the tutorial; KRS/correlation results also use different structures.
- **K-RBP:** only continuous mode is implemented reliably in the supplied routine. The signature advertises binary-mode parameters, but the function has no implemented binary analysis branch. Empty query lists with `k_out` fail before validation. Custom background/query overlaps are not removed automatically.
- **Effect fields:** `zscore_diff` is Cliff's delta; `effect_size` has the opposite sign for MWU. Plot thresholds and return-table filtering are separate.
- **Word clouds:** the inherited code takes top rows from results sorted by FDR before computing display weights; the printed claim that it selects “top by Cliff's delta” is misleading. Some manuscript descriptions refer to alternate ranking rules. Explicitly select profiles for final figures in a versioned analysis notebook.
- **Correlation:** supply `species="input"` and initialize `KEA_results` to a dictionary if no previous extraction has populated it.
- **External helpers:** interval intersection contains a hard-coded BEDTools executable path. Genome mapping references laboratory paths and scripts absent here. Some helpers use shell commands and remove temporary directories; do not invoke them with valuable existing directories. These paths were not rewritten speculatively.
- **Cross-species and plotting helpers:** not covered by the tutorial's validation. Their availability in the API inventory does not establish portability or biological validity.

See [CHANGELOG](../CHANGELOG.md) for the limited changes made to the original code and [Data and reproducibility](Data-and-reproducibility.md) for manuscript-specific open items.

## Kmap and CDS domains

Kmap catches metric errors and returns zero, reuses cached external tables by key without checking changed input files, and breaks tied Sankey ranks by row order. Validate inputs as described in [Kmap](Kmap.md).

Domain tests use dependent occurrence counts and are exploratory. Genomic UniProt overlaps are not isoform-specific domain assignments. Matching FASTA/GTF/BED resources are required; missing/invalid models are skipped or raise errors. See [Domain enrichment](Domain-enrichment.md).
