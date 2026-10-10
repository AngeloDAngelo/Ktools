# Known limitations

[Documentation](Documentation) · [Methodology](Methodology)

The tutorial exercises the core FASTA workflow. The following constraints matter when changing inputs or using additional utilities.

- **Counting boundary:** the reference-class counter omits the last possible window. Sequences of length k yield no counts; effects are proportionally larger for short sequences.
- **Sparse inputs:** delta-median scaling divides by quartile distances without a zero-spread guard. Inspect for NaN/Inf. All-skipped transcript tests can cause an empty-DataFrame error.
- **Input handling:** alphabet normalization, duplicate FASTA IDs, empty groups and short/ambiguous sequences need explicit validation by the caller. The demo checks its inputs.
- **Constant distributions:** the test routine skips a k-mer when both groups are constant, even if constants differ.
- **Intersections:** `IntersectEnrichedKmersAndPlotVenn` omits empty method selections and includes all populated methods. It is not a safe substitute for the explicit two-method intersection shown in the tutorial; KRS/correlation results also use different structures.
- **K-RBP:** only continuous mode is implemented reliably in the supplied routine. The signature advertises binary-mode parameters, but the function has no implemented binary analysis branch. Empty query lists with `k_out` fail before validation. Custom background/query overlaps are not removed automatically.
- **Effect fields:** `zscore_diff` is Cliff's delta; `effect_size` has the opposite sign for MWU. Plot thresholds and return-table filtering are separate.
- **Word clouds:** `wordcloud_top_n` limits passing profiles in FDR order; `wordcloud_score` controls displayed weights. See [parameters](Module-parameters).
- **Correlation:** supply `species="input"` and initialize `KEA_results` to a dictionary if no previous extraction has populated it.
- **External helpers:** BEDTools must be on PATH. Genome mapping requires separately configured mapper resources and an external script.
- **Cross-species and plotting helpers:** not covered by the tutorial's validation. Their availability in the API inventory does not establish portability or biological validity.

## Kmap and CDS domains

Kmap catches metric errors and returns zero, reuses cached external tables by key without checking changed input files, and breaks tied Sankey ranks by row order. Validate inputs as described in [Kmap](Kmap).

Domain tests use dependent occurrence counts and are exploratory. Genomic UniProt overlaps are not isoform-specific domain assignments. Matching FASTA/GTF/BED resources are required; missing/invalid models are skipped or raise errors. See [K-CDS](K-CDS).
