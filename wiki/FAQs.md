# Frequently asked questions

[Documentation](Documentation) · [Tutorial](Tutorial)

### How do I choose between KEA and KRS?

Use KEA to compare two RNA populations. Use KRS to characterize one RNA relative to a population. Either signature can be passed to K-RBP.

### Is K-tools a web server?

K-tools is a Python library.

### Does it accept RNA sequences containing U?

Convert U to T and uppercase before FASTA counting. The enumerated count universe is A/C/G/T. K-RBP can convert T to U to match RNA-formatted score matrices.

### Does it count overlapping occurrences?

Yes, but the reference-class counter omits the final possible window: `AAAA`, k=2 yields two `AA` counts. Use sequences longer than k.

### Is a k=7 analysis always appropriate?

No. Choose k according to sequence length, data size and sparsity. A 4ᵏ universe grows exponentially. The bundled tutorial uses k=3 for a small, fast demonstration. Sparse k=7 populations can have zero quartile spread and undefined delta-median scores.

### Is `sign_emp_pval` a P-value cutoff?

For delta selections it is a tail fraction of sorted scores. Statistical significance is controlled separately by `fdr_alpha` and `correction_method` in `stat_log2fc`.

### Why is my integrated signature empty?

The group-level and transcript-level selections may not overlap. This is a valid outcome. Inspect the full tables, group design, sparsity and thresholds. Do not drop the empty method and report the other method's list as the intersection.

### Why is a k-mer missing from the statistical table?

The code skips k-mers for which both groups have constant normalized frequencies, even if the two constant values differ. This is an implementation limitation. If every k-mer is skipped, the current method can fail rather than return an empty table.

### Does KRS `top_pct=5` select 5% of k-mer types?

No. It selects k-mers for which the target ranks at or above the 95th percentile among transcripts containing that k-mer. Signature size varies. The target is included in the ranking.

### Is `zscore_diff` a Z score?

In K-RBP it stores Cliff's delta. Positive values indicate higher query scores. The `effect_size` field has the opposite sign under the MWU branch; use `zscore_diff` for the documented biological filters.

### Why is my K-RBP plot empty?

Check that profiles pass both the significance and effect-size filters. Set `fc_thresh` explicitly; the default 1.0 can exclude all points under a strict filter because Cliff's delta cannot exceed 1. A word cloud raises an error when no profiles pass. The full returned table is separate from plot selection.

### Is the RBP score matrix included?

Two PEKA matrices are included in [data/peka](https://github.com/AngeloDAngelo/Ktools/blob/main/data/peka/README.md). Choose one input file. The synthetic demo uses a fictional matrix.

### Are results deterministic?

Set NumPy's global seed before KEA control resampling. K-RBP uses seed 42 for background subsampling but set-derived ordering can change the sampled background between processes. Pin input ordering and package versions, and record the code commit.

### Do the advanced mapping helpers run on my computer?

BEDTools helpers require `bedtools` on PATH. Genome mapping also requires the external mapper script and resources configured with `KTOOLS_MAPPER_SCRIPT` and `KTOOLS_MAPPER_DIR`.

### What should I cite before the preprint is published?

Cite the repository URL and the exact release or commit used, with access date. The preprint DOI will be added when available.

### Is uploading the `wiki/` folder enough to enable GitHub Wiki?

No. GitHub Wiki is a separate Git repository. The folder contains prepared page sources; publish them through the Wiki editor or export them to a clone of the Wiki repository. The main `docs/` pages work immediately after uploading the repository.
