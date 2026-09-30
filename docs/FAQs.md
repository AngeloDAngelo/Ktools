# Frequently asked questions

[Documentation](Documentation.md) · [Tutorial](Tutorial.md)

### How do I choose between KEA and KRS?

Use KEA to compare two RNA populations. Use KRS to characterize one RNA relative to a population. Either signature can be passed to K-RBP.

### Is K-tools a web server?

This distribution is a Python library. The catRAPID web pages linked as design references describe another tool; their server limits, scores and retention policies do not apply to K-tools.

### Does it accept RNA sequences containing U?

Convert U to T and uppercase before FASTA counting. The enumerated count universe is A/C/G/T. K-RBP can convert T to U to match RNA-formatted score matrices.

### Does it count overlapping occurrences?

Yes. `AAAA`, k=2 contains three `AA` occurrences. The supplied source missed the final window; this repository corrects it. Rerun the paper analyses before attaching the benchmark claims to this corrected software release.

### Is a k=7 analysis always appropriate?

No. The manuscript commonly uses k=7, but data size and sparsity matter. A 4ᵏ universe grows exponentially. The bundled tutorial uses k=3 for a small, fast demonstration. Sparse k=7 populations can have zero quartile spread and undefined delta-median scores.

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

The real matrix used in the paper is not in the supplied source folder. The demo generates a fictional matrix for testing. Obtain the PEKA resource and record the exact preprocessing and profile selection before biological analysis.

### Are results deterministic?

Set NumPy's global seed before KEA control resampling. K-RBP uses seed 42 for background subsampling and this snapshot sorts set-derived lists before sampling. Pin input ordering and package versions, and record the code commit.

### Can I reproduce the paper figures from this repository alone?

Not yet. The available material contains code, manuscript text and figure artwork, but not the raw benchmark inputs, notebooks or curated matrices. The repository provides methods, reported summaries and a runnable demonstration. See the [reproducibility inventory](Data-and-reproducibility.md).

### Do the advanced mapping helpers run on my computer?

Several retain laboratory-specific filesystem paths and external scripts. They are outside the portable tutorial and require configuration and separate validation.

### What should I cite before the preprint is published?

Cite the repository URL and the exact release or commit used, with access date. Add the approved author list and preprint DOI when available. Do not cite an invented bioRxiv record.

### Is uploading the `wiki/` folder enough to enable GitHub Wiki?

No. GitHub Wiki is a separate Git repository. The folder contains prepared page sources; publish them through the Wiki editor or export them to a clone of the Wiki repository. The main `docs/` pages work immediately after uploading the repository.
