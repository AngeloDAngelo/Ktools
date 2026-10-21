# Co-occurrence — WORK IN PROGRESS

**Experimental, pending validation on biological datasets.** Technical tests
cover counting, randomization, distances and reproducibility. They do not
establish biological validity or suitability for publication conclusions.
Results record `development_status: work_in_progress` in `settings.json`.

The entire implementation lives in `KEA.cooccurrence` in `KEA.py`. Import KEA
as usual; no support module or additional import is required.

```python
from KEA import KEA

analysis = KEA("results", "ref.fa", "ctrl.fa")
# Supply the enriched k-mers already selected in your analysis.
cooccurrence_table, distance_table = analysis.cooccurrence(
    kmers=KEA_motifs,
    mode="cooccurrence_distance",
    n_resamples=999,
    n_reference=99,
    seed=42,
    plot=True,
)
display(cooccurrence_table)
display(distance_table)
print(analysis.cooccurrence_results["output_dir"])
```

`mode="cooccurrence"` returns one pandas DataFrame and tests all unordered distinct motif pairs using
`sum_t count(A,t) * count(B,t)`. It preserves global A/B occurrence totals and
randomly assigns distinct eligible k-mer windows across the reference FASTA,
weighting transcripts by their available windows. Fixing both counts within
each transcript would make this co-occurrence statistic constant.

`mode="cooccurrence_distance"` returns a tuple of two DataFrames:
`(cooccurrence_table, distance_table)`. It first runs the co-occurrence test,
retains pairs with `q_cooccurrence <= alpha`, and uses
every transcript containing both motifs. Within each such transcript it
preserves the A and B counts and randomizes their positions. It computes all
Cartesian A×B distances: one A and two B occurrences give two distances.
Distances are **`abs(start_A - start_B)`**, in nucleotides, using 0-based start
coordinates only. AB and BA are the same motif pair and appear only once.

The co-occurrence table reports motif identities, occurrence counts, number
of co-occurring transcripts, observed site pairs, analytical/random expected
pairs, random mean/median/sample SD and 2.5–97.5 percentiles, excess pairs,
observed/random enrichment, descriptive z-score, empirical p-value, BH
q-value and test status. Null percentiles are reference-distribution bounds,
not confidence intervals for the biological effect.

The distance table includes the identities and co-occurrence selection
statistics, number of transcripts and distances, observed and random
**mean, median, SD, minimum, maximum, Q25 and Q75**, the CDF-difference test
statistic (`distance_KS_statistic`), empirical `p_distance`, BH `q_distance`
and status. Observed summaries are computed for all pairs with occurrences,
even when unselected or non-testable. Random summaries and distance p/q are
missing for unselected/non-testable pairs. Sample SD uses ddof=1 and is
missing for a single distance. Random summaries describe the pooled
histogram from `n_reference` replicates. For observed distances 20 and
120 nt: mean=median=70, SD=70.7107, min=20, max=120.

For compatibility, `short` and `long` are accepted as aliases, as is the
spelling `cooccurence_distance`. The return contract is the new table/tuple
contract for all aliases.

Random windows have the same k and contain only A/C/G/T. Windows equal to
either motif of the current pair are excluded. Sampling is **without
replacement**, including between random A and B sites in a replicate. Sites
can recur across replicates. Windows at different start coordinates may
overlap. Observed overlapping occurrences are also counted. Scanning follows
the orientation of the supplied transcript FASTA; reverse complements are not
added. Motifs must have one common length (1–31 nt); U is normalized to T.

Long compares the observed distance CDF to an independent random reference
CDF. A maximum-CDF-difference statistic is calibrated by complete random
replicates, preserving dependence among distances that share sites; it does
not use an iid two-sample KS p-value. Short tests excess pairs; long tests
differences anywhere in the distance distribution. Empirical p-values use
`(1 + number of null statistics at least as extreme)/(n_resamples + 1)`.
BH correction is applied separately to the short family and to the selected
long family. Long q-values do not claim joint FDR control for both stages.

If a required pool lacks enough distinct random sites, that pair receives an
explicit non-testable status and missing p/q values. Long does not silently
discard transcripts with insufficient sites. `distance_filter=False` skips
short selection for diagnostics or an explicitly different analysis.

With `kmers=None`, selection is the strict intersection of stored enriched
sets for `delta_median` and `stat_log2fc` (configurable via `methods`), under
`self.KEA_results[species or "input"]`. Empty sets remain empty. Supplying
`kmers` reuses exactly your selected motifs without rerunning enrichment.
The random pool uses the reference FASTA; the control FASTA is not used in
this stage. Associations are relative to this uniform-site null, which does
not automatically control GC, RNA structure or expression.

Each call saves a fresh `<dir_out>/Cooccurrence/<species>/<mode>/run_*`
directory. Returned tables are also saved as `cooccurrence.tsv` and, in
distance mode, `cooccurrence_distance.tsv`. Outputs include the full
`pair_tests.tsv`, occurrence coordinates and
transcript counts, exact distance histograms/CDFs in NPZ, and parameters,
seed, FASTA hash and development status in JSON. `store_null=True` saves null
statistics. `save=False` keeps results in memory. `plot=True` displays up to
`top_n` figures (default 6). Results are also stored in
`self.cooccurrence_results`, whose `cooccurrence` and `distances` keys contain
the returned tables. DataFrame `.attrs` records the output directory, mode
and WORK IN PROGRESS status.

**Monte Carlo resolution matters:** B=999 permits a minimum p of 0.001.
For 163 motifs there are 13,203 hypotheses; isolated signals may require many
more replicates to reach BH significance. The method never silently reduces
the number of pairs or replicates. Full long analyses can be expensive.
