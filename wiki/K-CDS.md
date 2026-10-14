# K-CDS

[Documentation](Documentation) · [K-map](K-map) · [API reference](API-reference)

`analysis.K_CDS(...)` performs the occurrence-level UniProt association analysis and, by default, directly creates the two final panels: domain enrichment and CDS positional enrichment with coding/shuffle frame bars and In Frame Index. No subsequent notebook cells or external R script preparation are needed. Occurrence-level Fisher tests remain exploratory; genomic overlaps do not establish isoform-specific protein function.

## Required resources and software

Supply a matching GRCh38 GTF with exon/CDS features and numeric phases, a directory of the selected UniProt BED tracks, and full spliced reference/control transcript FASTAs in transcript orientation. UTRs must be retained. GTF/FASTA IDs are matched after stripping terminal version suffixes for the CDS projection; duplicate normalized IDs are rejected. GTF and FASTA inputs may be gzip-compressed. UniProt inputs are uncompressed genomic BED files. GTF coordinates are converted to zero-based half-open coordinates and chromosome names prefixed with `chr`; `chrMT` becomes `chrM`. BED names are used unchanged.

The final panels require `Rscript` on PATH, R packages **ComplexHeatmap** and **circlize**, and a graphics-capable R installation. Install these separately from the Python package:

```r
install.packages(c("BiocManager", "circlize"))
BiocManager::install("ComplexHeatmap")
```

BEDTools on PATH is needed only when missing positional profiles must be computed through `EnrichmentAnalysis`. It is not used when complete cached profiles are available. The embedded R script is written to the output directory and run by Python. Dependency failures stop before the expensive positional analysis. `plot=False` disables rendering by default and requires neither R nor positional inputs; domain/frame tables are still calculated.

## Parameters

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

## Positional input reuse and automatic calculation

With `recompute_positional=False`, the method first looks for per-word `EnrichmentAnalysis` entries in `self.enrichment`, identified by exact `<kmer>_Enr_VS_Ctrl` analysis labels. Repeated stored entries use the last row per bin. Pooled signature analyses such as `ZipCodes_Enr_VS_Ctrl` cannot provide individual word profiles. If no reusable per-word entries exist, the previous `self.k_cds` positional table/output directory and then `dir_out/K_CDS/positional_enrichment.tsv` and `dir_out/positional_enrichment.tsv` are considered. An explicit `positional_enrichment` overrides discovery.

If profiles are missing, the method loads both full-transcript FASTAs and calls `EnrichmentAnalysis` separately for missing words. It uses `self.metatx` when available; otherwise it looks for one compatible `metatx*.bed` file under `dir_out`, or reconstructs a metatranscript from the matching FASTA/GTF using mean regional length proportions. Reconstructed regional bin counts are three times rounded percentage proportions, with a minimum of one bin per region; zero-width transcript bins are omitted. The reconstructed coordinates can differ from a separately supplied annotation table: retain the original metatranscript/cache to reproduce an earlier analysis exactly.

Occurrence coordinates come from `bed_feature`, then a saved `dir_out/kmers_bed.bed`, or an automatic scan for overlapping exact matches in both FASTA populations. Automatic positional counting uses the existing KEA helper's **k=7** convention; for other word lengths, supply complete positional estimates. All computed and cached profiles must share the same CDS bin set; incompatible partial caches are recalculated. `self.enrichment` is restored after internal calls, and the assembled profiles are stored in the returned result and saved for reuse. Original FASTA IDs are used for transcript-coordinate enrichment and must match `self.metatx`/BED IDs.

## Clusters, ordering and panels

CDS bins are sorted numerically. Complete unmasked finite log2 OR estimates are clipped to ±`color_cap`, then each row is centered on its own mean. Euclidean distance and Ward.D2 are applied to these profiles and the tree is cut into two groups. The cluster with the larger mean contrast between the first `start_fraction` of bins and the remaining bins is called **CDS start**; the other is **CDS body**. The contrast uses clipped, uncentered values. These are relative labels, not independent significance tests, and body rows may still have positive initial contrast. The numbers of words are data-dependent and are not fixed at 22/52.

Both panels use body first, then start, with decreasing **In Frame Index** inside each group and DNA sequence as a tie-breaker. In Frame Index is 100 times the largest coding-frame fraction. The positional heatmap displays clipped values without centering. Its grey/white/black bars represent first/second/third starting codon positions for coding and shuffle groups. Missing frame groups have unavailable fractions and explicit missing-data labels; finite groups sum to 100%. K-mer labels use U. Domain rows below `frame_threshold` are marked `N.C.`; stars mark FDR ≤`alpha` in either direction.

Frame 0/1/2 is `(start - anchor) % 3`, with anchor equal to CDS start plus first GTF phase. Domain shuffle relocates intervals; it does not shuffle nucleotides or require randomized intervals to contain the word. Random starts can repeat. The frame-distribution background draws one unrestricted random start per reference CDS occurrence, independently of domain filtering and background mode. No frame significance test is performed.

## Supported tracks

`unipDomain`, `unipInterest`, `unipStruct`, `unipLocTransMemb`, `unipLocExtra`, `unipLocSignal`, `unipLocCytopl`, `unipRepeat`. For `unipDomain` and `unipInterest`, use BED column 27 when present, otherwise column 4. `unipInterest` retains only `Disordered`; `unipLocSignal` only `Signal peptide`. Other tracks use column 4.

## Statistical interpretation

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

## Example

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

## Outputs and quality control

The returned dictionary is stored as both `analysis.k_cds` and `analysis.domain_enrichment`.

- `statistics`, `occurrences`, `qc`, `log2_or`, `fdr`, `frame_distribution`, `settings`, `output_dir`: original tables and metadata, saved as their matching TSV/JSON files.
- `positional_enrichment`: assembled positional DataFrame; also `positional_enrichment.tsv` when final panels are enabled.
- `cluster_assignments`: per-word DataFrame with group, maximum frame fraction, preferred frame, capped/raw initial contrasts, In Frame Index and RNA label. None without final panels.
- `final_figures`: dictionaries of `pdf`, `png`, `svg` paths for `main_domains` and `supplement_positions_frames`; empty without final panels.
- `figures_dir`: final-panel directory, present when rendering succeeds. Contains `cluster_assignments.tsv`, both panels in all three formats, and `plot_settings.txt` with the R session information.
- `heatmap_script.R` and `heatmap_R.log`: generated renderer and its stdout/stderr in the run directory.
- `figure` and `frame_figure`: compatibility keys set to None. The old `heatmap.png/pdf` and `frame_distribution.png/pdf` are no longer generated.

QC exclusions include missing/noncoding models, missing exons, sequence/GTF length mismatch, noncontiguous CDS and inconsistent phase. Conflicting loci, overlapping exons, invalid BED blocks and CDS outside an exon raise errors. Empty eligible occurrence sets or absent same-strand UniProt overlaps also raise errors. GTF CDS excludes the stop codon. On an R rendering failure the existing domain/frame tables remain saved, with the R log for diagnosis.

## What makes another BED compatible?

Compatible files are not limited to the original downloads. They must satisfy all of these requirements:

- **Format:** uncompressed, tab-separated genomic BED text with at least six columns: `chrom`, `chromStart`, `chromEnd`, `name`, `score`, `strand`. Start is zero-based and end is excluded, so `chr1\t100\t110\tFeature\t0\t+` represents ten nucleotides. Intervals must have start < end and strand `+` or `-`; unknown strand will not match the GTF. With BED12 and `split_blocks=True`, column 10 is blockCount, column 11 contains comma-separated blockSizes, and column 12 contains comma-separated blockStarts relative to chromStart. Counts must agree and all blocks must lie within the outer interval. Files use the selected track name, e.g. `unipLocCytopl.bed`.
- **Labels:** feature names come from column 4, except `unipDomain` and `unipInterest`, which use column 27 when present. Names must be nonempty. `unipLocSignal` retains only the exact, case-sensitive label `Signal peptide`; `unipInterest` retains only `Disordered`. Other tracks retain their labels as supplied. With multiple tracks, labels are prefixed with the track name. A different label convention must be adapted before running the method.
- **Assembly:** genomic coordinates must use the same reference assembly as the GTF (this workflow uses human GRCh38/hg38). Renaming an hg19 file does not convert its coordinates to hg38. BED chromosome names must match the projected GTF names: `chr1`, `chr2`, …, `chrX`, `chrY`, `chrM`. The method prefixes GTF names with `chr` and converts GTF `chrMT` to `chrM`, but it does not normalize BED chromosome names or lift coordinates between assemblies. Full-transcript FASTA sequences and IDs must also match the GTF annotation release.

Upload your working BED files to `data/uniprot`, preserving their names and annotations, and record source, assembly, release and download date.
