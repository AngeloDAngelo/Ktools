# K-CDS

[Documentation](Documentation.md) · [K-map](K-map.md) · [API reference](API-reference.md)

`K_CDS` asks whether occurrences of selected sequence words overlap genomic UniProt domain/feature annotations more often in reference transcripts than in a chosen background. This is an exploratory occurrence-level association, not a test of independent biological replicates or isoform-specific domain function.

The repository reserves [data/uniprot](../data/uniprot/README.md) for BED resources uploaded by the authors. Pass that folder as `uniprot_dir`; the files are not downloaded automatically.

## Required resources

Supply `gtf_file` and `uniprot_dir` explicitly; there are no bundled annotation datasets or laboratory-path defaults.

- Full spliced reference/control transcript FASTAs, in transcript orientation, including UTRs. CDS-only FASTAs will generally fail the length check.
- A matching GRCh38 GTF containing exon and CDS features and numeric CDS phases. Use the same transcript release as the FASTAs.
- A directory containing the selected UniProt genomic BED files, for example `unipDomain.bed`. Coordinate assembly and chromosome names must agree with the GTF projection.

FASTA and GTF transcript IDs have terminal version suffixes stripped. Duplicate IDs after stripping are rejected; control/reference sets must remain disjoint. FASTAs and GTF can be gzip-compressed. UniProt BED files are read as uncompressed text.

GTF coordinates are converted from 1-based inclusive to 0-based half-open. GTF chromosomes are prefixed with `chr` if needed; `chrMT` becomes `chrM`. BED chromosome names are used as supplied. Domains are assigned on the same strand with at least one nucleotide overlap.

## Parameters

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
| `seed` | 42 | NumPy generator seed for shuffled positions, including the independent frame-plot background. |
| `output_dir` | None | None creates a new `K_CDS_*` directory under `dir_out`. An explicit destination must not already exist, even if empty. |
| `plot` | True | Generate and save the domain heatmap and CODING/SHUFFLE frame plot; False still calculates and saves all tables. |
| `max_domains` | 60 | Maximum heatmap columns, prioritized by each domain's smallest FDR. None shows all. This does not reduce the tested domains or saved matrices. |

### Supported tracks

`unipDomain`, `unipInterest`, `unipStruct`, `unipLocTransMemb`, `unipLocExtra`, `unipLocSignal`, `unipLocCytopl`, `unipRepeat`.

For `unipDomain` and `unipInterest`, labels use BED column 27 when available, otherwise column 4. `unipInterest` retains only `Disordered`; `unipLocSignal` retains only `Signal peptide`. Other tracks use column 4. A generic BED with a different label convention must be adapted before use.

### Frame and shuffled background

The codon anchor is the start of the spliced CDS plus its first phase. Frame is `(occurrence_start - anchor) % 3`. Frame filtering does not translate the word and does not test amino-acid motifs.

Shuffle mode relocates intervals, not sequence words: the random interval need not contain the original k-mer. Positions are sampled independently and can repeat or overlap real hits. With frame filtering, shuffled starts also respect the codon anchor. Control FASTA is not read in shuffle mode.

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

## Outputs and quality control

The returned dictionary is also stored as `analysis.k_cds` (also `analysis.domain_enrichment` for compatibility).

- `statistics.tsv`: A/B/C/D, raw odds ratio, corrected display log2 OR, raw P, BH FDR and selection flags.
- `occurrences.tsv`: transcript, word, spliced start/end, frame, group and overlapping labels. Shuffle-control rows are random intervals, not observed sequence matches.
- `qc.tsv`: transcript status.
- `log2_or.tsv` and `fdr.tsv`: full word-by-feature matrices.
- `settings.json`: input paths, selected words and analysis/plot settings.
- `heatmap.png` and `heatmap.pdf` when `plot=True`.

QC exclusions include missing/noncoding GTF models, missing exons, FASTA/GTF length mismatch, noncontiguous spliced CDS and inconsistent phase. Overlapping exons, conflicting loci, invalid BED blocks or CDS outside an exon raise errors. GTF CDS does not include the stop codon.

An empty eligible occurrence set or absence of same-strand overlaps raises an error. Inspect assembly, IDs, full-transcript sequences and tracks before changing thresholds.

The returned `figure` can be displayed or closed by the caller. Plot generation saves figures but does not call `plt.show()`.

## CDS frame distribution

K-CDS is called as `analysis.K_CDS(...)` in Python. The Python method is `K_CDS`; the former `DomainEnrichment` name is no longer exposed.

For each k-mer, two horizontal stacked bars compare **CODING** and **SHUFFLE**. Blue, orange and green show the percentage of starts at the first, second and third codon positions (frames 0, 1 and 2), respecting the GTF phase. CODING pools all CDS-contained reference occurrences, before `in_frame` filtering. SHUFFLE draws one unrestricted random CDS start per reference occurrence, in the same transcript and with the same interval length. Its independent seeded generator does not alter the domain-enrichment background. Even with `in_frame=True`, this plot includes all three frames. The shuffle is generated for this plot in both background modes.

Percentages sum to 100 within each k-mer/group with occurrences. Missing groups have NaN percentages and are marked “No occurrences”. Short CDSs and finite sampling can produce unequal shuffle proportions; no frame significance test is performed.

- `frame_distribution.tsv`: k-mer, group (`coding`, `shuffle`, or `control`), frame, count and percent; control counts are included when control transcripts are loaded.
- `frame_distribution.png` and `frame_distribution.pdf`: CODING/SHUFFLE plot when `plot=True`.
- `result["frame_distribution"]`: the count/percentage table.
- `result["frame_figure"]`: the Matplotlib figure (`None` with `plot=False`).

This module still requires GTF and UniProt resources and an evaluable domain overlap; it is not a standalone frame-only analysis.

## What makes another BED compatible?

Compatible files are not limited to the original downloads. They must satisfy all of these requirements:

- **Format:** uncompressed, tab-separated genomic BED text with at least six columns: `chrom`, `chromStart`, `chromEnd`, `name`, `score`, `strand`. Start is zero-based and end is excluded, so `chr1\t100\t110\tFeature\t0\t+` represents ten nucleotides. Intervals must have start < end and strand `+` or `-`; unknown strand will not match the GTF. With BED12 and `split_blocks=True`, column 10 is blockCount, column 11 contains comma-separated blockSizes, and column 12 contains comma-separated blockStarts relative to chromStart. Counts must agree and all blocks must lie within the outer interval. Files use the selected track name, e.g. `unipLocCytopl.bed`.
- **Labels:** feature names come from column 4, except `unipDomain` and `unipInterest`, which use column 27 when present. Names must be nonempty. `unipLocSignal` retains only the exact, case-sensitive label `Signal peptide`; `unipInterest` retains only `Disordered`. Other tracks retain their labels as supplied. With multiple tracks, labels are prefixed with the track name. A different label convention must be adapted before running the method.
- **Assembly:** genomic coordinates must use the same reference assembly as the GTF (this workflow uses human GRCh38/hg38). Renaming an hg19 file does not convert its coordinates to hg38. BED chromosome names must match the projected GTF names: `chr1`, `chr2`, …, `chrX`, `chrY`, `chrM`. The method prefixes GTF names with `chr` and converts GTF `chrMT` to `chrM`, but it does not normalize BED chromosome names or lift coordinates between assemblies. Full-transcript FASTA sequences and IDs must also match the GTF annotation release.

Upload your working BED files to `data/uniprot`, preserving their names and annotations, and record source, assembly, release and download date.

## Separate figures

With `plot=True`, the domain heatmap and CODING/SHUFFLE frame distribution are separate figures. Display `result["figure"]` for the heatmap and `result["frame_figure"]` for the frame plot. The frame legend is outside the SHUFFLE panel. Figures are saved as `heatmap.png`/`.pdf` and `frame_distribution.png`/`.pdf`.
