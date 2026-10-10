# Documentation

[Home](Home) · [Tutorial](Tutorial) · [FAQs](FAQs)

K-tools links sequence enrichment to candidate post-transcriptional regulators. Select RNA populations from your experiment, count overlapping k-mers, extract a signature, and compare that signature with RBP-associated sequence preferences.

## Start here

1. [Install the package and its dependencies](Installation).
2. Read [the parameter guide](Module-parameters): exact thresholds, numerical examples, background construction and output interpretation.
3. Run [the simulated tutorial](Tutorial) and compare the included expected output.
4. Use [the API reference](API-reference) for complete call signatures and [methodology](Methodology) for the statistical definitions.

## Modules

### KEA: comparisons between RNA populations

Use KEA when you have a reference group (for example, RNAs sharing an experimental property) and a meaningful control group. At k=7 there are 4⁷ = 16,384 sequence types. Two complementary analyses compare pooled group frequencies and per-transcript frequency distributions. Their intersection defines the integrated signature.

The code also exposes a rank-difference approach and a continuous score-correlation analysis. These are distinct selections, not interchangeable statistical significance estimates. [Methodology](Methodology) describes their definitions.

### KRS: one RNA against a population

KRS ranks the target's normalized frequency separately for each k-mer among transcripts containing that k-mer. A k-mer at percentile ≥0.95 enters the default signature. The target participates in the ranking. This is not a selection of the top 5% of k-mer types within the target.

![KRS module schematic](https://github.com/AngeloDAngelo/Ktools/raw/refs/heads/main/docs/assets/krs-workflow.png)

### K-RBP: candidate RNA-binding proteins

K-RBP compares signature and background scores for each RBP–cell-line profile. A one-sided Mann–Whitney U test asks whether signature scores are higher; multiple-testing corrections and Cliff's delta support interpretation. When working from 7-mers against a 5-mer matrix, overlapping 5-mers are deduplicated before analysis.

![K-RBP module schematic](https://github.com/AngeloDAngelo/Ktools/raw/refs/heads/main/docs/assets/krbp-workflow.png)

RBP association depends on the reference matrix and experimental context. It is not a calibrated interaction probability, proof of direct binding, or evidence that other RBPs are absent.

## Inputs and resources

| Input | Format | Needed for |
| --- | --- | --- |
| Reference/control sequences | FASTA, uppercase DNA alphabet A/C/G/T, transcript orientation | KEA |
| Target plus background sequences | FASTA, unique identifiers; target included exactly once | KRS |
| RBP reference matrix | Numeric TSV, k-mers as rows, RBP profiles as columns | Biological K-RBP |
| Transcript scores | Two-column DataFrame: transcript identifier and numeric score | Correlation analysis |
| Transcript annotation or FASTA | Annotation with `ID` or `ensembl_transcript_id` and `cdna`, or FASTA | Positional profiling |

FASTA identifiers are the first whitespace-delimited token after `>`. Convert U to T for counting, resolve ambiguous nucleotides, and keep identifiers disjoint across reference/control files. The current counter does not automatically standardize the alphabet. Non-ACGT windows do not match the enumerated k-mer universe. Use sequences at least k nucleotides long.

Use matching annotation releases for sequences, region boundaries and mapped experimental observations.

## Outputs

KEA writes `all_kmers.tsv`, `enriched.txt` and `depleted.txt` for each requested selection. KRS writes `signature.txt` and `kmer_ranks.tsv`; the full rank matrix is optional and can be large. K-RBP returns a DataFrame; save it explicitly with `to_csv`. The [API reference](API-reference) records parameters, return values and output locations.

## Scope of this distribution

The supported tutorial covers FASTA counting, KEA, KRS and continuous-mode K-RBP. Advanced genomic mapping and BEDTools helpers require separately configured tools and resources. See [Known limitations](Known-limitations) before using them.

## Navigation

- [Installation](Installation)
- [Tutorial](Tutorial)
- [Methodology](Methodology)
- [Detailed parameter guide](Module-parameters)
- [API reference](API-reference)
- [Data and reproducibility](Data-and-reproducibility)
- [FAQs](FAQs)
- [Citation](Citation)

## Kmap and CDS domain analysis

- [Kmap: similarity metrics, parameters and Sankey visualization](Kmap)
- [K-CDS: inputs, frame, backgrounds and Fisher tests](K-CDS)

### Kmap: compare individual RNA profiles

Rank candidate transcripts against a query using complete normalized k-mer profiles. Visualize selected pairs against a background using binned Sankey-style rank alignment. See [Kmap](Kmap).

![Kmap: compare the target RNA k-mer frequency profile with matching and low-similarity profiles](https://github.com/AngeloDAngelo/Ktools/raw/refs/heads/main/docs/assets/kmap-workflow.png)

### K-CDS

Test selected CDS occurrence overlaps with UniProt genomic features using control transcripts or a shuffled-position null. Inputs, frame filtering, QC and statistical interpretation are described in [K-CDS](K-CDS).

## K_miR

[miRNA family annotation: orientation, parameters, scores and plots](MiRNA-seeds) · [CDR1as worked example](CDR1as-miRNA-example)
