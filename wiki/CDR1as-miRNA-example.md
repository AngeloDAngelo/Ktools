# K-miR example: CDR1as / ciRS-7

[K-miR](K-miR) · [Tutorial](Tutorial)

CDR1as is a published miR-7-binding circular RNA ([Hansen et al., 2013](https://doi.org/10.1038/nature11993)). This example annotates its selected seven-nucleotide words with complementary human miRNA seeds.

## Sequence

The [circBase record hsa_circ_0001946](https://www.circbase.org/cgi-bin/singlerecord.cgi?id=hsa_circ_0001946) defines a 1,485-nt, positive-strand circular interval at hg19 chrX:139865339–139866824 (zero-based, half-open). Ensembl assembly mapping places it at GRCh38 X:140783175–140784659 (one-based, inclusive). The sequence was retrieved from the Ensembl genomic sequence endpoint. [Coordinates, retrieval URL and checksum](https://github.com/AngeloDAngelo/Ktools/blob/main/examples/data/CDR1as_provenance.json).

## Run

Supply a background FASTA and use a new output directory:

```bash
python examples/run_cdr1as_mirna.py --background /path/to/lncRNA.fa --out results/cdr1as_mirna --top-pct 5
```

The example script streams transcripts rather than creating a dense count matrix. It normalizes canonical counts per transcript, excludes zero-frequency observations for each word, includes the target, and averages tied ranks. It then calls `AnnotateMiRNASeeds` on words with KRS percentile ≥0.95. Duplicate identifiers are rejected; identical target sequences are excluded from the background.

The script uses the reference-class counting boundary: 1,478 seven-nucleotide windows from the 1,485-nt linearized sequence. It omits the final possible window and does not count junction-spanning windows.

## Demo plot

![K-miR CDR1as demo: matched miRNA families and target k-mer percentile scores](https://github.com/AngeloDAngelo/Ktools/raw/refs/heads/main/examples/cdr1as_demo/seed_scores.png)

`AnnotateMiRNASeeds(source="KRS", plot=True)` produces this bar-plot format. Each row identifies a matched miRNA family and its complementary target word; bar length shows the target's KRS percentile score. The demo plots up to 20 matched rows, ordered by percentile. Results depend on the comparison FASTA and selection parameters.

The miR-7-5p seed `GGAAGAC` matches `GTCTTCC` (DNA) / `GUCUUCC` (RNA). The target word occurs 67 times in the included CDR1as sequence.

Exact seed matches nominate compatible miRNAs. They do not establish expression, binding or sponge activity; KRS percentiles are not P values.
