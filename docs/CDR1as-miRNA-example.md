# CDR1as / ciRS-7: miR-7 seed annotation

[miRNA module](MiRNA-seeds.md) · [Tutorial](Tutorial.md)

This is a separate real-data case study, not the synthetic core tutorial. CDR1as/ciRS-7 is a published miR-7-binding circRNA ([Hansen et al., 2013](https://doi.org/10.1038/nature11993)).

## Sequence provenance

The [circBase record hsa_circ_0001946](https://www.circbase.org/cgi-bin/singlerecord.cgi?id=hsa_circ_0001946) defines a single-exon 1,485-nt circRNA on the positive strand, hg19 BED interval chrX:139865339-139866824. Ensembl's assembly-mapping endpoint maps this intact interval to GRCh38, 1-based X:140783175-140784659 (+). Sequence was retrieved from Ensembl's GRCh38 genomic sequence endpoint. The response, coordinate conventions and sequence checksum are in [CDR1as_provenance.json](../examples/data/CDR1as_provenance.json).

This is a genomic extraction of the mapped circRNA boundaries, **not an Ensembl 99 linear-transcript sequence**. The supplied `lncRNA.fa` was designated by the author as the comparison population; its headers alone do not establish the annotation release. Its checksum is recorded in the output summary. It is not bundled; provide the same file to reproduce these ranks.

## Run

```bash
python examples/run_cdr1as_mirna.py --background /path/to/lncRNA.fa --out results/cdr1as_mirna --top-pct 5
```

Use a new output directory. The script calculates target-present 7-mer KRS ranks progressively, without a dense transcriptome matrix. It normalizes canonical counts per transcript, excludes zero-frequency observations for each word, includes the target, and uses average ranks for ties: `(number_lower + (number_equal_including_target + 1)/2) / number_nonzero_including_target`.

It then calls `AnnotateMiRNASeeds` on the selected KRS signature. Duplicate identifiers are rejected and identical target sequences in the background are excluded. Counting uses a linearized single traversal: six possible junction-crossing 7-mer windows are not included. Neither this example nor the core counter claims circular junction coverage.

## Observed results against the supplied background

- 17,781 usable background transcripts.
- 77 target words selected at percentile >=0.95.
- miR-7-5p seed: `GGAAGAC`.
- Target word: `GTCTTCC` (DNA) / `GUCUUCC` (RNA).
- Target word count: **67**.
- Frequency: **67/1479 = 0.04530**.
- KRS percentile: **1.0**, among 1,823 nonzero observations including the target.

These 67 exact 7mer-m8 occurrences are not a reproduction of the publication's total binding-site count, which uses its own site definitions and sequence material. A rank of 1 is not a binding probability or statistical P value. The seed match cannot distinguish expression/activity of family members.

![Matched family KRS scores in CDR1as](../examples/cdr1as_expected/seed_scores.png)

[Full annotations](../examples/cdr1as_expected/seed_matches.tsv), [miR-7 word score](../examples/cdr1as_expected/mir7_target_score.tsv), [all target-present word ranks](../examples/cdr1as_expected/target_kmer_ranks.tsv) and [run summary](../examples/cdr1as_expected/summary.json) are included.
