# Changelog

## Current development version

- Integrated autonomous K-CDS domain/positional panels with per-k-mer cache reuse, automatic missing-profile calculation, Ward.D2 on capped centered profiles and In Frame Index ordering. Removed legacy Matplotlib plots. Default background is now shuffle and default tracks are the four localization tracks.
- Added `kmer_protein_scores` to K-RBP by default, preserving original signature words and averaging their overlapping pentamer scores without changing enrichment results.
- Updated parameter/API guides, installation requirements and Wiki source pages.

- Updated KRS and K-RBP workflow diagrams and the K-miR CDR1as example figure.
- Aligned `AnnotateMiRNASeeds` with the reference class: the family table is configured through `seed_family_file`.
- Included the complete miRNA family table in the installed package and removed laboratory-specific default paths.
- Matched reference-class counting and K-RBP background ordering; documented the terminal-window boundary.
- Updated documentation, runnable examples and citation guidance.

## Modules

- KEA, KRS and K-RBP analysis and visualization.
- K-map transcript similarity and Sankey-style rank alignment.
- K-CDS frame distributions and UniProt feature overlap analysis.
- K-miR exact reverse-complement miRNA seed annotation.
