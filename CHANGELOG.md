# Changelog

## 0.1.0.dev0 — repository preparation

- Added packaging, documentation, tutorials, wiki source pages, original figure extracts, synthetic demonstration inputs and regression tests.
- Corrected `kmers_counter` to count the final valid start position (`range(len(seq)-k+1)`). A sequence of length k now yields one occurrence, and `AAAA`, k=2 yields three `AA` occurrences.
- Added missing `venn2` and `venn3` imports from matplotlib-venn.
- Sorted set-derived query/background lists in K-RBP before selection and seeded background sampling, making input ordering stable across processes.

The counting fix and stable background ordering can change scientific results. This snapshot is not certified to reproduce manuscript figures. Preserve the original analysis environment and rerun relevant analyses before choosing the paper-associated release. The exact changes to the supplied source are recorded in `provenance/source-changes.patch`; source checksums are recorded in `provenance/source-manifest.json`.

## Added: CDS domains and Kmap documentation

- Added supplied DomainEnrichment method with explicit resource paths instead of laboratory defaults; additional analysis/plot settings recorded.
- Documented existing transcript similarity and Sankey rank alignment without changing their calculations.
