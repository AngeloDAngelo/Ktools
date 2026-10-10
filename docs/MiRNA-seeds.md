# K_miR: miRNA family annotation from 7-mer seeds

[Documentation](Documentation.md) · [Module parameters](Module-parameters.md) · [CDR1as example](CDR1as-miRNA-example.md)

`AnnotateMiRNASeeds` merges RNA sequence signatures with family seed annotations. It reports candidate seed-compatible miRNAs, not a new miRNA-binding probability or a family-level enrichment test.

## Orientation and site definition

The supplied `Seed+m8` is positions 2-8 of the mature miRNA. All 2,606 human rows were checked against their mature sequences. [TargetScan defines families using these positions](https://www.targetscan.org/docs/seed.html).

To match a transcript-oriented query word, reverse complement the seed. Example: miR-7-5p seed `GGAAGAC` gives target `GUCUUCC` in RNA, or `GTCTTCC` in DNA.

This module annotates exact **7mer-m8** sites only. It does not include wobble, mismatches, 6mer, 7mer-A1, 8mer/context classification, accessibility, conservation scoring or miRNA expression. An 8mer includes a matching 7mer-m8 segment, but this join does not classify it as an 8mer. See [TargetScan site definitions](https://www.targetscan.org/docs/7mer.html).

## Parameters

```python
AnnotateMiRNASeeds(self, kmers=None, source="KEA",
    reference=None, species=None, species_id=9606,
    plot=True, top_n=20, output_dir=None)
```

| Parameter | Behavior |
| --- | --- |
| `analysis.seed_family_file` | Path to the tab-separated family table. Defaults to the complete table installed with K-tools; set this attribute to use a custom file. Required columns: `miR family`, `Seed+m8`, `Species ID`, `MiRBase ID`. If `Mature sequence` is present, positions 2–8 are validated. |
| `species_id` | Taxonomy ID for row filtering, default human 9606. The bundled human subset cannot annotate other species. |
| `source` | `KEA`: score join from stored delta-median and stat_log2fc tables; `KRS`: score join from stored KRS target results; `manual`: sequence annotation without analysis scores. |
| `kmers` | Explicit canonical 7-mers. If omitted in KEA, use the intersection of enriched delta-median and stat_log2fc lists. If omitted in KRS, use the stored selected signature. Required for manual mode. DNA/RNA alphabets are accepted and normalized to DNA. Empty signatures return empty matches. Other lengths raise an error rather than silently decomposing or transferring parent scores. |
| `reference` | Optional validation of the stored KRS target ID. To analyze another target, rerun KRS first. |
| `species` | Count/results key, default `input`; distinct from taxonomy `species_id`. |
| `plot` | Draw/save scored annotation plots, default True. Manual mode has no analysis scores and does not plot. |
| `top_n` | Positive integer, default 20; limit plotted matched word/family rows, not returned annotations. KEA sorts by log2FoldChange; KRS by target percentile. |
| `output_dir` | New output directory; None creates a unique `MiRNASeeds_*` subdirectory of the analysis output. Existing directories are refused. |

## Analysis scores and interpretation

KEA has an integrated **selection** but no defined scalar combined score. The method retains both `delta_median` and `log2FoldChange` and plots them in separate panels for matched words from the integrated signature. Original `PValue`, `FDR` and `bonf` are retained when available. They are k-mer test values, not new miRNA-family P values.

KRS retains `target_frequency`, `target_rank_percentile` and `selected_KRS`. The plot shows target percentile, not a calibrated P value or fold enrichment.

Family rows sharing a target word are preserved. Multiple mature miRNAs within a family/seed are grouped into semicolon-separated IDs/accessions. This prevents treating each mature ID as an independent discovered site; family annotations sharing a word still share the same evidence. Plots label each family and word explicitly and do not sum scores across miRNAs.

Unmatched words remain in `query_kmers.tsv`; the annotation does not discard them silently. Supplied explicit words may be outside the selected signature; their scores are joined where available, and missing score values remain missing.

## Launch examples

The bundled table is selected automatically. To use another table, set `analysis.seed_family_file = "path/to/miR_Family_Info.txt"` before calling the method.

```python
# After k=7 counting and both KEA extraction branches:
mirna = analysis.AnnotateMiRNASeeds(
    kmers=signature,
    source="KEA",
    plot=True,
)

# After KRS(reference="my_target", ...):
mirna = analysis.AnnotateMiRNASeeds(
    source="KRS",
    reference="my_target",
    plot=True,
)

# Orientation check, independent of any analysis:
mirna = analysis.AnnotateMiRNASeeds(
    kmers=["GTCTTCC"], source="manual", plot=False,
)
print(mirna["matches"])
```

## Outputs

Returns `matches`, `query`, `seed_map`, `figures`, `settings`, `output_dir`; also stores `analysis.mirna_seed_annotation`.

- `seed_matches.tsv`: matched words, family/seed, mature miRNA IDs and original scores.
- `query_kmers.tsv`: all unique supplied words and matched-family count.
- `family_seed_map.tsv`: filtered family-to-target-word mapping.
- `settings.json`: inputs, site type, selected words and parameters.
- `seed_scores.png` / `seed_scores.pdf` when scored matches are plotted.

Figure objects remain available for the caller to display or close. No extra runtime dependencies are needed.
