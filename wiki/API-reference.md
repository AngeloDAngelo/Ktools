# API reference

[Documentation](Documentation) · [Tutorial](Tutorial) · [Known limitations](Known-limitations)

Import with `from KEA import KEA`. The class exposes the analysis modules. Signatures below are extracted from the distributed code; scientific interpretation is described in [Methodology](Methodology).

For the behavior and practical effect of the parameters below, see the [detailed parameter guide](Module-parameters).

## Core workflow

### `__init__`

```python
__init__(self, dir_out, ref_fasta, ctrl_fasta)
```

Creates the output directory and stores reference/control FASTA paths. It does not load or count sequences.

### `CreateCombinedFasta`

```python
CreateCombinedFasta(self, species=None)
```

Concatenates input FASTA files and stores `combined_fasta`. Call before counting. Inputs must end with a newline and use distinct IDs. Returns None.

### `KmersCountsTable`

```python
KmersCountsTable(self, k, cores=5, subset_by_subseq=None, species=None)
```

Enumerates A/C/G/T k-mers, counts with a multiprocessing pool, and stores the table in `kmers_count_table_by_species["input"]` for default inputs. The `k` column is metadata. Returns None. Use a guarded Python script.

### `ExtractKmers`

```python
ExtractKmers(self, k_selected, species=None, extraction_methods=('delta_rank_percentile', 'delta_median', 'stat_log2fc'), n_iter=1000, sign_emp_pval=0.0025, log2fc_threshold=1.0, fdr_alpha=0.05, stat_test='mwu', correction_method='bonferroni', median=True)
```

Returns `KEA_results`, nested by species key then method. Each method holds `enriched`, `depleted` and `all_kmers`. Saves tables in `ExtractKmers/<species>/<method_parameters>/`. Request the two integrated methods explicitly. The default tail fraction is 0.0025.

### `KRS`

```python
KRS(self, reference, species=None, top_pct=5, store=True, save=True, save_full_rank_table=False, verbose=True)
```

Returns a list of signature k-mers. Stores a KRS result dictionary and, by default, saves `ExtractKmers/<species>/KRS/<safe_target_id>/signature.txt` and `kmer_ranks.tsv`. `top_pct=5` means rank ≥0.95.

### `KmersProteinEnrichment`

```python
KmersProteinEnrichment(self, kmers, df_zscore, k_out=None, deconvolute_unique=True, convert_to_rna=True, mode='continuous', test='mwu', threshold=None, top_n=None, fdr_alpha=0.05, plot=True, p_col='fdr', p_thresh=0.05, fc_thresh=1.0, highlight=None, highlight_top_n=10, save=None, bg_kmers=None, plot_type='lollipop', wordcloud_score='cliff_delta', wordcloud_top_n=10, wordcloud_font_path=None)
```

Returns a DataFrame containing all evaluated profiles. Core fields: `protein`, `mean_group1`, `mean_group2`, `median_group1`, `median_group2`, `pvalue`, `fdr`, `bonf`, `zscore_diff` (Cliff's delta), `n1`, `n2`. `save` saves figures, not the result table: call `to_csv` yourself. Use continuous mode and explicit plot thresholds. Historical binary-mode parameters are not implemented.

### `KmerCorrelationWithScore`

```python
KmerCorrelationWithScore(self, txids_and_score: pd.DataFrame, species=None, method='spearman', sign_emp_pval=0.05)
```

Accepts a two-column DataFrame (transcript ID, numeric score). Returns a correlation table and writes `KmerCorrelationWithScore/<species>/corr_analysis.<method>.txt` plus tail lists. For default inputs pass `species="input"`; ensure `KEA_results` is a dictionary.

### `KmerSlidingWindowProfile`

```python
KmerSlidingWindowProfile(self, kmers, tx_id, ann_df, fasta_file=None, kmer_groups=None, window_size=50, step=25, strand='+', n_clusters='auto', cluster_range=(2, 10), convert_to_rna=True, plot_line=True, plot_heatmap=True, plot_cluster_heatmap=True, plot_silhouette=True, save_prefix=None, verbose=True, improvement_threshold=0.1)
```

Accepts transcript sequence via `fasta_file` or annotation (`ID`/`ensembl_transcript_id` and `cdna`). Counts and clusters positional signatures. Advanced functionality not exercised by the minimal tutorial.

## Advanced method inventory

Genomic interval and mapping helpers require external tools and configured resources.

- `get_kmer_table_for_species`
- `LogoByList`
- `CompareTranscriptKmerProfiles`
- `plot_kmer_rank_alignment`
- `FetchOrthologsAndSaveFastaOneToOne`
- `IntersectEnrichedKmersAndPlotVenn`
- `ImportMetaTx`
- `ImportMetaTx_Mapped`
- `ImportMetaTx_MappedScore`
- `ImportMetaTx_Enrichment`
- `GenerateMetaTx`
- `MapperTxToGenome`
- `ScoreByBigWig`
- `EnrichmentAnalysis`
- `ScoreMatrix_Enrichment`
- `KmerBinFrequencyMatrix`
- `PlotKmerMatrixFreq`

## Standalone helpers

`kmers_counter(seq, k, all_kmers_df=None, tx_id=None)` counts overlapping windows using the reference-class boundary, which omits the final possible start. If an enumerated k-mer universe is supplied, absent entries are filled with zero.

`run_stat_test(x, y, test)` supports `mwu`, `ks` and `ttest`; it skips empty or jointly constant groups. Other helper functions support interval handling, annotations and plotting and are outside the portable workflow.

## Kmap and CDS domain analysis

- [Kmap: similarity metrics, parameters and Sankey visualization](Kmap)
- [K-CDS: inputs, frame, backgrounds and Fisher tests](K-CDS)

## `K_CDS`

```python
K_CDS(self, gtf_file=None, uniprot_dir=None,
    kmers=None, method="stat_log2fc", species=None,
    in_frame=False, background="control", tracks=("unipDomain",),
    alpha=0.05, fdr_scope="global", split_blocks=True,
    seed=42, output_dir=None, plot=True, max_domains=60,
    ref_fasta=None, ctrl_fasta=None)
```

Returns a dictionary of statistics, occurrences, QC, matrices, figure, settings and output directory; also stores `self.domain_enrichment`. Supply GTF and UniProt paths explicitly. Optional FASTA overrides apply only to this call. See [complete parameter explanations](K-CDS).

## Kmap calls

```python
CompareTranscriptKmerProfiles(self, target_tx, other_txs,
    target_species=None, other_species=None, metric="spearman",
    pseudocount=1e-9, k=None, fasta_input=False,
    fasta_species_key="_external_db")

plot_kmer_rank_alignment(self, ref_tx, ctrl_tx,
    ref_species="input", ctrl_species="input", bg_species="input",
    kmers=None, exclude_ctrl_from_bg=True, top_n=None,
    min_freq=None, bin_size=100, alpha_band=0.4, color_by="gc",
    kmer_categories=None, figsize=(10, 10), save=None)
```

See [Kmap](Kmap) for metric definitions, cached external FASTAs, ranked output and Sankey-style visualization.

## K_miR

```python
AnnotateMiRNASeeds(self, kmers=None, source="KEA", reference=None,
    species=None, species_id=9606, plot=True, top_n=20, output_dir=None)
```

`analysis.seed_family_file` selects the family table. The default is the resource installed with K-tools.

[miRNA family annotation: orientation, parameters, scores and plots](MiRNA-seeds) · [CDR1as worked example](CDR1as-miRNA-example)

K-CDS uses `K_CDS` as its Python method name; The former `DomainEnrichment` name is no longer exposed. Results also include `frame_distribution` and `frame_figure`, and are stored as `self.k_cds`. See [frame plot and saved outputs](K-CDS#cds-frame-distribution).

## Separate figures

With `plot=True`, the domain heatmap and CODING/SHUFFLE frame distribution are separate figures. Display `result["figure"]` for the heatmap and `result["frame_figure"]` for the frame plot. The frame legend is outside the SHUFFLE panel. Figures are saved as `heatmap.png`/`.pdf` and `frame_distribution.png`/`.pdf`.
