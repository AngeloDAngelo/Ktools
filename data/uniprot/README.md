# UniProt genomic feature tracks

Upload the UniProt BED tracks used by `KEA.DomainEnrichment` here. No BED resources are included yet.

## Filenames

The function loads `<track>.bed` from the directory passed as `uniprot_dir`:

- `unipDomain.bed` (default domain track)
- `unipInterest.bed` (Disordered features only)
- `unipStruct.bed`
- `unipLocTransMemb.bed`
- `unipLocExtra.bed`
- `unipLocSignal.bed` (Signal peptide features only)
- `unipLocCytopl.bed`
- `unipRepeat.bed`

Only upload tracks you use; unselected tracks are not required. Choose them through `tracks`. Preserve exact filenames or rename your files to this convention.

## Format and compatibility

Use uncompressed genomic BED text with at least six tab-separated fields: chromosome, zero-based start, half-open end, feature name, score and strand. BED12 blocks are respected with `split_blocks=True`. Standard UCSC extended UniProt files may have feature labels in column 27: this is used for `unipDomain`/`unipInterest` when present, otherwise column 4 is used. Review the label columns before substituting another BED format.

Resources must match the assembly used by the transcript/GTF inputs (documented workflow: GRCh38/hg38). BED chromosome names must match the GTF projection (`chr1`, ..., `chrM`). Preserve + or - strand. GTF coordinates are converted automatically; do not convert a BED to 1-based coordinates.

```python
result = analysis.DomainEnrichment(
    gtf_file="my_annotation.gtf",
    uniprot_dir="data/uniprot",
    ref_fasta="my_reference_transcripts.fa",
    ctrl_fasta="my_control_transcripts.fa",
    kmers=signature,
    tracks=("unipDomain",),
    background="control",
    plot=True,
)
```

The FASTA overrides are optional: omit them to use the object's reference/control files. These must contain full spliced transcripts matching the GTF, not only CDS sequences. All paths are relative to the working directory unless absolute.

For each uploaded track, record source URL, download date, genome assembly, annotation version, label definition and reuse/citation terms. See [domain parameters](../../docs/Domain-enrichment.md).
