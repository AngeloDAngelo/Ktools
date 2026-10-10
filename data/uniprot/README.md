# UniProt genomic feature tracks

Upload the UniProt BED tracks used by `KEA.K_CDS` here. No BED resources are included yet.

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
result = analysis.K_CDS(
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

For each uploaded track, record source URL, download date, genome assembly, annotation version, label definition and reuse/citation terms. See [domain parameters](../../docs/K-CDS.md).

## BED files for the localization analysis

Upload your original `unipLocCytopl.bed`, `unipLocExtra.bed`, `unipLocSignal.bed` and `unipLocTransMemb.bed` into this directory, then select these four names with `tracks`. No BED files are bundled yet.

The function is not restricted to a particular physical copy of a BED file: it requires the filenames, columns, labels, strand and genome assembly described above. In particular, `unipLocSignal` keeps rows whose feature label is exactly `Signal peptide`. Preserve the annotation format of your working files.

`gtf_file` defaults to `None` and must be supplied for a matching GTF. Pass this directory explicitly with `uniprot_dir="data/uniprot"`.

## What makes another BED compatible?

Compatible files are not limited to the original downloads. They must satisfy all of these requirements:

- **Format:** uncompressed, tab-separated genomic BED text with at least six columns: `chrom`, `chromStart`, `chromEnd`, `name`, `score`, `strand`. Start is zero-based and end is excluded, so `chr1\t100\t110\tFeature\t0\t+` represents ten nucleotides. Intervals must have start < end and strand `+` or `-`; unknown strand will not match the GTF. With BED12 and `split_blocks=True`, column 10 is blockCount, column 11 contains comma-separated blockSizes, and column 12 contains comma-separated blockStarts relative to chromStart. Counts must agree and all blocks must lie within the outer interval. Files use the selected track name, e.g. `unipLocCytopl.bed`.
- **Labels:** feature names come from column 4, except `unipDomain` and `unipInterest`, which use column 27 when present. Names must be nonempty. `unipLocSignal` retains only the exact, case-sensitive label `Signal peptide`; `unipInterest` retains only `Disordered`. Other tracks retain their labels as supplied. With multiple tracks, labels are prefixed with the track name. A different label convention must be adapted before running the method.
- **Assembly:** genomic coordinates must use the same reference assembly as the GTF (this workflow uses human GRCh38/hg38). Renaming an hg19 file does not convert its coordinates to hg38. BED chromosome names must match the projected GTF names: `chr1`, `chr2`, …, `chrX`, `chrY`, `chrM`. The method prefixes GTF names with `chr` and converts GTF `chrMT` to `chrM`, but it does not normalize BED chromosome names or lift coordinates between assemblies. Full-transcript FASTA sequences and IDs must also match the GTF annotation release.

Upload your working BED files to `data/uniprot`, preserving their names and annotations, and record source, assembly, release and download date.
