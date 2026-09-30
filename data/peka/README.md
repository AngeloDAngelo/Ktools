# PEKA Z-score matrices for K-RBP

This folder contains two PEKA k-mer × RBP-profile Z-score matrices.

Choose the matrix appropriate for your analysis. The K-RBP pipeline is the same for both files: only the input filename changes. No changes to `KEA.py` are required.

## Available matrices

| File | Description |
| --- | --- |
| `FIRST_MATRIX.tsv` | <!-- Describe the dataset and experimental context. --> |
| `SECOND_MATRIX.tsv` | <!-- Describe the dataset and experimental context. --> |

## Matrix format

Matrices are tab-separated files with:

- K-mers in the first column, used as row identifiers.
- Individual RBP profiles in the remaining columns.
- Numeric PEKA Z scores as values.

K-mer identifiers must be unique, uppercase and of consistent length. Preserve RBP and experimental context, such as cell line, in the profile names.

## Choose and load a matrix

Run the analysis from the repository root. Set `matrix_path` to either file:

```python
import pandas as pd

matrix_path = "data/peka/FIRST_MATRIX.tsv"
# Alternatively:
# matrix_path = "data/peka/SECOND_MATRIX.tsv"

matrix = pd.read_csv(
    matrix_path,
    sep="\t",
    index_col=0,
)
```

The loaded DataFrame is passed directly to K-RBP as `df_zscore`. K-RBP uses the supplied scores; it does not calculate PEKA scores or standardize the matrix.

## Run K-RBP

After obtaining a KEA or KRS signature, use the same analysis call for either matrix:

```python
results = analysis.KmersProteinEnrichment(
    kmers=signature,
    df_zscore=matrix,
    k_out=5,
    convert_to_rna=True,
    mode="continuous",
    test="mwu",
    fdr_alpha=0.05,
    plot=False,
)

results.to_csv(
    "results/rbp_associations.tsv",
    sep="\t",
    index=False,
)
```

This example assumes 5-mer matrices and a query signature containing words of length at least 5.

- `k_out=5` decomposes longer signature words into overlapping 5-mers.
- `convert_to_rna=True` converts uppercase T to U in the query and matrix row labels.
- The default background contains matrix k-mers outside the query.
- `test="mwu"` tests whether query scores are higher than background scores.
- `plot=False` returns the statistical results without generating figures.

If both matrices share the same k-mer length and alphabet, switching between them requires changing only `matrix_path`.

## Interpret the output

The returned table includes all evaluated RBP profiles.

- `pvalue`: raw one-sided test P value.
- `fdr`: Benjamini–Hochberg-adjusted P value.
- `bonf`: Bonferroni-adjusted P value.
- `zscore_diff`: Cliff's delta, despite its historical name.
- `n1`, `n2`: valid query and background observations for each profile.

Positive Cliff's delta indicates higher query scores. Apply your chosen adjusted-P and positive-effect thresholds explicitly to the returned table.

Missing scores are excluded independently for each profile. Profiles with fewer than two valid observations in either group are skipped.

## Matrix provenance

For each file, document:

- Source publication or dataset accession.
- RBP profiles and experimental context.
- PEKA version and processing parameters.
- K-mer length and alphabet.
- Any filtering or preprocessing.
- Applicable citation and reuse terms.

## Documentation

- [Tutorial](../../docs/Tutorial.md)
- [Module parameters](../../docs/Module-parameters.md)
- [Methodology](../../docs/Methodology.md)
