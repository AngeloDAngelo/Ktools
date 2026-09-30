# PEKA Z-score matrices for K-RBP

Place your experimental PEKA k-mer × RBP-profile matrices here. This folder currently contains no experimental matrix. The tutorial uses only a synthetic matrix.

Recommended format: a tab-separated `.tsv` file whose first column contains unique, equal-length k-mers and whose remaining columns contain numeric Z scores for individual RBP/cell-line profiles. Preserve profile identity in column names. Use uppercase RNA k-mers (A/C/G/U) with `convert_to_rna=True`, or uppercase DNA k-mers (A/C/G/T) with `convert_to_rna=False`. Missing values can be represented as empty cells; profiles with too few valid observations are skipped.

```python
import pandas as pd
matrix = pd.read_csv("data/peka/your_matrix.tsv", sep="\t", index_col=0)
results = analysis.KmersProteinEnrichment(
    signature, matrix, k_out=5, convert_to_rna=True,
    mode="continuous", test="mwu", plot=False,
)
results.to_csv("results/rbp_associations.tsv", sep="\t", index=False)
```

Replace the filename and `k_out` with your actual matrix name and k-mer length. Record the source dataset/accession, PEKA version, processing parameters, score definition, alphabet, k-mer length and reuse terms alongside every matrix. No real score values or dataset provenance have been invented.
