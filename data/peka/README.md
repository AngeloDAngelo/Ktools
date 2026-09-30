# PEKA matrices for K-RBP

This folder provides two PEKA score matrices for K-RBP analysis.

Choose the matrix you want to use. Both files share the same format, k-mer length and alphabet, so the analysis pipeline remains unchanged: only the input filename changes.

## Available matrices

| File | K-mers | RBP profiles |
| --- | ---: | ---: |
| [eclip_clippy.txt](eclip_clippy.txt) | 1,024 RNA 5-mers | 223 |
| [eclip_nr.txt](eclip_nr.txt) | 1,024 RNA 5-mers | 215 |

Each column represents an RBP/experimental profile. Preserve the complete column names, including cell-line information, when reporting results.

The files have different profile sets and missing-value patterns and may therefore produce different associations.

## File format

Although their extension is `.txt`, both files are tab-separated tables:

- First column: unique RNA 5-mer identifiers.
- Remaining columns: numeric scores for individual profiles.
- Alphabet: uppercase A/C/G/U.
- Missing scores: read as `NaN` by pandas.

Both matrices contain the complete set of 1,024 possible RNA 5-mers.

## Choose and load a matrix

Run from the repository root and choose either filename:

```python
import pandas as pd

matrix_path = "data/peka/eclip_clippy.txt"
# Alternatively:
# matrix_path = "data/peka/eclip_nr.txt"

matrix = pd.read_csv(
    matrix_path,
    sep="\t",
    index_col=0,
)
```

Do not replace missing scores with zero: missing values and measured zero scores have different meanings.

K-RBP uses the supplied scores directly. It does not calculate PEKA scores or standardize matrix columns.

## Run K-RBP

After generating a KEA or KRS signature, pass the chosen matrix as `df_zscore`:

```python
results = analysis.KmersProteinEnrichment(
    kmers=signature,
    df_zscore=matrix,
    k_out=5,
    deconvolute_unique=True,
    convert_to_rna=True,
    mode="continuous",
    test="mwu",
    fdr_alpha=0.05,
    plot=False,
)
```

This call is identical for both matrices.

- `k_out=5` decomposes longer signature words into overlapping 5-mers. Query words must have a consistent length of at least 5.
- `convert_to_rna=True` converts uppercase T to U, matching the matrices' RNA alphabet.
- The default background consists of matrix words outside the query. Large backgrounds are subsampled by the implementation.
- `test="mwu"` tests whether query scores are higher than background scores for each profile.
- Missing scores are excluded independently for each profile.
- Profiles with fewer than two valid observations in either group are skipped.
- `plot=False` returns the statistics without generating figures.

No changes to `KEA.py` are required when switching between these files.

## Save results

Save the complete result table and record which matrix was used:

```python
from pathlib import Path

output_dir = Path("results")
output_dir.mkdir(parents=True, exist_ok=True)

matrix_name = Path(matrix_path).stem

results.to_csv(
    output_dir / f"krbp_{matrix_name}_results.tsv",
    sep="\t",
    index=False,
)
```

Using the matrix name in the output filename avoids overwriting results when comparing both inputs.

## Interpret the output

The returned table contains all evaluated profiles.

| Column | Meaning |
| --- | --- |
| `protein` | Original matrix profile name |
| `pvalue` | Raw one-sided test P value |
| `fdr` | Benjamini–Hochberg-adjusted P value |
| `bonf` | Bonferroni-adjusted P value |
| `zscore_diff` | Cliff's delta, despite the historical column name |
| `n1` | Valid query observations |
| `n2` | Valid background observations |

Positive Cliff's delta indicates higher query scores.

`fdr_alpha` controls significance flags; it does not remove rows from the returned table. Apply your chosen adjusted-P and positive-effect thresholds explicitly.

## Provenance

When reporting an analysis, identify the exact matrix filename and document its source, PEKA processing version, profile selection and any preprocessing.

The source and preprocessing differences between these two files should be documented explicitly rather than inferred from their filenames.

## Documentation

- [Tutorial](../../docs/Tutorial.md)
- [Module parameters](../../docs/Module-parameters.md)
- [Methodology](../../docs/Methodology.md)
