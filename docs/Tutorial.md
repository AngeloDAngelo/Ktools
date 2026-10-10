# Tutorial: simulated RNA sequences

[Documentation](Documentation.md) · [Methodology](Methodology.md) · [FAQs](FAQs.md)

This tutorial uses **only simulated data**. It demonstrates KEA, KRS and K-RBP with fictional sequences and RBP profiles. All commands run from the repository root after [installation](Installation.md).

## Run everything

```bash
python examples/run_demo.py --out results/demo
```

The script uses 24 reference and 48 control sequences, each 600 nt long. Reference sequences carry an injected `ACG` repeat. It counts 3-mers so the example runs quickly, then creates a fictional RBP score matrix. Profile names beginning with `SIMULATED_` have no biological meaning.

To regenerate the inputs:

```bash
python examples/generate_data.py
```

Sequence generation uses seed 42. See [demonstration checks](Benchmarks.md) for expected outputs.

## 1. Load sequences and count k-mers

The FASTA inputs are `examples/data/reference.fa` and `examples/data/control.fa`. They use unique IDs and uppercase A/C/G/T. Save the following steps together in a Python script; keep analysis calls inside the `__main__` block because counting uses multiprocessing.

```python
from pathlib import Path
import numpy as np
import pandas as pd
from KEA import KEA

if __name__ == "__main__":
    np.random.seed(42)
    out = Path("results/tutorial")
    analysis = KEA(str(out), "examples/data/reference.fa", "examples/data/control.fa")
    analysis.CreateCombinedFasta()
    analysis.KmersCountsTable(k=3, cores=1)
```

The count matrix has 4³=64 k-mer rows, one column per sequence, and a `k` metadata column. It is available as `analysis.kmers_count_table_by_species["input"]`. Each 600-nt input has 597 counted overlapping 3-mer starts (the reference counter omits the last possible window).

## 2. Compare the two groups with KEA

Continue inside the same block:

```python
    selected = analysis.ExtractKmers(
        k_selected=3,
        extraction_methods=("delta_median", "stat_log2fc"),
        n_iter=100,
        sign_emp_pval=0.10,
        log2fc_threshold=0.5,
        correction_method="fdr",
        fdr_alpha=0.05,
    )["input"]
    signature = sorted(
        set(selected["delta_median"]["enriched"])
        & set(selected["stat_log2fc"]["enriched"])
    )
    print(signature)
    assert "ACG" in signature
```

The delta-median tail is 10% here because there are only 64 sequence types. `sign_emp_pval` is a selected fraction, not a P value. The transcript-level test uses BH-adjusted P≤0.05 and log₂FC≥0.5. The integrated signature is their intersection. An empty method selection must produce an empty intersection.

The full script saves `kea_signature.txt` and writes separate method tables under `ExtractKmers/input/`. Inspect `all_kmers.tsv` before interpreting filtered lists. A larger k or much smaller dataset can produce undefined delta-median scores when a quartile spread is zero.

## 3. Characterize one simulated RNA with KRS

```python
    krs_signature = analysis.KRS(
        reference="reference_006",
        top_pct=10,
        save=True,
        save_full_rank_table=False,
    )
```

This selects k-mers whose normalized frequency in `reference_006` reaches percentile ≥0.90 among the sequences containing that k-mer. The target is part of the ranked population. `top_pct=10` does not mean selecting 10% of its k-mer types.

The output directory is `ExtractKmers/input/KRS/reference_006/`. Open `kmer_ranks.tsv` to inspect `target_frequency`, `target_rank_percentile`, and `selected_KRS`; `signature.txt` contains one selected k-mer per line.

## 4. Build a fictional RBP matrix

The following deliberately creates a profile that favors the recovered signature, a neutral profile, and a profile that disfavors it. This tests the API and effect direction; it is not evidence of predictive performance.

```python
    kmers = sorted(analysis.kmers_count_table_by_species["input"].index)
    rng = np.random.default_rng(123)
    scores = pd.DataFrame({
        "SIMULATED_query_favoring": [
            8 + rng.normal(0, 0.2) if k in signature else rng.normal(0, 0.2)
            for k in kmers
        ],
        "SIMULATED_neutral": rng.normal(0, 1, len(kmers)),
        "SIMULATED_query_disfavoring": [
            -8 + rng.normal(0, 0.2) if k in signature else rng.normal(0, 0.2)
            for k in kmers
        ],
    }, index=pd.Index(kmers, name="kmer"))
    scores.to_csv(out / "synthetic_rbp_matrix.tsv", sep="\t")
```

The table has k-mers as rows and profiles as columns. For a real analysis, the matrix must come from an appropriate experimental resource. Here, it remains completely synthetic.

## 5. Run K-RBP

```python
    background = [word for word in kmers if word not in signature][:10 * len(signature)]
    rbp = analysis.KmersProteinEnrichment(
        kmers=signature,
        df_zscore=scores,
        bg_kmers=background,
        convert_to_rna=False,
        mode="continuous",
        test="mwu",
        plot=False,
    )
    rbp.to_csv(out / "synthetic_rbp_results.tsv", sep="\t", index=False)
    candidates = rbp[(rbp["fdr"] < 0.05) & (rbp["zscore_diff"] >= 0.47)]
    print(candidates[["protein", "fdr", "zscore_diff"]])
```

The demo uses the first 30 non-query words in alphabetical order as an explicit background, avoiding process-dependent background subsampling. The query-favoring profile should have positive Cliff's delta and the disfavoring profile negative delta. **`zscore_diff` contains Cliff's delta**, despite its historical name. The one-sided test asks whether query scores exceed background scores. The full result table is returned independently of plotting filters.

The tutorial uses a 3-mer matrix, so no decomposition is required. For signatures longer than the matrix k-mer length, set `k_out` to the matrix length and retain `deconvolute_unique=True`. Do not request `k_out` larger than the query length.

## 6. Optional plot

```python
    analysis.KmersProteinEnrichment(
        signature, scores, bg_kmers=background,
        convert_to_rna=False,
        plot=True, plot_type="volcano",
        p_col="fdr", p_thresh=0.05, fc_thresh=0.47,
        save=str(out / "synthetic_rbp_volcano.png"),
    )
```

Set `fc_thresh` explicitly. The inherited default 1.0 excludes every point under a strict Cliff's-delta filter. Plot filtering does not remove rows from the returned DataFrame.

## Expected result

The verified demonstration recovers the KEA signature `ACG`, `CGA`, `GAC`, produces 8 KRS k-mers for `reference_006`, and evaluates 3 fictional RBP profiles. Reference outputs are saved in [examples/expected](../examples/expected/). Floating-point statistics can vary slightly with dependency versions.

## Output checklist

| File or folder | Meaning |
| --- | --- |
| `kea_signature.txt` | Integrated simulated-group signature, saved by `run_demo.py` |
| `krs_signature.txt` | Single simulated transcript's signature, saved by `run_demo.py` |
| `ExtractKmers/input/` | Full tables and selection lists |
| `synthetic_rbp_matrix.tsv` | Fictional input scores |
| `synthetic_rbp_results.tsv` | Complete association statistics |
| `summary.json` | Seeds, dimensions and signature summary from `run_demo.py` |

The demonstration includes assertions for the injected motif and expected effect direction. For your own runs, record input checksums, parameters, seeds and environment versions. Consult [Known limitations](Known-limitations.md) when changing input sizes or using advanced helpers.

## Kmap and CDS domain analysis

- [Kmap: similarity metrics, parameters and Sankey visualization](Kmap.md)
- [K-CDS: inputs, frame, backgrounds and Fisher tests](K-CDS.md)

## K_miR

[miRNA family annotation: orientation, parameters, scores and plots](MiRNA-seeds.md) · [CDR1as worked example](CDR1as-miRNA-example.md)
