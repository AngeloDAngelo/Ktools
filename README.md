# K-tools

### From RNA sequence to candidate regulatory elements and proteins

K-tools identifies k-mer signatures in transcript populations and individual RNAs, then associates those signatures with experimentally derived RNA-binding protein (RBP) profiles. It works with transcript-oriented sequences and does not require preselected binding regions.

**[Documentation](docs/Documentation.md) · [Tutorial](docs/Tutorial.md) · [FAQs](docs/FAQs.md) · [Benchmarks](docs/Benchmarks.md) · [Methodology](docs/Methodology.md) · [Wiki](https://github.com/AngeloDAngelo/Ktools/wiki)**

![KEA workflow: compare group and transcript frequencies, then intersect selected k-mers](docs/assets/kea-workflow.png)

![KRS workflow: rank each target k-mer against the RNA population](docs/assets/krs-workflow.png)

![K-RBP workflow: compare signature k-mers with experimental RBP score profiles](docs/assets/krbp-workflow.png)

The KRS artwork uses the historical label “KPE” for the downstream step now called K-RBP.

## Choose a module

| Module | Biological question | Input | Main output |
| :--- | :--- | :--- | :--- |
| **KEA** — K-mer enrichment analysis | Which sequence elements distinguish two RNA populations? | Reference and control FASTA files | Enriched and depleted k-mers |
| **KRS** — K-mer RNA signature | Which k-mers are unusually abundant in one RNA? | Target RNA and a background population | Transcript-specific percentile-rank signature |
| **K-RBP** — K-mer RBP association | Which RBP profiles favor the signature? | Signature and k-mer × RBP score matrix | Association statistics and candidate proteins |

All three modules are exposed through the Python `KEA` class. Positional profiling can further locate signature elements along a transcript.

## Install

Download this repository using **Code → Download ZIP**, extract it, and open a terminal inside the extracted directory.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
python examples/run_demo.py --out results/demo
```

On Windows, activate with `.venv\Scripts\activate` instead. The demonstration uses small synthetic inputs and a synthetic RBP score matrix, clearly labeled as such. It exercises KEA, KRS and K-RBP; it does not reproduce the paper benchmarks. See [installation details](docs/Installation.md).

## Dependencies

Python ≥3.10 is required. `python -m pip install .` installs the runtime libraries declared in `pyproject.toml`: NumPy, pandas, SciPy, scikit-learn, statsmodels, Biopython, matplotlib, seaborn, Logomaker, pybiomart, WordCloud and matplotlib-venn. See [installation and library roles](docs/Installation.md) for version constraints and development installation.

## How the modules work

**KEA** counts overlapping k-mers in reference and control RNA populations. It combines a population-level frequency comparison with a statistical comparison of normalized transcript frequencies. Intersecting their selected lists yields a signature supported by both analyses.

**KRS** uses the combined RNA population to rank each k-mer's frequency in a target transcript against its nonzero frequencies in other transcripts. The chosen upper percentile defines the target's signature; background composition and the percentile threshold determine its interpretation.

**K-RBP** compares scores for signature k-mers against background k-mers in each RBP profile. Longer signature elements can be decomposed into shorter matrix k-mers. The output includes association statistics, multiple-testing corrections and effect sizes. These associations nominate candidate proteins; they are not direct binding measurements for the target RNA.

See [module parameters](docs/Module-parameters.md), [methodology](docs/Methodology.md) and the [simulated tutorial](docs/Tutorial.md). Place experimental PEKA matrices in [data/peka](data/peka/README.md).

## Manuscript and future DOI

- Paper/preprint DOI: <!-- Insert the DOI once available. -->
- bioRxiv preprint URL: <!-- Insert the public preprint URL once available. -->
- Manuscript DOCX: <!-- Add a link to manuscript/your-file.docx when ready. -->

[Manuscript placeholder folder](manuscript/README.md). No manuscript DOCX is included in this package.

## Minimal analysis

Save this code as a Python script and run it after installation:

```python
import numpy as np
from KEA import KEA

if __name__ == "__main__":
    np.random.seed(42)
    analysis = KEA("results/my_analysis", "reference.fa", "control.fa")
    analysis.CreateCombinedFasta()
    analysis.KmersCountsTable(k=7, cores=1)
    results = analysis.ExtractKmers(
        k_selected=7,
        extraction_methods=("delta_median", "stat_log2fc"),
        sign_emp_pval=0.01,
        log2fc_threshold=1.0,
        correction_method="bonferroni",
        fdr_alpha=0.05,
    )
    selected = results["input"]
    signature = sorted(
        set(selected["delta_median"]["enriched"])
        & set(selected["stat_log2fc"]["enriched"])
    )
    print(signature)
```

Use uppercase A/C/G/T FASTA sequences, with unique identifiers across both files. The `__main__` guard is required for multiprocessing on platforms using spawn. The [tutorial](docs/Tutorial.md) explains inputs, outputs, thresholds, KRS and K-RBP in detail.

## Benchmarks

K-tools has been evaluated using synthetic motif-injection designs and full-length RNA groups defined from eCLIP data. See [Benchmarks](docs/Benchmarks.md) for the evaluation setup, reported results and available reproducibility material.

**Status:** development snapshot; associated manuscript in preparation. See [CHANGELOG](CHANGELOG.md) for code changes and [Known limitations](docs/Known-limitations.md) for implementation caveats.

## Repository guide

| Location | Contents |
| :--- | :--- |
| `KEA.py` | Python implementation |
| `docs/` | Documentation, tutorial, FAQs, methodology, benchmarks and figures |
| `examples/` | Runnable demonstration and synthetic FASTA inputs |
| `benchmarks/` | Reported summary values and reproducibility inventory |
| `tests/` | Counting, rank and statistical regression tests |
| `wiki/` | Markdown source pages for the GitHub Wiki |
| `scripts/` | Wiki export helper |
| `provenance/` | Source checksums and code-change patch |

The `wiki/` folder supplies content for the separate GitHub Wiki; uploading the folder alone does not publish the Wiki. Repository documentation is readable immediately after upload.

## Citation and reuse

Until a public manuscript and archived software release are available, cite **K-tools**, this repository's URL, the exact commit or release used, and the access date. Add the manuscript citation once the preprint is published. No author list, DOI or release archive has been invented. See [Citation](docs/Citation.md).

A software reuse license has not yet been selected by the authors. See [LICENSE-STATUS](LICENSE-STATUS.md) before reusing or distributing the code.

## Support

Please open a GitHub issue with your Python version, input dimensions, parameters and a small reproducible example. See [CONTRIBUTING](CONTRIBUTING.md).
