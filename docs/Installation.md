# Installation

[Documentation](Documentation.md) · [Tutorial](Tutorial.md)

Use a dedicated Python environment. Package metadata requires Python ≥3.10 and pandas ≥2.2,<3 because the implementation sets `future.no_silent_downcasting`. The distribution name is `ktools-rna`, while the import remains `from KEA import KEA`. This name is local packaging metadata; no PyPI publication is claimed.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
python examples/run_demo.py --out results/demo
```

Run these commands from the repository root. On Windows use `.venv\Scripts\activate`. For development use `python -m pip install -e '.[dev]'` and run `python -m pytest`.

Dependencies include NumPy, pandas, SciPy, scikit-learn, statsmodels, Biopython, matplotlib, seaborn, Logomaker, pybiomart, WordCloud and matplotlib-venn. Internet access is needed to install packages. The included demonstration runs offline after installation.

## Resources

At k=7 the count table has 16,384 rows. A dense float64 matrix with 20,000 transcript columns alone occupies approximately 2.44 GiB; DataFrame overhead, multiprocessing copies and intermediate rank tables add substantially to this. Start with `cores=1` and small inputs. Larger k increases the k-mer universe exponentially.

`KmersCountsTable` creates a process pool even with one worker. Run analysis in a `.py` script with an `if __name__ == "__main__":` guard. If notebook multiprocessing fails, run the provided script in a terminal.

## Optional advanced utilities

The core tutorial does not need BEDTools, genome mapping scripts, BigWig utilities or BioMart network access. BEDTools helpers require `bedtools` on PATH. For genome mapping, set `KTOOLS_MAPPER_DIR` to the mapper resource directory and `KTOOLS_MAPPER_SCRIPT` to the external mapper script. See [Known limitations](Known-limitations.md).

## Runtime library roles and versions

These constraints come from the package metadata; installation resolves compatible versions automatically.

| Library | Constraint | Role |
| :--- | :--- | :--- |
| NumPy | ≥1.24, <3 | Arrays, numerical operations and sampling |
| pandas | ≥2.2, <3 | Count, rank and result tables |
| SciPy | ≥1.11, <2 | Statistical tests |
| scikit-learn | ≥1.3, <2 | Clustering and analysis utilities |
| statsmodels | ≥0.14, <1 | Multiple-testing correction |
| Biopython | ≥1.81, <2 | FASTA and sequence handling |
| matplotlib | ≥3.7, <4 | Figures |
| seaborn | ≥0.13, <1 | Statistical plots |
| Logomaker | ≥0.8, <1 | Sequence logos |
| pybiomart | ≥0.2, <1 | Advanced BioMart annotation queries |
| WordCloud | ≥1.9, <2 | Word-cloud visualization |
| matplotlib-venn | ≥0.11, <2 | Venn diagrams |

For contributors, `python -m pip install -e '.[dev]'` adds pytest (≥8,<10) and build (≥1,<2). [Tested environment](../examples/tested-environment.txt) records versions from the earlier local validation; it is not a lockfile or a guarantee for every platform.

## Installation versus source execution

Run `python -m pip install .` from the extracted repository root. This builds and installs the `KEA` module plus its declared dependencies into the active environment; it is not merely adding the current folder to Python's path. `requirements.txt` contains `-e .`, so `python -m pip install -r requirements.txt` installs this same project in editable mode using the dependency declarations in `pyproject.toml`.

To check the installed copy rather than accidentally importing local `KEA.py`, run this from another directory:

```bash
python -c "import KEA; print(KEA.__file__)"
python -m pip check
```

The import path should point inside the active environment's `site-packages` for a regular install. The package name is `ktools-rna`; the Python module name is `KEA`. This project is not published to PyPI: use the local repository install command rather than `pip install ktools-rna`.

After upload, another user can download/extract the repository or clone it, create a virtual environment, and run the same installation commands. The repository's `docs`, figures and demo inputs remain in the downloaded repository; the miRNA family table is installed as package data, while documentation and example inputs remain in the repository.
