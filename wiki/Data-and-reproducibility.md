# Data and reproducibility

[Documentation](Documentation.md) · [Benchmarks](Benchmarks.md)

## Tutorial data

All tutorial inputs are simulated. The repository contains 24 reference and 48 control FASTA sequences generated with seed 42, a script to regenerate them, and a demonstration that creates a fictional RBP score matrix with seed 123. No biological case-study data are needed. Run:

```bash
python examples/run_demo.py --out results/demo
```

The demonstration is not the study's in silico benchmark. It is a small exercise covering the public API.

## Benchmark data

`benchmarks/reported_summary.tsv` contains approximate summary values reported in the supplied manuscript and explicit eCLIP/PEKA counts from the presentation. Original per-dataset measurements and benchmark pipelines are not included in the available source material.

To reproduce the in silico and eCLIP evaluations, the following still need to be deposited:

| Resource | Required content |
| --- | --- |
| Synthetic sequences | FASTA groups, injected motifs, interval ground truth and seeds |
| Comparator workflows | Scripts, commands, exact environments and logs |
| Performance tables | Per-dataset scores with valid, failed, missing and undefined status |
| Transcript resources | Exact Ensembl 99 sequence files, annotation and selection rules |
| eCLIP resources | ENCODE file accessions, transcript groups, peak and control intervals |
| RBP score matrices | PEKA source sheets, chosen profiles, preprocessing and checksums |
| Analysis notebooks | Complete benchmark evaluation and figure-generation code |

The PEKA resource cited in the benchmark methods is Additional file 5, Table S4 of [Kuret et al. (2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9461102/). The exact curated matrix used in the evaluation is not bundled.

## Run records

Record the code commit, package versions, platform, input checksums, annotation release, group membership, k, test, correction, thresholds and seeds. Preserve complete result tables as well as filtered signatures. Capture your environment with `python -m pip freeze > environment.txt`.

The current snapshot fixes a final-window counting error and stabilizes K-RBP background ordering. These changes can affect results; benchmark summaries have not been recomputed against this snapshot. See [CHANGELOG](../CHANGELOG.md).

## Provenance

`provenance/source-manifest.json` records checksums of the supplied source code and documentation inputs. `provenance/source-changes.patch` records the code edits. The draft manuscript and full working presentation are not distributed in the repository.
