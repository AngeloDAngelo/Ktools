# Examples

## Synthetic KEA, KRS and K-RBP demo

From the repository root after installation:

```bash
python examples/run_demo.py --out results/demo
```

The demo uses simulated sequences and fictional RBP profiles. It writes count/selection tables, signatures, association statistics and `summary.json`. See the [tutorial](../docs/Tutorial.md) and [expected outputs](expected/).

To regenerate the inputs:

```bash
python examples/generate_data.py
```

Generation uses seed 42: 24 reference and 48 control sequences, each 600 nt long. Reference sequences contain a variable-length `ACG` repeat. The demo counts 3-mers and generates fictional RBP scores with seed 123.

## K_miR: CDR1as / miR-7

The [CDR1as example](../docs/CDR1as-miRNA-example.md) includes a public circular RNA sequence, seed annotations, KRS scores and a figure. Supply a comparison FASTA to run it:

```bash
python examples/run_cdr1as_mirna.py --background /path/to/lncRNA.fa --out results/cdr1as_mirna
```
