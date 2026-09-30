# Synthetic tutorials

All inputs in this directory are simulated. No transcript, signature, RBP profile or experimental dataset from the manuscript is used.

From the repository root after installation:

```bash
python examples/run_demo.py --out results/demo
```

The script runs KEA, KRS and K-RBP. It writes full tables and a compact `summary.json`. The [step-by-step tutorial](../docs/Tutorial.md) explains each call.

To regenerate the bundled FASTA files deterministically:

```bash
python examples/generate_data.py
```

Generation uses Python's `random.Random(42)`. Each group contains sequences of length 600 nt: 24 reference and 48 control. A variable-length central segment of each reference sequence is replaced with an `ACG` repeat. This simple composition shift is intentionally easy to detect at k=3 and is not the manuscript's motif-injection benchmark.

The demonstration constructs its RBP score matrix from a separate seed (123), with fictional profiles deliberately favoring, ignoring or disfavoring the recovered query. It demonstrates software behavior and does not validate RBP prediction accuracy.
