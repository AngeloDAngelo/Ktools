# Demonstration checks

[Documentation](Documentation.md) · [Tutorial](Tutorial.md)

The bundled synthetic example demonstrates KEA, KRS and K-RBP with small inputs and fixed seeds. It is a software check, not a measure of biological prediction accuracy.

```bash
python examples/run_demo.py --out results/demo
```

Demo outputs are in [examples/demo_output](../examples/demo_output/):

- KEA recovers the injected `ACG` repeat signature (`ACG`, `CGA`, `GAC`).
- KRS selects eight words for the example target at `top_pct=10`.
- K-RBP evaluates three fictional profiles and assigns a positive Cliff's delta to the query-favoring profile.

The regression tests cover reference counting boundaries and overlaps, KRS ranks and ties, K-RBP effect direction, K-map, K-CDS and miRNA seed orientation:

```bash
python -m pip install -e '.[dev]'
python -m pytest
```
