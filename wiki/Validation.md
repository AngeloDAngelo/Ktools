# Installation and validation

Verified on 30 September 2026 with Python 3.12 on macOS arm64, in a newly created virtual environment without inherited site-packages.

1. `python -m pip install ./outputs/ktools` downloaded and installed the declared runtime dependencies and built/installed the package wheel successfully.
2. `python -m pip check` reported `No broken requirements found`.
3. Import was checked from outside the repository: `KEA.__file__` resolved to the new environment's `site-packages/KEA.py`, confirming that the installed copy was exercised.
4. The complete synthetic demo ran with that installed copy: KEA recovered `ACG`, `CGA`, `GAC`; KRS selected 8 words; K-RBP evaluated 3 simulated profiles.
5. After installing pytest in the same environment, all four regression tests passed (terminal and overlapping counts, total counts, KRS absence/ties, K-RBP direction/correction).
6. Local Markdown links, Wiki export and ZIP integrity were checked.

[Dependency versions](../examples/tested-environment.txt) record this environment. [Expected outputs](../examples/expected/) provide a comparison for the demo. These tests cover the portable core workflow; manuscript benchmarks were not recomputed and advanced genomic utilities were not validated. Python 3.10/Linux checks are configured in GitHub Actions but were not executed by this local check. The original Desktop source was not edited.
