import numpy as np

import pandas as pd

from KEA import KEA, kmers_counter

def test_counts_match_reference_boundary_and_overlaps():
    assert kmers_counter("AAAA", 2).loc["AA", 1] == 2
    assert kmers_counter("ACG", 3).empty
    result = kmers_counter("ACGT", 3)
    assert result.loc["ACG", 1] == 1
    assert "CGT" not in result.index
    assert kmers_counter("AC", 3).empty

def test_complete_counts_match_number_of_windows():
    from itertools import product
    universe = ["".join(k) for k in product("ATCG", repeat=2)]
    counts = kmers_counter("ACGTACGT", 2, universe, "tx")
    assert counts["tx"].sum() == 6
    assert len(counts) == 16
    assert counts.loc["AA", "tx"] == 0

def test_krbp_effect_direction_and_correction(tmp_path):
    from itertools import product
    analysis = KEA(str(tmp_path), "unused.fa", "unused.fa")
    kmers = ["".join(k) for k in product("ACGT", repeat=2)]
    query = kmers[:4]
    scores = pd.DataFrame({"favored": [10] * 4 + [0] * 12,
                           "disfavored": [0] * 4 + [10] * 12}, index=kmers)
    result = analysis.KmersProteinEnrichment(query, scores, plot=False).set_index("protein")
    assert result.loc["favored", "zscore_diff"] == 1
    assert result.loc["disfavored", "zscore_diff"] == -1
    assert result.loc["favored", "pvalue"] < 0.01
    assert result.loc["disfavored", "pvalue"] > 0.95
    assert (result["fdr"] >= result["pvalue"]).all()
