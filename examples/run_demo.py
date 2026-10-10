"""Offline KEA, KRS and K-RBP tutorial using only synthetic inputs."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO
from KEA import KEA


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="results/demo")
    args = parser.parse_args()
    out = Path(args.out)
    data = Path(__file__).resolve().parent / "data"
    reference, control = data / "reference.fa", data / "control.fa"
    records = list(SeqIO.parse(reference, "fasta")) + list(SeqIO.parse(control, "fasta"))
    assert len({r.id for r in records}) == len(records), "Duplicate IDs"
    assert all(set(str(r.seq)) <= set("ACGT") and len(r.seq) >= 3 for r in records)
    np.random.seed(42)
    analysis = KEA(str(out), str(reference), str(control))
    analysis.CreateCombinedFasta()
    analysis.KmersCountsTable(k=3, cores=1)
    selected = analysis.ExtractKmers(
        k_selected=3,
        extraction_methods=("delta_median", "stat_log2fc"),
        n_iter=100,
        sign_emp_pval=0.10,
        log2fc_threshold=0.5,
        correction_method="fdr",
        fdr_alpha=0.05,
    )["input"]
    assert np.isfinite(selected["delta_median"]["all_kmers"]["delta_median"]).all()
    signature = sorted(set(selected["delta_median"]["enriched"])
                       & set(selected["stat_log2fc"]["enriched"]))
    assert "ACG" in signature, "The injected motif should be recovered"
    (out / "kea_signature.txt").write_text("".join(k + "\n" for k in signature))

    # KRS is independent of KEA: characterize one fictional transcript.
    krs = analysis.KRS(reference="reference_006", top_pct=10)
    (out / "krs_signature.txt").write_text("".join(k + "\n" for k in krs))

    # Fictional 3-mer scores. The query-favoring profile is intentionally
    # constructed to demonstrate the API, not discover a biological RBP.
    kmers = sorted(analysis.kmers_count_table_by_species["input"].index)
    rng = np.random.default_rng(123)
    matrix = pd.DataFrame({
        "SIMULATED_query_favoring": [8 + rng.normal(0, 0.2) if k in signature
                                     else rng.normal(0, 0.2) for k in kmers],
        "SIMULATED_neutral": rng.normal(0, 1, len(kmers)),
        "SIMULATED_query_disfavoring": [-8 + rng.normal(0, 0.2) if k in signature
                                        else rng.normal(0, 0.2) for k in kmers],
    }, index=pd.Index(kmers, name="kmer"))
    matrix.to_csv(out / "synthetic_rbp_matrix.tsv", sep="\t")
    background = [word for word in kmers if word not in signature][:10 * len(signature)]
    rbp = analysis.KmersProteinEnrichment(signature, matrix, bg_kmers=background, plot=False,
                                          convert_to_rna=False, mode="continuous")
    rbp.to_csv(out / "synthetic_rbp_results.tsv", sep="\t", index=False)
    favored = rbp.set_index("protein").loc["SIMULATED_query_favoring"]
    assert favored["zscore_diff"] > 0.9
    summary = {"data": "synthetic demonstration only", "seed_sequences": 42,
               "seed_kea": 42, "seed_matrix": 123, "k": 3,
               "reference_transcripts": 24, "control_transcripts": 48,
               "kea_signature": signature, "krs_size": len(krs),
               "rbp_profiles_tested": len(rbp)}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"Outputs: {out.resolve()}")


if __name__ == "__main__":
    main()
