"""Technical checks for the WORK IN PROGRESS method, not biological validation."""
import numpy as np
import pandas as pd
import pytest
from KEA import KEA

A, B = 'AACCGGT', 'TTGGCCA'

def make_analysis(tmp_path, records, name='ref.fa'):
    fasta = tmp_path/name
    fasta.write_text(''.join(f'>{key}\n{seq}\n' for key, seq in records))
    return KEA(str(tmp_path/'results'), str(fasta), str(fasta))

def clustered_records():
    records = []
    for i in range(30):
        seq = list('CG'*200)
        if i < 6:
            seq[40:47], seq[70:77] = A, B
        records.append((f'tx{i:02d}', ''.join(seq)))
    return records

def test_short_counts_randomization_and_reproducibility(tmp_path):
    records = clustered_records()
    instance = make_analysis(tmp_path, records)
    first = instance.cooccurrence([A,B], n_resamples=199, seed=123,
                                 save=False, store_null=True, verbose=False)
    row = first['pairs'].iloc[0]
    assert row.n_a == row.n_b == row.observed_pairs == row.n_cooccurring_transcripts == 6
    assert row.cooccurrence_significant and row.p_cooccurrence >= 1/200
    assert first['settings']['development_status'] == 'work_in_progress'
    reverse = make_analysis(tmp_path, list(reversed(records)), 'reverse.fa')
    repeated = reverse.cooccurrence([B,A], n_resamples=199, seed=123, save=False, verbose=False)
    pd.testing.assert_frame_equal(first['pairs'], repeated['pairs'])

def test_long_selects_pair_and_calibrates_all_distances(tmp_path):
    instance = make_analysis(tmp_path, clustered_records())
    result = instance.cooccurrence([A,B], mode='long', n_resamples=199,
                                  n_reference=29, seed=123, save=False, verbose=False)
    row = result['pairs'].iloc[0]
    assert row.distance_status == 'complete' and row.distance_significant
    dist = result['distributions'][(A,B)]
    assert dist.observed_count.sum() == 6
    assert dist.loc[dist.distance_nt.eq(30), 'observed_count'].iloc[0] == 6

def test_one_a_two_b_retains_both_distances(tmp_path):
    seq = list('CG'*150)
    seq[40:47], seq[60:67], seq[160:167] = A, B, B
    instance = make_analysis(tmp_path, [('tx',''.join(seq))])
    result = instance.cooccurrence([A,B], mode='long', distance_filter=False,
                                  n_resamples=29, n_reference=19, save=False, verbose=False)
    dist = result['distributions'][(A,B)].set_index('distance_nt')
    assert dist.observed_count.sum() == 2
    assert dist.loc[20,'observed_count'] == dist.loc[120,'observed_count'] == 1

def test_insufficient_sites_stays_missing(tmp_path):
    instance = make_analysis(tmp_path, [('small','AAAAAT'),('background','CG'*40)])
    result = instance.cooccurrence(['AAA','AAT'], mode='long', distance_filter=False,
                                  n_resamples=19, n_reference=9, save=False, verbose=False)
    row = result['pairs'].iloc[0]
    assert row.distance_status == 'insufficient_background_in_transcript'
    assert np.isnan(row.p_distance) and np.isnan(row.q_distance)

def test_saved_outputs_and_no_enrichment_rerun(tmp_path):
    instance = make_analysis(tmp_path, clustered_records())
    def unexpected(*args, **kwargs):
        raise AssertionError('Enrichment must not be rerun')
    instance.ExtractKmers = unexpected
    result = instance.cooccurrence([A,B], mode='long', n_resamples=49,
                                  n_reference=9, seed=123, save=True, verbose=False)
    from pathlib import Path
    import json
    out = Path(result['output_dir'])
    assert (out/'pair_tests.tsv').is_file() and (out/'positions.tsv.gz').is_file()
    settings = json.loads((out/'settings.json').read_text())
    assert settings['state'] == 'complete' and settings['development_status'] == 'work_in_progress'
    assert settings['source_unchanged']
    with np.load(result['distributions'][(A,B)]) as dist:
        assert dist['observed_count'].sum() == 6

def test_strict_enrichment_intersection_and_duplicate_ids(tmp_path):
    instance = make_analysis(tmp_path, clustered_records())
    instance.KEA_results = {'input': {'delta_median': {'enriched':[A,B]},
                                    'stat_log2fc': {'enriched':[]}}}
    with pytest.raises(ValueError, match='at least two'):
        instance.cooccurrence(save=False, verbose=False)
    duplicate = make_analysis(tmp_path, [('dup','ACGT'*20),('dup','ACGT'*20)], 'dup.fa')
    with pytest.raises(ValueError, match='Duplicate FASTA'):
        duplicate.cooccurrence([A,B], save=False, verbose=False)
