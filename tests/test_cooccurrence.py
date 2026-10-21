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
    returned = instance.cooccurrence([A,B], n_resamples=199, seed=123,
                                    save=False, store_null=True, verbose=False)
    assert isinstance(returned, pd.DataFrame)
    first = instance.cooccurrence_results
    row = first['pairs'].iloc[0]
    assert row.n_a == row.n_b == row.observed_pairs == row.n_cooccurring_transcripts == 6
    assert row.cooccurrence_significant and row.p_cooccurrence >= 1/200
    assert first['settings']['development_status'] == 'work_in_progress'
    reverse = make_analysis(tmp_path, list(reversed(records)), 'reverse.fa')
    repeated = reverse.cooccurrence([B,A], n_resamples=199, seed=123, save=False, verbose=False)
    pd.testing.assert_frame_equal(returned, repeated)

def test_long_selects_pair_and_calibrates_all_distances(tmp_path):
    instance = make_analysis(tmp_path, clustered_records())
    enrichment, distances = instance.cooccurrence([A,B], mode='cooccurrence_distance', n_resamples=199,
                                               n_reference=29, seed=123, save=False, verbose=False)
    assert isinstance(enrichment, pd.DataFrame) and isinstance(distances, pd.DataFrame)
    result = instance.cooccurrence_results
    row = result['pairs'].iloc[0]
    assert row.distance_status == 'complete' and row.distance_significant
    dist = result['distributions'][(A,B)]
    assert dist.observed_count.sum() == 6
    assert dist.loc[dist.distance_nt.eq(30), 'observed_count'].iloc[0] == 6

def test_one_a_two_b_retains_both_distances(tmp_path):
    seq = list('CG'*150)
    seq[40:47], seq[60:67], seq[160:167] = A, B, B
    instance = make_analysis(tmp_path, [('tx',''.join(seq))])
    enrichment, distances = instance.cooccurrence([A,B], mode='cooccurrence_distance', distance_filter=False,
                                               n_resamples=29, n_reference=19, save=False, verbose=False)
    row = distances.iloc[0]
    assert row.mean_distance_observed_nt == row.median_distance_observed_nt == 70
    assert np.isclose(row.sd_distance_observed_nt, np.std([20,120], ddof=1))
    result = instance.cooccurrence_results
    dist = result['distributions'][(A,B)].set_index('distance_nt')
    assert dist.observed_count.sum() == 2
    assert dist.loc[20,'observed_count'] == dist.loc[120,'observed_count'] == 1

def test_insufficient_sites_stays_missing(tmp_path):
    instance = make_analysis(tmp_path, [('small','AAAAAT'),('background','CG'*40)])
    instance.cooccurrence(['AAA','AAT'], mode='cooccurrence_distance', distance_filter=False,
                          n_resamples=19, n_reference=9, save=False, verbose=False)
    result = instance.cooccurrence_results
    row = result['pairs'].iloc[0]
    assert row.distance_status == 'insufficient_background_in_transcript'
    assert np.isnan(row.p_distance) and np.isnan(row.q_distance)

def test_saved_outputs_and_no_enrichment_rerun(tmp_path):
    instance = make_analysis(tmp_path, clustered_records())
    def unexpected(*args, **kwargs):
        raise AssertionError('Enrichment must not be rerun')
    instance.ExtractKmers = unexpected
    instance.cooccurrence([A,B], mode='cooccurrence_distance', n_resamples=49,
                          n_reference=9, seed=123, save=True, verbose=False)
    result = instance.cooccurrence_results
    from pathlib import Path
    import json
    out = Path(result['output_dir'])
    assert (out/'pair_tests.tsv').is_file() and (out/'positions.tsv.gz').is_file()
    pd.testing.assert_frame_equal(pd.read_csv(out/'cooccurrence.tsv', sep='\t'),result['cooccurrence'])
    assert (out/'cooccurrence_distance.tsv').is_file()
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

def test_unordered_pair_table_has_no_self_or_reversed_pairs(tmp_path):
    instance = make_analysis(tmp_path, clustered_records())
    third = 'GGATCCA'
    table = instance.cooccurrence([B,A,A,third], mode='cooccurrence',
                                  n_resamples=19, save=False, verbose=False)
    assert isinstance(table, pd.DataFrame) and len(table) == 3
    assert (table.kmer_a < table.kmer_b).all()
    assert not table[['kmer_a','kmer_b']].duplicated().any()
    assert {'observed_pairs','null_mean_pairs','null_median_pairs','null_sd_pairs',
            'pair_enrichment','p_cooccurrence','q_cooccurrence','cooccurrence_status'} <= set(table)

def test_distance_summaries_exist_even_when_pair_is_not_selected(tmp_path):
    instance = make_analysis(tmp_path, clustered_records())
    cooc, distance = instance.cooccurrence([A,B], mode='cooccurrence_distance',
                                         n_resamples=19, n_reference=9, alpha=.0001,
                                         save=False, verbose=False)
    row = distance.iloc[0]
    assert row.distance_status == 'not_selected_by_short_filter'
    assert row.distance_pairs == 6
    assert row.mean_distance_observed_nt == row.median_distance_observed_nt == 30
    assert row.sd_distance_observed_nt == 0
    assert np.isnan(row.p_distance) and np.isnan(row.mean_distance_null_nt)
    assert not instance.cooccurrence_results['distributions']

def test_pooled_random_summaries_match_all_saved_start_distances(tmp_path):
    seq = list('CG'*150)
    seq[40:47], seq[60:67], seq[160:167] = A, B, B
    instance = make_analysis(tmp_path, [('tx',''.join(seq))])
    cooc, distance = instance.cooccurrence([B,A], mode='cooccurence_distance',
                                         distance_filter=False, n_resamples=19,
                                         n_reference=9, save=False, verbose=False)
    hist = instance.cooccurrence_results['distributions'][(A,B)]
    expanded = np.repeat(hist.distance_nt, np.rint(hist.random_expected_count*9).astype(int))
    row = distance.iloc[0]
    assert len(expanded) == row.null_distance_pairs == 18
    assert np.isclose(row.mean_distance_null_nt, expanded.mean())
    assert row.median_distance_null_nt == np.median(expanded)
    assert np.isclose(row.sd_distance_null_nt, expanded.std(ddof=1))
    assert row.min_distance_observed_nt == 20 and row.max_distance_observed_nt == 120
    assert row.q25_distance_observed_nt == 45 and row.q75_distance_observed_nt == 95
    assert 0 < row.p_distance <= 1 and 0 <= row.distance_KS_statistic <= 1

def test_single_distance_has_missing_sample_sd_and_zero_pairs_have_stable_columns(tmp_path):
    seq = list('CG'*100)
    seq[20:27], seq[100:107] = A, B
    instance = make_analysis(tmp_path, [('tx',''.join(seq))])
    cooc, distance = instance.cooccurrence([A,B], mode='cooccurrence_distance',
                                         distance_filter=False, n_resamples=19,
                                         n_reference=9, save=False, verbose=False)
    row = distance.iloc[0]
    assert row.distance_pairs == 1 and row.mean_distance_observed_nt == 80
    assert np.isnan(row.sd_distance_observed_nt)
    absent = make_analysis(tmp_path, [('none','CG'*100)], 'none.fa')
    zero_cooc, zero_distance = absent.cooccurrence([A,B], mode='cooccurrence_distance',
                                                 n_resamples=19, save=False, verbose=False)
    assert list(distance.columns) == list(zero_distance.columns)
    row = zero_distance.iloc[0]
    assert row.distance_pairs == 0 and row.distance_status == 'zero_observed_pairs'
    assert np.isnan(row.mean_distance_observed_nt) and np.isnan(row.p_distance)
