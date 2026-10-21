import ast
import json
from pathlib import Path
import textwrap

import pandas as pd
import pytest
import KEA as module


def analysis(tmp_path, suffix=''):
    ref, ctrl = tmp_path/'ref.fa', tmp_path/'ctrl.fa'
    if not ref.exists():
        ref.write_text('>ref\nAAAAATCGATCG\n')
        ctrl.write_text('>ctrl\nATCGATATTTTT\n')
    obj = module.KEA(str(tmp_path/('results'+suffix)), str(ref), str(ctrl))
    obj.CreateCombinedFasta()
    return obj


def saved(obj, key='input'):
    return Path(obj.kmers_count_cache[key]['table'])


def test_original_counts_and_default_complete_table_saved(tmp_path):
    obj = analysis(tmp_path)
    from itertools import product
    words = [''.join(word) for word in product('ATCG', repeat=2)]
    expected = pd.DataFrame(0., index=pd.Index(words, name=0), columns=['ref','ctrl'])
    for word,ref,ctrl in [('AA',4,0),('AT',2,3),('TC',2,1),('CG',1,1),
                          ('GA',1,1),('TA',0,1),('TT',0,3)]:
        expected.loc[word] = [ref,ctrl]
    expected['k'] = 2
    assert obj.KmersCountsTable(k=2, cores=1) is None
    table = obj.kmers_count_table_by_species['input']
    pd.testing.assert_frame_equal(table, expected)
    assert table.shape == (16,3)
    assert saved(obj).parent == Path(obj.dir_out)
    pd.testing.assert_frame_equal(pd.read_pickle(saved(obj)), table)


def test_fresh_object_reuses_cache_without_pool_even_with_different_cores(tmp_path, monkeypatch):
    obj = analysis(tmp_path)
    obj.KmersCountsTable(k=2, cores=1)
    fresh = module.KEA(obj.dir_out, obj.reference_fasta, obj.control_fasta)
    def unexpected(*args, **kwargs):
        raise AssertionError('Reused cache must not create a worker pool')
    monkeypatch.setattr(module.mp, 'Pool', unexpected)
    fresh.KmersCountsTable(k=2, cores=15)
    assert fresh.kmers_count_cache['input']['state'] == 'reused'
    pd.testing.assert_frame_equal(obj.kmers_count_table_by_species['input'],
                                  fresh.kmers_count_table_by_species['input'])


def test_changed_fasta_same_size_invalidates_cache(tmp_path, monkeypatch):
    obj = analysis(tmp_path)
    obj.KmersCountsTable(k=2, cores=1)
    old = obj.kmers_count_table_by_species['input'].copy()
    p = Path(obj.combined_fasta)
    p.write_text(p.read_text().replace('AAAAATCGATCG','TTTTTTCGATCG'))
    real_pool = module.mp.Pool
    calls = []
    def counted(*args, **kwargs):
        calls.append(1)
        return real_pool(*args, **kwargs)
    monkeypatch.setattr(module.mp, 'Pool', counted)
    obj.KmersCountsTable(k=2, cores=1)
    assert calls == [1] and obj.kmers_count_cache['input']['state'] == 'computed_saved'
    assert not old.equals(obj.kmers_count_table_by_species['input'])


def test_k_filter_and_species_use_separate_caches(tmp_path):
    obj = analysis(tmp_path)
    obj.KmersCountsTable(k=2, cores=1)
    full = saved(obj)
    obj.KmersCountsTable(k=3, cores=1)
    other_k = saved(obj)
    obj.KmersCountsTable(k=2, cores=1, subset_by_subseq=['AA'])
    subset = saved(obj)
    assert obj.kmers_count_table_by_species['input'].index.tolist() == ['AA']
    obj.combined_fasta_by_species['mouse'] = obj.combined_fasta
    obj.KmersCountsTable(k=2, cores=1, species='mouse')
    assert len({full, other_k, subset, saved(obj,'mouse')}) == 4
    assert all(p.is_file() for p in [full, other_k, subset, saved(obj,'mouse')])


def test_force_and_cache_false_recount(tmp_path, monkeypatch):
    obj = analysis(tmp_path)
    obj.KmersCountsTable(k=2, cores=1)
    real_pool = module.mp.Pool
    calls = []
    def counted(*args, **kwargs):
        calls.append(1)
        return real_pool(*args, **kwargs)
    monkeypatch.setattr(module.mp, 'Pool', counted)
    obj.KmersCountsTable(k=2, cores=1, force=True)
    cache = saved(obj)
    stamp = cache.stat().st_mtime_ns
    obj.KmersCountsTable(k=2, cores=1, cache=False)
    assert calls == [1,1] and cache.stat().st_mtime_ns == stamp
    assert obj.kmers_count_cache['input']['state'] == 'computed_cache_disabled'


def test_corrupt_table_and_manifest_recover_by_recounting(tmp_path):
    obj = analysis(tmp_path)
    obj.KmersCountsTable(k=2, cores=1)
    expected = obj.kmers_count_table_by_species['input'].copy()
    saved(obj).write_bytes(b'broken')
    obj.KmersCountsTable(k=2, cores=1)
    pd.testing.assert_frame_equal(expected, obj.kmers_count_table_by_species['input'])
    Path(obj.kmers_count_cache['input']['metadata']).write_text('invalid JSON')
    obj.KmersCountsTable(k=2, cores=1)
    pd.testing.assert_frame_equal(expected, obj.kmers_count_table_by_species['input'])


def test_input_changed_during_counting_does_not_save_cache(tmp_path, monkeypatch):
    obj = analysis(tmp_path)
    real_pool = module.mp.Pool
    class ChangingPool:
        def __enter__(self):
            self.pool = real_pool(1)
            return self
        def __exit__(self,*args):
            return self.pool.__exit__(*args)
        def starmap(self, func, args):
            result = self.pool.starmap(func,args)
            p = Path(obj.combined_fasta)
            p.write_text(p.read_text()+'\n')
            return result
    monkeypatch.setattr(module.mp,'Pool',lambda *args, **kwargs: ChangingPool())
    with pytest.raises(RuntimeError,match='changed during counting'):
        obj.KmersCountsTable(k=2,cores=1)
    assert not list(Path(obj.dir_out).glob('*.pkl.gz'))


def test_recreated_identical_combined_fasta_keeps_cache_valid(tmp_path, monkeypatch):
    obj = analysis(tmp_path)
    obj.KmersCountsTable(k=2,cores=1)
    obj.CreateCombinedFasta()
    def unexpected(*args, **kwargs):
        raise AssertionError('Identical content remains cacheable after rewriting FASTA')
    monkeypatch.setattr(module.mp,'Pool',unexpected)
    obj.KmersCountsTable(k=2,cores=1)
    assert obj.kmers_count_cache['input']['state']=='reused'
