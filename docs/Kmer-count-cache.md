# Automatic k-mer count-table cache

```python
analysis.KmersCountsTable(cores=15, k=7)
```

Disk caching is enabled by default. The complete pandas count table, with
all selected k-mer rows, transcript columns and the existing `k` column,
is saved inside `analysis.dir_out`, the directory supplied to the KEA
constructor. For input k=7, the files are:

- `kmers_counts_input_k7.pkl.gz`: compressed DataFrame; readable with
  `pandas.read_pickle(path)`.
- `kmers_counts_input_k7.pkl.gz.json`: cache identity and integrity manifest.

The next identical call loads the saved table without starting a counting
worker pool. Reuse also works with a new KEA object and a different core
count. Table content, index/column order and data types are preserved.

A cache is valid only for the same combined FASTA **content**, k, species,
ordered `subset_by_subseq` filter, counter implementation and numpy/pandas
versions. The compressed file checksum and table schema are checked.
Changed inputs or invalid/incomplete caches trigger recounting. Rewriting
an identical combined FASTA remains cacheable. Input changes during counting
or loading raise an error instead of saving/reusing inconsistent results.
Species, k values and subsequence filters use separate cache filenames.

If a new object's combined FASTA attribute is empty, the method recognizes
the standard `combined_sequences.fa` (or species-specific equivalent) already
present in the KEA directory. If it is absent, run `CreateCombinedFasta` first
as before. The actual combined FASTA used by the counter is the input checked
for cache identity.

```python
# Explicitly recount and replace this cache.
analysis.KmersCountsTable(cores=15, k=7, force=True)

# Recount without reading or writing a cache.
analysis.KmersCountsTable(cores=15, k=7, cache=False)

table = analysis.kmers_count_table_by_species['input']
print(analysis.kmers_count_cache['input'])
```

The existing `None` return value and counting algorithm are unchanged.
Only `KmersCountsTable` changes; other KEA methods remain unchanged. Eight technical tests cover exact counts, fresh-object reuse,
content-based invalidation, separate cache keys, override flags, corrupt
caches, input mutation and identical FASTA rewrites.
