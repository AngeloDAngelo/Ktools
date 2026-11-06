def test_spacing_fdr_and_input_preservation(tmp_path):
    from pathlib import Path
    import importlib.util,hashlib,inspect
    import numpy as np,pandas as pd
    from statsmodels.stats.multitest import multipletests
    source=Path(__file__).resolve().parents[1]/'KEA.py'
    spec=importlib.util.spec_from_file_location('spacing_v11_test',source);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    root=tmp_path
    ref=root/'ref.fa';ref.write_text('>tx1\nAAATGCCCTAACCGGATCNNNAAACCCGTACGTACGT\n>tx2\nCCCTACGATCGAAACGTCGATCCCGATCGTAAA\n')
    sha=hashlib.sha256(ref.read_bytes()).hexdigest()
    obj=m.KEA(str(root),str(ref),str(ref))
    kwargs=dict(kmers=['AAA','CCC','GGG'],n_resamples=29,n_background=10000,window_width=5,min_distance=4,min_enrichment=1.,pvalue_threshold=.99,min_support_fraction=0.,seed=9,verbose=False)
    base=obj.cooccurrence(**kwargs,plot=False)
    assert 'pair_FDR' not in base['pair_summary'].columns
    corrected=obj.cooccurrence(**kwargs,fdr_threshold=.99,plot=True)
    p=corrected['pair_summary'];valid=p.pair_scan_pvalue.notna()
    q=multipletests(p.pair_scan_pvalue.fillna(1),method='fdr_bh')[1]
    assert np.allclose(p.loc[valid,'pair_FDR'],q[valid])
    assert p.loc[~valid,'pair_FDR'].isna().all()
    wb=base['branches']['unmasked']['windows'];wc=corrected['branches']['unmasked']['windows']
    assert np.array_equal(wb.pvalue_scan,wc.pvalue_scan)
    assert np.array_equal(wc.selected,wc.selected_pre_fdr&wc.passes_pair_fdr)
    assert np.array_equal(wc.selected_pre_fdr,wb.selected)
    assert corrected['settings']['fdr_computed'] and corrected['settings']['fdr_pair_family_size']==3
    assert not (wc.selected&~wb.selected).any()
    print('BH test PASS: full family, NaN for untestable pairs, identical local p-values, and monotone filtering.',flush=True)
    # Strict FDR rejects the complete family at attainable MC resolution.
    strict=obj.cooccurrence(**kwargs,fdr_threshold=.0001,plot=False)
    assert not strict['branches']['unmasked']['windows'].selected.any()
    assert strict['branches']['unmasked']['regions'].empty
    assert strict['branches']['unmasked']['support_transcripts'].empty
    print('Strict FDR test PASS: rejected pairs clear regions and supporting memberships.',flush=True)
    for overrides in [dict(pairs=[('AAA','CCC')],fdr_threshold=.05),dict(kmers=['AAA','CCC'],fdr_threshold=.05),dict(fdr_threshold=0),dict(fdr_threshold=1),dict(fdr_method='invalid')]:
     bad=kwargs|overrides
     try:obj.cooccurrence(**bad,plot=False)
     except ValueError:pass
     else:raise AssertionError('Invalid FDR input accepted: '+str(overrides))
    print('Validation test PASS: partial family, single-pair FDR, invalid threshold/method rejected.',flush=True)
    by=obj.cooccurrence(**kwargs,fdr_threshold=.99,fdr_method='fdr_by',plot=False)
    pby=by['pair_summary'];valid=pby.pair_scan_pvalue.notna();expected=multipletests(pby.pair_scan_pvalue.fillna(1),method='fdr_by')[1]
    assert np.allclose(pby.loc[valid,'pair_FDR'],expected[valid])
    assert hashlib.sha256(ref.read_bytes()).hexdigest()==sha
    print('BY and input preservation PASS.',flush=True)
    print('Signature:',inspect.signature(m.KEA.cooccurrence),flush=True)
