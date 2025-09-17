"""
general test for different cutoff methods for interpreting HMM model output
"""

import argparse
import sys
import os
from pathlib import Path

import isoform
import isohint

parser = argparse.ArgumentParser(
    description='Test geniso cutoff methods performance with basis sites')
parser.add_argument('model', type=str, metavar='<splice model file>',
    help='splice model file for isoform (e.g., worm.splicemodel)')
parser.add_argument('hmm', type=str, metavar='<hmm executable>',
    help='path to hmm executable')
parser.add_argument('--hmm_model', type=str, metavar='<hmm model dir>',
    default=None,
    help='hmm model directory containing PWM/MM files (optional, hmm has defaults)')
parser.add_argument('--test_dir', type=str, metavar='<dir>',
    default='/Users/gongchen/Code/isoforms/gong/test_data/smallgenes',
    help='directory with test fasta files [%(default)s]')
parser.add_argument('--max_files', type=int, default=0,
    metavar='<int>', help='max files to process (0=all) [%(default)i]')
parser.add_argument('--min_intron', type=int, default=35,
    metavar='<int>', help='minimum intron length [%(default)i]')
parser.add_argument('--min_exon', type=int, default=25,
    metavar='<int>', help='minimum exon length [%(default)i]')
parser.add_argument('--flank', type=int, default=99,
    metavar='<int>', help='genomic flank [%(default)i]')
parser.add_argument('--limit', type=int, default=100,
    metavar='<int>', help='limit isoforms [%(default)i]')

# cutoff method parameters
parser.add_argument('--k', type=int, default=5,
    help='k for smooth/topk [%(default)i]')
parser.add_argument('--p', type=int, default=25,
    help='percentile for perc [%(default)i]')
parser.add_argument('--t', type=float, default=-5.0,
    help='threshold for thresh [%(default).1f]')
parser.add_argument('--min_sites', type=int, default=5,
    help='min sites for adapt [%(default)i]')

# basis sites option
parser.add_argument('--include_basis', action='store_true',
    help='include out-of-bag GT/AG sites not found by HMM')

arg = parser.parse_args()

def get_isoform_sigs(locus):
    sigs = set()
    for iso in locus.isoforms:
        sig = (tuple(iso.dons), tuple(iso.accs))
        sigs.add(sig)
    
    return sigs

def calc_overlap(base_sigs, test_sigs):    
    if not base_sigs:
        return 0.0
    
    common = len(base_sigs & test_sigs)
    total  = len(base_sigs)
    
    return (common / total) * 100

def test_fasta(fasta, model):
    fname      = Path(fasta).stem
    methods    = ['gap', 'smooth', 'perc', 'topk', 'thresh', 'adapt']
    
    constraints = {
        'min_intron': arg.min_intron,
        'min_exon':   arg.min_exon,
        'flank':      arg.flank
    }
    
    name, seq = next(isoform.read_fasta(fasta))
    
    base_locus = isoform.Locus(
        name, seq, model,
        constraints=constraints,
        limit=arg.limit
    )
    
    if not base_locus or not base_locus.isoforms:
        return None
    
    base_sigs = get_isoform_sigs(base_locus)
    
    results = {
        'name':        fname,
        'length':      len(base_locus.seq),
        'base_dons':   len(base_locus.dons),
        'base_accs':   len(base_locus.accs),
        'base_iso':    len(base_locus.isoforms),
        'base_reject': base_locus.rejected,
        'methods':     {}
    }
    
    hmm_output          = isohint.run_hmm(arg.hmm, fasta, arg.hmm_model)
    raw_dons, raw_accs  = isohint.parse_hint(hmm_output)
    
    basis_dons = []
    basis_accs = []
    if arg.include_basis:
        basis_dons, basis_accs = isohint.get_basis_sites(
            seq, raw_dons, raw_accs, 
            flank=arg.flank, minex=arg.min_exon
        )
        results['basis_dons'] = len(basis_dons)
        results['basis_accs'] = len(basis_accs)
    
    for method in methods:
        method_kwargs = {}
        if   method in ['smooth', 'topk']:
            method_kwargs['k'] = arg.k
        elif method == 'perc':
            method_kwargs['p'] = arg.p
        elif method == 'thresh':
            method_kwargs['t'] = arg.t
        elif method == 'adapt':
            method_kwargs['min_sites'] = arg.min_sites
        
        cutoff_methods = {
            'gap':    isohint.gapstats,
            'smooth': lambda s: isohint.smoothed_gapstats(s, method_kwargs.get('k', 5)),
            'perc':   lambda s: isohint.percentile(s, method_kwargs.get('p', 25)),
            'topk':   lambda s: isohint.topk(s, method_kwargs.get('k', 10)),
            'thresh': lambda s: isohint.threshold(s, method_kwargs.get('t', -5.0)),
            'adapt':  lambda s: isohint.adaptive_cutoff(s, method_kwargs.get('min_sites', 5)),
        }
        
        cutfn = cutoff_methods[method]
        dons_cut = cutfn(raw_dons) if raw_dons else []
        accs_cut = cutfn(raw_accs) if raw_accs else []
        
        if arg.include_basis:
            dons_cut = isohint.combine_with_basis(dons_cut, basis_dons)
            accs_cut = isohint.combine_with_basis(accs_cut, basis_accs)
        
        test_locus = isoform.Locus(
            name, seq, model,
            constraints=constraints,
            limit=arg.limit,
            dons=dons_cut,
            accs=accs_cut
        )
        
        if test_locus and test_locus.isoforms:
            test_sigs = get_isoform_sigs(test_locus)
            overlap   = calc_overlap(base_sigs, test_sigs)
            
            results['methods'][method] = {
                'dons':     len(test_locus.dons),
                'accs':     len(test_locus.accs),
                'isoforms': len(test_locus.isoforms),
                'rejected': test_locus.rejected,
                'overlap':  overlap,
                'maxprob':  test_locus.isoforms[0].prob,
                'minprob':  test_locus.isoforms[-1].prob
            }
        else:
            results['methods'][method] = None
    
    return results

def print_report(all_results):
    
    print('# Geniso cutoff methods performance test')
    if arg.include_basis:
        print('# INCLUDING out-of-bag GT/AG sites (basis sites)')
    print(f'# Parameters: min_intron={arg.min_intron} min_exon={arg.min_exon} flank={arg.flank} limit={arg.limit}')
    print(f'# Method params: k={arg.k} p={arg.p}% t={arg.t} min_sites={arg.min_sites}')
    print(f'# Files processed: {len([r for r in all_results if r])}')
    
    if arg.include_basis:
        basis_don_counts = [r['basis_dons'] for r in all_results if r and 'basis_dons' in r]
        basis_acc_counts = [r['basis_accs'] for r in all_results if r and 'basis_accs' in r]
        if basis_don_counts:
            avg_basis_dons = sum(basis_don_counts) / len(basis_don_counts)
            avg_basis_accs = sum(basis_acc_counts) / len(basis_acc_counts)
            print(f'# Average basis sites (not in HMM): donors={avg_basis_dons:.1f}, acceptors={avg_basis_accs:.1f}')
    
    print('#')
    
    methods = ['gap', 'smooth', 'perc', 'topk', 'thresh', 'adapt']
    
    print('# Average performance across all files:')
    print('# Method   Overlap%  Don_reduce%  Acc_reduce%  Iso_ratio%')
    print('# -------  --------  -----------  -----------  ----------')
    
    for m in methods:
        overlaps  = []
        don_reduc = []
        acc_reduc = []
        iso_ratio = []
        
        for res in all_results:
            if res and res['methods'].get(m):
                md = res['methods'][m]
                overlaps.append(md['overlap'])
                
                if res['base_dons'] > 0:
                    don_reduc.append(100 * (res['base_dons'] - md['dons']) / res['base_dons'])
                if res['base_accs'] > 0:
                    acc_reduc.append(100 * (res['base_accs'] - md['accs']) / res['base_accs'])
                if res['base_iso'] > 0:
                    iso_ratio.append(100 * md['isoforms'] / res['base_iso'])
        
        if overlaps:
            avg_ovr = sum(overlaps) / len(overlaps)
            avg_dr  = sum(don_reduc) / len(don_reduc) if don_reduc else 0
            avg_ar  = sum(acc_reduc) / len(acc_reduc) if acc_reduc else 0
            avg_ir  = sum(iso_ratio) / len(iso_ratio) if iso_ratio else 0
            
            print(f'# {m:8s}   {avg_ovr:6.1f}     {avg_dr:7.1f}      {avg_ar:7.1f}      {avg_ir:7.1f}')

def main():
    test_dir = Path(arg.test_dir)
    fastas   = sorted(test_dir.glob('*.fa'))
    
    if not fastas:
        for subdir in test_dir.iterdir():
            if subdir.is_dir():
                fa_files = sorted(subdir.glob('*.fa'))
                if fa_files:
                    fastas.extend(fa_files)
    
    if arg.max_files > 0:
        fastas = fastas[:arg.max_files]
    
    if os.path.isdir(arg.model):
        model_files = list(Path(arg.model).glob('*.splicemodel'))
        model_file = str(model_files[0])
        model = isoform.read_splicemodel(model_file)
    else:
        model = isoform.read_splicemodel(arg.model)
    
    all_results = []
    
    for i, fasta in enumerate(fastas, 1):
        result = test_fasta(str(fasta), model)
        if result:
            all_results.append(result)
    
    if all_results:
        print_report(all_results)

if __name__ == '__main__':
    main()