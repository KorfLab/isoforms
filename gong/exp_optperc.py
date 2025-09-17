"""
Test what's the optimal percentile range for gene
"""

import argparse
import sys
import os
from pathlib import Path
import time
import math
from collections import defaultdict

import isoform
import isohint

parser = argparse.ArgumentParser(
    description='Analyze gene length vs optimal percentile cutoff relationships')
parser.add_argument('model', type=str, metavar='<splice model file>',
    help='splice model file for isoform (e.g., worm.splicemodel)')
parser.add_argument('hmm', type=str, metavar='<hmm executable>',
    help='path to hmm executable')
parser.add_argument('--hmm_model', type=str, metavar='<hmm model dir>',
    default=None,
    help='hmm model directory containing PWM/MM files (optional)')
parser.add_argument('--test_dir', type=str, metavar='<dir>',
    default='/Users/gongchen/Code/isoforms/gong/test_data/smallgenes',
    help='directory with test fasta files [%(default)s]')
parser.add_argument('--max_files', type=int, default=0,
    metavar='<int>', help='max files to process (0=all) [%(default)i]')
parser.add_argument('--target_quality', type=float, default=90.0,
    metavar='<float>', help='target quality retention percentage [%(default).1f]')
parser.add_argument('--percentile_step', type=int, default=5,
    metavar='<int>', help='step size for percentile search [%(default)i]')
parser.add_argument('--min_percentile', type=int, default=0,
    metavar='<int>', help='minimum percentile to start testing [%(default)i]')
parser.add_argument('--max_percentile', type=int, default=50,
    metavar='<int>', help='maximum percentile to test [%(default)i]')
parser.add_argument('--use_oob', action='store_true',
    help='add out-of-bag sites (baseline sites not in HMM) to HMM hints')
parser.add_argument('--min_intron', type=int, default=35,
    metavar='<int>', help='minimum intron length [%(default)i]')
parser.add_argument('--min_exon', type=int, default=25,
    metavar='<int>', help='minimum exon length [%(default)i]')
parser.add_argument('--flank', type=int, default=99,
    metavar='<int>', help='genomic flank [%(default)i]')
parser.add_argument('--limit', type=int, default=100,
    metavar='<int>', help='limit isoforms [%(default)i]')
parser.add_argument('--score_percentile', type=float, default=75.0,
    metavar='<float>', help='percentile threshold for high-scoring isoforms [%(default).1f]')
parser.add_argument('--min_isoforms', type=int, default=10000,
    metavar='<int>', help='minimum potential isoforms required to test gene [%(default)i]')
parser.add_argument('--summary', type=str, default='summary.txt',
    metavar='<file>', help='summary output file path [%(default)s]')

arg = parser.parse_args()

def calc_complexity(locus):
    if not locus or not locus.isoforms:
        return 0.0
    return isoform.complexity(locus.isoforms)

def get_high_score_isoforms(locus, percentile=75):
    if not locus or not locus.isoforms:
        return set()
    
    scores = [iso.score for iso in locus.isoforms]
    if not scores:
        return set()
    
    # Calculate percentile threshold
    sorted_scores = sorted(scores)
    idx = int(len(sorted_scores) * (percentile / 100))
    if idx >= len(sorted_scores):
        idx = len(sorted_scores) - 1
    threshold = sorted_scores[idx]
    
    # Get signatures of high-scoring isoforms
    high_score_sigs = set()
    for iso in locus.isoforms:
        if iso.score >= threshold:
            sig = (tuple(iso.dons), tuple(iso.accs))
            high_score_sigs.add(sig)
    
    return high_score_sigs

def get_isoform_sigs(locus):
    sigs = set()
    for iso in locus.isoforms:
        sig = (tuple(iso.dons), tuple(iso.accs))
        sigs.add(sig)
    return sigs

def calc_retention(base_sigs, test_sigs):
    if not base_sigs:
        return 0.0
    
    common = len(base_sigs & test_sigs)
    total = len(base_sigs)
    
    return (common / total) * 100

def get_out_of_bag_sites(gtag_dons, gtag_accs, hmm_dons, hmm_accs):
    oob_dons = [d for d in gtag_dons if d not in hmm_dons] if hmm_dons else []
    oob_accs = [a for a in gtag_accs if a not in hmm_accs] if hmm_accs else []
    return oob_dons, oob_accs

def merge_sites(hmm_sites, oob_sites):
    combined = list(set(hmm_sites + oob_sites))
    return sorted(combined)

def find_optimal_percentile(fasta, model, target_quality):    
    fname = Path(fasta).stem
    constraints = {
        'min_intron': arg.min_intron,
        'min_exon':   arg.min_exon,
        'flank':      arg.flank
    }
    
    name, seq   = next(isoform.read_fasta(fasta))
    gene_length = len(seq)
    
    # Get GT-AG sites
    gtag_dons, gtag_accs = isoform.gtag_sites(seq, arg.flank, arg.min_exon)
    
    # Check if gene has enough potential isoforms before processing
    potential_isoforms = isohint.countiso(
        gtag_dons, gtag_accs,
        arg.min_intron, arg.min_exon,
        limit=arg.min_isoforms + 1  # Check if > min_isoforms
    )
    
    if potential_isoforms <= arg.min_isoforms:
        return {
            'name': fname,
            'length': gene_length,
            'status': 'skipped',
            'potential_isoforms': potential_isoforms,
            'reason': f'Too few isoforms ({potential_isoforms} < {arg.min_isoforms})'
        }
    
    # Run baseline geniso with no hints
    try:
        base_locus = isoform.Locus(
            name, seq, model,
            constraints=constraints,
            limit=arg.limit
        )
    except Exception as e:
        return {
            'name': fname,
            'length': gene_length,
            'status': 'failed_baseline',
            'reason': str(e)
        }
    
    if not base_locus or not base_locus.isoforms:
        return {
            'name': fname,
            'length': gene_length,
            'status': 'no_baseline_isoforms',
            'reason': 'No isoforms generated in baseline'
        }
    
    # Get baseline metrics
    base_high_score_sigs    = get_high_score_isoforms(base_locus, arg.score_percentile)
    base_complexity         = calc_complexity(base_locus)
    base_isoform_sigs       = get_isoform_sigs(base_locus)
    
    # Run HMM once
    try:
        hmm_output = isohint.run_hmm(arg.hmm, fasta, arg.hmm_model)
        raw_hmm_dons, raw_hmm_accs = isohint.parse_hint(hmm_output)
        
        # Check if HMM produced valid output
        if not raw_hmm_dons and not raw_hmm_accs:
            return {
                'name': fname,
                'length': gene_length,
                'optimal_percentile': 'N/A',
                'high_score_retention': 100.0,
                'all_retention': 100.0,
                'complexity_retained': 100.0,
                'base_donors': len(base_locus.dons),
                'base_acceptors': len(base_locus.accs),
                'base_isoforms': len(base_locus.isoforms),
                'base_complexity': base_complexity,
                'test_donors': len(base_locus.dons),
                'test_acceptors': len(base_locus.accs),
                'test_isoforms': len(base_locus.isoforms),
                'donor_reduction': 0.0,
                'acceptor_reduction': 0.0,
                'status': 'no_hmm'
            }
            
    except Exception as e:
        return {
            'name': fname,
            'length': gene_length,
            'optimal_percentile': 'N/A',
            'high_score_retention': 100.0,
            'all_retention': 100.0,
            'complexity_retained': 100.0,
            'base_donors': len(base_locus.dons),
            'base_acceptors': len(base_locus.accs),
            'base_isoforms': len(base_locus.isoforms),
            'base_complexity': base_complexity,
            'test_donors': len(base_locus.dons),
            'test_acceptors': len(base_locus.accs),
            'test_isoforms': len(base_locus.isoforms),
            'donor_reduction': 0.0,
            'acceptor_reduction': 0.0,
            'status': 'hmm_failed'
        }
    
    # Get out-of-bag sites if requested
    oob_dons = []
    oob_accs = []
    if arg.use_oob:
        oob_dons, oob_accs = get_out_of_bag_sites(
            gtag_dons, gtag_accs,
            raw_hmm_dons, raw_hmm_accs
        )
    
    # Test percentiles from min to max
    optimal_percentile = None
    optimal_metrics = None
    
    for p in range(arg.min_percentile, arg.max_percentile + 1, arg.percentile_step):
        try:
            # Apply percentile cutoff to HMM sites
            if p == 0:
                # Use all HMM sites (no cutoff)
                hmm_dons_cut = raw_hmm_dons
                hmm_accs_cut = raw_hmm_accs
            else:
                hmm_dons_cut = isohint.percentile(raw_hmm_dons, p) if raw_hmm_dons else []
                hmm_accs_cut = isohint.percentile(raw_hmm_accs, p) if raw_hmm_accs else []
            
            # Merge with out-of-bag sites if using OOB
            if arg.use_oob and (oob_dons or oob_accs):
                final_dons = merge_sites(hmm_dons_cut, oob_dons)
                final_accs = merge_sites(hmm_accs_cut, oob_accs)
            else:
                final_dons = hmm_dons_cut
                final_accs = hmm_accs_cut
            
            # Create locus with hints
            test_locus = isoform.Locus(
                name, seq, model,
                constraints=constraints,
                limit=arg.limit,
                dons=final_dons,
                accs=final_accs
            )
            
            if test_locus and test_locus.isoforms:
                test_high_score_sigs = get_high_score_isoforms(test_locus, arg.score_percentile)
                test_complexity = calc_complexity(test_locus)
                test_isoform_sigs = get_isoform_sigs(test_locus)
                
                # Calculate metrics
                high_score_retention = calc_retention(base_high_score_sigs, test_high_score_sigs)
                all_isoform_retention = calc_retention(base_isoform_sigs, test_isoform_sigs)
                complexity_ratio = (test_complexity / base_complexity * 100) if base_complexity > 0 else 0
                
                # Site reduction
                don_reduction = 100 * (len(base_locus.dons) - len(test_locus.dons)) / len(base_locus.dons) if len(base_locus.dons) > 0 else 0
                acc_reduction = 100 * (len(base_locus.accs) - len(test_locus.accs)) / len(base_locus.accs) if len(base_locus.accs) > 0 else 0
                
                # Check if we meet the target quality
                if high_score_retention >= target_quality:
                    # This percentile is good, save it as optimal
                    optimal_percentile = p
                    optimal_metrics = {
                        'high_score_retention': high_score_retention,
                        'all_retention': all_isoform_retention,
                        'complexity_ratio': complexity_ratio,
                        'don_reduction': don_reduction,
                        'acc_reduction': acc_reduction,
                        'test_dons': len(test_locus.dons),
                        'test_accs': len(test_locus.accs),
                        'test_isoforms': len(test_locus.isoforms)
                    }
                else:
                    # Quality dropped below target, stop here
                    break
                    
        except Exception as e:
            continue
    
    # If we never found a good percentile, report the best we got
    if optimal_percentile is None:
        optimal_percentile = 'None'
        optimal_metrics = {
            'high_score_retention': 0,
            'all_retention': 0,
            'complexity_ratio': 0,
            'don_reduction': 0,
            'acc_reduction': 0,
            'test_dons': len(base_locus.dons),
            'test_accs': len(base_locus.accs),
            'test_isoforms': len(base_locus.isoforms)
        }
    
    return {
        'name': fname,
        'length': gene_length,
        'optimal_percentile': optimal_percentile,
        'high_score_retention': optimal_metrics['high_score_retention'],
        'all_retention': optimal_metrics['all_retention'],
        'complexity_retained': optimal_metrics['complexity_ratio'],
        'base_donors': len(base_locus.dons),
        'base_acceptors': len(base_locus.accs),
        'base_isoforms': len(base_locus.isoforms),
        'base_complexity': base_complexity,
        'test_donors': optimal_metrics['test_dons'],
        'test_acceptors': optimal_metrics['test_accs'],
        'test_isoforms': optimal_metrics['test_isoforms'],
        'donor_reduction': optimal_metrics['don_reduction'],
        'acceptor_reduction': optimal_metrics['acc_reduction'],
        'status': 'valid'
    }

def write_results_header():
    print(f'# Gene Length vs Optimal Percentile Analysis')
    print(f'# Date: {time.strftime("%Y-%m-%d %H:%M:%S")}')
    print(f'# Target quality: {arg.target_quality}% high-score isoform retention')
    print(f'# Percentile search range: {arg.min_percentile}-{arg.max_percentile}% (step={arg.percentile_step})')
    print(f'# High-score threshold: {arg.score_percentile}th percentile')
    print(f'# Minimum isoforms required: {arg.min_isoforms}')
    if arg.use_oob:
        print(f'# Out-of-bag sites: ENABLED')
    else:
        print(f'# Out-of-bag sites: DISABLED')
    print('#')
    print('# Column Definitions:')
    print('#   Col 1:  gene       - Name of the gene/sequence')
    print('#   Col 2:  length     - Gene sequence length (bp)')
    print('#   Col 3:  opt_pct    - Optimal percentile cutoff')
    print('#   Col 4:  hi_ret%    - High-score isoform retention %')
    print('#   Col 5:  all_ret%   - All isoform retention %')
    print('#   Col 6:  cplx_ret%  - Complexity retained %')
    print('#   Col 7:  base_don   - Baseline donor sites')
    print('#   Col 8:  base_acc   - Baseline acceptor sites')
    print('#   Col 9:  base_iso   - Baseline isoforms')
    print('#   Col 10: base_cplx  - Baseline complexity')
    print('#   Col 11: test_don   - Test donor sites')
    print('#   Col 12: test_acc   - Test acceptor sites')
    print('#   Col 13: test_iso   - Test isoforms')
    print('#   Col 14: don_red%   - Donor reduction %')
    print('#   Col 15: acc_red%   - Acceptor reduction %')
    print('#   Col 16: status     - Analysis status')
    print('#')

def write_results_row(result):
    if not result:
        return
    
    # Format with fixed widths for alignment
    row = []
    row.append(f"{result['name']:<20}")
    row.append(f"{result['length']:>8}")
    row.append(f"{str(result['optimal_percentile']):>8}")
    row.append(f"{result['high_score_retention']:>8.1f}")
    row.append(f"{result['all_retention']:>9.1f}")
    row.append(f"{result['complexity_retained']:>10.1f}")
    row.append(f"{result['base_donors']:>9}")
    row.append(f"{result['base_acceptors']:>9}")
    row.append(f"{result['base_isoforms']:>9}")
    row.append(f"{result['base_complexity']:>10.3f}")
    row.append(f"{result['test_donors']:>9}")
    row.append(f"{result['test_acceptors']:>9}")
    row.append(f"{result['test_isoforms']:>9}")
    row.append(f"{result['donor_reduction']:>9.1f}")
    row.append(f"{result['acceptor_reduction']:>9.1f}")
    row.append(f"{result.get('status', 'unknown'):>9}")
    
    print(' '.join(row))

def write_summary(results, filename):
    
    valid_results   = [r for r in results if r and r.get('status') == 'valid']
    successful      = [r for r in valid_results if r['high_score_retention'] >= arg.target_quality]
    no_hmm          = [r for r in results if r and r.get('status') in ['no_hmm', 'hmm_failed']]
    failed          = [r for r in results if r and r.get('status') in ['failed_baseline', 'no_baseline_isoforms']]
    skipped         = [r for r in results if r and r.get('status') == 'skipped']
    
    with open(filename, 'w') as f:
        f.write('='*70 + '\n')
        f.write('GENE LENGTH VS OPTIMAL PERCENTILE ANALYSIS SUMMARY\n')
        f.write('='*70 + '\n\n')
        
        f.write(f'Analysis Date: {time.strftime("%Y-%m-%d %H:%M:%S")}\n')
        f.write(f'Parameters:\n')
        f.write(f'  Target quality: {arg.target_quality}%\n')
        f.write(f'  Percentile range: {arg.min_percentile}-{arg.max_percentile}% (step={arg.percentile_step})\n')
        f.write(f'  Minimum isoforms: {arg.min_isoforms}\n')
        f.write(f'  Out-of-bag sites: {"ENABLED" if arg.use_oob else "DISABLED"}\n')
        f.write(f'  Min intron: {arg.min_intron}, Min exon: {arg.min_exon}\n\n')
        
        f.write('-'*70 + '\n')
        f.write('PROCESSING SUMMARY\n')
        f.write('-'*70 + '\n')
        f.write(f'Total genes analyzed: {len(results)}\n')
        f.write(f'  Skipped (< {arg.min_isoforms} isoforms): {len(skipped)}\n')
        if skipped and len(skipped) <= 5:
            for r in skipped:
                f.write(f'    {r["name"]}: {r.get("potential_isoforms", "?")} isoforms\n')
        elif skipped:
            for r in skipped[:3]:
                f.write(f'    {r["name"]}: {r.get("potential_isoforms", "?")} isoforms\n')
            f.write(f'    ... and {len(skipped) - 3} more\n')
        f.write(f'  Valid results: {len(valid_results)}\n')
        f.write(f'  No HMM output: {len(no_hmm)}\n')
        f.write(f'  Failed baseline: {len(failed)}\n')
        
        if valid_results:
            f.write(f'\nQuality Target Achievement:\n')
            f.write(f'  Meeting {arg.target_quality}% target: {len(successful)}/{len(valid_results)} ')
            f.write(f'({100*len(successful)/len(valid_results):.1f}%)\n')
        
        # Linear regression analysis
        if len(successful) > 1:
            f.write('\n' + '-'*70 + '\n')
            f.write('LINEAR REGRESSION ANALYSIS\n')
            f.write('-'*70 + '\n')
            
            lengths = [r['length'] for r in successful]
            percentiles = [r['optimal_percentile'] for r in successful if isinstance(r['optimal_percentile'], (int, float))]
            
            if len(percentiles) == len(lengths):
                n = len(lengths)
                mean_x = sum(lengths) / n
                mean_y = sum(percentiles) / n
                
                num = sum((x - mean_x) * (y - mean_y) for x, y in zip(lengths, percentiles))
                den = sum((x - mean_x) ** 2 for x in lengths)
                
                if den > 0:
                    slope = num / den
                    intercept = mean_y - slope * mean_x
                    
                    # Calculate R-squared
                    ss_tot = sum((y - mean_y) ** 2 for y in percentiles)
                    ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(lengths, percentiles))
                    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
                    
                    f.write(f'Relationship: optimal_percentile = {slope:.6f} * length + {intercept:.2f}\n')
                    f.write(f'  R-squared: {r_squared:.3f}\n')
                    f.write(f'  Sample size: {n} genes\n\n')
                    
                    # Predicted values
                    f.write(f'Predicted optimal percentiles:\n')
                    for length in [500, 1000, 2000, 5000, 10000]:
                        pred = slope * length + intercept
                        pred = max(arg.min_percentile, min(arg.max_percentile, pred))
                        f.write(f'  {length:5d} bp: {pred:.1f}%\n')
        
        # Performance metrics
        if successful:
            f.write('\n' + '-'*70 + '\n')
            f.write('PERFORMANCE METRICS\n')
            f.write('-'*70 + '\n')
            
            avg_retention = sum(r['high_score_retention'] for r in successful) / len(successful)
            avg_don_red = sum(r['donor_reduction'] for r in successful) / len(successful)
            avg_acc_red = sum(r['acceptor_reduction'] for r in successful) / len(successful)
            
            f.write(f'Average metrics for genes meeting target:\n')
            f.write(f'  High-score retention: {avg_retention:.1f}%\n')
            f.write(f'  Donor site reduction: {avg_don_red:.1f}%\n')
            f.write(f'  Acceptor site reduction: {avg_acc_red:.1f}%\n')
        
        f.write('\n' + '='*70 + '\n')
        f.write('END OF SUMMARY\n')
        f.write('='*70 + '\n')

def main():
    # Find test files
    test_dir = Path(arg.test_dir)
    if not test_dir.exists():
        print(f"Error: test directory {test_dir} not found", file=sys.stderr)
        sys.exit(1)
    
    # Look for fasta files
    fastas = []
    fastas.extend(sorted(test_dir.glob('*.fa')))
    if not fastas:
        for subdir in test_dir.iterdir():
            if subdir.is_dir():
                fa_files = sorted(subdir.glob('*.fa'))
                if fa_files:
                    fastas.extend(fa_files)
    if not fastas:
        fastas = sorted(test_dir.rglob('*.fa'))
    
    if not fastas:
        print(f"Error: no fasta files found in {test_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Limit files if requested
    if arg.max_files > 0:
        fastas = fastas[:arg.max_files]
    
    # Load model once
    try:
        if os.path.isdir(arg.model):
            model_files = list(Path(arg.model).glob('*.splicemodel'))
            if not model_files:
                print(f"Error: No .splicemodel files found in {arg.model}", file=sys.stderr)
                sys.exit(1)
            model_file = str(model_files[0])
            model = isoform.read_splicemodel(model_file)
        else:
            if not os.path.exists(arg.model):
                print(f"Error: Model file not found at {arg.model}", file=sys.stderr)
                sys.exit(1)
            model = isoform.read_splicemodel(arg.model)
    except Exception as e:
        print(f"Error loading model: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Verify HMM executable
    if not os.path.exists(arg.hmm):
        print(f"Error: HMM executable not found at {arg.hmm}", file=sys.stderr)
        sys.exit(1)
    
    # Write header
    write_results_header()
    
    # Analyze each file
    all_results = []
    skipped_count = 0
    processed_count = 0
    
    for i, fasta in enumerate(fastas, 1):
        result = find_optimal_percentile(
            str(fasta), 
            model, 
            arg.target_quality
        )
        
        if result:
            all_results.append(result)
            # Only write row if not skipped
            if result.get('status') != 'skipped':
                write_results_row(result)
                processed_count += 1
            else:
                skipped_count += 1
    
    # Write summary to file
    if all_results:
        print(f'#')
        print(f'# Processed: {processed_count} genes, Skipped: {skipped_count} genes')
        print(f'# Summary written to: {arg.summary}')
        write_summary(all_results, arg.summary)

if __name__ == '__main__':
    main()