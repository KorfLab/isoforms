"""
general test for testing what percentile range shall be for hints cutoff
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
    description='Test geniso percentile cutoff performance with complexity and score metrics')
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

# percentile values to test
parser.add_argument('--percentiles', type=int, nargs='+',
    default=[30, 20, 10],
    help='percentile cutoff values to test [%(default)s]')

# baseline addition control
parser.add_argument('--add_baseline', action='store_true',
    help='add GT-AG baseline sites to HMM hints')
parser.add_argument('--baseline_isoform_threshold', type=int, default=5000,
    metavar='<int>', help='isoform count threshold for adding baseline sites [%(default)i]')

# constraint parameters
parser.add_argument('--min_intron', type=int, default=35,
    metavar='<int>', help='minimum intron length [%(default)i]')
parser.add_argument('--min_exon', type=int, default=25,
    metavar='<int>', help='minimum exon length [%(default)i]')
parser.add_argument('--flank', type=int, default=99,
    metavar='<int>', help='genomic flank [%(default)i]')
parser.add_argument('--limit', type=int, default=100,
    metavar='<int>', help='limit isoforms [%(default)i]')

# Score threshold for high-quality isoforms
parser.add_argument('--score_percentile', type=float, default=75.0,
    metavar='<float>', help='percentile threshold for high-scoring isoforms [%(default).1f]')

# Output file for summary
parser.add_argument('--summary', type=str, default='summary.txt',
    metavar='<file>', help='summary output file [%(default)s]')

arg = parser.parse_args()

def calc_complexity(locus):
    if not locus or not locus.isoforms:
        return 0.0
    return isoform.complexity(locus.isoforms)

def get_isoform_sigs(locus):
    sigs = set()
    for iso in locus.isoforms:
        sig = (tuple(iso.dons), tuple(iso.accs))
        sigs.add(sig)
    return sigs

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

def calc_retention(base_sigs, test_sigs):
    if not base_sigs:
        return 0.0
    
    common = len(base_sigs & test_sigs)
    total = len(base_sigs)
    
    return (common / total) * 100

def calc_score_distribution_similarity(base_locus, test_locus):
    if not base_locus or not test_locus:
        return 0.0
    if not base_locus.isoforms or not test_locus.isoforms:
        return 0.0
    
    # Get probability distributions
    base_probs = sorted([iso.prob for iso in base_locus.isoforms], reverse=True)
    test_probs = sorted([iso.prob for iso in test_locus.isoforms], reverse=True)
    
    # Normalize to compare top N isoforms (where N is min of both)
    n = min(len(base_probs), len(test_probs))
    if n == 0:
        return 0.0
    
    # Calculate similarity using sum of min probabilities
    similarity = 0
    for i in range(n):
        similarity += min(base_probs[i], test_probs[i])
    
    return similarity * 100  # Convert to percentage

def get_out_of_bag_sites(base_dons, base_accs, hmm_dons, hmm_accs):
    """
    Get sites that are in baseline but not in HMM predictions (out-of-bag sites).
    These should be preserved when applying HMM hints.
    """
    oob_dons = [d for d in base_dons if d not in hmm_dons] if hmm_dons else base_dons
    oob_accs = [a for a in base_accs if a not in hmm_accs] if hmm_accs else base_accs
    return oob_dons, oob_accs

def merge_sites(hmm_sites, oob_sites):
    """
    Merge HMM sites with out-of-bag sites and sort.
    """
    combined = list(set(hmm_sites + oob_sites))
    return sorted(combined)

def test_fasta(fasta, model, percentiles):
    fname = Path(fasta).stem
    
    constraints = {
        'min_intron': arg.min_intron,
        'min_exon':   arg.min_exon,
        'flank':      arg.flank
    }
    
    name, seq = next(isoform.read_fasta(fasta))

    # Run baseline (no hints) with timing
    base_start = time.perf_counter()
    try:
        base_locus = isoform.Locus(
            name, seq, model,
            constraints=constraints,
            limit=arg.limit
        )
    except Exception as e:
        return None
    base_runtime = time.perf_counter() - base_start
    
    if not base_locus or not base_locus.isoforms:
        return None
    
    # Calculate baseline metrics
    base_sigs = get_isoform_sigs(base_locus)
    base_complexity = calc_complexity(base_locus)
    base_high_score_sigs = get_high_score_isoforms(base_locus, arg.score_percentile)
    
    results = {
        'name':              fname,
        'length':            len(base_locus.seq),
        'base_dons':         len(base_locus.dons),
        'base_accs':         len(base_locus.accs),
        'base_iso':          len(base_locus.isoforms),
        'base_complexity':   base_complexity,
        'base_high_score':   len(base_high_score_sigs),
        'base_runtime':      base_runtime,
        'percentiles':       {}
    }
    
    # Run HMM only once (and time it separately)
    hmm_start = time.perf_counter()
    try:
        hmm_output = isohint.run_hmm(arg.hmm, fasta, arg.hmm_model)
        raw_dons, raw_accs = isohint.parse_hint(hmm_output)
    except Exception as e:
        return results  # Return baseline results only
    hmm_runtime = time.perf_counter() - hmm_start
    results['hmm_runtime'] = hmm_runtime
    
    # Determine if we should add baseline sites
    oob_dons = []
    oob_accs = []
    if arg.add_baseline:
        # Count potential isoforms with HMM sites only
        hmm_iso_count = isohint.countiso(
            raw_dons, raw_accs,
            arg.min_intron, arg.min_exon,
            limit=arg.baseline_isoform_threshold + 1
        )
        
        # If below threshold, add baseline sites
        if hmm_iso_count < arg.baseline_isoform_threshold:
            oob_dons, oob_accs = get_out_of_bag_sites(
                base_locus.dons, base_locus.accs,
                raw_dons, raw_accs
            )
    
    # Test each percentile cutoff
    for p in percentiles:
        try:
            # Apply percentile cutoff (this is very fast, but we can time it)
            cutoff_start = time.perf_counter()
            dons_cut = isohint.percentile(raw_dons, p) if raw_dons else []
            accs_cut = isohint.percentile(raw_accs, p) if raw_accs else []
            
            # Merge with out-of-bag sites if using baseline addition
            if arg.add_baseline and (oob_dons or oob_accs):
                dons_cut = merge_sites(dons_cut, oob_dons)
                accs_cut = merge_sites(accs_cut, oob_accs)
            
            cutoff_time = time.perf_counter() - cutoff_start
            
            # Create locus with the cutoff hints (this is the main computation)
            geniso_start = time.perf_counter()
            test_locus = isoform.Locus(
                name, seq, model,
                constraints=constraints,
                limit=arg.limit,
                dons=dons_cut,
                accs=accs_cut
            )
            geniso_runtime = time.perf_counter() - geniso_start
            
            if test_locus and test_locus.isoforms:
                test_sigs = get_isoform_sigs(test_locus)
                test_complexity = calc_complexity(test_locus)
                test_high_score_sigs = get_high_score_isoforms(test_locus, arg.score_percentile)
                
                # Calculate all metrics
                retention = calc_retention(base_sigs, test_sigs)
                high_score_retention = calc_retention(base_high_score_sigs, test_high_score_sigs)
                score_similarity = calc_score_distribution_similarity(base_locus, test_locus)
                
                results['percentiles'][p] = {
                    'dons':                len(test_locus.dons),
                    'accs':                len(test_locus.accs),
                    'isoforms':            len(test_locus.isoforms),
                    'retention':           retention,
                    'high_score_retention': high_score_retention,
                    'score_similarity':    score_similarity,
                    'complexity':          test_complexity,
                    'complexity_ratio':    test_complexity / base_complexity if base_complexity > 0 else 0,
                    'don_reduction':       100 * (len(base_locus.dons) - len(test_locus.dons)) / len(base_locus.dons) if len(base_locus.dons) > 0 else 0,
                    'acc_reduction':       100 * (len(base_locus.accs) - len(test_locus.accs)) / len(base_locus.accs) if len(base_locus.accs) > 0 else 0,
                    'geniso_runtime':      geniso_runtime,
                    'runtime_ratio':       geniso_runtime / base_runtime if base_runtime > 0 else 0,
                    'cutoff_time':         cutoff_time,
                    'speedup':             base_runtime / geniso_runtime if geniso_runtime > 0 else 0,
                }
            else:
                results['percentiles'][p] = None
                
        except Exception as e:
            results['percentiles'][p] = None
    
    return results

def write_data_stream_header(percentiles):
    """Write header information for the data stream output"""
    print('# Geniso Percentile Cutoff Performance Test - Raw Data Stream')
    print(f'# Date: {time.strftime("%Y-%m-%d %H:%M:%S")}')
    print(f'# Parameters: min_intron={arg.min_intron} min_exon={arg.min_exon} flank={arg.flank} limit={arg.limit}')
    print(f'# High-score threshold: {arg.score_percentile}th percentile')
    print(f'# Percentiles tested: {percentiles}')
    if arg.add_baseline:
        print(f'# Baseline addition: ENABLED (threshold: {arg.baseline_isoform_threshold} isoforms)')
    else:
        print(f'# Baseline addition: DISABLED')
    print('#')
    print('# Column Definitions:')
    print('#   Col 1: gene        - Name of the gene/sequence')
    print('#   Col 2: length      - Gene sequence length (bp)')
    print('#   Col 3: base_don    - Baseline donor sites')
    print('#   Col 4: base_acc    - Baseline acceptor sites')
    print('#   Col 5: base_iso    - Baseline isoforms')
    print('#   Col 6: base_cplx   - Baseline complexity score')
    print('#   Col 7: base_rt     - Baseline runtime (s)')
    print('#   Col 8: hmm_rt      - HMM runtime (s)')
    
    # Add columns for each percentile
    col_num = 9
    for p in percentiles:
        print(f'# --- Percentile {p}% columns ---')
        print(f'#   Col {col_num}: p{p}_don     - Donor sites at {p}%')
        print(f'#   Col {col_num+1}: p{p}_acc     - Acceptor sites at {p}%')
        print(f'#   Col {col_num+2}: p{p}_iso     - Isoforms at {p}%')
        print(f'#   Col {col_num+3}: p{p}_ret     - Retention % at {p}%')
        print(f'#   Col {col_num+4}: p{p}_hsr     - High-score retention % at {p}%')
        print(f'#   Col {col_num+5}: p{p}_cplx    - Complexity ratio at {p}%')
        print(f'#   Col {col_num+6}: p{p}_speed   - Speedup factor at {p}%')
        print(f'#   Col {col_num+7}: p{p}_rt      - Runtime at {p}% (s)')
        col_num += 8
    
    print('#')
    # NO HEADER ROW - removed the header line printing

def write_data_stream_row(result, percentiles):
    """Write a single row of data to the stream with fixed width formatting"""
    if not result:
        return
    
    # Start with basic data - use fixed width formatting
    row = []
    row.append(f"{result['name']:<20}")  # 20 chars for gene name
    row.append(f"{result['length']:>8}")  # 8 chars for length
    row.append(f"{result['base_dons']:>8}")
    row.append(f"{result['base_accs']:>8}")
    row.append(f"{result['base_iso']:>8}")
    row.append(f"{result['base_complexity']:>10.4f}")
    row.append(f"{result['base_runtime']:>10.4f}")
    row.append(f"{result.get('hmm_runtime', 0):>10.4f}")
    
    # Add data for each percentile with fixed width
    for p in percentiles:
        if p in result['percentiles'] and result['percentiles'][p]:
            pd = result['percentiles'][p]
            row.append(f"{pd['dons']:>8}")
            row.append(f"{pd['accs']:>8}")
            row.append(f"{pd['isoforms']:>8}")
            row.append(f"{pd['retention']:>8.2f}")
            row.append(f"{pd['high_score_retention']:>8.2f}")
            row.append(f"{pd['complexity_ratio']:>8.4f}")
            row.append(f"{pd['speedup']:>8.2f}")
            row.append(f"{pd['geniso_runtime']:>10.4f}")
        else:
            # No data for this percentile - use consistent width
            row.extend([f"{'NA':>8}"] * 7 + [f"{'NA':>10}"])
    
    print(' '.join(row))

def write_summary_to_file(all_results, percentiles, filename):
    """Write the summary report to a file"""
    
    with open(filename, 'w') as f:
        f.write('='*70 + '\n')
        f.write('GENISO PERCENTILE CUTOFF PERFORMANCE TEST - SUMMARY REPORT\n')
        f.write('='*70 + '\n\n')
        
        f.write(f'Analysis Date: {time.strftime("%Y-%m-%d %H:%M:%S")}\n')
        f.write(f'Parameters:\n')
        f.write(f'  Min intron: {arg.min_intron}\n')
        f.write(f'  Min exon: {arg.min_exon}\n')
        f.write(f'  Flank: {arg.flank}\n')
        f.write(f'  Limit: {arg.limit}\n')
        f.write(f'  High-score threshold: {arg.score_percentile}th percentile\n')
        f.write(f'  Percentiles tested: {percentiles}\n')
        if arg.add_baseline:
            f.write(f'  Baseline addition: ENABLED (threshold: {arg.baseline_isoform_threshold})\n')
        else:
            f.write(f'  Baseline addition: DISABLED\n')
        f.write(f'  Files processed: {len([r for r in all_results if r])}\n\n')
        
        # Calculate average baseline and HMM runtimes
        base_runtimes = [r['base_runtime'] for r in all_results if r and 'base_runtime' in r]
        hmm_runtimes = [r['hmm_runtime'] for r in all_results if r and 'hmm_runtime' in r]
        
        if base_runtimes:
            avg_base_runtime = sum(base_runtimes) / len(base_runtimes)
            f.write(f'Average baseline geniso runtime: {avg_base_runtime:.4f}s\n')
        if hmm_runtimes:
            avg_hmm_runtime = sum(hmm_runtimes) / len(hmm_runtimes)
            f.write(f'Average HMM computation time: {avg_hmm_runtime:.4f}s (one-time cost)\n')
        
        f.write('\n' + '-'*70 + '\n')
        f.write('AVERAGE PERFORMANCE ACROSS ALL FILES\n')
        f.write('-'*70 + '\n\n')
        f.write('Percentile | Retention% | HighScore% | Complex% | Runtime(s) | Speedup | Don_red% | Acc_red%\n')
        f.write('-----------|------------|------------|----------|------------|---------|----------|----------\n')
        
        # Calculate averages for each percentile
        for p in percentiles:
            retentions = []
            high_score_rets = []
            complex_ratios = []
            don_reducs = []
            acc_reducs = []
            runtimes = []
            speedups = []
            
            for res in all_results:
                if res and res['percentiles'].get(p):
                    pd = res['percentiles'][p]
                    retentions.append(pd['retention'])
                    high_score_rets.append(pd['high_score_retention'])
                    complex_ratios.append(pd['complexity_ratio'] * 100)
                    don_reducs.append(pd['don_reduction'])
                    acc_reducs.append(pd['acc_reduction'])
                    runtimes.append(pd['geniso_runtime'])
                    speedups.append(pd['speedup'])
            
            if retentions:
                avg_ret = sum(retentions) / len(retentions)
                avg_hsr = sum(high_score_rets) / len(high_score_rets) if high_score_rets else 0
                avg_cr = sum(complex_ratios) / len(complex_ratios) if complex_ratios else 0
                avg_dr = sum(don_reducs) / len(don_reducs) if don_reducs else 0
                avg_ar = sum(acc_reducs) / len(acc_reducs) if acc_reducs else 0
                avg_rt = sum(runtimes) / len(runtimes) if runtimes else 0
                avg_sp = sum(speedups) / len(speedups) if speedups else 0
                
                f.write(f'    {p:3d}%    |   {avg_ret:6.1f}   |   {avg_hsr:6.1f}   |  {avg_cr:6.1f}  |   {avg_rt:7.4f}  |  {avg_sp:5.2f}x  |  {avg_dr:6.1f}  |  {avg_ar:6.1f}\n')
        
        f.write('\n' + '='*70 + '\n')
        f.write('END OF SUMMARY REPORT\n')
        f.write('='*70 + '\n')

def main():
    # Find test files
    test_dir = Path(arg.test_dir)
    if not test_dir.exists():
        print(f"Error: test directory {test_dir} not found", file=sys.stderr)
        sys.exit(1)
    
    # Look for fasta files
    fastas = []
    
    # First check the main directory
    fastas.extend(sorted(test_dir.glob('*.fa')))
    
    # If no files found, check subdirectories
    if not fastas:
        for subdir in test_dir.iterdir():
            if subdir.is_dir():
                fa_files = sorted(subdir.glob('*.fa'))
                if fa_files:
                    fastas.extend(fa_files)
    
    # If still no files, try recursive search
    if not fastas:
        fastas = sorted(test_dir.rglob('*.fa'))
    
    if not fastas:
        print(f"Error: no fasta files found in {test_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Limit files if requested
    if arg.max_files > 0:
        fastas = fastas[:arg.max_files]
    
    print(f"# Testing {len(fastas)} files with percentiles: {arg.percentiles}", file=sys.stderr)
    print(f"# High-score threshold: {arg.score_percentile}th percentile", file=sys.stderr)
    if arg.add_baseline:
        print(f"# Baseline addition: ENABLED (threshold: {arg.baseline_isoform_threshold})", file=sys.stderr)
    else:
        print(f"# Baseline addition: DISABLED", file=sys.stderr)
    print(f"# Summary will be written to: {arg.summary}", file=sys.stderr)
    print(f"#", file=sys.stderr)
    
    # Load model once
    try:
        if os.path.isdir(arg.model):
            model_files = list(Path(arg.model).glob('*.splicemodel'))
            if not model_files:
                print(f"Error: No .splicemodel files found in {arg.model}", file=sys.stderr)
                sys.exit(1)
            model_file = str(model_files[0])
            print(f"# Using splice model: {model_file}", file=sys.stderr)
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
    
    # Write data stream header
    write_data_stream_header(arg.percentiles)
    
    # Test each file
    all_results = []
    start_time = time.time()
    
    for i, fasta in enumerate(fastas, 1):
        # Progress indicator to stderr
        if i % 10 == 0 or i == len(fastas):
            elapsed = time.time() - start_time
            rate = i / elapsed if elapsed > 0 else 0
            print(f"# Progress: {i}/{len(fastas)} ({100*i/len(fastas):.0f}%) - {rate:.1f} files/sec", file=sys.stderr)
        
        result = test_fasta(str(fasta), model, arg.percentiles)
        if result:
            all_results.append(result)
            # Write raw data to stdout immediately
            write_data_stream_row(result, arg.percentiles)
    
    # Write summary to file
    if all_results:
        print(f"\n# Analysis complete. Writing summary to {arg.summary}", file=sys.stderr)
        write_summary_to_file(all_results, arg.percentiles, arg.summary)
        print(f"# Summary written successfully", file=sys.stderr)
        
        # Also print brief stats to stderr
        valid_results = len([r for r in all_results if r])
        print(f"# Processed {valid_results} files successfully", file=sys.stderr)

if __name__ == '__main__':
    main()