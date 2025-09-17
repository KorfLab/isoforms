#!/usr/bin/env python3

import argparse
import os
import random

parser = argparse.ArgumentParser(description="Test HMM with fasta input")
parser.add_argument("fasta", help="Path to input fasta file")
parser.add_argument("--hmm", default="edhmm2", help="Path to HMM executable (default: edhmm2)")
parser.add_argument("--models", help="Path to models directory", default="../models")
parser.add_argument('--min_intron', required=False, type=int, default=35,
    metavar='<int>', help='minimum length of intron [%(default)i]')
parser.add_argument('--min_exon', required=False, type=int, default=25,
    metavar='<int>', help='minimum length exon [%(default)i]')
parser.add_argument('--flank', required=False, type=int, default=99,
    metavar='<int>', help='genomic flank on each side [%(default)i]')

args = parser.parse_args()

import isohint
import randomf
import isoform

from randomf import IsoformTree

def randf_rnaseq_trail(trails: int = 1, lower_pct: float = 0.1, upper_pct: float = 0.3):
    ''' running geniso + random forest '''
    
    models_dir = os.path.abspath(args.models)
    
    if not os.path.exists(args.fasta):
        print(f"Error: Fasta file '{args.fasta}' not found")
        return None, None
    
    if not os.path.exists(models_dir):
        print(f"Error: Models directory '{models_dir}' not found")
        return None, None
    
    output          = isohint.run_hmm(args.hmm, args.fasta, models_dir)
    hints           = isohint.parse_hint_rf(output)
    name, seq       = next(isoform.read_fasta(args.fasta))
    bdons, baccs    = randomf.basis_sample(seq, hints, args.flank, args.min_exon)
    
    all_dons, all_accs = isoform.gtag_sites(seq, args.flank, args.min_exon)
    total_dons = len(all_dons)
    total_accs = len(all_accs)
    
    dons, accs, pos2info = randomf.prepare_rf_input(hints)

    isotree = IsoformTree(dons, accs, pos2info)
    
    pred_dons, pred_accs = randomf.collect_oob_predictions(
        isotree.oob_samples, 
        isotree.rules, 
        isotree.output,
        total_dons,
        total_accs,
        lower_pct,
        upper_pct
    )
    
    final_dons = sorted(list(set(bdons + pred_dons)))
    final_accs = sorted(list(set(baccs + pred_accs)))
    
    return final_dons, final_accs

def main():
    '''
    Main function to run the pipeline
    '''
    
    # Run the random forest pipeline
    final_dons, final_accs = randf_rnaseq_trail(
        trails=1, 
        lower_bound=5, 
        upper_bound=50
    )
    
    if final_dons is not None and final_accs is not None:
        print(f"Pipeline completed successfully!")
        print(f"Total donor sites: {len(final_dons)}")
        print(f"Total acceptor sites: {len(final_accs)}")
        
        # You can now use final_dons and final_accs with geniso
        # Example: geniso_output = isohint.run_geniso2(geniso_path, args.fasta, model, hints=True)
        
    else:
        print("Pipeline failed!")

if __name__ == "__main__":
    exit(main())