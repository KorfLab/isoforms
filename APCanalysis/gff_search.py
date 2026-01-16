import argparse
import gzip

parser = argparse.ArgumentParser(
	description='get overlapping introns in genomic region of interest'
parser.add_argument('annotation', help='gff3 genome annotation')
parser.add_argument('coordinates', help='chromosome:loc i.e. X:99:999; '
	'chromosomes: I, II, III, IV, V, X, MtDNA')
parser.add_argument('--seq_desc', help='add description for sequence '
	'in first line of fasta and file name')

args = parser.parse_args()

'''
dyn-1 whole gene sequence is X:15568921:15573021
for APC, only use 15571571 to 15573021
includes only last 3 exons
WBGene00001130
'''

read_arg = args.coordinates.split(':')
selected_chrom = read_arg[0]
gen_coors = [int(read_arg[1]), int(read_arg[2])]

open_type = gzip.open if args.genome.endswith('.gz') else open
