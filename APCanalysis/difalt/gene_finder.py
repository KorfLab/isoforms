import argparse
import gzip

parser = argparse.ArgumentParser(
	description='find gene sequences in the C. elegans genome')
parser.add_argument('genome', help='')
parser.add_argument('coordinates', help='chromosome:loc i.e. X:99:999; '
	'chromosomes: I, II, III, IV, V, X, MtDNA')
parser.add_argument('--seq_desc', required=True, help='add description for sequence '
	'in first line of fasta file')
parser.add_argument('--fname', required=True, help='add file name')
parser.add_argument('--rev', required=False, action='store_true', 
	help='is the sequence on the reverse strand?')

args = parser.parse_args()


'''
dyn-1 whole gene sequence is X:15568921:15573021
for APC, only use 15571571 to 15573021
includes only last 3 exons
may be too long
1kb: 15572021 15573021
1066 bp: 15571955 15573021
1kb is too long to run
don't forget to set flank to 0
test 900 bp: 15572121 15573021 *too long
test 800 bp: 15572221 15573021 *too long
test 700 bp: 15572321 15573021 *need longer first exon
test 750 bp: 15572271 15573021
test 500 bp: 15572521 15573021

unc-16 to test exon 16
unc-16 900 bp: III:9553541:9554441

mod-1 test exon 0
1082 bp V:8909913:8910995


'''

read_arg = args.coordinates.split(':')
selected_chrom = read_arg[0]
gen_coors = [int(read_arg[1]), int(read_arg[2])]

open_type = gzip.open if args.genome.endswith('.gz') else open

with open_type(args.genome, 'rt') as fp:
	current_chrom = None
	gen_seq = []
	n_counts = 0
	for line in fp:
		line = line.rstrip()
				
		# check if loop is currently on the chromosome of interest
		if line.startswith('>'):
			current_chrom = line.split('>')[1]
			continue
		if current_chrom != selected_chrom: continue
		#print(current_chrom, line)
	
		# build seq string
		for n in line:
			n_counts += 1
			if n_counts >= gen_coors[0] and n_counts <= gen_coors[1]:
				#print(n_counts, n)
				gen_seq.append(n)

if args.rev:
	rev_seq = []
	comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	for n in gen_seq:
		rev_seq.append(comp[n])
	gen_seq = rev_seq

# organize sequences into 80 nt lines
seq_lines = []
seq_line = []
for n in gen_seq:
	if len(seq_line) < 80:
		seq_line.append(n)
	else:
		seq_lines.append(''.join(seq_line))
		seq_line = []
		seq_line.append(n)

if len(seq_line) != 0:
	seq_lines.append(''.join(seq_line))
	
with open(f'{args.fname}.fa', 'wt') as fp:
	fp.write(f'>{args.seq_desc}\n')
	for seq in seq_lines:
		fp.write(seq)


