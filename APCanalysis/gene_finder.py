import argparse
import gzip

parser = argparse.ArgumentParser(
	description='find gene sequences in the C. elegans genome')
parser.add_argument('genome', help='')
parser.add_argument('coordinates', help='chromosome:loc i.e. X:99:999; '
	'chromosomes: I, II, III, IV, V, X, MtDNA')
parser.add_argument('--seq_desc', help='add description for sequence '
	'in first line of fasta and file name')

args = parser.parse_args()

'''
dyn-1 whole gene sequence is X:15568921:15573021
for APC, only use 15571571 to 15573021
includes only last 3 exons
may be too long
1kb: 15572021 15573021
1066 bp: 15571955 15573021
test 500 bp: 15572521 15573021

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

# write output file
if args.seq_desc:
	seq_desc = args.seq_desc
else:
	seq_desc = '_'.join(read_arg)	
	
with open(f'{seq_desc}.fa', 'wt') as fp:
	fp.write(f'>{'_'.join(read_arg)}\n')
	for seq in seq_lines:
		fp.write(seq)


