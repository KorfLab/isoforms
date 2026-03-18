import argparse
import gzip

parser = argparse.ArgumentParser(
	description='find gene sequences in the C. elegans genome')
parser.add_argument('genome', help='genomic fasta file')
parser.add_argument('WBGenes', help='csv with list of WBGene IDs and '
	'region of interest i.e. '
	'WBGene, gene name, chr, left bound, right bound, '
	'first exon start, sense | '
	'WBGene00003386,mod1,V,8910090,8910840,8913992,-')

args = parser.parse_args()

gene_info = {}
with open(args.WBGenes, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		if line.startswith('#'): continue
		line = line.split(',')
		gene_info[line[0]] = line[1:]

open_type = gzip.open if args.genome.endswith('.gz') else open

# gather seq strings for each entry
with open_type(args.genome, 'rt') as fp:
	current_chrom = None
	gen_seqs = {}
	n_counts = 0
	# go line by line, don't save whole chr seq
	for line in fp:
		line = line.rstrip()
		if line.startswith('>'):
			current_chrom = line.split(' ')[0][1:]
			n_counts = 0
			continue
		# match line index to region of interest 
		for info in gene_info.items():
			# gene name and strandedness
			gns = f'{info[1][1]}{info[1][5]}'
			if gns not in gen_seqs:
				gen_seqs[gns] = []
			if info[1][1] != current_chrom: continue
			beg = int(info[1][2])
			end = int(info[1][3])
			for n in line:
				n_counts += 1
				if n_counts >= beg and n_counts <= end:
					gen_seqs[gns].append(n)
				
# flip - strands
for gen_seq in gen_seqs.items():
	if gen_seq[0].endswith('-'):
		comp_seq = [] 
		comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
		for n in gen_seq[1]:
			comp_seq.append(comp[n])
		rev_seq = []
		for i in range(len(comp_seq)):
			i = i+1 
			rev_seq.append(comp_seq[-i])
		gen_seqs[gen_seq[0]] = rev_seq

# organize sequences into 80 nt lines
gen_seqs_80 = {}
for gen_seq in gen_seqs.items():
	seq_line = []
	seq_lines = []
	for n in gen_seq[1]:
		if len(seq_line) < 80:
			seq_line.append(n)
		else:
			seq_lines.append(''.join(seq_line))
			seq_line = []
			seq_line.append(n)
	if len(seq_line) != 0:
		seq_lines.append(''.join(seq_line))

	gen_seqs_80[gen_seq[0]] = seq_lines
	
print(gen_seqs_80)

'''
	
with open(f'{args.fname}.fa', 'wt') as fp:
	fp.write(f'>{args.seq_desc}\n')
	for seq in seq_lines:
		fp.write(f'{seq}\n')
'''

