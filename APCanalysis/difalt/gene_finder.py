import argparse
import gzip
import os
import glob

parser = argparse.ArgumentParser(
	description='find gene sequences in the C. elegans genome')
parser.add_argument('genome', help='genomic fasta file')
parser.add_argument('gffs', help='directory of gffs with genes of '
	'interest, output from search_gff.py')
parser.add_argument('out_dir', help='name of directory to store '
	'fa and gff files')
parser.add_argument('--flank', required=False, type=int, default=99, 
	metavar='<int>', help='add n nucleotides of genomic '
	'flank to sequence')

args = parser.parse_args()

if not args.gffs.endswith('/'):
	args.gffs = f'{args.gffs}/'

# collect info for sequence from gffs
genes = {}	
for file in glob.glob(f'{args.gffs}*'):
	wbg_gene = None
	gene_name = file.split('/')[-1].split('.')[0]
	chrom = None
	sense = None
	start = None
	end = None
	n_cds = 0
	gff_lines = []
	with open(file, 'rt') as fp:
		for line in fp:
			line = line.rstrip().split('\t')
			gff_lines.append(line)
		
		n_cds = 0
		for line in gff_lines:
			if line[2] == 'CDS': n_cds += 1
		
		n_cds2 = 0
		for line in gff_lines:
			if line[2] == 'gene':
				wbg_gene = line[8].split(';')[0].split(':')[1]
				chrom = line[0]
				sense = line[6]
			if line[2] == 'CDS':
				n_cds2 += 1
				if n_cds2 == 1:
					start = int(line[3]) - args.flank
				if n_cds2 == n_cds:
					end = int(line[4]) + args.flank

	g_info = f'>{gene_name} {chrom}:{start}-{end} {sense} {wbg_gene}'
	genes[g_info] = []
	
open_type = gzip.open if args.genome.endswith('.gz') else open

# gather seq strings for each entry
with open_type(args.genome, 'rt') as fp:
	current_chrom = None
	gen_seqs = {}
	n_counts = {}
	# go line by line, don't save whole chr seq
	for line in fp:
		line = line.rstrip()
		if line.startswith('>'):
			current_chrom = line.split(' ')[0][1:]
			# reset counts between different genes, not just chroms
			n_counts = {}
			continue
			
		# match line index to region of interest
		for gene in genes.items():
			chrom = gene[0].split(' ')[1].split(':')[0]
			if chrom != current_chrom: continue
			coors = gene[0].split(' ')[1].split(':')[1].split('-')
			coors = [int(coors[0]), int(coors[1])]
			for n in line:
				if gene[0] not in n_counts:
					n_counts[gene[0]] = 0
				n_counts[gene[0]] += 1
				if coors[0] <= n_counts[gene[0]] <= coors[1]:
					if gene[0] not in gen_seqs:
						gen_seqs[gene[0]] = []
					gen_seqs[gene[0]].append(n)

# not sure this is necessary
# i think the fasta should not be translated
# forgot why i did this
'''
# flip - strands
for info_seq in gen_seqs.items():
	if info_seq[0].split(' ')[2] == '-':
		comp_seq = []
		comps = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
		for n in info_seq[1]:
			comp_seq.append(comps[n])
		rev_seq = []
		for i in range(len(comp_seq)):
			i = i+1
			rev_seq.append(comp_seq[-i])
		gen_seqs[info_seq[0]] = rev_seq
'''

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

if args.out_dir.endswith('/'): 
	out = args.out_dir
else:
	out = f'{args.out_dir}/'
	
if not os.path.exists(out):
	os.makedirs(out)
	
for item in gen_seqs_80.items():
	info = item[0].split(' ')
	seq_desc = f'{info[0]} {info[1]} {info[2]} {info[3]}'
	fname = info[0].split('>')[1]
	with open(f'{out}{fname}.fa', 'wt') as fp:
		fp.write(f'{seq_desc}\n')
		for seq in item[1]:
			fp.write(f'{seq}\n')
			
# need to adjust gff coors to 0
for item in genes.items():
	coors = item[0].split(' ')[1].split(':')[1].split('-')
	start = int(coors[0])
	end = int(coors[1])
	
for item in gen_seqs_80.items():
	print(item[0])
