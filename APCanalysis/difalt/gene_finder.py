import argparse
import gzip
import os
import glob

parser = argparse.ArgumentParser(
	description='find gene sequences in the C. elegans genome')
parser.add_argument('genome', help='genomic fasta file')
parser.add_argument('gffs', help='directory of gffs with genes of '
	'interest, output from search_gff.py')
	
'''
parser.add_argument('WBGenes', help='csv with list of WBGene IDs and '
	'region of interest i.e. '
	'WBGene, gene name, chr, left bound, right bound, '
	'first exon start, sense | '
	'WBGene00003386,mod1,V,8910090,8910840,8913992,-')
parser.add_argument('out_dir', help='name of directory to store '
			'fa and gff files')
'''

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
					start = line[3]
				if n_cds2 == n_cds:
					end = line[4]

	g_info = f'>{gene_name} {chrom}:{start}-{end} {sense} {wbg_gene}'
	genes[g_info] = []

for gene in genes.items():
	print(gene)
	
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
			n_counts = {}
			coors = gene[0].split(' ')[1].split(':')[1].split('-')
			coors = [int(coors[0]), int(coors[1])]
			for n in line:
				if gene[0] not in n_counts:
					n_counts[gene[0]] = 0
				n_counts[gene[0]] += 1
				print(n_counts)
				if coors[0] <= n_counts[gene[0]] <= end:
					gen_seqs[item[0]].append(n)

print(gen_seqs)

'''
		# match line index to region of interest 
		for info in gene_info.items():
			if info[1][1] != current_chrom: continue
			# save sequence identifiers and info
			sq_ids = [info[0]]
			info_list = info[1]
			sq_ids.extend(info_list)
			sq_ids = tuple(sq_ids)
			if sq_ids not in gen_seqs:
				gen_seqs[sq_ids] = []
			beg = int(info[1][2])
			end = int(info[1][3])
			for n in line:
				# add counts based on gene, not chrom
				if info[1][0] not in n_counts:
					n_counts[info[1][0]] = 0
				n_counts[info[1][0]] += 1
				if n_counts[info[1][0]] >= beg and n_counts[info[1][0]] <= end:
					gen_seqs[sq_ids].append(n)

# flip - strands
for info_seq in gen_seqs.items():
	if info_seq[0][6] == '-':
		comp_seq = []
		comps = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
		for n in info_seq[1]:
			comp_seq.append(comps[n])
		rev_seq = []
		for i in range(len(comp_seq)):
			i = i+1
			rev_seq.append(comp_seq[-i])
		gen_seqs[info_seq[0]] = rev_seq

# organize sequences into 80 nt lines
gen_seqs_80 = {}
for gen_seq in gen_seqs.items():
	#####################################print(len(gen_seq[1]))
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
	
for items in gen_seqs_80.items():
	seq_desc = (
		f'>{items[0][1]} {items[0][2]}:{items[0][3]}-{items[0][4]} '
		f'{items[0][7]} Gene:{items[0][0]} Gene_Start={items[0][5]} '
		f'Gene_End={items[0][6]}'
	)
	with open(f'{out}{items[0][1]}.fa', 'wt') as fp:
		fp.write(f'{seq_desc}\n')
		for seq in items[1]:
			fp.write(f'{seq}\n')
'''
