import argparse
import gzip
import os

parser = argparse.ArgumentParser(description='create gffs for WBGenes '
	'of interest')
parser.add_argument('WBGenes', help='csv with list of WBGene IDs and '
	'region of interest i.e. '
	'WBGene, gene name, chr, region start, region end, '
	'gene start, gene end, sense | '
	'WBGene00003386,mod1,V,8910090,8910840,8913992,-')
parser.add_argument('annotation', help='annotation gff')
parser.add_argument('out_dir', help='name of directory to store '
			'fa and gff files')

args = parser.parse_args()

gene_info = {}
with open(args.WBGenes, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		if line.startswith('#'): continue
		line = line.split(',')
		gene_info[line[0]] = line[1:]

open_type = gzip.open if args.annotation.endswith('.gz') else open

gene_lines = {}
with open_type(args.annotation, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		line = line.split('\t')
		for info in gene_info.items():
			gid_info = (info[0], info[1][0])
			if gid_info not in gene_lines:
				gene_lines[gid_info] = []
			if len(line) == 9:
				f_start, f_end = int(line[3]), int(line[4])
				g_start, g_end = int(info[1][2]), int(info[1][3])
				
				# match WBGene
				if line[2] == 'gene':
					wbg_gff = line[8].split(';')[0].split(':')[1]
					if info[0] == wbg_gff:
						gene_lines[gid_info].append(line)
				
				# match WBGene
				if line[2] == 'mRNA':
					wbg_gff = line[8].split(';')[1].split(':')[1]
					if info[0] == wbg_gff:
						# only add first entry of WBGene
						if len(gene_lines[gid_info]) == 0:
							gene_lines[gid_info].append(line)
						
				'''
				# get any CDS regions that overlap region of interest
				if (line[0] == info[1][1] and line[1] == 'WormBase' and
					line[2] == 'CDS'):
					if f_start <= g_end and f_end >= g_start:
						gene_lines[gid_info].append(line)
				'''
				
				# get any CDS regions that overlap entire gene
				if (line[0] == info[1][1] and line[1] == 'WormBase' and
					line[2] == 'CDS'):
					gene_lines[gid_info].append(line)
				
				# get any introns only within 
				if line[0] == info[1][1] and line[1] == 'RNASeq_splice': 
					if f_start >= g_start and f_end <= g_end:		
						gene_lines[gid_info].append(line)
	
if args.out_dir.endswith('/'): 
	out = args.out_dir
else:
	out = f'{args.out_dir}/'
	
if not os.path.exists(out):
	os.makedirs(out)
						
for item in gene_lines.items():
	with open(f'{out}{item[0][1]}.gff3', 'wt') as fp:
		for line in item[1]:
			fp.write('\t'.join(line)+'\n')
				
			
