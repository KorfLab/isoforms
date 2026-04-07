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
			if line[0] != info[1][1]: continue
			if gid_info not in gene_lines:
				gene_lines[gid_info] = []
			if len(line) == 9:
				gff_start, gff_end = int(line[3]), int(line[4])
				subreg_start, subreg_end = int(info[1][2]), int(info[1][3])
				gen_start, gen_end = int(info[1][4]), int(info[1][5])
				
				# match WBGene
				if line[2] == 'gene':
					wbg_gff = line[8].split(';')[0].split(':')[1]
					if info[0] == wbg_gff:
						gene_lines[gid_info].append(line)
				
				# match WBGene
				if line[2] == 'mRNA':
					wbg_gff = line[8].split(';')[1].split(':')[1]
					if info[0] == wbg_gff:
						# only add first mRNA entry
						# not all mRNA ID transcripts start with 'a'
						if 'mRNA' not in [x[2] for x in gene_lines[gid_info]]:
							gene_lines[gid_info].append(line)
				
				# get any CDS regions that overlaps entire gene
				if (line[0] == info[1][1] and line[1] == 'WormBase' and
					line[2] == 'CDS'): 
					if gff_start <= gen_end and gff_end >= gen_start:
						gene_lines[gid_info].append(line)
						
				# get any CDS regions that overlap region of interest
				if (line[0] == info[1][1] and line[1] == 'WormBase' and
					line[2] == 'CDS'): 
					if gff_start <= subreg_end and gff_end >= subreg_start:
						gene_lines[gid_info].append(line)
				
				# get any introns only within region
				if line[0] == info[1][1] and line[1] == 'RNASeq_splice': 
					#print(subreg_start, subreg_end)
					if gff_start >= subreg_start and gff_end <= subreg_end:
						# remove low count introns for testing/readability
						if int(line[5]) < 20000: continue
						#print(line)
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
				
			
