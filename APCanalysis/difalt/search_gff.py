import argparse
import gzip
import glob
import os
import shutil

parser = argparse.ArgumentParser(description='create gffs for WBGenes '
	'of interest')
parser.add_argument('WBGenes', help='csv with list of WBGene IDs and '
	'region of interest i.e. '
	'WBGene, gene name, chr, left bound, right bound, '
	'first exon start, sense | '
	'WBGene00003386,mod1,V,8910090,8910840,8913992,-')
parser.add_argument('annotation', help='annotation gff')
parser.add_argument('out_dir', help='name of directory to store '
			'fa and gff files')

args = parser.parse_args()

ginfo = {}
with open(args.WBGenes, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		if line.startswith('#'): continue
		line = line.split(',')
		ginfo[line[0]] = line[1:]
		
print(ginfo)

open_type = gzip.open if args.annotation.endswith('.gz') else open

count = 0
with open_type(args.annotation, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		line = line.split('\t')
		#if len(line) == 9:
			#print(line)
			'''
			if (line[0] == chrom and line[1] == 'WormBase' and
				line[2] == 'CDS'):
				# get any CDS regions that overlap
				if ((int(line[3]) <= gen_coors[1] and 
					int(line[3]) >= gen_coors[0]) or 
					(int(line[4]) <= gen_coors[1] and
					int(line[4]) >= gen_coors[0])):
						# negative CDS coors needed for make_svg.py
						line[3] = str(int(line[3]) - gen_coors[0] +1)
						line[4] = str(int(line[4]) - gen_coors[0] +1)
						print('\t'.join(line))
			if line[0] == chrom and line[1] == 'RNASeq_splice':
				# only get introns within coors
				if (int(line[3]) <= gen_coors[1] and 
					int(line[3]) >= gen_coors[0] and
					int(line[4]) <= gen_coors[1] and
					int(line[4]) >= gen_coors[0]):
						line[3] = str(int(line[3]) - gen_coors[0] +1)
						line[4] = str(int(line[4]) - gen_coors[0] +1)
						print('\t'.join(line))
			'''

