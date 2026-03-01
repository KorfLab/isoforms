import argparse
import gzip
import glob
import os
import shutil

parser = argparse.ArgumentParser(description='find genes of interest '
			'from haman genes')
parser.add_argument('WBGenes', help='text with list of WBGene IDs and '
			'region of interest i.e. '
			'mod1,WBGene00003386,V,8910090,8910840,-')
parser.add_argument('gene_dir', help='directory with gff and fa files')
parser.add_argument('new_dir', help='name of directory to store '
			'fa and gff files')

args = parser.parse_args()

gfiles = {}
with open(args.WBGenes, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		line = line.split(':')
		gfiles[line[1]] = [line[0]]

for gff_file in glob.glob(args.gene_dir + '*.gff3'):
	with open(gff_file, 'rt') as fp:
		for line in fp:
			line = line.rstrip().split('\t')
			if line[1] == 'WormBase' and line[2] == 'gene':
				wbg = line[8].split(':')[1]
				if wbg in gfiles:
					gfiles[wbg].append(gff_file)
					base = gff_file.split('/')[-1].split('.')[:-1]
					gfiles[wbg].append(f'{args.gene_dir}'
										f'{'.'.join(base)}.fa')

if not os.path.isdir(f'{args.new_dir}'):
	os.mkdir(f'{args.new_dir}/')

for item in gfiles.items():
	if len(item[1]) == 1: continue
	src1 = item[1][1]
	src2 = item[1][2]
	dest = args.new_dir
	shutil.copy(src1, dest)
	shutil.copy(src2, dest)




# old code to search gff
# now using haman gene builds
'''
parser = argparse.ArgumentParser(
	description='get overlapping introns in genomic region of interest')
parser.add_argument('annotation', help='gff3 genome annotation')
parser.add_argument('coordinates', help='chromosome:loc i.e. X:99:999; '
	'chromosomes: I, II, III, IV, V, X, MtDNA')
#parser.add_argument('--fname', help='name file')

args = parser.parse_args()

read_arg = args.coordinates.split(':')
chrom = read_arg[0]
gen_coors = [int(read_arg[1]), int(read_arg[2])]

open_type = gzip.open if args.annotation.endswith('.gz') else open

count = 0
with open_type(args.annotation, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		line = line.split('\t')
		if len(line) == 9:
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
			
			
			
			
			
			
			
			
			
			
