import argparse
import gzip

parser = argparse.ArgumentParser(
	description='get overlapping introns in genomic region of interest')
parser.add_argument('annotation', help='gff3 genome annotation')
parser.add_argument('coordinates', help='chromosome:loc i.e. X:99:999; '
	'chromosomes: I, II, III, IV, V, X, MtDNA')
#parser.add_argument('--fname', help='name file')

args = parser.parse_args()

'''
dyn-1 whole gene sequence is X:15568921:15573021
for APC, only use 15571571 to 15573021
includes only last 3 exons
WBGene00001130
'''

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
			
			
			
			
			
			
			
			
			
			
			
