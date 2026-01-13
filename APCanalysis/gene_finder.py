import argparse
import gzip

parser = argparse.ArgumentParser(
	description='find gene sequences in the C. elegans genome')
parser.add_argument('genome', help='')
parser.add_argument('coordinates', help='chromosome:loc i.e. X:99:999; ' 
	'I, II, III, IV, V, X, MtDNA')

args = parser.parse_args()

# dyn-1 X:15568921:15573031

# put chromosome sequences into dictionary
# storing in dictionary is slow
# thow away unneeded text instead
read_arg = args.coordinates.split(':')
chrom = read_arg[0]
coors = [read_arg[1], read_arg[2]]

print(chrom, coors)

open_type = gzip.open if args.genome.endswith('.gz') else open
with open_type(args.genome, 'rt') as fp:
	count = 0
	for line in fp:
		line = line.rstrip()
		if line.startswith('>'):
			if chrom == line.split('>')[1]:
				print(line)
				for n in line:
					count += 1
				
#print(chr_seqs)
#coors = args.coordinates.split(':')

# use 1 based indexing for coordinates

#print(chr_seqs[coors[0]][int(coors[1])-1:int(coors[2])])
