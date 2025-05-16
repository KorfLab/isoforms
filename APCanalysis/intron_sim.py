# number of introns vs sequencing depth
# number of isoforms vs sequence length

# generate random sequences, 300 - 2000 base pairs

import random
import argparse

parser = argparse.ArgumentParser(description='sequencing simulator')
parser.add_argument('gff', help='gff file for example gene')

args = parser.parse_args()

GC = 0.36
slen = 300

#for i in range(low, 
rseq = ''.join(random.choices(['A', 'C', 'G', 'T'], 
				weights = [(1-GC)/2, GC/2, GC/2, (1-GC)/2], k = slen))
				
#print(rseq)

# resample a gene with a bunch of introns, not necessarily a small gene
# number of different introns depends on how deeply you sequence
# RNA seq values of score, use as weights to random sampler
# give 1000 of these reads based on the weights
# better to do with real data

with open(args.gff, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		print(line)
