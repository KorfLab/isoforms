# number of introns vs sequencing depth
# number of isoforms vs sequence length

# generate random sequences, 300 - 2000 base pairs

import random
import argparse

parser = argparse.ArgumentParser(description='sequencing simulator')
parser.add_argument('gff', help='gff file for example gene')
parser.add_argument('big_gff', help='genomic gff')

args = parser.parse_args()

GC = 0.36
slen = 300
depth = 10000
inc = 10

#for i in range(low, 
rseq = ''.join(random.choices(['A', 'C', 'G', 'T'], 
				weights = [(1-GC)/2, GC/2, GC/2, (1-GC)/2], k = slen))
				
#print(rseq)

# resample a gene with a bunch of introns, not necessarily a small gene
# number of different introns depends on how deeply you sequence
# RNA seq values of score, use as weights to random sampler
# give 1000 of these reads based on the weights
# better to do with real data

introns = []
scores = []
with open(args.gff, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		line = line.split('\t')
		if line[1] == 'RNASeq_splice':
			introns.append((line[3], line[4]))
			scores.append(float(line[5]))

total = 0
for s in scores:
	total += s
	
wt = [x/total for x in scores]

int_sub = [x for x in range(len(introns))]

res = {}
for i in range(inc, depth+inc, inc):
	int_samp = random.choices(int_sub, weights = wt, k = i)
	seen = []
	for intron in int_samp:
		if intron not in seen:
			seen.append(intron)
	res[i] = len(seen)
	
#for r in res:
#	print(r, res[r])
	
# need a gene with a ton of introns
# 1pc just goe by chromosome
# how do i find a gene with a lot of introns
with open(args.big_gff, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		line = line.split()
		if 'intron' in line:
			print(line)

			
		
	


