import argparse
import glob
import math
import random
import os
import statistics
import sys

from grimoire.genome import Reader

#######
# CLI #
#######

parser = argparse.ArgumentParser(description='introns vs. read depth')
parser.add_argument('dir', help='path to smallgenes directory')
arg = parser.parse_args()

print('introns\tcount\n')
d = {}
genes = 0
for ff in glob.glob(f'{arg.dir}/*.fa'):
	gf = ff[:-2] + 'gff3'
	genome = Reader(gff=gf, fasta=ff)
	chrom = next(genome)
	genes += 1

	introns = {}
	for f in chrom.ftable.features:
		if f.source == 'RNASeq_splice':
			introns[(f.beg, f.end)] = f.score
	
	print(ff, end=': ', file=sys.stderr)
	weights = list(introns.values())
	sigs = list(introns.keys())
	for reads in range(1000, 10001, 1000):
		obs = []
		for i in range(100):
			seen = set()
			for sig in random.choices(sigs, weights=weights, k=reads):
				seen.add(sig)
			obs.append(len(seen))
		print(statistics.mean(obs), end=' ', file=sys.stderr)
		if reads not in d: d[reads] = 0
		d[reads] += statistics.mean(obs)
	print(file=sys.stderr, flush=True)
	
for n, x in d.items():
	print(n, x/genes, sep='\t')

			
	

"""
	
	
	for exp in range(5):
		seen = set()
		for i in range(10, 1000, 10):
			for intron in random.choices(introns, weights=scores, k=i):
				seen.add(intron)
			n = len(seen) # number of introns seen so far
			if i not in icount: icount[i] = []
			icount[i].append(n)
	
	for samples, data in icount.items():
		print(samples, statistics.mean(data), flush=True)
"""