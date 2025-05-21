#!/usr/bin/env python3

import argparse
import glob
import math
import random
import os
import statistics
import sys

from grimoire.genome import Reader
import isoform

#######
# CLI #
#######

parser = argparse.ArgumentParser(description='model builder for apc set')
parser.add_argument('apc', help='path to apc directory')
arg = parser.parse_args()

for ff in glob.glob(f'{arg.apc}/*.fa'):
	print(ff)
	icount = {}
	gf = ff[:-2] + 'gff3'
	genome = Reader(gff=gf, fasta=ff)
	chrom = next(genome)
	introns = []
	scores = []
	for f in chrom.ftable.features:
		if f.source == 'RNASeq_splice':
			introns.append( (f.beg, f.end) )
			scores.append(f.score)
	
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
	