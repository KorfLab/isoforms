#!/usr/bin/env python3

import argparse
import random
from statistics import mean, stdev
import sys

import isoform
from isoform import Locus

parser = argparse.ArgumentParser(description='Random sequence experiemnt')
parser.add_argument('samples', type=int, help='number of runs at each length')
parser.add_argument('--init', required=False, type=int, default=300,
	metavar='<int>', help='starting size [%(default)i]')
parser.add_argument('--term', required=False, type=int, default=600,
	metavar='<int>', help='ending size [%(default)i]')
parser.add_argument('--step', required=False, type=int, default=100,
	metavar='<int>', help='step size [%(default)i]')
parser.add_argument('--min_intron', required=False, type=int, default=35,
	metavar='<int>', help='minimum length of intron [%(default)i]')
parser.add_argument('--min_exon', required=False, type=int, default=25,
	metavar='<int>', help='minimum length exon [%(default)i]')
parser.add_argument('--flank', required=False, type=int, default=99,
	metavar='<int>', help='genomic flank on each side [%(default)i]')
arg = parser.parse_args()

flank = 'A' * arg.flank
models = (None, None, None, None, None, None)
weights = (None, None, None, None, None, None)
icost = 0
for size in range(arg.init, arg.term+1, arg.step):
	count = []
	for _ in range(arg.samples):
		seq = flank + ''.join(random.choices('ACGT', k=size)) + flank
		locus = Locus('x', seq, arg.min_intron, arg.min_exon, arg.flank,
			models, weights, icost, countonly=True)
		count.append(locus.isocount)
	print(f'{size}\t{mean(count):.1f}\t{stdev(count):.1f}\t{min(count)}\t{max(count)}')
