#!/usr/bin/env python3

import argparse
import glob
import math
import os
import sys

from grimoire.genome import Reader

parser = argparse.ArgumentParser(description='model builder for apc set')
parser.add_argument('apc', help='path to apc directory')
arg = parser.parse_args()

for ff in glob.glob(f'{arg.apc}/*.fa'):
	gf = ff[:-2] + 'gff3'
	genome = Reader(gff=gf, fasta=ff)
	chrom = next(genome)
	genes = chrom.ftable.build_genes()
	if len(genes) != 1: sys.exit('multiple genes')
	gene = genes[0]
	
	introns = {}
	for tx in gene.transcripts():
		for intron in tx.introns:
			sig = (intron.beg, intron.end)
			if sig not in introns: introns[sig] = True

	maxexp = 0
	for f in chrom.ftable.features:
		sig = (f.beg, f.end)
		if f.source == 'RNASeq_splice' and sig in introns and f.score > maxexp:
			maxexp = int(f.score)
	print(maxexp)