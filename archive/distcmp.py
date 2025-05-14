import argparse
import glob
import math
import os
import random
import re
import sys

from grimoire.genome import Reader

def harvest_introns(fasta, gff, rna=False):
	# get all introns
	introns = []
	genome = Reader(gff=gff, fasta=fasta)
	chrom = next(genome)
	for f in chrom.ftable.features:
		if f.type != 'intron': continue
		if rna:
			if f.source != 'RNASeq_splice': continue
			if f.strand != '+': continue
		introns.append(f)

	# aggregate introns by position
	ipos = {}
	total = 0
	for intron in introns:
		sig = (intron.beg, intron.end)
		if sig not in ipos: ipos[sig] = 0
		ipos[sig] += intron.score
		total += intron.score
	for sig in ipos: ipos[sig] /= total # convert to probs

	return ipos

def freq_dist(vals, pseudo=0):
	assert(pseudo >= 0)
	total = sum(vals) + pseudo * len(vals)
	return [(x+pseudo) / total for x in vals]


def d1(P, Q): # manhattan
	return sum([abs(p - q) for p, q in zip(P, Q)])

def d2(P, Q): # cartesian
	return math.sqrt(sum([(p-q)**2 for p, q in zip(P, Q)]))

def dB(P, Q): # Bhattacharyya
	return -math.log(sum([math.sqrt(p * q) for p, q in zip(P, Q)]))

def dKL(P, Q): # symmetrical version
	return sum([p * math.log2(p/q) + q * math.log2(q/p) for p, q in zip(P, Q)])


parser = argparse.ArgumentParser(description='some testing thing')
parser.add_argument('apc', help='path to apc directory')
parser.add_argument('build', help='path to build directory')
arg = parser.parse_args()

for ff in glob.glob(f'{arg.apc}/*.fa'):
	m = re.search('ch\.\d+', ff)
	ch = m.group()
	obs = f'{arg.apc}/{ch}.gff3'
	exp = f'{arg.build}/isoformer/{ch}.gff'

	obsi = harvest_introns(ff, obs, rna=True)
	expi = harvest_introns(ff, exp)
	for sig in obsi:
		if sig not in expi: expi[sig] = 0
	for sig in expi:
		if sig not in obsi: obsi[sig] = 0
	P = []
	Q = []
	for sig in obsi:
		P.append(obsi[sig])
		Q.append(expi[sig])
	P = freq_dist(P, pseudo=1e-6)
	Q = freq_dist(Q, pseudo=1e-6)
	print(ch, d1(P,Q), d2(P,Q), dB(P,Q), dKL(P,Q), sep='\t')

