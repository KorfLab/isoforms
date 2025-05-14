import argparse
import glob
import math
import os
import random
import re
import sys

parser = argparse.ArgumentParser(description='more testing')
parser.add_argument('apc', help='path to apc directory')
parser.add_argument('build', help='path to build directory')
arg = parser.parse_args()

for ff in glob.glob(f'{arg.apc}/*.fa'):
	m = re.search('ch\.\d+', ff)
	ch = m.group()
	
	## RNA
	rna = {}
	total_reads = 0
	with open(f'{arg.apc}/{ch}.gff3') as fp:
		for line in fp:
			f = line.split()
			if f[1] != 'RNASeq_splice': continue
			beg = int(f[3])
			end = int(f[4])
			score = float(f[5])
			rna[(beg,end)] = int(score)
			total_reads += int(score)
	
	## APC
	apc = {}
	total_apc = 0
	with open(f'{arg.build}/isoformer/{ch}.gff') as fp:
		for line in fp:
			f = line.split()
			if len(f) < 8 or f[2] != 'intron': continue
			beg = int(f[3])
			end = int(f[4])
			score = float(f[5])
			sig = (beg, end)
			if sig not in apc: apc[sig] = 0
			apc[sig] += score
			total_apc += score
	for sig in apc: apc[sig] /= total_apc
	
	## Sample
	probs = list(apc.values())
	sigs = list(apc.keys())
	sample = {}
	for sig in random.choices(sigs, probs, k=total_reads):
		if sig not in sample: sample[sig] = 0
		sample[sig] += 1
	
	## Reporting
	all_sigs = list(rna.keys())
	for sig in apc:
		if sig not in all_sigs: all_sigs.append(sig)
	print(ch, total_reads)
	for sig in all_sigs:
		if sig not in sample: sample[sig] = 0
		if sig not in rna: rna[sig] = 0
		print(sig, rna[sig], sample[sig], sep='\t')
	print()
		