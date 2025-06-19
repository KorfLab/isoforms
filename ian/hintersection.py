import argparse
import glob
import re
import sys

import isoform

def intersect(rna, apc):
	shared = 0
	for sig in rna:
		if sig in apc: shared += min(rna[sig], apc[sig])
	return 1 - shared


parser = argparse.ArgumentParser()
parser.add_argument('smallgenes', type=str, metavar='<dir>',
	help='smallgenes directory')
parser.add_argument('apc', type=str, metavar='<dir>',
	help='apc isos directory')
arg = parser.parse_args()

cutoffs = [0] * 101
total = 0
for gff3 in glob.glob(f'{arg.smallgenes}/*.gff3'):
	m = re.search(r'(ce\.\d+\.\d+)\.gff3', gff3)
	locus = m.group(1)
	apcfile = glob.glob(f'{arg.apc}/{locus}.*')[0]
	rna_introns = isoform.get_introns(gff3)
	apc_introns = isoform.get_introns(apcfile)	
	d = int(intersect(rna_introns, apc_introns) * 100)
	cutoffs[d] += 1
	total += 1

running = 0
for i in range(len(cutoffs)):
	same = 100 - i
	running += cutoffs[i]
	print(running/total)

