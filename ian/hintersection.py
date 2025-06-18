import argparse
import glob
import re
import sys

import isoform

# need to compare different distance methods
# intersections, manhattan, pruning vs. keeping



parser = argparse.ArgumentParser()
parser.add_argument('smallgenes', type=str, metavar='<dir>',
	help='smallgenes directory')
parser.add_argument('apc', type=str, metavar='<dir>',
	help='apc isos directory')
arg = parser.parse_args()

for gff3 in glob.glob(f'{arg.smallgenes}/*.gff3'):
	m = re.search(r'(ce\.\d+\.\d+)\.gff3', gff3)
	locus = m.group(1)
	apcfile = glob.glob(f'{arg.apc}/{locus}.*')[0]
	rna_introns = isoform.get_introns(gff3)
	apc_introns = isoform.get_introns(apcfile)
	
