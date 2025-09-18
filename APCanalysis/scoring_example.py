# get data for figure showing isoform scoring
# using the simplest gene

import sys
import isoform2.py

fasta = sys.argv[1]
gff = sys.argv[2]

import argparse

parser = argparse.ArgumentParser(description='get data for figure '
	'showing how isoforms are scored')
parser.add_argument('fasta')
parser.add_argument('gff')
parser.add_argument('model')
parser.add_argument('--weights', required=False)

arg = parser.parse_args()


