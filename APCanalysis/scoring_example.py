# get data for figure showing isoform scoring
# using the simplest gene

import sys
import isoform2
from isoform2 import Locus

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

model = isoform2.read_splicemodel(arg.model)

weights = {
	'wacc': 1.0,
	'wdon': 1.0,
	'wexs': 1.0,
	'wins': 1.0,
	'wexl': 1.0,
	'winl': 1.0,
	'winf': 1.0,
} 

if arg.weights:
	with open(arg.weights, 'r') as fp:
		for line in fp.readlines():
			line = line.rstrip()
			line = line.split(',')
			weights['wacc'] = float(line[2])
			weights['wdon'] = float(line[3])
			weights['wexs'] = float(line[4])
			weights['wins'] = float(line[5])
			weights['wexl'] = float(line[6])
			weights['winl'] = float(line[7])
			weights['winf'] = float(line[8])
						
constraints = {
	'min_intron': 35,
	'min_exon': 25,
	'flank': 99
}

name, seq = next(isoform2.read_fasta(arg.fasta))

locus = Locus(name, seq, model, constraints=constraints, 
	weights=weights, limit=100)
















