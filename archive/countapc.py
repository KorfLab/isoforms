#!/usr/bin/env python3

import argparse
import random
from statistics import mean, stdev
import sys

import isoform
from isoform import Locus

parser = argparse.ArgumentParser(description='Count isoforms in fasta files')
parser.add_argument('files', nargs='+', help='fasta files')
parser.add_argument('--limit', required=False, type=int, default=100,
	metavar='<int>', help='limit number of transcripts [%(default)i]')
parser.add_argument('--min_intron', required=False, type=int, default=35,
	metavar='<int>', help='minimum length of intron [%(default)i]')
parser.add_argument('--min_exon', required=False, type=int, default=25,
	metavar='<int>', help='minimum length exon [%(default)i]')
parser.add_argument('--flank', required=False, type=int, default=99,
	metavar='<int>', help='genomic flank on each side [%(default)i]')
arg = parser.parse_args()

models = (None, None, None, None, None, None)
weights = (None, None, None, None, None, None)
icost = 0

for path in arg.files:
	name, seq = next(isoform.read_fasta(path))
	locus = Locus(name, seq, arg.min_intron, arg.min_exon, arg.flank,
		models, weights, icost, limit=arg.limit, countonly=True)
	print(locus.name, len(seq), locus.isocount, sep='\t', flush=True)
