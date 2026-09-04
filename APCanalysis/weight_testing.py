import argparse

parser = argparse.ArgumentParser(description='get genomic intron and '
	'exon length stats in csv format')
parser.add_argument('genome')
parser.add_argument('annotation')

args = parser.parse_args()
