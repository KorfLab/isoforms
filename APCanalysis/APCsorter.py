import argparse
import glob
import csv

parser = argparse.ArgumentParser(description=
	'code to interpret the APC results')
parser.add_argument('distances', type=str, metavar='<file>', 
	help='text file with distances')
parser.add_argument('smallgenes', type=str, metavar='<directory>',
	help='directory with smallgenes fasta and gff files')
	
args = parser.parse_args()

with open(args.distances, 'r') as csvfile:
	
