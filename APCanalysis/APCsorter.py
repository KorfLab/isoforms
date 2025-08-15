import argparse
import glob
import csv

parser = argparse.ArgumentParser(description=
	'code to interpret the APC results')
parser.add_argument('res', type=str, metavar='<file>', 
	help='text file with distances')
parser.add_argument('smallgenes', type=str, metavar='<directory>',
	help='directory with smallgenes fasta and gff files')
	
args = parser.parse_args()

genes = {}
count = 0
with open(args.res, 'rt') as cfp:
	reader = csv.reader(cfp, delimiter=',')
	for row in reader:
		if row[0].startswith('c'):
			genes[row[0]] = float(row[1])
			count += 1

genes = sorted(genes.items(), key=lambda item: item[1])

bins = {x/10: 0 for x in range(0, 10, 1)}
for gene in genes:
	for i, b in enumerate(bins):
		if gene[1] > b and gene[1] < b+0.1:
			bins[i/10] += 1
'''
print(bins)
print(count)
print(sum([item[1] for item in bins.items()]))
'''

for b in bins.items():
	print(f'{b[0]},{b[1]}')

'''
for file in glob.glob(f'{args.smallgenes}/*.gff3'):
	print(file)
'''	
