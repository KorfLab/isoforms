import argparse
import glob

parser = argparse.ArgumentParser(description='get gene symbols from '
	'smallgenes WBGene IDs')
parser.add_argument('smallgenes')
parser.add_argument('annotation')

args = parser.parse_args()

wbgenes = []
for gff in glob.glob(f'{args.smallgenes}*gff3'):
	with open(gff, 'rt') as fp:
		for line in fp:
			line = line.rstrip()
			line = line.split('\t')
			if line[2] == 'gene':
				wbgene = line[8].split(':')[1]
				wbgenes.append(wbgene)
				break

symbols = {}
with open(args.annotation, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		line = line.split()
		if line[2] == 'gene':
			info = {x.split('=')[0]: x.split('=')[1] for x in line[8].split(';')}
			wbg = info['Name']
			sym = info['sequence_name']
			symbols[wbg] = sym
			
print(symbols)

