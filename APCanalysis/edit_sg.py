# edit smallgenes gff for use with geniso_dif
import sys

gff = sys.argv[1]

with open(gff, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		line = line.split('\t')
		if line[2] in ['intron', 'CDS', 'exon']:
			print('\t'.join(line))
