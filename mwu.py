import argparse
import glob
from scipy.stats import mannwhitneyu

from grimoire.genome import Reader
import isoform2

def list2prob(vec):
	s = sum(vec)
	p = []
	for v in vec: p.append(v/s)
	return p



parser = argparse.ArgumentParser()
parser.add_argument('genesdir', help='directory of genes (e.g. smallgenes)')
parser.add_argument('model', help='model file (e.g. models/worm.splicemodel')
parser.add_argument('--limit', required=False, type=int, default=100,
	metavar='<int>', help='limit number of transcripts [%(default)i]')
parser.add_argument('--something', action='store_true')
arg = parser.parse_args()

model = isoform2.read_splicemodel(arg.model)

for ff in glob.glob(f'{arg.genesdir}/*.fa'):
	gff = ff[:-2] + 'gff3'
	genome = Reader(fasta=ff, gff=gff)
	chrom = next(genome)
	#print(chrom.name, flush=True)

	# RNA counts
	rna = {}
	for f in chrom.ftable.features:
		if f.source != 'RNASeq_splice': continue
		sig = (f.beg-1, f.end-1)
		rna[sig] = f.score
	#print('rna', rna)

	# APC counts
	locus = isoform2.Locus(chrom.name, chrom.seq, model, limit=arg.limit)
	apc = {}
	for iso in locus.isoforms:
		for sig in iso.introns:
			if sig not in apc: apc[sig] = 0
			apc[sig] += iso.prob
	#print('apc', apc)

	# create intersection of APC and rna as probabilities
	rna2 = []
	apc2 = []
	for sig in rna.keys() & apc.keys():
		rna2.append(rna[sig])
		apc2.append(apc[sig])
	rna2 = list2prob(rna2)
	apc2 = list2prob(apc2)
	#print(rna2)
	#print(apc2)

	# do the test
	U1, p = mannwhitneyu(apc2, rna2)
	print(chrom.name, p, sep='\t', flush=True)
