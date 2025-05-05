import argparse
import glob

from grimoire.genome import Reader
import isoform2

def convert_dict2prob(d):
	s = sum(d.values())
	for k, v in d.items():	d[k] /= s

def compare(d1, d2):
	p = []
	q = []
	for sig in d1.keys() | d2.keys():
		if sig in d1: p.append(d1[sig])
		else:         p.append(0)
		if sig in d2: q.append(d2[sig])
		else:         q.append(0)
	return isoform2.intersection(p, q), isoform2.manhattan(p, q)

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

	# RNA counts
	rna = {}
	for f in chrom.ftable.features:
		if f.source != 'RNASeq_splice': continue
		# remove introns that are in the flanks (there can be some)
		if f.beg < 100 or f.end > len(chrom.seq) -100: continue
		sig = (f.beg-1, f.end-1)
		rna[sig] = f.score
	convert_dict2prob(rna)

	# create intron counts from annotation data
	gene = chrom.ftable.build_genes()[0] # there is only ever one
	ann = {} # all annoated transcripts equally weighted
	for tx in gene.transcripts():
		for f in tx.introns:
			sig = (f.beg-1, f.end-1)
			if sig not in ann: ann[sig] = 0
			ann[sig] += 1
	convert_dict2prob(ann)

	# APC counts
	locus = isoform2.Locus(chrom.name, chrom.seq, model, limit=arg.limit)
	apc = {}
	for iso in locus.isoforms:
		for sig in iso.introns:
			if sig not in apc: apc[sig] = 0
			apc[sig] += iso.prob
	convert_dict2prob(apc)

	print(chrom.name, gene.id, compare(rna, ann), compare(rna, apc), sep='\t', flush=True)



