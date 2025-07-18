# data set generator for distance from start vs. expression

import argparse
import glob
import sys

from grimoire.genome import Reader

parser = argparse.ArgumentParser()
parser.add_argument('fasta')
parser.add_argument('gff3')
parser.add_argument('--exon', type=int, default=30, metavar='<int>',
	help='length of exon flanks [%(default)i]')
parser.add_argument('--minscore', type=int, default=10000, metavar='<int>',
	help='minimum intron count [%(default)i]')
parser.add_argument('--maxdiff', type=float, default=2, metavar='<float>',
	help='maximum difference between adjacent introns [%(default).1f]')
arg = parser.parse_args()

for chrom in Reader(fasta=arg.fasta, gff=arg.gff3):
	# get all the RNA-seq data
	rna = {}
	for f in chrom.ftable.features:
		if f.source != 'RNASeq_splice': continue
		sig = f.beg, f.end
		rna[sig] = f.score

	for gene in chrom.ftable.build_genes():
		txs = gene.transcripts()
		if len(txs) == 0: continue
		tx = txs[0]

		# make sure all exons of transcript are in RNA-seq data and > arg.min
		all_good = True
		for f in tx.introns:
			sig = f.beg, f.end
			if sig not in rna or rna[sig] < arg.minscore:
				all_good = False
				break

		if not all_good: continue

		# create list of (dist, exp) tuple
		if gene.strand == '+':
			exons = tx.exons
			introns = tx.introns
		else:
			exons = tx.exons[::-1]
			introns = tx.introns[::-1]

		distexp = []
		tx_dist = 0
		for exon, intron in zip(tx.exons, tx.introns):
			sig = intron.beg, intron.end
			tx_dist += exon.end - exon.beg + 1
			distexp.append((tx_dist, rna[sig]))

		# remove too different
		kill = False
		for i in range(1, len(distexp)):
			d = distexp[i-1][1] / distexp[i][1]
			if d < 1/arg.maxdiff or d > arg.maxdiff:
				kill = True
				break
		if kill: continue

		# out
		for dist, exp in distexp:
			print(gene.id, dist, exp, sep='\t')
