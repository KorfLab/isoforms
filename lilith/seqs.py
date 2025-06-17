# sequencing bias for spliced sequences
# why do some introns get sequenced more often than others
# what flanking sequences affect sequencing/fragmentation/alignment bias
# what kmers are associated with high vs low
# grab sequence of exons neighboring introns
# are they higher or lower
# filters on bad genes that are too different

import argparse
import glob

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
		rna[sig] = round(float(f.score))

	# get all the uniqe intron pairs and assign scores
	pair = {}
	for gene in chrom.ftable.build_genes():

		for tx in gene.transcripts():
			short = False
			txscos = []

			if len(tx.introns) < 2:
				continue
			
			for i in range(len(tx.exons) - 1):
				dist = tx.exons[i].end - tx.exons[i].beg
				if dist < (2 * arg.exon):
					short = True

			if short == True:
				continue


			for i in range(len(tx.introns)):
				#b1 = tx.introns[i-1].beg
				#e1 = tx.introns[i-1].end
				b2 = tx.introns[i].beg
				e2 = tx.introns[i].end
				#sig1 = b1, e1
				sig2 = b2, e2
				# require intron is validate by rna-seq
				#if sig1 not in rna: continue
				if sig2 not in rna: continue
				# require some level of rna expression
				#if rna[sig1] < arg.minscore: continue
				if rna[sig2] < arg.minscore: continue
				# remove really different values (probably errors)
				#if rna[sig1] / rna[sig2] > arg.maxdiff: continue
				#if rna[sig2] / rna[sig1] > arg.maxdiff: continue
				#pair[(gene.strand, sig1, sig2)] = (rna[sig1], rna[sig2])
				txscos.append(rna[sig2])
			print(txscos)

'''

	# get all the flanks
	for (strand, (b1, e1), (b2, e2)), (s1, s2) in pair.items():
		e1a = chrom.seq[b1-arg.exon-1:b1+1]
		e1b = chrom.seq[e1-2:e1+arg.exon]
		e2a = chrom.seq[b2-arg.exon-1:b2+1]
		e2b = chrom.seq[e2-2:e2+arg.exon]
		if len(e1a) != arg.exon +2: continue
		if len(e1b) != arg.exon +2: continue
		if len(e2a) != arg.exon +2: continue
		if len(e2b) != arg.exon +2: continue

		#print(f'{strand} {e1a}..{e1b} {s1} {e2a}..{e2b} {s2}')
		#print(f'{strand} {b1}-{e1}:{s1} {b2}-{e2}:{s2}')

'''