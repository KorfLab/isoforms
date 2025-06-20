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
import sys


def find_median(l):
	m = len(l) % 2
	m_index = len(l) // 2

	if len(l) == 0:
		return None
	
	if m != 0:
		return l[m_index]
	else:
		return (l[m_index-1] + l[m_index]) / 2




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
	for gene in chrom.ftable.build_genes():
		txs = []
		for tx in gene.transcripts():
			short = False
			txscos = []
			junc = {}

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
				b = tx.introns[i].beg
				e = tx.introns[i].end
				#sig1 = b1, e1
				sig = b, e
				# require intron is validate by rna-seq
				#if sig1 not in rna: continue
				if sig not in rna: continue
				# require some level of rna expression
				if rna[sig] < arg.minscore: continue
				junc[(gene.strand, sig)] = rna[sig]
				txscos.append(rna[sig])

			med = find_median(txscos)

			for k, v in junc.items():
				if junc[k] < med:
					junc[k] = (v, 'LOW')
				else:
					junc[k] = (v, 'HIGH')

			txs.append(junc)

		# get all the flanks
		for tx in txs:
			for (strand, (b, e)), (s, typ) in tx.items():
				ebeg = chrom.seq[b-arg.exon-1:b+1]
				eend = chrom.seq[e-2:e+arg.exon]
				if len(ebeg) != arg.exon +2: continue
				if len(eend) != arg.exon +2: continue

				print(f'{strand} {ebeg}..{eend} {s} {typ}')

