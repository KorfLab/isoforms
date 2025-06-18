import argparse
import glob
import statistics

from grimoire.genome import Reader

parser = argparse.ArgumentParser()
parser.add_argument('fasta')
parser.add_argument('gff3')
parser.add_argument('--resolution', type=int, default=500, metavar='<int>',
	help='histogram bin size [%(default)i]')
parser.add_argument('--minscore', type=int, default=10000, metavar='<int>',
	help='minimum intron count [%(default)i]')
parser.add_argument('--maxdiff', type=float, default=2, metavar='<float>',
	help='maximum difference between adjacent introns [%(default).1f]')
arg = parser.parse_args()

distance = {}
for chrom in Reader(fasta=arg.fasta, gff=arg.gff3):
	# get all the RNA-seq data
	rna = {}
	for f in chrom.ftable.features:
		if f.source != 'RNASeq_splice': continue
		sig = f.beg, f.end
		rna[sig] = f.score

	for gene in chrom.ftable.build_genes():
		tx = gene.transcripts()
		if len(tx) == 0: continue
		tx = tx[0]
		if tx.issues: continue
		if len(tx.introns) < 2: continue
		
		# all introns must be in RNA seq data
		confirmed = True
		total_rna = 0
		for intron in tx.introns:
			sig = intron.beg, intron.end
			if sig not in rna:
				confirmed = False
				break
			total_rna += rna[sig]
		if not confirmed: continue
		
		# average intron count must meet threshold
		ave_rna = total_rna / len(tx.introns)
		if ave_rna < arg.minscore: continue
		
		# remove genes that might be fusions
		fusion_suspect = False
		for i in range(1, len(tx.introns)):
			prev = rna[tx.introns[i-1].beg, tx.introns[i-1].end]
			this = rna[tx.introns[i].beg, tx.introns[i].end]
			if prev / this > arg.maxdiff or this/prev > arg.maxdiff:
				fusion_suspect = True
				break
		if fusion_suspect: continue
		
		# count distances
		for intron in tx.introns:
			intron_mid = (intron.beg + intron.end) // 2
			if tx.strand == '+':
				d = (tx.end - intron_mid) // arg.resolution
			else:
				d = (intron_mid - tx.beg) // arg.resolution
			if d not in distance: distance[d] = []
			distance[d].append(rna[intron.beg, intron.end] / ave_rna)
			
for d in sorted(distance):
	print(d * arg.resolution, len(distance[d]), statistics.mean(distance[d]))