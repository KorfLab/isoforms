import argparse
import glob
import re
import statistics

from grimoire.genome import Reader

parser = argparse.ArgumentParser()
parser.add_argument('smallgenes')
arg = parser.parse_args()

all_distances = []
for gff3 in glob.glob(f'{arg.smallgenes}/*.gff3'):
	m = re.search(r'(ce\.\d+\.\d+)\.gff3', gff3)
	locus = m.group(1)
	fasta = f'{arg.smallgenes}/{locus}.fa'

	for chrom in Reader(fasta=fasta, gff=gff3):
		# get all the RNA-seq data
		rna = {}
		for f in chrom.ftable.features:
			if f.source != 'RNASeq_splice': continue
			sig = f.beg, f.end
			rna[sig] = f.score
	
		
		for gene in chrom.ftable.build_genes():
			tx = gene.transcripts()
			tx = tx[0]
			if len(tx.introns) < 2: continue
			intron_counts = []
			total_counts = 0
			for intron in tx.introns:
				intron_counts.append(rna[intron.beg, intron.end])
				total_counts += rna[intron.beg, intron.end]
			intron_freqs = []
			for n in intron_counts: intron_freqs.append(n/total_counts)
			distances = []
			for i in range(len(intron_freqs)):
				for j in range(i+1, len(intron_freqs)):
					distances.append(abs(intron_freqs[i] - intron_freqs[j]))
			all_distances.extend(distances)
print(statistics.mean(all_distances), statistics.stdev(all_distances))			
