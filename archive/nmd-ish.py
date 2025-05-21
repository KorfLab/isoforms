import argparse
import re
import sys

import isoform2
from grimoire.genome import Reader

def display_isoform(gseq, tx, cds_beg, wrap=80, flank=99):
	# used for debugging
	print(iso.aaseq)
	rna = [' '] * len(gseq)
	for beg, end in tx.exons:
		for i in range(beg+1, end): rna[i] = 'X'
		rna[beg] = '['
		rna[end] = ']'
	for beg, end in tx.introns:
		for i in range(beg, end+1): rna[i] = '-'
	tseq = ''.join(rna)
	txcoor = 0
	tpos = []
	for nt in tseq:
		if nt != '-': txcoor += 1
		tpos.append(txcoor)

	cds = [' '] * len(gseq)
	x = cds_beg - tx.exons[0][0] # first atg
	while True:
		c0 = tx.rnaidx_to_dnaidx(x)
		c1 = tx.rnaidx_to_dnaidx(x+1)
		c2 = tx.rnaidx_to_dnaidx(x+2)
		if c0 is None: break
		if c1 is None: break
		if c2 is None: break
		codon = gseq[c0] + gseq[c1] + gseq[c2]
		aa = isoform2.GCODE[codon]
		cds[c1] = aa
		if codon == 'TAA' or codon == 'TAG' or codon == 'TGA': break
		x += 3
	cseq = ''.join(cds)

	for i in range(0, len(gseq), wrap):
		rbeg = i
		rend = i + wrap

		# genome coor
		for j in range(0, wrap, 10):
			print(f'{i+j+10:>10}', end='')
		print()
		for j in range(0, wrap, 10):
			print(' ' * 9, end='')
			print('|', end='')
		print()

		# sequences
		print(gseq[i:i+wrap])
		has_cds = False
		has_exon = False
		if re.search(r'\S', cseq[i:i+wrap]): has_cds = True
		if re.search(r'\S', tseq[i:i+wrap]): has_exon = True

		if has_cds: print(cseq[i:i+wrap])

		if has_exon:
			print(tseq[i:i+wrap])
			for j in range(0, wrap, 10):
				print(' ' * 9, end='')
				print('|', end='')
			print()

			for j in range(0, wrap, 10):
				p = i + j + 9
				if p >= len(tseq): break
				x = tpos[p]

				if tseq[p] == '.': print(' ' * 10, end='')
				else: print(f'{x:>10}', end='')
			print()
		print()
	sys.exit()


parser = argparse.ArgumentParser()
parser.add_argument('model', help='splicing model file (from smallgenes)')
parser.add_argument('fasta', help='fasta file (from smallgenes)')
parser.add_argument('gff', help='gff file (from smallgenes)')
parser.add_argument('--min-orf', type=int, default=25,
	help='minimum distance from start to stop [%(default)i]')
parser.add_argument('--nonstop', type=float, default=1e-3,
	help='degredation from no stop codon [%(default)g]')
parser.add_argument('--nmd', type=float, default=1e-3,
	help='degredation from NMD [%(default)g]')
parser.add_argument('--min-ejc', type=int, default=10,
	help='NMD minimum distance from stop to ejc [%(default)i]')
parser.add_argument('--utr', type=int, default=300,
	help='NMD long tail trigger length [%(default)i]')
parser.add_argument('--limit', type=int, default=100,
	help='maximum number of isoforms [%(default)i]')
parser.add_argument('--flank', type=int, default=99,
	help='flanking, non-coding sequence [%(default)i]')
parser.add_argument('--debug', action='store_true')
arg = parser.parse_args()

# create the locus
model = isoform2.read_splicemodel(arg.model)
reader = Reader(fasta=arg.fasta, gff=arg.gff)
region = next(reader)
locus = isoform2.Locus(region.name, region.seq, model, limit=arg.limit)

# find the canonical start codon (5' atg if there are more than one listed)
gene = region.ftable.build_genes()[0]
txs = gene.transcripts()
atgs = set()
for tx in txs:
	cdss = sorted(tx.cdss, key=lambda x: x.beg)
	atgs.add(cdss[0].beg -1)
cds_beg = sorted(list(atgs))[0]

# examine the isoforms to determine if they are NMD targets
prevs = [iso.prob for iso in locus.isoforms]
posts = []
rtypes = []
for iso in locus.isoforms:
	iso.translate(cds_beg)
	if arg.debug: display_isoform(region.seq, iso, cds_beg)
	if   iso.rnatype == 'non-stop': iso.prob *= arg.nonstop
	elif iso.rnatype == 'nmd-target': iso.prob *= arg.nmd
	posts.append(iso.prob)
	rtypes.append(iso.rnatype)

# recompute probabilities
total = sum([iso.prob for iso in locus.isoforms])
for iso in locus.isoforms:
	iso.prob /= total

print('rna-type', 'original', 'reduced', 'normalized', sep='\t')
for rtype, prev, post, iso in zip(rtypes, prevs, posts, locus.isoforms):
	print(f'{rtype}\t{prev:.3g}\t{post:.3g}\t{iso.prob:.3g}')
