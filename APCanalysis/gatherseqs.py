import argparse
import glob

from grimoire.genome import Reader

parser = argparse.ArgumentParser()
#parser.add_argument('genome')
#parser.add_argument('gff')

parser.add_argument('dir')
parser.add_argument('annot')

args = parser.parse_args()

# this method is only getting the cannonical intron
for ff in glob.glob(f'{args.dir}/*.fa'):
	gf = ff[:-2] + 'gff3'
	genome = Reader(gff=gf, fasta=ff)
	tx = next(genome).ftable.build_genes()[0].transcripts()[0]
	for f in tx.exons:
		print(f)
	print('#####')
	for f in tx.introns:
		print(f)
	break
	
'''
with open(args.annot, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		print(line)
'''

'''
introns = []
exons = []
with open(args.gff, 'r') as gp:
	for line in gp.readlines():
		line = line.rstrip()1
		line = line.split('\t')
		if line[2] == 'intron':
			intron = (line[0], int(line[3]), int(line[4]), line[6])
			introns.append(intron)
		if line[2] == 'exon':
			exon = (line[0], int(line[3]), int(line[4]), line[6])
			exons.append(exon)

chroms = {}
with open(args.genome, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			chrom = line.split(' ')[0][1:]
		else:
			if chrom not in chroms:
				chroms[chrom] = line
			else:
				chroms[chrom] += line

def revcomp(seq):

	comps = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

	rev = ''
	for i in range(1, len(seq) + 1):
		rev += comps[seq[-i]]

	return rev

accs = []
dons = []
eseqs = []
iseqs = []
exonsum = 0

for i in introns:
	beg = i[1]-1
	end = i[2]
	if i[3] == '-':
		miseq = chroms[i[0]][beg:end]
		iseq = revcomp(miseq)
	else:
		iseq = chroms[i[0]][beg:end]
	aseq = iseq[-6:]
	dseq = iseq[:5]
	if not aseq.endswith('AG'): continue
	accs.append(aseq)
	if not dseq.startswith('GT'): continue
	dons.append(dseq)
	iseqs.append(iseq)

for e in exons:
	beg = e[1]-1
	end = e[2]
	if e[3] == '-':
		meseq = chroms[e[0]][beg:end]
		eseq = revcomp(meseq)
	else:
		eseq = chroms[i[0]][beg:end]
	exonsum += len(eseq)
	eseqs.append(eseq)
'''
# i need to ge only intron and exon seqs from annotations.gff3







	
# need to weight the introns, so highly expressed ones are represented more
'''
for e in exons:
	beg = e[1]-1
	end = e[2]
	eseq = chroms[i[0]][beg:end]
	exonsum += len(eseq)
	eseqs.append(eseq)
'''
