import argparse

parser = argparse.ArgumentParser()
parser.add_argument('genome')
parser.add_argument('gff')

args = parser.parse_args()

def revcomp(seq):

	comps = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

	rev = ''
	for i in range(1, len(seq) + 1):
		rev += comps[seq[-i]]

	return rev

## collect subsequences
introns = []
exons = []
with open(args.gff, 'r') as gp:
	for line in gp.readlines():
		line = line.rstrip()
		line = line.split('\t')
		if len(line) < 9: continue
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

accs = []
dons = []
iseqs = []

for i in introns:
	beg = i[1]-1
	end = i[2]
	if i[3] == '-':
		miseq = chroms[i[0]][beg-5:end+5]
		iseq = revcomp(miseq)
	else:
		iseq = chroms[i[0]][beg-5:end+5]
	aseq = iseq[-6:]
	dseq = iseq[:5]
	if aseq.endswith('AG'):
		accs.append(aseq)
	if dseq.startswith('GT'):
		dons.append(dseq)
	iseqs.append(iseq)

with open('introns.txt', 'w') as ifp:
	for i in iseqs:
		ifp.write(f'{i}\n')





