import argparse

parser = argparse.ArgumentParser(description='get splice site '
	'sequences for each isoform')
parser.add_argument('gff')
parser.add_argument('fasta')

args = parser.parse_args()

iso_count = 0
isos = {}
with open(args.gff, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		if line.startswith('#'): continue
		if line == '': continue
		line = line.split('\t')
		if line[2] == 'mRNA':
			iso_count += 1
		if line[2] == 'intron':
			if iso_count not in isos:
				isos[iso_count] = [(int(line[3]), int(line[4]))]
			else:
				isos[iso_count].append((int(line[3]), int(line[4])))
				
seq = ''
with open(args.fasta, 'rt') as fp:
	for line in fp:	
		line = line.rstrip()
		if line.startswith('>'): continue
		seq += line
		
ss_seqs = {}
for item in isos.items():
	ss_seqs[item[0]] = {}
	for intron in item[1]:
		don = intron[0]
		acc = intron[1]
		dseq = seq[don-1:don+4]
		aseq = seq[acc-6:acc]
		ss_seqs[item[0]][intron] = (dseq, aseq)

print(ss_seqs[1])








