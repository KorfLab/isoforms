import argparse

parser = argparse.ArgumentParser(
	description='gets donor and acceptor sequences for modelbuilder')
parser.add_argument('introns', type=str, help='list of introns')

args = parser.parse_args()

# hard coded based on pwm of all intron sequences in genome
# donor length 8, 2 bp upstream of GT and 4 bp downstream
# acceptor length 8, 6 bp upsream of AG

introns = []
with open(args.introns) as fp:
	for line in fp:
		line = line.rstrip()
		introns.append(line)
		
d_seqs = []
a_seqs = []
for iseq in introns:
	d_seq = iseq[3:11]
	a_seq = iseq[-13:-5]
	d_seqs.append(d_seq)
	a_seqs.append(a_seq)
	
with open('donors.txt', 'w') as fp:
	for seq in d_seqs:
		fp.write(f'{seq}\n')

with open('acceptors.txt', 'w') as fp:
	for seq in a_seqs:
		fp.write(f'{seq}\n')
		
with open('splice_sites.txt', 'w') as fp:
	for donor, acceptor in zip(d_seqs, a_seqs):
		fp.write(f'{donor},{acceptor}\n')
