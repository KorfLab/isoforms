# get data for figure showing isoform scoring
# using the simplest gene

import sys
import isoform2
#from isoform2 import Locus
import isoform
from isoform import Locus

fasta = sys.argv[1]
gff = sys.argv[2]

import argparse

parser = argparse.ArgumentParser(description='get data for figure '
	'showing how isoforms are scored')
parser.add_argument('fasta')
parser.add_argument('gff')
parser.add_argument('model')
parser.add_argument('--weights', required=False)

arg = parser.parse_args()

model = isoform2.read_splicemodel(arg.model)

weights = {
	'wacc': 1.0,
	'wdon': 1.0,
	'wexs': 1.0,
	'wins': 1.0,
	'wexl': 1.0,
	'winl': 1.0,
	'winf': 1.0,
} 

if arg.weights:
	with open(arg.weights, 'r') as fp:
		for line in fp.readlines():
			line = line.rstrip()
			line = line.split(',')
			weights['wacc'] = float(line[2])
			weights['wdon'] = float(line[3])
			weights['wexs'] = float(line[4])
			weights['wins'] = float(line[5])
			weights['wexl'] = float(line[6])
			weights['winl'] = float(line[7])
			weights['winf'] = float(line[8])
						
constraints = {
	'min_intron': 35,
	'min_exon': 25,
	'flank': 99
}

name, seq = next(isoform2.read_fasta(arg.fasta))

locus = Locus(name, seq, model, constraints=constraints, 
	weights=weights, limit=100)
	
acc = model['acc']	
don = model['don']
exl = model['exl']
inl = model['inl']
exs = model['exs']
ins = model['ins']

# APC considers the UTRs part of the exon!!!!!!!
# next modification: scan for ATG...


for iso in locus.isoforms:
	if iso.introns != [(232, 278)]: continue
	for i in iso.accs:
		s = isoform2.score_pwm(acc, iso.seq, i -len(acc)+1, memo=None)
		a_seq = iso.seq[i-len(acc)+1:i+1]
		print('Acceptor: ', a_seq, i, s)
	for i in iso.dons:
		s = isoform2.score_pwm(don, iso.seq, i, memo=None)
		d_seq = iso.seq[i:i+5]
		print('Donor: ', d_seq, i, s)
	for b, e in iso.exons:
		print(b, e, '###')
		elen = e - b -2 + 1
		slen = isoform2.score_len(exl, elen)
		smm = isoform2.score_markov(exs, iso.seq, b, e, memo=None)
		
		isoform.score_markov(exs, iso.seq, b, e, memo=None)
		
		print(iso.seq[b:e])
		
		# first kmer score is not ATGN, but 3 bp upstream of exon start
		print(f'{iso.seq[b+3:b+3+6]}...{iso.seq[e-5:e+1]}')
		
		print('Ex len: ', elen, slen)
		print('Ex mm: ', smm)
	for b, e in iso.introns:
		ilen = e - b + 1
		s = isoform2.score_len(inl, ilen)
		print('In len: ', ilen, s)
	print(iso.score)

	
seq = ''
with open(arg.fasta) as fp:
	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'): continue
		seq += line
		















