import isoform2
# isoform2 scoring is incorrect
#from isoform2 import Locus
import isoform
from isoform import Locus
from grimoire.genome import Reader

import argparse

parser = argparse.ArgumentParser(description='get data for figure '
	'showing how isoforms are scored')
parser.add_argument('fasta')
parser.add_argument('gff')
parser.add_argument('model')
parser.add_argument('apc_gff')
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
	total = 0
	# for ce.1.5
	#if iso.introns != [(190, 235)]: continue
	# for ce.3.1
	#if iso.introns != [(182, 232), (308, 364), (462, 515)]: continue
	# for ce.1.1
	#if iso.introns != [(232, 278)]: continue
	# for ce.2.1
	if iso.introns != [(187, 241), (286, 333)]: continue
	for i in iso.accs:
		sacc = isoform2.score_pwm(acc, iso.seq, i -len(acc)+1, memo=None)
		total += sacc
		a_seq = iso.seq[i-len(acc)+1:i+1]
		print('Acceptor: ', a_seq, i, sacc)
	for i in iso.dons:
		sdon = isoform2.score_pwm(don, iso.seq, i, memo=None)
		total += sdon
		d_seq = iso.seq[i:i+5]
		print('Donor: ', d_seq, i, sdon)
	for b, e in iso.exons:
		elen = e - b + 1
		selen = isoform2.score_len(exl, elen)
		total += selen
		semm = isoform2.score_markov(exs, iso.seq, b, e, memo=None)
		total += semm
		print('Exon: ', f'{iso.seq[b:b+6]}...{iso.seq[e-5:e+1]}', f'{b, e}')
		print('Ex len: ', elen, selen)
		print('Ex mm: ', semm, '$$')
	for b, e in iso.introns:
		ilen = e - b + 1
		silen = isoform2.score_len(inl, ilen)
		total += silen
		# don/acc seq not part of intron mm 
		simm = isoform2.score_markov(ins, iso.seq, b + len(don), 
			e - len(acc), memo=None)
		total += simm
		print('Intron: ', f'{iso.seq[b+5:b+11]}...{iso.seq[e-11:e-5]}')
		print('In len: ', ilen, silen)
		print('In mm: ', simm)
	total += model['inf'] * len(iso.introns)
	print('iso score: ', iso.score, 'total score: ', total)

print('#####')

weight = []
total = 0
count = 0
for tx in locus.isoforms:
	w = 2 ** tx.score
	weight.append(w)
	total += w
	if count <= 10:
		print(tx.score)
	count += 1
prob = [w / total for w in weight]

print('#####')

count = 0 
for p in prob:
	if count <= 10:
		print(p)
	count += 1



'''
10.907457181792099
5.663082135803814
2.381151251216277
2.281113512556373
9.877905709654959
14.08864071944679
2.2309530782449323
9.856043635846964
-6.252390479016578
'''


print('##### Distance example #####')



i1 = isoform.get_introns(arg.apc_gff)
i2 = isoform.get_introns(arg.gff)

d1 = 0
for intron in i1.keys() | i2.keys():
	if intron in i1 and intron in i2: 
		d1 += abs(i1[intron] - i2[intron])
		print(intron, f'{i1[intron]:.2e}', f'{i2[intron]:.2e}', 
				f'{i1[intron] - i2[intron]:.2e}')
	elif intron in i1: 
		d1 += i1[intron]
		print(intron, f'{i1[intron]:.2e}', 0, f'{i1[intron]:.2e}')
	else: 
		d1 += i2[intron]
		print(intron, 0, f'{i2[intron]:.2e}', f'{i1[intron]:.2e}')
	
dist = d1 / 2

print(d1)
print(dist)

'''
for i in i1.items():
	print(i)

for tx in locus.isoforms:
	for intron in tx.introns:
		print(intron, tx.prob)
'''

'''
seq = ''
with open(arg.fasta) as fp:
	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'): continue
		seq += line
		
print('##########')

print(model['exs']['mm']['AAAA'])
print(model['exs']['mm']['AAAT'])
print(model['exs']['mm']['AATG'])

print(model['exl']['tail'])


reader = Reader(fasta=arg.fasta, gff=arg.gff)
region = next(reader)

# find the canonical start codon
gene = region.ftable.build_genes()[0]
txs = gene.transcripts()
atgs = set()
for tx in txs:
	cdss = sorted(tx.cdss, key=lambda x: x.beg)
	atgs.add(cdss[0].beg - 1)
cds_beg = sorted(atgs)[0]

nonstop = 1e-3
nmd = 1e-3

prevs = [iso.prob for iso in locus.isoforms]
for iso in locus.isoforms:
	iso.translate(cds_beg)
	if iso.rnatype == 'non-stop':
		iso.prob *= nonstop
	elif iso.rnatype == 'nmd-target':
		iso.prob *= nmd

# recompute probabilities
total = sum(iso.prob for iso in locus.isoforms)
if total > 0:
	for iso in locus.isoforms:
		iso.prob /= total


# x is the start codon coordinate
# x += 3 finds the stop codon coordinate

# isoforms can also be labeled as non-stop

# what does splice_after_stop = -1 mean?
# this value is continuously updated
'''
'''
ce.2.1
APC.base first isoform
prob 0.9231
gene/mRNA 100 425
exon 100 187
exon 243 286
exon 335 425
intron 188 242
intron 287 334
'''



