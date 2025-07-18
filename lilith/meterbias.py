import files
import json
import math
import os
import sys
import random


TYPES = ['HIGH', 'LOW']

EMPTY_PDIC = {
	'A': 0,
	'C': 0,
	'T': 0,
	'G': 0
}

class GeneJunction:
	def __init__(self, strand, seq, count, typ, flank = 10):
		self.strand = strand
		self.seq = ''.join(seq.split('..'))
		self.count = count
		self.type = typ
		self.kcounts = {}

	def __str__(self):
		return (
			f'Strand: {self.strand}\n'
			f'\tSeq: {self.seq}\n\tCount: {self.count}\n\tType: {self.type}'
		)

	def ksearch(self, k = 5):
		ktot = 0

		for i in range(0, len(self.seq) - k + 1):
			kmer = self.seq[i:i+k]
			if kmer not in self.kcounts:
				#print(kmer)
				self.kcounts[kmer] = 0
			self.kcounts[kmer] += 1

def make_genejunctions(l):
	fields = l.split()
	s = fields[0]
	j = fields[1]
	num = int(fields[2].rstrip('.0'))
	t = fields[3]

	
	return GeneJunction(s, j, num, t)

def make_pwms(trlist, ntf, flank = 10):
	pluspwm = [EMPTY_PDIC.copy() for i in range(len(trlist[0].seq))]
	minuspwm = [EMPTY_PDIC.copy() for i in range(len(trlist[0].seq))]

	for trg in trlist:
		if trg.strand == '+':
			for i in range(len(trg.seq)):
				pluspwm[i][trg.seq[i]] += 1
		else:
			for i in range(len(trg.seq)):
				minuspwm[i][trg.seq[i]] += 1

	for posp, posn in zip(pluspwm, minuspwm):
		for k1, k2 in zip(posp.keys(), posn.keys()):
			if posp[k1] == 0:
				posp[k1] += 1
			if posn[k2] == 0:
				posn[k2] += 1

		ptot = sum(list(posp.values()))
		ntot = sum(list(posn.values()))

		if ptot == 0:
			ptot += 1

		if ntot == 0:
			ntot += 1

		for nt in posp:
			posp[nt] /= ptot
			posn[nt] /= ntot



	for i in range(len(pluspwm)):
		for nt in pluspwm[i].keys():

			pluspwm[i][nt] = math.log2(pluspwm[i][nt] / ntf[nt])
			minuspwm[i][nt] = math.log2(minuspwm[i][nt] / ntf[nt])


	return pluspwm, minuspwm


def score_juncs(telist, tmdic, k):
	'''
	hit = {
		'HIGH': 0,
		'LOW': 0
	}

	mis = {
		'HIGH': 0,
		'LOW': 0
	}
	'''
	results = {
		'+': {'hits': 0, 'misses': 0},
		'-': {'hits': 0, 'misses': 0}
	}

	for tejunc in telist:
		if tejunc.type not in TYPES:
			continue

		logprob = 0 # HIGH / LOW
		for i in range(len(tejunc.seq) - k + 1):
			kme = tejunc.seq[i:i+k]
			#print(kme)
			if kme not in tmdic[tejunc.strand]['HIGH'] or kme not in tmdic[tejunc.strand]['LOW']:
				continue

			logprob += math.log2(tmdic[tejunc.strand]['HIGH'][kme] / tmdic[tejunc.strand]['LOW'][kme])



		if logprob > 0:
			result = 'HIGH'
		else:
			result = 'LOW'

		if result == tejunc.type:
			results[tejunc.strand]['hits'] += 1
		else:
			results[tejunc.strand]['misses'] += 1
	return results




def make_sets(fname, select = 0.10):
	trains = []
	tests = [] 

	with files.getfp(fname) as fp:
		for line in fp:
			gjunc = make_genejunctions(line)
			if random.random() >= 1 - select:
					trains.append(gjunc)
			else:
				tests.append(gjunc)

	return trains, tests

def make_trainmeter(trlist, k):
	typetots = {
		'+': {
			'HIGH': 0,
			'LOW': 0
		},
		'-': {
			'HIGH': 0,
			'LOW': 0
		}
	}
	tmdic = {
		'+': {
			'HIGH': {},
			'LOW': {}
		},
		'-': {
			'HIGH': {},
			'LOW': {}
		}
	}

	for trj in trlist:
		if trj.type not in TYPES:
			continue

		for i in range(len(trj.seq) - k + 1):
			kme = trj.seq[i:i+k]

			if kme not in tmdic[trj.strand][trj.type]:
				tmdic[trj.strand][trj.type][kme] = 0

			tmdic[trj.strand][trj.type][kme] += 1

			typetots[trj.strand][trj.type] += 1

	for typ, ks in tmdic['+'].items():
		for k, v in ks.items():
			ks[k] = v / typetots['+'][typ]

	for typ, ks in tmdic['-'].items():
		for k, v in ks.items():
			ks[k] = v / typetots['-'][typ]
	
	return tmdic




def find_ntfreq(infile):
	ntdic = {
		'A': 0,
		'C': 0,
		'G': 0,
		'T': 0
	}

	total = 0

	with files.getfp(infile) as fp:
		for line in fp:
			line = line.strip()
			for nt in line:
				if nt not in ntdic.keys(): continue
				ntdic[nt] += 1
				total += 1

	for nt in ntdic.keys():
		ntdic[nt] /= total

	return ntdic








txtfile = sys.argv[1]
exonfile = sys.argv[2]
ksize = int(sys.argv[3])
trials = int(sys.argv[4])
seed = sys.argv[5]


if seed.isdigit():
	random.seed(int(seed))

strrestults = {
	'+': 0.0,
	'-': 0.0
}

#meter = json.load(files.getfp(exonfile))

#ntfreqs = find_ntfreq(exonfile)

#print(ntfreqs)

for trial in range(trials):
	trainers, testers = make_sets(txtfile)
	print('wedidit')

	for junc in trainers:
		junc.ksearch()


	for junc in testers:
		junc.ksearch()

	print('wediditAGAIN')


	trainmeter = make_trainmeter(trainers, ksize)
	print('wediditAGAINAGAIN')

	#print(json.dumps(trainmeter, sort_keys=True, indent=4))
	#ppwm, npwm = make_pwms(trainers, ntfreqs)
	#sys.exit()

	res = score_juncs(testers, trainmeter, k=ksize)

	for key in res.keys():
		acc = res[key]['hits'] / (res[key]['hits'] + res[key]['misses'])
		print(f'Trial {trial}, strand {key} percentage of correct junction IDs: {acc}')


		strrestults[key] += acc


print(f'Average accuracy across {trials}:')
print(f'+: {strrestults['+'] / trials}')
print(f'-: {strrestults['-'] / trials}')
'''
data struct:

{
	strand: -
	high: {
		seq: seq,
		count: count
	},
	low: {
		seq: seq,
		count: count
	}
}
'''