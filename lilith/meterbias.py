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
				self.kcounts[kmer] = 0
			self.kcounts[kmer] += 1

def make_genejunctions(l):
	fields = l.split()
	s = fields[0]
	j = fields[1]
	num = int(fields[2].rstrip('.0'))
	t = fields[3]

	
	return GeneJunction(s, j, num, t)


def make_pwms(trlist, flank = 10):
	hi5pwm = []
	lo5pwm = []
	hi3pwm = []
	low3pwm = []
'''
	for trgj in trlist:
		if trgj.type == 'high':
	return hi5pwm, lo5pwm, hi3pwm, lo3pwm
'''

def score_juncs(meter, telist, tmdic, k):
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

	hit = 0
	mis = 0

	for tejunc in telist:
		if tejunc.type not in TYPES:
			continue

		logprob = 0 # HIGH / LOW
		for i in range(len(tejunc.seq) - k + 1):
			kme = tejunc.seq[i:i+k]

			if kme not in tmdic['HIGH'] or kme not in tmdic['LOW']:
				continue

			logprob += math.log2(tmdic['HIGH'][kme] / tmdic['LOW'][kme])



		if logprob > 0:
			result = 'HIGH'
		else:
			result = 'LOW'

		if result == tejunc.type:
			hit += 1
		else:
			mis += 1

	return hit, mis




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

	for tmstrand, tmkmers in tmdic.items():
		for typ, ks in tmkmers.items():
			for km in ks.keys():
				ks[km] = ks[km] / typetots[tmstrand][typ]
	return tmdic



txtfile = sys.argv[1]
exonfile = sys.argv[2]
ksize = int(sys.argv[3])
trials = int(sys.argv[4])
seed = sys.argv[5]


if seed.isdigit():
	random.seed(int(seed))

results = 0

for trial in range(trials):
	trainers, testers = make_sets(txtfile)


	for junc in trainers:
		junc.ksearch()


	for junc in testers:
		junc.ksearch()



	meter = json.load(files.getfp(exonfile))

	trainmeter = make_trainmeter(trainers, ksize)

	print(json.dumps(trainmeter, sort_keys=True, indent=4))

	sys.exit()
	hits, misses = score_juncs(meter, testers, trainmeter, k=ksize)

	print(f'Trial {trial} percentage of correct junction IDs: {hits / (hits + misses)}')

	results += hits / (hits + misses)


print(f'Average accuracy across {trials} trials: {results/trials}')
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