import files
import os
import sys
import random


TYPES = ['high', 'low']

EMPTY_PDIC = {
	'A': 0,
	'C': 0,
	'T': 0,
	'G': 0
}

class GeneJunction:
	def __init__(self, strand, seq, typ, count, flank = 10):
		self.strand = strand
		self.seq = seq
		self.count = count
		self.type = typ
		self.kfreqs = {}

	def __str__(self):
		return (
			f'Strand: {self.strand}\n'
			f'\tSeq: {self.seq}\n\tCount: {self.count}\n\tType: {self.type}'
		)


	def find_kmerfreqs(self, k = 5):
		kmer_tot = self._ksearch(k)

		for kmer, count in self.kfreqs.items():
			self.kfreqs[kmer] = count / kmer_tot


	def _ksearch(self, k = 5):
		ktot = 0
		fullseq = ''.join(self.seq.split('..'))

		for i in range(0, len(fullseq) - k + 1):
			kmer = fullseq[i:i+k]
			if kmer not in self.kfreqs:
				self.kfreqs[kmer] = 0
			self.kfreqs[kmer] += 1
			ktot += 1

		return ktot

def make_genejunctions(l):
	fields = l.split()
	s = fields[0]
	j1 = fields[1]
	num1 = int(fields[2].rstrip('.0'))
	j2 = fields[3]
	num2 = int(fields[4].rstrip('.0'))

	if num1 > num2:
		hiseq, hicount = j1, num1
		loseq, locount = j2, num2
	else:
		hiseq, hicount = j2, num2
		loseq, locount = j1, num1

	for typ, seq, count in zip(TYPES, [hiseq, loseq], [hicount, locount]):
		yield GeneJunction(s, seq, typ, count)


def make_pwms(trlist, flank = 10):
	hi5pwm = []
	lo5pwm = []
	hi3pwm = []
	low3pwm = []

	for trgj in trlist:
		if trgj.type == 'high':

	return hi5pwm, lo5pwm, hi3pwm, lo3pwm



def make_sets(fname, select = 0.10):
	trains = []
	tests = []

	with files.getfp(fname) as fp:
		for line in fp:
			for gjunc in make_genejunctions(line):
				if random.random() >= 1 - select:
						trains.append(gjunc)
				else:
					tests.append(gjunc)

	return trains, tests


txtfile = sys.argv[1]
seed = sys.argv[2]

if seed.isdigit():
	random.seed(int(seed))

trainers, testers = make_sets(txtfile)


for junc in trainers:
	junc.find_kmerfreqs()


for junc in testers:
	junc.find_kmerfreqs()



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