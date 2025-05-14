import glob
import os
import statistics
import sys

print('\t'.join( ('name', 'gene', 'loc', 'str', 'base', 'opt', 'introns',
	'don', 'acc', 'emm', 'imm', 'elen', 'ilen', 'icost') ))

data = []
for ff in (glob.glob(f'apc/*.fa')):

	# fasta
	with open(ff) as fp:
		line = fp.readline()
		name, loc, strand, gene = line.split()
		name = name[1:]

	# gff
	introns = 0
	with open(f'apc/{name}.gff3') as fp:
		for line in fp:
			f = line.split()
			if f[1] == 'RNASeq_splice': introns += 1

	# isoformer base performance
	with open(f'build/isoformer/{name}.txt') as fp:
		base = fp.readline().rstrip()

	# optiso optimized performance
	with open(f'build/optiso/{name}.txt') as fp:
		line = fp.readline()
		opt, don, acc, emm, imm, elen, ilen, icost, n = line.split();

	# report
	print('\t'.join( (name, gene, loc, strand, base, opt,
		str(introns), don, acc, emm, imm, elen, ilen, icost)))
