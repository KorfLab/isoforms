import argparse
import glob
import math
import random
import os
import statistics
import sys

from grimoire.genome import Reader
import isoform

def pwm(size):
	pwm = []
	for _ in range(size):
		pwm.append({'A':0, 'C':0, 'G':0, 'T':0})
	return pwm

def print_pwm(name, pwm):
	print(f'% PWM {name} {len(pwm)}')
	for pos in pwm:
		total = sum(pos.values())
		vals = []
		for nt in pos:
			print(f'{pos[nt]/total:.6f} ', end='')
			if pos[nt] != 0: vals.append(pos[nt]/total)
		print(f'--> {2 - isoform.entropy(vals)}')
	print()

#######
# CLI #
#######

parser = argparse.ArgumentParser(description='make pwms several ways')
parser.add_argument('apc', help='path to smallgenes directory')
parser.add_argument('--don', type=int, default=5)
parser.add_argument('--acc', type=int, default=6)
arg = parser.parse_args()

don_gene = pwm(arg.don)
acc_gene = pwm(arg.acc)

don_rna1 = pwm(arg.don)
acc_rna1 = pwm(arg.acc)

don_rnaa = pwm(arg.don)
acc_rnaa = pwm(arg.acc)

don_rnar = pwm(arg.don)
acc_rnar = pwm(arg.acc)

for ff in glob.glob(f'{arg.apc}/*.fa'):
	gf = ff[:-2] + 'gff3'
	genome = Reader(gff=gf, fasta=ff)
	chrom = next(genome)
	gene = chrom.ftable.build_genes()[0]
	tx = gene.transcripts()[0]

	# gene-based counting
	for intron in tx.introns:
		iseq = intron.seq_str()
		donor = ''
		for i in range(arg.don):
			nt = iseq[i]
			donor += nt
			don_gene[i][nt] += 1
		print(donor)
		for i in range(arg.acc):
			nt = iseq[-arg.acc+i]
			acc_gene[i][nt] += 1

	# transcript-based counting
	introns = []
	total_score = 0
	for f in chrom.ftable.features:
		if f.source != 'RNASeq_splice': continue
		#print(f)
		introns.append(f)
		total_score += f.score

	for intron in introns:
		iseq = intron.seq_str()

		# one-based counting
		#donor = ''
		for i in range(arg.don):
			nt = iseq[i]
			#donor += nt
			don_rna1[i][nt] += 1
		for i in range(arg.acc):
			nt = iseq[-arg.acc+i]
			acc_rna1[i][nt] += 1
		#print(donor)
		# absolute counting
		for i in range(arg.don):
			nt = iseq[i]
			don_rnaa[i][nt] += intron.score
		for i in range(arg.acc):
			nt = iseq[-arg.acc+i]
			acc_rnaa[i][nt] += intron.score

		# relative counting
		for i in range(arg.don):
			nt = iseq[i]
			don_rnar[i][nt] += intron.score / total_score
		for i in range(arg.acc):
			nt = iseq[-arg.acc+i]
			acc_rnar[i][nt] += intron.score / total_score

print_pwm('don_gene', don_gene)
#print_pwm('don_rna_one', don_rna1)
'''
print_pwm('don_rna_abs', don_rnaa)
print_pwm('don_rna_rel', don_rnar)
'''
