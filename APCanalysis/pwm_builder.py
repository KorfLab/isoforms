import argparse
import csv
import math

parser = argparse.ArgumentParser(description='create pwm for figures')

parser.add_argument('introns', type=str, 
	help='text file with intron sequences')
parser.add_argument('--pwm_size', required=False, type=int, default=30, 
	help='length of pwm %(default)i')
	
args = parser.parse_args()

donor_sides = []
acceptor_sides = []
with open(args.introns, 'r') as fp:
	for line in fp.readlines():
		seq = line.rstrip()
		donor_side = seq[:args.pwm_size]
		acceptor_side = seq[-args.pwm_size:]
		intron_seq = seq[5:len(seq)-5]
		
		donor_sides.append(donor_side)
		acceptor_sides.append(acceptor_side)	
		
# introns/exons are NOT weighted using this method
# each one counts the same as everyone else
# WormBase does not have counts
# RNASeq_splice does have counts
def build_pwm(seqs, pwm_size):
	
	counts = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for x in range(pwm_size)]
	
	for seq in seqs:
		for i, nt in enumerate(seq):
			counts[i][nt] += 1
			
	ppm = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for x in range(pwm_size)]
	
	for i, site in enumerate(counts):
		for nt in site:
			ppm[i][nt] = site[nt]/len(seqs)
			
	pwm = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for x in range(pwm_size)]
			
	for i, site in enumerate(ppm):
		uncertainty = 0
		for item in site.items():
			uncertainty += item[1] * math.log2(item[1])
		uncertainty = -uncertainty
		info_content = 2 - uncertainty
		for nt in site:
			pwm[i][nt] = site[nt] * info_content

	return pwm

dpwm = build_pwm(donor_sides, args.pwm_size)
apwm = build_pwm(acceptor_sides, args.pwm_size)

with open('donor_side_pwm.csv', 'w') as csvfile:
	dwriter = csv.writer(csvfile)
	for dsite in dpwm:
		row = [dsite[x] for x in dsite]
		dwriter.writerow(row)
	
with open('acceptor_side_pwm.csv', 'w') as csvfile:
	awriter = csv.writer(csvfile)
	for asite in apwm:
		row = [asite[x] for x in asite]
		awriter.writerow(row)

	
	
	
	
	
	
	
	
	
	
