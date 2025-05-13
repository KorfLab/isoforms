import argparse
import csv

parser = argparse.ArgumentParser(description='create pwm for figures')

parser.add_argument('introns', type=str, 
	help='text file with intron sequences')
parser.add_argument('--pwm_size', required=False, type=int, default=30, 
	help='length of pwm %(default)i')
	
args = parser.parse_args()

dsides = []
asides = []
with open(args.introns, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		dside = line[:args.pwm_size]
		aside = line[-args.pwm_size:]
		dsides.append(dside)
		asides.append(aside)

def build_pwm(seqs, pwm_size):
	
	counts = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for x in range(pwm_size)]
	
	for seq in seqs:
		for i, nt in enumerate(seq):
			counts[i][nt] += 1
			
	pwm = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for x in range(pwm_size)]
	
	for i, site in enumerate(counts):
		for nt in site:
			pwm[i][nt] = site[nt]/len(seqs)
			
	return pwm
	
dpwm = build_pwm(dsides, args.pwm_size)
apwm = build_pwm(asides, args.pwm_size)

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
	
	
	
	
	
	
	
	
	
	
	
