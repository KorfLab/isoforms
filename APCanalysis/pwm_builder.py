import argparse
import csv
import math
import apc_analysis as aa

parser = argparse.ArgumentParser(description='create pwm for figures')

parser.add_argument('introns', type=str, 
	help='text file with intron sequences')
parser.add_argument('--weighted', action='store_true', 
	help='used weighted introns')
parser.add_argument('--pwm_size', required=False, type=int, default=30, 
	help='length of pwm %(default)i')
	
args = parser.parse_args()

# input intron seqs need to include bases upstream of splice sites
# if you want the pwm to include those bases
if args.weighted:
	
	donor_sides = []
	acceptor_sides = []
	with open(args.introns, 'rt') as fp:
		for line in fp:
			line = line.rstrip()
			line = line.split(',')
			seq = line[0]
			weight = line[1]
			donor_side = seq[:args.pwm_size]
			acceptor_side = seq[-args.pwm_size:]
			donor_sides.append((donor_side, float(weight)))
			acceptor_sides.append((acceptor_side, float(weight)))
			
'''
if not args.weighted:
	
	donor_sides = []
	acceptor_sides = []
	with open(args.introns, 'rt') as fp:
		for line in fp.readlines():
			seq = line.rstrip()
			donor_side = seq[:args.pwm_size]
			acceptor_side = seq[-args.pwm_size:]			
			donor_sides.append(donor_side)
			acceptor_sides.append(acceptor_side)	
'''

# includes addition weights for each sequence, normalized to region
def build_weighted_pwm(seqs):
	
	pwm_size = len(seqs[0][0])
	
	lengths = []
	for seq in seqs:
		lengths.append(len(seq[0]))
		
	# make sure all sequences are the same length
	assert len(set(lengths)) == 1, 'bad sequence length'
	
	counts = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} 
				for x in range(pwm_size)]
	
	for seq in seqs:
		for i, nt in enumerate(seq[0]):
			# add weighted counts
			counts[i][nt] += seq[1]

	ppm = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for x in range(pwm_size)]
	
	for i, site in enumerate(counts):
		for nt in site:
			ppm[i][nt] = site[nt]/len(seqs)
			
	pwm = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for x in range(pwm_size)]
			
	for i, site in enumerate(ppm):
		uncertainty = 0
		for item in site.items():
			if item[1] == 0:
				uncertainty += -1e-100
			else:
				# Shannon entropy = uncertainty
				uncertainty += item[1] * math.log2(item[1])
		uncertainty = -uncertainty
		info_content = 2 - uncertainty
		for nt in site:
			pwm[i][nt] = site[nt] * info_content

	return pwm	

dpwm = build_weighted_pwm(donor_sides)

with open('weighted_donor_side_pwm.csv', 'w') as csvfile:
	dwriter = csv.writer(csvfile)
	for dsite in dpwm:
		row = [dsite[x] for x in dsite]
		dwriter.writerow(row)
	
'''
dpwm = aa.build_pwm(donor_sides, args.pwm_size)
apwm = aa.build_pwm(acceptor_sides, args.pwm_size)

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
'''
	
	
	
	
	
	
	
	
	
	
