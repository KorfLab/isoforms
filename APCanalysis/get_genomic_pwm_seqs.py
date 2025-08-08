import argparse
import apc_analysis as aa
import csv

parser = argparse.ArgumentParser(description='gathers intron sequences '
	'from genome annotation for PWM training')
parser.add_argument('fa', help='fasta file with genomic sequences')
parser.add_argument('gff', help='gff file with genome annotation')
parser.add_argument('--c', action='store_true', help='gather only '
	'canonical GT/AG splice sites')
parser.add_argument('--flank', required=False, type=int, default=0,
	metavar='<integer>', help='get extra bases on each side of intron')
parser.add_argument('--pwm_size', required=False, type=int, default=30, 
	help='length of pwm %(default)i')

args = parser.parse_args()

# get only lines with information of interest
'''
zless c_elegans.PRJNA13758.WS282.annotations.gff3.gz | awk -F'\t' '($2 == "WormBase" && $3 == "gene") || ($2 == "RNASeq_splice" && $3 == "intron")'
'''

# get regions based on WormBase genes
# this does not get WormBase annotated introns
genes = {}
introns = {}
with open(args.gff, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		line = line.split('\t')
		if line[1] == 'WormBase' and line[2] == 'gene':
			
			# skip pseudogenes
			pseudo = False
			for info in line[8].split(';'):
				if info == 'biotype=pseudogene':
					pseudo = True
			if pseudo == True: continue
			
			wbgene = line[8].split(';')[0].split(':')[1]
			region = (wbgene, int(line[3]), int(line[4]), line[6])
			if line[0] not in genes:
				genes[line[0]] = [region]
			else:
				genes[line[0]].append(region)
		if line[1] == 'RNASeq_splice' and line[2] == 'intron':
			intron = (line[3], line[4], line[5], line[6])
			if line[0] not in introns:
				introns[line[0]] =[intron]
			else:	
				introns[line[0]].append(intron)			
				
# gather overlapping genes for removal
overlapping_genes = {}
for item in genes.items():
	overlapping_genes[item[0]] = []
	for i, g1 in enumerate(item[1]):
		for j, g2 in enumerate(item[1]):
			if i == j: continue
			if aa.overlap(g1, g2):
				if g1 not in overlapping_genes[item[0]]:
					overlapping_genes[item[0]].append(g1)
				if g2 not in overlapping_genes[item[0]]:
					overlapping_genes[item[0]].append(g2)	

# remove overlapping genes	
for item in overlapping_genes.items():
	for gene in item[1]:
		genes[item[0]].remove(gene)

# get chromosomal sequences
sequences = {}
current_chrom = None
with open(args.fa, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		if line.startswith('>'):
			chrom = line.split('>')[1].split(' ')[0]
			current_chrom = chrom
			sequences[chrom] = ''
		else:
			sequences[current_chrom] += line

# gather introns that overlap with gene regions
assigned_introns = {}
for item in genes.items():
	for gene_region in item[1]:
		region_info = (item[0],) + gene_region
		if item[0] in introns:
			for intron in introns[item[0]]:
				int_beg = int(intron[0])
				int_end = int(intron[1])
				reg_beg = int(gene_region[1])
				reg_end = int(gene_region[2])
				if intron[3] == gene_region[3]:
					if ((int_beg >= reg_beg and int_beg <= reg_end) or
						(int_end >= reg_beg and int_end <= reg_end) or
						(reg_beg >= int_beg and reg_beg <= int_end) or
						(reg_end >= int_beg and reg_end <= int_end)):
							if region_info not in assigned_introns:
								assigned_introns[region_info] = []
								assigned_introns[region_info].append(intron)
							else:
								assigned_introns[region_info].append(intron)
		
F = args.flank

# get intron seqences with weights
weighted_introns = {}
for item in assigned_introns.items():
	total_score = sum([float(x[2]) for x in item[1]])
	weighted_introns[item[0]] = []
	for intron in item[1]:
		int_beg = int(intron[0])
		int_end = int(intron[1])
		int_seq = sequences[item[0][0]][int_beg-1-F:int_end+F]
		if intron[3] == '-':
			int_seq = aa.revcomp(int_seq)
		if args.c:
			if int_seq[F:F+2] == 'GT' and int_seq[-F-2:-F] == 'AG':
				weighted_introns[item[0]].append([int_seq, 
								float(intron[2])/total_score])
		else:
			weighted_introns[item[0]].append([int_seq, 
							float(intron[2])/total_score])

'''
with open('weighted_genomic_introns.txt', 'wt') as fp:
	for region in weighted_introns:
		for intron in weighted_introns[region]:
			fp.write(f'{intron[0]},{intron[1]}\n')
'''

introns = []
for item in weighted_introns.items():
	for intron in item[1]:
		introns.append(intron)
		
w_dons = []
w_accs = []
dons = []
accs = []
for intron in introns:
	donor_side = intron[0][:args.pwm_size]
	acceptor_side = intron[0][-args.pwm_size:]
	w_dons.append([donor_side, intron[1]])
	w_accs.append([acceptor_side, intron[1]])
	dons.append(donor_side)
	accs.append(acceptor_side)
	
print(w_dons)
	
dpwm = aa.build_pwm(dons, args.pwm_size)
apwm = aa.build_pwm(accs, args.pwm_size)

wdpwm = aa.build_weighted_pwm(w_dons, args.pwm_size)
wapwm = aa.build_weighted_pwm(w_accs, args.pwm_size)

if args.c:
	with open('canon_weighted_genomic_donor_pwm.csv', 'w') as csvfile:
		writer = csv.writer(csvfile)
		for site in wdpwm:
			row = [dsite[x] for x in dsite]
			writer.writerow(row)

	with open('canon_weighted_genomic_acceptor_pwm.csv', 'w') as csvfile:
		writer = csv.writer(csvfile)
		for asite in wapwm:
			row = [asite[x] for x in asite]
			dwriter.writerow(row)
			
	with open('canon_genomic_donor_pwm.csv', 'w') as csvfile:
		writer = csv.writer(csvfile)
		for dsite in dpwm:
			row = [dsite[x] for x in dsite]
			dwriter.writerow(row)

	with open('canon_genomic_acceptor_pwm.csv', 'w') as csvfile:
		writer = csv.writer(csvfile)
		for asite in apwm:
			row = [asite[x] for x in asite]
			dwriter.writerow(row)

else:
	with open('weighted_genomic_donor_pwm.csv', 'w') as csvfile:
		writer = csv.writer(csvfile)
		for dsite in wdpwm:
			row = [dsite[x] for x in dsite]
			dwriter.writerow(row)

	with open('weighted_genomic_acceptor_pwm.csv', 'w') as csvfile:
		writer = csv.writer(csvfile)
		for asite in wapwm:
			row = [asite[x] for x in asite]
			dwriter.writerow(row)
			
	with open('genomic_donor_pwm.csv', 'w') as csvfile:
		writer = csv.writer(csvfile)
		for dsite in dpwm:
			row = [dsite[x] for x in dsite]
			dwriter.writerow(row)

	with open('genomic_acceptor_pwm.csv', 'w') as csvfile:
		writer = csv.writer(csvfile)
		for asite in apwm:
			row = [asite[x] for x in asite]
			dwriter.writerow(row)



# notes
# need to search for either WormBase or RNAseq_splice
# this just gets either, bad
# coordinate could appear twice from different transcripts
# then you accidently double count
# the counts on RNAseq_splice are associated with the PROMOTER
# not the splice site
# in region you are looking at, what is the one count?
# normalize the counts to region you are looking at
# small genes are already regionized
# highest count gets one count, rest gets partial count
# weight from 0-1
# compare all different methods in smallgenes set instead of 
# whole genome. much easier


