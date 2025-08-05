import argparse
import glob
import apc_analysis as aa
from apc_analysis import SpliceSites
import copy

parser = argparse.ArgumentParser(
	description='creates PWMs from smallgenes')
parser.add_argument('smallgenes', type=str, metavar='<directory>', 
	help='path to smallgenes directory')
parser.add_argument('--don_len', required=False, type=int, default=5, 
	metavar='<integer>', help='donor site length')
parser.add_argument('--acc_len', required=False, type=int, default=5, 
	metavar='<integer>', help='acceptor site length')
parser.add_argument('--don_left', required=False, type=int, default=0, 
	metavar='<integer>', 
	help='get bases on donor site left flank [%(default)i]')
parser.add_argument('--acc_right', required=False, type=int, default=0, 
	metavar='<integer>', 
	help='get bases on acceptor site right flank [%(default)i]')

args = parser.parse_args()	

# shorten arguments for indexing
DN = args.don_len
AN = args.acc_len
DL = args.don_left 
AR = args.acc_right 

# make sure parameters work
#print(f'DN: {DN} AN: {AN} DL: {DL} AR: {AR}')
	
if args.smallgenes.endswith('/'):
	pass 
else:
	args.smallgenes = args.smallgenes + '/'
	
annotated_splice_sites = []
rnaseq_splice_sites = {}
for file in glob.glob(f'{args.smallgenes}*.fa'):
	ff = file
	gf = ff[:-2] + 'gff3'
	seq = ''
	with open(ff, 'rt') as fp:
		for line in fp:
			line = line.rstrip()
			if line.startswith('>'): continue
			else:
				seq += line
	# there are - strand features mixed in for 13 genes
	# only get + strand features
	with open(gf, 'rt') as fp:
		splice_site_set = {}
		for line in fp:
			line = line.rstrip()
			line = line.split('\t')
			
			if line[6] == '-': continue
			
			if line[1] == 'WormBase' and line[2] == 'intron':
				donor = seq[int(line[3])-DL-1:int(line[3])+DN-1]
				acceptor = seq[int(line[4])-AN:int(line[4])+AR]
				print(donor, acceptor)
				annotated_splice_sites.append((donor, acceptor))
				
			if line[1] == 'RNASeq_splice' and line[2] == 'intron':
				donor = seq[int(line[3])-3:int(line[3])+5]
				acceptor = seq[int(line[4])-8:int(line[4])]
				splice_site_set[line[5]] = (donor, acceptor) 
				
		# assign splice sites to score
		rnaseq_splice_sites[line[0]] = splice_site_set
		
#print(annotated_splice_sites)

#print(rnaseq_splice_sites)

for file in glob.glob(f'{args.smallgenes}*.fa'):
	ff = file
	gf = ff[:-2] + 'gff3'
	wb_pwm = SpliceSites(ff, gf, DN, DL, AN, AR, source='WormBase')
	rna_pwm = SpliceSites(ff, gf, DN, DL, AN, AR, source='RNASeq')
	for d, a in wb_pwm.splice_sites():
		print(d, a)
	for d, a, s in rna_pwm.splice_sites():
		print(d, a, s)
		
#genome_pwm = SpliceSites('1pct.fa', '1pct.gff3', DN, DL, AN, AR, source='WormBase')

#for i, j in genome_pwm.splice_sites():
#	print(i, j)
'''
# each chromosome restarts coordinates
#for d, a in genome_pwm.splice_sites():
#	print(d, a)
'''

'''
with open(self.gff, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		line = line.split('\t')
		if line[6] == '-': continue
		if line[1] == 'WormBase' and line[2] == 'intron':
			print('wow')
			donor = self.seq[int(line[3])-self.DL-1:
							int(line[3])+self.DN-1]
			acceptor = self.seq[int(line[4])-self.AN:
								int(line[4])+self.AR]
			print(donor, acceptor)
'''
print('########')

# get regions based on WormBase genes
genes = {}
introns = {}
with open('1pct.gff3', 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		line = line.split('\t')
		if line[1] == 'WormBase' and line[2] == 'gene':
			wbgene = line[8].split(';')[0].split(':')[1]
			region = (wbgene, int(line[3]), int(line[4]), line[6])
			if line[0] not in genes:
				genes[line[0]] = [region]
			elif len(genes[line[0]]) < 100:
				genes[line[0]].append(region)
		if line[1] == 'RNASeq_splice' and line[2] == 'intron':
			intron = (line[3], line[4], line[5], line[6])
			if line[0] not in introns:
				introns[line[0]] =[intron]
			elif len(introns[line[0]]) < 100:
				introns[line[0]].append(intron)			
				
# gather overlapping genes
overlapping_genes = {}
for item in genes.items():
	overlapping_genes[item[0]] = []
	for i, g1 in enumerate(item[1]):
		for j, g2 in enumerate(item[1]):
			if i == j: continue
			#print(item[0], i, j, g1, g2, overlap(g1, g2))
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
with open('1pct.fa', 'rt') as fp:
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
		
# need to remove pseudogenes			
# get intron seqences with weights
weighted_introns = {}
for item in assigned_introns.items():
	total_score = sum([int(x[2]) for x in item[1]])
	weighted_introns[item[0]] = []
	for intron in item[1]:
		int_beg = int(intron[0])
		int_end = int(intron[1])
		int_seq = sequences[item[0][0]][int_beg-1-DL:int_end+AR]
		if intron[3] == '-':
			int_seq = aa.revcomp(int_seq)
		#print(int_seq[:10], int_seq[-10:], item[0][1], item[0][2], 
		#		item[0][3], intron)
		weighted_introns[item[0]].append([int_seq, 
							int(intron[2])/total_score])
	
with open('weighted_introns.txt', 'wt') as fp:
	for region in weighted_introns:
		for intron in weighted_introns[region]:
			if len(intron[0]) < 100:
				print(intron[0])
			fp.write(f'{intron[0]},{intron[1]}\n')
		#print('####')


							
							

		
		
			



	
	
	
	
	
	
	
	
	
