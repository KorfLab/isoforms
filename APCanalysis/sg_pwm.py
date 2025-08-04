import argparse
import glob
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
			elif len(genes[line[0]]) < 10:
				genes[line[0]].append(region)
		if line[1] == 'RNASeq_splice' and line[2] == 'intron':
			intron = (line[3], line[4], line[5], line[6])
			if line[0] not in introns:
				introns[line[0]] =[intron]
			elif len(introns[line[0]]) < 10:
				introns[line[0]].append(intron)

i1 = ('w', 100, 200, '+')
i2 = ('w', 150, 250, '-')
i3 = ('w', 300, 400, '+')				
				
def overlap(g1, g2):
	
	if g1[1] >= g2[1] and g1[1] <= g2[2]:
		return True
	elif g2[1] >= g1[1] and g2[1] <= g1[2]:
		return True
	else:
		return False
	
	
print(overlap(i2, i3)) 
# remove overlapping genes
# should also remove things like ncRNA_gene?
'''
gene_coors = {}
for item in genes.items():
	for gene in item[1]:
		current_intron = (int(gene[1]), int(gene[2]))
		if item[0] not in gene_coors:
			gene_coors[item[0]] = [current_intron]
		else:
			for intron in gene_coors[item[0]]:
				print(item[0], intron, current_intron, overlap(intron, current_intron))
			gene_coors[item[0]].append(current_intron)
			
#print(gene_coors['X'])
'''
'''
genes_no_overlap = {}
for item in genes.items():
	print(item[0])
	for current_gene in item[1]:
		if item[0] not in genes_no_overlap:
			genes_no_overlap[item[0]] = [current_gene]
		else:
			for gene in genes_no_overlap[item[0]]:
				print(gene, current_gene, overlap(gene, current_gene))
				if overlap(gene, current_gene):
					genes_no_overlap[item[0]].remove(current_gene)
'''

overlapping_genes = {}
for item in genes.items():
	overlapping_genes[item[0]] = []
	for i, g1 in enumerate(item[1]):
		for j, g2 in enumerate(item[1]):
			if i == j: continue
			#print(item[0], i, j, g1, g2, overlap(g1, g2))
			if overlap(g1, g2):
				if g1 not in overlapping_genes[item[0]]:
					overlapping_genes[item[0]].append(g1)
				if g2 not in overlapping_genes[item[0]]:
					overlapping_genes[item[0]].append(g2)	
			
for item in overlapping_genes.items():
	for gene in item[1]:
		genes[item[0]].remove(gene)
		
for item in genes.items():
	print(item[0])
	for g in item[1]:
		print(g)
		
	
				
'''

4116...10230	11495...16837	17484...26781	22882...23600



'''

				
			
				
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
		
# there is weird stuff here, like a 21 bp 'gene'
# need to filter?
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
						
def revcomp(seq):

	comps = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

	rev = ''
	for i in range(1, len(seq) + 1):
		rev += comps[seq[-i]]

	return rev
	
								
adjusted_intron_counts = {}
for item in assigned_introns.items():
	#print(item)
	for intron in item[1]:
		#print(intron)
		int_beg = int(intron[0])
		int_end = int(intron[1])
		int_seq = sequences[item[0][0]][int_beg-1:int_end]
		if intron[3] == '-':
			int_seq = revcomp(int_seq)
		#print(int_seq)
	break
		

							
							

		
		
			



	
	
	
	
	
	
	
	
	
