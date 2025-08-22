import math

def overlap(g1, g2):
	
	if g1[1] >= g2[1] and g1[1] <= g2[2]:
		return True
	elif g2[1] >= g1[1] and g2[1] <= g1[2]:
		return True
	else:
		return False
		
def revcomp(seq):

	comps = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

	rev = ''
	for i in range(1, len(seq) + 1):
		rev += comps[seq[-i]]

	return rev
	
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
			if item[1] == 0:
				uncertainty += -1e-100
			else:
				uncertainty += item[1] * math.log2(item[1])
		uncertainty = -uncertainty
		info_content = 2 - uncertainty
		for nt in site:
			pwm[i][nt] = site[nt] * info_content

	return ppm	

# includes addition of weights for each sequence, normalized to region
def build_weighted_pwm(seqs, pwm_size):
	
	counts = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} 
				for x in range(pwm_size)]
	
	for seq in seqs:
		# discard sequences that are too short
		if len(seq[0]) < pwm_size: continue
		for i, nt in enumerate(seq[0]):
			# add weighted counts
			counts[i][nt] += float(seq[1])

	ppm = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for x in range(pwm_size)]
	
	for i, site in enumerate(counts):
		site_sum = sum(site.values())
		for nt in site:
			ppm[i][nt] = site[nt]/site_sum
				
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

	return ppm

def print_pwm(pwm, pwm_name):

	for row in pwm:
		site = [x[1] for x in row.items()]
		ent = sum(site)
		#print(ent, '***')
		print(f'{pwm_name},{site[0]},{site[1]},{site[2]},{site[3]}')

class SpliceSites:

	def __init__(self, fasta, gff, don_len, don_left, acc_len, 
					acc_right, source=None):

		self.fasta = fasta
		self.gff = gff
		self.seq = None
		self.source = source
		self.DN = don_len
		self.DL = don_left
		self.AN = acc_len
		self.AR = acc_right
		
	def _intron_seq(self):
		
		self.seq = ''
		with open(self.fasta, 'rt') as fp:
			for line in fp:
				line = line.rstrip()
				if line.startswith('>'): continue
				else:
					self.seq += line
					
	def _revcomp(self):

		comps = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

		rev = ''
		for i in range(1, len(self.seq) + 1):
			rev += comps[self.seq[-i]]

		return rev
	
	def _annotated_splice_sites(self):
		
		with open(self.gff, 'rt') as fp:
			for line in fp:
				line = line.rstrip()
				line = line.split('\t')
				# smallgenes coverts - to +?
				if line[6] == '-': continue
				if line[1] == 'WormBase' and line[2] == 'intron':
					donor = self.seq[int(line[3])-self.DL-1:
										int(line[3])+self.DN-1]
					acceptor = self.seq[int(line[4])-self.AN:
										int(line[4])+self.AR]
						
					yield donor, acceptor
	
	def _rnaseq_splice_sites(self):
		
		with open(self.gff, 'rt') as fp:
			for line in fp:
				line = line.rstrip()
				line = line.split('\t')
				if line[6] == '-': continue
				if line[1] == 'RNASeq_splice' and line[2] == 'intron':
					donor = self.seq[int(line[3])-self.DL-1:
										int(line[3])+self.DN-1]
					acceptor = self.seq[int(line[4])-self.AN:
										int(line[4])+self.AR]				
										
					yield donor, acceptor, line[5]
	
	def splice_sites(self):

		self._intron_seq()
		splice_sites = []
		
		if self.source == 'WormBase': 
			for don, acc in self._annotated_splice_sites():
				splice_sites.append((don, acc))
			
			return splice_sites
			
		if self.source == 'RNASeq':
			for don, acc, score in self._rnaseq_splice_sites():
				splice_sites.append((don, acc, score))
			
			return splice_sites
		
'''
class GenomicSpliceSites:

	def __init__(self, genome, annotation, flank=0, pwm_size=30, 
				canonical=False):

		self.fasta = genome
		self.gff = annotation
		self.flank = flank
		self.pwm_size = pwm_size
		self.canonical = canonical
		self.genes = {}
		self.introns = {}
		self.chr_seqs = {}
	
	# get regions based on WormBase genes
	# this does not get WormBase annotated introns
	def _gather_regions(self):

		with open(self.gff, 'rt') as fp:
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
					region = (wbgene, int(line[3]), int(line[4]), 
								line[6])
					if line[0] not in self.genes:
						self.genes[line[0]] = [region]
					else:
						self.genes[line[0]].append(region)
				if line[1] == 'RNASeq_splice' and line[2] == 'intron':
					intron = (line[3], line[4], line[5], line[6])
					if line[0] not in self.introns:
						self.introns[line[0]] =[intron]
					else:	
						self.introns[line[0]].append(intron)
	
	def _overlapping_genes(self):
		
		# gather overlapping genes
		overlapping_genes = {}
		for item in self.genes.items():
			overlapping_genes[item[0]] = []
			for i, g1 in enumerate(item[1]):
				for j, g2 in enumerate(item[1]):
					if i == j: continue
					if overlap(g1, g2):
						if g1 not in overlapping_genes[item[0]]:
							overlapping_genes[item[0]].append(g1)
						if g2 not in overlapping_genes[item[0]]:
							overlapping_genes[item[0]].append(g2)
		
		# remove overlapping genes
		for item in overlapping_genes.items():
			for gene in item[1]:
				self.genes[item[0]].remove(gene)

	def _chromosomal_sequences(self):

		# get chromosomal sequences
		current_chrom = None
		with open(self.fasta, 'rt') as fp:
			for line in fp:
				line = line.rstrip()
				if line.startswith('>'):
					chrom = line.split('>')[1].split(' ')[0]
					current_chrom = chrom
					self.chr_seqs[chrom] = ''
				else:
					self.chr_seqs[current_chrom] += line

	def res(self):

		self._gather_regions
		self._overlapping_genes
		self._chromosomal_sequences()

		self.chr_seqs

	
'''
class PWM:
	
	def __init__(self):
		print('pwm')
	
				

			

				
				
