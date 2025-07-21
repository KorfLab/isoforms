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
		
	def _annotated_splice_sites(self):
		
		with open(self.gff, 'rt') as fp:
			for line in fp:
				line = line.rstrip()
				line = line.split('\t')
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
		
	
				

			

				
				
