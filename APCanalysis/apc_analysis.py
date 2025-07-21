class PWM:

	def __init__(self, fasta, gff, don_len, don_left, acc_len, 
					acc_right):
		self.fasta = fasta
		self.gff = gff
		self.seq = None
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
						
		return self.seq
		
	def splice_sites(self):
		
		annotated_splice_sites = []
		with open(self.gff, 'rt') as fp:
			splice_site_set = {}
			for line in fp:
				line = line.rstrip()
				line = line.split('\t')
			
				if line[6] == '-': continue
			
				if line[1] == 'WormBase' and line[2] == 'intron':
					acceptor = self.seq[int(line[4])-self.AN:int(line[4])+self.AR]
					annotated_splice_sites.append((acceptor))
			
			return annotated_splice_sites
				
				
