class PWM:

	def __init__(self, fasta, gff):
		self.fasta = fasta
		self.gff = gff
		self.seq = None
		
	def intron_seq(self):
		
		self.seq = ''
		with open(self.fasta, 'rt') as fp:
			for line in fp:
				line = line.rstrip()
				if line.startswith('>'): continue
				else:
					self.seq += line
						
		return self.seq
		
	def get_ss(self, seq, gff):
		
		print(wow)
	
