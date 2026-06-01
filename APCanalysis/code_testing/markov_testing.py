seqs = [
	'GGACGTCATG',
	'ACTGCGATTG',
	'CGTAGCTGCG',
	'CTAGCTGGTC',
]

order = 3
beg = 0
end = 0

count = {}
for seq in seqs:
	for i in range(beg+order, len(seq) - end):
		ctx = seq[i-order:i]
		nt = seq[i]
		if ctx not in count: 
			count[ctx] = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
		count[ctx][nt] += 1
		
mm = {}
for kmer in count:
	mm[kmer] = {}
	total = 0
	for nt in count[kmer]: total += count[kmer][nt]
	for nt in count[kmer]: mm[kmer][nt] = count[kmer][nt] / total


for item in count.items():
	print(item)
	
for item in mm.items():
	print(item)
