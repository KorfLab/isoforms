seqs = [
	'GGACGTCATG',
	'ACTGCGATTG',
	'CGTAGCTGCG',
	'CTAGCTGGTC',
]

order = 3
beg = 0
end = 0

for seq in seqs:
	
	# + 1 after end is not in isoform.py
	# adding this includes the last kmer
	# maybe is a bug in isoform.py?
	
	for i in range(beg+order, len(seq) - end + 1):
		ctx = seq[i-order:i]
		print(ctx, i)
	break
