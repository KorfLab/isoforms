# example markov chain for dissertation

import argparse

parser = argparse.ArgumentParser(description='generate example Markov '
	'chain data for dissertation/paper')
parser.add_argument('fasta')
parser.add_argument('--seq_len', type=int, required=False, default=20)
parser.add_argument('--order', type=int, required=False, default=1)

args = parser.parse_args()

seq = []
with open(args.fasta) as fp:
	n_count = 0
	for line in fp:
		line = line.rstrip()
		if line.startswith('>'): continue
		for n in line:
			if n_count < args.seq_len:
				seq.append(n)
			n_count += 1
			
print(seq, len(seq))

seq = ''.join(seq)
order = args.order

# A, C, G, T	
counts = [[0, 0, 0, 0] for x in range(4**order)]
n_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
for i in range(order, len(seq)):
	
	prev = seq[i-order:i]
	prev_idx = [n_idx[x] for x in prev]
	
	curr = seq[i]
	curr_idx = n_idx[curr]
	
	row_idx = 0
	for j in range(len(prev_idx)):
		row_idx += prev_idx[j]*(4**(order-j-1))
	
	counts[row_idx][curr_idx] += 1
	
print(counts)

# match output
'''
count = {}
for i in range(order, len(seq)):
	ctx = seq[i-order:i]
	nt = seq[i]
	if ctx not in count: count[ctx] = {'A':0, 'C':0, 'G':0, 'T':0}
	count[ctx][nt] += 1
	
print(count)
'''
# output matches

	
	
	
	
	
	
	
	
			
