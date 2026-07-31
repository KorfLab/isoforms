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
	
#print(counts)
#print('##########')

freqs = []
for i in range(len(counts)):
	row_sum = 0
	for ct in counts[i]:
		row_sum += ct
	fq_row = []
	for ct in counts[i]:
		if row_sum != 0:
			freq = round(ct/row_sum, 2)
		else:
			freq = 0
		fq_row.append(freq)
	if sum(fq_row) != 1:
		fq_row[3] = fq_row[3] + (1-sum(fq_row))
		fq_row[3] = f'{fq_row[3]:.2f}'
	freqs.append(fq_row)

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

# recursion from copilot
def markov_rows(order, nts=('A', 'C', 'G', 'T')):
	
	if order == 0:
		return ['']
	prev = markov_rows(order-1, nts)
	return [p + n for p in prev for n in nts]
	
# cartesian product
'''
l1 = [1, 2]
l2 = ['A', 'B']

res = [(x, y) for x in l1 for y in l2]

print(res)
'''

res = markov_rows(order)

# format to latex
'''
$$
\begin{array}{c|ccc}
     & c_1 & c_2 & c_3 \\ \hline
r_1  & 1   & 0   & 0   \\
r_2  & 0   & 1   & 0   \\
r_3  & 0   & 0   & 1
\end{array}
$$
'''

print(seq)

print('$$')
print(r'\begin{array}{c|cccc}')
print(r'    & A   & C   & G   & T   \\ \hline')
for i in range(len(res)):
	print(f'{res[i]}  {''.join([f'& {str(x)}   ' for x in counts[i]])}'
			r'\\')
print(r'\end{array}')
print('$$')

print('$$')
print(r'\begin{array}{c|cccc}')
print(r'    & A   & C   & G   & T   \\ \hline')
for i in range(len(res)):
	print(f'{res[i]}  {''.join([f'& {str(x)}   ' for x in freqs[i]])}'
			r'\\')
print(r'\end{array}')
print('$$')




'''
for i in range(len(res)):
	print(f'{res[i]},{','.join([str(x) for x in counts[i]])}')
	
	
	
for i in range(len(res)):
	print(f'{res[i]},{','.join([str(x) for x in freqs[i]])}')
'''
	

	
	
	
	
	
			
