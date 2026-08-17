import argparse

parser = argparse.ArgumentParser(
	description='replace high gc intron with low gc intron')
parser.add_argument('high_fasta', help='fasta with high gc content intron')
parser.add_argument('high_int', help='high gc content intron coordinates')
parser.add_argument('low_fasta', help='fasta with low gc content intron')
parser.add_argument('low_int', help='low gc content intron coordinates')

args = parser.parse_args()

# high gc ce.1.542 261,633
# low gc ce.1.150 238,292

def get_seq_info(fasta, int_coor):
	
	header = ''
	seq = []
	with open(fasta, 'rt') as fp:
		for line in fp:
			line = line.rstrip()
			if line.startswith('>'): 
				header = line
				continue
			for n in line:
				seq.append(n)
				
	coor = [int(x) for x in int_coor.split(',')]
	
	sub_seq = seq[coor[0]-1:coor[1]]
	
	return header, coor, seq, sub_seq
	
(h_header, h_int_coor,
h_seq, h_int_seq) = get_seq_info(args.high_fasta, args.high_int)

(l_header, l_int_coor,
l_seq, l_int_seq) = get_seq_info(args.low_fasta, args.low_int)

def remove_ss(int_seq):
	
	rm_seq = []
	for i in range(len(int_seq)):
		if ''.join(int_seq[i:i+2]) == 'GT':
			rm_seq.append('C')
			continue
		if ''.join(int_seq[i:i+2]) == 'AG':
			rm_seq.append('T')
			continue
				
	
'''		
low_int_seq = low_seq[low_coor[0]-1:low_coor[1]]
print(low_int_seq, len(low_int_seq))

# remove GT and AG sites
low_seq = []
for i in range(len(low_int_seq)):
	if ''.join(low_int_seq[i:i+2]) == 'GT':
		low_seq.append('C')
		continue
	if ''.join(low_int_seq[i:i+2]) == 'AG':
		low_seq.append('T')
		continue
	else:
		low_seq.append(low_int_seq[i])
		
print(low_seq, len(low_seq))
	
print('######')

high_seq = []
header = None
with open(args.fasta, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		if line.startswith('>'): 
			header = line
			continue
		for n in line:
			high_seq.append(n)
			
high_coor = [261, 633]
high_int_seq = high_seq[high_coor[0]-1:high_coor[1]]

print(high_int_seq)

# find GT and AG site coordinates
gt_sites = []
ag_sites = []
for i in range(len(high_int_seq)):
	#print(i, high_int_seq[i:i+2])
	if ''.join(high_int_seq[i:i+2]) == 'GT':
		gt_sites.append(i)
	if ''.join(high_int_seq[i:i+2]) == 'AG':
		ag_sites.append(i)

print(gt_sites, ag_sites)

rep_int_seq = []
len_ct = 0
for i in range(len(high_int_seq)):
	if  len_ct == len(low_seq):
		len_ct = 0
		rep_int_seq.append(low_seq[len_ct])
	else:
		rep_int_seq.append(low_seq[len_ct])
	len_ct += 1
	
print(rep_int_seq)

print(len(rep_int_seq), len(high_int_seq))
print('########')
# put gt and ag sites back
rep_int = []
skip = False
for i in range(len(rep_int_seq)):
	if skip == True:
		skip = False
		continue
	if i in gt_sites:
		rep_int.append('G')
		rep_int.append('T')
		skip = True
		continue
	if i in ag_sites:
		rep_int.append('A')
		rep_int.append('G')
		skip = True
		continue
	else:
		rep_int.append(rep_int_seq[i])
		


print(rep_int)
	
print(len(rep_int), len(rep_int_seq))

print('###')

def gc_calc(seq):
	
	nts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
	for nt in seq:
		nts[nt] += 1
		
	gc_calc = (
		((nts['G'] + nts['C']) / 
		(nts['A'] + nts['C'] + nts['G'] + nts['T'])) * 100
	)	
	
	return round(gc_calc, 2)
	
print(gc_calc(rep_int), gc_calc(high_int_seq), gc_calc(low_int_seq))


# create new fasta with replaced intron
replacement = []
rep_idx = 0
for i in range(len(high_seq)):
	if i >= low_coor[0]-1 and i <= low_coor[1]-1:
		print(i, rep_idx, low_coor[0]-1)
		replacement.append(rep_int[rep_idx])
		rep_idx += 1
	else:
		replacement.append(high_seq[i])
		
print('######')

print(replacement, len(replacement), len(high_seq))
print(gc_calc(replacement), gc_calc(high_seq))

# print in fasta format
fa_lines = {}
n_ct = 0
ln_ct = 0
for i in range(len(replacement)):
	if n_ct <= 80:
		if ln_ct not in fa_lines:
			fa_lines[ln_ct] = [replacement[i]]
		else:
			fa_lines[ln_ct].append(replacement[i])
		n_ct += 1
	if n_ct == 80:
		n_ct = 0
		ln_ct += 1
		
print(len(fa_lines[11]))

print(header)
for item in fa_lines.items():
	print(''.join(item[1]))
	
print('########')
print(''.join(replacement[261:633]))

print(''.join(high_seq[261:633]))
		
'''
	
	
	
	
	
			
