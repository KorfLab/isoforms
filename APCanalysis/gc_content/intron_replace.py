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

# high test 7,18
# low test 6,17

def gc_calc(seq):
	
	nts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
	for nt in seq:
		nts[nt] += 1
		
	gc_calc = (
		((nts['G'] + nts['C']) / 
		(nts['A'] + nts['C'] + nts['G'] + nts['T'])) * 100
	)	
	
	return round(gc_calc, 2)

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

# remove splice sites from low gc intron
l_rm_seq = []
skip = False
for i in range(len(l_int_seq)):
	if skip:
		skip = False
		continue
	# don't change gc calc
	if ''.join(l_int_seq[i:i+2]) == 'GT':
		l_rm_seq.append('C')
		continue
	# don't change gc calc and don't add stop codons
	if ''.join(l_int_seq[i:i+2]) == 'AG':
		l_rm_seq.append('C')
		l_rm_seq.append('T')
		skip = True
		continue
	else:
		l_rm_seq.append(l_int_seq[i])
			
#print(l_int_seq)

#print(l_rm_seq)
			
#print(gc_calc(l_int_seq), 'l_int_seq')
#print(gc_calc(l_rm_seq), 'l_rm_seq')

# find intron-based splice site coordinates in high gc intron
gt_sites = []
ag_sites = []
for i in range(len(h_int_seq)):
	if ''.join(h_int_seq[i:i+2]) == 'GT':
		gt_sites.append(i) #+h_int_coor[0])
	if ''.join(h_int_seq[i:i+2]) == 'AG':
		ag_sites.append(i) #+h_int_coor[0])

#print(gt_sites, ag_sites)		

# create replacement intron
rep_int_seq = []
len_ct = 0
for i in range(len(h_int_seq)):
	if  len_ct == len(l_int_seq):
		len_ct = 0
		rep_int_seq.append(l_rm_seq[len_ct])
	else:
		rep_int_seq.append(l_rm_seq[len_ct])
	len_ct += 1
	
#print(rep_int_seq, 'rep_int_seq')
		
#print(gt_sites, ag_sites)
		
# restore splice sites
new_int_seq = []
skip = False
for i in range(len(rep_int_seq)):
	if skip == True:
		skip = False
		continue
	if i in gt_sites:
		new_int_seq.append('G')
		new_int_seq.append('T')
		skip = True
		continue
	if i in ag_sites:
		new_int_seq.append('A')
		new_int_seq.append('G')
		skip = True
		continue
	else:
		new_int_seq.append(rep_int_seq[i])
		
#print(new_int_seq, 'new_int_seq')

#print(gc_calc(rep_int_seq), gc_calc(l_int_seq))
#print(gc_calc(new_int_seq), gc_calc(h_int_seq))

# replace high gc intron with low gc intron in the gene seq
new_seq = []
rep_idx = 0
for i in range(len(h_seq)):
	if i >= h_int_coor[0]-1 and i <= h_int_coor[1]-1:
		new_seq.append(new_int_seq[rep_idx])
		rep_idx += 1
	else:
		new_seq.append(h_seq[i])
		
#print(h_seq, 'h_seq')
#print(new_seq, 'new_seq')
	
line_len = 80
	
# print in fasta format
fa_lines = {}
n_ct = 0
ln_ct = 0
for i in range(len(new_seq)):
	if n_ct <= line_len:
		if ln_ct not in fa_lines:
			fa_lines[ln_ct] = [new_seq[i]]
		else:
			fa_lines[ln_ct].append(new_seq[i])
		n_ct += 1
	if n_ct == line_len:
		n_ct = 0
		ln_ct += 1
		
print(h_header)
for item in fa_lines.items():
	print(''.join(item[1]))
	
	
	
	
	
			
