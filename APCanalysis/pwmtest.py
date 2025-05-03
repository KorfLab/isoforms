import sys
import csv

introns = sys.argv[1]
genome = sys.argv[2]

# potential donors/acceptors
pds = []
pas = []
nt_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
total = 0
with open(introns, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		pd = line[:25]
		if pd[5:7] == 'GT':
			pds.append(pd)
		pa = line[-25:]
		if pa[-7:].startswith('AG'):
			pas.append(pa)
		for nt in line:
			nt_counts[nt] += 1
			total += 1	

# nt frequencies
f = [nt_counts[x]/total for x in nt_counts]
nseqs = len(pds)
# nt background counts
bg = [nseqs * f[0], nseqs * f[1], nseqs * f[2], nseqs * f[3]]
print(f)
print(sum(f))

# get nt freq for entire genome
ntc = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
nt_total = 0
with open(genome, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'): continue
		for nt in line:
			ntc[nt] += 1
			nt_total += 1
			
fg = [ntc[x]/nt_total for x in nt_counts]
print(fg)
print(sum(fg))

# instead just make a table that can be used in R	
'''
def chi_test(seqs, background):
	
	count = []
	total = 0
	for seq in seqs:
		print(seq, '#')
		
chi_test(pds, bg)
'''
		
# get donor/acceptor nt counts

def counter(seqs):
	
	counts = []
	total = 0
	for seq in seqs:
		for i, nt in enumerate(seq):
			if len(counts) <= i:
				counts.append({'A': 0, 'C': 0, 'G': 0, 'T': 0})
			counts[i][nt] += 1
			total += 1
			
	return counts
		
dcounts = counter(pds)
acounts = counter(pas)

with open('donor_counts.csv', 'w') as csvfile:
	iwriter = csv.writer(csvfile)
	iwriter.writerow([x for x in fg])
	for d in dcounts:
		iwriter.writerow([d[x] for x in d])
	

with open('acceptor_counts.csv', 'w') as csvfile:
	iwriter = csv.writer(csvfile)
	iwriter.writerow([x for x in fg])
	for a in acounts:
		iwriter.writerow([a[x] for x in d])

# just do this in R
# include some bases from before the GT and after the AG
'''
chis = []
for site in dcounts:
	chi = 0
	for i, nt in enumerate(site):
		# even if GT is 100%, background for G and T is different
		c = ((site[nt] - bg[i]) ** 2) / bg[i]
		chi += c
	chis.append(chi)
	
for x in chis:
	print(x)
'''


'''
seqs = [
	'GTAAG',
	'GTACG',
	'GTCCG',
	'GTAAG',
	'GTCAC'
]

count = []
for seq in seqs:
	for i, nt in enumerate(seq):
		if len(count) <= i:
			count.append({'A': 0, 'C': 0, 'G': 0, 'T': 0})
		count[i][nt] += 1
	
print(count)

# chi square test
counts = {'A': 83, 'C': 30, 'G': 49, 'T': 84}
freqs = {}
total = 0
for n in counts:
	total += counts[n]
	
print(total)

for n in counts:
	if n not in freqs:
		freqs[n] = counts[n]/total
		
print(freqs)
		
bg = {'A': 0.3279, 'C': 0.1915, 'G': 0.2043, 'T': 0.2763}
ecounts = {'A': 0.3279 * 246, 'C': 0.1915 * 246, 'G': 0.2043 * 246, 'T': 0.2763 * 246}
print(ecounts)
tot = 0
for n in counts:
	tot += ((counts[n]-ecounts[n])**2)/(ecounts[n])

print(tot)

# chi square 2
counts = {'A': 103, 'C': 44, 'G': 46, 'T': 53}
ecounts = {'A': 0.3279 * 246, 'C': 0.1915 * 246, 'G': 0.2043 * 246, 'T': 0.2763 * 246}
tot = 0
for n in counts:
	tot += ((counts[n]-ecounts[n])**2)/(ecounts[n])

print(tot)

# chi square 2
counts = {'A': 121, 'C': 36, 'G': 38, 'T': 51}
ecounts = {'A': 0.3279 * 246, 'C': 0.1915 * 246, 'G': 0.2043 * 246, 'T': 0.2763 * 246}
tot = 0
for n in counts:
	tot += ((counts[n]-ecounts[n])**2)/(ecounts[n])

print(tot)

# alpha for rejection
a = 0.003012705

# site 1 P is 0.0177
# 3 degrees of freedom? 4 - 1 = 3
'''



















