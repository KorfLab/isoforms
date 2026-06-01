pwm_size = 5

seqs1 = [
		['GTCAT', 0.97],
		['CTGAT', 0.01],
		['ACGGG', 0.005],
		['ATCTT', 0.005]
	]

counts1 = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for x in range(pwm_size)]

wt_total = 0
for seq in seqs1:
	for i, nt in enumerate(seq[0]):
		counts1[i][nt] += seq[1]	
		wt_total += seq[1]

print(counts1)

ppm1 = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for x in range(pwm_size)]

for i, site in enumerate(counts1):
	site_sum = sum(site.values())
	for nt in site:
		ppm1[i][nt] = site[nt]/site_sum
		
print(ppm1)
for site in ppm1:
	print(site)
	print(sum(site.values()))

print('##########')
		
seqs2 = ['GTCAT', 'CTGAT', 'ACGGG', 'ATCTT']

counts2 = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for x in range(pwm_size)]

for seq in seqs2:
	for i, nt in enumerate(seq):
		counts2[i][nt] += 1
		
print(counts2)

ppm2 = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for x in range(pwm_size)]

for i, site in enumerate(counts2):
	for nt in site:
		ppm2[i][nt] = site[nt]/len(seqs2)
		
print(ppm2)
