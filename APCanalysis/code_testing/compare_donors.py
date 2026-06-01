import sys

d1 = sys.argv[1]
d2 = sys.argv[2]

d1_counts = {}
d1_total = 0
with open(d1, 'rt') as fp:
	for line in fp:
		dseq1 = line.rstrip()
		if dseq1 not in d1_counts:
			d1_counts[dseq1] = 1
			d1_total += 1
		else:
			d1_counts[dseq1] += 1
			d1_total += 1
		
dlist1 = sorted(d1_counts.items(), key=lambda item: item[0])

d2_counts = {}
d2_total = 0
with open(d2, 'rt') as fp:
	for line in fp:
		dseq2 = line.rstrip()
		if dseq2 not in d2_counts:
			d2_counts[dseq2] = 1
			d2_total += 1
		else:
			d2_counts[dseq2] += 1
			d2_total += 1
			
dlist2 = sorted(d2_counts.items(), key=lambda item: item[0])

print(len(dlist1), len(dlist2))
print(d1_total, d2_total)

print(dlist1)
print('###')
print(dlist2)
