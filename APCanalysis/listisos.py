import argparse
import os
import isoform
import numpy as np

parser = argparse.ArgumentParser(
	description='create list of ranked APC genes')
parser.add_argument('APC_gffs')
parser.add_argument('WB_gffs')

args = parser.parse_args()

'''
fyle = None
for fname in os.listdir(args.APC_gffs):
	with open(f'{args.APC_gffs}{fname}') as fp:
		fyle = fp
		for line in fp.readlines():
			line = line.rstrip()
			line = line.split('\t')
			#if line[2] == 'intron':
				#print(line)
		break

	
print(fyle.name)
ints = isoform.get_introns(fyle.name)
	


print(ints)

print('#####')

introns = {}
icounts = {}
totsc = 0
with open(fyle.name) as fp:
	for line in fp.readlines():
		line = line.rstrip()
		line = line.split('\t')
		if len(line) == 9:
			if line[2] == 'intron':
				intron = (line[3], line[4])
				prob = float(line[5])
				if intron not in introns:
					introns[intron] = prob
					totsc += prob
					icounts[intron] = 1
				else:
					introns[intron] += prob
					icounts[intron] += 1
					totsc += prob
		
total = 0
for p, c in zip(introns.items(), icounts.items()):
	print(p, c)
	total += c[1]

print(total)
print(totsc)

def intron_hist(gff):
	
	introns = {}
	icounts = {}
	totsc = 0
	with open(gff) as fp:
		for line in fp.readlines():
			line = line.rstrip()
			line = line.split('\t')
			if len(line) < 8: continue
			if line[2] == 'intron':
				intron = (line[3], line[4])
				prob = float(line[5])
				if intron not in introns:
					introns[intron] = prob
					totsc += prob
					icounts[intron] = 1
				else:
					introns[intron] += prob
					icounts[intron] += 1
					totsc += prob
						
	print(totsc, 'c')
	

for fname in os.listdir(args.APC_gffs):
	intron_hist(f'{args.APC_gffs}{fname}')
	
	
	
i2 = isoform.get_introns('smallgenes/ce.1.182.gff3')

	
print(i2)


for fname in os.listdir(args.WB_gffs):
	with open(f'{args.WB_gffs}{fname}') as fp:
		for line in fp.readlines():
			line = line.rstrip()
			line = line.split('\t')
			if len(line) < 8: continue
			if line[6] == '-':
				print(fname, line[5])
				
i3 = isoform.get_introns('smallgenes/ce.1.84.gff3')

print(i3)
'''

i4 = isoform.get_introns('smallgenes/ce.1.84.gff3')
i5 = isoform.get_introns('APCisos.base/ce.1.84.APC.base.gff')

hist1 = [0.06, 0.3, 0.5, 0.05]
hist2 = [0.08, 0.2, 0.6, 0.04]

min_sum = 0
for p1, p2 in zip(hist1, hist2):
	min_sum += min(p1, p2)
	
dist = min_sum/sum(hist2)

print(dist, min_sum, 'intersection')

def intersection(hist1, hist2):
	minima = np.minimum(hist1, hist2)
	inters = np.true_divide(np.sum(minima), np.sum(hist2))
	
	print(minima)
	print(inters)
	print(np.sum(minima))
	print(np.sum(hist2))
	
intersection(hist1, hist2)

d = 0
for pi, qi in zip(hist1, hist2):
	d += abs(pi - qi)

print(d, d/2)

mu_1 = -4
mu_2 = 4
data_1 = np.random.normal(mu_1, 2.0, 1000)
data_2 = np.random.normal(mu_2, 2.0, 1000)
hist_1, _ = np.histogram(data_1, bins=100, range=[-15, 15])
hist_2, _ = np.histogram(data_2, bins=100, range=[-15, 15])


def return_intersection(hist_1, hist_2):
    minima = np.minimum(hist_1, hist_2)
    intersection = np.true_divide(np.sum(minima), np.sum(hist_2))
    return intersection
			
i = return_intersection(hist_1, hist_1)

print(i)

print('#######')

h1 = [0.7, 0.2, 0.03, 0.025, 0.02, 0.015, 0.01]
h2 = [0.8, 0.1, 0.1, 0, 0, 0, 0]
h3 = [0.7, 0.2, 0.03]
h4 = [0.8, 0.1, 0.1]

min_sum = 0
for p1, p2 in zip(h1, h2):
	min_sum += min(p1, p2)
	
dist = min_sum/sum(h2)
print(dist)

min_sum = 0
for p1, p2 in zip(h3, h4):
	min_sum += min(p1, p2)
	
dist = min_sum/sum(h4)
print(dist)















	
