import argparse
import os
import locallib as ll
import isoform

parser = argparse.ArgumentParser(
	description='create list of ranked APC genes')
parser.add_argument('APC_gffs')
parser.add_argument('WB_gffs')

args = parser.parse_args()

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
			
	
