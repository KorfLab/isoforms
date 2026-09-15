import argparse

parser = argparse.ArgumentParser(description='find canonical intron '
	'in APC gff')
parser.add_argument('gff')
parser.add_argument('--int_coors', help='example: 188,244:406,453')

args = parser.parse_args()

iso_count = 0
isos = {}
with open(args.gff, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		if line.startswith('#'): continue
		if line == '': continue
		line = line.split('\t')
		if line[2] == 'mRNA':
			iso_count += 1
		if line[2] == 'intron':
			if iso_count not in isos:
				isos[iso_count] = [[int(line[3]), int(line[4])]]
			else:
				isos[iso_count].append([int(line[3]), int(line[4])])

	

int_coors = args.int_coors.split(':')
canon_int = []
for ic in int_coors:
	intc = [int(x) for x in ic.split(',')]
	canon_int.append(intc)
	


for iso in isos:
	if isos[iso] == canon_int:
		print('exact match found')
		print('iso #:', iso)
		print('introns:', isos[iso])

print('#########')

for iso in canon_int:
	print(iso)
	for i in isos:
		if iso in isos[i]:
			print(i, isos[i])
			break
			
print(len(isos))
	












