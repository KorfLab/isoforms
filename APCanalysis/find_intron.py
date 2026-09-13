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
		print(line)
		if line[2] == 'mRNA':
			iso_count += 1
		if line[2] == 'intron':
			if iso_count not in isos:
				isos[iso_count] = [(line[3], line[4])]
			else:
				isos[iso_count].append((line[3], line[4]))

print(isos[2])			
		
		
