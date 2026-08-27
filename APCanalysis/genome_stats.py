import argparse

parser = argparse.ArgumentParser(description='get genomic intron and '
	'exon length stats')
parser.add_argument('genome')
parser.add_argument('annotation')

args = parser.parse_args()

exls = {}
inls = {}
with open(args.annotation, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		line = line.split('\t')
		if line[1] == 'WormBase':
			if line[2] == 'exon':
				exl = int(line[4]) - int(line[3]) + 1
				if exl not in exls:
					exls[exl] = 1
				else:
					exls[exl] += 1
			if line[2] == 'intron':
				inl = int(line[4]) - int(line[3]) + 1
				if inl not in inls:
					inls[inl] = 1
				else:
					inls[inl] += 1
					
print(exls)
print(inls)

print(dict(sorted(exls.items())))

print(f'exin,length,count')
for ln in sorted(exls.items()):
	print(f'exon,{ln[0]},{ln[1]}')
	
for ln in sorted(inls.items()):
	print(f'intron,{ln[0]},{ln[1]}')









