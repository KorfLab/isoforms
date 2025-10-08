import argparse

parser = argparse.ArgumentParser(
	description='create svg file from APC isoforms')
parser.add_argument('apc_gff')
parser.add_argument('wb_gff')

arg = parser.parse_args()


c = 0 
isoforms = {}
with open(arg.apc_gff, 'r') as fp:
	for line in fp:
		line = line.rstrip()
		if line.startswith('#'): continue
		line = line.split('\t')
		if len(line) != 9: continue
		if line[2] == 'mRNA': continue
		iso_parent = line[8].split('=')[1]
		print(iso_parent)
		print(line)
		isoforms[iso_parent] = []
		break

		c += 1
		if c == 20: break
		
print('#########')		

rna_introns = {}
with open(arg.wb_gff, 'r') as fp:
	for line in fp:
		line = line.rstrip()
		line = line.split('\t')
		if line[1] == 'RNASeq_splice' and line[2] == 'intron':
			intron = (line[3], line[4])
			score = line[5]
			rna_introns[intron] = score
			
print(rna_introns)

with open('isoforms.svg', 'w') as fp:
	fp.write(f'<svg width="320" height="400">\n')
	fp.write(f'<rect width="300" height="100" x="10" '
			'style="fill:rgb(0,0,255);stroke-width:3;stroke:red" />\n')
	fp.write(f'<rect width="300" height="100" x="110" y="220" '
			'style="fill:rgb(0,0,255);stroke-width:3;stroke:red" />\n')
	fp.write(f'</svg>\n')
