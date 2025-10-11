import argparse

parser = argparse.ArgumentParser(
	description='create svg file from APC isoforms')
parser.add_argument('apc_gff')
parser.add_argument('wb_gff')
parser.add_argument('--limit', type=int, required=False, default=10)

arg = parser.parse_args()

isoforms = {}
with open(arg.apc_gff, 'r') as fp:
	for line in fp:
		line = line.rstrip()
		line = line.split('\t')
		if len(line) != 9: continue
		if line[2] == 'mRNA' or line[2] == 'gene': continue
		feature = [line[2], int(line[3]), int(line[4])]
		iso_parent = line[8].split('=')[1]
		if iso_parent not in isoforms:
			isoforms[iso_parent] = [feature]
		else:
			isoforms[iso_parent].append(feature)		

sorted_isos = []
for item in isoforms.items():
	sorted_iso = sorted(item[1], key=lambda index : index[2])
	sorted_isos.append(sorted_iso)
	
def draw_rect(width, height, x, y, fill):
	
	rect = f'<rect width="{width}" height="{height}" x="{x}" y="{y}" fill="{fill}" />\n'
	
	return rect
	
with open('isoforms.svg', 'w') as fp:
	fp.write(f'<svg width="1000" height="1000">\n')
	y = 0
	for iso in sorted_isos[:arg.limit]:
		print(iso)
		for exin in iso:
			if exin[0] == 'exon':
				height = 20
				width = exin[2] - exin[1] + 1
				rect = draw_rect(width, height, exin[1], y, 'blue')
				fp.write(rect)
			if exin[0] == 'intron':
				height = 6
				width = exin[2] - exin[1] + 1
				rect = draw_rect(width, height, exin[1], y+7, 'pink')
				fp.write(rect)
		y += 30
	fp.write(f'</svg>')
	
	

	

'''
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
'''
