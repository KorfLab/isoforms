import argparse

parser = argparse.ArgumentParser(
	description='create svg file from APC isoforms')
parser.add_argument('apc_gff')
parser.add_argument('wb_gff')
parser.add_argument('--limit', type=int, required=False, default=10)
parser.add_argument('--out_name', type=str, required=False, 
					default='isoforms.svg')
parser.add_argument('--width', type=int, required=False, default=1000)
parser.add_argument('--height', type=int, required=False, default=1000)

arg = parser.parse_args()

isoforms = {}
with open(arg.apc_gff, 'r') as fp:
	for line in fp:
		line = line.rstrip()
		line = line.split('\t')
		if len(line) != 9: continue
		if line[2] == 'mRNA' or line[2] == 'gene': continue
		feature = [line[2], int(line[3]), int(line[4]), line[5]]
		iso_parent = line[8].split('=')[1]
		if iso_parent not in isoforms:
			isoforms[iso_parent] = [feature]
		else:
			isoforms[iso_parent].append(feature)		

sorted_isos = []
for item in isoforms.items():
	sorted_iso = sorted(item[1], key=lambda index : index[2])
	sorted_isos.append(sorted_iso)
	
rna_introns = {}
wb_cds = {}
with open(arg.wb_gff, 'r') as fp:
	for line in fp:
		line = line.rstrip()
		line = line.split('\t')
		if line[1] == 'RNASeq_splice' and line[2] == 'intron':
			intron = (int(line[3]), int(line[4]))
			score = line[5]
			rna_introns[intron] = score
		if line[1] == 'WormBase' and line[2] == 'CDS':
			cds = [int(line[3]), int(line[4])]
			parent_tx = line[8].split(';')[1].split(':')[1]
			if parent_tx not in wb_cds:
				wb_cds[parent_tx] = [cds]
			else:
				wb_cds[parent_tx].append(cds)

sorted_cdss = []
for item in wb_cds.items():
	sorted_cds = sorted(item[1], key=lambda index : index[1])
	sorted_cdss.append(sorted_cds)	
	
def draw_rect(width, height, x, y, fill):
	
	rect = f'<rect width="{width}" height="{height}" x="{x}" y="{y}" fill="{fill}" />\n'
	
	return rect
	
def draw_text(text, x, y):
	
	text = f'<text x="{x}" y="{y}" class="sm">{text}</text>\n'
	
	return text
	
with open(arg.out_name, 'w') as fp:
	fp.write(f'<svg width="{arg.width}" height="{arg.height}">\n')
	
	y = 100
	x_offset = 0
	# draw WormBase gene
	for cdss in sorted_cdss:
		
		# this may not work for every gff3 (check CDS starts)
		if cdss[0][0] < 0: 
			x_offset = abs(cdss[0][0])
			
		intron_coors = []
		# first intron
		# 2 CDS must share same parent transcript
		first_beg = cdss[0][1] + 1
		first_end = cdss[1][0] - 1
		intron_coors.append([first_beg, first_end])
	
		# inner introns
		for i, cds in enumerate(cdss[1:len(cdss)-2]):
			inner_int = [cds[1]+1, cdss[i+2][0]-1]
			intron_coors.append(inner_int)
	
		# last intron
		last_beg = cdss[-2][1] + 1
		last_end = cdss[-1][0] - 1
		if [last_beg, last_end] not in intron_coors:
			intron_coors.append([last_beg, last_end])
	
		for cds in cdss:
			height = 20
			width = cds[1] - cds[0] + 1
			rect = draw_rect(width, height, cds[0]+x_offset, y, 'blue')
			fp.write(rect)
			
		for int_c in intron_coors:
			height = 6
			width = int_c[1] - int_c[0] + 1 
			rect = draw_rect(width, height, int_c[0]+x_offset, y+7, 'black')
			fp.write(rect)
			
		int_string = '|'.join([f'{x[0]},{x[1]}' for x in intron_coors])
		text = draw_text(int_string, x_offset+cdss[-1][0]+100, y+14)
		fp.write(text)
		
		y += 30
	
	# draw RNA-seq introns
	start_pt = min(rna_introns.items(), key=lambda index : index[0][0])
	end_pt = max(rna_introns.items(), key=lambda index : index[0][0])
	for intron in rna_introns.items():
		score = intron[1]
		height = 6
		width = intron[0][1] - intron[0][0] + 1
		rect = draw_rect(width, height, intron[0][0]+x_offset, y, 'green')
		text1 = draw_text(score, start_pt[0][0]-50+x_offset, y+7)
		text2 = draw_text(f'{intron[0][0]},{intron[0][1]}', 
							end_pt[0][0]+100+x_offset, y+7)
		fp.write(rect)
		fp.write(text1)
		fp.write(text2)
		y += 20
	
	# draw APC isoforms
	y += 10
	for iso in sorted_isos[:arg.limit]:
		prob = iso[0][3]
		int_def = []
		for exin in iso:
			if exin[0] == 'exon':
				prob = exin[3]
				height = 20
				width = exin[2] - exin[1] + 1
				rect = draw_rect(width, height, exin[1]+x_offset, y, 'blue')
				fp.write(rect)
			if exin[0] == 'intron':
				int_def.append([exin[1], exin[2]])
				prob = exin[3]
				height = 6
				width = exin[2] - exin[1] + 1
				rect = draw_rect(width, height, exin[1]+x_offset, y+7, 'black')
				fp.write(rect)
		int_text = '|'.join([f'{x[0]},{x[1]}' for x in int_def])
		text1 = draw_text(prob, iso[0][1]+x_offset-55, y+15)
		text2 = draw_text(int_text, iso[-1][2]+x_offset+10, y+15)
		fp.write(text1)
		fp.write(text2)
		y += 30
	fp.write(f'</svg>')
	
	

# ce.4.1 has 2 canonical isoforms
















