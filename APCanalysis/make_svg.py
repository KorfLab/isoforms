import argparse
import json

parser = argparse.ArgumentParser(
	description='create svg file from APC isoforms')
parser.add_argument('apc_gff')
parser.add_argument('wb_gff')
parser.add_argument('--limit', type=int, required=False, default=10)
parser.add_argument('--out_name', type=str, required=False, 
					default='isoforms.svg')
parser.add_argument('--width', type=int, required=False, default=1000)
parser.add_argument('--height', type=int, required=False, default=1000)
parser.add_argument('--gcc', type=str, required=False, nargs=2, 
					help='Input json files with GC content data; must '
					'be in order, with the smallgenes/rnaseq data '
					'first and apc data second. '
					'Example: --gcc gcc_res/rnaseq_gcc/ce.1.1.gcc.json '
					'gcc_res/APC_base_gcc/ce.1.1.gcc.json')

arg = parser.parse_args()

isoforms = {}
with open(arg.apc_gff, 'r') as agfp:
	for line in agfp:
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
with open(arg.wb_gff, 'r') as wgfp:
	for line in wgfp:
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
	
gcy_offset = 0

# read in json data
if arg.gcc:
	
	gcy_offset = 20
	
	# smallgenes/rnaseq 
	with open(arg.gcc[0], 'r') as gc0fp:
		rnaseq_gc = json.load(gc0fp)
	rna_gc_cds = {}
	rna_gc_introns = {}
	for iso in rnaseq_gc.items():
		for ft in iso[1].items():
			
			if ft[0] == 'cds':
				for k_v in ft[1].items():
					k = k_v[0]
					v = k_v[1]
					rna_gc_cds[k] = v
					
			if ft[0] == 'intron':
				if not ft[1]: continue
				for k_v in ft[1].items():
					k = k_v[0]
					v = k_v[1]
					rna_gc_introns[k] = v
	
	# apc results
	with open(arg.gcc[1], 'r') as gc1fp:
		apc_gc = json.load(gc1fp)
		
	apc_gc_cds = {}
	apc_gc_introns = {}
	for iso in apc_gc.items():
		for ft in iso[1].items():
			for k_v in ft[1].items():
				k = k_v[0]
				v = k_v[1]
				if ft[0] == 'cds':
					apc_gc_cds[k] = v
				if ft[0] == 'intron':
					apc_gc_introns[k] = v

with open(arg.out_name, 'w') as onfp:
	onfp.write(f'<svg width="{arg.width}" height="{arg.height}">\n')
	
	y = 100
	x_offset = 0
	# draw WormBase gene
	for cdss in sorted_cdss:
		
		# gc value text overlaps if cds is too short
		# attempt to build a single string for all gc value data
		wb_gc_vals = {}
		
		# condition if there are no introns
		if len(cdss) == 1:
			height = 20
			width = cds[1] - cds[0] + 1
			rect = draw_rect(width, height, cds[0]+x_offset, y, 'blue')
			onfp.write(rect)
			
			if arg.gcc:
				cds_mid = int(round((cds[0] + cds[1])/2, 0))
				cds_key = ','.join(map(str, map(lambda x: x-1, cds)))
				gc_val = rna_gc_cds[cds_key]
				if 'cds' not in wb_gc_vals:
					wb_gc_vals['cds'] = [gc_val]
				else:
					wb_gc_vals['cds'].append(gc_val)
				gc_text = draw_text(gc_val, cds_mid-10, y-5)
				onfp.write(gc_text)

			y += 30 + gcy_offset
			continue
			
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
			onfp.write(rect)
			
			if arg.gcc:
				cds_mid = int(round((cds[0] + cds[1])/2, 0))
				cds_key = ','.join(map(str, map(lambda x: x-1, cds)))
				gc_val = rna_gc_cds[cds_key]
				if 'cds' not in wb_gc_vals:
					wb_gc_vals['cds'] = [gc_val]
				else:
					wb_gc_vals['cds'].append(gc_val)
				gc_text = draw_text(gc_val, cds_mid-10, y-5)
				onfp.write(gc_text)
				
		for int_c in intron_coors:
			height = 6
			width = int_c[1] - int_c[0] + 1 
			rect = draw_rect(width, height, int_c[0]+x_offset, y+7, 'black')
			onfp.write(rect)
			
			if arg.gcc:
				int_mid = int(round((int_c[0] + int_c[1])/2, 0))
				intron_key = ','.join(map(str, map(lambda x: x-1, int_c)))
				gc_val = rna_gc_introns[intron_key]
				if 'intron' not in wb_gc_vals:
					wb_gc_vals['intron'] = [gc_val]
				else:
					wb_gc_vals['intron'].append(gc_val)
				gc_text = draw_text(gc_val, int_mid-10, y+32)
				onfp.write(gc_text)
		
		int_string = '|'.join([f'{x[0]},{x[1]}' for x in intron_coors])
		text = draw_text(int_string, x_offset+cdss[-1][1]+10, y+14)
		onfp.write(text)
		
		y += 30 + gcy_offset
	
	# draw RNA-seq introns
	start_pt = min(rna_introns.items(), key=lambda index : index[0][0])
	for intron in rna_introns.items():
		score = intron[1]
		height = 6
		width = intron[0][1] - intron[0][0] + 1
		rect = draw_rect(width, height, intron[0][0]+x_offset, y, 'green')
		text1 = draw_text(score, start_pt[0][0]-50+x_offset, y+7)
		text2 = draw_text(f'{intron[0][0]},{intron[0][1]}', 
							intron[0][1]+10+x_offset, y+7)
		onfp.write(rect)
		onfp.write(text1)
		onfp.write(text2)
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
				onfp.write(rect)
				
				if arg.gcc:
					exon_mid = int(round((exin[1] + exin[2])/2, 0))
					exon_key = ','.join(map(str, [exin[1]-1, exin[2]-1]))
					gc_val = apc_gc_cds[exon_key]
					gc_text = draw_text(gc_val, exon_mid-10, y-5)
					onfp.write(gc_text)
				
			if exin[0] == 'intron':
				int_def.append([exin[1], exin[2]])
				prob = exin[3]
				height = 6
				width = exin[2] - exin[1] + 1
				rect = draw_rect(width, height, exin[1]+x_offset, y+7, 'black')
				onfp.write(rect)
				
				if arg.gcc:
					i_mid = int(round((exin[1] + exin[2])/2, 0))
					i_key = ','.join(map(str, [exin[1]-1, exin[2]-1]))
					gc_val = apc_gc_introns[i_key]
					gc_text = draw_text(gc_val, i_mid-10, y-5)
					onfp.write(gc_text)
				
		int_text = '|'.join([f'{x[0]},{x[1]}' for x in int_def])
		text1 = draw_text(prob, iso[0][1]+x_offset-55, y+15)
		text2 = draw_text(int_text, iso[-1][2]+x_offset+10, y+15)
		onfp.write(text1)
		onfp.write(text2)
		y += 30 + gcy_offset
	onfp.write(f'</svg>')
	
	

# ce.4.1 has 2 canonical isoforms
















