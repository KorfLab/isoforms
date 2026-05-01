import argparse
import gzip
import os

'''
parser = argparse.ArgumentParser(description='create gffs for WBGenes '
	'of interest with in frame start/stop coordinates')
parser.add_argument('WBGenes', help='csv with list of WBGene IDs and '
	'region of interest i.e. '
	'WBGene, gene name, chr, region start, region end, '
	'gene start, gene end, sense | '
	'WBGene00003386,mod1,V,8910090,8910840,8913992,-')
parser.add_argument('annotation', help='annotation gff')
parser.add_argument('out_dir', help='name of directory to store '
			'fa and gff files')

args = parser.parse_args()

gene_info = {}
with open(args.WBGenes, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		if line.startswith('#'): continue
		line = line.split(',')
		gene_info[line[0]] = line[1:]

open_type = gzip.open if args.annotation.endswith('.gz') else open

gene_lines = {}
with open_type(args.annotation, 'rt') as fp:
	for line in fp:
		line = line.rstrip()
		line = line.split('\t')
		for info in gene_info.items():
			gid_info = (info[0], info[1][0], info[1][2], info[1][3])
			if line[0] != info[1][1]: continue
			if gid_info not in gene_lines:
				gene_lines[gid_info] = []
			if len(line) == 9:
				gff_start, gff_end = int(line[3]), int(line[4])
				subreg_start, subreg_end = int(info[1][2]), int(info[1][3])
				gen_start, gen_end = int(info[1][4]), int(info[1][5])
				
				# match WBGene
				if line[2] == 'gene':
					wbg_gff = line[8].split(';')[0].split(':')[1]
					if info[0] == wbg_gff:
						gene_lines[gid_info].append(line)
				
				# match WBGene
				if line[2] == 'mRNA':
					wbg_gff = line[8].split(';')[1].split(':')[1]
					if info[0] == wbg_gff:
						idt = line[8].split(';')[0].split(':')[1].split('.')
						if idt[1].endswith('a') or idt[1] == '1':
							gene_lines[gid_info].append(line)
						
				# get any CDS regions that overlap region of interest
				if (line[0] == info[1][1] and line[1] == 'WormBase' and
					line[2] == 'CDS'): 
					if gff_start <= subreg_end and gff_end >= subreg_start:
						gene_lines[gid_info].append(line)
				
				# get any introns only within region
				if line[0] == info[1][1] and line[1] == 'RNASeq_splice': 
					if gff_start >= subreg_start and gff_end <= subreg_end:
						
			############# filter low count introns for testing/readability
						#if int(line[5]) < 20000: continue
						
						gene_lines[gid_info].append(line)

# remove CDS regions not matching mRNA ID=Transcript
gene_lines_2 = {}
for item in gene_lines.items():
	gene_lines_2[item[0]] = []
	
	# get transcript ID
	m_tid = None
	for line in item[1]:
		if line[2] == 'mRNA':
			m_tid = line[8].split(';')[0].split(':')[1]
			break
	
	# add to new dictionary
	for line in item[1]:
		if line[2] == 'CDS':
			c_tid = line[8].split(';')[0].split(':')[1]
			# only add 'a' or '1' CDS entries
			# that match mRNA tx id
			if c_tid.split('.')[0] == m_tid.split('.')[0]:
				if c_tid.split('.')[1].endswith('a'):
					gene_lines_2[item[0]].append(line)
				if c_tid.split('.')[1] == '1':
					gene_lines_2[item[0]].append(line)
		# keep everything else
		else:
			gene_lines_2[item[0]].append(line)

'''

'''
for item in gene_lines_2.items():
	print(','.join(item[0]))
	for line in item[1]:
		print(','.join(line))
'''

	

import sys

gl2 = sys.argv[1]

gene_lines_2 = {}
with open(gl2, 'rt') as fp:
	current_lid = None
	for line in fp:
		line = line.rstrip()
		if line.startswith('WBGene'):
			line_id = tuple(line.split(','))
			current_lid = line_id
			gene_lines_2[line_id] = []
		else:
			gene_lines_2[current_lid].append(line.split(','))

		
# assume CDS regions appear in order without sorting
# all CDS should be on the same strand
# add in frame start or stop codon
# same for + and - strand
gene_lines_3 = {}
for item in gene_lines_2.items():
	print(item[0])
	print('########')
	# count number of CDS
	n_cds = 0
	for line in item[1]:
		if line[2] == 'CDS':
			n_cds += 1

	n_cds2 = 0
	total_len = 0
	for line in item[1]:
		if line[2] == 'CDS':
			new_line = line.copy()
			#print(line, '###')
			print(n_cds2, n_cds)
			# add first cds length
			if n_cds2 == 0:
				total_len += int(item[0][2]) - int(line[3]) + 1
				new_line[3] = item[0][2]
				gene_lines_3[item[0]] = [new_line]
				n_cds2 += 1
				continue
			
			# add inner CDS lengths
			if 1 <= n_cds2 < (n_cds - 1):
				total_len += int(line[4]) - int(line[3]) + 1
				gene_lines_3[item[0]].append(new_line)
				n_cds2 += 1
				continue
				
			# adjust last CDS coor for in frame stop or start codon
			else:
				cds_len = int(item[0][3]) - int(line[3]) + 1
				total_len += cds_len
				print(total_len, item[0][3], total_len%3)
				if total_len%3 == 0:
					new_end = int(item[0][3])
				if total_len%3 == 1:
					new_end = int(item[0][3]) + 1
				if total_len%3 == 2:
					new_end = int(item[0][3]) + 1
				print(new_end, int(item[0][3]))
				new_line[4] = str(new_end)
				gene_lines_3[item[0]].append(new_line)
				continue
				
		# keep all other lines
		else:
			if item[0] not in gene_lines_3:
				gene_lines_3[item[0]] = [line]
			else:
				gene_lines_3[item[0]].append(line)
				
print('#########')
				
for item in gene_lines_3.items():
	total = 0
	for line in item[1]:
		if line[2] == 'CDS':
			print(int(line[4]), int(line[3]))
			cds_len = int(line[4]) - int(line[3]) + 1
			total += cds_len
			print(cds_len)
	print(total/3)

'''
if args.out_dir.endswith('/'): 
	out = args.out_dir
else:
	out = f'{args.out_dir}/'
	
if not os.path.exists(out):
	os.makedirs(out)
						
for item in gene_lines_3.items():
	with open(f'{out}{item[0][1]}.gff3', 'wt') as fp:
		for line in item[1]:
			fp.write('\t'.join(line)+'\n')
'''


		
			
