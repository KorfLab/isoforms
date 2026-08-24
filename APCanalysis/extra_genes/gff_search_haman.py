import argparse
import gzip
import glob
import os
import shutil

parser = argparse.ArgumentParser(description='find genes of interest '
			'from haman genes')
parser.add_argument('WBGenes', help='text with list of WBGene IDs and '
			'region of interest i.e. '
			'WBGene, gene name, chr, left bound, right bound, '
			'first exon start, sense '
			'WBGene00003386,mod1,V,8910090,8910840,8913992,-')
parser.add_argument('gene_dir', help='directory with gff and fa files')
parser.add_argument('out_dir', help='name of directory to store '
			'fa and gff files')

args = parser.parse_args()

ginfo = {}
with open(args.WBGenes, 'rt') as fp:
	for line in fp:
		# to test a single WBGene
		if line.startswith('#'): continue
		line = line.rstrip()
		line = line.split(',')
		ginfo[line[0]] = line[1:]

# some haman genes have more than one WBGene in a gff

gene_chs = {}
for gff_file in glob.glob(args.gene_dir + '*.gff3'):
	with open(gff_file, 'rt') as fp:
		for line in fp:
			line = line.rstrip().split('\t')
			if line[1] == 'WormBase' and line[2] == 'gene':
				print(line[0], line[8])


# get any fasta/gff files that mention a WBGene
for gff_file in glob.glob(args.gene_dir + '*.gff3'):
	with open(gff_file, 'rt') as fp:
		for line in fp:
			line = line.rstrip().split('\t')
			if line[1] == 'WormBase' and line[2] == 'gene':
				wbg = line[8].split(':')[1]
				if wbg in ginfo:
					ginfo[wbg].append(gff_file)
					base = gff_file.split('/')[-1].split('.')[:-1]
					ginfo[wbg].append(f'{args.gene_dir}'
										f'{'.'.join(base)}.fa')							

# i don't think this method will work with optiso
# need to rework make_svg instead
'''
coors1 = [15, 25, 30, '-']
coors2 = [15, 25, 10, '+']
flank = 5
if coors1[3] == '-':
	coors1 = [abs(coors1[1]-coors1[2])+flank+1, abs(coors1[0]-coors1[2])+flank+1]
if coors2[3] == '+':
	coors2 = [coors2[0]-coors2[2]+flank+1, coors2[1]-coors2[2]+flank+1]
'''

lines_to_write = {}
# hard coded, use same flank as haman genes
flank = 100
for item in ginfo.items():
	if len(item[1]) < 7: continue
	# adjust genomic coordinates to gene 
	if item[1][5] == '-':
		start = abs(int(item[1][3]) - int(item[1][4])) + flank+1
		end = abs(int(item[1][2]) - int(item[1][4])) + flank+1
		gen_coors = [start, end]
	else:
		start = abs(int(item[1][2]) - int(item[1][4])) + flank+1
		end = abs(int(item[1][3]) - int(item[1][4])) + flank+1
		gen_coors = [start, end]
		
	keep_lines = []	
	with open(item[1][6], 'rt') as fp:
		for line in fp:
			line = line.rstrip()
			line = line.split('\t')
			if len(line) >= 8:
				if (line[1] == 'WormBase' and line[2] 
					in ['exon', 'gene', 'mRNA', 'CDS', 'intron']):
						#print(line)
					
						
						# get any regions that overlap
						f_start, f_end = int(line[3]), int(line[4])
						g_start, g_end = gen_coors
						if f_start <= gen_coors[1] and f_end >= gen_coors[0]:
							
								# negative CDS coors needed for make_svg.py
								line[3] = str(int(line[3]) - gen_coors[0] +1)
								line[4] = str(int(line[4]) - gen_coors[0] +1)
								keep_lines.append(line)
				
				if line[1] == 'RNASeq_splice':
					# only get introns within coors
					if (int(line[3]) <= gen_coors[1] and 
						int(line[3]) >= gen_coors[0] and
						int(line[4]) <= gen_coors[1] and
						int(line[4]) >= gen_coors[0]):
							line[3] = str(int(line[3]) - gen_coors[0] +1)
							line[4] = str(int(line[4]) - gen_coors[0] +1)
							keep_lines.append(line)

	lines_to_write[item[1][0]] = keep_lines
	
if not os.path.isdir(f'{args.out_dir}'):
	os.mkdir(f'{args.new_dir}/')
				
for item in lines_to_write.items():
	with open(f'gpt{item[0]}.gff3', 'wt') as fp:
		for line in item[1]:
			fp.write(f'{'\t'.join(line)}\n')


'''
for item in gfiles.items():
	if len(item[1]) == 1: continue
	src1 = item[1][1]
	src2 = item[1][2]
	dest = args.new_dir
	shutil.copy(src1, dest)
	shutil.copy(src2, dest)
'''

			
			
			
			
			
			
			
			
			
			
