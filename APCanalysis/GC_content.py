import argparse
import glob
import json

parser = argparse.ArgumentParser(description='get GC content data from '
	'smallgenes and APC')
parser.add_argument('smallgenes')
parser.add_argument('--apc_results', nargs='+')
parser.add_argument('--out_file', default='GC_con', 
					help='name of output file [%(default)s]')
	
args = parser.parse_args()

if not args.smallgenes.endswith('/'):
	args.smallgenes = f'{args.smallgenes}/'
	
# gather file paths from each directory
gene_file_paths = {}
for gff_path in glob.glob(f'{args.smallgenes}*gff3'):
	gid = gff_path.split('/')[-1]
	gid = '.'.join(gid.split('.')[:-1])
	fa_path = f'{args.smallgenes}{gid}.fa'
	gene_file_paths[gid] = [fa_path, gff_path]

for apc_dir in args.apc_results:
	
	if not apc_dir.endswith('/'):
		apc_dir = f'{apc_dir}/'

	for fpath in glob.glob(f'{apc_dir}*'):
		gid = fpath.split('/')[-1]
		gid = '.'.join(gid.split('.')[:-2])
		gene_file_paths[gid].append(fpath)
			
def gc_calc(seq):
	
	nts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
	for nt in seq:
		nts[nt] += 1
		
	gc_calc = (
		((nts['G'] + nts['C']) / 
		(nts['A'] + nts['C'] + nts['G'] + nts['T'])) * 100
	)	
	
	return round(gc_calc, 2)
	
# get matching cds and transcript regions from smallgenes gff
def get_wb_cds(gff):
	
	gff_lines = []
	with open(gff, 'rt') as fp:
		
		for line in fp:
			line = line.rstrip()
			line = line.split('\t')  
			gff_lines.append(line)
			
	# get transcript IDs
	mrna_cds = {}
	for line in gff_lines:	
		if line[2] == 'mRNA':
			m_tid = line[8].split(';')[0].split(':')[1]
			m_tid = '.'.join(m_tid.split('.')[:-1])
			mrna_cds[m_tid] = []
	
	# match CDS to transcript IDs
	for line in gff_lines:
		if line[2] == 'CDS':
			c_tid = line[8].split(';')[0].split(':')[1]
			if c_tid in mrna_cds:
				cds_coor = (int(line[3]), int(line[4]))
				mrna_cds[c_tid].append(cds_coor)
			
	return mrna_cds
	
# get matching cds and isoforms from apc gff results
def get_apc_cds(gff):
	
	with open(gff, 'rt') as fp:
		
		# collect cds coors for each isoform
		iso_n = 0
		iso_exons = {}
		for line in fp:
			line = line.rstrip()
			if line.startswith('#'): continue
			if line == '': continue
			line = line.split('\t')
			if line[2] == 'gene': continue
			if line[2] == 'mRNA': 
				iso_n += 1
				iso_exons[iso_n] = []
			if line[2] == 'exon':
				exon_coor = (int(line[3]), int(line[4]))
				iso_exons[iso_n].append(exon_coor)

		return iso_exons
		
def gc_content(fasta, exons):
	
	with open(fasta, 'rt') as fp:
		
		seq = []
		for line in fp:
			line = line.rstrip()
			if line.startswith('>'): continue
			for n in line:
				seq.append(n)

	for iso in exons.items():
		cds_res = {}
		intron_res = {}
		intron_bounds = []
		cds_count = 0
		for cds in iso[1]:
			# swap to 0 based indexing
			# get cds seq/gc
			l_cds = cds[0] - 1
			r_cds = cds[1] - 1
			cds = (l_cds, r_cds)
			cds_seq = seq[l_cds:r_cds+1]
			cres = gc_calc(cds_seq)
			cds_key = f'{l_cds},{r_cds}'
			cds_res[cds_key] = cres
			
			# get left and right intron boundaries
			if cds_count == 0:
				l_int = r_cds + 1
				intron_bounds.append(l_int)
				
			if cds_count > 0:
				
				if cds_count < len(iso[1]) - 1:
					l_int = l_cds - 1
					r_int = r_cds + 1
					intron_bounds.append(l_int)
					intron_bounds.append(r_int)
					
				if cds_count == len(iso[1]) - 1:
					r_int = l_cds - 1
					intron_bounds.append(r_int)
					
			cds_count += 1
			
		# ok to sort after 
		cds_res = dict(sorted(cds_res.items()))
		intron_bounds = sorted(intron_bounds)

		# some isoforms have no introns
		if len(intron_bounds) > 1:
			# get intron seq/gc
			for i in range(0, len(intron_bounds), 2):
				intron = tuple(intron_bounds[i:i+2])
				intron_seq = seq[intron[0]:intron[1]+1]
				ires = gc_calc(intron_seq)
				intron_key = f'{intron[0]},{intron[1]}'
				intron_res[intron_key] = ires
		else:
			intron_res = None
					
		yield {'cds': cds_res, 'intron': intron_res}
	
data = {}
for gene in gene_file_paths.items():
	fasta = gene[1][0]
	wormbase = gene[1][1]
	apc_files = gene[1][2:]
	
	wb_cds = get_wb_cds(wormbase)
	wb_gc = gc_content(fasta, wb_cds)
	
	wb_data = {}
	iso_n = 0
	for gc in wb_gc:
		if gene[1][1] not in wb_data:
			wb_data[gene[1][1]] = {iso_n: gc}
		else:
			wb_data[gene[1][1]][iso_n] = gc		
		iso_n += 1

	apc_data = {}
	for apcf in apc_files:
		apc_cds = get_apc_cds(apcf)
		apc_gc = gc_content(fasta, apc_cds)
		iso_n = 0
		for gc in apc_gc:
			if apcf not in apc_data:
				apc_data[apcf] = {iso_n: gc}
			else:
				apc_data[apcf][iso_n] = gc
			iso_n += 1
	
	data[gene[0]] = {'wb_data': wb_data, 'apc_data': apc_data}

with open(f'{args.out_file}.json', 'w') as jfile:
	json.dump(data, jfile, indent=4)



########## code testing ##########

# 1.370 WBGene00008792	isoform with no introns
#test_f = '../data/smallgenes/ce.1.370.fa', '../data/smallgenes/ce.1.370.gff3'

#cds = get_wb_cds(test_f[1])
#gc = gc_content(test_f[0], cds)

#for c in cds.items():
#	print(c)
	
#for g in gc:
#	print(g)


# ce.2.273 has an isoform with no introns (one CDS)		
		
# ce.2.243 has 7 cds
# 2 annotated gene structures
# WBGene00045461

# ce.1.350		


'''
for gene in gene_file_paths.items():
	
	if gene[0] == 'ce.1.542':
	
		fasta = gene[1][0]
		wormbase = gene[1][1]
		apc_res = gene[1][2:]
	
		wb_cds = get_wb_cds(wormbase)
		gc = gc_content(fasta, wb_cds)
		for g in gc:
			print(g)
	
		print('##########')
	
		for ar in apc_res:
			apc_cds = get_apc_cds(ar)
			gc = gc_content(fasta, apc_cds)
			g_count = 0
			for g in gc:
				if g_count < 5:
					print(g)
				g_count += 1
'''




	


