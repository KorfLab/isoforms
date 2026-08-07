import argparse
import glob

parser = argparse.ArgumentParser(description='get GC content data from '
	'smallgenes and APC')
parser.add_argument('smallgenes')
parser.add_argument('--apc_results', nargs='+')
	
args = parser.parse_args()

if not args.smallgenes.endswith('/'):
	args.smallgenes = f'{args.smallgenes}/'
	
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
		
def gc_content(fasta, gff):
	
	with open(fasta) as fp:
		
		seq = []
		for line in fp:
			line = line.rstrip()
			if line.startswith('>'): continue
			for n in line:
				seq.append(n)
				
	with open(gff) as fp:
	
		total_cds = 0		
		for line in fp:
			line = line.rstrip()
			line = line.split('\t')
			if line[2] == 'CDS': total_cds += 1
			

	with open(gff) as fp:
		
		cds_res = {}
		intron_res = {}
		intron_bounds = []
		cds_count = 1
		for line in fp:
			line = line.rstrip()
			line = line.split('\t')
			if line[2] == 'CDS':
				l_cds = int(line[3]) - 1
				r_cds = int(line[4]) - 1
				cds = (l_cds, r_cds)
				cds_seq = seq[l_cds:r_cds+1]
				cres = gc_calc(cds_seq)
				cds_res[cds] = cres
				
				if cds_count == 1:
					l_int = r_cds + 1  
					intron_bounds.append(l_int)
					
				if cds_count > 1:
					
					if cds_count < total_cds:
						l_int = l_cds - 1
						r_int = r_cds + 1
						intron_bounds.append(l_int)
						intron_bounds.append(r_int)
						
					if cds_count == total_cds:
						r_int = l_cds - 1
						intron_bounds.append(r_int)
					
				cds_count += 1
				
		cds_res = dict(sorted(cds_res.items()))
		intron_bounds = sorted(intron_bounds)
				
		for i in range(0, len(intron_bounds), 2):
			intron = tuple(intron_bounds[i:i+2])
			intron_seq = seq[intron[0]:intron[1]+1]
			ires = gc_calc(intron_seq)
			intron_res[intron] = ires
				
		return {'cds': cds_res, 'intron': intron_res}		
		
def gc_content(fasta, exons):
	
	with open(fasta, 'rt') as fp:
		
		seq = []
		for line in fp:
			line = line.rstrip()
			if line.startswith('>'): continue
			for n in line:
				seq.append(n)
				
	for item in exons.items():
		print(item)
				
		
		

			
			

for gene in gene_file_paths.items():
	
	fasta = gene[1][0]
	wormbase = gene[1][1]
	apc_res = gene[1][2:]
	
	#g = gc_content(fasta, wormbase)
	
	wb_cds = get_wb_cds(wormbase)
	print('#################')
	for r in apc_res:
		apc_cds = get_apc_cds(r)
		break
		
	gc_content(fasta, wb_cds)
	
	#gc_content(fasta, apc_cds)
	
	break

	
# ce.2.273 has an isoform with no introns (one CDS)		
		
# ce.2.243 has 7 cds
# 2 annotated gene structures
# WBGene00045461

# ce.1.350
	


	
	
	

	
	
	
'''
gc_contents = {}
for gene in gene_file_paths.items():
	
	if not gene[0] == 'ce.2.58': continue
	
	print(gene[0])
	print(gene[1][0], gene[1][1])
	

	
	seq = []
	with open(gene[1][1]) as fp:
		
		for line in fp:
			line = line.rstrip()
			if line.startswith('>'): continue
			for n in line:
				seq.append(n)
			
	print(seq)
	print('#######')
			
	with open(gene[1][0]) as fp:
		
		total_cds = 0
		for line in fp:
			line = line.rstrip()
			line = line.split('\t')
			if line[2] == 'CDS': total_cds += 1
				

	cds_res = {}
	intron_res = {}

	with open(gene[1][0]) as fp:
		
		intron_bounds = []
		cds_count = 1
		for line in fp:
			line = line.rstrip()
			line = line.split('\t')
			if line[2] == 'CDS':
				
				l_cds = int(line[3]) - 1
				r_cds = int(line[4]) - 1
				cds = (l_cds, r_cds)
				cds_seq = seq[l_cds:r_cds+1]
				print(cds)
				print(cds_seq)
				cres = gc_calc(cds_seq)
				cds_res[cds] = cres
				
				if cds_count == 1:
					l_int = r_cds + 1  
					intron_bounds.append(l_int)
				if cds_count > 1:
					if cds_count < total_cds:
						l_int = l_cds - 1
						r_int = r_cds + 1
						intron_bounds.append(l_int)
						intron_bounds.append(r_int)
					if cds_count == total_cds:
						r_int = l_cds - 1
						intron_bounds.append(r_int)
					
				cds_count += 1
				

		cds_res = dict(sorted(cds_res.items()))
		
		print(cds_res)
		print(sorted(intron_bounds))
		
		intron_bounds = sorted(intron_bounds)
		
		for item in cds_res.items():
			print(item)
				
		for i in range(0, len(intron_bounds), 2):
			intron = tuple(intron_bounds[i:i+2])
			intron_seq = seq[intron[0]:intron[1]+1]
			print(intron)
			print(intron_seq)
			#ires = gc_calc(intron_seq)
			#intron_res[intron] = ires
	
	#gc_contents[gene[0]] = {'cds': cds_res, 'intron': intron_res}
	
	#print(gc_contents)
			

	#break
'''



	





'''
	
if not args.APC.endswith('/'):
	args.APC = f'{args.APC}/'

gene_paths = {}
for fasta in glob.glob(f'{args.smallgenes}*.fa'):
		gid = '.'.join(fasta.split('/')[-1].split('.')[:-1])
		gff = f'{args.smallgenes}{gid}.gff3'
		gene_paths[gid] = [fasta, gff]
		

for fasta in glob.glob(f'{args.APC}*'):
	gid = '.'.join(fasta.split('/')[-1].split('.')[:-2])
	print(gid)
	print(fasta)
	
'''
