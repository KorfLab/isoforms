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
	gene_file_paths[gid] = [gff_path, fa_path]

for apc_dir in args.apc_results:
	
	if not apc_dir.endswith('/'):
		apc_dir = f'{apc_dir}/'

	for fpath in glob.glob(f'{apc_dir}*'):
		gid = fpath.split('/')[-1]
		gid = '.'.join(gid.split('.')[:-2])
		gene_file_paths[gid].append(fpath)
			
def GC(seq):
	
	nts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
	for nt in seq:
		nts[nt] += 1
		
	gc_con = (
		((nts['G'] + nts['C']) / 
		(nts['A'] + nts['C'] + nts['G'] + nts['T'])) * 100
	)	
	
	return round(gc_con, 2)
	
# get cds and intron sequences from annotation
gc_contents = {}
for gene in gene_file_paths.items():
	print(gene[0])
	print(gene[1][0], gene[1][1])
	
	seq = []
	with open(gene[1][1]) as fp:
		
		for line in fp:
			line = line.rstrip()
			if line.startswith('>'): continue
			for n in line:
				seq.append(n)
			
	with open(gene[1][0]) as fp:
		
		total_cds = 0
		for line in fp:
			line = line.rstrip()
			line = line.split('\t')
			if line[2] == 'CDS': total_cds += 1

	with open(gene[1][0]) as fp:
		
		intron_bounds = []
		cds_count = 1
		for line in fp:
			line = line.rstrip()
			line = line.split('\t')
			if line[2] == 'CDS':
				
				l_cds = int(line[3]) - 1
				r_cds = int(line[4]) - 1
				cds = [l_cds, r_cds]
				cds_seq = seq[l_cds:r_cds+1]
				res = GC(cds_seq)
				print(cds, res)
				
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
				
		print('###############')
		for i in range(0, len(intron_bounds), 2):
			intron = intron_bounds[i:i+2]
			intron_seq = seq[intron[0]:intron[1]+1]
			res = GC(intron_seq)
			print(intron, res)
			

	break
	



	





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
