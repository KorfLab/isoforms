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
			
# get intron sequences from annotation
for gene in gene_file_paths.items():
	print(gene[1][0], gene[1][1])

	with open(gene[1][0]) as fp:
		for line in fp:
			line = line.rstrip()
			line = line.split('\t')
			if line[2] == 'CDS':
				print(line)
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
