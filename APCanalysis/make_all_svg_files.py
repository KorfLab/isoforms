import argparse
import glob
import os
import subprocess

parser = argparse.ArgumentParser(
	description='make svg files for all genes using make_svg.py')
parser.add_argument('smallgenes', help='directory with smallgenes '
					'fasta and gff files')
parser.add_argument('apc_results', help='directory with APC gff file '
					'results')
parser.add_argument('--out_dir', help='name of output directory') 
parser.add_argument('--program', required=False, 
					default='make_svg.py', help='path to make_svg.py '
					'[%(default)s]')
parser.add_argument('--gc', type=str, required=False, nargs=2, 
					help='directories of json files with GC content '
					'data; must be in order, with the '
					'smallgenes/rnaseq data first and apc data second. '
					'Example: --gc gcc_res/rnaseq_gcc/ '
					'gcc_res/APC_base_gcc/')

args = parser.parse_args()

if not args.smallgenes.endswith('/'):
	args.smallgenes = f'{args.smallgenes}/'
	
if not args.apc_results.endswith('/'):
	args.apc_results = f'{args.apc_results}/'
	
if not args.out_dir.endswith('/'):
	args.out_dir = f'{args.out_dir}/'
	
if not os.path.isdir(args.out_dir):
	os.mkdir(args.out_dir)
	
if args.gc:
	if not args.gc[0].endswith('/'):
		args.gc[0] = f'{args.gc[0]}/'	
	if not args.gc[1].endswith('/'):
		args.gc[1] = f'{args.gc[1]}/'		
	
# add file paths in order
gene_files = {}
for apc_path in glob.glob(f'{args.apc_results}*'):
	gid = '.'.join(apc_path.split('/')[-1].split('.')[:-2])
	gene_files[gid] = [apc_path]

for rna_path in glob.glob(f'{args.smallgenes}*.gff3'):
	gid = '.'.join(rna_path.split('/')[-1].split('.')[:-1])
	gene_files[gid].append(rna_path)
	
#for rna_gc_path in glob.glob(f'{args.gc[0]}')

#for apc_gc_path in glob.glob(f'{args.gc[0]}')

if args.gc:
	for gc_path in zip(glob.glob(f'{args.gc[0]}*'), 
						glob.glob(f'{args.gc[1]}*')):
		gid = '.'.join(gc_path[0].split('/')[-1].split('.')[:-2])
		gene_files[gid].append(gc_path[0])
		gene_files[gid].append(gc_path[1])
	
def build_cmd(prog, apc_gff, rna_gff, out_name, gc=None):
	
	if gc:
		cmd = (
			f'python3 {prog} {apc_gff} {rna_gff} --out_name {out_name} '
			f'--gc {gc[0]} {gc[1]}'
		)		
	else:
		cmd = (
			f'python3 {prog} {apc_gff} {rna_gff} --out_name {out_name}'
		)
		
	return cmd
	
for gpaths in gene_files.items():
	print(gpaths[0])
	print(len(gpaths[1]))
	if len(gpaths[1]) == 4:
		gc_paths = [gpaths[1][2], gpaths[1][3]]
	else:
		gc_paths = None
	cmd = build_cmd(
			args.program, gpaths[1][0], gpaths[1][1], 
			f'{args.out_dir}{gpaths[0]}.svg', gc=gc_paths
		)
	print(cmd)
	break




#subprocess.run('ls')
'''
python3 make_svg.py ../data/APCisos_optiso/ce.2.243.APC_optiso.gff 
../data/smallgenes/ce.2.243.gff3 --out_name ce.2.243.svg 
--gc gcc_res/rnaseq_gcc/ce.2.243.gcc.json 
gcc_res/APC_optiso_gcc/ce.2.243.gcc.json
'''

